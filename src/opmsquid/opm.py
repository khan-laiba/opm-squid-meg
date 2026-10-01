"""Explicit finite OPM arrays (G2, G3).

Sensor model (study assumptions, labelled in docs/provenance_register.md)
------------------------------------------------------------------------
* Sensing volume: a 10-mm cubic vapour cell (Jas et al. 2026, experimental system, as stated in
  GOAL.md), represented by a custom MNE coil definition ``OPM_COIL_TYPE`` with 1, 8 (2x2x2
  Gauss-Legendre) or 27 (3x3x3) integration points averaging the field component along the
  sensitive axis. This is a single-axis, point-sampling-free magnetometer model; the cell size
  is a declared assumption, not a vendor specification.
* Placement: the sensing centre lies ``standoff`` = 7 mm from the helmet's inner surface
  (Jas et al. 2026 experimental geometry, as stated in GOAL.md) plus a separate, declared
  ``scalp_gap`` between the helmet's inner surface and the MRI scalp (default 0 mm; swept).
* Sites are defined on the MRI scalp; the sensitive axis is the normal of a smooth reference
  surface averaged within ``normal_radius`` (G2: the BEM head surface within 15 mm; without
  ``axis_surface``, the MRI scalp itself). Clearance, coverage and packing rules are documented
  in the builders below (A-OPM-CLEAR, A-OPM-COVER, A-OPM-PACK); the clearance rules check the
  sensing centre and every integration point of the cell, as oriented in the forward model, with
  exact point-to-triangle distances to the BEM head surface.

In MNE 1.13.2, ``mne.use_coil_def(fname)`` *adds* the coils defined in ``fname`` to the standard
``coil_def.dat`` (``_read_coil_defs``), and its parser rejects blank lines. ``coil_def_file``
therefore writes only the study OPM coil; wrap forward computations for OPM arrays in
``with mne.use_coil_def(coil_def_file()):``.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import mne
from mne.io.constants import FIFF
from scipy.spatial import cKDTree

from . import paths

OPM_COIL_TYPE = 9901  # study-specific id (unused by MNE 1.13.2's coil_def.dat)
CELL_SIZE = 0.010  # m, cubic vapour cell (assumption A-OPM-CELL)
STANDOFF = 0.007  # m, sensing centre to helmet inner surface (assumption A-OPM-STANDOFF)


def _gauss_cube(n: int, size: float):
    x, w = np.polynomial.legendre.leggauss(n)
    x = x * size / 2.0
    w = w / w.sum()
    g = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3)
    wg = np.einsum("i,j,k->ijk", w, w, w).ravel()
    return wg, g


def _coil_block(coil_type: int, size: float, description: str) -> str:
    lines = []
    for accuracy, n in ((0, 1), (1, 2), (2, 3)):
        if n == 1:
            w, g = np.array([1.0]), np.zeros((1, 3))
        else:
            w, g = _gauss_cube(n, size)
        lines.append(f'1   {coil_type}    {accuracy}  {len(w):>2d}  {size:.3e}  0.000e+00\t"{description}"')
        lines += [f"  {wi:.12e}  {x:.12e}  {y:.12e}  {z:.12e}  0.000  0.000  1.000" for wi, (x, y, z) in zip(w, g)]
    return "\n".join(lines) + "\n"


FOUR_POINT_COILS = {3014: 9014, 3024: 9024}  # T3 coils with MNE's 4-point 'normal' rule at every level


def _four_point_blocks() -> str:
    """Study coils 9014/9024: MNE's own 'normal' (4-point) integration of the T3 planar
    gradiometer 3014 and magnetometer 3024, repeated at every accuracy level so that MNE
    1.13.2's forward (which always requests 'accurate') uses the 4-point rule stated by
    Hunold et al. (2016)."""
    from mne.forward._make_forward import _read_coil_defs

    defs = {(int(c["coil_type"]), int(c["accuracy"])): c for c in _read_coil_defs()}
    out = []
    for src_type, new_type in FOUR_POINT_COILS.items():
        c = defs[(src_type, 1)]
        cls = 3 if src_type == 3014 else 1
        for accuracy in (0, 1, 2):
            out.append(f'{cls}   {new_type}    {accuracy}  {len(c["w"]):>2d}  {c["size"]:.3e}  {c["base"]:.3e}\t'
                       f'"{c["desc"]} [4-point normal rule at all levels]"')
            out += [f"  {w:.12e}  {r[0]:.12e}  {r[1]:.12e}  {r[2]:.12e}  {n[0]:.6f}  {n[1]:.6f}  {n[2]:.6f}"
                    for w, r, n in zip(c["w"], c["rmag"], c["cosmag"])]
    return "\n".join(out) + "\n"


def coil_def_file(cell_size: float = CELL_SIZE) -> Path:
    """Path of the study's extra coil-definition file (for use_coil_def): the OPM cell coil
    and the 4-point Neuromag emulation coils."""
    extra = _coil_block(OPM_COIL_TYPE, cell_size, f"Study OPM, {cell_size * 1e3:.1f}-mm cubic cell (assumption)")
    text = ("# opmsquid study coils (src/opmsquid/opm.py); added to MNE's coil_def.dat\n" + extra
            + _four_point_blocks())
    digest = hashlib.sha1(text.encode()).hexdigest()[:10]
    out = paths.CACHE / "coils" / f"coil_def_{digest}.dat"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    return out


@dataclass
class OPMArray:
    """Sensing centres and sensitive axes in the head frame, and bookkeeping."""

    pos: np.ndarray  # (n, 3) head frame [m]
    axis: np.ndarray  # (n, 3) unit sensitive axes
    scalp_point: np.ndarray  # (n, 3) head frame, where the site meets the scalp
    site_id: np.ndarray  # (n,) index of the originating site (e.g. Neuromag location)
    standoff: float
    scalp_gap: float
    label: str

    @property
    def n_sites(self) -> int:
        return len(np.unique(self.site_id))


def smoothed_normals(rr: np.ndarray, nn: np.ndarray, points: np.ndarray, radius: float = 0.010) -> np.ndarray:
    """Average vertex normal of a surface within ``radius`` of each point (helmet-scale)."""
    tree = cKDTree(rr)
    out = np.empty((len(points), 3))
    for i, idx in enumerate(tree.query_ball_point(points, radius)):
        v = nn[idx].sum(axis=0) if len(idx) else nn[tree.query(points[i])[1]]
        out[i] = v / np.linalg.norm(v)
    return out


def place_on_scalp(scalp_points: np.ndarray, scalp_normals: np.ndarray, site_id, standoff: float = STANDOFF,
                   scalp_gap: float = 0.0, label: str = "opm") -> OPMArray:
    """Sensing centres at scalp + (scalp_gap + standoff) along the (outward) scalp normal."""
    pos = scalp_points + (scalp_gap + standoff) * scalp_normals
    return OPMArray(pos, scalp_normals.copy(), scalp_points.copy(), np.asarray(site_id), standoff, scalp_gap, label)


def cell_frame(n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """In-plane axes (ex, ey) of the cell whose sensitive axis is ``n`` (as written to ch['loc'])."""
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    ex = np.cross(a, n)
    ex /= np.linalg.norm(ex)
    return ex, np.cross(n, ex)


def cell_points(pos: np.ndarray, n: np.ndarray, cell_size: float = CELL_SIZE, to_head: np.ndarray | None = None) -> np.ndarray:
    """The 27 integration points (MNE 'accurate' level) of a cell centred at ``pos`` with axis ``n``.
    The cell's in-plane orientation is the one ``make_info`` writes, which ``cell_frame`` derives
    from the head-frame axis: for ``pos``/``n`` in another frame (e.g. MRI), pass ``to_head``, the
    rotation taking that frame's vectors to the head frame."""
    if to_head is None:
        ex, ey = cell_frame(n)
    else:
        ex_h, ey_h = cell_frame(to_head @ n)
        ex, ey = to_head.T @ ex_h, to_head.T @ ey_h
    g = _gauss_cube(3, cell_size)[1]
    return pos + g[:, :1] * ex + g[:, 1:2] * ey + g[:, 2:] * n


def make_info(array: OPMArray, sfreq: float = 1000.0) -> mne.Info:
    """MNE info for a single-axis OPM array; device frame = head frame (head-mounted array)."""
    names = [f"OPM{i:03d}" for i in range(len(array.pos))]
    info = mne.create_info(names, sfreq, "mag")
    with info._unlock():
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")
        for ch, p, n in zip(info["chs"], array.pos, array.axis):
            ex, ey = cell_frame(n)
            ch["loc"][:3], ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = p, ex, ey, n
            ch["coil_type"] = OPM_COIL_TYPE
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
            ch["unit"] = FIFF.FIFF_UNIT_T
    return info


def signed_distance(points: np.ndarray, surf: dict) -> np.ndarray:
    """Distance [m] from each point to the closed MNE surface ``surf`` (same frame), negative
    inside (nearest vertex of the surface subdivided three times, ~1-mm vertex spacing)."""
    from mne.surface import _CheckInside
    from .anatomy import Surface

    fine = Surface(surf["rr"], surf["tris"], surf["nn"]).subdivided(3)
    d = cKDTree(fine.rr).query(points)[0]
    return np.where(_CheckInside(surf)(points), -d, d)


def farthest_point_subset(pos: np.ndarray, n: int, start: int | None = None) -> np.ndarray:
    """Indices of ``n`` points chosen by farthest-point sampling (Euclidean), starting from
    ``start`` (default: the highest point, max z). Spreads a subset evenly over an array."""
    k = int(np.argmax(pos[:, 2])) if start is None else int(start)
    chosen = [k]
    dmin = np.linalg.norm(pos - pos[k], axis=1)
    while len(chosen) < n:
        k = int(np.argmax(dmin))
        chosen.append(k)
        dmin = np.minimum(dmin, np.linalg.norm(pos - pos[k], axis=1))
    return np.sort(np.array(chosen))


def subset(arr: "OPMArray", idx: np.ndarray, label: str) -> "OPMArray":
    return OPMArray(pos=arr.pos[idx], axis=arr.axis[idx], scalp_point=arr.scalp_point[idx], site_id=arr.site_id[idx],
                    standoff=arr.standoff, scalp_gap=arr.scalp_gap, label=label)


def min_spacing(pos: np.ndarray) -> np.ndarray:
    """Distance from each sensing centre to its nearest neighbour [m]."""
    d, _ = cKDTree(pos).query(pos, k=2)
    return d[:, 1]


def fiducials_head(info: mne.Info) -> dict:
    """LPA, nasion and RPA [m] in the head frame from the digitisation in ``info``."""
    idents = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}
    out = {idents[d["ident"]]: np.asarray(d["r"], float) for d in info["dig"] or []
           if d["kind"] == FIFF.FIFFV_POINT_CARDINAL and d["ident"] in idents}
    if len(out) != 3:
        raise ValueError("info has no digitised LPA/nasion/RPA")
    return out


def above_brow_plane(points_head: np.ndarray, fids: dict, brow_offset: float = 0.030) -> np.ndarray:
    """True for points above the plane through LPA, RPA and a point ``brow_offset`` above the
    nasion (head frame, +z up): the scalp region an OPM helmet can cover (no face, no neck)."""
    brow = fids["nasion"] + np.array([0.0, 0.0, brow_offset])
    normal = np.cross(fids["rpa"] - fids["lpa"], brow - fids["lpa"])
    normal *= np.sign(normal[2]) / np.linalg.norm(normal)
    return (points_head - fids["lpa"]) @ normal > 0


def resolve_clearance(pos: np.ndarray, axis: np.ndarray, scalp, min_clearance: float, step: float = 0.0005,
                      max_extra: float = 0.015, model_surface: dict | None = None,
                      model_clearance: float = 0.0, cell_clearance: float | None = None,
                      to_head: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Move sensing centres outward along their axes until every point is at least
    ``min_clearance`` from the scalp (e.g. over the ear pinna or brow ridge, where the offset
    along a smoothed normal comes close to other parts of the head) and, with ``model_surface``
    (the BEM head surface, same frame), the centre is at least ``model_clearance`` outside it and,
    with ``cell_clearance``, every integration point of the cell at least that far outside it.
    Distances to the model surface are exact point-to-triangle distances with the inside test;
    the cell has the orientation the forward model uses (``to_head``: rotation from this frame to
    the head frame, see ``cell_points``). Returns (positions, extra outward shift per site,
    feasible mask: False if more than ``max_extra`` would be needed)."""
    from .anatomy import MeshDistance

    tree = cKDTree(scalp.rr)
    model = MeshDistance(model_surface) if model_surface is not None else None

    def too_close(p, n):
        if tree.query(p)[0] < min_clearance:
            return True
        if model is not None and model.signed(p)[0] < model_clearance:
            return True
        return (model is not None and cell_clearance is not None
                and model.signed(cell_points(p, n, to_head=to_head)).min() < cell_clearance)

    pos = pos.copy()
    extra = np.zeros(len(pos))
    ok = np.ones(len(pos), bool)
    for i in range(len(pos)):
        while too_close(pos[i], axis[i]):
            if extra[i] >= max_extra:
                ok[i] = False
                break
            pos[i] += step * axis[i]
            extra[i] += step
    return pos, extra, ok


MIN_CENTER_SPACING = 0.017  # m, packing assumption A-OPM-PACK (10-mm cell in a ~12-17 mm package)
MODEL_CLEARANCE = 0.004  # m, A-OPM-CLEAR: sensing centre to the BEM head surface (the forward model's outer boundary)
CELL_CLEARANCE = 0.001  # m, A-OPM-CLEAR: every integration point of the cell outside the BEM head surface by this much
MAX_EXTRA_SHIFT = 0.005  # m, A-OPM-CLEAR: a package may sit up to 5 mm beyond its nominal standoff to clear the
# scalp (pinna, brow, occipital curvature); a site needing more (e.g. in the occipito-cervical crease) is infeasible


def prune_to_spacing(centres: np.ndarray, min_dist: float, candidates: np.ndarray | None = None) -> np.ndarray:
    """Greedy pruning so that no two kept sensing centres are closer than ``min_dist``: repeatedly
    drop the site with the most violations (ties: the later one). Returns a keep mask."""
    keep = np.ones(len(centres), bool) if candidates is None else candidates.copy()
    while True:
        idx = np.flatnonzero(keep)
        if len(idx) < 2:
            return keep
        d = np.linalg.norm(centres[idx][:, None] - centres[idx][None], axis=-1) + np.eye(len(idx)) * 1e9
        viol = (d < min_dist).sum(axis=1)
        if viol.max() == 0:
            return keep
        worst = np.flatnonzero(viol == viol.max())[-1]
        keep[idx[worst]] = False


def dense_array(scalp, trans: mne.transforms.Transform, digitisation: mne.Info, min_spacing_m: float,
                max_sites: int | None = None, standoff: float = STANDOFF, scalp_gap: float = 0.0,
                normal_radius: float = 0.010, brow_offset: float = 0.030, seed_point=None,
                min_center_spacing: float = MIN_CENTER_SPACING, outer_skin: dict | None = None,
                max_depth_below_skin: float = 0.002, max_height_above_skin: float = 0.004, axis_surface=None,
                ear_clearance: float = 0.0, fov_margin: float = 0.020,
                max_extra_shift: float = MAX_EXTRA_SHIFT) -> tuple[OPMArray, dict]:
    """Scalp-normal single-axis OPM array filling the coverage region (above the brow plane) by
    farthest-point sampling of the MRI scalp: sites are added in order of largest distance to the
    already placed ones until that distance falls below ``min_spacing_m`` (or ``max_sites`` is
    reached). Sensing centres at scalp_gap + standoff along the smoothed normal, with clearance
    resolution as for the matched array, then pruned so that sensing centres are at least
    ``min_center_spacing`` apart (physical packing, A-OPM-PACK).

    With ``outer_skin`` (the MNE BEM head surface, MRI frame), scalp points lying more than
    ``max_depth_below_skin`` inside that smooth surface (ear canals, pinna folds: not reachable
    by a sensor package) or more than ``max_height_above_skin`` outside it (the ear pinna, and the
    flat cap where the MRI head surface is cut at the edge of the field of view, which lies above
    the brow plane at the back of a pitched head) are not used as sites, nor points within
    ``ear_clearance`` of the preauricular points, nor points within ``fov_margin`` (MRI z) of the
    lowest point of the MRI head surface: next to the field-of-view cut the smoothed normals and
    the clearance check see the artificial cap instead of anatomy (A-OPM-COVER). With ``axis_surface`` (a smooth closed surface, e.g. the
    BEM head surface with outward normals) the sensitive axes are its normals averaged within
    ``normal_radius`` (A-OPM-AXIS), instead of the dense scalp's."""
    mri_head = np.linalg.inv(trans["trans"])
    rr_head = scalp.rr @ mri_head[:3, :3].T + mri_head[:3, 3]
    fids = fiducials_head(digitisation)
    cover = above_brow_plane(rr_head, fids, brow_offset) & (scalp.rr[:, 2] >= scalp.rr[:, 2].min() + fov_margin)
    if ear_clearance > 0:
        cover &= np.minimum(np.linalg.norm(rr_head - fids["lpa"], axis=1), np.linalg.norm(rr_head - fids["rpa"], axis=1)) >= ear_clearance
    if outer_skin is not None:
        h = signed_distance(scalp.rr, outer_skin)
        cover &= (h >= -max_depth_below_skin) & (h <= max_height_above_skin)
    cand = np.flatnonzero(cover)
    pts = scalp.rr[cand]
    start = int(np.argmax(pts[:, 2])) if seed_point is None else int(cKDTree(pts).query(seed_point)[1])
    chosen = [start]
    dmin = np.linalg.norm(pts - pts[start], axis=1)
    while True:
        k = int(np.argmax(dmin))
        if dmin[k] < min_spacing_m or (max_sites is not None and len(chosen) >= max_sites):
            break
        chosen.append(k)
        dmin = np.minimum(dmin, np.linalg.norm(pts - pts[k], axis=1))
    sp_mri = pts[chosen]
    ref = axis_surface if axis_surface is not None else scalp
    nrm = smoothed_normals(ref.rr, ref.nn, sp_mri, normal_radius)
    centres = sp_mri + (scalp_gap + standoff) * nrm
    centres, extra, feasible = resolve_clearance(centres, nrm, scalp, standoff - 0.001, max_extra=max_extra_shift,
                                                 model_surface=outer_skin, model_clearance=MODEL_CLEARANCE,
                                                 cell_clearance=CELL_CLEARANCE, to_head=mri_head[:3, :3])
    keep = feasible & prune_to_spacing(centres, min_center_spacing, feasible)
    arr = OPMArray(pos=centres[keep] @ mri_head[:3, :3].T + mri_head[:3, 3], axis=nrm[keep] @ mri_head[:3, :3].T,
                   scalp_point=sp_mri[keep] @ mri_head[:3, :3].T + mri_head[:3, 3], site_id=np.arange(int(keep.sum())),
                   standoff=standoff, scalp_gap=scalp_gap, label=f"dense OPM ({keep.sum()} sites, >= {min_spacing_m * 1e3:g} mm)")
    report = dict(n_sites=int(keep.sum()), min_spacing_mm=float(min_spacing(arr.pos).min() * 1e3),
                  median_spacing_mm=float(np.median(min_spacing(arr.pos)) * 1e3), n_moved_out=int(np.sum(extra[keep] > 0)),
                  max_extra_shift_mm=round(float(extra[keep].max() * 1e3), 3) if keep.any() else 0.0)
    return arr, report


def matched_to_neuromag(squid_info: mne.Info, trans: mne.transforms.Transform, scalp, digitisation: mne.Info,
                        standoff: float = STANDOFF, scalp_gap: float = 0.0, normal_radius: float = 0.010,
                        brow_offset: float = 0.030, min_clearance: float | None = None, axis_surface=None,
                        ear_clearance: float = 0.0, max_extra_shift: float = MAX_EXTRA_SHIFT,
                        outer_skin: dict | None = None) -> tuple[OPMArray, dict]:
    """Matched-site OPM array: every Neuromag sensor location (magnetometer coil centre and
    normal, MRI frame) is projected along its inward normal onto the scalp. A site is kept if
    that scalp point lies above the brow plane; its sensing centre is placed at
    scalp_gap + standoff along the smoothed scalp normal and, where needed, moved further out
    until it is ``min_clearance`` (default standoff - 1 mm) from every scalp point and, with
    ``outer_skin`` (BEM head surface, MRI frame), ``MODEL_CLEARANCE`` from that surface; a site
    needing more than ``max_extra_shift`` is dropped. Returns the array (head frame) and a report
    with the per-site extra shift."""
    from . import neuromag  # local import to avoid a cycle at module import

    geo = neuromag.sensor_geometry(squid_info, frame="mri", trans=trans)
    mags = np.flatnonzero(geo.kind == "mag")
    dist = neuromag.ray_mesh_distance(geo.pos[mags], -geo.normal[mags], scalp.rr, scalp.tris)
    hit = np.isfinite(dist)
    scalp_pts = geo.pos[mags] - np.nan_to_num(dist)[:, None] * geo.normal[mags]
    ref = axis_surface if axis_surface is not None else scalp
    nrm = smoothed_normals(ref.rr, ref.nn, scalp_pts, normal_radius)
    mri_head = np.linalg.inv(trans["trans"])
    pts_head = scalp_pts @ mri_head[:3, :3].T + mri_head[:3, 3]
    fids = fiducials_head(digitisation)
    keep = hit & above_brow_plane(pts_head, fids, brow_offset)
    if ear_clearance > 0:
        keep &= np.minimum(np.linalg.norm(pts_head - fids["lpa"], axis=1), np.linalg.norm(pts_head - fids["rpa"], axis=1)) >= ear_clearance
    centres = scalp_pts + (scalp_gap + standoff) * nrm
    min_clearance = standoff - 0.001 if min_clearance is None else min_clearance
    centres, extra, feasible = resolve_clearance(centres, nrm, scalp, min_clearance, max_extra=max_extra_shift,
                                                 model_surface=outer_skin, model_clearance=MODEL_CLEARANCE,
                                                 cell_clearance=CELL_CLEARANCE, to_head=mri_head[:3, :3])
    keep &= feasible
    arr = OPMArray(pos=centres[keep] @ mri_head[:3, :3].T + mri_head[:3, 3], axis=nrm[keep] @ mri_head[:3, :3].T,
                   scalp_point=pts_head[keep], site_id=geo.site[mags][keep], standoff=standoff, scalp_gap=scalp_gap,
                   label=f"OPM matched to Neuromag ({keep.sum()} sites)")
    report = dict(n_neuromag_sites=len(mags), n_ray_hits=int(hit.sum()), n_kept=int(keep.sum()),
                  squid_coil_to_scalp_mm=dist * 1e3, excluded_sites=geo.site[mags][~keep].tolist(),
                  extra_shift_mm={int(s): round(float(e) * 1e3, 1) for s, e, k in zip(geo.site[mags], extra, keep)
                                  if k and e > 0},
                  min_spacing_mm=min_spacing(arr.pos) * 1e3)
    return arr, report
