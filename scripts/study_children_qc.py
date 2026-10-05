#!/usr/bin/env python3
"""Quality control of the school-aged children's anatomy against their MRIs (G3B, G4; referee round 1).

Requests (round-1 referee reports):
* Fable, single biggest weakness (i): the children's "white surfaces lie 3.7-5.7 mm from the 'scalp'
  (impossible for an 8-year-old, whose scalp plus skull is >= 8 mm) ... child B has 247 targets
  shallower than 10 mm"; to fix: "verify the children's scalp/white-surface geometry against the MRIs
  (scalp-to-cortex distance maps; compare to the dataset's T1) and either repair or drop them".
* Codex, major issue 2: "Provide MRI-based overlays, coordinate-transform checks and quantitative
  surface/registration quality measures. Correct and rerun affected analyses, or remove those
  anatomies from substantive conclusions."

Anatomies, with the surfaces exactly as G3B and G4 load them (``opmsquid.anatomy``): children A-C
(sub-Z213, sub-Z209, sub-Z226; configs/g3b_pediatric.toml) and, as references, the adult (MNE sample;
its dense scalp was made by the same FreeSurfer tool as the children's, mkheadsurf) and the 2-year
template (the closest in age; ``--anatomies`` selects others). Checks:

1. Coordinate frames: the volume geometry (vox2ras: dimensions, voxel size, direction cosines, c_ras)
   stored in every FreeSurfer surface used (white, dense MRI scalp, watershed BEM) against each other
   and against the T1 header (largest displacement over the volume's corners); for the children, the
   prepared MNE files against the FreeSurfer surfaces they were written from. Image-based: the T1 edge
   across the white surface (steepest intensity fall along its outward normal within +-2.5 mm),
   decomposed by least squares into a uniform offset, a translation and a small rotation of the edge
   relative to the surface; and the rigid shift of the white surface that maximises the white/grey
   contrast (T1 1 mm inside minus 1 mm outside the surface).
2. The cortex near the scalp: exact distance of every white-surface vertex to the dense scalp used by
   G3B; area and Desikan-Killiany parcels closer than 8 mm (the referee's scalp-plus-skull bound, also
   the modelled skull depth of A-BEM-CHILD) and 10 mm (the referee's target count), with the G3B
   targets below both (results/g3b/g3b_targets_<key>.csv); by head sector; and, where the T1 exists,
   the same distances to the MRI head boundary.
3. The scalp against the MRI head boundary: T1 profiles along the scalp's outward normals (smoothed
   over 10 neighbour rings, as anatomy.model_skull). On each profile the boundary (primary) is where
   the intensity falls, going outward, through the volume's air/tissue level (Otsu's threshold of the
   whole T1): the outermost such fall, then (second pass) the fall nearest to the clipped and smoothed
   first-pass offsets of the neighbourhood, so that the pinna's edge is not taken for the head's.
   Sensitivity: the half-maximum edge (the outermost crossing of the level half-way between the air and
   the scalp's peak within 6 mm inside that fall; it moves inward where the subcutaneous fat is much
   brighter than the skin) and the steepest fall. Signed offset, positive where the MRI boundary lies
   outside the surface (the surface is inside the head), summarised over the cap (head frame z >= 0, as
   pediatric.head_size) and the whole head, and decomposed as in 1. Where the head mask that mkheadsurf
   tessellated (mri/seghead.mgz) is available: the surface against the mask's edge, and the T1's edge
   against the mask's.
4. Head circumference (occipitofrontal, ``pediatric.head_size``) on the surface used, on the MRI
   boundary points of 3 and on the outer boundary of a T1 head mask (at the Otsu level, holes filled,
   largest component).
5. Overlays: T1 slices (axial, coronal, sagittal) through the white vertex closest to the scalp, with
   the dense scalp, the modelled outer and inner skull (the BEM as used), the white surface and the
   MRI boundary points; for the children also their watershed outer skin and inner skull.
6. Verdict per anatomy (rules in the JSON: tolerance one voxel of the 1-mm T1; the adult, whose scalp
   was made the same way, is the reference for a systematic offset, which must hold for the Otsu and
   the half-maximum edge; misregistration is judged on the frames, the white surface and the head mask)
   and its consequence for G3B/G4.
7. Correction (``corrected_scalp``): the dense scalp moved along its normals to the MRI boundary (the
   cap's offsets cleaned and smoothed over the mesh, extended smoothly below the cap and tapered to none
   at the neck's cut edge; local folds relaxed); then measured again (closure). Validated on the adult:
   its scalp displaced inward by 4 mm, and separately translated by 2 mm (and its white surface
   translated by 2 mm), then recovered from its T1. For a child found inside (or outside) its MRI scalp,
   the corrected scalp is written (corrected_scalp_<key>.fif) for a rerun of the preparation, G3B and G4.
   Every anatomy with a T1 also gets its distances to this corrected surface (the MRI head boundary).

Surface-only checks that need no T1 (all anatomies): the watershed outer skin, outer skull and inner
skull against the dense scalp (signed, over the cap; the adult's own segmented BEM as reference).

Inputs: the T1 volumes and head masks. Adult: MNE-sample-data/subjects/sample/mri/T1.mgz and
seghead.mgz; templates: their mri/T1.mgz; children: data/external/school_subjects/<subject>/mri/T1.mgz
and seghead.mgz (configs/school_subjects_qc_manifest.json: S3 object versions, sizes and SHA-256; not
among the 45 files of configs/school_subjects_manifest.json). An anatomy without its T1 gets the
surface-only checks and the verdict 'undetermined'.
Outputs (default results/g3b_children_qc/; ``--out`` for tests): children_qc.json (derived quantities
only), Figure_QC_overlay_<key>.png, Figure_QC_summary.png, Figure_QC_near_scalp.png and, for a child
found inside its MRI scalp, corrected_scalp_<key>.fif.

  .venv/bin/python scripts/study_children_qc.py                    # default anatomies -> results/g3b_children_qc/
  .venv/bin/python scripts/study_children_qc.py --out DIR --anatomies adult childB
"""
from __future__ import annotations

import argparse
import gzip
import re
import struct
import sys
import time
import tomllib
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from scipy import ndimage  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from opmsquid import anatomy, fsio, io, paths, pediatric, plotting  # noqa: E402

OUT = ROOT / "results" / "g3b_children_qc"
STATUS = ("NEW check (referee round 1): the school-aged children's anatomy against their MRIs "
          "(D-G3-ANAT, A-BEM-CHILD, A-G3-FID): frames, cortex near the scalp, scalp vs MRI head boundary, "
          "head circumference, overlays, verdict, correction")
TEMPLATES = {"infant2yr": "ANTS2-0Years3T", "infant18mo": "ANTS18-0Months3T", "infant12mo": "ANTS12-0Months3T"}  # as G3B
DEFAULT_ANATOMIES = ("adult", "infant2yr", "childA", "childB", "childC")

# Declared analysis parameters (choices of this check, not results; written to the JSON)
PROFILE_MM = (-20.0, 20.0, 0.25)  # samples along the scalp normal: first, last, step [mm]
AIR_FROM_MM = 12.0  # air level: median of the cap's profile samples at least this far outside the surface
PEAK_WINDOW_MM = 6.0  # scalp peak: maximum within this distance inside the outermost Otsu crossing (within the scalp:
#                       a brighter tissue deeper in, within the window, would bias the half level and the edge inward)
GRAD_WINDOW_MM = 3.0  # steepest fall (sensitivity) searched within this distance of the Otsu crossing
WHITE_WINDOW_MM = 2.5  # white-surface edge searched within +- this distance along the white normal
WHITE_CONTRAST_MM = 1.0  # white/grey contrast: T1 this far inside minus this far outside the white surface
WHITE_SAMPLE = 30000  # white vertices sampled (evenly) for the image-based registration checks
SHIFT_COARSE_MM = (6.0, 2.0)  # rigid shift search: range and step of the coarse grid [mm]
SHIFT_FINE_MM = (1.5, 0.5)  # and of the fine grid around the coarse optimum
SMOOTH_RINGS = 10  # normal and offset smoothing (neighbour averages), anatomy.model_skull's default
TOLERANCE_MM = 1.0  # verdict tolerance: one voxel of the 1-mm T1
PRIMARY = "otsu"  # the MRI head boundary: where the T1 falls through the volume's air/tissue (Otsu) level; 'half_max' and
#                   'steepest' are sensitivity definitions (a half-maximum of the local peak moves inward where the
#                   subcutaneous fat is much brighter than the skin, as in child C)
HEADER_TOLERANCE_MM = 0.01  # volume geometries agree if no corner of the volume moves more than this
NEAR_MM = (8.0, 10.0)  # Fable: scalp plus skull >= 8 mm (= A-BEM-CHILD depth); targets shallower than 10 mm
TEST_INWARD_MM = 4.0  # correction test (adult): uniform inward displacement, the order of child B's apparent deficit
TEST_SHIFT_MM = (2.0, 0.0, 0.0)  # registration test (adult): rigid translation of the scalp and of the white surface
SECTOR_MM = 30.0  # head sectors (head frame): left x < -30, right x > 30, front y > 30, back y < -30, top z > 60 mm
TOP_MM = 60.0


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pct(x, q=(5, 10, 25, 50, 75, 90, 95)) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return {f"p{p:g}": float(np.percentile(x, p)) for p in q} if len(x) else {f"p{p:g}": None for p in q}


# ----------------------------------------------------------------------------------------------
# volumes and frames
MGH_DTYPES = {0: ">u1", 1: ">i4", 3: ">f4", 4: ">i2"}


def read_mgh(path: Path) -> dict:
    """First frame of a FreeSurfer MGH/MGZ volume (float32, voxel axes i, j, k) with its scanner
    vox2ras and its tkregister vox2ras (surface RAS = MNE's MRI frame) [mm], as nibabel's
    MGHHeader.get_affine and get_vox2ras_tkr compute them (the study's environment has no nibabel)."""
    opener = gzip.open if str(path).endswith((".mgz", ".gz")) else open
    with opener(path, "rb") as f:
        raw = f.read()
    version, w, h, d, nframes, code, _ = struct.unpack(">7i", raw[:28])
    if version != 1 or code not in MGH_DTYPES:
        raise ValueError(f"{path}: unsupported MGH file (version {version}, data type {code})")
    good_ras = int(struct.unpack(">h", raw[28:30])[0])
    spacing = np.array(struct.unpack(">3f", raw[30:42]), float)
    mdc = np.array(struct.unpack(">9f", raw[42:78]), float).reshape(3, 3).T  # columns: x, y, z direction cosines
    c_ras = np.array(struct.unpack(">3f", raw[78:90]), float)
    dims = np.array([w, h, d])
    data = np.frombuffer(raw, np.dtype(MGH_DTYPES[code]), count=w * h * d, offset=284).reshape((w, h, d), order="F")
    m = mdc * spacing
    vox2ras = np.eye(4)
    vox2ras[:3, :3] = m
    vox2ras[:3, 3] = c_ras - m @ (dims / 2.0)
    tkr = np.array([[-spacing[0], 0, 0, spacing[0] * w / 2], [0, 0, spacing[2], -spacing[2] * d / 2],
                    [0, -spacing[1], 0, spacing[1] * h / 2], [0, 0, 0, 1]], float)
    tail = raw[284 + w * h * d * nframes * np.dtype(MGH_DTYPES[code]).itemsize:]
    history = []  # FreeSurfer's processing record after the voxel data: program, arguments (base names), version, time
    for rec in re.finditer(rb"ProgramName: (\S+)\s+ProgramArguments: (.*?)\s+ProgramVersion: (\S+)\s+TimeStamp: (\S+)", tail):
        args = " ".join(Path(a).name for a in rec.group(2).decode(errors="replace").split())
        history.append(f"{args} (FreeSurfer {rec.group(3).decode(errors='replace')}, {rec.group(4).decode(errors='replace')})")
    return dict(data=data.astype(np.float32), dims=dims, spacing=spacing, mdc=mdc, c_ras=c_ras, good_ras=good_ras,
                nframes=int(nframes), vox2ras=vox2ras, vox2ras_tkr=tkr, path=Path(path), history=history)


def info_vox2ras(info: dict) -> np.ndarray:
    """Scanner vox2ras [mm] of the volume a FreeSurfer surface was made from (its volume information)."""
    dims, vs = np.asarray(info["volume"], float), np.asarray(info["voxelsize"], float)
    m = np.column_stack([info["xras"], info["yras"], info["zras"]]).astype(float) * vs
    out = np.eye(4)
    out[:3, :3] = m
    out[:3, 3] = np.asarray(info["cras"], float) - m @ (dims / 2.0)
    return out


def corner_displacement(a: np.ndarray, b: np.ndarray, dims) -> float:
    """Largest distance [mm] between the scanner positions that two vox2ras give the volume's corners."""
    c = np.array([[i, j, k, 1.0] for i in (0, dims[0]) for j in (0, dims[1]) for k in (0, dims[2])])
    return float(np.max(np.linalg.norm((c @ a.T - c @ b.T)[:, :3], axis=1)))


def source_volume(info: dict) -> str:
    """The last four components of the path of the volume a surface was made from (which run made it)."""
    return "/".join(Path(str(info.get("filename", ""))).parts[-4:])


def frame_checks(files: dict, t1: dict | None) -> dict:
    """Volume geometry stored in each FreeSurfer surface file against the first file's and the T1's."""
    geo = {}
    for name, f in files.items():
        _, _, info = fsio.read_geometry(f, read_metadata=True)
        geo[name] = (info_vox2ras(info), np.asarray(info["volume"], int), source_volume(info)) if info else None
    ref = next(k for k, v in geo.items() if v is not None)
    out = dict(reference_file=ref, files={})
    for name, g in geo.items():
        if g is None:
            out["files"][name] = dict(volume_info=False)
            continue
        v2r, dims, src = g
        rec = dict(volume_info=True, source_volume=src, dims=dims.tolist(),
                   corner_displacement_vs_reference_mm=corner_displacement(v2r, geo[ref][0], dims))
        if t1 is not None:
            rec["same_dims_as_t1"] = bool(np.array_equal(dims, t1["dims"]))
            rec["corner_displacement_vs_t1_mm"] = corner_displacement(v2r, t1["vox2ras"], t1["dims"])
        out["files"][name] = rec
    vals = [r["corner_displacement_vs_reference_mm"] for r in out["files"].values() if r.get("volume_info")]
    out["max_corner_displacement_between_files_mm"] = float(max(vals))
    out["source_volumes"] = sorted({r["source_volume"] for r in out["files"].values() if r.get("volume_info")})
    if t1 is not None:
        vt = [r["corner_displacement_vs_t1_mm"] for r in out["files"].values() if r.get("volume_info")]
        out["max_corner_displacement_vs_t1_mm"] = float(max(vt))
        out["t1"] = dict(file="/".join(t1["path"].parts[-4:]), dims=t1["dims"].tolist(), voxel_mm=t1["spacing"].tolist(),
                         c_ras=t1["c_ras"].tolist(), good_ras_flag=t1["good_ras"], frames=t1["nframes"], history=t1["history"])
        out["headers_match_t1"] = bool(max(vt) <= HEADER_TOLERANCE_MM and all(
            r.get("same_dims_as_t1", True) for r in out["files"].values()))
    return out


# ----------------------------------------------------------------------------------------------
# meshes
def adjacency(tris: np.ndarray, n: int) -> tuple[sp.csr_matrix, np.ndarray]:
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    adj = sp.coo_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()
    adj.data[:] = 1.0
    return adj, np.asarray(adj.sum(axis=1)).ravel()


def smooth_vectors(nn: np.ndarray, adj, deg, rings: int = SMOOTH_RINGS) -> np.ndarray:
    """Unit vectors averaged over ``rings`` neighbour rings (the scheme of anatomy.model_skull)."""
    n = np.asarray(nn, float).copy()
    for _ in range(int(rings)):
        n = 0.5 * n + 0.5 * (adj @ n) / deg[:, None]
        n /= np.linalg.norm(n, axis=1, keepdims=True)
    return n


def smooth_field(x: np.ndarray, adj, deg, rings: int = SMOOTH_RINGS) -> np.ndarray:
    x = np.asarray(x, float).copy()
    for _ in range(int(rings)):
        x = 0.5 * x + 0.5 * (adj @ x) / deg
    return x


def outward_normals(rr: np.ndarray, tris: np.ndarray) -> np.ndarray:
    surf = mne.surface.complete_surface_info(dict(rr=np.asarray(rr, float), tris=np.asarray(tris), np=len(rr), ntri=len(tris)),
                                             copy=False, verbose=False)
    return anatomy._outward(surf).nn


class SignedDistance:
    """Exact point-to-mesh distance [unit of the mesh] (closest point on the triangles around the ``k``
    nearest vertices, as anatomy.MeshDistance), with the sign of the closest triangle's outward normal
    (positive outside) when vertex normals ``nn`` are given. The sign is meant for points near the
    surface (BEM layers against the scalp), where it is unambiguous; MeshDistance's inside test is
    impractical on a 300,000-triangle scalp."""

    def __init__(self, rr: np.ndarray, tris: np.ndarray, nn: np.ndarray | None = None, k: int = 10):
        self.rr, self.tris, self.k = np.asarray(rr, float), np.asarray(tris), k
        self.tree = cKDTree(self.rr)
        order = np.argsort(self.tris.ravel(), kind="stable")
        tri_of = order // 3
        counts = np.bincount(self.tris.ravel(), minlength=len(self.rr))
        start = np.r_[0, np.cumsum(counts)[:-1]]
        width = max(int(counts.max()), 1)
        cols = np.minimum(start[:, None] + np.arange(width)[None, :], (start + np.maximum(counts, 1) - 1)[:, None])
        self.incident = tri_of[np.clip(cols, 0, len(tri_of) - 1)]  # padded by repeats (the minimum is unchanged)
        self.fn = None
        if nn is not None:
            fn = np.cross(self.rr[self.tris[:, 1]] - self.rr[self.tris[:, 0]], self.rr[self.tris[:, 2]] - self.rr[self.tris[:, 0]])
            fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-300
            s = np.sign(np.sum(fn * np.asarray(nn, float)[self.tris].mean(axis=1), axis=1))
            self.fn = fn * np.where(s == 0, 1.0, s)[:, None]

    def __call__(self, points: np.ndarray, chunk: int = 1000) -> np.ndarray:
        points = np.atleast_2d(points)
        _, near = self.tree.query(points, k=self.k)
        out = np.empty(len(points))
        for lo in range(0, len(points), chunk):
            p = points[lo:lo + chunk, None, :]
            t = self.incident[near[lo:lo + chunk]].reshape(len(p), -1)
            a, b, c = (self.rr[self.tris[t, j]] for j in range(3))
            cp = anatomy.closest_point_on_triangles(p, a, b, c)
            d = np.linalg.norm(cp - p, axis=-1)
            i = np.argmin(d, axis=1)
            r = np.arange(len(p))
            if self.fn is None:
                out[lo:lo + chunk] = d[r, i]
                continue
            s = np.sign(np.sum((p[:, 0, :] - cp[r, i]) * self.fn[t[r, i]], axis=1))
            out[lo:lo + chunk] = np.where(s == 0, 1.0, s) * d[r, i]
        return out


def plane_segments(rr: np.ndarray, tris: np.ndarray, axis: int, value: float) -> np.ndarray:
    """Segments (n, 2, 3) where a triangle mesh crosses the plane rr[:, axis] = value."""
    d = rr[:, axis] - (value + 1e-6)
    sd = np.sign(d[tris])
    t = tris[(sd.min(axis=1) < 0) & (sd.max(axis=1) > 0)]
    if not len(t):
        return np.zeros((0, 2, 3))
    pts = []
    for i0, i1 in ((0, 1), (1, 2), (2, 0)):
        a, b = t[:, i0], t[:, i1]
        da, db = d[a], d[b]
        with np.errstate(divide="ignore", invalid="ignore"):
            p = rr[a] + (da / (da - db))[:, None] * (rr[b] - rr[a])
        pts.append(np.where((da * db < 0)[:, None], p, np.nan))
    p = np.stack(pts, axis=1)
    keep = ~np.isnan(p[..., 0])
    return p[keep].reshape(-1, 2, 3)


# ----------------------------------------------------------------------------------------------
# profiles and edges
def sample_along(vol: dict, rr_m: np.ndarray, nn: np.ndarray, t_mm: np.ndarray, chunk: int = 10000) -> np.ndarray:
    """T1 intensity (trilinear) at rr + t * nn for every point (rows) and offset t [mm] (columns)."""
    inv = np.linalg.inv(vol["vox2ras_tkr"])
    out = np.empty((len(rr_m), len(t_mm)), np.float32)
    for lo in range(0, len(rr_m), chunk):
        p = np.asarray(rr_m[lo:lo + chunk], float)[:, None, :] * 1e3 + t_mm[None, :, None] * nn[lo:lo + chunk, None, :]
        ijk = p.reshape(-1, 3) @ inv[:3, :3].T + inv[:3, 3]
        out[lo:lo + chunk] = ndimage.map_coordinates(vol["data"], ijk.T, order=1, mode="constant", cval=0.0).reshape(len(p), -1)
    return out


def otsu_threshold(values: np.ndarray, nbins: int = 256) -> float:
    """Otsu's threshold of the intensities (maximum between-class variance over a 256-bin histogram)."""
    h, edges = np.histogram(values, bins=nbins)
    c = 0.5 * (edges[1:] + edges[:-1])
    w0 = np.cumsum(h).astype(float)
    w1 = w0[-1] - w0
    s0 = np.cumsum(h * c)
    m0, m1 = s0 / np.maximum(w0, 1), (s0[-1] - s0) / np.maximum(w1, 1)
    return float(c[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def head_boundary(prof: np.ndarray, t: np.ndarray, otsu: float, air: float, reference: np.ndarray | None = None) -> dict:
    """Offsets [mm] of the MRI head boundary along each profile (rows: samples at offsets t along the
    outward normal). The boundary is a fall of the intensity through the volume's Otsu level (tissue
    to air, going outward): the outermost one, or with ``reference`` (one offset per profile) the one
    nearest to it, so that the edge of a structure outside the head (the pinna, beside the head above
    the ears) is not taken for the head's. Returned: 'otsu' (primary, PRIMARY: that fall), 'half_max'
    (the outermost crossing of the level half-way between the air level and the scalp's peak, within
    PEAK_WINDOW_MM inside and GRAD_WINDOW_MM outside that fall; the peak is the maximum within
    PEAK_WINDOW_MM inside it) and 'steepest' (the steepest fall within GRAD_WINDOW_MM of it). NaN where
    undetermined: no fall through the Otsu level, or tissue at the profile's outer end."""
    n, m = prof.shape
    dt = float(t[1] - t[0])
    rows = np.arange(n)
    above = prof >= otsu
    falls = above[:, :-1] & ~above[:, 1:]  # between samples j and j + 1
    a, b = prof[:, :-1], prof[:, 1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.clip(np.where(a > b, (a - otsu) / (a - b), 0.0), 0.0, 1.0)
    t_falls = np.where(falls, t[:-1][None, :] + frac * dt, np.nan)
    ok = falls.any(axis=1) & ~above[:, -1]
    if reference is None:
        j = (m - 2) - np.argmax(falls[:, ::-1], axis=1)  # the outermost fall
    else:
        j = np.argmin(np.where(falls, np.abs(t_falls - np.asarray(reference, float)[:, None]), np.inf), axis=1)
    j = np.where(ok, j, 0)
    t_otsu = t_falls[rows, j]
    w, g = int(round(PEAK_WINDOW_MM / dt)), int(round(GRAD_WINDOW_MM / dt))
    idx = np.arange(m)[None, :]
    inside = (idx >= (j - w)[:, None]) & (idx <= j[:, None])
    peak = np.where(inside, prof, -np.inf).max(axis=1)
    half = air + 0.5 * (peak - air)
    # the outermost sample at or above the half level near this fall (the half level can lie below the Otsu
    # level, its crossing then just outside); searched from the outside, so that a dark layer inside the scalp
    # (skull, CSF) cannot be taken for the boundary
    near = (idx >= (j - w)[:, None]) & (idx <= np.minimum(j + g, m - 2)[:, None])
    at = (prof >= half[:, None]) & near
    jh = m - 1 - np.argmax(at[:, ::-1], axis=1)
    jh = np.minimum(jh, m - 2)
    ok &= at.any(axis=1) & (prof[rows, jh + 1] < half)
    ah, bh = prof[rows, jh], prof[rows, jh + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_half = t[jh] + np.clip(np.where(ah > bh, (ah - half) / (ah - bh), 0.0), 0.0, 1.0) * dt
    grad = np.diff(prof, axis=1) / dt
    tm = t[:-1] + dt / 2
    close = np.abs(tm[None, :] - np.nan_to_num(t_otsu)[:, None]) <= GRAD_WINDOW_MM
    t_grad = tm[np.argmin(np.where(close, grad, np.inf), axis=1)]
    nan = np.full(n, np.nan)
    return dict(half_max=np.where(ok, t_half, nan), otsu=np.where(ok, t_otsu, nan), steepest=np.where(ok, t_grad, nan),
                peak=np.where(ok, peak, nan), n_falls=falls.sum(axis=1))


def consensus(offset_mm: np.ndarray, adj, deg, rings: int = SMOOTH_RINGS) -> np.ndarray:
    """A smooth reference offset per vertex: the offsets clipped to 3 robust SD about their median
    (undetermined ones set to it), smoothed over ``rings`` neighbour rings."""
    x = np.asarray(offset_mm, float)
    med = float(np.nanmedian(x))
    s = 1.4826 * float(np.nanmedian(np.abs(x - med)))
    x = np.where(np.isfinite(x), np.clip(x, med - 3.0 * s, med + 3.0 * s), med)
    return smooth_field(x, adj, deg, rings)


def rigid_decomposition(offset_mm: np.ndarray, rr_m: np.ndarray, nn: np.ndarray, n_iter: int = 5) -> dict:
    """offset ~ a + n.b + (r x n).w by least squares, trimmed at 3 robust SD (r about the points'
    centroid): the uniform offset a [mm], translation b [mm] and small rotation w of the edge relative
    to the surface; the rotation is also given as the largest displacement it causes on these points.
    On a nearly spherical cap b and w trade off (a rotation about the sphere's centre moves nothing
    along the normals), so the verdict uses the rigid part's normal offset n.b + (r x n).w itself
    ('rigid_normal_offset': its 95th percentile of absolute values), which the fit determines."""
    ok = np.isfinite(offset_mm)
    if ok.sum() < 100:
        return dict(n_points=int(ok.sum()))
    r = np.asarray(rr_m, float)[ok] * 1e3
    r = r - r.mean(axis=0)
    n = np.asarray(nn, float)[ok]
    y = np.asarray(offset_mm, float)[ok]
    x = np.column_stack([np.ones(len(y)), n, np.cross(r, n)])
    keep = np.ones(len(y), bool)
    for _ in range(int(n_iter)):
        coef = np.linalg.lstsq(x[keep], y[keep], rcond=None)[0]
        res = y - x @ coef
        s = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        keep = np.abs(res - np.median(res[keep])) <= 3.0 * max(s, 1e-3)
    a, b, w = float(coef[0]), coef[1:4], coef[4:7]
    rot_disp = np.linalg.norm(np.cross(w, r), axis=1)
    rigid_normal = np.abs(x[:, 1:] @ coef[1:])
    return dict(n_points=int(ok.sum()), kept_share=float(keep.mean()), uniform_mm=a, translation_mm=b.tolist(),
                translation_norm_mm=float(np.linalg.norm(b)), rotation_deg=np.degrees(w).tolist(),
                rotation_norm_deg=float(np.degrees(np.linalg.norm(w))), rotation_max_displacement_mm=float(rot_disp.max()),
                rigid_normal_offset_p50_mm=float(np.median(rigid_normal)),
                rigid_normal_offset_p95_mm=float(np.percentile(rigid_normal, 95)),
                residual_rms_mm=float(np.sqrt(np.mean(res[keep] ** 2))),
                explained_share=float(1.0 - np.var(res[keep]) / max(np.var(y[keep]), 1e-12)))


def scalp_vs_mri(vol: dict, rr_m: np.ndarray, nn: np.ndarray, cap: np.ndarray, tris: np.ndarray) -> dict:
    """Profiles along the scalp normals and the boundary offset of every vertex (``head_boundary``) in
    two passes: the outermost fall first, then the fall nearest to the first pass's ``consensus`` (its
    neighbours' clipped and smoothed offsets). The air level and the Otsu threshold are the volume's
    (air: over the cap's profiles)."""
    t = np.arange(PROFILE_MM[0], PROFILE_MM[1] + PROFILE_MM[2] / 2, PROFILE_MM[2])
    otsu = otsu_threshold(vol["data"].ravel())  # every voxel: air against tissue
    prof = sample_along(vol, rr_m, nn, t)
    pc = prof[cap]
    air = float(np.median(pc[:, t >= AIR_FROM_MM]))
    adj, deg = adjacency(np.asarray(tris), len(rr_m))
    first = head_boundary(prof, t, otsu, air)
    edge = head_boundary(prof, t, otsu, air, reference=consensus(first[PRIMARY], adj, deg))
    edge["first_pass"] = first[PRIMARY]
    median_profile = np.median(pc, axis=0)
    q25, q75 = np.percentile(pc, [25, 75], axis=0)
    del pc
    return dict(t=t, edge=edge, otsu=otsu, air=air, median_profile=median_profile, q25=q25, q75=q75)


def summarise_offsets(edge: dict, cap: np.ndarray) -> dict:
    out = {}
    for region, m in (("cap", cap), ("whole_head", np.ones_like(cap))):
        rec = {}
        for k in ("half_max", "otsu", "steepest"):
            x = edge[k][m]
            v = x[np.isfinite(x)]
            rec[k] = dict(n=int(m.sum()), determined_share=float(np.isfinite(x).mean()),
                          median_mm=float(np.median(v)) if len(v) else None, **{f"{q}_mm": val for q, val in pct(v).items()})
            if len(v):
                rec[k].update({f"share_above_{s:g}mm": float(np.mean(v > s)) for s in (1.0, 2.0, 4.0)},
                              share_below_minus_1mm=float(np.mean(v < -1.0)))
        out[region] = rec
    return out


def white_alignment(vol: dict, rr_m: np.ndarray, nn: np.ndarray, shift_mm=(0.0, 0.0, 0.0)) -> dict:
    """Image-based registration of the white surface (optionally displaced by ``shift_mm``): the edge
    (steepest fall outward within +-WHITE_WINDOW_MM; parabolic refinement) decomposed by
    ``rigid_decomposition``; the white/grey contrast at zero shift; the rigid shift (coarse then fine
    grid) that maximises the mean contrast, i.e. the shift that would realign the surface."""
    step = max(1, len(rr_m) // WHITE_SAMPLE)
    rr = np.asarray(rr_m[::step], float) + np.asarray(shift_mm, float) * 1e-3
    n = np.asarray(nn[::step], float)
    dt = 0.25
    t = np.arange(-WHITE_WINDOW_MM - dt, WHITE_WINDOW_MM + 1.5 * dt, dt)
    prof = sample_along(vol, rr, n, t)
    g = np.diff(prof, axis=1) / dt
    tm = t[:-1] + dt / 2
    inside = np.abs(tm) <= WHITE_WINDOW_MM
    gi = np.where(inside[None, :], g, np.inf)
    jg = np.argmin(gi, axis=1)
    jl, jr = np.clip(jg - 1, 0, len(tm) - 1), np.clip(jg + 1, 0, len(tm) - 1)
    rows = np.arange(len(g))
    gl, g0, gr = g[rows, jl], g[rows, jg], g[rows, jr]
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where((gl - 2 * g0 + gr) > 0, 0.5 * (gl - gr) / (gl - 2 * g0 + gr), 0.0)
    edge = tm[jg] + np.clip(np.nan_to_num(delta), -0.5, 0.5) * dt
    edge[(g0 >= 0) | (jg == 0) | (jg == len(tm) - 1)] = np.nan  # no fall inside the window
    inv = np.linalg.inv(vol["vox2ras_tkr"])

    def contrast(s_mm):
        p = rr * 1e3 + np.asarray(s_mm, float)
        pin, pout = p - WHITE_CONTRAST_MM * n, p + WHITE_CONTRAST_MM * n
        vals = [ndimage.map_coordinates(vol["data"], (q @ inv[:3, :3].T + inv[:3, 3]).T, order=1, mode="constant") for q in (pin, pout)]
        return float(np.mean(vals[0] - vals[1]))

    def search(centre, half, step_mm):
        g1 = np.arange(-half, half + step_mm / 2, step_mm)
        best, arg = -np.inf, np.zeros(3)
        for dx in g1:
            for dy in g1:
                for dz in g1:
                    s = np.asarray(centre) + (dx, dy, dz)
                    c = contrast(s)
                    if c > best:
                        best, arg = c, s
        return arg, best

    c0 = contrast((0.0, 0.0, 0.0))
    coarse, _ = search((0.0, 0.0, 0.0), *SHIFT_COARSE_MM)
    fine, cbest = search(coarse, *SHIFT_FINE_MM)
    return dict(n_vertices=int(len(rr)), edge_median_mm=float(np.nanmedian(edge)), edge=pct(edge),
                edge_determined_share=float(np.isfinite(edge).mean()), decomposition=rigid_decomposition(edge, rr, n),
                contrast_at_zero=c0, best_shift_mm=fine.tolist(), best_shift_norm_mm=float(np.linalg.norm(fine)),
                contrast_at_best=cbest, contrast_gain_share=float(cbest / c0 - 1.0) if c0 > 0 else None)


# ----------------------------------------------------------------------------------------------
# correction
def fill_from_neighbours(x: np.ndarray, adj) -> np.ndarray:
    """NaN entries filled ring by ring with the mean of their determined neighbours (a smooth
    extension of the determined values over the mesh)."""
    x = np.asarray(x, float).copy()
    for _ in range(100000):
        miss = np.isnan(x)
        if not miss.any():
            break
        num, den = adj @ np.nan_to_num(x), adj @ (~miss).astype(float)
        fill = miss & (den > 0)
        if not fill.any():
            x[miss] = np.nanmedian(x) if np.isfinite(x).any() else 0.0
            break
        x[fill] = num[fill] / den[fill]
    return x


def clean_offsets(offset_mm: np.ndarray, use: np.ndarray, adj, deg, rings: int = SMOOTH_RINGS) -> tuple[np.ndarray, dict]:
    """Displacement field from the MRI-boundary offsets: only the vertices in ``use`` (the cap, where the
    profiles are clean) are taken, clipped to 3 robust SD about their median (a spatially coherent
    misdetection, e.g. along the pinna, which reaches into the cap, is not an outlier of its
    neighbourhood and would otherwise be spread below the cap); undetermined ones and outliers (more
    than 3 robust SD from the field smoothed over ``rings`` neighbour rings) are refilled from their
    neighbours; the field is extended smoothly to the rest of the head (``fill_from_neighbours``) and
    smoothed."""
    x = np.where(use, np.asarray(offset_mm, float), np.nan)
    n_used = int(np.isfinite(x).sum())
    med = float(np.nanmedian(x))
    spread = 1.4826 * float(np.nanmedian(np.abs(x - med)))
    n_clipped = int(np.sum(np.abs(x - med) > 3.0 * spread))
    x = np.where(np.isfinite(x), np.clip(x, med - 3.0 * spread, med + 3.0 * spread), np.nan)
    filled = fill_from_neighbours(x, adj)
    sm = smooth_field(filled, adj, deg, rings)
    res = filled - sm
    s = 1.4826 * np.median(np.abs(res[use] - np.median(res[use])))
    out = use & (np.abs(res - np.median(res[use])) > 3.0 * max(s, 1e-3))
    x[out] = np.nan
    final = smooth_field(fill_from_neighbours(x, adj), adj, deg, rings)
    return final, dict(n_vertices=len(x), offsets_used=n_used, clipped=n_clipped, clip_range_mm=[med - 3.0 * spread, med + 3.0 * spread],
                       outliers_replaced=int(out.sum()), extended_to=int((~use).sum()))


def untangle(rr: np.ndarray, tris: np.ndarray, ref_rr: np.ndarray, adj, deg, max_iter: int = 200,
             bound_mm: float = TOLERANCE_MM) -> tuple[np.ndarray, dict]:
    """Triangles whose orientation flipped relative to ``ref_rr`` (an offset surface folds where a
    concave feature, e.g. behind the ear, is tighter than the displacement; a sliver triangle flips
    under a tiny displacement difference) relaxed locally by shrink-free (Taubin) smoothing steps of
    their vertices and one ring around them; no vertex moves more than ``bound_mm`` from its displaced
    position (it is then frozen), so the surface stays within one voxel of the measured boundary. The
    state with the fewest flipped triangles is returned with a report (flips before and after and
    their area in the reference surface: the folds left are reported, not hidden)."""
    def face(x):
        return np.cross(x[tris[:, 1]] - x[tris[:, 0]], x[tris[:, 2]] - x[tris[:, 0]])

    ref = face(ref_rr)
    area = 0.5 * np.linalg.norm(ref, axis=1) * 1e6  # mm^2
    rr = np.asarray(rr, float).copy()
    start = rr.copy()

    def flipped(x):
        return np.einsum("ij,ij->i", ref, face(x)) <= 0

    bad = flipped(rr)
    before = (int(bad.sum()), float(area[bad].sum()))
    best_n, best = before[0], rr.copy()
    frozen = np.zeros(len(rr), bool)
    it = 0
    for it in range(1, max_iter + 1):
        if not bad.any():
            break
        v = np.zeros(len(rr), bool)
        v[np.unique(tris[bad])] = True
        v |= (adj @ v.astype(float)) > 0
        v &= ~frozen
        if not v.any():
            break
        for f in (0.5, -0.53):
            lap = (adj @ rr) / deg[:, None] - rr
            rr[v] += f * lap[v]
        moved = np.linalg.norm(rr - start, axis=1) * 1e3
        over = moved > bound_mm
        rr[over] = start[over] + (rr[over] - start[over]) * (bound_mm / moved[over])[:, None]
        frozen |= over
        bad = flipped(rr)
        if bad.sum() < best_n:
            best_n, best = int(bad.sum()), rr.copy()
    bad = flipped(best)
    moved = np.linalg.norm(best - start, axis=1) * 1e3
    return best, dict(flipped_before=before[0], flipped_area_before_mm2=before[1], flipped_after=int(bad.sum()),
                      flipped_area_after_mm2=float(area[bad].sum()), surface_area_mm2=float(area.sum()), iterations=it,
                      relaxed_vertices=int(np.sum(moved > 1e-6)), relaxation_max_mm=float(moved.max()))


def taper_weight(height_m: np.ndarray) -> np.ndarray:
    """1 from the top of the head down to half-way between the fiducial plane (head frame z = 0) and
    the bottom of the mesh (the field-of-view cut of the neck), then linearly to 0 at the bottom: the
    correction leaves the cut edge, where an offset surface folds, as it is."""
    z = np.asarray(height_m, float)
    lo = float(z.min())
    return np.clip((z - lo) / max(-0.5 * lo, 1e-6), 0.0, 1.0) if lo < 0 else np.ones_like(z)


def corrected_scalp(scalp: anatomy.Surface, offset_mm: np.ndarray, normals: np.ndarray, use: np.ndarray,
                    weight: np.ndarray | None = None) -> tuple[anatomy.Surface, dict]:
    """The dense scalp moved along ``normals`` (its smoothed outward unit normals) to the MRI head
    boundary: the displacement is the offset measured on the vertices in ``use`` (the cap), cleaned,
    extended smoothly over the rest of the head and smoothed (``clean_offsets``), times ``weight``
    (``taper_weight`` of the head-frame height: none at the neck's cut edge); folds are relaxed locally
    within one voxel (``untangle``). Same vertices and triangles, so the steps that use the scalp's
    vertices (A-BEM-CONFORM, the fiducial projection, target depths, head size) can take it as the MRI
    scalp directly; triangles still flipped (small concave features and slivers of the dense
    FreeSurfer mesh fold under an offset of a few mm) are reported with their area, in and outside
    the cap, not hidden."""
    adj, deg = adjacency(np.asarray(scalp.tris), len(scalp.rr))
    disp, rep = clean_offsets(offset_mm, use, adj, deg)
    if weight is not None:
        disp = disp * np.asarray(weight, float)
    tris = np.asarray(scalp.tris)
    rr0 = np.asarray(scalp.rr, float)
    rr = rr0 + disp[:, None] * 1e-3 * np.asarray(normals, float)
    rr, unt = untangle(rr, tris, rr0, adj, deg)
    f0 = np.cross(rr0[tris[:, 1]] - rr0[tris[:, 0]], rr0[tris[:, 2]] - rr0[tris[:, 0]])
    f1 = np.cross(rr[tris[:, 1]] - rr[tris[:, 0]], rr[tris[:, 2]] - rr[tris[:, 0]])
    bad = np.einsum("ij,ij->i", f0, f1) <= 0
    in_cap = use[tris].all(axis=1)
    area = 0.5 * np.linalg.norm(f0, axis=1) * 1e6
    unt.update(flipped_in_cap=int(np.sum(bad & in_cap)), flipped_area_in_cap_mm2=float(area[bad & in_cap].sum()),
               cap_area_mm2=float(area[in_cap].sum()))
    rep.update(displacement_mm=pct(disp, (1, 5, 50, 95, 99)), displacement_cap_mm=pct(disp[use], (1, 5, 50, 95, 99)), folds=unt)
    return anatomy.Surface(rr, tris, outward_normals(rr, tris)), rep


# ----------------------------------------------------------------------------------------------
# anatomies
def load(key: str, cfg: dict):
    """(subject as G3B loads it, T1 path (may not exist), white-surface files, the segmentation's BEM
    surface files as distributed, the child's configuration entry or None)."""
    school = {c["key"]: c for c in cfg["anatomy"]["school"]}
    if key == "adult":
        s = anatomy.load_sample()
        sd = paths.SUBJECTS_DIR / "sample"
        bem = {f: sd / "bem" / f"{f}.surf" for f in ("outer_skin", "outer_skull", "inner_skull")}
        return s, sd / "mri" / "T1.mgz", [sd / "surf" / f"{h}.white" for h in ("lh", "rh")], bem, None
    if key in TEMPLATES:
        name = TEMPLATES[key]
        s = anatomy.load_template(name)
        td = paths.EXTERNAL / anatomy.INFANT_SUBJECTS / name
        bem = {f: td / "bem" / f"{f}.surf" for f in ("outer_skin", "outer_skull", "inner_skull")}
        return s, td / "mri" / "T1.mgz", [td / "surf" / f"{h}.white" for h in ("lh", "rh")], bem, None
    if key in school:
        c = school[key]
        root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS
        s = anatomy.load_school(c["subject"])
        bem = {f: root / c["bem_folder"] / "surf" / f"{f}.surf" for f in ("outer_skin", "outer_skull", "inner_skull")}
        whites = [root / c["subject"] / "surf" / f"{h}.white" for h in ("lh", "rh")]
        return s, root / c["subject"] / "mri" / "T1.mgz", whites, bem, c
    raise ValueError(f"unknown anatomy {key}")


def head_frame(subject, rr_m: np.ndarray) -> np.ndarray:
    mh = np.linalg.inv(subject.trans["trans"])
    return np.asarray(rr_m, float) @ mh[:3, :3].T + mh[:3, 3]


def ofc_of_points(subject, points_m: np.ndarray) -> dict:
    """pediatric.head_size (occipitofrontal circumference and dimensions) of a point cloud in the
    subject's MRI frame, with the subject's head frame."""
    pseudo = SimpleNamespace(scalp=SimpleNamespace(rr=np.asarray(points_m, float)), trans=subject.trans,
                             fiducials=subject.fiducials)
    return pediatric.head_size(pseudo)


def mask_boundary_points(vol: dict, level: float) -> tuple[np.ndarray, dict]:
    """Outer boundary voxels [m, MRI frame] of the T1 head mask (intensity >= level, holes filled,
    largest connected component), and how many mask voxels lie on each face of the volume (a head cut
    by the field of view, whose true scalp no surface made from this T1 could reach)."""
    mask = ndimage.binary_fill_holes(vol["data"] >= level)
    lab, n = ndimage.label(mask)
    if n > 1:
        mask = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    edge = mask & ~ndimage.binary_erosion(mask)
    ijk = np.argwhere(edge).astype(float)
    tkr = vol["vox2ras_tkr"]
    faces = {}
    for axis in range(3):
        for side, sl in (("first", 0), ("last", -1)):
            count = int(np.take(mask, sl, axis=axis).sum())
            direction = tkr[:3, axis] * (-1 if side == "first" else 1)  # outward direction of this face, MRI frame
            name = "RAS"[int(np.argmax(np.abs(direction)))] + ("+" if direction[np.argmax(np.abs(direction))] > 0 else "-")
            faces[name] = count
    return (ijk @ tkr[:3, :3].T + tkr[:3, 3]) * 1e-3, faces


def region_names(subject, cortex) -> np.ndarray:
    names = np.concatenate([plotting.read_freesurfer_annot(subject.labels / f"{h}.aparc.annot") for h in ("lh", "rh")])
    if len(names) != cortex.n:
        raise ValueError(f"{subject.name}: aparc has {len(names)} vertices, the white surface {cortex.n}")
    return np.array([f"{'lh' if h == 0 else 'rh'}.{nm}" for h, nm in zip(cortex.hemi, names)], dtype=object)


def near_surface(d_mm: np.ndarray, cortex, regions: np.ndarray, hf: np.ndarray) -> dict:
    """Distances [mm] of the white vertices to a head surface: percentiles, and per threshold the
    vertices, area and parcels closer than it; the 5th percentile per head sector."""
    area = cortex.area * 1e4  # cm^2
    out = dict(min_mm=float(d_mm.min()), **{f"{k}_mm": v for k, v in pct(d_mm, (0.1, 1, 5, 10, 50)).items()})
    for thr in NEAR_MM:
        m = d_mm < thr
        parcels = {}
        for r in np.unique(regions[m]):
            sel = m & (regions == r)
            parcels[r] = dict(area_cm2=float(area[sel].sum()), n_vertices=int(sel.sum()), min_mm=float(d_mm[sel].min()))
        out[f"below_{thr:g}mm"] = dict(n_vertices=int(m.sum()), area_cm2=float(area[m].sum()),
                                       n_usable_vertices=int((m & cortex.usable).sum()),
                                       parcels=dict(sorted(parcels.items(), key=lambda kv: -kv[1]["area_cm2"])))
    sectors = dict(left=hf[:, 0] < -SECTOR_MM * 1e-3, right=hf[:, 0] > SECTOR_MM * 1e-3, front=hf[:, 1] > SECTOR_MM * 1e-3,
                   back=hf[:, 1] < -SECTOR_MM * 1e-3, top=hf[:, 2] > TOP_MM * 1e-3)
    out["p5_by_sector_mm"] = {k: float(np.percentile(d_mm[m], 5)) for k, m in sectors.items() if m.any()}
    return out


def g3b_targets(key: str, d_used: np.ndarray | None = None, d_mri: np.ndarray | None = None, cortex=None) -> dict | None:
    """G3B's targets of this anatomy below NEAR_MM (their stored depth: nearest dense-scalp vertex),
    and where the MRI boundary is known, the same targets' depth change."""
    f = paths.RESULTS / "g3b" / f"g3b_targets_{key}.csv"
    if not f.exists():
        return None
    rows = io.read_csv(f)
    depth = np.array([float(r["depth_mm"]) for r in rows])
    area = np.array([float(r["area_mm2"]) for r in rows])
    out = dict(n_targets=len(rows), depth=pct(depth, (1, 5, 50)))
    for thr in NEAR_MM:
        m = depth < thr
        out[f"below_{thr:g}mm"] = dict(n=int(m.sum()), area_cm2=float(area[m].sum() / 100.0))
    if d_mri is not None and cortex is not None:
        n_lh = int(np.sum(cortex.hemi == 0))
        g = np.array([int(r["vertno"]) + (n_lh if int(r["hemi"]) == 1 else 0) for r in rows])
        change = d_mri[g] - d_used[g]
        out["depth_change_to_mri_boundary_mm"] = pct(change, (5, 50, 95))
        for thr in NEAR_MM:
            out[f"below_{thr:g}mm_to_mri_boundary"] = int(np.sum(d_mri[g] < thr))
    return out


def watershed_vs_scalp(subject, bem_files: dict, sd: SignedDistance, cap_of) -> dict:
    """Signed distance [mm] of the segmentation's BEM surfaces (as distributed, before A-BEM-CONFORM
    and A-BEM-CHILD) to the dense scalp used, over the cap (positive: outside the scalp)."""
    out = {}
    for name, f in bem_files.items():
        if not f.exists():
            continue
        rr, _ = fsio.read_geometry(f)
        rr = rr / 1e3
        d = sd(rr) * 1e3
        cap = cap_of(rr)
        out[name] = dict(cap=dict(median_mm=float(np.median(d[cap])), **{f"{k}_mm": v for k, v in pct(d[cap], (10, 90)).items()}),
                         whole=dict(median_mm=float(np.median(d))))
    return out


def prepared_file_checks(key: str, c: dict, subject, white_files) -> dict:
    """The prepared MNE files of a child against the FreeSurfer surfaces they were written from."""
    root = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS / c["subject"]
    seg, _ = fsio.read_geometry(root / "surf" / "lh.seghead")
    head = mne.read_bem_surfaces(root / "bem" / f"{c['subject']}-head.fif", verbose=False)[0]
    whites = [fsio.read_geometry(f)[0] for f in white_files]
    src_dev = max(float(np.max(np.abs(s["rr"] - w / 1e3))) for s, w in zip(subject.src, whites))
    bem_head = next(b for b in subject.bem_surfaces if b["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    on_scalp = cKDTree(subject.scalp.rr).query(bem_head["rr"])[0]
    return dict(head_fif_vs_seghead_max_um=float(np.max(np.abs(head["rr"] - seg / 1e3)) * 1e6),
                source_space_vs_white_max_um=src_dev * 1e6, bem_head_on_scalp_max_um=float(on_scalp.max() * 1e6),
                head_frame_from="fiducials transferred from the adult (A-G3-FID)")


# ----------------------------------------------------------------------------------------------
# figures
SURF_STYLE = dict(scalp=("gold", 1.2, "solid", "dense scalp used"), outer_skull=("deepskyblue", 1.0, "solid", "outer skull (BEM used)"),
                  inner_skull=("limegreen", 1.0, "solid", "inner skull (BEM used)"), white=("red", 0.6, "solid", "white surface"),
                  watershed_skin=("darkviolet", 0.9, "dashed", "segmentation's outer skin (as distributed)"),
                  watershed_inner=("darkorange", 0.9, "dashed", "segmentation's inner skull (as distributed)"))
PLANES = ((2, "axial", (0, 1), ("R", "A")), (1, "coronal", (0, 2), ("R", "S")), (0, "sagittal", (1, 2), ("A", "S")))


def overlay_figure(path: Path, title: str, centre_mm: np.ndarray, meshes: dict, vol: dict | None, boundary_mm: np.ndarray | None,
                   mask: dict | None = None):
    """2 x 3 panels: axial, coronal and sagittal planes through ``centre_mm`` (MRI frame), full head and
    a 50-mm zoom: the T1 (if any) with every mesh's intersection, the MRI boundary points near the plane
    and the outline of mkheadsurf's head mask (``mask``, if any)."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10.4))
    vmax = float(np.percentile(vol["data"][vol["data"] > 0], 99.5)) if vol is not None else None
    head = meshes["scalp"][0]
    cap = head[head[:, 2] >= np.percentile(head[:, 2], 25)]  # the head above its lowest quarter (not the neck)
    centre_head = 0.5 * (cap.min(axis=0) + cap.max(axis=0))
    half_head = float(0.5 * np.max(cap.max(axis=0) - cap.min(axis=0)) + 15.0)
    for col, (axis, name, (u, v), (lu, lv)) in enumerate(PLANES):
        value = float(centre_mm[axis])
        for row, half in ((0, half_head), (1, 25.0)):
            ax = axes[row, col]
            cu, cv = (float(centre_head[u]), float(centre_head[v])) if row == 0 else (float(centre_mm[u]), float(centre_mm[v]))
            if vol is not None:
                res = 0.5 if row == 0 else 0.2
                gu = np.arange(cu - half, cu + half + res / 2, res)
                gv = np.arange(cv - half, cv + half + res / 2, res)
                uu, vv = np.meshgrid(gu, gv)
                p = np.zeros(uu.shape + (3,))
                p[..., u], p[..., v], p[..., axis] = uu, vv, value
                inv = np.linalg.inv(vol["vox2ras_tkr"])
                ijk = p.reshape(-1, 3) @ inv[:3, :3].T + inv[:3, 3]
                img = ndimage.map_coordinates(vol["data"], ijk.T, order=1, mode="constant").reshape(uu.shape)
                ax.imshow(img, cmap="gray", origin="lower", extent=(gu[0], gu[-1], gv[0], gv[-1]), vmin=0, vmax=vmax)
                if mask is not None:
                    inv_m = np.linalg.inv(mask["vox2ras_tkr"])
                    mimg = ndimage.map_coordinates(mask["data"], (p.reshape(-1, 3) @ inv_m[:3, :3].T + inv_m[:3, 3]).T, order=1,
                                                   mode="constant").reshape(uu.shape)
                    ax.contour(uu, vv, mimg, levels=[0.5 * float(mask["data"].max())], colors="cyan", linewidths=0.8, linestyles="dotted")
                    ax.plot([], [], ":", color="cyan", lw=0.8, label="mkheadsurf head mask (seghead.mgz)")
            else:
                ax.set_facecolor("0.93")
            for k, (rr, tris) in meshes.items():
                seg = plane_segments(rr, tris, axis, value)
                if len(seg):
                    color, lw, ls, label = SURF_STYLE[k]
                    ax.add_collection(LineCollection(seg[:, :, [u, v]], colors=color, linewidths=lw, linestyles=ls, label=label))
            if boundary_mm is not None:
                b = boundary_mm[np.abs(boundary_mm[:, axis] - value) < 0.5]
                ax.plot(b[:, u], b[:, v], ".", color="magenta", ms=1.5, label="MRI head boundary (profiles)")
            ax.plot(centre_mm[u], centre_mm[v], "+", color="cyan", ms=14, mew=2)
            ax.set_xlim(cu - half, cu + half)
            ax.set_ylim(cv - half, cv + half)
            ax.set_aspect("equal")
            ax.set_xlabel(f"{lu} [mm, MRI frame]")
            ax.set_ylabel(f"{lv} [mm]")
            ax.set_title(f"{name} {'xyz'[axis]} = {value:.1f} mm" + (" (zoom)" if row else ""), fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    for ax in axes.ravel()[1:]:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96), h_pad=2.5)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def summary_figure(path: Path, recs: dict, distances: dict, profiles: dict):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    colors = {"adult": "k", "infant2yr": "tab:red", "infant18mo": "tab:orange", "infant12mo": "tab:purple",
              "childA": "tab:brown", "childB": "tab:pink", "childC": "tab:olive"}
    ax = axes[0]
    for key, d in distances.items():
        for which, ls in (("used", "-"), ("mri", "--")):
            if d.get(which) is not None:
                x = np.sort(d[which])
                ax.plot(x, np.arange(1, len(x) + 1) / len(x), ls, color=colors.get(key, "gray"),
                        label=f"{key} ({'scalp used' if which == 'used' else 'MRI boundary'})")
    for thr in NEAR_MM:
        ax.axvline(thr, color="0.6", lw=0.8, ls=":")
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 0.25)
    ax.set_xlabel("white-surface vertex to head surface [mm]")
    ax.set_ylabel("cumulative share of white vertices")
    ax.set_title("cortex near the scalp (lowest quarter)")
    ax.legend(fontsize=7)
    ax = axes[1]
    for key, p in profiles.items():
        ax.plot(p["t"], p["median_profile"], color=colors.get(key, "gray"), label=key)
        ax.fill_between(p["t"], p["q25"], p["q75"], color=colors.get(key, "gray"), alpha=0.15)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlabel("offset along the outward scalp normal [mm] (0 = surface used)")
    ax.set_ylabel("T1 intensity (median, quartiles; cap)")
    ax.set_title("T1 across the scalp surface used")
    ax.legend(fontsize=8)
    ax = axes[2]
    keys = list(recs)
    first_mri, first_ws = True, True
    for i, key in enumerate(keys):
        sm = recs[key].get("scalp_vs_mri")
        if sm:
            q = np.nanpercentile(sm["_cap_offsets"], [5, 25, 50, 75, 95])
            ax.plot([i, i], [q[0], q[4]], color="k", lw=1, label="MRI boundary: 5-95 % (cap)" if first_mri else None)
            ax.plot([i, i], [q[1], q[3]], color="k", lw=6, solid_capstyle="butt", alpha=0.5,
                    label="MRI boundary: quartiles" if first_mri else None)
            ax.plot(i, q[2], "o", color="white", mec="k", ms=7, label="MRI boundary: median" if first_mri else None)
            first_mri = False
        else:
            ax.text(i, 0.3, "no T1", ha="center", fontsize=8, color="0.3")
        ws = recs[key].get("watershed_vs_scalp", {}).get("outer_skin")
        if ws:
            ax.plot(i + 0.15, ws["cap"]["median_mm"], "D", color="tab:blue", ms=6,
                    label="segmentation's outer skin: median" if first_ws else None)
            first_ws = False
        mk = recs[key].get("scalp_vs_seghead_mask")
        if mk:
            ax.plot(i - 0.15, mk["mask_edge_minus_surface"]["median_mm"], "^", color="tab:cyan", mec="k", ms=7,
                    label="mkheadsurf mask edge: median" if not any(h.get_label().startswith("mkheadsurf") for h in ax.lines) else None)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(range(len(keys)), keys)
    ax.set_ylim(-5, 8)
    ax.set_ylabel("minus the surface used [mm] (cap; + = surface inside)")
    ax.set_title("MRI head boundary and segmentation outer skin\nagainst the dense scalp used", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def near_scalp_figure(path: Path, maps: dict):
    """White vertices closer than 10 mm to the scalp used (colour: distance), over the cortex outline."""
    keys = list(maps)
    views = (("left", 1, 2, 0, -1), ("top", 0, 1, 2, 1), ("right", 1, 2, 0, 1), ("back", 0, 2, 1, -1))
    fig, axes = plt.subplots(len(keys), len(views), figsize=(4 * len(views), 3.4 * len(keys)), squeeze=False)
    for i, key in enumerate(keys):
        hf, d = maps[key]
        for j, (name, a, b, depth_axis, sign) in enumerate(views):
            ax = axes[i, j]
            side = (sign * hf[:, depth_axis]) > (0 if name in ("left", "right", "back") else 0.03)
            ax.plot(hf[side, a][::7] * 1e3, hf[side, b][::7] * 1e3, ",", color="0.8")
            m = side & (d < NEAR_MM[1])
            sc = ax.scatter(hf[m, a] * 1e3, hf[m, b] * 1e3, c=d[m], s=1.5, cmap="inferno", vmin=3, vmax=NEAR_MM[1])
            ax.set_aspect("equal")
            ax.set_title(f"{key}: {name}", fontsize=9)
            ax.tick_params(labelsize=7)
        fig.colorbar(sc, ax=axes[i, -1], fraction=0.05, label="white to scalp [mm]")
    fig.suptitle("White-surface vertices closer than 10 mm to the scalp used (head frame, mm)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=100)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# analysis of one anatomy
def analyse(key: str, cfg: dict, out_dir: Path, ref: dict) -> tuple[dict, dict]:
    t0 = time.time()
    subject, t1_path, white_files, bem_files, child = load(key, cfg)
    cortex = anatomy.full_resolution(subject)
    regions = region_names(subject, cortex)
    rec = dict(subject=subject.name, description=subject.description)
    scalp = subject.scalp
    adj, deg = adjacency(np.asarray(scalp.tris), len(scalp.rr))
    nn_s = smooth_vectors(scalp.nn, adj, deg)
    height = head_frame(subject, scalp.rr)[:, 2]
    cap = height >= 0
    rec["scalp"] = dict(n_vertices=int(len(scalp.rr)), cap_vertices=int(cap.sum()))
    t1 = read_mgh(t1_path) if t1_path.exists() else None
    rec["t1_available"] = t1 is not None
    rec["t1_expected_at"] = str(t1_path.relative_to(ROOT)) if t1_path.is_relative_to(ROOT) else str(t1_path)
    mask_path = t1_path.with_name("seghead.mgz")  # mkheadsurf's head mask, which the dense scalp was tessellated from
    mask = read_mgh(mask_path) if (t1 is not None and mask_path.exists()) else None

    # 1. frames
    files = {f.name: f for f in white_files}
    if child is not None:
        files["lh.seghead"] = paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS / child["subject"] / "surf" / "lh.seghead"
    files.update({f"{k}.surf (BEM as distributed)": f for k, f in bem_files.items() if f.exists()})
    rec["frames"] = frame_checks(files, t1)
    if mask is not None:
        rec["frames"]["seghead_mask"] = dict(
            file="/".join(mask["path"].parts[-4:]), dims=mask["dims"].tolist(), values=np.unique(mask["data"])[:5].tolist(),
            corner_displacement_vs_t1_mm=corner_displacement(mask["vox2ras"], t1["vox2ras"], t1["dims"]), history=mask["history"])
    if child is not None:
        rec["prepared_files"] = prepared_file_checks(key, child, subject, white_files)
        rec["surfaces_used"] = dict(
            frame="MNE's MRI frame = FreeSurfer surface RAS (tkregister RAS) of the conformed 256^3 1-mm T1, in m",
            dense_scalp=f"{child['subject']}/surf/lh.seghead (FreeSurfer mkheadsurf, from the seghead.mgz mask), written unchanged "
                        f"to bem/{child['subject']}-head.fif: the scalp of G3B/G4 (OPM placement, helmet contact, depths)",
            bem_head="the watershed outer_skin.surf (stored in the preceding subject's folder: "
                     f"{child['bem_folder']}) with every vertex moved to the nearest dense-scalp vertex (A-BEM-CONFORM)",
            skull="modelled (A-BEM-CHILD): inner skull 8 mm below the dense scalp where the watershed's was shallower, outer "
                  "skull half-way to the scalp",
            white=f"{child['subject']}/surf/lh.white, rh.white; source volumes in frames.files")
    log(f"{key}: frames done ({time.time() - t0:.0f} s); T1 {'found' if t1 is not None else 'NOT available'}")

    # 2. cortex near the scalp used
    d_used = SignedDistance(scalp.rr, scalp.tris)(cortex.rr) * 1e3
    hf_cortex = head_frame(subject, cortex.rr)
    rec["white_to_scalp_used"] = near_surface(d_used, cortex, regions, hf_cortex)
    closest = int(np.argmin(d_used))
    rec["closest_white_vertex"] = dict(region=str(regions[closest]), distance_mm=float(d_used[closest]),
                                       mri_frame_mm=(cortex.rr[closest] * 1e3).tolist())
    log(f"{key}: white-to-scalp done ({time.time() - t0:.0f} s): min {d_used.min():.2f} mm in {regions[closest]}")

    # surface-only: segmentation BEM against the scalp
    sd = SignedDistance(scalp.rr, scalp.tris, scalp.nn)
    rec["watershed_vs_scalp"] = watershed_vs_scalp(subject, bem_files, sd, lambda rr: head_frame(subject, rr)[:, 2] >= 0)

    d_mri, boundary_mm, corrected = None, None, None
    extra = dict(d_used=d_used, hf=hf_cortex)
    if t1 is not None:
        # 3. scalp against the MRI head boundary
        res = scalp_vs_mri(t1, scalp.rr, nn_s, cap, scalp.tris)
        edge = res["edge"]
        changed = np.abs(edge["first_pass"] - edge[PRIMARY]) > TOLERANCE_MM
        rec["scalp_vs_mri"] = dict(primary=PRIMARY, otsu_threshold=res["otsu"], air_level=res["air"],
                                   scalp_peak_median=float(np.nanmedian(edge["peak"][cap])),
                                   cap_profiles_with_several_falls_share=float(np.mean(edge["n_falls"][cap] > 1)),
                                   cap_second_pass_changed_share=float(np.mean(changed[cap])),
                                   offsets=summarise_offsets(edge, cap),
                                   decomposition_cap=rigid_decomposition(np.where(cap, edge[PRIMARY], np.nan), scalp.rr, nn_s),
                                   decomposition_whole_head=rigid_decomposition(edge[PRIMARY], scalp.rr, nn_s),
                                   _cap_offsets=edge[PRIMARY][cap], _cap_height=height[cap])
        extra["profile"] = res
        log(f"{key}: scalp vs MRI done ({time.time() - t0:.0f} s): median cap offset "
            f"{rec['scalp_vs_mri']['offsets']['cap'][PRIMARY]['median_mm']:+.2f} mm ({PRIMARY}; half-max "
            f"{rec['scalp_vs_mri']['offsets']['cap']['half_max']['median_mm']:+.2f} mm)")
        if mask is not None:
            # the dense scalp against the mask it was tessellated from (registration of surface and mask), and the
            # mask's edge against the T1's (is the mask itself inside the head?)
            # the mask's edge is its half-maximum (0.5 between 0 and 1); Otsu's threshold of a two-valued volume is
            # degenerate (any level between the values separates them equally well)
            medge = scalp_vs_mri(mask, scalp.rr, nn_s, cap, scalp.tris)["edge"]["half_max"]
            diff = edge[PRIMARY] - medge
            rec["scalp_vs_seghead_mask"] = dict(
                mask_edge_minus_surface=summarise_offsets(dict(half_max=medge, otsu=medge, steepest=medge), cap)["cap"]["half_max"],
                t1_edge_minus_mask_edge_cap=dict(median_mm=float(np.nanmedian(diff[cap])), **{f"{k}_mm": v for k, v in pct(diff[cap]).items()}),
                t1_edge_minus_mask_edge_whole_head_median_mm=float(np.nanmedian(diff)))
            extra["mask_edge_cap"] = medge[cap]
            log(f"{key}: mask edge {rec['scalp_vs_seghead_mask']['mask_edge_minus_surface']['median_mm']:+.2f} mm from the surface, "
                f"T1 edge {rec['scalp_vs_seghead_mask']['t1_edge_minus_mask_edge_cap']['median_mm']:+.2f} mm from the mask edge (cap medians)")
        # MRI head boundary as a surface (the correction applied), distances to it
        try:
            corrected, crep = corrected_scalp(scalp, edge[PRIMARY], nn_s, cap, taper_weight(height))
        except ValueError as err:
            corrected, crep = None, dict(error=str(err))
        rec["mri_boundary_surface"] = crep
        if corrected is not None:
            nn_c = smooth_vectors(corrected.nn, adj, deg)  # closure: the corrected surface measured again
            again = scalp_vs_mri(t1, corrected.rr, nn_c, cap, corrected.tris)["edge"][PRIMARY]
            crep["remeasured_cap_offset_mm"] = pct(again[cap], (5, 25, 50, 75, 95))
            d_mri = SignedDistance(corrected.rr, corrected.tris)(cortex.rr) * 1e3
            rec["white_to_mri_boundary"] = near_surface(d_mri, cortex, regions, hf_cortex)
            extra["d_mri"] = d_mri
        ok = np.isfinite(edge[PRIMARY])
        boundary_mm = (scalp.rr[ok] * 1e3 + edge[PRIMARY][ok, None] * nn_s[ok])
        # 4. head circumference (the boundary points of the cap: the planes searched lie 0-60 mm above the fiducials)
        level = res["otsu"]  # the head mask at the boundary's own level
        cap_ok = ok & cap
        cap_points = (scalp.rr[cap_ok] * 1e3 + edge[PRIMARY][cap_ok, None] * nn_s[cap_ok]) * 1e-3
        mask_points, fov_faces = mask_boundary_points(t1, level)
        rec["head_circumference"] = dict(
            surface_used=pediatric.head_size(subject), mri_boundary_points=ofc_of_points(subject, cap_points),
            t1_mask=ofc_of_points(subject, mask_points), t1_mask_level=level)
        rec["frames"]["head_mask_voxels_on_volume_faces"] = fov_faces  # R+/R-/A+/A-/S+/S-: a head cut by the field of view
        # 1b. image-based registration of the white surface
        rec["white_vs_t1"] = white_alignment(t1, cortex.rr, cortex.nn)
        # is the cortex closest to the scalp genuine (white matter inside, grey outside)? T1 at -1 and +1 mm
        tt = np.array([-WHITE_CONTRAST_MM, WHITE_CONTRAST_MM])
        every = np.arange(0, cortex.n, max(1, cortex.n // WHITE_SAMPLE))
        near = np.flatnonzero(d_used < NEAR_MM[0])
        cmp = {}
        for name, idx in (("all", every), (f"below_{NEAR_MM[0]:g}mm", near)):
            if len(idx):
                p = sample_along(t1, cortex.rr[idx], cortex.nn[idx], tt)
                cmp[name] = dict(n=int(len(idx)), inside_median=float(np.median(p[:, 0])), outside_median=float(np.median(p[:, 1])),
                                 contrast_median=float(np.median(p[:, 0] - p[:, 1])), share_without_contrast=float(np.mean(p[:, 0] <= p[:, 1])))
        rec["white_contrast_near_scalp"] = cmp
        log(f"{key}: white vs T1 done ({time.time() - t0:.0f} s): edge median {rec['white_vs_t1']['edge_median_mm']:+.2f} mm, "
            f"best shift {np.round(rec['white_vs_t1']['best_shift_mm'], 2).tolist()} mm")
    else:
        rec["head_circumference"] = dict(surface_used=pediatric.head_size(subject))
    white_area = float(cortex.area.sum() * 1e4)
    rec["head_circumference"]["white_surface_area_cm2"] = white_area
    rec["head_circumference"]["ofc_per_sqrt_white_area"] = float(rec["head_circumference"]["surface_used"]["ofc_mm"] / 10.0
                                                                 / np.sqrt(white_area))
    rec["g3b_targets"] = g3b_targets(key, d_used, d_mri, cortex)

    # 5. overlay
    surfs = {b["id"]: b for b in subject.bem_surfaces}
    meshes = dict(scalp=(scalp.rr * 1e3, scalp.tris),
                  outer_skull=(surfs[FIFF.FIFFV_BEM_SURF_ID_SKULL]["rr"] * 1e3, surfs[FIFF.FIFFV_BEM_SURF_ID_SKULL]["tris"]),
                  inner_skull=(surfs[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["rr"] * 1e3, surfs[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["tris"]),
                  white=(cortex.rr * 1e3, cortex.tris))
    if child is not None:
        for name, f in (("watershed_skin", "outer_skin"), ("watershed_inner", "inner_skull")):
            meshes[name] = fsio.read_geometry(bem_files[f])  # [mm]
    title = (f"{key} ({subject.name}): white vertex closest to the scalp used ({regions[closest]}, {d_used[closest]:.1f} mm)"
             + ("" if t1 is not None else " -- T1 NOT AVAILABLE LOCALLY: surfaces only"))
    overlay_figure(out_dir / f"Figure_QC_overlay_{key}.png", title, cortex.rr[closest] * 1e3, meshes, t1, boundary_mm, mask)

    # 6. verdict
    rec["verdict"] = verdict(key, rec, ref)
    if child is not None and corrected is not None and rec["verdict"]["class"] in ("inside", "outside"):
        surf = dict(id=FIFF.FIFFV_BEM_SURF_ID_HEAD, sigma=1.0, coord_frame=FIFF.FIFFV_COORD_MRI, rr=corrected.rr,
                    tris=corrected.tris, np=len(corrected.rr), ntri=len(corrected.tris))
        surf = mne.surface.complete_surface_info(surf, copy=False, verbose=False)
        mne.write_bem_surfaces(out_dir / f"corrected_scalp_{key}.fif", [surf], overwrite=True, verbose=False)
        rec["verdict"]["corrected_scalp_file"] = f"corrected_scalp_{key}.fif"
    rec["runtime_s"] = time.time() - t0
    log(f"{key}: {rec['verdict']['class']} ({rec['runtime_s']:.0f} s)")
    return rec, extra


def verdict(key: str, rec: dict, ref: dict) -> dict:
    """usable / inside (the surface lies inside the MRI head boundary) / outside / misregistered /
    undetermined (no T1); rules in the returned dict. The systematic offset is taken relative to the
    adult (``ref``: its scalp was made by the same FreeSurfer tool, mkheadsurf) and must hold for two
    edge definitions (the primary Otsu crossing and the half-maximum), since a dim outer skin layer over
    bright fat moves the half-maximum edge inward by the layer's thickness while the Otsu crossing stays
    outside it."""
    rules = (f"misregistered: the T1 header differs from the volume geometry of a surface or of the head mask by more than "
             f"{HEADER_TOLERANCE_MM} mm, or the white surface's realigning shift, edge translation or rotation displacement "
             f"exceeds {TOLERANCE_MM} mm, or the dense scalp lies more than {TOLERANCE_MM} mm (cap median) from the edge of "
             f"the head mask it was tessellated from; inside/outside: the cap median of the MRI boundary's offset exceeds "
             f"(falls short of) the adult reference's by more than {TOLERANCE_MM} mm for both the Otsu-crossing (primary) and "
             f"the half-maximum edge; usable: otherwise; undetermined: no T1. The rigid pattern of the scalp offsets "
             f"(n.b + (r x n).w) and the share of the cap beyond the tolerance are reported, not used: the dense scalp was made "
             f"from this very T1, so a pattern reflects the edge definition or the tissue, not a misregistration")
    w = rec["white_to_scalp_used"]
    indications = dict(white_to_scalp_p5_mm=w["p5_mm"], white_to_scalp_min_mm=w["min_mm"],
                       area_below_8mm_cm2=w["below_8mm"]["area_cm2"], area_below_10mm_cm2=w["below_10mm"]["area_cm2"],
                       p5_by_sector_mm=w["p5_by_sector_mm"])
    ws = rec.get("watershed_vs_scalp", {}).get("outer_skin")
    if ws:
        indications["segmentation_outer_skin_minus_scalp_cap_median_mm"] = ws["cap"]["median_mm"]
    hc = rec["head_circumference"]
    indications.update(ofc_surface_used_mm=hc["surface_used"]["ofc_mm"], white_surface_area_cm2=hc["white_surface_area_cm2"],
                       ofc_cm_per_sqrt_white_area_cm=hc["ofc_per_sqrt_white_area"])
    refs = {k: v for k, v in ref.get("white_to_scalp_p5_mm", {}).items() if k != key and not k.startswith("child")}
    if refs:  # the uniform outward move of the scalp that would give this anatomy the reference's 5th percentile
        indications["p5_shortfall_vs_references_mm"] = {k: float(v - w["p5_mm"]) for k, v in refs.items()}
    if not rec["t1_available"]:
        return dict(**{"class": "undetermined"}, rules=rules, reason=f"T1 not available at {rec['t1_expected_at']}",
                    surface_only_indications=indications)
    fr, wv, sm = rec["frames"], rec["white_vs_t1"], rec["scalp_vs_mri"]
    off = sm["offsets"]["cap"][PRIMARY]["median_mm"]
    off_h = sm["offsets"]["cap"]["half_max"]["median_mm"]
    ref_off = ref.get("adult_offset_mm", off) if key != "adult" else off
    ref_h = ref.get("adult_offset_half_max_mm", off_h) if key != "adult" else off_h
    excess, excess_h = off - ref_off, off_h - ref_h
    dec, wd = sm["decomposition_cap"], wv["decomposition"]
    cap_off = np.asarray(sm["_cap_offsets"], float)
    cap_h = np.asarray(sm["_cap_height"], float)
    det = np.isfinite(cap_off)
    out_of = det & (cap_off < ref_off - TOLERANCE_MM)  # the surface outside the boundary by more than the adult's position allows
    in_of = det & (cap_off > ref_off + TOLERANCE_MM)
    mk = rec.get("scalp_vs_seghead_mask")
    mask_frame = fr.get("seghead_mask", {}).get("corner_displacement_vs_t1_mm")
    checks = dict(headers_match_t1=fr["headers_match_t1"], mask_header_vs_t1_mm=mask_frame,
                  white_best_shift_mm=wv["best_shift_norm_mm"], white_edge_translation_mm=wd.get("translation_norm_mm"),
                  white_edge_rotation_max_displacement_mm=wd.get("rotation_max_displacement_mm"),
                  scalp_vs_mask_edge_cap_median_mm=mk["mask_edge_minus_surface"]["median_mm"] if mk else None,
                  scalp_offset_cap_median_mm=dict(otsu=off, half_max=off_h),
                  adult_reference_mm=dict(otsu=ref_off, half_max=ref_h),
                  scalp_offset_excess_over_adult_mm=dict(otsu=excess, half_max=excess_h),
                  cap_share_beyond_tolerance_of_adult_median=dict(outside=float(out_of.sum() / max(det.sum(), 1)),
                                                                   inside=float(in_of.sum() / max(det.sum(), 1))),
                  cap_height_where_outside_mm=pct(cap_h[out_of] * 1e3, (10, 50, 90)) if out_of.any() else None,
                  cap_height_where_inside_mm=pct(cap_h[in_of] * 1e3, (10, 50, 90)) if in_of.any() else None,
                  scalp_rigid_normal_offset_p95_mm=dec.get("rigid_normal_offset_p95_mm"))
    mis = (not fr["headers_match_t1"]) or (mask_frame is not None and mask_frame > HEADER_TOLERANCE_MM) or \
        wv["best_shift_norm_mm"] > TOLERANCE_MM or (wd.get("translation_norm_mm") or 0) > TOLERANCE_MM or \
        (wd.get("rotation_max_displacement_mm") or 0) > TOLERANCE_MM or \
        (mk is not None and abs(mk["mask_edge_minus_surface"]["median_mm"]) > TOLERANCE_MM)
    if mis:
        cls = "misregistered"
    elif excess > TOLERANCE_MM and excess_h > TOLERANCE_MM:
        cls = "inside"
    elif excess < -TOLERANCE_MM and excess_h < -TOLERANCE_MM:
        cls = "outside"
    else:
        cls = "usable"
    text = {"misregistered": "surface and T1 misregistered (see checks)",
            "inside": f"scalp surface inside the MRI head boundary by {off:.1f} mm (cap median, Otsu edge; half-maximum edge "
                      f"{off_h:.1f} mm), {excess:.1f} mm more than the adult's",
            "outside": f"scalp surface outside the MRI head boundary by {-off:.1f} mm (cap median, Otsu edge; half-maximum edge "
                       f"{-off_h:.1f} mm)",
            "usable": f"scalp surface on the MRI head boundary as the adult's is (cap median offset {off:+.1f} mm Otsu edge, "
                      f"{off_h:+.1f} mm half-maximum; adult {ref_off:+.1f} / {ref_h:+.1f} mm); white surface registered to the "
                      f"T1"}[cls]
    notes = []
    if abs(off - off_h) > TOLERANCE_MM:
        notes.append(f"the Otsu and half-maximum edges differ by {abs(off - off_h):.1f} mm (cap medians): where a dim outer skin "
                     f"layer lies over brighter fat, the half-maximum edge falls at the fat's edge")
    if out_of.mean() > 0 and det.any():
        h = checks["cap_height_where_outside_mm"]
        notes.append(f"the surface lies more than {TOLERANCE_MM:g} mm further out than the adult's position over "
                     f"{100 * out_of.sum() / det.sum():.0f} % of the cap (median height {h['p50']:.0f} mm above the fiducial plane)"
                     + (f"; more than {TOLERANCE_MM:g} mm further in over {100 * in_of.sum() / det.sum():.0f} %" if in_of.any() else ""))
    if (dec.get("rigid_normal_offset_p95_mm") or 0) > TOLERANCE_MM:
        notes.append(f"the offsets vary over the cap (rigid part p95 {dec['rigid_normal_offset_p95_mm']:.1f} mm, translation "
                     f"{np.round(dec['translation_mm'], 1).tolist()} mm); the surface follows its head mask, made from this T1, so "
                     f"this is a property of the edge definition or of the tissue, not a misregistration")
    out = dict(**{"class": cls}, text=text, rules=rules, checks=checks, notes=notes, surface_only_indications=indications)
    if rec.get("g3b_targets") and "depth_change_to_mri_boundary_mm" in rec["g3b_targets"]:
        gt = rec["g3b_targets"]
        out["consequence"] = dict(
            g3b_target_depth_change_to_mri_boundary_mm=gt["depth_change_to_mri_boundary_mm"],
            g3b_targets_below_10mm_used_vs_mri_boundary=[gt["below_10mm"]["n"], gt["below_10mm_to_mri_boundary"]],
            g3b_targets_below_8mm_used_vs_mri_boundary=[gt["below_8mm"]["n"], gt["below_8mm_to_mri_boundary"]],
            white_area_below_8mm_cm2_used_vs_mri_boundary=[w["below_8mm"]["area_cm2"],
                                                           rec["white_to_mri_boundary"]["below_8mm"]["area_cm2"]],
            white_to_scalp_min_mm_used_vs_mri_boundary=[w["min_mm"], rec["white_to_mri_boundary"]["min_mm"]],
            opm_cell_inner_face_to_mri_boundary_mm=float(2.0 - off), neuromag_top_contact_to_mri_boundary_mm=float(20.0 - off),
            excess_over_adult_mm=excess,
            note=("OPM sensing centres sit 7 mm (A-OPM-STANDOFF) outside the surface used, so the 10-mm cell's inner face "
                  "is 2 mm from it; the Neuromag top contact puts the nearest magnetometer 20 mm from it (A-G3-PLACE; Dewar "
                  "spacing 18 mm). Both are measured from the surface used; with the MRI boundary further out (positive "
                  "offset) both sit closer to the head than intended by the offset, and the targets are deeper by the depth "
                  "change. The child-adult contrasts (Delta) depend on the excess over the adult, whose scalp is offset the "
                  "same way."))
    return out


# ----------------------------------------------------------------------------------------------
# correction and registration tests on the adult
def correction_tests(out_dir: Path) -> dict:
    """The adult's scalp displaced inward by TEST_INWARD_MM along its smoothed normals, and separately
    translated by TEST_SHIFT_MM, then measured against its T1 and (inward case) corrected; the adult's
    white surface translated by TEST_SHIFT_MM and measured as in white_alignment."""
    s = anatomy.load_sample()
    t1 = read_mgh(paths.SUBJECTS_DIR / "sample" / "mri" / "T1.mgz")
    adj, deg = adjacency(np.asarray(s.scalp.tris), len(s.scalp.rr))
    nn_s = smooth_vectors(s.scalp.nn, adj, deg)
    height = head_frame(s, s.scalp.rr)[:, 2]
    cap, weight = height >= 0, taper_weight(height)
    base = scalp_vs_mri(t1, s.scalp.rr, nn_s, cap, s.scalp.tris)["edge"][PRIMARY]
    out = dict(declared_test_values=dict(inward_mm=TEST_INWARD_MM, translation_mm=list(TEST_SHIFT_MM)))
    # inward displacement
    deg_rr = s.scalp.rr - TEST_INWARD_MM * 1e-3 * nn_s
    deg_surf = anatomy.Surface(deg_rr, s.scalp.tris, outward_normals(deg_rr, s.scalp.tris))
    nn_d = smooth_vectors(deg_surf.nn, adj, deg)
    e = scalp_vs_mri(t1, deg_rr, nn_d, cap, s.scalp.tris)["edge"][PRIMARY]
    out["inward"] = dict(measured_offset_cap=pct(e[cap], (5, 50, 95)), adult_own_offset_cap=pct(base[cap], (5, 50, 95)),
                         recovered_displacement_cap_median_mm=float(np.nanmedian(e[cap] - base[cap])),
                         decomposition_cap=rigid_decomposition(np.where(cap, e, np.nan), deg_rr, nn_d))
    try:
        fixed, rep = corrected_scalp(deg_surf, e, nn_d, cap, weight)
        base_fixed, _ = corrected_scalp(s.scalp, base, nn_s, cap, weight)
        resid = SignedDistance(base_fixed.rr, base_fixed.tris, base_fixed.nn)(fixed.rr[cap]) * 1e3
        to_original = SignedDistance(s.scalp.rr, s.scalp.tris, s.scalp.nn)(fixed.rr[cap]) * 1e3
        out["inward"].update(correction=rep, corrected_minus_adult_corrected_cap=pct(resid, (5, 50, 95)),
                             corrected_minus_original_scalp_cap=pct(to_original, (5, 50, 95)))
    except ValueError as err:
        out["inward"]["correction"] = dict(error=str(err))
    # rigid translation of the scalp
    tr_rr = s.scalp.rr + np.asarray(TEST_SHIFT_MM) * 1e-3
    e2 = scalp_vs_mri(t1, tr_rr, nn_s, cap, s.scalp.tris)["edge"][PRIMARY]
    out["translation"] = dict(decomposition_cap=rigid_decomposition(np.where(cap, e2, np.nan), tr_rr, nn_s),
                              decomposition_whole_head=rigid_decomposition(e2, tr_rr, nn_s),
                              decomposition_cap_untranslated=rigid_decomposition(np.where(cap, base, np.nan), s.scalp.rr, nn_s),
                              expected_translation_mm=(-np.asarray(TEST_SHIFT_MM)).tolist())
    # white surface translated
    cortex = anatomy.full_resolution(s)
    out["white_translation"] = dict(measured=white_alignment(t1, cortex.rr, cortex.nn, shift_mm=TEST_SHIFT_MM),
                                    expected_realigning_shift_mm=(-np.asarray(TEST_SHIFT_MM)).tolist())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=OUT, help="output folder (default results/g3b_children_qc)")
    ap.add_argument("--anatomies", nargs="+", default=list(DEFAULT_ANATOMIES),
                    help="adult, infant2yr/infant18mo/infant12mo, childA/childB/childC")
    ap.add_argument("--skip-tests", action="store_true", help="skip the correction and registration tests on the adult")
    args = ap.parse_args(argv)
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    mne.set_log_level("WARNING")
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    keys = list(dict.fromkeys(["adult"] + list(args.anatomies)))  # the adult first: it is the verdict's reference
    out = dict(status=STATUS, requests=dict(
        fable="verify the children's scalp/white-surface geometry against the MRIs (scalp-to-cortex distance maps; compare to "
              "the dataset's T1) and either repair or drop them",
        codex="Provide MRI-based overlays, coordinate-transform checks and quantitative surface/registration quality measures. "
              "Correct and rerun affected analyses, or remove those anatomies from substantive conclusions."),
        parameters=dict(profile_mm=PROFILE_MM, air_from_mm=AIR_FROM_MM, peak_window_mm=PEAK_WINDOW_MM,
                        grad_window_mm=GRAD_WINDOW_MM, white_window_mm=WHITE_WINDOW_MM, white_contrast_mm=WHITE_CONTRAST_MM,
                        white_sample=WHITE_SAMPLE, shift_coarse_mm=SHIFT_COARSE_MM, shift_fine_mm=SHIFT_FINE_MM,
                        smooth_rings=SMOOTH_RINGS, tolerance_mm=TOLERANCE_MM, header_tolerance_mm=HEADER_TOLERANCE_MM,
                        near_mm=NEAR_MM, cap="head frame z >= 0 (plane of LPA, RPA, nasion), as pediatric.head_size",
                        sectors_mm=dict(lateral=SECTOR_MM, top=TOP_MM),
                        correction_taper="full above half the depth of the mesh below the fiducial plane, linear to none "
                                         "at the mesh's bottom (the neck's cut edge)",
                        tests_on_adult=dict(inward_mm=TEST_INWARD_MM, translation_mm=TEST_SHIFT_MM),
                        note="declared analysis choices of this check; thresholds 8 and 10 mm from the referee's request "
                             "(8 mm also configs/g3b_pediatric.toml child_skull_depth_mm); the adult test displacements are "
                             "declared test values"),
        anatomies={})
    ref, distances, profiles, maps = {}, {}, {}, {}
    for key in keys:
        rec, extra = analyse(key, cfg, args.out, ref)
        if key == "adult" and rec["t1_available"]:
            ref["adult_offset_mm"] = rec["scalp_vs_mri"]["offsets"]["cap"][PRIMARY]["median_mm"]
            ref["adult_offset_half_max_mm"] = rec["scalp_vs_mri"]["offsets"]["cap"]["half_max"]["median_mm"]
        ref.setdefault("white_to_scalp_p5_mm", {})[key] = rec["white_to_scalp_used"]["p5_mm"]
        distances[key] = dict(used=extra["d_used"], mri=extra.get("d_mri"))
        if "profile" in extra:
            profiles[key] = extra["profile"]
        maps[key] = (extra["hf"], extra["d_used"])
        out["anatomies"][key] = rec
    if not args.skip_tests:
        log("correction and registration tests on the adult ...")
        out["tests_on_adult"] = correction_tests(args.out)
    summary_figure(args.out / "Figure_QC_summary.png", out["anatomies"], distances, profiles)
    near_scalp_figure(args.out / "Figure_QC_near_scalp.png", maps)
    for rec in out["anatomies"].values():
        rec.get("scalp_vs_mri", {}).pop("_cap_offsets", None)
        rec.get("scalp_vs_mri", {}).pop("_cap_height", None)
    out["verdicts"] = {k: dict(**{"class": v["verdict"]["class"]}, text=v["verdict"].get("text", v["verdict"].get("reason")))
                       for k, v in out["anatomies"].items()}
    out["runtime_s"] = time.time() - t0
    io.write_json(out, args.out / "children_qc.json")
    log(f"done in {out['runtime_s']:.0f} s -> {args.out}")
    return out


if __name__ == "__main__":
    main()
