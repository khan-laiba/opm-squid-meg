"""Subject anatomy: source spaces, BEM surfaces, scalp, transforms and source descriptors.

Primary adult anatomy: the MNE ``sample`` subject (individual adult MRI with FreeSurfer surfaces,
3-layer BEM surfaces, oct-6 source space with distances and patch statistics, dense scalp, and
the head-MRI transform of its Neuromag recording). Coordinates are in the MRI frame [m] unless
stated otherwise.

Pediatric anatomies (G3B): ``load_template`` reads an infant template of O'Reilly et al. (2021)
as packaged by ``mne.datasets.fetch_infant_template`` (FreeSurfer surfaces, 3-layer BEM, dense
head surface, oct-6 source space, MRI-frame fiducials; native dimensions, head frame from its
fiducials); ``scaled`` makes a size-only control by scaling every coordinate of a subject about
its MRI origin (a scaled adult is not pediatric anatomy).

Source descriptors (Hunold et al. 2016 definitions, see docs/literature/hunold2016.md):
* depth: distance from the source to the scalp (nearest point of the dense scalp surface);
* orientation: angle between the source's cortical normal and the normal of the inner skull at
  the nearby inner-skull point (0 deg = normal to the skull, i.e. radial; 90 deg = tangential).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import mne
from mne.io.constants import FIFF
from scipy.spatial import cKDTree

from . import paths


@dataclass
class Surface:
    rr: np.ndarray  # (n, 3) MRI frame [m]
    tris: np.ndarray
    nn: np.ndarray  # outward unit vertex normals

    def subdivided(self, times: int = 1) -> "Surface":
        rr, nn, tris = self.rr, self.nn, self.tris
        for _ in range(times):
            edges = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
            edges, inv = np.unique(edges, axis=0, return_inverse=True)
            mid_nn = nn[edges[:, 0]] + nn[edges[:, 1]]
            mid_nn /= np.linalg.norm(mid_nn, axis=1, keepdims=True)
            m = inv.reshape(3, -1).T + len(rr)
            tris = np.concatenate([np.c_[tris[:, 0], m[:, 0], m[:, 2]], np.c_[tris[:, 1], m[:, 1], m[:, 0]],
                                   np.c_[tris[:, 2], m[:, 2], m[:, 1]], m])
            rr = np.vstack([rr, 0.5 * (rr[edges[:, 0]] + rr[edges[:, 1]])])
            nn = np.vstack([nn, mid_nn])
        return Surface(rr, tris, nn)


def _outward(surf: dict) -> Surface:
    rr, nn = surf["rr"], surf["nn"]
    sign = 1.0 if np.median(np.sum((rr - rr.mean(axis=0)) * nn, axis=1)) > 0 else -1.0
    return Surface(rr, surf["tris"], sign * nn)


@dataclass
class Subject:
    name: str
    subjects_dir: Path
    src: mne.SourceSpaces
    bem_surfaces: list  # MNE BEM surface dicts, outermost first (as stored)
    scalp: Surface  # dense head surface
    inner_skull: Surface
    trans: mne.transforms.Transform  # head -> MRI
    label_dir: Path | None = None  # FreeSurfer label directory (aparc annotations); default subjects_dir/name/label
    surface_src: Path | None = None  # source space holding the full white surface; default bem/{name}-all-src.fif
    fiducials: dict | None = None  # LPA, nasion, RPA [m], head frame (None: those of the sample recording)
    parent: "Subject | None" = None  # the subject a scaled control was made from
    scale: float = 1.0  # linear scale factor relative to ``parent``
    description: str = ""

    @property
    def labels(self) -> Path:
        return self.label_dir if self.label_dir is not None else self.subjects_dir / self.name / "label"

    def digitisation(self) -> mne.Info:
        """Measurement info carrying only the head-frame fiducials (what the OPM coverage rules
        use); for the sample subject, the digitisation of its recording."""
        if self.fiducials is None:
            return mne.io.read_info(paths.SAMPLE_MEG / "sample_audvis_raw.fif", verbose=False)
        return fiducial_info(self.fiducials)

    @property
    def src_rr(self) -> np.ndarray:
        return np.concatenate([s["rr"][s["vertno"]] for s in self.src])

    @property
    def src_nn(self) -> np.ndarray:
        return np.concatenate([s["nn"][s["vertno"]] for s in self.src])

    def bem_model(self, conductivity: tuple | None = None, head_refine: int | None = None) -> list:
        """BEM surfaces for ``mne.make_bem_solution``. One value -> inner skull only (standard
        MEG model); three values -> scalp, skull, brain conductivities as given (S/m). In a
        3-layer model the head surface is subdivided ``head_refine`` times (default
        ``HEAD_SURFACE_REFINE``; A-BEM-SKIN): on-scalp sensors sit a few mm from it, where the
        field of the 5,120-triangle surface is not converged (scripts/study_opm_near_mesh.py)."""
        surfs = [dict(s) for s in self.bem_surfaces]
        if conductivity is None or len(conductivity) == 1:
            inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
            inner["sigma"] = 0.3 if conductivity is None else float(conductivity[0])
            return [inner]
        if len(conductivity) != 3:
            raise ValueError("conductivity must have 1 or 3 values (scalp, skull, brain)")
        order = [FIFF.FIFFV_BEM_SURF_ID_HEAD, FIFF.FIFFV_BEM_SURF_ID_SKULL, FIFF.FIFFV_BEM_SURF_ID_BRAIN]
        refine = HEAD_SURFACE_REFINE if head_refine is None else int(head_refine)
        out = []
        for sid, sigma in zip(order, conductivity):
            s = next(s for s in surfs if s["id"] == sid)
            if sid == FIFF.FIFFV_BEM_SURF_ID_HEAD and refine > 0:
                s = dict(_refined_surface(self.name, s, refine))
            s["sigma"] = float(sigma)
            out.append(s)
        return out


HEAD_SURFACE_REFINE = 1  # A-BEM-SKIN: 3-layer head surface 5,120 -> 20,480 triangles (same flat geometry)
_REFINED: dict = {}


def _refined_surface(name: str, surf: dict, times: int) -> dict:
    """A BEM surface subdivided ``times`` by flat midpoints (same shape, 4^times triangles), cached
    by subject name, surface id, refinement and a hash of the geometry (a scaled control shares
    the name and mesh size of nothing else, but the hash makes the key safe regardless)."""
    import hashlib

    from mne.surface import complete_surface_info

    geometry = hashlib.sha1(np.ascontiguousarray(surf["rr"], np.float64).tobytes()
                            + np.ascontiguousarray(surf["tris"], np.int64).tobytes()).hexdigest()
    key = (name, int(surf["id"]), times, geometry)
    if key not in _REFINED:
        fine = Surface(surf["rr"], surf["tris"], surf["nn"]).subdivided(times)
        new = dict(id=surf["id"], sigma=surf.get("sigma", 1.0), coord_frame=surf["coord_frame"], rr=fine.rr.copy(),
                   tris=fine.tris.copy(), np=len(fine.rr), ntri=len(fine.tris))
        _REFINED[key] = complete_surface_info(new, copy=False, verbose=False)
    return _REFINED[key]


def refined_inner_skull(subject: "Subject", times: int = 1, sigma: float = 0.3) -> list:
    """Single-layer BEM surface list with the inner skull subdivided ``times`` (same shape,
    4^times more triangles): a discretisation check of the 5120-triangle model."""
    from mne.surface import complete_surface_info

    inner = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    fine = Surface(inner["rr"], inner["tris"], inner["nn"]).subdivided(times)
    surf = dict(id=FIFF.FIFFV_BEM_SURF_ID_BRAIN, sigma=float(sigma), coord_frame=inner["coord_frame"], rr=fine.rr,
                tris=fine.tris, np=len(fine.rr), ntri=len(fine.tris))
    return [complete_surface_info(surf, copy=False, verbose=False)]


def load_sample(spacing: str = "oct6") -> Subject:
    """The MNE sample subject. ``spacing``: 'oct6' (stored source space, 2 x 4098) or 'all'."""
    sd = paths.require(paths.SUBJECTS_DIR, "MNE sample subjects directory")
    bem_dir = sd / "sample" / "bem"
    src_file = {"oct6": "sample-oct-6-src.fif", "all": "sample-all-src.fif"}[spacing]
    src = mne.read_source_spaces(bem_dir / src_file, verbose=False)
    surfs = mne.read_bem_surfaces(bem_dir / "sample-5120-5120-5120-bem.fif", verbose=False)
    head = mne.read_bem_surfaces(bem_dir / "sample-head.fif", verbose=False)[0]
    inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    trans = mne.read_trans(paths.require(paths.SAMPLE_MEG / "sample_audvis_raw-trans.fif", "sample trans"))
    return Subject("sample", sd, src, surfs, _outward(head), _outward(inner), trans,
                   description="MNE sample subject (adult; individual MRI)")


def fiducial_info(fiducials: dict) -> mne.Info:
    """An empty measurement info whose digitisation is LPA, nasion and RPA (head frame [m])."""
    info = mne.create_info(["dummy"], 1000.0, "misc")
    montage = mne.channels.make_dig_montage(nasion=np.asarray(fiducials["nasion"], float), lpa=np.asarray(fiducials["lpa"], float),
                                            rpa=np.asarray(fiducials["rpa"], float), coord_frame="head")
    info.set_montage(montage, on_missing="ignore")
    return info


INFANT_SUBJECTS = "infant_subjects"  # data/external/infant_subjects (mne.datasets.fetch_infant_template)


def load_template(name: str = "ANTS2-0Years3T") -> Subject:
    """An infant template (O'Reilly et al. 2021, via ``mne.datasets.fetch_infant_template``) in its
    native dimensions: 3-layer BEM (5,120 triangles per surface), dense head surface, oct-6 source
    space (whose full white surface also provides the full-resolution cortex), aparc labels. The
    head frame is defined by the template's MRI-frame fiducials (Neuromag convention), so the
    head-to-MRI transform follows from them."""
    sd = paths.require(paths.EXTERNAL / INFANT_SUBJECTS, "infant template directory")
    bem_dir = paths.require(sd / name / "bem", f"{name} BEM directory")
    src_file = bem_dir / f"{name}-oct-6-src.fif"
    src = mne.read_source_spaces(src_file, verbose=False)
    surfs = mne.read_bem_surfaces(bem_dir / f"{name}-5120-5120-5120-bem.fif", verbose=False)
    head = mne.read_bem_surfaces(bem_dir / f"{name}-head.fif", verbose=False)[0]
    inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    fids, frame = mne.io.read_fiducials(bem_dir / f"{name}-fiducials.fif", verbose=False)
    if frame != FIFF.FIFFV_COORD_MRI:
        raise ValueError(f"{name} fiducials are not in the MRI frame")
    key = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}
    mri = {key[f["ident"]]: np.asarray(f["r"], float) for f in fids if f["ident"] in key}
    mri_head = mne.transforms.get_ras_to_neuromag_trans(mri["nasion"], mri["lpa"], mri["rpa"])
    trans = mne.transforms.Transform("head", "mri", np.linalg.inv(mri_head))
    fid_head = {k: mne.transforms.apply_trans(mri_head, v) for k, v in mri.items()}
    return Subject(name, sd, src, surfs, _outward(head), _outward(inner), trans, label_dir=sd / name / "label",
                   surface_src=src_file, fiducials=fid_head,
                   description=f"infant template {name} (O'Reilly et al. 2021; native dimensions)")


def _scaled_surface(surf: dict, s: float) -> dict:
    from mne.surface import complete_surface_info

    new = dict(id=surf["id"], sigma=surf.get("sigma", 1.0), coord_frame=surf["coord_frame"], rr=np.asarray(surf["rr"], float) * s,
               tris=np.asarray(surf["tris"]).copy(), np=len(surf["rr"]), ntri=len(surf["tris"]))
    return complete_surface_info(new, copy=False, verbose=False)


_SRC_LENGTHS = ("rr", "tri_cent", "use_tri_cent", "nearest_dist", "dist_limit")
_SRC_AREAS = ("tri_area", "use_tri_area")


def scaled(subject: Subject, factor: float, name: str | None = None) -> Subject:
    """Size-only control: every MRI-frame coordinate of ``subject`` multiplied by ``factor`` (about
    the MRI origin). Directions, the mesh topology and the vertex correspondence are unchanged, so
    every vertex is homologous to the original; areas scale by factor^2. The head frame (defined by
    the scaled fiducials) is the original one scaled about its origin, so the head-to-MRI transform
    keeps its rotation and scales its translation. Not pediatric anatomy (GOAL G3)."""
    s = float(factor)
    src = subject.src.copy()
    for hemi in src:
        for k in _SRC_LENGTHS:
            if hemi.get(k) is not None:
                hemi[k] = np.asarray(hemi[k], float) * s
        for k in _SRC_AREAS:
            if hemi.get(k) is not None:
                hemi[k] = np.asarray(hemi[k], float) * s**2
        if hemi.get("dist") is not None:
            hemi["dist"] = hemi["dist"] * s
    t = np.array(subject.trans["trans"], float)
    t[:3, 3] *= s
    fids = subject.fiducials if subject.fiducials is not None else _sample_fiducials()
    return Subject(name or f"{subject.name}_x{s:.3f}", subject.subjects_dir, src, [_scaled_surface(b, s) for b in subject.bem_surfaces],
                   Surface(subject.scalp.rr * s, subject.scalp.tris, subject.scalp.nn),
                   Surface(subject.inner_skull.rr * s, subject.inner_skull.tris, subject.inner_skull.nn),
                   mne.transforms.Transform("head", "mri", t), label_dir=subject.labels, surface_src=subject.surface_src,
                   fiducials={k: np.asarray(v, float) * s for k, v in fids.items()}, parent=subject, scale=s,
                   description=f"{subject.name} scaled by {s:.3f} (size-only control)")


def _sample_fiducials() -> dict:
    from .opm import fiducials_head

    return fiducials_head(mne.io.read_info(paths.SAMPLE_MEG / "sample_audvis_raw.fif", verbose=False))


def closest_point_on_triangles(p: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Closest points to ``p`` (3,) on triangles (a, b, c) (each (m, 3)); Ericson, Real-Time
    Collision Detection, 5.1.5 (Voronoi regions of vertices, edges and face)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = np.einsum("ij,ij->i", ab, ap), np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3, d4 = np.einsum("ij,ij->i", ab, bp), np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5, d6 = np.einsum("ij,ij->i", ab, cp), np.einsum("ij,ij->i", ac, cp)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 / (va + vb + vc)
        out = a + ab * (vb * denom)[:, None] + ac * (vc * denom)[:, None]  # inside the face
        rbc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        out = np.where(rbc[:, None], b + w[:, None] * (c - b), out)
        rac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        out = np.where(rac[:, None], a + (d2 / (d2 - d6))[:, None] * ac, out)
        rab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        out = np.where(rab[:, None], a + (d1 / (d1 - d3))[:, None] * ab, out)
    out = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, out)
    out = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, out)
    out = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, out)
    return out


class MeshDistance:
    """Exact point-to-mesh distance (closest point on the triangles around the ``k`` nearest
    vertices) with the sign of the inside test (negative inside the closed surface ``surf``)."""

    def __init__(self, surf: dict, k: int = 10):
        from mne.surface import _CheckInside

        self.rr, self.tris, self.k = np.asarray(surf["rr"], float), np.asarray(surf["tris"]), k
        self.tree = cKDTree(self.rr)
        inc = [[] for _ in range(len(self.rr))]
        for t, tri in enumerate(self.tris):
            for v in tri:
                inc[v].append(t)
        self.incident = [np.array(x, int) for x in inc]
        self.inside = _CheckInside(surf)

    def unsigned(self, points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(points)
        _, nearest = self.tree.query(points, k=self.k)
        out = np.empty(len(points))
        for i, (p, near) in enumerate(zip(points, nearest)):
            t = np.unique(np.concatenate([self.incident[v] for v in near]))
            a, b, c = (self.rr[self.tris[t, j]] for j in range(3))
            out[i] = np.min(np.linalg.norm(closest_point_on_triangles(p, a, b, c) - p, axis=1))
        return out

    def signed(self, points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(points)
        d = self.unsigned(points)
        return np.where(self.inside(points), -d, d)


def depth_to_surface(points: np.ndarray, surface: Surface) -> np.ndarray:
    """Distance [m] from each point to the nearest vertex of a (dense) surface."""
    return cKDTree(surface.rr).query(points)[0]


def orientation_angle(points: np.ndarray, normals: np.ndarray, inner_skull: Surface,
                      radius: float = 0.010) -> np.ndarray:
    """Angle [deg] between each source normal and the local inner-skull normal (average of the
    inner-skull vertex normals within ``radius`` of the nearest inner-skull point). 0 deg =
    perpendicular to the skull (radial), 90 deg = parallel to it (tangential). The sign of
    the source normal is ignored (angles folded into 0-90 deg)."""
    tree = cKDTree(inner_skull.rr)
    _, idx = tree.query(points)
    out = np.empty(len(points))
    for i, (j, n) in enumerate(zip(idx, normals)):
        near = tree.query_ball_point(inner_skull.rr[j], radius)
        m = inner_skull.nn[near].sum(axis=0)
        m /= np.linalg.norm(m)
        out[i] = np.degrees(np.arccos(min(1.0, abs(float(n @ m)))))
    return out


MIN_BEM_DISTANCE = 0.004  # [m] A-BEM-DIST: sources closer to the 5120-triangle inner skull are not used


@dataclass
class FullResCortex:
    """Full-resolution white surface of both hemispheres (sources at every vertex, as in Hunold
    et al. 2016 and Goldenholz et al. 2009). Vertices outside the inner-skull BEM surface are
    marked invalid (MNE cannot place BEM sources there). ``usable`` further drops vertices within
    MIN_BEM_DISTANCE (4 mm) of the 5120-triangle inner-skull mesh, where the linear-collocation
    lead fields are not converged: refining the mesh to 20,480 triangles changes the Neuromag gains
    by a median 14 % (90th percentile 67 %) at 2-3 mm, 2.5 % (13 %) at 3-4 mm and 0.8 % (3.6 %) at
    4-5 mm (OPM about half). This drops 8.7 % of the cortex, mostly gyral crowns nearest the skull.
    Forward matrices keep a column for every valid vertex; sources are chosen among usable ones."""

    rr: np.ndarray  # (n, 3) MRI [m]
    nn: np.ndarray  # (n, 3) unit normals (outward from white matter)
    area: np.ndarray  # (n,) vertex area [m^2] (one third of adjacent triangle areas)
    hemi: np.ndarray  # (n,) 0 = lh, 1 = rh
    vertno: np.ndarray  # (n,) vertex number within its hemisphere
    tris: np.ndarray  # (m, 3) triangles in global indices
    valid: np.ndarray  # (n,) inside the inner skull
    adjacency: object  # scipy.sparse csr (n, n) edge lengths [m]
    dist_inner_skull: np.ndarray | None = None  # (n,) distance to the inner-skull mesh [m]

    @property
    def n(self) -> int:
        return len(self.rr)

    @property
    def usable(self) -> np.ndarray:
        return self.valid & (self.dist_inner_skull >= MIN_BEM_DISTANCE)


def inner_skull_distance(subject: "Subject", points: np.ndarray) -> np.ndarray:
    """Distance [m] from each point to the subject's 5120-triangle inner-skull BEM surface (to the
    surface subdivided three times: nearest-vertex distance, which exceeds the exact point-to-surface
    distance by at most ~0.4 mm)."""
    inner = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    fine = Surface(inner["rr"], inner["tris"], inner["nn"]).subdivided(3)
    return cKDTree(fine.rr).query(points)[0]


def full_resolution(subject: Subject) -> FullResCortex:
    """Build (or load from ``cache/anatomy``) the full-resolution cortex of ``subject``. A scaled
    control is its parent's cortex scaled (same vertices; validity and the distance to the inner
    skull recomputed on the scaled meshes)."""
    import scipy.sparse as sp
    from mne.surface import _CheckInside

    if subject.parent is not None:
        p = full_resolution(subject.parent)
        s = subject.scale
        rr = p.rr * s
        inner = next(b for b in subject.bem_surfaces if b["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
        return FullResCortex(rr, p.nn.copy(), p.area * s**2, p.hemi, p.vertno, p.tris, _CheckInside(inner)(rr),
                             (p.adjacency * s).tocsr(), inner_skull_distance(subject, rr))
    cache = paths.CACHE / "anatomy" / f"{subject.name}_fullres.npz"
    if cache.exists():
        z = np.load(cache)
        adj = sp.csr_matrix((z["adj_data"], z["adj_indices"], z["adj_indptr"]), shape=(len(z["rr"]),) * 2)
        return FullResCortex(z["rr"], z["nn"], z["area"], z["hemi"], z["vertno"], z["tris"], z["valid"], adj,
                             inner_skull_distance(subject, z["rr"]))

    src_file = subject.surface_src or subject.subjects_dir / subject.name / "bem" / f"{subject.name}-all-src.fif"
    src = mne.read_source_spaces(src_file, verbose=False)  # every MNE surface source space holds the full surface
    rr, nn, area, hemi, vertno, tris = [], [], [], [], [], []
    offset = 0
    for h, s in enumerate(src):
        r, t = s["rr"], s["tris"]
        tri_area = 0.5 * np.linalg.norm(np.cross(r[t[:, 1]] - r[t[:, 0]], r[t[:, 2]] - r[t[:, 0]]), axis=1)
        rr.append(r)
        nn.append(s["nn"])
        area.append(np.bincount(t.ravel(), np.repeat(tri_area / 3.0, 3), minlength=len(r)))
        hemi.append(np.full(len(r), h))
        vertno.append(np.arange(len(r)))
        tris.append(t + offset)
        offset += len(r)
    rr, nn, tris = np.concatenate(rr), np.concatenate(nn), np.concatenate(tris)
    inner = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    valid = _CheckInside(inner)(rr)
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    w = np.linalg.norm(rr[e[:, 0]] - rr[e[:, 1]], axis=1)
    adj = sp.coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(len(rr),) * 2)
    out = FullResCortex(rr, nn, np.concatenate(area), np.concatenate(hemi), np.concatenate(vertno), tris, valid,
                        adj.tocsr(), inner_skull_distance(subject, rr))
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, rr=out.rr, nn=out.nn, area=out.area, hemi=out.hemi, vertno=out.vertno, tris=out.tris,
             valid=out.valid, adj_data=out.adjacency.data, adj_indices=out.adjacency.indices,
             adj_indptr=out.adjacency.indptr)
    return out


def vertex_areas(src: mne.SourceSpaces) -> list[np.ndarray]:
    """Cortical area [m^2] represented by each source (sum over its patch of the full-resolution
    surface, from MNE's patch statistics ``pinfo``)."""
    out = []
    for s in src:
        if s.get("pinfo") is None:
            raise ValueError("source space has no patch statistics (use add_source_space_distances)")
        rr, tris = s["rr"], s["tris"]
        tri_area = 0.5 * np.linalg.norm(np.cross(rr[tris[:, 1]] - rr[tris[:, 0]], rr[tris[:, 2]] - rr[tris[:, 0]]),
                                        axis=1)
        v_area = np.bincount(tris.ravel(), np.repeat(tri_area / 3.0, 3), minlength=len(rr))
        out.append(np.array([v_area[s["pinfo"][k]].sum() for k in s["patch_inds"]]))
    return out
