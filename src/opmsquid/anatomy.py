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
    head_conform: dict | None = None  # what ``head_on_scalp`` did to the BEM head surface at loading (A-BEM-CONFORM)

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


def crossing_edges(rr: np.ndarray, tris: np.ndarray) -> int:
    """Number of mesh edges that cross a triangle they do not belong to (0 for a surface without
    self-intersections). Candidates: triangles whose centroid lies within half the edge plus the
    largest triangle circumradius of the edge's midpoint."""
    rr, tris = np.asarray(rr, float), np.asarray(tris)
    tri = rr[tris]
    cen = tri.mean(axis=1)
    r_tri = float(np.max(np.linalg.norm(tri - cen[:, None, :], axis=2)))
    tree = cKDTree(cen)
    edges = np.unique(np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1), axis=0)
    n = 0
    for e in edges:
        a, b = rr[e[0]], rr[e[1]]
        cand = np.asarray(tree.query_ball_point(0.5 * (a + b), 0.5 * np.linalg.norm(b - a) + r_tri), int)
        cand = cand[~np.isin(tris[cand], e).any(axis=1)]
        if len(cand) and _segments_hit_triangles(a[None], b[None], tri[cand]).any():
            n += 1
    return n


def head_on_scalp(surfaces: list, scalp: Surface) -> tuple[list, dict]:
    """A-BEM-CONFORM: the BEM surfaces with every vertex of the head surface moved to the nearest
    vertex of the MRI scalp ``scalp`` (same mesh, same triangles; the other surfaces unchanged).
    The infant templates' head surfaces are built this way (each of their 2,562 vertices is a vertex
    of the template's MRI head surface), so for them this is the identity and the surfaces are
    returned as they are. The sample subject's (its own segmentation's outer skin) lies a median
    0.8 mm outside its MRI scalp over the OPM coverage region (up to 2.6 mm), which made the
    whole-cell OPM clearance (A-OPM-CLEAR) move the adult's sensors farther from the scalp than the
    children's (goal review, 2026-10-02). Refused if two vertices would merge, a triangle would flip,
    the surface would fold (an edge crossing a triangle) or the outer skull would not stay inside.
    Returns (surfaces, report)."""
    from mne.surface import _CheckInside, complete_surface_info

    k = next(i for i, s in enumerate(surfaces) if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    head = surfaces[k]
    old = np.asarray(head["rr"], float)
    tree = cKDTree(scalp.rr)
    _, idx = tree.query(old)
    n_shared = 0
    if len(np.unique(idx)) != len(idx):  # two vertices nearest to one scalp vertex: the farther one takes its next free one
        _, cand = tree.query(old, k=10)
        taken = set()
        for v in np.argsort(np.linalg.norm(np.asarray(scalp.rr)[idx] - old, axis=1)):
            free = [c for c in cand[v] if c not in taken]
            if not free:
                raise ValueError("A-BEM-CONFORM: no free scalp vertex among the 10 nearest")
            n_shared += int(free[0] != idx[v])
            idx[v] = free[0]
            taken.add(free[0])
    new_rr = np.asarray(scalp.rr, float)[idx]
    moved = np.linalg.norm(new_rr - old, axis=1)
    report = dict(n_vertices=len(old), moved_median_mm=float(np.median(moved) * 1e3), moved_p90_mm=float(np.percentile(moved, 90) * 1e3),
                  moved_max_mm=float(moved.max() * 1e3), reassigned=n_shared)
    if not np.any(moved > 0):
        return list(surfaces), dict(report, identity=True)
    if len(np.unique(idx)) != len(idx):
        raise ValueError("A-BEM-CONFORM: two head-surface vertices map to one scalp vertex")
    tris = np.asarray(head["tris"])

    def face_normals(rr):
        return np.cross(rr[tris[:, 1]] - rr[tris[:, 0]], rr[tris[:, 2]] - rr[tris[:, 0]])

    n_repaired = 0
    nbrs = None
    for _ in range(50):  # a vertex whose triangle would flip goes to the free scalp vertex nearest its neighbours' mean
        bad = np.einsum("ij,ij->i", face_normals(old), face_normals(new_rr)) <= 0
        if not bad.any():
            break
        if nbrs is None:
            nbrs = [set() for _ in range(len(old))]
            for a, b, c in tris:
                nbrs[a] |= {b, c}
                nbrs[b] |= {a, c}
                nbrs[c] |= {a, b}
        used = set(idx.tolist())
        for v in np.unique(tris[bad]):
            _, cand = tree.query(new_rr[list(nbrs[v])].mean(axis=0), k=10)
            free = [c for c in cand if c not in used or c == idx[v]]
            if free and free[0] != idx[v]:
                used.discard(idx[v])
                idx[v] = free[0]
                used.add(free[0])
                new_rr[v] = np.asarray(scalp.rr[free[0]], float)
                n_repaired += 1
    report.update(repaired=n_repaired, moved_max_mm=float(np.linalg.norm(new_rr - old, axis=1).max() * 1e3))
    if np.any(np.einsum("ij,ij->i", face_normals(old), face_normals(new_rr)) <= 0):
        raise ValueError("A-BEM-CONFORM: a head-surface triangle would flip")
    if crossing_edges(new_rr, tris):
        raise ValueError("A-BEM-CONFORM: the head surface would fold (an edge crosses a triangle)")
    new = complete_surface_info(dict(id=head["id"], sigma=head.get("sigma", 1.0), coord_frame=head["coord_frame"], rr=new_rr,
                                     tris=tris.copy(), np=len(new_rr), ntri=len(tris)), copy=False, verbose=False)
    skull = next(s for s in surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_SKULL)
    if not _CheckInside(new)(skull["rr"], verbose=False).all():
        raise ValueError("A-BEM-CONFORM: the outer skull would leave the head surface")
    report.update(identity=False, outer_skull_min_mm=float(MeshDistance(new).unsigned(skull["rr"]).min() * 1e3))
    out = list(surfaces)
    out[k] = new
    return out, report


def model_skull(inner: dict, scalp: Surface, cortex_rr: np.ndarray, depth: float = 0.008, min_cortex: float = 0.002,
                smooth: int = 10, max_skull: float = 0.005) -> tuple[dict, dict, dict]:
    """A-BEM-CHILD: an inner skull that lies too close to the scalp (a failed watershed segmentation,
    as in the school-aged children) moved inward along its normals to ``depth`` [m] below the MRI
    scalp wherever it is shallower, no vertex closer than ``min_cortex`` to a cortex vertex
    (``cortex_rr``, white-surface vertices; the surface between its vertices can come closer: 0.9-1.3
    mm in the children, scripts/study_school_anatomy.py); the displacement is smoothed over the mesh
    (``smooth`` neighbour averages) and the cortex limit applied again. Parts already deeper than
    ``depth`` (the skull base) stay where they are. The outer skull is the modelled inner skull moved
    outward by half its local distance to the scalp (at least 1 mm, leaving 2 mm of scalp where that
    distance allows; at most ``max_skull``, where the scalp is far: the skull base). Both keep the
    input's triangles. Refused if a triangle flips, the surface folds, a vertex of an even sample of
    the cortex (about 20,000) is left outside or the outer skull leaves the scalp; the report's cortex
    clearance is over the same sample. Returns (inner, outer skull, report); surfaces in the MRI frame
    [m]."""
    import scipy.sparse as sp
    from mne.surface import _CheckInside, complete_surface_info

    tris = np.asarray(inner["tris"])
    rr = np.asarray(inner["rr"], float)
    scalp_surf = dict(rr=scalp.rr, tris=scalp.tris, nn=scalp.nn, np=len(scalp.rr), ntri=len(scalp.tris))
    to_scalp = MeshDistance(scalp_surf)
    cortex = cKDTree(np.asarray(cortex_rr, float))
    d_s = to_scalp.unsigned(rr)
    d_c = cortex.query(rr)[0]
    limit = np.maximum(d_c - min_cortex, 0.0)
    delta = np.minimum(np.clip(depth - d_s, 0.0, None), limit)
    e = np.unique(np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1), axis=0)
    adj = sp.coo_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(len(rr),) * 2).tocsr()
    deg = np.asarray(adj.sum(axis=1)).ravel()

    def smoothed_normals(x):  # vertex normals averaged over ``smooth`` neighbour rings: offsets along them do not cross
        n = _outward(complete_surface_info(dict(rr=x, tris=tris, np=len(x), ntri=len(tris)), copy=False, verbose=False)).nn
        for _ in range(int(smooth)):
            n = 0.5 * n + 0.5 * (adj @ n) / deg[:, None]
            n /= np.linalg.norm(n, axis=1, keepdims=True)
        return n

    nn = smoothed_normals(rr)
    for _ in range(int(smooth)):
        delta = np.minimum(0.5 * delta + 0.5 * (adj @ delta) / deg, limit)
    new_rr = rr - delta[:, None] * nn

    def face_normals(x):
        return np.cross(x[tris[:, 1]] - x[tris[:, 0]], x[tris[:, 2]] - x[tris[:, 0]])

    def check(x, ref, what):
        flips = int(np.sum(np.einsum("ij,ij->i", face_normals(ref), face_normals(x)) <= 0))
        if flips:
            raise ValueError(f"A-BEM-CHILD: {flips} {what} triangles would flip")
        if crossing_edges(x, tris):
            raise ValueError(f"A-BEM-CHILD: the {what} would fold (an edge crosses a triangle)")

    check(new_rr, rr, "inner-skull")
    new_inner = complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_BRAIN, coord_frame=inner.get("coord_frame", FIFF.FIFFV_COORD_MRI),
                                           rr=new_rr, tris=tris.copy(), np=len(new_rr), ntri=len(tris)), copy=False, verbose=False)
    sample = np.asarray(cortex_rr, float)[:: max(1, len(cortex_rr) // 20000)]
    if not _CheckInside(new_inner)(sample, verbose=False).all():
        raise ValueError("A-BEM-CHILD: cortex outside the modelled inner skull")
    nn_new = smoothed_normals(new_rr)
    gap = to_scalp.unsigned(new_rr)
    lo, hi = 0.001, np.clip(gap - 0.002, 0.001, max_skull)
    t = np.clip(0.5 * gap, lo, hi)
    for _ in range(int(smooth)):
        t = np.clip(0.5 * t + 0.5 * (adj @ t) / deg, lo, hi)
    os_rr = new_rr + t[:, None] * nn_new
    for _ in range(30):  # where an offset triangle would flip (a crease of the modelled inner skull), thin the skull there
        bad = np.einsum("ij,ij->i", face_normals(new_rr), face_normals(os_rr)) <= 0
        if not bad.any():
            break
        v = np.unique(tris[bad])
        t[v] = np.maximum(0.7 * t[v], lo)
        os_rr = new_rr + t[:, None] * nn_new
    check(os_rr, new_rr, "outer-skull")
    outer = complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_SKULL, coord_frame=new_inner["coord_frame"], rr=os_rr, tris=tris.copy(),
                                       np=len(os_rr), ntri=len(tris)), copy=False, verbose=False)
    if not to_scalp.inside(os_rr, verbose=False).all():
        raise ValueError("A-BEM-CHILD: the outer skull would leave the scalp")
    d_new = to_scalp.unsigned(new_rr)
    report = dict(depth_mm=depth * 1e3, min_cortex_mm=min_cortex * 1e3, smooth=int(smooth), n_vertices=len(rr),
                  moved_share=float(np.mean(delta > 1e-4)), moved_median_mm=float(np.median(delta[delta > 1e-4]) * 1e3) if np.any(delta > 1e-4) else 0.0,
                  moved_max_mm=float(delta.max() * 1e3),
                  scalp_depth_before_mm=[float(x) for x in np.percentile(d_s * 1e3, [10, 50, 90])],
                  scalp_depth_after_mm=[float(x) for x in np.percentile(d_new * 1e3, [10, 50, 90])],
                  cortex_clearance_min_mm=float(MeshDistance(new_inner).unsigned(sample).min() * 1e3),
                  skull_thickness_mm=[float(x) for x in np.percentile(t * 1e3, [10, 50, 90])])
    return new_inner, outer, report


def similarity_fit(src: np.ndarray, tgt: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Least-squares similarity tgt ~ s R src + t of corresponding points (Umeyama 1991)."""
    mu_s, mu_t = src.mean(axis=0), tgt.mean(axis=0)
    a, b = src - mu_s, tgt - mu_t
    u, d, vt = np.linalg.svd(b.T @ a / len(src))
    e = np.ones(3)
    if np.linalg.det(u) * np.linalg.det(vt) < 0:
        e[-1] = -1.0
    r = u @ np.diag(e) @ vt
    s = float(np.sum(d * e) / np.mean(np.sum(a**2, axis=1)))
    return s, r, mu_t - s * r @ mu_s


def transfer_fiducials(adult_cortex: np.ndarray, adult_fids: dict, cortex: np.ndarray, scalp: np.ndarray, n_iter: int = 60,
                       trim: float = 90.0) -> tuple[dict, dict]:
    """A-G3-FID: fiducials of a head that has none of its own: the adult's (LPA, nasion, RPA; MRI
    frame [m]) mapped by a similarity transform fitted from the adult's cortex to ``cortex`` (white
    surface vertices) by trimmed iterative closest points (start: centroids and sizes aligned; the
    ``trim`` percent closest pairs refit each time), then moved to the nearest vertex of ``scalp``.
    The cortex is fitted rather than the scalp because faces scale differently from crania (on the
    infant templates the same fit between the scalps puts the nasion 34-45 mm from theirs; this one
    reproduces their fiducials to 1.6-14.4 mm and their head frame to 3.2-5.9 deg:
    scripts/study_school_anatomy.py). Returns (fiducials, report)."""
    src = np.asarray(adult_cortex, float)
    src = src[:: max(1, len(src) // 8000)]
    tgt = np.asarray(cortex, float)
    tgt = tgt[:: max(1, len(tgt) // 8000)]
    size = lambda x: np.linalg.norm(x.max(axis=0) - x.min(axis=0))  # noqa: E731
    s = size(tgt) / size(src)
    r = np.eye(3)
    t = tgt.mean(axis=0) - s * src.mean(axis=0)
    tree = cKDTree(tgt)
    for _ in range(int(n_iter)):
        d, idx = tree.query(s * src @ r.T + t)
        keep = d <= np.percentile(d, trim)
        s, r, t = similarity_fit(src[keep], tgt[idx[keep]])
    d, _ = tree.query(s * src @ r.T + t)
    mapped = {k: s * r @ np.asarray(v, float) + t for k, v in adult_fids.items()}
    stree = cKDTree(scalp)
    out = {k: np.asarray(scalp[stree.query(v)[1]], float) for k, v in mapped.items()}
    report = dict(scale=float(s), rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1)))),
                  cortex_fit_median_mm=float(np.median(d) * 1e3), projection_mm={k: float(np.linalg.norm(out[k] - mapped[k]) * 1e3) for k in out})
    return out, report


def load_sample(spacing: str = "oct6") -> Subject:
    """The MNE sample subject, its BEM head surface on its MRI scalp (``head_on_scalp``).
    ``spacing``: 'oct6' (stored source space, 2 x 4098) or 'all'."""
    sd = paths.require(paths.SUBJECTS_DIR, "MNE sample subjects directory")
    bem_dir = sd / "sample" / "bem"
    src_file = {"oct6": "sample-oct-6-src.fif", "all": "sample-all-src.fif"}[spacing]
    src = mne.read_source_spaces(bem_dir / src_file, verbose=False)
    surfs = mne.read_bem_surfaces(bem_dir / "sample-5120-5120-5120-bem.fif", verbose=False)
    head = mne.read_bem_surfaces(bem_dir / "sample-head.fif", verbose=False)[0]
    inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    trans = mne.read_trans(paths.require(paths.SAMPLE_MEG / "sample_audvis_raw-trans.fif", "sample trans"))
    scalp = _outward(head)
    surfs, conform = head_on_scalp(surfs, scalp)
    return Subject("sample", sd, src, surfs, scalp, _outward(inner), trans,
                   description="MNE sample subject (adult; individual MRI)", head_conform=conform)


def fiducial_info(fiducials: dict) -> mne.Info:
    """An empty measurement info whose digitisation is LPA, nasion and RPA (head frame [m])."""
    info = mne.create_info(["dummy"], 1000.0, "misc")
    montage = mne.channels.make_dig_montage(nasion=np.asarray(fiducials["nasion"], float), lpa=np.asarray(fiducials["lpa"], float),
                                            rpa=np.asarray(fiducials["rpa"], float), coord_frame="head")
    info.set_montage(montage, on_missing="ignore")
    return info


INFANT_SUBJECTS = "infant_subjects"  # data/external/infant_subjects (mne.datasets.fetch_infant_template)
SCHOOL_SUBJECTS = "school_subjects"  # data/external/school_subjects (scripts/prepare_school_subjects.py)


def _load_prepared(root: Path, name: str, description: str) -> Subject:
    """A subject laid out as MNE packages the infant templates: bem/{name}-5120-5120-5120-bem.fif,
    -head.fif (dense MRI scalp), -oct-6-src.fif, -fiducials.fif (MRI frame), label/*.annot. The head
    frame follows from the fiducials; the BEM head surface is put on the scalp (A-BEM-CONFORM)."""
    bem_dir = paths.require(root / name / "bem", f"{name} BEM directory")
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
    scalp = _outward(head)
    surfs, conform = head_on_scalp(surfs, scalp)
    return Subject(name, root, src, surfs, scalp, _outward(inner), trans, label_dir=root / name / "label",
                   surface_src=src_file, fiducials=fid_head, description=description, head_conform=conform)


def load_template(name: str = "ANTS2-0Years3T") -> Subject:
    """An infant template (O'Reilly et al. 2021, via ``mne.datasets.fetch_infant_template``) in its
    native dimensions: 3-layer BEM (5,120 triangles per surface), dense head surface, oct-6 source
    space (whose full white surface also provides the full-resolution cortex), aparc labels. The
    head frame is defined by the template's MRI-frame fiducials (Neuromag convention), so the
    head-to-MRI transform follows from them. Its head surface is already on its scalp (identity)."""
    sd = paths.require(paths.EXTERNAL / INFANT_SUBJECTS, "infant template directory")
    return _load_prepared(sd, name, f"infant template {name} (O'Reilly et al. 2021; native dimensions)")


def load_school(name: str) -> Subject:
    """A school-aged child of OpenNeuro ds005234 (Fadeev et al. 2024) prepared by
    ``scripts/prepare_school_subjects.py``: its own white surfaces, aparc labels and dense MRI scalp;
    its watershed outer skin as the BEM head surface (A-BEM-CONFORM: put on the scalp at loading);
    a modelled skull (A-BEM-CHILD) and fiducials transferred from the adult (A-G3-FID)."""
    sd = paths.require(paths.EXTERNAL / SCHOOL_SUBJECTS, "school-aged subjects directory")
    return _load_prepared(sd, name, f"school-aged child {name} (OpenNeuro ds005234; individual MRI)")


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
    """Closest points to ``p`` on triangles (a, b, c): ``p`` (3,) with triangles (m, 3), or any
    shapes that broadcast, e.g. points (n, 1, 3) with triangles (n, m, 3); Ericson, Real-Time
    Collision Detection, 5.1.5 (Voronoi regions of vertices, edges and face)."""
    def dot(x, y):
        return np.einsum("...j,...j->...", x, y)

    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 / (va + vb + vc)
        out = a + ab * (vb * denom)[..., None] + ac * (vc * denom)[..., None]  # inside the face
        rbc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        out = np.where(rbc[..., None], b + w[..., None] * (c - b), out)
        rac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        out = np.where(rac[..., None], a + (d2 / (d2 - d6))[..., None] * ac, out)
        rab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        out = np.where(rab[..., None], a + (d1 / (d1 - d3))[..., None] * ab, out)
    out = np.where(((d6 >= 0) & (d5 <= d6))[..., None], c, out)
    out = np.where(((d3 >= 0) & (d4 <= d3))[..., None], b, out)
    out = np.where(((d1 <= 0) & (d2 <= 0))[..., None], a, out)
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
        # the same as a table, each row padded with a repeat of its first triangle (or triangle 0 for
        # an unused vertex): repeats leave the minimum distance unchanged
        width = max(len(x) for x in inc)
        self.incident_table = np.array([x + [x[0] if x else 0] * (width - len(x)) for x in inc], int)
        self.inside = _CheckInside(surf)

    def unsigned(self, points: np.ndarray, chunk: int = 2000) -> np.ndarray:
        points = np.atleast_2d(points)
        _, nearest = self.tree.query(points, k=self.k)
        out = np.empty(len(points))
        for lo in range(0, len(points), chunk):
            p = points[lo:lo + chunk, None, :]
            t = self.incident_table[nearest[lo:lo + chunk]].reshape(len(p), -1)  # triangles around the k nearest vertices
            a, b, c = (self.rr[self.tris[t, j]] for j in range(3))
            out[lo:lo + chunk] = np.linalg.norm(closest_point_on_triangles(p, a, b, c) - p, axis=-1).min(axis=1)
        return out

    def signed(self, points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(points)
        d = self.unsigned(points)
        return np.where(self.inside(points), -d, d)


def segment_segment_distance(p1, q1, p2, q2) -> np.ndarray:
    """Distances between segments p1-q1 and p2-q2 (arrays (..., 3) that broadcast); Ericson,
    Real-Time Collision Detection, 5.1.9."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2

    def dot(x, y):
        return np.einsum("...j,...j->...", x, y)

    a, e, f, c, b = dot(d1, d1), dot(d2, d2), dot(d2, r), dot(d1, r), dot(d1, d2)
    denom = a * e - b * b
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(denom > 1e-30, np.clip((b * f - c * e) / denom, 0.0, 1.0), 0.0)
        t = (b * s + f) / e
        s = np.where(t < 0.0, np.clip(-c / a, 0.0, 1.0), np.where(t > 1.0, np.clip((b - c) / a, 0.0, 1.0), s))
        t = np.clip(t, 0.0, 1.0)
    return np.linalg.norm(p1 + s[..., None] * d1 - (p2 + t[..., None] * d2), axis=-1)


CUBE_CORNERS = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float)
CUBE_EDGES = np.array([(a, b) for a in range(8) for b in range(a + 1, 8) if np.sum(CUBE_CORNERS[a] != CUBE_CORNERS[b]) == 1])


class CubeMeshDistance:
    """Exact distance between a cube and a triangle mesh. Two convex sets that do not intersect
    have their closest points at a vertex of one and a face of the other or on an edge of each, so
    the distance is the minimum over cube corners to triangles, triangle vertices to the cube and
    cube edges to triangle edges; an intersection (a triangle edge through the cube or a cube edge
    through a triangle) gives 0. Candidate triangles are all those whose centroid lies within
    ``bound`` + the cube's circumradius + the largest triangle circumradius of the cube centre,
    ``bound`` being an upper bound on the distance (e.g. that of a point of the cube): any closer
    triangle lies within that radius, so the search is certified."""

    def __init__(self, surf: dict):
        self.rr = np.asarray(surf["rr"], float)
        self.tris = np.asarray(surf["tris"], int)
        tri = self.rr[self.tris]
        self.centroid = tri.mean(axis=1)
        self.r_tri = float(np.max(np.linalg.norm(tri - self.centroid[:, None, :], axis=2)))
        self.tree = cKDTree(self.centroid)

    def distance(self, centre, frame, half: float, bound: float) -> float:
        """``frame``: rows are the cube's unit axes; ``half``: half its edge [m]."""
        centre, frame = np.asarray(centre, float), np.asarray(frame, float)
        cand = self.tree.query_ball_point(centre, bound + half * np.sqrt(3) + self.r_tri)
        if not cand:
            return float("inf")
        tri = self.rr[self.tris[np.asarray(cand)]]  # (m, 3, 3)
        corners = centre + half * CUBE_CORNERS @ frame
        t0, t1 = tri.reshape(-1, 3), tri[:, [1, 2, 0]].reshape(-1, 3)  # triangle edges
        if _segments_hit_box((t0 - centre) @ frame.T, (t1 - centre) @ frame.T, half).any() or \
                _segments_hit_triangles(corners[CUBE_EDGES[:, 0]], corners[CUBE_EDGES[:, 1]], tri).any():
            return 0.0
        cp = closest_point_on_triangles(corners[:, None, :], tri[None, :, 0], tri[None, :, 1], tri[None, :, 2])
        d_corner = np.linalg.norm(cp - corners[:, None, :], axis=-1).min()
        q = (tri.reshape(-1, 3) - centre) @ frame.T
        d_vertex = np.linalg.norm(np.maximum(np.abs(q) - half, 0.0), axis=1).min()
        d_edge = segment_segment_distance(corners[CUBE_EDGES[:, 0]][:, None], corners[CUBE_EDGES[:, 1]][:, None],
                                          t0[None], t1[None]).min()
        return float(min(d_corner, d_vertex, d_edge))


def _segments_hit_box(p, q, half) -> np.ndarray:
    """Whether segments p-q (cube frame, (n, 3)) meet the box [-half, half]^3 (slab method)."""
    d = q - p
    with np.errstate(divide="ignore", invalid="ignore"):
        t1, t2 = (-half - p) / d, (half - p) / d
    lo, hi = np.minimum(t1, t2), np.maximum(t1, t2)
    flat = d == 0  # parallel to a slab: inside it or never
    lo = np.where(flat, np.where(np.abs(p) <= half, -np.inf, np.inf), lo)
    hi = np.where(flat, np.where(np.abs(p) <= half, np.inf, -np.inf), hi)
    return np.maximum(lo.max(axis=1), 0.0) <= np.minimum(hi.min(axis=1), 1.0)


def _segments_hit_triangles(a, b, tri) -> np.ndarray:
    """Whether any of the segments a-b ((k, 3)) crosses any triangle ((m, 3, 3)); Moller-Trumbore."""
    d = (b - a)[:, None, :]
    e1, e2 = tri[None, :, 1] - tri[None, :, 0], tri[None, :, 2] - tri[None, :, 0]
    h = np.cross(d, e2)
    det = np.einsum("...j,...j->...", e1, h)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / det
        s = a[:, None, :] - tri[None, :, 0]
        u = inv * np.einsum("...j,...j->...", s, h)
        qv = np.cross(s, e1)
        v = inv * np.einsum("...j,...j->...", d, qv)
        t = inv * np.einsum("...j,...j->...", e2, qv)
        return (np.abs(det) > 1e-30) & (u >= 0) & (v >= 0) & (u + v <= 1) & (t >= 0) & (t <= 1)


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
