"""Subject anatomy: source spaces, BEM surfaces, scalp, transforms and source descriptors.

Primary adult anatomy: the MNE ``sample`` subject (individual adult MRI with FreeSurfer surfaces,
3-layer BEM surfaces, oct-6 source space with distances and patch statistics, dense scalp, and
the head-MRI transform of its Neuromag recording). Coordinates are in the MRI frame [m] unless
stated otherwise.

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

    @property
    def src_rr(self) -> np.ndarray:
        return np.concatenate([s["rr"][s["vertno"]] for s in self.src])

    @property
    def src_nn(self) -> np.ndarray:
        return np.concatenate([s["nn"][s["vertno"]] for s in self.src])

    def bem_model(self, conductivity: tuple | None = None) -> list:
        """BEM surfaces for ``mne.make_bem_solution``. One value -> inner skull only (standard
        MEG model); three values -> scalp, skull, brain conductivities as given (S/m)."""
        surfs = [dict(s) for s in self.bem_surfaces]
        if conductivity is None or len(conductivity) == 1:
            inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
            inner["sigma"] = 0.3 if conductivity is None else float(conductivity[0])
            return [inner]
        if len(conductivity) != 3:
            raise ValueError("conductivity must have 1 or 3 values (scalp, skull, brain)")
        order = [FIFF.FIFFV_BEM_SURF_ID_HEAD, FIFF.FIFFV_BEM_SURF_ID_SKULL, FIFF.FIFFV_BEM_SURF_ID_BRAIN]
        out = []
        for sid, sigma in zip(order, conductivity):
            s = next(s for s in surfs if s["id"] == sid)
            s["sigma"] = float(sigma)
            out.append(s)
        return out


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
    return Subject("sample", sd, src, surfs, _outward(head), _outward(inner), trans)


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


MIN_BEM_DISTANCE = 0.002  # [m] A-BEM-DIST: sources closer to the 5120-triangle inner skull are not used


@dataclass
class FullResCortex:
    """Full-resolution white surface of both hemispheres (sources at every vertex, as in Hunold
    et al. 2016 and Goldenholz et al. 2009). Vertices outside the inner-skull BEM surface are
    marked invalid (MNE cannot place BEM sources there). ``usable`` further drops vertices within
    MIN_BEM_DISTANCE of the inner-skull mesh: there the linear-collocation BEM lead fields are
    numerical artefacts (up to ~4000x the energy of neighbouring vertices within 0.5 mm; no
    anomaly beyond 2 mm). Forward matrices keep a column for every valid vertex; sources are
    chosen among usable vertices."""

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
    surface subdivided three times, ~1-mm vertex spacing, so the error is below ~0.1 mm)."""
    inner = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    fine = Surface(inner["rr"], inner["tris"], inner["nn"]).subdivided(3)
    return cKDTree(fine.rr).query(points)[0]


def full_resolution(subject: Subject) -> FullResCortex:
    """Build (or load from ``cache/anatomy``) the full-resolution cortex of ``subject``."""
    import scipy.sparse as sp
    from mne.surface import _CheckInside

    cache = paths.CACHE / "anatomy" / f"{subject.name}_fullres.npz"
    if cache.exists():
        z = np.load(cache)
        adj = sp.csr_matrix((z["adj_data"], z["adj_indices"], z["adj_indptr"]), shape=(len(z["rr"]),) * 2)
        return FullResCortex(z["rr"], z["nn"], z["area"], z["hemi"], z["vertno"], z["tris"], z["valid"], adj,
                             inner_skull_distance(subject, z["rr"]))

    src = mne.read_source_spaces(subject.subjects_dir / subject.name / "bem" / f"{subject.name}-all-src.fif",
                                 verbose=False)
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
