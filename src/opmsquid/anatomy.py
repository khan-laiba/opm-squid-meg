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
