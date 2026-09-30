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
* Sites are defined on the scalp and oriented along a smoothed scalp normal (helmet-scale
  smoothing radius, default 10 mm), which is also the sensitive axis.

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


def coil_def_file(cell_size: float = CELL_SIZE) -> Path:
    """Path of an extra coil-definition file holding the study OPM coil (for use_coil_def)."""
    extra = _coil_block(OPM_COIL_TYPE, cell_size, f"Study OPM, {cell_size * 1e3:.1f}-mm cubic cell (assumption)")
    text = "# opmsquid study OPM coil (src/opmsquid/opm.py); added to MNE's coil_def.dat\n" + extra
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


def make_info(array: OPMArray, sfreq: float = 1000.0) -> mne.Info:
    """MNE info for a single-axis OPM array; device frame = head frame (head-mounted array)."""
    names = [f"OPM{i:03d}" for i in range(len(array.pos))]
    info = mne.create_info(names, sfreq, "mag")
    with info._unlock():
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")
        for ch, p, n in zip(info["chs"], array.pos, array.axis):
            a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ex = np.cross(a, n)
            ex /= np.linalg.norm(ex)
            ch["loc"][:3], ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = p, ex, np.cross(n, ex), n
            ch["coil_type"] = OPM_COIL_TYPE
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
            ch["unit"] = FIFF.FIFF_UNIT_T
    return info


def min_spacing(pos: np.ndarray) -> np.ndarray:
    """Distance from each sensing centre to its nearest neighbour [m]."""
    d, _ = cKDTree(pos).query(pos, k=2)
    return d[:, 1]
