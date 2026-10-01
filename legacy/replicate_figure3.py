#!/usr/bin/env python3
"""
End-to-end replication of Figure 3 of

    Jas M, Matsubara T, Sohrabpour A, Sundaram P, Mody M, Ahlfors SP (2026).
    "Signal-to-noise ratio of event-related fields in on-scalp and off-scalp MEG."
    bioRxiv 2026.08.17.744953.

Figure 3: simulated signal strength (A) and SNR (B) of on-scalp (OPM, xi = 0 mm) and
off-scalp (SQUID, xi = 18 mm) magnetometers measuring the radial magnetic field of a
tangential current dipole (Q = 30 nAm) at different depths d in a spherical head model
(head radius h = 95 mm, brain radius b = 80 mm), relative noise level
eta = sigma_OPM / sigma_SQUID = 3.  The dotted line marks the equal-SNR depth d_eq.

Pipeline
--------
1. Forward simulation with MNE-Python.  A ``mne.Dipole`` (one y-oriented 30 nAm dipole
   per depth, placed on the z-axis at r_Q = h - d) is converted with
   ``mne.make_forward_dipole`` in a spherical conductor (``mne.make_sphere_model``) and
   projected to sensor space with ``mne.simulation.simulate_evoked``.  The sensors are
   radial point magnetometers (``FIFFV_COIL_POINT_MAGNETOMETER``; the paper ignores the
   finite coil / vapour-cell size) densely covering the great circle of the xz-plane at
   radius s = h (OPM, on-scalp) and s = h + 18 mm (SQUID, Triux neo).  The signal is the
   peak of the radial-field topography, B_r^max(r_Q, s) (paper Eq. 1).
2. Checks: the MNE field equals the Biot-Savart radial field at every sensor, and its
   interpolated peak equals the closed-form Eq. (1).
3. SNR (Eqs. 2a,b): SNR_OPM = B_r^max(r_Q, h) / sigma_OPM and
   SNR_SQUID = B_r^max(r_Q, h + xi) / sigma_SQUID, with sigma_OPM = eta * sigma_SQUID.
4. Equal-SNR depth (Eq. 3): B_r^max(h - d_eq, h) / B_r^max(h - d_eq, h + xi) = eta; a
   crossing inside the brain exists for (113/95)^3 = 1.6829 < eta < 5.3086.
5. Figure 3, drawn with the geometry, fonts and annotations of the published figure.

Details recovered from the published figure (not stated in the text)
-------------------------------------------------------------------
The preprint embeds Fig. 3 as a raster whose strokes (0.8 pt = 13.3 px, 1.5 pt = 25 px)
and dash periods give 1200 px per figure inch; it was measured with ``verify_figure3.py``.

* Absolute noise level.  Only eta is given in the paper.  The SNR axis of Fig. 3B is
  reproduced when the SQUID noise equals the SQUID signal of the same 30 nAm dipole at
  r_Q = 0.8 b = 64 mm (d = 31 mm), i.e. the superficial reference source used as "noise" in
  Fig. 6 -- equivalently, SNR_SQUID = 1 for a source at d = 31 mm:
  sigma_SQUID = B_r^max(0.8 b, h + xi) = 0.3546 pT, sigma_OPM = eta * sigma_SQUID = 1.064 pT.
  The tick scales of Fig. 4B-D agree with this value (to ~0.1 %); the tick-label spacing of
  Fig. 3B suggests 0.3547-0.3548 pT (0.05 %, about 1 px at 1200 ppi).  Panels are autoscaled,
  so sigma only sets the SNR tick values, not the pixel shape of the curves.
* Curves: r_Q = linspace(0, b, 101) (0.8-mm steps; the grid of Fig. 4, whose d_eq markers at
  35.0 and 19.8 mm and staircase in Fig. 4E follow it).  Eq. (1) is undefined at r_Q = 0
  (gamma -> infinity), so the published curves end at r_Q = 0.8 mm (d = 94.2 mm); they are
  drawn from deep to superficial sources (this sets the phase of the SQUID dash pattern) and
  the OPM curve lies on top of the SQUID curve.
* Equal-SNR marker.  The published dotted line was measured at 27.53 mm, i.e. 0.14 mm
  shallower than the root of Eq. (3), 27.665 mm (the "28 mm" of the caption), and not at a
  point of the 0.8-mm grid (27.0 mm).  It is reproduced by the deepest depth with
  SNR_OPM > SNR_SQUID on d = linspace(h - b, h, 250) (27.530 mm); a few other grid sizes
  also land within the measurement precision, so the authors' exact rule is not
  identifiable.  Both values are reported; ``--deq-marker exact`` draws the root of Eq. (3).
* Artwork: curves 1.5 pt (blue solid / red dashed / black dotted), spines and ticks
  0.8 pt, tick labels DejaVu Sans 8 pt, all other text Myriad Pro 10 pt (Times New Roman
  italic for the variable d), subscripts at 58.3 % size / 33.3 % baseline shift, arrows
  with concave-sided heads, matplotlib's default legend box.  Positions of all elements
  were measured on the published raster (section 4).
* Raster registration (PNG only).  The figure was finished in Illustrator and rasterized
  there: axis-parallel strokes sit on the pixel grid, while the curves and the (unhinted,
  slightly heavier) tick labels keep sub-pixel positions.  The fitted offsets of those
  elements relative to the pixel-aligned frame (<= 1.4 px, 0.03 mm) are applied to the
  1200-dpi PNG only, so that it can be compared pixel by pixel with the published raster
  (``verify_figure3.py``); the PDF/SVG keep the undistorted geometry.
* Figure-specific layout (annotations, arrows, legend) is that of the published eta = 3
  figure; other ``--eta`` values only change the curves and the d_eq marker.

Usage
-----
    python replicate_figure3.py                 # writes figure3_output/ next to this file
    python replicate_figure3.py --outdir out --dtheta 0.05 --deq-marker exact
    python verify_figure3.py                    # quantitative comparison with the preprint

Outputs: Figure3_replicated.{png,pdf,svg}, figure3_data.csv, figure3_summary.json
Requires: mne, numpy, scipy, matplotlib >= 3.10, Pillow.  Text is set in open fonts shipped
with matplotlib (DejaVu Sans; STIX italic for d) unless ``--artwork-fonts`` asks for the published
artwork's Myriad Pro / Times New Roman (found e.g. inside Adobe Acrobat): that rendering is the
one checked pixel by pixel (verify_figure3.py), but its outputs embed licensed fonts and are not
for redistribution.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import matplotlib

if tuple(int(v) for v in matplotlib.__version__.split(".")[:2]) < (3, 10):
    raise ImportError("matplotlib >= 3.10 is required (FT2Font LoadFlags, text layout)")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager, patheffects  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.ft2font import LoadFlags  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.transforms import Affine2D, ScaledTranslation  # noqa: E402
from PIL import Image  # noqa: E402

import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402

# simulate_evoked warns that the (positive) source amplitudes are not signed currents of an
# inverse solution -- irrelevant for forward simulation
warnings.filterwarnings("ignore", message="Source estimate only contains currents with positive values")

# =========================================================================================
# 1. Model parameters (paper Sec. 2.1, Table 1 "Adult", Fig. 3 caption)
# =========================================================================================
MU0_OVER_4PI = 1e-7  # T m / A
Q_DIPOLE = 30e-9  # dipole moment [A m] (30 nAm), tangential (+y)
HEAD_RADIUS = 0.095  # h [m]
BRAIN_RADIUS = 0.080  # b [m]
XI_OPM = 0.0  # scalp-to-sensor distance, on-scalp OPM [m]
XI_SQUID = 0.018  # scalp-to-sensor distance, off-scalp SQUID (Triux neo) [m]
ETA = 3.0  # relative noise level sigma_OPM / sigma_SQUID
NOISE_REF_RADIUS = 0.8 * BRAIN_RADIUS  # see module docstring ("Absolute noise level")
N_CURVE = 101  # curves: r_Q = linspace(0, b, N_CURVE) without the undefined r_Q = 0 sample
N_DEQ_GRID = 250  # equal-SNR marker grid: d = linspace(h - b, h, N_DEQ_GRID)
DTHETA_DEG = 0.05  # angular spacing of the simulated sensor arcs [deg]


# =========================================================================================
# 2. Forward simulation with MNE-Python
# =========================================================================================
def make_arc_info(radius: float, dtheta_deg: float = DTHETA_DEG) -> tuple[mne.Info, np.ndarray]:
    """Radial point magnetometers, evenly spaced on the great circle y = 0 of a sphere."""
    alpha = np.linspace(0.0, 2.0 * np.pi, int(round(360.0 / dtheta_deg)), endpoint=False)
    names = [f"MAG{i:05d}" for i in range(alpha.size)]
    info = mne.create_info(names, sfreq=1000.0, ch_types="mag")
    e_y = np.array([0.0, 1.0, 0.0])
    with info._unlock():  # MNE keeps Info locked; the sensor geometry is set explicitly here
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")  # device == head
        for ch, a in zip(info["chs"], alpha):
            e_r = np.array([np.sin(a), 0.0, np.cos(a)])  # radial unit vector = coil normal
            e_x = np.cross(e_y, e_r)  # (e_x, e_y, e_r) right-handed coil frame
            ch["loc"][:3] = radius * e_r
            ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = e_x, e_y, e_r
            ch["coil_type"] = FIFF.FIFFV_COIL_POINT_MAGNETOMETER
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
    return info, alpha


def _parabolic_peak(values: np.ndarray) -> np.ndarray:
    """Peak of each column, refined by a parabola through the maximum sample and its two
    neighbours (evenly spaced samples on a closed circle, so neighbours wrap around)."""
    n = values.shape[0]
    i = values.argmax(axis=0)
    cols = np.arange(values.shape[1])
    y0, y1, y2 = values[(i - 1) % n, cols], values[i, cols], values[(i + 1) % n, cols]
    denom = y0 - 2.0 * y1 + y2
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(denom < 0, 0.5 * (y0 - y2) / denom, 0.0)
    return y1 - 0.25 * (y0 - y2) * delta


def radial_field_biot_savart(r_q: np.ndarray, sensor_radius: float, alpha: np.ndarray,
                             q: float = Q_DIPOLE) -> np.ndarray:
    """Radial field [T] of y-oriented dipoles at (0, 0, r_q) at the arc sensors (n_sensors,
    n_dipoles): B_r = mu0/4pi (Q x (r - r_Q)) . e_r / |r - r_Q|^3 = -mu0/4pi Q r_Q sin(alpha)
    / |r - r_Q|^3.  In a spherical conductor the volume currents add no radial field."""
    s, z = sensor_radius, np.asarray(r_q, float)[None, :]
    dist3 = (s**2 + z**2 - 2.0 * s * z * np.cos(alpha)[:, None]) ** 1.5
    return -MU0_OVER_4PI * q * z * np.sin(alpha)[:, None] / dist3


def make_head_sphere() -> mne.bem.ConductorModel:
    """Spherical head (origin at the sphere centre).  Only the origin matters for MEG
    (Sarvas); the two shells (brain b, scalp h) just satisfy MNE's geometry checks -- the
    inner radius must contain sources at r_Q = b exactly."""
    return mne.make_sphere_model(
        r0=(0.0, 0.0, 0.0),
        head_radius=HEAD_RADIUS,
        relative_radii=(np.nextafter(BRAIN_RADIUS / HEAD_RADIUS, 1.0), 1.0),
        sigmas=(0.33, 0.33),
        verbose=False,
    )


def simulate_peak_field_mne(r_q: np.ndarray, sensor_radius: float,
                            dtheta_deg: float = DTHETA_DEG) -> dict:
    """Simulate with MNE-Python the radial field of y-oriented Q_DIPOLE dipoles at
    (0, 0, r_q) on an arc of radial point magnetometers of the given radius.

    Returns a dict with ``peak`` (interpolated peak of |B_r| per dipole, i.e. paper Eq. 1),
    ``peak_pos`` / ``peak_neg`` (peaks of the two lobes) and ``forward_rel_err`` (largest
    sensor-wise deviation from the Biot-Savart field, relative to each dipole's peak)."""
    info, alpha = make_arc_info(sensor_radius, dtheta_deg)
    sphere = make_head_sphere()
    n = r_q.size
    dip = mne.Dipole(
        times=np.arange(n) * 1e-3,  # one dipole (depth) per "time sample"
        pos=np.column_stack([np.zeros(n), np.zeros(n), r_q]),
        amplitude=np.full(n, Q_DIPOLE),
        ori=np.tile([0.0, 1.0, 0.0], (n, 1)),
        gof=np.full(n, 100.0),
        verbose=False,
    )
    # on-scalp sensors lie on the scalp sphere itself (s = h)
    on_inside = "ignore" if sensor_radius <= HEAD_RADIUS else "raise"
    fwd, stc = mne.make_forward_dipole(dip, sphere, info, on_inside=on_inside, verbose=False)
    evoked = mne.simulation.simulate_evoked(fwd, stc, info, cov=None, nave=np.inf, verbose=False)
    # (n_sensors, n_dipoles) [T].  MNE follows Biot-Savart, B = mu0/4pi Q x (r - r_Q)/|r - r_Q|^3,
    # so the lobe at x > 0 is negative; the printed form of the paper has the opposite sign.
    # Only the peak magnitude enters Eq. (1), hence |B_r| below.
    b_r = evoked.data
    reference = radial_field_biot_savart(r_q, sensor_radius, alpha)
    scale = np.abs(reference).max(axis=0)
    ok = scale > 0
    forward_err = float(np.max(np.abs(b_r - reference)[:, ok] / scale[ok])) if ok.any() else 0.0
    return dict(
        peak=_parabolic_peak(np.abs(b_r)),
        peak_pos=_parabolic_peak(np.clip(b_r, 0.0, None)),
        peak_neg=_parabolic_peak(np.clip(-b_r, 0.0, None)),
        forward_rel_err=forward_err,
    )


def bmax_analytic(r_q, s, q: float = Q_DIPOLE) -> np.ndarray:
    """Closed-form peak radial field (paper Eq. 1; Boto et al., 2017); 0 at r_Q = 0."""
    r_q = np.atleast_1d(np.asarray(r_q, dtype=float))
    out = np.zeros_like(r_q)
    m = r_q > 0
    rq = r_q[m]
    gamma = (s**2 + rq**2) / (2.0 * s * rq)
    cos0 = 3.0 / (gamma + np.sqrt(gamma**2 + 3.0))  # = -gamma + sqrt(gamma^2 + 3), stable
    out[m] = (
        MU0_OVER_4PI * q * rq * np.sqrt(1.0 - cos0**2)
        / (s**2 + rq**2 - 2.0 * s * rq * cos0) ** 1.5
    )
    return out


# =========================================================================================
# 3. SNR and equal-SNR depth
# =========================================================================================
@dataclass
class Results:
    depth_mm: np.ndarray  # curve samples, deep -> superficial (d = 94.2 ... 15 mm)
    signal_opm_pT: np.ndarray
    signal_squid_pT: np.ndarray
    snr_opm: np.ndarray
    snr_squid: np.ndarray
    sigma_opm_pT: float
    sigma_squid_pT: float
    d_eq_grid_mm: float
    d_eq_exact_mm: float
    dtheta_deg: float
    max_rel_err_forward: float  # MNE vs Biot-Savart, sensor by sensor
    max_rel_err_peak: float  # interpolated peak vs Eq. (1)
    max_lobe_asymmetry: float


def run_simulation(n_curve: int = N_CURVE, dtheta_deg: float = DTHETA_DEG, eta: float = ETA,
                   n_deq_grid: int = N_DEQ_GRID) -> Results:
    h, b = HEAD_RADIUS, BRAIN_RADIUS
    # curve samples, deep -> superficial (r_Q = 0 is excluded: Eq. (1) is undefined there)
    r_curve = np.linspace(0.0, b, n_curve)[1:]
    # depth grid of the equal-SNR marker (brain surface -> sphere centre, r_Q = 0 included)
    depth_grid = np.linspace(h - b, h, n_deq_grid)
    r_grid = h - depth_grid
    # all dipoles (curves, marker grid, noise reference) are simulated in one run per array
    r_all = np.concatenate([r_curve, r_grid, [NOISE_REF_RADIUS]])
    sims = {
        "OPM": simulate_peak_field_mne(r_all, h + XI_OPM, dtheta_deg),
        "SQUID": simulate_peak_field_mne(r_all, h + XI_SQUID, dtheta_deg),
    }
    fwd_err = max(sim["forward_rel_err"] for sim in sims.values())
    peak_err, asym = 0.0, 0.0
    for key, xi in (("OPM", XI_OPM), ("SQUID", XI_SQUID)):
        ana = bmax_analytic(r_all, h + xi)
        pk = sims[key]["peak"]
        ok = ana > 0
        peak_err = max(peak_err, float(np.max(np.abs(pk[ok] - ana[ok]) / ana[ok])))
        pos, neg = sims[key]["peak_pos"][ok], sims[key]["peak_neg"][ok]
        asym = max(asym, float(np.max(np.abs(pos - neg) / np.maximum(pos, neg))))
        if not np.allclose(pk[~ok], 0.0, atol=1e-20):
            raise RuntimeError("non-zero field for a dipole at the sphere centre")
    if fwd_err > 1e-5:  # MNE stores sensor positions in float32 (~1e-7)
        raise RuntimeError(f"MNE forward deviates from Biot-Savart: rel. err {fwd_err:.2e}")
    tol = min(1e-6 + 1e-3 * dtheta_deg**3, 1e-4)  # interpolation error ~ dtheta^3, capped
    if peak_err > tol:
        raise RuntimeError(f"MNE peak deviates from Eq. (1): rel. err {peak_err:.2e} > {tol:.1e}")

    n1, n2 = r_curve.size, r_grid.size
    sigma_squid = sims["SQUID"]["peak"][-1]  # B_r^max(0.8 b, h + xi)
    sigma_opm = eta * sigma_squid
    b_opm, b_squid = sims["OPM"]["peak"][:n1], sims["SQUID"]["peak"][:n1]
    g_opm = sims["OPM"]["peak"][n1:n1 + n2] / sigma_opm
    g_squid = sims["SQUID"]["peak"][n1:n1 + n2] / sigma_squid

    # equal-SNR marker: deepest grid depth (from the brain surface down) with
    # SNR_OPM > SNR_SQUID.  The last grid point (r_Q = 0) has zero field for both sensors, so
    # "not better" always occurs there; a crossing inside the brain requires an earlier one.
    first_not_better = np.flatnonzero(~(g_opm > g_squid))
    if first_not_better.size and 0 < first_not_better[0] < n2 - 1:
        d_eq_grid = depth_grid[first_not_better[0] - 1]
    else:
        d_eq_grid = np.nan  # no crossing inside the brain (eta <= 1.6829 or eta >= 5.3086)

    def f(d):
        return bmax_analytic(h - d, h + XI_OPM)[0] / bmax_analytic(h - d, h + XI_SQUID)[0] - eta

    lo, hi = h - b, h - 1e-9
    d_eq_exact = brentq(f, lo, hi, xtol=1e-12) if f(lo) > 0 > f(hi) else np.nan

    return Results(
        depth_mm=(h - r_curve) * 1e3,
        signal_opm_pT=b_opm * 1e12,
        signal_squid_pT=b_squid * 1e12,
        snr_opm=b_opm / sigma_opm,
        snr_squid=b_squid / sigma_squid,
        sigma_opm_pT=float(sigma_opm * 1e12),
        sigma_squid_pT=float(sigma_squid * 1e12),
        d_eq_grid_mm=float(d_eq_grid * 1e3),
        d_eq_exact_mm=float(d_eq_exact * 1e3),
        dtheta_deg=float(dtheta_deg),
        max_rel_err_forward=fwd_err,
        max_rel_err_peak=peak_err,
        max_lobe_asymmetry=asym,
    )


# =========================================================================================
# 4. Figure 3 layout, measured on the published 1200-ppi raster (see verify_figure3.py)
#    Coordinates are pixels of that raster, origin top-left; 1 pt = 1200/72 px.
# =========================================================================================
PUB_PPI = 1200.0
PX_PER_PT = PUB_PPI / 72.0
CANVAS_PX = (4693.0, 3207.0)  # width, height of the published figure
AX_X0_PX = {"A": 194.80, "B": 2633.00}  # x pixel of d = 0 mm
AX_PX_PER_MM = 21.615  # 20-mm tick spacing = 432.3 px
AX_Y_BOTTOM_PX = 2681.0  # the x-axis (value 0)
AX_Y_TOP_PX = 233.2  # top of the axes box
XLIM_MM = (0.0, 95.0)

FONT_SIZE = 10.0  # pt, all Myriad Pro text
TICK_LABEL_SIZE = 8.0  # pt, DejaVu Sans (matplotlib default font)
LINE_WIDTH = 1.5  # pt, curves and d_eq line
SPINE_WIDTH = 0.8  # pt
SUB_SCALE, SUB_SHIFT = 0.583, 0.333  # Illustrator subscript: size and baseline shift

# Raster registration.  The published raster was produced from vector art (matplotlib
# output finished in Illustrator): axis-parallel strokes (spines, ticks, d_eq line, legend)
# sit on the pixel grid, the curves keep their sub-pixel positions, and the DejaVu tick
# labels are rendered unhinted with slightly heavier outlines.  Fitted on that raster; the
# offsets are applied to the PNG only (see save_outputs):
CURVE_OFFSET_PX = {"A": (0.31, 0.89), "B": (0.63, 1.42)}  # curves & tick labels vs. frame
TICK_LABEL_STROKE_PT = 0.07  # outline weight of the tick-label glyphs
TEXT_HINTING = "no_hinting"  # outline (unhinted) glyph rendering
FIGURE_RC = {"text.hinting": TEXT_HINTING}

# pen origins (x, baseline) of the text elements [px]
PANEL_LETTERS = {"A": (8.2, 139.0), "B": (2454.0, 139.0)}
PANEL_TITLES = {"A": ((298.75, 139.0), "Signal [pT]"), "B": ((2769.75, 139.0), "Signal-to-Noise Ratio")}
XLABEL_ORIGINS = {"A": (579.8, 3155.0), "B": (3003.55, 3155.0)}
XLABEL_D_NUDGE_PX = 3.0  # the italic d sits 3 px right of its pen position (DTP nudge)

# annotations: list of (origin, runs); run = (text, subscript?).  The font's own kerning
# (e.g. the Myriad Pro "PM" pair inside the OPM subscript) comes from matplotlib's layout.
ANNOTATIONS = {
    "A": [((1044.25, 1151.5), [("S", False), ("OPM", True), (" > ", False), ("S", False), ("SQUID", True)])],
    "B": [
        ((3351.5, 1081.4), [("SNR", False), ("OPM", True), (" > ", False), ("SNR", False), ("SQUID", True)]),
        ((3381.5, 1691.0), [("SNR", False), ("SQUID", True), (" > ", False), ("SNR", False), ("OPM", True)]),
    ],
}

# arrows: polyline vertices [px] (arrowheads at the flagged ends)
ARROW_WIDTH_PT = 0.785
ARROW_HEAD_LENGTH_PX = 106.0
ARROW_HEAD_WIDTH_PX = 64.0
ARROW_HEAD_CONTROL = (0.5, 0.30)  # Bezier control point: fraction of length, of half-width
ARROWS = [
    dict(xy=[(619.0, 1261.4), (981.3, 1195.8), (1078.9, 2102.9)], start=True, end=True),
    dict(xy=[(3335.1, 1161.6), (3070.2, 1280.1)], start=False, end=True),
    dict(xy=[(3973.2, 1877.0), (3823.1, 2315.7)], start=False, end=True),
]

# legend (matplotlib default legend box of the original, labels re-set in Myriad Pro)
LEGEND = dict(
    box=(1028.0, 316.0, 2165.0, 855.0),  # stroke centre-line: x0, y0 (top), x1, y1 [px]
    handles=[(1095.5, 1428.3, 451.0), (1095.5, 1428.3, 696.0)],  # x0, x1, y [px]
    labels=[(1538.5, 496.25), (1538.0, 739.25)],  # pen origins [px]
)

MYRIAD_FILES = ("MyriadPro-Regular.otf", "MyriadPro-Bold.otf")
MYRIAD_SEARCH = [
    "/Applications/Adobe Acrobat DC/Adobe Acrobat.app/Contents/Resources/Resource/Font",
    "/Applications/Adobe Acrobat Reader DC.app/Contents/Resources/Resource/Font",
    "/Applications/Adobe Acrobat Reader.app/Contents/Resources/Resource/Font",
    "/Library/Application Support/Adobe/Acrobat/DC/WebResources/Resource0/app1/fonts",
    "~/Library/Fonts",
    "/Library/Fonts",
    "C:/Program Files/Adobe/Acrobat DC/Resource/Font",
    "C:/Program Files (x86)/Adobe/Acrobat Reader DC/Reader/../Resource/Font",
    "C:/Windows/Fonts",
    "/usr/share/fonts",
    "~/.local/share/fonts",
    "~/.fonts",
]
FALLBACK_FAMILIES = ("Source Sans 3", "Source Sans Pro", "Segoe UI", "Helvetica", "Arial")


@dataclass
class Fonts:
    regular: FontProperties
    bold: FontProperties
    italic_serif: FontProperties
    tick: FontProperties
    note: str
    artwork: bool = False


OPEN_SANS, OPEN_SERIF = "DejaVu Sans", "STIXGeneral"  # Bitstream Vera licence / SIL OFL, shipped with matplotlib


def setup_fonts(font_dir: str | None = None, artwork: bool = False) -> Fonts:
    """Open fonts by default (DejaVu Sans; STIX italic serif). With ``artwork`` (or a font folder,
    or FIG3_FONT_DIR): Myriad Pro, the published figure font, from an installed family, a font
    folder or Adobe Acrobat (otherwise a similar sans-serif with a warning), and Times New Roman
    italic. Those fonts are licensed: outputs that embed them are not for redistribution."""
    tick = FontProperties(family="DejaVu Sans", size=TICK_LABEL_SIZE)
    if not (artwork or font_dir or os.environ.get("FIG3_FONT_DIR")):
        return Fonts(FontProperties(family=OPEN_SANS, size=FONT_SIZE), FontProperties(family=OPEN_SANS, weight="bold", size=FONT_SIZE),
                     FontProperties(family=OPEN_SERIF, style="italic", size=FONT_SIZE), tick,
                     f"open fonts ({OPEN_SANS}; {OPEN_SERIF} italic)", False)
    installed = {f.name for f in font_manager.fontManager.ttflist}

    def search(dirs):
        found = {}
        for d in dirs:
            for name in MYRIAD_FILES:
                if name not in found:
                    hits = sorted(glob.glob(os.path.join(os.path.expanduser(d), "**", name), recursive=True))
                    if hits:
                        found[name] = hits[0]
        return found if len(found) == len(MYRIAD_FILES) else None

    found = search([font_dir]) if font_dir else None
    if found is None and "Myriad Pro" in installed:
        regular = FontProperties(family="Myriad Pro", size=FONT_SIZE)
        bold = FontProperties(family="Myriad Pro", weight="bold", size=FONT_SIZE)
        note = "Myriad Pro (installed)"
    else:
        found = found or search([d for d in [os.environ.get("FIG3_FONT_DIR")] + MYRIAD_SEARCH if d])
        if found:
            for path in found.values():
                font_manager.fontManager.addfont(path)
            regular = FontProperties(fname=found["MyriadPro-Regular.otf"], size=FONT_SIZE)
            bold = FontProperties(fname=found["MyriadPro-Bold.otf"], size=FONT_SIZE)
            note = f"Myriad Pro ({os.path.dirname(found['MyriadPro-Regular.otf'])})"
        else:
            family = next((f for f in FALLBACK_FAMILIES if f in installed), "DejaVu Sans")
            regular = FontProperties(family=family, size=FONT_SIZE)
            bold = FontProperties(family=family, weight="bold", size=FONT_SIZE)
            note = f"Myriad Pro not found; using {family}"
            warnings.warn(note + " (pass --font-dir to use Myriad Pro)", stacklevel=2)
    serif = "Times New Roman" if "Times New Roman" in installed else "STIXGeneral"
    italic = FontProperties(family=serif, style="italic", size=FONT_SIZE)
    return Fonts(regular, bold, italic, tick, note, True)


def _advance_px(text: str, prop: FontProperties) -> float:
    """Pen advance of a string in published pixels, from matplotlib's own text layout
    (kerning included): the position of a space appended to the string."""
    font = font_manager.get_font(font_manager.findfont(prop))  # honours .ttc face indices
    font.set_size(prop.get_size_in_points(), PUB_PPI)
    return float(font.set_text(text + " ", 0.0, flags=LoadFlags.NO_HINTING)[-1, 0]) / 64.0


class Canvas:
    """Draw on the figure in published-pixel coordinates (origin top-left)."""

    def __init__(self, fig):
        self.fig = fig
        w, h = CANVAS_PX
        self.transform = Affine2D().scale(1.0 / w, -1.0 / h).translate(0.0, 1.0) + fig.transFigure

    def text(self, x, y, s, prop, **kw):
        return self.fig.text(x, y, s, fontproperties=prop, ha="left", va="baseline",
                             transform=self.transform, **kw)

    def typeset(self, x, y, runs):
        """Set consecutive runs like a DTP application.  A run is
        (text, FontProperties, baseline shift px[, nudge px[, kerning px]]): the kerning
        moves the pen before the run (and everything after it), the nudge moves only the
        run's own glyphs."""
        for run in runs:
            s, prop, shift = run[:3]
            nudge = run[3] if len(run) > 3 else 0.0
            x += run[4] if len(run) > 4 else 0.0
            self.text(x + nudge, y + shift, s, prop)
            x += _advance_px(s, prop)


class RasterRegistration:
    """Sub-pixel offsets of artists relative to the pixel-aligned axes frame of the published
    raster.  The offset is composed with each artist's transform when switched on, so it
    stays valid if the artists' own transforms change beforehand (e.g. tick_params)."""

    def __init__(self):
        self._shifts, self._plain = [], {}

    def add(self, artist, shift):
        self._shifts.append((artist, shift))

    def __call__(self, enabled: bool):
        for artist, shift in self._shifts:
            if enabled and artist not in self._plain:
                self._plain[artist] = artist.get_transform()
                artist.set_transform(self._plain[artist] + shift)
            elif not enabled and artist in self._plain:
                artist.set_transform(self._plain.pop(artist))


def _arrow_head_path(tip, direction, length, width, control):
    """Arrowhead with a straight base and concave (quadratic Bezier) sides."""
    u = np.asarray(direction, float) / np.linalg.norm(direction)
    n = np.array([-u[1], u[0]])
    tip = np.asarray(tip, float)
    base = tip - length * u
    cu, cw = control
    pts = [
        tip,
        tip - cu * length * u + cw * 0.5 * width * n,
        base + 0.5 * width * n,
        base - 0.5 * width * n,
        tip - cu * length * u - cw * 0.5 * width * n,
        tip,
        tip,
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3, MplPath.LINETO,
             MplPath.CURVE3, MplPath.CURVE3, MplPath.CLOSEPOLY]
    return MplPath(pts, codes)


def draw_arrow(canvas: Canvas, xy, start: bool, end: bool):
    """Polyline (butt caps, mitre joins) with concave-sided arrowheads at the flagged ends."""
    xy = np.asarray(xy, float)
    shaft = xy.copy()
    inset = 0.6 * ARROW_HEAD_LENGTH_PX  # the shaft ends inside the head
    heads = []
    if start:
        u = xy[0] - xy[1]
        heads.append((xy[0], u))
        shaft[0] = xy[0] - inset * u / np.linalg.norm(u)
    if end:
        u = xy[-1] - xy[-2]
        heads.append((xy[-1], u))
        shaft[-1] = xy[-1] - inset * u / np.linalg.norm(u)
    canvas.fig.add_artist(Line2D(shaft[:, 0], shaft[:, 1], transform=canvas.transform, color="k",
                                 lw=ARROW_WIDTH_PT, solid_capstyle="butt", solid_joinstyle="miter"))
    for tip, u in heads:
        path = _arrow_head_path(tip, u, ARROW_HEAD_LENGTH_PX, ARROW_HEAD_WIDTH_PX, ARROW_HEAD_CONTROL)
        canvas.fig.add_artist(PathPatch(path, transform=canvas.transform, facecolor="k",
                                        edgecolor="none", lw=0))


def draw_legend(canvas: Canvas, fonts: Fonts):
    """matplotlib's default legend box (fontsize 10: rounded frame, 0.8 grey edge, framealpha
    0.8, handlelength 2) with the entry labels set in Myriad Pro as in the published figure."""
    x0, y0, x1, y1 = LEGEND["box"]
    frame = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                           boxstyle=f"round,pad=0,rounding_size={0.2 * FONT_SIZE * PX_PER_PT}",
                           transform=canvas.transform, facecolor="white", edgecolor="0.8",
                           lw=1.0, alpha=0.8, zorder=5)
    canvas.fig.add_artist(frame)
    styles = [dict(color="b", ls="-"), dict(color="r", ls="--")]
    for (hx0, hx1, hy), st in zip(LEGEND["handles"], styles):
        canvas.fig.add_artist(Line2D([hx0, hx1], [hy, hy], transform=canvas.transform,
                                     lw=LINE_WIDTH, zorder=6, **st))
    for (tx, ty), label in zip(LEGEND["labels"], ("OPM", "SQUID")):
        canvas.text(tx, ty, label, fonts.regular, zorder=6)


def make_figure(res: Results, fonts: Fonts, deq_marker: str = "published", annotations: bool = True,
                source_word: str = "Dipole"):
    """Build Figure 3.  Returns (fig, axes, registration); calling ``registration(True)``
    applies the sub-pixel raster registration of the curves and tick labels.
    ``annotations=False`` omits the arrows and "S_OPM > S_SQUID"-type labels, which belong to
    the published eta = 3 dipole curves (used for the extended-source variants)."""
    w_px, h_px = CANVAS_PX
    fig = plt.figure(figsize=(w_px / PUB_PPI, h_px / PUB_PPI), dpi=PUB_PPI)
    canvas = Canvas(fig)
    registration = RasterRegistration()
    width_px = (XLIM_MM[1] - XLIM_MM[0]) * AX_PX_PER_MM
    series = {
        "A": (res.signal_opm_pT, res.signal_squid_pT),
        "B": (res.snr_opm, res.snr_squid),
    }
    axes = {}
    for key, (y_opm, y_squid) in series.items():
        rect = [AX_X0_PX[key] / w_px, 1.0 - AX_Y_BOTTOM_PX / h_px,
                width_px / w_px, (AX_Y_BOTTOM_PX - AX_Y_TOP_PX) / h_px]
        ax = fig.add_axes(rect)
        # drawing order of the published figure: SQUID, OPM on top, then the d_eq marker
        curves = [ax.plot(res.depth_mm, y_squid, "r--", lw=LINE_WIDTH, label="SQUID")[0],
                  ax.plot(res.depth_mm, y_opm, "b-", lw=LINE_WIDTH, label="OPM")[0]]
        d_eq = res.d_eq_grid_mm if deq_marker == "published" else res.d_eq_exact_mm
        if np.isfinite(d_eq):
            ax.axvline(d_eq, color="k", ls=":", lw=LINE_WIDTH)
        ax.set_xlim(*XLIM_MM)
        ax.set_ylim(bottom=0.0)  # top: matplotlib auto-scaling (5 % margin)
        ax.set_xticks([0, 20, 40, 60, 80])
        ax.set_yticks(np.arange(0, int(ax.get_ylim()[1]) + 1))  # fixed ticks: fixed labels
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_linewidth(SPINE_WIDTH)
        ax.tick_params(width=SPINE_WIDTH, length=3.5, pad=3.5)
        stroke = [patheffects.withStroke(linewidth=TICK_LABEL_STROKE_PT, foreground="k")]
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        for lab in labels:
            lab.set_fontproperties(fonts.tick)
            lab.set_path_effects(stroke)
        dx, dy = CURVE_OFFSET_PX[key]
        shift = ScaledTranslation(dx / PUB_PPI, -dy / PUB_PPI, fig.dpi_scale_trans)
        for artist in curves + labels:
            registration.add(artist, shift)
        axes[key] = ax

    sub_prop = fonts.regular.copy()
    sub_prop.set_size(FONT_SIZE * SUB_SCALE)
    sub_shift = SUB_SHIFT * FONT_SIZE * PX_PER_PT
    for key in ("A", "B"):
        canvas.text(*PANEL_LETTERS[key], key, fonts.bold)
        (tx, ty), title = PANEL_TITLES[key]
        canvas.text(tx, ty, title, fonts.regular)
        x, y = XLABEL_ORIGINS[key]
        canvas.typeset(x, y, [(f"{source_word} depth (", fonts.regular, 0.0),
                              ("d", fonts.italic_serif, 0.0, XLABEL_D_NUDGE_PX),
                              (") [mm]", fonts.regular, 0.0)])
        for (ox, oy), runs in (ANNOTATIONS[key] if annotations else []):
            canvas.typeset(ox, oy, [(s, sub_prop, sub_shift) if sub else (s, fonts.regular, 0.0)
                                    for s, sub in runs])
    draw_legend(canvas, fonts)
    for arrow in (ARROWS if annotations else []):
        draw_arrow(canvas, **arrow)
    return fig, axes, registration


# =========================================================================================
# 5. Driver
# =========================================================================================
def _json_safe(value):
    """NaN / inf -> None, recursively (strict JSON)."""
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_outputs(res: Results, fig, registration, outdir: Path, fonts: Fonts, eta: float,
                 deq_marker: str, raster_registration: bool = True) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext in ("png", "pdf", "svg"):
        # the sub-pixel registration only reproduces the published *raster*
        registration(raster_registration and ext == "png")
        p = outdir / f"Figure3_replicated.{ext}"
        fig.savefig(p, dpi=PUB_PPI, facecolor="white")
        paths[ext] = str(p)
    registration(False)
    with Image.open(paths["png"]) as im:  # Agg truncates the canvas size to whole pixels
        if im.size != tuple(int(v) for v in CANVAS_PX):
            raise RuntimeError(f"PNG size {im.size} differs from the published canvas {CANVAS_PX}")
    with open(outdir / "figure3_data.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["depth_mm", "signal_OPM_pT", "signal_SQUID_pT", "SNR_OPM", "SNR_SQUID"])
        order = np.argsort(res.depth_mm)
        cols = (res.depth_mm, res.signal_opm_pT, res.signal_squid_pT, res.snr_opm, res.snr_squid)
        for row in zip(*(c[order] for c in cols)):
            wr.writerow([f"{v:.10g}" for v in row])
    summary = {k: v for k, v in asdict(res).items() if not isinstance(v, np.ndarray)}
    surface = int(np.argmin(res.depth_mm))  # brain surface, d = h - b
    summary.update(
        eta=eta, head_radius_mm=HEAD_RADIUS * 1e3, brain_radius_mm=BRAIN_RADIUS * 1e3,
        xi_squid_mm=XI_SQUID * 1e3, dipole_nAm=round(Q_DIPOLE * 1e9, 9), n_curve_samples=int(res.depth_mm.size),
        noise_reference_radius_mm=NOISE_REF_RADIUS * 1e3, deq_marker=deq_marker,
        raster_registration_png=raster_registration,
        signal_at_brain_surface_pT={"OPM": float(res.signal_opm_pT[surface]),
                                    "SQUID": float(res.signal_squid_pT[surface])},
        snr_at_brain_surface={"OPM": float(res.snr_opm[surface]), "SQUID": float(res.snr_squid[surface])},
        font=fonts.note, artwork_fonts=fonts.artwork, mne_version=mne.__version__, matplotlib_version=matplotlib.__version__,
        outputs={k: Path(v).name for k, v in paths.items()},  # relative to the summary (no local absolute paths)
    )
    if eta == ETA:
        summary["d_eq_paper_caption_mm"] = 28.0
    with open(outdir / "figure3_summary.json", "w") as fh:
        json.dump(_json_safe(summary), fh, indent=2, allow_nan=False)
    return summary


def render_png(path, eta: float = ETA, deq_marker: str = "published", raster_registration: bool = True,
               dtheta_deg: float = DTHETA_DEG, font_dir: str | None = None, artwork_fonts: bool = False) -> Results:
    """Simulate and write only the 1200-dpi PNG (used by verify_figure3.py)."""
    res = run_simulation(dtheta_deg=dtheta_deg, eta=eta)
    fonts = setup_fonts(font_dir, artwork_fonts)
    with plt.rc_context(FIGURE_RC):
        fig, _, registration = make_figure(res, fonts, deq_marker)
        registration(raster_registration)
        fig.savefig(path, dpi=PUB_PPI, facecolor="white")
    plt.close(fig)
    return res


def _positive_float(text: str) -> float:
    value = float(text)
    if not (np.isfinite(value) and value > 0):
        raise argparse.ArgumentTypeError(f"must be a positive number, got {text}")
    return value


def _sensor_spacing(text: str) -> float:
    value = _positive_float(text)
    if value > 0.5:
        raise argparse.ArgumentTypeError(f"sensor spacing must be <= 0.5 deg, got {text}")
    return value


def main(argv=None) -> Results:
    ap = argparse.ArgumentParser(
        description="End-to-end replication of Figure 3 of Jas et al. (2026), bioRxiv "
                    "2026.08.17.744953: MNE-Python simulation of on-scalp (OPM) vs off-scalp "
                    "(SQUID) signal and SNR as a function of dipole depth.")
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent / "figure3_output"))
    ap.add_argument("--dtheta", type=_sensor_spacing, default=DTHETA_DEG,
                    help="sensor spacing on the arcs [deg], at most 0.5 (default: %(default)s)")
    ap.add_argument("--eta", type=_positive_float, default=ETA,
                    help="relative noise level (the annotations are those of the published eta = 3 figure)")
    ap.add_argument("--deq-marker", choices=("published", "exact"), default="published",
                    help="dotted line at the published position (250-point grid) or at the root of Eq. (3)")
    ap.add_argument("--no-raster-registration", action="store_true",
                    help="do not apply the sub-pixel raster registration to the PNG")
    ap.add_argument("--artwork-fonts", action="store_true",
                    help="set text in the published artwork's licensed fonts (Myriad Pro, Times New Roman), as checked by "
                         "verify_figure3.py; not for redistribution (default: open fonts)")
    ap.add_argument("--font-dir", default=None, help="folder containing MyriadPro-*.otf (implies --artwork-fonts)")
    args = ap.parse_args(argv)

    mne.set_log_level("WARNING")
    res = run_simulation(dtheta_deg=args.dtheta, eta=args.eta)
    fonts = setup_fonts(args.font_dir, args.artwork_fonts)
    with plt.rc_context(FIGURE_RC):
        fig, _, registration = make_figure(res, fonts, args.deq_marker)
        summary = save_outputs(res, fig, registration, Path(args.outdir), fonts, args.eta,
                               args.deq_marker, raster_registration=not args.no_raster_registration)
    plt.close(fig)

    surface = int(np.argmin(res.depth_mm))
    print(f"MNE-Python {mne.__version__}: forward vs. Biot-Savart max rel. error {res.max_rel_err_forward:.1e}; "
          f"peak vs. Eq. (1) {res.max_rel_err_peak:.1e}; lobe asymmetry {res.max_lobe_asymmetry:.1e}")
    print(f"Signal at brain surface (d = {res.depth_mm[surface]:.0f} mm): OPM {res.signal_opm_pT[surface]:.3f} pT, "
          f"SQUID {res.signal_squid_pT[surface]:.3f} pT")
    print(f"Noise: sigma_SQUID = {res.sigma_squid_pT:.4f} pT, sigma_OPM = {res.sigma_opm_pT:.4f} pT "
          f"(eta = {args.eta:g})")
    marker = "published grid rule" if args.deq_marker == "published" else "root of Eq. (3)"
    print(f"Equal-SNR depth: grid {res.d_eq_grid_mm:.2f} mm, exact {res.d_eq_exact_mm:.2f} mm "
          f"(dotted line: {marker})")
    print("Wrote:", ", ".join(str(Path(args.outdir) / v) for v in summary["outputs"].values()))
    return res


if __name__ == "__main__":
    main()
