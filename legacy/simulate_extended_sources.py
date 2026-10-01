#!/usr/bin/env python3
"""
Figure 3 with extended sources: cortical patches of different sizes instead of a current
dipole (extension of Jas et al., 2026, bioRxiv 2026.08.17.744953).

Everything except the source is taken from ``replicate_figure3.py``: spherical head
(h = 95 mm, b = 80 mm), radial point magnetometers on-scalp (OPM, s = h) and off-scalp (SQUID,
s = h + 18 mm), noise sigma_SQUID = 0.3546 pT and sigma_OPM = eta * sigma_SQUID with eta = 3
(the noise does not depend on the source of interest), SNR = peak |B_r| / sigma, and the
Figure 3 layout.

Source model
------------
A patch of cortex at depth d is the spherical cap of geodesic radius rho on the sphere
r = r_Q = h - d (the whole patch lies at depth d), centred on the z-axis, carrying a uniform
+y current density (tangential at the patch centre).  By default the total moment equals the
paper's dipole, Q = 30 nAm, so rho -> 0 recovers Figure 3 and only the spatial extent changes;
with ``--moment-density`` the patch instead has a fixed current density (e.g. 0.5 nAm/mm^2) and
its moment grows with its area.  (In a spherical conductor radial currents produce no external
field, so the uniform +y current is equivalent to its tangential projection, whose
MEG-effective moment is Q (1 - <sin^2 psi>/2): 98 % of Q for rho = 20 mm near the brain
surface, 2/3 of Q for a hemispherical cap.)  A cap is limited to a hemisphere
(rho <= pi/2 * r_Q): large patches do not exist for the deepest sources.  The equal-SNR depth
only depends on the OPM/SQUID signal ratio and is therefore independent of the moment
normalization.

Simulation (MNE-Python)
-----------------------
Each patch is discretised with 8 (Gauss-Legendre in cos psi) x 32 (uniform in phi) = 256
dipoles (peak-field error < 2e-10 for rho <= 20 mm and < 1e-6 at 40 mm; sensor-wise error
relative to the peak ~1e-8 and ~1e-4; it grows to ~1e-3 at 100 mm, and the script warns above
40 mm).  All element dipoles of a patch share one "time sample" of an ``mne.Dipole``;
``mne.make_forward_dipole`` + ``mne.simulation.simulate_evoked`` then sum their fields at the
sensor arcs.  By symmetry the patch topography peaks in the xz-plane (confirmed on a 2-D sensor
grid during validation), so the peak over the great circle is the patch signal.  The run stops
if the MNE fields deviate from the Biot-Savart superposition of the elements (sensor by sensor)
by more than 1e-5, or the interpolated peaks from refined Biot-Savart peaks by more than the
interpolation-tolerance rule of ``replicate_figure3.run_simulation`` (1.7e-5 at the 0.25-degree
sensor spacing).

Outputs (``extended_sources_output/``)
--------------------------------------
* ``Figure3_patch_rho_XXmm.{png,pdf}``: the Figure 3 layout for each patch radius (the dotted
  line is the equal-SNR depth of that patch; rho = 0 is the dipole).
* ``Figure_extended_sources.{png,pdf,svg}``: (A) signal and (B) SNR for all radii,
  (C) SNR_OPM / SNR_SQUID vs depth, (D) equal-SNR depth vs eta (cf. Fig. 4E); circles mark
  the equal-SNR depths at the chosen eta; legend areas are cap areas at d = 15 mm.
* ``extended_sources_data.csv`` and ``extended_sources_summary.json``.

Usage
-----
    python simulate_extended_sources.py
    python simulate_extended_sources.py --radii 0 5 10 15 20 25 --eta 3
    python simulate_extended_sources.py --moment-density 0.5      # nAm/mm^2
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

import replicate_figure3 as R  # noqa: E402  (sets the matplotlib backend)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402
import mne  # noqa: E402

PATCH_RADII_MM = (0.0, 5.0, 10.0, 15.0, 20.0)
N_PSI, N_PHI = 8, 32  # quadrature nodes per patch
DTHETA_DEG = 0.25  # sensor spacing on the arcs [deg]
MAX_DIPOLES_PER_RUN = 6000  # MNE forward batch size
ETA_GRID = np.round(np.linspace(1.0, 6.0, 201), 6)
FINE_DEPTH_STEP_MM = 0.25  # grid for locating equal-SNR crossings before refinement
H, B, XI = R.HEAD_RADIUS, R.BRAIN_RADIUS, R.XI_SQUID
SENSOR_RADII = {"OPM": H + R.XI_OPM, "SQUID": H + XI}


# =========================================================================================
# Patch geometry
# =========================================================================================
def patch_defined(r_q: float, radius: float) -> bool:
    """A cap of geodesic radius `radius` fits on the sphere r_Q as at most a hemisphere."""
    return radius == 0.0 or (r_q > 0.0 and radius <= 0.5 * np.pi * r_q)


def patch_area(r_q: float, radius: float) -> float:
    """Area [m^2] of the spherical cap."""
    return 0.0 if radius == 0.0 else 2.0 * np.pi * r_q**2 * (1.0 - np.cos(radius / r_q))


def cap_quadrature(r_q: float, radius: float, n_psi: int = N_PSI, n_phi: int = N_PHI):
    """Nodes (n, 3) and area weights (n,) summing to 1 for the cap of geodesic radius
    `radius` on the sphere of radius r_q centred on +z: Gauss-Legendre in u = cos(psi)
    (the area element is r^2 du dphi) times the midpoint rule in phi."""
    if radius == 0.0:
        return np.array([[0.0, 0.0, r_q]]), np.array([1.0])
    x, w = np.polynomial.legendre.leggauss(n_psi)
    c0 = np.cos(radius / r_q)
    u = 0.5 * (1.0 - c0) * x + 0.5 * (1.0 + c0)
    wu = 0.5 * (1.0 - c0) * w
    phi = (np.arange(n_phi) + 0.5) * 2.0 * np.pi / n_phi
    sin_psi = np.sqrt(1.0 - u**2)
    pos = r_q * np.stack([np.outer(sin_psi, np.cos(phi)), np.outer(sin_psi, np.sin(phi)),
                          np.outer(u, np.ones(n_phi))], axis=-1).reshape(-1, 3)
    weights = np.repeat(wu, n_phi) / wu.sum() / n_phi
    return pos, weights


def patch_moment(r_q: float, radius: float, density_nAm_mm2: float | None) -> float:
    """Total moment [A m]: 30 nAm, or current density x area."""
    if density_nAm_mm2 is None or radius == 0.0:
        return R.Q_DIPOLE
    return density_nAm_mm2 * 1e-9 * patch_area(r_q, radius) * 1e6


def radial_field(sensors: np.ndarray, pos: np.ndarray, moments: np.ndarray) -> np.ndarray:
    """Radial field [T] at sensors (m, 3) of +y dipoles at pos (n, 3): Biot-Savart,
    B_r = mu0/4pi sum_k m_k (y x (r - p_k)) . e_r / |r - p_k|^3 (volume currents add no radial
    field in a sphere).  Returns (m, n) per-dipole fields weighted by the moments."""
    e_r = sensors / np.linalg.norm(sensors, axis=1, keepdims=True)
    a = sensors[:, None, :] - pos[None, :, :]
    num = a[..., 2] * e_r[:, None, 0] - a[..., 0] * e_r[:, None, 2]
    return R.MU0_OVER_4PI * num / np.linalg.norm(a, axis=2) ** 3 * moments[None, :]


def peak_field_analytic(r_q: float, radius: float, s: float, moment: float = R.Q_DIPOLE) -> float:
    """Peak |B_r| of a patch on the great circle y = 0 of radius s (Biot-Savart superposition,
    grid search + bounded refinement).  Used to refine equal-SNR crossings."""
    pos, w = cap_quadrature(r_q, radius)
    m = moment * w

    def neg_abs_br(theta):
        th = np.atleast_1d(theta)
        sens = s * np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], axis=1)
        return -np.abs(radial_field(sens, pos, m).sum(axis=1))

    grid = np.deg2rad(np.arange(0.25, 180.0, 0.5))
    k = int(np.argmin(neg_abs_br(grid)))
    step = np.deg2rad(0.5)
    res = minimize_scalar(lambda t: neg_abs_br(t)[0], bounds=(grid[k] - step, grid[k] + step),
                          method="bounded", options={"xatol": 1e-12})
    return float(-res.fun)


# =========================================================================================
# MNE-Python simulation
# =========================================================================================
def simulate_patches_mne(r_q: np.ndarray, radius: float, sensor_radius: float,
                         density: float | None = None, dtheta_deg: float = DTHETA_DEG) -> dict:
    """Peak |B_r| [T] of the patch of geodesic radius `radius` at each r_Q (NaN where the
    patch is undefined), simulated with MNE-Python; plus the largest sensor-wise deviation
    of the MNE field from the Biot-Savart superposition."""
    info, alpha = R.make_arc_info(sensor_radius, dtheta_deg)
    sensors = sensor_radius * np.column_stack([np.sin(alpha), np.zeros_like(alpha), np.cos(alpha)])
    sphere = R.make_head_sphere()
    on_inside = "ignore" if sensor_radius <= R.HEAD_RADIUS else "raise"
    peaks = np.full(r_q.size, np.nan)
    defined = [i for i, r in enumerate(r_q) if patch_defined(r, radius)]
    patches = {i: cap_quadrature(r_q[i], radius) for i in defined}
    moments = {i: patch_moment(r_q[i], radius, density) for i in defined}
    forward_err = 0.0
    batch: list[int] = []
    for n, i in enumerate(defined):
        batch.append(i)
        n_dip = sum(patches[j][0].shape[0] for j in batch)
        if n == len(defined) - 1 or n_dip + patches[defined[n + 1]][0].shape[0] > MAX_DIPOLES_PER_RUN:
            # shrink by 1e-12 (field change ~1e-12): nodes at r_Q = b can round 1 ulp outside
            # MNE's inner sphere, which make_forward_dipole rejects
            pos = np.concatenate([patches[j][0] for j in batch]) * (1.0 - 1e-12)
            amp = np.concatenate([moments[j] * patches[j][1] for j in batch])
            times = np.concatenate([np.full(patches[j][0].shape[0], k * 1e-3) for k, j in enumerate(batch)])
            dip = mne.Dipole(times=times, pos=pos, amplitude=amp, ori=np.tile([0.0, 1.0, 0.0], (len(amp), 1)),
                             gof=np.full(len(amp), 100.0), verbose=False)
            fwd, stc = mne.make_forward_dipole(dip, sphere, info, on_inside=on_inside, verbose=False)
            b_r = mne.simulation.simulate_evoked(fwd, stc, info, cov=None, nave=np.inf, verbose=False).data
            # Biot-Savart reference: per-element fields summed per patch
            member = np.zeros((len(amp), len(batch)))
            member[np.arange(len(amp)), np.round(times * 1e3).astype(int)] = 1.0
            reference = radial_field(sensors, pos, amp) @ member
            scale = np.abs(reference).max(axis=0)
            forward_err = max(forward_err, float(np.max(np.abs(b_r - reference) / scale)))
            peaks[batch] = R._parabolic_peak(np.abs(b_r))
            batch = []
    return dict(peak=peaks, forward_rel_err=forward_err)


# =========================================================================================
# Study
# =========================================================================================
@dataclass
class PatchResult:
    radius_mm: float
    depth_mm: np.ndarray  # curve samples (as in Fig. 3), deep -> superficial
    signal_opm_pT: np.ndarray
    signal_squid_pT: np.ndarray
    snr_opm: np.ndarray
    snr_squid: np.ndarray
    moment_at_surface_nAm: float
    area_at_surface_cm2: float
    deepest_depth_mm: float  # deepest depth at which the patch exists
    d_eq_mm: float  # equal-SNR depth at the study's eta (shallowest crossing)
    crossings_mm: list  # all crossings at the study's eta
    d_eq_vs_eta_mm: np.ndarray  # on ETA_GRID; 15 mm = no OPM advantage, deepest = everywhere
    forward_rel_err: float
    peak_rel_err: float  # MNE (interpolated) peak vs refined Biot-Savart peak
    extras: dict = field(default_factory=dict)


def _ratio_function(radius: float, density: float | None):
    def ratio(d_mm: float) -> float:
        r_q = H - d_mm * 1e-3
        if not patch_defined(r_q, radius):
            return np.nan
        m = patch_moment(r_q, radius, density)
        return (peak_field_analytic(r_q, radius, SENSOR_RADII["OPM"], m)
                / peak_field_analytic(r_q, radius, SENSOR_RADII["SQUID"], m))
    return ratio


def _crossings(depths: np.ndarray, ratios: np.ndarray, eta: float, ratio) -> list[float]:
    """Depths where the OPM/SQUID signal ratio equals eta (sign changes of ratio - eta,
    refined with brentq)."""
    out = []
    g = ratios - eta
    for k in np.flatnonzero(np.isfinite(g[:-1]) & np.isfinite(g[1:]) & (np.sign(g[:-1]) != np.sign(g[1:]))):
        root = brentq(lambda d: ratio(d) - eta, depths[k], depths[k + 1], xtol=1e-9)
        if not out or abs(root - out[-1]) > 1e-6:  # a root on a grid node is bracketed twice
            out.append(root)
    return out


def run_study(radii=PATCH_RADII_MM, eta: float = R.ETA, density: float | None = None) -> tuple[list, dict]:
    r_curve = np.linspace(0.0, B, R.N_CURVE)[1:]  # as in Fig. 3: 0.8 ... 80 mm
    depth_mm = (H - r_curve) * 1e3
    # noise: identical to Fig. 3 (SQUID signal of the 30 nAm reference dipole at 0.8 b)
    ref = R.simulate_peak_field_mne(np.array([R.NOISE_REF_RADIUS]), SENSOR_RADII["SQUID"], R.DTHETA_DEG)
    sigma_squid = float(ref["peak"][0])
    sigma_opm = eta * sigma_squid
    fine_depths = np.arange((H - B) * 1e3, H * 1e3 - 1e-9, FINE_DEPTH_STEP_MM)
    results = []
    for radius_mm in radii:
        radius = radius_mm * 1e-3
        sims = {k: simulate_patches_mne(r_curve, radius, s, density) for k, s in SENSOR_RADII.items()}
        b_opm, b_squid = sims["OPM"]["peak"], sims["SQUID"]["peak"]
        ratio = _ratio_function(radius, density)
        # MNE (interpolated) peaks vs refined analytic peaks on a subset of depths
        idx = [i for i in range(0, r_curve.size, 9) if np.isfinite(b_opm[i])]
        peak_err = max((abs(b_opm[i] / peak_field_analytic(r_curve[i], radius, SENSOR_RADII["OPM"],
                                                           patch_moment(r_curve[i], radius, density)) - 1)
                        for i in idx), default=0.0)
        fwd_err = max(s["forward_rel_err"] for s in sims.values())
        tol = min(1e-6 + 1e-3 * DTHETA_DEG**3, 1e-4)  # as replicate_figure3.run_simulation
        if fwd_err > 1e-5 or peak_err > tol:
            raise RuntimeError(f"rho = {radius_mm:g} mm: MNE vs Biot-Savart field {fwd_err:.1e}, "
                               f"peak {peak_err:.1e} (tolerances 1e-5, {tol:.1e})")
        # equal-SNR depth(s): crossings of the signal ratio with eta, from the brain surface down
        fine_ratio = np.array([ratio(d) for d in fine_depths])
        crossings = _crossings(fine_depths, fine_ratio, eta, ratio)
        defined = np.isfinite(fine_ratio)
        deepest = float(fine_depths[defined][-1])
        # d_eq(eta): shallowest crossing (linear interpolation on the 0.25-mm grid); 15 mm if
        # there is no OPM advantage at all, the deepest depth if OPM wins wherever the patch exists
        d_vs_eta = np.empty(ETA_GRID.size)
        for k, e in enumerate(ETA_GRID):
            below = np.flatnonzero(defined & (fine_ratio < e))
            if fine_ratio[0] <= e:
                d_vs_eta[k] = fine_depths[0]
            elif below.size == 0:
                d_vs_eta[k] = deepest
            else:
                j = below[0]
                g0, g1 = fine_ratio[j - 1] - e, fine_ratio[j] - e
                d_vs_eta[k] = fine_depths[j - 1] + (fine_depths[j] - fine_depths[j - 1]) * g0 / (g0 - g1)
        r_surface = H - fine_depths[0] * 1e-3
        results.append(PatchResult(
            radius_mm=float(radius_mm), depth_mm=depth_mm,
            signal_opm_pT=b_opm * 1e12, signal_squid_pT=b_squid * 1e12,
            snr_opm=b_opm / sigma_opm, snr_squid=b_squid / sigma_squid,
            moment_at_surface_nAm=patch_moment(r_surface, radius, density) * 1e9,
            area_at_surface_cm2=patch_area(r_surface, radius) * 1e4,
            deepest_depth_mm=deepest,
            d_eq_mm=float(crossings[0]) if crossings else np.nan, crossings_mm=[float(c) for c in crossings],
            d_eq_vs_eta_mm=d_vs_eta,
            forward_rel_err=fwd_err, peak_rel_err=float(peak_err),
            extras=dict(signal_ratio_at_surface=float(fine_ratio[0]),
                        snr_ratio_at_surface=float(fine_ratio[0] / eta)),
        ))
    info = dict(sigma_squid_pT=sigma_squid * 1e12, sigma_opm_pT=sigma_opm * 1e12, eta=eta,
                moment_density_nAm_per_mm2=density, fine_depths_mm=fine_depths)
    return results, info


# =========================================================================================
# Figures
# =========================================================================================
def _as_fig3_results(p: PatchResult, info: dict) -> R.Results:
    return R.Results(
        depth_mm=p.depth_mm, signal_opm_pT=p.signal_opm_pT, signal_squid_pT=p.signal_squid_pT,
        snr_opm=p.snr_opm, snr_squid=p.snr_squid, sigma_opm_pT=info["sigma_opm_pT"],
        sigma_squid_pT=info["sigma_squid_pT"], d_eq_grid_mm=np.nan, d_eq_exact_mm=p.d_eq_mm,
        dtheta_deg=DTHETA_DEG, max_rel_err_forward=p.forward_rel_err, max_rel_err_peak=p.peak_rel_err,
        max_lobe_asymmetry=0.0)


def _label(p: PatchResult) -> str:
    if p.radius_mm == 0:
        return "dipole"
    return f"ρ = {p.radius_mm:g} mm ({p.area_at_surface_cm2:.1f} cm²)"


def figure3_per_patch(p: PatchResult, info: dict, fonts, outdir: Path) -> list[str]:
    """The Figure 3 layout for one patch radius; the dotted line is its equal-SNR depth."""
    with plt.rc_context(R.FIGURE_RC):
        fig, axes, _ = R.make_figure(_as_fig3_results(p, info), fonts, deq_marker="exact", annotations=False,
                                     source_word="Dipole" if p.radius_mm == 0 else "Patch")
        for ax in axes.values():  # large moments (--moment-density): at most ~7 integer ticks
            top = ax.get_ylim()[1]
            if top >= 8:
                ticks = MaxNLocator(nbins=6, integer=True).tick_values(0.0, top)
                ax.set_yticks(ticks[ticks <= top])
                for lab in ax.get_yticklabels():
                    lab.set_fontproperties(fonts.tick)
                    lab.set_path_effects([patheffects.withStroke(linewidth=R.TICK_LABEL_STROKE_PT,
                                                                 foreground="k")])
        # the published layout leaves room left of panel A for one-digit tick labels only:
        # widen the canvas to the left when longer labels stick out
        fig.canvas.draw()
        x0 = min(lab.get_window_extent().x0 for lab in axes["A"].get_yticklabels() if lab.get_text())
        extra_left_in = -x0 / fig.dpi + 0.02 if x0 < -20 else 0.0
        canvas = R.Canvas(fig)
        x0 = 3330.0  # free area of panel B (published layout coordinates, px)
        if p.radius_mm == 0:
            lines = ["Current dipole", f"Q = {p.moment_at_surface_nAm:.0f} nAm"]
        else:
            lines = [f"Patch ρ = {p.radius_mm:g} mm",
                     f"{p.area_at_surface_cm2:.1f} cm², {p.moment_at_surface_nAm:.0f} nAm"]
        for k, line in enumerate(lines):
            canvas.text(x0, 700.0 + 200.0 * k, line, fonts.regular)
        sub = fonts.regular.copy()
        sub.set_size(R.FONT_SIZE * R.SUB_SCALE)
        value = f" = {p.d_eq_mm:.1f} mm" if np.isfinite(p.d_eq_mm) else ": none"
        canvas.typeset(x0, 1100.0, [("d", fonts.italic_serif, 0.0),
                                    ("eq", sub, R.SUB_SHIFT * R.FONT_SIZE * R.PX_PER_PT),
                                    (value, fonts.regular, 0.0)])
        paths = []
        for ext in ("png", "pdf"):
            path = outdir / f"Figure3_patch_rho_{p.radius_mm:02g}mm.{ext}"
            bbox = Bbox.from_extents(-extra_left_in, 0.0, *fig.get_size_inches()) if extra_left_in else None
            fig.savefig(path, dpi=R.PUB_PPI if ext == "png" else 300, facecolor="white", bbox_inches=bbox)
            paths.append(str(path))
        plt.close(fig)
    return paths


def summary_figure(results: list, info: dict, fonts, outdir: Path) -> list[str]:
    """(A) signal and (B) SNR vs depth for all patch radii (OPM blue solid, SQUID red dashed;
    darker = smaller patch), (C) SNR_OPM / SNR_SQUID vs depth, (D) equal-SNR depth vs eta."""
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D

    n = len(results)
    shade = np.linspace(0.95, 0.45, n)
    blues, reds, greys = plt.cm.Blues(shade), plt.cm.Reds(shade), plt.cm.Greys(shade)
    fp = fonts.regular.copy()
    small = fonts.regular.copy()
    small.set_size(7)
    family = fonts.regular.get_name()
    serif = fonts.italic_serif.get_name()
    math_rc = {"mathtext.fontset": "custom", "mathtext.rm": family, "mathtext.bf": f"{family}:bold",
               "mathtext.it": f"{serif}:italic", "mathtext.default": "it"}
    with plt.rc_context(dict(R.FIGURE_RC, **math_rc, **{"axes.linewidth": 0.8, "xtick.major.width": 0.8,
                                                         "ytick.major.width": 0.8, "xtick.labelsize": 8,
                                                         "ytick.labelsize": 8})):
        fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.2), dpi=300)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.08, wspace=0.28, hspace=0.40)
        (ax_a, ax_b), (ax_c, ax_d) = axs
        handles = []
        for k, p in enumerate(results):
            h_opm, = ax_a.plot(p.depth_mm, p.signal_opm_pT, "-", color=blues[k], lw=1.3)
            h_sq, = ax_a.plot(p.depth_mm, p.signal_squid_pT, "--", color=reds[k], lw=1.3)
            handles.append((h_opm, h_sq))
            ax_b.plot(p.depth_mm, p.snr_opm, "-", color=blues[k], lw=1.3)
            ax_b.plot(p.depth_mm, p.snr_squid, "--", color=reds[k], lw=1.3)
            ax_c.plot(p.depth_mm, p.snr_opm / p.snr_squid, "-", color=greys[k], lw=1.3)
            ax_d.plot(ETA_GRID, p.d_eq_vs_eta_mm, "-", color=greys[k], lw=1.3)
            if np.isfinite(p.d_eq_mm):
                ok = np.isfinite(p.snr_opm)
                y = np.interp(p.d_eq_mm, p.depth_mm[ok][::-1], p.snr_opm[ok][::-1])
                ax_b.plot([p.d_eq_mm], [y], "o", mfc="white", mec="k", ms=3.5, mew=0.8, zorder=5)
                ax_c.plot([p.d_eq_mm], [1.0], "o", mfc="white", mec="k", ms=3.5, mew=0.8, zorder=5)
                ax_d.plot([info["eta"]], [p.d_eq_mm], "o", mfc="white", mec="k", ms=3.5, mew=0.8, zorder=5)
        for ax in (ax_a, ax_b, ax_c):
            ax.set_xlim(0, 95)
            ax.set_xticks([0, 20, 40, 60, 80])
            ax.set_xlabel("Source depth ($d$) [mm]", fontproperties=fp)
        ax_a.set_ylim(bottom=0)
        ax_b.set_ylim(bottom=0)
        ax_c.set_yscale("log")
        ax_c.axhline(1.0, color="k", ls=":", lw=1.0)
        ratios = np.concatenate([p.snr_opm / p.snr_squid for p in results])
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        c_lo, c_hi = min(0.3, ratios.min() / 1.1), max(2.5, ratios.max() * 1.1)
        ax_c.set_ylim(c_lo, c_hi)
        c_ticks = [t for t in (0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10) if c_lo <= t <= c_hi]
        ax_c.set_yticks(c_ticks)
        ax_c.set_yticklabels([f"{t:g}" for t in c_ticks])
        ax_c.minorticks_off()
        ax_d.set_xlim(1.0, 6.0)
        ax_d.set_ylim(95, 0)
        ax_d.axvline(info["eta"], color="k", ls=":", lw=1.0)
        ax_d.axhline((H - B) * 1e3, color="k", ls="-.", lw=0.8)
        ax_d.set_xlabel(r"Relative noise level ($\eta$)", fontproperties=fp)
        ax_d.set_ylabel(r"$d_\mathrm{eq}$ [mm]", fontproperties=fp)
        k0 = int(np.argmin(np.abs(ETA_GRID - 1.35)))
        plateaus = [p for p in results if np.isclose(p.d_eq_vs_eta_mm[k0], p.deepest_depth_mm)]
        if plateaus:  # point at the middle curve's "patch exists only down to here" plateau
            y = plateaus[len(plateaus) // 2].deepest_depth_mm
            ax_d.annotate("OPM better\nat all depths", xy=(ETA_GRID[k0], y), xytext=(2.1, y), va="center",
                          fontproperties=small, color="0.35",
                          arrowprops=dict(arrowstyle="->", lw=0.6, color="0.35", shrinkA=1, shrinkB=0))
        ax_d.text(5.95, 12.5, "no OPM advantage", fontproperties=small, ha="right", va="bottom", color="0.35")
        titles = {ax_a: ("A", "Signal [pT]"), ax_b: ("B", "Signal-to-Noise Ratio"),
                  ax_c: ("C", r"SNR$_\mathrm{OPM}$ / SNR$_\mathrm{SQUID}$  ($\eta$ = " + f"{info['eta']:g})"),
                  ax_d: ("D", "Equal-SNR depth")}
        for ax, (letter, title) in titles.items():
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.text(-0.13, 1.04, letter, transform=ax.transAxes, fontproperties=fonts.bold, va="bottom")
            ax.text(0.0, 1.04, title, transform=ax.transAxes, fontproperties=fp, va="bottom")
        ax_c.text(0.98, 0.96, "OPM better", transform=ax_c.transAxes, ha="right", va="top", fontproperties=small)
        ax_c.text(0.98, 0.04, "SQUID better", transform=ax_c.transAxes, ha="right", va="bottom",
                  fontproperties=small)
        ax_a.legend(handles, [_label(p) for p in results], prop=small, frameon=False, loc="upper right",
                    handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)}, handlelength=3.4,
                    title="OPM (solid) / SQUID (dashed)\narea: at d = 15 mm", title_fontproperties=small)
        ax_d.legend([Line2D([], [], color=greys[k], lw=1.3) for k in range(n)], [_label(p) for p in results],
                    prop=small, frameon=False, loc="center right")
        moment = ("Q = 30 nAm per patch" if info["moment_density_nAm_per_mm2"] is None
                  else f"{info['moment_density_nAm_per_mm2']:g} nAm/mm\u00b2 (dipole: 30 nAm)")
        circles = (f"; circles: equal-SNR depth at \u03b7 = {info['eta']:g}"
                   if any(np.isfinite(p.d_eq_mm) for p in results) else "")
        fig.text(0.99, 0.995, f"Spherical-cap patches, uniform current (tangential at the centre), {moment}{circles}",
                 ha="right", va="top", fontproperties=small, color="0.35")
        paths = []
        for ext in ("png", "pdf", "svg"):
            path = outdir / f"Figure_extended_sources.{ext}"
            fig.savefig(path, dpi=300, facecolor="white")
            paths.append(str(path))
        plt.close(fig)
    return paths


# =========================================================================================
# Driver
# =========================================================================================
def _json_safe(value):
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main(argv=None):
    ap = argparse.ArgumentParser(description="Figure 3 of Jas et al. (2026) with extended sources "
                                             "(cortical patches of different sizes), simulated with MNE-Python.")
    ap.add_argument("--radii", type=float, nargs="+", default=list(PATCH_RADII_MM),
                    help="patch radii [mm] (geodesic, on the sphere at the source depth); 0 = dipole")
    ap.add_argument("--eta", type=R._positive_float, default=R.ETA)
    ap.add_argument("--moment-density", type=R._positive_float, default=None,
                    help="current density [nAm/mm^2]; default: every patch has Q = 30 nAm")
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent / "extended_sources_output"))
    ap.add_argument("--artwork-fonts", action="store_true",
                    help="licensed artwork fonts (Myriad Pro; see replicate_figure3.py); default: open fonts")
    ap.add_argument("--font-dir", default=None, help="folder containing MyriadPro-*.otf (implies --artwork-fonts)")
    args = ap.parse_args(argv)
    max_radius_mm = 0.5 * np.pi * B * 1e3  # a hemispherical cap on the brain surface
    if any(not 0 <= r < max_radius_mm for r in args.radii):
        ap.error(f"patch radii must be in [0, {max_radius_mm:.1f}) mm (hemisphere of the brain surface)")
    if any(r > 40 for r in args.radii):
        print(f"warning: the {N_PSI} x {N_PHI} patch quadrature loses accuracy above 40 mm "
              "(peak error ~4e-5 at 60 mm, ~1e-3 at 100 mm)")

    mne.set_log_level("WARNING")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results, info = run_study(sorted(set(args.radii)), args.eta, args.moment_density)
    fonts = R.setup_fonts(args.font_dir, args.artwork_fonts)
    files = []
    for p in results:
        files += figure3_per_patch(p, info, fonts, outdir)
    files += summary_figure(results, info, fonts, outdir)

    with open(outdir / "extended_sources_data.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["patch_radius_mm", "depth_mm", "signal_OPM_pT", "signal_SQUID_pT", "SNR_OPM", "SNR_SQUID"])
        for p in results:
            order = np.argsort(p.depth_mm)
            for row in zip(p.depth_mm[order], p.signal_opm_pT[order], p.signal_squid_pT[order],
                           p.snr_opm[order], p.snr_squid[order]):
                wr.writerow([f"{p.radius_mm:g}"] + ["" if not np.isfinite(v) else f"{v:.10g}" for v in row])
    summary = dict(
        model="spherical-cap patches, uniform +y current (tangential at the cap centre), depth = cap depth",
        deepest_depth_note="deepest depth on the 0.25-mm crossing grid at which the patch exists",
        sigma_squid_pT=info["sigma_squid_pT"], sigma_opm_pT=info["sigma_opm_pT"], eta=info["eta"],
        moment=("30 nAm per patch" if args.moment_density is None
                else f"{args.moment_density} nAm/mm^2 x area (dipole: 30 nAm)"),
        quadrature=f"{N_PSI} x {N_PHI} dipoles per patch", sensor_spacing_deg=DTHETA_DEG,
        mne_version=mne.__version__,
        patches=[dict(radius_mm=p.radius_mm, area_at_surface_cm2=p.area_at_surface_cm2,
                      moment_at_surface_nAm=round(p.moment_at_surface_nAm, 9), deepest_depth_mm=p.deepest_depth_mm,
                      d_eq_mm=p.d_eq_mm, crossings_mm=p.crossings_mm,
                      signal_OPM_pT_at_surface=float(p.signal_opm_pT[-1]),
                      signal_SQUID_pT_at_surface=float(p.signal_squid_pT[-1]),
                      snr_OPM_at_surface=float(p.snr_opm[-1]), snr_SQUID_at_surface=float(p.snr_squid[-1]),
                      eta_range_with_crossing=[float(np.nanmin(p.snr_opm / p.snr_squid) * info["eta"]),
                                               p.extras["signal_ratio_at_surface"]],
                      forward_rel_err=p.forward_rel_err, peak_rel_err=p.peak_rel_err,
                      d_eq_vs_eta_mm=dict(zip(ETA_GRID.tolist(), p.d_eq_vs_eta_mm.tolist())))
                 for p in results],
        outputs=files)
    with open(outdir / "extended_sources_summary.json", "w") as fh:
        json.dump(_json_safe(summary), fh, indent=2, allow_nan=False)

    print(f"MNE-Python {mne.__version__}; noise sigma_SQUID = {info['sigma_squid_pT']:.4f} pT, "
          f"sigma_OPM = {info['sigma_opm_pT']:.4f} pT (eta = {info['eta']:g})")
    print(f"{'patch':>28s} | {'B_OPM/B_SQUID at d=15':>21s} | {'SNR_OPM/SNR_SQUID':>17s} | "
          f"{'d_eq [mm]':>9s} | {'exists to d':>11s} | MNE vs Biot-Savart")
    for p in results:
        print(f"{_label(p):>28s} | {p.extras['signal_ratio_at_surface']:21.3f} | "
              f"{p.extras['snr_ratio_at_surface']:17.3f} | {p.d_eq_mm:9.2f} | {p.deepest_depth_mm:9.1f} mm | "
              f"{max(p.forward_rel_err, p.peak_rel_err):.1e}")
    print("Wrote:", outdir)
    return results, info


if __name__ == "__main__":
    main()
