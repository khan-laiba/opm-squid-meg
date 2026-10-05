#!/usr/bin/env python3
"""Report figures in the manuscript's own terms (round-1 referee reports: no internal analysis
labels in the panels, one dB colour scale for the cortical maps, a legible geometry figure), drawn
from stored outputs only: nothing is simulated or re-analysed.

  R0   the spherical benchmark: (a) equal-SNR depth against the noise ratio eta in the sphere of
       Jas et al. (2026) (our reimplementation, results/g1a), (b) the same quantity on the realistic
       adult head for the matched and dense OPM arrays and the sphere with the real standoffs
       (results/g2/g2_summary.json :: bridge_to_sphere)
  R11  D (dB) on the adult's inflated cortex: dense and matched OPM arrays against Neuromag's 306
       channels, sensor + brain noise (results/g2/g2_targets.csv)
  R12  geometry: sagittal and coronal scalp sections of the adult, the 12-month template and child B,
       each with the fixed adult helmet at top contact, the helmet scaled with the head (laterally
       centred) and the dense OPM sites (results/g3b/g3b_geometry_sections.json)
  R13  D (dB) on the inflated cortex of the 24- and 12-month templates, dense array, top contact
       (results/g3b/g3b_targets_<anatomy>.csv), on the colour scale of R11
  R14  D (dB) on the adult's inflated cortex for the adult and the two scaled adults (school-age and
       2-year size), dense array, top contact (results/g3b/g3b_targets_<anatomy>.csv), same scale

D = 20 log10(detectability OPM / detectability Neuromag) per cortical target (the G2 and G3B
definition; scripts/g3b_pediatric_helmet.py d_db). Surfaces are read from the external data for
drawing the cortical maps only (OPMSQUID_DATA, src/opmsquid/paths.py): inflated white surfaces and
sulcal depth (shading) at the vertices of the oct-6 source spaces, of which the targets are a subset.
R12 is drawn from results/ only: its scalp and white-surface sections, coil centres and OPM sites are
those exported from the stored state of the G3B run by scripts/export_g3b_geometry.py, and the script
stops unless that export is of the run of results/g3b/g3b_summary.json and reproduces its moves,
scale factors, coil-to-scalp gaps and OPM site counts. A missing input stops the script before
anything is drawn: nothing is recomputed in its place.
Every number in the labels, captions and descriptions is read from these files or computed from them
as described in figures_clean.json; colour limits, section planes and the 20-mm sensor slab are
display choices (the slab is that of the G3B geometry figure).

Outputs (results/report/): Figure_R0_sphere.png, Figure_R11_maps_adult.png, Figure_R12_geometry.png,
Figure_R13_maps_heads.png, Figure_R14_maps_scaled.png and figures_clean.json (inputs, description,
alt text, draft caption and plotted values of every figure).

Usage: OPMSQUID_DATA=<data dir> PYTHONPATH=src .venv/bin/python scripts/report_figures_clean.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend before pyplot is used)
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from opmsquid import anatomy, io, paths, pediatric as P, plotting  # noqa: E402

G1A = "results/g1a/g1a_benchmark.json"
G1A_CURVES = "results/g1a/g1a_curves.csv"
G2 = "results/g2/g2_summary.json"
G2_TARGETS = "results/g2/g2_targets.csv"
G3B = "results/g3b/g3b_summary.json"
G3B_TARGETS = "results/g3b/g3b_targets_{}.csv"
GEOMETRY = "results/g3b/g3b_geometry_sections.json"  # scripts/export_g3b_geometry.py
OUT_JSON = "figures_clean.json"

COND = "intrinsic+brain"  # sensor (intrinsic) + brain noise: the primary condition
DB_PER_LOG2 = 20.0 * np.log10(2.0)
TEMPLATES = {"infant2yr": "ANTS2-0Years3T", "infant12mo": "ANTS12-0Months3T"}  # as scripts/g3b_pediatric_helmet.py
SCALED_HEADS = ("adult", "school", "size2yr")  # R14: the scaled adults keep the adult's vertices
GEOMETRY_HEADS = {"adult": None, "infant12mo": "ANTS12-0Months3T", "childB": "sub-Z209"}
GEOMETRY_NAMES = {"adult": "adult", "infant12mo": "12-month template", "childB": "child B"}

# cortical maps: one diverging dB scale for R11 and R13 (display choice: +-6 dB, about a factor of 2 in
# detectability, the range of the earlier dB maps), two greys for what carries no value
LIM_DB = 6.0
CMAP = plt.get_cmap("RdBu_r")
DATA, MEDIAL, EXCLUDED = 0, 1, 2
GREY = {MEDIAL: 0.64, EXCLUDED: 0.36}  # both darker than the palest colours of the scale, so not mistaken for 0 dB
SHADE = 0.12  # sulcal shading of the coloured cortex (opmsquid.plotting.render_view uses 0.2)
NEAR_SKULL_MM = anatomy.MIN_BEM_DISTANCE * 1e3  # targets keep this distance from the inner-skull mesh
VIEWS = ((0, -1.0, "Left, lateral"), (0, 1.0, "Left, medial"), (1, -1.0, "Right, medial"), (1, 1.0, "Right, lateral"))

# geometry: sections through the head origin, sensors within the slab (as the G3B geometry figure)
SLAB_MM = 20.0
ROUND_MM = 0.01  # rounding of the exported coordinates (scripts/export_g3b_geometry.py)
COIL_TOL_MM = 0.02  # drawn coil centre vs stored transform applied to the device-frame coil centre (both rounded)
SCALED_COLOR = "#D55E00"
CORTEX_GREY = "0.62"
SECTION_LIM = {"sagittal": ((-138.0, 152.0), (-68.0, 182.0)), "coronal": ((-145.0, 145.0), (-68.0, 182.0))}

OPM_FILL = "#56B4E9"
PRINTED_COLOR = "#D55E00"
MARKER = {"opm_dense": "o", "opm_matched": "s"}


def load(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text())


def col(rows: list[dict], key: str) -> np.ndarray:
    return np.array([float(r[key]) for r in rows])


def fmt_db(x: float) -> str:
    return f"{x:+.2f}".replace("-", "\u2212")


# ----------------------------------------------------------------------------------------------
# R0: the spherical benchmark and the realistic head
def figure_sphere() -> dict:
    g, s = load(G1A), load(G2)
    p, br = g["parameters"], s["bridge_to_sphere"]
    rows = io.read_csv(ROOT / G1A_CURVES)
    depth = col(rows, "depth_mm")
    ratio = col(rows, "B_OPM_pT") / col(rows, "B_SQUID_pT")  # SNR_OPM = SNR_SQUID <=> B_OPM / B_SQUID = eta
    if not np.all(np.diff(ratio) < 0):
        raise SystemExit(f"{G1A_CURVES}: the field ratio does not fall monotonically with depth")
    eta0, eta1 = g["eta_range_with_crossing"]
    exact = {float(e): float(v) for e, v in g["d_eq_mm"].items() if v is not None}

    def d_eq(e):  # depth where the stored field ratio equals eta (log-linear between grid depths)
        return float(np.interp(np.log(e), np.log(ratio[::-1]), depth[::-1]))

    err = max(abs(d_eq(e) - v) for e, v in exact.items())
    if err > 0.01:
        raise SystemExit(f"the stored field ratio does not reproduce the stored exact roots ({err:.3f} mm)")
    printed = next(x for x in g["printed_vs_exact"] if x["item"] == "Fig. 3 d_eq (eta 3)")
    h, b = float(p["h_mm"]), float(p["b_mm"])
    brain_top = h - b
    eta_c = np.r_[eta0, ratio[::-1]]  # from the sphere centre (eta0, depth h) to the brain surface
    d_c = np.r_[h, depth[::-1]]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(style.FULL_W, 3.15), sharey=True, layout="constrained")
    ax.fill_between(np.r_[1.0, eta_c], brain_top, np.r_[h, d_c], color=OPM_FILL, alpha=0.28, lw=0, zorder=1)
    ax.plot(eta_c, d_c, "-", color="k", lw=1.4, zorder=3)
    ax.plot(list(exact), list(exact.values()), "o", color="k", ms=3.4, zorder=4)
    ax.plot([3.0], [printed["printed_mm"]], "D", ms=7.5, mfc="none", mec=PRINTED_COLOR, mew=1.3, zorder=5)
    ax.axhline(brain_top, color="0.5", lw=0.7, ls=(0, (1, 2)), zorder=2)
    ax.text(5.97, brain_top - 1.2, "brain surface", ha="right", va="bottom", fontsize=7, color="0.35")
    ax.text(1.55, 40, "OPM SNR\nhigher", ha="center", va="center", fontsize=8, color="#0b4f75")
    ax.text(4.7, 52, "SQUID SNR\nhigher", ha="center", va="center", fontsize=8, color="0.25")
    ax.set_title("(a) Spherical head", loc="left")
    ax.set_xlabel(r"Noise ratio $\eta = \sigma_{\mathrm{OPM}}\,/\,\sigma_{\mathrm{SQUID}}$")
    ax.set_ylabel("Equal-SNR depth below the scalp (mm)")
    h_a = [Line2D([], [], color="k", lw=1.4, marker="o", ms=3.4,
                  label=f"OPM on the scalp, SQUID {p['xi_squid_mm']:g} mm above it\n(dots: exact roots)"),
           Line2D([], [], ls="none", marker="D", ms=7.5, mfc="none", mec=PRINTED_COLOR, mew=1.3,
                  label=f"Value printed in the preprint\n({printed['printed_mm']:g} mm at \u03b7 = 3)")]
    ax.legend(handles=h_a, loc="lower right", handlelength=1.8)

    sph = {float(e): v for e, v in br["sphere_d_eq_mm"]["realistic_standoffs"].items() if v is not None}
    xi_opm, xi_sq = br["sensor_distance_mm"]["opm_matched"]["median"], br["sensor_distance_mm"]["squid"]["median"]
    bx.plot(list(sph), list(sph.values()), "--", color="k", lw=1.3, zorder=3,
            label=f"Sphere with the real standoffs\n(OPM {xi_opm:.1f} mm, Neuromag {xi_sq:.1f} mm)")
    arr = {}
    for a in ("opm_matched", "opm_dense"):
        v = {float(e): x for e, x in br[f"{a}_d_eq_mm"].items() if x is not None}
        arr[a] = v
        bx.plot(list(v), list(v.values()), "-", color=style.ARRAY_COLOR[a], marker=MARKER[a], ms=3.6, lw=1.4, zorder=4,
                label=f"{style.ARRAY_LABEL[a]} ({s['arrays'][a]['n_sites']} sites), realistic head")
    ahead = {a: br[f"{a}_eta_opm_ahead_at_all_depths"] for a in arr}
    behind = {a: br[f"{a}_eta_squid_ahead_at_all_depths"] for a in arr}
    if ahead["opm_dense"] != ahead["opm_matched"]:
        raise SystemExit("the arrays differ in the eta range where the OPM leads at every depth: reword the annotation")
    bx.text(1.08, 22, f"\u03b7 \u2264 {max(ahead['opm_dense']):g}:\nOPM ahead at every\ndepth (both arrays)", ha="left",
            va="center", fontsize=7, color="0.3")
    bx.text(5.95, 30, f"Neuromag ahead at every depth:\nmatched \u03b7 \u2265 {min(behind['opm_matched']):g}, "
            f"dense \u03b7 \u2265 {min(behind['opm_dense']):g}", ha="right", va="center", fontsize=7, color="0.3")
    bx.set_title("(b) Realistic adult head", loc="left")
    bx.set_xlabel(r"Noise ratio $\eta = \sigma_{\mathrm{OPM}}\,/\,\sigma_{\mathrm{magnetometer}}$")
    bx.legend(loc="lower right", handlelength=2.2)
    bx.tick_params(labelleft=True)
    for z in (ax, bx):
        z.set_xlim(1.0, 6.0)
        z.set_ylim(97.0, 10.0)
        z.set_xticks(np.arange(1, 7))
    style.save(fig, "Figure_R0_sphere")

    n = {a: s["arrays"][a]["n_sites"] for a in arr}
    rng_txt = {a: f"eta {min(v):g}-{max(v):g}: {min(v.values()):.1f}-{max(v.values()):.1f} mm" for a, v in arr.items()}
    values = dict(sphere_exact_roots_mm=exact, sphere_eta_range_with_crossing=[eta0, eta1],
                  sphere_curve_max_abs_err_vs_exact_mm=err, printed_eta3_mm=printed["printed_mm"],
                  exact_eta3_mm=exact[3.0], brain_surface_depth_mm=brain_top,
                  sphere_real_standoffs_d_eq_mm=sph, real_standoffs_mm=dict(opm=xi_opm, squid_mag=xi_sq),
                  arrays_d_eq_mm=arr, eta_opm_ahead_at_all_depths=ahead, eta_squid_ahead_at_all_depths=behind, n_sites=n)
    q = s["config"]["sources"]["focal_nAm"]
    bins = br["opm_dense_ratio_vs_depth"]
    bin_w = {r["hi"] - r["lo"] for r in bins}
    etas = sorted(float(e) for e in br["sphere_d_eq_mm"]["realistic_standoffs"])
    steps = set(np.round(np.diff(etas), 6))
    if len(bin_w) != 1 or len(steps) != 1:
        raise SystemExit(f"{G2}: unequal depth bins or eta steps in bridge_to_sphere")
    bin_w, step = bin_w.pop(), steps.pop()
    cross_lo = min(min(v.values()) for v in arr.values())
    cross_hi = max(max(v.values()) for v in arr.values())
    eta_lo, eta_hi = min(min(v) for v in arr.values()), max(max(v) for v in arr.values())
    return dict(
        inputs=[f"{G1A_CURVES} :: depth_mm, B_OPM_pT, B_SQUID_pT (peak radial fields on the preprint's grid)",
                f"{G1A} :: d_eq_mm (exact roots of the preprint's Eq. 3), eta_range_with_crossing, "
                "printed_vs_exact['Fig. 3 d_eq (eta 3)'], parameters (h_mm, b_mm, Q_nAm, xi_opm_mm, xi_squid_mm)",
                f"{G2} :: bridge_to_sphere.sphere_d_eq_mm.realistic_standoffs, bridge_to_sphere.<array>_d_eq_mm, "
                "bridge_to_sphere.<array>_eta_opm_ahead_at_all_depths, "
                "bridge_to_sphere.<array>_eta_squid_ahead_at_all_depths, bridge_to_sphere.<array>_ratio_vs_depth (bin "
                "edges), bridge_to_sphere.sensor_distance_mm (opm_matched.median, squid.median), arrays.<array>.n_sites, "
                "config.sources.focal_nAm"],
        description=(
            "Adapted from Jas et al. (2026), CC BY 4.0 (bioRxiv 2026.08.17.744953; the preprint's Fig. 4E, redrawn from "
            f"our reimplementation of its spherical model). (a) Sphere of the preprint: head radius {h:g} mm, brain radius "
            f"{b:g} mm, {p['Q_nAm']:g}-nAm tangential dipole, OPM on the scalp (standoff {p['xi_opm_mm']:g} mm), SQUID "
            f"{p['xi_squid_mm']:g} mm above it, single-channel SNR = peak radial field / noise. SNR_OPM = SNR_SQUID where "
            f"B_OPM / B_SQUID = eta, so the curve is the depth at which the stored field ratio of {G1A_CURVES} equals eta "
            f"(log-linear interpolation between its {len(depth)} grid depths, {depth.min():g}-{depth.max():g} mm), from the "
            f"sphere centre (eta0 = {eta0:.4f}, depth {h:g} mm) to the brain surface (eta1 = {eta1:.4f}, depth "
            f"{brain_top:g} mm); dots: the stored exact roots of Eq. 3 at eta " + ", ".join(f"{e:g}" for e in exact)
            + f" (curve within {err:.4f} mm of them). Shaded: depths where the OPM's SNR is higher (every depth for eta < "
            "eta0); dotted line: the brain surface. Open diamond: the equal-SNR depth the preprint prints for eta = 3 "
            f"({printed['printed_mm']:g} mm, with its Fig. 3; drawn there at {printed['drawn_mm']:g} mm; exact root "
            f"{exact[3.0]:.2f} mm). (b) Realistic adult head (MNE sample subject, 3-layer BEM, {q:g}-nAm cortical-normal "
            "dipoles), the G2 bridge: per target the peak |B| of the OPM array over that of the best Neuromag "
            f"magnetometer; medians over targets in {bin_w:g}-mm bins of depth below the scalp (distance to the nearest MRI "
            "scalp point); d_eq(eta) = the depth where that median falls below eta, linearly interpolated between bin "
            "centres (sensor white noise only, peak-channel SNR, magnetometers as the comparator, as in the sphere); eta "
            f"grid {etas[0]:g}-{etas[-1]:g} in steps of {step:g}. Points are drawn where a crossing exists: "
            + "; ".join(f"{style.ARRAY_LABEL[a]} {rng_txt[a]}" for a in arr)
            + f"; the OPM leads in every depth bin for eta {min(ahead['opm_dense']):g}-{max(ahead['opm_dense']):g} (both "
            f"arrays) and Neuromag already in the shallowest bin for eta >= {min(behind['opm_matched']):g} (matched) and "
            f">= {min(behind['opm_dense']):g} (dense). Dashed: the sphere of (a) with the real median standoffs (OPM "
            f"matched array {xi_opm:.2f} mm, Neuromag magnetometer coils {xi_sq:.2f} mm from the scalp). Depth axis "
            "inverted (deeper down) and shared by both panels."),
        alt=("Two line charts of the equal-SNR depth below the scalp (vertical axis inverted, deeper down) against the "
             f"noise ratio eta from {etas[0]:g} to {etas[-1]:g}. Left, the sphere: the depth falls from the sphere centre "
             f"at eta {eta0:.2f} to the brain surface at eta {eta1:.2f}; above the curve the OPM's SNR is higher; an open "
             f"diamond marks the preprint's printed {printed['printed_mm']:g} mm at eta 3, next to the exact "
             f"{exact[3.0]:.1f} mm. Right, the realistic adult head: the matched and dense OPM arrays cross between "
             f"{cross_lo:.0f} and {cross_hi:.0f} mm for eta {eta_lo:g} to {eta_hi:g}, close to the dashed sphere with the "
             "real standoffs."),
        caption_draft=(
            "The spherical benchmark and the realistic head. (a) Depth below the scalp at which an OPM on the scalp and a "
            f"SQUID {p['xi_squid_mm']:g} mm above it have equal single-channel SNR, against their noise ratio eta, in the "
            f"adult sphere of Jas et al. (2026) (head radius {h:g} mm, brain radius {b:g} mm, {p['Q_nAm']:g}-nAm tangential "
            "dipole; our reimplementation; dots: exact roots; shaded: OPM SNR higher). A crossing exists only for eta "
            f"between {eta0:.2f} and {eta1:.2f}; the diamond marks the {printed['printed_mm']:g} mm printed in the preprint "
            f"at eta = 3 (exact {exact[3.0]:.1f} mm). (b) The same quantity on the realistic adult head: depth at which the "
            "median over targets of the peak field of the OPM array over that of the best Neuromag magnetometer equals eta "
            f"(sensor noise only), for the matched ({n['opm_matched']} sites) and dense ({n['opm_dense']} sites) arrays, "
            f"and the sphere of (a) with the real median standoffs ({xi_opm:.1f} and {xi_sq:.1f} mm). Adapted from Jas et "
            "al. (2026), CC BY 4.0."),
        values=values)


# ----------------------------------------------------------------------------------------------
# R11 and R13: D on the inflated cortex
def cortex_values(subjects_dir: Path, subject: str, src_file: Path, rows: list[dict], values: np.ndarray) -> list:
    """Per hemisphere: (inflated vertices, triangles, gyral mask, values, category) at the vertices of
    the oct-6 source space. Category: DATA (a target), MEDIAL (a target on the medial wall, region
    '*unknown': not cortex, kept in the G2 statistics, left out of the G3B ones) or EXCLUDED (a source
    space vertex that is not a target: closer than 4 mm to the inner-skull mesh or outside it)."""
    src = mne.read_source_spaces(src_file, verbose=False)
    views = plotting.inflated_views(subjects_dir, subject, src)
    hemi = np.array([int(r["hemi"]) for r in rows])
    vno = np.array([int(r["vertno"]) for r in rows])
    medial = np.char.endswith(np.array([r["region"] for r in rows]).astype(str), "unknown")
    out = []
    for h, (verts, tris, gyral) in enumerate(views):
        vert = np.asarray(src[h]["vertno"])
        m = hemi == h
        idx = np.searchsorted(vert, vno[m])
        if not np.array_equal(vert[np.minimum(idx, len(vert) - 1)], vno[m]):
            raise SystemExit(f"{subject}: stored targets are not vertices of {src_file.name}")
        val = np.full(len(vert), np.nan)
        cat = np.full(len(vert), EXCLUDED)
        cat[idx] = np.where(medial[m], MEDIAL, DATA)
        val[idx] = np.where(medial[m], np.nan, values[m])
        out.append((verts, tris, gyral, val, cat))
    return out


def render(ax, verts, tris, gyral, val, cat, camera_x: float, norm) -> None:
    """One orthographic view (opmsquid.plotting.render_view, with two grey categories): back faces
    culled, triangles painted back to front; a triangle takes the category that at least two of its
    vertices share and, if that is DATA, the mean of its targets' values."""
    d = np.array([camera_x, 0.0, 0.0])
    screen = np.column_stack([verts[:, 1] * camera_x, verts[:, 2]])
    tri_n = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]], verts[tris[:, 2]] - verts[tris[:, 0]])
    t = tris[(tri_n @ d) / np.linalg.norm(tri_n, axis=1) > 0]
    t = t[np.argsort(verts[t].mean(axis=1) @ d)]
    c = plotting.triangle_mode(cat[t])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN triangles (grey ones) give NaN
        v = np.nanmean(val[t], axis=1)
    rgba = CMAP(norm(np.nan_to_num(v)))
    for k, grey in GREY.items():
        rgba[c == k, :3] = grey
    lit = c == DATA
    rgba[lit, :3] *= (1.0 - SHADE + SHADE * gyral[t[lit]].mean(axis=1))[:, None]
    rgba[~lit, :3] *= (1.0 - SHADE / 2 + SHADE / 2 * gyral[t[~lit]].mean(axis=1))[:, None]
    ax.add_collection(PolyCollection(screen[t], facecolors=rgba, edgecolors="none", antialiased=False, rasterized=True))
    ax.set_xlim(screen[:, 0].min() - 1.5, screen[:, 0].max() + 1.5)
    ax.set_ylim(screen[:, 1].min() - 1.5, screen[:, 1].max() + 1.5)
    ax.set_aspect("equal")
    ax.axis("off")


def map_figure(name: str, rows: list[tuple[str, list]], note: str, extend: str) -> None:
    """Rows of four views (left lateral, left medial, right medial, right lateral) with a title above
    each row, one horizontal colour bar and the two greys below."""
    norm = Normalize(-LIM_DB, LIM_DB)
    asp = max(np.ptp(v[:, 2]) / np.ptp(v[:, 1]) for _, hem in rows for v, *_ in hem)
    W, lm, rm, gap = style.FULL_W, 0.04, 0.04, 0.06
    cw = (W - lm - rm - 3 * gap) / 4
    rh, th, top = cw * asp + 0.04, 0.27, 0.22
    H = top + len(rows) * (th + rh) + 0.62
    fig = plt.figure(figsize=(W, H))
    for i, (title, hem) in enumerate(rows):
        y0 = H - top - (i + 1) * (th + rh)
        fig.text(lm / W, (y0 + rh + th - 0.04) / H, title, ha="left", va="top", fontsize=9)
        for j, (h, cam, label) in enumerate(VIEWS):
            ax = fig.add_axes([(lm + j * (cw + gap)) / W, y0 / H, cw / W, rh / H])
            render(ax, *hem[h], cam, norm)
            if i == 0:
                fig.text((lm + j * (cw + gap) + cw / 2) / W, (H - 0.04) / H, label, ha="center", va="top",
                         fontsize=8, color="0.3")
    cw_bar = 2.8
    cax = fig.add_axes([lm / W, 0.30 / H, cw_bar / W, 0.11 / H])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal", extend=extend)
    ticks = np.arange(-LIM_DB, LIM_DB + 0.1, 2.0)
    cb.set_ticks(ticks)
    cb.ax.set_xticklabels([f"{x:+.0f}".replace("-", "\u2212") if x else "0" for x in ticks])
    cb.ax.tick_params(labelsize=7.5, length=2)
    cb.outline.set_linewidth(0.5)
    fig.text(lm / W, 0.47 / H, "Neuromag higher", ha="left", va="bottom", fontsize=7.5, color="#2166AC")
    fig.text((lm + cw_bar) / W, 0.47 / H, "OPM higher", ha="right", va="bottom", fontsize=7.5, color="#B2182B")
    fig.text((lm + cw_bar / 2) / W, 0.47 / H, "D (dB)", ha="center", va="bottom", fontsize=8)
    greys = [Patch(facecolor=str(GREY[MEDIAL]), edgecolor="none", label="medial wall (not cortex)"),
             Patch(facecolor=str(GREY[EXCLUDED]), edgecolor="none",
                   label=f"within {NEAR_SKULL_MM:g} mm of the inner skull (not simulated)")]
    x_leg = lm + cw_bar + 0.36
    fig.legend(handles=greys, loc="lower left", bbox_to_anchor=(x_leg / W, 0.24 / H), ncol=1, handlelength=1.2,
               handleheight=0.9, fontsize=7.5, borderaxespad=0.0)
    fig.text(x_leg / W, 0.04 / H, note, ha="left", va="bottom", fontsize=7.5, color="0.3")
    style.save(fig, name)


def extend_of(arrays: list[np.ndarray]) -> str:
    lo = any(np.nanmin(a) < -LIM_DB for a in arrays)
    hi = any(np.nanmax(a) > LIM_DB for a in arrays)
    return {(False, False): "neither", (True, False): "min", (False, True): "max", (True, True): "both"}[(lo, hi)]


def adult_maps() -> tuple[list, dict]:
    """R11 data: D per G2 target for the dense and matched arrays against Neuromag's 306 channels."""
    s = load(G2)
    rows = io.read_csv(ROOT / G2_TARGETS)
    sq = col(rows, f"detect_squid_combined_{COND}")
    medial = np.char.endswith(np.array([r["region"] for r in rows]).astype(str), "unknown")
    src_file = paths.SUBJECTS_DIR / "sample" / "bem" / "sample-oct-6-src.fif"
    out, vals = [], {}
    for letter, a in zip("ab", ("opm_dense", "opm_matched")):
        d = 20.0 * np.log10(col(rows, f"detect_{a}_opm_{COND}") / sq)
        stored = s["primary"]["oracle"][f"{a}/combined/{COND}"]["median_log2"] * DB_PER_LOG2
        if abs(np.median(d) - stored) > 0.005:
            raise SystemExit(f"{G2_TARGETS} does not reproduce the stored median for {a} "
                             f"({np.median(d):.4f} vs {stored:.4f} dB)")
        hem = cortex_values(paths.SUBJECTS_DIR, "sample", src_file, rows, d)
        out.append((f"({letter}) {style.ARRAY_LABEL[a]} ({s['arrays'][a]['n_sites']} sites) vs Neuromag "
                    f"({s['arrays']['squid']['channels']} channels)", hem, d))
        vals[a] = dict(median_all_targets_dB=float(np.median(d)), stored_median_dB=float(stored),
                       median_without_medial_wall_dB=float(np.median(d[~medial])), share_positive=float(np.mean(d > 0)),
                       min_dB=float(d.min()), max_dB=float(d.max()), share_above_limit=float(np.mean(d > LIM_DB)),
                       share_below_minus_limit=float(np.mean(d < -LIM_DB)))
    vals["targets"] = dict(n=len(rows), medial_wall=int(medial.sum()),
                           source_space_vertices=sum(len(hm[4]) for hm in hem),
                           not_targets=sum(int(np.sum(hm[4] == EXCLUDED)) for hm in hem))
    return out, vals


def template_maps() -> tuple[list, dict]:
    """R13 data: D per G3B target (dense array, Neuromag 306, top contact) for the 24- and 12-month templates."""
    g = load(G3B)
    sd = paths.EXTERNAL / anatomy.INFANT_SUBJECTS
    out, vals = [], {}
    for letter, (k, name) in zip("ab", TEMPLATES.items()):
        rows = io.read_csv(ROOT / G3B_TARGETS.format(k))
        d = 20.0 * np.log10(col(rows, f"detect_opm_dense_opm_{COND}") / col(rows, f"detect_squid_top_combined_{COND}"))
        medial = np.char.endswith(np.array([r["region"] for r in rows]).astype(str), "unknown")
        area = col(rows, "area_mm2")
        wmed = P.weighted_median(d[~medial], area[~medial])
        stored = g["placement_D"][f"{k}/top/combined/{COND}"]["median"]
        if abs(wmed - stored) > 0.005:
            raise SystemExit(f"{G3B_TARGETS.format(k)} does not reproduce the stored D ({wmed:.4f} vs {stored:.4f} dB)")
        hem = cortex_values(sd, name, sd / name / "bem" / f"{name}-oct-6-src.fif", rows, d)
        n_sites, n_sq = g["arrays"][k]["opm_dense"]["n_sites"], g["arrays"][k]["squid:top"]["n"]
        out.append((f"({letter}) {style.ANAT_LABEL[k]}: {style.ARRAY_LABEL['opm_dense']} ({n_sites} sites) vs "
                    f"Neuromag ({n_sq} channels)", hem, d))
        vals[k] = dict(template=name, n_sites=n_sites, n_squid_channels=n_sq,
                       area_weighted_median_without_medial_wall_dB=wmed, stored_D_dB=stored,
                       share_positive_area=float(np.sum(area[~medial][d[~medial] > 0]) / np.sum(area[~medial])),
                       min_dB=float(d.min()), max_dB=float(d.max()), share_above_limit=float(np.mean(d > LIM_DB)),
                       n_targets=len(rows), medial_wall=int(medial.sum()),
                       source_space_vertices=sum(len(hm[4]) for hm in hem),
                       not_targets=sum(int(np.sum(hm[4] == EXCLUDED)) for hm in hem))
    return out, vals


def scale_factors(g: dict) -> dict:
    """Linear scale factors of the adult and the two scaled adults as scripts/g3b_pediatric_helmet.py sets them
    (config.anatomy.school_age_scale; the 24-month template's head circumference over the adult's), checked against
    the scale notes stored with the anatomies."""
    a = g["anatomies"]
    f = {"adult": 1.0, "school": float(g["config"]["anatomy"]["school_age_scale"]),
         "size2yr": a["infant2yr"]["head_size"]["ofc_mm"] / a["adult"]["head_size"]["ofc_mm"]}
    for k in ("school", "size2yr"):
        if not a[k]["scale_note"].endswith(f": {f[k]:.4f}"):
            raise SystemExit(f"{G3B}: the scale note of {k} ('{a[k]['scale_note']}') does not give {f[k]:.4f}")
    return f


def scaled_maps() -> tuple[list, dict]:
    """R14 data: D per G3B target (dense array, Neuromag 306, top contact) for the adult and the two scaled adults,
    all on the adult's inflated cortex: a scaled adult's targets are adult targets (the same vertices, those still at
    least 4 mm inside its scaled inner skull)."""
    g = load(G3B)
    src_file = paths.SUBJECTS_DIR / "sample" / "bem" / "sample-oct-6-src.fif"
    adult = {(r["hemi"], r["vertno"]) for r in io.read_csv(ROOT / G3B_TARGETS.format("adult"))}
    f = scale_factors(g)
    out, vals = [], {}
    for letter, k in zip("abc", SCALED_HEADS):
        rows = io.read_csv(ROOT / G3B_TARGETS.format(k))
        own = {(r["hemi"], r["vertno"]) for r in rows}
        if not own <= adult:
            raise SystemExit(f"{G3B_TARGETS.format(k)}: targets that are not targets of the adult")
        d = 20.0 * np.log10(col(rows, f"detect_opm_dense_opm_{COND}") / col(rows, f"detect_squid_top_combined_{COND}"))
        medial = np.char.endswith(np.array([r["region"] for r in rows]).astype(str), "unknown")
        area = col(rows, "area_mm2")
        wmed = P.weighted_median(d[~medial], area[~medial])
        stored = g["placement_D"][f"{k}/top/combined/{COND}"]["median"]
        if abs(wmed - stored) > 0.005:
            raise SystemExit(f"{G3B_TARGETS.format(k)} does not reproduce the stored D ({wmed:.4f} vs {stored:.4f} dB)")
        hem = cortex_values(paths.SUBJECTS_DIR, "sample", src_file, rows, d)
        n_sites, n_sq = g["arrays"][k]["opm_dense"]["n_sites"], g["arrays"][k]["squid:top"]["n"]
        out.append((f"({letter}) {style.ANAT_LABEL[k]}: {style.ARRAY_LABEL['opm_dense']} ({n_sites} sites) vs "
                    f"Neuromag ({n_sq} channels)", hem, d))
        vals[k] = dict(scale=f[k], n_sites=n_sites, n_squid_channels=n_sq,
                       area_weighted_median_without_medial_wall_dB=wmed, stored_D_dB=stored,
                       share_positive_area=float(np.sum(area[~medial][d[~medial] > 0]) / np.sum(area[~medial])),
                       min_dB=float(d.min()), max_dB=float(d.max()), share_above_limit=float(np.mean(d > LIM_DB)),
                       share_below_minus_limit=float(np.mean(d < -LIM_DB)), n_targets=len(rows),
                       adult_targets_not_targets_here=len(adult - own), medial_wall=int(medial.sum()),
                       source_space_vertices=sum(len(hm[4]) for hm in hem),
                       not_targets=sum(int(np.sum(hm[4] == EXCLUDED)) for hm in hem))
    return out, vals


def figures_maps() -> dict:
    adult, va = adult_maps()
    heads, vh = template_maps()
    scaled, vs = scaled_maps()
    ext = extend_of([d for *_, d in adult] + [d for *_, d in heads] + [d for *_, d in scaled])
    g2s, g3 = load(G2), load(G3B)
    q = g2s["config"]["sources"]["focal_nAm"]
    map_figure("Figure_R11_maps_adult", [(t, h) for t, h, _ in adult],
               f"Adult at its measured head position; sensor + brain noise; {q:g}-nAm dipoles", ext)
    map_figure("Figure_R13_maps_heads", [(t, h) for t, h, _ in heads],
               f"Top contact in the adult helmet; sensor + brain noise; {q:g}-nAm dipoles", ext)
    map_figure("Figure_R14_maps_scaled", [(t, h) for t, h, _ in scaled],
               f"Top contact in the adult helmet; sensor + brain noise; {q:g}-nAm dipoles", ext)
    n_sq = g2s["arrays"]["squid"]["channels"]
    vd, vm = va["opm_dense"], va["opm_matched"]
    v2, v1 = vh["infant2yr"], vh["infant12mo"]
    # the adult under the pediatric convention (area-weighted, medial wall excluded): the G3B run, whose adult at its
    # measured pose must agree with G2 target by target
    adult_aw = g3["placement_D"][f"adult/centred/combined/{COND}"]["median"]
    adult_top = g3["placement_D"][f"adult/top/combined/{COND}"]["median"]
    rows2, rows3 = io.read_csv(ROOT / G2_TARGETS), io.read_csv(ROOT / G3B_TARGETS.format("adult"))
    if [(r["hemi"], r["vertno"]) for r in rows2] != [(r["hemi"], r["vertno"]) for r in rows3]:
        raise SystemExit("the G2 and G3B adult targets differ")
    d2 = 20.0 * np.log10(col(rows2, f"detect_opm_dense_opm_{COND}") / col(rows2, f"detect_squid_combined_{COND}"))
    d3 = 20.0 * np.log10(col(rows3, f"detect_opm_dense_opm_{COND}") / col(rows3, f"detect_squid_centred_combined_{COND}"))
    agree = float(np.max(np.abs(d2 - d3)))
    if agree > 0.01:
        raise SystemExit(f"the G3B adult at its measured pose differs from G2 by up to {agree:.4f} dB per target")
    vd.update(g3b_area_weighted_without_medial_wall_measured_pose_dB=adult_aw,
              g3b_measured_pose_max_abs_diff_per_target_dB=agree)
    vh["adult_top_contact_D_dB"] = adult_top
    clearance = g3["config"]["placement"]["clearance_mm"]
    grey_txt = ("Light grey: medial wall (targets in FreeSurfer's 'unknown' region: the cut through the corpus callosum "
                "and midbrain, not cortex); dark grey: source-space vertices that are not targets, i.e. within "
                f"{NEAR_SKULL_MM:g} mm of the 5,120-triangle inner-skull mesh or outside it (A-BEM-DIST; the shallowest "
                "cortex, which the noise model also leaves out of the brain background). A triangle takes the category "
                "of at least two of its vertices.")
    scale_txt = (f"Diverging colour scale (red: OPM higher), -{LIM_DB:g} to +{LIM_DB:g} dB, shared by R11, R13 and R14 "
                 f"(display choice); values beyond it take the end colour (colour bar extension: '{ext}'). Views: "
                 "orthographic, left lateral, left medial, right medial, right lateral, with sulcal shading "
                 "(opmsquid.plotting); each hemisphere fills its panel (not to scale).")
    r11 = dict(
        inputs=[f"{G2_TARGETS} :: hemi, vertno, region, detect_opm_dense_opm_{COND}, detect_opm_matched_opm_{COND}, "
                f"detect_squid_combined_{COND}",
                f"{G2} :: primary.oracle['<array>/combined/{COND}'].median_log2 (check), arrays.<array>.n_sites, "
                "arrays.squid.channels, config.sources.focal_nAm",
                f"{G3B} :: placement_D['adult/centred/combined/{COND}'].median (the adult under the pediatric "
                "convention)",
                f"{G3B_TARGETS.format('adult')} :: detect_opm_dense_opm_{COND}, detect_squid_centred_combined_{COND} "
                "(check against G2)",
                "MNE-sample-data/subjects/sample: bem/sample-oct-6-src.fif, surf/?h.inflated, surf/?h.sulc (drawing only)"],
        description=(
            "D = 20 log10(detect_<array>_opm / detect_squid_combined) per cortical target of the adult G2 run (MNE "
            f"sample subject at its measured head position, {va['targets']['n']:,} targets, {q:g}-nAm cortical-normal "
            "dipoles, sensor white noise + cortical background, oracle covariance), drawn at the target's vertex of the "
            "oct-6 source space on the inflated white surface. The median over all targets is "
            f"{vd['median_all_targets_dB']:.3f} dB (dense) and {vm['median_all_targets_dB']:.3f} dB (matched), equal to "
            f"the stored median log2 ratios times 20 log10 2 ({vd['stored_median_dB']:.3f} and "
            f"{vm['stored_median_dB']:.3f} dB; checked). Under the pediatric convention (area-weighted median, medial "
            f"wall excluded) the dense array's D at the same pose is {adult_aw:.3f} dB ({G3B} "
            f"placement_D['adult/centred/combined/{COND}']; that run's per-target D at this pose agrees with G2's within "
            f"{agree:.4f} dB). " + grey_txt + " " + scale_txt + f" Medial wall: {va['targets']['medial_wall']} targets; "
            f"not targets: {va['targets']['not_targets']} of {va['targets']['source_space_vertices']:,} source-space "
            "vertices."),
        alt=("Inflated left and right cortical hemispheres in lateral and medial views, two rows. Top, the dense OPM "
             f"array against Neuromag's {n_sq} channels: red almost everywhere ({vd['share_positive']:.1%} of targets "
             "above 0 dB), strongest on lateral convexities, paler on medial and ventral surfaces. Bottom, the matched "
             f"array: a mix of pale red and blue ({vm['share_positive']:.0%} above 0 dB), blue on the ventral and polar "
             "temporal and frontal surfaces. Grey: medial wall (light) and scattered vertices near the inner skull "
             "(dark)."),
        caption_draft=(
            "D on the adult's inflated cortex (MNE sample subject at its measured head position; sensor plus brain "
            f"noise; {q:g}-nAm dipoles). (a) Dense array ({g2s['arrays']['opm_dense']['n_sites']} sites) and (b) matched "
            f"array ({g2s['arrays']['opm_matched']['n_sites']} sites) against Neuromag's {n_sq} channels; red: the "
            f"OPM's detectability is higher. Median D {fmt_db(vd['median_all_targets_dB'])} dB (dense) and "
            f"{fmt_db(vm['median_all_targets_dB'])} dB (matched) over all {va['targets']['n']:,} targets (dense array, "
            f"area-weighted without the medial wall: {fmt_db(adult_aw)} dB). Light grey: medial wall; dark grey: cortex "
            f"within {NEAR_SKULL_MM:g} mm of the inner skull, which is not simulated as a target and is also left out of "
            "the brain background. The colour scale is that of Figure R13."),
        values=dict(va, colour_limit_dB=LIM_DB, extend=ext))
    w2 = v2["area_weighted_median_without_medial_wall_dB"]
    w1 = v1["area_weighted_median_without_medial_wall_dB"]
    r13 = dict(
        inputs=[f"{G3B_TARGETS.format(k)} :: hemi, vertno, region, area_mm2, detect_opm_dense_opm_{COND}, "
                f"detect_squid_top_combined_{COND}" for k in TEMPLATES]
               + [f"{G3B} :: placement_D['<template>/top/combined/{COND}'].median (check), "
                  f"placement_D['adult/top/combined/{COND}'].median (comparison), arrays.<template>.opm_dense.n_sites, "
                  "arrays.<template>['squid:top'].n, config.placement.clearance_mm"]
               + [f"infant_subjects/{n}: bem/{n}-oct-6-src.fif, surf/?h.inflated, surf/?h.sulc (drawing only)"
                  for n in TEMPLATES.values()],
        description=(
            "D = 20 log10(detect_opm_dense_opm / detect_squid_top_combined) per cortical target of the G3B run for the "
            "24-month (ANTS2-0Years3T) and 12-month (ANTS12-0Months3T) infant templates in their native dimensions, the "
            f"dense OPM array refitted to each head, Neuromag's {v2['n_squid_channels']} channels with the head at top "
            "contact in the fixed adult helmet (the adult's measured pose, then raised until the nearest magnetometer "
            f"coil is {clearance:g} mm from the scalp), sensor + brain noise, {q:g}-nAm cortical-normal dipoles, drawn "
            "on each template's own inflated cortex (no vertex correspondence with the adult). Area-weighted median D "
            f"without the medial wall (opmsquid.pediatric.weighted_median, weights area_mm2): {w2:.3f} dB (24 months) "
            f"and {w1:.3f} dB (12 months), equal to the stored values ({v2['stored_D_dB']:.3f} and "
            f"{v1['stored_D_dB']:.3f} dB; checked); the adult at top contact: {adult_top:.3f} dB "
            f"(placement_D['adult/top/combined/{COND}']). " + grey_txt + " " + scale_txt
            + f" Targets above +{LIM_DB:g} dB: {v2['share_above_limit']:.1%} (24 months), "
            f"{v1['share_above_limit']:.1%} (12 months). Infant templates: O'Reilly et al. (2021), built from the "
            "Neurodevelopmental MRI Database (Richards et al., 2016)."),
        alt=("Inflated cortical hemispheres of two infant templates in lateral and medial views, one row each (24 and "
             f"12 months). Both are red almost everywhere (D positive on {v2['share_positive_area']:.0%} and "
             f"{v1['share_positive_area']:.0%} of the cortical area), darkest on the lateral convexities, where "
             f"{v2['share_above_limit']:.0%} and {v1['share_above_limit']:.0%} of the targets lie beyond the end of the "
             "colour scale, and paler on the medial surfaces."),
        caption_draft=(
            "D on the inflated cortex of the (a) 24-month and (b) 12-month infant templates (O'Reilly et al., 2021; "
            f"Richards et al., 2016): dense OPM array ({v2['n_sites']} and {v1['n_sites']} sites) against Neuromag's "
            f"{v2['n_squid_channels']} channels with the head at top contact in the adult helmet; sensor plus brain "
            f"noise; {q:g}-nAm dipoles. Area-weighted median D {fmt_db(w2)} and {fmt_db(w1)} dB (medial wall excluded; "
            f"the adult at top contact {fmt_db(adult_top)} dB). Colour scale and greys as in Figure R11."),
        values=dict(vh, colour_limit_dB=LIM_DB, extend=ext))
    sa, ss, s2 = vs["adult"], vs["school"], vs["size2yr"]
    wa, ws, w2s = (v["area_weighted_median_without_medial_wall_dB"] for v in (sa, ss, s2))
    if abs(wa - adult_top) > 0.005:
        raise SystemExit(f"{G3B_TARGETS.format('adult')}: top-contact D {wa:.4f} dB, stored {adult_top:.4f} dB")
    r14 = dict(
        inputs=[f"{G3B_TARGETS.format(k)} :: hemi, vertno, region, area_mm2, detect_opm_dense_opm_{COND}, "
                f"detect_squid_top_combined_{COND}" for k in SCALED_HEADS]
               + [f"{G3B} :: placement_D['<head>/top/combined/{COND}'].median (check), arrays.<head>.opm_dense.n_sites, "
                  "arrays.<head>['squid:top'].n, config.anatomy.school_age_scale, anatomies.<adult|infant2yr>.head_size.ofc_mm "
                  "(the 2-year size factor), anatomies.<school|size2yr>.scale_note (check), config.placement.clearance_mm",
                  "MNE-sample-data/subjects/sample: bem/sample-oct-6-src.fif, surf/?h.inflated, surf/?h.sulc (drawing only; "
                  "all three rows on the adult's cortex)"],
        description=(
            "D = 20 log10(detect_opm_dense_opm / detect_squid_top_combined) per cortical target of the G3B run for the adult "
            "(MNE sample subject) and its two size-only controls, the adult scaled about its MRI origin by "
            f"{ss['scale']:.4f} (school-age size) and {s2['scale']:.4f} (2-year size; the scale notes stored with the run: "
            f"'{g3['anatomies']['school']['scale_note']}' and '{g3['anatomies']['size2yr']['scale_note']}'): the dense OPM "
            f"array refitted to each head ({sa['n_sites']}, {ss['n_sites']} and {s2['n_sites']} sites), Neuromag's "
            f"{sa['n_squid_channels']} channels with the head at top contact in the fixed adult helmet (the adult's measured "
            f"pose, then raised until the nearest magnetometer coil is {clearance:g} mm from the scalp), sensor + brain "
            f"noise, {q:g}-nAm cortical-normal dipoles. A scaled adult keeps the adult's vertices, so all three rows are "
            "drawn on the adult's inflated cortex and compare vertex by vertex; its targets are the adult's that stay "
            f"usable on the scaled meshes (inside the scaled inner skull and at least {NEAR_SKULL_MM:g} mm from it), so "
            f"{ss['adult_targets_not_targets_here']} and {s2['adult_targets_not_targets_here']} of the adult's "
            f"{sa['n_targets']:,} targets are not targets there and are drawn dark grey. Area-weighted median D without the "
            f"medial wall (opmsquid.pediatric.weighted_median, weights area_mm2): {wa:.3f} dB (adult), {ws:.3f} dB "
            f"(school-age size) and {w2s:.3f} dB (2-year size), equal to the stored values ({sa['stored_D_dB']:.3f}, "
            f"{ss['stored_D_dB']:.3f} and {s2['stored_D_dB']:.3f} dB; checked). The adult at its measured position is "
            "Figure R11a. " + grey_txt + " " + scale_txt + f" Targets above +{LIM_DB:g} dB: {sa['share_above_limit']:.1%} "
            f"(adult), {ss['share_above_limit']:.1%} (school-age size), {s2['share_above_limit']:.1%} (2-year size); below "
            f"-{LIM_DB:g} dB: " + (", ".join(f"{v['share_below_minus_limit']:.1%}" for v in (sa, ss, s2))
                                   if any(v["share_below_minus_limit"] for v in (sa, ss, s2)) else "none") + "."),
        alt=("Inflated left and right cortical hemispheres of the adult in lateral and medial views, three rows: the adult "
             "and the adult scaled to school-age size and to 2-year size, all at top contact. All three are red over most of "
             f"the cortex (D positive on {sa['share_positive_area']:.1%}, {ss['share_positive_area']:.1%} and "
             f"{s2['share_positive_area']:.1%} of the cortical area) and redder with each smaller size, most on the lateral "
             "convexities, while the medial surfaces stay pale. Grey: medial wall (light) and vertices near the inner skull "
             "(dark), more of them in the scaled heads."),
        caption_draft=(
            "D on the adult's inflated cortex for (a) the adult and the adult scaled (b) to school-age size "
            f"(\u00d7{ss['scale']:.3f}) and (c) to the 24-month template's head circumference (\u00d7{s2['scale']:.3f}), each "
            f"at top contact in the adult helmet: dense OPM array ({sa['n_sites']}, {ss['n_sites']} and {s2['n_sites']} "
            f"sites) against Neuromag's {sa['n_squid_channels']} channels; sensor plus brain noise; {q:g}-nAm dipoles. The "
            "scaled adults keep the adult's vertices, so the rows compare vertex by vertex. Area-weighted median D "
            f"{fmt_db(wa)}, {fmt_db(ws)} and {fmt_db(w2s)} dB (medial wall excluded). Dark grey also marks the adult's "
            f"targets that come within {NEAR_SKULL_MM:g} mm of the scaled inner skull. Colour scale and greys as in "
            "Figure R11."),
        values=dict(vs, colour_limit_dB=LIM_DB, extend=ext))
    return {"Figure_R11_maps_adult": r11, "Figure_R13_maps_heads": r13, "Figure_R14_maps_scaled": r14}


# ----------------------------------------------------------------------------------------------
# R12: geometry
def segments(polylines: list) -> np.ndarray:
    """Segments (n, 2, 2) between consecutive points of each polyline of the geometry export."""
    return np.concatenate([np.stack([p[:-1], p[1:]], axis=1) for p in map(np.asarray, polylines)])


def geometry_sections(g: dict) -> tuple[dict, dict]:
    """Sections, magnetometer coil centres and dense OPM sites of the three heads, head frame [mm]
    (results/g3b/g3b_geometry_sections.json, exported by scripts/export_g3b_geometry.py from the stored
    state of the G3B run), checked against results/g3b/g3b_summary.json: the export must be of the same
    run and give its moves, scale factors, coil-to-scalp gaps (recomputed by the export from the
    unrounded geometry) and OPM site counts, and the drawn coil centres must be the exported transforms
    applied to the device-frame coil centres (to the rounding)."""
    x = load(GEOMETRY)
    run = x["provenance"]["source_state"]["commit"]
    if run != g["provenance"]["commit"]:
        raise SystemExit(f"{GEOMETRY}: exported from the state of run {run}, {G3B} is run {g['provenance']['commit']}")
    dev = np.asarray(x["neuromag"]["coil_centres_device_mm"], float)
    out = {}
    for k, name in GEOMETRY_HEADS.items():
        a, stored = x["anatomies"][k], g["placements"][k]
        if name is not None and name not in a["description"]:  # the description below names the subjects
            raise SystemExit(f"{GEOMETRY}: {k} is '{a['description']}', not {name}")
        tr, cf = a["transforms"], a["transforms"]["counterfactual_x-centred"]
        checks = [("top moved_mm", tr["top"]["moved_mm"], stored["top"]["moved_mm"], 1e-9),
                  ("k", cf["k"], stored["counterfactual_x-centred"]["k"], 1e-9),
                  ("k_nominal", cf["k_nominal"], stored["counterfactual_x-centred"]["k_nominal"], 1e-9)]
        checks += [(f"{name} {key}", a["gaps_mm"][name][key], stored[name][key], 1e-6)
                   for name in ("top", "counterfactual_x-centred") for key in ("min_dist_mm", "median_dist_mm", "max_dist_mm")]
        for what, have, want, tol in checks:
            if abs(have - want) > tol:
                raise SystemExit(f"{GEOMETRY}: {k} {what} {have} differs from {G3B} ({want})")
        top = np.asarray(tr["top"]["device_to_head"], float)
        cfx = np.asarray(cf["device_to_head"], float)
        centre = np.linalg.inv(cfx)[:3, 3] * 1e3  # head origin, device frame [mm]
        expect = {"fixed_top": dev @ top[:3, :3].T + top[:3, 3] * 1e3,
                  "scaled_x_centred": (centre + cf["k"] * (dev - centre)) @ cfx[:3, :3].T + cfx[:3, 3] * 1e3}
        coils = {n: np.asarray(a["coils_mm"][n], float) for n in expect}
        coil_err = max(float(np.abs(coils[n] - expect[n]).max()) for n in expect)
        if coil_err > COIL_TOL_MM:
            raise SystemExit(f"{GEOMETRY}: {k} coil centres differ from the exported transforms by {coil_err:.4f} mm")
        opm = np.asarray(a["opm_dense_mm"], float)
        if len(opm) != g["arrays"][k]["opm_dense"]["n_sites"]:
            raise SystemExit(f"{GEOMETRY}: {k} has {len(opm)} dense OPM sites, {G3B} {g['arrays'][k]['opm_dense']['n_sites']}")
        sec = a["sections"]
        out[k] = dict(scalp={plane: segments(sec[plane]["scalp"]) for plane in ("sagittal", "coronal")},
                      white=[segments(sec["coronal"]["white"][h]) for h in ("lh", "rh")],
                      fixed=coils["fixed_top"], scaled=coils["scaled_x_centred"], opm=opm, k=float(cf["k"]),
                      coil_err=coil_err)
    return out, x


def figure_geometry() -> dict:
    g = load(G3B)
    heads, geo = geometry_sections(g)
    n_coils = int(geo["neuromag"]["magnetometers"])

    # panels of equal size and scale: the sagittal and coronal sections share the horizontal span
    xs, zl = SECTION_LIM["sagittal"]
    span_x, span_z = xs[1] - xs[0], zl[1] - zl[0]
    W, lm, rm, gap = style.FULL_W, 0.66, 0.04, 0.10
    pw = (W - lm - rm - 2 * gap) / 3
    ph = pw * span_z / span_x
    t_top, t_mid, t_bot = 0.62, 0.36, 0.92  # title block, sagittal x label, coronal x label + legend
    H = t_top + 2 * ph + t_mid + t_bot
    fig = plt.figure(figsize=(W, H))
    markers = ((dict(marker="o", ms=3.0, mfc="none", mec=style.ARRAY_COLOR["opm_dense"], mew=0.8), "opm"),
               (dict(marker="s", ms=3.6, mfc="none", mec=SCALED_COLOR, mew=0.9), "scaled"),
               (dict(marker="s", ms=3.0, mfc="k", mec="k", mew=0.0), "fixed"))
    vals = {}
    for j, k in enumerate(GEOMETRY_HEADS):
        hd, pl = heads[k], g["placements"][k]
        cf = pl["counterfactual_x-centred"]
        x0 = lm + j * (pw + gap)
        for i, (plane, axis) in enumerate((("sagittal", 0), ("coronal", 1))):
            ax = fig.add_axes([x0 / W, (H - t_top - ph if i == 0 else t_bot) / H, pw / W, ph / H])
            other = 1 - axis  # horizontal coordinate: y (anterior) in the sagittal, x (right) in the coronal section
            # section points are in-plane, (y, z) or (x, z); no pixel snapping: a segment made exactly vertical or
            # horizontal by the 0.01-mm rounding would otherwise be snapped, which the unrounded sections never were
            if plane == "coronal":  # the midsagittal plane passes between the hemispheres
                for seg in hd["white"]:
                    ax.add_collection(LineCollection(seg, colors=CORTEX_GREY, linewidths=0.4, zorder=1, snap=False))
            ax.add_collection(LineCollection(hd["scalp"][plane], colors="k", linewidths=0.9, zorder=2, snap=False))
            for kw, key in markers:
                pts = hd[key]
                off = np.abs(pts[:, axis])
                if np.any(np.abs(off - SLAB_MM) <= ROUND_MM / 2):  # rounding could move a sensor across the slab edge
                    raise SystemExit(f"{GEOMETRY}: {k} has a {key} position within the rounding of the {SLAB_MM:g}-mm slab edge")
                sel = off < SLAB_MM
                ax.plot(pts[sel, other], pts[sel, 2], ls="none", zorder=3, **kw)
            ax.set_xlim(*SECTION_LIM[plane][0])
            ax.set_ylim(*zl)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([-100, 0, 100])
            ax.set_yticks([-50, 0, 50, 100, 150])
            ax.tick_params(labelsize=7, length=2, pad=1.5)
            if j:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(f"{plane.capitalize()} section\nsuperior (mm)", fontsize=8)
            ax.set_xlabel("anterior (mm)" if plane == "sagittal" else "right (mm)", fontsize=7.5, labelpad=1)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        cx = (x0 + pw / 2) / W
        fig.text(cx, (H - 0.04) / H, style.ANAT_LABEL[k], ha="center", va="top", fontsize=9)
        fig.text(cx, (H - 0.24) / H,
                 f"median gap: fixed {pl['top']['median_dist_mm']:.1f}, scaled {cf['median_dist_mm']:.1f} mm\n"
                 f"helmet scale \u00d7{hd['k']:.3f}; {g['arrays'][k]['opm_dense']['n_sites']} OPM sites",
                 ha="center", va="top", fontsize=7, color="0.25", linespacing=1.3)
        vals[k] = dict(top=dict(moved_mm=pl["top"]["moved_mm"], median_gap_mm=pl["top"]["median_dist_mm"],
                                min_gap_mm=pl["top"]["min_dist_mm"], max_gap_mm=pl["top"]["max_dist_mm"]),
                       scaled_laterally_centred=dict(k=hd["k"], k_nominal=cf["k_nominal"],
                                                     median_gap_mm=cf["median_dist_mm"],
                                                     min_gap_mm=cf["min_dist_mm"], max_gap_mm=cf["max_dist_mm"],
                                                     shift_x_mm=pl["x-centred"]["shift_x_mm"]),
                       opm_dense_sites=g["arrays"][k]["opm_dense"]["n_sites"],
                       head_circumference_cm=g["anatomies"][k]["head_size"]["ofc_mm"] / 10.0)
    labels = ("dense OPM sites", "helmet scaled with the head, laterally centred",
              "fixed adult helmet (Neuromag), head at top contact")
    handles = [Line2D([], [], ls="none", label=lab, **kw) for (kw, _), lab in zip(markers, labels)][::-1]
    handles += [Line2D([], [], color="k", lw=0.9, label="scalp (MRI surface)"),
                Line2D([], [], color=CORTEX_GREY, lw=0.8, label="white-matter surface: the sources (coronal sections)")]
    fig.legend(handles=handles, loc="lower left", ncol=2, bbox_to_anchor=(lm / W, 0.0), handlelength=1.6,
               columnspacing=1.6, fontsize=7.5, borderaxespad=0.1)
    style.save(fig, "Figure_R12_geometry")

    g_fix = [vals[k]["top"]["median_gap_mm"] for k in vals]
    g_sc = [vals[k]["scaled_laterally_centred"]["median_gap_mm"] for k in vals]
    pc = g["config"]["placement"]
    standoffs = {g["arrays"][k]["opm_dense"]["standoff_mm"] for k in vals}
    if len(standoffs) != 1:
        raise SystemExit(f"the dense arrays' standoffs differ between the heads: {standoffs}")
    standoff = standoffs.pop()
    sc = {k: vals[k]["scaled_laterally_centred"] for k in vals}
    return dict(
        inputs=[f"{GEOMETRY} (written by scripts/export_g3b_geometry.py at {geo['provenance']['commit']} from the stored "
                f"state of the G3B run {geo['provenance']['source_state']['commit']}) :: anatomies[<head>].sections "
                "(sagittal.scalp, coronal.scalp, coronal.white.lh/rh: polylines, head frame, mm), coils_mm (fixed_top, "
                "scaled_x_centred), opm_dense_mm, transforms (top: device_to_head, moved_mm; counterfactual_x-centred: "
                "device_to_head, k, k_nominal), gaps_mm (checks), neuromag (magnetometers, coil_centres_device_mm: "
                "check), provenance (commit, source_state.commit)",
                f"{G3B} :: placements[<head>]['top'|'counterfactual_x-centred'|'x-centred'] (moved_mm, k, k_nominal, "
                "shift_x_mm, min/median/max_dist_mm: checks and labels), arrays[<head>].opm_dense (n_sites, "
                "standoff_mm), anatomies[<head>].head_size.ofc_mm, config.placement (clearance_mm, dewar_spacing_mm), "
                "provenance.commit"],
        description=(
            "Head-frame sections (Neuromag convention of each head's own fiducials; x right, y anterior, z up; mm) "
            "through the head origin: sagittal x = 0 and coronal y = 0. Lines: the intersection of the MRI scalp mesh "
            "(black) and, in the coronal sections, of the white-matter surfaces (grey; the cortical sources lie on them) "
            "with the plane, one segment per crossed triangle (the midsagittal plane passes between the hemispheres, so "
            "the sagittal sections show the scalp only), as polylines rounded to 0.01 mm; the sections are cropped below "
            f"z = {zl[0]:g} mm (face and neck). Markers: sensors within "
            f"{SLAB_MM:g} mm of the plane, projected onto it (as in the G3B geometry figure). Filled black squares: the "
            f"{n_coils} Neuromag magnetometer coil centres of the fixed adult helmet (VectorView geometry of the sample "
            "recording) with the head at top contact (the adult's measured pose, then raised along device +z until the "
            f"nearest coil is {pc['clearance_mm']:g} mm from the scalp; moved "
            + ", ".join(f"{GEOMETRY_NAMES[k]} {vals[k]['top']['moved_mm']:g} mm" for k in vals)
            + "). Open vermillion squares: the counterfactual helmet scaled with the head and laterally centred "
            "('counterfactual_x-centred': coil centres scaled by k about the head origin with the head at the measured "
            "pose shifted along device x to equal left/right median gaps, before any top contact; k = head-circumference "
            "ratio, raised in steps of 0.005 (opmsquid.pediatric.counterfactual_helmet) until no coil is closer than "
            f"{pc['dewar_spacing_mm']:g} mm to the scalp): k = " + ", ".join(f"{sc[k]['k']:.4f}" for k in vals)
            + " (nominal " + ", ".join(f"{sc[k]['k_nominal']:.4f}" for k in vals) + "), x shift "
            + ", ".join(f"{sc[k]['shift_x_mm']:+.1f} mm" for k in vals)
            + f". Open blue circles: the dense OPM sites (sensing centres, {standoff:g} mm from the scalp) refitted to "
            "each head (" + ", ".join(str(vals[k]["opm_dense_sites"]) for k in vals) + " sites). Order: adult, 12-month "
            "template (ANTS12-0Months3T), child B (sub-Z209, OpenNeuro ds005234). Median magnetometer-to-scalp gap "
            f"(nearest MRI scalp vertex, all {n_coils} coils), stored, and recomputed by the export from the unrounded "
            "geometry (equal within 1e-6 mm): fixed helmet at top contact " + ", ".join(f"{x:.2f}" for x in g_fix)
            + " mm; scaled, laterally centred helmet " + ", ".join(f"{x:.2f}" for x in g_sc) + " mm. Sections, coil "
            f"centres and OPM sites are read from {GEOMETRY}, exported by scripts/export_g3b_geometry.py from the stored "
            f"state of the G3B run of {G3B} (checked: run, moves, scale factors, gaps, site counts; the drawn coil centres "
            "equal the exported transforms applied to the device-frame coil centres within "
            f"{max(heads[k]['coil_err'] for k in heads):.4f} mm, the effect of the rounding); the coil centres are drawn, not "
            "the Dewar surface. Infant template: O'Reilly et al. (2021), built from the Neurodevelopmental MRI Database "
            "(Richards et al., 2016); child B: OpenNeuro ds005234 (Fadeev et al., 2024, 2025)."),
        alt=("Six panels: sagittal (top) and coronal (bottom) sections of three heads, adult, 12-month template and "
             "child B, drawn at the same scale. Each shows the scalp outline (the coronal ones with the folded "
             "white-matter outline inside), open blue circles of the dense OPM sites hugging the scalp, filled black "
             "squares of the fixed adult helmet and open vermillion squares of the helmet scaled with the head. Around "
             "the adult the two helmets nearly coincide; around the two smaller heads the fixed helmet touches at the top "
             f"and leaves a wide gap at the sides, front and back (median {g_fix[1]:.0f} and {g_fix[2]:.0f} mm against "
             f"the adult's {g_fix[0]:.0f} mm), while the scaled helmet follows the head ({g_sc[1]:.0f} and "
             f"{g_sc[2]:.0f} mm)."),
        caption_draft=(
            "Helmet geometry. Sagittal (top) and coronal (bottom) sections through the head origin of the adult, the "
            "12-month template and child B, at the same scale: scalp (black) and, in the coronal sections, white-matter "
            "surface (grey); filled squares, the Neuromag magnetometers of the fixed adult helmet with the head at top "
            "contact; open squares, the helmet scaled with the head and laterally centred (scale factor above each "
            f"column); circles, the dense OPM sites; only sensors within {SLAB_MM:g} mm of the section are drawn. Median "
            "magnetometer-to-scalp gap in the fixed helmet " + ", ".join(f"{x:.1f}" for x in g_fix) + " mm and in the "
            "scaled helmet " + ", ".join(f"{x:.1f}" for x in g_sc) + " mm (adult, 12-month template, child B)."),
        values=vals)


# ----------------------------------------------------------------------------------------------
def check_inputs() -> None:
    """Stop before drawing anything if an input is missing (nothing is recomputed in its place)."""
    paths.require(paths.SUBJECTS_DIR / "sample" / "bem" / "sample-oct-6-src.fif",
                  "MNE sample subject (set OPMSQUID_DATA)")
    for name in TEMPLATES.values():
        paths.require(paths.EXTERNAL / anatomy.INFANT_SUBJECTS / name, f"infant template {name} (set OPMSQUID_DATA)")
    paths.require(ROOT / GEOMETRY, "G3B geometry export (scripts/export_g3b_geometry.py)")


def main():
    style.apply()
    mne.set_log_level("WARNING")
    check_inputs()
    entries = {"Figure_R0_sphere": figure_sphere()}
    entries.update(figures_maps())
    entries["Figure_R12_geometry"] = figure_geometry()
    entries = {k: entries[k] for k in sorted(entries, key=lambda n: int(n.split("_")[1][1:]))}
    style.write_provenance(OUT_JSON, entries)
    for name in entries:
        print(f"results/report/{name}.png")
    print(f"results/report/{OUT_JSON}")


if __name__ == "__main__":
    main()
