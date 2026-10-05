#!/usr/bin/env python3
"""Report figures R5-R10 (children in the adult helmet, interictal spikes), drawn from stored
outputs only: nothing is re-analysed. Definitions are those of scripts/g3b_pediatric_helmet.py
(D, Delta, area-weighted medians without the medial wall, share_positive) and
scripts/g4_epilepsy_adult.py (paired S50 ratio, location sign-flip test).

  R5  D by lobe and for four parcels or parcel groups on the nine heads; fixed and scaled helmet
  R6  D_head and Delta = D_head - D_adult with 95 % intervals, both headline noise conditions; panel (b)
      names each class's Delta estimator (vertex-wise for the scaled adults, parcel-matched otherwise)
  R7  helmet fit vs head size: D at five helmet conditions, Neuromag gap, share of cortex with the
      OPM ahead, OPM site count
  R8  Delta per depth stratum
  R9  D against the OPM white-noise level
  R10 paired S50 ratio Neuromag / OPM per anatomy and depth band, Holm over the nine anatomies

Values come from results/g3b/g3b_summary.json and results/g4/g4_<anatomy>_summary.json; only
the four parcel rows of R5 are computed, from results/g3b/g3b_targets_<anatomy>.csv by the G3B
rule (area-weighted median: sort, cumulative area, first value reaching half the total), after
checking that the same rule on the same files reproduces the stored lobe medians.
Outputs: results/report/Figure_R5..R10_*.png and figures_pediatric.json (inputs, description,
alt text, caption draft and the plotted values of every figure).

Usage: PYTHONPATH=src .venv/bin/python scripts/report_figures_pediatric.py
"""
from __future__ import annotations

import csv
import json
import sys
import textwrap
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import Normalize, to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

from opmsquid import pediatric as P, plotting  # noqa: E402

G3B = "results/g3b/g3b_summary.json"
TARGETS = "results/g3b/g3b_targets_{}.csv"
G3A = "results/g3a/g3a_size_benchmark.json"
G4 = "results/g4/g4_{}_summary.json"
G4_REPORT = "results/g4/G4_pediatric_report.md"
G2_CFG = "configs/g2_adult.toml"
G2_SUMMARY = "results/g2/g2_summary.json"  # the projection's term count (the same noise model)
JAS = "docs/literature/jas2026.md"
OUT_JSON = "figures_pediatric.json"

CLASSES = ("adult", "scaled", "template", "child")
COND = "intrinsic+brain"
GREY = "0.45"
# R5 rows below the lobes: Desikan-Killiany parcels (both hemispheres) and the mesial temporal group
PARCELS = {"Precentral": ("precentral",), "Superior temporal": ("superiortemporal",),
           "Parahippocampal": ("parahippocampal",),
           "Mesial temporal\n(parahippocampal + entorhinal)": ("parahippocampal", "entorhinal")}
# the 12 source-blind placements in the fixed helmet (docs/methods.md: top contact and its variants)
FAMILY = ("top", "back", "x+5mm", "x-5mm", "y+5mm", "y-5mm", "pitch+10deg", "pitch-10deg", "roll+5deg", "roll-5deg",
          "yaw+10deg", "yaw-10deg")
HELMETS = (("centred", "centred\n(adult's\npose)"), ("top", "top\ncontact\n(primary)"),
           ("x-centred", "laterally\ncentred,\ntop contact"), ("counterfactual", "scaled\nwith the\nhead"),
           ("counterfactual_x-centred", "scaled,\nlaterally\ncentred"))
# approximate OPM empty-room floor of the Jas et al. (2026) recording: docs/literature/jas2026.md section 5.6
# ("about 28-35 fT/sqrt(Hz)", read from the preprint's Fig. 9 spectra; a broadband value, not 1-40 Hz)
FLOOR_FT = 30.0
G4_KEY = "{}/opm_vs_squid/combined/practical@1"
CM: dict = {}  # head circumference [cm] per anatomy, filled in main()


# ----------------------------------------------------------------------------------------------
# helpers
def load(rel: str):
    return json.loads((ROOT / rel).read_text())


def projection_terms() -> int:
    """Terms of the room-field projection of the 'projected' condition (the noise model's external basis), as the adult
    comparison stores them: retained rank with the room field minus retained rank after the projection."""
    r = load(G2_SUMMARY)["retained_rank"]
    return r["squid/intrinsic+brain+env"] - r["squid/projected"]


def order(s: dict) -> list[str]:
    """Anatomies by class (adult, scaled adult, infant template, individual child), then head
    circumference, largest first."""
    ofc = {k: s["anatomies"][k]["head_size"]["ofc_mm"] for k in style.ANAT_ORDER}
    return sorted(style.ANAT_ORDER, key=lambda k: (CLASSES.index(style.ANAT_CLASS[k]), -ofc[k]))


def label(k: str, short: bool = False, sep: str = ", ") -> str:
    """Anatomy name with its head circumference in cm."""
    return f"{(style.ANAT_SHORT if short else style.ANAT_LABEL)[k]}{sep}{CM[k]:.1f} cm"


def slots(keys, gap: float) -> dict:
    """Axis positions, one unit per anatomy, with ``gap`` added between classes."""
    pos, x, prev = {}, 0.0, None
    for k in keys:
        if prev is not None and style.ANAT_CLASS[k] != prev:
            x += gap
        pos[k], x, prev = x, x + 1.0, style.ANAT_CLASS[k]
    return pos


def groups(keys, pos) -> list[tuple[str, float, float]]:
    """(class, first, last position) of each run of one class."""
    out = []
    for k in keys:
        c = style.ANAT_CLASS[k]
        if out and out[-1][0] == c:
            out[-1] = (c, out[-1][1], pos[k])
        else:
            out.append((c, pos[k], pos[k]))
    return out


def mstyle(k: str, filled: bool = True, ms: float = 5.0, alpha: float = 1.0, color=None) -> dict:
    """Marker of an anatomy: class shape, anatomy colour; open = white face."""
    m = style.CLASS_MARKER[style.ANAT_CLASS[k]]
    c = style.ANAT_COLOR[k] if color is None else color
    size = ms * {"D": 0.85, "^": 1.1}.get(m, 1.0)
    return dict(marker=m, ms=size, mec=c, mfc=c if filled else "white", mew=0.9, alpha=alpha, ls="none")


def tint(c, f: float = 0.55):
    """Colour blended towards white by ``f``."""
    r = np.array(to_rgb(c))
    return tuple(r + (1 - r) * f)


def db(v: float) -> str:
    return f"{v:+.2f}".replace("-", "−")


def fmt_p(p: float) -> str:
    if p >= 0.001:
        return f"{p:.3g}"
    m, e = f"{p:.1e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m}×10{str(int(e)).translate(sup)}"


def holm(p) -> np.ndarray:
    """Holm step-down adjusted p values over the family given, in input order."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.maximum.accumulate(np.minimum(1.0, (len(p) - np.arange(len(p))) * p[o]))
    out = np.empty(len(p))
    out[o] = adj
    return out


def read_targets(k: str) -> dict:
    """Columns of results/g3b/g3b_targets_<k>.csv (the status line skipped)."""
    with open(ROOT / TARGETS.format(k), newline="") as fh:
        rows = list(csv.reader(line for line in fh if not line.startswith("#")))
    head, body = rows[0], rows[1:]
    return {h: np.array([r[i] for r in body], dtype=str if h in ("region", "lobe") else float) for i, h in enumerate(head)}


def anat_legend(keys, ms: float = 5.0) -> list:
    return [Line2D([], [], color=style.ANAT_COLOR[k], lw=0.9, label=label(k), **dict(mstyle(k, ms=ms), ls="-")) for k in keys]


def hline_label(ax, x, y, text, **kw):
    ax.text(x, y, text, fontsize=6.5, color=GREY, **kw)


def footnote(fig, text: str, y: float = 0.0) -> None:
    """Small grey note under the panels, wrapped to the figure's design width."""
    fig.text(0.0, y, textwrap.fill(text, 135), ha="left", va="top", fontsize=6.5, color="0.3")


def n_cortical(s: dict, keys) -> dict:
    """Cortical targets per head (all targets minus the medial wall), the n of every G3B median."""
    return {k: s["anatomies"][k]["n_targets"] - s["medial_wall_targets"][k] for k in keys}


def n_text(s: dict, keys) -> str:
    n = n_cortical(s, keys)
    kids = [n[k] for k in keys if k != "adult"]
    return f"{min(kids):,}\u2013{max(kids):,} per smaller head, {n['adult']:,} adult"


# ----------------------------------------------------------------------------------------------
# R5 regions x head size
def region_values(s: dict, keys) -> tuple[dict, dict, float]:
    """Area-weighted median D (dense OPM vs Neuromag combined, intrinsic + brain, top contact) per
    parcel group from the target files, targets per cell (parcels and lobes), and the largest
    difference between the same rule's lobe medians and the stored ones (a check of the columns)."""
    val, n, check = {}, {}, 0.0
    for k in keys:
        t = read_targets(k)
        cortical = ~np.char.endswith(t["region"], "unknown")  # the medial wall is left out (G3B)
        d = np.where(cortical, 20 * np.log10(t[f"detect_opm_dense_opm_{COND}"] / t[f"detect_squid_top_combined_{COND}"]), np.nan)
        parcel = np.array([r.split(".", 1)[1] for r in t["region"]])
        for name, members in PARCELS.items():
            m = cortical & np.isin(parcel, members)
            val[name, k], n[name, k] = P.weighted_median(d[m], t["area_mm2"][m]), int(m.sum())
        for lb, v in s["placement_D"][f"{k}/top/combined/{COND}"]["by_lobe"].items():
            m = cortical & (t["lobe"] == lb)
            check = max(check, abs(P.weighted_median(d[m], t["area_mm2"][m]) - v))
            n[lb.capitalize(), k] = int(m.sum())
    return val, n, check


def column_header(ax, keys, xs) -> None:
    """Above a heat map: class brackets, the anatomy names with head circumference, class markers."""
    tr = ax.get_xaxis_transform()
    ax.set_xticks([xs[k] for k in keys], [label(k, short=True, sep="\n") for k in keys], fontsize=6.8)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=13)
    glyph = offset_copy(tr, fig=ax.figure, y=7, units="points")
    bracket = offset_copy(tr, fig=ax.figure, y=40, units="points")
    for k in keys:
        ax.plot([xs[k]], [1], transform=glyph, clip_on=False, **mstyle(k, ms=4.5))
    for c, a, b in groups(keys, xs):
        ax.plot([a - 0.45, b + 0.45], [1, 1], transform=bracket, clip_on=False, color="0.65", lw=0.7)
        ax.text(0.5 * (a + b), 1, style.CLASS_LABEL[c], transform=offset_copy(bracket, fig=ax.figure, y=2, units="points"),
                ha="center", va="bottom", fontsize=7, color="0.3")


def heat_cells(ax, cells: dict, xs: dict, ys: dict, norm, cmap) -> None:
    """Coloured cells with their values (dB); white text on dark cells."""
    for (r, k), v in cells.items():
        rgba = cmap(norm(v))
        ax.add_patch(Rectangle((xs[k] - 0.5, ys[r] - 0.5), 1, 1, fc=rgba, ec="white", lw=1.2))
        lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        ax.text(xs[k], ys[r], db(v), ha="center", va="center", fontsize=7, color="white" if lum < 0.5 else "black")


def fig_r5(s: dict, keys, reg: dict, n: dict, check: float) -> dict:
    lobes = [lb.capitalize() for lb in plotting.DK_LOBES]
    rows = lobes + list(PARCELS)
    xs = slots(keys, gap=0.22)
    ys = {r: i + (0.4 if i >= len(lobes) else 0.0) for i, r in enumerate(rows)}
    panels = (("top", "(a)  Fixed adult helmet, top contact (primary placement)"),
              ("counterfactual_x-centred", "(b)  Helmet scaled with the head, laterally centred"))
    vals = {}
    for pl, _ in panels:
        for k in keys:
            by = s["placement_D"][f"{k}/{pl}/combined/{COND}"]["by_lobe"]
            vals.update({(pl, lb.capitalize(), k): by[lb] for lb in plotting.DK_LOBES})
    vals.update({("top", r, k): reg[r, k] for r in PARCELS for k in keys})
    vmax = float(np.ceil(max(abs(v) for v in vals.values()) * 2) / 2)
    norm, cmap = Normalize(-vmax, vmax), plt.get_cmap("RdBu_r")
    h_a, h_b = ys[rows[-1]] + 1, len(lobes) + 1.4  # panel (b): the lobes and one strip for the parcels not stored
    fig, axs = plt.subplots(2, 1, figsize=(style.FULL_W, 7.3), gridspec_kw=dict(hspace=0.5, height_ratios=[h_a, h_b], left=0.215,
                                                                                right=0.9, top=0.9, bottom=0.1))
    x_lo, x_hi = min(xs.values()) - 0.5, max(xs.values()) + 0.5
    for ax, (pl, title), h in zip(axs, panels, (h_a, h_b)):
        heat_cells(ax, {(r, k): v for (p_, r, k), v in vals.items() if p_ == pl}, xs, ys, norm, cmap)
        shown = rows if pl == "top" else lobes
        if pl != "top":  # the per-target files hold top, centred and the unshifted scaled helmet only
            y0 = len(lobes) - 0.5 + 0.4
            ax.add_patch(Rectangle((x_lo, y0), x_hi - x_lo, 1.0, fc="#f4f4f4", ec="#c8c8c8", hatch="////", lw=0.6))
            ax.text(0.5 * (x_lo + x_hi), y0 + 0.5, "precentral, superior temporal, parahippocampal, mesial temporal: "
                    "not stored for this helmet", ha="center", va="center", fontsize=7, color="0.3",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
            ax.text(x_lo - 0.12, y0 + 0.5, "Parcels", ha="right", va="center", fontsize=7.2)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(h - 0.5, -0.5)
        ax.set_yticks([ys[r] for r in shown], shown, fontsize=7.2)
        ax.tick_params(axis="y", length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        column_header(ax, keys, xs)
        ax.set_title(title, loc="left", fontsize=9, pad=60)
        for r in shown:  # targets per cell (the same target set for both helmets)
            nn = [n[r, k] for k in keys]
            ax.text(x_hi + 0.15, ys[r], f"{min(nn):,}\u2013{max(nn):,}", ha="left", va="center", fontsize=6.5, color=GREY)
        ax.text(x_hi + 0.15, -0.6, "targets\nper cell", ha="left", va="bottom", fontsize=6.5, color=GREY)
    cax = fig.add_axes([0.30, 0.04, 0.45, 0.013])
    cb = fig.colorbar(ScalarMappable(norm, cmap), cax=cax, orientation="horizontal")
    cb.set_label("D = dense OPM minus Neuromag (dB)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    cax.text(-0.02, 0.5, "Neuromag ahead", transform=cax.transAxes, ha="right", va="center", fontsize=7, color="0.3")
    cax.text(1.02, 0.5, "OPM ahead", transform=cax.transAxes, ha="left", va="center", fontsize=7, color="0.3")
    style.save(fig, "Figure_R5_regions_heads")
    allv = list(vals.values())
    cells = {pl: {r: {k: round(vals[pl, r, k], 4) for k in keys if (pl, r, k) in vals} for r in rows} for pl, _ in panels}
    return dict(order=keys, head_circumference_cm={k: round(CM[k], 2) for k in keys}, D_dB=cells,
                targets_per_cell={r: {k: n[r, k] for k in keys} for r in rows}, colour_limit_dB=vmax,
                min_cell_dB=round(min(allv), 4), max_cell_dB=round(max(allv), 4),
                lobe_check_max_abs_diff_dB=round(check, 6))


# ----------------------------------------------------------------------------------------------
# R6 children in the adult helmet
def adult_reference(s: dict, kids, arr: str, cond: str, field: str = "median", part: str = "d_adult") -> float:
    """The adult's value in every child comparison (identical by construction; checked)."""
    v = {s["comparisons"][f"{k}/{arr}/combined/{cond}/detect"][part][field] for k in kids}
    if len(v) != 1:
        raise ValueError(f"adult {field} differs between comparisons ({arr}, {cond}): {v}")
    return v.pop()


def row_axes(ax, keys, pos, labels=True, short=False):
    """Anatomies as rows (top to bottom), labelled on the axes that carries the labels, with a
    hairline between classes."""
    if labels:
        ax.set_yticks([-pos[k] for k in keys], [label(k, short=short) for k in keys])
    ax.tick_params(axis="y", length=0, labelleft=labels)
    ax.set_ylim(-max(pos.values()) - 0.7, 0.7)
    for (_, _, b), (_, a, _) in zip(groups(keys, pos)[:-1], groups(keys, pos)[1:]):
        ax.axhline(-0.5 * (b + a), color="0.88", lw=0.7, zorder=0)


def fig_r6(s: dict, keys) -> dict:
    kids = [k for k in keys if k != "adult"]
    pos = slots(kids, gap=0.5)
    ref = {c: adult_reference(s, kids, "opm_dense", c) for c in (COND, "projected")}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(style.FULL_W, 3.5), sharey=True,
                                 gridspec_kw=dict(wspace=0.07, left=0.315, right=0.985, top=0.8, bottom=0.16))
    out = {}
    for k in kids:
        for cond, dy, filled in ((COND, 0.16, True), ("projected", -0.16, False)):
            e = s["comparisons"][f"{k}/opm_dense/combined/{cond}/detect"]
            y = -pos[k] + dy
            for ax, q in ((a1, e["d_child"]), (a2, e["delta"])):
                ax.plot(q["ci95"], [y, y], color=style.ANAT_COLOR[k], lw=1.0, solid_capstyle="butt")
                ax.plot(q["median"], y, **mstyle(k, filled))
            out[f"{k}/{cond}"] = dict(d_child=e["d_child"]["median"], d_child_ci95=e["d_child"]["ci95"], delta=e["delta"]["median"],
                                      delta_ci95=e["delta"]["ci95"])
    for ax in (a1, a2):
        row_axes(ax, kids, pos, labels=ax is a1)
    a1.axvline(ref[COND], color="k", lw=0.9, zorder=0)
    a1.axvline(ref["projected"], color="k", lw=0.9, ls=(0, (4, 2)), zorder=0)
    a1.axvline(0, color="0.75", lw=0.7, zorder=0)
    a2.axvline(0, color="0.5", lw=0.8, zorder=0)
    a1.set_xlim(-0.15, 2.65)
    a2.set_xlim(-0.55, 1.65)
    a1.set_xlabel("D$_\\mathrm{head}$ = dense OPM minus Neuromag (dB)")
    a2.set_xlabel("Δ = D$_\\mathrm{head}$ − D$_\\mathrm{adult}$ (dB)")
    a1.set_title("(a)  OPM advantage per head", loc="left")
    a2.set_title("(b)  Change from the adult", loc="left")
    est = {}  # (b): the estimator of each class's Delta, named in the figure (Referee 1, minor 14)
    for c, a, b in groups(kids, pos):
        a1.text(-0.01, -a + 0.5, style.CLASS_LABEL[c], transform=a1.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=6.8, color=GREY, style="italic")
        d = [s["comparisons"][f"{k}/opm_dense/combined/{COND}/detect"]["delta"] for k in kids if style.ANAT_CLASS[k] == c]
        if all("n_parcels" in x for x in d):  # no vertex correspondence: parcel differences (stored 'method')
            n_p = sorted({x["n_parcels"] for x in d})
            rng = str(n_p[0]) if len(n_p) == 1 else f"{n_p[0]}\u2013{n_p[-1]}"
            est[c] = f"\u0394 parcel-matched ({rng} parcels)"
        elif not any("n_parcels" in x for x in d) and c == "scaled":  # the adult's own vertices
            est[c] = "\u0394 vertex-wise (the same vertices)"
        else:
            raise ValueError(f"{c}: mixed or unexpected Delta estimators")
        a2.text(0.985, -a + 0.5, est[c], transform=a2.get_yaxis_transform(), ha="right", va="bottom", fontsize=6.8,
                color=GREY, style="italic")
    h = [Line2D([], [], **mstyle("childC", True, color="0.3")), Line2D([], [], **mstyle("childC", False, color="0.3")),
         Line2D([], [], color="k", lw=0.9), Line2D([], [], color="k", lw=0.9, ls=(0, (4, 2)))]
    fig.legend(h, ["sensor plus brain noise (filled)", f"after the {projection_terms()}-term projection (open)",
                   f"adult D, sensor plus brain noise ({db(ref[COND])} dB)", f"adult D, after the projection ({db(ref['projected'])} dB)"],
               loc="upper center", ncol=2, bbox_to_anchor=(0.62, 1.0), handlelength=1.8, columnspacing=1.5)
    primary = {k: s["comparisons"][f"{k}/opm_dense/combined/{COND}/detect"] for k in kids}
    if any(primary[k]["d_child"]["n"] != n_cortical(s, kids)[k] for k in kids):
        raise ValueError("D_child n differs from the cortical target count")
    parc = [primary[k]["delta"]["n_parcels"] for k in kids if "n_parcels" in primary[k]["delta"]]
    footnote(fig, f"D: area-weighted median over the cortical targets ({n_text(s, keys)}). Lines: 95 % intervals, parcel "
                  f"bootstrap ({s['config']['strata']['n_boot']:,} resamples); \u0394 vertex-wise for the scaled adults, over "
                  f"{min(parc)}\u2013{max(parc)} parcels for the templates and children.", y=0.035)
    style.save(fig, "Figure_R6_pediatric_D")
    return dict(order=kids, adult_D_dB={c: ref[c] for c in ref}, values=out, delta_estimator_labels=est)


# ----------------------------------------------------------------------------------------------
# R7 helmet fit vs head size
def fig_r7(s: dict, keys) -> dict:
    kids = [k for k in keys if k != "adult"]
    pd_, sd, comp, cc = s["placement_D"], s["sensor_distances"], s["comparisons"], s["channel_count_control"]
    xi = {v["xi_squid_mm"] for v in load(G3A)["size_following"].values()}
    if len(xi) != 1:
        raise ValueError(f"{G3A}: the sphere model's SQUID gap differs between heads: {xi}")
    xi = xi.pop()
    fam = [pd_[f"adult/{n}/combined/{COND}"]["median"] for n in FAMILY if n not in s["infeasible_placements"]["adult"]]
    a_top, a_cfx = pd_[f"adult/top/combined/{COND}"]["median"], pd_[f"adult/counterfactual_x-centred/combined/{COND}"]["median"]
    if abs(a_top - adult_reference(s, kids, "opm_dense", COND)) > 1e-12:
        raise ValueError("the adult's top-contact D differs between placement_D and comparisons")
    fig = plt.figure(figsize=(style.FULL_W, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.11)
    ax_a, ax_b, ax_c, ax_d = (fig.add_subplot(gs[i, j]) for i in (0, 1) for j in (0, 1))
    out = dict(order=keys, adult_family_D_dB=dict(n=len(fam), min=min(fam), max=max(fam)), adult_top_D_dB=a_top,
               adult_counterfactual_x_centred_D_dB=a_cfx, sphere_gap_mm=xi)

    # (a) D at five helmet conditions
    xs = [0.0, 1.0, 2.0, 3.45, 4.45]
    ax_a.fill_between([-0.35, 2.35], min(fam), max(fam), color="0.9", lw=0, zorder=0)
    ax_a.hlines(a_top, -0.35, 2.35, color="k", lw=0.8, ls=(0, (4, 2)), zorder=1)
    ax_a.hlines(a_cfx, 3.1, 4.8, color="k", lw=0.8, ls=(0, (4, 2)), zorder=1)
    for k in keys:
        v = [pd_[f"{k}/{pl}/combined/{COND}"]["median"] for pl, _ in HELMETS]
        for sl in (slice(0, 3), slice(3, 5)):
            ax_a.plot(xs[sl], v[sl], color=style.ANAT_COLOR[k], lw=0.9, alpha=0.85, zorder=2)
        ax_a.plot(xs, v, zorder=3, **mstyle(k, ms=4.5))
        out[f"{k}/D_by_helmet_dB"] = dict(zip([pl for pl, _ in HELMETS], v))
    ax_a.set_xticks(xs, [t for _, t in HELMETS], fontsize=6.5)
    ax_a.set_xlim(-0.45, 4.9)
    ax_a.set_ylim(0, 3.9)
    ax_a.axhline(0, color="0.75", lw=0.7, zorder=0)
    ax_a.set_ylabel("D, dense OPM minus Neuromag (dB)")
    ax_a.text(1.0, 1.02, "fixed adult helmet", transform=ax_a.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color="0.3")
    ax_a.text(3.95, 1.02, "helmet scaled with the head", transform=ax_a.get_xaxis_transform(), ha="center", va="bottom", fontsize=7,
              color="0.3")
    hline_label(ax_a, -0.3, 0.12, f"grey band: adult over its {len(fam)} standard\nplacements ({db(min(fam))} to {db(max(fam))} dB)",
                ha="left", va="bottom")
    hline_label(ax_a, 2.33, a_top - 0.05, f"adult {db(a_top)}", ha="right", va="top")
    hline_label(ax_a, 4.85, a_cfx, f"adult\n{db(a_cfx)}", ha="left", va="center")
    ax_a.set_title("(a)  OPM advantage at five helmet conditions", loc="left", pad=14)

    # (b) Neuromag gap at top contact and in the scaled, laterally centred helmet
    pos = slots(keys, gap=0.5)
    for k in keys:
        g0, g1 = sd[k]["squid:top"]["median_mm"], sd[k]["squid:counterfactual_x-centred"]["median_mm"]
        ax_b.plot([pos[k]] * 2, [g0, g1], color=style.ANAT_COLOR[k], lw=0.8, alpha=0.7)
        ax_b.plot(pos[k], g0, **mstyle(k, True))
        ax_b.plot(pos[k], g1, **mstyle(k, False))
        out[f"{k}/gap_mm"] = dict(top=g0, counterfactual_x_centred=g1)
    ax_b.axhline(xi, color="0.35", lw=0.8, ls=(0, (1, 1.5)))
    hline_label(ax_b, max(pos.values()) + 0.4, xi + 0.6, f"sphere model gap, {xi:g} mm (Jas et al. 2026)", ha="right", va="bottom")
    ax_b.set_ylim(14, 44)
    ax_b.set_ylabel("median coil-to-scalp gap (mm)")
    ax_b.set_title("(b)  Neuromag gap", loc="left", pad=14)
    ax_b.legend([Line2D([], [], **mstyle("childC", True, color="0.3")), Line2D([], [], **mstyle("childC", False, color="0.3"))],
                ["top contact (primary)", "scaled, laterally centred"], loc="upper left", ncol=2, handlelength=1.0,
                columnspacing=1.0, borderaxespad=0.2)

    # (c) share of cortex with the OPM ahead
    for arr, alpha in (("opm_dense", 0.3), ("opm_matched", 1.0)):
        for cond, dx, filled in ((COND, -0.14, True), ("projected", 0.14, False)):
            for k in keys:
                v = (adult_reference(s, kids, arr, cond, "share_positive") if k == "adult"
                     else comp[f"{k}/{arr}/combined/{cond}/detect"]["d_child"]["share_positive"])
                ax_c.plot(pos[k] + dx, v, zorder=3 if alpha == 1 else 2, **mstyle(k, filled, ms=4.5, alpha=alpha))
                out[f"{k}/share_positive/{arr}/{cond}"] = v
    ax_c.axhline(0.5, color="0.8", lw=0.7, zorder=0)
    ax_c.set_ylim(0, 1.05)
    ax_c.set_ylabel("share of cortical area with D > 0")
    ax_c.set_title("(c)  Share of cortex where the OPM is ahead", loc="left", pad=14)
    ax_c.legend([Line2D([], [], **mstyle("childC", True, color="0.3")), Line2D([], [], **mstyle("childC", False, color="0.3")),
                 Line2D([], [], **mstyle("childC", True, color="0.3", alpha=0.3))],
                ["site-matched OPM, sensor plus brain noise", f"site-matched OPM, after the {projection_terms()}-term projection",
                 "dense OPM (faint)"], loc="lower left", ncol=1, handlelength=1.0, borderaxespad=0.2, fontsize=6.8)

    # (d) OPM site count and the adult's loss at that count
    a_n = s["arrays"]["adult"]["opm_dense"]["n"]
    ax_d.axhline(a_top, color="k", lw=0.8, ls=(0, (4, 2)), zorder=0)
    ax_d.plot(pos["adult"], a_top, **mstyle("adult"))
    ax_d.text(pos["adult"], a_top - 0.05, f"{a_n}", ha="center", va="top", fontsize=6.5, color="0.25")
    for k in kids:
        e = cc[f"{k}/combined"]
        if e["n_sites"] != s["arrays"][k]["opm_dense"]["n"]:
            raise ValueError(f"{k}: channel_count_control site count differs from the array")
        v = e["d_adult_subsampled"]["median"]
        ax_d.plot([pos[k]] * 2, [a_top, v], color=style.ANAT_COLOR[k], lw=0.8, alpha=0.7)
        ax_d.plot(pos[k], v, **mstyle(k))
        ax_d.text(pos[k], v - 0.05, f"{e['n_sites']}", ha="center", va="top", fontsize=6.5, color="0.25")
        out[f"{k}/channel_count"] = dict(n_sites=e["n_sites"], adult_subsampled_D_dB=v, adult_loss_dB=a_top - v)
    hline_label(ax_d, max(pos.values()) + 0.4, a_top + 0.03, f"adult, all {a_n} sites ({db(a_top)} dB)", ha="right", va="bottom")
    ax_d.set_ylim(0, 1.25)
    ax_d.set_ylabel("adult D, subsampled (dB)")
    ax_d.set_title("(d)  Dense OPM sites (numbers), the adult's loss", loc="left", pad=14)
    for ax in (ax_b, ax_c, ax_d):
        ax.set_xticks([pos[k] for k in keys], [style.ANAT_SHORT[k] for k in keys], rotation=40, ha="right", fontsize=6.8)
        ax.set_xlim(-0.6, max(pos.values()) + 0.6)
        ax.tick_params(axis="x", length=2)
    fig.legend(handles=anat_legend(keys, ms=4.5), loc="upper center", ncol=3, bbox_to_anchor=(0.53, 1.0), handlelength=2.0,
               columnspacing=1.2, fontsize=7)
    footnote(fig, f"D and shares: area-weighted over the cortical targets ({n_text(s, keys)}). Gap: median over the "
                  "magnetometer coil centres of the distance to the nearest scalp point.", y=0.012)
    style.save(fig, "Figure_R7_helmet_fit")
    return out


# ----------------------------------------------------------------------------------------------
# R8 depth-matched Delta
def fig_r8(s: dict, keys) -> dict:
    kids = [k for k in keys if k != "adult"]
    left = [k for k in kids if style.ANAT_CLASS[k] != "child"]
    right = [k for k in kids if style.ANAT_CLASS[k] == "child"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(style.FULL_W, 3.3), sharey=True, gridspec_kw=dict(wspace=0.06))
    out, n_child, n_adult = {}, [], []
    for ax, ks, ci in ((a1, left, False), (a2, right, True)):
        for i, k in enumerate(ks):
            rows = [r for r in s["comparisons"][f"{k}/opm_dense/combined/{COND}/detect"]["delta_by_depth"] if "delta" in r]
            x = np.array([0.5 * (r["lo"] + r["hi"]) for r in rows]) + (i - (len(ks) - 1) / 2) * (0.7 if ci else 0.0)
            y = np.array([r["delta"] for r in rows])
            ax.plot(x, y, color=style.ANAT_COLOR[k], lw=1.0, label=label(k), **dict(mstyle(k, ms=4.5), ls="-"))
            if ci:
                lo, hi = np.array([r["ci95"] for r in rows]).T
                ax.vlines(x, lo, hi, color=style.ANAT_COLOR[k], lw=0.9)
            out[k] = [dict(lo=r["lo"], hi=r["hi"], n_child=r["n_child"], n_adult=r["n_adult"], delta=r["delta"], ci95=r["ci95"])
                      for r in rows]
            n_child += [r["n_child"] for r in rows]
            n_adult += [r["n_adult"] for r in rows]
        ax.axhline(0, color="0.5", lw=0.8, zorder=0)
        ax.set_xlim(7, 92)
        ax.set_xticks(range(10, 100, 10))
        ax.legend(loc="upper center" if ax is a1 else "upper left", handlelength=2.0, fontsize=7)
    edges = s["config"]["strata"]["depth_edges_mm"]
    for ax in (a1, a2):  # stratum edges between the labelled ticks
        ax.xaxis.set_minor_locator(FixedLocator([e for e in edges if e % 10]))
        ax.tick_params(axis="x", which="minor", length=2)
    a1.set_ylim(-0.75, 2.55)
    fig.supxlabel("depth below the scalp (mm; points at stratum centres)", fontsize=8.5, y=-0.02)
    a1.set_ylabel("Δ in the stratum: D$_\\mathrm{head}$ − D$_\\mathrm{adult}$\nof the stratum medians (dB)")
    a1.set_title("(a)  Scaled adults and infant templates", loc="left")
    a2.set_title("(b)  Individual children (95 % parcel-bootstrap intervals)", loc="left")
    footnote(fig, f"Targets per stratum: {min(n_child):,}\u2013{max(n_child):,} (smaller heads), {min(n_adult):,}\u2013"
                  f"{max(n_adult):,} (adult); strata with fewer than {s['config']['strata']['min_n']} in either head omitted. "
                  f"Intervals: 95 %, parcels resampled in each anatomy ({s['config']['strata']['n_boot']:,} resamples).", y=-0.07)
    style.save(fig, "Figure_R8_depth_matched")
    return dict(order=kids, depth_edges_mm=edges, min_n=s["config"]["strata"]["min_n"], strata=out,
                n_child_range=[min(n_child), max(n_child)], n_adult_range=[min(n_adult), max(n_adult)])


# ----------------------------------------------------------------------------------------------
# R9 noise floor
def fig_r9(s: dict, keys) -> dict:
    sens = s["sensitivity_median_D_dB"]
    g2 = tomllib.loads((ROOT / G2_CFG).read_text())["sensors"]
    levels, primary = [float(a) for a in g2["opm_asd_fT_per_rtHz"]], float(g2["opm_asd_primary_fT_per_rtHz"])
    fig, axs = plt.subplots(1, 2, figsize=(style.FULL_W, 3.9), sharey=True, gridspec_kw=dict(wspace=0.06))
    out = dict(order=keys, levels_fT_per_rtHz=levels, primary_fT_per_rtHz=primary, floor_fT_per_rtHz=FLOOR_FT)
    n_sq = {s["arrays"][k]["squid:top"]["n"] for k in keys}
    if len(n_sq) != 1:
        raise ValueError(f"Neuromag channel counts differ between heads: {n_sq}")
    n_sq = n_sq.pop()
    for ax, arr, title in ((axs[0], "opm_dense", "(a)  Dense OPM array"), (axs[1], "opm_matched", "(b)  Site-matched OPM array")):
        for k in keys:
            v = [sens[f"{k}/opm_asd_{a:g}fT/{arr}/combined/{COND}"] for a in levels]
            prim = s["D_median_dB"][f"{k}/{arr}/combined/{COND}/detect"]
            if abs(v[levels.index(primary)] - prim) > 1e-9:
                raise ValueError(f"{k}/{arr}: the stored {primary:g}-fT sensitivity differs from the primary D")
            ax.plot(levels, v, color=style.ANAT_COLOR[k], lw=1.0, **dict(mstyle(k, ms=4.5), ls="-"))
            out[f"{k}/{arr}"] = dict(zip([f"{a:g}" for a in levels], v))
        ax.axhline(0, color="0.5", lw=0.8, zorder=0)
        ax.axvline(primary, color="0.7", lw=0.8, zorder=0)
        ax.axvline(FLOOR_FT, color="0.35", lw=0.8, ls=(0, (1, 1.5)), zorder=0)
        ax.set_xticks(levels, [f"{a:g}" for a in levels])
        ax.set_xlim(5.5, 31.5)
        ax.set_xlabel("OPM white noise (fT/√Hz)")
        ax.set_title(title, loc="left", pad=22)
        ax.text(primary, 1.01, "primary", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5, color=GREY)
        ax.text(FLOOR_FT, 1.01, "≈ OPM floor,\nJas et al. 2026", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=6.5, color=GREY)
    axs[0].set_ylabel("D, OPM minus Neuromag (dB)")
    axs[0].set_ylim(-1.0, 3.4)
    fig.legend(handles=anat_legend(keys, ms=4.5), loc="upper center", ncol=3, bbox_to_anchor=(0.52, 0.0), handlelength=2.0,
               columnspacing=1.2, fontsize=7)
    footnote(fig, f"D: OPM array minus Neuromag ({n_sq} channels), area-weighted median over the cortical targets "
                  f"({n_text(s, keys)}); Neuromag noise unchanged.", y=-0.16)
    style.save(fig, "Figure_R9_noise_floor")
    return out


# ----------------------------------------------------------------------------------------------
# R10 interictal spikes
def g4_bands() -> list[str]:
    """Depth-band labels depth0..depth3, from the table header of results/g4/G4_pediatric_report.md."""
    for line in (ROOT / G4_REPORT).read_text().splitlines():
        if line.startswith("| anatomy | detector |"):
            return [c.strip().replace(" mm", "") for c in line.strip("|").split("|")[2:]]
    raise ValueError(f"{G4_REPORT}: no S50 table header")


def check_bands(g: dict, bands: list[str], k: str) -> None:
    """Every stored location's depth lies in the band its stratum names."""
    for loc in g["locations"]:
        lo, hi = (float(x) for x in bands[loc["stratum"][0]].split("-"))
        if not lo <= loc["depth_mm"] < hi:
            raise ValueError(f"{k}: location at {loc['depth_mm']:.1f} mm outside band {bands[loc['stratum'][0]]}")


def arrow(ax, x0, x1, y, color, ls="-"):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0, ls=ls, mutation_scale=7,
                                                                shrinkA=0, shrinkB=0))


def draw_ratio(ax, y, sr: dict, color, marker: dict, lims) -> str:
    """One S50 ratio with its 95 % interval. An open interval end (resamples outside the tested
    strengths) is an arrow to the axis edge, an interval open at both ends a dotted line; a censored
    point estimate (only bounded) is an open marker at its bound with a dashed arrow beyond it;
    nothing is drawn when neither system reaches 50 % (the caller writes it)."""
    lo, hi = sr["ci95"] if sr["ci95"] else (None, None)
    v = sr["value"]
    if v is not None:
        if lo is None and hi is None:
            ax.plot(lims, [y, y], color=color, lw=0.8, ls=(0, (1, 1.5)))
            arrow(ax, lims[0] * 1.06, lims[0], y, color)
            arrow(ax, lims[1] / 1.06, lims[1], y, color)
        else:
            a, b = lo if lo is not None else lims[0], hi if hi is not None else lims[1]
            ax.plot([a, b], [y, y], color=color, lw=1.0, solid_capstyle="butt")
            if lo is None:
                arrow(ax, a * 1.06, lims[0], y, color)
            if hi is None:
                arrow(ax, b / 1.06, lims[1], y, color)
        ax.plot(v, y, **marker)
        return "value" if lo is not None and hi is not None else "value, interval open"
    b_lo, b_hi = sr["value_bounds"] or (None, None)
    if b_lo is not None or b_hi is not None:
        start, edge = (b_lo, lims[1]) if b_lo is not None else (b_hi, lims[0])
        arrow(ax, start, edge, y, color, ls=(0, (2, 1.5)))
        ax.plot(start, y, **dict(marker, mfc="white"))
        return "lower bound" if b_lo is not None else "upper bound"
    return "neither reaches 50 %"


def fig_r10(s: dict, keys) -> dict:
    bands = g4_bands()
    g4 = {k: load(G4.format(k)) for k in keys}
    for k in keys:
        check_bands(g4[k], bands, k)
    strengths = {tuple(g4[k]["config"]["events"]["strengths_nAm"]) for k in keys}
    if len(strengths) != 1:
        raise ValueError(f"G4 strengths differ between anatomies: {strengths}")
    strengths = strengths.pop()
    n_loc = {r["n_locations"] for k in keys for arr in ("opm_dense", "opm_matched")
             for r in g4[k]["paired"][G4_KEY.format(arr)].values()}
    n_ev = {r["n"] for k in keys for arr in ("opm_dense", "opm_matched") for r in g4[k]["paired"][G4_KEY.format(arr)].values()}
    if len(n_loc) != 1 or len(n_ev) != 1:
        raise ValueError(f"G4 locations or events per band differ: {n_loc}, {n_ev}")
    n_loc, n_ev = n_loc.pop(), n_ev.pop()
    pairs = {(k, arr): g4[k]["paired"][G4_KEY.format(arr)] for k in keys for arr in ("opm_dense", "opm_matched")}
    p_holm = {arr: dict(zip(keys, holm([pairs[k, arr]["depth0"]["location_sign_flip_p"] for k in keys])))
              for arr in ("opm_dense", "opm_matched")}
    pos = slots(keys, gap=0.5)
    fig = plt.figure(figsize=(style.FULL_W, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.36, left=0.305, right=0.99, top=0.9, bottom=0.1)
    top = gs[0].subgridspec(1, 2, width_ratios=[2.1, 0.9], wspace=0.03)
    ax, ax_t = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    bot = gs[1].subgridspec(1, len(bands) - 1, wspace=0.08)
    axs_b = [fig.add_subplot(bot[i]) for i in range(len(bands) - 1)]
    lims_a, lims_b = (0.7, 2.0), (0.7, 1.6)
    out = dict(order=keys, bands_mm=bands, strengths_nAm=list(strengths), ratios={},
               family="Holm over the nine anatomies at 10-20 mm (depth0), each OPM array separately",
               anatomy_notes={k: g4[k]["anatomy"] for k in keys})
    for k in keys:
        for arr, dy, light in (("opm_dense", 0.2, False), ("opm_matched", -0.2, True)):
            c = tint(style.ANAT_COLOR[k]) if light else style.ANAT_COLOR[k]
            tc = "0.5" if light else "0.1"
            for b, a in enumerate([ax] + axs_b):
                r = pairs[k, arr][f"depth{b}"]
                lims = lims_a if b == 0 else lims_b
                kind = draw_ratio(a, -pos[k] + dy, r["s50_ratio_squid_over_opm"], c, mstyle(k, color=c, ms=4.5), lims)
                if kind == "neither reaches 50 %":
                    a.text(np.sqrt(lims[0] * lims[1]), -pos[k] + dy, kind, ha="center", va="center", fontsize=6.4, color=tc,
                           style="italic", bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
                out["ratios"][f"{k}/{arr}/depth{b}"] = dict(r["s50_ratio_squid_over_opm"], drawn_as=kind,
                                                            locations_favouring_opm=r["locations_favouring_opm"],
                                                            locations_favouring_squid=r["locations_favouring_squid"],
                                                            location_sign_flip_p=r["location_sign_flip_p"],
                                                            p_holm=p_holm[arr][k] if b == 0 else None)
            r0 = pairs[k, arr]["depth0"]
            ax_t.text(0.25, -pos[k] + dy, f"{r0['locations_favouring_opm']} : {r0['locations_favouring_squid']}", ha="center",
                      va="center", fontsize=6.8, color=tc)
            ax_t.text(0.75, -pos[k] + dy, fmt_p(p_holm[arr][k]), ha="center", va="center", fontsize=6.8, color=tc)
    ax_t.text(0.25, 1.0, "locations\nOPM : Neuromag", transform=ax_t.transAxes, ha="center", va="bottom", fontsize=6.6, color="0.3")
    ax_t.text(0.75, 1.0, "Holm p,\n9 anatomies", transform=ax_t.transAxes, ha="center", va="bottom", fontsize=6.6, color="0.3")
    ax_t.set_xlim(0, 1)
    ax_t.axis("off")
    row_axes(ax, keys, pos)
    ax.set_yticks([-pos[k] for k in keys], [label(k) + (" \u2020" if k == "adult" else "") for k in keys])
    for a in [ax_t, ax] + axs_b:
        a.set_ylim(-max(pos.values()) - 0.65, 0.95)
    for i, a in enumerate(axs_b):
        row_axes(a, keys, pos, labels=i == 0, short=True)
        a.set_ylim(-max(pos.values()) - 0.65, 0.65)
        a.set_title(f"{bands[i + 1]} mm", fontsize=8.5)
    axs_b[0].set_yticks([-pos[k] for k in keys], [label(k, short=True) + (" \u2020" if k == "adult" else "") for k in keys])
    for a, lims, ticks in [(ax, lims_a, (0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0))] + [(a, lims_b, (0.8, 1.0, 1.25, 1.5)) for a in axs_b]:
        a.set_xscale("log")
        a.set_xlim(*lims)
        a.axvline(1.0, color="0.5", lw=0.8, zorder=0)
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_minor_locator(NullLocator())
        a.set_xticklabels([f"{x:g}" for x in ticks])
    ax.text(0.995, 0.995, "favours OPM \u2192", transform=ax.transAxes, ha="right", va="top", fontsize=6.6, color=GREY)
    ax.text(0.005, 0.995, "\u2190 favours Neuromag", transform=ax.transAxes, ha="left", va="top", fontsize=6.6, color=GREY)
    ax.set_xlabel("S50 ratio, Neuromag / OPM (log scale; > 1: the OPM detects at a lower strength)")
    axs_b[1].set_xlabel("S50 ratio, Neuromag / OPM (log scale; > 1 favours the OPM)")
    ax.set_title(f"(a)  Spikes {bands[0]} mm below the scalp", loc="left")
    axs_b[0].annotate("(b)  Deeper bands", xy=(0, 1), xycoords="axes fraction", xytext=(0, 20), textcoords="offset points",
                      ha="left", va="bottom", fontsize=9)
    footnote(fig, f"Arrow: interval end open (bootstrap resamples outside the tested {min(strengths):g}\u2013{max(strengths):g} nAm); "
                  "dotted: open at both ends; open marker with dashed arrow: Neuromag does not reach 50 % detection, the ratio "
                  f"lies above the bound. \u2020 Adult at its measured head position, the other heads at top contact in the fixed "
                  f"adult helmet. {n_loc} "
                  f"locations ({n_ev} focal events) per band and anatomy.", y=0.045)
    h = [Line2D([], [], **mstyle("childC", color="0.2", ms=4.5)), Line2D([], [], **mstyle("childC", color=tint("0.2"), ms=4.5)),
         Line2D([], [], **mstyle("childC", False, color="0.2", ms=4.5))]
    fig.legend(h, ["dense OPM (dark)", "site-matched OPM (light)", "open: censored point estimate (a bound)"], loc="upper center",
               ncol=3, bbox_to_anchor=(0.6, 0.995), handlelength=1.0, columnspacing=1.4, fontsize=7,
               title="Practical detector, thresholds frozen at 1 false event/min; paired on identical simulated spikes",
               title_fontsize=7)
    bounds = {r["value_censored"] for r in out["ratios"].values() if r["drawn_as"] in ("lower bound", "upper bound")}
    if bounds - {"Neuromag does not reach 50 %"}:  # the footnote describes this case only
        raise ValueError(f"censored S50 ratios of another kind: {bounds}")
    style.save(fig, "Figure_R10_spikes")
    out["p_holm_depth0"] = {arr: {k: p_holm[arr][k] for k in keys} for arr in p_holm}
    return out


# ----------------------------------------------------------------------------------------------
def names(ks) -> str:
    return ", ".join(style.ANAT_SHORT[k] for k in ks) if ks else "none"


def provenance(s: dict, keys, v: dict) -> dict:
    """Inputs, description, alt text and caption draft per figure; every number and every
    qualitative statement in them is computed from the plotted values."""
    r5, r6, r7, r8, r9, r10 = (v[f"R{i}"] for i in range(5, 11))
    kids = [k for k in keys if k != "adult"]
    child = [k for k in kids if style.ANAT_CLASS[k] == "child"]
    rng = lambda xs: f"{db(min(xs))} to {db(max(xs))}"  # noqa: E731
    order_txt = ("Heads ordered by class (adult; scaled adult, size only; infant template; individual child), then head "
                 "circumference, largest first: " + ", ".join(f"{style.ANAT_SHORT[k]} {CM[k]:.1f} cm" for k in keys) + ".")
    cm_rng = f"{min(CM.values()):.1f}-{max(CM.values()):.1f} cm"

    # R5
    a5, b5 = r5["D_dB"]["top"], {r: c for r, c in r5["D_dB"]["counterfactual_x-centred"].items() if c}
    above_a = sum(a5[r][k] > a5[r]["adult"] for r in a5 for k in kids)
    below_b = sum(b5[r][k] <= b5[r]["adult"] for r in b5 for k in kids)
    # R6
    ib = {k: r6["values"][f"{k}/{COND}"] for k in kids}
    pj = {k: r6["values"][f"{k}/projected"] for k in kids}
    ib_pos = [k for k in kids if ib[k]["delta"] > 0]
    ib_excl0 = [k for k in kids if ib[k]["delta_ci95"][0] > 0]
    pj_incl0 = [k for k in kids if pj[k]["delta_ci95"][0] <= 0 <= pj[k]["delta_ci95"][1]]
    n_units = {k: (s["comparisons"][f"{k}/opm_dense/combined/{COND}/detect"]["delta"].get("n_parcels")
                   or s["comparisons"][f"{k}/opm_dense/combined/{COND}/detect"]["delta"].get("n")) for k in kids}
    vert = [n_units[k] for k in kids if style.ANAT_CLASS[k] == "scaled"]
    parc = [n_units[k] for k in kids if style.ANAT_CLASS[k] != "scaled"]
    # R7
    cfx = {k: r7[f"{k}/D_by_helmet_dB"]["counterfactual_x-centred"] for k in kids}
    cfx_above = [k for k in kids if cfx[k] > r7["adult_counterfactual_x_centred_D_dB"]]
    g_top = {k: r7[f"{k}/gap_mm"]["top"] for k in keys}
    g_cfx = {k: r7[f"{k}/gap_mm"]["counterfactual_x_centred"] for k in keys}
    wider_top = [k for k in kids if g_top[k] > g_top["adult"]]
    cfx_wider = [k for k in kids if g_cfx[k] > g_cfx["adult"]]
    sp = {(k, arr, c): r7[f"{k}/share_positive/{arr}/{c}"] for k in keys for arr in ("opm_dense", "opm_matched")
          for c in (COND, "projected")}
    m_higher = [k for k in kids if all(sp[k, "opm_matched", c] > sp["adult", "opm_matched", c] for c in (COND, "projected"))]
    dense_ib = [sp[k, "opm_dense", COND] for k in keys]
    loss = [r7[f"{k}/channel_count"]["adult_loss_dB"] for k in kids]
    sites = [r7[f"{k}/channel_count"]["n_sites"] for k in kids]
    a_sites = s["arrays"]["adult"]["opm_dense"]["n"]
    # R8
    st = r8["strata"]
    left = [k for k in kids if style.ANAT_CLASS[k] != "child"]
    left_pos = all(r["delta"] > 0 for k in left for r in st[k])
    shallow = [r["delta"] for k in child for r in st[k] if r["hi"] <= 30]
    below_2025 = [k for k in child for r in st[k] if (r["lo"], r["hi"]) == (20.0, 25.0) and r["ci95"][1] < 0]
    deep_pos = all(r["delta"] > 0 for k in child for r in st[k] if r["lo"] >= 30)
    # R9
    lv = [f"{a:g}" for a in r9["levels_fT_per_rtHz"]]
    falls = all(all(np.diff([r9[f"{k}/{arr}"][x] for x in lv]) < 0) for k in keys for arr in ("opm_dense", "opm_matched"))
    d30 = {arr: {k: r9[f"{k}/{arr}"][f"{FLOOR_FT:g}"] for k in keys} for arr in ("opm_dense", "opm_matched")}
    neg30 = {arr: [k for k in keys if d30[arr][k] < 0] for arr in d30}
    # R10
    rt = r10["ratios"]
    dense0 = {k: rt[f"{k}/opm_dense/depth0"] for k in keys}
    match0 = {k: rt[f"{k}/opm_matched/depth0"] for k in keys}
    if any(r["value"] is None for r in list(dense0.values()) + list(match0.values())):
        raise ValueError("a 10-20 mm S50 ratio is censored: the R10 texts assume point estimates there")
    incl1 = lambda r: (r["ci95"][0] is None or r["ci95"][0] <= 1) and (r["ci95"][1] is None or r["ci95"][1] >= 1)  # noqa: E731
    d_incl1, m_incl1 = [k for k in keys if incl1(dense0[k])], [k for k in keys if incl1(match0[k])]
    censored = {b: [k for k in keys for arr in ("opm_dense", "opm_matched")
                    if rt[f"{k}/{arr}/depth{b}"]["drawn_as"] in ("lower bound", "upper bound", "neither reaches 50 %")]
                for b in range(1, len(r10["bands_mm"]))}
    holm_max = max(r10["p_holm_depth0"]["opm_dense"].values())
    holm_m = {k: p for k, p in r10["p_holm_depth0"]["opm_matched"].items() if p < 0.05}
    n_boot = s["config"]["strata"]["n_boot"]
    n_ext = projection_terms()
    return {
        "Figure_R5_regions_heads": dict(
            inputs=[G3B] + [TARGETS.format(k) for k in keys],
            description=("Area-weighted median D (dB; dense OPM array vs Neuromag (306 channels), sensor plus brain noise, 10-nAm "
                         "cortical-normal dipoles; medial wall excluded) per lobe from g3b_summary.json placement_D["
                         "'<anatomy>/<placement>/combined/intrinsic+brain'].by_lobe, placement 'top' in (a) and "
                         "'counterfactual_x-centred' in (b). The parcel rows of (a) are computed here from "
                         "g3b_targets_<anatomy>.csv: D = 20 log10(detect_opm_dense_opm_intrinsic+brain / "
                         "detect_squid_top_combined_intrinsic+brain) per target, targets with region '*unknown' (medial wall) "
                         "left out, both hemispheres pooled, area-weighted median with weights area_mm2 (opmsquid.pediatric."
                         "weighted_median: values sorted, the first reaching half the cumulative area); the same rule on the same "
                         "files reproduces the stored lobe medians to within "
                         f"{r5['lobe_check_max_abs_diff_dB']:.1e} dB (CSV rounding). Mesial temporal = parahippocampal + "
                         "entorhinal. The target files hold no detectability for the laterally centred scaled helmet "
                         "(squid:counterfactual_x-centred), so (b) has lobes only. Diverging colour scale centred at 0 dB, "
                         f"symmetric (+/-{r5['colour_limit_dB']:g} dB), shared by both panels. " + order_txt),
            alt=("Two heat maps of D in decibels: rows are six lobes and four temporal or precentral parcel groups, columns the "
                 f"nine heads grouped by class and ordered by head circumference ({cm_rng}). Every cell is positive "
                 f"({db(r5['min_cell_dB'])} to {db(r5['max_cell_dB'])} dB), so the dense OPM array is ahead of Neuromag in "
                 f"every region of every head. In the fixed adult helmet {above_a} of {len(a5) * len(kids)} cells of the "
                 "smaller heads exceed the adult's value in the same row; in the helmet scaled with the head "
                 f"{below_b} of {len(b5) * len(kids)} lobe cells are at or below it. The parcel rows of the second map are "
                 "hatched as not stored."),
            caption_draft=("Figure R5. Which brain regions favour the OPM at which head size. Area-weighted median D, dense OPM "
                           "array minus Neuromag (306 channels), sensor plus brain noise, by lobe and for precentral, superior temporal, "
                           "parahippocampal and mesial temporal (parahippocampal + entorhinal) cortex, both hemispheres, medial "
                           "wall excluded; right: cortical targets per cell (range over the heads). (a) The fixed adult helmet "
                           "at top contact. (b) The helmet scaled with the head, laterally centred (lobes "
                           f"only: parcel values were not stored for this helmet). Head circumference {cm_rng}; the infant heads "
                           "are average templates, the children individual MRIs."),
            values=r5),
        "Figure_R6_pediatric_D": dict(
            inputs=[G3B, f"{G2_SUMMARY} :: retained_rank (the projection's term count)"],
            description=("g3b_summary.json comparisons['<anatomy>/opm_dense/combined/<condition>/detect'] for condition "
                         "'intrinsic+brain' (filled) and 'projected' (open): d_child.median with ci95 in (a), delta.median "
                         "with ci95 in (b). The adult's D is d_adult.median, identical in every comparison (checked): "
                         f"{db(r6['adult_D_dB'][COND])} dB (sensor plus brain noise), {db(r6['adult_D_dB']['projected'])} dB "
                         f"(after the {n_ext}-term projection). Intervals: parcel bootstrap, 95 %, {n_boot:,} resamples. Delta is vertex-wise for the "
                         f"scaled copies (the same cortical vertex; {min(vert):,}-{max(vert):,} homologous targets) and the "
                         "area-weighted median of parcel differences for the templates and children (no vertex correspondence; "
                         f"{min(parc)}-{max(parc)} Desikan-Killiany parcels). " + order_txt),
            alt=("Two dot plots with one row per smaller head: left, each head's D with vertical lines at the adult's D; "
                 "right, Delta with 95 % intervals, each class labelled with its estimator. With sensor plus brain noise Delta is "
                 "positive in "
                 f"{len(ib_pos)} of {len(kids)} heads, with the interval above zero in {len(ib_excl0)} "
                 f"(not: {names([k for k in kids if k not in ib_excl0])}); after the projection the interval includes "
                 f"zero for {names(pj_incl0)}."),
            caption_draft=("Figure R6. Smaller heads in the fixed adult helmet. (a) D_head, dense OPM array minus Neuromag (306 "
                           "channels), at top contact; vertical lines: the adult's D (solid: sensor plus brain noise; dashed: the "
                           f"room field added and removed by the {n_ext}-term projection). (b) Delta = D_head - D_adult, vertex-wise "
                           "for the scaled adults and parcel-matched for the templates and children (labelled in the panel). "
                           f"Filled: sensor plus brain noise; open: after the {n_ext}-term projection. Lines: 95 % parcel-bootstrap "
                           "intervals (spatial, within one anatomy; no between-subject "
                           f"variability). With sensor plus brain noise Delta is {rng([ib[k]['delta'] for k in kids])} dB over "
                           f"the eight heads and {rng([ib[k]['delta'] for k in child])} dB for the individual children."),
            values=r6),
        "Figure_R7_helmet_fit": dict(
            inputs=[G3B, G3A, JAS, "docs/methods.md", f"{G2_SUMMARY} :: retained_rank (the projection's term count)"],
            description=("(a) g3b_summary.json placement_D['<anatomy>/<placement>/combined/intrinsic+brain'].median for "
                         "centred, top, x-centred, counterfactual and counterfactual_x-centred; dashed: the adult at top "
                         f"({db(r7['adult_top_D_dB'])} dB) and counterfactual_x-centred ({db(r7['adult_counterfactual_x_centred_D_dB'])}"
                         f" dB); grey band: the adult over its {r7['adult_family_D_dB']['n']} source-blind placements in the "
                         "fixed helmet (top, back and the ten +-5 mm x/y, pitch +-10, roll +-5 and yaw +-10 deg variants; "
                         f"docs/methods.md), {db(r7['adult_family_D_dB']['min'])} to {db(r7['adult_family_D_dB']['max'])} dB. "
                         "(b) sensor_distances['<anatomy>']['squid:top' | 'squid:counterfactual_x-centred'].median_mm (median "
                         "magnetometer coil-to-scalp distance); dotted: the sphere model's "
                         f"{r7['sphere_gap_mm']:g} mm (g3a_size_benchmark.json size_following[*].xi_squid_mm; Jas et al. 2026, "
                         "docs/literature/jas2026.md section 2.1). (c) comparisons['<anatomy>/<array>/combined/<condition>/"
                         "detect'].d_child.share_positive (adult: d_adult.share_positive, identical in every comparison), matched "
                         "array solid, dense array faint; filled sensor plus brain noise, open after the projection. (d) channel_count_control["
                         "'<anatomy>/combined'].n_sites (numbers) and .d_adult_subsampled.median (the adult's dense array "
                         "subsampled by farthest-point sampling to that count; same placement and noise) against the adult's full "
                         "array (dashed). " + order_txt),
            alt=("Four panels. (a) D of nine heads at five helmet conditions: in the fixed helmet the smaller heads' D is "
                 f"{rng([r7[f'{k}/D_by_helmet_dB']['centred'] for k in kids])} dB at the adult's pose and "
                 f"{rng([r7[f'{k}/D_by_helmet_dB']['top'] for k in kids])} dB at top contact; in the scaled, laterally centred "
                 f"helmet {rng(list(cfx.values()))} dB, above the adult's {db(r7['adult_counterfactual_x_centred_D_dB'])} dB "
                 f"for {names(cfx_above)} only. (b) The Neuromag gap at top contact is wider than the adult's for "
                 f"{len(wider_top)} of {len(kids)} smaller heads; in the scaled helmet it is narrower than the adult's except "
                 f"for {names(cfx_wider)}; every gap exceeds the sphere model's {r7['sphere_gap_mm']:g} mm. (c) The dense "
                 f"array's share of cortex with the OPM ahead is {min(dense_ib):.2f}-{max(dense_ib):.2f} with sensor plus "
                 f"brain noise; the site-matched array's is higher than the adult's in {len(m_higher)} of {len(kids)} smaller heads "
                 "under both conditions. (d) The adult's D falls by "
                 f"{min(loss):.2f}-{max(loss):.2f} dB when its {a_sites}-site array is cut to {min(sites)}-{max(sites)} sites."),
            caption_draft=("Figure R7. Helmet fit and head size, each panel against one fixed reference. (a) D (dense OPM array minus "
                           "Neuromag (306 channels), sensor plus brain noise) at five helmet conditions; dashed: the adult at top contact and "
                           "in its own scaled, laterally centred helmet; grey: the adult over its "
                           f"{r7['adult_family_D_dB']['n']} standard placements. (b) Median Neuromag coil-to-scalp gap: "
                           f"{min(g_top[k] for k in kids):.1f}-{max(g_top[k] for k in kids):.1f} mm at top contact (adult "
                           f"{g_top['adult']:.1f}), {min(g_cfx[k] for k in kids):.1f}-{max(g_cfx[k] for k in kids):.1f} mm in the "
                           f"scaled, laterally centred helmet (adult {g_cfx['adult']:.1f}); dotted: the "
                           f"{r7['sphere_gap_mm']:g}-mm gap the sphere model keeps for every head size. (c) Share of cortical area "
                           "with D > 0 (the preprint's metric), site-matched array (solid) and dense array (faint). (d) Dense OPM "
                           f"sites per head ({min(sites)}-{max(sites)}; adult {a_sites}) and the adult's D with its array "
                           "subsampled to each count. The helmet scaled with the head is an array construction (coil centres scaled), "
                           "not a device."),
            values=r7),
        "Figure_R8_depth_matched": dict(
            inputs=[G3B],
            description=("g3b_summary.json comparisons['<anatomy>/opm_dense/combined/intrinsic+brain/detect'].delta_by_depth: "
                         "per depth stratum (lo-hi below the scalp, each head's own depth), the area-weighted median D of the "
                         "head's cortical targets minus the adult's (delta), with a 95 % interval from resampling parcels in "
                         "each anatomy (ci95; drawn for the individual children only). Strata with fewer than "
                         f"{r8['min_n']} targets in either head have no delta and are skipped (0-10 mm: the adult has none). "
                         f"Targets per stratum: {r8['n_child_range'][0]:,}-{r8['n_child_range'][1]:,} (smaller heads), "
                         f"{r8['n_adult_range'][0]:,}-{r8['n_adult_range'][1]:,} (adult). Points at stratum centres (60-90 mm at "
                         "75 mm); the children's points are offset by 0.7 mm to keep the intervals apart."),
            alt=("Two line plots of Delta against depth below the scalp. Scaled adults and infant templates: Delta is "
                 f"{'positive in every stratum' if left_pos else 'not positive in every stratum'}. Individual children: Delta "
                 f"is {rng(shallow)} dB in the strata shallower than 30 mm, with the interval entirely below zero at 20-25 mm "
                 f"for {names(below_2025)}, and {'positive' if deep_pos else 'mixed'} from 30 mm deeper."),
            caption_draft=("Figure R8. Depth-matched comparison with the adult. Delta per depth stratum (dense OPM array vs Neuromag "
                           "(306 channels), sensor plus brain noise, top contact): the area-weighted median D of a head's cortical targets in the "
                           "stratum minus the adult's in the same stratum (a difference of the two heads' medians). (a) Scaled "
                           "adults and infant templates. (b) Individual children with 95 % parcel-bootstrap intervals. Strata with "
                           "fewer than "
                           f"{r8['min_n']} targets in either head are omitted; targets per stratum "
                           f"{r8['n_child_range'][0]:,}-{r8['n_child_range'][1]:,} (smaller heads), "
                           f"{r8['n_adult_range'][0]:,}-{r8['n_adult_range'][1]:,} (adult)."),
            values=r8),
        "Figure_R9_noise_floor": dict(
            inputs=[G3B, G2_CFG, JAS],
            description=("g3b_summary.json sensitivity_median_D_dB['<anatomy>/opm_asd_<n>fT/<array>/combined/intrinsic+brain'] "
                         f"for the OPM white-noise levels of configs/g2_adult.toml ({', '.join(lv)} fT/sqrt(Hz)); the "
                         f"{r9['primary_fT_per_rtHz']:g}-fT points are the stored sensitivity keys, checked equal to the "
                         "primary D_median_dB['<anatomy>/<array>/combined/intrinsic+brain/detect']. Neuromag noise unchanged. "
                         "Vertical lines: the primary level and about 30 fT/sqrt(Hz), the approximate OPM empty-room floor of the "
                         "Jas et al. (2026) recording (docs/literature/jas2026.md section 5.6: about 28-35 fT/sqrt(Hz), read from "
                         "the preprint's Fig. 9; a broadband value, not matched to the 1-40 Hz band). " + order_txt),
            alt=("Two line plots of D against the OPM white noise, one line per head; every line "
                 f"{'falls' if falls else 'does not fall monotonically'} as the noise rises. At "
                 f"{FLOOR_FT:g} fT per root hertz the dense array's D is {rng(list(d30['opm_dense'].values()))} dB "
                 f"(negative for {names(neg30['opm_dense'])}); the site-matched array's is "
                 f"{rng(list(d30['opm_matched'].values()))} dB (negative for {len(neg30['opm_matched'])} of {len(keys)} heads)."),
            caption_draft=("Figure R9. The noise condition. D (OPM minus Neuromag (306 channels), sensor plus brain noise, top "
                           "contact) against the OPM's white sensor noise, Neuromag unchanged. (a) Dense OPM array; (b) "
                           "site-matched OPM array. "
                           f"Grey line: the primary {r9['primary_fT_per_rtHz']:g} fT/sqrt(Hz); dotted: about {FLOOR_FT:g} "
                           "fT/sqrt(Hz), the approximate OPM empty-room floor in the Jas et al. (2026) recording (broadband). At "
                           f"{FLOOR_FT:g} fT/sqrt(Hz) the dense array's D is {rng(list(d30['opm_dense'].values()))} dB (adult "
                           f"{db(d30['opm_dense']['adult'])})."),
            values=r9),
        "Figure_R10_spikes": dict(
            inputs=[G3B, G4_REPORT] + [G4.format(k) for k in keys],
            description=("g4_<anatomy>_summary.json paired['<array>/opm_vs_squid/combined/practical@1']['depth<b>']"
                         ".s50_ratio_squid_over_opm (value, value_bounds, value_censored, ci95) for opm_dense (dark) and "
                         "opm_matched (light); depth bands from the table header of G4_pediatric_report.md "
                         f"({', '.join(r10['bands_mm'])} mm; every stored location checked to lie in its band). Annotations at "
                         "10-20 mm: locations_favouring_opm : locations_favouring_squid, and location_sign_flip_p Holm-adjusted "
                         "(step-down) over the nine anatomies, computed here for each array separately; the stored p values are "
                         "uncorrected, and the 10-20 mm band was taken as the single endpoint after seeing the data (post hoc). "
                         "An open interval end (bootstrap resamples outside the tested "
                         f"{min(r10['strengths_nAm']):g}-{max(r10['strengths_nAm']):g} nAm) is an arrow to the axis edge, an "
                         "interval open at both ends a dotted line; a censored point estimate (value_bounds) an open marker with a "
                         "dashed arrow; 'neither reaches 50 %': no estimate. The deeper panel shows all three deeper bands "
                         "(depth1-depth3). Adult (dagger) at its measured head position (G2), the other heads at top contact. 18 "
                         "locations per band and anatomy. " + order_txt),
            alt=("Forest plots of the S50 ratio on a log axis with a line at 1. At 10-20 mm the dense array's ratio is "
                 f"{min(r['value'] for r in dense0.values()):.2f}-{max(r['value'] for r in dense0.values()):.2f} in the nine "
                 f"heads, its interval including 1 for {names(d_incl1)}; the site-matched array's ratio is "
                 f"{min(r['value'] for r in match0.values()):.2f}-{max(r['value'] for r in match0.values()):.2f}, its interval "
                 f"including 1 for {len(m_incl1)} heads. In the deeper bands the ratios lie near 1 with wide or open "
                 f"intervals; censored or missing estimates: {len(censored[1])}, {len(censored[2])} and {len(censored[3])} of "
                 f"{2 * len(keys)} at {', '.join(r10['bands_mm'][1:])} mm."),
            caption_draft=("Figure R10. Interictal spikes. Strength for 50 % detection, Neuromag (306 channels) over OPM, paired on identical "
                           "simulated spikes (practical detector, thresholds frozen at 1 false event per minute); > 1: the OPM "
                           "detects at a lower strength. Dark: dense OPM; light: site-matched OPM. Lines: 95 % bootstrap intervals "
                           "over the 18 locations per band; arrows: open ends; open markers: censored estimates. (a) 10-20 mm, "
                           f"dense-array ratios {min(r['value'] for r in dense0.values()):.2f}-"
                           f"{max(r['value'] for r in dense0.values()):.2f}; right: locations favouring each system and the "
                           "location sign-flip p, Holm-adjusted over the nine anatomies (each array separately; dense array "
                           f"largest {fmt_p(holm_max)}; site-matched array below 0.05 for {names(list(holm_m))}); the endpoint was chosen "
                           "post hoc. (b) Deeper bands. Dagger: the adult sits at its measured head position, the smaller heads at "
                           "top contact in the fixed adult helmet."),
            values=r10),
    }


def main():
    style.apply()
    s = load(G3B)
    keys = order(s)
    CM.update({k: s["anatomies"][k]["head_size"]["ofc_mm"] / 10 for k in keys})
    reg, n, check = region_values(s, keys)
    if check > 0.005:
        raise SystemExit(f"target files do not reproduce the stored lobe medians (max difference {check:.4f} dB)")
    v = {"R5": fig_r5(s, keys, reg, n, check), "R6": fig_r6(s, keys), "R7": fig_r7(s, keys), "R8": fig_r8(s, keys),
         "R9": fig_r9(s, keys), "R10": fig_r10(s, keys)}
    style.write_provenance(OUT_JSON, provenance(s, keys, v))
    print(f"order: {', '.join(f'{k} ({CM[k]:.1f} cm)' for k in keys)}")
    print(f"R5 lobe check (target files vs stored lobe medians): max |diff| {check:.2e} dB")
    print("R10 Holm (dense, 10-20 mm): " + ", ".join(f"{k} {fmt_p(p)}" for k, p in v["R10"]["p_holm_depth0"]["opm_dense"].items()))
    print(f"wrote {len(v)} figures and results/report/{OUT_JSON}")


if __name__ == "__main__":
    main()
