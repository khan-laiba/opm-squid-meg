#!/usr/bin/env python3
"""Report figure R16: the confirmatory spike run in the manuscript's terms, drawn from stored outputs only (nothing is
simulated or re-analysed).

  R16 (a) The declared endpoint, per anatomy: the paired S50 ratio Neuromag (306 channels) / dense OPM for focal spikes
      in the declared depth band, practical scanning detector with thresholds frozen at the declared false-event rate,
      noise replicate 0, with its 95 % location-bootstrap interval (dark); beside it, the exploratory run's ratio for the
      same comparison (light; results/g4, endpoint chosen after the analyses). Right: locations favouring each array and
      the location sign-flip p, Holm-adjusted over the anatomies (a tick: below the declared alpha, the declared test passed).
      (b) The oracle detector (known topography, waveform and time) and (c) the declared detector-mismatch variant
      (thresholds frozen at the same rate), dense OPM, replicate 0, each with its own Holm-adjusted p (secondary analyses,
      not confirmatory); a grey bar marks the endpoint's ratio of panel (a) for reference.

Inputs: <results-dir>/g4_confirm_summary.json (the combine step of scripts/g4_confirmatory.py: per anatomy the endpoint,
the oracle and mismatch comparisons, Holm-adjusted p per family) and <results-dir>/g4c_<anatomy>_summary.json (design
numbers for the footnote), results/g4/g4_pediatric_comparison.json (the exploratory ratios, checked equal to the copies in
the combined summary). The exploratory Holm p is the combined summary's, over the anatomies present. Censoring is drawn
as in Figure R10 (scripts/report_figures_pediatric.py): an open interval end (resamples outside the tested strengths) is
an arrow to the axis edge, an interval open at both ends a dotted line, a censored point estimate (only bounded) an open
marker at its bound with a dashed arrow, and 'neither reaches 50 %' is written in the row. p-values are printed as the
report's facts print them (scripts/report_facts_g4.py pval). A run that is not confirmatory (incomplete, test settings,
several or uncommitted commits, several configurations: the combine step's rule) is drawn with a banner saying so.

Outputs: <out-dir>/Figure_R16_confirm.png and <out-dir>/figures_confirm.json (inputs, description, alt text, caption
draft and the plotted values, with the commit that drew them).

Usage: PYTHONPATH=src .venv/bin/python scripts/report_figures_confirm.py [--results-dir results/g4_confirm]
       [--out-dir results/report]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

NAME = "Figure_R16_confirm"
OUT_JSON = "figures_confirm.json"
EXPLO = "results/g4/g4_pediatric_comparison.json"
CLASSES = ("adult", "scaled", "template", "child")
GREY = "0.45"
PASS = "#009E73"  # Okabe-Ito bluish green: the tick of a passed test


def _facts_g4():
    spec = importlib.util.spec_from_file_location("report_facts_g4", ROOT / "scripts" / "report_facts_g4.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FG4 = _facts_g4()  # the report's number formats (p-values, ratios), so the figure prints what the text prints


def rel(p: Path) -> str:
    p = Path(p).resolve()
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)


# ----------------------------------------------------------------------------------------------
# drawing helpers (as scripts/report_figures_pediatric.py, Figure R10)
def slots(keys, gap: float) -> dict:
    """Axis positions, one unit per anatomy, with ``gap`` added between classes."""
    pos, x, prev = {}, 0.0, None
    for k in keys:
        if prev is not None and style.ANAT_CLASS[k] != prev:
            x += gap
        pos[k], x, prev = x, x + 1.0, style.ANAT_CLASS[k]
    return pos


def class_breaks(keys, pos) -> list[float]:
    """Positions half-way between two runs of different classes."""
    return [-0.5 * (pos[a] + pos[b]) for a, b in zip(keys[:-1], keys[1:]) if style.ANAT_CLASS[a] != style.ANAT_CLASS[b]]


def mstyle(k: str, filled: bool = True, ms: float = 5.0, color=None) -> dict:
    """Marker of an anatomy: class shape, anatomy colour; open = white face."""
    m = style.CLASS_MARKER[style.ANAT_CLASS[k]]
    c = style.ANAT_COLOR[k] if color is None else color
    size = ms * {"D": 0.85, "^": 1.1}.get(m, 1.0)
    return dict(marker=m, ms=size, mec=c, mfc=c if filled else "white", mew=0.9, ls="none")


def tint(c, f: float = 0.55):
    """Colour blended towards white by ``f``."""
    r = np.array(to_rgb(c))
    return tuple(r + (1 - r) * f)


def arrow(ax, x0, x1, y, color, ls="-"):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.0, ls=ls, mutation_scale=7,
                                                                shrinkA=0, shrinkB=0))


def draw_ratio(ax, y, sr: dict, color, marker: dict, lims, lw: float = 1.0) -> str:
    """One S50 ratio with its 95 % interval; returns how it was drawn (Figure R10's rules)."""
    lo, hi = sr["ci95"] if sr.get("ci95") else (None, None)
    v = sr.get("value")
    if v is not None:
        if lo is None and hi is None:
            ax.plot(lims, [y, y], color=color, lw=0.8 * lw, ls=(0, (1, 1.5)))
            arrow(ax, lims[0] * 1.06, lims[0], y, color)
            arrow(ax, lims[1] / 1.06, lims[1], y, color)
        else:
            a, b = lo if lo is not None else lims[0], hi if hi is not None else lims[1]
            ax.plot([a, b], [y, y], color=color, lw=lw, solid_capstyle="butt")
            if lo is None:
                arrow(ax, a * 1.06, lims[0], y, color)
            if hi is None:
                arrow(ax, b / 1.06, lims[1], y, color)
        ax.plot(v, y, **marker)
        return "value" if lo is not None and hi is not None else "value, interval open"
    b_lo, b_hi = sr.get("value_bounds") or (None, None)
    if b_lo is not None or b_hi is not None:
        start, edge = (b_lo, lims[1]) if b_lo is not None else (b_hi, lims[0])
        arrow(ax, start, edge, y, color, ls=(0, (2, 1.5)))
        ax.plot(start, y, **dict(marker, mfc="white"))
        return "lower bound" if b_lo is not None else "upper bound"
    ax.text(np.sqrt(lims[0] * lims[1]), y, "neither reaches 50 %", ha="center", va="center", fontsize=6.4, color=color,
            style="italic", bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
    return "neither reaches 50 %"


def limits(srs) -> tuple[float, float]:
    """Shared x limits (log axis): the drawn values and finite interval ends with a margin, at least 0.8-1.6, from a fixed
    set of round limits."""
    vals = [x for sr in srs for x in [sr.get("value"), *(sr.get("ci95") or []), *(sr.get("value_bounds") or [])]
            if x is not None and np.isfinite(x) and x > 0]
    lo, hi = min(vals + [0.8]), max(vals + [1.6])
    lo = max([x for x in (0.25, 0.33, 0.5, 0.6, 0.7, 0.8) if x <= lo / 1.04] or [0.25])
    hi = min([x for x in (1.6, 2.0, 2.5, 3.0, 4.0, 5.0) if x >= hi * 1.04] or [5.0])
    return lo, hi


def ticks_in(lims, narrow: bool = False) -> list[float]:
    full = (0.25, 0.33, 0.5, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
    return [t for t in ((0.25, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0) if narrow else full) if lims[0] <= t <= lims[1]]


def log_axis(ax, lims, narrow: bool = False):
    ax.set_xscale("log")
    ax.set_xlim(*lims)
    ax.axvline(1.0, color="0.5", lw=0.8, zorder=0)
    t = ticks_in(lims, narrow)
    ax.xaxis.set_major_locator(FixedLocator(t))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{x:g}" for x in t])


def holm_text(p: float, alpha: float) -> tuple[str, bool]:
    return FG4.pval(p) + (" ✓" if p < alpha else ""), p < alpha


def common(values: list, f=lambda x: f"{x:g}") -> str:
    """One value when all agree, else 'lo–hi' (a test run mixing settings)."""
    u = sorted(set(values))
    return f(u[0]) if len(u) == 1 else f"{f(u[0])}–{f(u[-1])}"


# ----------------------------------------------------------------------------------------------
def not_confirmatory_reasons(cf: dict, per: dict) -> list[str]:
    """The combine step's rule (scripts/g4_confirmatory.py combine), in words."""
    reasons = []
    if not cf["complete"]:
        reasons.append(f"{len(cf['anatomies'])} of {len(style.ANAT_ORDER)} anatomies")
    tests = [k for k in per if not per[k]["confirmatory"]]
    if tests:
        reasons.append("test settings")
    if len(cf["simulated_at_commits"]) > 1:
        reasons.append("several commits")
    if any("+dirty" in c for c in cf["simulated_at_commits"]):
        reasons.append("uncommitted changes")
    if len(cf["config_digests"]) > 1:
        reasons.append("several configurations")
    if bool(cf["confirmatory"]) != (not reasons):
        raise ValueError(f"combined summary: confirmatory = {cf['confirmatory']} disagrees with its own flags ({reasons})")
    return reasons


def figure(results_dir: Path, out_dir: Path) -> dict:
    cf_path = results_dir / "g4_confirm_summary.json"
    if not cf_path.is_file():
        raise FileNotFoundError(f"{cf_path} not found: the combine step of scripts/g4_confirmatory.py has not run")
    cf = json.loads(cf_path.read_text())
    labs = [k for k in style.ANAT_ORDER if k in cf["anatomy"]]
    per = {k: json.loads((results_dir / f"g4c_{k}_summary.json").read_text()) for k in labs}
    explo = json.loads((ROOT / EXPLO).read_text())["comparison"]
    ep = cf["endpoint"]["definition"]
    alpha, rate = ep["alpha"], ep["false_events_per_min"]
    rt = f"{rate:g}"
    comp = cf["endpoint"]["comparison"]  # '<dense>_vs_<comparator>/practical@<rate>/replicate0'
    pair = comp.split("/practical@")[0]
    if pair != "opm_dense/opm_vs_squid/combined" or not comp.endswith("/replicate0"):
        raise ValueError(f"unexpected endpoint comparison {comp!r}")
    fam = {"endpoint": cf["families"]["endpoint"]["holm_p"],
           "oracle": cf["families"][f"secondary/{pair}/oracle/replicate0"]["holm_p"],
           "mismatch": cf["families"][f"secondary/{pair}/mismatch@{rt}/replicate0"]["holm_p"]}
    reasons = not_confirmatory_reasons(cf, per)

    rows = {}
    for k in labs:
        a = cf["anatomy"][k]
        e = explo[f"{k}/paired/{pair}/practical@{rt}/depth0"]
        for key in ("s50_ratio_squid_over_opm", "location_sign_flip_p", "locations_favouring_opm", "locations_favouring_squid"):
            if a["exploratory"][key] != e[key]:
                raise ValueError(f"{k}: the combined summary's exploratory {key} differs from {EXPLO}")
        if a["endpoint"]["holm_p"] != fam["endpoint"][k]:
            raise ValueError(f"{k}: anatomy endpoint holm_p differs from families['endpoint']")
        rows[k] = dict(confirmatory=dict(a["endpoint"]["s50_ratio_squid_over_opm"], locations_favouring_opm=a["endpoint"]["locations_favouring_opm"],
                                         locations_favouring_squid=a["endpoint"]["locations_favouring_squid"],
                                         location_sign_flip_p=a["endpoint"]["location_sign_flip_p"], p_holm=fam["endpoint"][k]),
                       exploratory=dict(e["s50_ratio_squid_over_opm"], locations_favouring_opm=e["locations_favouring_opm"],
                                        locations_favouring_squid=e["locations_favouring_squid"],
                                        location_sign_flip_p=e["location_sign_flip_p"],
                                        p_holm_over_present_anatomies=a["exploratory"]["holm_p_over_present_anatomies"]),
                       oracle=dict(a["oracle"]["s50_ratio_squid_over_opm"], locations_favouring_opm=a["oracle"]["locations_favouring_opm"],
                                   locations_favouring_squid=a["oracle"]["locations_favouring_squid"],
                                   location_sign_flip_p=a["oracle"]["location_sign_flip_p"], p_holm=fam["oracle"][k]),
                       mismatch=dict(a["mismatch"]["s50_ratio_squid_over_opm"], locations_favouring_opm=a["mismatch"]["locations_favouring_opm"],
                                     locations_favouring_squid=a["mismatch"]["locations_favouring_squid"],
                                     location_sign_flip_p=a["mismatch"]["location_sign_flip_p"], p_holm=fam["mismatch"][k]))

    # design numbers for the footnote (every anatomy's summary; a test run may mix them)
    n_loc = common([per[k]["n_locations"] for k in labs])
    n_ev = common([per[k]["n_events_per_replicate"] for k in labs])
    strengths = per[labs[0]]["config"]["inherited"]["events"]["strengths_nAm"]
    n_boot = common([per[k]["config"]["confirmatory"]["design"]["bootstrap_resamples"] for k in labs], lambda x: f"{x:,}")
    tests = {per[k]["declared_choices"].get("sign_flip_test", "") for k in labs} - {""}
    method = " / ".join(sorted({t.split(" (")[0] for t in tests})) or "exact"
    o_alpha = common([o["alpha"] for k in labs for o in per[k]["oracle"].values()])
    n_loc_x = common([explo[f"{k}/n_locations_per_depth_band"]["depth0"] for k in labs])
    band = f"{ep['depth_band_mm'][0]:g}–{ep['depth_band_mm'][1]:g}"
    var = {x["name"]: x for x in per[labs[0]]["config"]["confirmatory"].get("variant", [])}["mismatch"]
    mismatch = (f"templates between the simulated morphologies and candidate fields from a "
                f"{len(var['dictionary_conductivity'])}-layer model with a {var['coreg_shift_mm']:g}-mm / "
                f"{var['coreg_angle_deg']:g}-deg coregistration error")

    # layout
    pos = slots(labs, gap=0.45)
    span = max(pos.values()) + 1.25  # rows in axis units
    h_a, h_b = 0.46 * span, 0.32 * span  # inches
    top_in = 1.0 + (0.2 if reasons else 0.0)  # banner, legend, panel title and column headers
    gap_in, bottom_in = 0.9, 1.3  # (a)'s axis label and (b, c)'s titles; (b, c)'s axis labels and the footnote
    H = top_in + h_a + gap_in + h_b + bottom_in
    fig = plt.figure(figsize=(style.FULL_W, H))
    gs = fig.add_gridspec(2, 1, height_ratios=[h_a, h_b], hspace=gap_in / ((h_a + h_b) / 2), left=0.255, right=0.995,
                          top=1 - top_in / H, bottom=bottom_in / H)
    ga = gs[0].subgridspec(1, 2, width_ratios=[2.15, 1.0], wspace=0.03)
    ax_a, ax_ta = fig.add_subplot(ga[0]), fig.add_subplot(ga[1])
    gb = gs[1].subgridspec(1, 4, width_ratios=[1.3, 0.34, 1.3, 0.34], wspace=0.04)
    ax_b, ax_tb, ax_c, ax_tc = (fig.add_subplot(gb[i]) for i in range(4))
    lims = limits([r[s] for r in rows.values() for s in ("confirmatory", "exploratory", "oracle", "mismatch")])
    drawn = {}
    ylim = (-max(pos.values()) - 0.62, 0.62)

    # (a) the declared endpoint, confirmatory (dark) and exploratory (light)
    for k in labs:
        y = -pos[k]
        c = style.ANAT_COLOR[k]
        drawn[f"{k}/confirmatory"] = draw_ratio(ax_a, y + 0.17, rows[k]["confirmatory"], c, mstyle(k, ms=5.0), lims)
        drawn[f"{k}/exploratory"] = draw_ratio(ax_a, y - 0.17, rows[k]["exploratory"], tint(c), mstyle(k, ms=4.2, color=tint(c)),
                                               lims, lw=0.9)
        r, x = rows[k]["confirmatory"], rows[k]["exploratory"]
        ax_ta.text(0.3, y + 0.17, f"{r['locations_favouring_opm']} : {r['locations_favouring_squid']}", ha="center", va="center",
                   fontsize=6.8, color="0.1")
        ax_ta.text(0.3, y - 0.17, f"{x['locations_favouring_opm']} : {x['locations_favouring_squid']}", ha="center",
                   va="center", fontsize=6.8, color="0.55")
        txt, ok = holm_text(r["p_holm"], alpha)
        ax_ta.text(0.8, y + 0.17, txt, ha="center", va="center", fontsize=6.8, color="0.1", weight="bold" if ok else "normal")
        ax_ta.text(0.8, y - 0.17, FG4.pval(x["p_holm_over_present_anatomies"]), ha="center", va="center", fontsize=6.8,
                   color="0.55")
    ax_ta.text(0.3, 1.0, "locations\nOPM : Neuromag", transform=ax_ta.transAxes, ha="center", va="bottom", fontsize=6.6, color="0.3")
    ax_ta.text(0.8, 1.0, f"Holm p\n({len(labs)} anatomies)", transform=ax_ta.transAxes, ha="center", va="bottom", fontsize=6.6,
               color="0.3")
    # (b) oracle and (c) mismatch, with the endpoint's ratio as a grey bar
    for ax, ax_t, key in ((ax_b, ax_tb, "oracle"), (ax_c, ax_tc, "mismatch")):
        for k in labs:
            y = -pos[k]
            ref = rows[k]["confirmatory"].get("value")
            if ref is not None:
                ax.plot(ref, y, marker="|", ms=8, mew=1.3, color="0.62", ls="none", zorder=1)
            drawn[f"{k}/{key}"] = draw_ratio(ax, y, rows[k][key], style.ANAT_COLOR[k], mstyle(k, ms=4.6), lims)
            txt, ok = holm_text(rows[k][key]["p_holm"], alpha)
            ax_t.text(0.5, y, txt, ha="center", va="center", fontsize=6.6, color="0.1", weight="bold" if ok else "normal")
        ax_t.text(0.5, 1.0, "Holm p", transform=ax_t.transAxes, ha="center", va="bottom", fontsize=6.6, color="0.3")
    for ax in (ax_a, ax_b, ax_c):
        log_axis(ax, lims, narrow=ax is not ax_a)
        ax.set_ylim(*ylim)
        ax.tick_params(axis="y", length=0)
        for b in class_breaks(labs, pos):
            ax.axhline(b, color="0.88", lw=0.7, zorder=0)
    for ax in (ax_ta, ax_tb, ax_tc):
        ax.set_xlim(0, 1)
        ax.set_ylim(*ylim)
        ax.axis("off")
    ax_a.set_yticks([-pos[k] for k in labs], [style.ANAT_LABEL[k] for k in labs])
    ax_b.set_yticks([-pos[k] for k in labs], [style.ANAT_LABEL[k] for k in labs])
    ax_c.set_yticks([-pos[k] for k in labs], [""] * len(labs))
    ax_a.text(0.995, 0.995, "favours OPM →", transform=ax_a.transAxes, ha="right", va="top", fontsize=6.6, color=GREY)
    ax_a.text(0.005, 0.995, "← favours Neuromag", transform=ax_a.transAxes, ha="left", va="top", fontsize=6.6, color=GREY)
    ax_a.set_xlabel("S50 ratio, Neuromag / OPM (log scale; > 1: the OPM detects at a lower strength)")
    ax_b.set_xlabel("S50 ratio, Neuromag / OPM")
    ax_c.set_xlabel("S50 ratio, Neuromag / OPM")
    ax_a.set_title(f"(a)  Declared endpoint, spikes {band} mm deep", loc="left")
    ax_b.set_title("(b)  Oracle detector", loc="left")
    ax_c.set_title("(c)  Mismatched detector", loc="left")
    h = [Line2D([], [], **mstyle("childC", color="0.15", ms=4.8)), Line2D([], [], **mstyle("childC", color=tint("0.15"), ms=4.2)),
         Line2D([], [], **mstyle("childC", False, color="0.15", ms=4.8)),
         Line2D([], [], marker="|", ms=8, mew=1.3, color="0.62", ls="none")]
    fig.legend(h, ["confirmatory run (dark)", "exploratory run (light; endpoint chosen after the analyses)",
                   "open: censored estimate (a bound)", "(b, c): endpoint ratio of (a)"],
               loc="upper center", ncol=2, bbox_to_anchor=(0.6, 1 - (0.04 + (0.2 if reasons else 0.0)) / H), handlelength=1.0,
               columnspacing=1.4, fontsize=7,
               title=f"Dense OPM vs Neuromag 306, thresholds frozen at {rt} false event/min; paired on identical simulated spikes",
               title_fontsize=7)
    if reasons:
        fig.text(0.5, 1 - 0.03 / H, "NOT CONFIRMATORY: " + "; ".join(reasons), ha="center", va="top", fontsize=8.5,
                 color="#D55E00", weight="bold")
    fig.text(0.0, 0.72 / H, textwrap.fill(
        f"Lines: 95 % location-bootstrap intervals ({n_boot} resamples) over the {n_loc} locations per anatomy ({n_ev} focal events, "
        f"noise replicate 0); arrows: open ends (resamples outside the tested {strengths[0]:g}–{strengths[-1]:g} nAm); open "
        "marker with dashed arrow: censored point estimate. p: two-sided sign-flip test on the per-location differences in "
        f"detection counts ({method}), Holm-adjusted over the anatomies within each panel's family; ✓ and bold: below "
        f"{alpha:g}. (b) known topography, waveform and time, per-trial false-positive probability {o_alpha}; (c) {mismatch}, "
        f"thresholds recalibrated to {rt} per minute. Exploratory run: {n_loc_x} locations per anatomy. The adult at its measured head "
        "position, the other heads at top contact.", 150), ha="left", va="top", fontsize=6.5, color="0.3")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{NAME}.png"
    fig.savefig(path, dpi=style.DPI)
    plt.close(fig)
    kinds = {v for v in drawn.values()}
    return dict(labs=labs, rows=rows, drawn=drawn, kinds=sorted(kinds), lims=list(lims), reasons=reasons, alpha=alpha, rate=rate,
                band=band, n_loc=n_loc, n_ev=n_ev, n_boot=n_boot, method=method, oracle_alpha=o_alpha, n_loc_exploratory=n_loc_x,
                mismatch=mismatch,
                strengths=[strengths[0], strengths[-1]], commits=cf["simulated_at_commits"], confirmatory=cf["confirmatory"],
                path=path)


# ----------------------------------------------------------------------------------------------
def provenance(v: dict, results_dir: Path) -> dict:
    """Inputs, description, alt text and caption draft; every number in them is computed from the plotted values."""
    labs, rows, alpha = v["labs"], v["rows"], v["alpha"]
    R = rel(results_dir)

    def rng(key):
        x = [rows[k][key]["value"] for k in labs if rows[k][key].get("value") is not None]
        return (f"{FG4.ratio(min(x))}–{FG4.ratio(max(x))}" if len(x) > 1 else (FG4.ratio(x[0]) if x else "not estimable"),
                len(labs) - len(x))

    def passed(key):
        return [k for k in labs if rows[k][key]["p_holm"] < alpha]

    names = lambda ks: ", ".join(style.ANAT_SHORT[k] for k in ks) if ks else "none"  # noqa: E731
    c_rng, c_cens = rng("confirmatory")
    x_rng, _ = rng("exploratory")
    o_rng, o_cens = rng("oracle")
    m_rng, m_cens = rng("mismatch")
    p_c, p_o, p_m = passed("confirmatory"), passed("oracle"), passed("mismatch")
    incl1 = [k for k in labs if (rows[k]["confirmatory"]["ci95"] or [None, None])[0] is None
             or rows[k]["confirmatory"]["ci95"][0] <= 1 <= (rows[k]["confirmatory"]["ci95"][1] or np.inf)]
    status = ("confirmatory (all nine anatomies, one commit, declared settings)" if v["confirmatory"]
              else "NOT confirmatory: " + "; ".join(v["reasons"]))
    cens_txt = lambda n: f" ({n} censored)" if n else ""  # noqa: E731
    return {NAME: dict(
        inputs=[f"{R}/g4_confirm_summary.json"] + [f"{R}/g4c_{k}_summary.json" for k in labs] + [EXPLO],
        description=(f"Run status: {status}; commits {', '.join(v['commits'])}. (a) g4_confirm_summary.json anatomy[<anatomy>]."
                     "endpoint.s50_ratio_squid_over_opm (value, value_bounds, value_censored, ci95), dark, and "
                     f"{EXPLO} comparison['<anatomy>/paired/opm_dense/opm_vs_squid/combined/practical@{v['rate']:g}/depth0']"
                     ".s50_ratio_squid_over_opm, light (checked equal to anatomy[<anatomy>].exploratory). Text: "
                     "locations_favouring_opm : locations_favouring_squid; Holm p from families['endpoint'].holm_p (dark) and "
                     "anatomy[<anatomy>].exploratory.holm_p_over_present_anatomies (light). (b) anatomy[<anatomy>].oracle and "
                     "(c) anatomy[<anatomy>].mismatch (.s50_ratio_squid_over_opm), with families['secondary/opm_dense/opm_vs_"
                     f"squid/combined/oracle/replicate0'] and ['.../mismatch@{v['rate']:g}/replicate0'].holm_p; grey bars: the "
                     "endpoint's point estimate. A tick and bold: Holm-adjusted p below the declared alpha "
                     f"({alpha:g}). Shared logarithmic x axis {v['lims'][0]:g}–{v['lims'][1]:g} (display choice: the drawn "
                     "values and finite interval ends with a margin); open interval ends are arrows to the axis edge. Anatomies "
                     "in the order of scripts/report_style.py ANAT_ORDER with its labels, classes separated by hairlines. "
                     f"Drawn as: {', '.join(v['kinds'])}."),
        alt=(f"Three forest plots of the S50 ratio Neuromag over OPM on a log axis with a line at 1, one row per anatomy "
             f"({len(labs)}). In the confirmatory run the dense array's ratio is {c_rng}{cens_txt(c_cens)}, against {x_rng} "
             f"in the exploratory run; the declared test passes after Holm correction in {len(p_c)} of {len(labs)} anatomies "
             f"({names(p_c)}), and the interval includes 1 for {names(incl1)}. The oracle detector's ratio is "
             f"{o_rng}{cens_txt(o_cens)} (Holm p below {alpha:g}: {names(p_o)}); the mismatched detector's {m_rng}"
             f"{cens_txt(m_cens)} (Holm p below {alpha:g}: {names(p_m)})."),
        caption_draft=(f"Figure R16. Confirmatory spike run. Strength for 50 % detection, Neuromag 306 over the dense OPM array, "
                       f"for focal spikes {v['band']} mm below the scalp, paired on identical simulated spikes at {v['n_loc']} "
                       f"newly drawn locations per anatomy (> 1: the OPM detects at a lower strength); lines: 95 % location-"
                       "bootstrap intervals; arrows: open ends; open markers: censored estimates. (a) The endpoint declared before "
                       f"the run (practical scanning detector, thresholds frozen at {v['rate']:g} false event per minute, first "
                       f"noise realization): ratios {c_rng}; the two-sided location sign-flip test, Holm-adjusted over the "
                       f"anatomies, is below {alpha:g} in {len(p_c)} of {len(labs)}. Light: the exploratory run "
                       f"({v['n_loc_exploratory']} locations per anatomy, endpoint chosen after the analyses), ratios {x_rng}. "
                       f"(b) Oracle detector (known topography, waveform and time): ratios {o_rng}. (c) Mismatched detector "
                       f"({v['mismatch']}), thresholds recalibrated to the same rate: ratios {m_rng}. Panels (b) and (c) are "
                       "secondary analyses, each Holm-adjusted within its own family; grey bars mark the endpoint's ratio."),
        values=dict(status=status, commits=v["commits"], alpha=alpha, false_events_per_min=v["rate"], band_mm=v["band"],
                    n_locations=v["n_loc"], n_events_per_replicate=v["n_ev"], bootstrap_resamples=v["n_boot"],
                    sign_flip_test=v["method"], oracle_alpha=v["oracle_alpha"], exploratory_n_locations=v["n_loc_exploratory"],
                    mismatch_variant=v["mismatch"],
                    x_limits=v["lims"], order=labs, rows=rows, drawn_as=v["drawn"]))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", type=Path, default=ROOT / "results" / "g4_confirm",
                    help="folder of the confirmatory run's outputs (default results/g4_confirm)")
    ap.add_argument("--out-dir", type=Path, default=style.OUT, help="folder for the figure and its provenance (default results/report)")
    args = ap.parse_args(argv)
    results_dir, out_dir = args.results_dir.resolve(), args.out_dir.resolve()  # a given path: relative to the working folder
    style.apply()
    v = figure(results_dir, out_dir)
    doc = {"provenance": {"commit": style.commit(), "results_dir": rel(results_dir), "confirmatory": v["confirmatory"]},
           "figures": provenance(v, results_dir)}
    (out_dir / OUT_JSON).write_text(json.dumps(doc, indent=1))
    print(f"{rel(v['path'])}\n{rel(out_dir / OUT_JSON)}" + ("" if v["confirmatory"] else f"\nNOT confirmatory: {'; '.join(v['reasons'])}"))


if __name__ == "__main__":
    main()
