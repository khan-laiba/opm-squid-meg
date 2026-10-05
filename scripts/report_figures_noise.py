#!/usr/bin/env python3
"""Report figure R17: how far the adult noise model can be trusted, drawn from stored outputs only
(nothing is simulated, re-fitted or re-tested).

  (a) Neuromag's known-topography detectability with its MEASURED noise covariance over that with the
      modelled covariance (median over the targets, 95 % parcel-bootstrap interval; finite-sample
      corrected), for the 305 good channels, the 102 magnetometers and the 203 gradiometers, with the
      heart's field removed, beside the range an exact model measured the same way would give (2.5th to
      97.5th percentile of the stored surrogate medians).
  (b) The implied OPM/Neuromag detectability ratio of the dense and the site-matched OPM arrays in the
      five scenarios of the covariance check, lettered A-E as in docs/methods.md and the register
      (A-E = S1-S5 of the result file; the lettering is checked against the methods text), and the
      published model's value.
  (c) The dense and site-matched ratios at the primary OPM white level with sensor plus brain noise and
      the omissions of the noise model added: near-skull cortex (as noise, as targets, both), coloured
      OPM noise (1/f corner 1, 3, 10 Hz, frequency-resolved, beside its white-noise reference), heart
      and eye sources at both stored levels, and the two joint runs; the published model as reference.
  (d) The ratio against the OPM white-noise level (the stored sweep, sensor plus brain noise and after
      the 8-term projection, the stored 95 % intervals as bands), the break-even levels with their
      intervals and the most adverse joint combination (dashed).
Every ratio axis carries a secondary axis in dB (20 log10 of the ratio; D for an OPM/Neuromag
detectability ratio), checked tick by tick against the ratio axis when drawn (as in R1, R3 and R4).

Inputs: results/g2_covariance_validation/covariance_validation.json, results/g2_noise_sensitivity/
noise_sensitivity_summary.json, results/g2/g2_summary.json (the published model) and docs/methods.md
(the scenario letters). Every plotted value is read from a named key and checked against the other keys
that store it (headline rows, the plain-text summary, the published summary); the script stops on any
mismatch. The checks and the plotted values are written with the figure's provenance.

Outputs (results/report/): Figure_R17_noise_checks.png (175 mm wide) and figures_noise.json (inputs,
description, alt text, caption draft, plotted values, checks, dB axes).

Usage: PYTHONPATH=src .venv/bin/python scripts/report_figures_noise.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend before pyplot is used)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402

from report_figures_adult import check_db_axes, db_axis, log_axis, to_db  # noqa: E402  (the dB axes of R1, R3 and R4)

COV = "results/g2_covariance_validation/covariance_validation.json"
NS = "results/g2_noise_sensitivity/noise_sensitivity_summary.json"
G2 = "results/g2/g2_summary.json"
METHODS = "docs/methods.md"
NAME = "Figure_R17_noise_checks"
OUT_JSON = "figures_noise.json"

# page width; the figure is drawn 0.5 mm narrower so that the saved PNG (tight bbox plus its 0.03-in pad) fits it (checked)
PAGE_MM, WIDTH_MM, HEIGHT_IN = 175.0, 174.5, 7.7
ARRAYS = ("opm_dense", "opm_matched")
ARRAY_TERM = {"opm_dense": "dense OPM array", "opm_matched": "site-matched OPM array"}  # the manuscript's terms
MARKER = {"opm_dense": "o", "opm_matched": "s"}
DODGE = {"opm_dense": -0.15, "opm_matched": 0.15}  # rows run downwards: the dense array above the site-matched
IB, PR = "combined/intrinsic+brain", "combined/projected"
C_MEAS, C_HEART, C_NULL = "#000000", "#D55E00", "0.72"  # Neuromag measured, heart removed (Okabe-Ito), exact-model range
FS_TITLE, FS_ROW, FS_HEAD, FS_NOTE = 9.0, 7.2, 7.0, 6.8
ROW_BAND = dict(color="0.91", lw=0, zorder=0)  # the published model's row in (b)
CI_LINE = dict(lw=1.3, solid_capstyle="butt", zorder=3)
# (b) the scenarios: letter, key in covariance_validation.json opm_implication.ratios (= its assumptions key), a short
# gloss (the caption defines each scenario in full)
SCENARIOS = (("A", "S1_same_relative_change", "common model error"),
             ("B", "S2_same_variance_excess_both_modelled", "both backgrounds × {k}"),
             ("C", "S3_neuromag_measured_opm_variance_excess", "OPM bears the excess (bound)"),
             ("D", "S4_neuromag_measured_opm_as_modelled", "Neuromag measured, OPM modelled"),
             ("E", "S5_neuromag_measured_empty_room", "measured empty room"))
# the phrases of docs/methods.md that define A-E, in this order (the letters must stay those of the methods)
METHODS_LETTERS = (("A", "the model's error common"), ("B", "both cortical backgrounds scaled"),
                   ("C", "a pessimistic bound, Neuromag as measured"), ("D", "Neuromag as measured, the OPM exactly as modelled"),
                   ("E", "Neuromag with its measured empty"))
CORNERS = (1, 3, 10)  # (c) 1/f corners, Hz
SWEEP_LS = {IB: "-", PR: (0, (1.0, 1.4))}  # (d) sensor plus brain noise solid, after the projection dotted
ADVERSE_LS = (0, (3.6, 2.0))  # (d) the most adverse joint combination, dashed
BE_DB_OFFSET = 0.13  # (d) break-even markers this many dB above (sensor plus brain noise) / below (projected) the ratio of 1


# ----------------------------------------------------------------------------------------------
def load():
    return tuple(json.loads((ROOT / p).read_text()) for p in (COV, NS, G2))


def close(a, b, tol=1e-12):
    return abs(float(a) - float(b)) <= tol


def same_entry(a, b, tol=1e-12):
    """Two stored comparisons with the same median ratio and interval."""
    return close(a["ratio"], b["ratio"], tol) and all(close(x, y, tol) for x, y in zip(a["ci95_ratio"], b["ci95_ratio"]))


def ri(e):
    """(ratio, lo, hi) of a stored comparison (ci95_ratio in the sensitivity run, ci95 in the covariance check)."""
    lo, hi = e["ci95_ratio"] if "ci95_ratio" in e else e["ci95"]
    return float(e["ratio"]), float(lo), float(hi)


def f2(v):
    return f"{v[0]:.2f}"


def fci(v, n=2):
    return f"{v[0]:.{n}f} [{v[1]:.{n}f}, {v[2]:.{n}f}]"


def fdb(r):
    return f"{to_db(r):+.2f} dB"


def published(g2):
    """The published model (the adult comparison's primary result): sensor plus brain noise, all 306 channels."""
    P = g2["primary"]["oracle"]
    return {a: (2.0 ** P[f"{a}/{IB}"]["median_log2"], *(2.0 ** x for x in P[f"{a}/{IB}"]["ci95"])) for a in ARRAYS}


def n_ext(ns):
    """The projection's term count, from the stored conventions ('8-term external subspace')."""
    return int(re.search(r"(\d+)-term", ns["conventions"]).group(1))


# ----------------------------------------------------------------------------------------------
def rows_layout(groups, gap=0.4):
    """y of every group header and row (rows run downwards from 0): [(header, y_header, [(row, y), ...]), ...]."""
    out, y = [], 0.0
    for head, rows in groups:
        hy = y if head else None
        y += 0.95 if head else 0.0
        ys = []
        for r in rows:
            ys.append((r, y))
            y += 1.0
        out.append((head, hy, ys))
        y += gap
    return out, y - gap - 1.0


def forest(ax, y, values, ms=4.0):
    """Both OPM arrays in one row: marker at the median, line over the 95 % interval."""
    for a in ARRAYS:
        v, lo, hi = values[a]
        col = style.ARRAY_COLOR[a]
        ax.plot([lo, hi], [y + DODGE[a]] * 2, "-", color=col, **CI_LINE)
        ax.plot(v, y + DODGE[a], MARKER[a], color=col, ms=ms, mec="white" if a == "opm_dense" else col, mew=0.4, zorder=4)


def reference(ax, pub):
    """The published model: a line at each array's median, its 95 % interval shaded."""
    for a in ARRAYS:
        ax.axvspan(pub[a][1], pub[a][2], color=style.ARRAY_COLOR[a], alpha=0.14, lw=0, zorder=0)
        ax.axvline(pub[a][0], color=style.ARRAY_COLOR[a], lw=0.8, alpha=0.9, zorder=1)


def headers(ax, layout):
    for head, hy, _ in layout:
        if head:
            ax.text(0.008, hy, head, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=FS_HEAD,
                    fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=0.4), zorder=6)


def row_axis(ax, layout, labels, last_y):
    ys = [y for _, _, rows in layout for _, y in rows]
    ax.set_yticks(ys, labels, fontsize=FS_ROW)
    for t in ax.get_yticklabels():
        t.set_multialignment("left")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(last_y + 0.62, -0.72)
    ax.spines["left"].set_visible(False)


def title(sf, text):
    sf.suptitle(text, x=0.0, ha="left", fontsize=FS_TITLE)


# ----------------------------------------------------------------------------------------------
def check_methods_letters():
    """docs/methods.md defines A-E in the order and with the phrases of METHODS_LETTERS (the covariance check)."""
    text = " ".join((ROOT / METHODS).read_text().split())
    anchor = "in five scenarios, A-E here (S1-S5 in the result file)"
    seg = text[text.index(anchor):text.index(anchor) + 900]
    pos = []
    for letter, phrase in METHODS_LETTERS:
        m = re.search(rf"\b{letter} {re.escape(phrase)}", seg)
        if m is None:
            raise ValueError(f"{METHODS}: scenario {letter} is no longer '{phrase}...'")
        pos.append(m.start())
    if pos != sorted(pos):
        raise ValueError(f"{METHODS}: the scenarios are no longer defined in the order A-E")
    return (f"{METHODS}: the covariance check ('{anchor}') defines " + "; ".join(f"{l} '{p} ...'" for l, p in METHODS_LETTERS)
            + ", in this order; A-E are drawn from S1-S5 of the result file in the same order")


def panel_a(sf, cov):
    """Neuromag's detectability with the measured covariance over the model, per channel set."""
    ax = sf.subplots()
    det, checks = cov["detectability"], []
    ch, comp, sur = det["channels"], det["comparisons"], det["surrogate_medians"]
    if ch["combined"] != ch["mag"] + ch["grad"]:
        raise ValueError(f"{COV}: detectability.channels {ch}")
    sets = (("combined", f"{ch['combined']} good channels"), ("mag", f"{ch['mag']} magnetometers"),
            ("grad", f"{ch['grad']} gradiometers"))
    vals = {}
    for j, (cs, lab) in enumerate(sets):
        em, eh = comp[f"measured_corrected_over_model/{cs}"], comp[f"without_heart_corrected_over_model/{cs}"]
        m, h, nl = ri(em), ri(eh), sur[f"null_corrected/{cs}"]
        assert em["n"] == eh["n"] == comp["measured_corrected_over_model/combined"]["n"]
        assert em["ci_method"] == eh["ci_method"]
        vals[cs] = dict(label=lab, measured=m, heart_removed=h, exact_model_range=[nl[0], nl[2]], exact_model_median=nl[1],
                        ci_method=em["ci_method"])
        line = next(x for x in cov["plain_summary"] if x.startswith(f"Neuromag {cs}: "))  # quotes the same values
        want = (f"{m[0]:.3f} [{m[1]:.3f}, {m[2]:.3f}]", f"would give {nl[0]:.3f}-{nl[2]:.3f}", f"heart removed {h[0]:.3f}")
        if not all(w in line for w in want):
            raise ValueError(f"{COV}: plain_summary does not quote {want}: {line}")
        y = float(j)
        ax.fill_betweenx([y - 0.38, y + 0.38], nl[0], nl[2], color=C_NULL, lw=0, zorder=1)
        for v, dy, col, mk in ((m, -0.14, C_MEAS, "o"), (h, 0.14, C_HEART, "D")):
            ax.plot([v[1], v[2]], [y + dy] * 2, "-", color=col, **CI_LINE)
            ax.plot(v[0], y + dy, mk, color=col, ms=4.0 if mk == "o" else 3.6, zorder=4)
    checks += [f"{COV}: detectability.channels combined = mag + grad ({ch['combined']} = {ch['mag']} + {ch['grad']})",
               f"{COV}: plain_summary quotes the plotted measured/model ratio, its interval, the exact-model range and the "
               "heart-removed ratio of all three channel sets (3 decimals)"]
    log_axis(ax, "x", [1, 1.05, 1.1, 1.15, 1.2], (0.98, 1.235))
    db_axis(ax, "x", 0.5, "dB")
    ax.set_yticks(range(len(sets)), [lab for _, lab in sets], fontsize=FS_ROW)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(5.3, -0.6)
    ax.spines["left"].set_visible(False)
    n_t = comp["measured_corrected_over_model/combined"]["n"]
    ax.set_xlabel(f"Neuromag detectability, measured over\nmodelled noise (median over {n_t:,} targets)")
    h = [Line2D([], [], color=C_MEAS, marker="o", ms=4, lw=1.3, label="measured noise\n(finite-sample corrected)"),
         Line2D([], [], color=C_HEART, marker="D", ms=3.6, lw=1.3, label="the same, heart's\nfield removed"),
         Patch(color=C_NULL, lw=0, label="an exact model, measured\nthe same way (2.5–97.5 %)")]
    ax.legend(handles=h, loc="lower right", fontsize=FS_NOTE, handlelength=1.7, borderaxespad=0.15, labelspacing=0.6,
              frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    title(sf, "(a) Neuromag: measured vs model")
    return ax, vals, checks


def panel_b(sf, cov, pub):
    """The implied OPM/Neuromag ratio in the five scenarios, and the published model."""
    ax = sf.subplots()
    imp, checks = cov["opm_implication"], []
    r, k = imp["ratios"], imp["checks"]["magnetometer_variance_excess_k"]
    if not close(k, cov["model"]["magnetometer_calibration_scale_ratio"]):
        raise ValueError(f"{COV}: magnetometer_variance_excess_k {k} != model.magnetometer_calibration_scale_ratio")
    rows = [("published", "Primary model\n(sensor plus brain noise)", pub)]
    for letter, key, lab in SCENARIOS:
        if key not in imp["assumptions"]:
            raise ValueError(f"{COV}: opm_implication.assumptions has no '{key}'")
        rows.append((letter, f"{letter}  " + lab.format(k=f"{k:.2f}").replace("\n", "\n      "),
                     {a: ri(r[f"{a}/combined/{key}"]) for a in ARRAYS}))
    for a in ARRAYS:
        s1, m305 = r[f"{a}/combined/S1_same_relative_change"], r[f"{a}/combined/model_305"]
        if not (close(s1["ratio"], m305["ratio"]) and s1["ci95"] == m305["ci95"]):
            raise ValueError(f"{COV}: scenario A of {a} is not the model's own ratio (model_305)")
        dev = abs(r[f"{a}/combined/model_306/intrinsic+brain"]["ratio"] / pub[a][0] - 1)
        if dev > 1e-4:
            raise ValueError(f"{COV}: model_306/intrinsic+brain of {a} does not reproduce {G2} ({dev:.1e})")
        checks.append(f"{COV}: opm_implication.ratios['{a}/combined/model_306/intrinsic+brain'].ratio reproduces the published "
                      f"median ({G2} primary.oracle) to {dev:.1e} (relative)")
    checks += [f"{COV}: scenario A (S1_same_relative_change) is identical to model_305, the model's own ratio, for both arrays",
               f"{COV}: opm_implication.checks.magnetometer_variance_excess_k = model.magnetometer_calibration_scale_ratio "
               f"= {k:.4f}"]
    n_quoted = 0
    for line in cov["plain_summary"]:  # the plain-text summary quotes the same scenario values
        if " / Neuromag (all channels, sensor + brain + room): " in line:
            a = line.split(" / ")[0]
            for letter, key, _ in SCENARIOS:
                v = ri(r[f"{a}/combined/{key}"])
                if f"{v[0]:.3f} [{v[1]:.3f}, {v[2]:.3f}]" not in line:
                    raise ValueError(f"{COV}: plain_summary does not quote scenario {letter} of {a}: {line}")
                n_quoted += 1
    assert n_quoted == len(ARRAYS) * len(SCENARIOS), n_quoted
    checks.append(f"{COV}: plain_summary quotes every plotted scenario value with its interval (3 decimals)")
    layout, last = rows_layout([(None, rows[:1]), ("A–E: sensor, brain and room noise", rows[1:])], gap=0.2)
    reference(ax, pub)
    for head, hy, ys in layout:
        for (id_, _, values), y in ys:
            if id_ == "published":
                ax.axhspan(y - 0.48, y + 0.48, **ROW_BAND)
            forest(ax, y, values)
    headers(ax, layout)
    log_axis(ax, "x", [0.7, 0.8, 0.9, 1, 1.1, 1.2], (0.625, 1.24))
    db_axis(ax, "x", 1.0, "D (dB)")
    row_axis(ax, layout, [lab for _, lab, _ in rows], last)
    ax.get_yticklabels()[0].set_fontweight("bold")
    ax.set_xlabel(f"OPM / Neuromag detectability\n(median over the {r['opm_dense/combined/model_305']['n']:,} targets)")
    title(sf, "(b) Implied OPM / Neuromag ratio, scenarios A–E")
    vals = {id_: dict(label=lab.replace("\n      ", " ").replace("\n", " "), **{a: list(v[a]) for a in ARRAYS})
            for id_, lab, v in rows}
    vals["k_magnetometer_variance_excess"] = k
    return ax, vals, checks


def panel_c(sf, ns, pub):
    """The omissions of the noise model at the primary white level (sensor plus brain noise)."""
    ax = sf.subplots()
    hl, checks = ns["headline"], []
    lvl = ns["config"]["adult"]["sensors"]["opm_asd_primary_fT_per_rtHz"]
    lk = f"{lvl:g}"
    sweep = ns["sweep"]["entries"][lk]["comparisons"]
    for a in ARRAYS:  # the reference: the published model, as stored in the sweep and in the published summary
        e = sweep[f"{a}/{IB}"]
        if not (close(e["ratio"], pub[a][0]) and close(e["stored_g2"]["ratio"], pub[a][0])
                and same_entry(e, hl[f"published model, OPM white {lk} fT/sqrt(Hz)"][f"{a}/{IB}"])):
            raise ValueError(f"{NS}: sweep.entries['{lk}'] of {a} is not the published model")
    checks.append(f"{NS}: sweep.entries['{lk}'] (= headline 'published model, OPM white {lk} fT/sqrt(Hz)') reproduces the "
                  f"published medians of {G2} exactly")
    fl = ns["near_skull"]["primary"]["floor_2mm"]
    rule = ns["config"]["sensitivity"]["near_skull"]["rule_mm"]
    co, ff, jr = ns["coloured"]["results"], ns["far_field"]["results"], ns["joint"]["results"]
    jcfg = ns["config"]["sensitivity"]["joint"]
    if not ns["joint"]["most_adverse"].startswith("joint_without_far_field:"):
        raise ValueError(f"{NS}: joint.most_adverse is no longer joint_without_far_field: {ns['joint']['most_adverse']}")
    assert jcfg["corner_hz"] == max(CORNERS) and jcfg["near_skull_floor_mm"] == fl["floor_mm"] == 2.0, jcfg
    assert jcfg["level"] == "magnetometer_shortfall", jcfg
    corner = f"{jcfg['corner_hz']:g}"
    # (label, stored comparisons, headline key with the same values, path of the stored comparisons)
    groups = [
        (f"Cortex {fl['floor_mm']:g}–{rule:g} mm from the inner skull", [
            ("as noise", fl["background"]["comparisons"], "near-skull cortex (>= 2 mm) in the background",
             "near_skull.primary.floor_2mm.background.comparisons"),
            ("as targets", fl["targets"]["all"], "near-skull cortex (>= 2 mm) as additional targets",
             "near_skull.primary.floor_2mm.targets.all"),
            ("as noise and targets", fl["both"]["all"], "near-skull cortex (>= 2 mm) in the background and as targets",
             "near_skull.primary.floor_2mm.both.all")]),
        ("Coloured OPM noise: 1/f corner", [
            ("none, frequency-resolved", co[f"frequency_resolved/flat/{lk}fT/0Hz"]["comparisons"],
             "OPM 1/f corner 0 Hz, frequency-resolved (flat signal)",
             f"coloured.results['frequency_resolved/flat/{lk}fT/0Hz'].comparisons")]
         + [(f"{c:g} Hz", co[f"frequency_resolved/flat/{lk}fT/{c}Hz"]["comparisons"],
             f"OPM 1/f corner {c} Hz, frequency-resolved (flat signal)",
             f"coloured.results['frequency_resolved/flat/{lk}fT/{c}Hz'].comparisons") for c in CORNERS]),
        ("Heart and eyes added to the noise", [
            ("at the room field's level", ff["cardiac+ocular/room_field/added"]["comparisons"],
             "heart + eyes, room field level, added", "far_field.results['cardiac+ocular/room_field/added'].comparisons"),
            ("at the magnetometer shortfall", ff["cardiac+ocular/magnetometer_shortfall/added"]["comparisons"],
             "heart + eyes, magnetometer shortfall level, added",
             "far_field.results['cardiac+ocular/magnetometer_shortfall/added'].comparisons")]),
        ("Joint, frequency-resolved", [
            (f"near-skull cortex + {corner}-Hz corner", jr[f"joint_without_far_field/{lk}fT"]["comparisons"],
             f"near-skull cortex + 1/f corner (most adverse combination), OPM white {lk} fT/sqrt(Hz)",
             f"joint.results['joint_without_far_field/{lk}fT'].comparisons"),
            ("the same + heart and eyes", jr[f"joint/{lk}fT"]["comparisons"],
             f"all alternatives together (far field filling the magnetometer shortfall), OPM white {lk} fT/sqrt(Hz)",
             f"joint.results['joint/{lk}fT'].comparisons")]),
    ]
    rows, vals = [], {}
    for head, items in groups:
        rr = []
        for lab, comps, hk, path in items:
            for a in ARRAYS:
                if not same_entry(comps[f"{a}/{IB}"], hl[hk][f"{a}/{IB}"]):
                    raise ValueError(f"{NS}: {path} differs from headline['{hk}'] ({a})")
            v = {a: ri(comps[f"{a}/{IB}"]) for a in ARRAYS}
            rr.append((lab, v))
            vals[f"{head} / {lab}"] = dict(source=f"{NS} :: {path}['<array>/{IB}']", headline=hk,
                                           **{a: list(v[a]) for a in ARRAYS})
        rows.append((head, rr))
    checks += [f"{NS}: every plotted row of (c) equals its headline row (ratio and interval) for both arrays",
               f"{NS}: joint.most_adverse names joint_without_far_field; config.sensitivity.joint: near-skull floor "
               f"{jcfg['near_skull_floor_mm']:g} mm, corner {corner} Hz, far field at the {jcfg['level']} level"]
    layout, last = rows_layout(rows)
    reference(ax, pub)
    for head, hy, ys in layout:
        for (lab, v), y in ys:
            forest(ax, y, v)
    headers(ax, layout)
    log_axis(ax, "x", [0.95, 1, 1.05, 1.1, 1.15], (0.935, 1.2))
    db_axis(ax, "x", 0.5, "D (dB)")
    row_axis(ax, layout, [lab for _, rr in rows for lab, _ in rr], last)
    ax.set_xlabel(f"OPM / Neuromag detectability ({ns['setup']['arrays']['squid']} channels),\n"
                  f"OPM {lk} fT/√Hz, sensor plus brain noise")
    title(sf, f"(c) Omissions of the noise model, at {lk} fT/√Hz")
    return ax, vals, checks


def panel_d(sf, ns, pub):
    """The ratio against the OPM white-noise level, the break-even levels and the most adverse joint combination."""
    ax = sf.subplots()
    sw, checks = ns["sweep"], []
    levels = [float(x) for x in sw["levels_fT"]]
    if levels != [float(x) for x in ns["config"]["adult"]["sensors"]["opm_asd_fT_per_rtHz"]] or \
            list(sw["entries"]) != [f"{x:g}" for x in levels]:
        raise ValueError(f"{NS}: sweep levels {levels}")
    hl, jr, be = ns["headline"], ns["joint"]["results"], sw["break_even"]
    vals = dict(levels_fT=levels, sweep={}, most_adverse={}, break_even={})
    x = np.array(levels)
    for a in ARRAYS:
        col = style.ARRAY_COLOR[a]
        for key in (IB, PR):
            e = [sw["entries"][f"{lv:g}"]["comparisons"][f"{a}/{key}"] for lv in levels]
            for lv, ee in zip(levels, e):
                if not same_entry(ee, hl[f"published model, OPM white {lv:g} fT/sqrt(Hz)"][f"{a}/{key}"]):
                    raise ValueError(f"{NS}: sweep.entries['{lv:g}'] of {a}/{key} differs from its headline row")
            v = np.array([ri(ee) for ee in e])
            vals["sweep"][f"{a}/{key}"] = v.tolist()
            ax.fill_between(x, v[:, 1], v[:, 2], color=col, alpha=0.17 if key == IB else 0.10, lw=0, zorder=1)
            ax.plot(x, v[:, 0], ls=SWEEP_LS[key], color=col, lw=1.5, zorder=3, marker=MARKER[a] if key == IB else None,
                    ms=3.4, mec="white" if a == "opm_dense" else col, mew=0.4)
            b = be["entries"][f"{a}/{key}"]  # the level at which the median ratio is 1
            bv, (blo, bhi) = b["break_even_fT"], b["ci95_fT"]
            if b["bound"] is not None or b["resamples_outside_grid"]:
                raise ValueError(f"{NS}: break-even of {a}/{key} is bounded or partly outside the grid: {b}")
            crossing = all(r > 1 for lv, r in zip(levels, v[:, 0]) if lv < bv) and \
                all(r < 1 for lv, r in zip(levels, v[:, 0]) if lv > bv)
            if not crossing:
                raise ValueError(f"{NS}: the break-even level of {a}/{key} ({bv:.2f}) does not match the sweep's crossing")
            vals["break_even"][f"{a}/{key}"] = [bv, blo, bhi]
        adv = [jr[f"joint_without_far_field/{lv:g}fT"]["comparisons"][f"{a}/{IB}"] for lv in levels]
        for lv, ee in zip(levels, adv):
            hk = f"near-skull cortex + 1/f corner (most adverse combination), OPM white {lv:g} fT/sqrt(Hz)"
            if not same_entry(ee, hl[hk][f"{a}/{IB}"]):
                raise ValueError(f"{NS}: joint.results['joint_without_far_field/{lv:g}fT'] differs from its headline row")
        va = np.array([ri(ee) for ee in adv])
        vals["most_adverse"][a] = va.tolist()
        ax.plot(x, va[:, 0], ls=ADVERSE_LS, color=col, lw=1.25, zorder=3)
    checks.append(f"{NS}: sweep.levels_fT = config.adult.sensors.opm_asd_fT_per_rtHz; every plotted sweep and most-adverse "
                  "value equals its headline row; every break-even level lies where the plotted sweep median crosses 1, "
                  "unbounded and with no resample outside its grid")
    for a in ARRAYS:  # break-even markers just above (sensor plus brain noise) or below (projected) the ratio of 1
        col = style.ARRAY_COLOR[a]
        for key, sgn in ((IB, 1), (PR, -1)):
            bv, blo, bhi = vals["break_even"][f"{a}/{key}"]
            yb = 10 ** (sgn * BE_DB_OFFSET / 20)
            ax.plot([blo, bhi], [yb, yb], "-", color=col, lw=2.4, alpha=0.6, solid_capstyle="butt", zorder=4)
            ax.plot(bv, yb, "D", color=col, ms=4.4, mfc=col if key == IB else "white", mec=col, mew=1.1, zorder=5)
            ax.annotate(f"{bv:.1f}", (bv, yb), xytext=(0, 4.5 * sgn), textcoords="offset points", ha="center",
                        va="bottom" if sgn > 0 else "top", fontsize=FS_NOTE, color=col, zorder=6,
                        bbox=dict(facecolor="white", edgecolor="none", pad=0.1, alpha=0.85))
    prim = ns["config"]["adult"]["sensors"]["opm_asd_primary_fT_per_rtHz"]
    ax.axvline(prim, color="0.7", lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(5.0, 40.0)
    ticks = [5, 7, 10, 15, 20, 30, 40]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    log_axis(ax, "y", [0.8, 0.9, 1, 1.1, 1.2, 1.3], (0.765, 1.36))
    db_axis(ax, "y", 1.0, "D (dB)")
    ax.set_xlabel("OPM white noise level (fT/√Hz)")
    ax.set_ylabel(f"OPM / Neuromag detectability\n({ns['setup']['arrays']['squid']} channels), log scale")
    h = [Line2D([], [], color="0.25", lw=1.5, ls=SWEEP_LS[IB], marker="o", ms=3.4, label="sensor plus brain noise"),
         Line2D([], [], color="0.25", lw=1.5, ls=SWEEP_LS[PR], label=f"after the {n_ext(ns)}-term projection"),
         Patch(color="0.5", alpha=0.28, lw=0, label="95 % intervals (bands)"),
         Line2D([], [], color="0.25", lw=1.25, ls=ADVERSE_LS,
                label=f"joint run: near-skull cortex\n+ {ns['config']['sensitivity']['joint']['corner_hz']:g}-Hz corner, "
                      "frequency-resolved"),
         Line2D([], [], color="0.25", lw=0, marker="D", ms=4.4, label="break-even level, 95 % interval (filled: sensor\n"
                "plus brain noise; open: after the projection)"),
         Line2D([], [], color="0.7", lw=0.8, label=f"primary model's level ({prim:g} fT/√Hz)")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=FS_NOTE, handlelength=2.4,
              borderaxespad=0.0, labelspacing=0.45, ncol=1)
    title(sf, "(d) OPM white-noise level")
    return ax, vals, checks


# ----------------------------------------------------------------------------------------------
def texts(v, cov, ns, g2, pub):
    """Inputs, description, alt text and caption draft, every number from the plotted values."""
    a_, b_, c_, d_ = v["a"], v["b"], v["c"], v["d"]
    dec = cov["declared_choices"]
    n_t = cov["detectability"]["comparisons"]["measured_corrected_over_model/combined"]["n"]
    lk = f"{ns['config']['adult']['sensors']['opm_asd_primary_fT_per_rtHz']:g}"
    nx, k = n_ext(ns), b_["k_magnetometer_variance_excess"]
    sq = g2["arrays"]["squid"]["channels"]
    sites = {a: g2["arrays"][a]["n_sites"] for a in ARRAYS}
    lv = d_["levels_fT"]
    be, sw, adv = d_["break_even"], d_["sweep"], d_["most_adverse"]
    cv = {key.split(" / ")[1]: val for key, val in c_.items()}
    grid = ns["sweep"]["break_even"]["grid_fT"]
    ff = ns["far_field"]["level_values_fT"]
    jw = cv[f"near-skull cortex + {ns['config']['sensitivity']['joint']['corner_hz']:g}-Hz corner"]
    jf = cv["the same + heart and eyes"]
    room, short = cv["at the room field's level"], cv["at the magnetometer shortfall"]
    near = [cv[x] for x in ("as noise", "as targets", "as noise and targets")]
    corner_max = f"{max(CORNERS):g}"
    corners = ", ".join(f"{c:g}" for c in CORNERS[:-1]) + f" and {CORNERS[-1]:g}"
    # the parcel labels of the bootstrap: the same in both result files
    methods = {a_["combined"]["ci_method"], ns["sweep"]["entries"][lk]["comparisons"][f"opm_dense/{IB}"]["ci_method"]}
    if len(methods) != 1:
        raise ValueError(f"different bootstrap groups: {methods}")
    n_lab = int(re.fullmatch(r"parcels \((\d+)\)", methods.pop()).group(1))
    bads = ", ".join(cov["data"]["bad_channels"])
    jr = ns["joint"]["results"]

    def joint_side(a):
        """Where the joint run with the heart and eyes lies against the most adverse combination, level by level."""
        low = [x for x in lv if jr[f"joint/{x:g}fT"]["comparisons"][f"{a}/{IB}"]["ratio"]
               < jr[f"joint_without_far_field/{x:g}fT"]["comparisons"][f"{a}/{IB}"]["ratio"]]
        if len(low) == len(lv):
            return "below it at every level"
        if not low:
            return "above it at every level"
        if len(low) <= len(lv) / 2:
            return "above it except at " + ", ".join(f"{x:g}" for x in low) + " fT/sqrt(Hz)"
        return "below it except at " + ", ".join(f"{x:g}" for x in lv if x not in low) + " fT/sqrt(Hz)"

    assert all(x["measured"][1] > 1 for x in a_.values())  # wording guard (alt text: 'intervals above 1')
    # (c): the largest changes against the primary model, per array (%)
    chg = {a: [100 * (x[a][0] / pub[a][0] - 1) for x in c_.values()] for a in ARRAYS}
    fr_white = {a: 100 * (cv["none, frequency-resolved"][a][0] / pub[a][0] - 1) for a in ARRAYS}

    def pct(x):
        return f"{x:+.1f}".replace("-", "\u2212")

    def both(val, n=2):
        return f"{val['opm_dense'][0]:.{n}f} / {val['opm_matched'][0]:.{n}f}"

    def bci(val, n=2):
        return f"{fci(val['opm_dense'], n)} / {fci(val['opm_matched'], n)}"

    def bev(key):
        return f"{be[key][0]:.1f} [{be[key][1]:.1f}, {be[key][2]:.1f}]"

    inputs = [
        f"{COV} :: detectability.comparisons['measured_corrected_over_model/<combined|mag|grad>'] (ratio, ci95, n, ci_method)",
        f"{COV} :: detectability.comparisons['without_heart_corrected_over_model/<combined|mag|grad>'] (ratio, ci95)",
        f"{COV} :: detectability.surrogate_medians['null_corrected/<combined|mag|grad>'] (2.5th, 50th, 97.5th percentile)",
        f"{COV} :: detectability.channels; declared_choices.n_surrogates, n_parcel_bootstrap",
        f"{COV} :: opm_implication.ratios['<array>/combined/<S1_same_relative_change|S2_same_variance_excess_both_modelled|"
        "S3_neuromag_measured_opm_variance_excess|S4_neuromag_measured_opm_as_modelled|S5_neuromag_measured_empty_room>'] "
        "(ratio, ci95) = scenarios A-E",
        f"{COV} :: opm_implication.ratios['<array>/combined/model_305'], ['<array>/combined/model_306/intrinsic+brain'] (checks)",
        f"{COV} :: opm_implication.assumptions, opm_implication.checks.magnetometer_variance_excess_k, "
        "model.magnetometer_calibration_scale_ratio",
        f"{COV} :: plain_summary (checks)",
        f"{G2} :: primary.oracle['<array>/combined/intrinsic+brain'] (median_log2, ci95) = the primary model",
        f"{G2} :: arrays.<array>.n_sites, arrays.squid.channels",
        f"{NS} :: sweep.levels_fT, sweep.entries['<level>'].comparisons['<array>/combined/<intrinsic+brain|projected>'] "
        "(ratio, ci95_ratio, stored_g2)",
        f"{NS} :: sweep.break_even.entries['<array>/combined/<intrinsic+brain|projected>'] (break_even_fT, ci95_fT, bound, "
        "resamples_outside_grid), sweep.break_even.grid_fT, sweep.break_even.method",
        f"{NS} :: near_skull.primary.floor_2mm.<background.comparisons|targets.all|both.all>['<array>/combined/intrinsic+brain']",
        f"{NS} :: coloured.results['frequency_resolved/flat/{lk}fT/<0|1|3|10>Hz'].comparisons"
        "['<array>/combined/intrinsic+brain']",
        f"{NS} :: far_field.results['cardiac+ocular/<room_field|magnetometer_shortfall>/added'].comparisons"
        "['<array>/combined/intrinsic+brain'], far_field.level_values_fT",
        f"{NS} :: joint.results['<joint_without_far_field|joint>/<level>fT'].comparisons['<array>/combined/intrinsic+brain'], "
        "joint.most_adverse, joint.definition",
        f"{NS} :: headline (every plotted row; checks)",
        f"{NS} :: config.adult.sensors.<opm_asd_fT_per_rtHz|opm_asd_primary_fT_per_rtHz>, config.sensitivity.<near_skull|joint>, "
        "setup.arrays.squid, setup.n_boot, conventions (projection terms)",
        f"{METHODS} :: section 8, the covariance check's scenarios A-E (lettering; checked)",
    ]
    description = (
        f"How far the adult noise model (MNE sample subject, {n_t:,} cortical targets, 1-40 Hz) can be trusted. Ratios are "
        "medians over the targets of a paired ratio of known-topography detectability (oracle covariance) with 95 % "
        f"parcel-bootstrap intervals ({dec['n_parcel_bootstrap']:,} resamples in the covariance check, "
        f"{ns['setup']['n_boot']:,} in the sensitivity run; {n_lab} labels: the Desikan-Killiany parcels and the "
        "medial-wall labels); every ratio axis is logarithmic with a line at 1 and carries the same ratio in dB (20 log10 "
        "of the ratio; "
        "for an OPM/Neuromag detectability ratio this is D) on its top or right-hand axis. "
        f"(a) Neuromag's detectability with the sample covariance of the measured noise (the task recording's pre-stimulus "
        f"windows; {bads} left out) over that with the modelled covariance, divided per target by the finite-sample bias "
        "of a covariance estimated from that many samples (detectability.comparisons['measured_corrected_over_model/<set>']), "
        "for the good channels, the magnetometers and the gradiometers; orange: the same with the heart's field "
        "(cardiac-locked average) removed (without_heart_corrected_over_model); grey: the 2.5th-97.5th percentile of the "
        f"medians that an exact model gives when 'measured' the same way ({dec['n_surrogates']} surrogates recoloured to "
        "the model and corrected as the data; surrogate_medians['null_corrected/<set>']), so a value outside the grey range "
        "is a real departure of the measured noise from the model. "
        "(b) The OPM/Neuromag ratio of the dense and site-matched arrays that this implies (sensor, brain and room noise, "
        "like for like with the measured covariance, which holds the room field; Neuromag's "
        f"{cov['detectability']['channels']['combined']} good channels), in the five scenarios of the covariance check "
        "(opm_implication.ratios; the result file's S1-S5 are the methods' A-E): A, the model's error common to both systems "
        "(the OPM's real covariance changes its detectability by the same factor as Neuromag's; identical to the model's own "
        f"ratio); B, both systems' cortical backgrounds scaled by the measured magnetometer variance excess k = {k:.2f}; C, a "
        "pessimistic bound, Neuromag with its measured covariance and the OPM bearing that whole excess as cortical noise; "
        "D, Neuromag as measured, the OPM exactly as modelled; E, Neuromag with its measured empty room and the modelled "
        f"brain noise, the OPM as modelled. Top row (grey): the primary model, sensor plus brain noise against all {sq} "
        f"channels ({G2} primary.oracle), also drawn as a line per array with its interval shaded in (b) and (c). "
        f"(c) The ratio at the primary OPM white level ({lk} fT/sqrt(Hz)) with sensor plus brain noise when what the model "
        "leaves out is added (the near-skull and coloured-noise variants recalibrate the cortical background on the "
        "measured gradiometer brain noise, as the primary model does): the cortex 2-4 mm from the inner skull (the valid "
        "vertices closer than 4 mm to it, "
        "left out of the primary model; below 2 mm the lead fields are not converged) in the background, as additional "
        f"targets or both (near_skull.primary.floor_2mm); coloured OPM noise with a 1/f corner at {corners} Hz under the "
        "frequency-resolved detector (coloured.results "
        "'frequency_resolved/flat'), whose white-noise value (corner 'none') is its reference and lies below the primary "
        "band-variance ratio; heart and eye sources (a cardiac current dipole 250 mm below the head and current dipoles at "
        f"the two eyes, each half the far-field variance) added on top of the primary model's cortical background at the room "
        "field's level "
        f"({ff['room_field']:.0f} fT at the median magnetometer) or filling the magnetometer brain-noise shortfall "
        f"({ff['magnetometer_shortfall']:.0f} fT; far_field.results['cardiac+ocular/<level>/added'], the mode adverse to the "
        f"OPM); and the joint runs (frequency-resolved): the near-skull background with the {corner_max}-Hz corner "
        "(joint_without_far_field, which the result file names the most adverse combination) and the same with the heart "
        "and eyes at the shortfall level (joint). "
        f"(d) The ratio against the OPM white level ({', '.join(f'{x:g}' for x in lv)} fT/sqrt(Hz); sweep.entries) with "
        f"sensor plus brain noise (solid, markers) and after the {nx}-term projection (dotted), the stored 95 % intervals as "
        "bands; the break-even level, the white level at which the median ratio is 1 (sweep.break_even: crossing by "
        f"log-linear interpolation on {len(grid)} log-spaced levels from {grid[0]:g} to {grid[-1]:g} fT/sqrt(Hz), interval "
        "from the same parcel resamples applied to the whole curve), drawn as a diamond on the ratio of 1 with its interval "
        f"(filled and {BE_DB_OFFSET:g} dB above the line: sensor plus brain noise; open and {BE_DB_OFFSET:g} dB below: after "
        "the projection; the offsets only keep the intervals apart); dashed: the most adverse combination "
        "(joint_without_far_field, sensor plus brain noise, frequency-resolved); grey vertical line: the primary model's "
        f"level. The joint run with the heart and eyes added (joint) lies, for the site-matched array, "
        f"{joint_side('opm_matched')} ({jf['opm_matched'][0]:.3f} against {jw['opm_matched'][0]:.3f} at {lk} fT/sqrt(Hz); "
        f"panel c) and, for the dense array, {joint_side('opm_dense')}. No OPM noise was measured: (a) checks the model for "
        "Neuromag only, "
        "and (b) bounds the OPM side under stated assumptions rather than estimating it.")
    alt = (
        "Four panels. (a) A small forest plot for Neuromag: with its measured noise, detectability is "
        f"{f2(a_['combined']['measured'])} times the modelled for all good channels, {f2(a_['mag']['measured'])} for the "
        f"magnetometers and {f2(a_['grad']['measured'])} for the gradiometers, with intervals above 1, while an exact model "
        "would give a narrow grey range at 1; removing the heart's field changes little. (b) Forest plot of the implied "
        f"OPM/Neuromag ratio for the dense and site-matched arrays: the primary model {both(b_['published'])}, scenarios A, "
        f"B and E about the same ({both(b_['A'])}, {both(b_['B'])}, {both(b_['E'])}), D lower ({both(b_['D'])}) and C, the "
        f"pessimistic bound, lowest ({both(b_['C'])}). (c) Forest plot of the ratios at {lk} fT/sqrt(Hz) with each omission "
        "of the noise model added: near-skull cortex, coloured OPM noise, heart and eye fields and the joint runs move the "
        f"dense ratio between {min(x['opm_dense'][0] for x in c_.values()):.2f} and "
        f"{max(x['opm_dense'][0] for x in c_.values()):.2f} and the site-matched between "
        f"{min(x['opm_matched'][0] for x in c_.values()):.2f} and {max(x['opm_matched'][0] for x in c_.values()):.2f}. "
        "(d) Line chart of the ratio against the OPM white noise level from "
        f"{lv[0]:g} to {lv[-1]:g} fT/sqrt(Hz): the dense ratio falls from {sw['opm_dense/' + IB][0][0]:.2f} to "
        f"{sw['opm_dense/' + IB][-1][0]:.2f} and reaches 1 at {be['opm_dense/' + IB][0]:.1f} fT/sqrt(Hz), the site-matched "
        f"from {sw['opm_matched/' + IB][0][0]:.2f} to {sw['opm_matched/' + IB][-1][0]:.2f}, reaching 1 at "
        f"{be['opm_matched/' + IB][0]:.1f}; after the projection both are lower, and the most adverse combination (dashed) "
        "lies below the sensor plus brain noise curves. All ratio axes also carry a dB scale.")
    caption = (
        "Figure R17. How far the noise model can be trusted. (a) Neuromag's detectability with its measured noise over "
        f"that with the modelled noise (median over the {n_t:,} targets, finite-sample corrected): "
        f"{fci(a_['combined']['measured'])} "
        f"({fdb(a_['combined']['measured'][0])}) for its {cov['detectability']['channels']['combined']} good channels, "
        f"{fci(a_['mag']['measured'])} for the {cov['detectability']['channels']['mag']} magnetometers and "
        f"{fci(a_['grad']['measured'])} for the {cov['detectability']['channels']['grad']} gradiometers; "
        f"{f2(a_['combined']['heart_removed'])}, {f2(a_['mag']['heart_removed'])} and {f2(a_['grad']['heart_removed'])} "
        f"with the heart's field removed, where an exact model measured the same way would give "
        f"{min(x['exact_model_range'][0] for x in a_.values()):.3f}-{max(x['exact_model_range'][1] for x in a_.values()):.3f} "
        "(grey). The measured noise is more favourable to Neuromag than the model, by about as much as the dense array's "
        f"modelled advantage ({fdb(pub['opm_dense'][0])}). (b) The OPM/Neuromag ratio this implies (dense / site-matched; "
        "sensor, brain and room noise) depends on how the OPM's real noise departs from its model, which no measurement here "
        f"constrains: A, the model's error common to both systems, {both(b_['A'])} (the model's own ratio); B, both cortical "
        f"backgrounds scaled by the magnetometer variance excess (× {k:.2f}), {both(b_['B'])}; C, the pessimistic bound, "
        f"Neuromag as measured and the OPM bearing the whole excess, {bci(b_['C'])}; D, Neuromag as measured and the OPM as "
        f"modelled, {bci(b_['D'])}; E, Neuromag's measured empty room with the modelled brain noise, {both(b_['E'])}; the "
        f"primary model (sensor plus brain noise, grey row and vertical lines) {bci(b_['published'])}. (c) At {lk} fT/√Hz "
        "with sensor plus brain noise, what the model leaves out moves the ratios by at most "
        f"{max(abs(x) for v in chg.values() for x in v):.0f} % against the primary model (dense "
        f"{pct(min(chg['opm_dense']))} to {pct(max(chg['opm_dense']))} %, site-matched {pct(min(chg['opm_matched']))} to "
        f"{pct(max(chg['opm_matched']))} %; the frequency-resolved detector alone, with white noise, gives "
        f"{pct(fr_white['opm_dense'])} and {pct(fr_white['opm_matched'])} %): the cortex 2-4 mm from the "
        f"inner skull as noise, as targets or both {near[0]['opm_dense'][0]:.3f}, {near[1]['opm_dense'][0]:.3f} and "
        f"{near[2]['opm_dense'][0]:.3f} (site-matched {min(x['opm_matched'][0] for x in near):.3f}-"
        f"{max(x['opm_matched'][0] for x in near):.3f}); coloured OPM noise "
        f"with a 1/f corner at {corners} Hz "
        + ", ".join(f"{cv[f'{c:g} Hz']['opm_dense'][0]:.3f}" for c in CORNERS)
        + f" against {cv['none, frequency-resolved']['opm_dense'][0]:.3f} with white noise under the same frequency-resolved "
        "detector (site-matched " + ", ".join(f"{cv[f'{c:g} Hz']['opm_matched'][0]:.3f}" for c in CORNERS)
        + f" against {cv['none, frequency-resolved']['opm_matched'][0]:.3f}); heart and eye fields added at the room "
        f"field's level or filling the magnetometer shortfall {both(room, 3)} and {both(short, 3)}; jointly, the near-skull "
        f"cortex with a {corner_max}-Hz corner "
        f"{both(jw, 3)}, with the heart and eyes added {both(jf, 3)}. (d) Against the OPM white-noise level (solid: sensor plus "
        f"brain noise; dotted: after the {nx}-term projection) the dense ratio falls from {sw['opm_dense/' + IB][0][0]:.2f} at "
        f"{lv[0]:g} fT/√Hz to {sw['opm_dense/' + IB][-1][0]:.2f} at {lv[-1]:g} fT/√Hz "
        f"({sw['opm_dense/' + PR][-1][0]:.2f} after the projection) and the site-matched from "
        f"{sw['opm_matched/' + IB][0][0]:.2f} to {sw['opm_matched/' + IB][-1][0]:.2f}; the ratio reaches 1 (diamonds) at "
        f"{bev('opm_dense/' + IB)} fT/√Hz for the dense array ({bev('opm_dense/' + PR)} after the projection) and "
        f"{bev('opm_matched/' + IB)} for the site-matched ({bev('opm_matched/' + PR)}); in the most adverse combination "
        f"(dashed) the dense ratio is {adv['opm_dense'][0][0]:.2f}, {adv['opm_dense'][lv.index(float(lk))][0]:.2f} and "
        f"{adv['opm_dense'][-1][0]:.2f} at {lv[0]:g}, {lk} and {lv[-1]:g} fT/√Hz. Medians over the targets with 95 % "
        "parcel-bootstrap intervals (lines and bands) within this one anatomy; top and right-hand axes: the same ratios in "
        f"dB (D for the OPM/Neuromag ratios); dense OPM array {sites['opm_dense']} sites, site-matched {sites['opm_matched']}, "
        f"Neuromag ({sq} channels). No OPM noise was measured: the model is checked against measured noise for Neuromag only.")
    return inputs, description, alt, caption


def main():
    style.apply()
    cov, ns, g2 = load()
    pub = published(g2)
    checks = [check_methods_letters()]
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_IN), layout="constrained")
    sf_top, sf_bot = fig.subfigures(2, 1, height_ratios=[1.0, 1.42])
    sf_a, sf_b = sf_top.subfigures(1, 2, width_ratios=[1.0, 1.55])
    sf_c, sf_d = sf_bot.subfigures(1, 2, width_ratios=[1.3, 1.0])
    _, va, ca = panel_a(sf_a, cov)
    _, vb, cb = panel_b(sf_b, cov, pub)
    _, vc, cc = panel_c(sf_c, ns, pub)
    _, vd, cd = panel_d(sf_d, ns, pub)
    h = [Line2D([], [], color=style.ARRAY_COLOR[a], marker=MARKER[a], ms=4.2, lw=1.3,
                mec="white" if a == "opm_dense" else style.ARRAY_COLOR[a], mew=0.4,
                label=f"{ARRAY_TERM[a]} ({g2['arrays'][a]['n_sites']} sites)") for a in ARRAYS]
    h += [Patch(facecolor="0.88", edgecolor="0.45", lw=0.6, label="primary model, 95 % interval shaded (b, c)")]
    fig.legend(handles=h, loc="outside upper center", ncol=3, fontsize=7.3, handlelength=1.8, columnspacing=1.1,
               title="Medians over the targets with their 95 % parcel-bootstrap intervals (lines, bands)", title_fontsize=7.3)
    db = check_db_axes(fig)
    path = style.save(fig, NAME)
    px = plt.imread(path).shape
    size_mm = [round(px[1] / style.DPI * 25.4, 1), round(px[0] / style.DPI * 25.4, 1)]
    if not PAGE_MM - 3.0 <= size_mm[0] <= PAGE_MM + 0.05:
        raise ValueError(f"{path.name} is {size_mm[0]} mm wide, not the page width ({PAGE_MM:g} mm)")
    v = dict(a=va, b=vb, c=vc, d=vd)
    inputs, description, alt, caption = texts(v, cov, ns, g2, pub)
    entry = dict(inputs=inputs, description=description, alt=alt, caption_draft=caption,
                 values=dict(a_neuromag_measured_over_model=va, b_scenarios=vb, c_omissions=vc, d_white_level=vd,
                             published_model={a: list(pub[a]) for a in ARRAYS}),
                 checks=checks + ca + cb + cc + cd, db_axes=db,
                 size=dict(page_width_mm=PAGE_MM, saved_mm=size_mm, saved_px=[px[1], px[0]], dpi=style.DPI))
    style.write_provenance(OUT_JSON, {NAME: entry})
    print(f"results/report/{NAME}.png")
    print(f"results/report/{OUT_JSON}")


if __name__ == "__main__":
    main()
