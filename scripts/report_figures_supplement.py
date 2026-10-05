#!/usr/bin/env python3
"""Supplementary figures of the simulated-spike study, drawn from stored outputs only: nothing is
simulated, re-fitted or re-tested.

  Figure_S_localization_effects     Paired localization-error differences, OPM array minus Neuromag
      (306 channels), with their stored 95 % bootstrap intervals over events, for the nine anatomies,
      both OPM arrays, focal and patch sources at the two strengths and the three inverse methods
      (dipole fit; dSPM peak with MNE-Python and with the study's implementation). Marked from the
      per-event tables: comparisons in which fewer than half of the paired events are detected by
      both systems (their errors are then mostly those of undetected, sub-threshold events), and the
      stored Wilcoxon p (nominal, and below a Bonferroni threshold within anatomy and OPM array).
  Figure_S_joint_detection_localization   Share of the injected events detected by the practical
      detector (thresholds frozen at 1 false event per minute) and share detected and localized
      within 10 mm, per anatomy, system, source, strength and inverse method, with the stored exact
      McNemar p of the paired joint-success endpoints.
  Figure_S_detection_curves_adult   Adult, exploratory run: detection probability against spike
      strength in the four depth bands, oracle and practical detectors, Wilson 95 % bands (event
      level, descriptive: the events of a point share 18 locations).

Inputs: results/g4/g4_localization[_<anatomy>]_summary.json and _events.csv (nine anatomies),
results/g4/g4_adult_summary.json, results/g2/g2_summary.json (the adult's array sizes),
configs/g4_epilepsy.toml (the patch radius; the stored localization configuration is checked against
it) and docs/methods.md (depth bands of section 9, the bootstrap count IC-BOOT-LOC of section 13).
The per-event localization tables are the only link between detection and localization; they are
re-summarised here, and the script stops unless they reproduce every stored summary value it uses
(detection shares, joint shares, discordant counts, median paired differences, shares of events
with the OPM closer). Every number in labels, descriptions and captions is read from these files or
derived from them as stated in figures_supplement.json.

Outputs (results/report/): Figure_S_localization_effects.png, Figure_S_joint_detection_localization.png,
Figure_S_detection_curves_adult.png and figures_supplement.json (inputs, description, alt text,
caption draft, plotted values and the checks of every figure).

Usage: PYTHONPATH=src .venv/bin/python scripts/report_figures_supplement.py
"""
from __future__ import annotations

import json
import math
import re
import sys
import textwrap
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend before pyplot is used)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

from opmsquid import io  # noqa: E402

LOC_SUMMARY = "results/g4/g4_localization{}_summary.json"
LOC_EVENTS = "results/g4/g4_localization{}_events.csv"
DET_ADULT = "results/g4/g4_adult_summary.json"
G2 = "results/g2/g2_summary.json"
CONFIG = "configs/g4_epilepsy.toml"
DOC = "docs/methods.md"
OUT_JSON = "figures_supplement.json"

SQUID = "squid"  # Neuromag, all channels: the primary comparator of the paired localization endpoints
OPM = ("opm_dense", "opm_matched")
SYSTEMS = (SQUID,) + OPM
FAMILIES = ("focal", "patch")
# inverse methods: per-event error column, joint-success share in results[...], paired joint endpoint (None: not
# stored as a paired endpoint), label
INVERSES = (
    dict(key="ecd", err="ecd_error_mm", joint="joint_detect_and_ecd_within_10mm", paired_joint="joint_ecd_10mm",
         label="Dipole fit (ECD)", short="dipole fit"),
    dict(key="dspm_mne", err="dspm_mne_error_mm", joint="joint_detect_and_dspm_mne_within_10mm", paired_joint=None,
         label="dSPM, MNE-Python", short="MNE-Python's dSPM"),
    dict(key="dspm", err="dspm_error_mm", joint="joint_detect_and_dspm_within_10mm", paired_joint="joint_dspm_10mm",
         label="dSPM, study's implementation", short="the study's dSPM"),
)
COLOR = {SQUID: style.ARRAY_COLOR["combined"], "opm_dense": style.ARRAY_COLOR["opm_dense"],
         "opm_matched": style.ARRAY_COLOR["opm_matched"]}
FAMILY_MARKER = {"focal": "o", "patch": "s"}
MINUS = "−"
EN = "–"
DAGGER = "†"
# the five detectors of the detection study, drawn in this order (primary systems last, on top)
DETECTORS = ("squid/grad", "squid/mag", "squid/combined", "opm_matched/opm", "opm_dense/opm")
PRIMARY_DET = ("squid/combined", "opm_matched/opm", "opm_dense/opm")
DASH = {"squid/grad": (0, (3.0, 1.5)), "squid/mag": (0, (1.0, 1.3))}  # the Neuromag channel subsets (no bands)


# ----------------------------------------------------------------------------------------------
# helpers
def rel(path: str, k: str) -> str:
    return path.format("" if k == "adult" else f"_{k}")


def load(path: str):
    return json.loads((ROOT / path).read_text())


def close(a, b, tol: float = 1e-9) -> bool:
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a)), abs(float(b)))


def need(ok: bool, what: str) -> None:
    if not ok:
        raise ValueError(what)


def ic_value(ident: str) -> str:
    """Value cell of an implementation-constant row in docs/methods.md section 13."""
    for line in (ROOT / DOC).read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 4 and cells[0] == ident:
            return cells[2]
    raise KeyError(ident)


def depth_bands() -> list[tuple[float, float]]:
    """The spike study's depth bands as written in docs/methods.md section 9."""
    m = re.search(r"stratified by depth \(([0-9, \-]+) mm\)", (ROOT / DOC).read_text())
    need(m is not None, f"{DOC}: no depth-band statement in section 9")
    return [tuple(float(x) for x in b.split("-")) for b in m.group(1).split(", ")]


def signed(x: float, nd: int = 1) -> str:
    s = f"{x:+.{nd}f}"
    return s.replace("-", MINUS) if s.strip("+-0.") else s.lstrip("+-")


def span(values, fmt: str = "{:.0f}") -> str:
    lo, hi = min(values), max(values)
    a, b = fmt.format(lo), fmt.format(hi)
    return a if a == b else f"{a}{EN}{b}"


def tint(c, f: float = 0.55):
    r = np.array(to_rgb(c))
    return tuple(r + (1 - r) * f)


def slots(keys, gap: float = 0.5) -> dict:
    """Row positions, one unit per anatomy, ``gap`` added between anatomy classes."""
    pos, y, prev = {}, 0.0, None
    for k in keys:
        if prev is not None and style.ANAT_CLASS[k] != prev:
            y += gap
        pos[k], y, prev = y, y + 1.0, style.ANAT_CLASS[k]
    return pos


def row_axes(ax, pos, labels: bool) -> None:
    """Anatomies as rows (top to bottom), labelled on the axes that carries the labels, a hairline
    between classes."""
    keys = list(pos)
    if labels:
        ax.set_yticks([-pos[k] for k in keys],
                      [style.ANAT_LABEL[k] + (f" {DAGGER}" if k == "adult" else "") for k in keys])
    else:
        ax.set_yticks([-pos[k] for k in keys])
    ax.tick_params(axis="y", length=0, labelleft=labels)
    ax.set_ylim(-max(pos.values()) - 0.62, 0.62)
    for a, b in zip(keys[:-1], keys[1:]):
        if style.ANAT_CLASS[a] != style.ANAT_CLASS[b]:
            ax.axhline(-0.5 * (pos[a] + pos[b]), color="0.88", lw=0.7, zorder=0)


def footnote(fig, text: str, y: float, width: int = 132) -> None:
    """Small grey note under the panels, wrapped to the figure's design width."""
    fig.text(0.01, y, textwrap.fill(text, width), ha="left", va="top", fontsize=6.6, color="0.25")


def axis_ticks(lo: float, hi: float, max_labels: int = 7) -> tuple[float, float, float]:
    """(start, end, step): the smallest step of 1, 2, 5, 10, 20 or 25 mm that labels [lo, hi] with at most
    ``max_labels`` ticks, and the axis limits rounded out to it."""
    for step in (1.0, 2.0, 5.0, 10.0, 20.0, 25.0):
        a, b = math.floor(lo / step) * step, math.ceil(hi / step) * step
        if round((b - a) / step) + 1 <= max_labels:
            return a, b, step
    raise ValueError(f"no tick step for {lo}-{hi}")


# ----------------------------------------------------------------------------------------------
# localization: stored summaries and the per-event link between detection and localization
def load_localization() -> dict:
    keys = style.ANAT_ORDER
    L = {k: load(rel(LOC_SUMMARY, k)) for k in keys}
    E = {k: io.read_csv(ROOT / rel(LOC_EVENTS, k)) for k in keys}
    cfg = L["adult"]["config"]
    for k in keys:
        need(L[k]["config"] == cfg, f"{k}: localization configuration differs from the adult's")
    toml = tomllib.loads((ROOT / CONFIG).read_text())
    need(all(toml["localization"].get(c) == v for c, v in cfg.items()),
         f"{CONFIG} [localization] differs from the stored localization configuration")
    strengths = sorted(cfg["strengths_nAm"], reverse=True)  # strongest first (top row)
    # the success radius and the false-event rate are written in the stored key names
    radius = {float(re.search(r"within_(\d+)mm", inv["joint"]).group(1)) for inv in INVERSES}
    radius |= {float(re.search(r"_(\d+)mm$", inv["paired_joint"]).group(1)) for inv in INVERSES if inv["paired_joint"]}
    need(len(radius) == 1, f"joint-success radii differ: {radius}")
    rates = {float(m.group(1)) for k in keys for m in [re.fullmatch(r"thresholds_(\d+(?:\.\d+)?)_per_min", c)
                                                        for c in L[k] if c.startswith("thresholds_") and "per_min" in c] if m}
    need(len(rates) == 1, f"detector operating points differ: {rates}")
    # the dagger notes: the adult at its measured head position, every other head at top contact (stored anatomy notes)
    need("measured head position" in L["adult"]["anatomy"] and all("'top' placement" in L[k]["anatomy"] for k in keys[1:]),
         "head placements differ from the dagger note")
    sites = {k: {v: L[k]["views"][v]["channels"] for v in SYSTEMS} for k in keys}
    need(len({sites[k][SQUID] for k in keys}) == 1, "Neuromag channel counts differ between anatomies")
    heldout = {k: {v: L[k]["thresholds_heldout"]["false_events"][v]["rate_per_min"] for v in SYSTEMS} for k in keys}
    n_cond = {r["n"] for k in keys for r in L[k]["results"].values()} | {r["n"] for k in keys for r in L[k]["paired"].values()}
    need(len(n_cond) == 1, f"localization conditions differ in their event counts: {n_cond}")
    # comparisons per anatomy and OPM array (conditions x endpoints): the Bonferroni family of the text
    conds = [c for c in L["adult"]["paired"] if c.startswith(f"{OPM[0]}_vs_{SQUID}/")]
    endpoints = [inv["err"] for inv in INVERSES] + [inv["paired_joint"] for inv in INVERSES if inv["paired_joint"]]
    need(all(e in L[k]["paired"][c] for k in keys for c in L[k]["paired"] for e in endpoints) and
         all(len(L[k]["paired"][c]) - 1 == len(endpoints) for k in keys for c in L[k]["paired"]),
         "the stored paired entries do not hold exactly the expected endpoints")
    fam_n = len(conds) * len(endpoints)
    patch = {toml["events"]["patch_radius_mm"], load(DET_ADULT)["config"]["events"]["patch_radius_mm"]}
    need(len(patch) == 1, f"patch radius differs between {CONFIG} and {DET_ADULT}: {patch}")
    D = dict(L=L, E=E, keys=keys, cfg=cfg, strengths=strengths, radius=radius.pop(), rate=rates.pop(), sites=sites,
             heldout=heldout, n=n_cond.pop(), family_per_array=fam_n, n_conditions_per_array=len(conds),
             n_endpoints=len(endpoints), patch_radius=patch.pop(), n_boot=int(ic_value("IC-BOOT-LOC")))
    D["link"] = link(D)
    return D


def link(D) -> dict:
    """Per anatomy, OPM array, source and strength, from the per-event table: which of the paired
    events each system detects, the tied dSPM errors, and the stored summaries recomputed (each must
    agree with the stored value)."""
    out, checks = {}, 0
    for k in D["keys"]:
        rows, L = D["E"][k], D["L"][k]
        need(len(rows) == L["n_events"] * len(L["views"]), f"{k}: event table has {len(rows)} rows")
        by = {(r["array"], int(r["event"])): r for r in rows}
        need(len(by) == len(rows), f"{k}: duplicate (array, event) rows")
        for fam in FAMILIES:
            for s in D["strengths"]:
                ev = sorted({int(r["event"]) for r in rows if r["family"] == fam and float(r["strength_nAm"]) == s})
                need(len(ev) == D["n"], f"{k} {fam} {s:g}: {len(ev)} events")
                det = {a: np.array([by[a, e]["detected"] == "True" for e in ev]) for a in SYSTEMS}
                need(all(by[a, e]["detected"] in ("True", "False") for a in SYSTEMS for e in ev), f"{k}: detected flags")
                err = {(a, inv["err"]): np.array([float(by[a, e][inv["err"]]) for e in ev]) for a in SYSTEMS for inv in INVERSES}
                ok = {(a, inv["key"]): det[a] & (err[a, inv["err"]] <= D["radius"]) for a in SYSTEMS for inv in INVERSES}
                cond = f"{fam}/{s:g}nAm"
                for a in SYSTEMS:
                    r = L["results"][f"{a}/{cond}"]
                    need(close(r["detected"], det[a].mean()), f"{k} {a} {cond}: detected share")
                    for inv in INVERSES:
                        need(close(r[inv["joint"]], ok[a, inv["key"]].mean()), f"{k} {a} {cond}: {inv['joint']}")
                    checks += 1 + len(INVERSES)
                for a in OPM:
                    p = L["paired"][f"{a}_vs_{SQUID}/{cond}"]
                    rec = dict(n=len(ev), both=int(np.sum(det[a] & det[SQUID])), opm_only=int(np.sum(det[a] & ~det[SQUID])),
                               squid_only=int(np.sum(~det[a] & det[SQUID])), neither=int(np.sum(~det[a] & ~det[SQUID])), ties={})
                    for inv in INVERSES:
                        diff = err[a, inv["err"]] - err[SQUID, inv["err"]]
                        need(np.all(np.isfinite(diff)), f"{k} {a} {cond}: non-finite {inv['err']}")
                        st = p[inv["err"]]
                        need(close(st["median_difference"], np.median(diff)), f"{k} {a} {cond}: median {inv['err']}")
                        need(close(st["share_opm_smaller"], np.mean(diff < 0)), f"{k} {a} {cond}: share {inv['err']}")
                        rec["ties"][inv["key"]] = int(np.sum(diff == 0))
                        checks += 2
                        if inv["paired_joint"]:
                            pj = p[inv["paired_joint"]]
                            oa, os_ = ok[a, inv["key"]], ok[SQUID, inv["key"]]
                            need(pj["only_opm"] == int(np.sum(oa & ~os_)) and pj["only_squid"] == int(np.sum(~oa & os_))
                                 and close(pj["opm"], oa.mean()) and close(pj["squid"], os_.mean()),
                                 f"{k} {a} {cond}: {inv['paired_joint']}")
                            checks += 1
                    out[k, a, fam, s] = rec
    D["n_checks"] = checks
    return out


# ----------------------------------------------------------------------------------------------
# Figure S: localization effect sizes
OFFSETS = (("opm_dense", "focal", 0.3), ("opm_dense", "patch", 0.1), ("opm_matched", "focal", -0.1),
           ("opm_matched", "patch", -0.3))


def stars(p: float, bonf: float) -> str:
    return "**" if p < bonf else ("*" if p < 0.05 else "")


def site_text(D, a: str) -> str:
    v = [D["sites"][k][a] for k in D["keys"] if k != "adult"]
    return f"{D['sites']['adult'][a]} sites in the adult, {span(v)} in the smaller heads"


def place_marks(ax, marks, markers, x0: float, x1: float) -> None:
    """Write each p mark just left of its interval, on the estimate's row, or right of it when a marker of a
    neighbouring row would sit under the mark (the asterisk glyph sits high, so the text is lowered a little)."""
    w = x1 - x0
    for y, lo, hi, mark in marks:
        width = 0.033 * w * len(mark)  # about 3.3 % of the panel per asterisk at this size and panel width
        left = (lo - 0.01 * w - width, lo - 0.01 * w)
        hit = any(abs(my - y) < 0.25 and left[0] - 0.02 * w <= mx <= left[1] + 0.02 * w for mx, my in markers)
        if hit and hi + 0.01 * w + width <= x1:
            ax.text(hi + 0.01 * w, y - 0.07, mark, ha="left", va="center", fontsize=6.6, color="0.1", zorder=4)
        else:
            ax.text(left[1], y - 0.07, mark, ha="right", va="center", fontsize=6.6, color="0.1", zorder=4)


def fig_effects(D) -> dict:
    keys, L, strengths, n = D["keys"], D["L"], D["strengths"], D["n"]
    pos = slots(keys)
    bonf = 0.05 / D["family_per_array"]
    half = n / 2
    fig = plt.figure(figsize=(style.FULL_W, 8.9))
    gs = fig.add_gridspec(2, 3, left=0.25, right=0.985, top=0.865, bottom=0.165, hspace=0.30, wspace=0.2)
    values, plotted = {}, []
    # x limits per strength: one for the dipole, one shared by the two dSPM implementations (room for the stars)
    lims = {}
    for s in strengths:
        for grp, invs in (("ecd", ("ecd",)), ("dspm", ("dspm_mne", "dspm"))):
            cis = [L[k]["paired"][f"{a}_vs_{SQUID}/{f}/{s:g}nAm"][inv["err"]]["ci95"] for k in keys for a in OPM for f in FAMILIES
                   for inv in INVERSES if inv["key"] in invs]
            lo, hi = min(c[0] for c in cis), max(c[1] for c in cis)
            # ticks from the data range; the left limit leaves room for the p marks left of the intervals
            t0, t1, step = axis_ticks(lo, hi)
            lims[s, grp] = (min(t0, lo - 0.09 * (hi - lo)), t1, step)
    letters = iter("abcdef")
    for row, s in enumerate(strengths):
        for col, inv in enumerate(INVERSES):
            ax = fig.add_subplot(gs[row, col])
            grp = "ecd" if inv["key"] == "ecd" else "dspm"
            x0, x1, step = lims[s, grp]
            marks = []  # (y, lo, hi, mark), placed after every marker is known
            for k in keys:
                for a, fam, dy in OFFSETS:
                    st = L[k]["paired"][f"{a}_vs_{SQUID}/{fam}/{s:g}nAm"][inv["err"]]
                    lk = D["link"][k, a, fam, s]
                    y = -pos[k] + dy
                    lo, hi = st["ci95"]
                    need(x0 + 0.06 * (x1 - x0) <= lo and hi <= x1, f"interval outside the axis: {k} {a} {fam} {s:g} {inv['key']}")
                    c = COLOR[a]
                    full = lk["both"] >= half
                    ax.plot([lo, hi], [y, y], color=c, lw=1.0, solid_capstyle="butt", zorder=2)
                    ax.plot(st["median_difference"], y, marker=FAMILY_MARKER[fam], ms=3.9 if fam == "focal" else 3.5, mec=c,
                            mfc=c if full else "white", mew=0.9, ls="none", zorder=3)
                    mark = stars(st["wilcoxon_p"], bonf)
                    if mark:
                        marks.append((y, lo, hi, mark))
                    tag = f"{k}/{a}/{fam}/{s:g}nAm/{inv['key']}"
                    values[tag] = dict(median_difference_mm=st["median_difference"], ci95_mm=st["ci95"], wilcoxon_p=st["wilcoxon_p"],
                                       share_opm_closer=st["share_opm_smaller"], n_events=lk["n"],
                                       n_tied=lk["ties"][inv["key"]], detected_both=lk["both"], detected_opm_only=lk["opm_only"],
                                       detected_neuromag_only=lk["squid_only"], detected_neither=lk["neither"],
                                       drawn_filled=bool(full), mark=mark)
                    plotted.append(dict(tag=tag, s=s, inv=inv["key"], k=k, a=a, fam=fam, full=full, **values[tag]))
            offset = {(a, f): dy for a, f, dy in OFFSETS}
            place_marks(ax, marks, [(p["median_difference_mm"], -pos[p["k"]] + offset[p["a"], p["fam"]])
                                    for p in plotted if p["s"] == s and p["inv"] == inv["key"]], x0, x1)
            ax.axvline(0.0, color="0.5", lw=0.8, zorder=1)
            row_axes(ax, pos, labels=col == 0)
            ax.set_xlim(x0, x1)
            ticks = np.arange(math.ceil(x0 / step) * step, x1 + 1e-9, step)
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xticklabels([signed(t, 0) if t else "0" for t in ticks], fontsize=7.0)
            ax.set_title(f"({next(letters)}) {inv['label']}", loc="center", fontsize=7.9)
            if row == 1:
                ax.set_xlabel(f"OPM {MINUS} Neuromag (mm)", fontsize=7.8)
    # row headers: strength and how many of the paired events both systems detect
    for row, s in enumerate(strengths):
        both = [D["link"][k, a, f, s]["both"] for k in keys for a in OPM for f in FAMILIES]
        y = 0.873 if row == 0 else 0.497
        extra = "; the errors are mostly those of undetected events" if max(both) < half else ""
        fig.text(0.01, y + 0.02, f"{s:g} nAm: both systems detect {span(both)} of the {n} paired events per condition{extra}",
                 ha="left", va="bottom", fontsize=8.2, color="0.1")
    fig.text(0.6175, 0.118, f"Median paired difference in localization error, OPM {MINUS} Neuromag (mm; negative: OPM closer)",
             ha="center", va="top", fontsize=7.8)
    h = [Line2D([], [], color=COLOR[a], lw=1.0, marker="o", ms=3.9, mfc=COLOR[a], mec=COLOR[a],
                label=f"{'dense' if a == 'opm_dense' else 'site-matched'} OPM array ({site_text(D, a)})") for a in OPM]
    h2 = [Line2D([], [], color="0.3", ls="none", marker="o", ms=3.9, label="focal dipole"),
          Line2D([], [], color="0.3", ls="none", marker="s", ms=3.5, label=f"{D['patch_radius']:g}-mm patch"),
          Line2D([], [], color="0.3", ls="none", marker="o", ms=3.9, mfc="0.3",
                 label=f"filled: both systems detect ≥ {half:g} of {n}"),
          Line2D([], [], color="0.3", ls="none", marker="o", ms=3.9, mfc="white", label=f"open: fewer than {half:g}")]
    leg = fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.005, 0.998), ncol=1, handlelength=2.0, fontsize=7.2)
    fig.legend(handles=h2, loc="upper left", bbox_to_anchor=(0.005, 0.955), ncol=4, handlelength=1.0, columnspacing=1.3, fontsize=7.2)
    fig.add_artist(leg)
    rate = f"{D['rate']:g} false event per minute"
    footnote(fig, (
        f"Each estimate: the median over the {n} events of one condition (one per location; both systems see the same event, "
        f"noise and coregistration draw) of the OPM array's localization error minus Neuromag's, with its 95 % bootstrap interval "
        f"over events ({D['n_boot']:,} resamples). Errors are read at the true peak sample for every event, detected or not, so the "
        f"comparison is not selected by detection; detection (practical detector, thresholds frozen at {rate}) only sets the marker "
        f"fill. dSPM peaks lie on a {D['cfg']['grid_spacing_mm']:g}-mm source grid, so paired errors often tie: medians and interval "
        "ends at exactly 0 are common. Beside an interval (left, or right where a marker is in the way): * p < 0.05, "
        f"** p < {bonf:.4g} "
        f"(Wilcoxon signed-rank, uncorrected; {bonf:.4g} = 0.05 / {D['family_per_array']}, the comparisons of one anatomy and OPM "
        f"array). {DAGGER} Adult at its measured head position; the other heads at top contact in the adult helmet."), y=0.098)
    style.save(fig, "Figure_S_localization_effects")
    return dict(values=values, plotted=plotted, lims={f"{s:g}nAm/{g}": list(v[:2]) for (s, g), v in lims.items()}, bonf=bonf)


# ----------------------------------------------------------------------------------------------
# Figure S: joint detection and localization
SYS_OFFSET = ((SQUID, 0.28), ("opm_dense", 0.0), ("opm_matched", -0.28))


def fig_joint(D) -> dict:
    keys, L, strengths, n = D["keys"], D["L"], D["strengths"], D["n"]
    pos = slots(keys)
    conds = [(f, s) for s in strengths for f in FAMILIES]
    fig, axs = plt.subplots(len(INVERSES), len(conds), figsize=(style.FULL_W, 9.0), sharex=True, sharey=True,
                            gridspec_kw=dict(left=0.25, right=0.985, top=0.885, bottom=0.125, hspace=0.24, wspace=0.14))
    values = {}
    dark = {SQUID: "0.12", "opm_dense": COLOR["opm_dense"], "opm_matched": COLOR["opm_matched"]}
    light = {SQUID: "0.86", "opm_dense": tint(COLOR["opm_dense"], 0.78), "opm_matched": tint(COLOR["opm_matched"], 0.72)}
    edge = {SQUID: "0.35", "opm_dense": COLOR["opm_dense"], "opm_matched": COLOR["opm_matched"]}
    for i, inv in enumerate(INVERSES):
        for j, (fam, s) in enumerate(conds):
            ax = axs[i, j]
            cond = f"{fam}/{s:g}nAm"
            for k in keys:
                for a, dy in SYS_OFFSET:
                    r = L[k]["results"][f"{a}/{cond}"]
                    det, joint = 100 * r["detected"], 100 * r[inv["joint"]]
                    y = -pos[k] + dy
                    ax.barh(y, det, height=0.25, color=light[a], edgecolor=edge[a], linewidth=0.45, zorder=2)
                    ax.barh(y, joint, height=0.25, color=dark[a], linewidth=0, zorder=3)
                    v = dict(detected_pct=det, detected_and_within_pct=joint)
                    if a != SQUID and inv["paired_joint"]:
                        pj = L[k]["paired"][f"{a}_vs_{SQUID}/{cond}"][inv["paired_joint"]]
                        v.update(only_opm=pj["only_opm"], only_neuromag=pj["only_squid"], mcnemar_exact_p=pj["mcnemar_exact_p"])
                        if pj["mcnemar_exact_p"] < 0.05:
                            ax.text(min(det + 1.5, 101.5), y, "*", ha="left", va="center", fontsize=7.0, color="0.1", zorder=4)
                    values[f"{k}/{a}/{cond}/{inv['key']}"] = v
            row_axes(ax, pos, labels=j == 0)
            ax.set_xlim(0, 100)
            ax.xaxis.set_major_locator(FixedLocator([0, 50, 100]))
            ax.xaxis.set_minor_locator(FixedLocator([25, 75]))
            ax.grid(axis="x", which="both", color="0.92", lw=0.6, zorder=0)
            ax.set_axisbelow(True)
            if i == len(INVERSES) - 1:  # end labels inside the panel, so neighbouring panels' 100 and 0 stay apart
                labels = ax.get_xticklabels()
                labels[0].set_ha("left")
                labels[-1].set_ha("right")
            if i == 0:
                source = "Focal dipole" if fam == "focal" else f"{D['patch_radius']:g}-mm patch"
                ax.set_title(f"{source}, {s:g} nAm", fontsize=8.0)
        y_top = axs[i, 0].get_position().y1
        fig.text(0.01, y_top + (0.028 if i == 0 else 0.006), f"({'abc'[i]})  {inv['label']}: detected and localized within "
                 f"{D['radius']:g} mm", ha="left", va="bottom", fontsize=8.2)
    fig.text(0.6175, 0.095, f"Share of the {n} injected events per condition (%)", ha="center", va="top", fontsize=7.8)
    systems = [Patch(facecolor=dark[a], edgecolor="none", label=lab) for a, lab in
               ((SQUID, f"Neuromag ({D['sites']['adult'][SQUID]} channels)"), ("opm_dense", "dense OPM array"),
                ("opm_matched", "site-matched OPM array"))]
    parts = [Patch(facecolor="0.12", edgecolor="none", label=f"dark: detected and localized within {D['radius']:g} mm"),
             Patch(facecolor="0.86", edgecolor="0.35", linewidth=0.45, label=f"pale: detected, error > {D['radius']:g} mm"),
             Patch(facecolor="white", edgecolor="none", label="bar end: share detected")]
    # legends fill column by column: interleave so that the systems form the first row
    fig.legend(handles=[x for pair in zip(systems, parts) for x in pair], loc="upper left", bbox_to_anchor=(0.005, 0.998),
               ncol=3, handlelength=1.4, columnspacing=1.5, fontsize=7.2)
    held = [D["heldout"][k][a] for k in keys for a in SYSTEMS]
    footnote(fig, (
        "Detection: the practical detector, idealized (templates: the simulated spike waveform; candidate sources from the true "
        f"forward model), thresholds set for {D['rate']:g} false event per minute on {D['cfg']['calibration_min']:g} min of null "
        f"data and frozen ({span(held, '{:.1f}')} per minute on {D['cfg']['holdout_min']:g} min of held-out null data, nine "
        "anatomies and three systems). Localization is read at the true peak sample, so joint success assumes the event time is "
        f"known. Dense OPM array: {site_text(D, 'opm_dense')}; site-matched: {D['sites']['adult']['opm_matched']} and "
        f"{span([D['sites'][k]['opm_matched'] for k in D['keys'][1:]])}. * paired OPM {MINUS} Neuromag difference in joint success "
        "with p < 0.05 (exact McNemar, uncorrected; stored for the dipole fit and the study's dSPM). "
        f"{DAGGER} Adult at its measured head position; the other heads at top contact in the adult helmet."), y=0.072)
    style.save(fig, "Figure_S_joint_detection_localization")
    return dict(values=values, conds=[f"{f}/{s:g}nAm" for f, s in conds], heldout_range=[min(held), max(held)])


# ----------------------------------------------------------------------------------------------
# Figure S: the adult's detection curves (exploratory run)
def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def fig_curves() -> dict:
    g, g2 = load(DET_ADULT), load(G2)
    cfg = g["config"]
    need(g["anatomy"].startswith("MNE sample subject") and "measured head position" in g["anatomy"],
         f"{DET_ADULT}: the adult is not at its measured head position")
    strengths = cfg["events"]["strengths_nAm"]
    bands = depth_bands()
    need(len(bands) == len({L["stratum"][0] for L in g["locations"]}), "depth bands: count differs from the stored strata")
    for L in g["locations"]:
        lo, hi = bands[L["stratum"][0]]
        need(lo <= L["depth_mm"] < hi, f"location at {L['depth_mm']:.1f} mm outside band {lo:g}-{hi:g}")
    n_loc = {b: sum(L["stratum"][0] == b for L in g["locations"]) for b in range(len(bands))}
    need(len(set(n_loc.values())) == 1, f"locations per band differ: {n_loc}")
    n_loc = n_loc[0]
    need(sorted(cfg["detector"]["arrays"]) == sorted(DETECTORS), "detector set differs from the stored one")
    op = cfg["detector"]["operating_points_per_min"][0]
    mode_p = f"practical@{op:g}"
    arr = g2["arrays"]
    n_ch = dict(combined=arr["squid"]["channels"], mag=arr["squid"]["sites"], grad=arr["squid"]["channels"] - arr["squid"]["sites"],
                dense=arr["opm_dense"]["n_sites"], matched=arr["opm_matched"]["n_sites"])
    loc_views = load(rel(LOC_SUMMARY, "adult"))["views"]
    need((loc_views["squid"]["channels"], loc_views["squid_mag"]["channels"], loc_views["squid_grad"]["channels"],
          loc_views["opm_dense"]["channels"], loc_views["opm_matched"]["channels"])
         == (n_ch["combined"], n_ch["mag"], n_ch["grad"], n_ch["dense"], n_ch["matched"]),
         "the adult's array sizes differ between the G2 summary and the spike study")
    label = {"squid/combined": f"Neuromag ({n_ch['combined']} channels)",
             "squid/grad": f"Neuromag gradiometers only ({n_ch['grad']})",
             "squid/mag": f"Neuromag magnetometers only ({n_ch['mag']})",
             "opm_matched/opm": f"site-matched OPM array ({n_ch['matched']} sites)",
             "opm_dense/opm": f"dense OPM array ({n_ch['dense']} sites)"}
    color = {"squid/combined": COLOR[SQUID], "squid/grad": "0.45", "squid/mag": "0.62", "opm_matched/opm": COLOR["opm_matched"],
             "opm_dense/opm": COLOR["opm_dense"]}
    fig, axs = plt.subplots(2, len(bands), figsize=(style.FULL_W, 5.35), sharex=True, sharey=True,
                            gridspec_kw=dict(left=0.085, right=0.99, top=0.80, bottom=0.215, hspace=0.30, wspace=0.08))
    values, n_pts, n_checked = {}, set(), 0
    for row, mode in enumerate(("oracle", mode_p)):
        for b, (lo, hi) in enumerate(bands):
            ax = axs[row, b]
            for key in DETECTORS:
                det = g["detectors"][key]
                pts = [det["curves"][f"{mode}/focal/depth{b}/{s:g}nAm"] for s in strengths]
                p = np.array([c["p"] for c in pts])
                for c in pts:
                    k_ = round(c["p"] * c["n"])
                    need(close(k_ / c["n"], c["p"]) and all(close(x, y) for x, y in zip(wilson(k_, c["n"]), c["ci"])),
                         f"{key} {mode} band {b}: stored Wilson interval not reproduced")
                    n_pts.add(c["n"])
                    n_checked += 1
                prim = key in PRIMARY_DET
                if prim:
                    ax.fill_between(strengths, [c["ci"][0] for c in pts], [c["ci"][1] for c in pts], color=color[key],
                                    alpha=0.13 if key != "squid/combined" else 0.09, lw=0, zorder=1)
                ax.plot(strengths, p, ls="-" if prim else DASH[key], color=color[key], lw=1.35 if prim else 0.9,
                        marker="o" if prim else None, ms=2.6, zorder=3 if prim else 2)
                values[f"{mode}/{key}/depth{b}"] = dict(n=[c["n"] for c in pts], p=[c["p"] for c in pts],
                                                       wilson95=[c["ci"] for c in pts],
                                                       s50_nAm=det["strength_for_50pct_nAm"][f"{mode}/depth{b}"])
            ax.axhline(0.5, color="0.82", lw=0.7, zorder=0)
            ax.set_xscale("log")
            ax.set_xlim(strengths[0] / 1.25, strengths[-1] * 1.25)
            ax.xaxis.set_major_locator(FixedLocator(strengths))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xticklabels([f"{s:g}" for s in strengths])
            ax.set_ylim(-0.02, 1.02)
            ax.yaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1.0]))
            ax.set_yticklabels(["0", "", "0.5", "", "1"])
            if row == 0:
                ax.set_title(f"{lo:g}{EN}{hi:g} mm below the scalp", fontsize=8.2)
    need(len(n_pts) == 1, f"events per point differ: {n_pts}")
    n_pt = n_pts.pop()
    stretches = cfg["events"]["stretches"]
    need(n_pt == n_loc * len(stretches), "events per point differ from locations x morphologies")
    axs[0, 0].set_ylabel("Detection probability")
    axs[1, 0].set_ylabel("Detection probability")
    pos0, pos1 = axs[0, 0].get_position(), axs[1, 0].get_position()
    fig.text(0.01, pos0.y1 + 0.05, "(a)  Oracle detector: knows source, waveform and peak time; false-positive probability "
             f"{cfg['detector']['oracle_alpha']:g} per trial", ha="left", va="bottom", fontsize=8.2)
    fig.text(0.01, pos1.y1 + 0.012, f"(b)  Practical detector: scans all times, {len(stretches)} templates, {g['n_dictionary']} "
             f"candidate sources; {op:g} false event per minute (frozen)", ha="left", va="bottom", fontsize=8.2)
    fig.text(0.54, 0.155, f"Focal dipole strength (nAm, log scale; {len(stretches)} spike-wave morphologies pooled)", ha="center",
             va="top", fontsize=7.8)
    h = [Line2D([], [], color=color[k], lw=1.35 if k in PRIMARY_DET else 0.9, ls="-" if k in PRIMARY_DET else DASH[k],
                marker="o" if k in PRIMARY_DET else None, ms=2.6, label=label[k])
         for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm", "squid/grad", "squid/mag")]
    h.append(Patch(facecolor="0.5", alpha=0.25, edgecolor="none", label="Wilson 95 % band (event level, descriptive)"))
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.005, 0.998), ncol=3, handlelength=2.2, columnspacing=1.3, fontsize=7.2)
    held = {k: g["detectors"][k]["heldout_false_per_min"][f"{op:g}"] for k in DETECTORS}
    names = {"squid/combined": f"Neuromag {n_ch['combined']}", "squid/grad": "gradiometers", "squid/mag": "magnetometers",
             "opm_dense/opm": "dense", "opm_matched/opm": "site-matched"}
    footnote(fig, (
        f"Adult (MNE sample subject at its measured head position), exploratory spike run: {n_loc} locations per depth band, "
        f"each with {len(stretches)} morphologies at every strength, so {n_pt} focal events per point; the noise (sensor, cortical "
        f"background and room field) is shared by all systems in each simulated segment. Bands: Wilson 95 % intervals that treat the "
        f"{n_pt} events of a point as independent; the events share {n_loc} locations, so the bands are descriptive, not "
        f"inferential (the location is the statistical unit of the paired comparisons). The practical detector is idealized: its "
        "templates are the simulated morphologies and its candidate sources come from the true forward model. Its thresholds were "
        "set on "
        f"{cfg['null']['calibration_min']:g} min of null data; on {cfg['null']['heldout_min']:g} min of held-out null data they gave "
        + ", ".join(f"{held[k]:g} ({names[k]})" for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm", "squid/grad",
                                                           "squid/mag"))
        + " false events per minute. Grey horizontal line: 50 % detection."), y=0.105)
    style.save(fig, "Figure_S_detection_curves_adult")
    return dict(values=values, bands=bands, n_locations_per_band=n_loc, n_events_per_point=n_pt, heldout_per_min=held,
                n_wilson_checked=n_checked, labels=label, mode_practical=mode_p, strengths=strengths, n_channels=n_ch,
                oracle_alpha=cfg["detector"]["oracle_alpha"], n_dictionary=g["n_dictionary"], stretches=stretches,
                calibration_min=cfg["null"]["calibration_min"], heldout_min=cfg["null"]["heldout_min"],
                simulated_at=g.get("simulated_at_commit"), seed=cfg["simulation"]["seed"])


# ----------------------------------------------------------------------------------------------
# provenance, alt text and caption drafts (every number from the plotted values)
def entry_effects(D, r) -> dict:
    keys, files = D["keys"], [rel(LOC_SUMMARY, k) for k in D["keys"]]
    P = r["plotted"]
    hi_s, lo_s = D["strengths"][0], D["strengths"][-1]

    def sel(**kw):
        return [p for p in P if all(p[a] == b for a, b in kw.items())]

    def rng(ps, f="median_difference_mm"):
        return f"{signed(min(p[f] for p in ps))} to {signed(max(p[f] for p in ps))} mm"

    short = {inv["key"]: inv["short"] for inv in INVERSES}

    n_all = len(P)
    nom = [p for p in P if p["wilcoxon_p"] < 0.05]
    bon = [p for p in P if p["wilcoxon_p"] < r["bonf"]]
    open_ = [p for p in P if not p["full"]]
    nom_open = [p for p in nom if not p["full"]]
    zero_med = [p for p in P if p["median_difference_mm"] == 0]
    pos_med = [p for p in P if p["median_difference_mm"] > 0]
    ecd_hi = sel(s=hi_s, inv="ecd")
    dspm_hi = [p for p in P if p["s"] == hi_s and p["inv"] != "ecd"]
    dspm_lo = [p for p in P if p["s"] == lo_s and p["inv"] != "ecd"]
    ecd_lo = sel(s=lo_s, inv="ecd")
    ecd_ci_hi = max(max(abs(p["ci95_mm"][0]), abs(p["ci95_mm"][1])) for p in ecd_hi)
    by_inv = {inv["key"]: len([p for p in nom if p["inv"] == inv["key"]]) for inv in INVERSES}
    bon_txt = "; ".join(f"{style.ANAT_LABEL[p['k']]}, {'dense' if p['a'] == 'opm_dense' else 'site-matched'} array, "
                        f"{'focal' if p['fam'] == 'focal' else 'patch'} {p['s']:g} nAm, {short[p['inv']]} "
                        f"{signed(p['median_difference_mm'])} mm (p = {p['wilcoxon_p']:.2g})" for p in bon)
    need(all(p["median_difference_mm"] < 0 for p in bon), "a Bonferroni survivor favours Neuromag (caption wording)")
    need(all(p["inv"] != "ecd" for p in zero_med), "an exactly zero dipole-fit median (caption wording: dSPM ties)")
    open_lo = [p for p in open_ if p["s"] == lo_s]
    layers = [int(x) for x in re.findall(r"(\d+)-layer", D["cfg"]["inverse_bem"])]
    need(len(layers) == 2, f"inverse BEM description not understood: {D['cfg']['inverse_bem']}")
    inputs = ([f"{f} :: paired['<opm_dense|opm_matched>_vs_squid/<focal|patch>/<80|320>nAm'].<ecd_error_mm|dspm_mne_error_mm|"
               "dspm_error_mm> (median_difference, ci95, wilcoxon_p, share_opm_smaller, n)" for f in files]
              + [f"{rel(LOC_EVENTS, k)} :: detected, ecd_error_mm, dspm_mne_error_mm, dspm_error_mm per (array, event) "
                 "(detection overlap and ties; the stored medians and shares recomputed and checked)" for k in keys]
              + [f"{files[0]} :: config (grid_spacing_mm, strengths_nAm; identical in all nine summaries and equal to "
                 f"{CONFIG} [localization]), views[*].channels, thresholds_1_per_min (operating point)",
                 f"{CONFIG} :: events.patch_radius_mm",
                 f"{DOC} :: section 13 IC-BOOT-LOC (event resamples of every paired median error difference)"])
    description = (
        "Exploratory localization study of the spike simulations (one event per location and condition: "
        f"{D['cfg']['n_locations']} locations, focal dipoles and {D['patch_radius']:g}-mm patches at "
        f"{' and '.join(f'{s:g}' for s in sorted(D['strengths']))} nAm; inverse with a {layers[0]}-layer BEM (the simulation "
        f"used {layers[1]} layers) and a {D['cfg']['coreg_shift_mm']:g}-mm / {D['cfg']['coreg_angle_deg']:g}-degree "
        "coregistration error shared by all systems). Each point is the stored paired statistic of g4_localization*_summary.json "
        "paired[...]: the median over the condition's events of the OPM array's error minus Neuromag's (all channels), with the "
        f"stored 95 % bootstrap interval over events ({D['n_boot']:,} resamples, docs/methods.md IC-BOOT-LOC) and the uncorrected "
        "Wilcoxon signed-rank p. Panels: rows = strength, columns = inverse method (equivalent current dipole; dSPM peak through "
        "MNE-Python; dSPM peak through the study's implementation); rows within a panel = the nine anatomies (style.ANAT_ORDER), "
        "four estimates each (dense array above the site-matched one; circle focal, square patch). The errors are read at the "
        "true peak sample for every event, detected or not; the marker fill is derived from the per-event tables: filled when "
        f"both systems detected at least half of the {D['n']} paired events (practical detector, thresholds frozen at "
        f"{D['rate']:g} false event per minute), open otherwise. Stars (beside the interval: left, or right where a "
        "neighbouring marker is in the way): * p < 0.05, ** p < "
        f"{r['bonf']:.4g} (0.05 / "
        f"{D['family_per_array']}: {D['n_conditions_per_array']} conditions x {D['n_endpoints']} endpoints per anatomy and OPM "
        "array, the within-anatomy-and-array Bonferroni family of the text). The x axis differs between strengths and between "
        "the dipole and the dSPM columns (limits: " + ", ".join(f"{k} {v[0]:.1f} to {v[1]:.1f} mm" for k, v in r["lims"].items())
        + "). Checks: the per-event tables reproduce every stored detection share, joint share, discordant count, median "
        f"difference and share with the OPM closer that the figure uses ({D['n_checks']} values). Selection and censoring: the "
        "paired statistics drawn include every event (no selection by detection), whereas the stored per-system medians among "
        "detected events (results[...].*_median_detected) are selected by each system's own detections and are not drawn; every "
        "event has a dipole and a dSPM error (docs/methods.md section 9: the dipole fit returned a dipole for every event), so no "
        "error is censored and gross errors are included. Neuromag's magnetometers or gradiometers alone (secondary comparators, "
        "paired_secondary) are not drawn. The localization summaries "
        "hold no depth-band split (6 locations per band would be too few), so depth is not drawn.")
    alt = (f"Six forest plots (two strengths by three inverse methods) of the median paired localization-error difference, OPM "
           f"array minus Neuromag, for nine anatomies and two OPM arrays. At {hi_s:g} nAm the dipole-fit differences range from "
           f"{rng(ecd_hi)} with every interval within {ecd_ci_hi:.0f} mm of zero, and the dSPM differences from {rng(dspm_hi)}. At "
           f"{lo_s:g} nAm, where most events go undetected (open markers), the dipole differences range from {rng(ecd_lo)} and the "
           f"dSPM differences from {rng(dspm_lo)}, with wide intervals.")
    lo_text = (f"all {len(open_lo)} estimates rest mostly on undetected events" if len(open_lo) == len(sel(s=lo_s))
               else f"{len(open_lo)} of {len(sel(s=lo_s))} estimates rest mostly on undetected events")
    caption = (
        f"Localization effect sizes, OPM array minus Neuromag ({D['sites']['adult'][SQUID]} channels), paired on "
        f"identical simulated spikes (median over the {D['n']} events of a condition, 95 % bootstrap interval over events; "
        f"negative: the OPM array localizes closer to the true source). Dark blue: dense OPM array; light blue: site-matched OPM "
        f"array; circles: focal dipoles; squares: {D['patch_radius']:g}-mm patches; open markers: fewer than half of the events "
        f"detected by both systems. (a{EN}c) {hi_s:g} nAm: dipole-fit differences {rng(ecd_hi)}; dSPM differences "
        f"{rng(dspm_hi)}. (d{EN}f) {lo_s:g} nAm: {lo_text}; dipole {rng(ecd_lo)}, dSPM {rng(dspm_lo)}. Of the "
        f"{n_all} error comparisons, {len(nom)} have p < 0.05 (Wilcoxon, uncorrected; dipole fit {by_inv['ecd']}, MNE-Python's "
        f"dSPM {by_inv['dspm_mne']}, the study's dSPM {by_inv['dspm']}), {len(nom_open)} of them with fewer than half "
        f"of the events detected by both systems, and {len(bon)} fall below {r['bonf']:.4g} (Bonferroni within anatomy and OPM "
        f"array): {bon_txt}. {len(zero_med)} medians are exactly zero (tied dSPM peaks on the shared grid) and {len(pos_med)} "
        "are positive. The errors include every event, detected or not; localization is read at the true peak sample.")
    values = dict(order=keys, row_offsets={f"{a}/{f}": dy for a, f, dy in OFFSETS}, x_limits_mm=r["lims"],
                  bonferroni_threshold=r["bonf"], family_per_anatomy_and_array=D["family_per_array"],
                  n_estimates=n_all, n_p_below_005=len(nom), n_p_below_bonferroni=len(bon), n_open=len(open_),
                  n_open_by_strength={f"{s:g}nAm": len([p for p in open_ if p["s"] == s]) for s in D["strengths"]},
                  n_p_below_005_by_method=by_inv, n_zero_median=len(zero_med), n_positive_median=len(pos_med),
                  estimates=r["values"])
    return dict(inputs=inputs, description=description, alt=alt, caption_draft=caption, values=values)


def entry_joint(D, r) -> dict:
    keys = D["keys"]
    files = [rel(LOC_SUMMARY, k) for k in keys]
    V = r["values"]
    hi_s, lo_s = D["strengths"][0], D["strengths"][-1]

    def vals(a=None, s=None, inv=None, f="detected_and_within_pct"):
        out = []
        for k in keys:
            for aa in SYSTEMS:
                for fam in FAMILIES:
                    for ss in D["strengths"]:
                        for ii in INVERSES:
                            if (a is None or aa == a) and (s is None or ss == s) and (inv is None or ii["key"] == inv):
                                out.append(V[f"{k}/{aa}/{fam}/{ss:g}nAm/{ii['key']}"][f])
        return out

    det_hi = vals(s=hi_s, inv="ecd", f="detected_pct")
    det_lo = vals(s=lo_s, inv="ecd", f="detected_pct")
    pairs = {inv["key"]: [v for t, v in V.items() if t.endswith("/" + inv["key"]) and "mcnemar_exact_p" in v]
             for inv in INVERSES if inv["paired_joint"]}
    sig = {k: [v for v in vs if v["mcnemar_exact_p"] < 0.05] for k, vs in pairs.items()}
    need(all(v["only_opm"] > v["only_neuromag"] for vs in sig.values() for v in vs), "a significant joint difference favours "
         "Neuromag (caption wording)")
    sig_lo = {k: len([t for t, v in V.items() if t.endswith("/" + k) and f"/{lo_s:g}nAm/" in t and v.get("mcnemar_exact_p", 1) < 0.05])
              for k in sig}

    def by_sys(s, inv):
        return "; ".join(f"{lab} {span(vals(a=a, s=s, inv=inv))} %" for a, lab in
                         ((SQUID, "Neuromag"), ("opm_dense", "dense"), ("opm_matched", "site-matched")))

    inputs = ([f"{f} :: results['<squid|opm_dense|opm_matched>/<focal|patch>/<80|320>nAm'].<detected|"
               "joint_detect_and_ecd_within_10mm|joint_detect_and_dspm_mne_within_10mm|joint_detect_and_dspm_within_10mm>"
               for f in files]
              + [f"{f} :: paired['<opm>_vs_squid/<condition>'].<joint_ecd_10mm|joint_dspm_10mm> (only_opm, only_squid, "
                 "mcnemar_exact_p)" for f in files]
              + [f"{f} :: thresholds_heldout.false_events[squid|opm_dense|opm_matched].rate_per_min, views[*].channels" for f in files]
              + [f"{rel(LOC_EVENTS, k)} :: detected and the three errors per (array, event): every share and discordant count "
                 "drawn here recomputed and checked" for k in keys]
              + [f"{files[0]} :: config.calibration_min, config.holdout_min", f"{CONFIG} :: events.patch_radius_mm"])
    description = (
        "Joint detection and localization in the exploratory localization study. For each anatomy (rows), system (Neuromag with "
        "all channels, dense OPM array, site-matched OPM array) and condition (columns: focal dipole or "
        f"{D['patch_radius']:g}-mm patch at {hi_s:g} and {lo_s:g} nAm), the bar length is the stored share of the {D['n']} "
        "injected events detected by the practical detector (results[...].detected) and its dark part the stored share detected "
        f"and localized within {D['radius']:g} mm (results[...].joint_detect_and_<method>_within_10mm; rows of panels: dipole "
        "fit, dSPM with MNE-Python, dSPM with the study's implementation). Both shares have the injected events as their "
        "denominator, so neither is conditioned on detection; the pale part is the share detected but localized more than "
        f"{D['radius']:g} mm away, the rest to 100 % the share missed. Stars: the stored exact McNemar p of the paired OPM minus "
        "Neuromag joint-success endpoint (dipole fit and study's dSPM only; none is stored for MNE-Python's dSPM) below 0.05, "
        "uncorrected. Localization is read at the true peak sample (oracle timing). Neuromag's magnetometers or gradiometers alone "
        "(secondary comparators) are not drawn. The per-event tables reproduce every share and discordant count drawn.")
    alt = (f"Twelve small bar charts: three inverse methods by four conditions, nine anatomies each with three bars (Neuromag, "
           f"dense OPM, site-matched OPM). At {hi_s:g} nAm the systems detect {span(det_hi)} % of the events; the share also "
           f"localized within {D['radius']:g} mm is, for the dipole fit, {by_sys(hi_s, 'ecd')}; for MNE-Python's dSPM "
           f"{by_sys(hi_s, 'dspm_mne')}. At {lo_s:g} nAm detection falls to {span(det_lo)} % and joint success with it.")
    caption = (
        f"Joint detection and localization of the simulated spikes. Bar length: share of the {D['n']} injected "
        "events per condition detected by the practical detector (thresholds frozen at "
        f"{D['rate']:g} false event per minute); dark part: detected and localized within {D['radius']:g} mm by (a) the dipole fit, "
        f"(b) dSPM through MNE-Python, (c) dSPM through the study's implementation. At {hi_s:g} nAm the systems detect "
        f"{span(det_hi)} % of the events; detected and localized within {D['radius']:g} mm: with the dipole fit "
        f"{by_sys(hi_s, 'ecd')}; with MNE-Python's dSPM {by_sys(hi_s, 'dspm_mne')}. At {lo_s:g} nAm they detect "
        f"{span(det_lo)} %, which bounds joint success. Of the "
        f"{len(pairs['ecd'])} paired joint-success comparisons per method, {len(sig['ecd'])} (dipole fit) and {len(sig['dspm'])} "
        f"(study's dSPM; {sig_lo['dspm']} of them at {lo_s:g} nAm) have p < 0.05 (exact McNemar, uncorrected), all favouring the "
        "OPM array. Localization is read at the true peak sample, so joint success assumes the event time is known.")
    values = dict(order=keys, conditions=r["conds"], radius_mm=D["radius"], heldout_false_per_min_range=r["heldout_range"],
                  n_mcnemar_below_005={k: len(v) for k, v in sig.items()}, n_mcnemar=dict((k, len(v)) for k, v in pairs.items()),
                  shares=V)
    return dict(inputs=inputs, description=description, alt=alt, caption_draft=caption, values=values)


def entry_curves(r) -> dict:
    V, bands = r["values"], r["bands"]
    mp = r["mode_practical"]
    s50 = {k: V[f"{mp}/{k}/depth0"]["s50_nAm"]["value"] for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm")}
    s50o = {k: V[f"oracle/{k}/depth0"]["s50_nAm"]["value"] for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm")}
    need(all(v is not None for v in list(s50.values()) + list(s50o.values())), "an S50 at 10-20 mm is censored (caption wording)")
    top = {k: V[f"{mp}/{k}/depth{len(bands) - 1}"]["p"][-1] for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm")}
    factor = {k: s50[k] / s50o[k] for k in s50}
    inputs = [f"{DET_ADULT} :: detectors['<detector>'].curves['<oracle|{mp}>/focal/depth<b>/<strength>nAm'] (n, p, ci = Wilson 95 %)",
              f"{DET_ADULT} :: detectors['<detector>'].heldout_false_per_min['{r['mode_practical'].split('@')[1]}'], "
              "detectors['<detector>'].strength_for_50pct_nAm (values only, in figures_supplement.json)",
              f"{DET_ADULT} :: config (events.strengths_nAm, events.stretches, detector.arrays, detector.operating_points_per_min, "
              "detector.oracle_alpha, null.calibration_min, null.heldout_min, simulation.seed), n_dictionary, locations[].stratum and "
              "locations[].depth_mm (locations per band; every location inside its band)",
              f"{G2} :: arrays.squid.channels, arrays.squid.sites, arrays.<opm_dense|opm_matched>.n_sites (checked against "
              f"{rel(LOC_SUMMARY, 'adult')} views[*].channels)",
              f"{DOC} :: section 9 'stratified by depth (10-20, 20-30, 30-45, 45-70 mm)'; 'The Wilson bands in the figure are "
              "event-level and descriptive'"]
    description = (
        f"Clean redraw of results/g4/Figure_G4_detection.png (exploratory run, simulated at {r['simulated_at']}, seed {r['seed']}): "
        "the adult's stored detection probability of focal events (three morphologies pooled) per strength, depth band and "
        "detector, with the stored Wilson 95 % intervals (recomputed here from the stored n and p for all "
        f"{r['n_wilson_checked']} points: identical). Columns: depth bands from docs/methods.md section 9, checked against the "
        f"stored location depths ({r['n_locations_per_band']} locations per band, {r['n_events_per_point']} events per point). "
        f"Rows: oracle detector (knows source, waveform and peak time; per-trial false-positive probability {r['oracle_alpha']:g}) "
        f"and practical detector at {mp.split('@')[1]} false event per minute (thresholds from {r['calibration_min']:g} min of "
        f"null data, frozen; held-out rates over {r['heldout_min']:g} min in the footnote). Solid lines with bands: Neuromag (all "
        "channels), dense and site-matched OPM arrays; grey lines without bands: Neuromag's gradiometers (dashed) or "
        "magnetometers (dotted) "
        "alone (secondary comparators). The bands treat each point's events as independent; the events share locations, so they "
        "are descriptive (docs/methods.md section 9). The stored S50 values (location bootstrap) are listed in values, not drawn.")
    alt = (f"Eight line charts of detection probability against focal spike strength from {r['strengths'][0]:g} to "
           f"{r['strengths'][-1]:g} nAm on a log axis, for four depth bands (columns) and two detectors (rows). Curves rise from 0 "
           f"to near 1 and shift to higher strengths with depth. At {bands[0][0]:g}{EN}{bands[0][1]:g} mm the practical detector "
           f"reaches 50 % at about {s50['opm_dense/opm']:.0f} nAm with the dense OPM array, {s50['opm_matched/opm']:.0f} nAm with "
           f"the site-matched array and {s50['squid/combined']:.0f} nAm with Neuromag, {span(factor.values(), '{:.1f}')} times "
           f"the oracle's strengths; in the deepest band it detects {min(top.values()):.2f} to {max(top.values()):.2f} of the "
           f"{r['strengths'][-1]:g}-nAm events.")
    caption = (
        "Detection of simulated spikes in the adult, exploratory run: probability of detecting a focal spike "
        f"({r['n_events_per_point']} events per point: {r['n_locations_per_band']} locations x {len(r['stretches'])} morphologies) "
        "against its strength, per depth band, for (a) the oracle detector and (b) the practical detector (thresholds frozen at "
        f"{mp.split('@')[1]} false event per minute). Shaded: Wilson 95 % intervals at the event level, descriptive only (the "
        "events of a point share locations; the paired comparisons treat the location as the unit). The strength for 50 % "
        f"detection at {bands[0][0]:g}{EN}{bands[0][1]:g} mm is {s50o['opm_dense/opm']:.0f}, {s50o['opm_matched/opm']:.0f} and "
        f"{s50o['squid/combined']:.0f} nAm with the oracle and {s50['opm_dense/opm']:.0f}, {s50['opm_matched/opm']:.0f} and "
        f"{s50['squid/combined']:.0f} nAm with the practical detector (dense OPM, site-matched OPM, Neuromag).")
    values = dict(bands_mm=[list(b) for b in bands], strengths_nAm=r["strengths"], n_locations_per_band=r["n_locations_per_band"],
                  n_events_per_point=r["n_events_per_point"], heldout_false_per_min=r["heldout_per_min"], labels=r["labels"],
                  n_channels=r["n_channels"], curves=V)
    return dict(inputs=inputs, description=description, alt=alt, caption_draft=caption, values=values)


def main():
    style.apply()
    D = load_localization()
    r_eff = fig_effects(D)
    r_joint = fig_joint(D)
    r_curves = fig_curves()
    entries = {"Figure_S_localization_effects": entry_effects(D, r_eff),
               "Figure_S_joint_detection_localization": entry_joint(D, r_joint),
               "Figure_S_detection_curves_adult": entry_curves(r_curves)}
    style.write_provenance(OUT_JSON, entries)
    for name in entries:
        print(f"results/report/{name}.png")
    print(f"results/report/{OUT_JSON}")


if __name__ == "__main__":
    main()
