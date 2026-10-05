#!/usr/bin/env python3
"""Report figures R1-R4: the adult comparison (G2), drawn from the stored G2 outputs.

Nothing is simulated or re-analysed. R1-R3 plot values stored in results/g2/g2_summary.json and
results/g2/g2_band_sensitivity.json (stored log2 medians and intervals are drawn as ratios 2**x);
R1 adds the per-depth-bin 95 % parcel-bootstrap intervals of results/g2/g2_depth_bins.json
(scripts/study_g2_depth_bins.py); R4 adds medians and quartiles of the stored per-target values in
results/g2/g2_targets.csv. Every number in the labels and texts is read from those files, the stored
config or docs/provenance_register.md. Every ratio axis of R1, R3 and R4 also carries a secondary axis
in dB (20 log10 of the ratio; for the detectability ratio this is D), checked tick by tick against the
ratio axis when drawn (db_axes in figures_adult.json).

Outputs (results/report/): Figure_R1_adult_depth.png, Figure_R2_noise_model.png,
Figure_R3_conditions.png, Figure_R4_regions_adult.png and figures_adult.json (inputs, description,
alt text and a draft caption for every figure).

Usage: PYTHONPATH=src .venv/bin/python scripts/report_figures_adult.py
"""
from __future__ import annotations

import json
import sys
import textwrap
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

from opmsquid import io  # noqa: E402

SUMMARY = "results/g2/g2_summary.json"
BANDS = "results/g2/g2_band_sensitivity.json"
TARGETS = "results/g2/g2_targets.csv"
REGISTER = "docs/provenance_register.md"
DEPTH_BINS = "results/g2/g2_depth_bins.json"  # R1: per-bin 95 % intervals (scripts/study_g2_depth_bins.py)
ARRAYS = ("opm_dense", "opm_matched")
MARKER = {"opm_dense": "o", "opm_matched": "s", "triaxial": "^"}
REF_LINE = dict(color="0.2", lw=0.8, zorder=1)
IQR_DASH = (0, (2.2, 1.3))
# R1: the three noise conditions, drawn per array (the projection's term count is filled in from the stored run)
CONDS = (("intrinsic", ":", "sensor noise only"), ("intrinsic+brain", "-", "sensor plus brain noise (primary)"),
         ("projected", "--", "+ room field, after the {n_ext}-term projection"))
R1_ARRAY = {"opm_dense": "dense OPM array", "opm_matched": "site-matched OPM array"}  # the manuscript's terms
R1_CI_CONDS = ("intrinsic+brain", "projected")  # conditions drawn with their 95 % intervals
R1_DODGE = {("opm_dense", "intrinsic+brain"): -0.75, ("opm_dense", "projected"): -0.25,
            ("opm_matched", "intrinsic+brain"): 0.25, ("opm_matched", "projected"): 0.75}  # mm, so the bars do not overlap
R1_DODGE_PEAK = {"opm_dense": -0.3, "opm_matched": 0.3}
R1_BAR = dict(fmt="none", elinewidth=0.9, capsize=1.6, capthick=0.9)
# secondary dB axes (20 log10 of the ratio; D for a detectability ratio): tick spacing in dB per panel
R1_DB_STEP = (2.0, 5.0)  # (a) detectability ratio, (b) peak-field ratio
R3_DB_STEP = 1.0
R4_DB_STEP = (1.0, 5.0)  # (a, b) detectability ratio, (d) peak-field ratio


def interval_bars(ax, x, lo, hi, color):
    """Bars from lo to hi at x (drawn about the interval's own centre, so they show the interval exactly)."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    ax.errorbar(x, 0.5 * (lo + hi), yerr=0.5 * (hi - lo), ecolor=color, zorder=4, **R1_BAR)
# R4: the six lobes (by_lobe of g2_summary.json) and three single parcels (both hemispheres pooled)
LOBES = (("frontal", "Frontal"), ("parietal", "Parietal"), ("temporal", "Temporal"), ("occipital", "Occipital"),
         ("cingulate", "Cingulate"), ("insula", "Insula"))
PARCELS = (("precentral", "Precentral", ""), ("superiortemporal", "Superior temporal", "lateral temporal"),
           ("parahippocampal", "Parahippocampal", ""))
# a lobe whose parcel bootstrap would rest on fewer parcels than this gets no interval (the insula: one parcel per
# hemisphere): no interval rather than a crude one
MIN_CI_PARCELS = 3


# ----------------------------------------------------------------------------------------------
def load():
    """The three G2 result files; the per-target CSV as numpy columns."""
    s = json.loads((ROOT / SUMMARY).read_text())
    b = json.loads((ROOT / BANDS).read_text())
    rows = io.read_csv(ROOT / TARGETS)
    t = {k: np.array([r[k] for r in rows]) for k in ("region", "lobe")}
    for k in rows[0]:
        if k not in t and k not in ("hemi", "vertno"):
            t[k] = np.array([float(r[k]) for r in rows])
    return s, b, t


def register(row_id: str) -> str:
    """Value cell of a docs/provenance_register.md row (e.g. 'J-h' -> '95 mm')."""
    for line in (ROOT / REGISTER).read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) > 2 and cells[0] == row_id:
            return cells[2]
    raise KeyError(row_id)


def study(s):
    """Labels and constants of the stored G2 run (config, array sizes, projection rank)."""
    cfg, sq = s["config"], s["arrays"]["squid"]
    lo, hi = cfg["band"]["l_freq_hz"], cfg["band"]["h_freq_hz"]
    return dict(q=cfg["sources"]["focal_nAm"], band=f"{lo:g}-{hi:g} Hz", width=hi - lo,
                asd=cfg["sensors"]["opm_asd_primary_fT_per_rtHz"], n_all=sq["channels"], n_mag=sq["sites"],
                n_grad=sq["channels"] - sq["sites"],
                n_ext=s["retained_rank"]["squid/intrinsic+brain+env"] - s["retained_rank"]["squid/projected"])


def ratio(entry):
    """(median ratio, CI low, CI high) of a stored comparison {median_log2, ci95}."""
    lo, hi = entry["ci95"]
    return 2.0 ** entry["median_log2"], 2.0 ** lo, 2.0 ** hi


def binned(rows):
    """Bin centres [mm] and the stored rows of the depth bins that have a median."""
    keep = [r for r in rows if r["median"] is not None]
    return np.array([0.5 * (r["lo"] + r["hi"]) for r in keep]), keep


def log_axis(ax, axis, ticks, lim):
    """Log scale with plain-number ticks and a reference line at a ratio of 1."""
    getattr(ax, f"set_{axis}scale")("log")
    getattr(ax, f"set_{axis}lim")(*lim)
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    a.set_minor_locator(NullLocator())
    (ax.axvline if axis == "x" else ax.axhline)(1.0, **REF_LINE)


def to_db(r):
    """20 log10 of a ratio (an amplitude ratio in dB; for the detectability ratio this is D)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20.0 * np.log10(r)


def from_db(x):
    return 10.0 ** (np.asarray(x, dtype=float) / 20.0)


def db_text(t):
    """Tick label of a dB value: '0', '+2', '−2'."""
    return "0" if t == 0 else f"{t:+g}".replace("-", "−")


def db_axis(ax, axis, step, label):
    """Secondary axis in dB (20 log10 of the ratio) beside a logarithmic ratio axis whose limits are final (set by
    log_axis): on the right of a ratio y axis, on top of a ratio x axis, ticks at the multiples of `step` dB. dB is
    linear in log(ratio), so a linear twin axis whose limits are the parent's limits in dB lines up exactly
    (matplotlib's secondary_*axis cannot show negative values beside a log axis). The ratio axis is unchanged;
    check_db_axes() verifies the alignment once drawn."""
    lo, hi = to_db(np.array(ax.get_ylim() if axis == "y" else ax.get_xlim()))
    ticks = [float(t) + 0.0 for t in np.arange(np.ceil(min(lo, hi) / step) * step, max(lo, hi) + 1e-9, step)]  # no -0.0
    tw = ax.twinx() if axis == "y" else ax.twiny()
    side = "right" if axis == "y" else "top"
    for name, sp in tw.spines.items():
        sp.set_visible(name == side)
    (tw.set_ylim if axis == "y" else tw.set_xlim)(lo, hi)
    a = tw.yaxis if axis == "y" else tw.xaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FixedFormatter([db_text(t) for t in ticks]))
    a.set_minor_locator(NullLocator())
    (tw.set_ylabel if axis == "y" else tw.set_xlabel)(label)
    DB_AXES.append((ax, tw, axis, ticks))
    return tw


DB_AXES = []  # (parent, dB twin, axis, dB ticks) of every dB axis drawn, for check_db_axes


def check_db_axes(fig):
    """Every dB tick sits where the ratio axis shows 10**(dB/20) (drawn positions, within 0.01 pixel), the parent's
    limits are those the dB axis was built for, and a few hand-checked conversions hold (ratio 1.144 = +1.17 dB,
    2 = +6.02 dB, 0.5 = -6.02 dB). Returns the checked ticks."""
    fig.canvas.draw()
    out = []
    for ax, tw, axis, ticks in [d for d in DB_AXES if d[0].get_figure(root=True) is fig]:  # also axes in subfigures
        i = 1 if axis == "y" else 0
        lim_p = np.array(ax.get_ylim() if axis == "y" else ax.get_xlim())
        lim_t = np.array(tw.get_ylim() if axis == "y" else tw.get_xlim())
        assert np.allclose(to_db(lim_p), lim_t, rtol=0, atol=1e-12), (axis, lim_p, lim_t)
        for t in ticks:
            pt, pr = [0.0, 0.0], [1.0, 1.0]
            pt[i], pr[i] = t, float(from_db(t))
            p_tw = tw.transData.transform([pt])[0][i]
            p_par = ax.transData.transform([pr])[0][i]
            assert abs(p_tw - p_par) < 0.01, (axis, t, p_tw, p_par)
        shown = [float(x) for x in (tw.get_yticks() if axis == "y" else tw.get_xticks())]
        assert shown == ticks, (shown, ticks)
        out.append(dict(axis={"y": "right", "x": "top"}[axis], ticks_dB=ticks,
                        at_ratio=[round(float(from_db(t)), 4) for t in ticks]))
    assert f"{to_db(1.144):+.2f}" == "+1.17" and f"{to_db(2.0):.2f}" == "6.02" and f"{to_db(0.5):.2f}" == "-6.02"
    return out


def sites(s, a):
    return f"{style.ARRAY_LABEL[a]} ({s['arrays'][a]['n_sites']} sites)"


def span(values, fmt="{:.2f}"):
    """'lo-hi' of a sequence."""
    return f"{fmt.format(min(values))}-{fmt.format(max(values))}"


def mm(r):
    return f"{r['lo']:g}-{r['hi']:g} mm"


def plain(label):
    """A figure label in plain text (for alt text and captions)."""
    return label.replace("\u221aHz", "sqrt(Hz)").replace("\u00d7", "x")


# ----------------------------------------------------------------------------------------------
def r1_sites(s, a):
    """An OPM array in the manuscript's terms: 'dense OPM array (208 sites)'."""
    return f"{R1_ARRAY[a]} ({s['arrays'][a]['n_sites']} sites)"


def depth_intervals(s):
    """The per-bin 95 % intervals of results/g2/g2_depth_bins.json, checked against the stored bins and medians:
    {key: {(lo, hi): (ratio, lo95, hi95, n_parcels)}} for '<array>/combined/<condition>' and 'peak_field/<array>'."""
    d = json.loads((ROOT / DEPTH_BINS).read_text())
    if d["inputs"]["summary_commit"] != s["provenance"]["commit"]:
        raise ValueError(f"{DEPTH_BINS} was computed from another {SUMMARY} (commit {d['inputs']['summary_commit']})")
    out = {}
    for key, c in d["comparisons"].items():
        stored = (s["log2_ratio_vs_depth"][key] if not key.startswith("peak_field/") else
                  s["bridge_to_sphere"][f"{key.split('/')[1]}_ratio_vs_depth"])
        assert [(b["lo"], b["hi"], b["n"]) for b in c["bins"]] == [(r["lo"], r["hi"], r["n"]) for r in stored], key
        out[key] = {(b["lo"], b["hi"]): (b["ratio"], *b["ci95_ratio"], b["n_parcels"]) for b in c["bins"] if "ratio" in b}
        assert set(out[key]) == {(r["lo"], r["hi"]) for r in stored if r["median"] is not None}, key
    return d, out


def figure_depth(s):
    """R1: (a) paired detectability ratio vs depth, three noise conditions, with per-bin 95 % parcel-bootstrap intervals;
    (b) peak-field ratio vs depth, with the same intervals."""
    c = study(s)
    L, br = s["log2_ratio_vs_depth"], s["bridge_to_sphere"]
    dbins, ci = depth_intervals(s)
    nb = dbins["method"]["n_boot"]
    conds = [(k, ls, lab.format(n_ext=c["n_ext"])) for k, ls, lab in CONDS]
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(style.FULL_W, 6.2), sharex=True, layout="constrained")
    med = {}
    for a in ARRAYS:
        col = style.ARRAY_COLOR[a]
        for cond, ls, _ in conds:
            x, rows = binned(L[f"{a}/combined/{cond}"])
            med[a, cond] = (x, 2.0 ** np.array([r["median"] for r in rows]), rows)
            primary = cond == "intrinsic+brain"
            xs = x + R1_DODGE.get((a, cond), 0.0)
            if primary:
                ax.fill_between(xs, 2.0 ** np.array([r["q25"] for r in rows]), 2.0 ** np.array([r["q75"] for r in rows]),
                                color=col, alpha=0.15, lw=0, zorder=2)
            ax.plot(xs, med[a, cond][1], ls, color=col, lw=1.5 if primary else 1.2, zorder=3)  # no markers: the bars mark the bins
            if cond in R1_CI_CONDS:
                iv = np.array([ci[f"{a}/combined/{cond}"][(r["lo"], r["hi"])][1:3] for r in rows]).T
                interval_bars(ax, xs, iv[0], iv[1], col)
    key0 = "opm_dense/combined/intrinsic+brain"
    x, rows = binned(L[key0])
    counts = [r["n"] for r in rows]
    assert all([r["n"] for r in binned(L[f"{a}/combined/{k}"])[1]] == counts for a in ARRAYS for k, _, _ in CONDS)
    width = rows[0]["hi"] - rows[0]["lo"]
    dropped = [r for r in L[key0] if r["median"] is None and r["n"] > 0]
    for xc, n in zip(x, counts):
        ax.text(xc, 0.025, f"{n:,}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color="0.3")
    ax.text(0.003, 0.11, "targets per bin", transform=ax.transAxes, ha="left", va="bottom", fontsize=7, color="0.3")
    log_axis(ax, "y", [0.5, 0.75, 1, 1.25, 1.5, 2], (0.33, 2.3))
    db_axis(ax, "y", R1_DB_STEP[0], "D (dB)")
    ax.set_ylabel(f"Detectability ratio, OPM over\nNeuromag ({c['n_all']} channels), log scale")
    ax.set_title(f"(a) Detectability ratio: median over the targets of each {width:g}-mm depth bin", loc="left")
    h = [Line2D([], [], color=style.ARRAY_COLOR[a], lw=1.5, label=r1_sites(s, a)) for a in ARRAYS]
    h += [Patch(color="0.5", alpha=0.3, label="IQR over targets (primary)"),
          Line2D([], [], color="0.15", lw=0, marker="|", ms=7, mew=0.9, label="95 % interval of the bin median")]
    h += [Line2D([], [], color="0.15", ls=ls, lw=1.3, label=lab) for _, ls, lab in conds]
    ax.legend(handles=h, loc="upper right", ncol=2, handlelength=2.6, columnspacing=1.2, fontsize=7)

    sph, centers = br["sphere_ratio_vs_depth"], np.array(br["depth_centers_mm"])
    lead = {}
    for a in ARRAYS:
        x, rows = binned(br[f"{a}_ratio_vs_depth"])
        lead[a] = np.array([r["median"] for r in rows])
        xs = x + R1_DODGE_PEAK[a]
        bx.fill_between(xs, [r["q25"] for r in rows], [r["q75"] for r in rows], color=style.ARRAY_COLOR[a], alpha=0.15, lw=0,
                        zorder=2)
        bx.plot(xs, lead[a], "-", color=style.ARRAY_COLOR[a], lw=1.5, zorder=3,
                label=f"{r1_sites(s, a)}: median, IQR (shaded) and 95 % interval")
        iv = np.array([ci[f"peak_field/{a}"][(r["lo"], r["hi"])][1:3] for r in rows]).T
        interval_bars(bx, xs, iv[0], iv[1], style.ARRAY_COLOR[a])
    dist = br["sensor_distance_mm"]
    xi_opm, xi_sq = dist["opm_matched"]["median"], dist["squid"]["median"]
    keep = centers <= max(r["hi"] for r in rows)  # the sphere over the depths that hold targets
    sph_r, sph_j = np.array(sph["realistic_standoffs"])[keep], np.array(sph["jas_xi0_18"])[keep]
    bx.plot(centers[keep], sph_r, "--", color="k", lw=1.1,
            label=f"sphere, standoffs {xi_opm:.1f} / {xi_sq:.1f} mm (real medians, OPM / Neuromag)")
    bx.plot(centers[keep], sph_j, ":", color="k", lw=1.4, label=f"sphere, standoffs {register('J-xi')} (Jas et al.)")
    log_axis(bx, "y", [1, 1.5, 2, 3, 4, 6], (0.92, 8.0))
    db_axis(bx, "y", R1_DB_STEP[1], "dB")
    bx.set_ylabel("Peak-field ratio, OPM over\nNeuromag magnetometer, log scale")
    bx.set_title("(b) Signal only: peak field of the OPM array over the best Neuromag magnetometer, median per bin", loc="left")
    bx.legend(loc="upper right", handlelength=2.6, fontsize=7)
    bx.set_xlim(10, 65)
    bx.set_xticks(np.arange(10, 66, 5))
    bx.set_xlabel(f"Depth below the scalp (mm; distance to the nearest scalp point; {c['q']:g}-nAm cortical-normal dipoles)")
    db_checked = check_db_axes(fig)
    style.save(fig, "Figure_R1_adult_depth")

    # provenance and caption, every number from the dicts above
    d_ib, m_ib = med["opm_dense", "intrinsic+brain"], med["opm_matched", "intrinsic+brain"]
    d_in, m_in = med["opm_dense", "intrinsic"], med["opm_matched", "intrinsic"]
    d_pr, m_pr = med["opm_dense", "projected"], med["opm_matched", "projected"]
    deep = d_ib[0] > 35.0
    deep_mm = f"{min(r['lo'] for r, k in zip(d_ib[2], deep) if k):g}-{max(r['hi'] for r, k in zip(d_ib[2], deep) if k):g} mm"

    def below(m):
        """First bin below a ratio of 1 (every deeper bin is below 1 too: checked)."""
        i = next(i for i, r in enumerate(m[2]) if r["median"] < 0)
        assert all(r["median"] < 0 for r in m[2][i:])
        return f"{mm(m[2][i])} ({2 ** m[2][i]['median']:.3f})"

    def bins_of(key, side):
        """'10-15, 15-20 mm' or 'none': the bins whose 95 % interval lies above / spans / lies below 1."""
        labs = dbins["comparisons"][key][f"bins_ci_{side}_1"]
        return "none" if not labs else ", ".join(labs) + " mm"

    def ci_at(key, r):
        v, lo, hi, _ = ci[key][(r["lo"], r["hi"])]
        return f"{v:.2f} [{lo:.2f}, {hi:.2f}]"

    first, last = d_ib[2][0], d_ib[2][-1]
    lo_d, lo_m = int(np.argmin(lead["opm_dense"])), int(np.argmin(lead["opm_matched"]))
    lead_rows = binned(br["opm_dense_ratio_vs_depth"])[1]
    # wording guards: the qualitative statements below hold for these values
    assert d_ib[1][0] > 1 and m_ib[1][0] > 1 and max(d_in[1]) < 1.005 and max(m_in[1]) < 1
    assert min(lead["opm_dense"]) > 1 and min(lead["opm_matched"]) > 1 and lo_d == lo_m
    assert all(v[1] > 1 for k in ("peak_field/opm_dense", "peak_field/opm_matched") for v in ci[k].values())
    assert not dbins["comparisons"]["opm_dense/combined/intrinsic+brain"]["bins_ci_spans_1"] + \
        dbins["comparisons"]["opm_dense/combined/intrinsic+brain"]["bins_ci_below_1"]  # 'above 1 in every bin'
    n_t = s["n_targets"]
    n_par = [v[3] for k, x in ci.items() if not k.startswith("peak_field/") for v in x.values()]
    k_ib, k_pr = "opm_dense/combined/intrinsic+brain", "opm_dense/combined/projected"
    m_kib, m_kpr = "opm_matched/combined/intrinsic+brain", "opm_matched/combined/projected"
    entry = dict(
        inputs=[f"{SUMMARY} :: log2_ratio_vs_depth['{a}/combined/{k}'] (lo, hi, n, median, q25, q75)" for a in ARRAYS
                for k, _, _ in CONDS]
        + [f"{SUMMARY} :: bridge_to_sphere.{a}_ratio_vs_depth (lo, hi, median, q25, q75)" for a in ARRAYS]
        + [f"{DEPTH_BINS} :: comparisons['{a}/combined/{k}'].bins[*].ci95_ratio (95 % parcel-bootstrap intervals)"
           for a in ARRAYS for k in R1_CI_CONDS]
        + [f"{DEPTH_BINS} :: comparisons['peak_field/{a}'].bins[*].ci95_ratio" for a in ARRAYS]
        + [f"{DEPTH_BINS} :: method.n_boot, comparisons[*].bins_ci_above_1, bins_ci_spans_1, bins_ci_below_1",
           f"{SUMMARY} :: bridge_to_sphere.sphere_ratio_vs_depth.realistic_standoffs",
           f"{SUMMARY} :: bridge_to_sphere.sphere_ratio_vs_depth.jas_xi0_18",
           f"{SUMMARY} :: bridge_to_sphere.depth_centers_mm",
           f"{SUMMARY} :: bridge_to_sphere.sensor_distance_mm (opm_matched.median, squid.median)",
           f"{SUMMARY} :: arrays.<array>.n_sites, arrays.squid.channels", f"{SUMMARY} :: n_targets",
           f"{SUMMARY} :: config.sources.focal_nAm, config.band, config.sensors.opm_asd_primary_fT_per_rtHz",
           f"{SUMMARY} :: retained_rank (squid/intrinsic+brain+env minus squid/projected = projection terms)",
           f"{REGISTER} :: J-h (sphere radius), J-xi (preprint standoffs)"],
        description=(
            f"Adult comparison (MNE sample subject, {n_t:,} cortical targets, {c['q']:g}-nAm cortical-normal dipoles, "
            f"{c['band']}, OPM {c['asd']:g} fT/sqrt(Hz)) against depth below the scalp (distance to the nearest MRI scalp "
            f"point) in {width:g}-mm bins; bins without a stored median are not drawn ("
            + ", ".join(f"{mm(r)}: {r['n']} targets" for r in dropped) + "). (a) The stored statistic is the median "
            f"over the bin's targets of log2(d_OPM / d_Neuromag) against all {c['n_all']} Neuromag channels, with its 25th "
            "and 75th percentiles (log2_ratio_vs_depth); it is drawn as 2**x, the median paired ratio of known-topography "
            "detectability with the oracle covariance. Line style = noise condition (sensor noise only; sensor plus brain "
            f"noise, the primary condition, with its interquartile range shaded; + room field, after the {c['n_ext']}-term "
            "projection), colour = OPM array. Bars: the 95 % interval of each bin's median in the primary and the projected "
            f"condition, from {nb:,} bootstrap resamples of the parcel labels present in the bin (Desikan-Killiany parcels "
            "and the medial-wall labels; their targets in the bin pooled), the adult analysis's own bootstrap "
            f"({DEPTH_BINS}; {min(n_par)} to {max(n_par)} parcel labels per bin, fewest in the deepest bins, whose intervals "
            "are therefore uneven); the sensor-noise-only medians are drawn without intervals (they are stored in the same "
            "file). The two conditions with intervals are offset by up to 0.75 mm within each bin so that their bars do not "
            "overlap. The grey numbers along the bottom of (a) are the targets per bin. (b) Signal only: per target the peak "
            "|B| of the OPM array over the peak |B| of the Neuromag magnetometers, median and interquartile range per bin "
            f"(bridge_to_sphere.<array>_ratio_vs_depth) with the 95 % interval of the median as in (a) ({DEPTH_BINS} "
            f"peak_field; the two arrays offset by {abs(R1_DODGE_PEAK['opm_dense']):g} mm), against the analytical sphere (radius "
            f"{register('J-h')}, Jas et al. Eq. 1, stored in bridge_to_sphere.sphere_ratio_vs_depth) at the bin centres with "
            "the real median standoffs (site-matched OPM sites and Neuromag magnetometers) and with the standoffs of Jas et al. "
            f"({register('J-xi')}). Both y axes are logarithmic with a line at a ratio of 1; the right-hand axes give the same "
            "ratios in dB, 20 log10 of the ratio (in (a) this is D, the OPM advantage in dB; ticks every "
            f"{R1_DB_STEP[0]:g} dB in (a) and every {R1_DB_STEP[1]:g} dB in (b)). The shading is the spread over "
            "targets, the bars are confidence intervals of the bin medians (within this one anatomy)."),
        alt=(f"Two stacked line charts sharing a depth axis from {first['lo']:g} to {last['hi']:g} mm. Top: the OPM / "
             "Neuromag detectability ratio for the dense and site-matched OPM arrays under three noise conditions, with "
             "short vertical bars giving each bin's 95 % interval; with sensor plus brain noise the dense array falls from "
             f"{d_ib[1][0]:.2f} near the scalp to about {d_ib[1][-1]:.2f}, its intervals above 1 in every bin, and the "
             f"site-matched from {m_ib[1][0]:.2f} to about {m_ib[1][-1]:.2f}; with sensor noise only both are at or below 1; "
             "after the projection both fall below 1 at depth. Bottom: the OPM / magnetometer peak-field ratio falls from "
             f"about {lead['opm_dense'][0]:.1f} to about {lead['opm_dense'][lo_d]:.1f} with depth, close to two analytical "
             "sphere curves. Both ratio axes also carry a dB scale on the right."),
        caption_draft=(
            f"Figure R1. The OPM advantage is largest near the scalp. (a) Median paired detectability ratio OPM / Neuromag "
            f"({c['n_all']} channels) per {width:g}-mm depth bin ({sum(counts):,} of {n_t:,} targets in {len(counts)} bins; "
            f"shaded: interquartile range over targets; bars: 95 % parcel-bootstrap intervals of the bin medians, {nb:,} "
            f"resamples; right-hand axis: D = 20 log10 of the ratio, in dB). With sensor plus brain noise the dense array's "
            f"ratio falls from {ci_at(k_ib, first)} at "
            f"{mm(first)} (n = {first['n']}) to {span(d_ib[1][deep])} at {deep_mm}, its interval above 1 in all "
            f"{len(dbins['comparisons'][k_ib]['bins_ci_above_1'])} bins; the site-matched array's falls from {ci_at(m_kib, first)} to "
            f"{span(m_ib[1][deep])}, its interval above 1 at {bins_of(m_kib, 'above')}, spanning 1 at "
            f"{bins_of(m_kib, 'spans')} and below 1 at {bins_of(m_kib, 'below')}. With sensor noise only the dense ratio is "
            f"{d_in[1][0]:.2f} at {mm(first)} and {span(d_in[1][1:])} deeper, the site-matched {span(m_in[1])}. After the "
            f"{c['n_ext']}-term projection the dense median is below 1 from {below(d_pr)} to {d_pr[1][-1]:.2f} at "
            f"{mm(last)}, its interval below 1 only at {bins_of(k_pr, 'below')} and spanning 1 at {bins_of(k_pr, 'spans')}; "
            f"the site-matched median is below 1 from {below(m_pr)}, its interval below 1 at {bins_of(m_kpr, 'below')}. (b) "
            f"The signal alone (right-hand axis in dB): the median peak-field ratio OPM / Neuromag magnetometer falls from "
            f"{ci_at('peak_field/opm_dense', first)} (dense) and {ci_at('peak_field/opm_matched', first)} (site-matched) at "
            f"{mm(first)} to {lead['opm_dense'][lo_d]:.2f} and {lead['opm_matched'][lo_m]:.2f} at {mm(lead_rows[lo_d])} "
            f"({lead['opm_dense'][-1]:.2f} and {lead['opm_matched'][-1]:.2f} at {mm(last)}); the sphere with the real median "
            f"standoffs ({xi_opm:.1f} and {xi_sq:.1f} mm) gives {sph_r[0]:.2f} to {sph_r[-1]:.2f} at the same bin centres, "
            f"with the standoffs of Jas et al. ({register('J-xi')}) {sph_j[0]:.2f} to {sph_j[-1]:.2f}. The OPM signal is "
            "larger at every depth, every interval above 1; whether its detectability is higher depends on the noise (a). "
            "The intervals resample parcels within this one anatomy."),
        db_axes=db_checked,
    )
    return entry


# ----------------------------------------------------------------------------------------------
def figure_noise(s, bands):
    """R2: (a) in-band noise per channel by component; (b) measured vs model; (c) sensor-noise share of variance."""
    c = study(s)
    comp, val = s["noise_composition"], s["noise_validation"]
    share = bands["bands"][c["band"].replace(" ", "")]["intrinsic_share_of_variance"]  # the primary band
    # (c) on the basis of (a): the sensor component's share of the summed median variances, per channel type
    vshare = {k: comp[k]["intrinsic_rms"] ** 2 / sum(comp[k][p] ** 2 for p in ("intrinsic_rms", "brain_rms", "env_rms"))
              for k in ("opm_dense", "opm_matched", "squid_mag", "squid_grad")}
    meas, model = val["measured"], val["model"]
    expl = val["environment_explained_fraction"]
    pc = val["per_channel_brain_ratio_model_over_measured"]
    sq_asd = s["config"]["sensors"]["squid_asd"]
    fig = plt.figure(figsize=(style.FULL_W, 6.4), layout="constrained")
    top, bot, foot = fig.add_gridspec(3, 1, height_ratios=[1, 1.1, 0.07])
    a1, a2 = (fig.add_subplot(g) for g in top.subgridspec(1, 2, width_ratios=[3, 1.05]))
    b1, b2, c1 = (fig.add_subplot(g) for g in bot.subgridspec(1, 3, width_ratios=[1.2, 1.2, 1.25]))
    parts = (("intrinsic_rms", "sensor"), ("brain_rms", "brain"), ("env_rms", "room field"))

    def bars(ax, names, scale, fmt, width):
        n = len(names)
        for j, (name, col, lab) in enumerate(names):
            xs = np.arange(len(parts)) + (j - (n - 1) / 2) * width
            v = np.array([comp[name][p] for p, _ in parts]) * scale
            ax.bar(xs, v, width * 0.92, color=col, label=lab, zorder=2)
            for x_, v_ in zip(xs, v):
                ax.text(x_, v_, fmt.format(v_), ha="center", va="bottom", fontsize=6.3)
        ax.set_xticks(np.arange(len(parts)), [lab for _, lab in parts])

    bars(a1, [("opm_dense", style.ARRAY_COLOR["opm_dense"], sites(s, "opm_dense")),
              ("opm_matched", style.ARRAY_COLOR["opm_matched"], sites(s, "opm_matched")),
              ("squid_mag", style.ARRAY_COLOR["mag"], style.ARRAY_LABEL["mag"])], 1e15, "{:.0f}", 0.27)
    a1.set_ylim(0, 760)
    a1.set_ylabel("RMS per channel (fT)")
    a1.set_title("(a) Modelled noise per channel by component (in-band RMS)", loc="left")
    a1.legend(loc="upper right")
    a1.text(0.01, 0.97, "No OPM noise was measured:\nthe OPM bars are model only", transform=a1.transAxes,
            ha="left", va="top", fontsize=6.5, style="italic", color="0.25")
    bars(a2, [("squid_grad", style.ARRAY_COLOR["grad"], style.ARRAY_LABEL["grad"])], 1e13, "{:.1f}", 0.55)
    a2.set_ylim(0, 45)
    a2.set_ylabel("RMS per channel (fT/cm)")
    a2.set_title(f"gradiometers ({c['n_grad']})", loc="left")

    def pair(ax, items, unit, ymax, colour, fmt, xlabel):
        for i, (lab, m, mo, status) in enumerate(items):
            ax.bar(i - 0.19, m, 0.36, color="white", edgecolor="0.1", hatch="////", lw=0.8, zorder=2)
            ax.bar(i + 0.19, mo, 0.36, color=colour, zorder=2)
            ax.text(i - 0.19, m, fmt.format(m), ha="center", va="bottom", fontsize=6.3)
            ax.text(i + 0.19, mo, fmt.format(mo), ha="center", va="bottom", fontsize=6.3)
            head, rest = status.split("\n", 1)
            ax.text(i, ymax * 0.985, head, ha="center", va="top", fontsize=6.6, fontweight="bold")
            ax.text(i, ymax * 0.91, rest, ha="center", va="top", fontsize=6.3, linespacing=1.1)
        ax.set_xticks(range(len(items)), [lab for lab, *_ in items])
        ax.set_xlim(-0.6, len(items) - 0.4)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(f"RMS per channel ({unit})")
        ax.set_xlabel(xlabel, fontsize=6.5)

    pc_txt = "{}; per channel, model /\nmeasured brain variance\n{:.2f}\u2013{:.2f} (p5\u2013p95), median {:.2f}"
    pair(b1, [("empty room", meas["empty_room_rms_mag_fT"], model["empty_room_rms_mag_fT"],
               f"FITTED\nroom field:\n{100 * expl['mag']:.0f} % of variance"),
              ("brain noise", meas["brain_rms_mag_fT"], model["brain_rms_mag_fT"],
               f"PREDICTED\nmodel =\n{model['brain_mag_model_over_measured']:.2f}\u00d7 measured")],
         "fT", 370, style.ARRAY_COLOR["mag"], "{:.0f}",
         pc_txt.format(f"magnetometers ({c['n_mag']})", pc["mag"]["p5"], pc["mag"]["p95"], pc["mag"]["median"]))
    b1.set_title("(b) Measured vs model (Neuromag)", loc="left")
    b1.legend(handles=[Patch(facecolor="white", edgecolor="0.1", hatch="////", label="measured"),
                       Patch(facecolor="0.45", label="model")], loc="upper left", bbox_to_anchor=(0.0, 0.7))
    pair(b2, [("empty room", meas["empty_room_rms_grad_fT_cm"], model["empty_room_rms_grad_fT_cm"],
               f"CHECK\ndatasheet sensor noise;\nroom field {100 * expl['grad']:.1f} %"),
              ("brain noise", meas["brain_rms_grad_fT_cm"], model["brain_rms_grad_fT_cm"], "CALIBRATED\nequal by\nconstruction")],
         "fT/cm", 53, style.ARRAY_COLOR["grad"], "{:.1f}",
         pc_txt.format(f"gradiometers ({c['n_grad']})", pc["grad"]["p5"], pc["grad"]["p95"], pc["grad"]["median"]))

    keys = (("opm_dense", "dense\nOPM", "opm_dense"), ("opm_matched", "site-matched\nOPM", "opm_matched"),
            ("squid_mag", "Neuromag\nmag.", "mag"), ("squid_grad", "Neuromag\ngrad.", "grad"))
    v = [100 * vshare[k] for k, _, _ in keys]
    c1.bar(range(len(keys)), v, 0.62, color=[style.ARRAY_COLOR[col] for _, _, col in keys], zorder=2)
    for i, v_ in enumerate(v):
        c1.text(i, v_, f"{v_:.1f} %", ha="center", va="bottom", fontsize=6.5)
    c1.set_xticks(range(len(keys)), [lab for _, lab, _ in keys], fontsize=6.3)
    c1.set_ylim(0, 28)
    c1.set_ylabel("Share of channel variance (%)")
    c1.set_xlabel("sensor-noise variance over the summed median\nvariances of the three components, as in (a)", fontsize=6.3)
    c1.set_title("(c) Sensor-noise share", loc="left")
    note = fig.add_subplot(foot)
    note.axis("off")
    note.text(0.0, 1.0, textwrap.fill(
        f"Note: no OPM noise is measured in this model. The OPM sensor noise is a declared white {c['asd']:g} fT/\u221aHz; the "
        "OPM brain and room noise are the Neuromag-calibrated cortical background and the room field fitted to the Neuromag "
        "empty-room recording, seen through each OPM array's own lead fields and coils. Band "
        f"{c['band'].replace('-', chr(0x2013))}; RMS: the square root of each component's median channel variance.", 128),
        transform=note.transAxes, ha="left", va="top", fontsize=6.5, color="0.2")
    style.save(fig, "Figure_R2_noise_model")

    fT = {k: {p: comp[k][p] * 1e15 for p, _ in parts} for k in ("opm_dense", "opm_matched", "squid_mag")}
    gr = {p: comp["squid_grad"][p] * 1e13 for p, _ in parts}
    ratio_mag = model["brain_mag_model_over_measured"]
    entry = dict(
        inputs=[f"{SUMMARY} :: noise_composition.{k}.{{intrinsic_rms, brain_rms, env_rms}} (T; T/m for squid_grad)"
                for k in ("opm_dense", "opm_matched", "squid_mag", "squid_grad")]
        + [f"{SUMMARY} :: noise_validation.measured.{{empty_room_rms_mag_fT, empty_room_rms_grad_fT_cm, brain_rms_mag_fT, "
           "brain_rms_grad_fT_cm}",
           f"{SUMMARY} :: noise_validation.model.{{empty_room_rms_mag_fT, empty_room_rms_grad_fT_cm, brain_rms_mag_fT, "
           "brain_rms_grad_fT_cm, brain_mag_model_over_measured}",
           f"{SUMMARY} :: noise_validation.environment_explained_fraction.{{mag, grad}}",
           f"{SUMMARY} :: noise_validation.per_channel_brain_ratio_model_over_measured.{{mag, grad}}.{{p5, median, p95}}",
           f"{SUMMARY} :: noise_composition.{{opm_dense, opm_matched, squid_mag, squid_grad}} (panel c: sensor share of the "
           "summed component variances)",
           f"{SUMMARY} :: arrays.<array>.n_sites, arrays.squid.{{sites, channels}}", f"{SUMMARY} :: enbw_hz",
           f"{SUMMARY} :: config.sensors.{{squid_asd, opm_asd_primary_fT_per_rtHz}}, config.band",
           f"{SUMMARY} :: retained_rank (projection terms)"],
        description=(
            f"The adult noise model in the {c['band']} band (G2, MNE sample recording). (a) In-band RMS per channel of each "
            "noise component, the square root of the median over channels of that component's variance "
            "(noise_composition, stored in T and T/m, drawn in fT and fT/cm): sensor (intrinsic) noise = ASD x "
            f"sqrt(ENBW {s['enbw_hz']:.1f} Hz) with {c['asd']:g} fT/sqrt(Hz) for the OPM (declared, not a device value) and "
            f"{sq_asd['mag_fT_per_rtHz']:g} fT/sqrt(Hz) and {sq_asd['grad_fT_per_cm_rtHz']:g} fT/(cm sqrt(Hz)) for the "
            "Neuromag magnetometers and gradiometers (the TRIUX datasheet's typical values); brain = the distributed cortical "
            "background, its "
            f"scale calibrated once on the Neuromag gradiometers; room field = the {c['n_ext']}-term external field fitted "
            "to the Neuromag empty-room recording. The components are drawn side by side because RMS values do not add "
            "(variances do). (b) The only measured quantities, all from the Neuromag recording, against the model: the "
            "empty-room RMS (the model's room field is FITTED to this recording; it explains most of the magnetometers' "
            "empty-room variance but only the stated small share of the gradiometers', whose agreement is therefore a CHECK "
            "of the datasheet sensor noise) and the brain noise (task baseline "
            "minus empty room): the gradiometer level is CALIBRATED, equal by construction, and the magnetometer level is "
            f"the one independent PREDICTION (model / measured RMS {ratio_mag:.2f}). The x-axis labels give the per-channel "
            "spread of the model / measured brain-noise variance ratio (5th and 95th percentiles over channels; a variance "
            "ratio, not an RMS ratio). (c) The sensor-noise share of the in-band variance with brain noise and the room "
            "field, on the basis of (a): the sensor component's median channel variance over the sum of the three "
            "components' median channel variances, for the OPM arrays and for Neuromag's magnetometers and gradiometers "
            "separately. Nothing about the OPM noise is measured: its sensor noise is a declared "
            "white level and its brain and room noise are the Neuromag-calibrated model seen through the OPM arrays' own "
            "lead fields and coil responses."),
        alt=(f"Three bar-chart panels. Top: modelled noise per channel by component; the OPM channels carry "
             f"{fT['opm_dense']['brain_rms']:.0f}-{fT['opm_matched']['brain_rms']:.0f} fT of brain noise against "
             f"{fT['squid_mag']['brain_rms']:.0f} fT at the Neuromag magnetometers, similar room-field noise (about "
             f"{fT['squid_mag']['env_rms']:.0f}-{fT['opm_matched']['env_rms']:.0f} fT) and {fT['opm_dense']['intrinsic_rms']:.0f} "
             f"fT of sensor noise against {fT['squid_mag']['intrinsic_rms']:.0f} fT; gradiometers on their own fT/cm axis. "
             "Bottom left: measured versus model empty-room and brain noise for the Neuromag magnetometers and gradiometers, "
             f"labelled fitted, check, predicted ({ratio_mag:.2f} times) and calibrated. Bottom right: sensor noise is "
             f"{100 * vshare['opm_dense']:.1f} % of the dense OPM channels' variance, {100 * vshare['squid_mag']:.1f} % of "
             f"the Neuromag magnetometers' and {100 * vshare['squid_grad']:.1f} % of its gradiometers'."),
        caption_draft=(
            f"Figure R2. The noise model and what anchors it. (a) Modelled in-band ({c['band']}) RMS noise per channel, "
            f"sensor / brain / room field: dense OPM array {fT['opm_dense']['intrinsic_rms']:.0f} / {fT['opm_dense']['brain_rms']:.0f} / "
            f"{fT['opm_dense']['env_rms']:.0f} fT, site-matched OPM array {fT['opm_matched']['intrinsic_rms']:.0f} / "
            f"{fT['opm_matched']['brain_rms']:.0f} / {fT['opm_matched']['env_rms']:.0f} fT, Neuromag magnetometers "
            f"{fT['squid_mag']['intrinsic_rms']:.0f} / {fT['squid_mag']['brain_rms']:.0f} / {fT['squid_mag']['env_rms']:.0f} fT, "
            f"gradiometers {gr['intrinsic_rms']:.1f} / {gr['brain_rms']:.1f} / {gr['env_rms']:.1f} fT/cm. (b) Measured vs "
            f"model (Neuromag only): empty room {meas['empty_room_rms_mag_fT']:.1f} vs {model['empty_room_rms_mag_fT']:.1f} fT "
            f"(magnetometers) and {meas['empty_room_rms_grad_fT_cm']:.1f} vs {model['empty_room_rms_grad_fT_cm']:.1f} fT/cm "
            f"(gradiometers), with the room field fitted to this recording (explaining {100 * expl['mag']:.0f} % and "
            f"{100 * expl['grad']:.1f} % of the variance, so the gradiometer agreement checks the datasheet sensor noise); "
            f"gradiometer brain noise {meas['brain_rms_grad_fT_cm']:.1f} fT/cm, "
            f"calibrated (equal by construction); magnetometer brain noise {meas['brain_rms_mag_fT']:.0f} fT measured vs "
            f"{model['brain_rms_mag_fT']:.0f} fT predicted ({ratio_mag:.2f}x); per channel the model / measured brain-noise "
            f"variance spans {pc['mag']['p5']:.2f}-{pc['mag']['p95']:.2f} (magnetometers) and {pc['grad']['p5']:.2f}-"
            f"{pc['grad']['p95']:.2f} (gradiometers; 5th-95th percentiles). (c) Sensor noise is {100 * vshare['opm_dense']:.1f} % "
            f"(dense) and {100 * vshare['opm_matched']:.1f} % (site-matched) of the OPM channels' in-band variance, "
            f"{100 * vshare['squid_mag']:.1f} % of the Neuromag magnetometers' and {100 * vshare['squid_grad']:.1f} % of its "
            "gradiometers' (each component's median channel variance). No OPM noise was measured: the OPM sensor "
            "noise is declared, its brain and room noise are the Neuromag-calibrated model."),
    )
    return entry


# ----------------------------------------------------------------------------------------------
def condition_rows(s):
    """R3 rows in drawing order: dict(id, group, label, values {array: (ratio, lo, hi)}, inputs)."""
    c, cfg = study(s), s["config"]
    P, S, J = s["primary"], s["sensitivity"], s["sensitivity_joint_asd_gap"]
    T, C = s["triaxial_control"], s["convergence"]["bem"]
    rows = []

    def add(id_, group, label, values, inputs):
        rows.append(dict(id=id_, group=group, label=label, values=values, inputs=inputs))

    def both(block, path, key):
        return {a: ratio(block[f"{a}/{key}"]) for a in ARRAYS}, [f"{SUMMARY} :: {path}['<array>/{key}']"]

    g = "Comparator (detectability, sensor plus brain noise)"
    add("primary", g, f"vs Neuromag ({c['n_all']} channels), primary", *both(P["oracle"], "primary.oracle", "combined/intrinsic+brain"))
    add("grad", g, f"vs {c['n_grad']} gradiometers", *both(P["oracle"], "primary.oracle", "grad/intrinsic+brain"))
    add("mag", g, f"vs {c['n_mag']} magnetometers", *both(P["oracle"], "primary.oracle", "mag/intrinsic+brain"))
    g = "Metric (sensor plus brain noise)"
    add("peak306", g, f"Peak-channel SNR vs the best of {c['n_all']} channels", *both(P["peak"], "primary.peak", "combined/intrinsic+brain"))
    add("peak102", g, f"Peak-channel SNR vs {c['n_mag']} magnetometers", *both(P["peak"], "primary.peak", "mag/intrinsic+brain"))
    add("meanpow", g, f"Mean-power SNR vs {c['n_all']} channels (as amplitude)",
        *both(P["meanpow_db"], "primary.meanpow_db", "combined/intrinsic+brain"))
    for sec in cfg["covariance"]["estimate_seconds"]:
        lab = f"T{sec:g}"
        n = s["n_estimate_samples"][lab]
        assert n == round(2 * c["width"] * sec)
        add(lab, g, f"Plug-in covariance, {n:,} samples (2BT of {sec:g} s)",
            *both(P[f"plugin_{lab}"], f"primary.plugin_{lab}", "combined/intrinsic+brain"))
    g = f"Noise and hardware (detectability vs {c['n_all']} channels)"
    add("intrinsic", g, "Sensor noise only", *both(P["oracle"], "primary.oracle", "combined/intrinsic"))
    add("env", g, "+ room field", *both(P["oracle"], "primary.oracle", "combined/intrinsic+brain+env"))
    add("projected", g, f"+ room field, after the {c['n_ext']}-term projection", *both(P["oracle"], "primary.oracle", "combined/projected"))
    for asd in cfg["sensors"]["opm_asd_fT_per_rtHz"]:
        if asd != c["asd"]:
            add(f"asd{asd:g}", g, f"OPM sensor noise {asd:g} fT/\u221aHz (primary {c['asd']:g})",
                *both(S[f"opm_asd_{asd:g}fT"], f"sensitivity.opm_asd_{asd:g}fT", "combined/intrinsic+brain"))
    for gap in cfg["sensors"]["opm_scalp_gap_mm"]:
        if gap:
            add(f"gap{gap:g}", g, f"OPM scalp gap {gap:g} mm",
                *both(S[f"gap_{gap:g}mm"], f"sensitivity.gap_{gap:g}mm", "combined/intrinsic+brain"))
    jk = "gap3mm/asd30fT"
    add("joint", g, "30 fT/\u221aHz and a 3-mm gap (joint)",
        {a: (J[jk][a]["ratio"], *J[jk][a]["ci95"]) for a in ARRAYS},
        [f"{SUMMARY} :: sensitivity_joint_asd_gap['{jk}'].<array>.{{ratio, ci95}}"])
    for id_, key, lab in (("tri", "equal_noise", f"Triaxial at the site-matched sites, {T['channels']} channels"),
                          ("tri2", "tangential_noise_x2", "Triaxial, tangential axes at 2\u00d7 noise")):
        add(id_, g, lab, {"triaxial": ratio(T[f"{key}/combined/intrinsic+brain"])},
            [f"{SUMMARY} :: triaxial_control['{key}/combined/intrinsic+brain']"])
    add("bem1", g, f"1-layer BEM ({cfg['convergence']['subset_targets']:,}-target subset; no interval)",
        {a: (2.0 ** C["bem1_5120"]["median_log2"][f"{a}/combined/intrinsic+brain"], None, None) for a in ARRAYS},
        [f"{SUMMARY} :: convergence.bem.bem1_5120.median_log2['<array>/combined/intrinsic+brain']",
         f"{SUMMARY} :: convergence.bem.subset_reference_median_log2['<array>/combined/intrinsic+brain']"])
    return rows


def figure_conditions(s, t):
    """R3: forest plot of the adult median ratio under every stored condition, both arrays."""
    c = study(s)
    rows = condition_rows(s)
    groups = list(dict.fromkeys(r["group"] for r in rows))
    y, ys, heads = 0.0, [], {}
    for gname in groups:
        heads[gname] = y
        y -= 1.0
        for r in rows:
            if r["group"] == gname:
                ys.append(y)
                y -= 1.0
        y -= 0.35
    fig = plt.figure(figsize=(style.FULL_W, 7.4), layout="constrained")
    ax, vx = fig.subplots(1, 2, sharey=True, gridspec_kw=dict(width_ratios=[4.0, 1.2]))
    get = {r["id"]: r["values"] for r in rows}
    for a in ARRAYS:
        ax.axvspan(get["primary"][a][1], get["primary"][a][2], color=style.ARRAY_COLOR[a], alpha=0.13, lw=0, zorder=0)
    off = {"opm_dense": 0.18, "opm_matched": -0.18, "triaxial": -0.18}
    for yy, r in zip(ys, rows):
        if r["id"] == "primary":
            ax.axhspan(yy - 0.48, yy + 0.48, color="0.9", lw=0, zorder=0)
        for a, (v, lo, hi) in r["values"].items():
            col = style.ARRAY_COLOR["opm_matched" if a == "triaxial" else a]
            if lo is not None:
                ax.plot([lo, hi], [yy + off[a]] * 2, "-", color=col, lw=1.3, zorder=3, solid_capstyle="butt")
            ax.plot(v, yy + off[a], MARKER[a], color=col, ms=4.2 if a != "triaxial" else 4.8, zorder=4,
                    mec="white" if a == "opm_dense" else col, mew=0.4)
            vx.text(0.06 if a == "opm_dense" else 0.5, yy, ("\u25b2 " if a == "triaxial" else "") + f"{v:.2f}", ha="left",
                    va="center", fontsize=7, color=col, transform=vx.get_yaxis_transform())
    for gname, yy in heads.items():
        ax.text(0.004, yy, gname, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.8, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=5)
    log_axis(ax, "x", [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4], (0.48, 1.45))
    db_axis(ax, "x", R3_DB_STEP, "dB")
    ax.set_yticks(ys, [r["label"] for r in rows])
    for lab, r in zip(ax.get_yticklabels(), rows):
        if r["id"] == "primary":
            lab.set_fontweight("bold")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(min(ys) - 0.8, 0.7)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"Median ratio OPM / Neuromag over the {s['n_targets']:,} targets (log scale; > 1: OPM higher)")
    vx.axis("off")
    for x_, a, lab in ((0.06, "opm_dense", "dense"), (0.5, "opm_matched", "site-matched")):
        vx.text(x_, 0.7, lab, transform=vx.get_yaxis_transform(), ha="left", va="center", fontsize=7.5, fontweight="bold",
                color=style.ARRAY_COLOR[a])
    n_par = int(s["primary"]["oracle"]["opm_dense/combined/intrinsic+brain"]["ci_method"].split("(")[1].rstrip(")"))
    medial = sorted({r for r in t["region"] if r.endswith("unknown")})
    assert n_par == len(set(t["region"])), "the bootstrap groups are the region labels of the targets"
    h = [Line2D([], [], color=style.ARRAY_COLOR[a], marker=MARKER[a], ms=4.5, lw=1.3, label=sites(s, a)) for a in ARRAYS]
    h += [Line2D([], [], color=style.ARRAY_COLOR["opm_matched"], marker="^", ms=4.8, lw=1.3,
                 label=f"triaxial OPM ({s['triaxial_control']['sites']} sites, {s['triaxial_control']['channels']} channels)"),
          Patch(color="0.5", alpha=0.25, label="primary 95 % interval, per array (shaded)")]
    fig.legend(handles=h, loc="outside lower center", ncol=2, title=f"Lines: 95 % parcel-bootstrap intervals ({n_par} labels)",
               title_fontsize=7.5)
    db_checked = check_db_axes(fig)
    style.save(fig, "Figure_R3_conditions")

    def v(id_, a="opm_dense"):
        return f"{get[id_][a][0]:.2f}"

    def vb(id_):
        return f"{v(id_)} / {v(id_, 'opm_matched')}"

    def ci(id_, a):
        return f"{get[id_][a][0]:.3f} [{get[id_][a][1]:.3f}-{get[id_][a][2]:.3f}]"

    under = {a: [f"{plain(r['label'])}: {r['values'][a][0]:.2f}" for r in rows if a in r["values"] and r["values"][a][0] < 1]
             for a in ARRAYS}
    n_rows = {a: sum(a in r["values"] for r in rows) for a in ARRAYS}
    sub = s["convergence"]["bem"]["subset_reference_median_log2"]
    secs = s["config"]["covariance"]["estimate_seconds"]
    n_est = [s["n_estimate_samples"][f"T{x:g}"] for x in secs]
    asds = [a for a in s["config"]["sensors"]["opm_asd_fT_per_rtHz"] if a != c["asd"]]
    gaps = [g for g in s["config"]["sensors"]["opm_scalp_gap_mm"] if g]
    entry = dict(
        inputs=sorted({k for r in rows for k in r["inputs"]})
        + [f"{SUMMARY} :: n_estimate_samples", f"{SUMMARY} :: triaxial_control.{{sites, channels}}",
           f"{SUMMARY} :: arrays.<array>.n_sites, arrays.squid.{{sites, channels}}", f"{SUMMARY} :: n_targets",
           f"{SUMMARY} :: primary.oracle[...].ci_method ('parcels (N)')",
           f"{SUMMARY} :: config.{{sensors, covariance, convergence.subset_targets, band}}",
           f"{TARGETS} :: region (the bootstrap groups: Desikan-Killiany parcels and medial-wall labels)"],
        description=(
            f"Forest plot of the adult (G2) median over {s['n_targets']:,} cortical targets of the paired ratio OPM / "
            "Neuromag under every stored condition, for the dense and the site-matched OPM arrays and the triaxial control at "
            "the site-matched sites. Every value except the joint row is stored as median_log2 = the median over targets of "
            f"log2(m_OPM / m_Neuromag) with a 95 % CI from a bootstrap over {n_par} parcels (the region labels of the "
            f"targets: {n_par - len(medial)} Desikan-Killiany parcels and {len(medial)} medial-wall labels); it is drawn "
            "as 2**x. The joint OPM-noise x scalp-gap row (sensitivity_joint_asd_gap) is stored as a ratio with its CI "
            "from the same parcel bootstrap and is drawn as stored; the sensitivity rows used fewer bootstrap resamples "
            "than the primary rows (scripts/g2_adult_comparison.py). m is known-topography detectability with the oracle "
            "covariance except in the metric group: peak-channel SNR (the best single channel of each system; against "
            f"the best of all {c['n_all']} channels or of the {c['n_mag']} magnetometers), mean-power SNR (stored per "
            "target in dB; the analysis converted it to an amplitude ratio 10**(dB/20) before the same log2 summary, so "
            "it is on the amplitude scale of detectability), and plug-in detectability (Ledoit-Wolf covariance estimated "
            f"from one common noise realization of {n_est[0]:,} or {n_est[1]:,} independent samples, the 2BT equivalent of "
            f"{secs[0]:g} or {secs[1]:g} s in the {c['width']:g}-Hz-wide band, not estimates from coloured time series). "
            "Conditions vary one factor at a time from the primary (oracle covariance, sensor plus brain noise, OPM "
            f"{c['asd']:g} fT/sqrt(Hz), no scalp gap, 3-layer BEM, measured head position) except the joint row. The "
            "1-layer BEM row is a convergence check on a random "
            f"{s['config']['convergence']['subset_targets']:,}-target subset without an interval; on the same subset the "
            f"3-layer reference is {2 ** sub['opm_dense/combined/intrinsic+brain']:.3f} (dense) and "
            f"{2 ** sub['opm_matched/combined/intrinsic+brain']:.3f} (site-matched). Grey row: the primary result; shaded "
            "vertical bands: each array's primary 95 % CI; logarithmic x axis with a line at 1 and, on top, the same ratio "
            f"in dB (20 log10 of the ratio, which for the detectability rows is D; ticks every {R3_DB_STEP:g} dB); the "
            "right-hand columns repeat the medians (triangle: triaxial)."),
        alt=(f"Forest plot with {len(rows)} rows in three groups (comparator, metric, noise and hardware). Each row shows "
             "the OPM / Neuromag median ratio with a 95 % interval for the dense array (dark blue circles) and the site-matched "
             "array (light blue squares) on a logarithmic axis with a line at 1 and a dB scale on top. The dense array is "
             "below 1 in "
             f"{len(under['opm_dense'])} of {n_rows['opm_dense']} rows (" + "; ".join(under["opm_dense"]) + "); the site-matched "
             f"array, {v('primary', 'opm_matched')} in the primary row, is below 1 in {len(under['opm_matched'])} of "
             f"{n_rows['opm_matched']} rows. Two triaxial rows show a triaxial array at the site-matched sites."),
        caption_draft=(
            "Figure R3. The adult ranking depends on the comparator, the metric and the OPM noise. Median ratio OPM / "
            f"Neuromag over {s['n_targets']:,} targets with parcel-bootstrap 95 % CIs; values dense / site-matched. Primary "
            f"(detectability against all {c['n_all']} channels, sensor plus brain noise): dense {ci('primary', 'opm_dense')}, "
            f"site-matched {ci('primary', 'opm_matched')}. Against the {c['n_grad']} gradiometers alone {vb('grad')}, the "
            f"{c['n_mag']} magnetometers alone {vb('mag')}. Peak-channel SNR {vb('peak306')} against the best of the "
            f"{c['n_all']} channels and {vb('peak102')} against the best magnetometer; mean-power SNR {vb('meanpow')}; "
            f"plug-in covariance {vb(f'T{secs[0]:g}')} ({n_est[0]:,} samples) and {vb(f'T{secs[1]:g}')} ({n_est[1]:,}). "
            f"Sensor noise only {vb('intrinsic')}; with the room field {vb('env')}; after the {c['n_ext']}-term projection {vb('projected')}. OPM "
            f"sensor noise {', '.join(f'{a:g}' for a in asds)} fT/sqrt(Hz): {', '.join(v(f'asd{a:g}') for a in asds)} "
            f"(dense) and {', '.join(v(f'asd{a:g}', 'opm_matched') for a in asds)} (site-matched; primary {c['asd']:g}); scalp "
            "gap " + ", ".join(f"{g:g} mm {vb(f'gap{g:g}')}" for g in gaps)
            + f"; 30 fT/sqrt(Hz) with a 3-mm gap {ci('joint', 'opm_dense')} (dense) and {ci('joint', 'opm_matched')} "
            f"(site-matched). Triaxial OPM at the site-matched sites {v('tri', 'triaxial')} ({v('tri2', 'triaxial')} with the "
            f"tangential axes at twice the noise); 1-layer BEM {vb('bem1')} "
            f"({s['config']['convergence']['subset_targets']:,}-target subset, no interval). One factor varied at a time "
            "except the joint row. Top axis: the same ratios in dB (20 log10 of the ratio)."),
        db_axes=db_checked,
    )
    return entry


# ----------------------------------------------------------------------------------------------
def parcel_stats(t, m, ref, a):
    """Median and quartiles over the targets in mask m of the per-target ratio OPM / Neuromag (CSV values)."""
    r = t[f"detect_{a}_opm_intrinsic+brain"][m] / t[f"detect_squid_{ref}_intrinsic+brain"][m]
    return float(np.median(r)), float(np.percentile(r, 25)), float(np.percentile(r, 75))


def region_values(s, t):
    """Per region: kind, label, n, (ratio, lo, hi) per array vs 306 and vs 102, parcels in a lobe, depth, peak-field ratio."""
    B = s["by_lobe"]
    peak = t["amp_opm_dense"] / t["amp_squid_mag"]
    out = []
    for lobe, name in LOBES:
        m = t["lobe"] == lobe
        n = int(m.sum())
        assert n == B["opm_dense/combined/intrinsic+brain"][lobe]["n"] == B["opm_matched/mag/intrinsic+brain"][lobe]["n"]
        for a in ARRAYS:  # the CSV (rounded values) reproduces the stored lobe medians
            for ref in ("combined", "mag"):
                assert abs(parcel_stats(t, m, ref, a)[0] / ratio(B[f"{a}/{ref}/intrinsic+brain"][lobe])[0] - 1) < 1e-3
        out.append(dict(kind="lobe", name=name, n=n, label=f"{name}\nn = {n:,}",
                        v306={a: ratio(B[f"{a}/combined/intrinsic+brain"][lobe]) for a in ARRAYS},
                        v102={a: ratio(B[f"{a}/mag/intrinsic+brain"][lobe]) for a in ARRAYS},
                        parcels=int(B["opm_dense/combined/intrinsic+brain"][lobe]["ci_method"].split("(")[1].rstrip(")")),
                        depth=float(np.median(t["depth_mm"][m])), peak=float(np.median(peak[m]))))
    for parcel, name, alias in PARCELS:
        m = np.isin(t["region"], [f"lh.{parcel}", f"rh.{parcel}"])
        n = int(m.sum())
        out.append(dict(kind="parcel", name=name, alias=alias, n=n,
                        label=f"{name}\n({alias})\nn = {n:,}" if alias else f"{name}\nn = {n:,}",
                        v306={a: parcel_stats(t, m, "combined", a) for a in ARRAYS},
                        v102={a: parcel_stats(t, m, "mag", a) for a in ARRAYS},
                        depth=float(np.median(t["depth_mm"][m])), peak=float(np.median(peak[m]))))
    return out


def figure_regions(s, t):
    """R4: median ratio per brain region (six lobes, three parcels) vs 306 and vs 102; depth and peak-field ratio."""
    c = study(s)
    regs = region_values(s, t)
    ys, y, head = [], 0.0, {}
    for kind in ("lobe", "parcel"):
        head[kind] = y
        y -= 0.8
        for r in regs:
            if r["kind"] == kind:
                ys.append(y)
                y -= 1.0
        y -= 0.25
    fig = plt.figure(figsize=(style.FULL_W, 6.2), layout="constrained")
    axs = fig.subplots(1, 4, sharey=True, gridspec_kw=dict(width_ratios=[1.65, 1.65, 1.05, 1.05]))
    off = {"opm_dense": 0.19, "opm_matched": -0.19}
    no_ci = [r["name"] for r in regs if r["kind"] == "lobe" and r["parcels"] < MIN_CI_PARCELS]
    assert no_ci == ["Insula"], f"lobes without an interval: {no_ci} (the caption names the insula only)"
    for ax, key, title in ((axs[0], "v306", f"(a) vs all {c['n_all']} channels"),
                           (axs[1], "v102", f"(b) vs {c['n_mag']} magnetometers")):
        for yy, r in zip(ys, regs):
            for a in ARRAYS:
                v, lo, hi = r[key][a]
                ax.barh(yy + off[a], v - 1.0, 0.36, left=1.0, color=style.ARRAY_COLOR[a], zorder=2)
                if r["name"] in no_ci:
                    continue
                eb = ax.errorbar(v, yy + off[a], xerr=[[v - lo], [hi - v]], fmt="none", ecolor="0.05", elinewidth=0.9,
                                 capsize=2.0 if r["kind"] == "lobe" else 0, zorder=3)
                if r["kind"] == "parcel":
                    eb[2][0].set_linestyle(IQR_DASH)
            if r["name"] in no_ci:  # beside the longer bar, inside the row
                x_end = max(max(r[key][a][0] for a in ARRAYS), 1.0)
                ax.text(x_end * 1.03, yy, f"no interval\n({r['parcels']} parcels)", ha="left", va="center", fontsize=6.5,
                        style="italic", color="0.25", linespacing=1.1, zorder=4)
        log_axis(ax, "x", [0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4], (0.78, 1.5))
        db_axis(ax, "x", R4_DB_STEP[0], "D (dB)")
        ax.set_title(title, loc="left")
        ax.set_xlabel("Detectability ratio\nOPM / Neuromag (log)")
    for kind, txt in (("lobe", "Lobes"), ("parcel", "Single parcels")):
        axs[0].text(0.015, head[kind] - 0.05, txt, transform=axs[0].get_yaxis_transform(), ha="left", va="center",
                    fontsize=7.8, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=5)
    ax = axs[2]
    ax.barh(ys, [r["depth"] for r in regs], 0.5, color="0.6", zorder=2)
    for yy, r in zip(ys, regs):
        ax.text(r["depth"] + 1.5, yy, f"{r['depth']:.0f}", ha="left", va="center", fontsize=6.8)
    ax.set_xlim(0, 64)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_title("(c) Depth", loc="left")
    ax.set_xlabel("median depth\nbelow scalp (mm)", fontsize=7.5)
    ax = axs[3]
    ax.barh(ys, [r["peak"] - 1.0 for r in regs], 0.5, left=1.0, color=style.ARRAY_COLOR["opm_dense"], zorder=2)
    for yy, r in zip(ys, regs):
        ax.text(r["peak"] * 1.05, yy, f"{r['peak']:.1f}", ha="left", va="center", fontsize=6.8)
    log_axis(ax, "x", [1, 2, 4], (0.9, 7.0))
    db_axis(ax, "x", R4_DB_STEP[1], "dB")
    ax.set_title("(d) Signal", loc="left")
    ax.set_xlabel("median peak-field ratio\ndense OPM / mag. (log)", fontsize=7.5)
    axs[0].set_yticks(ys, [r["label"] for r in regs])
    for ax in axs:
        ax.tick_params(axis="y", length=0)
    axs[0].set_ylim(min(ys) - 0.65, 0.5)
    h = [Patch(color=style.ARRAY_COLOR[a], label=sites(s, a)) for a in ARRAYS]
    h += [Line2D([], [], color="0.05", lw=0.9, marker="|", ms=5, label="lobes: 95 % parcel-bootstrap interval"),
          Line2D([], [], color="0.05", lw=0.9, ls=IQR_DASH, label="single parcels: IQR over targets (a spread)")]
    fig.legend(handles=h, loc="outside upper center", ncol=2,
               title="Sensor plus brain noise; both hemispheres pooled; bars start at a ratio of 1", title_fontsize=7.5)
    fig.align_titles()  # (c) has no dB axis on top: its title level with the others
    db_checked = check_db_axes(fig)
    style.save(fig, "Figure_R4_regions_adult")

    lobes = [r for r in regs if r["kind"] == "lobe"]
    parcels = [r for r in regs if r["kind"] == "parcel"]

    def nm(r):
        return r["name"].lower() + (f" ('{r['alias']}')" if r.get("alias") else "")

    medial = s["medial_wall"]["n_targets"]
    assert medial == int(np.sum(t["lobe"] == "other"))
    assert all(r[k]["opm_dense"][0] > 1 for r in regs for k in ("v306", "v102"))  # wording guard (alt text)
    entry = dict(
        inputs=[f"{SUMMARY} :: by_lobe['<array>/combined/intrinsic+brain'][<lobe>] (median_log2, ci95, n, ci_method)",
                f"{SUMMARY} :: by_lobe['<array>/mag/intrinsic+brain'][<lobe>] (median_log2, ci95, n, ci_method)",
                f"{TARGETS} :: region, lobe, depth_mm, amp_opm_dense, amp_squid_mag, detect_opm_dense_opm_intrinsic+brain, "
                "detect_opm_matched_opm_intrinsic+brain, detect_squid_combined_intrinsic+brain, "
                "detect_squid_mag_intrinsic+brain",
                f"{SUMMARY} :: medial_wall.n_targets", f"{SUMMARY} :: arrays.<array>.n_sites, arrays.squid.{{sites, channels}}"],
        description=(
            "Adult (G2) comparison by brain region with sensor plus brain noise, oracle covariance and OPM "
            f"{c['asd']:g} fT/sqrt(Hz); both hemispheres pooled. Regions: the six Desikan-Killiany lobes of the study's "
            f"grouping (the {medial} medial-wall targets belong to no lobe) and three single parcels emphasised by the "
            "epilepsy literature: precentral, superior temporal ('lateral temporal') and parahippocampal (one of the two "
            "parcels of the report's mesial temporal group). (a, b) Bars start at a ratio of 1 on a logarithmic axis (the "
            f"top axes give the same ratio in dB, 20 log10 of the ratio, i.e. D; ticks every {R4_DB_STEP[0]:g} dB) and "
            "end at the median over the region's "
            f"targets of the paired detectability ratio OPM / Neuromag, against all {c['n_all']} channels (a) or the "
            f"{c['n_mag']} magnetometers (b). Lobes: median and 95 % parcel-bootstrap interval from g2_summary.json by_lobe "
            "(2**median_log2; a "
            "bootstrap over the lobe's own parcels: " + ", ".join(f"{r['name'].lower()} {r['parcels']}" for r in lobes)
            + f" parcels; a lobe with fewer than {MIN_CI_PARCELS} parcels, the insula, is drawn without an interval and "
            "marked 'no interval' with its parcel count, since a bootstrap over so few parcels would be crude). Single "
            "parcels: the per-target ratio "
            "detect_<array>_opm_intrinsic+brain / detect_squid_<combined|mag>_intrinsic+brain from g2_targets.csv (values "
            "as stored, rounded), unweighted median and interquartile range over the parcel's targets (a spread, not a "
            "confidence interval; a parcel bootstrap is not possible within one parcel). (c) Median depth below the scalp "
            "(distance to the nearest scalp point) of the region's targets; (d) median over the region's targets of the "
            "peak-field ratio amp_opm_dense / amp_squid_mag (peak |B| of the dense OPM array over the best Neuromag "
            "magnetometer; signal only), both from g2_targets.csv, with the same ratio in dB on top (ticks every "
            f"{R4_DB_STEP[1]:g} dB). The n under each label is the number of targets; the "
            "lobe medians recomputed from the CSV agree with by_lobe."),
        alt=("Horizontal bar chart of nine brain regions: six lobes, then the precentral, superior temporal and "
             "parahippocampal parcels. For each region two bars give the OPM / Neuromag detectability ratio of the dense "
             f"and site-matched arrays, against all {c['n_all']} channels and against the {c['n_mag']} magnetometers; the "
             "dense array is above 1 in every region; against all channels the site-matched array is below 1 in "
             + ", ".join(nm(r) for r in regs if r["v306"]["opm_matched"][0] < 1)
             + f". Two narrow panels give the median depth ({span([r['depth'] for r in regs], '{:.0f}')} mm) and the median "
             f"peak-field ratio ({span([r['peak'] for r in regs], '{:.1f}')}). The ratio axes carry dB scales on top."),
        caption_draft=(
            f"Figure R4. The adult result by brain region (sensor plus brain noise; both hemispheres). (a) Median "
            f"detectability ratio OPM / Neuromag ({c['n_all']} channels), dense / site-matched: "
            + ", ".join(f"{r['name'].lower()} {r['v306']['opm_dense'][0]:.2f} / {r['v306']['opm_matched'][0]:.2f}" for r in lobes)
            + " (lobes, 95 % parcel-bootstrap intervals; none for the "
            + ", ".join(f"{r['name'].lower()}, which has {r['parcels']} parcels" for r in lobes if r["name"] in no_ci) + "); "
            + ", ".join(f"{nm(r)} {r['v306']['opm_dense'][0]:.2f} / {r['v306']['opm_matched'][0]:.2f}" for r in parcels)
            + f" (single parcels, interquartile ranges over targets). (b) Against the {c['n_mag']} magnetometers alone: "
            + ", ".join(f"{r['name'].lower()} {r['v102']['opm_dense'][0]:.2f} / {r['v102']['opm_matched'][0]:.2f}" for r in regs)
            + ". (c, d) Median depth below the scalp and median peak-field ratio dense OPM / magnetometer: "
            + ", ".join(f"{r['name'].lower()} {r['depth']:.1f} mm / {r['peak']:.2f}" for r in regs)
            + ". Targets: " + ", ".join(f"{r['name'].lower()} {r['n']:,}" for r in regs) + ". Top axes of (a), (b) and "
            "(d): the same ratios in dB (20 log10 of the ratio)."),
        db_axes=db_checked,
    )
    return entry


def main():
    style.apply()
    s, b, t = load()
    entries = {"Figure_R1_adult_depth": figure_depth(s), "Figure_R2_noise_model": figure_noise(s, b),
               "Figure_R3_conditions": figure_conditions(s, t), "Figure_R4_regions_adult": figure_regions(s, t)}
    style.write_provenance("figures_adult.json", entries)
    for name in entries:
        print(f"results/report/{name}.png")
    print("results/report/figures_adult.json")


if __name__ == "__main__":
    main()
