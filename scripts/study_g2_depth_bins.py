#!/usr/bin/env python3
"""Per-depth-bin intervals of the adult comparison (G2): a 95 % parcel-bootstrap interval for the median ratio in every
5-mm depth bin of results/g2/g2_summary.json (log2_ratio_vs_depth), from the stored per-target values.

Revision: the depth-bin crossings of the adult depth curve were descriptive because no per-bin intervals were stored.
Nothing is simulated here. Inputs:
  results/g2/g2_targets.csv        per-target detectabilities (4 decimals) and peak fields (5 significant digits)
  results/g2/g2_targets_depth.csv  the unrounded depth of every target and its 5-mm bin (column depth_bin_5mm; written by
                                   scripts/export_target_precision.py). g2_targets.csv prints the depth to 0.01 mm, which
                                   puts 14 targets on a bin edge, so its rounded depth does not give the stored bin counts.
  results/g2/g2_summary.json       the bins, their target counts and the stored medians (log2_ratio_vs_depth;
                                   bridge_to_sphere.<array>_ratio_vs_depth for the peak-field ratio)
  configs/g2_adult.toml            sources.seed (the bootstrap streams)

Statistic and interval: those of the adult analysis, through its own function (scripts/g2_adult_comparison.py, compare(),
imported): per bin, the median over the bin's targets of log2(m_OPM / m_Neuromag), and a 95 % interval from the 2.5th
and 97.5th percentiles of the medians of <n_boot> bootstrap resamples of the parcel labels present in the bin
(Desikan-Killiany parcels and the two medial-wall labels, drawn with replacement, their targets in the bin pooled);
n_boot is compare()'s default, the count of the adult analysis's primary intervals. Reported as ratios 2**x. Bins with
fewer targets than binned()'s min_n (10) get no median, as in g2_summary.json. Comparisons: the dense (208 sites) and
site-matched (98 sites) OPM arrays against Neuromag (306 channels), known-topography detectability with the oracle
covariance, in three noise conditions (sensor noise only; sensor plus brain noise, the primary; + room field after the
8-term projection); and the peak-field ratio (peak |B| of the OPM array over the best Neuromag magnetometer, signal only;
g2_summary.json bridge_to_sphere stores its linear median, here the median of the log ratio and its interval).

Checks (nothing is written if one fails): the two CSVs have the same keys (hemi, vertno) in the same order; every printed
depth is the unrounded one rounded to 0.01 mm; every bin label is the bin of the unrounded depth; each bin's target count
is the stored one; every median reproduces the stored one within the precision of the CSV (for each bin, the largest
change that the rounding of the bin's printed values can make to the median: the median moves by no more than the
largest change of any one value).

Random streams: one generator per comparison and bin, numpy default_rng(SeedSequence([seed, comparison index, bin
index])), seed = sources.seed of configs/g2_adult.toml, comparisons indexed in the order of the output.

Output: results/g2/g2_depth_bins.json (provenance by opmsquid.io). Light: seconds, well under 1 GB.
Usage: PYTHONPATH=src .venv/bin/python scripts/study_g2_depth_bins.py
"""
from __future__ import annotations

import inspect
import json
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

import g2_adult_comparison as G2  # noqa: E402  (the adult analysis: compare(), binned(), DEPTH_EDGES)
from opmsquid import io  # noqa: E402

TARGETS = "results/g2/g2_targets.csv"
DEPTH = "results/g2/g2_targets_depth.csv"
SUMMARY = "results/g2/g2_summary.json"
CONFIG = "configs/g2_adult.toml"
OUT = "results/g2/g2_depth_bins.json"
STATUS = ("NEW (revision: 95 % parcel-bootstrap intervals of the adult depth curves per 5-mm depth bin, "
          "from the stored per-target values; nothing simulated)")
N_BOOT = inspect.signature(G2.compare).parameters["n_boot"].default  # the adult analysis's primary resample count
MIN_N = inspect.signature(G2.binned).parameters["min_n"].default  # bins with fewer targets have no stored median
ARRAYS = ("opm_dense", "opm_matched")
CONDITIONS = ("intrinsic", "intrinsic+brain", "projected")
ARRAY_TEXT = {"opm_dense": "dense OPM array (208 sites)", "opm_matched": "site-matched OPM array (98 sites)"}
COND_TEXT = {"intrinsic": "sensor noise only", "intrinsic+brain": "sensor plus brain noise (the primary condition)",
             "projected": "sensor plus brain noise and the room field, after the 8-term projection"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------------------------
def status_line(rel: str) -> str:
    """The '#' status line of a result CSV (its first line), without the '# '."""
    with open(ROOT / rel) as fh:
        first = fh.readline().rstrip("\n")
    if not first.startswith("#"):
        raise ValueError(f"{rel}: no '#' status line")
    return first.lstrip("# ")


def half_unit(text: str) -> float:
    """Half a unit in the last printed digit of a number as printed ('1.3834' -> 5e-5; '8.6787e-14' -> 5e-19)."""
    mant, _, exp = text.lower().partition("e")
    decimals = len(mant.split(".")[1]) if "." in mant else 0
    return 0.5 * 10.0 ** ((int(exp) if exp else 0) - decimals)


def label(lo: float, hi: float) -> str:
    return f"{lo:g}-{hi:g}"


def load_targets() -> dict:
    """The per-target columns, the exact bin label of every target and the printing precision of the value columns."""
    rows = io.read_csv(ROOT / TARGETS)
    comp = io.read_csv(ROOT / DEPTH)
    if len(rows) != len(comp) or any((r["hemi"], r["vertno"]) != (c["hemi"], c["vertno"]) for r, c in zip(rows, comp)):
        raise ValueError(f"{DEPTH}: keys or row order differ from {TARGETS}")
    exact = np.array([float(c["depth_mm"]) for c in comp])
    printed = [r["depth_mm"] for r in rows]
    bad = [i for i, (p, x) in enumerate(zip(printed, exact)) if f"{x:.2f}" != p]
    if bad:
        raise ValueError(f"{DEPTH}: the unrounded depth of {len(bad)} targets does not print as {TARGETS} (first row {bad[0]})")
    edges = G2.DEPTH_EDGES
    rule = np.full(len(exact), "", dtype=object)
    for lo, hi in zip(edges[:-1], edges[1:]):
        rule[(exact >= lo) & (exact < hi)] = label(lo, hi)
    bins = np.array([c["depth_bin_5mm"] for c in comp], dtype=object)
    if not np.array_equal(bins, rule):
        raise ValueError(f"{DEPTH}: depth_bin_5mm is not the 5-mm bin (lo <= depth < hi) of the unrounded depth")
    t = dict(region=np.array([r["region"] for r in rows]), bin=bins, depth=exact,
             depth_printed=np.array([float(p) for p in printed]))
    cols = [f"detect_{a}_opm_{c}" for a in ARRAYS for c in CONDITIONS] + [f"detect_squid_combined_{c}" for c in CONDITIONS]
    cols += [f"amp_{a}" for a in ARRAYS] + ["amp_squid_mag"]
    for k in cols:
        t[k] = np.array([float(r[k]) for r in rows])
        t[f"{k}/half_unit"] = np.array([half_unit(r[k]) for r in rows])
    return t


def precision_bound_log2(a, ea, b, eb) -> float:
    """Largest change of log2(a / b) over the targets when a and b may each be off by half a printed unit."""
    lr = np.log2(a / b)
    hi = np.log2((a + ea) / (b - eb)) - lr
    lo = lr - np.log2((a - ea) / (b + eb))
    return float(np.max(np.maximum(hi, lo)))


def precision_bound_ratio(a, ea, b, eb) -> float:
    """Largest change of a / b over the targets when a and b may each be off by half a printed unit."""
    r = a / b
    return float(np.max(np.maximum((a + ea) / (b - eb) - r, r - (a - ea) / (b + eb))))


def rng_for(seed: int, i_comparison: int, i_bin: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([int(seed), int(i_comparison), int(i_bin)]))


def side_of_one(ci_log2) -> str:
    """Where a 95 % interval (log2) lies relative to a ratio of 1."""
    lo, hi = ci_log2
    return "above" if lo > 0 else ("below" if hi < 0 else "spans")


# ----------------------------------------------------------------------------------------------
def comparisons(summary: dict) -> list[dict]:
    """The comparisons in output order: key, numerator and denominator columns, the stored bins and how they store it."""
    out = []
    L, B = summary["log2_ratio_vs_depth"], summary["bridge_to_sphere"]
    for a in ARRAYS:
        for c in CONDITIONS:
            key = f"{a}/combined/{c}"
            out.append(dict(key=key, kind="detectability", num=f"detect_{a}_opm_{c}", den=f"detect_squid_combined_{c}",
                            stored=L[key], stored_path=f"{SUMMARY} :: log2_ratio_vs_depth['{key}']",
                            description=(f"{ARRAY_TEXT[a]} over Neuromag (306 channels), known-topography detectability with "
                                         f"the oracle covariance, {COND_TEXT[c]}")))
    for a in ARRAYS:
        out.append(dict(key=f"peak_field/{a}", kind="peak_field", num=f"amp_{a}", den="amp_squid_mag",
                        stored=B[f"{a}_ratio_vs_depth"], stored_path=f"{SUMMARY} :: bridge_to_sphere.{a}_ratio_vs_depth",
                        description=(f"peak |B| of the {ARRAY_TEXT[a]} over the peak |B| of the best Neuromag magnetometer "
                                     "(signal only; 10-nAm cortical-normal dipoles)")))
    return out


def analyse(t: dict, summary: dict, seed: int) -> dict:
    """Every comparison and bin: counts, the median and its interval, and the comparison with the stored median."""
    edges = G2.DEPTH_EDGES
    stored_bins = [(b["lo"], b["hi"]) for b in summary["log2_ratio_vs_depth"]["opm_dense/combined/intrinsic+brain"]]
    if stored_bins != [(float(lo), float(hi)) for lo, hi in zip(edges[:-1], edges[1:])]:
        raise ValueError(f"{SUMMARY}: the stored depth bins are not g2_adult_comparison.DEPTH_EDGES")
    out, worst = {}, dict(log2=(0.0, ""), log2_over_bound=(0.0, ""), linear=(0.0, ""), linear_over_bound=(0.0, ""),
                          geometric_vs_linear=(0.0, ""))
    for ic, cmp_ in enumerate(comparisons(summary)):
        key, a, b = cmp_["key"], t[cmp_["num"]], t[cmp_["den"]]
        ea, eb = t[f"{cmp_['num']}/half_unit"], t[f"{cmp_['den']}/half_unit"]
        if [(x["lo"], x["hi"]) for x in cmp_["stored"]] != stored_bins:
            raise ValueError(f"{cmp_['stored_path']}: bins differ from log2_ratio_vs_depth")
        rows = []
        for ib, st in enumerate(cmp_["stored"]):
            m = t["bin"] == label(st["lo"], st["hi"])
            n = int(m.sum())
            if n != st["n"]:
                raise ValueError(f"{cmp_['stored_path']}[{ib}]: {n} targets in the bin, {st['n']} stored")
            row = dict(lo=st["lo"], hi=st["hi"], n=n, n_parcels=int(len(np.unique(t["region"][m]))))
            if (st["median"] is None) != (n < MIN_N):
                raise ValueError(f"{cmp_['stored_path']}[{ib}]: a stored median where n < {MIN_N}, or none where n >= {MIN_N}")
            if st["median"] is not None:
                c = G2.compare(a[m], b[m], rng_for(seed, ic, ib), N_BOOT, groups=t["region"][m])
                if c["n"] != n or c["ci_method"] != f"parcels ({row['n_parcels']})":
                    raise ValueError(f"{key} bin {ib}: compare() used {c['n']} targets and {c['ci_method']}")
                where = f"{key} {label(st['lo'], st['hi'])} mm"
                # the bin's ratio is the stored median (as published); the CSV's own median is its check, the
                # interval comes from resampling the CSV's values
                if cmp_["kind"] == "detectability":
                    bound = precision_bound_log2(a[m], ea[m], b[m], eb[m])
                    diff = abs(c["median_log2"] - st["median"])
                    row.update(ratio=2.0 ** st["median"], median_log2=st["median"], median_log2_csv=c["median_log2"],
                               abs_diff_log2=diff, precision_bound_log2=bound)
                    if diff > bound + 1e-12:
                        raise ValueError(f"{where}: median {c['median_log2']!r} differs from the stored {st['median']!r} by "
                                         f"more than the CSV's precision allows ({bound:.2e})")
                    worst["log2"] = max(worst["log2"], (diff, where))
                    worst["log2_over_bound"] = max(worst["log2_over_bound"], (diff / bound, where))
                else:  # the stored statistic is the median of the linear ratio
                    lin = float(np.median(a[m] / b[m]))
                    bound = precision_bound_ratio(a[m], ea[m], b[m], eb[m])
                    diff = abs(lin - st["median"])
                    gl = abs(2.0 ** c["median_log2"] / lin - 1.0)
                    row.update(ratio=st["median"], median_ratio_linear_csv=lin, median_log2_csv=c["median_log2"],
                               abs_diff_ratio=diff, precision_bound_ratio=bound, rel_diff_log_vs_linear_median=gl)
                    if diff > bound + 1e-12:
                        raise ValueError(f"{where}: linear median {lin!r} differs from the stored {st['median']!r} by more "
                                         f"than the CSV's precision allows ({bound:.2e})")
                    worst["linear"] = max(worst["linear"], (diff, where))
                    worst["linear_over_bound"] = max(worst["linear_over_bound"], (diff / bound, where))
                    worst["geometric_vs_linear"] = max(worst["geometric_vs_linear"], (gl, where))
                row.update(ci95=c["ci95"], ci95_ratio=[2.0 ** x for x in c["ci95"]], ci_vs_1=side_of_one(c["ci95"]))
                row["ratio_inside_ci"] = bool(row["ci95_ratio"][0] <= row["ratio"] <= row["ci95_ratio"][1])
            rows.append(row)
        pop = [r for r in rows if "ratio" in r]
        out[key] = dict(description=cmp_["description"], stored=cmp_["stored_path"],
                        numerator=f"{TARGETS} :: {cmp_['num']}", denominator=f"{TARGETS} :: {cmp_['den']}", bins=rows,
                        bins_ci_above_1=[label(r["lo"], r["hi"]) for r in pop if r["ci_vs_1"] == "above"],
                        bins_ci_spans_1=[label(r["lo"], r["hi"]) for r in pop if r["ci_vs_1"] == "spans"],
                        bins_ci_below_1=[label(r["lo"], r["hi"]) for r in pop if r["ci_vs_1"] == "below"])
        log(f"{key}: " + ", ".join(f"{label(r['lo'], r['hi'])} {r['ratio']:.3f} [{r['ci95_ratio'][0]:.3f}, "
                                   f"{r['ci95_ratio'][1]:.3f}]" for r in pop))
    return dict(comparisons=out, worst=worst)


def main():
    t0 = time.time()
    cfg = tomllib.loads((ROOT / CONFIG).read_text())
    seed = cfg["sources"]["seed"]
    summary = json.loads((ROOT / SUMMARY).read_text())
    t = load_targets()
    regions = np.unique(t["region"])
    n_labels = int(summary["primary"]["oracle"]["opm_dense/combined/intrinsic+brain"]["ci_method"].split("(")[1].rstrip(")"))
    if len(regions) != n_labels or len(t["region"]) != summary["n_targets"]:
        raise ValueError(f"{TARGETS}: {len(regions)} region labels and {len(t['region'])} targets; {SUMMARY} has {n_labels} "
                         f"parcel labels and {summary['n_targets']} targets")
    moved = int(np.sum(np.array([label(lo, hi) for lo, hi in zip(G2.DEPTH_EDGES[:-1], G2.DEPTH_EDGES[1:])], dtype=object)
                       [np.clip(np.searchsorted(G2.DEPTH_EDGES, t["depth_printed"], side="right") - 1, 0,
                                len(G2.DEPTH_EDGES) - 2)] != t["bin"]))
    res = analyse(t, summary, seed)
    w = res["worst"]
    out = dict(
        status=STATUS,
        inputs=dict(targets=TARGETS, targets_status=status_line(TARGETS), depth=DEPTH, depth_status=status_line(DEPTH),
                    summary=SUMMARY, summary_commit=summary["provenance"]["commit"], config=CONFIG),
        method=dict(
            statistic=("per bin, the median over the bin's targets of log2(m_OPM / m_Neuromag) (m: detectability or peak "
                       f"|B|). 'ratio' is the median stored in {SUMMARY} (2**median_log2 of log2_ratio_vs_depth; the linear "
                       "median of bridge_to_sphere for the peak field), as published; median_log2_csv (and "
                       "median_ratio_linear_csv) is the same median recomputed from the CSV, the check; the interval is "
                       "the bootstrap of the CSV's values, ci95 in log2, ci95_ratio = 2**ci95"),
            interval=(f"95 %: 2.5th and 97.5th percentiles of the medians of {N_BOOT} bootstrap resamples of the parcel labels "
                      "present in the bin (drawn with replacement, their targets in the bin pooled); "
                      "scripts/g2_adult_comparison.py compare(), the adult analysis's own function and count"),
            function="scripts/g2_adult_comparison.py :: compare", n_boot=N_BOOT, unit="parcel labels present in the bin",
            n_parcel_labels=n_labels, seed=seed,
            seed_rule=(f"numpy default_rng(SeedSequence([{seed} ({CONFIG} sources.seed), comparison index, bin index])), the "
                       "comparisons indexed in the order of 'comparisons'"),
            bins=dict(edges_mm=[float(x) for x in G2.DEPTH_EDGES], rule=(f"lo <= depth < hi on the unrounded depth ({DEPTH} "
                                                                          "depth_bin_5mm)"),
                      min_n=MIN_N, note=f"bins with fewer than {MIN_N} targets have no median, as in {SUMMARY}"),
            precision=(f"{TARGETS} prints detectability with 4 decimals and peak fields with 5 significant digits; each "
                       "median is checked against the stored one within the largest change that rounding the bin's "
                       "printed values can make (precision_bound_*)"),
            peak_field=("bridge_to_sphere stores the median of the linear ratio, the bootstrap resamples the median of "
                        "log2 (the two differ only when a bin has an even count: rel_diff_log_vs_linear_median)")),
        comparisons=res["comparisons"],
        checks=dict(
            keys_and_order_equal=True, printed_depth_equals_rounded_unrounded=True, bin_labels_follow_unrounded_depth=True,
            bin_counts_equal_stored=True, medians_within_csv_precision=True,
            targets_whose_printed_depth_falls_in_another_bin=moved,
            max_abs_diff_median_log2=w["log2"][0], max_abs_diff_median_log2_at=w["log2"][1],
            max_diff_over_precision_bound_log2=w["log2_over_bound"][0], max_diff_over_precision_bound_log2_at=w["log2_over_bound"][1],
            max_abs_diff_peak_field_median=w["linear"][0], max_abs_diff_peak_field_median_at=w["linear"][1],
            max_diff_over_precision_bound_peak_field=w["linear_over_bound"][0],
            max_rel_diff_log_vs_linear_median_peak_field=w["geometric_vs_linear"][0],
            max_rel_diff_log_vs_linear_median_peak_field_at=w["geometric_vs_linear"][1]),
        runtime_s=time.time() - t0,
    )
    io.write_json(out, ROOT / OUT)
    log(f"largest deviation of a median from the stored one: {w['log2'][0]:.2e} (log2, at {w['log2'][1]}); peak field "
        f"{w['linear'][0]:.2e} (at {w['linear'][1]}); {moved} targets print into another bin; wrote {OUT}")


if __name__ == "__main__":
    main()
