#!/usr/bin/env python3
"""Report facts for G4 (merged by scripts/report_facts.py): simulated interictal-spike detection
(g4_), bounded localization (loc_) and head motion (mot_), read from the committed results/g4
files and the configuration the runs recorded.

facts(root) returns name -> {"value": text as printed, "raw": unrounded number(s) or text,
"source": "results/<file> :: <key path>" or "results/<file> :: derived: <how>"}. Derived values
are ranges, counts and Holm adjustments over stored entries and tallies of the stored per-event
tables; nothing is simulated or re-estimated.

Formats: ratios 2 decimals (plus a _3dp twin when 0.95-1.05); dB and paired differences signed,
U+2212 for minus; mm 1 decimal; nAm, fT, pT and counts integer with thousands separators; shares,
sensitivities and false-event rates 2 decimals; percentages integer (1 decimal below 10 %);
p-values 2 significant digits down to 0.0001, "<0.0001" below; rotation thresholds 2 significant
digits; configuration constants as declared. Intervals "[lo, hi]" ("open" for an open end), ranges
"lo to hi" from unrounded values (one value when both ends print alike). Rounding is Python's
(round half to even on exact binary ties, e.g. 0.125 -> 0.12), as in the repository's reports.
Names: band tags 10_20mm, 20_30mm, 30_45mm, 45_70mm (depth below the scalp); Neuromag 306 channels
"combined", OPM arrays "dense" and "matched"; "practical" = scanning detector at 1 false event/min.

Usage: .venv/bin/python scripts/report_facts_g4.py [--grep TEXT]   (prints name, value, source)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MINUS = "\u2212"
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
TEMPLATES_CONTROLS = ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")
CHILDREN = ("childA", "childB", "childC")
# the referees' split of the smaller heads (round 3): the principal pediatric evidence (the scaled adults and the 24- and
# 12-month templates; with the adult where a statement includes it) and the provisional heads reported separately (the
# 18-month template, classed misregistered by the MRI check, and children A-C)
PRINCIPAL = ("school", "size2yr", "infant2yr", "infant12mo")
PROVISIONAL = ("infant18mo", "childA", "childB", "childC")
SPLIT = (("principal", PRINCIPAL, "the four principal heads"),
         ("adult_principal", ("adult",) + PRINCIPAL, "the adult and the four principal heads"),
         ("provisional", PROVISIONAL, "the four provisional heads"))
BANDS = ("10_20mm", "20_30mm", "30_45mm", "45_70mm")  # depth0..depth3
BAND_EDGES = (10, 20, 30, 45, 70)  # docs/methods.md section 9
DETECTORS = {"squid/combined": "combined", "squid/grad": "grad", "squid/mag": "mag", "opm_matched/opm": "matched",
             "opm_dense/opm": "dense"}
OPM = {"opm_dense/opm": "dense", "opm_matched/opm": "matched"}
MODES = {"practical@1": "practical", "oracle": "oracle"}
ARRAYS = {"squid": "combined", "opm_matched": "matched", "opm_dense": "dense"}  # localization views (primary)
CONDS = ("focal/80nAm", "focal/320nAm", "patch/80nAm", "patch/320nAm")
STRENGTH_WORDS = {10.0: "10", 20.0: "20", 40.0: "40", 80.0: "80", 160.0: "160", 320.0: "320"}


# ----------------------------------------------------------------------------------------------
# formatting
def _fmt(x, nd: int, sep: bool = False, sign: bool = False) -> str:
    s = format(x, f"{',' if sep else ''}.{nd}f")
    if s.lstrip("-").strip("0.,") == "":  # rounds to zero: no sign
        return s.lstrip("-")
    if sign and not s.startswith("-"):
        s = "+" + s
    return s.replace("-", MINUS)


def ratio(x) -> str:
    return _fmt(x, 2)


def ratio3(x) -> str:
    return _fmt(x, 3)


def db(x) -> str:
    return _fmt(x, 2, sign=True)


def dbu(x) -> str:
    """An unsigned dB quantity (a loss: the change with its sign reversed)."""
    return _fmt(x, 2)


def mm(x) -> str:
    return _fmt(x, 1)


def dmm(x) -> str:
    return _fmt(x, 1, sign=True)


def integer(x) -> str:
    return _fmt(x, 0, sep=True)


def count(n) -> str:
    return f"{int(n):,d}"


def share(x) -> str:
    return _fmt(x, 2)


def pct(frac) -> str:
    v = 100.0 * frac
    return (_fmt(v, 1) if float(format(abs(v), ".1f")) < 10 else _fmt(v, 0, sep=True)) + "%"


def sig(x, n: int = 2) -> str:
    """n significant digits in plain decimals (integers keep all their digits)."""
    if x == 0:
        return "0"
    e = int(format(abs(x), f".{n - 1}e").split("e")[1])
    return _fmt(x, max(n - 1 - e, 0), sep=True)


def pval(p) -> str:
    return "<0.0001" if p < 1e-4 else sig(p, 2)


def _join(items: list) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def const(x) -> str:
    """A configuration constant as declared (no added decimals); a list as 'a, b and c'."""
    if isinstance(x, (list, tuple)):
        return _join([const(v) for v in x])
    if isinstance(x, float) and x.is_integer():
        x = int(x)
    return (f"{x:,d}" if isinstance(x, int) else format(x, "g")).replace("-", MINUS)


def iv(lo, hi, f) -> str:
    return "[" + ("open" if lo is None else f(lo)) + ", " + ("open" if hi is None else f(hi)) + "]"


def span(values, f) -> str:
    a, b = f(min(values)), f(max(values))
    return a if a == b else f"{a} to {b}"


# ----------------------------------------------------------------------------------------------
def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def _csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def _add(out: dict, name: str, value: str, raw, source: str) -> None:
    if name in out:
        raise ValueError(f"fact {name!r} defined twice")
    if name != name.lower() or not name.startswith(("g4_", "loc_", "mot_")) or " " in name:
        raise ValueError(f"bad fact name {name!r}")
    if str(value).strip() == "":
        raise ValueError(f"fact {name!r} has an empty value")
    out[name] = dict(value=value, raw=raw, source=source)


def _spread(out: dict, name: str, values: list, f, source: str) -> None:
    """name_min, name_max and name_range (from unrounded values)."""
    _add(out, f"{name}_min", f(min(values)), min(values), f"{source} (minimum)")
    _add(out, f"{name}_max", f(max(values)), max(values), f"{source} (maximum)")
    _add(out, f"{name}_range", span(values, f), [min(values), max(values)], f"{source} (minimum to maximum)")


def _tag(lab: str) -> str:
    return lab.lower()


def _garwood(k: int, minutes: float) -> list:
    """Exact Poisson 95 % interval of a rate (as opmsquid.detection.rate_with_ci)."""
    from scipy.stats import chi2

    lo = 0.0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2.0
    return [lo / minutes, chi2.ppf(0.975, 2 * k + 2) / 2.0 / minutes]


def _holm(ps: dict) -> dict:
    run, out = 0.0, {}
    items = sorted(ps.items(), key=lambda kv: kv[1])
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (len(items) - i) * p))
        out[k] = run
    return out


def _ratio_fact(out, name, sr: dict, source: str) -> float | None:
    """Paired S50 ratio Neuromag / OPM (+ _ci, + _3dp twins near 1); returns the point value."""
    v, ci = sr.get("value"), sr.get("ci95") or [None, None]
    if v is not None:
        _add(out, name, ratio(v), v, f"{source}.value")
        _add(out, f"{name}_ci", iv(ci[0], ci[1], ratio), ci, f"{source}.ci95 (location bootstrap; open = censored end)")
        if 0.95 <= v <= 1.05:
            _add(out, f"{name}_3dp", ratio3(v), v, f"{source}.value")
            _add(out, f"{name}_ci_3dp", iv(ci[0], ci[1], ratio3), ci, f"{source}.ci95")
        return v
    lo, hi = sr.get("value_bounds") or [None, None]
    why = sr.get("value_censored") or "censored"
    if lo is not None and hi is None:
        txt, raw = "> " + ratio(lo), [lo, None]
    elif hi is not None and lo is None:
        txt, raw = "< " + ratio(hi), [None, hi]
    else:
        txt, raw = "not estimable", "not estimable"
    _add(out, name, txt, raw, f"{source}.value_bounds (point estimate censored: {why})")
    _add(out, f"{name}_ci", iv(ci[0], ci[1], ratio), ci, f"{source}.ci95 (location bootstrap; open = censored end)")
    return None


# ----------------------------------------------------------------------------------------------
# detection
def _detection(root: Path, out: dict) -> None:
    det_file = {lab: f"results/g4/g4_{lab}_summary.json" for lab in LABELS}
    S = {lab: _json(root / det_file[lab]) for lab in LABELS}
    CMPF = "results/g4/g4_pediatric_comparison.json"
    C = _json(root / CMPF)["comparison"]
    MRF = "results/g4/g4_matched_rate.json"
    MR = _json(root / MRF)["anatomies"]
    CONSF = "results/g4/g4_g2_consistency.json"
    CONS = _json(root / CONSF)["bands"]
    A = "results/g4/g4_adult_summary.json"
    cfg = S["adult"]["config"]
    for lab in LABELS:  # the same detection configuration in every anatomy
        for sec in ("simulation", "null", "events", "detector"):
            if S[lab]["config"][sec] != cfg[sec]:
                raise ValueError(f"{lab}: detection configuration [{sec}] differs from the adult's")
    same = " (identical in all nine g4_<anatomy>_summary.json)"
    sim, nul, ev, dc = cfg["simulation"], cfg["null"], cfg["events"], cfg["detector"]

    # configuration constants
    for name, val, key in (
            ("g4_band_lo_hz", sim["band_hz"][0], "simulation.band_hz[0]"),
            ("g4_band_hi_hz", sim["band_hz"][1], "simulation.band_hz[1]"),
            ("g4_decimation_factor", sim["decimate"], "simulation.decimate"),
            ("g4_segment_s", sim["segment_s"], "simulation.segment_s"),
            ("g4_event_spacing_s", sim["event_spacing_s"], "simulation.event_spacing_s"),
            ("g4_opm_noise_ft_per_rthz", sim["opm_asd_fT_per_rtHz"], "simulation.opm_asd_fT_per_rtHz (white, fT/sqrt(Hz))"),
            ("g4_simulation_seed", sim["seed"], "simulation.seed"),
            ("g4_null_baseline_minutes", nul["baseline_min"], "null.baseline_min (whitener)"),
            ("g4_null_calibration_minutes", nul["calibration_min"], "null.calibration_min (thresholds, then frozen)"),
            ("g4_null_heldout_minutes", nul["heldout_min"], "null.heldout_min (false-event check)"),
            ("g4_n_locations", ev["n_locations"], "events.n_locations"),
            ("g4_strengths_nam", ev["strengths_nAm"], "events.strengths_nAm"),
            ("g4_strength_min_nam", ev["strengths_nAm"][0], "events.strengths_nAm[0]"),
            ("g4_strength_max_nam", ev["strengths_nAm"][-1], "events.strengths_nAm[-1]"),
            ("g4_n_strengths", len(ev["strengths_nAm"]), "events.strengths_nAm (count)"),
            ("g4_stretches", ev["stretches"], "events.stretches (spike-wave duration scale)"),
            ("g4_n_morphologies", len(ev["stretches"]), "events.stretches (count)"),
            ("g4_patch_radius_mm", ev["patch_radius_mm"], "events.patch_radius_mm"),
            ("g4_dictionary_spacing_mm", dc["dictionary_spacing_mm"], "detector.dictionary_spacing_mm"),
            ("g4_refractory_s", dc["refractory_s"], "detector.refractory_s"),
            ("g4_operating_point_per_minute", dc["operating_points_per_min"][0], "detector.operating_points_per_min[0]"),
            ("g4_operating_point_low_per_minute", dc["operating_points_per_min"][1], "detector.operating_points_per_min[1]"),
            ("g4_oracle_alpha", dc["oracle_alpha"], "detector.oracle_alpha (per-trial false-positive probability)"),
            ("g4_n_detectors", len(dc["arrays"]), "detector.arrays (count: Neuromag combined, grad, mag; OPM matched, dense)")):
        _add(out, name, str(val) if name.endswith("_seed") else const(val), val, f"{A} :: config.{key}{same}")
    _add(out, "g4_strength_range", f"{const(ev['strengths_nAm'][0])} to {const(ev['strengths_nAm'][-1])}",
         [ev["strengths_nAm"][0], ev["strengths_nAm"][-1]], f"{A} :: config.events.strengths_nAm (first to last)")
    for i, s_nam in enumerate(ev["strengths_nAm"]):  # each tested strength, as printed beside per-strength counts
        _add(out, f"g4_strength_{const(s_nam)}_nam", const(s_nam), s_nam, f"{A} :: config.events.strengths_nAm[{i}]{same}")
    steps = {b / a for a, b in zip(ev["strengths_nAm"][:-1], ev["strengths_nAm"][1:])}
    if len(steps) == 1:
        _add(out, "g4_strength_step_factor", const(steps.pop()), ev["strengths_nAm"],
             f"{A} :: derived: ratio of successive config.events.strengths_nAm")
    _add(out, "g4_hit_tolerance_ms", const(round(dc["hit_tolerance_s"] * 1000, 6)), dc["hit_tolerance_s"],
         f"{A} :: config.detector.hit_tolerance_s x 1000 (+/- around the true peak)")
    _add(out, "g4_sampling_rate_hz", integer(S["adult"]["fs_out"]), S["adult"]["fs_out"], f"{A} :: fs_out (after decimation)")
    per_seg = int((sim["segment_s"] - 1.5) // sim["event_spacing_s"])
    _add(out, "g4_events_per_segment", count(per_seg), per_seg,
         f"{A} :: derived: floor((config.simulation.segment_s - 1.5 s first event) / config.simulation.event_spacing_s) "
         "(configs/g4_epilepsy.toml: 'first at 1.5 s')")
    n_ev = S["adult"]["n_events"]
    if {S[lab]["n_events"] for lab in LABELS} != {n_ev}:
        raise ValueError("n_events differs between anatomies")
    n_focal = ev["n_locations"] * len(ev["strengths_nAm"]) * len(ev["stretches"])
    _add(out, "g4_n_events_per_anatomy", count(n_ev), n_ev, f"{A} :: n_events{same}")
    _add(out, "g4_n_focal_events", count(n_focal), n_focal,
         f"{A} :: derived: config.events n_locations x strengths x stretches (focal dipoles)")
    _add(out, "g4_n_patch_events", count(n_ev - n_focal), n_ev - n_focal, f"{A} :: derived: n_events - focal events (10-mm patches)")
    _add(out, "g4_n_focal_events_per_location", count(len(ev["strengths_nAm"]) * len(ev["stretches"])),
         len(ev["strengths_nAm"]) * len(ev["stretches"]), f"{A} :: derived: strengths x stretches")
    p0 = S["adult"]["paired"]["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]
    if {(r["n"], r["n_locations"]) for lab in LABELS for blk in S[lab]["paired"].values() for r in blk.values()} != {(p0["n"], p0["n_locations"])}:
        raise ValueError("paired comparisons differ in their event or location counts")
    _add(out, "g4_n_focal_events_per_band", count(p0["n"]), p0["n"],
         f"{A} :: paired['opm_dense/opm_vs_squid/combined/practical@1'].depth0.n (every band and anatomy)")
    per_band = {v for lab in LABELS for v in C[f"{lab}/n_locations_per_depth_band"].values()}
    if len(per_band) != 1:
        raise ValueError(f"locations per depth band differ between anatomies: {per_band}")
    n_pb = per_band.pop()
    _add(out, "g4_locations_per_band", count(n_pb), n_pb,
         f"{CMPF} :: comparison['<anatomy>/n_locations_per_depth_band'] (the same in every band of all nine anatomies)")
    strata = {tuple(L["stratum"]) for L in S["adult"]["locations"]}
    n_or = len({s[1] for s in strata})
    _add(out, "g4_n_depth_bands", count(len({s[0] for s in strata})), len({s[0] for s in strata}),
         f"{A} :: derived: distinct locations[].stratum[0]")
    _add(out, "g4_n_orientation_bands", count(n_or), n_or, f"{A} :: derived: distinct locations[].stratum[1] "
         "(orientation 0-30, 30-60, 60-90 deg to the scalp normal: scripts/g4_epilepsy_adult.py ORIENT_BANDS)")
    _add(out, "g4_n_strata", count(len(strata)), len(strata), f"{A} :: derived: distinct locations[].stratum")
    _add(out, "g4_locations_per_stratum", count(ev["n_locations"] // len(strata)), ev["n_locations"] // len(strata),
         f"{A} :: derived: config.events.n_locations / strata")
    for b, (lo, hi) in enumerate(zip(BAND_EDGES[:-1], BAND_EDGES[1:])):
        depths = [L["depth_mm"] for L in S["adult"]["locations"] if L["stratum"][0] == b]
        if not all(lo <= d < hi for d in depths):
            raise ValueError(f"depth band {b} edges do not match the stored locations")
        _add(out, f"g4_band_{BANDS[b]}", f"{lo}\u2013{hi}", [lo, hi],
             "docs/methods.md :: section 9 'stratified by depth (10-20, 20-30, 30-45, 45-70 mm)' (depth below the scalp; "
             "consistent with results/g4/g4_adult_summary.json locations[].depth_mm)")
        _add(out, f"g4_band_{BANDS[b]}_lo", const(lo), lo, "docs/methods.md :: section 9 depth bands")
        _add(out, f"g4_band_{BANDS[b]}_hi", const(hi), hi, "docs/methods.md :: section 9 depth bands")
    _add(out, "g4_s50_level", "50%", 0.5, f"{A} :: strength_for_50pct_nAm (S50 = strength for 50 % detection, log interpolation)")
    _add(out, "g4_bootstrap_resamples", "1,000", 1000, "docs/methods.md :: section 9 'a bootstrap over locations (1,000 resamples)'")
    for op, k in ((dc["operating_points_per_min"][0], "1pm"), (dc["operating_points_per_min"][1], "02pm")):
        n_cal = int(math.floor(op * nul["calibration_min"]))
        _add(out, f"g4_threshold_calibration_events_{k}", count(n_cal), n_cal,
             f"{A} :: derived: floor(rate x config.null.calibration_min) calibration events above each frozen threshold "
             "(opmsquid.detection.threshold_for_rate)")
        ci = _garwood(n_cal, nul["calibration_min"])
        _add(out, f"g4_threshold_rate_poisson_ci_{k}", iv(ci[0], ci[1], share), ci,
             f"{A} :: derived: exact (Garwood) Poisson 95 % interval of {n_cal} events in config.null.calibration_min, per minute")
    # the paired family per anatomy
    pk = [k for k in S["adult"]["paired"] if k.startswith("opm_dense/")]
    n_fam = len(pk) * len(BANDS)
    _add(out, "g4_family_per_opm_array", count(n_fam), n_fam,
         f"{A} :: derived: paired keys per OPM array (3 Neuromag comparators x 2 detectors) x 4 depth bands")
    _add(out, "g4_family_per_anatomy", count(2 * n_fam), 2 * n_fam, f"{A} :: derived: both OPM arrays x {n_fam}")
    _add(out, "g4_bonferroni_24_threshold", pval(0.05 / n_fam), 0.05 / n_fam, f"{A} :: derived: 0.05 / {n_fam} comparisons")
    _add(out, "g4_bonferroni_48_threshold", pval(0.05 / (2 * n_fam)), 0.05 / (2 * n_fam), f"{A} :: derived: 0.05 / {2 * n_fam}")
    _add(out, "g4_n_anatomies", count(len(LABELS)), len(LABELS), f"{CMPF} :: labels (count)")
    _add(out, "g4_n_smaller_heads", count(len(LABELS) - 1), len(LABELS) - 1, f"{CMPF} :: labels (count without the adult)")
    nd = {lab: S[lab]["n_dictionary"] for lab in LABELS}
    for lab in LABELS:
        _add(out, f"g4_{_tag(lab)}_n_dictionary", count(nd[lab]), nd[lab], f"{det_file[lab]} :: n_dictionary (scanner candidates)")
    _spread(out, "g4_n_dictionary", list(nd.values()), count, f"{', '.join(det_file.values())} :: n_dictionary over the nine anatomies")

    # adult pose (G3B placements): the G4 adult sits at its measured position, the smaller heads at top contact
    G3F = "results/g3b/g3b_summary.json"
    g3 = _json(root / G3F)
    d_meas = g3["placement_D"]["adult/centred/combined/intrinsic+brain"]["median"]
    d_top = g3["placement_D"]["adult/top/combined/intrinsic+brain"]["median"]
    _add(out, "g4_adult_d_measured_position_db", db(d_meas), d_meas,
         f"{G3F} :: placement_D['adult/centred/combined/intrinsic+brain'].median (placements.adult.centred = the measured "
         "device-to-head transform, as in G4)")
    _add(out, "g4_adult_d_top_contact_db", db(d_top), d_top, f"{G3F} :: placement_D['adult/top/combined/intrinsic+brain'].median")
    _add(out, "g4_adult_pose_d_difference_db", db(d_meas - d_top), d_meas - d_top,
         f"{G3F} :: derived: placement_D adult centred (measured) minus top, combined, intrinsic+brain medians")
    mv = {lab: g3["placements"][lab]["top"]["moved_mm"] for lab in LABELS}
    _add(out, "g4_adult_top_contact_shift_mm", mm(mv["adult"]), mv["adult"], f"{G3F} :: placements.adult.top.moved_mm")
    _spread(out, "g4_smaller_heads_top_contact_shift_mm", [mv[lab] for lab in LABELS[1:]], mm,
            f"{G3F} :: placements.<anatomy>.top.moved_mm over the eight smaller heads")

    # per anatomy: S50, paired comparisons, held-out rates, sensitivities, locations
    for lab in LABELS:
        s, f, t = S[lab], det_file[lab], _tag(lab)
        for d, dt in DETECTORS.items():
            for mode, mt in MODES.items():
                for b, bt in enumerate(BANDS):
                    v = s["detectors"][d]["strength_for_50pct_nAm"][f"{mode}/depth{b}"]
                    key = f"{f} :: detectors['{d}'].strength_for_50pct_nAm['{mode}/depth{b}']"
                    name = f"g4_{t}_s50_{dt}_{mt}_{bt}"
                    if v["value"] is None:
                        _add(out, name, "not reached", "not reached", f"{key}.value (null: 50 % not reached by 320 nAm)")
                    else:
                        _add(out, name, integer(v["value"]), v["value"], f"{key}.value")
                    _add(out, f"{name}_ci", iv(v["ci95"][0], v["ci95"][1], integer), v["ci95"],
                         f"{key}.ci95 (location bootstrap; open = beyond the tested strengths)")
            hr = s["detectors"][d]["heldout_false_per_min"]
            _add(out, f"g4_{t}_heldout_rate_{dt}", share(hr["1"]), hr["1"],
                 f"{f} :: detectors['{d}'].heldout_false_per_min['1'] (per minute, frozen 1-per-minute threshold)")
            _add(out, f"g4_{t}_heldout_rate_{dt}_op02", share(hr["0.2"]), hr["0.2"],
                 f"{f} :: detectors['{d}'].heldout_false_per_min['0.2'] (frozen 0.2-per-minute threshold)")
            sv = C[f"{lab}/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min"]
            _add(out, f"g4_{t}_sens40_matchedrate_{dt}", share(sv), sv,
                 f"{CMPF} :: comparison['{lab}/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min'] (focal 40 nAm at "
                 "10-30 mm, threshold set on the held-out null to 1 false event/min)")
        for a, at in OPM.items():
            for mode, mt in MODES.items():
                for b, bt in enumerate(BANDS):
                    r = s["paired"][f"{a}_vs_squid/combined/{mode}"][f"depth{b}"]
                    key = f"{f} :: paired['{a}_vs_squid/combined/{mode}'].depth{b}"
                    base = f"g4_{t}_{at}_vs_combined_{mt}_{bt}"
                    _add(out, f"{base}_locs_opm", count(r["locations_favouring_opm"]), r["locations_favouring_opm"],
                         f"{key}.locations_favouring_opm (of {r['n_locations']})")
                    _add(out, f"{base}_locs_squid", count(r["locations_favouring_squid"]), r["locations_favouring_squid"],
                         f"{key}.locations_favouring_squid (of {r['n_locations']})")
                    _add(out, f"{base}_p", pval(r["location_sign_flip_p"]), r["location_sign_flip_p"],
                         f"{key}.location_sign_flip_p (exact sign-flip over locations, uncorrected)")
                    _add(out, f"{base}_events_opm_only", count(r["detected_only_opm"]), r["detected_only_opm"],
                         f"{key}.detected_only_opm (of {r['n']} focal events)")
                    _add(out, f"{base}_events_squid_only", count(r["detected_only_squid"]), r["detected_only_squid"],
                         f"{key}.detected_only_squid (of {r['n']} focal events)")
                    rv = _ratio_fact(out, f"{base}_ratio", r["s50_ratio_squid_over_opm"], f"{key}.s50_ratio_squid_over_opm")
                    if rv is not None:
                        _add(out, f"{base}_reduction_pct", pct(1 - 1 / rv), 1 - 1 / rv,
                             f"{key}.s50_ratio_squid_over_opm.value; derived: 1 - 1/ratio (the OPM's S50 is lower by this share)")
        # sampled locations
        for b, bt in enumerate(BANDS):
            dep = sorted(L["depth_mm"] for L in s["locations"] if L["stratum"][0] == b)
            med = _median(dep)
            _add(out, f"g4_{t}_median_depth_{bt}", mm(med), med, f"{f} :: derived: median locations[].depth_mm with stratum[0] = {b}")
        lobes = {}
        for L in s["locations"]:
            lobes[L["lobe"]] = lobes.get(L["lobe"], 0) + 1
        for lobe in ("frontal", "temporal", "parietal", "occipital", "cingulate", "insula"):
            _add(out, f"g4_{t}_locations_{lobe}", count(lobes.get(lobe, 0)), lobes.get(lobe, 0),
                 f"{f} :: derived: count of locations[].lobe == '{lobe}'")

    # frozen-threshold sensitivities from the per-event tables (identical events in every array)
    n_sens = set()
    for lab in LABELS:
        f = f"results/g4/g4_{lab}_events.csv"
        rows = _csv(root / f)
        thr = {d: S[lab]["detectors"][d]["thresholds"]["1"] for d in DETECTORS}
        sup = [r for r in rows if r["family"] == "focal" and int(r["depth_band"]) <= 1 and float(r["strength_nAm"]) == 40.0]
        deep = [r for r in rows if r["family"] == "focal" and int(r["depth_band"]) >= 2 and float(r["strength_nAm"]) == 160.0]
        n_sens |= {len(sup), len(deep)}
        for d, dt in DETECTORS.items():
            for sel, nm, what in ((sup, "sens40_frozen", "focal 40 nAm at 10-30 mm"), (deep, "sens160deep_frozen", "focal 160 nAm at 30-70 mm")):
                v = sum(float(r[f"{d}_event_height"]) > thr[d] for r in sel) / len(sel)
                _add(out, f"g4_{_tag(lab)}_{nm}_{dt}", share(v), v,
                     f"{f} :: derived: share of {len(sel)} {what} events with {d}_event_height > detectors['{d}'].thresholds['1'] "
                     f"of {det_file[lab]}")
    if len(n_sens) != 1:
        raise ValueError(f"sensitivity event counts differ: {n_sens}")
    n_sens = n_sens.pop()
    _add(out, "g4_sens_n_events", count(n_sens), n_sens, "results/g4/g4_adult_events.csv :: derived: focal 40-nAm events at "
         "10-30 mm, and 160-nAm events at 30-70 mm (the same count in every anatomy)")

    # adult superficial events: who detects what, by strength (practical detector, dense vs combined)
    f = "results/g4/g4_adult_events.csv"
    rows = _csv(root / f)
    thr = {d: S["adult"]["detectors"][d]["thresholds"]["1"] for d in DETECTORS}
    sel = [r for r in rows if r["family"] == "focal" and int(r["depth_band"]) == 0]
    det = {d: [float(r[f"{d}_event_height"]) > thr[d] for r in sel] for d in ("opm_dense/opm", "squid/combined")}
    only_o = [float(r["strength_nAm"]) for r, x, y in zip(sel, det["opm_dense/opm"], det["squid/combined"]) if x and not y]
    only_s = [float(r["strength_nAm"]) for r, x, y in zip(sel, det["opm_dense/opm"], det["squid/combined"]) if y and not x]
    jp = S["adult"]["paired"]["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]
    if (len(only_o), len(only_s)) != (jp["detected_only_opm"], jp["detected_only_squid"]):
        raise ValueError("adult event table does not reproduce the stored paired counts")
    how = (f"{f} :: derived: focal events at 10-20 mm (depth_band 0), detected = event_height > detectors[...].thresholds['1'] "
           f"of {A}; reproduces paired['opm_dense/opm_vs_squid/combined/practical@1'].depth0 detected_only_opm/squid")
    _add(out, "g4_adult_events_10_20mm_n", count(len(sel)), len(sel), how)
    _add(out, "g4_adult_events_10_20mm_dense_only", count(len(only_o)), len(only_o), how)
    _add(out, "g4_adult_events_10_20mm_combined_only", count(len(only_s)), len(only_s), how)
    net = (len(only_o) - len(only_s)) / len(sel)
    _add(out, "g4_adult_events_10_20mm_net_pp", _fmt(100 * net, 1), 100 * net,
         f"{how}; (dense-only - combined-only) / n x 100, percentage points (the prose adds the unit)")
    for st, w in STRENGTH_WORDS.items():
        _add(out, f"g4_adult_events_10_20mm_dense_only_{w}nam", count(only_o.count(st)), only_o.count(st), f"{how}; strength {w} nAm")
        _add(out, f"g4_adult_events_10_20mm_combined_only_{w}nam", count(only_s.count(st)), only_s.count(st), f"{how}; strength {w} nAm")
        n_st = sum(float(r["strength_nAm"]) == st for r in sel)
        for d, dt in (("opm_dense/opm", "dense"), ("squid/combined", "combined")):
            k = sum(x for r, x in zip(sel, det[d]) if float(r["strength_nAm"]) == st)
            _add(out, f"g4_adult_events_10_20mm_detected_{dt}_{w}nam", count(k), k, f"{how}; strength {w} nAm, of {n_st}")
    _add(out, "g4_adult_events_10_20mm_per_strength", count(n_st), n_st, f"{how}; events per strength (18 locations x 3 morphologies)")
    low = [x for x in only_o if x <= 40.0]
    _add(out, "g4_adult_events_10_20mm_dense_only_20_40nam", count(len(low)), len(low), f"{how}; dense-only events at 20-40 nAm")
    _add(out, "g4_adult_events_10_20mm_dense_only_ge80nam", count(len(only_o) - len(low)), len(only_o) - len(low),
         f"{how}; dense-only events at 80 nAm or more")
    _add(out, "g4_adult_events_10_20mm_dense_only_strengths", const(sorted(set(only_o))), sorted(set(only_o)),
         f"{how}; strengths (nAm) of the dense-only events")
    _add(out, "g4_adult_events_10_20mm_combined_only_strengths", const(sorted(set(only_s))), sorted(set(only_s)),
         f"{how}; strengths (nAm) of the combined-only events")

    # adult: practical vs oracle S50 at 10-20 mm (the cost of searching time and source)
    pr = []
    for d, dt in (("squid/combined", "combined"), ("opm_matched/opm", "matched"), ("opm_dense/opm", "dense")):
        sp = S["adult"]["detectors"][d]["strength_for_50pct_nAm"]
        x = sp["practical@1/depth0"]["value"] / sp["oracle/depth0"]["value"]
        pr.append(x)
        _add(out, f"g4_adult_practical_over_oracle_s50_{dt}_10_20mm", ratio(x), x,
             f"{A} :: derived: detectors['{d}'].strength_for_50pct_nAm practical@1/depth0 / oracle/depth0 values")
    _spread(out, "g4_adult_practical_over_oracle_s50_10_20mm", pr, ratio,
            f"{A} :: derived: practical / oracle S50 at 10-20 mm, Neuromag combined, OPM matched and dense")

    # held-out rates over detectors and anatomies
    src_all = f"{CMPF} :: comparison['<anatomy>/<detector>/heldout_false_per_min']"
    for op, nm in (("1", "g4_heldout_rate"), ("0.2", "g4_heldout_rate_op02")):
        _spread(out, nm, [C[f"{lab}/{d}/heldout_false_per_min"][op] for lab in LABELS for d in DETECTORS], share,
                f"{src_all}['{op}'] over the five detectors and nine anatomies")
    _spread(out, "g4_adult_heldout_rate", [S["adult"]["detectors"][d]["heldout_false_per_min"]["1"] for d in DETECTORS], share,
            f"{A} :: detectors[*].heldout_false_per_min['1'] over the five detectors")
    _spread(out, "g4_adult_heldout_rate_op02", [S["adult"]["detectors"][d]["heldout_false_per_min"]["0.2"] for d in DETECTORS],
            share, f"{A} :: detectors[*].heldout_false_per_min['0.2'] over the five detectors")
    sens = lambda labs, d: [C[f"{lab}/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min"] for lab in labs]  # noqa: E731
    for d, dt in (("opm_dense/opm", "dense"), ("squid/combined", "combined"), ("opm_matched/opm", "matched")):
        _spread(out, f"g4_smaller_heads_sens40_matchedrate_{dt}", sens(LABELS[1:], d), share,
                f"{CMPF} :: comparison['<anatomy>/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min'] over the eight smaller heads")

    # the 10-20 mm endpoint across anatomies (dense vs combined, practical)
    k1020 = "paired/opm_dense/opm_vs_squid/combined/practical@1/depth0"
    P = {lab: C[f"{lab}/{k1020}"] for lab in LABELS}
    src = f"{CMPF} :: comparison['<anatomy>/{k1020}']"
    ps = {lab: P[lab]["location_sign_flip_p"] for lab in LABELS}
    adj = _holm(ps)
    for lab in LABELS:
        _add(out, f"g4_holm_10_20mm_{_tag(lab)}_p_adj", pval(adj[lab]), adj[lab],
             f"{src.replace('<anatomy>', lab)}.location_sign_flip_p; derived: Holm adjustment over the nine anatomies "
             "(the 10-20 mm comparison taken as a single endpoint, post hoc)")
    _add(out, "g4_holm_10_20mm_max_p_adj", pval(max(adj.values())), max(adj.values()),
         f"{src}.location_sign_flip_p; derived: largest Holm-adjusted p over the nine anatomies")
    n_pass = sum(v < 0.05 for v in adj.values())
    _add(out, "g4_holm_10_20mm_n_pass", count(n_pass), n_pass, f"{src}.location_sign_flip_p; derived: anatomies with Holm-adjusted p < 0.05")
    for group, labs in (("all", LABELS), ("smaller_heads", LABELS[1:]), ("templates_controls", TEMPLATES_CONTROLS),
                        ("childrenabc", CHILDREN)) + tuple((g, labs_) for g, labs_, _ in SPLIT):
        r = [P[lab]["s50_ratio_squid_over_opm"]["value"] for lab in labs]
        p = [P[lab]["location_sign_flip_p"] for lab in labs]
        nm = f"g4_{group}_dense_vs_combined_practical_10_20mm"
        _spread(out, f"{nm}_ratio", r, ratio, f"{src}.s50_ratio_squid_over_opm.value over {len(labs)} anatomies")
        _spread(out, f"{nm}_reduction_pct", [1 - 1 / x for x in r], pct,
                f"{src}.s50_ratio_squid_over_opm.value; derived: 1 - 1/ratio over {len(labs)} anatomies")
        _spread(out, f"{nm}_p", p, pval, f"{src}.location_sign_flip_p over {len(labs)} anatomies")
        lo = [P[lab]["locations_favouring_opm"] for lab in labs]
        _spread(out, f"{nm}_locs_opm", lo, count, f"{src}.locations_favouring_opm over {len(labs)} anatomies")
    n_bonf = sum(P[lab]["location_sign_flip_p"] < 0.05 / n_fam for lab in LABELS)
    _add(out, "g4_dense_vs_combined_practical_10_20mm_n_below_bonferroni_24", count(n_bonf), n_bonf,
         f"{src}.location_sign_flip_p; derived: anatomies with p < 0.05/{n_fam}")
    n_bonf = sum(P[lab]["location_sign_flip_p"] < 0.05 / (2 * n_fam) for lab in LABELS)
    _add(out, "g4_dense_vs_combined_practical_10_20mm_n_below_bonferroni_48", count(n_bonf), n_bonf,
         f"{src}.location_sign_flip_p; derived: anatomies with p < 0.05/{2 * n_fam}")
    n_sig = sum(P[lab]["location_sign_flip_p"] < 0.05 for lab in LABELS)
    _add(out, "g4_dense_vs_combined_practical_10_20mm_n_sig", count(n_sig), n_sig, f"{src}.location_sign_flip_p; derived: anatomies with p < 0.05")

    # deeper bands, oracle and the matched-site array, summarised over anatomies
    def cell(lab, a, mode, b):
        return C[f"{lab}/paired/{a}_vs_squid/combined/{mode}/depth{b}"]

    deep = [(lab, b) for lab in LABELS for b in (1, 2, 3)]
    sig_deep = [(lab, b) for lab, b in deep if cell(lab, "opm_dense/opm", "practical@1", b)["location_sign_flip_p"] < 0.05]
    srcd = f"{CMPF} :: comparison['<anatomy>/paired/opm_dense/opm_vs_squid/combined/practical@1/depth1..3'].location_sign_flip_p"
    _add(out, "g4_deeper_dense_vs_combined_practical_n_sig", count(len(sig_deep)), len(sig_deep),
         f"{srcd}; derived: cells (anatomy x band 20-70 mm) with p < 0.05: " + ", ".join(f"{lab} {BANDS[b]}" for lab, b in sig_deep))
    sd_small = [x for x in sig_deep if x[0] != "adult"]
    _add(out, "g4_deeper_dense_vs_combined_practical_n_sig_smaller_heads", count(len(sd_small)), len(sd_small),
         f"{srcd}; derived: as above without the adult")
    _add(out, "g4_deeper_dense_vs_combined_practical_n_cells", count(len(deep)), len(deep), f"{srcd}; derived: 9 anatomies x 3 bands")
    min_deep = min(cell(lab, "opm_dense/opm", "practical@1", b)["location_sign_flip_p"] for lab, b in deep if lab != "adult")
    _add(out, "g4_deeper_dense_vs_combined_practical_min_p_smaller_heads", pval(min_deep), min_deep, f"{srcd}; derived: smallest p, eight smaller heads")
    # direction of the deeper cells (locations favouring each system), whatever their p
    lean = {}
    for lab, b in deep:
        e = cell(lab, "opm_dense/opm", "practical@1", b)
        lean[(lab, b)] = (e["locations_favouring_opm"], e["locations_favouring_squid"])
    srcl = srcd.replace(".location_sign_flip_p", ".locations_favouring_opm, .locations_favouring_squid")
    for nm, test in (("favour_opm", lambda o, q: o > q), ("tie", lambda o, q: o == q), ("favour_squid", lambda o, q: o < q)):
        cells_ = [f"{lab} {BANDS[b]}" for (lab, b), (o, q) in lean.items() if test(o, q)]
        _add(out, f"g4_deeper_dense_vs_combined_practical_n_{nm}", count(len(cells_)), len(cells_),
             f"{srcl}; derived: cells (anatomy x band 20-70 mm) in which more locations favour "
             + {"favour_opm": "the dense OPM", "tie": "neither system (equal counts)", "favour_squid": "Neuromag"}[nm]
             + (": " + ", ".join(cells_) if nm != "favour_opm" else ""))
    n = sum(lean[x][0] > lean[x][1] for x in sig_deep)
    _add(out, "g4_deeper_dense_vs_combined_practical_n_sig_favour_opm", count(n), n,
         f"{srcl}; derived: of the cells with p < 0.05, those in which more locations favour the dense OPM")
    r2030 = [cell(lab, "opm_dense/opm", "practical@1", 1)["s50_ratio_squid_over_opm"]["value"] for lab in LABELS]
    n = sum(v is not None and v > 1 for v in r2030)
    _add(out, "g4_dense_vs_combined_practical_20_30mm_n_ratio_above1", count(n), n,
         f"{CMPF} :: comparison['<anatomy>/paired/opm_dense/opm_vs_squid/combined/practical@1/depth1'].s50_ratio_squid_over_opm.value; "
         "derived: anatomies with a point estimate above 1")
    # Holm over the nine anatomies, band by band, for the deeper bands (the rule of the 10-20 mm endpoint applied per band)
    passed = []
    for b in (1, 2, 3):
        adj_b = _holm({lab: cell(lab, "opm_dense/opm", "practical@1", b)["location_sign_flip_p"] for lab in LABELS})
        passed += [f"{lab} {BANDS[b]}" for lab, v in adj_b.items() if v < 0.05]
        if b == 3:
            _add(out, "g4_holm_45_70mm_adult_p_adj", pval(adj_b["adult"]), adj_b["adult"],
                 f"{srcd.replace('depth1..3', 'depth3')}; derived: Holm adjustment over the nine anatomies (the adult's value)")
    _add(out, "g4_holm_deeper_n_pass", count(len(passed)), len(passed),
         f"{srcd}; derived: cells with Holm-adjusted p < 0.05 when each deeper band is Holm-adjusted over the nine anatomies: "
         + (", ".join(passed) or "none"))
    m_ci = [lab for lab in LABELS
            if (cell(lab, "opm_matched/opm", "practical@1", 0)["s50_ratio_squid_over_opm"]["ci95"] or [None, None])[0] is not None
            and cell(lab, "opm_matched/opm", "practical@1", 0)["s50_ratio_squid_over_opm"]["ci95"][0] > 1]
    _add(out, "g4_matched_vs_combined_practical_10_20mm_n_ci_above1", count(len(m_ci)), len(m_ci),
         f"{CMPF} :: comparison['<anatomy>/paired/opm_matched/opm_vs_squid/combined/practical@1/depth0'].s50_ratio_squid_over_opm.ci95; "
         "derived: anatomies whose interval lies above 1: " + ", ".join(m_ci))
    for group, labs in (("templates_controls", TEMPLATES_CONTROLS), ("childrenabc", CHILDREN)):
        rest = [(lab, b) for lab in labs for b in (1, 2, 3) if (lab, b) not in sig_deep]
        pr_ = [cell(lab, "opm_dense/opm", "practical@1", b)["location_sign_flip_p"] for lab, b in rest]
        rr_ = [cell(lab, "opm_dense/opm", "practical@1", b)["s50_ratio_squid_over_opm"]["value"] for lab, b in rest]
        rr_ = [x for x in rr_ if x is not None]
        _add(out, f"g4_{group}_deeper_nonsig_min_p", pval(min(pr_)), min(pr_),
             f"{srcd}; derived: smallest p among the 20-70 mm cells with p >= 0.05")
        _spread(out, f"g4_{group}_deeper_nonsig_ratio", rr_, ratio,
                f"{srcd.replace('.location_sign_flip_p', '.s50_ratio_squid_over_opm.value')}; derived: uncensored ratios of the "
                "20-70 mm cells with p >= 0.05")
    or0 = [cell(lab, "opm_dense/opm", "oracle", 0) for lab in LABELS]
    srco = f"{CMPF} :: comparison['<anatomy>/paired/opm_dense/opm_vs_squid/combined/oracle/depth0']"
    _spread(out, "g4_oracle_dense_vs_combined_10_20mm_locs_opm", [x["locations_favouring_opm"] for x in or0], count,
            f"{srco}.locations_favouring_opm over the nine anatomies")
    _add(out, "g4_oracle_dense_vs_combined_10_20mm_locs_squid_max", count(max(x["locations_favouring_squid"] for x in or0)),
         max(x["locations_favouring_squid"] for x in or0), f"{srco}.locations_favouring_squid, largest over the nine anatomies")
    for b in (1, 2, 3):
        n = [lab for lab in LABELS if cell(lab, "opm_dense/opm", "oracle", b)["location_sign_flip_p"] < 0.05]
        _add(out, f"g4_oracle_dense_vs_combined_{BANDS[b]}_n_sig", count(len(n)), len(n),
             f"{srco.replace('depth0', f'depth{b}')}.location_sign_flip_p; derived: anatomies with p < 0.05: " + ", ".join(n))
    m0 = [lab for lab in LABELS if cell(lab, "opm_matched/opm", "practical@1", 0)["location_sign_flip_p"] < 0.05]
    _add(out, "g4_matched_vs_combined_practical_10_20mm_n_sig", count(len(m0)), len(m0),
         f"{CMPF} :: comparison['<anatomy>/paired/opm_matched/opm_vs_squid/combined/practical@1/depth0'].location_sign_flip_p; "
         "derived: anatomies with p < 0.05: " + ", ".join(m0))
    for lab in LABELS:
        pm = [cell(lab, "opm_matched/opm", "practical@1", b)["location_sign_flip_p"] for b in range(len(BANDS))]
        _add(out, f"g4_{_tag(lab)}_matched_vs_combined_practical_p_min", pval(min(pm)), min(pm),
             f"{CMPF} :: comparison['{lab}/paired/opm_matched/opm_vs_squid/combined/practical@1/depth0..3'].location_sign_flip_p; "
             "derived: smallest over the four bands")
    m3 = [cell(lab, "opm_matched/opm", "practical@1", 3) for lab in LABELS]
    n_def = sum(x["locations_favouring_squid"] > x["locations_favouring_opm"] and x["location_sign_flip_p"] < 0.05 for x in m3)
    _add(out, "g4_matched_vs_combined_practical_45_70mm_n_deficit", count(n_def), n_def,
         f"{CMPF} :: comparison['<anatomy>/paired/opm_matched/opm_vs_squid/combined/practical@1/depth3']; derived: anatomies "
         "where Neuromag is ahead with p < 0.05")

    # sampled location depths over anatomies
    for b, bt in enumerate(BANDS):
        meds = [_median([L["depth_mm"] for L in S[lab]["locations"] if L["stratum"][0] == b]) for lab in LABELS]
        _spread(out, f"g4_median_depth_{bt}", meds, mm,
                ", ".join(det_file.values()) + f" :: derived: median locations[].depth_mm with stratum[0] = {b}, over the nine anatomies")

    # matched held-out rates (in-sample check); its 'frozen' block must be the committed summaries
    for lab in LABELS:
        t = _tag(lab)
        for a in OPM:
            for b in range(len(BANDS)):
                fr = MR[lab]["frozen"][f"{a}_vs_squid/combined/depth{b}"]
                sm = S[lab]["paired"][f"{a}_vs_squid/combined/practical@1"][f"depth{b}"]
                if any(fr[k] != sm[k] for k in fr):
                    raise ValueError(f"{MRF}: frozen entry of {lab} {a} depth{b} differs from {det_file[lab]}")
        blk = MR[lab]["matched"]
        for a, at in OPM.items():
            for b, bt in enumerate(BANDS):
                r = blk[f"{a}_vs_squid/combined/depth{b}"]
                key = f"{MRF} :: anatomies['{lab}'].matched['{a}_vs_squid/combined/depth{b}']"
                base = f"g4_{t}_matchedrate_{at}_vs_combined_{bt}"
                _add(out, f"{base}_locs_opm", count(r["locations_favouring_opm"]), r["locations_favouring_opm"], f"{key}.locations_favouring_opm")
                _add(out, f"{base}_locs_squid", count(r["locations_favouring_squid"]), r["locations_favouring_squid"],
                     f"{key}.locations_favouring_squid")
                _add(out, f"{base}_p", pval(r["location_sign_flip_p"]), r["location_sign_flip_p"], f"{key}.location_sign_flip_p")
                rv = _ratio_fact(out, f"{base}_ratio", r["s50_ratio_squid_over_opm"], f"{key}.s50_ratio_squid_over_opm")
                if rv is not None:
                    _add(out, f"{base}_reduction_pct", pct(1 - 1 / rv), 1 - 1 / rv,
                         f"{key}.s50_ratio_squid_over_opm.value; derived: 1 - 1/ratio")
        for d, dt in (("squid/combined", "combined"), ("opm_matched/opm", "matched"), ("opm_dense/opm", "dense")):
            for b, bt in enumerate(BANDS):
                v = blk[f"s50/{d}/depth{b}"]
                key = f"{MRF} :: anatomies['{lab}'].matched['s50/{d}/depth{b}']"
                name = f"g4_{t}_matchedrate_s50_{dt}_{bt}"
                _add(out, name, "not reached" if v["value"] is None else integer(v["value"]),
                     "not reached" if v["value"] is None else v["value"], f"{key}.value")
                _add(out, f"{name}_ci", iv(v["ci95"][0], v["ci95"][1], integer), v["ci95"], f"{key}.ci95")
    mr0 = [MR[lab]["matched"]["opm_dense/opm_vs_squid/combined/depth0"] for lab in LABELS]
    srcm = f"{MRF} :: anatomies['<anatomy>'].matched['opm_dense/opm_vs_squid/combined/depth0']"
    _spread(out, "g4_matchedrate_dense_vs_combined_10_20mm_ratio", [x["s50_ratio_squid_over_opm"]["value"] for x in mr0], ratio,
            f"{srcm}.s50_ratio_squid_over_opm.value over the nine anatomies")
    _spread(out, "g4_matchedrate_dense_vs_combined_10_20mm_p", [x["location_sign_flip_p"] for x in mr0], pval,
            f"{srcm}.location_sign_flip_p over the nine anatomies")
    n_same = sum((x["location_sign_flip_p"] < 0.05) == (P[lab]["location_sign_flip_p"] < 0.05) for lab, x in zip(LABELS, mr0))
    _add(out, "g4_matchedrate_dense_vs_combined_10_20mm_n_same_conclusion", count(n_same), n_same,
         f"{srcm}.location_sign_flip_p vs {CMPF} frozen-threshold p; derived: anatomies with the same p < 0.05 verdict")

    # consistency with G2 at the G4 locations
    for b, bt in enumerate(BANDS):
        blk = CONS[f"depth{b}"]
        for k, nm in (("g2_median_ratio_dense_over_combined", "g4_g2_ratio_dense_vs_combined"),
                      ("oracle_s50_ratio_combined_over_dense", "g4_g2_oracle_s50_ratio")):
            v = blk[k]
            _add(out, f"{nm}_{bt}", ratio(v), v, f"{CONSF} :: bands.depth{b}.{k}")
            if 0.95 <= v <= 1.05:
                _add(out, f"{nm}_{bt}_3dp", ratio3(v), v, f"{CONSF} :: bands.depth{b}.{k}")

    # earlier runs (superseded result files; numbers recorded in the methods text only)
    M = "docs/methods.md :: section 9 'Scoring and run-to-run variability' (earlier run; its result files were superseded)"
    for name, txt, raw, what in (
            ("g4_earlier_v3_matched_vs_combined_45_70mm", "1/9", [1, 9], "v3 matched OPM vs Neuromag 45-70 mm, locations OPM/Neuromag"),
            ("g4_earlier_v3_matched_vs_combined_45_70mm_p", "0.016", 0.016, "v3 matched OPM 45-70 mm sign-flip p"),
            ("g4_earlier_v2_matched_vs_combined_45_70mm", "0/7", [0, 7], "v2 matched OPM vs Neuromag 45-70 mm, locations OPM/Neuromag"),
            ("g4_earlier_v3_dense_vs_combined_10_20mm", "8/2", [8, 2], "v3 dense OPM vs Neuromag 10-20 mm, locations OPM/Neuromag "
             "('v3 ... gave dense 8/2, 5/1, 4/0 and 0/2')"),
            ("g4_earlier_v3_dense_vs_combined_10_20mm_p", "0.037", 0.037, "v3 dense OPM 10-20 mm sign-flip p ('p = 0.037, ...')"),
            ("g4_earlier_v3_dense_vs_combined_10_20mm_ratio", "1.29", 1.29, "v3 dense OPM 10-20 mm S50 ratio ('S50 ratio 1.29 "
             "[1.03-1.51] at 10-20 mm')")):
        _add(out, name, txt, raw, f"{M}: {what}")


def _median(v: list) -> float:
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


# ----------------------------------------------------------------------------------------------
# localization
def _locfile(lab: str) -> str:
    return "results/g4/g4_localization_summary.json" if lab == "adult" else f"results/g4/g4_localization_{lab}_summary.json"


def _locevents(lab: str) -> str:
    return "results/g4/g4_localization_events.csv" if lab == "adult" else f"results/g4/g4_localization_{lab}_events.csv"


CELL = (("detected", "detected", share), ("ecd_det_mm", "ecd_error_mm_median_detected", mm),
        ("ecd_all_mm", "ecd_error_mm_median_all", mm), ("ecd_sensor_det_mm", "ecd_error_sensor_frame_mm_median_detected", mm),
        ("dspm_det_mm", "dspm_error_mm_median_detected", mm), ("dspm_all_mm", "dspm_error_mm_median_all", mm),
        ("dspm_mne_det_mm", "dspm_mne_error_mm_median_detected", mm), ("dspm_mne_all_mm", "dspm_mne_error_mm_median_all", mm),
        ("joint_dspm", "joint_detect_and_dspm_within_10mm", share), ("joint_dspm_mne", "joint_detect_and_dspm_mne_within_10mm", share),
        ("joint_ecd", "joint_detect_and_ecd_within_10mm", share), ("gof_det_pct", "ecd_gof_median_detected", lambda x: pct(x / 100)),
        ("chi2_det", "ecd_khi2_per_dof_median_detected", ratio), ("support", "support_recovery_median", share))
ERRORS = (("dspm_error_mm", "dspm"), ("dspm_mne_error_mm", "dspm_mne"), ("ecd_error_mm", "ecd"))
JOINT = (("joint_dspm_10mm", "joint_dspm"), ("joint_ecd_10mm", "joint_ecd"))


def _cond_tag(cond: str) -> str:
    fam, s = cond.split("/")
    return f"{fam}{s.replace('nAm', '')}"


def _localization(root: Path, out: dict) -> None:
    files = {lab: _locfile(lab) for lab in LABELS}
    L = {lab: _json(root / files[lab]) for lab in LABELS}
    ALL = ", ".join(files.values())
    A = files["adult"]
    cfg = L["adult"]["config"]
    for lab in LABELS:
        if L[lab]["config"] != cfg:
            raise ValueError(f"{lab}: localization configuration differs from the adult's")
    same = " (identical in all nine localization summaries)"
    for name, val, key in (
            ("loc_n_locations", cfg["n_locations"], "n_locations"), ("loc_strengths_nam", cfg["strengths_nAm"], "strengths_nAm"),
            ("loc_strength_low_nam", cfg["strengths_nAm"][0], "strengths_nAm[0]"),
            ("loc_strength_high_nam", cfg["strengths_nAm"][1], "strengths_nAm[1]"),
            ("loc_grid_spacing_mm", cfg["grid_spacing_mm"], "grid_spacing_mm"), ("loc_coreg_shift_mm", cfg["coreg_shift_mm"], "coreg_shift_mm"),
            ("loc_coreg_angle_deg", cfg["coreg_angle_deg"], "coreg_angle_deg"), ("loc_coreg_draws", cfg["coreg_draws"], "coreg_draws"),
            ("loc_baseline_minutes", cfg["baseline_min"], "baseline_min (noise covariance)"),
            ("loc_calibration_minutes", cfg["calibration_min"], "calibration_min (detector thresholds, 1 false event/min)"),
            ("loc_holdout_minutes", cfg["holdout_min"], "holdout_min (independent held-out null)"),
            ("loc_snr", cfg["snr"], "snr (lambda^2 = 1/snr^2)"), ("loc_depth_weighting", cfg["depth"], "depth"),
            ("loc_ecd_min_dist_mm", cfg["ecd_min_dist_mm"], "ecd_min_dist_mm (fit_dipole min_dist to the inner skull)"),
            ("loc_seed", cfg["seed"], "seed")):
        _add(out, name, str(val) if name.endswith("_seed") else const(val), val, f"{A} :: config.{key}{same}")
    layers = [int(x) for x in re.findall(r"(\d+)-layer", cfg["inverse_bem"])]
    _add(out, "loc_inverse_bem_layers", count(layers[0]), layers[0], f"{A} :: config.inverse_bem ('{cfg['inverse_bem']}'){same}")
    _add(out, "loc_truth_bem_layers", count(layers[1]), layers[1], f"{A} :: config.inverse_bem ('{cfg['inverse_bem']}'){same}")
    per_draw = cfg["n_locations"] // cfg["coreg_draws"]
    _add(out, "loc_locations_per_coreg_draw", count(per_draw), per_draw,
         f"{A} :: derived: config.n_locations / config.coreg_draws (docs/methods.md section 9: location i uses draw i mod 8)")
    _add(out, "loc_coreg_draws_per_location", "1", 1, "docs/methods.md :: section 9 'location i uses draw i mod 8 in every condition'")
    n_ev = L["adult"]["n_events"]
    _add(out, "loc_n_events_per_anatomy", count(n_ev), n_ev, f"{A} :: n_events{same}")
    per = cfg["n_locations"] // 12
    _add(out, "loc_locations_per_stratum", count(per), per, f"{A} :: derived: config.n_locations / 12 depth x orientation strata")
    n_cond = {r["n"] for lab in LABELS for r in L[lab]["results"].values()}
    if len(n_cond) != 1:
        raise ValueError("localization conditions differ in their event counts")
    n_cond = n_cond.pop()
    _add(out, "loc_events_per_condition", count(n_cond), n_cond,
         f"{A} :: results[*].n (one event per location and condition; the same in all nine localization summaries)")
    fam_n = len(CONDS) * (len(ERRORS) + len(JOINT))
    _add(out, "loc_family_per_opm_array", count(fam_n), fam_n, f"{A} :: derived: paired conditions (4) x endpoints (5) per OPM array")
    _add(out, "loc_family_per_anatomy", count(2 * fam_n), 2 * fam_n, f"{A} :: derived: two OPM arrays x {fam_n}")
    for k, d in (("20", fam_n), ("40", 2 * fam_n), ("320", 8 * 2 * fam_n), ("360", 9 * 2 * fam_n)):
        _add(out, f"loc_bonferroni_{k}_threshold", pval(0.05 / d), 0.05 / d, f"{A} :: derived: 0.05 / {d} comparisons")

    grids = {lab: L[lab]["inverse_grid"] for lab in LABELS}
    coreg = {lab: L[lab]["coreg_displacement_mm_median"] for lab in LABELS}
    held_all, held_primary = [], []
    for lab in LABELS:
        s, f, t = L[lab], files[lab], _tag(lab)
        _add(out, f"loc_{t}_inverse_grid", count(grids[lab]), grids[lab], f"{f} :: inverse_grid (5-mm grid sources, true vertices excluded)")
        _add(out, f"loc_{t}_coreg_displacement_mm", mm(coreg[lab]), coreg[lab], f"{f} :: coreg_displacement_mm_median")
        for v, info in s["views"].items():
            _add(out, f"loc_{t}_n_channels_{v.replace('squid_', '').replace('opm_', '').replace('squid', 'combined')}",
                 count(info["channels"]), info["channels"], f"{f} :: views['{v}'].channels")
        for v, h in s["thresholds_heldout"]["false_events"].items():
            vt = v.replace("squid_", "").replace("opm_", "").replace("squid", "combined")
            _add(out, f"loc_{t}_heldout_rate_{vt}", share(h["rate_per_min"]), h["rate_per_min"],
                 f"{f} :: thresholds_heldout.false_events['{v}'].rate_per_min ({h['count']} events in {h['minutes']:g} min)")
            _add(out, f"loc_{t}_heldout_rate_{vt}_ci", iv(h["ci95"][0], h["ci95"][1], share), h["ci95"],
                 f"{f} :: thresholds_heldout.false_events['{v}'].ci95 (exact Poisson)")
            held_all.append(h["rate_per_min"])
            if v in ARRAYS:
                held_primary.append(h["rate_per_min"])
        # per array and condition
        for arr, at in ARRAYS.items():
            for cond in CONDS:
                r = s["results"][f"{arr}/{cond}"]
                for nm, k, fn in CELL:
                    if r.get(k) is None:
                        continue
                    _add(out, f"loc_{t}_{at}_{_cond_tag(cond)}_{nm}", fn(r[k]), r[k], f"{f} :: results['{arr}/{cond}'].{k}")
                    if lab == "adult" and nm == "joint_dspm_mne":  # shares of 24 events (0.125, 0.333, ...) also at 3 decimals
                        _add(out, f"loc_{t}_{at}_{_cond_tag(cond)}_{nm}_3dp", ratio3(r[k]), r[k], f"{f} :: results['{arr}/{cond}'].{k}")
        # paired, OPM minus Neuromag combined, with untied counts from the event table
        ev = _csv(root / _locevents(lab))
        by = {(r["array"], r["event"]): r for r in ev}
        n_sig = 0
        for arr in ("opm_dense", "opm_matched"):
            for cond in CONDS:
                fam, st = cond.split("/")
                v = s["paired"][f"{arr}_vs_squid/{cond}"]
                key = f"{f} :: paired['{arr}_vs_squid/{cond}']"
                base = f"loc_{t}_{ARRAYS[arr]}_vs_combined_{_cond_tag(cond)}"
                evs = sorted({r["event"] for r in ev if r["family"] == fam and float(r["strength_nAm"]) == float(st[:-3])}, key=int)
                for k, et in ERRORS:
                    r = v[k]
                    _add(out, f"{base}_{et}_diff_mm", dmm(r["median_difference"]), r["median_difference"],
                         f"{key}.{k}.median_difference (OPM minus Neuromag; negative = OPM closer)")
                    _add(out, f"{base}_{et}_diff_mm_ci", iv(r["ci95"][0], r["ci95"][1], dmm), r["ci95"], f"{key}.{k}.ci95 (bootstrap over events)")
                    _add(out, f"{base}_{et}_p", pval(r["wilcoxon_p"]), r["wilcoxon_p"], f"{key}.{k}.wilcoxon_p (uncorrected)")
                    _add(out, f"{base}_{et}_share_opm_closer", share(r["share_opm_smaller"]), r["share_opm_smaller"],
                         f"{key}.{k}.share_opm_smaller")
                    if et != "ecd":  # dSPM peaks on the shared grid can tie; dipole errors do not
                        diff = [float(by[arr, e][k]) - float(by["squid", e][k]) for e in evs]
                        for nm, n in (("n_opm_closer", sum(x < 0 for x in diff)), ("n_squid_closer", sum(x > 0 for x in diff)),
                                      ("n_tied", sum(x == 0 for x in diff))):
                            _add(out, f"{base}_{et}_{nm}", count(n), n,
                                 f"{_locevents(lab)} :: derived: {k} of {arr} minus squid per event, {fam} {st} ({len(evs)} events)")
                    n_sig += r["wilcoxon_p"] < 0.05
                for k, et in JOINT:
                    r = v[k]
                    _add(out, f"{base}_{et}_opm", share(r["opm"]), r["opm"], f"{key}.{k}.opm (detected and within 10 mm)")
                    _add(out, f"{base}_{et}_squid", share(r["squid"]), r["squid"], f"{key}.{k}.squid")
                    _add(out, f"{base}_{et}_only_opm", count(r["only_opm"]), r["only_opm"], f"{key}.{k}.only_opm (locations)")
                    _add(out, f"{base}_{et}_only_squid", count(r["only_squid"]), r["only_squid"], f"{key}.{k}.only_squid (locations)")
                    _add(out, f"{base}_{et}_p", pval(r["mcnemar_exact_p"]), r["mcnemar_exact_p"], f"{key}.{k}.mcnemar_exact_p (uncorrected)")
                    n_sig += r["mcnemar_exact_p"] < 0.05
        _add(out, f"loc_{t}_n_sig_p05", count(n_sig), n_sig, f"{f} :: paired (all 40 entries: 2 OPM arrays x 4 conditions x 5 endpoints); "
             "derived: count with p < 0.05 (Wilcoxon for errors, exact McNemar for joint success)")

    # the adult's secondary comparators (Neuromag magnetometers or gradiometers alone)
    s = L["adult"]
    for ref, rt in (("squid_mag", "mag"), ("squid_grad", "grad")):
        n_sig, n_all = 0, 0
        for arr in ("opm_dense", "opm_matched"):
            for cond in CONDS:
                v = s["paired_secondary"][f"{arr}_vs_{ref}/{cond}"]
                for k, et in ERRORS:
                    r = v[k]
                    base = f"loc_adult_{ARRAYS[arr]}_vs_{rt}_{_cond_tag(cond)}_{et}"
                    key = f"{A} :: paired_secondary['{arr}_vs_{ref}/{cond}'].{k}"
                    _add(out, f"{base}_diff_mm", dmm(r["median_difference"]), r["median_difference"], f"{key}.median_difference")
                    _add(out, f"{base}_p", pval(r["wilcoxon_p"]), r["wilcoxon_p"], f"{key}.wilcoxon_p (uncorrected)")
                    n_sig += r["wilcoxon_p"] < 0.05
                    n_all += 1
                for k, et in JOINT:
                    n_sig += v[k]["mcnemar_exact_p"] < 0.05
                    n_all += 1
        _add(out, f"loc_adult_vs_{rt}_n_sig_p05", count(n_sig), n_sig,
             f"{A} :: paired_secondary['*_vs_{ref}/*'] ({n_all} entries); derived: count with p < 0.05")

    # summaries over arrays and anatomies
    _spread(out, "loc_inverse_grid", list(grids.values()), count, f"{ALL} :: inverse_grid over the nine anatomies")
    _spread(out, "loc_coreg_displacement_mm", list(coreg.values()), mm, f"{ALL} :: coreg_displacement_mm_median over the nine anatomies")
    _spread(out, "loc_heldout_rate_primary", held_primary, share,
            f"{ALL} :: thresholds_heldout.false_events[squid|opm_matched|opm_dense].rate_per_min over the nine anatomies")
    _spread(out, "loc_heldout_rate_all_views", held_all, share,
            f"{ALL} :: thresholds_heldout.false_events[*].rate_per_min, five views (Neuromag combined, mag, grad; OPM matched, dense)")
    _spread(out, "loc_adult_heldout_rate", [h["rate_per_min"] for h in L["adult"]["thresholds_heldout"]["false_events"].values()], share,
            f"{A} :: thresholds_heldout.false_events[*].rate_per_min over the five views")
    res = lambda labs, arrs, cond, k: [L[lab]["results"][f"{a}/{cond}"][k] for lab in labs for a in arrs]  # noqa: E731
    _spread(out, "loc_ecd_det_focal320_mm", res(LABELS, ARRAYS, "focal/320nAm", "ecd_error_mm_median_detected"), mm,
            f"{ALL} :: results['<array>/focal/320nAm'].ecd_error_mm_median_detected, nine anatomies x three arrays")
    _spread(out, "loc_childrenabc_ecd_det_focal320_mm", res(CHILDREN, ARRAYS, "focal/320nAm", "ecd_error_mm_median_detected"), mm,
            f"{ALL} :: results['<array>/focal/320nAm'].ecd_error_mm_median_detected, children A-C x three arrays")
    for k, labs, what in SPLIT:
        _spread(out, f"loc_{k}_ecd_det_focal320_mm", res(labs, ARRAYS, "focal/320nAm", "ecd_error_mm_median_detected"), mm,
                f"{ALL} :: results['<array>/focal/320nAm'].ecd_error_mm_median_detected, {what} x three arrays")
        _spread(out, f"loc_{k}_coreg_displacement_mm", [coreg[lab] for lab in labs], mm,
                f"{ALL} :: coreg_displacement_mm_median over {what}")
    _spread(out, "loc_adult_patch80_detected", res(["adult"], ARRAYS, "patch/80nAm", "detected"), share,
            f"{A} :: results['<array>/patch/80nAm'].detected over the three arrays")
    _spread(out, "loc_adult_chi2_det", [L["adult"]["results"][f"{a}/{c}"]["ecd_khi2_per_dof_median_detected"] for a in ARRAYS for c in CONDS],
            ratio, f"{A} :: results['<array>/<condition>'].ecd_khi2_per_dof_median_detected over three arrays and four conditions")
    adiff = [abs(L["adult"]["paired"][f"{a}_vs_squid/{c}"]["ecd_error_mm"]["median_difference"]) for a in ("opm_dense", "opm_matched") for c in CONDS]
    ap = [L["adult"]["paired"][f"{a}_vs_squid/{c}"]["ecd_error_mm"]["wilcoxon_p"] for a in ("opm_dense", "opm_matched") for c in CONDS]
    _add(out, "loc_adult_ecd_abs_diff_max_mm", mm(max(adiff)), max(adiff), f"{A} :: paired[*].ecd_error_mm.median_difference; derived: largest |difference|")
    _add(out, "loc_adult_ecd_p_min", pval(min(ap)), min(ap), f"{A} :: paired[*].ecd_error_mm.wilcoxon_p; derived: smallest")

    # the family of paired comparisons over the nine anatomies
    rows = []
    for lab in LABELS:
        for key, v in L[lab]["paired"].items():
            arr = key.split("_vs_")[0]
            fam, st = key.split("/")[1:]
            for k, et in ERRORS:
                m = v[k]["median_difference"]
                rows.append(dict(lab=lab, arr=arr, fam=fam, st=st, ep=et, p=v[k]["wilcoxon_p"], dir="opm" if m < 0 else ("squid" if m > 0 else "zero"),
                                 eff=m))
            for k, et in JOINT:
                d = v[k]["only_opm"] - v[k]["only_squid"]
                rows.append(dict(lab=lab, arr=arr, fam=fam, st=st, ep=et, p=v[k]["mcnemar_exact_p"],
                                 dir="opm" if d > 0 else ("squid" if d < 0 else "zero"), eff=d))
    src = f"{ALL} :: paired (all entries); derived"
    sig_ = [r for r in rows if r["p"] < 0.05]
    _add(out, "loc_family_n", count(len(rows)), len(rows), f"{src}: comparisons (9 anatomies x 2 OPM arrays x 4 conditions x 5 endpoints)")
    _add(out, "loc_family_per_endpoint", count(len(rows) // 5), len(rows) // 5, f"{src}: comparisons per endpoint")
    _add(out, "loc_family_n_sig", count(len(sig_)), len(sig_), f"{src}: p < 0.05 (uncorrected)")
    _add(out, "loc_family_expected_by_chance", count(round(0.05 * len(rows))), 0.05 * len(rows),
         f"{src}: 0.05 x comparisons (if independent, which they are not)")
    for k, et in ERRORS + JOINT:
        n = sum(r["ep"] == et for r in sig_)
        _add(out, f"loc_family_n_sig_{et}", count(n), n, f"{src}: p < 0.05 among the {et} comparisons")
    for dname in ("opm", "zero", "squid"):
        n = sum(r["dir"] == dname for r in sig_)
        _add(out, f"loc_family_n_sig_favour_{dname}", count(n), n,
             f"{src}: p < 0.05 and {'favouring an OPM array' if dname == 'opm' else ('zero median difference' if dname == 'zero' else 'favouring Neuromag')} "
             "(errors: sign of the median OPM-minus-Neuromag difference; joint: more locations successful only with the OPM)")
    n = sum(r["dir"] in ("opm", "zero") for r in sig_)
    _add(out, "loc_family_n_sig_favour_opm_or_zero", count(n), n, f"{src}: p < 0.05 and not favouring Neuromag (all endpoints)")
    fav_err = sum(r["dir"] in ("opm", "zero") and r["ep"] in ("dspm", "dspm_mne") for r in sig_)
    _add(out, "loc_family_n_sig_favour_opm_or_zero_dspm_errors", count(fav_err), fav_err, f"{src}: p < 0.05, study or MNE dSPM errors, not favouring Neuromag")
    for st in ("80nAm", "320nAm"):
        n = sum(r["st"] == st for r in sig_)
        _add(out, f"loc_family_n_sig_{st.lower()}", count(n), n, f"{src}: p < 0.05 at {st[:-3]} nAm")
    abc = [r for r in rows if r["lab"] in CHILDREN]
    _add(out, "loc_childrenabc_n", count(len(abc)), len(abc), f"{src}: comparisons in children A-C")
    _add(out, "loc_childrenabc_n_sig", count(sum(r["p"] < 0.05 for r in abc)), sum(r["p"] < 0.05 for r in abc), f"{src}: p < 0.05, children A-C")
    _add(out, "loc_childrenabc_n_sig_favour_opm", count(sum(r["p"] < 0.05 and r["dir"] == "opm" for r in abc)),
         sum(r["p"] < 0.05 and r["dir"] == "opm" for r in abc), f"{src}: p < 0.05 favouring an OPM array, children A-C")
    pmin = min(r["p"] for r in abc)
    _add(out, "loc_childrenabc_min_p", pval(pmin), pmin, f"{src}: smallest p in children A-C")
    _add(out, "loc_childrenabc_min_p_3sf", sig(pmin, 3), pmin,
         f"{src}: smallest p in children A-C, 3 significant digits (2 would print it equal to the 0.0025 threshold it misses)")
    for k, labs, what in SPLIT:
        sub = [r for r in rows if r["lab"] in labs]
        sg = [r for r in sub if r["p"] < 0.05]
        _add(out, f"loc_{k}_n", count(len(sub)), len(sub), f"{src}: comparisons in {what}")
        _add(out, f"loc_{k}_n_sig", count(len(sg)), len(sg), f"{src}: p < 0.05, {what}")
        _add(out, f"loc_{k}_expected_by_chance", count(round(0.05 * len(sub))), 0.05 * len(sub),
             f"{src}: 0.05 x comparisons in {what} (if independent, which they are not)")
        _add(out, f"loc_{k}_n_sig_favour_opm", count(sum(r["dir"] == "opm" for r in sg)), sum(r["dir"] == "opm" for r in sg),
             f"{src}: p < 0.05 favouring an OPM array, {what}")
        n = sum(r["p"] < 0.05 / fam_n for r in sub)
        _add(out, f"loc_survivors_bonferroni_20_n_{k}", count(n), n,
             f"{src}: p < 0.05/{fam_n} (within anatomy and OPM array), survivors in {what}")
    small = [r for r in rows if r["lab"] != "adult"]
    _add(out, "loc_smaller_heads_min_p", pval(min(r["p"] for r in small)), min(r["p"] for r in small), f"{src}: smallest p in the eight smaller heads")
    _add(out, "loc_family_min_p", pval(min(r["p"] for r in rows)), min(r["p"] for r in rows), f"{src}: smallest p over the nine anatomies")
    for k, d in (("20", fam_n), ("40", 2 * fam_n)):
        sv = [r for r in rows if r["p"] < 0.05 / d]
        _add(out, f"loc_survivors_bonferroni_{k}", count(len(sv)), len(sv),
             f"{src}: p < 0.05/{d} (within anatomy{' and OPM array' if k == '20' else ', both arrays'}): "
             + "; ".join(f"{r['lab']} {r['arr']} {r['fam']} {r['st']} {r['ep']} p={r['p']:.2g}" for r in sv))
        if k == "20":
            for st in ("80nAm", "320nAm"):
                n = sum(r["st"] == st for r in sv)
                _add(out, f"loc_survivors_bonferroni_20_n_{st.lower()}", count(n), n, f"{src}: survivors at {st[:-3]} nAm")
            for lab in LABELS:
                n = sum(r["lab"] == lab for r in sv)
                if n:
                    _add(out, f"loc_survivors_bonferroni_20_n_{_tag(lab)}", count(n), n, f"{src}: survivors in {lab}")
                    for st in ("80nAm", "320nAm"):
                        n = sum(r["lab"] == lab and r["st"] == st for r in sv)
                        _add(out, f"loc_survivors_bonferroni_20_n_{_tag(lab)}_{st.lower()}", count(n), n,
                             f"{src}: survivors in {lab} at {st[:-3]} nAm")
            n = sum(r["dir"] == "opm" for r in sv)
            _add(out, "loc_survivors_bonferroni_20_n_favour_opm", count(n), n, f"{src}: survivors favouring an OPM array")
    for k, labs, d in (("smaller_heads", LABELS[1:], 8 * 2 * fam_n), ("all", LABELS, 9 * 2 * fam_n)) + tuple(
            (g, labs_, len(labs_) * 2 * fam_n) for g, labs_, _ in SPLIT):
        n = sum(r["p"] < 0.05 / d for r in rows if r["lab"] in labs)
        _add(out, f"loc_survivors_bonferroni_{k}", count(n), n, f"{src}: p < 0.05/{d} over {len(labs)} anatomies")

    # failed fits under declared criteria
    FF = "results/g4/g4_fit_failures.json"
    ff = _json(root / FF)
    _add(out, "loc_gross_error_threshold_mm", const(ff["criteria"]["gross_error_mm"]), ff["criteria"]["gross_error_mm"], f"{FF} :: criteria.gross_error_mm")
    _add(out, "loc_joint_threshold_mm", "10", 10.0, "results/g4/g4_localization_summary.json :: results[*].joint_detect_and_dspm_within_10mm, "
         "joint_detect_and_ecd_within_10mm (the success radius in the stored keys; docs/methods.md section 9 'Detected and "
         "localized within 10 mm')")
    _add(out, "loc_unconstrained_volume_threshold_cm3", const(ff["criteria"]["unconstrained_volume_cm3"]),
         ff["criteria"]["unconstrained_volume_cm3"], f"{FF} :: criteria.unconstrained_volume_cm3")
    tot = ff["totals"]
    rates_ecd, rates_dspm = [], []
    for arr, sfx, at in (("squid", "", "combined"), ("opm_matched", "_matched", "matched"), ("opm_dense", "_dense", "dense")):
        n = tot["n" if not sfx else f"n{sfx}"]
        ndet = sum(v["n_detected"] for a in ff["anatomies"].values() for k2, v in a.items() if k2.startswith(arr + "/"))
        _add(out, f"loc_failures_{at}_n_events", count(n), n, f"{FF} :: totals.n{sfx} (nine anatomies)")
        _add(out, f"loc_failures_{at}_n_detected", count(ndet), ndet, f"{FF} :: derived: sum of anatomies[*]['{arr}/*'].n_detected")
        for k in ("ecd_gross", "ecd_gross_detected", "dspm_gross", "dspm_gross_detected"):
            _add(out, f"loc_failures_{at}_{k}", count(tot[k + sfx]), tot[k + sfx], f"{FF} :: totals.{k}{sfx}")
        uc = sum(v["ecd_unconstrained"] for a in ff["anatomies"].values() for k2, v in a.items() if k2.startswith(arr + "/"))
        ucd = sum(v["ecd_unconstrained_detected"] for a in ff["anatomies"].values() for k2, v in a.items() if k2.startswith(arr + "/"))
        _add(out, f"loc_failures_{at}_ecd_unconstrained", count(uc), uc, f"{FF} :: derived: sum of anatomies[*]['{arr}/*'].ecd_unconstrained")
        _add(out, f"loc_failures_{at}_ecd_unconstrained_detected", count(ucd), ucd,
             f"{FF} :: derived: sum of anatomies[*]['{arr}/*'].ecd_unconstrained_detected")
        re_, rd_ = tot["ecd_gross_detected" + sfx] / ndet, tot["dspm_gross_detected" + sfx] / ndet
        rates_ecd.append(re_)
        rates_dspm.append(rd_)
        _add(out, f"loc_failures_{at}_ecd_gross_rate_detected", pct(re_), re_, f"{FF} :: derived: totals.ecd_gross_detected{sfx} / detected events")
        _add(out, f"loc_failures_{at}_dspm_gross_rate_detected", pct(rd_), rd_, f"{FF} :: derived: totals.dspm_gross_detected{sfx} / detected events")
    _spread(out, "loc_failures_ecd_gross_rate_detected", rates_ecd, pct, f"{FF} :: derived: gross ECD errors among detected events, three arrays")
    _add(out, "loc_failures_gross_undetected_share", pct(tot["gross_undetected_share"]), tot["gross_undetected_share"],
         f"{FF} :: totals.gross_undetected_share")
    _add(out, "loc_failures_gross_80nam_share", pct(tot["gross_80nAm_share"]), tot["gross_80nAm_share"], f"{FF} :: totals.gross_80nAm_share")
    ev = [r for r in _csv(root / _locevents("adult")) if r["array"] in ARRAYS]
    g = [r for r in ev if float(r["ecd_error_mm"]) > ff["criteria"]["gross_error_mm"]]
    how = f"{_locevents('adult')} :: derived: three primary arrays, ecd_error_mm > {ff['criteria']['gross_error_mm']:g}"
    _add(out, "loc_adult_ecd_n", count(len(ev)), len(ev), f"{how.split(', ecd')[0]}: dipole fits")
    _add(out, "loc_adult_ecd_gross_n", count(len(g)), len(g), how)
    gu = [r for r in g if r["detected"] != "True"]
    _add(out, "loc_adult_ecd_gross_undetected_n", count(len(gu)), len(gu), f"{how}, undetected")
    gu80 = sum(float(r["strength_nAm"]) == 80.0 for r in gu)
    _add(out, "loc_adult_ecd_gross_undetected_80nam_n", count(gu80), gu80, f"{how}, undetected, 80 nAm")

    # earlier runs (superseded; numbers recorded in the methods text only)
    M = ("docs/methods.md :: section 9 'Earlier runs: v3 0.0 mm (p = 0.50 and 0.81), v2 -4.8 and -5.3 mm (p = 0.14 and 0.046)' "
         "(320-nAm patch dSPM, OPM minus Neuromag; dense then matched, the order of the sentence before it; earlier runs whose "
         "result files were superseded)")
    M11 = ("docs/methods.md :: section 11 'v3 had four survivors under 16 comparisons, two of them among these; v2 seven' "
           "(earlier runs whose result files were superseded)")
    for name, txt, raw, src_ in (
            ("loc_earlier_v3_patch320_dspm_diff_mm", "0.0", 0.0, M),
            ("loc_earlier_v3_dense_patch320_dspm_p", "0.50", 0.50, M),
            ("loc_earlier_v3_matched_patch320_dspm_p", "0.81", 0.81, M),
            ("loc_earlier_v2_dense_patch320_dspm_diff_mm", MINUS + "4.8", -4.8, M),
            ("loc_earlier_v2_matched_patch320_dspm_diff_mm", MINUS + "5.3", -5.3, M),
            ("loc_earlier_v2_dense_patch320_dspm_p", "0.14", 0.14, M),
            ("loc_earlier_v2_matched_patch320_dspm_p", "0.046", 0.046, M),
            ("loc_earlier_v3_survivors", "4", 4, M11),
            ("loc_earlier_v3_family_per_opm_array", "16", 16, M11),
            ("loc_earlier_v2_survivors", "7", 7, M11)):
        _add(out, name, txt, raw, src_)


# ----------------------------------------------------------------------------------------------
# motion
def _case(case: str) -> str:
    """'down 2 mm' -> 'down2mm', 'x+2 mm' -> 'xplus2mm', 'pitch -10 deg' -> 'pitchminus10deg', 'slip x +1 deg' -> 'slipxplus1deg'."""
    return case.replace("+", "plus").replace("-", "minus").replace(" ", "")


def _amount(case: str) -> tuple:
    """(magnitude, 'mm' or 'deg') of a displacement or slip label."""
    num, unit = case.split()[-2:]
    return float(num.lstrip("+-xyz")), unit


CORR = {"none": "none", "homogeneous": "homog", "homogeneous+gradient": "eightterm"}
CAL = {"tilt0deg_gain0pct": "cal0", "tilt1deg_gain1pct": "cal1", "tilt3deg_gain3pct": "cal3"}
LEVEL = {"loss_1dB": "loss1db", "loss_3dB": "loss3db", "D_0dB": "d0"}


def _motion(root: Path, out: dict) -> None:
    F = "results/g4/g4_motion_summary.json"
    m = _json(root / F)
    cfg = m["config"]
    g, cb, tc = cfg["geometry"], cfg["coupling"], cfg["timecourse"]
    for name, val, key in (
            ("mot_n_anatomies", len(m["anatomies"]), "anatomies (count: adult, 24- and 12-month templates)"),
            ("mot_translations_mm", g["translations_mm"], "geometry.translations_mm (down, +-x, +-y)"),
            ("mot_rotations_deg", g["rotations_deg"], "geometry.rotations_deg (pitch, roll, yaw, both signs)"),
            ("mot_slip_deg", g["slip_deg"], "geometry.slip_deg (about x, y, z, both signs)"),
            ("mot_wrong_polarity_floor_db", g["floor_db"], "geometry.floor_db"),
            ("mot_pivot_below_origin_mm", round(-cb["pivot_head_m"][2] * 1000, 6), "coupling.pivot_head_m[2] x -1000"),
            ("mot_n_draws", cb["n_draws"], "coupling.n_draws (common random numbers)"),
            ("mot_rotation_grid_min_deg", cb["rotation_rms_deg"][0], "coupling.rotation_rms_deg[0]"),
            ("mot_rotation_grid_max_deg", cb["rotation_rms_deg"][-1], "coupling.rotation_rms_deg[-1]"),
            ("mot_calibration_tilt_deg", [c[0] for c in cb["calibration"]], "coupling.calibration[*][0] (RMS sensitive-axis tilt)"),
            ("mot_loss_levels_db", cb["loss_db"], "coupling.loss_db"), ("mot_coupling_seed", cb["seed"], "coupling.seed"),
            ("mot_tc_duration_s", tc["duration_s"], "timecourse.duration_s"), ("mot_tc_fs_hz", tc["fs_hz"], "timecourse.fs_hz"),
            ("mot_tc_drift_deg", tc["drift_deg"], "timecourse.drift_deg (below 0.1 Hz, peak per axis)"),
            ("mot_tc_inband_rotation_deg", tc["inband_rotation_rms_deg"], "timecourse.inband_rotation_rms_deg (RMS per axis)"),
            ("mot_tc_b0_nt", tc["b0_nT"], "timecourse.b0_nT"), ("mot_tc_gradient_nt_per_m", tc["gradient_nT_per_m"], "timecourse.gradient_nT_per_m"),
            ("mot_tc_calibration_tilt_deg", tc["calibration"][0], "timecourse.calibration[0]"),
            ("mot_tc_seed", tc["seed"], "timecourse.seed")):
        _add(out, name, str(val) if name.endswith("_seed") else const(val), val, f"{F} :: config.{key} (= configs/g4_motion.toml)")
    _add(out, "mot_translation_largest_mm", const(max(g["translations_mm"])), max(g["translations_mm"]),
         f"{F} :: config.geometry.translations_mm (the largest displacement; = configs/g4_motion.toml)")
    _add(out, "mot_slip_largest_deg", const(max(g["slip_deg"])), max(g["slip_deg"]),
         f"{F} :: config.geometry.slip_deg (the largest cap slip; = configs/g4_motion.toml)")
    gains = [round(100 * c[1], 6) for c in cb["calibration"]]
    _add(out, "mot_calibration_gain_levels_pct", _join([const(x) + "%" for x in gains]), gains,
         f"{F} :: config.coupling.calibration[*][1] x 100 (RMS gain error; = configs/g4_motion.toml)")
    _add(out, "mot_tc_calibration_gain_pct", const(round(100 * tc["calibration"][1], 6)) + "%", tc["calibration"][1],
         f"{F} :: config.timecourse.calibration[1] x 100 (= configs/g4_motion.toml)")
    _add(out, "mot_infeasible_contact_mm", "18", 18.0, "docs/provenance_register.md :: A-MOT-GEOM 'infeasible if a magnetometer coil "
         "centre comes within 18 mm of the scalp'")
    G3F = "results/g3b/g3b_summary.json"
    enbw = _json(root / G3F)["enbw_hz"]
    asd = _json(root / "results/g4/g4_adult_summary.json")["config"]["simulation"]["opm_asd_fT_per_rtHz"]
    _add(out, "mot_opm_intrinsic_inband_ft", integer(asd * math.sqrt(enbw)), asd * math.sqrt(enbw),
         f"{G3F} :: derived: OPM white noise 15 fT/sqrt(Hz) (results/g4/g4_adult_summary.json config.simulation.opm_asd_fT_per_rtHz) "
         "x sqrt(enbw_hz), RMS in the 1-40 Hz band")

    # A. sustained displacement
    geo = m["geometry"]
    sq_rows, opm_rows = [], []
    for lab in m["anatomies"]:
        t = _tag(lab)
        for case, r in geo[lab]["squid"].items():
            base, key = f"mot_{t}_squid_{_case(case)}", f"{F} :: geometry.{lab}.squid['{case}']"
            _add(out, f"{base}_min_dist_mm", mm(r["min_dist_mm"]), r["min_dist_mm"], f"{key}.min_dist_mm (nearest magnetometer to the scalp)")
            sq_rows.append((lab, case, r))
            if not r["feasible"]:
                continue
            for k in ("known", "mismatched"):
                _add(out, f"{base}_{k}_db", db(r[k]["median_db"]), r[k]["median_db"], f"{key}.{k}.median_db (median change in detectability)")
            _add(out, f"{base}_mismatched_share_gt3db", pct(r["mismatched"]["share_loss_gt_3db"]), r["mismatched"]["share_loss_gt_3db"],
                 f"{key}.mismatched.share_loss_gt_3db (area-weighted share of the cortex)")
        for case, r in geo[lab]["opm"].items():
            base, key = f"mot_{t}_opm_{_case(case)}", f"{F} :: geometry.{lab}.opm['{case}']"
            opm_rows.append((lab, case, r))
            for k in ("known", "mismatched"):
                _add(out, f"{base}_{k}_db", db(r[k]["median_db"]), r[k]["median_db"], f"{key}.{k}.median_db")
            _add(out, f"{base}_shift_mm", mm(r["median_sensor_shift_mm"]), r["median_sensor_shift_mm"], f"{key}.median_sensor_shift_mm")
            _add(out, f"{base}_n_lifted", count(r["n_lifted"]), r["n_lifted"], f"{key}.n_lifted")
            _add(out, f"{base}_max_lift_mm", mm(r["max_lift_mm"]), r["max_lift_mm"], f"{key}.max_lift_mm")
    n_inf = sum(not r["feasible"] for _, _, r in sq_rows)
    _add(out, "mot_squid_n_cases", count(len(sq_rows)), len(sq_rows), f"{F} :: geometry.<anatomy>.squid (count over three anatomies)")
    _add(out, "mot_squid_n_infeasible", count(n_inf), n_inf, f"{F} :: geometry.<anatomy>.squid[*].feasible == false (count)")
    feas = [(lab, c, r) for lab, c, r in sq_rows if r["feasible"]]
    src = f"{F} :: geometry.<anatomy>.squid['<case>']"

    def sq(sel, k="mismatched"):
        return [r[k]["median_db"] for lab, c, r in feas if sel(lab, c)]

    def add_change(name, vals, what):
        _spread(out, f"{name}_db", vals, db, f"{src}.{what}")
        _spread(out, f"{name}_loss_db", [-v for v in vals], dbu, f"{src}.{what}; derived: loss = -median_db")

    for d in g["translations_mm"]:
        dd = const(d)
        add_change(f"mot_squid_down{dd}mm_mismatched", sq(lambda lab, c: c == f"down {dd} mm"), f"mismatched.median_db, 'down {dd} mm'")
        add_change(f"mot_squid_down{dd}mm_known", sq(lambda lab, c: c == f"down {dd} mm", "known"), f"known.median_db, 'down {dd} mm'")
        add_change(f"mot_squid_{dd}mm_all_mismatched", sq(lambda lab, c: _amount(c) == (d, "mm")),
                   f"mismatched.median_db, every feasible {dd}-mm translation (down, x, y)")
    add_change("mot_squid_10mm_lateral_mismatched", sq(lambda lab, c: _amount(c) == (10.0, "mm") and not c.startswith("down")),
               "mismatched.median_db, feasible 10-mm translations along x or y (templates only)")
    _spread(out, "mot_squid_down10mm_mismatched_share_gt3db", [r["mismatched"]["share_loss_gt_3db"] for lab, c, r in feas if c == "down 10 mm"],
            pct, f"{src}.mismatched.share_loss_gt_3db, 'down 10 mm'")
    for deg in g["rotations_deg"]:
        dg = const(deg)
        add_change(f"mot_squid_rot{dg}deg_mismatched", sq(lambda lab, c: _amount(c) == (deg, "deg")),
                   f"mismatched.median_db, every feasible {dg}-deg rotation (pitch, roll, yaw)")
    tp = [(lab, c, r) for lab, c, r in feas if lab != "adult" and c in ("pitch +10 deg", "pitch -10 deg", "roll +10 deg", "roll -10 deg")]
    add_change("mot_squid_templates_pitch_roll10deg_mismatched", [r["mismatched"]["median_db"] for _, _, r in tp],
               "mismatched.median_db, feasible 10-deg pitch and roll of the templates")
    _spread(out, "mot_squid_templates_pitch_roll10deg_share_gt3db", [r["mismatched"]["share_loss_gt_3db"] for _, _, r in tp], pct,
            f"{src}.mismatched.share_loss_gt_3db, feasible 10-deg pitch and roll of the templates")
    kn = [r["known"]["median_db"] for _, _, r in feas]
    _add(out, "mot_squid_known_max_loss_db", dbu(-min(kn)), -min(kn), f"{src}.known.median_db; derived: largest loss over every feasible case")
    src_o = f"{F} :: geometry.<anatomy>.opm['<case>']"
    for deg in g["slip_deg"]:
        dg = const(deg)
        sel = [r for _, c, r in opm_rows if _amount(c) == (deg, "deg")]
        if len(sel) != 6 * len(m["anatomies"]):
            raise ValueError(f"expected six {dg}-deg slips per anatomy, found {len(sel)}")
        mis = [r["mismatched"]["median_db"] for r in sel]
        _spread(out, f"mot_opm_slip{dg}deg_mismatched_db", mis, db, f"{src_o}.mismatched.median_db, {dg}-deg slips")
        _spread(out, f"mot_opm_slip{dg}deg_mismatched_loss_db", [-v for v in mis], dbu,
                f"{src_o}.mismatched.median_db, {dg}-deg slips; derived: loss = -median_db")
        _spread(out, f"mot_opm_slip{dg}deg_known_loss_db", [-r["known"]["median_db"] for r in sel], dbu,
                f"{src_o}.known.median_db, {dg}-deg slips; derived: loss = -median_db")
        _spread(out, f"mot_opm_slip{dg}deg_shift_mm", [r["median_sensor_shift_mm"] for r in sel], mm, f"{src_o}.median_sensor_shift_mm, {dg}-deg slips")
        _spread(out, f"mot_opm_slip{dg}deg_n_lifted", [r["n_lifted"] for r in sel], count, f"{src_o}.n_lifted, {dg}-deg slips")
        _spread(out, f"mot_opm_slip{dg}deg_max_lift_mm", [r["max_lift_mm"] for r in sel], mm, f"{src_o}.max_lift_mm, {dg}-deg slips")
    _spread(out, "mot_adult_opm_slip3deg_shift_mm", [r["median_sensor_shift_mm"] for lab, c, r in opm_rows if lab == "adult" and _amount(c) == (3.0, "deg")],
            mm, f"{F} :: geometry.adult.opm['slip <axis> +-3 deg'].median_sensor_shift_mm")
    kn_o = [r["known"]["median_db"] for _, _, r in opm_rows]
    _add(out, "mot_opm_slip_known_max_loss_db", dbu(-min(kn_o)), -min(kn_o), f"{src_o}.known.median_db; derived: largest loss over every slip")

    # B. in-band rotation in a static residual field
    cp = m["coupling"]
    for lab in m["anatomies"]:
        t = _tag(lab)
        for corr, ct in CORR.items():
            st = cp[lab][f"static/{corr}"]
            _add(out, f"mot_{t}_static_d_{ct}_db", db(st["D_db"]), st["D_db"], f"{F} :: coupling.{lab}['static/{corr}'].D_db (dense OPM vs Neuromag, no motion)")
            _add(out, f"mot_{t}_static_opm_change_{ct}_db", db(st["opm_change_db"]), st["opm_change_db"],
                 f"{F} :: coupling.{lab}['static/{corr}'].opm_change_db (the correction alone)")
        for key, c in cp[lab]["cases"].items():
            corr, cal, field, pivot = key.split("/")
            base = f"mot_{t}_{CORR[corr]}_{CAL[cal]}_{field}_{pivot}"
            ksrc = f"{F} :: coupling.{lab}.cases['{key}']"
            for lv, lt in LEVEL.items():
                v = c["unmodelled"]["thresholds_deg_unit_field"][lv]["median_curve"]
                _add(out, f"{base}_thr_{lt}_deg", "> " + const(cb["rotation_rms_deg"][-1]) if v is None else sig(v), "not reached" if v is None else v,
                     f"{ksrc}.unmodelled.thresholds_deg_unit_field.{lv}.median_curve (deg RMS per axis in a unit field; null = not reached "
                     "within the tested rotations)")
            ol = c["oracle"]["opm_change_db_median"][-1]
            _add(out, f"{base}_oracle_change5_db", db(ol), ol, f"{ksrc}.oracle.opm_change_db_median[-1] (at {const(cb['rotation_rms_deg'][-1])} deg)")
            cf = cp[lab]["coupling_fT"][key]["median"]
            _add(out, f"{base}_artefact_ft", integer(cf), cf, f"{F} :: coupling.{lab}.coupling_fT['{key}'].median (per channel, 1 deg RMS, unit field)")

    def thr(corr, cal, field, lv, pivot="neck"):
        return [cp[lab]["cases"][f"{corr}/{cal}/{field}/{pivot}"]["unmodelled"]["thresholds_deg_unit_field"][lv]["median_curve"]
                for lab in m["anatomies"]]

    src_c = f"{F} :: coupling.<anatomy>.cases['<correction>/<calibration>/<field>/neck'].unmodelled.thresholds_deg_unit_field"
    for corr, ct in CORR.items():
        for cal, kt in CAL.items():
            for field in ("uniform", "gradient"):
                for lv, lt in LEVEL.items():
                    v = thr(corr, cal, field, lv)
                    if all(x is not None for x in v):
                        _spread(out, f"mot_{ct}_{kt}_{field}_thr_{lt}_deg", v, sig,
                                f"{src_c}.{lv}.median_curve, {corr}, {cal}, {field}, over the three anatomies")
    for field in ("uniform", "gradient"):  # no correction: calibration does not matter; every calibration level pooled
        v = [x for cal in CAL for x in thr("none", cal, field, "loss_1dB")]
        _spread(out, f"mot_none_allcal_{field}_thr_loss1db_deg", v, sig, f"{src_c}.loss_1dB.median_curve, none, every calibration level")
    v = [x / 10 for x in thr("homogeneous+gradient", "tilt1deg_gain1pct", "uniform", "loss_1dB")]
    _spread(out, "mot_eightterm_cal1_uniform_thr_loss1db_deg_at_10nt", v, sig,
            f"{src_c}.loss_1dB.median_curve, homogeneous+gradient, tilt1deg_gain1pct, uniform; derived: divided by 10 (a 10-nT field)")
    ratios, below, above, piv_h, piv_8 = [], [], [], [], []
    for lab in m["anatomies"]:
        for key, c in cp[lab]["cases"].items():
            t_ = c["unmodelled"]["thresholds_deg_unit_field"]
            l1, d0 = t_["loss_1dB"]["median_curve"], t_["D_0dB"]["median_curve"]
            if key.endswith("/neck") and l1 and d0:
                ratios.append(d0 / l1)
            mc, q = t_["loss_1dB"]["median_curve"], t_["loss_1dB"]["per_draw_p10_p90"]
            if mc:
                if q[0] is not None:
                    below.append(1 - q[0] / mc)
                if q[1] is not None:
                    above.append(q[1] / mc - 1)
            corr, cal, field, pivot = key.split("/")
            if pivot == "origin" and cal != "tilt0deg_gain0pct":
                neck = cp[lab]["cases"][f"{corr}/{cal}/{field}/neck"]["unmodelled"]["thresholds_deg_unit_field"]["loss_1dB"]["median_curve"]
                if l1 and neck:
                    (piv_h if corr == "homogeneous" else piv_8).append(l1 / neck - 1)
    _spread(out, "mot_d0_over_loss1db_ratio", ratios, ratio, f"{src_c}; derived: D_0dB / loss_1dB median-curve thresholds (neck pivot)")
    _spread(out, "mot_thr_p10_below_median_pct", below, pct,
            f"{src_c}.loss_1dB.per_draw_p10_p90[0]; derived: 1 - p10 / median_curve, every case and anatomy (draws beyond the grid excluded)")
    _spread(out, "mot_thr_p90_above_median_pct", above, pct,
            f"{src_c}.loss_1dB.per_draw_p10_p90[1]; derived: p90 / median_curve - 1, every case and anatomy (draws beyond the grid excluded)")
    _add(out, "mot_pivot_origin_change_homog_max_pct", pct(max(abs(x) for x in piv_h)), max(abs(x) for x in piv_h),
         f"{src_c}.loss_1dB.median_curve; derived: largest |origin / neck - 1|, homogeneous projection with calibration errors")
    _spread(out, "mot_pivot_origin_change_eightterm_pct", piv_8, pct,
            f"{src_c}.loss_1dB.median_curve; derived: origin / neck - 1, 8-term projection with calibration errors (reached cases)")
    orc = [(c["oracle"]["opm_change_db_median"][-1], lab, key) for lab in m["anatomies"] for key, c in cp[lab]["cases"].items()]
    worst = min(orc)
    _add(out, "mot_oracle_change5_worst_db", db(worst[0]), worst[0],
         f"{F} :: coupling.{worst[1]}.cases['{worst[2]}'].oracle.opm_change_db_median[-1] (smallest over every case and anatomy)")
    _add(out, "mot_oracle_loss5_max_db", dbu(-worst[0]), -worst[0], f"{F} :: coupling.*.cases[*].oracle.opm_change_db_median[-1]; derived: largest loss")
    for nm, key in (("mot_artefact_none_uniform", "none/tilt0deg_gain0pct/uniform/neck"),
                    ("mot_artefact_homog_cal1_uniform", "homogeneous/tilt1deg_gain1pct/uniform/neck"),
                    ("mot_artefact_eightterm_cal1_gradient", "homogeneous+gradient/tilt1deg_gain1pct/gradient/neck"),
                    ("mot_artefact_eightterm_cal1_uniform", "homogeneous+gradient/tilt1deg_gain1pct/uniform/neck")):
        v = [cp[lab]["coupling_fT"][key]["median"] for lab in m["anatomies"]]
        _spread(out, f"{nm}_ft", v, integer, f"{F} :: coupling.<anatomy>.coupling_fT['{key}'].median over the three anatomies")
        _spread(out, f"{nm}_pt", [x / 1000 for x in v], sig, f"{F} :: coupling.<anatomy>.coupling_fT['{key}'].median / 1000, three anatomies")
    for corr, ct in (("homogeneous", "homog"), ("homogeneous+gradient", "eightterm")):
        _spread(out, f"mot_static_d_{ct}_db", [cp[lab][f"static/{corr}"]["D_db"] for lab in m["anatomies"]], db,
                f"{F} :: coupling.<anatomy>['static/{corr}'].D_db over the three anatomies")
        _spread(out, f"mot_templates_static_d_after_projections_{ct}_db",
                [cp[lab][f"static/{corr}"]["D_db"] for lab in m["anatomies"] if lab != "adult"], db,
                f"{F} :: coupling.<template>['static/{corr}'].D_db, 24- and 12-month templates")
    _spread(out, "mot_adult_static_d_after_projections_db", [cp["adult"][f"static/{c}"]["D_db"] for c in ("homogeneous", "homogeneous+gradient")],
            db, f"{F} :: coupling.adult['static/homogeneous'|'static/homogeneous+gradient'].D_db")
    _spread(out, "mot_templates_static_d_after_projections_db",
            [cp[lab][f"static/{c}"]["D_db"] for lab in m["anatomies"] if lab != "adult" for c in ("homogeneous", "homogeneous+gradient")],
            db, f"{F} :: coupling.<template>['static/homogeneous'|'static/homogeneous+gradient'].D_db")

    # B'. exact rigid motion over 60 s
    T = m["timecourse"]
    for corr, ct in CORR.items():
        r = T["corrections"][corr]
        _add(out, f"mot_tc_artefact_{ct}_ft", integer(r["inband_rms_fT_median"]), r["inband_rms_fT_median"],
             f"{F} :: timecourse.corrections['{corr}'].inband_rms_fT_median (exact, median over channels)")
        _add(out, f"mot_tc_linear_{ct}_ft", integer(r["linear_prediction_fT_median"]), r["linear_prediction_fT_median"],
             f"{F} :: timecourse.corrections['{corr}'].linear_prediction_fT_median")
        _add(out, f"mot_tc_linear_max_deviation_{ct}_pct", pct(r["exact_over_linear_max_abs_deviation"]), r["exact_over_linear_max_abs_deviation"],
             f"{F} :: timecourse.corrections['{corr}'].exact_over_linear_max_abs_deviation")
    dev = [T["corrections"][c]["exact_over_linear_max_abs_deviation"] for c in CORR]
    _spread(out, "mot_tc_linear_max_deviation_pct", dev, pct, f"{F} :: timecourse.corrections[*].exact_over_linear_max_abs_deviation")
    p90 = max(max(abs(T["corrections"][c]["exact_over_linear_percentiles"][q] - 1) for q in ("p5", "p95")) for c in CORR)
    _add(out, "mot_tc_linear_p5_p95_deviation_pct", pct(p90), p90,
         f"{F} :: timecourse.corrections[*].exact_over_linear_percentiles p5/p95; derived: largest |ratio - 1| (90 % of channels within)")
    pk = T["peak_field_change_pT"]
    _add(out, "mot_tc_peak_field_change_median_pt", integer(pk["median"]), pk["median"], f"{F} :: timecourse.peak_field_change_pT.median")
    _add(out, "mot_tc_peak_field_change_max_pt", integer(pk["max"]), pk["max"], f"{F} :: timecourse.peak_field_change_pT.max")
    _spread(out, "mot_tc_inband_rotation_measured_deg", T["inband_rotation_rms_deg"], sig,
            f"{F} :: timecourse.inband_rotation_rms_deg (per axis, measured in band)")


# ----------------------------------------------------------------------------------------------
def facts(root: Path = ROOT) -> dict:
    """Every G4 fact: name -> {"value", "raw", "source"}."""
    root = Path(root)
    out: dict = {}
    _detection(root, out)
    _localization(root, out)
    _motion(root, out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grep", default="", help="only facts whose name contains this text")
    args = ap.parse_args(argv)
    f = facts(ROOT)
    for k, v in f.items():
        if args.grep in k:
            print(f"{k}\t{v['value']}\t{v['source']}")
    print(f"{len(f)} facts")


if __name__ == "__main__":
    main()
