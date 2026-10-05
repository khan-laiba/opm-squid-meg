#!/usr/bin/env python3
"""Report facts for the referee round-1 revision analyses (prefix rev_), merged by scripts/report_facts.py.

facts(root) -> {name: {"value": text as printed, "raw": unrounded number(s) or text, "source": "results/<file> :: <key path>"
or "results/<file> :: derived: <how>"}}. Read only; nothing is re-run:
  rev_ns_   results/g2_noise_sensitivity/noise_sensitivity_summary.json: the adult noise model's sensitivity (Fable major
            issue 1 (a)-(d), Codex major issue 1). sweep: the OPM white-noise sweep with 1,000-resample parcel-bootstrap
            intervals and the break-even white levels; near: the cortex within 4 mm of the inner skull in the background and
            as targets, in three forward models; col: coloured (1/f + white) OPM noise with swept corners; ff: far-field
            cardiac and ocular sources and the energy the projection leaves of them; joint: all alternatives together; alt:
            the range over the run's own headline table. Its stored configuration must equal configs/g2_noise_sensitivity.toml.
  rev_cov_  results/g2_covariance_validation/covariance_validation.json: the noise model against the measured Neuromag
            covariance of the MNE sample recording (Codex major issue 1, Fable major issue 1): what is fitted and what is
            predicted, channel variances (held-out halves), channel patterns, correlation against distance, eigenstructure,
            the heart, sub-bands, Neuromag detectability with the measured covariance, the implied OPM/Neuromag scenarios.
  rev_qc_   results/g3b_children_qc/children_qc.json: the school-aged children's surfaces against their MRIs, with the adult
            and the 2-year template as references (Fable weakness (i), Codex major issue 2); configs/school_subjects_qc_manifest.json.
  Reserved for analyses whose results are not yet committed (marked section at the end): rev_cgap_ (the gap-matched
  counterfactual helmet) and rev_conf_ (the confirmatory spike run).

Name tokens. Arrays dense (208 sites), matched (98), neuromag (the 306-channel array where a quantity is per array);
Neuromag comparators combined (all channels), mag, grad; conditions ib (sensor + cortical background), proj (+ room field,
8-term projection), ibenv (sensor + background + room field, not projected: the covariance validation's like-for-like
condition). A comparison is <array>_vs_<comparator>_<condition> with _ratio (median over targets of d_OPM/d_Neuromag), _ci
(95 % interval), _share (share of the targets with the OPM ahead). OPM white levels asd<level>ft (fT/sqrt(Hz)); 1/f corners
c<f>hz; sub-bands sb<lo>_<hi> (Hz) and b<lo>_<hi> (the covariance validation's bands); depth bins <lo>_<hi> (mm below the
scalp); distance bands from the inner skull band_<lo>_<hi> (mm); anatomies adult, infant2yr, child_a, child_b, child_c.
The covariance validation's OPM/Neuromag scenarios are rev_cov_opm_<array>_vs_<comparator>_<scenario> (s1 to s5,
s4_uncorrected, model306, model305, corr5mm_model, ...), in the condition ibenv unless the token says ib (model306_ib) or
proj (proj_...); rev_cov_det_<comparator>_<what> are Neuromag's detectability ratios, measured (or hybrid) over modelled.
Formats as scripts/report_facts_g12.py, whose helpers are imported: ratios 2 decimals, with _3dp twins for every
comparison against combined Neuromag (quoted at 3 decimals) and for any ratio near 1 (0.95-1.05); intervals "[lo, hi]" in
_ci facts; percentages integer from 10 %, 1 decimal below, 2 significant digits below 0.1 %; changes signed with 1 decimal
(+0.7%, −1.3%); mm 1 decimal (_2dp twins where finer differences matter), signed where a direction is meant (+1.3: the MRI
head boundary lies outside the surface); fT levels 1 decimal, field RMS in fT as integers; correlations 2 decimals; counts
with thousands separators; diagnostic ratios 2 significant digits (integers from 100); ranges "lo to hi" from unrounded
values with _min and _max; U+2212 for negatives; text values as stored. "raw" holds the unrounded value in the printed
unit (percent for % facts).

Usage: .venv/bin/python scripts/report_facts_rev.py [--grep TEXT]   (prints name, value, source)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_g12():
    spec = importlib.util.spec_from_file_location("report_facts_g12", Path(__file__).with_name("report_facts_g12.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G12 = _load_g12()  # the formatting helpers and the Facts container of the G1/G2 module
MINUS, num, signed, sig2, pct, count, sci, r2, r3, tag = (G12.MINUS, G12.num, G12.signed, G12.sig2, G12.pct, G12.count,
                                                         G12.sci, G12.r2, G12.r3, G12.tag)

NS = "results/g2_noise_sensitivity/noise_sensitivity_summary.json"
COV = "results/g2_covariance_validation/covariance_validation.json"
QC = "results/g3b_children_qc/children_qc.json"
CFG_NS = "configs/g2_noise_sensitivity.toml"
QC_MANIFEST = "configs/school_subjects_qc_manifest.json"

ARR = {"opm_dense": "dense", "opm_matched": "matched"}
ARR3 = {"squid": "neuromag", "opm_matched": "matched", "opm_dense": "dense"}
COND = {"intrinsic+brain": "ib", "projected": "proj", "intrinsic+brain+env": "ibenv"}
COND_DESC = {"intrinsic+brain": "sensor + cortical background", "projected": "+ room field, 8-term projection",
             "intrinsic+brain+env": "sensor + background + room field, not projected"}


# ------------------------------------------------------------------------------------------------
# formatting (beyond the imported helpers)
def g(x) -> str:
    """A declared value as written ('15', '0.5'), U+2212 for a negative."""
    return f"{x:g}".replace("-", MINUS)


def spct(p) -> str:
    """A signed percentage change, 1 decimal: '+0.7%', '−1.3%', '0.0%' (no sign when it rounds to zero)."""
    s = f"{p:+.1f}"
    return "0.0%" if float(s) == 0 else s.replace("-", MINUS) + "%"


def mm(x) -> str:
    return num(x, 1)


def mm2(x) -> str:
    return num(x, 2)


def smm(x) -> str:
    return signed(x, 1)


def smm2(x) -> str:
    return signed(x, 2)


def pctk(p) -> str:
    """A percentage as pct(), with thousands separators from 1,000 %."""
    return f"{p:,.0f}%".replace("-", MINUS) if abs(p) >= 1000 else pct(p)


def small(x) -> str:
    """A tiny check value: '0' or scientific notation ('2.7 × 10⁻⁵')."""
    return "0" if x == 0 else sci(x)


def diag(x) -> str:
    """A diagnostic ratio: 2 significant digits, integers with separators from 100."""
    return count(x) if abs(x) >= 100 else sig2(x)


def yes(b) -> str:
    return "yes" if b else "no"


def listing(values, f=g) -> str:
    """'7, 10, 15 and 20'; 'none' for an empty list."""
    items = [f(v) for v in values]
    if not items:
        return "none"
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def ci_text(lo, hi, f=r2) -> str:
    return f"[{f(lo)}, {f(hi)}]"


def ctok(key: str) -> str:
    """'opm_dense/combined/intrinsic+brain' -> 'dense_vs_combined_ib'."""
    arr, comp, cond = key.split("/")
    return f"{ARR[arr]}_vs_{comp}_{COND[cond]}"


def band_tok(lo, hi) -> str:
    return f"{tag(lo)}_{tag(hi)}"


def parcels_of(ci_method: str) -> str:
    m = re.fullmatch(r"parcels \((\d+)\)", ci_method)
    if not m:
        raise ValueError(f"unexpected ci_method {ci_method!r}")
    return m.group(1)


class Facts(G12.Facts):
    """G12's container, with the prefix and the source form of this module enforced."""

    def add(self, name, value, raw, source):
        if not name.startswith("rev_"):
            raise ValueError(f"fact {name!r} lacks the rev_ prefix")
        if str(value).strip() == "":
            raise ValueError(f"fact {name!r} has an empty value")
        if " :: " not in source:
            raise ValueError(f"fact {name!r}: source lacks ' :: <key path>'")
        super().add(name, value, raw, source)


def load(root, rel):
    return json.loads((Path(root) / rel).read_text())


# ------------------------------------------------------------------------------------------------
# noise-model sensitivity (rev_ns_)
def ns_cmp(F, d, base, e, path, key, parcels=False, note=""):
    """One stored comparison (ratio, ci95_ratio as stored, checked to be 2**median_log2 and 2**ci95) and its share(s)."""
    r, ci = e["ratio"], e["ci95_ratio"]
    if abs(math.log2(r) - e["median_log2"]) > 1e-9 or any(abs(math.log2(c) - lg) > 1e-9 for c, lg in zip(ci, e["ci95"])):
        raise ValueError(f"{NS} :: {path}: ratio or ci95_ratio is not 2**median_log2 or 2**ci95")
    cond = key.split("/")[2]
    src = (f"{NS} :: {path} (ratio, ci95_ratio: median over {e['n']:,} targets of d_OPM/d_Neuromag, {COND_DESC[cond]}"
           f"{'; ' + note if note else ''}; 95 % CI from {d['setup']['n_boot']:,} parcel-bootstrap resamples over "
           f"{parcels_of(e['ci_method'])} parcel labels)")
    F.ratio(base, r, src, ci, three="/combined/" in key)
    share = 100 * e["share_opm_better"]
    F.add(f"{base}_share", pct(share), share, f"{NS} :: {path}.share_opm_better (share of the targets with the OPM ahead)")
    if pct(share) == "100%" and share < 100:  # 'all but k': the count the rounded share hides
        k = round((1 - e["share_opm_better"]) * e["n"])
        F.add(f"{base}_n_neuromag_ahead", count(k), k, f"{NS} :: derived: (1 - {path}.share_opm_better) x n ({e['n']:,} targets)")
    if parcels:
        F.add(f"{base}_share_parcels", pct(100 * e["share_parcels_opm_better"]), 100 * e["share_parcels_opm_better"],
              f"{NS} :: {path}.share_parcels_opm_better (share of the parcel labels whose median favours the OPM)")


def depth_summary(F, base, bins, path):
    """A stored depth profile: the shallowest populated bin (ratio, CI, n), the span of the populated bins and the bin
    from which every deeper populated bin has a ratio below 1."""
    pop = [b for b in bins if b.get("ratio") is not None]
    if not pop:
        return
    b0 = pop[0]
    i0 = bins.index(b0)
    bt = f"{b0['lo']:.0f}_{b0['hi']:.0f}"
    src = (f"{NS} :: {path}[{i0}] (bin {b0['lo']:g}-{b0['hi']:g} mm below the scalp, the shallowest populated: ratio, "
           f"ci95_ratio of its {b0['n']:,} targets, parcel bootstrap)")
    F.ratio(f"{base}_depth_{bt}", b0["ratio"], src, b0["ci95_ratio"], three=True)
    F.add(f"{base}_depth_{bt}_n", count(b0["n"]), b0["n"], f"{NS} :: {path}[{i0}].n")
    src = f"{NS} :: {path}[*].ratio (median ratio per 5-mm depth bin, the populated bins)"
    F.span(f"{base}_depth_bins", [b["ratio"] for b in pop], r2, src)
    below = [i for i, b in enumerate(pop) if b["ratio"] < 1]
    F.add(f"{base}_depth_n_bins_below1", count(len(below)), len(below), f"{src}; bins with ratio < 1 of {len(pop)}")
    if below and below == list(range(below[0], len(pop))):
        b = pop[below[0]]
        src = f"{NS} :: derived: {path}, the shallowest bin from which every deeper populated bin has ratio < 1"
        F.add(f"{base}_depth_below1_from_mm", f"{b['lo']:.0f} to {b['hi']:.0f}", [b["lo"], b["hi"]], src)
        F.add(f"{base}_depth_below1_from_bin", f"{b['lo']:.0f}–{b['hi']:.0f}", [b["lo"], b["hi"]], src + " (bin label, mm)")


def ns_setup_facts(F, root, d):
    cfg = d["config"]["sensitivity"]
    if tomllib.loads((Path(root) / CFG_NS).read_text()) != cfg:
        raise ValueError(f"{NS}: config.sensitivity differs from {CFG_NS}: the stored run is not the declared one")
    s, p = d["setup"], f"{NS} :: setup"
    if s["n_boot"] != cfg["general"]["n_boot"]:
        raise ValueError(f"{NS}: setup.n_boot differs from the declared general.n_boot")
    F.add("rev_ns_n_targets", count(s["n_targets"]), s["n_targets"], f"{p}.n_targets (cortical targets, those of G2)")
    F.add("rev_ns_n_background", count(s["n_background_grid"]), s["n_background_grid"],
          f"{p}.n_background_grid (background sources on the 7-mm grid)")
    F.add("rev_ns_n_boot", count(s["n_boot"]), s["n_boot"], f"{p}.n_boot (parcel-bootstrap resamples; {CFG_NS} general.n_boot)")
    F.add("rev_ns_boot_seed", str(cfg["general"]["bootstrap_seed"]), cfg["general"]["bootstrap_seed"],
          f"{NS} :: config.sensitivity.general.bootstrap_seed (the run's own generator)")
    methods = {e["ci_method"] for x in d["sweep"]["entries"].values() for e in x["comparisons"].values()}
    if len(methods) != 1:
        raise ValueError(f"{NS}: the sweep's comparisons use several ci_methods {methods}")
    n = parcels_of(methods.pop())
    F.add("rev_ns_n_parcels", n, int(n), f"{NS} :: sweep.entries[*].comparisons[*].ci_method ('parcels ({n})': the "
          "Desikan-Killiany parcels and the two medial-wall labels resampled by the bootstrap)")
    F.add("rev_ns_enbw_hz", num(s["enbw_hz"], 1), s["enbw_hz"], f"{p}.enbw_hz (equivalent noise bandwidth of the 1-40 Hz filter)")
    for arr, k in s["arrays"].items():
        F.add(f"rev_ns_n_channels_{ARR3[arr]}", count(k), k, f"{p}.arrays.{arr} (channels)")
    prim = d["config"]["adult"]["sensors"]["opm_asd_primary_fT_per_rtHz"]
    F.add("rev_ns_primary_asd", g(prim), prim, f"{NS} :: config.adult.sensors.opm_asd_primary_fT_per_rtHz (fT/sqrt(Hz); "
          "configs/g2_adult.toml)")
    F.add("rev_ns_commit", d["provenance"]["commit"], d["provenance"]["commit"], f"{NS} :: provenance.commit (code of the run)")


def ns_sweep_facts(F, d):
    s = d["sweep"]
    levels = s["levels_fT"]
    if levels != d["config"]["adult"]["sensors"]["opm_asd_fT_per_rtHz"]:
        raise ValueError(f"{NS}: sweep.levels_fT differ from config.adult.sensors.opm_asd_fT_per_rtHz")
    F.add("rev_ns_sweep_levels", listing(levels), levels, f"{NS} :: sweep.levels_fT (OPM white levels, fT/sqrt(Hz))")
    F.add("rev_ns_sweep_levels_range", f"{g(min(levels))} to {g(max(levels))}", [min(levels), max(levels)], f"{NS} :: sweep.levels_fT")
    x = s["reproduces_stored_g2_max_abs_diff_median_log2"]
    F.add("rev_ns_sweep_g2_max_diff_log2", g(x), x, f"{NS} :: sweep.reproduces_stored_g2_max_abs_diff_median_log2 (largest "
          f"difference of a median log2 ratio from results/g2/g2_summary.json, stored at commit {s['stored_g2_commit']})")
    for lvl in levels:
        lk = f"{lvl:g}"
        e = s["entries"][lk]
        for key, c in e["comparisons"].items():
            ns_cmp(F, d, f"rev_ns_sweep_asd{lk}ft_{ctok(key)}", c, f"sweep.entries['{lk}'].comparisons['{key}']", key, parcels=True)
        for key, bins in e["depth"].items():
            depth_summary(F, f"rev_ns_sweep_asd{lk}ft_{ctok(key)}", bins, f"sweep.entries['{lk}'].depth['{key}']")
    prim = f"{d['config']['adult']['sensors']['opm_asd_primary_fT_per_rtHz']:g}"
    for key, rg in s["range"].items():
        tok = ctok(key)
        ratios = [s["entries"][f"{lv:g}"]["comparisons"][key]["ratio"] for lv in levels]
        stored = (rg["ratio_at_min_noise"], rg["ratio_at_max_noise"], rg["ratio_at_primary"])
        if stored != (ratios[0], ratios[-1], s["entries"][prim]["comparisons"][key]["ratio"]):
            raise ValueError(f"{NS}: sweep.range['{key}'] does not match the sweep entries")
        src = f"{NS} :: sweep.entries[*].comparisons['{key}'].ratio (over the white levels {listing(levels)} fT/sqrt(Hz))"
        F.span(f"rev_ns_sweep_{tok}", ratios, r2, src, extra=(("3dp", r3),))
        lo, hi = rg["ci95_envelope"]
        src = f"{NS} :: sweep.range['{key}'].ci95_envelope (lowest lower and highest upper 95 % bound over the sweep)"
        F.add(f"rev_ns_sweep_{tok}_envelope_ci", ci_text(lo, hi), [lo, hi], src)
        F.add(f"rev_ns_sweep_{tok}_envelope_ci_3dp", ci_text(lo, hi, r3), [lo, hi], src)
        for side, k in (("above1", "levels_with_ci_above_1"), ("below1", "levels_with_ci_below_1")):
            src = f"{NS} :: sweep.range['{key}'].{k} (white levels, fT/sqrt(Hz), whose 95 % CI lies {side[:5]} 1)"
            F.add(f"rev_ns_sweep_{tok}_levels_ci_{side}", listing(rg[k]), rg[k], src)
            F.add(f"rev_ns_sweep_{tok}_n_levels_ci_{side}", count(len(rg[k])), len(rg[k]), src + f"; of {len(levels)}")
    # break-even OPM white levels (median ratio = 1)
    be = s["break_even"]
    grid = be["grid_fT"]
    src = f"{NS} :: sweep.break_even.grid_fT (log-spaced OPM white levels, fT/sqrt(Hz))"
    F.add("rev_ns_breakeven_grid_range", f"{g(grid[0])} to {g(grid[-1])}", [grid[0], grid[-1]], src)
    F.add("rev_ns_breakeven_grid_points", count(len(grid)), len(grid), src)
    outside, by_arr = [], {}
    for key, e in be["entries"].items():
        tok, v, ci = ctok(key), e["break_even_fT"], e["ci95_fT"]
        if v is None or e["bound"] is not None or ci is None:
            raise ValueError(f"{NS}: sweep.break_even.entries['{key}'] does not cross 1 inside the grid")
        src = (f"{NS} :: sweep.break_even.entries['{key}'] (OPM white level, fT/sqrt(Hz), at which the median ratio is 1: "
               "log-linear interpolation on the grid; 95 % CI from the same parcel resamples applied to the whole curve)")
        F.add(f"rev_ns_breakeven_{tok}_asd", num(v, 1), v, src)
        F.add(f"rev_ns_breakeven_{tok}_asd_2dp", num(v, 2), v, src)
        F.add(f"rev_ns_breakeven_{tok}_asd_ci", ci_text(ci[0], ci[1], lambda x: num(x, 1)), list(ci), src)
        F.add(f"rev_ns_breakeven_{tok}_asd_ci_2dp", ci_text(ci[0], ci[1], lambda x: num(x, 2)), list(ci), src)
        outside.append(e["resamples_outside_grid"])
        by_arr.setdefault(key.split("/")[0], []).append(v)
    F.add("rev_ns_breakeven_resamples_outside_grid_max", count(max(outside)), max(outside),
          f"{NS} :: sweep.break_even.entries[*].resamples_outside_grid (largest over the comparisons)")
    for arr, vals in by_arr.items():
        F.span(f"rev_ns_breakeven_{ARR[arr]}_asd", vals, lambda x: num(x, 1),
               f"{NS} :: sweep.break_even.entries['{arr}/*'].break_even_fT (three comparators, two conditions; fT/sqrt(Hz))",
               extra=(("int", lambda x: num(x, 0)),))


NEAR_MODELS = {"primary": "primary", "one_layer_coarse": "coarse", "one_layer_refined": "refined"}
NEAR_VARIANTS = {"background": "bg", "targets": "targets", "both": "both"}


def floor_tok(floor_mm) -> str:
    return f"f{floor_mm:g}mm".replace(".", "p")


def ns_near_facts(F, d):
    n, cfg, p = d["near_skull"], d["config"]["sensitivity"]["near_skull"], f"{NS} :: near_skull"
    F.add("rev_ns_near_rule_mm", g(cfg["rule_mm"]), cfg["rule_mm"], f"{NS} :: config.sensitivity.near_skull.rule_mm "
          f"(A-BEM-DIST: {n['rule']})")
    if n["n_near_vertices_used"] != n["n_near_vertices"] or n["n_near_targets_used"] != n["n_near_oct6_targets"]:
        raise ValueError(f"{NS}: not every near-skull vertex or target was used")
    F.add("rev_ns_near_n_vertices", count(n["n_near_vertices"]), n["n_near_vertices"], f"{p}.n_near_vertices")
    F.add("rev_ns_near_share_valid_pct", pct(100 * n["share_of_valid_vertices"]), 100 * n["share_of_valid_vertices"],
          f"{p}.share_of_valid_vertices (share of the valid white-surface vertices)")
    F.add("rev_ns_near_area_cm2", count(n["near_area_cm2"]), n["near_area_cm2"], f"{p}.near_area_cm2")
    F.add("rev_ns_near_area_cm2_1dp", num(n["near_area_cm2"], 1), n["near_area_cm2"], f"{p}.near_area_cm2")
    F.add("rev_ns_near_valid_area_cm2", count(n["valid_area_cm2"]), n["valid_area_cm2"], f"{p}.valid_area_cm2 (all valid cortex)")
    share = 100 * n["near_area_cm2"] / n["valid_area_cm2"]
    F.add("rev_ns_near_area_pct", pct(share), share, f"{NS} :: derived: near_skull.near_area_cm2 / near_skull.valid_area_cm2")
    F.add("rev_ns_near_n_targets", count(n["n_near_oct6_targets"]), n["n_near_oct6_targets"],
          f"{p}.n_near_oct6_targets (near-skull vertices of the G2 target grid)")
    for i, b in enumerate(n["vertices_by_band"]):
        bt = band_tok(b["lo_mm"], b["hi_mm"])
        src = f"{p}.vertices_by_band[{i}] (vertices {b['lo_mm']:g}-{b['hi_mm']:g} mm from the inner skull)"
        F.add(f"rev_ns_near_band_{bt}_label", f"{g(b['lo_mm'])}–{g(b['hi_mm'])}", [b["lo_mm"], b["hi_mm"]], src + " (label, mm)")
        F.add(f"rev_ns_near_band_{bt}_n_vertices", count(b["n"]), b["n"], src + ".n")
        F.add(f"rev_ns_near_band_{bt}_area_cm2", num(b["area_cm2"], 1), b["area_cm2"], src + ".area_cm2")
    sub2 = [b for b in n["vertices_by_band"] if b["hi_mm"] <= 2.0]
    src = f"{NS} :: derived: near_skull.vertices_by_band, the bands closer than 2 mm to the inner skull (summed)"
    F.add("rev_ns_near_below2mm_n_vertices", count(sum(b["n"] for b in sub2)), sum(b["n"] for b in sub2), src)
    F.add("rev_ns_near_below2mm_area_cm2", num(sum(b["area_cm2"] for b in sub2), 1), sum(b["area_cm2"] for b in sub2), src)
    floors = {}
    for mk, m in NEAR_MODELS.items():
        blk = n[mk]
        F.add(f"rev_ns_near_{m}_model", blk["model"], blk["model"], f"{p}.{mk}.model")
        for key, e in blk["baseline"]["comparisons"].items():
            if m == "primary" or "/combined/" in key:
                ns_cmp(F, d, f"rev_ns_near_{m}_baseline_{ctok(key)}", e, f"near_skull.{mk}.baseline.comparisons['{key}']", key)
        base_rms = None
        for fk in [k for k in blk if k.startswith("floor_")]:
            f = blk[fk]
            ft = floor_tok(f["floor_mm"])
            geo = (f["near_vertices"], f["near_area_cm2"], f["near_targets"])
            if floors.setdefault(ft, (geo, f"{mk}.{fk}"))[0] != geo:
                raise ValueError(f"{NS}: near_skull floors of {ft} differ between the models")
            bg = f["background"]
            if base_rms is None:
                base_rms = bg["brain_rms_baseline_fT"]
                for arr, v in base_rms.items():
                    F.add(f"rev_ns_near_{m}_baseline_brain_rms_{ARR3[arr]}_ft", num(v, 0), v, f"{p}.{mk}.{fk}.background."
                          f"brain_rms_baseline_fT.{arr} (median in-band brain-noise RMS without the near-skull cortex, fT)")
            elif bg["brain_rms_baseline_fT"] != base_rms:
                raise ValueError(f"{NS}: near_skull.{mk} baselines differ between floors")
            F.ratio(f"rev_ns_near_{m}_{ft}_bg_scale", bg["brain_scale_over_baseline"], f"{p}.{mk}.{fk}.background."
                    "brain_scale_over_baseline (cortical background scale recalibrated on the gradiometers, over the baseline's)",
                    three=True)
            for arr, v in bg["brain_rms_fT"].items():
                F.add(f"rev_ns_near_{m}_{ft}_brain_rms_{ARR3[arr]}_ft", num(v, 0), v,
                      f"{p}.{mk}.{fk}.background.brain_rms_fT.{arr} (median in-band brain-noise RMS, fT)")
            variants = (("bg", bg["comparisons"], "background.comparisons"), ("targets", f["targets"]["all"], "targets.all"),
                        ("targets_near", f["targets"]["near_only"], "targets.near_only"), ("both", f["both"]["all"], "both.all"),
                        ("both_near", f["both"]["near_only"], "both.near_only"))
            for vt, comps, path in variants:
                for key, e in comps.items():
                    if m == "primary" or "/combined/" in key:
                        ns_cmp(F, d, f"rev_ns_near_{m}_{ft}_{vt}_{ctok(key)}", e, f"near_skull.{mk}.{fk}.{path}['{key}']", key)
            for dv, profiles in f.get("depth", {}).items():
                for key, bins in profiles.items():
                    depth_summary(F, f"rev_ns_near_{m}_{ft}_{NEAR_VARIANTS.get(dv, dv)}_{ctok(key)}", bins,
                                  f"near_skull.{mk}.{fk}.depth.{dv}['{key}']")
    for ft, ((nv, area, nt), where) in floors.items():
        src = f"{p}.{where} (near-skull cortex used from this distance to the inner skull on)"
        F.add(f"rev_ns_near_{ft}_n_vertices", count(nv), nv, src + ".near_vertices")
        F.add(f"rev_ns_near_{ft}_area_cm2", count(area), area, src + ".near_area_cm2")
        F.add(f"rev_ns_near_{ft}_area_cm2_1dp", num(area, 1), area, src + ".near_area_cm2")
        F.add(f"rev_ns_near_{ft}_n_targets", count(nt), nt, src + ".near_targets")
    # inclusion effects: change of the median log2 ratio by the near-skull cortex, every model and floor
    largest = (0.0, "")
    for key, x in n["inclusion_effect_log2"].items():
        tok, groups = ctok(key), {}
        for k, v in x.items():
            mk, fk, var = k.split("/")
            ft, vt = floor_tok(float(fk.split("_")[1].removesuffix("mm"))), NEAR_VARIANTS[var]
            ch = 100 * (2 ** v - 1)
            F.add(f"rev_ns_near_effect_{NEAR_MODELS[mk]}_{ft}_{vt}_{tok}_pct", spct(ch), ch,
                  f"{p}.inclusion_effect_log2['{key}']['{k}'] (median log2 ratio with minus without the near-skull cortex; "
                  "% = 100 (2**x - 1), the change of the median ratio)")
            groups.setdefault(vt, []).append(ch)
            if abs(ch) > abs(largest[0]):
                largest = (ch, f"{key} / {k}")
        for vt, vals in groups.items():
            F.span(f"rev_ns_near_effect_{vt}_{tok}_pct", vals, spct, f"{p}.inclusion_effect_log2['{key}'][*/{vt}] (% change "
                   "of the median ratio, over the three models and their floors)")
    F.add("rev_ns_near_effect_max_abs_pct", pct(abs(largest[0])), abs(largest[0]), f"{NS} :: derived: largest |100 (2**x - 1)| "
          f"over near_skull.inclusion_effect_log2 (every comparison, model, floor and variant; at {largest[1]})")
    # lead-field accuracy near the inner skull
    F.add("rev_ns_near_nbr_flag", g(cfg["neighbour_energy_flag"]), cfg["neighbour_energy_flag"],
          f"{NS} :: config.sensitivity.near_skull.neighbour_energy_flag (a vertex's lead-field energy over the median of its mesh "
          "neighbours' above which it counts as anomalous)")
    nbr = {"primary": n["primary_neighbour_energy"], "coarse": n["one_layer_coarse"]["neighbour_energy"],
           "refined": n["one_layer_refined"]["neighbour_energy"]}
    paths = {"primary": "primary_neighbour_energy", "coarse": "one_layer_coarse.neighbour_energy",
             "refined": "one_layer_refined.neighbour_energy"}
    for m, per_arr in nbr.items():
        for arr, rows in per_arr.items():
            for i, b in enumerate(rows):
                base, src = f"rev_ns_near_nbr_{m}_{ARR3[arr]}_band_{band_tok(b['lo_mm'], b['hi_mm'])}", f"{p}.{paths[m]}.{arr}[{i}]"
                F.add(f"{base}_n_flagged", count(b["n_above_flag"]), b["n_above_flag"], f"{src}.n_above_flag (of {b['n']:,} vertices)")
                if m == "primary":
                    F.add(f"{base}_n", count(b["n"]), b["n"], f"{src}.n")
                    for st in ("median", "p99", "max"):
                        F.add(f"{base}_{st}", diag(b[st]), b[st], f"{src}.{st} (lead-field energy over the neighbours' median)")
    for arr, rows in n["gain_change_coarse_to_refined"].items():
        for i, b in enumerate(rows):
            base = f"rev_ns_near_gainchange_{ARR3[arr]}_band_{band_tok(b['lo_mm'], b['hi_mm'])}"
            src = (f"{p}.gain_change_coarse_to_refined.{arr}[{i}] (relative change of the noise-normalised lead field from the "
                   "5,120- to the 20,480-triangle inner skull, 1-layer BEM, %)")
            F.add(f"{base}_n", count(b["n"]), b["n"], src + ".n")
            for st in ("median", "p90", "max"):
                F.add(f"{base}_{st}_pct", pctk(100 * b[st]), 100 * b[st], f"{src}.{st}")


COL_VARIANTS = {"frequency_resolved/flat": "fr", "frequency_resolved/spike-wave": "frsw", "band_variance": "bv",
                "frequency_resolved_squid_measured/flat": "frsq"}
COL_KEY = {vt: k for k, vt in COL_VARIANTS.items()}
COL_DESC = {"fr": "frequency-resolved whitening, flat signal spectrum", "frsw": "frequency-resolved, the G4 spike-wave spectrum",
            "bv": "band variance (spatial whitening only)", "frsq": "frequency-resolved, flat signal, measured Neuromag sensor spectrum"}


def sb_tok(lo, hi) -> str:
    return f"sb{tag(lo)}_up" if hi > 100 else f"sb{tag(lo)}_{tag(hi)}"


def ns_col_facts(F, d):
    c, cfg, p = d["coloured"], d["config"]["sensitivity"]["coloured_noise"], f"{NS} :: coloured"
    if c["corners_hz"] != cfg["corner_hz"]:
        raise ValueError(f"{NS}: coloured.corners_hz differ from the declared corner_hz")
    F.add("rev_ns_col_corners_hz", listing(c["corners_hz"]), c["corners_hz"], f"{p}.corners_hz (1/f corners of the OPM noise, Hz)")
    F.add("rev_ns_col_exponent", g(cfg["exponent"]), cfg["exponent"], f"{NS} :: config.sensitivity.coloured_noise.exponent "
          "(OPM noise PSD w^2 (1 + (f_c / f)^exponent))")
    F.add("rev_ns_col_n_subbands", count(len(c["subbands"])), len(c["subbands"]), f"{p}.subbands (sub-bands of the 1-40 Hz filter)")
    pc = c["partition_check"]
    F.add("rev_ns_col_partition_n_subbands", count(pc["n_subbands"]), pc["n_subbands"],
          f"{p}.partition_check.n_subbands (every sub-band split in {pc['split']})")
    mx = max(x["abs_diff_log2"] for r in pc["results"].values() for x in r.values())
    F.add("rev_ns_col_partition_max_change_pct", pct(100 * (2 ** mx - 1)), 100 * (2 ** mx - 1),
          f"{NS} :: derived: largest coloured.partition_check.results[*][*].abs_diff_log2 as % (finer partition vs the primary)")
    F.ratio("rev_ns_col_brain_grad_var_sum", c["checks"]["brain_grad_variance_sum_over_primary"],
            f"{p}.checks.brain_grad_variance_sum_over_primary (sub-band gradiometer brain variances summed, over the band's)", three=True)
    for corner, fct in c["band_variance_factor"].items():
        if float(corner) == 0:
            continue
        ct = f"c{corner}hz"
        F.ratio(f"rev_ns_col_bv_factor_{ct}", fct, f"{p}.band_variance_factor['{corner}'] (band-averaged OPM noise variance over "
                f"the white level's, {corner}-Hz corner)")
        F.add(f"rev_ns_col_bv_factor_{ct}_pct", pct(100 * (fct - 1)), 100 * (fct - 1), f"{NS} :: derived: coloured."
              f"band_variance_factor['{corner}'] - 1 (the band-averaged OPM noise variance increase of the {corner}-Hz corner, %)")
    # sub-bands
    meas = c["squid_measured_rms"]
    share = {"flat": {}, "spike-wave": {}}
    for i, sb in enumerate(c["subbands"]):
        st, src = sb_tok(sb["lo_hz"], sb["hi_hz"]), f"{p}.subbands[{i}]"
        hi = f"{sb['hi_hz']:g}" if sb["hi_hz"] <= 100 else f"{sb['hi_hz']:.0f}"
        F.add(f"rev_ns_col_{st}_label", f"{g(sb['lo_hz'])}–{hi}", [sb["lo_hz"], sb["hi_hz"]], src + ".lo_hz, .hi_hz (Hz)")
        for spec, tk in (("flat", "flat"), ("spike-wave", "spikewave")):
            v = 100 * sb["signal_share"][spec]
            F.add(f"rev_ns_col_{st}_signal_{tk}_pct", pct(v), v, f"{src}.signal_share['{spec}'] (signal energy share)")
            share[spec][(sb["lo_hz"], sb["hi_hz"])] = v
        for corner, fct in sb["opm_variance_factor"].items():
            if float(corner) != 0:
                F.ratio(f"rev_ns_col_{st}_opm_var_factor_c{corner}hz", fct, f"{src}.opm_variance_factor['{corner}'] (OPM noise "
                        "variance over the white level's in this sub-band)")
        F.add(f"rev_ns_col_{st}_brain_mag_ft", num(sb["brain_rms_mag_fT"], 1), sb["brain_rms_mag_fT"], f"{src}.brain_rms_mag_fT")
        F.add(f"rev_ns_col_{st}_brain_grad_ftcm", num(sb["brain_rms_grad_fT_cm"], 1), sb["brain_rms_grad_fT_cm"],
              f"{src}.brain_rms_grad_fT_cm")
        for kind in ("mag", "grad"):
            v = 100 * sb["env_explained"][kind]
            F.add(f"rev_ns_col_{st}_room_explained_{kind}_pct", pct(v), v, f"{src}.env_explained.{kind} (room-field share)")
        q = meas[i]
        if q["lo_hz"] != sb["lo_hz"]:
            raise ValueError(f"{NS}: coloured.squid_measured_rms[{i}] is not sub-band {i}")
        src = f"{p}.squid_measured_rms[{i}] (median Neuromag magnetometer sensor-noise RMS in the sub-band, fT)"
        F.add(f"rev_ns_col_{st}_squid_mag_measured_ft", num(q["mag_fT"], 1), q["mag_fT"], src + ".mag_fT (empty room)")
        F.add(f"rev_ns_col_{st}_squid_mag_brochure_ft", num(q["brochure_mag_fT"], 1), q["brochure_mag_fT"], src + ".brochure_mag_fT")
        F.ratio(f"rev_ns_col_{st}_squid_measured_over_brochure", q["mag_fT"] / q["brochure_mag_fT"],
                f"{NS} :: derived: coloured.squid_measured_rms[{i}].mag_fT / .brochure_mag_fT")
    for spec, tk in (("flat", "flat"), ("spike-wave", "spikewave")):
        for lo_, hi_, nm in ((0.0, 4.0, "below4hz"), (8.0, 24.0, "8_24hz")):
            v = sum(x for (lo, hi), x in share[spec].items() if lo >= lo_ and hi <= hi_)
            F.add(f"rev_ns_col_signal_{nm}_{tk}_pct", pct(v), v, f"{NS} :: derived: sum of coloured.subbands[*].signal_share['{spec}'] "
                  f"over the sub-bands within {lo_:g}-{hi_:g} Hz")
    # detectability ratios
    res, groups = c["results"], {}
    published = d["sweep"]["entries"]
    for rk, r in res.items():
        parts = rk.split("/")
        vt, lvl, corner = COL_VARIANTS["/".join(parts[:-2])], parts[-2].removesuffix("fT"), parts[-1].removesuffix("Hz")
        for key, e in r["comparisons"].items():
            tok = ctok(key)
            base = f"rev_ns_col_{vt}_asd{lvl}ft_c{corner}hz_{tok}"
            ns_cmp(F, d, base, e, f"coloured.results['{rk}'].comparisons['{key}']", key, note=COL_DESC[vt])
            groups.setdefault((vt, key), {}).setdefault(lvl, {})[corner] = e["ratio"]
        for key, bins in (r.get("depth") or {}).items():
            depth_summary(F, f"rev_ns_col_{vt}_asd{lvl}ft_c{corner}hz_{ctok(key)}", bins, f"coloured.results['{rk}'].depth['{key}']")
    for (vt, key), by_lvl in groups.items():  # every comparator at the primary level, the combined one at every level
        allv, vk, tok = [], COL_KEY[vt], ctok(key)
        for lvl, by_corner in by_lvl.items():
            white = by_corner["0"]
            for corner, rc in by_corner.items():
                if corner != "0":
                    ch = 100 * (rc / white - 1)
                    F.add(f"rev_ns_col_{vt}_asd{lvl}ft_c{corner}hz_{tok}_vs_white_pct", spct(ch), ch,
                          f"{NS} :: derived: coloured.results['{vk}/{lvl}fT/{corner}Hz'] over coloured.results['{vk}/{lvl}fT/0Hz'], "
                          f"comparisons['{key}'].ratio (% change by the {corner}-Hz corner)")
            F.span(f"rev_ns_col_{vt}_asd{lvl}ft_{tok}", list(by_corner.values()), r2, f"{p}.results['{vk}/{lvl}fT/*Hz']."
                   f"comparisons['{key}'].ratio (corners {listing(sorted(by_corner, key=float), str)} Hz)", extra=(("3dp", r3),))
            allv += list(by_corner.values())
            if vt == "fr":
                ch = 100 * (white / published[lvl]["comparisons"][key]["ratio"] - 1)
                F.add(f"rev_ns_col_fr_asd{lvl}ft_{tok}_white_vs_published_pct", spct(ch), ch,
                      f"{NS} :: derived: coloured.results['{vk}/{lvl}fT/0Hz'] over sweep.entries['{lvl}'].comparisons['{key}'], "
                      "ratio (the frequency-resolved detector with white noise against the band detector, % change)")
        if len(by_lvl) > 1:
            F.span(f"rev_ns_col_{vt}_{tok}", allv, r2, f"{p}.results['{vk}/*fT/*Hz'].comparisons['{key}'].ratio (every white level "
                   "and corner)", extra=(("3dp", r3),))


FF_SOURCES = {"cardiac/current_dipole/200mm": "heart_cd200", "cardiac/current_dipole/250mm": "heart_cd250",
              "cardiac/magnetic_dipole/250mm": "heart_md250", "cardiac/current_dipole/300mm": "heart_cd300",
              "ocular/head_bem": "eye_bem", "ocular/current_dipole": "eye_cd", "cardiac+ocular": "heart_eye"}
FF_LEVELS = {"room_field": "room", "magnetometer_shortfall": "short"}
FF_MODES = {"added": "add", "recalibrated": "recal"}
FF_SETS = {"squid/combined": "combined", "squid/mag": "mag", "squid/grad": "grad", "opm_matched/all": "matched",
           "opm_dense/all": "dense"}


def ns_ff_facts(F, d):
    f, cfg, p = d["far_field"], d["config"]["sensitivity"]["far_field"], f"{NS} :: far_field"
    cp = f"{NS} :: config.sensitivity.far_field"
    F.add("rev_ns_ff_heart_depths_mm", listing(cfg["heart_depth_mm"]), cfg["heart_depth_mm"], f"{cp}.heart_depth_mm (below the "
          "head origin)")
    F.add("rev_ns_ff_heart_x_mm", g(cfg["heart_head_mm"][0]), cfg["heart_head_mm"][0], f"{cp}.heart_head_mm[0] (head frame, left -)")
    F.add("rev_ns_ff_heart_y_mm", g(cfg["heart_head_mm"][1]), cfg["heart_head_mm"][1], f"{cp}.heart_head_mm[1] (anterior +)")
    F.add("rev_ns_ff_combined_share_pct", pct(100 * cfg["combined_share"]), 100 * cfg["combined_share"],
          f"{cp}.combined_share (each source's share of the level's variance in the combined heart + eyes term)")
    F.add("rev_ns_ff_eye_distance_mm", mm(f["eyes"]["interocular_distance_mm"]), f["eyes"]["interocular_distance_mm"],
          f"{p}.eyes.interocular_distance_mm (eye centres located in the sample T1)")
    mx = max(f["coil_integration_check_max_rel"].values())
    F.add("rev_ns_ff_coil_check_max_rel", small(mx), mx, f"{p}.coil_integration_check_max_rel (largest over the arrays)")
    lv = f["level_values_fT"]
    for k, nm in (("room_field", "room"), ("magnetometer_shortfall", "short"), ("measured_brain_mag", "measured_brain_mag"),
                  ("modelled_brain_mag", "modelled_brain_mag")):
        F.add(f"rev_ns_ff_level_{nm}_ft", num(lv[k], 0), lv[k], f"{p}.level_values_fT.{k} (median RMS at the Neuromag "
              "magnetometers, fT)")
    if abs(lv["magnetometer_shortfall"] ** 2 - (lv["measured_brain_mag"] ** 2 - lv["modelled_brain_mag"] ** 2)) > 1e-6 * lv["measured_brain_mag"] ** 2:
        raise ValueError(f"{NS}: far_field.level_values_fT.magnetometer_shortfall is not sqrt(measured^2 - modelled^2)")
    for sk, st in FF_SOURCES.items():
        if sk not in f["projection_energy_left"]:
            continue
        for setk, x in f["projection_energy_left"][sk].items():
            base = f"rev_ns_ff_left_{st}_{FF_SETS[setk]}"
            src = (f"{p}.projection_energy_left['{sk}']['{setk}'] (share of the source's field energy left after the 8-term "
                   "projection)")
            F.add(f"{base}_pct", pct(100 * x["total"]), 100 * x["total"], src + ".total")
            for ax, v in x["per_axis"].items():
                F.add(f"{base}_{ax}_pct", pct(100 * v), 100 * v, f"{src}.per_axis.{ax} (moment along {ax})")
    for arr, x in f["projection_context"].items():
        base, src = f"rev_ns_ff_context_{ARR3[arr]}", f"{p}.projection_context.{arr} (field energy share left after the projection)"
        F.add(f"{base}_targets_median_pct", pct(100 * x["cortical_targets_median"]), 100 * x["cortical_targets_median"],
              src + ".cortical_targets_median")
        lo, hi = x["cortical_targets_p5_p95"]
        F.add(f"{base}_targets_p5_p95_range", f"{pct(100 * lo)} to {pct(100 * hi)}", [100 * lo, 100 * hi], src + ".cortical_targets_p5_p95")
        F.add(f"{base}_background_pct", pct(100 * x["cortical_background"]), 100 * x["cortical_background"], src + ".cortical_background")
        F.add(f"{base}_room_fraction", small(x["room_field"]), x["room_field"], src + ".room_field (fraction)")
    prim = f"{d['config']['adult']['sensors']['opm_asd_primary_fT_per_rtHz']:g}"
    pub = d["sweep"]["entries"][prim]["comparisons"]
    ratios, infeasible = {}, []
    for rk, r in f["results"].items():
        sk, lk, mk = rk.rsplit("/", 2)
        base, cal, src = f"rev_ns_ff_{FF_SOURCES[sk]}_{FF_LEVELS[lk]}_{FF_MODES[mk]}", r["calibration"], f"{p}.results['{rk}']"
        F.add(f"{base}_feasible", yes(cal["feasible"]), cal["feasible"], src + ".calibration.feasible" +
              (f" ({cal['reason']})" if not cal["feasible"] else ""))
        F.ratio(f"{base}_grad_model_over_measured", cal["grad_rms_model_over_measured"], f"{src}.calibration."
                "grad_rms_model_over_measured (median gradiometer brain-noise RMS, model over measured)", three=True)
        F.ratio(f"{base}_mag_model_over_measured", cal["mag_rms_model_over_measured"], f"{src}.calibration."
                "mag_rms_model_over_measured (median magnetometer brain-noise RMS, model over measured)", three=True)
        if not cal["feasible"]:
            infeasible.append(rk)
            continue
        F.ratio(f"{base}_bg_scale", cal["scale_over_gradiometer_calibration"], f"{src}.calibration.scale_over_gradiometer_calibration "
                "(cortical background scale over the published one)", three=True)
        for k, nm, nd in (("squid_mag_fT", "mag_ft", 0), ("squid_grad_fT_cm", "grad_ftcm", 1), ("opm_matched_fT", "matched_ft", 0),
                          ("opm_dense_fT", "dense_ft", 0)):
            F.add(f"{base}_rms_{nm}", num(r["rms"][k], nd), r["rms"][k], f"{src}.rms.{k} (far-field median in-band RMS)")
        for key, e in r["comparisons"].items():
            ns_cmp(F, d, f"{base}_{ctok(key)}", e, f"far_field.results['{rk}'].comparisons['{key}']", key)
            ratios.setdefault(key, []).append(e["ratio"])
        for key, bins in (r.get("depth") or {}).items():
            depth_summary(F, f"{base}_{ctok(key)}", bins, f"far_field.results['{rk}'].depth['{key}']")
    F.add("rev_ns_ff_n_results", count(len(f["results"])), len(f["results"]), f"{p}.results (sources x levels x calibrations)")
    F.add("rev_ns_ff_n_infeasible", count(len(infeasible)), len(infeasible), f"{p}.results[*].calibration.feasible (false: "
          f"{listing(infeasible, str)})")
    for key, vals in ratios.items():
        tok = ctok(key)
        F.span(f"rev_ns_ff_{tok}", vals, r2, f"{p}.results[*].comparisons['{key}'].ratio (every feasible source, level and "
               "calibration)", extra=(("3dp", r3),))
        ch = [100 * (v / pub[key]["ratio"] - 1) for v in vals]
        F.span(f"rev_ns_ff_{tok}_change_pct", ch, spct, f"{NS} :: derived: far_field.results[*].comparisons['{key}'].ratio over "
               f"sweep.entries['{prim}'].comparisons['{key}'].ratio (the published model), % change (every feasible far-field result)")


JOINT = {"baseline_white": "white", "joint": "all", "joint_room_field_level": "allroom", "joint_without_near_skull": "nonear",
         "joint_without_far_field": "nofar", "joint_without_coloured": "nocol"}


def ns_joint_facts(F, d):
    j, cfg, p = d["joint"], d["config"]["sensitivity"]["joint"], f"{NS} :: joint"
    cp = f"{NS} :: config.sensitivity.joint"
    F.add("rev_ns_joint_floor_mm", g(cfg["near_skull_floor_mm"]), cfg["near_skull_floor_mm"], f"{cp}.near_skull_floor_mm")
    F.add("rev_ns_joint_corner_hz", g(cfg["corner_hz"]), cfg["corner_hz"], f"{cp}.corner_hz")
    by_lvl = {}
    for rk, r in j["results"].items():
        vk, lv = rk.rsplit("/", 1)
        vt, lvl = JOINT[vk], lv.removesuffix("fT")
        for key, e in r["comparisons"].items():
            if vt in ("all", "nofar") or "/combined/" in key:
                ns_cmp(F, d, f"rev_ns_joint_{vt}_asd{lvl}ft_{ctok(key)}", e, f"joint.results['{rk}'].comparisons['{key}']", key)
            if vt != "white" and "/combined/" in key:
                by_lvl.setdefault((lvl, key), []).append((e["ratio"], e["ci95_ratio"], rk))
        for key, bins in (r.get("depth") or {}).items():
            depth_summary(F, f"rev_ns_joint_{vt}_asd{lvl}ft_{ctok(key)}", bins, f"joint.results['{rk}'].depth['{key}']")
    ff = {}
    for rk, r in j["results"].items():
        if r.get("far_field_rms_mag_fT") is not None:
            vt = JOINT[rk.rsplit("/", 1)[0]]
            if ff.setdefault(vt, (r["far_field_rms_mag_fT"], rk))[0] != r["far_field_rms_mag_fT"]:
                raise ValueError(f"{NS}: joint far-field RMS of {vt} differs between white levels")
    for vt, (vals, rk) in ff.items():
        F.span(f"rev_ns_joint_{vt}_ff_rms_mag_ft", vals, lambda x: num(x, 0), f"{p}.results['{rk}'].far_field_rms_mag_fT (far-field "
               "median RMS at the magnetometers per sub-band, fT; the same at every white level)")
    for (lvl, key), rows in by_lvl.items():
        r, ci, rk = min(rows)
        base = f"rev_ns_joint_worst_asd{lvl}ft_{ctok(key)}"
        F.ratio(base, r, f"{NS} :: derived: the lowest joint.results['*/{lvl}fT'].comparisons['{key}'] (ratio, ci95_ratio) over the "
                f"joint runs at this white level: {rk}", ci, three=True)
        white = j["results"][f"baseline_white/{lvl}fT"]["comparisons"][key]["ratio"]
        ch = 100 * (r / white - 1)
        F.add(f"{base}_vs_white_pct", spct(ch), ch, f"{NS} :: derived: {rk} over joint.results['baseline_white/{lvl}fT'], "
              f"comparisons['{key}'] ratio (% change against the frequency-resolved white-noise reference)")
    for vk, vals in j["gradiometer_rms_model_over_measured"].items():
        F.span(f"rev_ns_joint_{JOINT[vk]}_grad_model_over_measured", vals, r2, f"{p}.gradiometer_rms_model_over_measured['{vk}'] "
               "(median gradiometer brain-noise RMS, model over measured, per sub-band)")


def headline_level(label: str, primary: float) -> float:
    """The white level of a headline row: the one it names, else the primary level."""
    m = re.search(r"(\d+(?:\.\d+)?) fT/sqrt\(Hz\)", label)
    return float(m.group(1)) if m else primary


def ns_alt_facts(F, d):
    """The range over the run's own headline table (every model alternative), per white level and over all."""
    h, prim = d["headline"], d["config"]["adult"]["sensors"]["opm_asd_primary_fT_per_rtHz"]
    # a row that names no white level must be a primary-level result of the near-skull, coloured or far-field runs
    known = []
    for fk in [k for k in d["near_skull"]["primary"] if k.startswith("floor_")]:
        f = d["near_skull"]["primary"][fk]
        known += [f["background"]["comparisons"], f["targets"]["all"], f["both"]["all"]]
    known += [r["comparisons"] for rk, r in d["coloured"]["results"].items() if f"/{prim:g}fT/" in rk]
    known += [r["comparisons"] for r in d["far_field"]["results"].values() if r.get("comparisons")]
    for label, row in h.items():
        if "fT/sqrt(Hz)" not in label and not any(all(c.get(k, {}).get("ratio") == v["ratio"] for k, v in row.items()) for c in known):
            raise ValueError(f"{NS}: headline row {label!r} names no white level and matches no primary-level result")
    F.add("rev_ns_alt_n_rows", count(len(h)), len(h), f"{NS} :: headline (rows: model alternatives x white levels)")
    keys = list(next(iter(h.values())).keys())
    levels = sorted({headline_level(lb, prim) for lb in h})
    for lvl in levels + [None]:
        rows = {lb: r for lb, r in h.items() if lvl is None or headline_level(lb, prim) == lvl}
        lt = "all" if lvl is None else f"asd{lvl:g}ft"
        where = "every row" if lvl is None else (f"the rows at {lvl:g} fT/sqrt(Hz) (rows naming no level are at the primary "
                                                 f"{prim:g})" if lvl == prim else f"the rows at {lvl:g} fT/sqrt(Hz)")
        F.add(f"rev_ns_alt_{lt}_n_rows", count(len(rows)), len(rows), f"{NS} :: headline, {where}")
        for key in keys:
            tok, vals = ctok(key), {lb: r[key] for lb, r in rows.items()}
            src = f"{NS} :: headline[*]['{key}'] ({where})"
            F.span(f"rev_ns_alt_{lt}_{tok}", [v["ratio"] for v in vals.values()], r2, src + ".ratio", extra=(("3dp", r3),))
            above = sum(v["ci95_ratio"][0] > 1 for v in vals.values())
            below = sum(v["ci95_ratio"][1] < 1 for v in vals.values())
            F.add(f"rev_ns_alt_{lt}_{tok}_n_ci_above1", count(above), above, src + ".ci95_ratio (rows whose 95 % CI lies above 1)")
            F.add(f"rev_ns_alt_{lt}_{tok}_n_ci_below1", count(below), below, src + ".ci95_ratio (rows whose 95 % CI lies below 1)")
            lo = min(vals, key=lambda lb: vals[lb]["ratio"])
            hi = max(vals, key=lambda lb: vals[lb]["ratio"])
            F.add(f"rev_ns_alt_{lt}_{tok}_min_row", lo, lo, src + " (the row with the lowest ratio)")
            F.add(f"rev_ns_alt_{lt}_{tok}_max_row", hi, hi, src + " (the row with the highest ratio)")


def ns_facts(F, root):
    d = load(root, NS)
    ns_setup_facts(F, root, d)
    ns_sweep_facts(F, d)
    ns_near_facts(F, d)
    ns_col_facts(F, d)
    ns_ff_facts(F, d)
    ns_joint_facts(F, d)
    ns_alt_facts(F, d)


# ------------------------------------------------------------------------------------------------
# noise model against the measured Neuromag covariance (rev_cov_)
COV_COMPS = {"brain": "brain", "total": "total", "empty_room": "emptyroom", "brain_without_heart": "noheart"}
SPLITS = {"fit_left_predict_right": "left_to_right", "fit_right_predict_left": "right_to_left",
          "fit_posterior_predict_anterior": "post_to_ant", "fit_anterior_predict_posterior": "ant_to_post",
          "fit_inferior_predict_superior": "inf_to_sup", "fit_superior_predict_inferior": "sup_to_inf"}
DET = {"measured_over_model": "measured", "measured_corrected_over_model": "corrected", "hybrid_over_model": "hybrid",
       "hybrid_corrected_over_model": "hybrid_corrected", "measured_corrected_over_model_intrinsic_brain": "corrected_ib",
       "half1_over_model": "half1", "half2_over_model": "half2",
       "time_locked_removed_corrected_over_model": "timelocked_corrected", "without_heart_corrected_over_model": "noheart_corrected",
       "measured_corrected_over_model_independent_mag_calibrated": "corrected_vs_magcal",
       "measured_corrected_over_model_correlated_5mm": "corrected_vs_corr5mm",
       "measured_corrected_over_model_correlated_10mm": "corrected_vs_corr10mm", "model305_over_model306": "model305_vs_306"}
SCEN = {"model_306": "model306", "model_306/intrinsic+brain": "model306_ib", "model_305": "model305",
        "S1_same_relative_change": "s1", "S2_same_variance_excess_both_modelled": "s2",
        "S3_neuromag_measured_opm_variance_excess": "s3", "S4_neuromag_measured_opm_as_modelled": "s4",
        "S4_uncorrected": "s4_uncorrected", "S5_neuromag_measured_empty_room": "s5",
        "variant_correlated_5mm/model": "corr5mm_model", "variant_correlated_5mm/neuromag_measured": "corr5mm_measured",
        "variant_correlated_10mm/model": "corr10mm_model", "variant_correlated_10mm/neuromag_measured": "corr10mm_measured",
        "projected/model_305": "proj_model305", "projected/S4_uncorrected": "proj_s4_uncorrected"}
VARIANTS = {"independent_mag_calibrated": "magcal", "correlated_5mm": "corr5mm", "correlated_10mm": "corr10mm"}


def dist_tok(label: str) -> str:
    return "same_site" if label == "same site" else "d" + label.replace("-", "_")


def cov_cmp(F, base, e, path, three=True, note=""):
    """A stored detectability ratio of the covariance validation (ratio, ci95 on the ratio scale, share above 1)."""
    src = (f"{COV} :: {path} (median over {e['n']:,} targets{'; ' + note if note else ''}; 95 % CI from the parcel bootstrap "
           f"over {parcels_of(e['ci_method'])} parcel labels)")
    F.ratio(base, e["ratio"], src, e["ci95"], three)
    F.add(f"{base}_share", pct(100 * e["share_above_1"]), 100 * e["share_above_1"], f"{COV} :: {path}.share_above_1 (share of "
          "the targets with a ratio above 1)")


def cov_setup_facts(F, c):
    dc, da, p = c["declared_choices"], c["data"], f"{COV} :: declared_choices"
    F.add("rev_cov_n_windows", count(da["n_windows"]), da["n_windows"], f"{COV} :: data.n_windows (task baseline windows)")
    F.add("rev_cov_window_s", f"{g(da['window_s'][0])} to {g(da['window_s'][1])}", da["window_s"], f"{COV} :: data.window_s "
          "(baseline window relative to the stimulus, s)")
    F.add("rev_cov_n_baseline_samples", count(da["n_baseline_samples"]), da["n_baseline_samples"], f"{COV} :: data.n_baseline_samples")
    F.add("rev_cov_n_empty_room_samples", count(da["n_empty_room_samples"]), da["n_empty_room_samples"],
          f"{COV} :: data.n_empty_room_samples")
    F.add("rev_cov_band_hz", f"{g(da['band_hz'][0])} to {g(da['band_hz'][1])}", da["band_hz"], f"{COV} :: data.band_hz")
    F.add("rev_cov_bad_channels", listing(da["bad_channels"], str), da["bad_channels"], f"{COV} :: data.bad_channels")
    F.add("rev_cov_subbands_hz", listing(dc["sub_bands_hz"], lambda b: f"{g(b[0])}–{g(b[1])}"), dc["sub_bands_hz"],
          f"{p}.sub_bands_hz ({dc['sub_band_note']})")
    bins = dc["distance_bins_mm"]
    steps = {b - a for a, b in zip(bins[:-1], bins[1:])}
    if len(steps) != 1:
        raise ValueError(f"{COV}: declared_choices.distance_bins_mm are not equally spaced")
    F.add("rev_cov_distance_bins_range", f"{g(bins[0])} to {g(bins[-1])}", [bins[0], bins[-1]], f"{p}.distance_bins_mm (mm)")
    F.add("rev_cov_distance_bin_mm", g(steps.pop()), bins[1] - bins[0], f"{p}.distance_bins_mm (bin width, mm)")
    F.add("rev_cov_subspace_dims", listing(dc["subspace_dimensions"]), dc["subspace_dimensions"], f"{p}.subspace_dimensions")
    for k, nm in (("n_surrogates", "n_surrogates"), ("n_band_surrogates", "n_band_surrogates"), ("n_random_splits", "n_random_splits"),
                  ("n_parcel_bootstrap", "n_boot")):
        F.add(f"rev_cov_{nm}", count(dc[k]), dc[k], f"{p}.{k}")
    F.add("rev_cov_seed", str(dc["seed"]), dc["seed"], f"{p}.seed")
    F.add("rev_cov_commit", c["provenance"]["commit"], c["provenance"]["commit"], f"{COV} :: provenance.commit (code of the run)")
    ck = c["checks"]
    F.add("rev_cov_check_brain_scale_dev", g(ck["brain_scale_rel_dev_from_g2"]), ck["brain_scale_rel_dev_from_g2"],
          f"{COV} :: checks.brain_scale_rel_dev_from_g2 (relative deviation of the cortical scale from G2's)")
    mx = max(ck["g2_detectability_max_abs_dev_vs_stored_csv"].values())
    F.add("rev_cov_check_g2_detect_max_dev", small(mx), mx, f"{COV} :: checks.g2_detectability_max_abs_dev_vs_stored_csv (largest; "
          "the CSV stores 4 decimals)")
    mx = max(ck["diag_baseline_vs_g2_max_rel_dev"], ck["diag_empty_room_vs_g2_max_rel_dev"])
    F.add("rev_cov_check_diag_max_rel_dev", small(mx), mx, f"{COV} :: checks.diag_baseline_vs_g2_max_rel_dev, "
          "diag_empty_room_vs_g2_max_rel_dev (channel variances against G2's, largest relative deviation)")
    # what is fitted and what is predicted
    rows = c["fitted_vs_predicted"]
    fitted = [r["quantity"] for r in rows if r["role"] == "fitted"]
    by_construction = [r["quantity"] for r in rows if r["role"].startswith("fitted (")]
    predicted = [r["quantity"] for r in rows if r["role"].startswith("independent prediction")]
    not_fitted = [r["quantity"] for r in rows if r["role"].startswith("not fitted")]
    if len(fitted) + len(by_construction) + len(predicted) + len(not_fitted) != len(rows):
        raise ValueError(f"{COV}: a fitted_vs_predicted role is not recognised")
    src = f"{COV} :: fitted_vs_predicted[*].role"
    F.add("rev_cov_n_quantities", count(len(rows)), len(rows), f"{COV} :: fitted_vs_predicted (rows)")
    F.add("rev_cov_n_fitted", count(len(fitted)), len(fitted), src + " ('fitted')")
    F.add("rev_cov_n_fitted_by_construction", count(len(by_construction)), len(by_construction), src + " ('fitted (equal by construction)')")
    F.add("rev_cov_n_predicted", count(len(predicted)), len(predicted), src + " ('independent prediction ...')")
    F.add("rev_cov_n_not_fitted", count(len(not_fitted)), len(not_fitted), src + " ('not fitted ...': the brochure sensor noise)")
    F.add("rev_cov_fitted_quantities", listing(fitted, str), fitted, f"{COV} :: fitted_vs_predicted[role == 'fitted'].quantity")
    F.add("rev_cov_predicted_quantities", "; ".join(predicted), predicted,
          f"{COV} :: fitted_vs_predicted[role == 'independent prediction ...'].quantity")
    m = c["model"]
    F.ratio("rev_cov_k", m["magnetometer_calibration_scale_ratio"], f"{COV} :: model.magnetometer_calibration_scale_ratio (k: the "
            "background scale that matches the magnetometers over the gradiometer-calibrated one; = opm_implication.checks."
            "magnetometer_variance_excess_k)", three=True)
    if m["magnetometer_calibration_scale_ratio"] != c["opm_implication"]["checks"]["magnetometer_variance_excess_k"]:
        raise ValueError(f"{COV}: the two stored k differ")
    F.add("rev_cov_enbw_hz", num(m["enbw_hz"], 1), m["enbw_hz"], f"{COV} :: model.enbw_hz")
    ind = m["variants"]["independent"]["scale"]
    for vk, vt in VARIANTS.items():
        v = m["variants"][vk]
        F.ratio(f"rev_cov_variant_{vt}_scale", v["scale"] / ind, f"{COV} :: derived: model.variants['{vk}'].scale / "
                f"model.variants['independent'].scale ({v['label']})")
    for kind, v in m["environment_explained_fraction"].items():
        F.add(f"rev_cov_room_explained_{kind}_pct", pct(100 * v), 100 * v, f"{COV} :: model.environment_explained_fraction.{kind} "
              "(empty-room variance explained by the fitted room field)")


def cov_variance_facts(F, c):
    for vk, v in c["variances"].items():
        comp, kind = vk.split("/")
        base, p = f"rev_cov_{COV_COMPS[comp]}_{kind}", f"{COV} :: variances['{vk}']"
        F.ratio(f"{base}_amp", v["amp_ratio_of_medians_model_over_measured"], f"{p}.amp_ratio_of_medians_model_over_measured "
                "(median channel RMS, model over measured)", three=True)
        lo, _, hi = v["null95"]
        src = f"{p}.null95 (2.5th and 97.5th percentiles of the same ratio for an exact model measured the same way)"
        F.add(f"{base}_amp_null_ci", ci_text(lo, hi), [lo, hi], src)
        F.add(f"{base}_amp_null_ci_3dp", ci_text(lo, hi, r3), [lo, hi], src)
        F.add(f"{base}_amp_null_range", f"{r2(lo)} to {r2(hi)}", [lo, hi], src)
        q = v["var_ratio_model_over_measured"]
        for st, x in q.items():
            F.add(f"{base}_var_ratio_{st}", r2(x), x, f"{p}.var_ratio_model_over_measured.{st} (per-channel variance, model over measured)")
        F.add(f"{base}_var_ratio_p5_p95_range", f"{r2(q['p5'])} to {r2(q['p95'])}", [q["p5"], q["p95"]],
              f"{p}.var_ratio_model_over_measured.p5, .p95")
        a, b = v["null_var_ratio_p5"][1], v["null_var_ratio_p95"][1]
        F.add(f"{base}_var_ratio_null_p5_p95_range", f"{r2(a)} to {r2(b)}", [a, b], f"{p}.null_var_ratio_p5[1], "
              "null_var_ratio_p95[1] (medians over the exact-model surrogates of the 5th and 95th percentiles)")
        F.add(f"{base}_log_var_r", r2(v["log_var_pearson"]), v["log_var_pearson"], f"{p}.log_var_pearson (Pearson r of the log "
              "channel variances, model vs measured)")
        F.add(f"{base}_log_var_rho", r2(v["log_var_spearman"]), v["log_var_spearman"], f"{p}.log_var_spearman")
        lo, _, hi = v["null_log_var_pearson"]
        F.add(f"{base}_log_var_r_null_ci", ci_text(lo, hi), [lo, hi], f"{p}.null_log_var_pearson (exact model, 2.5th-97.5th percentiles)")
        F.add(f"{base}_log_var_r_split_half", r2(v["split_half_log_var_pearson"]), v["split_half_log_var_pearson"],
              f"{p}.split_half_log_var_pearson (the measurement's own reproducibility: one half of the windows vs the other)")
        F.add(f"{base}_n_nonpositive", count(v["n_nonpositive_measured"]), v["n_nonpositive_measured"],
              f"{p}.n_nonpositive_measured (channels whose measured variance is not positive)")
    # per channel: where the model is below the measurement
    ch = c["channels"]
    for kind in ("mag", "grad"):
        idx = [i for i, k in enumerate(ch["kind"]) if k == kind and ch["measured_brain_var"][i] > 0]
        below = sum(ch["model_brain_var"][i] < ch["measured_brain_var"][i] for i in idx)
        src = (f"{COV} :: derived: channels.model_brain_var < channels.measured_brain_var over the good {kind} channels with a "
               "positive measured brain variance")
        F.add(f"rev_cov_brain_{kind}_n_model_below", count(below), below, src)
        F.add(f"rev_cov_brain_{kind}_n_channels_positive", count(len(idx)), len(idx), src + " (their number)")
    h = c["held_out"]
    rs, p = h["random_site_splits"], f"{COV} :: held_out.random_site_splits"
    F.add("rev_cov_heldout_n_splits", count(rs["n_splits"]), rs["n_splits"], f"{p}.n_splits (random halves of the gradiometer sites)")
    F.add("rev_cov_heldout_n_sites", count(rs["n_sites"]), rs["n_sites"], f"{p}.n_sites (gradiometer sites)")
    for k, nm in (("heldout_amp_ratio_of_medians", "grad_amp"), ("mag_amp_ratio_of_medians", "mag_amp")):
        q = rs[k]
        F.ratio(f"rev_cov_heldout_random_{nm}", q["p50"], f"{p}.{k} (median over the splits; interval: 2.5th-97.5th percentiles "
                "over the splits)", [q["p2.5"], q["p97.5"]], three=True)
    for k, nm, what in (("heldout_log_var_pearson", "log_var_r", "Pearson r of the held-out log channel variances, model vs "
                         "measured"), ("heldout_median_abs_log10_ratio", "abs_log10", "median |log10 model/measured| of the "
                                       "held-out channel variances")):
        q, src = rs[k], f"{p}.{k} ({what})"
        F.add(f"rev_cov_heldout_random_{nm}", r2(q["p50"]), q["p50"], src + ".p50 (median over the splits)")
        F.add(f"rev_cov_heldout_random_{nm}_ci", ci_text(q["p2.5"], q["p97.5"]), [q["p2.5"], q["p97.5"]],
              src + ".p2.5, .p97.5 (over the splits)")
    x = 10 ** rs["heldout_median_abs_log10_ratio"]["p50"]
    F.add("rev_cov_heldout_random_typical_factor", r2(x), x, f"{COV} :: derived: 10**held_out.random_site_splits."
          "heldout_median_abs_log10_ratio.p50 (the typical factor between a held-out channel's model and measured variance)")
    amps = []
    for sk, st in SPLITS.items():
        x, src = h["spatial_splits"][sk], f"{COV} :: held_out.spatial_splits.{sk} ({h['spatial_split_rule']})"
        F.ratio(f"rev_cov_heldout_{st}_grad_amp", x["heldout_amp_ratio_of_medians"], src + ".heldout_amp_ratio_of_medians", three=True)
        F.ratio(f"rev_cov_heldout_{st}_mag_amp", x["mag_amp_ratio_of_medians"], src + ".mag_amp_ratio_of_medians", three=True)
        F.add(f"rev_cov_heldout_{st}_log_var_r", r2(x["heldout_log_var_pearson"]), x["heldout_log_var_pearson"], src + ".heldout_log_var_pearson")
        F.add(f"rev_cov_heldout_{st}_abs_log10", r2(x["heldout_median_abs_log10_ratio"]), x["heldout_median_abs_log10_ratio"],
              src + ".heldout_median_abs_log10_ratio")
        F.add(f"rev_cov_heldout_{st}_n_cal", count(x["n_calibration"]), x["n_calibration"], src + ".n_calibration (gradiometers)")
        F.add(f"rev_cov_heldout_{st}_n_heldout", count(x["n_heldout"]), x["n_heldout"], src + ".n_heldout (gradiometers)")
        amps.append(x["heldout_amp_ratio_of_medians"])
    F.span("rev_cov_heldout_spatial_grad_amp", amps, r2, f"{COV} :: held_out.spatial_splits[*].heldout_amp_ratio_of_medians "
           "(six spatial halves)", extra=(("3dp", r3),))
    F.add("rev_cov_heldout_all_n_grad", count(h["all_gradiometers"]["n_calibration"]), h["all_gradiometers"]["n_calibration"],
          f"{COV} :: held_out.all_gradiometers.n_calibration (good gradiometers)")


def cov_heart_facts(F, c):
    ca, p = c["cardiac_check"], f"{COV} :: cardiac_check"
    F.add("rev_cov_heart_n_beats", count(ca["n_beats"]), ca["n_beats"], f"{p}.n_beats ({ca['method'][:70]}...)")
    F.add("rev_cov_heart_rr_s", num(ca["rr_s"], 2), ca["rr_s"], f"{p}.rr_s (median R-R interval, s)")
    F.add("rev_cov_heart_rate_bpm", num(ca["heart_rate_bpm"], 0), ca["heart_rate_bpm"], f"{p}.heart_rate_bpm")
    for kind, unit, nd in (("mag", "fT", 0), ("grad", "fT/cm", 1)):
        x, base = ca[kind], f"rev_cov_heart_{kind}"
        u = unit.replace("/", "").lower()
        for k, nm in (("median_rms_in_baseline", "rms"), ("max_rms_in_baseline", "rms_max"), ("median_rms_average_waveform", "rms_waveform")):
            F.add(f"{base}_{nm}_{u}", num(x[k], nd), x[k], f"{p}.{kind}.{k} ({unit})")
        F.add(f"{base}_share_brain_pct", pct(100 * x["median_share_of_measured_brain_variance"]),
              100 * x["median_share_of_measured_brain_variance"], f"{p}.{kind}.median_share_of_measured_brain_variance")
        if x["share_of_model_shortfall_at_the_median"] is not None:
            F.add(f"{base}_share_shortfall_pct", pct(100 * x["share_of_model_shortfall_at_the_median"]),
                  100 * x["share_of_model_shortfall_at_the_median"], f"{p}.{kind}.share_of_model_shortfall_at_the_median")
        F.ratio(f"{base}_amp_without", x["brain_amp_ratio_model_over_measured_without_heart"],
                f"{p}.{kind}.brain_amp_ratio_model_over_measured_without_heart (median RMS, model over measured, heart removed)", three=True)
        F.add(f"{base}_leading_power_pct", pct(100 * x["leading_component_power_share"]), 100 * x["leading_component_power_share"],
              f"{p}.{kind}.leading_component_power_share (of the cardiac field)")
        F.add(f"{base}_leading_external_pct", pct(100 * x["leading_component_external_share"]),
              100 * x["leading_component_external_share"], f"{p}.{kind}.leading_component_external_share (inside the 8-term subspace)")
        F.add(f"{base}_leading_cos2_excess", r2(x["leading_component_cos2_with_leading_excess"]),
              x["leading_component_cos2_with_leading_excess"], f"{p}.{kind}.leading_component_cos2_with_leading_excess")
    t = c["time_locked_check"]
    for kind in ("mag", "grad"):
        x, p = t[kind], f"{COV} :: time_locked_check.{kind}"
        F.ratio(f"rev_cov_timelocked_{kind}_brain_change", x["brain_amp_ratio_time_locked_removed_over_primary"],
                f"{p}.brain_amp_ratio_time_locked_removed_over_primary (measured brain RMS with the time-locked average removed, "
                "over without)", three=True)
        F.ratio(f"rev_cov_timelocked_{kind}_amp", x["model_over_measured_time_locked_removed"],
                f"{p}.model_over_measured_time_locked_removed", three=True)


def cov_structure_facts(F, c):
    s, null = c["structure"], c["structure_null"]
    for ck, ct in COV_COMPS.items():
        for kind in ("mag", "grad"):
            x, p = s[ck][kind], f"{COV} :: structure.{ck}.{kind}"
            F.add(f"rev_cov_{ct}_{kind}_corr_pattern_r", r2(x["corr_pattern_pearson"]), x["corr_pattern_pearson"],
                  f"{p}.corr_pattern_pearson (Pearson r of the pairwise channel correlations, model vs measured)")
            if ck == "brain_without_heart":
                F.ratio(f"rev_cov_noheart_{kind}_amp", x["amp_ratio_of_medians"], f"{p}.amp_ratio_of_medians (model over measured, "
                        "the cardiac-locked average removed)", three=True)
                F.add(f"rev_cov_noheart_{kind}_log_var_r", r2(x["log_var_pearson"]), x["log_var_pearson"], f"{p}.log_var_pearson")
                q = x["var_ratio"]
                F.add(f"rev_cov_noheart_{kind}_var_ratio_p5_p95_range", f"{r2(q['p5'])} to {r2(q['p95'])}", [q["p5"], q["p95"]],
                      f"{p}.var_ratio.p5, .p95")
    for kind in ("mag", "grad"):
        x, p = s["brain"][kind], f"{COV} :: structure.brain.{kind}"
        F.add(f"rev_cov_brain_{kind}_pattern_reliability", r2(x["measured_pattern_reliability"]), x["measured_pattern_reliability"],
              f"{p}.measured_pattern_reliability (Spearman-Brown of the split-half log-variance correlation)")
        F.add(f"rev_cov_brain_{kind}_log_var_r_disattenuated", r2(x["log_var_pearson_disattenuated"]),
              x["log_var_pearson_disattenuated"], f"{p}.log_var_pearson_disattenuated")
        for i, b in enumerate(x["corr_vs_distance"]):
            if "agreement" not in b:
                continue
            base, src = f"rev_cov_dist_{kind}_{dist_tok(b['bin'])}", f"{p}.corr_vs_distance[{i}] (channel pairs {b['bin']} mm apart)"
            F.add(f"{base}_n_pairs", count(b["n_pairs"]), b["n_pairs"], src + ".n_pairs")
            F.add(f"{base}_measured_r", r2(b["median_r_a"]), b["median_r_a"], src + ".median_r_a (median correlation, measured)")
            F.add(f"{base}_model_r", r2(b["median_r_b"]), b["median_r_b"], src + ".median_r_b (median correlation, model)")
            F.add(f"{base}_agreement", r2(b["agreement"]), b["agreement"], src + ".agreement (Pearson r of the pairs' correlations, "
                  "model vs measured)")
            nk = f"brain.{kind}.corr_vs_distance.{b['bin']}.agreement"
            if nk in null:
                lo, _, hi = null[nk]
                F.add(f"{base}_agreement_null_ci", ci_text(lo, hi), [lo, hi], f"{COV} :: structure_null['{nk}'] (exact model, "
                      "2.5th-97.5th percentiles)")
    # eigenstructure and dominant subspaces (normalised channels)
    for sk in ("mag", "grad", "combined_eigen"):
        x, p = s["brain"][sk], f"{COV} :: structure.brain.{sk}"
        st = sk.replace("_eigen", "")
        for side, who in (("eig_a", "measured"), ("eig_b", "model")):
            e, src = x[side], f"{p}.{side} (eigenvalue spectrum, {who})"
            for k in ("top1", "top5", "top10"):
                F.add(f"rev_cov_eig_{st}_{who}_{k}_pct", pct(100 * e[k]), 100 * e[k], f"{src}.{k} (variance share)")
            for k in ("n_for_50pct", "n_for_90pct", "n_for_95pct"):
                F.add(f"rev_cov_eig_{st}_{who}_{k}", count(e[k]), e[k], f"{src}.{k} (components for this share)")
            F.add(f"rev_cov_eig_{st}_{who}_participation", num(e["participation_ratio"], 1), e["participation_ratio"],
                  f"{src}.participation_ratio")
        for k, o in x["subspace"].items():
            base, src = f"rev_cov_sub_{st}_k{k}", f"{p}.subspace['{k}'] (top-{k} subspaces, model vs measured)"
            F.add(f"{base}_overlap", r2(o["overlap"]), o["overlap"], src + ".overlap (mean cos^2 of the principal angles)")
            F.add(f"{base}_random", sig2(o["random_overlap"]), o["random_overlap"], src + ".random_overlap (k/n: a random subspace)")
            F.add(f"{base}_captured_pct", pct(100 * o["captured_variance"]), 100 * o["captured_variance"], src + ".captured_variance")
            F.add(f"{base}_largest_angle_deg", num(o["largest_angle_deg"], 0), o["largest_angle_deg"], src + ".largest_angle_deg")
            sh = s["brain"]["split_half"][sk]["subspace"][k]["overlap"]
            F.add(f"{base}_split_half_overlap", r2(sh), sh, f"{COV} :: structure.brain.split_half.{sk}.subspace['{k}'].overlap "
                  "(one half of the windows vs the other)")
            nk = f"brain.{sk}.subspace.{k}.overlap"
            if nk in null:
                lo, _, hi = null[nk]
                F.add(f"{base}_overlap_null_ci", ci_text(lo, hi), [lo, hi], f"{COV} :: structure_null['{nk}'] (exact model)")
    # the measured excess over the model
    for ck in ("brain", "total", "brain_without_heart"):
        for sk, x in s[ck]["excess"].items():
            base, p = f"rev_cov_excess_{COV_COMPS[ck]}_{sk}", f"{COV} :: structure.{ck}.excess.{sk}"
            F.ratio(f"{base}_pos_over_model", x["positive_part_over_model_trace"], f"{p}.positive_part_over_model_trace "
                    "(positive part of measured - model over the model's trace)")
            for k, v in x["top_k_share_of_positive_part"].items():
                F.add(f"{base}_top{k}_pct", pct(100 * v), 100 * v, f"{p}.top_k_share_of_positive_part['{k}']")
            for j_, v in x["leading_component_external_share"].items():
                F.add(f"{base}_lead{j_}_external_pct", pct(100 * v), 100 * v, f"{p}.leading_component_external_share['{j_}'] "
                      "(share inside the 8-dimensional external subspace)")
            F.add(f"{base}_external_random_pct", pct(100 * x["external_share_random"]), 100 * x["external_share_random"],
                  f"{p}.external_share_random (a random direction)")
            if ck == "brain":  # what an exact model, measured the same way, gives for the same statistics
                for nk, nm, f in ((f"brain.excess.{sk}.positive_part_over_model_trace", "pos_over_model", r2),
                                  (f"brain.excess.{sk}.top_k_share_of_positive_part.3", "top3", None),
                                  (f"brain.excess.{sk}.leading_component_external_share.1", "lead1_external", None)):
                    if nk in null:
                        lo, _, hi = null[nk]
                        val, raw = ((f"{f(lo)} to {f(hi)}", [lo, hi]) if f else (f"{pct(100 * lo)} to {pct(100 * hi)}",
                                                                                [100 * lo, 100 * hi]))
                        F.add(f"{base}_{nm}_null_range", val, raw, f"{COV} :: structure_null['{nk}'] (exact model, 2.5th-97.5th "
                              "percentiles)")
    for vk, vt in VARIANTS.items():
        for kind, x in c["structure_variants"][vk].items():
            base, p = f"rev_cov_variant_{vt}_{kind}", f"{COV} :: structure_variants['{vk}'].{kind}"
            F.ratio(f"{base}_amp", x["amp_ratio_of_medians"], f"{p}.amp_ratio_of_medians (model over measured)", three=True)
            for st in ("p5", "p50", "p95"):
                F.add(f"{base}_var_ratio_{st}", r2(x["var_ratio"][st]), x["var_ratio"][st], f"{p}.var_ratio.{st}")
            F.add(f"{base}_log_var_r", r2(x["log_var_pearson"]), x["log_var_pearson"], f"{p}.log_var_pearson")
            F.add(f"{base}_corr_pattern_r", r2(x["corr_pattern_pearson"]), x["corr_pattern_pearson"], f"{p}.corr_pattern_pearson")


def cov_detect_facts(F, c):
    dt, p = c["detectability"], f"{COV} :: detectability"
    for sk, n in dt["channels"].items():
        F.add(f"rev_cov_det_n_channels_{sk}", count(n), n, f"{p}.channels.{sk}")
    for key, e in dt["comparisons"].items():
        parts = key.split("/")
        what, st = parts[0], "_".join(parts[1:]).replace("projected", "proj")
        if what == "finite_sample_bias":
            src = f"{p}.comparisons['{key}'] (finite-sample bias of detectability with a covariance estimated from these samples)"
            F.ratio(f"rev_cov_det_{st}_bias", e["median"], src + ".median", three=True)
            F.add(f"rev_cov_det_{st}_bias_p5_p95_range", f"{r3(e['p5'])} to {r3(e['p95'])}", [e["p5"], e["p95"]], src + ".p5, .p95")
            F.ratio(f"rev_cov_det_{st}_bias_hybrid", e["hybrid_median"], src + ".hybrid_median", three=True)
            continue
        cov_cmp(F, f"rev_cov_det_{st}_{DET[what]}", e, f"detectability.comparisons['{key}']",
                note="Neuromag detectability with the measured (or hybrid) covariance over that with the model's")
    for key, v in dt["surrogate_medians"].items():
        what, sk = key.split("/")
        base = f"rev_cov_surr_{what}_{sk}"
        src = (f"{p}.surrogate_medians['{key}'] (median ratio over the targets per surrogate recording: the 50th percentile over "
               "the surrogates, and the 2.5th to 97.5th as interval and range)")
        F.ratio(base, v[1], src, [v[0], v[2]], three=True)
        F.add(f"{base}_range", f"{r2(v[0])} to {r2(v[2])}", [v[0], v[2]], src)
        F.add(f"{base}_range_3dp", f"{r3(v[0])} to {r3(v[2])}", [v[0], v[2]], src)
    meds = []
    for i, b in enumerate(dt["by_depth_combined"]):
        if b["median"] is None:
            continue
        base = f"rev_cov_det_depth_{b['lo']:.0f}_{b['hi']:.0f}"
        src = f"{p}.by_depth_combined[{i}] (targets {b['lo']:g}-{b['hi']:g} mm deep: measured-corrected over model, combined)"
        F.add(f"{base}_median", r2(b["median"]), b["median"], src + ".median")
        F.add(f"{base}_iqr_range", f"{r2(b['q25'])} to {r2(b['q75'])}", [b["q25"], b["q75"]], src + ".q25, .q75")
        F.add(f"{base}_n", count(b["n"]), b["n"], src + ".n")
        meds.append(b["median"])
    F.span("rev_cov_det_depth_median", meds, r2, f"{p}.by_depth_combined[*].median (populated 5-mm bins)")
    lobes = {}
    for lb, x in dt["by_lobe_combined"].items():
        F.add(f"rev_cov_det_lobe_{lb}_n_targets", count(x["n_targets"]), x["n_targets"], f"{p}.by_lobe_combined['{lb}'].n_targets")
        for k, nm in (("measured_corrected", "corrected"), ("without_heart_corrected", "noheart_corrected")):
            cov_cmp(F, f"rev_cov_det_lobe_{lb}_{nm}", x[k], f"detectability.by_lobe_combined['{lb}'].{k}",
                    note="Neuromag combined, measured covariance (finite-sample corrected) over the model's")
        lobes[lb] = x["measured_corrected"]["ratio"]
    F.span("rev_cov_det_lobe_corrected", list(lobes.values()), r2, f"{p}.by_lobe_combined[*].measured_corrected.ratio (lobes)",
           extra=(("3dp", r3),))


def cov_opm_facts(F, c):
    o, p = c["opm_implication"], f"{COV} :: opm_implication"
    mx = max(o["checks"]["g2_opm_detectability_max_abs_dev_vs_stored_csv"].values())
    F.add("rev_cov_check_g2_opm_detect_max_dev", small(mx), mx, f"{p}.checks.g2_opm_detectability_max_abs_dev_vs_stored_csv (largest)")
    scen = {}
    for key, e in o["ratios"].items():
        arr, sk, rest = key.split("/", 2)
        tok = SCEN[rest]
        cond = ("sensor + cortical background (the G2 headline condition)" if rest.endswith("intrinsic+brain") else
                "projected; Neuromag measured rows not finite-sample corrected" if rest.startswith("projected/") else
                "sensor + background + room field, not projected")
        cov_cmp(F, f"rev_cov_opm_{ARR[arr]}_vs_{sk}_{tok}", e, f"opm_implication.ratios['{key}']",
                note=f"d_OPM/d_Neuromag, {cond}")
        if tok in ("s1", "s2", "s3", "s4", "s5"):
            scen.setdefault((ARR[arr], sk), []).append(e["ratio"])
    for (arr, sk), vals in scen.items():
        F.span(f"rev_cov_opm_{arr}_vs_{sk}_scenarios", vals, r2, f"{p}.ratios['opm_{arr}/{sk}/S1...S5'].ratio (the five scenarios)",
               extra=(("3dp", r3),))


def cov_band_facts(F, c):
    sb, amps, excess = c["sub_bands"], {}, {}
    for bk, x in sb.items():
        if not bk.endswith("Hz"):
            continue
        lo, hi = bk.removesuffix("Hz").split("-")
        bt, p = f"b{lo}_{hi}", f"{COV} :: sub_bands['{bk}']"
        F.add(f"rev_cov_band_{bt}_label", f"{lo}–{hi}", [float(lo), float(hi)], f"{p} (band, Hz)")
        F.add(f"rev_cov_band_{bt}_enbw_hz", num(x["enbw_hz"], 1), x["enbw_hz"], f"{p}.enbw_hz")
        m, gr = x["mag"], x["grad"]
        F.ratio(f"rev_cov_band_{bt}_mag_amp", m["brain_amp_ratio_model_over_measured"], f"{p}.mag.brain_amp_ratio_model_over_measured "
                "(median magnetometer brain RMS, model over measured; model scale refitted per band on the gradiometers)", three=True)
        F.ratio(f"rev_cov_band_{bt}_mag_emptyroom_amp", m["empty_room_amp_ratio_model_over_measured"],
                f"{p}.mag.empty_room_amp_ratio_model_over_measured", three=True)
        F.ratio(f"rev_cov_band_{bt}_grad_emptyroom_amp", gr["empty_room_amp_ratio_model_over_measured"],
                f"{p}.grad.empty_room_amp_ratio_model_over_measured", three=True)
        for kind, y in (("mag", m), ("grad", gr)):
            for k, nm in (("brain_log_var_pearson", "log_var_r"), ("brain_corr_pattern_pearson", "corr_pattern_r"),
                          ("brain_subspace_overlap_k5", "overlap_k5"), ("split_half_log_var_pearson", "split_half_log_var_r")):
                F.add(f"rev_cov_band_{bt}_{kind}_{nm}", r2(y[k]), y[k], f"{p}.{kind}.{k}")
            F.add(f"rev_cov_band_{bt}_{kind}_n_nonpositive", count(y["n_nonpositive_brain"]), y["n_nonpositive_brain"],
                  f"{p}.{kind}.n_nonpositive_brain")
        for side in ("measured", "model"):
            v = m[f"brain_{side}"]
            F.add(f"rev_cov_band_{bt}_mag_brain_{side}_ft", num(v["median_amplitude"], 0), v["median_amplitude"],
                  f"{p}.mag.brain_{side}.median_amplitude (median magnetometer RMS, fT)")
            F.add(f"rev_cov_band_{bt}_mag_brain_{side}_asd", num(v["median_asd"], 1), v["median_asd"],
                  f"{p}.mag.brain_{side}.median_asd (fT/sqrt(Hz))")
        F.add(f"rev_cov_band_{bt}_grad_brain_ftcm", num(gr["brain_measured"]["median_amplitude"], 1),
              gr["brain_measured"]["median_amplitude"], f"{p}.grad.brain_measured.median_amplitude (fT/cm; the model equals it)")
        F.ratio(f"rev_cov_band_{bt}_mag_over_grad_excess", x["mag_over_grad_brain_variance_measured_over_model"],
                f"{p}.mag_over_grad_brain_variance_measured_over_model (measured magnetometer/gradiometer brain-variance ratio over "
                "the model's)")
        lo_, _, hi_ = x["null95"]["mag.brain_amp_ratio_model_over_measured"]
        F.add(f"rev_cov_band_{bt}_mag_amp_null_ci", ci_text(lo_, hi_), [lo_, hi_],
              f"{p}.null95['mag.brain_amp_ratio_model_over_measured'] (exact model, 2.5th-97.5th percentiles)")
        if bk != "1-40Hz":
            amps[bt] = m["brain_amp_ratio_model_over_measured"]
            excess[bt] = x["mag_over_grad_brain_variance_measured_over_model"]
    src = f"{COV} :: sub_bands[*].mag.brain_amp_ratio_model_over_measured (the six sub-bands)"
    F.span("rev_cov_band_mag_amp", list(amps.values()), r2, src, extra=(("3dp", r3),))
    F.span("rev_cov_band_mag_over_grad_excess", list(excess.values()), r2,
           f"{COV} :: sub_bands[*].mag_over_grad_brain_variance_measured_over_model (the six sub-bands)")
    for k, v in sb["closure_median"].items():
        rec, kind = k.split("/")
        F.ratio(f"rev_cov_band_closure_{rec.replace('_', '')}_{kind}", v, f"{COV} :: sub_bands.closure_median['{k}'] (per channel, sum "
                "of the sub-band variances over the full-band variance; median over channels)", three=True)


def cov_facts(F, root):
    c = load(root, COV)
    cov_setup_facts(F, c)
    cov_variance_facts(F, c)
    cov_heart_facts(F, c)
    cov_structure_facts(F, c)
    cov_detect_facts(F, c)
    cov_opm_facts(F, c)
    cov_band_facts(F, c)


# ------------------------------------------------------------------------------------------------
# the children's anatomy against their MRIs (rev_qc_)
QC_TOK = {"adult": "adult", "infant2yr": "infant2yr", "childA": "child_a", "childB": "child_b", "childC": "child_c"}
CHILDREN = ("childA", "childB", "childC")
EDGES = {"otsu": "otsu", "half_max": "halfmax", "steepest": "steepest"}
HC = {"surface_used": "surface", "mri_boundary_points": "mri", "t1_mask": "mask"}


def qc_distance_facts(F, base, x, path, kind):
    """Distances of the white surface to a head surface (white_to_scalp_used or white_to_mri_boundary)."""
    src = f"{QC} :: {path} (white-surface vertex to {kind}, mm)"
    F.add(f"{base}_min_mm", mm(x["min_mm"]), x["min_mm"], src + ".min_mm (the closest vertex)")
    F.add(f"{base}_min_mm_2dp", mm2(x["min_mm"]), x["min_mm"], src + ".min_mm")
    for k in ("p0.1", "p1", "p5", "p10", "p50"):
        F.add(f"{base}_{k.replace('.', '')}_mm", mm(x[f"{k}_mm"]), x[f"{k}_mm"], f"{src}.{k}_mm")
    for k, nm in (("below_8mm", "lt8mm"), ("below_10mm", "lt10mm")):
        y = x[k]
        F.add(f"{base}_{nm}_n_vertices", count(y["n_vertices"]), y["n_vertices"], f"{src}.{k}.n_vertices")
        F.add(f"{base}_{nm}_area_cm2", num(y["area_cm2"], 1), y["area_cm2"], f"{src}.{k}.area_cm2")
        F.add(f"{base}_{nm}_n_usable", count(y["n_usable_vertices"]), y["n_usable_vertices"],
              f"{src}.{k}.n_usable_vertices (vertices of the usable cortex: at least 4 mm inside the inner skull)")
        if y.get("parcels"):
            lab = {}
            for pk, pv in y["parcels"].items():
                lab[pk.split(".", 1)[1]] = lab.get(pk.split(".", 1)[1], 0.0) + pv["area_cm2"]
            top = max(lab, key=lab.get)
            share = 100 * lab[top] / y["area_cm2"]
            F.add(f"{base}_{nm}_top_parcel", top, top, f"{QC} :: derived: {path}.{k}.parcels, the Desikan-Killiany label (both "
                  "hemispheres) with the largest area")
            F.add(f"{base}_{nm}_top_parcel_pct", pct(share), share, f"{QC} :: derived: {path}.{k}.parcels['?h.{top}'].area_cm2 "
                  f"summed over hemispheres / {path}.{k}.area_cm2")
    for sec, v in x["p5_by_sector_mm"].items():
        F.add(f"{base}_p5_{sec}_mm", mm(v), v, f"{src}.p5_by_sector_mm.{sec} (5th percentile in the {sec} sector)")


def qc_anatomy_facts(F, key, a):
    t, p = QC_TOK[key], f"{QC} :: anatomies['{key}']"
    F.add(f"rev_qc_{t}_subject", a["subject"], a["subject"], f"{p}.subject ({a['description']})")
    F.add(f"rev_qc_{t}_t1_available", yes(a["t1_available"]), a["t1_available"], f"{p}.t1_available")
    F.add(f"rev_qc_{t}_scalp_n_vertices", count(a["scalp"]["n_vertices"]), a["scalp"]["n_vertices"], f"{p}.scalp.n_vertices")
    F.add(f"rev_qc_{t}_scalp_cap_vertices", count(a["scalp"]["cap_vertices"]), a["scalp"]["cap_vertices"],
          f"{p}.scalp.cap_vertices (above the fiducial plane)")
    # frames
    fr = a["frames"]
    F.add(f"rev_qc_{t}_headers_match_t1", yes(fr["headers_match_t1"]), fr["headers_match_t1"], f"{p}.frames.headers_match_t1")
    F.add(f"rev_qc_{t}_frame_max_disp_mm", small(fr["max_corner_displacement_vs_t1_mm"]), fr["max_corner_displacement_vs_t1_mm"],
          f"{p}.frames.max_corner_displacement_vs_t1_mm (largest displacement of a volume corner, surface geometry vs T1)")
    F.add(f"rev_qc_{t}_frame_max_disp_files_mm", small(fr["max_corner_displacement_between_files_mm"]),
          fr["max_corner_displacement_between_files_mm"], f"{p}.frames.max_corner_displacement_between_files_mm")
    dims = fr["t1"]["dims"]
    F.add(f"rev_qc_{t}_t1_dims", f"{dims[0]}³" if len(set(dims)) == 1 else " × ".join(map(str, dims)), dims, f"{p}.frames.t1.dims")
    vox = fr["t1"]["voxel_mm"]
    F.add(f"rev_qc_{t}_t1_voxel_mm", g(vox[0]) if len(set(vox)) == 1 else " × ".join(map(g, vox)), vox, f"{p}.frames.t1.voxel_mm")
    faces = fr["head_mask_voxels_on_volume_faces"]
    F.add(f"rev_qc_{t}_mask_voxels_inferior_face", count(faces["S-"]), faces["S-"], f"{p}.frames.head_mask_voxels_on_volume_faces['S-'] "
          "(head-mask voxels on the volume's inferior face)")
    other = sum(v for k, v in faces.items() if k != "S-")
    F.add(f"rev_qc_{t}_mask_voxels_other_faces", count(other), other, f"{p}.frames.head_mask_voxels_on_volume_faces (the five other faces)")
    pf = a.get("prepared_files")
    if pf:
        for k, nm in (("head_fif_vs_seghead_max_um", "head_vs_seghead"), ("source_space_vs_white_max_um", "source_space_vs_white"),
                      ("bem_head_on_scalp_max_um", "bem_head_on_scalp")):
            F.add(f"rev_qc_{t}_prepared_{nm}_um", small(pf[k]), pf[k], f"{p}.prepared_files.{k} (largest vertex distance, µm)")
    # cortex near the scalp and near the MRI head boundary
    qc_distance_facts(F, f"rev_qc_{t}_white_scalp", a["white_to_scalp_used"], f"anatomies['{key}'].white_to_scalp_used",
                      "the dense scalp used by G3B/G4")
    qc_distance_facts(F, f"rev_qc_{t}_white_mri", a["white_to_mri_boundary"], f"anatomies['{key}'].white_to_mri_boundary",
                      "the MRI head boundary")
    cw = a["closest_white_vertex"]
    if abs(cw["distance_mm"] - a["white_to_scalp_used"]["min_mm"]) > 1e-9:
        raise ValueError(f"{QC}: anatomies['{key}'].closest_white_vertex is not the minimum of white_to_scalp_used")
    F.add(f"rev_qc_{t}_closest_region", cw["region"], cw["region"], f"{p}.closest_white_vertex.region (hemisphere.parcel of the "
          "white vertex closest to the dense scalp)")
    wm = a["white_to_mri_boundary"]
    for k in ("below_8mm", "below_10mm"):  # the parcel holding the white vertex closest to the MRI head boundary, where stored
        hit = [pk for pk, pv in (wm[k].get("parcels") or {}).items() if abs(pv["min_mm"] - wm["min_mm"]) < 1e-9]
        if hit:
            F.add(f"rev_qc_{t}_closest_mri_region", hit[0], hit[0], f"{QC} :: derived: anatomies['{key}'].white_to_mri_boundary.{k}."
                  "parcels, the parcel whose min_mm equals white_to_mri_boundary.min_mm (the white vertex closest to the MRI boundary)")
            break
    for surf, x in a["watershed_vs_scalp"].items():
        F.add(f"rev_qc_{t}_watershed_{surf}_cap_mm", smm(x["cap"]["median_mm"]), x["cap"]["median_mm"],
              f"{p}.watershed_vs_scalp.{surf}.cap.median_mm (watershed BEM surface minus the dense scalp, signed, mm)")
    # the scalp against the MRI head boundary
    sv = a["scalp_vs_mri"]
    for region, rt in (("cap", "cap"), ("whole_head", "head")):
        for edge, et in EDGES.items():
            x, src = sv["offsets"][region][edge], f"{p}.scalp_vs_mri.offsets.{region}.{edge}"
            stats = ("p5", "p25", "p50", "p75", "p95") if region == "cap" else ("p5", "p50", "p95")
            what = " (MRI head boundary minus the surface along its outward normal, mm: + where the boundary lies outside)"
            F.add(f"rev_qc_{t}_offset_{rt}_{et}_median_mm", smm(x["median_mm"]), x["median_mm"], src + ".median_mm" + what)
            F.add(f"rev_qc_{t}_offset_{rt}_{et}_median_mm_2dp", smm2(x["median_mm"]), x["median_mm"], src + ".median_mm" + what)
            for st in stats:
                if st != "p50":
                    F.add(f"rev_qc_{t}_offset_{rt}_{et}_{st}_mm", smm(x[f"{st}_mm"]), x[f"{st}_mm"], f"{src}.{st}_mm")
            F.add(f"rev_qc_{t}_offset_{rt}_{et}_p5_p95_range", f"{smm(x['p5_mm'])} to {smm(x['p95_mm'])}", [x["p5_mm"], x["p95_mm"]],
                  f"{src}.p5_mm, .p95_mm")
            if region == "cap":
                for k, nm in (("share_above_1mm", "above1mm"), ("share_above_2mm", "above2mm"), ("share_above_4mm", "above4mm"),
                              ("share_below_minus_1mm", "below_minus1mm"), ("determined_share", "determined")):
                    F.add(f"rev_qc_{t}_offset_cap_{et}_{nm}_pct", pct(100 * x[k]), 100 * x[k], f"{src}.{k}")
    dc, src = sv["decomposition_cap"], f"{p}.scalp_vs_mri.decomposition_cap (cap offsets as uniform + translation + rotation)"
    F.add(f"rev_qc_{t}_offset_cap_uniform_mm", smm(dc["uniform_mm"]), dc["uniform_mm"], src + ".uniform_mm")
    F.add(f"rev_qc_{t}_offset_cap_uniform_mm_2dp", smm2(dc["uniform_mm"]), dc["uniform_mm"], src + ".uniform_mm")
    F.add(f"rev_qc_{t}_offset_cap_translation_mm", mm2(dc["translation_norm_mm"]), dc["translation_norm_mm"], src + ".translation_norm_mm")
    F.add(f"rev_qc_{t}_offset_cap_rotation_deg", mm2(dc["rotation_norm_deg"]), dc["rotation_norm_deg"], src + ".rotation_norm_deg")
    F.add(f"rev_qc_{t}_offset_cap_rotation_disp_mm", mm2(dc["rotation_max_displacement_mm"]), dc["rotation_max_displacement_mm"],
          src + ".rotation_max_displacement_mm")
    F.add(f"rev_qc_{t}_offset_cap_residual_mm", mm2(dc["residual_rms_mm"]), dc["residual_rms_mm"], src + ".residual_rms_mm")
    F.add(f"rev_qc_{t}_offset_cap_explained_pct", pct(100 * dc["explained_share"]), 100 * dc["explained_share"], src + ".explained_share")
    sm = a.get("scalp_vs_seghead_mask")
    if sm:
        x, src = sm["mask_edge_minus_surface"], f"{p}.scalp_vs_seghead_mask.mask_edge_minus_surface (head-mask edge minus the surface, cap, mm)"
        F.add(f"rev_qc_{t}_mask_vs_scalp_median_mm", smm2(x["median_mm"]), x["median_mm"], src + ".median_mm")
        F.add(f"rev_qc_{t}_mask_vs_scalp_p5_p95_range", f"{smm2(x['p5_mm'])} to {smm2(x['p95_mm'])}", [x["p5_mm"], x["p95_mm"]],
              src + ".p5_mm, .p95_mm")
        y = sm["t1_edge_minus_mask_edge_cap"]
        F.add(f"rev_qc_{t}_t1_vs_mask_edge_median_mm", smm(y["median_mm"]), y["median_mm"],
              f"{p}.scalp_vs_seghead_mask.t1_edge_minus_mask_edge_cap.median_mm (T1 Otsu edge minus the mask's edge, cap, mm)")
    mb, src = a["mri_boundary_surface"], f"{p}.mri_boundary_surface (the dense scalp moved to the MRI head boundary)"
    for st in ("p5", "p50", "p95"):
        F.add(f"rev_qc_{t}_boundary_disp_cap_{st}_mm", smm(mb["displacement_cap_mm"][st]), mb["displacement_cap_mm"][st],
              f"{src}.displacement_cap_mm.{st}")
        F.add(f"rev_qc_{t}_boundary_remeasured_{st}_mm", smm2(mb["remeasured_cap_offset_mm"][st]), mb["remeasured_cap_offset_mm"][st],
              f"{src}.remeasured_cap_offset_mm.{st} (offsets measured again on the moved surface)")
    F.add(f"rev_qc_{t}_boundary_clipped", count(mb["clipped"]), mb["clipped"], f"{src}.clipped")
    F.add(f"rev_qc_{t}_boundary_outliers", count(mb["outliers_replaced"]), mb["outliers_replaced"], f"{src}.outliers_replaced")
    F.add(f"rev_qc_{t}_boundary_flipped_cap", count(mb["folds"]["flipped_in_cap"]), mb["folds"]["flipped_in_cap"],
          f"{src}.folds.flipped_in_cap (faces folded over in the cap after relaxation)")
    F.add(f"rev_qc_{t}_boundary_flipped_cap_mm2", num(mb["folds"]["flipped_area_in_cap_mm2"], 1), mb["folds"]["flipped_area_in_cap_mm2"],
          f"{src}.folds.flipped_area_in_cap_mm2")
    # head circumference on the surface used, the MRI boundary and the T1 head mask
    hc = a["head_circumference"]
    for k, ht in HC.items():
        x, src = hc[k], f"{p}.head_circumference.{k}"
        F.add(f"rev_qc_{t}_ofc_{ht}_cm", num(x["ofc_mm"] / 10, 1), x["ofc_mm"] / 10, f"{src}.ofc_mm / 10 (occipitofrontal, cm)")
        F.add(f"rev_qc_{t}_ofc_{ht}_mm", mm(x["ofc_mm"]), x["ofc_mm"], f"{src}.ofc_mm")
        for q in ("breadth_mm", "length_mm", "vertex_height_mm", "ofc_plane_height_mm"):
            F.add(f"rev_qc_{t}_{q.removesuffix('_mm')}_{ht}_mm", mm(x[q]), x[q], f"{src}.{q}")
        F.add(f"rev_qc_{t}_cap_volume_{ht}_cm3", count(x["cap_volume_cm3"]), x["cap_volume_cm3"], f"{src}.cap_volume_cm3")
    for ht, k in (("mri", "mri_boundary_points"), ("mask", "t1_mask")):
        dv = hc[k]["ofc_mm"] - hc["surface_used"]["ofc_mm"]
        F.add(f"rev_qc_{t}_ofc_{ht}_minus_surface_mm", smm(dv), dv, f"{QC} :: derived: anatomies['{key}'].head_circumference.{k}.ofc_mm "
              "- .surface_used.ofc_mm")
    F.add(f"rev_qc_{t}_white_area_cm2", count(hc["white_surface_area_cm2"]), hc["white_surface_area_cm2"],
          f"{p}.head_circumference.white_surface_area_cm2")
    F.add(f"rev_qc_{t}_ofc_per_sqrt_white_area", r2(hc["ofc_per_sqrt_white_area"]), hc["ofc_per_sqrt_white_area"],
          f"{p}.head_circumference.ofc_per_sqrt_white_area (cm per cm)")
    # the white surface against the T1
    wt, src = a["white_vs_t1"], f"{p}.white_vs_t1 (T1 edge across the white surface)"
    F.add(f"rev_qc_{t}_white_edge_median_mm", smm2(wt["edge_median_mm"]), wt["edge_median_mm"], src + ".edge_median_mm")
    F.add(f"rev_qc_{t}_white_edge_p5_p95_range", f"{smm(wt['edge']['p5'])} to {smm(wt['edge']['p95'])}",
          [wt["edge"]["p5"], wt["edge"]["p95"]], src + ".edge.p5, .p95 (mm)")
    F.add(f"rev_qc_{t}_white_translation_mm", mm2(wt["decomposition"]["translation_norm_mm"]), wt["decomposition"]["translation_norm_mm"],
          src + ".decomposition.translation_norm_mm")
    F.add(f"rev_qc_{t}_white_rotation_disp_mm", mm2(wt["decomposition"]["rotation_max_displacement_mm"]),
          wt["decomposition"]["rotation_max_displacement_mm"], src + ".decomposition.rotation_max_displacement_mm")
    F.add(f"rev_qc_{t}_white_best_shift_mm", mm(wt["best_shift_norm_mm"]), wt["best_shift_norm_mm"], src + ".best_shift_norm_mm "
          "(rigid shift of the white surface that maximises its contrast)")
    F.add(f"rev_qc_{t}_white_contrast_gain_pct", pct(100 * wt["contrast_gain_share"]), 100 * wt["contrast_gain_share"],
          src + ".contrast_gain_share")
    for k, nm in (("all", ""), ("below_8mm", "_lt8mm")):
        x = a["white_contrast_near_scalp"].get(k)
        if x:
            src = f"{p}.white_contrast_near_scalp.{k} (T1 1 mm inside minus 1 mm outside the white surface)"
            F.add(f"rev_qc_{t}_white_contrast{nm}_median", num(x["contrast_median"], 1), x["contrast_median"], src + ".contrast_median")
            F.add(f"rev_qc_{t}_white_no_contrast{nm}_pct", pct(100 * x["share_without_contrast"]), 100 * x["share_without_contrast"],
                  src + ".share_without_contrast")
            F.add(f"rev_qc_{t}_white_contrast{nm}_n", count(x["n"]), x["n"], src + ".n")
    # the G3B targets
    tg, src = a["g3b_targets"], f"{p}.g3b_targets"
    F.add(f"rev_qc_{t}_targets_n", count(tg["n_targets"]), tg["n_targets"], src + ".n_targets")
    for st, v in tg["depth"].items():
        F.add(f"rev_qc_{t}_targets_depth_{st}_mm", mm(v), v, f"{src}.depth.{st} (below the scalp used, mm)")
    for k, nm in (("below_8mm", "lt8mm"), ("below_10mm", "lt10mm")):
        F.add(f"rev_qc_{t}_targets_{nm}_n", count(tg[k]["n"]), tg[k]["n"], f"{src}.{k}.n (below the scalp used)")
        F.add(f"rev_qc_{t}_targets_{nm}_area_cm2", num(tg[k]["area_cm2"], 1), tg[k]["area_cm2"], f"{src}.{k}.area_cm2")
        F.add(f"rev_qc_{t}_targets_{nm}_mri_n", count(tg[f"{k}_to_mri_boundary"]), tg[f"{k}_to_mri_boundary"],
              f"{src}.{k}_to_mri_boundary (below the MRI head boundary)")
    for st, v in tg["depth_change_to_mri_boundary_mm"].items():
        F.add(f"rev_qc_{t}_targets_depth_change_{st}_mm", smm(v), v, f"{src}.depth_change_to_mri_boundary_mm.{st} (depth below the "
              "MRI boundary minus below the scalp used)")
    # verdict and its checks
    v, src = a["verdict"], f"{p}.verdict"
    ck = v["checks"]
    for edge, et in (("otsu", "otsu"), ("half_max", "halfmax")):
        x = ck["scalp_offset_excess_over_adult_mm"][edge]
        F.add(f"rev_qc_{t}_excess_{et}_mm", smm(x), x, f"{src}.checks.scalp_offset_excess_over_adult_mm.{edge} (cap median offset "
              "minus the adult's, mm)")
        F.add(f"rev_qc_{t}_excess_{et}_mm_2dp", smm2(x), x, f"{src}.checks.scalp_offset_excess_over_adult_mm.{edge}")
    for side in ("outside", "inside"):
        x = ck["cap_share_beyond_tolerance_of_adult_median"][side]
        F.add(f"rev_qc_{t}_cap_{side}_tolerance_pct", pct(100 * x), 100 * x, f"{src}.checks.cap_share_beyond_tolerance_of_adult_median."
              f"{side} (share of the cap where the surface lies further {side[:-4]} than the adult's position by more than the "
              "tolerance)")
    for ref, x in (v["surface_only_indications"].get("p5_shortfall_vs_references_mm") or {}).items():
        F.add(f"rev_qc_{t}_p5_shortfall_vs_{QC_TOK[ref]}_mm", smm(x), x, f"{src}.surface_only_indications.p5_shortfall_vs_references_mm."
              f"{ref} (5th-percentile white-to-scalp distance of the reference minus this anatomy's, mm)")
    cq = v["consequence"]
    F.add(f"rev_qc_{t}_opm_face_to_mri_mm", mm(cq["opm_cell_inner_face_to_mri_boundary_mm"]), cq["opm_cell_inner_face_to_mri_boundary_mm"],
          f"{src}.consequence.opm_cell_inner_face_to_mri_boundary_mm (2 mm intended)")
    F.add(f"rev_qc_{t}_neuromag_contact_to_mri_mm", mm(cq["neuromag_top_contact_to_mri_boundary_mm"]),
          cq["neuromag_top_contact_to_mri_boundary_mm"], f"{src}.consequence.neuromag_top_contact_to_mri_boundary_mm (20 mm intended)")


def qc_facts(F, root):
    q = load(root, QC)
    p = f"{QC} :: parameters"
    pa = q["parameters"]
    F.add("rev_qc_tolerance_mm", g(pa["tolerance_mm"]), pa["tolerance_mm"], f"{p}.tolerance_mm (one voxel of the 1-mm T1)")
    F.add("rev_qc_header_tolerance_mm", g(pa["header_tolerance_mm"]), pa["header_tolerance_mm"], f"{p}.header_tolerance_mm")
    F.add("rev_qc_near_mm", listing(pa["near_mm"]), pa["near_mm"], f"{p}.near_mm ({pa['note'][:80]}...)")
    F.add("rev_qc_smooth_rings", count(pa["smooth_rings"]), pa["smooth_rings"], f"{p}.smooth_rings (profiles smoothed over the "
          "mesh neighbourhood)")
    lo, hi, step = pa["profile_mm"]
    F.add("rev_qc_profile_range_mm", f"{g(lo)} to {g(hi)}", [lo, hi], f"{p}.profile_mm (T1 profile along the outward normal)")
    F.add("rev_qc_profile_step_mm", g(step), step, f"{p}.profile_mm[2]")
    F.add("rev_qc_commit", q["provenance"]["commit"], q["provenance"]["commit"], f"{QC} :: provenance.commit (code of the run)")
    man = load(root, QC_MANIFEST)
    F.add("rev_qc_manifest_n_files", count(len(man["files"])), len(man["files"]), f"{QC_MANIFEST} :: files (T1 and head-mask volumes)")
    F.add("rev_qc_manifest_dataset", man["dataset"], man["dataset"], f"{QC_MANIFEST} :: dataset")
    usable = []
    for key, a in q["anatomies"].items():
        qc_anatomy_facts(F, key, a)
        vd = q["verdicts"][key]
        if (vd["class"], vd["text"]) != (a["verdict"]["class"], a["verdict"]["text"]):
            raise ValueError(f"{QC}: verdicts['{key}'] differs from anatomies['{key}'].verdict")
        F.add(f"rev_qc_{QC_TOK[key]}_verdict", vd["class"], vd["class"], f"{QC} :: verdicts['{key}'].class")
        F.add(f"rev_qc_{QC_TOK[key]}_verdict_text", vd["text"], vd["text"], f"{QC} :: verdicts['{key}'].text")
        if key in CHILDREN and vd["class"] == "usable":
            usable.append(key)
    F.add("rev_qc_n_anatomies", count(len(q["anatomies"])), len(q["anatomies"]), f"{QC} :: anatomies")
    F.add("rev_qc_n_children_usable", count(len(usable)), len(usable), f"{QC} :: verdicts[childA, childB, childC].class == "
          f"'usable' ({listing(usable, str)})")
    an = q["anatomies"]
    # the children as a group
    groups = (("white_scalp_min_mm", lambda a: a["white_to_scalp_used"]["min_mm"], mm, "white_to_scalp_used.min_mm"),
              ("white_scalp_p5_mm", lambda a: a["white_to_scalp_used"]["p5_mm"], mm, "white_to_scalp_used.p5_mm"),
              ("white_mri_min_mm", lambda a: a["white_to_mri_boundary"]["min_mm"], mm, "white_to_mri_boundary.min_mm"),
              ("white_mri_p5_mm", lambda a: a["white_to_mri_boundary"]["p5_mm"], mm, "white_to_mri_boundary.p5_mm"),
              ("white_scalp_lt8mm_area_cm2", lambda a: a["white_to_scalp_used"]["below_8mm"]["area_cm2"], lambda x: num(x, 1),
               "white_to_scalp_used.below_8mm.area_cm2"),
              ("white_mri_lt8mm_area_cm2", lambda a: a["white_to_mri_boundary"]["below_8mm"]["area_cm2"], lambda x: num(x, 1),
               "white_to_mri_boundary.below_8mm.area_cm2"),
              ("offset_cap_otsu_median_mm", lambda a: a["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"], smm,
               "scalp_vs_mri.offsets.cap.otsu.median_mm"),
              ("offset_cap_halfmax_median_mm", lambda a: a["scalp_vs_mri"]["offsets"]["cap"]["half_max"]["median_mm"], smm,
               "scalp_vs_mri.offsets.cap.half_max.median_mm"),
              ("excess_otsu_mm", lambda a: a["verdict"]["checks"]["scalp_offset_excess_over_adult_mm"]["otsu"], smm2,
               "verdict.checks.scalp_offset_excess_over_adult_mm.otsu"),
              ("excess_halfmax_mm", lambda a: a["verdict"]["checks"]["scalp_offset_excess_over_adult_mm"]["half_max"], smm2,
               "verdict.checks.scalp_offset_excess_over_adult_mm.half_max"),
              ("ofc_surface_cm", lambda a: a["head_circumference"]["surface_used"]["ofc_mm"] / 10, lambda x: num(x, 1),
               "head_circumference.surface_used.ofc_mm / 10"),
              ("ofc_mri_cm", lambda a: a["head_circumference"]["mri_boundary_points"]["ofc_mm"] / 10, lambda x: num(x, 1),
               "head_circumference.mri_boundary_points.ofc_mm / 10"),
              ("ofc_mri_minus_surface_mm", lambda a: (a["head_circumference"]["mri_boundary_points"]["ofc_mm"]
                                                      - a["head_circumference"]["surface_used"]["ofc_mm"]), smm,
               "head_circumference.mri_boundary_points.ofc_mm - .surface_used.ofc_mm"),
              ("targets_lt10mm_n", lambda a: a["g3b_targets"]["below_10mm"]["n"], count, "g3b_targets.below_10mm.n"),
              ("targets_lt10mm_mri_n", lambda a: a["g3b_targets"]["below_10mm_to_mri_boundary"], count,
               "g3b_targets.below_10mm_to_mri_boundary"),
              ("targets_depth_change_p50_mm", lambda a: a["g3b_targets"]["depth_change_to_mri_boundary_mm"]["p50"], smm,
               "g3b_targets.depth_change_to_mri_boundary_mm.p50"),
              ("opm_face_to_mri_mm", lambda a: a["verdict"]["consequence"]["opm_cell_inner_face_to_mri_boundary_mm"], mm,
               "verdict.consequence.opm_cell_inner_face_to_mri_boundary_mm"))
    for nm, get, fmt, path in groups:
        F.span(f"rev_qc_children_{nm}", [get(an[k]) for k in CHILDREN], fmt, f"{QC} :: anatomies[childA, childB, childC].{path}")
    mx = max(an[k]["frames"]["max_corner_displacement_vs_t1_mm"] for k in CHILDREN)
    F.add("rev_qc_children_frame_max_disp_mm", small(mx), mx, f"{QC} :: anatomies[childA, childB, childC].frames."
          "max_corner_displacement_vs_t1_mm (largest)")
    for k, nm in (("translation_norm_mm", "white_translation_max_mm"), ("rotation_max_displacement_mm", "white_rotation_disp_max_mm")):
        mx = max(an[c]["white_vs_t1"]["decomposition"][k] for c in CHILDREN)
        F.add(f"rev_qc_children_{nm}", mm2(mx), mx, f"{QC} :: anatomies[childA, childB, childC].white_vs_t1.decomposition.{k} (largest)")
    mx = max(an[c]["white_vs_t1"]["best_shift_norm_mm"] for c in CHILDREN)
    F.add("rev_qc_children_white_best_shift_max_mm", mm(mx), mx, f"{QC} :: anatomies[childA, childB, childC].white_vs_t1."
          "best_shift_norm_mm (largest)")
    # the 2-year template: its scalp surface lies outside its MRI head boundary
    a, src = an["infant2yr"], f"{QC} :: anatomies['infant2yr']"
    for edge, et in (("otsu", ""), ("half_max", "_halfmax")):
        x = a["scalp_vs_mri"]["offsets"]["cap"][edge]["median_mm"]
        F.add(f"rev_qc_infant2yr_outward{et}_mm", mm(-x), -x, f"{src}.scalp_vs_mri.offsets.cap.{edge}.median_mm, sign reversed (the "
              "surface lies outside the MRI head boundary by this much, cap median, mm)")
        y = a["verdict"]["checks"]["scalp_offset_excess_over_adult_mm"][edge]
        F.add(f"rev_qc_infant2yr_outward_excess{et}_mm", mm(-y), -y, f"{src}.verdict.checks.scalp_offset_excess_over_adult_mm.{edge}, "
              "sign reversed (outward relative to the adult's own offset, mm)")
    # the correction method on the adult (declared displacements recovered from its T1)
    ta, p = q["tests_on_adult"], f"{QC} :: tests_on_adult"
    F.add("rev_qc_test_inward_declared_mm", g(ta["declared_test_values"]["inward_mm"]), ta["declared_test_values"]["inward_mm"],
          f"{p}.declared_test_values.inward_mm (the adult's scalp moved inward)")
    x = ta["inward"]["recovered_displacement_cap_median_mm"]
    F.add("rev_qc_test_inward_recovered_mm", mm(x), x, f"{p}.inward.recovered_displacement_cap_median_mm")
    F.add("rev_qc_test_inward_recovered_mm_2dp", mm2(x), x, f"{p}.inward.recovered_displacement_cap_median_mm")
    x = ta["inward"]["decomposition_cap"]["translation_norm_mm"]
    F.add("rev_qc_test_inward_translation_mm", mm2(x), x, f"{p}.inward.decomposition_cap.translation_norm_mm (spurious translation)")
    x = ta["inward"]["corrected_minus_adult_corrected_cap"]["p50"]
    F.add("rev_qc_test_inward_corrected_vs_adult_mm", smm2(x), x, f"{p}.inward.corrected_minus_adult_corrected_cap.p50 (the "
          "corrected displaced scalp minus the adult's corrected scalp, cap median, mm)")
    tv = ta["declared_test_values"]["translation_mm"]
    F.add("rev_qc_test_translation_declared_mm", g(math.hypot(*tv)), math.hypot(*tv), f"{p}.declared_test_values.translation_mm "
          f"({listing(tv)} mm; its length)")
    tr = ta["translation"]
    F.add("rev_qc_test_translation_recovered_mm", mm2(tr["decomposition_cap"]["translation_norm_mm"]),
          tr["decomposition_cap"]["translation_norm_mm"], f"{p}.translation.decomposition_cap.translation_norm_mm")
    F.add("rev_qc_test_translation_explained_pct", pct(100 * tr["decomposition_cap"]["explained_share"]),
          100 * tr["decomposition_cap"]["explained_share"], f"{p}.translation.decomposition_cap.explained_share")
    F.add("rev_qc_test_untranslated_mm", mm2(tr["decomposition_cap_untranslated"]["translation_norm_mm"]),
          tr["decomposition_cap_untranslated"]["translation_norm_mm"], f"{p}.translation.decomposition_cap_untranslated."
          "translation_norm_mm (the same fit without the displacement)")
    wm = ta["white_translation"]
    F.add("rev_qc_test_white_shift_declared_mm", g(math.hypot(*wm["expected_realigning_shift_mm"])),
          math.hypot(*wm["expected_realigning_shift_mm"]), f"{p}.white_translation.expected_realigning_shift_mm (its length)")
    F.add("rev_qc_test_white_shift_recovered_mm", mm(wm["measured"]["best_shift_norm_mm"]), wm["measured"]["best_shift_norm_mm"],
          f"{p}.white_translation.measured.best_shift_norm_mm")
    F.add("rev_qc_test_white_contrast_gain_pct", pct(100 * wm["measured"]["contrast_gain_share"]),
          100 * wm["measured"]["contrast_gain_share"], f"{p}.white_translation.measured.contrast_gain_share")


# ------------------------------------------------------------------------------------------------
# RESERVED: the facts of the two revision analyses whose results are not committed yet. Write them here when the result
# files land, call them from facts() below, and name them with these prefixes and sources:
#   rev_cgap_  the counterfactual helmet fitted at the adult's gap: scripts/study_g3b_constant_gap.py ->
#              results/g3b_constant_gap/g3b_constant_gap_summary.json
#   rev_conf_  the confirmatory spike run (endpoint, seeds and sample size declared in configs/g4_confirmatory.toml):
#              scripts/g4_confirmatory.py -> results/g4_confirm/g4_confirm_summary.json and g4c_<anatomy>_summary.json
# Until then no fact carries these prefixes, so a report reference to one fails the site build instead of printing a
# number that does not exist yet; tests/test_report_facts_rev.py checks that every such fact cites its own result file.
RESERVED = {"rev_cgap_": "results/g3b_constant_gap/g3b_constant_gap_summary.json",
            "rev_conf_": "results/g4_confirm/g4_confirm_summary.json"}


# ------------------------------------------------------------------------------------------------
def facts(root: Path = ROOT) -> dict:
    """Every revision fact, name -> {"value", "raw", "source"}."""
    root = Path(root)
    F = Facts()
    ns_facts(F, root)
    cov_facts(F, root)
    qc_facts(F, root)
    return dict(F)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grep", default="", help="print only facts whose name contains this text")
    args = ap.parse_args(argv)
    for name, f in facts(ROOT).items():
        if args.grep in name:
            print(f"{name} = {f['value']}    [{f['source']}]")


if __name__ == "__main__":
    main()
