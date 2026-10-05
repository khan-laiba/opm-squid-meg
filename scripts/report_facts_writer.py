#!/usr/bin/env python3
"""Facts added by the writers of the revised report (prefix wr_), merged last by scripts/report_facts.py.

facts(root) -> {name: {"value": text as printed, "raw": unrounded number(s) or text, "source": "<file> :: <key path>"}}.
Read only; nothing is re-run. Every fact is computed from unrounded stored values:
  wr_adult_<array>_<cond>_db      the adult's primary comparison as D = 20 log10(d_OPM / d_Neuromag) in dB, from the stored
                                  median log2 ratio of results/g2/g2_summary.json (the report's one dB convention, applied
                                  to the adult's stored ratios); _ci twins from the stored log2 interval
  wr_opm_cell_peak_effect_pct     the largest effect of the 10-mm OPM averaging volume on the field at its peak, read from the
                                  register row A-OPM-CELL of docs/provenance_register.md (its magnitude, as a percentage)
  wr_db_per_percent, wr_ratio_per_db  constants of the dB convention D = 20 log10(ratio) (docs/methods.md, A-G3-METRIC): the dB
                                  change for a 1 % change of the ratio, and the ratio that corresponds to 1 dB
  wr_locfig_*, wr_jointfig_*      the summary counts stored with the supplementary localization figures
                                  (results/report/figures_supplement.json: the family of paired error comparisons drawn, how
                                  many reach p < 0.05 by inverse method, how many fall below the within-family Bonferroni
                                  threshold, how many medians are zero or positive, how many estimates rest on mostly
                                  undetected events; the exact McNemar tests of joint success and the held-out false-event
                                  rates of the localization study)
  wr_qc_*_offset_*_mm            the MRI quality check's cap-median offsets between the MRI head boundary and the scalp
                                  used, as unsigned magnitudes with the direction in the name (inside: the scalp lies inside
                                  the boundary, a positive stored offset; outside: a negative one), from
                                  results/g3b_children_qc/children_qc.json
  wr_cov_band_mag_emptyroom_from8hz_range  the magnetometers' empty-room prediction (model over measured amplitude) as the
                                  range over the four sub-bands from 8 Hz up (results/g2_covariance_validation/...)
  wr_adult_dense_ib_ratio_area_weighted_cortical_3dp  the adult's dense ratio area-weighted without the medial wall, the
                                  pediatric convention, from the stored dB value of results/g3b/g3b_summary.json link_to_g2
  wr_cgap_adult_fixed_minus_fittedtop_combined_ib(_ci)  the adult's within-head contrast against the helmet fitted at its own
                                  top-contact gap (results/g3b_constant_gap/g3b_constant_gap_summary.json headline)
  wr_seed_*                       seeds the supplementary seed table lacked (the depth-bin bootstrap; the motion study)
  wr_fig_geometry_slab_mm         the half-thickness of the slab of sensors drawn in the geometry figure
  wr_g4_oracle_*                  how many anatomies have the oracle's exploratory S50 ratio above the practical detector's
  wr_confirm_*                    declarations of configs/g4_confirmatory.toml: the mismatch variant's coregistration error
                                  and the pilot runs recorded in the declaration of the location count
  wr_n_*_words, wr_*_words        small counts in words for running prose, each from the file that holds the count
  wr_*_share_pos_pct              shares of cortical area where the OPM is ahead, as percentages
  wr_cov_s4_dense_*               scenario D of the covariance check (Neuromag with its measured covariance, the OPM as
                                  modelled; sensor, brain and room noise) by lobe and by 5-mm depth bin: medians over the
                                  per-target values of results/g2_covariance_validation/covariance_validation_targets.csv
                                  (point estimates, no intervals); the overall median reproduces the stored scenario ratio
  wr_cov_dist_*_label             the distance-bin labels of the covariance check's correlation-against-distance comparison
  wr_cgap_depth_*_label           the depth-stratum labels of the constant-gap run's per-stratum Delta
  wr_cgap_scaled_templates_fitted_delta_ci_bounds_range  the lowest and highest interval bound of Delta at the adult's gap
                                  over the scaled adults and templates
  wr_n_children_below_adult_fitted_words  how many children have a negative Delta at the adult's gap (point estimate)
  wr_p_exploratory_threshold      the conventional uncorrected threshold of the exploratory tests, recovered from the stored
                                  Bonferroni threshold of the localization figure and its family size
  wr_qc_near_lo_mm, wr_qc_near_hi_mm  the two near-scalp thresholds of the children's MRI check
  wr_ns_depth_bin_mm              the depth-bin width of the noise-sensitivity run's own depth results
Formats as scripts/report_facts_g12.py (whose helpers are imported): dB signed with 2 decimals, intervals "[lo, hi]",
U+2212 for negatives, counts with thousands separators.

Usage: .venv/bin/python scripts/report_facts_writer.py [--grep TEXT]
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
G2 = "results/g2/g2_summary.json"
REGISTER = "docs/provenance_register.md"
METHODS = "docs/methods.md"
FIGS_SUPP = "results/report/figures_supplement.json"
QC = "results/g3b_children_qc/children_qc.json"
COVVAL = "results/g2_covariance_validation/covariance_validation.json"
G3B = "results/g3b/g3b_summary.json"
CGAP = "results/g3b_constant_gap/g3b_constant_gap_summary.json"
DEPTH_BINS = "results/g2/g2_depth_bins.json"
MOTION = "results/g4/g4_motion_summary.json"
FIGS_CLEAN = "results/report/figures_clean.json"
G4_PED = "results/g4/g4_pediatric_comparison.json"
G4_SUMMARY = "results/g4/g4_{anatomy}_summary.json"
CONFIRM_CFG = "configs/g4_confirmatory.toml"
HEAD_LABELS = {"adult": "adult", "school": "school-age size", "size2yr": "2-year size", "infant2yr": "24-month template",
               "infant18mo": "18-month template", "infant12mo": "12-month template", "childA": "child A",
               "childB": "child B", "childC": "child C"}
WORDS = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
         10: "ten", 11: "eleven", 12: "twelve"}


def _load_g12():
    spec = importlib.util.spec_from_file_location("report_facts_g12", Path(__file__).with_name("report_facts_g12.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G12 = _load_g12()
MINUS, signed, count = G12.MINUS, G12.signed, G12.count


class Facts(G12.Facts):
    """G12's container, with the prefix and the source form of this module enforced."""

    def add(self, name, value, raw, source):
        if not name.startswith("wr_"):
            raise ValueError(f"fact {name!r} lacks the wr_ prefix")
        if str(value).strip() == "":
            raise ValueError(f"fact {name!r} has an empty value")
        if " :: " not in source:
            raise ValueError(f"fact {name!r}: source lacks ' :: <key path>'")
        super().add(name, value, raw, source)


def db_of_log2(x: float) -> float:
    """20 log10 of a ratio stored as log2: 20 * x * log10(2)."""
    return 20.0 * x * math.log10(2.0)


ARRAYS = {"opm_dense": "dense", "opm_matched": "matched"}
CONDS = {"intrinsic+brain": "ib", "projected": "proj", "intrinsic": "int"}


def adult_db_facts(F: Facts, root: Path) -> None:
    d = json.loads((root / G2).read_text())
    oracle = d["primary"]["oracle"]
    for arr, atok in ARRAYS.items():
        for cond, ctok in CONDS.items():
            key = f"{arr}/combined/{cond}"
            e = oracle[key]
            med, ci = float(e["median_log2"]), [float(c) for c in e["ci95"]]
            src = (f"{G2} :: primary.oracle['{key}'].median_log2 (median over {e['n']:,} targets of log2(d_OPM/d_Neuromag), "
                   f"all 306 Neuromag channels; derived: D = 20 log10(2**median_log2) = 20 log10(2) x median_log2 dB)")
            F.add(f"wr_adult_{atok}_{ctok}_db", signed(db_of_log2(med), 2), db_of_log2(med), src)
            lo, hi = db_of_log2(ci[0]), db_of_log2(ci[1])
            F.add(f"wr_adult_{atok}_{ctok}_db_ci", f"[{signed(lo, 2)}, {signed(hi, 2)}]", [lo, hi],
                  f"{G2} :: primary.oracle['{key}'].ci95 (95 % parcel-bootstrap interval of the median log2 ratio; "
                  f"derived: each end as 20 log10(2) x log2 value, dB)")


_CELL_ROW = re.compile(r"^\|\s*A-OPM-CELL\s*\|.*Effect at the field peak:\s*(-?[0-9.]+)\s*%\s*\((\d+) mm source distance\)",
                       re.M)


def register_facts(F: Facts, root: Path) -> None:
    text = (root / REGISTER).read_text()
    m = _CELL_ROW.search(text)
    if m is None:
        raise ValueError(f"{REGISTER}: row A-OPM-CELL with 'Effect at the field peak' not found")
    effect, dist = float(m.group(1)), int(m.group(2))
    F.add("wr_opm_cell_peak_effect_pct", f"{abs(effect):.1f}%", abs(effect),
          f"{REGISTER} :: row A-OPM-CELL, 'Effect at the field peak' at {dist} mm source distance (magnitude of the "
          f"largest, nearest-source value; the register lists smaller effects at larger distances)")


def convention_facts(F: Facts, root: Path) -> None:
    """Two constants of the dB convention D = 20 log10(ratio), for the estimand map of the supplementary text."""
    if "20 log10" not in (root / METHODS).read_text():
        raise ValueError(f"{METHODS}: the dB convention '20 log10' is not stated")
    src = f"{METHODS} :: section 10, G3B metric (A-G3-METRIC): D = 20 log10 of the detectability ratio; derived: "
    per_pct = 20.0 * math.log10(1.01)
    F.add("wr_db_per_percent", f"{per_pct:.3f}", per_pct, src + "20 log10(1.01), the dB change for a 1 % change of the ratio")
    per_db = 10.0 ** (1.0 / 20.0)
    F.add("wr_ratio_per_db", f"{per_db:.3f}", per_db, src + "10**(1/20), the ratio that corresponds to 1 dB")


def figure_summary_facts(F: Facts, root: Path) -> None:
    """The summary counts stored with the supplementary localization figures (results/report/figures_supplement.json)."""
    d = json.loads((root / FIGS_SUPP).read_text())["figures"]
    loc = d["Figure_S_localization_effects"]["values"]
    src = f"{FIGS_SUPP} :: figures.Figure_S_localization_effects.values."
    F.add("wr_locfig_n_error_comparisons", count(loc["n_estimates"]), loc["n_estimates"],
          src + "n_estimates (paired error comparisons drawn: 9 anatomies x 2 OPM arrays x 4 conditions x 3 inverse methods)")
    F.add("wr_locfig_n_p_below_005", count(loc["n_p_below_005"]), loc["n_p_below_005"],
          src + "n_p_below_005 (Wilcoxon signed-rank, uncorrected)")
    for method, tok in (("ecd", "ecd"), ("dspm_mne", "dspm_mne"), ("dspm", "dspm")):
        n = loc["n_p_below_005_by_method"][method]
        F.add(f"wr_locfig_n_p_below_005_{tok}", count(n), n, src + f"n_p_below_005_by_method['{method}']")
    F.add("wr_locfig_n_p_below_bonferroni", count(loc["n_p_below_bonferroni"]), loc["n_p_below_bonferroni"],
          src + "n_p_below_bonferroni (p below 0.05 / 20, the comparisons of one anatomy and OPM array)")
    thr = float(loc["bonferroni_threshold"])
    F.add("wr_locfig_bonferroni_threshold", f"{thr:.4f}", thr, src + "bonferroni_threshold (0.05 / family_per_anatomy_and_array)")
    F.add("wr_locfig_family_per_anatomy_and_array", count(loc["family_per_anatomy_and_array"]),
          loc["family_per_anatomy_and_array"], src + "family_per_anatomy_and_array")
    F.add("wr_locfig_n_zero_median", count(loc["n_zero_median"]), loc["n_zero_median"],
          src + "n_zero_median (median paired difference exactly zero: tied dSPM peaks on the shared grid)")
    F.add("wr_locfig_n_positive_median", count(loc["n_positive_median"]), loc["n_positive_median"],
          src + "n_positive_median (median paired difference above zero: Neuromag closer)")
    for strength in ("80nAm", "320nAm"):
        n = loc["n_open_by_strength"][strength]
        F.add(f"wr_locfig_n_open_{strength.lower()}", count(n), n,
              src + f"n_open_by_strength['{strength}'] (estimates with fewer than half of the events detected by both systems)")
    joint = d["Figure_S_joint_detection_localization"]["values"]
    src = f"{FIGS_SUPP} :: figures.Figure_S_joint_detection_localization.values."
    n_ecd, n_dspm = joint["n_mcnemar"]["ecd"], joint["n_mcnemar"]["dspm"]
    if n_ecd != n_dspm:
        raise ValueError(f"{FIGS_SUPP}: the McNemar families of the dipole fit and the study's dSPM differ ({n_ecd}, {n_dspm})")
    F.add("wr_jointfig_n_mcnemar_per_method", count(n_ecd), n_ecd,
          src + "n_mcnemar['ecd'] (= n_mcnemar['dspm']: paired joint-success comparisons per inverse method)")
    for method in ("ecd", "dspm"):
        n = joint["n_mcnemar_below_005"][method]
        F.add(f"wr_jointfig_n_mcnemar_below_005_{method}", count(n), n,
              src + f"n_mcnemar_below_005['{method}'] (exact McNemar, uncorrected)")
    lo, hi = (float(x) for x in joint["heldout_false_per_min_range"])
    F.add("wr_jointfig_heldout_false_per_min_range", f"{lo:.1f} to {hi:.1f}", [lo, hi],
          src + "heldout_false_per_min_range (false events per minute on the localization study's held-out null data, "
          "over the nine anatomies and three systems)")


def words(n: int) -> str:
    """A small count in words, for running prose."""
    n = int(n)
    if n not in WORDS:
        raise ValueError(f"no word form for {n}")
    return WORDS[n]


def _name_list(names):
    """'a, b and c'."""
    names = list(names)
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def qc_offset_facts(F: Facts, root: Path) -> None:
    """Unsigned cap-median offsets of the children's MRI check, with the direction in the fact name."""
    d = json.loads((root / QC).read_text())["anatomies"]
    src = f"{QC} :: anatomies[{{}}].scalp_vs_mri.offsets.cap.{{}}.median_mm (MRI head boundary minus the scalp used, along " \
          "its outward normal; derived: the magnitude, the direction being in the fact name)"
    children = ("childA", "childB", "childC")
    for edge, tok in (("otsu", "otsu"), ("half_max", "halfmax")):
        vals = [float(d[c]["scalp_vs_mri"]["offsets"]["cap"][edge]["median_mm"]) for c in children]
        if min(vals) <= 0:
            raise ValueError(f"{QC}: a child's {edge} offset is not positive (scalp inside the boundary)")
        lo, hi = min(vals), max(vals)
        F.add(f"wr_qc_children_offset_inside_{tok}_mm_range", f"{lo:.1f} to {hi:.1f}", [lo, hi],
              src.format("childA, childB, childC", edge))
        a = float(d["adult"]["scalp_vs_mri"]["offsets"]["cap"][edge]["median_mm"])
        if a <= 0:
            raise ValueError(f"{QC}: the adult's {edge} offset is not positive")
        F.add(f"wr_qc_adult_offset_inside_{tok}_mm", f"{a:.1f}", a, src.format("'adult'", edge))
    t = float(d["infant2yr"]["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"])
    if t >= 0:
        raise ValueError(f"{QC}: the 24-month template's otsu offset is not negative (scalp outside the boundary)")
    F.add("wr_qc_infant2yr_offset_outside_otsu_mm", f"{abs(t):.1f}", abs(t), src.format("'infant2yr'", "otsu"))


def covariance_band_facts(F: Facts, root: Path) -> None:
    """The magnetometers' empty-room prediction over the sub-bands from 8 Hz up, and the number of sub-bands."""
    sb = json.loads((root / COVVAL).read_text())["sub_bands"]
    bands = ["8-13Hz", "13-20Hz", "20-30Hz", "30-40Hz"]
    vals = [float(sb[b]["mag"]["empty_room_amp_ratio_model_over_measured"]) for b in bands]
    lo, hi = min(vals), max(vals)
    F.add("wr_cov_band_mag_emptyroom_from8hz_range", f"{lo:.2f} to {hi:.2f}", [lo, hi],
          f"{COVVAL} :: sub_bands['8-13Hz', '13-20Hz', '20-30Hz', '30-40Hz'].mag.empty_room_amp_ratio_model_over_measured "
          "(derived: min and max over the four sub-bands from 8 Hz up)")
    n = len([k for k in sb if k.endswith("Hz") and k != "1-40Hz"])
    F.add("wr_n_subbands", count(n), n, f"{COVVAL} :: sub_bands (count of the sub-bands, the whole 1-40 Hz band excluded)")
    F.add("wr_n_subbands_words", words(n), n, f"{COVVAL} :: sub_bands (count of the sub-bands, in words)")


def link_facts(F: Facts, root: Path) -> None:
    """The adult's dense ratio under the pediatric convention (area-weighted, medial wall excluded)."""
    link = json.loads((root / G3B).read_text())["link_to_g2"]
    db = float(link["adult_centred_area_weighted_cortical_dB"])
    r = 10.0 ** (db / 20.0)
    F.add("wr_adult_dense_ib_ratio_area_weighted_cortical_3dp", f"{r:.3f}", r,
          f"{G3B} :: link_to_g2.adult_centred_area_weighted_cortical_dB (derived: 10**(D/20), the ratio under the "
          "pediatric convention at the adult's measured head position)")


def cgap_facts(F: Facts, root: Path) -> None:
    """The adult's within-head contrast against the helmet fitted at its own top-contact gap."""
    w = json.loads((root / CGAP).read_text())["headline"]["by_anatomy"]["adult"]["within_head"]["gap_matched_top"]
    med, ci = float(w["median"]), [float(c) for c in w["ci95"]]
    src = f"{CGAP} :: headline.by_anatomy.adult.within_head.gap_matched_top (D at top contact in the fixed helmet minus " \
          "D in the helmet fitted at the adult's top-contact gap, dense OPM array against Neuromag's 306 channels, sensor " \
          "plus brain noise, dB)"
    F.add("wr_cgap_adult_fixed_minus_fittedtop_combined_ib", signed(med, 2), med, src + ".median")
    F.add("wr_cgap_adult_fixed_minus_fittedtop_combined_ib_ci", f"[{signed(ci[0], 2)}, {signed(ci[1], 2)}]", ci,
          src + ".ci95 (95 % parcel-bootstrap interval)")


def seed_facts(F: Facts, root: Path) -> None:
    """Seeds that the supplementary seed table lacked."""
    s = json.loads((root / DEPTH_BINS).read_text())["method"]["seed"]
    F.add("wr_seed_g2_depth_bins", str(int(s)), int(s), f"{DEPTH_BINS} :: method.seed (parcel bootstrap of the depth bins)")
    cfg = json.loads((root / MOTION).read_text())["config"]
    for part in ("coupling", "timecourse"):
        s = cfg[part]["seed"]
        F.add(f"wr_seed_motion_{part}", str(int(s)), int(s), f"{MOTION} :: config.{part}.seed")


_SLAB = re.compile(r"within (\d+) mm of the plane")
_QC_SECTION = re.compile(r"each (\d+) mm wide")
FIGS_QC = "results/report/figures_qc.json"
NOISE_SENS = "results/g2_noise_sensitivity/noise_sensitivity_summary.json"


def figure_constant_facts(F: Facts, root: Path) -> None:
    """Plot constants of the report figures that their captions state."""
    desc = json.loads((root / FIGS_CLEAN).read_text())["figures"]["Figure_R12_geometry"]["description"]
    m = _SLAB.search(desc)
    if m is None:
        raise ValueError(f"{FIGS_CLEAN}: the geometry figure's slab half-thickness ('within N mm of the plane') not found")
    F.add("wr_fig_geometry_slab_mm", m.group(1), int(m.group(1)),
          f"{FIGS_CLEAN} :: figures.Figure_R12_geometry.description ('sensors within {m.group(1)} mm of the plane')")
    fig = json.loads((root / FIGS_QC).read_text())["figures"]["Figure_S_children_qc"]
    for field in ("alt", "description", "caption_draft"):
        m = _QC_SECTION.search(str(fig.get(field, "")))
        if m is not None:
            break
    else:
        raise ValueError(f"{FIGS_QC}: the MRI check figure's section width ('each N mm wide') not found")
    F.add("wr_fig_qc_section_mm", m.group(1), int(m.group(1)),
          f"{FIGS_QC} :: figures.Figure_S_children_qc.{field} ('MRI sections, each {m.group(1)} mm wide')")


def far_field_facts(F: Facts, root: Path) -> None:
    """The depths of the cardiac current dipoles of the far-field sensitivity analysis."""
    heart = json.loads((root / NOISE_SENS).read_text())["far_field"]["heart_head_mm"]
    depths = sorted(abs(float(p[2])) for p in heart)
    if len(depths) != 3 or any(d <= 0 for d in depths):
        raise ValueError(f"{NOISE_SENS}: expected three cardiac dipole positions below the head, found {heart}")
    F.add("wr_ns_ff_heart_depths_mm", _name_list(f"{d:g}" for d in depths), depths,
          f"{NOISE_SENS} :: far_field.heart_head_mm (the cardiac current dipoles' head-frame positions; derived: their depths "
          "below the head origin, mm)")
    F.add("wr_ns_ff_heart_mid_depth_mm", f"{depths[1]:g}", depths[1],
          f"{NOISE_SENS} :: far_field.heart_head_mm (derived: the middle of the three depths, the one the magnetic dipole shares)")


def oracle_facts(F: Facts, root: Path) -> None:
    """How many anatomies have the oracle's exploratory S50 ratio (10-20 mm, dense vs Neuromag 306) above the practical's."""
    labels = json.loads((root / G4_PED).read_text())["labels"]
    above, below = [], []
    for a in labels:
        p = json.loads((root / G4_SUMMARY.format(anatomy=a)).read_text())["paired"]
        o = float(p["opm_dense/opm_vs_squid/combined/oracle"]["depth0"]["s50_ratio_squid_over_opm"]["value"])
        q = float(p["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]["s50_ratio_squid_over_opm"]["value"])
        (above if o > q else below).append(a)
    src = (f"{G4_PED} :: labels (the {len(labels)} anatomies), each read from {G4_SUMMARY.format(anatomy='<anatomy>')} "
           "paired['opm_dense/opm_vs_squid/combined/oracle'].depth0.s50_ratio_squid_over_opm.value against "
           "paired['opm_dense/opm_vs_squid/combined/practical@1'].depth0.s50_ratio_squid_over_opm.value (derived)")
    F.add("wr_g4_oracle_above_practical_n", count(len(above)), len(above), src + ": count with oracle > practical")
    F.add("wr_g4_oracle_below_practical_heads", _name_list(HEAD_LABELS[a] for a in below), below,
          src + ": the anatomies with oracle <= practical")


_PILOT = re.compile(r"after (\w+) pilot runs on test seeds \(child A, (\d+) locations: ([\d/]+) and ([\d/]+)\s+#?\s*"
                    r"locations, p = ([0-9.]+) and ([0-9.]+)\)", re.S)


def confirmatory_facts(F: Facts, root: Path) -> None:
    """Declarations of the confirmatory configuration: the mismatch variant's coregistration error and the pilot runs."""
    text = (root / CONFIRM_CFG).read_text()
    cfg = tomllib.loads(text)
    variants = cfg["variant"]
    if len(variants) != 1:
        raise ValueError(f"{CONFIRM_CFG}: expected one declared variant, found {len(variants)}")
    v = variants[0]
    F.add("wr_confirm_coreg_shift_mm", f"{float(v['coreg_shift_mm']):g}", float(v["coreg_shift_mm"]),
          f"{CONFIRM_CFG} :: [[variant]] coreg_shift_mm (the mismatch variant's coregistration error, shared by all arrays)")
    F.add("wr_confirm_coreg_angle_deg", f"{float(v['coreg_angle_deg']):g}", float(v["coreg_angle_deg"]),
          f"{CONFIRM_CFG} :: [[variant]] coreg_angle_deg")
    n_rep = int(cfg["design"]["noise_replicates"])
    F.add("wr_confirm_replicates_words", words(n_rep), n_rep, f"{CONFIRM_CFG} :: [design] noise_replicates (in words)")
    comment = re.sub(r"\n\s*#\s*", " ", text)
    m = _PILOT.search(comment)
    if m is None:
        raise ValueError(f"{CONFIRM_CFG}: the pilot-run record in the comment on locations_per_stratum not found")
    n_words, n_locs, locs1, locs2, p1, p2 = m.groups()
    src = f"{CONFIRM_CFG} :: [design] locations_per_stratum, the declaration's comment (pilot runs on test seeds, child A)"
    F.add("wr_confirm_pilot_runs_words", n_words, n_words, src + ": number of pilot runs")
    F.add("wr_confirm_pilot_locs", f"{locs1} and {locs2}", [locs1, locs2],
          src + ": locations favouring the OPM / Neuromag in each pilot run")
    F.add("wr_confirm_pilot_p", f"{p1} and {p2}", [float(p1), float(p2)], src + ": sign-flip p of each pilot run")
    F.add("wr_confirm_pilot_locs_per_band", n_locs, int(n_locs), src + ": locations per band in the pilot runs")


def count_word_facts(F: Facts, root: Path) -> None:
    """Small counts in words, for running prose."""
    labels = json.loads((root / G4_PED).read_text())["labels"]
    smaller = [a for a in labels if a != "adult"]
    F.add("wr_n_smaller_heads_words", words(len(smaller)), len(smaller), f"{G4_PED} :: labels (count without the adult, in words)")
    F.add("wr_n_anatomies_words", words(len(labels)), len(labels), f"{G4_PED} :: labels (count, in words)")
    children = [a for a in labels if a.startswith("child")]
    templates = [a for a in labels if a.startswith("infant")]
    scaled = [a for a in labels if a in ("school", "size2yr")]
    F.add("wr_n_children_words", words(len(children)), len(children), f"{G4_PED} :: labels (the school-aged children, in words)")
    F.add("wr_n_templates_words", words(len(templates)), len(templates), f"{G4_PED} :: labels (the infant templates, in words)")
    F.add("wr_n_scaled_words", words(len(scaled)), len(scaled), f"{G4_PED} :: labels (the scaled adults, in words)")
    stretches = json.loads((root / G4_SUMMARY.format(anatomy="adult")).read_text())["config"]["events"]["stretches"]
    F.add("wr_n_morphologies_words", words(len(stretches)), len(stretches),
          f"{G4_SUMMARY.format(anatomy='adult')} :: config.events.stretches (count, in words)")


def share_facts(F: Facts, root: Path) -> None:
    """Shares of cortical area with D > 0 (fixed adult helmet, top contact, sensor plus brain noise), as percentages."""
    comp = json.loads((root / G3B).read_text())["comparisons"]
    dense = float(comp["school/opm_dense/combined/intrinsic+brain/detect"]["d_adult"]["share_positive"]) * 100
    F.add("wr_adult_dense_ib_share_pos_pct", f"{dense:.1f}%", dense,
          f"{G3B} :: comparisons['school/opm_dense/combined/intrinsic+brain/detect'].d_adult.share_positive (the adult at "
          "top contact; the same in every comparison; as a percentage)")
    matched = float(comp["school/opm_matched/combined/intrinsic+brain/detect"]["d_adult"]["share_positive"]) * 100
    F.add("wr_adult_matched_ib_share_pos_pct", f"{matched:.0f}%", matched,
          f"{G3B} :: comparisons['school/opm_matched/combined/intrinsic+brain/detect'].d_adult.share_positive (as a percentage)")
    heads = [a for a in HEAD_LABELS if a != "adult"]
    vals = [float(comp[f"{a}/opm_matched/combined/intrinsic+brain/detect"]["d_child"]["share_positive"]) * 100 for a in heads]
    lo, hi = min(vals), max(vals)
    F.add("wr_smaller_matched_ib_share_pos_pct_range", f"{lo:.0f}% to {hi:.0f}%", [lo, hi],
          f"{G3B} :: comparisons['<head>/opm_matched/combined/intrinsic+brain/detect'].d_child.share_positive (derived: min "
          "and max over the eight smaller heads, as percentages)")


COVVAL_TARGETS = "results/g2_covariance_validation/covariance_validation_targets.csv"
S4_NUM = "detect_opm_dense_independent_intrinsic_brain_env"
S4_DEN = "detect_neuromag_combined_measured_corrected"
S4_KEY = "opm_dense/combined/S4_neuromag_measured_opm_as_modelled"
LOBES = ("frontal", "parietal", "temporal", "occipital", "cingulate", "insula")
S4_DEPTH_EDGES = [float(x) for x in range(10, 70, 5)]
S4_MIN_N = 10


def _read_targets_csv(path: Path) -> list:
    """The per-target table, whose first line is a provenance comment."""
    import csv
    with path.open() as f:
        first = f.readline()
        if not first.startswith("#"):
            f.seek(0)
        return list(csv.DictReader(f))


def _median(xs) -> float:
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def scenario_d_facts(F: Facts, root: Path) -> None:
    """Scenario D (Neuromag measured, finite-sample corrected; the OPM as modelled) by lobe and by depth bin, as medians
    over the per-target values of the covariance check; the overall median must reproduce the stored scenario ratio."""
    rows = _read_targets_csv(root / COVVAL_TARGETS)
    stored = json.loads((root / COVVAL).read_text())["opm_implication"]["ratios"][S4_KEY]
    ratios = [float(r[S4_NUM]) / float(r[S4_DEN]) for r in rows]
    overall = _median(ratios)
    if len(rows) != int(stored["n"]) or abs(overall - float(stored["ratio"])) > 5e-4:
        raise ValueError(f"{COVVAL_TARGETS}: the per-target scenario-D median ({overall:.4f} over {len(rows)} targets) does "
                         f"not reproduce the stored ratio {stored['ratio']:.4f} (n = {stored['n']})")
    src = (f"{COVVAL_TARGETS} :: median over targets of {S4_NUM} / {S4_DEN} (scenario D: Neuromag's measured covariance, "
           "finite-sample corrected, 305 good channels; the dense OPM array as modelled; sensor, brain and room noise), "
           "derived from the per-target values; a point estimate without interval; the median over all targets reproduces "
           f"{COVVAL} :: opm_implication.ratios['{S4_KEY}'].ratio")
    F.add("wr_cov_s4_dense_ratio_3dp", f"{overall:.3f}", overall, src + " (all targets)")
    for lobe in LOBES:
        vals = [x for r, x in zip(rows, ratios) if r["lobe"] == lobe]
        if len(vals) < S4_MIN_N:
            raise ValueError(f"{COVVAL_TARGETS}: lobe {lobe!r} has {len(vals)} targets")
        m = _median(vals)
        F.add(f"wr_cov_s4_dense_lobe_{lobe}_ratio", f"{m:.2f}", m, src + f" (lobe '{lobe}', {len(vals)} targets)")
    first_below = None
    for lo in S4_DEPTH_EDGES:
        hi = lo + 5.0
        vals = [x for r, x in zip(rows, ratios) if lo <= float(r["depth_mm"]) < hi]
        if len(vals) < S4_MIN_N:
            continue
        m = _median(vals)
        tok = f"{lo:.0f}_{hi:.0f}"
        F.add(f"wr_cov_s4_dense_depth_{tok}_ratio", f"{m:.2f}", m,
              src + f" (targets {lo:.0f} <= depth_mm < {hi:.0f}, {len(vals)} targets)")
        if first_below is None and m < 1.0:
            first_below = (lo, hi)
    if first_below is None:
        raise ValueError(f"{COVVAL_TARGETS}: the scenario-D dense median is at or above 1 in every depth bin")
    F.add("wr_cov_s4_dense_depth_first_below1_label", f"{first_below[0]:.0f}\u2013{first_below[1]:.0f}", list(first_below),
          src + " (the shallowest 5-mm depth bin, from 10 mm, whose median lies below 1; label in mm)")


def distance_label_facts(F: Facts, root: Path) -> None:
    """The distance-bin labels of the covariance check's correlation-against-distance comparison."""
    bins = json.loads((root / COVVAL).read_text())["structure"]["brain"]["mag"]["corr_vs_distance"]
    wanted = {"25-50", "50-75", "100-125", "200-225"}
    found = set()
    for i, b in enumerate(bins):
        label = str(b.get("bin", ""))
        if label in wanted:
            lo, hi = label.split("-")
            F.add(f"wr_cov_dist_d{lo}_{hi}_label", f"{lo}\u2013{hi}", [float(lo), float(hi)],
                  f"{COVVAL} :: structure.brain.mag.corr_vs_distance[{i}].bin (inter-sensor distance bin, mm)")
            found.add(label)
    if found != wanted:
        raise ValueError(f"{COVVAL}: distance bins {sorted(wanted - found)} not found")


def cgap_depth_label_facts(F: Facts, root: Path) -> None:
    """The depth-stratum labels of the constant-gap run's per-stratum Delta (child A, the strata the text names)."""
    key = "childA/gap_matched/opm_dense/combined/intrinsic+brain"
    strata = json.loads((root / CGAP).read_text())["delta_same_rule"][key]["delta_by_depth"]
    wanted = {(10.0, 15.0), (20.0, 25.0), (30.0, 40.0), (40.0, 50.0)}
    found = set()
    for i, b in enumerate(strata):
        lo, hi = float(b["lo"]), float(b["hi"])
        if (lo, hi) in wanted:
            F.add(f"wr_cgap_depth_{lo:.0f}_{hi:.0f}_label", f"{lo:.0f}\u2013{hi:.0f}", [lo, hi],
                  f"{CGAP} :: delta_same_rule['{key}'].delta_by_depth[{i}].lo, .hi (depth stratum below the scalp, mm)")
            found.add((lo, hi))
    if found != wanted:
        raise ValueError(f"{CGAP}: depth strata {sorted(wanted - found)} not found")


def cgap_interval_facts(F: Facts, root: Path) -> None:
    """The interval bounds of Delta at the adult's gap over the scaled adults and templates, and how many children lie
    below the adult there (point estimates)."""
    d = json.loads((root / CGAP).read_text())["delta_same_rule"]
    heads = ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")
    cis = [[float(c) for c in d[f"{h}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["ci95"]] for h in heads]
    lo, hi = min(c[0] for c in cis), max(c[1] for c in cis)
    F.add("wr_cgap_scaled_templates_fitted_delta_ci_bounds_range", f"{signed(lo, 2)} to {signed(hi, 2)}", [lo, hi],
          f"{CGAP} :: delta_same_rule['<head>/gap_matched/opm_dense/combined/intrinsic+brain'].delta.ci95 (derived: the "
          "lowest lower bound and the highest upper bound over the scaled adults and the infant templates, dB)")
    children = ("childA", "childB", "childC")
    below = [c for c in children
             if float(d[f"{c}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["median"]) < 0]
    F.add("wr_n_children_below_adult_fitted_words", words(len(below)), len(below),
          f"{CGAP} :: delta_same_rule['child<X>/gap_matched/opm_dense/combined/intrinsic+brain'].delta.median (derived: "
          f"the children with a negative median, {_name_list(HEAD_LABELS[c] for c in below)}; in words)")


def threshold_facts(F: Facts, root: Path) -> None:
    """The conventional uncorrected threshold of the exploratory tests, the MRI check's near-scalp thresholds, the
    rounding shortfall of the check's count for child B, and the noise-sensitivity run's depth-bin width."""
    loc = json.loads((root / FIGS_SUPP).read_text())["figures"]["Figure_S_localization_effects"]["values"]
    p = float(loc["bonferroni_threshold"]) * int(loc["family_per_anatomy_and_array"])
    if abs(p - 0.05) > 1e-9:
        raise ValueError(f"{FIGS_SUPP}: bonferroni_threshold x family_per_anatomy_and_array is {p}, not 0.05")
    F.add("wr_p_exploratory_threshold", f"{p:.2f}", p,
          f"{FIGS_SUPP} :: figures.Figure_S_localization_effects.values.bonferroni_threshold x family_per_anatomy_and_array "
          "(derived: the uncorrected threshold the family-wise threshold was built from)")
    near = [float(x) for x in json.loads((root / QC).read_text())["parameters"]["near_mm"]]
    if len(near) != 2 or near[0] >= near[1]:
        raise ValueError(f"{QC}: parameters.near_mm should hold two increasing thresholds, found {near}")
    F.add("wr_qc_near_lo_mm", f"{near[0]:g}", near[0], f"{QC} :: parameters.near_mm[0] (the lower near-scalp threshold, mm)")
    F.add("wr_qc_near_hi_mm", f"{near[1]:g}", near[1], f"{QC} :: parameters.near_mm[1] (the upper near-scalp threshold, mm)")
    entries = json.loads((root / NOISE_SENS).read_text())["sweep"]["entries"]
    key = next(k for k in entries if float(k) == 15.0)
    bins = entries[key]["depth"]["opm_dense/combined/intrinsic+brain"]
    widths = {float(b["hi"]) - float(b["lo"]) for b in bins}
    if len(widths) != 1:
        raise ValueError(f"{NOISE_SENS}: the depth bins of the sweep are not of one width ({sorted(widths)})")
    w = widths.pop()
    F.add("wr_ns_depth_bin_mm", f"{w:g}", w,
          f"{NOISE_SENS} :: sweep.entries['{key}'].depth['opm_dense/combined/intrinsic+brain'][*].lo, .hi (bin width, mm)")


def facts(root: Path = ROOT) -> dict:
    """Every writer fact, name -> {"value", "raw", "source"}."""
    root = Path(root)
    F = Facts()
    adult_db_facts(F, root)
    register_facts(F, root)
    convention_facts(F, root)
    figure_summary_facts(F, root)
    qc_offset_facts(F, root)
    covariance_band_facts(F, root)
    link_facts(F, root)
    cgap_facts(F, root)
    seed_facts(F, root)
    figure_constant_facts(F, root)
    far_field_facts(F, root)
    oracle_facts(F, root)
    confirmatory_facts(F, root)
    count_word_facts(F, root)
    share_facts(F, root)
    scenario_d_facts(F, root)
    distance_label_facts(F, root)
    cgap_depth_label_facts(F, root)
    cgap_interval_facts(F, root)
    threshold_facts(F, root)
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
