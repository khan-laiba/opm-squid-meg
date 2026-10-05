"""Facts added by the writers of the revised report (scripts/report_facts_writer.py): the prefix, the source form, the
existence of every source file, and every value recomputed from the stored result files and the register."""
import csv
import importlib.util
import json
import math
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
G2 = "results/g2/g2_summary.json"
REGISTER = "docs/provenance_register.md"
METHODS = "docs/methods.md"
FIGS_SUPP = "results/report/figures_supplement.json"
HAVE = all((ROOT / p).exists() for p in (G2, REGISTER, METHODS, FIGS_SUPP))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


W = _load("report_facts_writer")


class TestContainer(unittest.TestCase):
    def test_prefix_and_source_form_are_enforced(self):
        F = W.Facts()
        with self.assertRaises(ValueError):
            F.add("rev_x", "1", 1, f"{G2} :: x")
        with self.assertRaises(ValueError):
            F.add("wr_x", "", 1, f"{G2} :: x")
        with self.assertRaises(ValueError):
            F.add("wr_x", "1", 1, G2)
        F.add("wr_x", "1", 1, f"{G2} :: x")
        with self.assertRaises(ValueError):
            F.add("wr_x", "1", 1, f"{G2} :: x")

    def test_db_conversion(self):
        self.assertAlmostEqual(W.db_of_log2(1.0), 20 * math.log10(2), places=12)
        self.assertAlmostEqual(W.db_of_log2(0.0), 0.0, places=12)


@unittest.skipUnless(HAVE, "result files not present")
class TestWriterFacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.F = W.facts(ROOT)
        cls.g2 = json.loads((ROOT / G2).read_text())

    def test_every_fact_has_the_prefix_a_value_and_an_existing_source(self):
        self.assertGreaterEqual(len(self.F), 31)
        for name, f in self.F.items():
            self.assertRegex(name, r"^wr_[a-z0-9_]+$")
            self.assertEqual(set(f), {"value", "raw", "source"}, name)
            self.assertTrue(str(f["value"]).strip(), name)
            src = f["source"].split(" :: ")[0]
            self.assertTrue((ROOT / src).exists(), f"{name}: source {src} missing")

    def test_adult_db_values_match_the_stored_log2_medians(self):
        oracle = self.g2["primary"]["oracle"]
        for arr, atok in (("opm_dense", "dense"), ("opm_matched", "matched")):
            for cond, ctok in (("intrinsic+brain", "ib"), ("projected", "proj"), ("intrinsic", "int")):
                e = oracle[f"{arr}/combined/{cond}"]
                db = 20 * math.log10(2 ** e["median_log2"])
                f = self.F[f"wr_adult_{atok}_{ctok}_db"]
                self.assertAlmostEqual(f["raw"], db, places=9)
                self.assertEqual(f["value"], f"{db:+.2f}".replace("-", "−"))
                lo, hi = (20 * math.log10(2 ** c) for c in e["ci95"])
                ci = self.F[f"wr_adult_{atok}_{ctok}_db_ci"]
                self.assertAlmostEqual(ci["raw"][0], lo, places=9)
                self.assertAlmostEqual(ci["raw"][1], hi, places=9)
                self.assertEqual(ci["value"], f"[{lo:+.2f}, {hi:+.2f}]".replace("-", "−"))

    def test_the_headline_values_as_printed(self):
        # the adult's primary dense ratio 1.144 is +1.17 dB; the matched 1.009 is +0.08 dB
        self.assertEqual(self.F["wr_adult_dense_ib_db"]["value"], "+1.17")
        self.assertEqual(self.F["wr_adult_matched_ib_db"]["value"], "+0.08")

    def test_cell_effect_matches_the_register_row(self):
        text = (ROOT / REGISTER).read_text()
        m = re.search(r"^\|\s*A-OPM-CELL\s*\|.*Effect at the field peak:\s*(-?[0-9.]+)\s*%", text, re.M)
        self.assertIsNotNone(m)
        f = self.F["wr_opm_cell_peak_effect_pct"]
        self.assertAlmostEqual(f["raw"], abs(float(m.group(1))))
        self.assertEqual(f["value"], f"{abs(float(m.group(1))):.1f}%")

    def test_db_convention_constants(self):
        per_pct, per_db = self.F["wr_db_per_percent"], self.F["wr_ratio_per_db"]
        self.assertAlmostEqual(per_pct["raw"], 20 * math.log10(1.01), places=12)
        self.assertEqual(per_pct["value"], "0.086")
        self.assertAlmostEqual(per_db["raw"], 10 ** 0.05, places=12)
        self.assertEqual(per_db["value"], "1.122")
        for f in (per_pct, per_db):
            self.assertTrue(f["source"].startswith(f"{METHODS} :: "))

    def test_figure_summary_counts_match_the_stored_figure_values(self):
        figs = json.loads((ROOT / FIGS_SUPP).read_text())["figures"]
        loc = figs["Figure_S_localization_effects"]["values"]
        joint = figs["Figure_S_joint_detection_localization"]["values"]
        self.assertEqual(self.F["wr_locfig_n_error_comparisons"]["raw"], loc["n_estimates"])
        self.assertEqual(self.F["wr_locfig_n_error_comparisons"]["value"], f"{loc['n_estimates']:,}")
        self.assertEqual(self.F["wr_locfig_n_p_below_005"]["raw"], loc["n_p_below_005"])
        for method in ("ecd", "dspm_mne", "dspm"):
            self.assertEqual(self.F[f"wr_locfig_n_p_below_005_{method}"]["raw"], loc["n_p_below_005_by_method"][method])
        self.assertEqual(sum(loc["n_p_below_005_by_method"].values()), loc["n_p_below_005"])
        self.assertEqual(self.F["wr_locfig_n_p_below_bonferroni"]["raw"], loc["n_p_below_bonferroni"])
        self.assertEqual(self.F["wr_locfig_bonferroni_threshold"]["value"], f"{loc['bonferroni_threshold']:.4f}")
        self.assertEqual(self.F["wr_locfig_family_per_anatomy_and_array"]["raw"], loc["family_per_anatomy_and_array"])
        self.assertEqual(self.F["wr_locfig_n_zero_median"]["raw"], loc["n_zero_median"])
        self.assertEqual(self.F["wr_locfig_n_positive_median"]["raw"], loc["n_positive_median"])
        self.assertEqual(self.F["wr_locfig_n_open_80nam"]["raw"], loc["n_open_by_strength"]["80nAm"])
        self.assertEqual(self.F["wr_locfig_n_open_320nam"]["raw"], loc["n_open_by_strength"]["320nAm"])
        self.assertEqual(self.F["wr_jointfig_n_mcnemar_per_method"]["raw"], joint["n_mcnemar"]["ecd"])
        for method in ("ecd", "dspm"):
            self.assertEqual(self.F[f"wr_jointfig_n_mcnemar_below_005_{method}"]["raw"], joint["n_mcnemar_below_005"][method])
        lo, hi = joint["heldout_false_per_min_range"]
        self.assertEqual(self.F["wr_jointfig_heldout_false_per_min_range"]["value"], f"{lo:.1f} to {hi:.1f}")
        # the figure's error-comparison family and the text's localization family agree: 216 error comparisons plus
        # 72 joint-success comparisons for each of the two methods that store one make the family of 360
        self.assertEqual(loc["n_estimates"] + 2 * joint["n_mcnemar"]["ecd"], 360)

    def test_qc_offsets_are_unsigned_magnitudes_of_the_stored_cap_medians(self):
        qc = json.loads((ROOT / "results/g3b_children_qc/children_qc.json").read_text())["anatomies"]
        cap = lambda a, edge: qc[a]["scalp_vs_mri"]["offsets"]["cap"][edge]["median_mm"]
        kids = [cap(a, "otsu") for a in ("childA", "childB", "childC")]
        f = self.F["wr_qc_children_offset_inside_otsu_mm_range"]
        self.assertEqual(f["value"], f"{min(kids):.1f} to {max(kids):.1f}")
        self.assertTrue(all(v > 0 for v in kids))
        self.assertEqual(self.F["wr_qc_adult_offset_inside_otsu_mm"]["value"], f"{cap('adult', 'otsu'):.1f}")
        self.assertEqual(self.F["wr_qc_adult_offset_inside_halfmax_mm"]["value"], f"{cap('adult', 'half_max'):.1f}")
        t = cap("infant2yr", "otsu")
        self.assertLess(t, 0)
        self.assertEqual(self.F["wr_qc_infant2yr_offset_outside_otsu_mm"]["value"], f"{abs(t):.1f}")
        self.assertNotIn("−", self.F["wr_qc_infant2yr_offset_outside_otsu_mm"]["value"])

    def test_covariance_sub_band_range_and_count(self):
        sb = json.loads((ROOT / "results/g2_covariance_validation/covariance_validation.json").read_text())["sub_bands"]
        vals = [sb[b]["mag"]["empty_room_amp_ratio_model_over_measured"] for b in ("8-13Hz", "13-20Hz", "20-30Hz", "30-40Hz")]
        self.assertEqual(self.F["wr_cov_band_mag_emptyroom_from8hz_range"]["value"], f"{min(vals):.2f} to {max(vals):.2f}")
        n = len([k for k in sb if k.endswith("Hz") and k != "1-40Hz"])
        self.assertEqual(self.F["wr_n_subbands"]["raw"], n)
        self.assertEqual(self.F["wr_n_subbands_words"]["value"], W.words(n))

    def test_area_weighted_adult_ratio_from_the_stored_db(self):
        link = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())["link_to_g2"]
        r = 10 ** (link["adult_centred_area_weighted_cortical_dB"] / 20)
        f = self.F["wr_adult_dense_ib_ratio_area_weighted_cortical_3dp"]
        self.assertAlmostEqual(f["raw"], r, places=9)
        self.assertEqual(f["value"], f"{r:.3f}")
        # the pediatric convention gives a larger ratio than the unweighted median over all targets (1.144)
        self.assertGreater(r, link["adult_centred_unweighted_all_targets"])

    def test_adult_within_head_contrast_at_the_top_contact_gap(self):
        w = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())
        w = w["headline"]["by_anatomy"]["adult"]["within_head"]["gap_matched_top"]
        f = self.F["wr_cgap_adult_fixed_minus_fittedtop_combined_ib"]
        self.assertAlmostEqual(f["raw"], w["median"], places=9)
        self.assertEqual(f["value"], f"{w['median']:+.2f}".replace("-", "−"))
        ci = self.F["wr_cgap_adult_fixed_minus_fittedtop_combined_ib_ci"]
        self.assertEqual(ci["value"], f"[{w['ci95'][0]:+.2f}, {w['ci95'][1]:+.2f}]".replace("-", "−"))

    def test_seeds_and_plot_constants(self):
        db = json.loads((ROOT / "results/g2/g2_depth_bins.json").read_text())["method"]["seed"]
        self.assertEqual(self.F["wr_seed_g2_depth_bins"]["raw"], db)
        cfg = json.loads((ROOT / "results/g4/g4_motion_summary.json").read_text())["config"]
        self.assertEqual(self.F["wr_seed_motion_coupling"]["raw"], cfg["coupling"]["seed"])
        self.assertEqual(self.F["wr_seed_motion_timecourse"]["raw"], cfg["timecourse"]["seed"])
        desc = json.loads((ROOT / "results/report/figures_clean.json").read_text())["figures"]["Figure_R12_geometry"]["description"]
        self.assertIn(f"within {self.F['wr_fig_geometry_slab_mm']['value']} mm of the plane", desc)
        fig = json.loads((ROOT / "results/report/figures_qc.json").read_text())["figures"]["Figure_S_children_qc"]
        self.assertIn(f"each {self.F['wr_fig_qc_section_mm']['value']} mm wide", json.dumps(fig))
        self.assertEqual(self.F["wr_fig_qc_section_mm"]["value"], "50")

    def test_far_field_cardiac_depths(self):
        heart = json.loads((ROOT / "results/g2_noise_sensitivity/noise_sensitivity_summary.json").read_text())["far_field"]["heart_head_mm"]
        depths = sorted(abs(p[2]) for p in heart)
        self.assertEqual(self.F["wr_ns_ff_heart_depths_mm"]["raw"], depths)
        self.assertEqual(self.F["wr_ns_ff_heart_depths_mm"]["value"], "200, 250 and 300")
        self.assertEqual(self.F["wr_ns_ff_heart_mid_depth_mm"]["value"], "250")

    def test_oracle_above_practical_count_recomputed(self):
        labels = json.loads((ROOT / "results/g4/g4_pediatric_comparison.json").read_text())["labels"]
        n = 0
        for a in labels:
            p = json.loads((ROOT / f"results/g4/g4_{a}_summary.json").read_text())["paired"]
            o = p["opm_dense/opm_vs_squid/combined/oracle"]["depth0"]["s50_ratio_squid_over_opm"]["value"]
            q = p["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]["s50_ratio_squid_over_opm"]["value"]
            n += o > q
        self.assertEqual(self.F["wr_g4_oracle_above_practical_n"]["raw"], n)
        self.assertEqual(len(self.F["wr_g4_oracle_below_practical_heads"]["raw"]), len(labels) - n)
        self.assertEqual(self.F["wr_g4_oracle_below_practical_heads"]["value"], "18-month template and child C")

    def test_confirmatory_declarations_and_pilot_record(self):
        import tomllib
        text = (ROOT / "configs/g4_confirmatory.toml").read_text()
        cfg = tomllib.loads(text)
        v = cfg["variant"][0]
        self.assertEqual(self.F["wr_confirm_coreg_shift_mm"]["raw"], v["coreg_shift_mm"])
        self.assertEqual(self.F["wr_confirm_coreg_angle_deg"]["raw"], v["coreg_angle_deg"])
        self.assertEqual(self.F["wr_confirm_replicates_words"]["value"], W.words(cfg["design"]["noise_replicates"]))
        self.assertIn("two pilot runs", text)
        self.assertEqual(self.F["wr_confirm_pilot_runs_words"]["value"], "two")
        self.assertEqual(self.F["wr_confirm_pilot_locs"]["value"], "8/1 and 9/4")
        self.assertEqual(self.F["wr_confirm_pilot_p"]["value"], "0.074 and 0.21")
        self.assertEqual(self.F["wr_confirm_pilot_locs_per_band"]["raw"], 18)

    def test_confirmatory_equal_rate_counts_and_rate_cells(self):
        labels = json.loads((ROOT / "results/g4_confirm/g4_confirm_summary.json").read_text())["anatomies"]
        per = {a: json.loads((ROOT / f"results/g4_confirm/g4c_{a}_summary.json").read_text()) for a in labels}
        for thresholds in ("frozen", "matched"):
            key = f"opm_dense/opm_vs_squid/combined/primary/{thresholds}"
            low = [a for a in labels
                   if per[a]["false_event_rate_equality_on_evaluation_null"][key]["conditional_binomial_p"] < 0.05]
            n = self.F[f"wr_confirm_rate_equality_endpoint_{thresholds}_n_p05"]
            heads = self.F[f"wr_confirm_rate_equality_endpoint_{thresholds}_heads_p05"]
            self.assertEqual(n["raw"], len(low))
            self.assertEqual(n["value"], str(len(low)))
            self.assertEqual(heads["raw"], low)
            self.assertEqual(heads["value"], W._name_list([W.HEAD_LABELS[a] for a in low]) if low else "none")
        self.assertEqual(self.F["wr_confirm_rate_equality_endpoint_frozen_heads_p05"]["value"], "24-month template and child B")
        self.assertEqual(self.F["wr_confirm_rate_equality_endpoint_matched_heads_p05"]["value"], "school-age size")
        cells = sum(1 for a in labels for k in per[a]["false_events"] if k.endswith("|primary"))
        self.assertEqual(self.F["wr_confirm_n_rate_cells"]["raw"], cells)
        self.assertEqual(cells, 3 * len(labels))
        self.assertEqual(self.F["wr_confirm_n_rate_cells"]["value"], "27")

    def test_count_words_match_the_label_lists(self):
        labels = json.loads((ROOT / "results/g4/g4_pediatric_comparison.json").read_text())["labels"]
        self.assertEqual(self.F["wr_n_smaller_heads_words"]["raw"], len(labels) - 1)
        self.assertEqual(self.F["wr_n_smaller_heads_words"]["value"], "eight")
        self.assertEqual(self.F["wr_n_anatomies_words"]["value"], "nine")
        self.assertEqual(self.F["wr_n_children_words"]["value"], "three")
        self.assertEqual(self.F["wr_n_templates_words"]["value"], "three")
        self.assertEqual(self.F["wr_n_scaled_words"]["value"], "two")
        stretches = json.loads((ROOT / "results/g4/g4_adult_summary.json").read_text())["config"]["events"]["stretches"]
        self.assertEqual(self.F["wr_n_morphologies_words"]["value"], W.words(len(stretches)))
        with self.assertRaises(ValueError):
            W.words(13)

    def test_share_percentages_from_the_stored_shares(self):
        comp = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())["comparisons"]
        dense = comp["school/opm_dense/combined/intrinsic+brain/detect"]["d_adult"]["share_positive"] * 100
        self.assertEqual(self.F["wr_adult_dense_ib_share_pos_pct"]["value"], f"{dense:.1f}%")
        matched = comp["school/opm_matched/combined/intrinsic+brain/detect"]["d_adult"]["share_positive"] * 100
        self.assertEqual(self.F["wr_adult_matched_ib_share_pos_pct"]["value"], f"{matched:.0f}%")
        heads = [a for a in W.HEAD_LABELS if a != "adult"]
        vals = [comp[f"{a}/opm_matched/combined/intrinsic+brain/detect"]["d_child"]["share_positive"] * 100 for a in heads]
        self.assertEqual(self.F["wr_smaller_matched_ib_share_pos_pct_range"]["value"], f"{min(vals):.0f}% to {max(vals):.0f}%")

    def test_scenario_d_medians_recomputed_from_the_per_target_table(self):
        import csv
        import statistics
        cv = json.loads((ROOT / "results/g2_covariance_validation/covariance_validation.json").read_text())
        stored = cv["opm_implication"]["ratios"]["opm_dense/combined/S4_neuromag_measured_opm_as_modelled"]
        with (ROOT / "results/g2_covariance_validation/covariance_validation_targets.csv").open() as f:
            f.readline()  # the provenance comment
            rows = list(csv.DictReader(f))
        ratios = [float(r["detect_opm_dense_independent_intrinsic_brain_env"])
                  / float(r["detect_neuromag_combined_measured_corrected"]) for r in rows]
        self.assertEqual(len(rows), stored["n"])
        self.assertAlmostEqual(statistics.median(ratios), stored["ratio"], places=4)
        self.assertEqual(self.F["wr_cov_s4_dense_ratio_3dp"]["value"], f"{stored['ratio']:.3f}")
        for lobe in W.LOBES:
            vals = [x for r, x in zip(rows, ratios) if r["lobe"] == lobe]
            f = self.F[f"wr_cov_s4_dense_lobe_{lobe}_ratio"]
            self.assertAlmostEqual(f["raw"], statistics.median(vals), places=9)
            self.assertEqual(f["value"], f"{statistics.median(vals):.2f}")
        # the direction the text states: the lead shrinks most in the cingulate and grows in the temporal lobe
        self.assertLess(self.F["wr_cov_s4_dense_lobe_cingulate_ratio"]["raw"], 1.0)
        self.assertGreater(self.F["wr_cov_s4_dense_lobe_temporal_ratio"]["raw"], 1.0)
        first = None
        for lo in W.S4_DEPTH_EDGES:
            vals = [x for r, x in zip(rows, ratios) if lo <= float(r["depth_mm"]) < lo + 5]
            if len(vals) < W.S4_MIN_N:
                continue
            m = statistics.median(vals)
            self.assertAlmostEqual(self.F[f"wr_cov_s4_dense_depth_{lo:.0f}_{lo + 5:.0f}_ratio"]["raw"], m, places=9)
            if first is None and m < 1:
                first = lo
        self.assertEqual(self.F["wr_cov_s4_dense_depth_first_below1_label"]["raw"], [first, first + 5])
        self.assertEqual(self.F["wr_cov_s4_dense_depth_first_below1_label"]["value"], f"{first:.0f}\u2013{first + 5:.0f}")

    def test_bin_and_stratum_labels_match_the_stored_edges(self):
        cv = json.loads((ROOT / "results/g2_covariance_validation/covariance_validation.json").read_text())
        bins = {b.get("bin") for b in cv["structure"]["brain"]["mag"]["corr_vs_distance"]}
        for lo, hi in ((25, 50), (50, 75), (100, 125), (200, 225)):
            self.assertIn(f"{lo}-{hi}", bins)
            self.assertEqual(self.F[f"wr_cov_dist_d{lo}_{hi}_label"]["value"], f"{lo}\u2013{hi}")
        cg = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())
        strata = cg["delta_same_rule"]["childA/gap_matched/opm_dense/combined/intrinsic+brain"]["delta_by_depth"]
        edges = {(b["lo"], b["hi"]) for b in strata}
        for lo, hi in ((10, 15), (20, 25), (30, 40), (40, 50)):
            self.assertIn((float(lo), float(hi)), edges)
            self.assertEqual(self.F[f"wr_cgap_depth_{lo}_{hi}_label"]["value"], f"{lo}\u2013{hi}")

    def test_fitted_helmet_interval_bounds_and_children_below_the_adult(self):
        d = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())["delta_same_rule"]
        cis = [d[f"{h}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["ci95"]
               for h in ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")]
        lo, hi = min(c[0] for c in cis), max(c[1] for c in cis)
        f = self.F["wr_cgap_scaled_templates_fitted_delta_ci_bounds_range"]
        self.assertAlmostEqual(f["raw"][0], lo, places=9)
        self.assertAlmostEqual(f["raw"][1], hi, places=9)
        self.assertEqual(f["value"], f"{lo:+.2f} to {hi:+.2f}".replace("-", "\u2212"))
        self.assertEqual(f["value"], "\u22120.09 to +0.33")
        n = sum(d[f"{c}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["median"] < 0
                for c in ("childA", "childB", "childC"))
        self.assertEqual(self.F["wr_n_children_below_adult_fitted_words"]["raw"], n)
        self.assertEqual(self.F["wr_n_children_below_adult_fitted_words"]["value"], "two")

    def test_thresholds_counts_and_depth_bin_width(self):
        loc = json.loads((ROOT / FIGS_SUPP).read_text())["figures"]["Figure_S_localization_effects"]["values"]
        self.assertAlmostEqual(self.F["wr_p_exploratory_threshold"]["raw"],
                               loc["bonferroni_threshold"] * loc["family_per_anatomy_and_array"], places=12)
        self.assertEqual(self.F["wr_p_exploratory_threshold"]["value"], "0.05")
        qc = json.loads((ROOT / "results/g3b_children_qc/children_qc.json").read_text())
        self.assertEqual([self.F["wr_qc_near_lo_mm"]["raw"], self.F["wr_qc_near_hi_mm"]["raw"]], qc["parameters"]["near_mm"])
        self.assertEqual(self.F["wr_qc_near_lo_mm"]["value"], "8")
        self.assertEqual(self.F["wr_qc_near_hi_mm"]["value"], "10")
        g3b = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())["comparisons"]
        full = g3b["childB/opm_dense/combined/intrinsic+brain/detect"]["delta_by_depth"][0]["n_child"]
        # the check now counts from the full-precision depths (results/g3b/g3b_targets_childB_depth.csv)
        self.assertEqual(qc["anatomies"]["childB"]["g3b_targets"]["below_10mm"]["n"], full)
        self.assertNotIn("wr_qc_child_b_targets_lt10mm_shortfall", self.F)
        ns = json.loads((ROOT / "results/g2_noise_sensitivity/noise_sensitivity_summary.json").read_text())
        bins = ns["sweep"]["entries"]["15"]["depth"]["opm_dense/combined/intrinsic+brain"]
        self.assertEqual({b["hi"] - b["lo"] for b in bins}, {self.F["wr_ns_depth_bin_mm"]["raw"]})
        self.assertEqual(self.F["wr_ns_depth_bin_mm"]["value"], "5")

    def test_sweep_list_and_depth_band_words(self):
        levels = self.g2["config"]["sensors"]["opm_asd_fT_per_rtHz"]
        self.assertEqual(self.F["wr_opm_asd_sweep_list"]["raw"], [float(x) for x in levels])
        self.assertEqual(self.F["wr_opm_asd_sweep_list"]["value"], "7, 10, 15, 20 and 30")
        locs = json.loads((ROOT / "results/g4/g4_adult_summary.json").read_text())["locations"]
        n = len({l["stratum"][0] for l in locs})
        self.assertEqual(self.F["wr_n_depth_bands_words"]["raw"], n)
        self.assertEqual(self.F["wr_n_depth_bands_words"]["value"], "four")

    def test_background_scaling_against_the_unscaled_adult(self):
        d = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())
        adult = d["D_median_dB"]["adult/opm_dense/combined/intrinsic+brain/detect"]
        heads = ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
        both = []
        for factor, tok in (("0.5", "x0p5"), ("2", "x2")):
            vals = [d["sensitivity_median_D_dB"][f"{h}/background_x{factor}/opm_dense/combined/intrinsic+brain"] - adult
                    for h in heads]
            f = self.F[f"wr_g3b_smaller_bgonly_{tok}_d_minus_adult_range"]
            self.assertAlmostEqual(f["raw"][0], min(vals), places=9)
            self.assertAlmostEqual(f["raw"][1], max(vals), places=9)
            self.assertEqual(f["value"], f"{min(vals):+.2f} to {max(vals):+.2f}".replace("-", "−"))
            # the stored 'delta/...' sensitivity scales the adult's background too, so it differs from this difference
            stored = [d["sensitivity_median_D_dB"][f"delta/{h}/background_x{factor}/opm_dense/combined/intrinsic+brain"]
                      for h in heads]
            self.assertNotAlmostEqual(min(stored), min(vals), places=3)
            self.assertTrue(all(v > 0 for v in vals))
            both += vals
        f = self.F["wr_g3b_smaller_bgonly_d_minus_adult_range"]
        self.assertAlmostEqual(f["raw"][0], min(both), places=9)
        self.assertAlmostEqual(f["raw"][1], max(both), places=9)
        self.assertEqual(f["value"], "+0.12 to +1.16")

    def test_maps_heads_colour_limit_and_sphere_depth(self):
        v = json.loads((ROOT / "results/report/figures_clean.json").read_text())["figures"]["Figure_R13_maps_heads"]["values"]
        self.assertEqual(self.F["wr_fig_maps_heads_colour_limit_db"]["raw"], v["colour_limit_dB"])
        self.assertEqual(self.F["wr_fig_maps_heads_colour_limit_db"]["value"], "6")
        for head in ("infant2yr", "infant12mo"):
            s = v[head]["share_above_limit"] * 100
            self.assertEqual(self.F[f"wr_fig_maps_heads_above_limit_pct_{head}"]["value"], f"{s:.1f}%")
        self.assertEqual(self.F["wr_fig_maps_heads_above_limit_pct_infant2yr"]["value"], "3.1%")
        self.assertEqual(self.F["wr_fig_maps_heads_above_limit_pct_infant12mo"]["value"], "5.7%")
        p = json.loads((ROOT / "results/g1a/g1a_benchmark.json").read_text())["parameters"]
        self.assertEqual(self.F["wr_sphere_brain_depth_mm"]["raw"], p["h_mm"] - p["b_mm"])
        self.assertEqual(self.F["wr_sphere_brain_depth_mm"]["value"], "15")

    def test_maps_heads_cortical_shares_recomputed_from_the_target_tables(self):
        import csv
        v = json.loads((ROOT / "results/report/figures_clean.json").read_text())["figures"]["Figure_R13_maps_heads"]["values"]
        lim = v["colour_limit_dB"]
        for head in ("infant2yr", "infant12mo"):
            with (ROOT / f"results/g3b/g3b_targets_{head}.csv").open() as f:
                f.readline()  # the provenance comment
                rows = list(csv.DictReader(f))
            d = [20 * math.log10(float(r["detect_opm_dense_opm_intrinsic+brain"])
                                 / float(r["detect_squid_top_combined_intrinsic+brain"])) for r in rows]
            self.assertEqual(len(rows), v[head]["n_targets"])
            self.assertAlmostEqual(sum(x > lim for x in d) / len(d), v[head]["share_above_limit"], places=12)
            cort = [x for r, x in zip(rows, d) if not r["region"].endswith(".unknown")]
            self.assertEqual(len(cort), v[head]["n_targets"] - v[head]["medial_wall"])
            share = sum(x > lim for x in cort) / len(cort) * 100
            f = self.F[f"wr_fig_maps_heads_above_limit_cortical_pct_{head}"]
            self.assertAlmostEqual(f["raw"], share, places=9)
            self.assertEqual(f["value"], f"{share:.1f}%")
            # the cortical share exceeds the all-target share, since no medial-wall target lies above the limit
            self.assertGreater(f["raw"], self.F[f"wr_fig_maps_heads_above_limit_pct_{head}"]["raw"])
        self.assertEqual(self.F["wr_fig_maps_heads_above_limit_cortical_pct_infant2yr"]["value"], "3.3%")
        self.assertEqual(self.F["wr_fig_maps_heads_above_limit_cortical_pct_infant12mo"]["value"], "6.1%")

    def test_projected_fitted_helmet_intervals_spanning_zero(self):
        d = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())["delta_same_rule"]
        n = 0
        for h in ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC"):
            lo, hi = d[f"{h}/gap_matched/opm_dense/combined/projected"]["delta"]["ci95"]
            n += lo <= 0 <= hi
        f = self.F["wr_cgap_fitted_proj_delta_n_ci_includes_zero_words"]
        self.assertEqual(f["raw"], n)
        self.assertEqual(f["value"], W.words(n))
        self.assertEqual(f["value"], "three")

    def test_anatomies_below_the_adult_depth_ratio_in_the_confirmatory_run(self):
        bins = json.loads((ROOT / "results/g2/g2_depth_bins.json").read_text())["comparisons"]
        ref = next(b["ratio"] for b in bins["opm_dense/combined/intrinsic+brain"]["bins"] if b["lo"] == 15 and b["hi"] == 20)
        labels = json.loads((ROOT / "results/g4_confirm/g4_confirm_summary.json").read_text())["anatomies"]
        below = []
        for a in labels:
            c = json.loads((ROOT / f"results/g4_confirm/g4c_{a}_summary.json").read_text())["comparisons"]
            if c["opm_dense/opm_vs_squid/combined/practical@1/replicate0"]["s50_ratio_squid_over_opm"]["value"] < ref:
                below.append(a)
        f = self.F["wr_cf_dense_practical_ratio_below_adult_depth_heads"]
        self.assertEqual(f["raw"], below)
        self.assertEqual(below, ["infant2yr", "childA", "childB", "childC"])
        self.assertEqual(f["value"], "the 24-month template, child A, child B and child C")
        # the adult's own ratio lies above its 15-20 mm detectability ratio
        self.assertNotIn("adult", below)

    def test_exact_sign_flip_enumeration_and_the_endpoint_bound(self):
        # the enumeration reproduces the stored exact values (20 or fewer non-zero differences) and simple cases
        self.assertEqual(W.exact_sign_flip_p([]), 1.0)
        self.assertEqual(W.exact_sign_flip_p([1, 1, 1]), 0.25)  # all three signs alike: 2 of 8 patterns
        self.assertEqual(W.exact_sign_flip_p([2, -1]), 1.0)  # |sum| = 1 is reached by every pattern
        labels = json.loads((ROOT / "results/g4_confirm/g4_confirm_summary.json").read_text())["anatomies"]
        fam = "opm_dense/opm_vs_squid/combined/practical@1/replicate0"
        exact, stored = {}, {}
        for a in labels:
            c = json.loads((ROOT / f"results/g4_confirm/g4c_{a}_summary.json").read_text())["comparisons"][fam]
            exact[a] = W.exact_sign_flip_p(c["location_differences"])
            stored[a] = c["location_sign_flip_p"]
            if sum(1 for v in c["location_differences"] if v != 0) <= 20:
                self.assertAlmostEqual(exact[a], stored[a], places=12)
        G4 = _load("report_facts_g4")
        adj = G4._holm(exact)
        zeros = [a for a in labels if stored[a] == 0]
        self.assertTrue(zeros)
        m = max(adj[a] for a in zeros)
        f = self.F["wr_cf_endpoint_p_holm_exact_max_mc_zero"]
        self.assertAlmostEqual(f["raw"], m, places=15)
        self.assertLess(m, 1e-4)  # every endpoint value printed '<0.0001' is below 0.0001 by exact enumeration
        self.assertEqual(f["value"], "2.2 × 10⁻⁵")
        self.assertEqual(W._sci(5e-5), "5.0 × 10⁻⁵")
        self.assertEqual(W._sci(0.0099), "9.9 × 10⁻³")

    def test_round3_facts_against_their_files(self):
        g2s = json.loads((ROOT / "results/g2/g2_summary.json").read_text())
        self.assertAlmostEqual(self.F["wr_g2_dense_ib_ratio_no_medial_wall"]["raw"],
                               g2s["medial_wall"]["ratios_without"]["opm_dense/combined/intrinsic+brain"], places=12)
        qc = json.loads((ROOT / "results/g3b_children_qc/children_qc.json").read_text())["anatomies"]["childB"]["g3b_targets"]
        self.assertEqual(self.F["wr_qc_child_b_targets_lt10mm_mri_n"]["raw"], qc["below_10mm_to_mri_boundary"])
        ex = json.loads((ROOT / "results/g3b/g3b_cortex_exclusion.json").read_text())["anatomies"]
        for key, tok in (("adult", "adult"), ("childA", "child_a"), ("childB", "child_b"), ("childC", "child_c")):
            w = ex[key]["whole_surface"]
            self.assertAlmostEqual(self.F[f"wr_excl_{tok}_pct"]["raw"],
                                   100 * (w["within_4mm_of_inner_skull_share"] + w["outside_inner_skull_share"]), places=9)
        self.assertEqual(self.F["wr_excl_child_b_lt8mm_usable_n"]["raw"], ex["childB"]["within_8mm_of_scalp"]["n_usable_vertices"])
        self.assertEqual(self.F["wr_confirm_calib_expected_events"]["value"], "20")
        self.assertEqual(self.F["wr_confirm_calib_poisson_pct"]["value"], "22%")

    def test_round3b_facts_against_their_files(self):
        ex = json.loads((ROOT / "results/g3b/g3b_cortex_exclusion.json").read_text())["anatomies"]
        self.assertEqual(set(ex), {"adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC"})
        left = {k: 100 * (v["whole_surface"]["within_4mm_of_inner_skull_share"] + v["whole_surface"]["outside_inner_skull_share"])
                for k, v in ex.items()}
        self.assertEqual(self.F["wr_excl_scaled_pct_range"]["raw"], [min(left["school"], left["size2yr"]),
                                                                     max(left["school"], left["size2yr"])])
        near = [100 * (1 - ex[k]["within_8mm_of_scalp"]["usable_share"]) for k in ("childA", "childB", "childC")]
        self.assertEqual(self.F["wr_excl_children_lt8mm_pct_range"]["raw"], [min(near), max(near)])
        cg = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())
        heads = ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")
        it = [cg["interaction"][f"{h}/gap_matched/combined/intrinsic+brain"]["interaction"] for h in heads]
        self.assertTrue(all(x["ci95"][0] > 0 for x in it))
        self.assertEqual(self.F["wr_cgap_scaled_templates_interaction_range"]["raw"],
                         [min(x["median"] for x in it), max(x["median"] for x in it)])
        top = [cg["delta_vs_adult_top"][f"{h}/gap_matched_top/opm_dense/combined/intrinsic+brain"]["delta"] for h in heads]
        self.assertEqual(self.F["wr_cgap_scaled_templates_fittedtop_vs_top_n_ci_above0"]["raw"],
                         sum(x["ci95"][0] > 0 for x in top))
        self.assertRegex(self.F["wr_build_commit"]["value"], r"^([0-9a-f]{7,}(\+dirty)?|\(not a git checkout\))$")

    def test_round3c_facts_against_their_files(self):
        qt = json.loads((ROOT / "results/g3b_templates_qc/children_qc.json").read_text())
        q24 = json.loads((ROOT / "results/g3b_children_qc/children_qc.json").read_text())
        self.assertEqual(set(qt["anatomies"]), {"adult", "infant18mo", "infant12mo"})
        self.assertEqual(qt["anatomies"]["adult"]["scalp_vs_mri"], q24["anatomies"]["adult"]["scalp_vs_mri"])
        out = {k: -q["anatomies"][k]["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"]
               for k, q in (("infant2yr", q24), ("infant18mo", qt), ("infant12mo", qt))}
        self.assertTrue(all(v > 0 for v in out.values()))
        self.assertEqual(self.F["wr_qc_templates_offset_outside_otsu_mm_range"]["raw"], [min(out.values()), max(out.values())])
        self.assertEqual(self.F["wr_qc_infant18mo_offset_outside_otsu_mm"]["raw"], out["infant18mo"])
        adult = q24["anatomies"]["adult"]["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"]
        self.assertEqual(self.F["wr_qc_template_adult_convention_mm_range"]["raw"],
                         [adult + min(out.values()), adult + max(out.values())])
        self.assertEqual(self.F["wr_qc_infant18mo_verdict"]["value"], qt["verdicts"]["infant18mo"]["class"])
        self.assertEqual(self.F["wr_qc_infant12mo_verdict"]["value"], qt["verdicts"]["infant12mo"]["class"])
        self.assertEqual(self.F["wr_qc_infant18mo_white_best_shift_mm"]["raw"],
                         qt["anatomies"]["infant18mo"]["white_vs_t1"]["best_shift_norm_mm"])
        with open(ROOT / "results/g2/g2_targets_metrics.csv", newline="") as fh:  # per-target flags, independent of the summary
            rows = [r for r in csv.DictReader(line for line in fh if not line.startswith("#"))]
        self.assertEqual(self.F["wr_g2_ahead_n"]["raw"], sum(int(r["ahead_dense_intrinsic+brain"]) for r in rows))
        self.assertEqual(self.F["wr_g2_ahead_n"]["value"], "7,657")
        self.assertEqual(self.F["wr_cover_brow_offset_mm"]["raw"], 30.0)
        ci = json.loads((ROOT / "results/g4_confirm/g4c_childC_summary.json").read_text())["comparisons"][
            "opm_matched/opm_vs_squid/combined/practical@1/replicate0"]["s50_ratio_squid_over_opm"]["ci95"]
        self.assertLess(self.F["wr_cf_childc_matched_practical_ci_lo_4dp"]["raw"], 1)
        self.assertEqual(self.F["wr_cf_childc_matched_practical_ci_lo_4dp"]["raw"], ci[0])
        cg = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())["interaction"]
        it = {h: cg[f"{h}/gap_matched/combined/intrinsic+brain"]["interaction"]["median"]
              for h in ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")}
        tm = [it[h] for h in ("infant2yr", "infant18mo", "infant12mo")]
        self.assertEqual(self.F["wr_cgap_templates_interaction_range"]["raw"], [min(tm), max(tm)])
        lo, hi = self.F["wr_cgap_templates_interaction_range"]["raw"]
        self.assertTrue(min(it["school"], it["size2yr"]) <= lo <= hi <= max(it["school"], it["size2yr"]))
        self.assertEqual(self.F["wr_cgap_templates_interaction_range"]["value"], "+0.54 to +0.70")
        # the remaining template-check facts, and printed values
        for key, q in (("infant18mo", qt), ("infant12mo", qt), ("infant2yr", q24)):
            cap = q["anatomies"][key]["scalp_vs_mri"]["offsets"]["cap"]
            self.assertEqual(self.F[f"wr_qc_{key}_offset_outside_halfmax_mm"]["raw"], -cap["half_max"]["median_mm"])
        self.assertEqual(self.F["wr_qc_infant12mo_offset_outside_otsu_mm"]["raw"],
                         -qt["anatomies"]["infant12mo"]["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"])
        w18 = qt["anatomies"]["infant18mo"]["white_vs_t1"]
        self.assertEqual(self.F["wr_qc_infant18mo_white_contrast_gain_pct"]["raw"], 100 * w18["contrast_gain_share"])
        self.assertEqual(self.F["wr_qc_infant18mo_white_contrast_gain_pct"]["value"], "41%")
        self.assertEqual(self.F["wr_qc_infant18mo_white_edge_translation_mm"]["raw"],
                         qt["anatomies"]["infant18mo"]["verdict"]["checks"]["white_edge_translation_mm"])
        self.assertGreater(w18["best_shift_norm_mm"], qt["parameters"]["tolerance_mm"])
        self.assertEqual(self.F["wr_qc_templates_other_white_best_shift_mm"]["raw"],
                         [q24["anatomies"]["infant2yr"]["white_vs_t1"]["best_shift_norm_mm"],
                          qt["anatomies"]["infant12mo"]["white_vs_t1"]["best_shift_norm_mm"]])
        dch = [q["anatomies"][k]["g3b_targets"]["depth_change_to_mri_boundary_mm"]["p50"]
               for k, q in (("infant2yr", q24), ("infant18mo", qt), ("infant12mo", qt))]
        self.assertEqual(self.F["wr_qc_templates_targets_depth_change_p50_mm_range"]["raw"], [min(dch), max(dch)])
        self.assertEqual(self.F["wr_cf_childc_matched_practical_ci_lo_4dp"]["value"], "0.9997")
        self.assertEqual(self.F["wr_qc_template_adult_convention_mm_range"]["value"], "2.6 to 3.0")

    def test_round4_facts_against_their_files(self):
        qt = json.loads((ROOT / "results/g3b_templates_qc/children_qc.json").read_text())["anatomies"]["infant18mo"]
        self.assertEqual(self.F["wr_qc_infant18mo_white_scalp_min_mm"]["raw"], qt["white_to_scalp_used"]["min_mm"])
        self.assertAlmostEqual(qt["white_to_scalp_used"]["min_mm"], qt["closest_white_vertex"]["distance_mm"], places=9)
        self.assertEqual(self.F["wr_qc_infant18mo_white_scalp_min_mm"]["value"], "7.3")

    def test_registered_last_in_the_merged_facts(self):
        R = _load("report_facts")
        self.assertEqual(R.MODULES[-1], "report_facts_writer")
        merged = R.build_facts(ROOT)
        for name, f in self.F.items():
            self.assertEqual(merged[name]["value"], f["value"])
            self.assertEqual(merged[name]["module"], "report_facts_writer")


if __name__ == "__main__":
    unittest.main()
