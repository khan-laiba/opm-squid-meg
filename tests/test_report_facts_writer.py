"""Facts added by the writers of the revised report (scripts/report_facts_writer.py): the prefix, the source form, the
existence of every source file, and every value recomputed from the stored result files and the register."""
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

    def test_registered_last_in_the_merged_facts(self):
        R = _load("report_facts")
        self.assertEqual(R.MODULES[-1], "report_facts_writer")
        merged = R.build_facts(ROOT)
        for name, f in self.F.items():
            self.assertEqual(merged[name]["value"], f["value"])
            self.assertEqual(merged[name]["module"], "report_facts_writer")


if __name__ == "__main__":
    unittest.main()
