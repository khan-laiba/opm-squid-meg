"""Report facts of the referee round-1 revision analyses (scripts/report_facts_rev.py) against the stored results:
formats, sources, registration, literal values read from the result files, the study's own printed summary, and the
guard on the noise-sensitivity configuration."""
import importlib.util
import json
import math
import re
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NS = "results/g2_noise_sensitivity/noise_sensitivity_summary.json"
COV = "results/g2_covariance_validation/covariance_validation.json"
QC = "results/g3b_children_qc/children_qc.json"
CFG_NS = "configs/g2_noise_sensitivity.toml"
CGAP = "results/g3b_constant_gap/g3b_constant_gap_summary.json"
DB = "results/g2/g2_depth_bins.json"
G2T = "results/g2/g2_targets.csv"
NEEDED = [NS, COV, QC, CFG_NS, "configs/school_subjects_qc_manifest.json", CGAP, DB]
HAVE = all((ROOT / p).exists() for p in NEEDED)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _load("report_facts_rev")


class TestFormats(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(R.spct(0.04), "0.0%")
        self.assertEqual(R.spct(-1.26), "−1.3%")
        self.assertEqual(R.spct(0.66), "+0.7%")
        self.assertEqual(R.pctk(1499.08), "1,499%")
        self.assertEqual(R.pctk(38.0), "38%")
        self.assertEqual(R.listing([]), "none")
        self.assertEqual(R.listing([7.0, 10.0, 15.0]), "7, 10 and 15")
        self.assertEqual(R.small(0), "0")
        self.assertEqual(R.small(2.74e-05), "2.7 × 10⁻⁵")
        self.assertEqual(R.diag(5092.7), "5,093")
        self.assertEqual(R.diag(19.88), "20")
        self.assertEqual(R.diag(0.9929), "0.99")
        self.assertEqual(R.g(-30.0), "−30")
        self.assertEqual(R.ctok("opm_matched/grad/projected"), "matched_vs_grad_proj")
        self.assertEqual(R.band_tok(0.5, 1.0), "0p5_1")
        self.assertEqual(R.sb_tok(48.0, 300.3), "sb48_up")

    def test_facts_container_enforces_the_conventions(self):
        F = R.Facts()
        with self.assertRaises(ValueError):
            F.add("meth_x", "1", 1, f"{NS} :: x")  # another module's prefix
        with self.assertRaises(ValueError):
            F.add("rev_x", " ", 1, f"{NS} :: x")  # empty value
        with self.assertRaises(ValueError):
            F.add("rev_x", "1", 1, NS)  # no key path
        F.add("rev_x", "1", 1, f"{NS} :: x")
        with self.assertRaises(ValueError):
            F.add("rev_x", "1", 1, f"{NS} :: x")  # defined twice


@unittest.skipUnless(HAVE, "revision result files not present")
class TestRevisionFacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.F = R.facts(ROOT)
        cls.v = {k: f["value"] for k, f in cls.F.items()}
        cls.ns, cls.cov, cls.qc = (json.loads((ROOT / p).read_text()) for p in (NS, COV, QC))

    def test_every_fact_has_a_value_and_an_existing_source(self):
        self.assertGreater(len(self.F), 5000)
        for name, f in self.F.items():
            self.assertRegex(name, r"^rev_(ns|cov|qc|cgap|db|g2_notahead|seed|boot)_[a-z0-9_]+$")
            self.assertEqual(set(f), {"value", "raw", "source"}, name)
            self.assertTrue(str(f["value"]).strip(), name)
            files, sep, how = f["source"].partition(" :: ")
            self.assertTrue(sep and how.strip(), f"{name}: source lacks ' :: <key path>'")
            paths = [p for p in files.replace(",", " ").split() if p.startswith(("results/", "configs/", "docs/"))]
            self.assertTrue(paths, f"{name}: no source file")
            for p in paths:
                self.assertTrue((ROOT / p).exists(), f"{name}: {p} not found")
            json.dumps(f["raw"])  # plain JSON types only
        self.assertEqual(_load("report_facts").check(self.F, ROOT), [])

    def test_number_formats(self):
        """U+2212 for minus, no signed zero, intervals '[lo, hi]', ranges 'lo to hi'."""
        for name, f in self.F.items():
            v = str(f["value"])
            self.assertNotRegex(v, r"(^|[\s\[(])-\d", f"{name}: ASCII minus")
            if re.search(r"\d", v) and not re.search(r"[A-Za-z]", v):
                self.assertNotIn("-", v, name)
            self.assertNotRegex(v, "[+−]0\\.0+(?![0-9])", f"{name}: signed zero")
            if name.endswith(("_ci", "_ci_3dp", "_ci_2dp")):
                self.assertRegex(v, r"^\[\S+, \S+\]$", name)
            if name.endswith("_range"):
                self.assertRegex(v, r"^\S+ to \S+$", name)

    def test_registered_after_the_methods_module_without_clashes(self):
        rf = _load("report_facts")
        self.assertEqual(rf.MODULES.index("report_facts_rev"), rf.MODULES.index("report_facts_methods") + 1)
        combined = rf.build_facts(ROOT)
        self.assertTrue(all(combined[k]["module"] == "report_facts_rev" for k in self.F))

    def test_reserved_prefixes_cite_their_own_results(self):
        """The constant-gap (rev_cgap_) and confirmatory (rev_conf_) facts, once added, must cite their result files."""
        self.assertEqual(set(R.RESERVED), {"rev_cgap_", "rev_conf_"})
        for prefix, rel in R.RESERVED.items():
            self.assertTrue(rel.startswith("results/"), rel)
            for name, f in self.F.items():
                if name.startswith(prefix):
                    self.assertIn(rel, f["source"].split(" :: ")[0], name)

    def test_key_values(self):
        """Literals read from the three result files (comments: the stored value)."""
        expect = {
            # noise sensitivity
            "rev_ns_n_boot": "1,000",
            "rev_ns_sweep_asd15ft_dense_vs_combined_ib_ratio_3dp": "1.144",      # sweep.entries['15'] ratio 1.14391
            "rev_ns_sweep_asd15ft_dense_vs_combined_ib_ci_3dp": "[1.124, 1.167]",  # ci95_ratio [1.12356, 1.16730]
            "rev_ns_sweep_asd15ft_dense_vs_combined_ib_n_neuromag_ahead": "4",     # (1 - 0.99948) x 7,661
            "rev_ns_sweep_asd30ft_matched_vs_combined_ib_ratio": "0.93",           # 0.93309
            "rev_ns_sweep_matched_vs_combined_proj_levels_ci_below1": "15, 20 and 30",
            "rev_ns_sweep_dense_vs_combined_ib_range_3dp": "1.008 to 1.298",
            "rev_ns_breakeven_dense_vs_combined_ib_asd": "31.5",                   # 31.5157
            "rev_ns_breakeven_dense_vs_combined_ib_asd_ci": "[28.5, 34.9]",
            "rev_ns_breakeven_matched_vs_combined_proj_asd": "8.9",                # 8.8809
            "rev_ns_near_n_vertices": "27,059",
            "rev_ns_near_share_valid_pct": "8.7%",                                 # 0.08673
            "rev_ns_near_f2mm_n_targets": "463",
            "rev_ns_near_primary_f2mm_bg_dense_vs_combined_ib_ratio_3dp": "1.137",  # 1.13664
            "rev_ns_near_primary_f2mm_targets_near_dense_vs_combined_ib_ratio": "1.40",  # 1.40082
            "rev_ns_near_effect_primary_f2mm_bg_dense_vs_combined_ib_pct": "−0.6%",  # 2**-0.009191 - 1
            "rev_ns_near_gainchange_dense_band_0_0p5_median_pct": "38%",           # 0.37523
            "rev_ns_col_fr_asd15ft_c10hz_dense_vs_combined_ib_ratio_3dp": "1.083",  # 1.08313
            "rev_ns_col_fr_asd15ft_c10hz_dense_vs_combined_ib_vs_white_pct": "−3.9%",  # 1.08313 / 1.12712 - 1
            "rev_ns_col_bv_factor_c10hz_ratio": "1.99",                            # 1.98622
            "rev_ns_ff_left_heart_cd250_grad_pct": "52%",                          # 0.51807
            "rev_ns_ff_left_eye_bem_dense_pct": "22%",                             # 0.22435
            "rev_ns_ff_level_short_ft": "178",                                     # 178.34
            "rev_ns_ff_n_infeasible": "1",
            "rev_ns_joint_all_asd15ft_matched_vs_combined_proj_ratio_3dp": "0.890",  # 0.89007
            "rev_ns_alt_asd15ft_n_rows": "20",
            # covariance validation
            "rev_cov_brain_mag_amp_ratio": "0.73",                                 # 0.73281
            "rev_cov_brain_mag_var_ratio_p5_p95_range": "0.13 to 0.78",
            "rev_cov_heldout_inf_to_sup_grad_amp_ratio": "1.34",                   # 1.33705
            "rev_cov_heart_mag_share_shortfall_pct": "51%",                        # 0.50682
            "rev_cov_n_fitted": "2", "rev_cov_n_predicted": "7",
            "rev_cov_det_combined_corrected_ratio_3dp": "1.142",                   # 1.14173
            "rev_cov_surr_null_corrected_combined_range_3dp": "0.997 to 1.006",
            "rev_cov_opm_dense_vs_combined_s4_ratio_3dp": "1.034",                 # 1.03385
            "rev_cov_opm_matched_vs_combined_s3_ci_3dp": "[0.647, 0.703]",
            "rev_cov_band_b30_40_mag_amp_ratio": "0.33",                           # 0.32833
            "rev_cov_det_lobe_temporal_corrected_ratio": "0.89",                   # 0.89277
            # children's MRI quality check
            "rev_qc_child_a_white_scalp_min_mm": "3.7",                            # 3.6962
            "rev_qc_child_c_white_scalp_min_mm": "5.7",                            # 5.7360
            "rev_qc_child_a_white_mri_min_mm": "5.2",                              # 5.2003
            "rev_qc_child_a_closest_region": "rh.lateraloccipital",
            "rev_qc_child_b_offset_cap_otsu_median_mm": "+1.5",                    # 1.47782
            "rev_qc_infant2yr_outward_mm": "1.6",                                  # -(-1.63892)
            "rev_qc_infant2yr_outward_excess_mm": "3.0",                           # -(-2.98547)
            "rev_qc_child_a_ofc_mri_cm": "52.9",                                   # 528.73 mm
            "rev_qc_child_b_targets_lt10mm_n": "247",  # full-precision depths since the re-run at bd70b45
            "rev_qc_child_b_targets_lt10mm_mri_n": "47",
            "rev_qc_child_a_verdict": "usable", "rev_qc_infant2yr_verdict": "outside",
            "rev_qc_n_children_usable": "3",
            "rev_qc_test_inward_recovered_mm": "4.0",                              # 3.99999
        }
        for name, value in expect.items():
            self.assertEqual(self.v[name], value, name)

    def test_raw_values_are_the_stored_ones(self):
        e = self.ns["sweep"]["entries"]["15"]["comparisons"]["opm_dense/combined/intrinsic+brain"]
        self.assertAlmostEqual(self.F["rev_ns_sweep_asd15ft_dense_vs_combined_ib_ratio"]["raw"], 2 ** e["median_log2"], places=12)
        be = self.ns["sweep"]["break_even"]["entries"]["opm_matched/combined/projected"]
        self.assertEqual(self.F["rev_ns_breakeven_matched_vs_combined_proj_asd"]["raw"], be["break_even_fT"])
        self.assertEqual(self.F["rev_ns_breakeven_matched_vs_combined_proj_asd_ci"]["raw"], be["ci95_fT"])
        x = self.ns["near_skull"]["inclusion_effect_log2"]["opm_matched/combined/projected"]["one_layer_refined/floor_0mm/both"]
        self.assertAlmostEqual(self.F["rev_ns_near_effect_refined_f0mm_both_matched_vs_combined_proj_pct"]["raw"], 100 * (2 ** x - 1))
        self.assertEqual(self.F["rev_cov_k_ratio"]["raw"], self.cov["model"]["magnetometer_calibration_scale_ratio"])
        off = self.qc["anatomies"]["infant2yr"]["scalp_vs_mri"]["offsets"]["cap"]["otsu"]["median_mm"]
        self.assertEqual(self.F["rev_qc_infant2yr_outward_mm"]["raw"], -off)
        hc = self.qc["anatomies"]["childA"]["head_circumference"]
        self.assertAlmostEqual(self.F["rev_qc_child_a_ofc_mri_minus_surface_mm"]["raw"],
                               hc["mri_boundary_points"]["ofc_mm"] - hc["surface_used"]["ofc_mm"])
        # the worst joint case at 15 fT/sqrt(Hz) is the lowest of the joint runs there
        key = "opm_matched/combined/projected"
        runs = [r["comparisons"][key]["ratio"] for k, r in self.ns["joint"]["results"].items()
                if k.endswith("/15fT") and not k.startswith("baseline_white")]
        self.assertEqual(self.F["rev_ns_joint_worst_asd15ft_matched_vs_combined_proj_ratio"]["raw"], min(runs))

    def test_the_study_summary_is_reproduced(self):
        """The covariance validation prints its own plain summary: the facts must read the same."""
        text = " ".join(self.cov["plain_summary"])
        v = self.v
        quoted = [
            f"model/measured {v['rev_cov_brain_mag_amp_ratio']}",
            f"would give {v['rev_cov_brain_mag_amp_null_range'].replace(' to ', '-')}",
            f"{v['rev_cov_heldout_random_grad_amp_ratio']} {v['rev_cov_heldout_random_grad_amp_ci']}",
            f"but {v['rev_cov_heldout_inf_to_sup_grad_amp_ratio']} (fitted below",
            f"and {v['rev_cov_heldout_sup_to_inf_grad_amp_ratio']} (fitted above",
            f"magnetometers {v['rev_cov_brain_mag_log_var_r']}, gradiometers {v['rev_cov_brain_grad_log_var_r']}",
            f"{v['rev_cov_brain_mag_var_ratio_p5_p95_range'].replace(' to ', '-')} (magnetometers)",
            f"{v['rev_cov_brain_grad_var_ratio_p5_p95_range'].replace(' to ', '-')} (gradiometers)",
            f"leading 3 components hold {v['rev_cov_excess_brain_mag_top3_pct']}",
            f"{v['rev_cov_heart_n_beats']} beats) adds {v['rev_cov_heart_mag_rms_ft']} fT RMS",
            f"{v['rev_cov_heart_mag_share_brain_pct']} of the measured magnetometer brain variance and "
            f"{v['rev_cov_heart_mag_share_shortfall_pct']} of the model's shortfall",
            f"rises to {v['rev_cov_noheart_mag_log_var_r']} (magnetometers) and {v['rev_cov_noheart_grad_log_var_r']}",
            f"30-40 Hz {v['rev_cov_band_b30_40_mag_amp_ratio']}",
        ]
        for st in ("combined", "mag", "grad"):
            quoted += [f"{v[f'rev_cov_det_{st}_corrected_ratio_3dp']} {v[f'rev_cov_det_{st}_corrected_ci_3dp']}",
                       f"uncorrected {v[f'rev_cov_det_{st}_measured_ratio_3dp']}",
                       f"would give {v[f'rev_cov_surr_null_corrected_{st}_range_3dp'].replace(' to ', '-')}",
                       f"heart removed {v[f'rev_cov_det_{st}_noheart_corrected_ratio_3dp']}"]
        for lobe in ("cingulate", "frontal", "temporal"):
            quoted.append(f"{lobe} {v[f'rev_cov_det_lobe_{lobe}_corrected_ratio']}")
        for arr in ("dense", "matched"):
            for tok, label in (("model306", "both modelled"), ("s4", "S4 Neuromag measured, OPM as modelled"),
                               ("s3", "S3 Neuromag measured, OPM background x k")):
                quoted.append(f"{label} {v[f'rev_cov_opm_{arr}_vs_combined_{tok}_ratio_3dp']} "
                              f"{v[f'rev_cov_opm_{arr}_vs_combined_{tok}_ci_3dp']}")
        for q in quoted:
            self.assertIn(q, text)

    def test_the_noise_sensitivity_configuration_is_the_declared_one(self):
        """A stored run whose configuration differs from configs/g2_noise_sensitivity.toml is refused."""
        with tempfile.TemporaryDirectory() as d:
            for rel in (NS, CFG_NS):
                (Path(d) / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(ROOT / rel, Path(d) / rel)
            cfg = Path(d) / CFG_NS
            cfg.write_text(cfg.read_text().replace("n_boot = 1000", "n_boot = 999"))
            with self.assertRaisesRegex(ValueError, "config.sensitivity differs from configs/g2_noise_sensitivity.toml"):
                R.facts(Path(d))

    def test_children_ranges_cover_the_three_children(self):
        mins = [self.qc["anatomies"][c]["white_to_scalp_used"]["min_mm"] for c in ("childA", "childB", "childC")]
        self.assertEqual(self.F["rev_qc_children_white_scalp_min_mm_range"]["raw"], [min(mins), max(mins)])
        self.assertTrue(math.isclose(self.F["rev_qc_children_white_scalp_min_mm_min"]["raw"], 3.6709052216543094))


# ================================================================================================
# STAGE B2 SECTIONS: rev_cgap_ (helmet fitted at the adult's gap), rev_db_ (per-depth-bin intervals),
# rev_g2_notahead_ (targets where the dense array is not ahead), rev_seed_ (seeds no other module states)
@unittest.skipUnless(HAVE, "revision result files not present")
class TestStageB2Facts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.F = R.facts(ROOT)
        cls.v = {k: f["value"] for k, f in cls.F.items()}
        cls.cg, cls.db = (json.loads((ROOT / p).read_text()) for p in (CGAP, DB))

    def test_cgap_key_values(self):
        """Literals read from results/g3b_constant_gap/g3b_constant_gap_summary.json (comments: the stored value)."""
        expect = {
            "rev_cgap_adult_gap_mm_2dp": "29.55",                                  # target_gaps.gap_matched.gap_mm 29.55017
            "rev_cgap_adult_top_gap_mm_2dp": "28.38",                              # target_gaps.gap_matched_top.gap_mm 28.37951
            "rev_cgap_adult_gap_minus_top_mm": "1.2",                              # 29.55017 - 28.37951
            "rev_cgap_child_c_fitted_gap_mm": "32.8",                              # helmets.childC.gap_matched.median_mm 32.8295
            "rev_cgap_child_c_fitted_binding": "yes",
            "rev_cgap_child_c_fitted_gap_above_target_mm": "3.3",                  # 32.8295 - 29.5502
            "rev_cgap_child_c_fitted_gap_nearest_mm": "18.0",                      # min_mm 18.000: the Dewar spacing
            "rev_cgap_fitted_binding_heads": "child C",
            "rev_cgap_fitted_n_binding": "1",
            "rev_cgap_school_fitted_k": "0.916",                                   # k 0.91617
            "rev_cgap_child_c_fitted_k_at_target": "0.920",                        # k_at_target 0.92025
            "rev_cgap_smaller_fitted_gap_mm_range": "29.6 to 32.8",
            "rev_cgap_all_fixed_gap_mm_range": "28.4 to 40.4",                     # top: adult 28.38 ... infant12mo 40.42
            "rev_cgap_adult_fixed_dense_vs_combined_ib_d": "+1.00",                # D adult/top 0.99739
            "rev_cgap_adult_fitted_dense_vs_combined_ib_d": "+1.28",               # D adult/gap_matched 1.28350
            "rev_cgap_adult_fitted_dense_vs_combined_ib_d_ci": "[+1.08, +1.47]",   # [1.08255, 1.47399]
            "rev_cgap_school_fitted_dense_vs_combined_ib_delta": "+0.03",          # delta_same_rule ... delta.median 0.02563
            "rev_cgap_school_fitted_dense_vs_combined_ib_delta_ci": "[0.00, +0.05]",  # [0.00014, 0.05387]
            "rev_cgap_size2yr_fitted_dense_vs_combined_ib_delta": "+0.02",         # 0.01922
            "rev_cgap_child_a_fitted_dense_vs_combined_ib_delta": "−0.30",         # -0.30461
            "rev_cgap_child_a_fitted_dense_vs_combined_ib_delta_ci": "[−0.48, −0.22]",
            "rev_cgap_child_b_fitted_dense_vs_combined_ib_delta": "−0.39",         # -0.39300
            "rev_cgap_child_b_fitted_dense_vs_combined_ib_delta_ci": "[−0.54, −0.17]",
            "rev_cgap_child_c_fitted_dense_vs_combined_ib_delta": "+0.01",         # 0.00593
            "rev_cgap_child_a_fitted_dense_vs_combined_ib_delta_n": "66",          # n_parcels
            "rev_cgap_school_fitted_dense_vs_combined_ib_delta_n": "6,949",        # common vertices
            "rev_cgap_templates_fitted_dense_vs_combined_ib_delta_range": "+0.04 to +0.15",
            "rev_cgap_child_a_fittedtop_vs_fixed_dense_vs_combined_ib_delta_ci": "[−0.40, +0.06]",
            "rev_cgap_fitted_dense_vs_combined_ib_delta_templates_n_ci_includes_zero": "3",
            "rev_cgap_fitted_dense_vs_combined_ib_delta_children_ci_below_zero_heads": "child A and child B",
            "rev_cgap_fitted_dense_vs_combined_ib_delta_smaller_n_near_zero": "6",
            "rev_cgap_near_zero_db": "0.15",
            "rev_cgap_templates_scaled_fixed_minus_fitted_combined_ib_range": "+0.17 to +0.76",
            "rev_cgap_child_a_fixed_minus_fitted_combined_ib": "+0.29",            # within_head childA/gap_matched 0.28507
            "rev_cgap_child_b_fixed_minus_fitted_combined_ib": "+0.42",            # 0.42436
            "rev_cgap_child_c_fixed_minus_fitted_combined_ib": "−0.12",            # -0.12159
            "rev_cgap_adult_fixed_minus_fitted_combined_ib": "−0.23",              # -0.22858
            "rev_cgap_fixed_minus_fitted_combined_ib_smaller_n_positive": "7",
            "rev_cgap_smaller_fitted_eqsites_dense_vs_combined_ib_delta_range": "+0.16 to +0.49",
            "rev_cgap_child_a_eqsites_n_sites": "155",
            "rev_cgap_school_band_dense_vs_combined_ib_range": "+0.20 to +1.12",   # placement_band min 0.20061, max 1.12202
            "rev_cgap_child_a_band_dense_vs_combined_ib_range": "+0.06 to +0.37",
            "rev_cgap_child_b_band_dense_vs_combined_ib_range": "+0.19 to +0.63",
            "rev_cgap_child_c_band_dense_vs_combined_ib_range": "−0.10 to +0.70",
            "rev_cgap_child_c_band_dense_vs_combined_ib_n": "11",
            "rev_cgap_child_c_band_dense_vs_combined_ib_excluded": "x-5mm",
            "rev_cgap_dense_vs_combined_ib_smaller_band_includes_zero_heads": "child C",
            "rev_cgap_n_boot_primary": "1,000",
            "rev_cgap_n_boot_secondary": "200",
            "rev_cgap_check_max_abs_diff": "7.1 × 10⁻¹⁵",
        }
        for name, value in expect.items():
            self.assertEqual(self.v[name], value, name)

    def test_cgap_raw_values_are_the_stored_ones(self):
        d, F = self.cg, self.F
        e = d["delta_same_rule"]["childA/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]
        self.assertEqual(F["rev_cgap_child_a_fitted_dense_vs_combined_ib_delta"]["raw"], e["median"])
        self.assertEqual(F["rev_cgap_child_a_fitted_dense_vs_combined_ib_delta_ci"]["raw"], e["ci95"])
        tpl = [d["delta_same_rule"][f"{a}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["median"]
               for a in ("infant2yr", "infant18mo", "infant12mo")]
        self.assertEqual(F["rev_cgap_templates_fitted_dense_vs_combined_ib_delta_range"]["raw"], [min(tpl), max(tpl)])
        wh = [d["within_head"][f"{a}/gap_matched/combined/intrinsic+brain"]["median"]
              for a in ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo")]
        self.assertEqual(F["rev_cgap_templates_scaled_fixed_minus_fitted_combined_ib_range"]["raw"], [min(wh), max(wh)])
        self.assertEqual(F["rev_cgap_child_c_fitted_gap_mm"]["raw"], d["helmets"]["childC"]["gap_matched"]["median_mm"])
        self.assertAlmostEqual(F["rev_cgap_adult_fitted_dense_vs_combined_ib_d_ratio"]["raw"],
                               10 ** (d["D"]["adult/gap_matched/opm_dense/combined/intrinsic+brain"]["median"] / 20))
        # the counts are recounted here from the stored values
        near = [a for a in ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
                if abs(d["delta_same_rule"][f"{a}/gap_matched/opm_dense/combined/intrinsic+brain"]["delta"]["median"]) <= 0.15]
        self.assertEqual(F["rev_cgap_fitted_dense_vs_combined_ib_delta_smaller_n_near_zero"]["raw"], len(near))

    def test_cgap_headline_is_the_detailed_sections(self):
        """A headline value that differs from its detailed section is refused (the facts read the sections)."""
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / CGAP).parent.mkdir(parents=True)
            x = json.loads((ROOT / CGAP).read_text())
            x["headline"]["by_anatomy"]["childA"]["delta"]["gap_matched"]["median"] += 0.01
            (Path(d) / CGAP).write_text(json.dumps(x))
            with self.assertRaisesRegex(ValueError, "headline value differs"):
                R.cgap_facts(R.Facts(), Path(d))

    def test_depth_bins(self):
        """Literals of results/g2/g2_depth_bins.json and its medians against results/g2/g2_summary.json."""
        db, v = self.db, self.v
        self.assertEqual(v["rev_db_n_boot"], "1,000")
        self.assertEqual(v["rev_db_seed"], "2026")
        self.assertEqual(v["rev_db_n_parcel_labels"], "70")
        self.assertEqual(v["rev_db_bin_mm"], "5")
        self.assertEqual(v["rev_db_dense_vs_combined_ib_n_bins"], "11")
        self.assertEqual(v["rev_db_dense_vs_combined_ib_10_15_n"], "222")
        summary = json.loads((ROOT / "results/g2/g2_summary.json").read_text())
        for key, c in db["comparisons"].items():
            stored = (summary["log2_ratio_vs_depth"][key] if not key.startswith("peak_field/") else
                      summary["bridge_to_sphere"][key.split("/")[1] + "_ratio_vs_depth"])
            for b, st in zip(c["bins"], stored):
                self.assertEqual((b["lo"], b["hi"], b["n"]), (st["lo"], st["hi"], st["n"]), key)
                if "ratio" not in b:
                    self.assertIsNone(st["median"])
                    continue
                self.assertLessEqual(b["ci95_ratio"][0], b["ratio"])  # the stored median lies inside its interval
                self.assertLessEqual(b["ratio"], b["ci95_ratio"][1])
                if key.startswith("peak_field/"):  # the bin's ratio is the stored median; the CSV reproduces it
                    self.assertEqual(b["ratio"], st["median"])
                    self.assertLessEqual(abs(b["median_ratio_linear_csv"] - st["median"]), b["precision_bound_ratio"] + 1e-12)
                else:
                    self.assertEqual((b["median_log2"], b["ratio"]), (st["median"], 2 ** st["median"]))
                    self.assertLessEqual(abs(b["median_log2_csv"] - st["median"]), b["precision_bound_log2"] + 1e-12)
        b0 = db["comparisons"]["opm_dense/combined/intrinsic+brain"]["bins"][0]
        self.assertEqual(self.F["rev_db_dense_vs_combined_ib_10_15_ci"]["raw"], b0["ci95_ratio"])
        self.assertEqual(v["rev_db_dense_vs_combined_ib_10_15_ratio"], f"{b0['ratio']:.2f}")
        # printed as the existing depth facts print the stored medians (the two must never disagree)
        g2f = _load("report_facts_g12").facts(ROOT)
        for arr in ("dense", "matched"):
            for cond in ("int", "ib", "proj"):
                for lo in range(10, 65, 5):
                    name = f"{arr}_vs_combined_{cond}_{lo}_{lo + 5}_ratio"
                    self.assertEqual(v[f"rev_db_{name}"], g2f[f"g2_depth_{name}"]["value"], name)
        self.assertEqual(v["rev_db_dense_vs_combined_proj_50_55_ratio_3dp"], "0.995")      # stored 0.99493, CSV 0.99541
        self.assertEqual(v["rev_db_dense_vs_combined_proj_bins_ci_below1"], "55–60")
        self.assertEqual(v["rev_db_matched_vs_combined_ib_bins_ci_spans1"], "30–35 and 35–40")
        self.assertEqual(v["rev_db_matched_vs_combined_ib_ci_below1_from_bin"], "40–45")
        self.assertEqual(v["rev_db_dense_vs_combined_ib_n_bins_ci_above1"], "11")
        self.assertEqual(v["rev_db_dense_vs_combined_ib_10_15_ci"], "[1.55, 1.65]")
        self.assertEqual(v["rev_db_check_n_moved"], "5")

    def test_not_ahead_targets(self):
        """Referee 1, minor 13: the 4 of 7,661 targets at which the dense array is not ahead (sensor + brain noise)."""
        v = self.v
        self.assertEqual(v["rev_g2_notahead_n"], "4")
        self.assertEqual(v["rev_g2_notahead_n_csv_undecided"], "1")       # equal at the CSV's 4 decimals
        self.assertEqual(v["rev_g2_notahead_regions"], "left medial orbitofrontal (3) and left medial wall (1)")
        self.assertEqual(v["rev_g2_notahead_depth_mm_range"], "26.5 to 35.1")  # depth_mm 26.53, 30.30, 33.67, 35.08
        self.assertEqual(v["rev_g2_notahead_ratio_range"], "0.956 to 0.997")   # 0.1816/0.19 ... 0.5081/0.5097
        self.assertEqual(v["rev_g2_notahead_1_region"], "lh.unknown")
        self.assertEqual(v["rev_g2_notahead_2_region"], "lh.medialorbitofrontal")
        self.assertEqual(v["rev_g2_notahead_2_depth_mm"], "30.3")
        self.assertEqual(v["rev_g2_notahead_hemispheres"], "left")

    def test_seeds(self):
        """Seeds no other module states (values set in code are read from their line when the facts are built)."""
        expect = {"rev_seed_g2_patch_boot": "1", "rev_seed_g2_convergence": "2", "rev_seed_g2_band": "2026",
                  "rev_seed_g2_hse": "2026", "rev_seed_g2_near_mesh": "11", "rev_seed_bem_sphere": "5", "rev_seed_g1b": "2016",
                  "rev_seed_g1b_noise": "7", "rev_seed_g1c": "2009", "rev_seed_g3b_patches": "2027", "rev_seed_g3b_useful": "0",
                  "rev_seed_cgap": "2026", "rev_seed_signflip_mc": "0", "rev_seed_confirm_replicates": "5"}
        for name, value in expect.items():
            self.assertEqual(self.v[name], value, name)
        self.assertFalse([k for k in self.F if k.startswith("rev_boot_")])  # every resample count is a fact elsewhere


class TestCodeConstants(unittest.TestCase):
    def test_code_int_reads_the_line_and_refuses_a_missing_one(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "x.py").write_text("def a():\n    rng = np.random.default_rng(1)\n\n\ndef b():\n    rng = np.random.default_rng(2)\n")
            self.assertEqual(R.code_int(d, "x.py", r"default_rng\((\d+)\)"), (1, 2))
            self.assertEqual(R.code_int(d, "x.py", r"default_rng\((\d+)\)", after=r"^def b\("), (2, 6))
            with self.assertRaises(ValueError):
                R.code_int(d, "x.py", r"default_rng\((\d+)\)", after=r"^def c\(")
            with self.assertRaises(ValueError):
                R.code_int(d, "x.py", r"SEED = (\d+)")


if __name__ == "__main__":
    unittest.main()
