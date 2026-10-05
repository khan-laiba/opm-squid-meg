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
NEEDED = [NS, COV, QC, CFG_NS, "configs/school_subjects_qc_manifest.json"]
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
            self.assertRegex(name, r"^rev_(ns|cov|qc)_[a-z0-9_]+$")
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
            "rev_qc_child_b_targets_lt10mm_n": "245",
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


if __name__ == "__main__":
    unittest.main()
