import importlib.util
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEEDED = ["results/g1a/g1a_benchmark.json", "results/g1a/fig3/figure3_summary.json", "results/g1b/g1b_summary.json",
          "results/g1c/g1c_summary.json", "results/g2/g2_summary.json", "results/g2/g2_band_sensitivity.json",
          "results/g2/g2_targets.csv", "results/g2/g2_patch_targets.csv", "results/g2/head_surface_effect.json",
          "results/g2/near_mesh_check.json", "results/g2/bem_skin_refinement.json",
          "results/g2/bem_sphere_check.json"]
PREFIXES = ("g1a_", "g1b_", "g1c_", "g2_", "reg_", "lit_")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(all((ROOT / p).exists() for p in NEEDED), "G1/G2 result files not available")
class TestReportFactsG12(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.facts = _load("report_facts_g12").facts(ROOT)

    def test_every_fact_has_a_value_and_an_existing_source(self):
        self.assertGreater(len(self.facts), 1000)
        for name, f in self.facts.items():
            self.assertEqual(set(f), {"value", "raw", "source"}, name)
            self.assertTrue(str(f["value"]).strip(), name)
            files = [p for p in f["source"].split("::")[0].replace(",", " ").split()
                     if p.startswith(("results/", "configs/", "docs/"))]
            self.assertTrue(files, f"{name}: no file in {f['source']!r}")
            for p in files:
                self.assertTrue((ROOT / p).exists(), f"{name}: {p} not found")
        self.assertEqual(_load("report_facts").check(self.facts, ROOT), [])

    def test_names_and_formats(self):
        for name, f in self.facts.items():
            self.assertRegex(name, r"^[a-z0-9_]+$")
            self.assertTrue(name.startswith(PREFIXES), name)
            self.assertIsNone(re.search(r"(^|[\s\[(])-\d", f["value"]), f"{name}: ASCII minus in {f['value']!r}")
            if name.endswith(("_ci", "_ci_3dp")):
                self.assertRegex(f["value"], r"^\[\S+, \S+\]$", name)
            if name.endswith("_range"):
                self.assertIn(" to ", f["value"], name)
            json.dumps(f["raw"])  # plain JSON types only

    def test_key_values(self):
        expect = {  # literals read from the result files (ratios = 2**median_log2 where stored as log2)
            "g1a_deq_eta3_mm_3dp": "27.665", "g1a_eta_crossing_range": "1.683 to 5.309", "g1a_fig4b_deq_drawn_mm": "35.02",
            "g1b_dipole_mag_r_3dp": "0.961", "g1b_dipole_mag_level_ratio_3dp": "0.736", "g1b_sign_agree_dipole_n": "53",
            "g1b_opm_vs_mag_row_20_25_ratio": "1.37", "g1b_dipole_mag_within_class": "9.7%",
            "g1c_focal_model_pooled_median_db": "−22.13", "g1c_opm30ft_minus_mag_focal_db": "−1.01",
            "g1c_patch_gain_model_pooled_db": "+5.41", "g1c_source_sd_4000_nam": "1.77",
            "g2_n_targets": "7,661", "g2_n_background": "1,755",
            "g2_dense_vs_combined_ib_ratio": "1.14", "g2_dense_vs_combined_ib_ci": "[1.12, 1.17]",
            "g2_dense_vs_combined_ib_ci_3dp": "[1.121, 1.168]", "g2_matched_vs_combined_ib_ratio_3dp": "1.009",
            "g2_matched_vs_combined_ib_share": "56%", "g2_dense_vs_combined_ib_n_squid_better": "4",
            "g2_matched_vs_combined_int_ratio_3dp": "0.531", "g2_dense_vs_combined_proj_ratio_3dp": "1.115",
            "g2_peak_dense_vs_combined_ib_ratio": "0.90", "g2_peak_dense_vs_mag_ib_ratio": "1.23",
            "g2_depth_dense_vs_combined_ib_10_15_ratio": "1.60", "g2_depth_10_15_n": "222",
            "g2_depth_dense_vs_combined_proj_below1_from_mm": "45 to 50", "g2_brainmult_dense_vs_combined_range": "1.51 to 1.68",
            "g2_breakeven_dense_vs_combined_asd": "11.13", "g2_bridge_deq_matched_eta3_mm": "29.6",
            "g2_bridge_deq_sphere_real_eta3_mm_2dp": "28.95", "g2_measured_brain_mag_ft": "262", "g2_model_brain_mag_ft": "192",
            "g2_noise_rms_dense_vs_mag_ratio": "2.30", "g2_band_1_40hz_intrinsic_share_dense_pct": "3.0%",
            "g2_band_1_40hz_matched_vs_combined_ib_ci_3dp": "[1.007, 1.011]", "g2_sens_asd30ft_dense_vs_combined_ib_ratio_3dp": "1.008",
            "g2_joint_gap3mm_asd30ft_dense_vs_combined_ib_ci_3dp": "[0.956, 0.978]", "g2_plugin10s_samples": "780",
            "g2_lobe_temporal_dense_vs_combined_ib_ratio": "1.16", "g2_hse_v3_dense_vs_combined_proj_ratio_3dp": "1.048",
            "reg_superiortemporal_n": "281", "reg_parahippocampal_depth_mm": "28.7",
            "reg_parahippocampal_matched_vs_combined_proj_ratio": "0.78", "reg_precentral_peakfield_dense_vs_mag_ratio": "3.77",
            "lit_sef_noise_ratio_closest": "4.6",
        }
        for name, value in expect.items():
            self.assertEqual(self.facts[name]["value"], value, name)

    def test_ranges_come_from_unrounded_values(self):
        f = self.facts["g2_depth_dense_vs_combined_ib_35_60_range_3dp"]
        lo, hi = f["raw"]
        self.assertLess(lo, hi)
        self.assertEqual(self.facts["g2_depth_dense_vs_combined_ib_35_60_min"]["raw"], lo)


if __name__ == "__main__":
    unittest.main()
