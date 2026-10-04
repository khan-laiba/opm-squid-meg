"""Report facts of G3A, G3B and the regions by head (scripts/report_facts_g3.py) against the stored results."""
import importlib.util
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEEDED = [ROOT / "results/g3a/g3a_size_benchmark.json", ROOT / "results/g3b/g3b_summary.json",
          ROOT / "results/g3b/school_anatomy_checks.json", ROOT / "results/g3b/child_bem_validation.json",
          ROOT / "results/g3b/school_subjects_preparation.json", ROOT / "configs/g2_adult.toml", ROOT / "docs/methods.md"] + [
    ROOT / f"results/g3b/g3b_targets_{a}.csv" for a in ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo",
                                                        "childA", "childB", "childC")]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(all(p.exists() for p in NEEDED), "G3A/G3B result files not present")
class TestReportFactsG3(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.F = _load("report_facts_g3").facts(ROOT)

    def test_every_fact_has_value_raw_and_existing_source(self):
        self.assertGreater(len(self.F), 1000)
        for name, f in self.F.items():
            self.assertRegex(name, r"^(g3a|g3b|regh)(_[a-z0-9]+)+$")
            self.assertEqual(set(f), {"value", "raw", "source"}, name)
            self.assertTrue(str(f["value"]).strip(), name)
            self.assertIsNotNone(f["raw"], name)
            files = [p for p in f["source"].split("::")[0].replace(",", " ").split()]
            self.assertTrue(files, name)
            for p in files:
                self.assertTrue(p.startswith(("results/", "configs/", "docs/")), f"{name}: {p}")
                self.assertTrue((ROOT / p).exists(), f"{name}: {p}")
        self.assertEqual(_load("report_facts").check({k: dict(f) for k, f in self.F.items()}, ROOT), [])

    def test_number_formats(self):
        """dB with U+2212 and no signed zero; no ASCII hyphen in any numeric value; ranges 'lo to hi'."""
        for name, f in self.F.items():
            v = f["value"]
            if re.search(r"\d", v) and not re.search(r"[A-Za-z]", v):
                self.assertNotIn("-", v, name)
            self.assertNotRegex(v, "[+\u2212]0\\.0+(?![0-9])", name)
            if name.endswith("_range"):
                self.assertRegex(v, r"^\S+ to \S+$", name)

    def test_key_values(self):
        """Literals read from the stored results (g3a_size_benchmark.json, g3b_summary.json, school_anatomy_checks.json, CSVs)."""
        expect = {
            "g3a_adult_follow_eta3_deq": "27.7",                    # size_following['Adult']['eta3'].d_eq_exact_mm 27.665
            "g3a_newborn_follow_eta3_deq": "30.8",                  # 30.814
            "g3a_newborn_follow_eta3_norm": "50%",                  # normalized_exact_pct 49.61
            "g3a_newborn_follow_eta3_norm_1dp": "49.6%",
            "g3a_adult_follow_eta3_volume_1dp": "40.4%",            # volume_fraction_pct 40.37
            "g3a_newborn_fixed_eta0_1dp": "8.7",                    # fixed_adult_shell['Infant (newborn)'].eta0 8.6726
            "g3a_child8yr_fixed_eta3_deq": "48.0",                  # 48.032
            "g3b_adult_dense_vs_combined_ib_d": "+1.00",            # D_median_dB 0.99739
            "g3b_adult_ofc_cm": "58.6",                             # anatomies.adult.head_size.ofc_mm 586.10
            "g3b_adult_n_targets": "7,661",
            "g3b_child_a_cortex_area_cm2": "1,852",                 # 1851.75
            "g3b_school_dense_vs_combined_ib_delta": "+0.44",       # comparisons[...].delta.median 0.44193
            "g3b_school_dense_vs_combined_ib_delta_ci": "[+0.34, +0.56]",
            "g3b_child_c_dense_vs_combined_ib_delta_ci": "[\u22120.07, +0.34]",
            "g3b_infant2yr_dense_vs_combined_ib_delta_n": "66",
            "g3b_adult_top_gap": "28.4",                            # sensor_distances.adult['squid:top'].median_mm 28.380
            "g3b_child_c_cfx_gap": "33.1",
            "g3b_child_c_cfx_k": "0.943",
            "g3b_adult_cfx_dense_vs_combined_ib_d": "+1.28",        # placement_D 1.28350
            "g3b_child_c_cfx_dense_vs_combined_ib_delta": "+0.04",  # delta_other_placements 0.03916
            "g3b_child_a_cfx_dense_vs_combined_ib_delta_ci": "[\u22120.84, \u22120.68]",
            "g3b_school_eqcount_dense_vs_combined_ib_d_adult_sub": "+0.73",
            "g3b_child_b_shallow_targets_lt10mm": "247",
            "g3b_child_b_pooled_diff": "+0.52",                     # template_depth_checks.childB 0.51730
            "g3b_child_b_reweighted_diff": "+0.05",
            "g3b_child_b_noise_dense_brain": "751",                 # noise_rms_median.childB['opm_dense/mag'].brain 751.31
            "g3b_adult_useful_ib_100nam_both": "0.660",
            "g3b_child_a_white_to_scalp_closest_2dp": "3.70",       # children_thinnest_layers.childA 3.6962
            "g3b_n_infeasible_placements": "1",
            "regh_adult_top_vs_combined_parahippocampal": "+0.74",
            "regh_adult_parahippocampal_n": "92",
        }
        for name, value in expect.items():
            self.assertEqual(self.F[name]["value"], value, name)

    def test_derived_counterfactual_and_penalty(self):
        s = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())
        d_c = s["delta_other_placements"]["childC/counterfactual_x-centred_vs_adult_counterfactual_x-centred/combined"]["d_child"]["median"]
        d_top = s["placement_D"]["adult/top/combined/intrinsic+brain"]["median"]
        self.assertAlmostEqual(self.F["g3b_child_c_cfx_dense_vs_combined_ib_d_minus_adult_top"]["raw"], d_c - d_top, places=12)
        pen = (s["comparisons"]["school/opm_dense/combined/intrinsic+brain/detect"]["d_adult"]["median"]
               - s["channel_count_control"]["school/combined"]["d_adult_subsampled"]["median"])
        self.assertAlmostEqual(self.F["g3b_school_site_penalty_vs_combined"]["raw"], pen, places=12)
        self.assertEqual(self.F["g3b_smaller_site_penalty_vs_combined_range"]["value"], "0.27 to 0.53")

    def test_csv_regions_reproduce_stored_lobes(self):
        """The CSV recomputation behind every lobe fact agrees with placement_D by_lobe within 0.001 dB."""
        n = 0
        for name, f in self.F.items():
            m = re.search(r"gives ([+-]\d+\.\d+) dB", f["source"]) if name.startswith("regh_") else None
            if m:
                n += 1
                self.assertAlmostEqual(float(m.group(1)), f["raw"], delta=0.001, msg=name)
        self.assertEqual(n, 9 * 2 * 3 * 6)  # anatomies x (top, counterfactual) x comparators x lobes


if __name__ == "__main__":
    unittest.main()
