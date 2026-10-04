"""Report facts for G4 (scripts/report_facts_g4.py): the formatting helpers, and every fact built
from the committed results/g4 files has a value and an existing source file; key values are
checked against literals read from the result files."""
import importlib.util
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
G4 = ROOT / "results" / "g4"
NEEDED = ["g4_pediatric_comparison.json", "g4_matched_rate.json", "g4_g2_consistency.json", "g4_fit_failures.json",
          "g4_motion_summary.json", "g4_localization_summary.json", "g4_adult_events.csv", "g4_localization_events.csv"]
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
NEEDED += [f"g4_{lab}_summary.json" for lab in LABELS] + [f"g4_{lab}_events.csv" for lab in LABELS]
NEEDED += [f"g4_localization_{lab}_{kind}" for lab in LABELS[1:] for kind in ("summary.json", "events.csv")]
HAVE = all((G4 / f).exists() for f in NEEDED) and (ROOT / "results" / "g3b" / "g3b_summary.json").exists()


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("report_facts_g4")


class TestFormats(unittest.TestCase):
    def test_numbers(self):
        self.assertEqual(M.pval(6.1e-5), "<0.0001")
        self.assertEqual(M.pval(0.00390625), "0.0039")
        self.assertEqual(M.pval(0.0224609375), "0.022")
        self.assertEqual(M.pval(0.5), "0.50")
        self.assertEqual(M.pval(1.0), "1.0")
        self.assertEqual(M.pval(0.05 / 24), "0.0021")
        self.assertEqual(M.db(0.44), "+0.44")
        self.assertEqual(M.db(-0.38), "−0.38")
        self.assertEqual(M.db(-0.004), "0.00")  # no signed zero
        self.assertEqual(M.ratio(1.3642996), "1.36")
        self.assertEqual(M.pct(27 / 456), "5.9%")
        self.assertEqual(M.pct(77 / 456), "17%")
        self.assertEqual(M.pct(0.0996), "10%")
        self.assertEqual(M.integer(1394.94), "1,395")
        self.assertEqual(M.count(1728), "1,728")
        self.assertEqual(M.dmm(-6.97), "−7.0")
        self.assertEqual(M.sig(0.016823), "0.017")
        self.assertEqual(M.sig(3.16597), "3.2")
        self.assertEqual(M.const([10.0, 20.0, 40.0]), "10, 20 and 40")
        self.assertEqual(M.const(-60.0), "−60")

    def test_intervals_and_ranges(self):
        self.assertEqual(M.iv(241.46, None, M.integer), "[241, open]")
        self.assertEqual(M.iv(1.1016, 1.5756, M.ratio), "[1.10, 1.58]")
        self.assertEqual(M.span([-1.7165, -1.6444, -1.6663], M.db), "−1.72 to −1.64")
        self.assertEqual(M.span([14.86, 15.06], lambda x: M.sig(x)), "15")  # equal ends print once

    def test_holm(self):
        adj = M._holm({"a": 0.01, "b": 0.04, "c": 0.03})
        self.assertAlmostEqual(adj["a"], 0.03)
        self.assertAlmostEqual(adj["c"], 0.06)
        self.assertAlmostEqual(adj["b"], 0.06)  # monotone


@unittest.skipUnless(HAVE, "G4 result files not available")
class TestG4Facts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f = M.facts(ROOT)

    def test_every_fact_has_a_value_and_an_existing_source(self):
        self.assertGreater(len(self.f), 1000)
        name_ok = re.compile(r"^(g4|loc|mot)_[a-z0-9_]+$")
        for name, fact in self.f.items():
            self.assertRegex(name, name_ok)
            self.assertEqual(set(fact), {"value", "raw", "source"}, name)
            self.assertNotEqual(str(fact["value"]).strip(), "", name)
            files, sep, how = fact["source"].partition(" :: ")
            self.assertTrue(sep and how.strip(), f"{name}: source lacks ' :: <key path>'")
            paths = [p for p in files.replace(",", " ").split() if p.startswith(("results/", "configs/", "docs/"))]
            self.assertTrue(paths, f"{name}: no source file")
            for p in paths:
                self.assertTrue((ROOT / p).exists(), f"{name}: {p} not found")
        # the combiner's own check agrees
        rf = _load("report_facts")
        self.assertEqual(rf.check(self.f, ROOT), [])

    def test_no_ascii_minus_in_numbers(self):
        for name, fact in self.f.items():
            self.assertNotRegex(str(fact["value"]), r"(^|[\s\[(])-\d", name)

    def test_key_values_against_the_result_files(self):
        f = {k: v["value"] for k, v in self.f.items()}
        adult = json.loads((G4 / "g4_adult_summary.json").read_text())
        s50 = adult["detectors"]["squid/combined"]["strength_for_50pct_nAm"]["practical@1/depth0"]
        self.assertAlmostEqual(s50["value"], 46.5051762592345)  # the stored literal behind the fact
        self.assertEqual(f["g4_adult_s50_combined_practical_10_20mm"], "47")
        self.assertEqual(f["g4_adult_s50_combined_practical_10_20mm_ci"], "[36, 63]")
        self.assertEqual(f["g4_adult_s50_combined_practical_45_70mm_ci"], "[241, open]")
        self.assertEqual(f["g4_adult_dense_vs_combined_practical_10_20mm_ratio"], "1.36")
        self.assertEqual(f["g4_adult_dense_vs_combined_practical_10_20mm_ratio_ci"], "[1.10, 1.58]")
        self.assertEqual(f["g4_adult_dense_vs_combined_practical_10_20mm_p"], "0.0039")
        self.assertEqual((f["g4_adult_dense_vs_combined_practical_10_20mm_locs_opm"],
                          f["g4_adult_dense_vs_combined_practical_10_20mm_locs_squid"]), ("11", "1"))
        self.assertEqual(f["g4_adult_dense_vs_combined_practical_10_20mm_reduction_pct"], "27%")
        self.assertEqual(f["g4_adult_dense_vs_combined_oracle_10_20mm_p"], "<0.0001")
        self.assertEqual(f["g4_childa_dense_vs_combined_practical_10_20mm_ratio_ci"], "[0.93, 1.63]")
        self.assertEqual(f["g4_infant18mo_dense_vs_combined_practical_45_70mm_ratio"], "> 1.13")
        self.assertEqual(f["g4_holm_10_20mm_max_p_adj"], "0.022")
        self.assertEqual(f["g4_holm_10_20mm_n_pass"], "9")
        self.assertEqual(f["g4_heldout_rate_range"], "0.40 to 1.55")
        self.assertEqual(f["g4_adult_sens40_matchedrate_dense"], "0.36")
        self.assertEqual(f["g4_adult_matchedrate_dense_vs_combined_10_20mm_ratio_ci"], "[1.12, 1.59]")
        self.assertEqual((f["g4_adult_events_10_20mm_dense_only"], f["g4_adult_events_10_20mm_combined_only"],
                          f["g4_adult_events_10_20mm_n"]), ("19", "1", "324"))
        self.assertEqual(f["g4_n_events_per_anatomy"], "1,728")
        self.assertEqual(f["g4_bonferroni_24_threshold"], "0.0021")
        self.assertEqual((f["loc_family_n"], f["loc_family_n_sig"], f["loc_survivors_bonferroni_20"]), ("360", "65", "8"))
        self.assertEqual([f[f"loc_family_n_sig_{e}"] for e in ("dspm", "dspm_mne", "ecd", "joint_dspm", "joint_ecd")],
                         ["21", "27", "6", "11", "0"])
        self.assertEqual(f["loc_adult_coreg_displacement_mm"], "2.9")
        self.assertEqual(f["loc_adult_dense_vs_combined_patch320_dspm_diff_mm"], "−5.6")
        self.assertEqual(f["loc_failures_combined_ecd_gross_rate_detected"], "5.9%")
        self.assertEqual(f["loc_inverse_grid_range"], "2,114 to 3,821")
        self.assertEqual(f["mot_squid_down10mm_mismatched_db_range"], "−1.72 to −1.64")
        self.assertEqual(f["mot_adult_static_d_eightterm_db"], "+0.67")
        self.assertEqual(f["mot_eightterm_cal1_gradient_thr_loss1db_deg_range"], "3.2 to 3.9")
        self.assertEqual(f["mot_tc_artefact_none_ft"], "1,395")
        self.assertEqual(f["mot_squid_n_infeasible"], "18")

    def test_raw_values_are_the_stored_ones(self):
        mr = json.loads((G4 / "g4_matched_rate.json").read_text())
        r = mr["anatomies"]["childA"]["matched"]["opm_dense/opm_vs_squid/combined/depth0"]["s50_ratio_squid_over_opm"]["value"]
        self.assertEqual(self.f["g4_childa_matchedrate_dense_vs_combined_10_20mm_ratio"]["raw"], r)
        ff = json.loads((G4 / "g4_fit_failures.json").read_text())
        self.assertEqual(self.f["loc_failures_dense_ecd_gross_detected"]["raw"], ff["totals"]["ecd_gross_detected_dense"])
        mo = json.loads((G4 / "g4_motion_summary.json").read_text())
        self.assertEqual(self.f["mot_tc_peak_field_change_max_pt"]["raw"], mo["timecourse"]["peak_field_change_pT"]["max"])


if __name__ == "__main__":
    unittest.main()
