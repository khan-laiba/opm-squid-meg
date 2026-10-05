"""Full-precision depth companions of the per-target tables (scripts/export_target_precision.py): the
half-open bin rule and the exact decimal form, then the exported companions against their stored tables
(keys, row order, printed depth, cortical flag, bins) and against every depth-bin and orientation-stratum
count stored in the summaries. Reads result files only (no anatomy is loaded)."""
import importlib.util
import tomllib
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


X = _load("export_target_precision")


class TestBins(unittest.TestCase):
    def test_half_open_bins(self):
        edges = (0.0, 10.0, 15.0, 90.0)
        x = np.array([-1.0, 0.0, np.nextafter(10.0, 0.0), 10.0, 14.99, 15.0, 89.99, 90.0, np.nan])
        want = ["", "0-10", "0-10", "10-15", "10-15", "15-90", "15-90", "", ""]
        self.assertEqual(list(X.bin_labels(x, edges)), want)
        self.assertEqual(list(X.labels_by_rule(x, edges)), want)

    def test_exact_form_reads_back(self):
        for x in (9.999935645831501, np.nextafter(10.0, 0.0), 0.1 + 0.2, 65.42377386317769, 7.0):
            self.assertEqual(float(X.exact(x)), float(x))
        self.assertEqual(f"{9.999935645831501:.2f}", "10.00")  # printed to 0.01 mm, a target below 10 mm reads as 10.00

    def test_contiguous_bands(self):
        self.assertEqual(X.edges_of(((10.0, 20.0), (20.0, 45.0))), (10.0, 20.0, 45.0))
        with self.assertRaises(ValueError):
            X.edges_of(((10.0, 20.0), (25.0, 45.0)))


SPECS = X.specs()


@unittest.skipUnless(all(s.companion.exists() for s in SPECS), "companions not exported (scripts/export_target_precision.py)")
class TestCompanions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.strata = X.strata_edges(tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text()))
        cls.tables = {s: X.read_table(s.table) for s in SPECS}
        cls.comps = {s: X.read_table(s.companion, s.columns) for s in SPECS}
        cls.summaries = X.load_summaries()
        cls.adult = next(s for s in SPECS if s.family == "adult")

    def test_header_lines(self):
        for s in SPECS:
            c = self.comps[s]
            self.assertEqual(c["_header"], list(s.columns), s.companion.name)
            self.assertRegex(c["_comments"][0], r"^# .+ \| commit [0-9a-f]{7,40}(\+dirty)?$")
            self.assertIn(f"{s.table.name} (commit {self.tables[s]['_commit']})", c["_comments"][0])
            self.assertEqual(len(c["_comments"]), 2)  # status and legend

    def test_rows_printed_depth_and_bins(self):
        chk = X.Checks()
        for s in SPECS:
            X.check_against_table(chk, s, self.comps[s], self.tables[s], self.strata)
        self.assertEqual(chk.failed, [])
        self.assertEqual(len(chk.n), len(SPECS))

    def test_patch_table_has_the_same_rows(self):
        patch = X.read_table(X.PATCH_TABLE, ("hemi", "vertno", "depth_mm"))
        t = self.tables[self.adult]
        self.assertEqual((patch["hemi"], patch["vertno"], patch["depth_mm"]), (t["hemi"], t["vertno"], t["depth_mm"]))

    def test_adult_counts(self):
        chk = X.Checks()
        for name in X.ADULT_SUMMARIES:
            X.check_adult_counts(chk, name, self.summaries[name], self.comps[self.adult])
        self.assertEqual(chk.failed, [])
        for name in X.ADULT_SUMMARIES:  # every summary had bin counts to check
            self.assertGreater(chk.n.get(X.SUMMARIES[name].name, 0), 0, name)

    def test_smaller_head_counts(self):
        chk = X.Checks()
        for fam, name in (("pediatric", "g3b_summary"), ("constant_gap", "constant_gap_summary")):
            X.check_strata_counts(chk, name, self.summaries[name], {s.anatomy: self.comps[s] for s in SPECS if s.family == fam},
                                  self.strata)
        self.assertEqual(chk.failed, [])
        for name in ("g3b_summary", "constant_gap_summary"):
            self.assertGreater(chk.n.get(X.SUMMARIES[name].name, 0), 1000, name)

    def test_child_b_below_10mm(self):
        """Referee 2's example: the count of g3b_summary.json follows from the companion (the rounded table prints
        some of these targets as 10.00)."""
        s = next(x for x in SPECS if x.family == "pediatric" and x.anatomy == "childB")
        c = self.comps[s]
        row = self.summaries["g3b_summary"]["comparisons"]["childB/opm_dense/combined/intrinsic+brain/detect"]["delta_by_depth"][0]
        self.assertEqual((row["lo"], row["hi"]), (0.0, 10.0))
        self.assertEqual(sum(k == "1" and b == "0-10" for k, b in zip(c["cortical"], c["depth_stratum"], strict=True)), row["n_child"])
        self.assertEqual(sum(k == "1" and float(d) < 10.0 for k, d in zip(c["cortical"], c["depth_mm"], strict=True)), row["n_child"])


if __name__ == "__main__":
    unittest.main()
