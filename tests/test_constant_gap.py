"""Constant-gap helmet control (scripts/study_g3b_constant_gap.py): the gap-matched helmet, the
placement band read from the stored G3B summary and the strata-free comparison."""
import importlib.util
import json
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import mne

from opmsquid import anatomy, neuromag, pediatric as P

ROOT = Path(__file__).resolve().parents[1]
G3B_SUMMARY = ROOT / "results/g3b/g3b_summary.json"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load("study_g3b_constant_gap")


class TestGapMatchedHelmet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.info = neuromag.load_info("T3")
        base = np.array(cls.info["dev_head_t"]["trans"])
        cls.adult = anatomy.load_sample()
        cls.pl = P.placements(cls.info, cls.adult, base)
        cls.pose = cls.pl["x-centred"]["pose"]
        cls.cfx = P.counterfactual_helmet(cls.info, cls.adult, cls.pose, 1.0)
        cls.child = anatomy.scaled(cls.adult, 0.85)
        cls.child_pose = P.placements(cls.info, cls.child, base)["x-centred"]["pose"]

    def test_adult_under_its_own_rule_is_its_own_helmet(self):
        g = S.gap_matched_helmet(self.info, self.adult, self.pose, self.cfx["median_dist_mm"] * 1e-3)
        self.assertEqual(g["k"], 1.0)
        self.assertFalse(g["clearance_binding"])
        for a, b in zip(g["info"]["chs"], self.cfx["info"]["chs"]):
            np.testing.assert_array_equal(a["loc"], b["loc"])  # bitwise: the cached lead fields are reused

    def test_target_gap_is_reached(self):
        target = self.pl["top"]["median_dist_mm"] * 1e-3  # the adult's top-contact gap
        for sub, pose in ((self.adult, self.pose), (self.child, self.child_pose)):
            g = S.gap_matched_helmet(self.info, sub, pose, target)
            self.assertFalse(g["clearance_binding"])
            self.assertTrue(g["feasible"])
            self.assertAlmostEqual(g["median_dist_mm"], target * 1e3, delta=1e-3)
            self.assertGreaterEqual(g["median_dist_mm"], target * 1e3)
            np.testing.assert_allclose(P.HelmetFit(g["info"], sub).distances(pose).min() * 1e3, g["min_dist_mm"])
        g = S.gap_matched_helmet(self.info, self.child, self.child_pose, target)
        self.assertTrue(0.85 < g["k"] < 1.0)  # between the helmet scaled with this head (0.85) and the adult's own (1)

    def test_clearance_binds_before_a_small_target(self):
        g = S.gap_matched_helmet(self.info, self.child, self.child_pose, 0.020)  # 20-mm median gap: some coil within 18 mm
        self.assertTrue(g["clearance_binding"])
        self.assertGreater(g["k"], g["k_at_target"])
        self.assertTrue(g["feasible"])
        self.assertGreaterEqual(g["min_dist_mm"], 18.0)
        self.assertGreater(g["median_dist_mm"], 20.0)
        centre = np.linalg.inv(self.child_pose)[:3, 3]
        smaller = P.HelmetFit(P.scaled_helmet(self.info, g["k"] - 1e-5, centre), self.child).describe(self.child_pose)
        self.assertFalse(smaller["feasible"])  # the smallest feasible factor
        other = S.gap_matched_helmet(self.info, self.child, self.child_pose, 0.022)  # another target below the same limit
        self.assertTrue(other["clearance_binding"])
        self.assertEqual(other["k"], g["k"])  # the same helmet (the same cached lead fields)


@unittest.skipUnless(G3B_SUMMARY.exists(), "G3B summary not present")
class TestPlacementBand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g3b = json.loads(G3B_SUMMARY.read_text())
        cls.band = S.placement_band(cls.g3b)

    def test_band_is_the_difference_of_the_stored_medians(self):
        plc = self.g3b["placement_D"]
        for c in S.CHILDREN:
            for ref in S.REFS:
                b = self.band[f"{c}/{ref}/intrinsic+brain"]
                bad = set(self.g3b["infeasible_placements"][c]) | set(self.g3b["infeasible_placements"]["adult"])
                self.assertEqual(b["n_placements"], len([p for p in S.FAMILY if p not in bad]))
                for p, v in b["delta_by_placement"].items():
                    self.assertAlmostEqual(v, plc[f"{c}/{p}/{ref}/intrinsic+brain"]["median"] - plc[f"adult/{p}/{ref}/intrinsic+brain"]["median"])
                self.assertLessEqual(b["min"], b["top"])
                self.assertLessEqual(b["top"], b["max"])
                self.assertLessEqual(b["crossed"]["lo"], b["min"])
                self.assertGreaterEqual(b["crossed"]["hi"], b["max"])

    def test_infeasible_placement_left_out(self):
        self.assertIn("x-5mm", self.g3b["infeasible_placements"]["childC"])
        b = self.band["childC/combined/intrinsic+brain"]
        self.assertNotIn("x-5mm", b["delta_by_placement"])
        self.assertEqual(b["excluded_infeasible"], ["x-5mm"])


class TestLeanComparison(unittest.TestCase):
    def test_lean_configuration_keeps_the_estimators_and_drops_the_strata(self):
        rng = np.random.default_rng(0)
        cfg = dict(strata=dict(depth_edges_mm=[0.0, 20.0, 90.0], orientation_edges_deg=[0.0, 45.0, 90.1], min_n=10, n_boot=50))

        def anat(key, n):
            region = np.array([f"lh.p{i % 12}" for i in range(n)])
            src = SimpleNamespace(target=np.arange(n), region=region, lobe=np.array(["frontal"] * n), depth_mm=rng.uniform(5, 80, n),
                                  orientation_deg=rng.uniform(0, 90, n))
            return SimpleNamespace(key=key, src=src, weights=rng.uniform(0.5, 1.5, n), cortical=np.ones(n, bool), nt=n)

        for key in ("school", "childA"):  # vertex-wise and parcel-matched estimators
            c, a = anat(key, 600), anat("adult", 600)
            xc, xa = rng.normal(1.0, 1.0, 600), rng.normal(0.5, 1.0, 600)
            full = S.G3B.compare(c, a, xc, xa, cfg, None, 0)
            lean = S.G3B.compare(c, a, xc, xa, S.lean(cfg), None, 0)
            for q in ("d_child", "d_adult", "delta"):
                self.assertEqual(lean[q]["median"], full[q]["median"])
            self.assertEqual(lean["delta_by_depth"], [])
            self.assertTrue(full["delta_by_depth"])
            self.assertNotIn("delta_by_depth", S.trim(lean, False))


if __name__ == "__main__":
    unittest.main()
