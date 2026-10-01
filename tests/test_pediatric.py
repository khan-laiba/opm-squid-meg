"""G3B geometry: scaled size-only controls, the infant template loader, helmet placements, the
counterfactual helmet and the area-weighted summaries."""
import unittest

import numpy as np
import mne

from opmsquid import anatomy, neuromag, paths, pediatric as P


class TestScaledSubject(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.adult = anatomy.load_sample()
        cls.s = 0.9
        cls.child = anatomy.scaled(cls.adult, cls.s)

    def test_head_frame_is_the_adult_one_scaled_about_its_origin(self):
        np.testing.assert_allclose(P.scalp_head_frame(self.child), self.s * P.scalp_head_frame(self.adult), atol=1e-12)
        fa = anatomy._sample_fiducials()
        for k, v in self.child.fiducials.items():
            np.testing.assert_allclose(v, self.s * fa[k], atol=1e-12)

    def test_transform_keeps_rotation_and_scales_translation(self):
        ta, tc = self.adult.trans["trans"], self.child.trans["trans"]
        np.testing.assert_allclose(tc[:3, :3], ta[:3, :3])
        np.testing.assert_allclose(tc[:3, 3], self.s * ta[:3, 3])

    def test_refined_head_surface_follows_the_geometry_not_the_name(self):
        same_name = anatomy.scaled(self.adult, self.s, name="sample")  # the cache key must include the geometry
        ba, bc = self.adult.bem_model((0.3, 0.006, 0.3)), same_name.bem_model((0.3, 0.006, 0.3))
        self.assertEqual(len(bc[0]["tris"]), 20480)
        np.testing.assert_allclose(bc[0]["rr"], self.s * ba[0]["rr"], atol=1e-12)
        np.testing.assert_allclose(bc[2]["rr"], self.s * ba[2]["rr"], atol=1e-12)

    def test_full_resolution_cortex_is_homologous(self):
        ca, cc = anatomy.full_resolution(self.adult), anatomy.full_resolution(self.child)
        np.testing.assert_allclose(cc.rr, self.s * ca.rr)
        np.testing.assert_allclose(cc.area, self.s**2 * ca.area)
        np.testing.assert_array_equal(cc.valid, ca.valid)
        np.testing.assert_allclose(cc.dist_inner_skull, self.s * ca.dist_inner_skull, rtol=1e-9, atol=1e-12)


@unittest.skipUnless((paths.EXTERNAL / anatomy.INFANT_SUBJECTS / "ANTS2-0Years3T").exists(), "infant template not downloaded")
class TestTemplate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.t = anatomy.load_template()

    def test_head_frame_from_fiducials(self):
        f = self.t.fiducials
        self.assertGreater(f["nasion"][1], 0.05)
        np.testing.assert_allclose(f["nasion"][[0, 2]], 0.0, atol=1e-9)
        for k in ("lpa", "rpa"):
            np.testing.assert_allclose(f[k][1:], 0.0, atol=1e-9)
        self.assertLess(f["lpa"][0], 0.0)
        self.assertGreater(f["rpa"][0], 0.0)

    def test_native_dimensions(self):
        size = P.head_size(self.t)
        self.assertTrue(450 < size["ofc_mm"] < 520, size)  # 2-year head circumference, ~48-50 cm
        self.assertTrue(110 < size["inter_auricular_mm"] < 130, size)


class TestPlacements(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.info = neuromag.load_info("T3")
        cls.base = cls.info["dev_head_t"]["trans"]
        cls.adult = anatomy.load_sample()
        cls.child = anatomy.scaled(cls.adult, 0.85)
        cls.pl = P.placements(cls.info, cls.child, cls.base)

    def test_centred_is_the_adult_measured_position(self):
        np.testing.assert_allclose(self.pl["centred"]["trans"], self.base)

    def test_top_contact(self):
        top = self.pl["top"]
        self.assertGreater(top["moved_mm"], 10.0)  # a smaller head is raised
        self.assertTrue(20.0 < top["min_dist_mm"] <= 20.6, top["min_dist_mm"])  # last 0.5-mm step before 20 mm
        fit = P.HelmetFit(self.info, self.child)
        d = fit.distances(P.moved(top["trans"], P.translate([0, 0, 0.0005])))
        self.assertLessEqual(d.min(), P.CLEARANCE)

    def test_all_child_placements_feasible(self):
        for name, v in self.pl.items():
            self.assertTrue(v["feasible"], name)
            self.assertGreaterEqual(v["min_dist_mm"], 18.0, name)

    def test_contact_never_moves_down(self):
        fit = P.HelmetFit(self.info, self.adult)
        close = P.moved(self.base, P.translate([0.0, 0.0, 0.008]))  # adult raised 8 mm: within 20 mm of a coil
        self.assertLessEqual(fit.distances(close).min(), P.CLEARANCE)
        t, d = fit.contact(close, (0, 0, 1))
        self.assertEqual(d, 0.0)
        np.testing.assert_allclose(t, close)

    def test_counterfactual_helmet_scales_every_gap(self):
        # helmet scaled by s about the head origin around a head scaled by s: the adult's fit, scaled
        cf = P.counterfactual_helmet(self.info, self.child, self.base, 0.85)
        self.assertAlmostEqual(cf["k"], 0.85)
        dc = P.HelmetFit(cf["info"], self.child).distances(self.base)
        da = P.HelmetFit(self.info, self.adult).distances(self.base)
        np.testing.assert_allclose(dc, 0.85 * da, rtol=1e-9)
        same = P.scaled_helmet(self.info, 1.0, [0.01, 0.02, 0.03])
        np.testing.assert_allclose([c["loc"] for c in same["chs"]], [c["loc"] for c in self.info["chs"]])

    def test_helmet_regions_partition_magnetometers(self):
        regions = P.helmet_regions(self.info)
        self.assertEqual(len(regions), 8)
        idx = np.concatenate(list(regions.values()))
        self.assertEqual(len(idx), len(np.unique(idx)))  # disjoint


class TestSummaries(unittest.TestCase):
    def test_weighted_median(self):
        self.assertEqual(P.weighted_median(np.array([3.0, 1.0, 2.0])), 2.0)
        self.assertEqual(P.weighted_median(np.array([1.0, 2.0, 3.0]), np.array([1.0, 1.0, 10.0])), 3.0)
        self.assertEqual(P.weighted_median(np.array([1.0, np.nan, 3.0, 2.0]), np.ones(4)), 2.0)

    def test_grouped_bootstrap_of_a_constant(self):
        rng = np.random.default_rng(0)
        ci = P.grouped_bootstrap(np.full(30, 1.5), np.repeat(["a", "b", "c"], 10), None, rng, 50)
        self.assertEqual(ci, [1.5, 1.5])

    def test_target_areas_partition_the_usable_cortex(self):
        adult = anatomy.load_sample()
        cortex = anatomy.full_resolution(adult)
        targets = np.flatnonzero(cortex.usable)[::50]
        a = P.target_areas(cortex, targets)
        self.assertAlmostEqual(a.sum(), cortex.area[cortex.usable].sum(), places=12)


if __name__ == "__main__":
    unittest.main()
