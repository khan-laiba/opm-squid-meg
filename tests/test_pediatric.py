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
    NAMES = ("ANTS2-0Years3T", "ANTS18-0Months3T", "ANTS12-0Months3T")  # configs/g3b_pediatric.toml

    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.t = {n: anatomy.load_template(n) for n in cls.NAMES}

    def test_head_frame_from_fiducials(self):
        for n, t in self.t.items():
            with self.subTest(template=n):
                f = t.fiducials
                self.assertGreater(f["nasion"][1], 0.05)
                np.testing.assert_allclose(f["nasion"][[0, 2]], 0.0, atol=1e-9)
                for k in ("lpa", "rpa"):
                    np.testing.assert_allclose(f[k][1:], 0.0, atol=1e-9)
                self.assertLess(f["lpa"][0], 0.0)
                self.assertGreater(f["rpa"][0], 0.0)

    def test_native_dimensions(self):
        size = {n: P.head_size(t) for n, t in self.t.items()}
        for n, s in size.items():
            with self.subTest(template=n):
                self.assertTrue(450 < s["ofc_mm"] < 520, s)  # 12-24-month head circumference, ~45-50 cm
                self.assertTrue(110 < s["inter_auricular_mm"] < 130, s)
        ofc = [size[n]["ofc_mm"] for n in self.NAMES]
        self.assertEqual(ofc, sorted(ofc, reverse=True))  # the older template has the larger head

    def test_head_surface_already_on_the_scalp(self):
        for n, t in self.t.items():  # A-BEM-CONFORM is the identity for the templates
            with self.subTest(template=n):
                self.assertTrue(t.head_conform["identity"])
                self.assertEqual(t.head_conform["moved_max_mm"], 0.0)

    def test_bem_surfaces_nested(self):
        from scipy.spatial import cKDTree
        for n, t in self.t.items():
            with self.subTest(template=n):
                rr = {s["id"]: s["rr"] for s in t.bem_surfaces}
                ids = sorted(rr)  # brain (inner skull), skull, head
                for a, b in zip(ids, ids[1:]):
                    self.assertGreater(cKDTree(rr[b]).query(rr[a])[0].min(), 1e-4)  # the 12-month skull is 0.25 mm at its thinnest


class TestPlacements(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        cls.info = neuromag.load_info("T3")
        cls.base = cls.info["dev_head_t"]["trans"]
        cls.adult = anatomy.load_sample()
        cls.child = anatomy.scaled(cls.adult, 0.85)
        cls.pl = P.placements(cls.info, cls.child, cls.base)

    def test_rotation_variants(self):
        for name in ("pitch+10deg", "pitch-10deg", "roll+5deg", "roll-5deg", "yaw+10deg", "yaw-10deg"):
            self.assertIn(name, self.pl)
        # a yaw turns the head about device z: the rotation part of the device-to-head transform changes by 10 deg (to the
        # single precision of the measured transform)
        r = self.pl["yaw+10deg"]["trans"][:3, :3] @ self.base[:3, :3].T
        self.assertAlmostEqual(np.degrees(np.arccos((np.trace(r) - 1) / 2)), 10.0, delta=0.01)

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

    def test_lateral_centring_balances_the_helmet_halves(self):
        pose, dx = P.lateral_centring(self.info, self.child, self.base)
        fit = P.HelmetFit(self.info, self.child)
        regions = P.helmet_regions(self.info)
        left = np.concatenate([v for k, v in regions.items() if k.startswith("left")])
        right = np.concatenate([v for k, v in regions.items() if k.startswith("right")])
        d = fit.distances(pose)
        self.assertLess(abs(np.median(d[left]) - np.median(d[right])), 0.001)
        self.assertLessEqual(abs(dx), 0.020)
        self.assertIn("x-centred", self.pl)
        self.assertTrue(18.0 < self.pl["top-18mm"]["min_dist_mm"] <= 18.6)

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

    def test_categorical_maps_are_not_averaged(self):
        from opmsquid import plotting

        v = np.array([[3.0, 0.0, 0.0], [1.0, 2.0, 3.0], [2.0, 2.0, 1.0], [0.0, 3.0, 3.0]])
        np.testing.assert_array_equal(plotting.triangle_mode(v), [0.0, 1.0, 2.0, 3.0])  # never the mean (1.0 for 3, 0, 0)

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


class TestG3BSummaries(unittest.TestCase):
    """The summary helpers of scripts/g3b_pediatric_helmet.py on synthetic values."""

    @classmethod
    def setUpClass(cls):
        import sys

        sys.path.insert(0, str(paths.ROOT / "scripts"))
        import g3b_pediatric_helmet as G3

        cls.G3 = G3
        cls.rng = np.random.default_rng(0)

    def test_delta_strata_of_shifted_values(self):
        n = 60
        vc, va = np.linspace(10, 40, n), np.linspace(10, 40, n)
        groups = np.repeat(np.array(["a", "b", "c", "d"]), n // 4)
        rows = self.G3.delta_strata(np.full(n, 2.0), np.ones(n), groups, vc, np.full(n, 0.5), np.ones(n), groups, va,
                                    np.array([10.0, 25.0, 40.1]), self.rng, 50, 10)
        for row in rows:
            self.assertAlmostEqual(row["delta"], 1.5)
            self.assertEqual(row["ci95"], [1.5, 1.5])

    def test_sparse_strata_are_reported_not_summarised(self):
        x = np.ones(12)
        rows = self.G3.strata_rows(x, np.ones(12), np.repeat(["a", "b"], 6), np.r_[np.full(11, 5.0), 50.0],
                                   np.array([0.0, 10.0, 90.0]), self.rng, 10, 10)
        self.assertIn("median", rows[0])
        self.assertNotIn("median", rows[1])  # one target: sparse
        self.assertEqual(rows[1]["n"], 1)

    def test_parcel_table_skips_the_medial_wall_and_sparse_parcels(self):
        gc = np.array(["lh.a"] * 12 + ["lh.unknown"] * 12 + ["lh.b"] * 3)
        xc = np.r_[np.full(12, 3.0), np.full(12, 9.0), np.full(3, 1.0)]
        rows = {r["parcel"]: r for r in self.G3.parcel_table(xc, np.ones(27), gc, xc - 1.0, np.ones(27), gc, 10)}
        self.assertNotIn("lh.unknown", rows)
        self.assertAlmostEqual(rows["lh.a"]["delta"], 1.0)
        self.assertNotIn("delta", rows["lh.b"])


if __name__ == "__main__":
    unittest.main()
