"""Integration tests on the MNE sample subject (skipped if the sample data are missing)."""
import unittest

import numpy as np
import mne

from opmsquid import anatomy, forward, neuromag, opm, paths

mne.set_log_level("WARNING")
HAVE_SAMPLE = (paths.SAMPLE_MEG / neuromag.RAW_FILE).exists()


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestForward(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subject = anatomy.load_sample("oct6")
        cls.info = neuromag.load_info("T3")
        cls.bem = cls.subject.bem_model()  # 1-layer, 0.3 S/m
        src = cls.subject.src
        rng = np.random.default_rng(0)
        cls.pick = [np.sort(rng.choice(s["nuse"], 10, replace=False)) for s in src]
        # a reduced surface source space with the picked vertices only
        cls.src_small = src.copy()
        for s, p in zip(cls.src_small, cls.pick):
            keep = s["vertno"][p]
            s["inuse"] = np.zeros(s["np"], int)
            s["inuse"][keep] = 1
            s["vertno"] = keep
            s["nuse"] = len(keep)
            for k in ("use_tris", "dist", "dist_limit", "patch_inds", "pinfo", "nearest", "nearest_dist"):
                s[k] = None
            s["nuse_tri"] = 0

    def test_discrete_equals_surface_forward(self):
        g_surf, _ = forward.fixed_gain(self.info, self.subject.trans, self.src_small, self.bem, use_cps=False,
                                       use_cache=False)
        rr = np.concatenate([s["rr"][s["vertno"]] for s in self.src_small])
        nn = np.concatenate([s["nn"][s["vertno"]] for s in self.src_small])
        g_disc, _ = forward.discrete_gain(self.info, self.subject.trans, rr, nn, self.bem, use_cache=False)
        rel = np.abs(g_disc - g_surf).max(axis=0) / np.abs(g_surf).max(axis=0)
        self.assertLess(rel.max(), 1e-5)

    def test_cache_roundtrip(self):
        rr = self.subject.src_rr[:5]
        nn = self.subject.src_nn[:5]
        g1, m1 = forward.discrete_gain(self.info, self.subject.trans, rr, nn, self.bem)
        g2, m2 = forward.discrete_gain(self.info, self.subject.trans, rr, nn, self.bem)
        np.testing.assert_array_equal(g1, g2)
        self.assertEqual(m1["key"], m2["key"])
        self.assertEqual(len(m1["ch_names"]), 306)

    def test_matched_opm_array_is_physically_feasible(self):
        dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE)
        arr, rep = opm.matched_to_neuromag(self.info, self.subject.trans, self.subject.scalp, dig)
        self.assertEqual(rep["n_ray_hits"], 102)
        self.assertEqual(arr.n_sites, rep["n_kept"])
        self.assertGreaterEqual(rep["n_kept"], 95)
        from scipy.spatial import cKDTree
        clearance = cKDTree(self.subject.scalp.rr).query(mne.transforms.apply_trans(self.subject.trans, arr.pos))[0]
        self.assertGreaterEqual(clearance.min(), opm.STANDOFF - 0.001 - 1e-9)
        self.assertGreater(rep["min_spacing_mm"].min(), 17.0)  # no two 10-mm cells closer than ~17 mm
        np.testing.assert_allclose(np.linalg.norm(arr.axis, axis=1), 1.0, atol=1e-12)

    def test_opm_coil_in_realistic_forward(self):
        """An OPM on the scalp above a source sees the same field with the 1-point and the
        27-point cell model to within the expected finite-cell effect (< 3 %)."""
        rr, nn = self.subject.src_rr[:1], self.subject.src_nn[:1]
        head_pts = mne.transforms.apply_trans(np.linalg.inv(self.subject.trans["trans"]), rr)
        d = head_pts[0] / np.linalg.norm(head_pts[0])
        center = head_pts[0] + 0.03 * d
        arr = opm.OPMArray(center[None], d[None], center[None], np.array([0]), 0.0, 0.0, "test")
        info = opm.make_info(arr)
        g_cell, _ = forward.discrete_gain(info, self.subject.trans, rr, nn, self.bem, coil_def=opm.coil_def_file(),
                                          use_cache=False)
        g_point, _ = forward.discrete_gain(info, self.subject.trans, rr, nn, self.bem,
                                           coil_def=opm.coil_def_file(cell_size=1e-6), use_cache=False)
        self.assertLess(abs(g_cell[0, 0] / g_point[0, 0] - 1), 0.03)


if __name__ == "__main__":
    unittest.main()
