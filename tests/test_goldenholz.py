import unittest

import numpy as np
import scipy.sparse as sp

from opmsquid import goldenholz


class TestGoldenholz(unittest.TestCase):
    def test_eq1_with_one_over_n(self):
        topo = np.array([[3.0], [4.0]])  # a b_k for two channels
        var = np.array([1.0, 4.0])
        self.assertAlmostEqual(float(goldenholz.eq1_snr_db(topo, var)[0]), 10 * np.log10((9 + 4) / 2))

    def test_calibration_rule(self):
        kinds = np.array(["grad"] * 4 + ["mag"] * 2)
        rec = np.array([2.0, 2.0, 2.0, 8.0, 3.0, 5.0])
        aat = np.ones(6)
        s2, per = goldenholz.calibrate_source_variance(rec, aat, kinds)
        self.assertEqual(per, {"grad": 2.0, "mag": 4.0})
        self.assertAlmostEqual(s2, (4 * 2.0 + 2 * 4.0) / 6)

    def test_poisson_disk_spacing(self):
        rng = np.random.default_rng(0)
        pts = rng.uniform(0, 0.05, size=(4000, 3))
        keep = goldenholz.poisson_disk(pts, 0.007, rng)
        d = np.linalg.norm(pts[keep][:, None] - pts[keep][None], axis=-1) + np.eye(len(keep))
        self.assertGreaterEqual(d.min(), 0.007)
        self.assertGreater(len(keep), 50)

    def test_patch_topographies_match_explicit_sum(self):
        rng = np.random.default_rng(0)
        g = rng.normal(size=(5, 40)).astype(np.float32)
        members = [np.array([1, 5, 7]), np.array([0]), np.arange(10, 30)]
        w = rng.uniform(1, 2, 40)
        ref = np.stack([g[:, m].astype(np.float64) @ w[m] for m in members], axis=1)
        np.testing.assert_allclose(goldenholz.patch_topographies(g, members, np.arange(40), w), ref, rtol=1e-12)
        # global vertex indices mapped to gain columns (as with the valid-vertex full-resolution matrices)
        col_of = np.full(100, -1)
        col_of[np.arange(0, 80, 2)] = np.arange(40)  # vertex 2k -> column k
        members_global = [2 * m for m in members]
        np.testing.assert_allclose(goldenholz.patch_topographies(g, members_global, col_of, np.repeat(w, 2)[:100]),
                                   np.stack([g[:, m].astype(np.float64) @ np.repeat(w, 2)[:100][2 * m] for m in members], axis=1),
                                   rtol=1e-12)

    def test_eq1_chunking_and_scale(self):
        rng = np.random.default_rng(1)
        t = rng.normal(size=(6, 50))
        v = rng.uniform(1, 2, 6)
        np.testing.assert_allclose(goldenholz.eq1_snr_db(t, v, scale=3.0, chunk=7), goldenholz.eq1_snr_db(3.0 * t, v))

    def test_geodesic_patches_on_grid(self):
        n, h = 41, 0.001
        idx = np.arange(n * n).reshape(n, n)
        r, c = [], []
        for a, b in ((idx[:, :-1], idx[:, 1:]), (idx[:-1, :], idx[1:, :])):
            r += [a.ravel(), b.ravel()]
            c += [b.ravel(), a.ravel()]
        r, c = np.concatenate(r), np.concatenate(c)
        adj = sp.csr_matrix((np.full(len(r), h), (r, c)), shape=(n * n, n * n))
        centre = int(idx[20, 20])
        (m,) = goldenholz.geodesic_patches(adj, np.array([centre]), 0.005, np.ones(n * n, bool))
        # Manhattan (edge-path) ball of radius 5 edges: 2*5*6 + 1 = 61 vertices
        self.assertEqual(len(m), 61)
        # invalid vertices are never members (they may still carry the path)
        valid = np.ones(n * n, bool)
        valid[idx[20, 21:23]] = False
        (m2,) = goldenholz.geodesic_patches(adj, np.array([centre]), 0.005, valid)
        self.assertEqual(len(m2), 59)
        self.assertFalse(np.any(np.isin(idx[20, 21:23], m2)))


if __name__ == "__main__":
    unittest.main()
