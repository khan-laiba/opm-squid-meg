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


if __name__ == "__main__":
    unittest.main()
