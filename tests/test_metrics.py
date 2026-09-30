import unittest

import numpy as np

from opmsquid import metrics


def _random_problem(rng, n=12, rank=None):
    a = rng.normal(size=(n, n if rank is None else rank))
    cov = a @ a.T + (0 if rank is not None else 0.1 * np.eye(n))
    s = rng.normal(size=(n, 3))
    return s, cov


class TestMetrics(unittest.TestCase):
    def test_definitions_on_diagonal_noise(self):
        s = np.array([3.0, -4.0, 1.0])
        var = np.array([1.0, 4.0, 1.0])
        self.assertAlmostEqual(metrics.peak_channel_snr(s, var), 3.0)
        self.assertAlmostEqual(metrics.mean_power_snr_db(s, var), 10 * np.log10((9 + 4 + 1) / 3))
        self.assertAlmostEqual(metrics.detectability(s, np.diag(var)), np.sqrt(9 + 4 + 1))

    def test_unit_invariance(self):
        """Rescaling channels (e.g. T/m -> fT/cm for gradiometers) leaves every metric unchanged."""
        rng = np.random.default_rng(1)
        s, cov = _random_problem(rng)
        scale = 10.0 ** rng.uniform(-13, 13, size=len(s))
        s2, cov2 = s * scale[:, None], cov * np.outer(scale, scale)
        for f in (metrics.peak_channel_snr, metrics.mean_power_snr_db, metrics.detectability):
            np.testing.assert_allclose(f(s2, cov2), f(s, cov), rtol=1e-9)

    def test_detectability_matches_inverse_for_well_conditioned_noise(self):
        rng = np.random.default_rng(2)
        s, cov = _random_problem(rng)
        direct = np.sqrt(np.einsum("ij,ij->j", s, np.linalg.solve(cov, s)))
        np.testing.assert_allclose(metrics.detectability(s, cov), direct, rtol=1e-9)

    def test_rank_deficient_whitening(self):
        rng = np.random.default_rng(3)
        n, r = 20, 7
        s, cov = _random_problem(rng, n, rank=r)
        w = metrics.whitener(cov)
        self.assertEqual(w.rank, r)
        np.testing.assert_allclose(w.matrix @ cov @ w.matrix.T, np.eye(r), atol=1e-8)
        # a topography inside the noise subspace keeps its full detectability
        v = cov @ rng.normal(size=n)
        np.testing.assert_allclose(metrics.detectability(v, w=w), np.sqrt(v @ np.linalg.pinv(cov) @ v), rtol=1e-6)

    def test_ledoit_wolf(self):
        rng = np.random.default_rng(4)
        true = np.diag(np.linspace(1, 3, 10)) * 1e-26
        x = rng.multivariate_normal(np.zeros(10), true, size=200).T
        est, shrink = metrics.ledoit_wolf_covariance(x)
        self.assertTrue(0.0 <= shrink <= 1.0)
        # shrinkage improves on the empirical estimate for this well-conditioned diagonal truth
        emp = metrics.empirical_covariance(x)
        self.assertLess(np.linalg.norm(est - true), np.linalg.norm(emp - true))
        # unit invariance of the shrinkage coefficient
        self.assertAlmostEqual(metrics.ledoit_wolf_covariance(x * 1e13)[1], shrink, places=10)


if __name__ == "__main__":
    unittest.main()
