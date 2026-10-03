import unittest

import numpy as np

from opmsquid import metrics, noisemodel


def _model(rng, n=12, units=None):
    units = np.ones(n) if units is None else units
    e = rng.normal(size=(n, 8)) * units[:, None]
    a = rng.normal(size=(n, 20)) * units[:, None]
    return noisemodel.ArrayNoise(intrinsic_var=rng.uniform(1, 2, n) * units**2, brain_cov=a @ a.T,
                                 env_cov=e @ np.diag(rng.uniform(1, 3, 8)) @ e.T, ext_basis=e)


class TestPairedNoise(unittest.TestCase):
    """G2 plug-in covariances from one common noise realization (review, 2026-10-02)."""

    def setUp(self):
        from opmsquid import background, environment, g2

        self.g2, self.background = g2, background
        rng = np.random.default_rng(1)
        self.areas = rng.uniform(0.5, 1.5, 30)
        self.gains = {"a": rng.normal(size=(6, 30)), "b": rng.normal(size=(4, 30))}
        coef_cov = np.diag(rng.uniform(1, 2, 8))
        self.env = environment.EnvironmentModel(coef_cov, None, 1.0, {}, np.zeros(3))
        self.noises = {}
        for k, g in self.gains.items():
            e = rng.normal(size=(g.shape[0], 8))
            self.noises[k] = noisemodel.ArrayNoise(rng.uniform(1, 2, g.shape[0]), 2.0 * background.sensor_covariance(g, self.areas),
                                                   e @ coef_cov @ e.T, e)

    def test_one_realization_seen_by_every_array(self):
        n = 200000
        s = self.g2.paired_noise_samples(self.noises, self.gains, self.areas, 2.0, self.env, n, np.random.default_rng(2))
        cross = s["a"]["brain"] @ s["b"]["brain"].T / n  # cross-array background covariance: exact pairing
        expected = 2.0 * (self.gains["a"] * self.areas) @ self.gains["b"].T
        np.testing.assert_allclose(cross, expected, atol=0.03 * np.abs(expected).max())
        for k, nz in self.noises.items():
            env_c = s[k]["env"] @ s[k]["env"].T / n
            np.testing.assert_allclose(env_c, nz.env_cov, atol=0.03 * np.abs(nz.env_cov).max())
            np.testing.assert_allclose(np.var(s[k]["intrinsic"], axis=1), nz.intrinsic_var, rtol=0.02)
        self.assertLess(np.abs(s["a"]["intrinsic"] @ s["b"]["intrinsic"].T / n).max(), 0.02)  # independent sensors

    def test_measured_squid_variance(self):
        import dataclasses

        kinds = np.array(["mag", "grad", "grad", "mag", "grad"])
        env = dataclasses.replace(self.env, residual_var=np.array([1.0, 4.0, np.nan, 3.0, 6.0]))
        np.testing.assert_allclose(self.g2.measured_squid_variance(env, kinds), [1.0, 4.0, 5.0, 3.0, 6.0])  # bad: type median
        with self.assertRaisesRegex(ValueError, "no good mag"):
            self.g2.measured_squid_variance(dataclasses.replace(env, residual_var=np.array([np.nan, 4.0, 5.0, np.nan, 6.0])), kinds)
        with self.assertRaisesRegex(ValueError, "refit"):
            self.g2.measured_squid_variance(self.env, kinds)  # an environment model without the residual

    def test_plugin_covariances_estimate_each_condition(self):
        n = 100000
        s = self.g2.paired_noise_samples(self.noises, self.gains, self.areas, 2.0, self.env, n, np.random.default_rng(3))
        for k, nz in self.noises.items():
            est = self.g2.plugin_covariances(s[k], nz, noisemodel.CONDITIONS)
            for cond in noisemodel.CONDITIONS:
                c = nz.covariance(cond)
                with self.subTest(array=k, condition=cond):
                    np.testing.assert_allclose(est[cond], c, atol=0.03 * np.abs(c).max())


class TestNoiseModel(unittest.TestCase):
    def test_projector_removes_external_subspace_and_is_idempotent(self):
        m = _model(np.random.default_rng(0))
        p = m.projector()
        np.testing.assert_allclose(p @ m.ext_basis, 0.0, atol=1e-10)
        np.testing.assert_allclose(p @ p, p, atol=1e-10)
        c = m.covariance("projected")
        self.assertEqual(np.linalg.matrix_rank(c, tol=1e-8 * np.abs(c).max()), 12 - 8)

    def test_conditions_are_cumulative(self):
        m = _model(np.random.default_rng(1))
        c1, c2, c3 = (m.covariance(k) for k in noisemodel.CONDITIONS[:3])
        np.testing.assert_allclose(c2 - c1, m.brain_cov)
        np.testing.assert_allclose(c3 - c2, m.env_cov)

    def test_projected_detectability_is_unit_invariant(self):
        """Mixed T and T/m channels: a consistent change of units leaves the projected metric unchanged."""
        rng = np.random.default_rng(2)
        units = np.r_[np.ones(4), np.full(8, 1e-2)]
        m1 = _model(np.random.default_rng(3))
        s = rng.normal(size=(12, 3))
        m2 = noisemodel.ArrayNoise(m1.intrinsic_var * units**2, m1.brain_cov * np.outer(units, units),
                                   m1.env_cov * np.outer(units, units), m1.ext_basis * units[:, None])
        d1 = metrics.detectability(m1.signal(s, "projected"), m1.covariance("projected"))
        d2 = metrics.detectability(m2.signal(s * units[:, None], "projected"), m2.covariance("projected"))
        np.testing.assert_allclose(d1, d2, rtol=1e-6)

    def test_grid_areas_conserve_total(self):
        rng = np.random.default_rng(4)
        v = rng.uniform(0, 0.1, (500, 3))
        a = rng.uniform(0.5, 1.5, 500) * 1e-6
        g = noisemodel.grid_areas(v[:40], v, a)
        self.assertAlmostEqual(g.sum(), a.sum(), places=15)


if __name__ == "__main__":
    unittest.main()
