import unittest

import numpy as np

from opmsquid import metrics, noisemodel


def _model(rng, n=12, units=None):
    units = np.ones(n) if units is None else units
    e = rng.normal(size=(n, 8)) * units[:, None]
    a = rng.normal(size=(n, 20)) * units[:, None]
    return noisemodel.ArrayNoise(intrinsic_var=rng.uniform(1, 2, n) * units**2, brain_cov=a @ a.T,
                                 env_cov=e @ np.diag(rng.uniform(1, 3, 8)) @ e.T, ext_basis=e)


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
