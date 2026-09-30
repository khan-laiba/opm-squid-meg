import unittest

import numpy as np

from opmsquid import background


def _plane(n):
    """n x n grid of sources on a 40 x 40 mm plane at z = 0, with their areas."""
    x = (np.arange(n) + 0.5) / n * 0.04 - 0.02
    xx, yy = np.meshgrid(x, x, indexing="ij")
    rr = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(n * n)])
    return rr, np.full(n * n, (0.04 / n) ** 2)


def _gain(rr, sensors):
    """A smooth, dipole-like lead field: 1 / |r - p|^3 for sensors p above the plane."""
    d = np.linalg.norm(sensors[:, None, :] - rr[None, :, :], axis=-1)
    return 1e-9 / d**3


SENSORS = np.array([[0.0, 0.0, 0.03], [0.01, 0.0, 0.03], [0.0, -0.015, 0.035]])


class TestBackground(unittest.TestCase):
    def test_independent_model_is_mesh_invariant(self):
        covs = []
        for n in (20, 40, 80):
            rr, a = _plane(n)
            covs.append(background.sensor_covariance(_gain(rr, SENSORS), background.moment_covariance(a)))
        np.testing.assert_allclose(covs[1], covs[2], rtol=2e-3)
        np.testing.assert_allclose(covs[0], covs[2], rtol=1e-2)

    def test_correlated_model_tends_to_independent_and_is_mesh_invariant(self):
        rr, a = _plane(60)
        g = _gain(rr, SENSORS)
        ind = background.sensor_covariance(g, background.moment_covariance(a))
        small = background.sensor_covariance(g, background.moment_covariance(a, rr, 0.0006))
        np.testing.assert_allclose(small, ind, rtol=0.1)
        rr2, a2 = _plane(40)
        c1 = background.sensor_covariance(g, background.moment_covariance(a, rr, 0.005))
        c2 = background.sensor_covariance(_gain(rr2, SENSORS), background.moment_covariance(a2, rr2, 0.005))
        np.testing.assert_allclose(c1, c2, rtol=0.02)

    def test_joint_factor_reproduces_joint_covariance(self):
        rng = np.random.default_rng(0)
        g1, g2 = rng.normal(size=(5, 50)), rng.normal(size=(7, 50)) * 1e3  # different units
        mc = rng.uniform(1, 2, 50)
        f, sl = background.joint_factor([g1, g2], mc)
        joint = background.sensor_covariance(np.vstack([g1, g2]), mc)
        np.testing.assert_allclose(f @ f.T, joint, rtol=1e-8, atol=1e-10 * np.abs(joint).max())
        self.assertEqual(sl[1], slice(5, 12))

    def test_calibration_and_pink_noise(self):
        c = np.diag([1.0, 4.0, 9.0])
        s2 = background.calibrate(c, [0, 1, 2], 8.0)
        self.assertAlmostEqual(float(np.median(np.diag(s2 * c))), 8.0)
        x = background.pink_noise(4, 20_000, 1000.0, np.random.default_rng(1))
        np.testing.assert_allclose(x.std(axis=1), 1.0)
        f = np.fft.rfftfreq(20_000, 1e-3)
        p = np.mean(np.abs(np.fft.rfft(x, axis=1)) ** 2, axis=0)
        lo, hi = p[(f > 2) & (f < 4)].mean(), p[(f > 20) & (f < 40)].mean()
        self.assertAlmostEqual(np.log10(lo / hi), 1.0, delta=0.15)  # ~1/f: a decade in power


if __name__ == "__main__":
    unittest.main()
