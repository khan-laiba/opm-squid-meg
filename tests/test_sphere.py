import unittest

import numpy as np

from opmsquid import sphere

H, B, XI = 0.095, 0.080, 0.018


class TestSphere(unittest.TestCase):
    def test_eq1_matches_independent_numerical_maximum(self):
        for r_q in (0.01, 0.03, 0.05, 0.064, 0.08):
            for s in (H, H + XI):
                closed = sphere.bmax_radial(r_q, s)[0]
                numerical, theta, phi = sphere.numerical_peak_radial(r_q, s)
                self.assertLess(abs(numerical / closed - 1.0), 1e-8, (r_q, s))
                # the peak lies in the plane perpendicular to the dipole (phi = 0 or pi)
                self.assertLess(min(abs(np.sin(phi)), 1.0), 1e-5)

    def test_sarvas_radial_equals_biot_savart(self):
        rng = np.random.default_rng(0)
        u = rng.normal(size=(50, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        r = 0.1 * u
        r_q, q = np.array([0.01, -0.02, 0.05]), np.array([10e-9, 20e-9, -5e-9])
        b = sphere.sarvas_field(r, r_q, q)
        a = r - r_q
        bs = sphere.MU0_OVER_4PI * np.cross(q, a) / np.linalg.norm(a, axis=1, keepdims=True) ** 3
        np.testing.assert_allclose(np.sum(b * u, axis=1), np.sum(bs * u, axis=1), rtol=1e-10, atol=1e-25)

    def test_radial_dipole_is_silent(self):
        r = 0.1 * np.array([[0.3, 0.4, np.sqrt(1 - 0.25)], [0.0, 0.6, 0.8]])
        r_q = np.array([0.0, 0.0, 0.05])
        b = sphere.sarvas_field(r, r_q, 30e-9 * r_q / np.linalg.norm(r_q))
        self.assertLess(np.max(np.abs(b)), 1e-25)

    def test_sphere_centre_is_excluded(self):
        self.assertEqual(sphere.bmax_radial(0.0, H)[0], 0.0)
        self.assertEqual(sphere.numerical_peak_radial(0.0, H)[0], 0.0)
        self.assertTrue(np.isnan(sphere.signal_ratio(H)[0]))

    def test_equal_snr_depth_eta3(self):
        d = sphere.equal_snr_depth(3.0, H, B, 0.0, XI)
        self.assertAlmostEqual(d * 1e3, 27.665025, places=5)  # paper: "approximately 28 mm"
        self.assertEqual(round(d * 1e3), 28)

    def test_equal_snr_depth_existence_range(self):
        lo = sphere.centre_limit_ratio(H, 0.0, XI)  # (113/95)^3
        hi = sphere.signal_ratio(H - B)[0]  # at the brain surface
        self.assertAlmostEqual(lo, (113 / 95) ** 3, places=12)
        self.assertAlmostEqual(hi, 5.308565, places=5)
        self.assertTrue(np.isnan(sphere.equal_snr_depth(0.99 * lo)))
        self.assertTrue(np.isnan(sphere.equal_snr_depth(1.01 * hi)))
        self.assertFalse(np.isnan(sphere.equal_snr_depth(1.01 * lo)))
        self.assertAlmostEqual(sphere.equal_snr_depth(hi * (1 - 1e-12)), H - B, places=6)


if __name__ == "__main__":
    unittest.main()
