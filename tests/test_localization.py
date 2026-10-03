import unittest

import numpy as np

from opmsquid import localization


class TestLocalization(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.n_ch, self.n_src = 40, 120
        # smooth, spatially structured gain: sources on a line, sensors see Gaussian profiles
        x_s = np.linspace(0, 1, self.n_src)
        x_c = np.linspace(0, 1, self.n_ch)
        self.gain = np.exp(-((x_c[:, None] - x_s[None, :]) ** 2) / (2 * 0.05**2)) * (1 + 0.2 * rng.random(self.n_src))
        a = rng.normal(size=(self.n_ch, self.n_ch)) * 0.1
        self.cov = np.eye(self.n_ch) + a @ a.T
        self.chol = np.linalg.cholesky(self.cov)
        self.src_rr = np.c_[x_s, np.zeros(self.n_src), np.zeros(self.n_src)]

    def test_dspm_is_unit_variance_on_noise(self):
        inv = localization.MNEInverse.make(self.gain, self.cov)
        noise = self.chol @ np.random.default_rng(1).standard_normal((self.n_ch, 40000))
        d = inv.apply(noise, "dSPM")
        np.testing.assert_allclose(d.var(axis=1), 1.0, rtol=0.05)

    def test_noiseless_source_is_localized(self):
        inv = localization.MNEInverse.make(self.gain, self.cov, snr=10.0)
        for j in (20, 60, 100):
            est = inv.apply(self.gain[:, [j]], "dSPM")[:, 0]
            err, i = localization.peak_error(est, self.src_rr, self.src_rr[j])
            self.assertLessEqual(abs(i - j), 2)
            self.assertLess(err, 0.02)

    def test_support_recovery(self):
        est = np.exp(-0.5 * ((np.arange(self.n_src) - 55) / 3.0) ** 2)  # peaked at source 55
        members = np.zeros(self.n_src, bool)
        members[53:58] = True  # the simulated patch: 5 sources around 55
        self.assertEqual(localization.support_recovery(est, members), 1.0)
        far = np.zeros(self.n_src, bool)
        far[3:8] = True
        self.assertEqual(localization.support_recovery(est, far), 0.0)
        # a neighbour outside the patch gets no credit, even if it is close to the centre
        lopsided = np.zeros(self.n_src, bool)
        lopsided[[55, 56, 57, 58, 59]] = True
        self.assertEqual(localization.support_recovery(est, lopsided), 0.6)  # the 5 strongest are 53-57

    def test_channel_views(self):
        from types import SimpleNamespace

        arrays = {"squid": SimpleNamespace(n=5, kinds=np.array(["mag", "grad", "grad", "mag", "grad"])), "opm": SimpleNamespace(n=3, kinds=None)}
        v = localization.channel_views(arrays, ("mag", "grad"))
        self.assertEqual(list(v), ["squid", "opm", "squid_mag", "squid_grad"])
        np.testing.assert_array_equal(v["squid_mag"][1], [0, 3])
        np.testing.assert_array_equal(v["squid_grad"][1], [1, 2, 4])
        self.assertEqual(v["opm"][0], "opm")
        np.testing.assert_array_equal(v["opm"][1], [0, 1, 2])
        with self.assertRaisesRegex(ValueError, "no Neuromag channels"):
            localization.channel_views(arrays, ("squid_mag",))

    def test_perturb_trans_bounds(self):
        rng = np.random.default_rng(3)
        t = np.eye(4)
        t[:3, 3] = [0.01, -0.02, 0.03]
        p = localization.perturb_trans(t, rng, shift=0.002, angle_deg=2.0)
        m = np.linalg.inv(t) @ p
        self.assertAlmostEqual(np.linalg.norm(m[:3, 3]), 0.002, places=12)
        self.assertAlmostEqual(np.degrees(np.arccos((np.trace(m[:3, :3]) - 1) / 2)), 2.0, places=9)


if __name__ == "__main__":
    unittest.main()
