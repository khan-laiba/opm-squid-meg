import unittest

import numpy as np
from scipy import signal

from opmsquid import ied, noise


class TestIED(unittest.TestCase):
    def test_waveform_and_stretch(self):
        w = ied.ied_waveform(1000.0)
        self.assertAlmostEqual(w.max(), 1.0)
        self.assertEqual(int(np.argmax(w)), 24)  # digitised peak at 24 ms
        w2 = ied.ied_waveform(1000.0, stretch=1.5)
        self.assertEqual(int(np.argmax(w2)), 36)
        self.assertAlmostEqual(len(w2) / len(w), 1.5, delta=0.02)

    def test_csd_synthesis_reproduces_coherence(self):
        """Two channels: x2 = x1 filtered + independent noise; the synthesised pair has the same
        normalised covariance and magnitude coherence at 10 Hz."""
        rng = np.random.default_rng(0)
        fs, n = 200.0, 400_000
        a = rng.standard_normal(n)
        b, a_ = signal.butter(2, [5, 20], btype="bandpass", fs=fs)
        x1 = signal.lfilter(b, a_, a)
        x2 = 0.8 * x1 + 0.3 * signal.lfilter(b, a_, rng.standard_normal(n))
        x = np.vstack([x1, x2])
        f, s = ied.csd_matrix(x, fs, nperseg=1024)
        y = ied.synthesize_from_csd(f, s, 200_000, fs, np.random.default_rng(1))
        cx, cy = np.corrcoef(x)[0, 1], np.corrcoef(y)[0, 1]
        self.assertAlmostEqual(cx, cy, delta=0.02)
        fx, coh_x = signal.coherence(x[0], x[1], fs=fs, nperseg=1024)
        _, coh_y = signal.coherence(y[0], y[1], fs=fs, nperseg=1024)
        k = np.argmin(np.abs(fx - 10))
        self.assertAlmostEqual(coh_x[k], coh_y[k], delta=0.03)

    def test_noise_generator_matches_band_covariance(self):
        """Per-channel variance of the generated noise = diag(G S G^T + B C B^T) + ASD^2 ENBW."""
        rng = np.random.default_rng(2)
        fs = 600.0
        filt = noise.AnalysisFilter(fs=fs, l_freq=1.0, h_freq=40.0, order=4)
        n_ch, n_grid = 6, 50
        g = rng.normal(size=(n_ch, n_grid))
        basis = rng.normal(size=(n_ch, 8))
        grid_var = rng.uniform(0.5, 1.5, n_grid)
        env_coef = filt.apply(rng.standard_normal((8, 60_000)))
        env_cov = np.cov(env_coef)
        asd = np.full(n_ch, 0.2)
        gen = ied.NoiseGenerator({"a": ied.ArraySpec(g, basis, asd)}, grid_var, env_coef, env_cov, filt, decimate=4)
        y = gen.segment(120.0, rng)["a"]
        expected = np.diag(g @ np.diag(grid_var) @ g.T + basis @ env_cov @ basis.T) + asd**2 * filt.enbw()
        np.testing.assert_allclose(y.var(axis=1), expected, rtol=0.05)
        self.assertEqual(y.shape[1], int(120 * fs) // 4)

    def test_injection_places_the_peak(self):
        filt = noise.AnalysisFilter(fs=600.0, l_freq=1.0, h_freq=40.0, order=4)
        tpl, k = ied.filtered_template(filt, decimate=4)
        data = np.zeros((2, 500))
        ied.inject(data, np.array([1.0, -2.0]), tpl, k, 3.0, 250)
        self.assertLessEqual(abs(int(np.argmax(data[0])) - 250), 1)
        np.testing.assert_allclose(data[1], -2.0 * data[0])
        ied.inject(data, np.array([1.0, 1.0]), tpl, k, 1.0, 2)  # partly outside: no error
        self.assertTrue(np.isfinite(data).all())


if __name__ == "__main__":
    unittest.main()
