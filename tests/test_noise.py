import unittest

import numpy as np

from opmsquid import noise


class TestNoise(unittest.TestCase):
    def test_unfiltered_white_noise_variance(self):
        rng = np.random.default_rng(0)
        x = noise.white_noise(10e-15, 1000.0, (4, 200_000), rng)
        np.testing.assert_allclose(np.var(x, axis=1), (10e-15) ** 2 * 500.0, rtol=0.02)

    def test_filtered_variance_matches_simulation(self):
        """Integrating PSD x |H_composite|^2 reproduces the variance of simulated filtered noise,
        for zero-phase (|H|^4 power) and causal (|H|^2) application."""
        rng = np.random.default_rng(1)
        fs, asd = 1000.0, 7e-15
        for zero_phase in (True, False):
            filt = noise.AnalysisFilter(fs=fs, l_freq=1.0, h_freq=40.0, order=4, zero_phase=zero_phase)
            x = filt.apply(noise.white_noise(asd, fs, (8, 400_000), rng))[:, 20_000:-20_000]
            simulated = np.mean(np.var(x, axis=1))
            self.assertLess(abs(simulated / noise.white_variance(asd, filt) - 1.0), 0.02, zero_phase)

    def test_coloured_psd(self):
        filt = noise.AnalysisFilter(fs=1000.0, l_freq=None, h_freq=100.0, order=8, zero_phase=False)
        # a flat PSD integrated through the filter equals asd^2 * ENBW
        self.assertAlmostEqual(filt.variance(4.0), 4.0 * filt.enbw(), places=6)
        # ENBW of a sharp 100-Hz low-pass is close to 100 Hz
        self.assertLess(abs(filt.enbw() / 100.0 - 1.0), 0.05)


if __name__ == "__main__":
    unittest.main()
