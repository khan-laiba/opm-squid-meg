import unittest

import numpy as np
import scipy.sparse as sp

from opmsquid import hunold


def _grid_mesh(n=30, spacing=0.001):
    """Planar n x n grid mesh (edge lengths `spacing`), vertex areas spacing^2."""
    idx = np.arange(n * n).reshape(n, n)
    rows, cols = [], []
    for a, b in ((idx[:, :-1], idx[:, 1:]), (idx[:-1, :], idx[1:, :])):
        rows += [a.ravel(), b.ravel()]
        cols += [b.ravel(), a.ravel()]
    r, c = np.concatenate(rows), np.concatenate(cols)
    adj = sp.csr_matrix((np.full(len(r), spacing), (r, c)), shape=(n * n, n * n))
    return adj, np.full(n * n, spacing**2), idx


class TestHunold(unittest.TestCase):
    def test_spike_waveform_landmarks(self):
        w = hunold.spike_waveform(1000.0)
        self.assertEqual(len(w), 200)
        self.assertAlmostEqual(w.max(), 600e-9)
        self.assertEqual(int(np.argmax(w)), 24)  # main positive peak at 24 ms
        self.assertEqual(int(np.argmin(w[:100])), 54)  # main negative trough at 54 ms
        self.assertAlmostEqual(w.min() / w.max(), -253 / 600, delta=0.01)
        self.assertEqual(w[199], 0.0)

    def test_grow_patch_area_and_orientation_window(self):
        adj, area, idx = _grid_mesh()
        orient = np.full(area.size, 45.0)
        orient[idx[:, 20:].ravel()] = 70.0  # outside the +/-10 deg window of a 45-deg seed
        valid = np.ones(area.size, bool)
        seed = int(idx[15, 10])
        target = 20.5e-6  # not a multiple of the 1-mm^2 vertex area (avoids exact ties)
        members = hunold.grow_patch(seed, adj, orient, area, valid, target_area=target)
        self.assertGreater(area[members].sum(), target)
        self.assertLessEqual(area[members].sum() - area[members[-1]], target)  # stops when first exceeding
        self.assertTrue(np.all(np.abs(orient[members] - 45.0) <= 10.0))
        # a seed confined to a too-small admissible region is rejected
        orient2 = np.full(area.size, 45.0)
        orient2[seed] = 5.0
        orient2[idx[15, 11]] = 3.0
        self.assertIsNone(hunold.grow_patch(seed, adj, orient2, area, valid, target_area=20e-6))

    def test_background_timecourses(self):
        x = hunold.background_timecourses(20, 6000, 1000.0, np.random.default_rng(0))
        np.testing.assert_allclose(np.abs(x).max(axis=1), 10e-9)
        f = np.fft.rfftfreq(6000, 1e-3)
        p = np.mean(np.abs(np.fft.rfft(x, axis=1)) ** 2, axis=0)
        self.assertLess(p[f > 60].sum() / p.sum(), 1e-3)  # band-limited to the EEG bands (<= 45 Hz)

    def test_background_amplitude_rayleigh(self):
        rng = np.random.default_rng(1)
        x = rng.standard_normal((5, 60000))
        a = hunold.background_amplitude(x, slice(10000, 50000))
        np.testing.assert_allclose(a, 2 * np.sqrt(np.pi / 2), rtol=0.02)  # 2 <e> = 2.51 sigma

    def test_spike_snr_selects_max_noise_free_channel(self):
        w = hunold.spike_waveform(1000.0)
        topo = np.array([[1.0, -3.0], [-2.0, 1.0], [0.5, 2.0]])
        bg = np.zeros((3, 1000))
        out = hunold.spike_snr(topo, w, bg, 500, a_bg=np.array([1.0, 2.0, 4.0]))
        np.testing.assert_array_equal(out["channel"], [1, 0])
        np.testing.assert_allclose(out["peak"], [2.0 * 600e-9 / 2.0, 3.0 * 600e-9 / 1.0])
        np.testing.assert_allclose(out["noisy_peak"], out["peak"])  # zero background

    def test_bins(self):
        r, c = hunold.bin_index(np.array([19.9, 20.0, 24.99, 59.9, 60.0]), np.array([0.0, 9.99, 10.0, 90.0, 45.0]))
        np.testing.assert_array_equal(r, [-1, 0, 0, 7, -1])
        np.testing.assert_array_equal(c, [-1, 0, 1, 8, -1])
        self.assertEqual(int(hunold.PAPER_DIPOLE_COUNTS.sum()), 3783)


if __name__ == "__main__":
    unittest.main()
