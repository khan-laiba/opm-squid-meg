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

    def test_background_is_stationary(self):
        """The per-source maximum (which sets the normalisation) is not held by filter edge
        transients: maxima fall within 300 ms of an edge about as often as uniform (10 %)."""
        x = hunold.background_timecourses(400, 6000, 1000.0, np.random.default_rng(3))
        at = np.abs(x).argmax(axis=1)
        self.assertLess(np.mean((at < 300) | (at >= 5700)), 0.2)
        np.testing.assert_allclose(x[:, :500].std(), x[:, 2750:3250].std(), rtol=0.1)

    def test_rendered_centroid_sd(self):
        """The Fig. 6 rendering emulation scales linearly with the drawing scale and reads a
        band-limited trace a few per cent low (7.65-ms pixel columns average fast components)."""
        x = hunold.background_timecourses(1, 6000, 1000.0, np.random.default_rng(5))[0]
        x = x / x.std()
        a = hunold.rendered_centroid_sd(x, 1000.0, px_per_unit=3.0)
        b = hunold.rendered_centroid_sd(x, 1000.0, px_per_unit=6.0)
        self.assertAlmostEqual(b / a, 2.0, delta=0.06)
        true = 3.0 * np.sqrt(np.mean([np.var(x[i:i + 900]) for i in range(0, 5101, 900)]))
        self.assertTrue(0.85 < a / true < 1.0, a / true)

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
        p2p = (w.max() - w.min()) / w.max()  # 853/600 for the digitised complex
        np.testing.assert_allclose(out["p2p"], out["peak"] * p2p)
        np.testing.assert_allclose(out["noisy_p2p"], out["p2p"])
        self.assertAlmostEqual(p2p, 853 / 600, delta=0.01)

    def test_bin_means_and_stratified_sampling(self):
        rng = np.random.default_rng(4)
        depth = rng.uniform(20, 60, 5000)
        orient = rng.uniform(0, 90, 5000)
        counts = np.full((8, 9), 3)
        counts[0, 0] = 10_000  # more than available: capped
        idx, achieved = hunold.sample_by_bins(depth, orient, np.arange(5000), counts, rng)
        r, c = hunold.bin_index(depth[idx], orient[idx])
        self.assertEqual(achieved[0, 0], int(np.sum((r == 0) & (c == 0))))
        self.assertTrue(np.all(achieved[1:, :] == 3) and len(np.unique(idx)) == len(idx))
        vals = np.arange(len(idx), dtype=float)
        mean, count = hunold.bin_means(vals, r, c, (8, 9))
        np.testing.assert_array_equal(count, achieved)
        np.testing.assert_allclose(mean[3, 4], vals[(r == 3) & (c == 4)].mean())

    def test_bem_node_descriptors(self):
        """Radial source under a flat 'scalp' plane: depth = distance to the nearest node,
        orientation = angle to the nearest inner-skull node normal, folded to 0-90 deg."""
        g = np.stack(np.meshgrid(np.arange(-5, 6) * 0.005, np.arange(-5, 6) * 0.005), -1).reshape(-1, 2)
        scalp = np.c_[g, np.full(len(g), 0.10)]
        skull = np.c_[g, np.full(len(g), 0.09)]
        skull_nn = np.tile([0.0, 0.0, 1.0], (len(g), 1))
        pts = np.array([[0.0, 0.0, 0.07], [0.0, 0.0, 0.07], [0.0, 0.0, 0.07]])
        nrm = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, np.sin(np.radians(30)), -np.cos(np.radians(30))]])
        d, o = hunold.bem_node_descriptors(pts, nrm, scalp, skull, skull_nn)
        np.testing.assert_allclose(d, 30.0)
        np.testing.assert_allclose(o, [0.0, 90.0, 30.0], atol=1e-9)

    def test_bins(self):
        r, c = hunold.bin_index(np.array([19.9, 20.0, 24.99, 59.9, 60.0]), np.array([0.0, 9.99, 10.0, 90.0, 45.0]))
        np.testing.assert_array_equal(r, [-1, 0, 0, 7, -1])
        np.testing.assert_array_equal(c, [-1, 0, 1, 8, -1])
        self.assertEqual(int(hunold.PAPER_DIPOLE_COUNTS.sum()), 3783)


if __name__ == "__main__":
    unittest.main()
