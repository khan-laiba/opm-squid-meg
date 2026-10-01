import unittest

import numpy as np

from opmsquid import detection


def _template():
    t = np.arange(40)
    tpl = np.exp(-0.5 * ((t - 12) / 3.0) ** 2) - 0.4 * np.exp(-0.5 * ((t - 24) / 5.0) ** 2)
    return tpl, int(np.argmax(tpl))


class TestDetection(unittest.TestCase):
    def test_matched_filter_alignment(self):
        tpl, pk = _template()
        y = np.zeros(300)
        y[100 - pk:100 - pk + len(tpl)] = tpl
        out = detection.matched_filter(y, tpl, pk)
        self.assertEqual(int(np.argmax(out)), 100)
        self.assertAlmostEqual(out[100], float(np.sum(tpl**2)))
        np.testing.assert_allclose(detection.matched_filter(np.vstack([y, 2 * y]), tpl, pk)[1], 2 * out, atol=1e-12)

    def test_oracle_null_is_standard_normal_and_detects(self):
        rng = np.random.default_rng(0)
        n_ch = 12
        mix = rng.normal(size=(n_ch, n_ch))
        null = mix @ rng.standard_normal((n_ch, 60000))
        w = detection.whitener_from_null(null)
        tpl, pk = _template()
        s = rng.normal(size=n_ch)
        oracle, samples = detection.Oracle.calibrate(null, w, {1.0: (tpl, pk)}, {"src": s}, rng=rng)
        z = samples[("src", 1.0)]
        self.assertAlmostEqual(z.std(), 1.0, places=6)
        self.assertLess(abs(np.mean(z)), 0.05)
        data = mix @ rng.standard_normal((n_ch, 2000))
        amp = 6.0 / (np.linalg.norm(w.apply(s)) * np.linalg.norm(tpl))  # expected z ~ 6
        data[:, 1000 - pk:1000 - pk + len(tpl)] += amp * np.outer(s, tpl)
        self.assertGreater(oracle.statistic(data, "src", 1.0, s, 1000), 3.5)

    def test_scan_detector_rates_and_detection(self):
        rng = np.random.default_rng(1)
        n_ch, n_cand = 10, 30
        mix = rng.normal(size=(n_ch, n_ch)) * 0.3 + np.eye(n_ch)
        cal = mix @ rng.standard_normal((n_ch, 60000))
        held = mix @ rng.standard_normal((n_ch, 60000))
        w = detection.whitener_from_null(cal)
        tpl, pk = _template()
        cands = rng.normal(size=(n_ch, n_cand))
        det = detection.ScanDetector.build(cal, w, {1.0: (tpl, pk)}, cands, refractory=50)
        # outputs are normalised on the calibration data
        part = next(det._outputs(cal, 60000))
        np.testing.assert_allclose(np.sqrt(np.mean(part[0] ** 2, axis=1)), 1.0, rtol=1e-6)  # RMS-normalised
        stat_cal, _ = det.statistic(cal)
        p_cal, h_cal = det.events(stat_cal)
        minutes = 60000 / 100.0 / 60.0  # nominal 100 Hz
        thr = detection.threshold_for_rate(h_cal, minutes, 2.0)
        self.assertLessEqual(np.sum(h_cal > thr), 2.0 * minutes)
        stat_h, _ = det.statistic(held)
        _, h_h = det.events(stat_h)
        self.assertLess(np.sum(h_h > thr), 4 * 2.0 * minutes)  # held-out rate of the same order
        # a strong event on candidate 7 is found at the right time and candidate
        data = mix @ rng.standard_normal((n_ch, 3000))
        amp = 12.0 / (np.linalg.norm(w.apply(cands[:, 7])) * np.linalg.norm(tpl))
        data[:, 1500 - pk:1500 - pk + len(tpl)] += amp * np.outer(cands[:, 7], tpl)
        stat, arg = det.statistic(data, chunk=700)
        peaks, heights = det.events(stat)
        hit = detection.match_events(peaks, heights, np.array([1500]), tol=3, threshold=thr)
        self.assertTrue(hit[0])
        self.assertEqual(int(arg[1500]), 7)

    def test_match_events(self):
        peaks = np.array([10, 50, 90])
        heights = np.array([5.0, 1.0, 5.0])
        hit = detection.match_events(peaks, heights, np.array([12, 50, 95, 200]), tol=3, threshold=2.0)
        np.testing.assert_array_equal(hit, [True, False, False, False])

    def test_event_height_scores_with_the_emitted_events(self):
        peaks = np.array([10, 13, 50, 90])
        heights = np.array([5.0, 7.0, 1.0, 5.0])
        truth = [12, 50, 95, 200]
        h = [detection.event_height(peaks, heights, t, 3) for t in truth]
        self.assertEqual(h, [7.0, 1.0, 0.0, 0.0])  # the highest emitted event within +/- 3 samples, 0 if none
        # the same decision as match_events at any threshold
        for thr in (0.5, 2.0, 6.0):
            np.testing.assert_array_equal(np.array(h) > thr, detection.match_events(peaks, heights, np.array(truth), 3, thr))
        # an injected peak the refractory period removed (a higher event 20 samples away) is not scored
        stat = np.zeros(200)
        stat[100], stat[120] = 4.0, 6.0
        det = detection.ScanDetector(None, {}, np.zeros((1, 1)), np.ones((1, 1)), refractory=50)
        self.assertEqual(detection.event_height(*det.events(stat), 100, 3), 0.0)
        self.assertEqual(detection.event_height(*det.events(stat), 118, 3), 6.0)


class TestSignFlip(unittest.TestCase):
    def test_exact_values(self):
        # 5 locations all favouring one side: only the all-same-sign patterns reach |sum| = 5
        self.assertAlmostEqual(detection.sign_flip_p([1, 1, 1, 1, 1]), 2 / 32)
        self.assertAlmostEqual(detection.sign_flip_p([2, 0, 0, 1]), 2 / 4)  # zeros carry no sign
        self.assertEqual(detection.sign_flip_p([1, -1]), 1.0)
        self.assertEqual(detection.sign_flip_p([0, 0]), 1.0)

    def test_monte_carlo_agrees_with_exact(self):
        x = np.random.default_rng(0).normal(0.3, 1.0, 22)  # 22 > 20 locations: Monte Carlo branch
        signs = 1 - 2 * ((np.arange(2 ** 22)[:, None] >> np.arange(22)) & 1).astype(np.int8)
        exact = np.mean(np.abs(signs @ x) >= abs(x.sum()) - 1e-9)
        self.assertAlmostEqual(detection.sign_flip_p(x, n_mc=50000), exact, delta=0.01)


if __name__ == "__main__":
    unittest.main()
