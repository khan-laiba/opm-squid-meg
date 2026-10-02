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


class TestCensoredS50(unittest.TestCase):
    strengths = np.array([10.0, 20.0, 40.0, 80.0, 160.0, 320.0])

    def test_s50_bounds(self):
        s = self.strengths
        self.assertEqual(detection.s50_bounds([0.1, 0.2, 0.4, 0.6, 0.9, 1.0], s), (detection.s50_from([0.1, 0.2, 0.4, 0.6, 0.9, 1.0], s),) * 2)
        self.assertEqual(detection.s50_bounds([0.6, 0.8, 1, 1, 1, 1], s), (0.0, 10.0))  # already at 50 % at the weakest strength
        self.assertEqual(detection.s50_bounds([0, 0, 0.1, 0.2, 0.3, 0.4], s), (320.0, np.inf))  # never reached

    def test_ratio_bounds_and_labels(self):
        exact_o, exact_s = (40.0, 40.0), (60.0, 60.0)
        never, weakest = (320.0, np.inf), (0.0, 10.0)
        self.assertEqual(detection.ratio_bounds(exact_o, exact_s), (1.5, 1.5))
        self.assertEqual(detection.ratio_bounds(exact_o, never), (8.0, np.inf))  # Neuromag beyond the range: ratio > 320 / 40
        self.assertEqual(detection.ratio_bounds(never, exact_s), (0.0, 60.0 / 320.0))  # the OPM beyond: ratio < 60 / 320, not 0
        self.assertEqual(detection.ratio_bounds(never, never), (0.0, np.inf))
        self.assertEqual(detection.ratio_bounds(weakest, exact_s), (6.0, np.inf))  # the OPM at 50 % already at 10 nAm
        self.assertIsNone(detection.censoring_label(exact_o, exact_s))
        self.assertEqual(detection.censoring_label(never, never), "neither reaches 50 %")
        self.assertEqual(detection.censoring_label(exact_o, never), "Neuromag does not reach 50 %")
        self.assertEqual(detection.censoring_label(weakest, exact_s), "the OPM reaches 50 % at the weakest strength")

    def test_censored_interval_keeps_every_resample(self):
        r = np.linspace(1.0, 2.0, 1000)
        np.testing.assert_allclose(detection.censored_interval(r, r),
                                   np.percentile(r, [2.5, 97.5], method="inverted_cdf"))  # no censoring: plain percentiles
        # resamples where neither system reaches 50 % are unrestricted and widen both ends (the review's example)
        lo = np.r_[np.zeros(900), np.full(100, 2.0)]
        hi = np.r_[np.full(900, np.inf), np.full(100, 2.0)]
        self.assertEqual(detection.censored_interval(lo, hi), [None, None])
        # a fully censored sample of OPM-not-reached resamples is bounded above, not at 0
        self.assertEqual(detection.censored_interval(np.zeros(1000), np.full(1000, 0.2)), [None, 0.2])
        # a few censored resamples (< 2.5 %) leave the ends finite
        lo = np.r_[np.zeros(20), r[20:]]
        hi = np.r_[np.full(20, np.inf), r[20:]]
        a, b = detection.censored_interval(lo, hi)
        self.assertIsNotNone(a)
        self.assertIsNotNone(b)

    def test_format(self):
        self.assertEqual(detection.format_s50_ratio(dict(value=1.5, ci95=[1.2, None])), "1.50 [1.20-open]")
        self.assertEqual(detection.format_s50_ratio(dict(value=None, value_bounds=[1.13, None], value_censored="Neuromag does not reach 50 %",
                                                         ci95=[1.05, None])), "> 1.13 (Neuromag does not reach 50 %) [1.05-open]")
        self.assertEqual(detection.format_s50_ratio(dict(value=None, value_bounds=[None, None], value_censored="neither reaches 50 %",
                                                         ci95=None)), "neither reaches 50 %")


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
