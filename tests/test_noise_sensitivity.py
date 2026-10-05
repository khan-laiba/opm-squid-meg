"""Adult noise-model sensitivity analyses (scripts/study_noise_sensitivity.py): the sub-band
partition of the analysis filter, the frequency-resolved detectability, the far-field calibration,
the fixed-weight projector, the far-field coil integration and the break-even interpolation. Light:
no BEM solution and no OPM array construction."""
import importlib.util
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import mne

from opmsquid import environment, g2, neuromag, noise, noisemodel, paths

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # registered before execution: its dataclasses look their module up there
    spec.loader.exec_module(mod)
    return mod


S = _load("study_noise_sensitivity")
FS = 600.614990234375  # the sample recording's sampling rate
EDGES = [0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0]


def _subbands(filt):
    edges = EDGES + [filt.fs / 2.0 + 1e-9]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        f = S.SubbandFilter(filt, lo, hi)
        out.append(dict(lo=lo, hi=hi, filt=f, enbw=f.enbw()))
    return out


class TestSubbands(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filt = noise.AnalysisFilter(fs=FS, l_freq=1.0, h_freq=40.0, order=4)
        cls.subs = _subbands(cls.filt)

    def test_subbands_partition_the_band(self):
        self.assertAlmostEqual(sum(s["enbw"] for s in self.subs) / self.filt.enbw(), 1.0, places=9)

    def test_fft_subbands_add_up_to_the_band_filter(self):
        x = noise.white_noise(1.0, FS, (3, int(60 * FS)), np.random.default_rng(0))
        total = sum(s["filt"].apply(x) for s in self.subs)
        whole = S.SubbandFilter(self.filt, 0.0, FS).apply(x)
        np.testing.assert_allclose(total, whole, atol=1e-9 * np.abs(whole).max())
        trim = int(5 * FS)  # the FFT filter equals the Butterworth one away from the record's ends
        ref = self.filt.apply(x)[:, trim:-trim]
        self.assertLess(np.abs(whole[:, trim:-trim] - ref).max() / np.abs(ref).max(), 1e-3)

    def test_one_over_f_factor(self):
        self.assertEqual(S.opm_shape(self.filt, 0.0, 1.0), 1.0)
        f1, f10 = S.opm_shape(self.filt, 1.0, 1.0), S.opm_shape(self.filt, 10.0, 1.0)
        self.assertTrue(1.0 < f1 < f10)
        self.assertAlmostEqual((f10 - 1.0) / (f1 - 1.0), 10.0, places=6)  # the 1/f excess is linear in the corner

    def test_signal_weights_sum_to_one(self):
        for kind in ("flat", "spike-wave"):
            w = S.signal_weights(self.subs, self.filt, kind)
            self.assertAlmostEqual(w.sum(), 1.0, places=12)
            self.assertTrue(np.all(w >= 0))


class TestFrequencyResolvedDetectability(unittest.TestCase):
    """With every noise term sharing one spectrum the frequency-resolved detectability equals the
    band detectability; a coloured OPM term puts it between the white sensor and the band-variance
    (spatial whitening only) value for every target."""

    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(3)
        ns, no, nt = 24, 16, 30
        kinds = np.array(["mag"] * 8 + ["grad"] * 16)
        cls.arrays = {"squid": g2.Array("squid", None, kinds, None, {}), "opm_x": g2.Array("opm_x", None, np.array(["mag"] * no), None, {})}
        cls.topo = {"squid": rng.standard_normal((ns, nt)), "opm_x": rng.standard_normal((no, nt))}

        def psd(n):
            a = rng.standard_normal((n, n + 5))
            return a @ a.T / (n + 5)

        cls.B = {"squid": psd(ns), "opm_x": psd(no)}
        cls.E = {"squid": 0.3 * psd(ns), "opm_x": 0.3 * psd(no)}
        cls.basis = {"squid": rng.standard_normal((ns, 8)), "opm_x": rng.standard_normal((no, 8))}
        cls.ivar_sq = np.where(kinds == "mag", 0.5, 2.0)
        cls.w_opm = 0.7
        cls.filt = noise.AnalysisFilter(fs=FS, l_freq=1.0, h_freq=40.0, order=4)
        cls.subs = _subbands(cls.filt)
        enbw = sum(s["enbw"] for s in cls.subs)
        cls.rows = []
        for s in cls.subs:
            r = s["enbw"] / enbw
            cls.rows.append(dict(noises={
                "squid": noisemodel.ArrayNoise(cls.ivar_sq * r, cls.B["squid"] * r, cls.E["squid"] * r, cls.basis["squid"]),
                "opm_x": noisemodel.ArrayNoise(np.full(no, np.nan), cls.B["opm_x"] * r, cls.E["opm_x"] * r, cls.basis["opm_x"])}))
        cls.ctx = SimpleNamespace(arrays=cls.arrays, topo=lambda n: cls.topo[n])
        cls.white = np.sqrt(cls.w_opm / enbw) * 1e15  # fr_detect takes fT/sqrt(Hz)
        cls.conds = ["intrinsic+brain", "projected"]

    def _band(self, name, ivar):
        return S.detect(self.topo[name], noisemodel.ArrayNoise(ivar, self.B[name], self.E[name], self.basis[name]), self.arrays[name], self.conds)

    def test_common_spectrum_reproduces_band_detectability(self):
        for kind in ("flat", "spike-wave"):
            u = S.signal_weights(self.subs, self.filt, kind)
            d = S.fr_detect(self.ctx, self.subs, self.rows, u, self.white, 0.0, 1.0, ["squid", "opm_x"], self.conds)
            for name, ivar in (("squid", self.ivar_sq), ("opm_x", np.full(16, self.w_opm))):
                ref = self._band(name, ivar)
                for k in ref:
                    np.testing.assert_allclose(d[name][k], ref[k], rtol=1e-9)

    def test_coloured_opm_noise_between_white_and_band_variance(self):
        u = S.signal_weights(self.subs, self.filt, "flat")
        for corner in (1.0, 3.0, 10.0):
            d = S.fr_detect(self.ctx, self.subs, self.rows, u, self.white, corner, 1.0, ["opm_x"], self.conds)["opm_x"]
            white = self._band("opm_x", np.full(16, self.w_opm))
            band = self._band("opm_x", np.full(16, self.w_opm * S.opm_shape(self.filt, corner, 1.0)))
            for k in white:
                self.assertTrue(np.all(d[k] <= white[k] * (1 + 1e-12)))
                self.assertTrue(np.all(d[k] >= band[k] * (1 - 1e-12)))


class TestFarFieldCalibration(unittest.TestCase):
    """At the study's real scales (brain scale ~3.5e-14): the gradiometer brain noise is matched and
    the far field takes its level on the magnetometers (a fixed RMS, or the magnetometer shortfall)."""

    def setUp(self):
        rng = np.random.default_rng(1)
        n = 306
        kinds = np.array((["mag"] + ["grad"] * 2) * 102)
        self.ctx = SimpleNamespace(grads=kinds == "grad", mags=kinds == "mag")
        s0 = 3.5356e-14
        self.du = np.where(self.ctx.grads, (3.71e-12) ** 2 / s0, (1.92e-13) ** 2 / s0) * rng.uniform(0.5, 1.5, n)
        f = rng.standard_normal((n, 3)) * np.where(self.ctx.grads, 0.05, 1.0)[:, None]
        self.ff = f @ f.T
        self.ff /= np.median(np.diag(self.ff)[self.ctx.mags])
        self.target = np.median(self.du[self.ctx.grads]) * s0
        self.s0 = s0

    def test_fixed_level(self):
        s, a2, info = S.calibrate_with_far_field(self.ctx, np.diag(self.du), self.ff, self.target, (2.62e-13) ** 2, 1.13e-13, "recalibrated")
        self.assertTrue(info["feasible"])
        self.assertGreater(s / self.s0, 0.99)  # the far field barely touches the gradiometers
        self.assertAlmostEqual(np.sqrt(a2 * np.median(np.diag(self.ff)[self.ctx.mags])) / 1.13e-13, 1.0, places=9)
        self.assertAlmostEqual(np.median(s * self.du[self.ctx.grads] + a2 * np.diag(self.ff)[self.ctx.grads]) / self.target, 1.0, places=6)

    def test_shortfall_level(self):
        meas = (2.62e-13) ** 2
        s, a2, info = S.calibrate_with_far_field(self.ctx, np.diag(self.du), self.ff, self.target, meas, "magnetometer_shortfall", "recalibrated")
        self.assertTrue(info["feasible"])
        self.assertGreater(s / self.s0, 0.99)
        self.assertAlmostEqual(a2 * np.median(np.diag(self.ff)[self.ctx.mags]) / (meas - np.median(s * self.du[self.ctx.mags])), 1.0, places=6)


    def test_added_mode_keeps_the_published_background(self):
        meas = (2.62e-13) ** 2
        s, a2, info = S.calibrate_with_far_field(self.ctx, np.diag(self.du), self.ff, self.target, meas, "magnetometer_shortfall", "added")
        self.assertEqual(s, self.target / np.median(self.du[self.ctx.grads]))
        self.assertAlmostEqual(a2 * np.median(np.diag(self.ff)[self.ctx.mags]) / (meas - np.median(s * self.du[self.ctx.mags])), 1.0, places=9)
        self.assertGreaterEqual(info["grad_rms_model_over_measured"], 1.0)  # the far field on top over-predicts the gradiometers

    def test_infeasible_when_the_far_field_exceeds_the_gradiometers(self):
        ff = self.ff + np.diag(np.where(self.ctx.grads, 1e3, 0.0))  # a far field with strong gradiometer content
        _, _, info = S.calibrate_with_far_field(self.ctx, np.diag(self.du), ff, self.target, (2.62e-13) ** 2, "magnetometer_shortfall", "recalibrated")
        self.assertFalse(info["feasible"])


class TestFixedProjector(unittest.TestCase):
    def test_band_weights_kept(self):
        rng = np.random.default_rng(0)
        e, iv = rng.standard_normal((20, 8)), rng.uniform(1, 2, 20)
        a = noisemodel.ArrayNoise(iv, np.eye(20), np.eye(20), e)
        b = S.FixedProjectorNoise(iv * rng.uniform(1, 3, 20), np.eye(20), np.eye(20), e, iv)
        np.testing.assert_allclose(a.projector(), b.projector(), atol=1e-12)


@unittest.skipUnless((paths.SAMPLE_MEG / neuromag.RAW_FILE).exists(), "MNE sample data not available")
class TestFarFieldFields(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mne.set_log_level("WARNING")
        info = neuromag.load_info("T3")
        cls.arr = g2.Array("squid", info, neuromag.channel_kinds(info), None, {})

    def test_coil_integration_matches_room_field_basis(self):
        basis = environment.external_basis(self.arr.info, (0.0, 0.0, 0.04), None)
        b0 = np.array([1.0, -2.0, 0.5])
        resp = np.array([c["w"] @ (c["cosmag"] @ b0) for c in environment.coil_integration(self.arr.info, None)])
        np.testing.assert_allclose(resp, basis[:, :3] @ b0, rtol=1e-12, atol=1e-12 * np.abs(resp).max())

    def test_dipole_far_field_falls_with_distance(self):
        mags = self.arr.kinds == "mag"
        g = [np.median(np.abs(S.dipole_field_gain(self.arr, np.array([0.0, 0.05, -d]), "current_dipole")[mags])) for d in (0.2, 0.3)]
        self.assertGreater(g[0], g[1])
        # a 1-A m current dipole's primary field: |B| <= mu0/(4 pi) / r^2 at every magnetometer integration point
        r0 = np.array([0.0, 0.0, -0.2])
        coils = environment.coil_integration(self.arr.info, None)
        r_min = min(np.linalg.norm(c["rmag"] - r0, axis=1).min() for c, m in zip(coils, mags) if m)
        self.assertLess(np.abs(S.dipole_field_gain(self.arr, r0, "current_dipole")[mags]).max(), 1e-7 / r_min ** 2)


class TestBreakEven(unittest.TestCase):
    def test_log_linear_crossing(self):
        grid = np.array([5.0, 10.0, 20.0, 40.0])
        med = np.log2(np.array([1.4, 1.2, 0.9, 0.7]))
        c = S.crossing(grid, med)
        self.assertIsNone(c["bound"])
        x = np.interp(0.0, [med[2], med[1]], [np.log(20.0), np.log(10.0)])
        self.assertAlmostEqual(c["value"], float(np.exp(x)))
        self.assertEqual(S.crossing(grid, np.log2(np.array([1.4, 1.3, 1.2, 1.1])))["bound"], ">40")
        self.assertEqual(S.crossing(grid, np.log2(np.array([0.9, 0.8, 0.7, 0.6])))["bound"], "<5")


if __name__ == "__main__":
    unittest.main()
