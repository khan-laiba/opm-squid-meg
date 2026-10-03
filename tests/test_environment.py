"""Physical response tests of the environmental basis (GOAL G5: uniform-field gradiometer
cancellation, known-gradient response)."""
import unittest

import numpy as np
import mne

from opmsquid import environment, neuromag, opm, paths

mne.set_log_level("WARNING")
HAVE_SAMPLE = (paths.SAMPLE_MEG / neuromag.RAW_FILE).exists()


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestEnvironmentBasis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.info = neuromag.load_info("T3")
        cls.basis = environment.external_basis(cls.info)
        cls.kind = np.array(["grad" if int(ch["coil_type"]) == 3014 else "mag" for ch in cls.info["chs"]])

    def test_planar_gradiometers_reject_uniform_fields(self):
        grads, mags = self.kind == "grad", self.kind == "mag"
        uniform = self.basis[:, :3]
        # a 1 T uniform field gives O(1) T on magnetometers and ~0 T/m on gradiometers
        self.assertGreater(np.median(np.linalg.norm(uniform[mags], axis=1)), 0.5)
        self.assertLess(np.max(np.abs(uniform[grads])), 1e-9)

    def test_known_gradient_response(self):
        """A planar gradiometer measures d(B . n)/du along its baseline u: for the field
        B = G (r - r0) with G symmetric, the reading is k n^T G u, where k = sum_i w_i x_i is the
        effective baseline response of MNE's integration rule. For coil 3014 the tabulated
        weights (+/-14.9858) give k = 0.99991 rather than 1 (rounding in coil_def.dat)."""
        from mne.forward._make_forward import _read_coil_defs
        coil = next(c for c in _read_coil_defs() if c["coil_type"] == 3014 and c["accuracy"] == 2)
        k_eff = float(coil["w"] @ coil["rmag"][:, 0])
        self.assertLess(abs(k_eff - 1.0), 1e-4)
        chs = self.info["chs"]
        dev_head = self.info["dev_head_t"]["trans"]
        for k in np.flatnonzero(self.kind == "grad")[:20]:
            rot = dev_head[:3, :3] @ chs[k]["loc"][3:12].reshape(3, 3).T  # columns ex, ey, ez (head frame)
            ex, ez = rot[:, 0], rot[:, 2]
            for gi, g in enumerate(environment.GRADIENT_BASIS):
                expected = k_eff * (ez @ g @ ex)  # baseline of Vectorview planar gradiometers lies along ex
                self.assertAlmostEqual(self.basis[k, 3 + gi], expected, places=12)

    def test_magnetometer_reads_field_along_normal(self):
        chs = self.info["chs"]
        R = self.info["dev_head_t"]["trans"][:3, :3]
        for k in np.flatnonzero(self.kind == "mag")[:10]:
            ez = R @ chs[k]["loc"][9:12]
            np.testing.assert_allclose(self.basis[k, :3], ez, atol=1e-12)

    def test_opm_cell_averages_linear_field_exactly(self):
        """For linear fields the cubature average equals the field at the cell centre."""
        arr = opm.OPMArray(np.array([[0.02, -0.01, 0.12]]), np.array([[0.0, 0.6, 0.8]]),
                           np.zeros((1, 3)), np.array([0]), 0.0, 0.0, "t")
        b = environment.external_basis(opm.make_info(arr), coil_def=opm.coil_def_file())
        rel = arr.pos[0] - np.array([0.0, 0.0, 0.04])
        expected = np.concatenate([arr.axis[0], [arr.axis[0] @ g @ rel for g in environment.GRADIENT_BASIS]])
        np.testing.assert_allclose(b[0], expected, atol=1e-12)

    def test_empty_room_fit_ignores_bad_channels(self):
        """Synthetic room field through the real coil responses plus white noise, one channel
        corrupted: with that channel marked bad the fit recovers the field (explained > 99 %)."""
        from opmsquid import noise

        rng = np.random.default_rng(3)
        fs, n = 600.0, 12000
        coef = rng.normal(size=(8, n)) * np.r_[np.full(3, 1e-12), np.full(5, 1e-11)][:, None]
        data = self.basis @ coef
        data += rng.normal(size=data.shape) * np.where(self.kind == "mag", 2e-15, 2e-13)[:, None]
        bad = self.info.ch_names[10]
        data[10] += rng.normal(size=n) * 1e-9
        raw = mne.io.RawArray(data, mne.create_info(self.info.ch_names, fs, ["mag" if k == "mag" else "grad" for k in self.kind]))
        filt = noise.AnalysisFilter(fs=fs, l_freq=1.0, h_freq=100.0, order=4)
        env = environment.fit_empty_room(raw, self.info, filt, trim_s=1.0, bads=[bad])
        self.assertGreater(env.explained_fraction["mag"], 0.99)
        self.assertGreater(env.explained_fraction["grad"], 0.9)
        # the per-channel residual (measured-spectrum SQUID scenario): about the white noise left, NaN for the bad channel
        rv = env.residual_var
        self.assertTrue(np.isnan(rv[10]) and np.isfinite(np.delete(rv, 10)).all())
        white = filt.enbw() / (fs / 2) * np.where(self.kind == "mag", 2e-15, 2e-13) ** 2
        ratio = np.delete(rv / white, 10)
        self.assertTrue(0.8 < np.median(ratio) < 1.05, np.median(ratio))


if __name__ == "__main__":
    unittest.main()
