import unittest

import mne
import numpy as np

from opmsquid import neuromag, paths


class TestNeuromag(unittest.TestCase):
    def test_ssp_projector_matches_mne(self):
        """Same projector as MNE's (private) make_projector for the sample recording's SSP."""
        from mne._fiff.proj import make_projector

        raw_info = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
        names = neuromag.load_info("T3").ch_names
        p, n = neuromag.ssp_projector(raw_info["projs"], names)
        p_mne, n_mne, _ = make_projector(raw_info["projs"], names)
        self.assertEqual(n, n_mne)
        np.testing.assert_allclose(p, p_mne, atol=1e-12)
        np.testing.assert_allclose(p @ p, p, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
