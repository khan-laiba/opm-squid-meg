import unittest

import mne
import numpy as np

from opmsquid import anatomy, g2, neuromag, paths

mne.set_log_level("WARNING")
HAVE_SAMPLE = (paths.SAMPLE_MEG / neuromag.RAW_FILE).exists()


class TestG2Helpers(unittest.TestCase):
    def test_sample_covariance_converges(self):
        rng = np.random.default_rng(0)
        a = rng.normal(size=(5, 5))
        cov = a @ a.T + np.eye(5)
        for shrink in (False, True):
            est = g2.sample_covariance(cov, 200_000, rng, shrink=shrink)
            self.assertLess(np.linalg.norm(est - cov) / np.linalg.norm(cov), 0.02)


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestHeadPositions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.info = neuromag.load_info("T3")
        cls.variants = g2.head_position_variants(cls.info, anatomy.load_sample())

    def test_translation_moves_coils_the_other_way_in_the_head_frame(self):
        base = self.variants["measured"]["trans"]
        t = self.variants["z+5mm"]["trans"]
        coil = np.array([0.01, 0.02, 0.10, 1.0])  # a device-frame point
        shift = (t @ coil - base @ coil)[:3]
        np.testing.assert_allclose(shift, -base[:3, :3] @ np.array([0, 0, 0.005]), atol=1e-12)

    def test_pitch_is_a_rotation_about_the_head_origin(self):
        base = self.variants["measured"]["trans"]
        t = self.variants["pitch+5deg"]["trans"]
        origin_dev = np.linalg.inv(base) @ np.array([0, 0, 0, 1.0])
        np.testing.assert_allclose(t @ origin_dev, base @ origin_dev, atol=1e-12)  # head origin fixed
        motion = np.linalg.inv(np.linalg.inv(base) @ t)  # the rigid head motion in the device frame
        np.testing.assert_allclose(np.degrees(np.arccos((np.trace(motion[:3, :3]) - 1) / 2)), 5.0, atol=1e-9)
        np.testing.assert_allclose(motion[:3, 0], [1, 0, 0], atol=1e-12)  # about the device x axis

    def test_feasibility_and_well_fitted_clearance(self):
        v = self.variants
        self.assertAlmostEqual(v["well_fitted"]["min_dist"], 0.020, delta=0.0006)
        self.assertLess(v["well_fitted"]["median_dist"], v["measured"]["median_dist"])
        self.assertTrue(all(x["feasible"] for x in v.values()))
        self.assertEqual(len(v), 10)


if __name__ == "__main__":
    unittest.main()
