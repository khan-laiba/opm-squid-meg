import unittest

import numpy as np
import mne

from opmsquid import opm, sphere

mne.set_log_level("WARNING")


def _cube_average_reference(center, axis, r_q, q_vec, size, n=10):
    """Independent reference: field component along `axis`, averaged over a cube of side
    `size` centred at `center` (cube axes aligned with the coil frame), 10x10x10 Gauss-Legendre."""
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = x * size / 2, w / w.sum()
    a = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    ex = np.cross(a, axis)
    ex /= np.linalg.norm(ex)
    ey = np.cross(axis, ex)
    g = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3)
    wg = np.einsum("i,j,k->ijk", w, w, w).ravel()
    pts = center + g[:, :1] * ex + g[:, 1:2] * ey + g[:, 2:3] * axis
    b = sphere.sarvas_field(pts, r_q, q_vec) @ axis
    return float(wg @ b)


class TestOPMCoil(unittest.TestCase):
    def setUp(self):
        self.fname = opm.coil_def_file()

    def test_coil_file_adds_opm_to_standard_coils(self):
        text = self.fname.read_text()
        self.assertIn(f"   {opm.OPM_COIL_TYPE}    ", text)
        self.assertNotIn("\n\n", text)  # MNE's parser rejects blank lines
        from mne.forward._make_forward import _read_coil_defs
        with mne.use_coil_def(self.fname):
            defs = {(c["coil_type"], c["accuracy"]): c for c in _read_coil_defs()}
        for coil in (3014, 3024, opm.OPM_COIL_TYPE):  # standard T3 coils stay available
            self.assertIn((coil, 2), defs)
        self.assertEqual(len(defs[(opm.OPM_COIL_TYPE, 2)]["w"]), 27)
        self.assertAlmostEqual(float(np.sum(defs[(opm.OPM_COIL_TYPE, 2)]["w"])), 1.0, places=10)
        # weights of each accuracy level sum to one
        for acc, n in ((0, 1), (1, 8), (2, 27)):
            w, g = (np.array([1.0]), np.zeros((1, 3))) if n == 1 else opm._gauss_cube(round(n ** (1 / 3)), 0.01)
            self.assertAlmostEqual(w.sum(), 1.0, places=12)
            self.assertEqual(len(w), n)

    def _forward(self, r_q, center, axis):
        """MNE forward (sphere model) of a 30 nAm +y dipole at (0, 0, r_q) for one OPM. MNE
        1.13.2's make_forward_solution always uses the 'accurate' coil level (27 points)."""
        arr = opm.OPMArray(center[None], axis[None], center[None], np.array([0]), 0.0, 0.0, "test")
        info = opm.make_info(arr)
        sph = mne.make_sphere_model(r0=(0.0, 0.0, 0.0), head_radius=0.09, relative_radii=(0.9, 1.0),
                                    sigmas=(0.33, 0.33))
        src = mne.setup_volume_source_space(pos=dict(rr=np.array([[0.0, 0.0, r_q]]), nn=np.array([[0.0, 1.0, 0.0]])))
        with mne.use_coil_def(self.fname):
            fwd = mne.make_forward_solution(info, None, src, sph)
        return float(fwd["sol"]["data"][0, 1]) * 30e-9  # column 1 = y-oriented unit dipole

    def test_mne_cube_integration_matches_independent_average(self):
        r_q, s = 0.07, 0.10
        theta = 0.25  # near the peak of the radial field
        center = s * np.array([np.sin(theta), 0.0, np.cos(theta)])
        axis = center / s
        q = np.array([0.0, 30e-9, 0.0])
        rq = np.array([0.0, 0.0, r_q])
        ref27 = _cube_average_reference(center, axis, rq, q, opm.CELL_SIZE, n=3)
        ref_dense = _cube_average_reference(center, axis, rq, q, opm.CELL_SIZE, n=10)
        mne_value = self._forward(r_q, center, axis)
        self.assertLess(abs(mne_value / ref27 - 1), 1e-6)  # MNE uses exactly the 27-point rule
        self.assertLess(abs(mne_value / ref_dense - 1), 1e-3)  # which approximates the cell average

    def test_finite_cell_effect_at_the_peak(self):
        """At the radial-field peak the 10-mm cell lowers the reading, less so for distant
        sources: about 1.9 % at 10 mm, 0.4 % at 15 mm and 0.04 % at 30 mm source-sensor distance."""
        q, s = np.array([0.0, 30e-9, 0.0]), 0.10
        effects = []
        for d in (0.010, 0.015, 0.030):
            r_q = s - d
            g = (s**2 + r_q**2) / (2 * s * r_q)
            th = np.arccos(3 / (g + np.sqrt(g * g + 3)))  # peak angle (paper Eq. 1)
            c = s * np.array([np.sin(th), 0.0, np.cos(th)])
            point = sphere.sarvas_field(c[None], np.array([0, 0, r_q]), q)[0] @ (c / s)
            cell = _cube_average_reference(c, c / s, np.array([0, 0, r_q]), q, opm.CELL_SIZE)
            effects.append(1 - cell / point)
        self.assertTrue(effects[0] > effects[1] > effects[2] > 0)
        self.assertAlmostEqual(effects[0], 0.0188, delta=0.001)



class TestTriaxial(unittest.TestCase):
    """A-OPM-TRIAX: three orthogonal axes per matched site over the same cell (channel-count control)."""

    def setUp(self):
        from opmsquid import g2

        rng = np.random.default_rng(4)
        pos = rng.normal(size=(4, 3)) * 0.01 + np.array([0.0, 0.0, 0.1])
        nn = rng.normal(size=(4, 3))
        nn /= np.linalg.norm(nn, axis=1, keepdims=True)
        arr = opm.OPMArray(pos, nn, pos.copy(), np.arange(4), opm.STANDOFF, 0.0, "test")
        self.matched = g2.Array("opm_matched", opm.make_info(arr), np.array(["mag"] * 4), opm.coil_def_file(), dict(standoff_mm=7.0))
        self.tri = g2.triaxial_opm(self.matched)

    def test_frames(self):
        self.assertEqual(self.tri.n, 12)
        for i, src in enumerate(self.matched.info["chs"]):
            np.testing.assert_array_equal(self.tri.info["chs"][3 * i]["loc"], src["loc"])  # the normal channel is the matched one
            for j in range(3):
                loc = self.tri.info["chs"][3 * i + j]["loc"]
                f = np.array([loc[3:6], loc[6:9], loc[9:12]])
                np.testing.assert_allclose(f @ f.T, np.eye(3), atol=1e-12)
                np.testing.assert_allclose(np.cross(f[0], f[1]), f[2], atol=1e-12)  # right-handed
                np.testing.assert_array_equal(loc[:3], src["loc"][:3])
        self.assertEqual(list(self.tri.meta["axis_role"][:3]), ["normal", "tangential_x", "tangential_y"])

    def test_uniform_field_component(self):
        from opmsquid import environment

        basis = environment.external_basis(self.tri.info, coil_def=self.tri.coil_def)[:, :3]  # unit homogeneous fields
        axes = np.array([ch["loc"][9:12] for ch in self.tri.info["chs"]])
        np.testing.assert_allclose(basis, axes, atol=1e-12)  # each channel reads the field along its own axis


class TestArrayHelpers(unittest.TestCase):
    def test_farthest_point_subset_spreads_and_keeps_spacing(self):
        rng = np.random.default_rng(0)
        pos = rng.normal(size=(300, 3))
        pos /= np.linalg.norm(pos, axis=1, keepdims=True)
        idx = opm.farthest_point_subset(pos, 120)
        self.assertEqual(len(np.unique(idx)), 120)
        self.assertEqual(int(idx[np.argmax(pos[idx, 2])]), int(np.argmax(pos[:, 2])))  # starts at the top
        self.assertGreaterEqual(opm.min_spacing(pos[idx]).min(), opm.min_spacing(pos).min())
        # every left-out point is closer to the subset than the subset's own typical spacing
        d = np.min(np.linalg.norm(pos[:, None] - pos[idx][None], axis=-1), axis=1)
        self.assertLess(d.max(), 1.5 * np.median(opm.min_spacing(pos[idx])))

if __name__ == "__main__":
    unittest.main()
