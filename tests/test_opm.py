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


if __name__ == "__main__":
    unittest.main()
