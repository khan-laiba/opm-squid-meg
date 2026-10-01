"""Motion module (G4 bounded extension): exact rigid-motion readings, linearisation, projections."""
import unittest

import numpy as np

from opmsquid import environment, motion, neuromag


def synthetic_sensors(n=60, radius=0.11, cell=True, seed=0):
    """Radial single-axis sensors on a sphere above z = 0 (8-point cells if ``cell``)."""
    rng = np.random.default_rng(seed)
    u = rng.standard_normal((n, 3))
    u[:, 2] = np.abs(u[:, 2]) + 0.3
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    centres = u * radius
    if cell:
        offs = 0.0025 * np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float)
    else:
        offs = np.zeros((1, 3))
    rmag = np.concatenate([c + offs for c in centres])
    cos = np.repeat(u, len(offs), axis=0)
    w = np.full(len(rmag), 1.0 / len(offs))
    ch = np.repeat(np.arange(n), len(offs))
    return motion.Sensors(rmag, cos, w, ch, n)


class TestMotion(unittest.TestCase):
    def setUp(self):
        self.s = synthetic_sensors()
        self.rng = np.random.default_rng(1)
        self.x_ref = np.array([0.0, 0.0, 0.04])
        self.pivot = np.zeros(3)
        self.b0, self.g = motion.random_field(self.rng, 2e-9, 5e-9)

    def test_random_field(self):
        b0, g = motion.random_field(self.rng, 3e-9, 7e-9)
        self.assertAlmostEqual(np.linalg.norm(b0), 3e-9)
        self.assertAlmostEqual(np.linalg.norm(g), 7e-9)
        np.testing.assert_allclose(g, g.T)
        self.assertAlmostEqual(np.trace(g), 0.0)

    def test_rotation_is_orthonormal(self):
        r = motion.rotation([0.1, -0.2, 0.3])
        np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(r), 1.0)
        np.testing.assert_allclose(motion.rotation([0.0, 0.0, np.pi / 2]) @ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], atol=1e-12)

    def test_jacobian_matches_finite_differences(self):
        jac = motion.jacobian(self.s, self.b0, self.g, self.x_ref, self.pivot)
        h = 1e-6
        for k in range(6):
            v = np.zeros(6)
            v[k] = h
            plus = motion.readings(self.s, self.b0, self.g, self.x_ref, motion.rotation(v[:3]), v[3:], self.pivot)
            minus = motion.readings(self.s, self.b0, self.g, self.x_ref, motion.rotation(-v[:3]), -v[3:], self.pivot)
            np.testing.assert_allclose((plus - minus) / (2 * h), jac[:, k], rtol=1e-5, atol=1e-6 * np.abs(jac).max())

    def test_jacobian_matches_finite_differences_with_pivot_and_errors(self):
        act, gains = motion.with_calibration_errors(self.s, np.random.default_rng(9), 2.0, 0.02)
        pivot = np.array([0.01, -0.02, -0.06])
        jac = motion.jacobian(act, self.b0, self.g, self.x_ref, pivot, gains)
        h = 1e-6
        for k in range(6):
            v = np.zeros(6)
            v[k] = h
            plus = motion.readings(act, self.b0, self.g, self.x_ref, motion.rotation(v[:3]), v[3:], pivot, gains)
            minus = motion.readings(act, self.b0, self.g, self.x_ref, motion.rotation(-v[:3]), -v[3:], pivot, gains)
            np.testing.assert_allclose((plus - minus) / (2 * h), jac[:, k], rtol=1e-5, atol=1e-6 * np.abs(jac).max())

    def test_any_rigid_motion_is_removed_exactly_by_the_8_term_projection(self):
        e = motion.external_basis(self.s, self.x_ref)
        p8 = motion.projector(e)
        pivot = np.array([0.0, -0.01, -0.06])
        rest = motion.readings(self.s, self.b0, self.g, self.x_ref, np.eye(3), np.zeros(3), pivot)
        for rv, tr in (([0.3, -0.5, 0.2], [0.01, -0.005, 0.02]), ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])):
            d = motion.readings(self.s, self.b0, self.g, self.x_ref, motion.rotation(rv), tr, pivot) - rest
            self.assertLess(np.linalg.norm(p8 @ d), 1e-9 * np.linalg.norm(d))
        act, gains = motion.with_calibration_errors(self.s, np.random.default_rng(4), 1.0, 0.01)
        rest = motion.readings(act, self.b0, self.g, self.x_ref, np.eye(3), np.zeros(3), pivot, gains)
        d = motion.readings(act, self.b0, self.g, self.x_ref, motion.rotation([0.01, 0.0, 0.0]), np.zeros(3), pivot, gains) - rest
        self.assertGreater(np.linalg.norm(p8 @ d), 1e-4 * np.linalg.norm(d))  # but not with calibration errors

    def test_unmodelled_and_oracle_detectability(self):
        from opmsquid import metrics

        rng = np.random.default_rng(6)
        n = 40
        a = rng.standard_normal((n, 2 * n))
        cs = a @ a.T / (2 * n) + 0.1 * np.eye(n)
        u = rng.standard_normal((n, 3))
        cu = u @ u.T
        s = rng.standard_normal((n, 25))
        thetas = [0.0, 0.1, 1.0, 10.0]
        ws = metrics.whitener(cs)
        static = np.linalg.norm(ws.apply(s), axis=0)
        un = motion.unmodelled_detectability(ws, s, cs, cu, thetas)
        orc = motion.oracle_detectability(s, cs, cu, thetas)
        np.testing.assert_allclose(un[0], static, rtol=1e-10)
        np.testing.assert_allclose(orc[0], static, rtol=1e-10)
        self.assertTrue(np.all(un <= orc + 1e-9))
        self.assertTrue(np.all(orc <= static[None] + 1e-9))
        self.assertTrue(np.all(np.diff(un, axis=0) <= 1e-12))  # more artefact, less detectability

    def test_crossing(self):
        x = [0.01, 0.1, 1.0, 10.0]
        self.assertAlmostEqual(motion.crossing(x, [0.0, -0.5, -1.5, -4.0], -1.0), 10 ** -0.5)
        self.assertIsNone(motion.crossing(x, [0.0, -0.1, -0.2, -0.3], -1.0))
        self.assertEqual(motion.crossing(x, [-2.0, -3.0, -4.0, -5.0], -1.0), 0.01)

    def test_basis_matches_environment(self):
        info = neuromag.load_info("T3")
        s = motion.sensors(info)
        np.testing.assert_allclose(motion.external_basis(s, self.x_ref), environment.external_basis(info, self.x_ref),
                                   rtol=1e-10, atol=1e-14)

    def test_uniform_field_rotation_removed_exactly_by_the_homogeneous_projection(self):
        rest = motion.readings(self.s, self.b0, 0 * self.g, self.x_ref, np.eye(3), np.zeros(3), self.pivot)
        moved = motion.readings(self.s, self.b0, 0 * self.g, self.x_ref, motion.rotation([0.2, -0.1, 0.15]), [0.01, 0, 0], self.pivot)
        p3 = motion.projector(motion.external_basis(self.s, self.x_ref)[:, :3])
        d = moved - rest
        self.assertLess(np.linalg.norm(p3 @ d), 1e-10 * np.linalg.norm(d))

    def test_gradient_translation_is_uniform_and_rotation_needs_the_gradient_terms(self):
        e = motion.external_basis(self.s, self.x_ref)
        p3, p8 = motion.projector(e[:, :3]), motion.projector(e)
        jac = motion.jacobian(self.s, 0 * self.b0, self.g, self.x_ref, self.pivot)
        for k in range(3, 6):  # translation in a gradient: a uniform change
            self.assertLess(np.linalg.norm(p3 @ jac[:, k]), 1e-10 * np.linalg.norm(jac[:, k]))
        for k in range(3):  # rotation in a gradient: a gradient change, outside the homogeneous span
            self.assertGreater(np.linalg.norm(p3 @ jac[:, k]), 0.05 * np.linalg.norm(jac[:, k]))
            self.assertLess(np.linalg.norm(p8 @ jac[:, k]), 1e-10 * np.linalg.norm(jac[:, k]))

    def test_calibration_errors_leave_a_residual_that_grows_with_the_error(self):
        e = motion.external_basis(self.s, self.x_ref)
        p8 = motion.projector(e)
        res = []
        for err in (0.5, 2.0):
            act, gains = motion.with_calibration_errors(self.s, np.random.default_rng(5), err, 0.0)
            jac = motion.jacobian(act, self.b0, self.g, self.x_ref, self.pivot, gains)
            res.append(np.linalg.norm(p8 @ jac[:, :3]) / np.linalg.norm(jac[:, :3]))
        self.assertGreater(res[0], 1e-4)
        self.assertAlmostEqual(res[1] / res[0], 4.0, delta=0.2)  # first order in the tilt

    def test_tilt_has_the_declared_rms(self):
        s = synthetic_sensors(n=2000, cell=False, seed=3)
        act, _ = motion.with_calibration_errors(s, np.random.default_rng(2), 2.0, 0.0)
        ang = np.degrees(np.arccos(np.clip(np.sum(act.cosmag * s.cosmag, axis=1), -1, 1)))
        self.assertAlmostEqual(float(np.sqrt(np.mean(ang**2))), 2.0, delta=0.1)

    def test_mismatched_detectability(self):
        rng = np.random.default_rng(3)
        a = rng.standard_normal((30, 50))
        b = a + 0.3 * rng.standard_normal((30, 50))
        np.testing.assert_allclose(motion.mismatched_detectability(a, a), np.linalg.norm(a, axis=0))
        self.assertTrue(np.all(motion.mismatched_detectability(a, b) <= np.linalg.norm(b, axis=0) + 1e-12))

    def test_motion_covariance_matches_simulation(self):
        jac = motion.jacobian(self.s, self.b0, self.g, self.x_ref, self.pivot)
        th, tr = np.radians(0.1), 1e-4
        cov = motion.motion_covariance(jac, th, tr)
        x = np.concatenate([th * self.rng.standard_normal((3, 20000)), tr * self.rng.standard_normal((3, 20000))])
        emp = np.cov(jac @ x)
        np.testing.assert_allclose(np.diag(emp), np.diag(cov), rtol=0.05)


if __name__ == "__main__":
    unittest.main()
