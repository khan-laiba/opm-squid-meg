import unittest

import mne
import numpy as np

from opmsquid import anatomy, neuromag, paths

mne.set_log_level("WARNING")

HAVE_SAMPLE = (paths.SAMPLE_MEG / neuromag.RAW_FILE).exists()


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestFullResCortex(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subject = anatomy.load_sample()
        cls.cortex = anatomy.full_resolution(cls.subject)

    def test_usable_sources_keep_clear_of_the_inner_skull_mesh(self):
        c = self.cortex
        self.assertTrue(np.all(c.valid[c.usable]))
        self.assertGreaterEqual(c.dist_inner_skull[c.usable].min(), anatomy.MIN_BEM_DISTANCE)
        dropped = np.mean(~c.usable[c.valid])
        self.assertTrue(0.05 < dropped < 0.12, dropped)  # ~8.7 % of the valid vertices (4 mm)

    def test_inner_skull_distance_is_accurate(self):
        """Distance to the subdivided mesh vs exact point-to-triangle distance on a few vertices."""
        from mne.io.constants import FIFF

        inner = next(s for s in self.subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
        rng = np.random.default_rng(0)
        idx = rng.choice(np.flatnonzero(self.cortex.dist_inner_skull < 0.006), 20, replace=False)
        pts = self.cortex.rr[idx]
        exact = np.array([_point_mesh_distance(p, inner["rr"], inner["tris"]) for p in pts])
        diff = self.cortex.dist_inner_skull[idx] - exact
        self.assertTrue(np.all(diff > -1e-6) and np.all(diff < 4e-4), diff)  # never below, at most 0.4 mm above


def _point_mesh_distance(p, rr, tris):
    """Exact distance from p to a triangle mesh (closest point on each triangle, Ericson 2005)."""
    a, b, c = rr[tris[:, 0]], rr[tris[:, 1]], rr[tris[:, 2]]
    best = np.inf
    for i in np.argsort(np.linalg.norm((a + b + c) / 3 - p, axis=1))[:200]:
        best = min(best, np.linalg.norm(p - _closest_on_triangle(p, a[i], b[i], c[i])))
    return best


def _closest_on_triangle(p, a, b, c):
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return a
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return b
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        return a + d1 / (d1 - d3) * ab
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return c
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        return a + d2 / (d2 - d6) * ac
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        return b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)
    denom = 1.0 / (va + vb + vc)
    return a + ab * vb * denom + ac * vc * denom


if __name__ == "__main__":
    unittest.main()


class TestMeshDistance(unittest.TestCase):
    def test_exact_distance_to_a_tetrahedron(self):
        from mne.io.constants import FIFF

        rr = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        tris = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])  # outward-facing
        surf = dict(rr=rr, tris=tris, np=4, ntri=4, coord_frame=FIFF.FIFFV_COORD_MRI, id=FIFF.FIFFV_BEM_SURF_ID_HEAD)
        surf = mne.surface.complete_surface_info(surf, copy=False, verbose=False)
        md = anatomy.MeshDistance(surf, k=4)
        pts = np.array([[-1.0, 0.2, 0.2], [0.1, 0.1, 0.1], [2.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        expected = [1.0, -0.1, 1.0, np.sqrt(3) * (1 - 1 / 3)]  # face x=0 (outside), inside, beyond vertex 1, above face 123
        np.testing.assert_allclose(md.signed(pts), expected, atol=1e-12)
