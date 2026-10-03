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



class TestHeadOnScalp(unittest.TestCase):
    """A-BEM-CONFORM: the BEM head surface's vertices moved to the nearest MRI-scalp vertices."""

    @staticmethod
    def _sphere(grade, radius, sid):
        from mne.io.constants import FIFF
        from mne.surface import _get_ico_surface

        ico = _get_ico_surface(grade)
        surf = dict(rr=ico["rr"] * radius, tris=ico["tris"], np=len(ico["rr"]), ntri=len(ico["tris"]),
                    coord_frame=FIFF.FIFFV_COORD_MRI, id=sid)
        return mne.surface.complete_surface_info(surf, copy=False, verbose=False)

    def _model(self, skull_r=0.085):
        from mne.io.constants import FIFF

        return [self._sphere(2, 0.0905, FIFF.FIFFV_BEM_SURF_ID_HEAD), self._sphere(2, skull_r, FIFF.FIFFV_BEM_SURF_ID_SKULL),
                self._sphere(2, 0.080, FIFF.FIFFV_BEM_SURF_ID_BRAIN)]

    def _scalp(self, grade=4, radius=0.090):
        s = self._sphere(grade, radius, 4)
        return anatomy.Surface(s["rr"], s["tris"], s["nn"])

    def test_vertices_move_onto_the_scalp(self):
        model = self._model()
        surfs, rep = anatomy.head_on_scalp(model, self._scalp())
        # the ico-2 directions are ico-4 vertices (MNE's icosahedra are nested; radii 1 to ~5e-5): every head vertex lands
        # on the scalp vertex in its own direction
        from scipy.spatial import cKDTree

        scalp = self._scalp()
        d, i = cKDTree(scalp.rr).query(surfs[0]["rr"])
        self.assertEqual(d.max(), 0.0)
        np.testing.assert_array_equal(i, np.arange(len(i)))
        np.testing.assert_array_equal(surfs[0]["tris"], model[0]["tris"])
        self.assertAlmostEqual(rep["moved_median_mm"], 0.5, places=4)
        self.assertFalse(rep["identity"])
        self.assertGreater(rep["outer_skull_min_mm"], 3.0)
        self.assertIs(surfs[1], model[1])  # skull and brain untouched
        self.assertIs(surfs[2], model[2])
        again, rep2 = anatomy.head_on_scalp(surfs, self._scalp())  # already on the scalp: returned as it is
        self.assertTrue(rep2["identity"])
        self.assertIs(again[0], surfs[0])

    def test_crossing_edges(self):
        s = self._sphere(2, 0.09, 4)
        self.assertEqual(anatomy.crossing_edges(s["rr"], s["tris"]), 0)
        folded = np.array(s["rr"], float)
        folded[0] = -1.2 * folded[0]  # one vertex pushed through to the far side: its edges cross the surface
        self.assertGreater(anatomy.crossing_edges(folded, s["tris"]), 0)

    def test_refusals(self):
        with self.assertRaisesRegex(ValueError, "one scalp vertex"):  # a scalp coarser than the head surface
            anatomy.head_on_scalp(self._model(), self._scalp(grade=1))
        with self.assertRaisesRegex(ValueError, "outer skull"):  # the skull would stick out of the conformed head
            anatomy.head_on_scalp(self._model(skull_r=0.0899), self._scalp(radius=0.0895))


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestSampleHeadSurface(unittest.TestCase):
    def test_sample_head_surface_on_its_scalp(self):
        from mne.io.constants import FIFF
        from scipy.spatial import cKDTree

        sub = anatomy.load_sample()
        stored = mne.read_bem_surfaces(sub.subjects_dir / "sample" / "bem" / "sample-5120-5120-5120-bem.fif", verbose=False)
        for s, f in zip(sub.bem_surfaces, stored):
            with self.subTest(surface=int(s["id"])):
                np.testing.assert_array_equal(s["tris"], f["tris"])
                if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD:  # every vertex a vertex of the MRI scalp (as the templates')
                    self.assertEqual(cKDTree(sub.scalp.rr).query(s["rr"])[0].max(), 0.0)
                else:
                    np.testing.assert_array_equal(s["rr"], f["rr"])
        rep = sub.head_conform
        self.assertFalse(rep["identity"])
        self.assertTrue(0.9 < rep["moved_median_mm"] < 1.2, rep)  # the stored outer skin lay ~1 mm off the scalp
        self.assertGreater(rep["outer_skull_min_mm"], 2.0)


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


class TestCubeMeshDistance(unittest.TestCase):
    @staticmethod
    def _plane(h=0.2, n=40):
        x = np.linspace(-h, h, n)
        gx, gy = np.meshgrid(x, x)
        rr = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
        tris = [(i * n + j, i * n + j + 1, (i + 1) * n + j + 1) for i in range(n - 1) for j in range(n - 1)]
        tris += [(i * n + j, (i + 1) * n + j + 1, (i + 1) * n + j) for i in range(n - 1) for j in range(n - 1)]
        return dict(rr=rr, tris=np.array(tris))

    def test_known_distances_to_a_plane(self):
        cm = anatomy.CubeMeshDistance(self._plane())
        half, h = 0.005, 0.02
        c, s = np.cos(np.pi / 4), np.sin(np.pi / 4)
        rot_x = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        self.assertAlmostEqual(cm.distance([0.001, 0.002, h], np.eye(3), half, 0.05), h - half, places=12)  # face down
        self.assertAlmostEqual(cm.distance([0.001, 0.002, h], rot_x.T, half, 0.05), h - half * np.sqrt(2), places=12)  # edge down
        a = np.arctan(np.sqrt(2))  # body diagonal (1, 1, 1) turned onto the z axis
        u = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
        k = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
        rot = np.eye(3) + np.sin(a) * k + (1 - np.cos(a)) * k @ k
        self.assertAlmostEqual(cm.distance([0.0, 0.0, h], rot.T, half, 0.05), h - half * np.sqrt(3), places=12)  # corner down
        self.assertEqual(cm.distance([0.0, 0.0, 0.003], np.eye(3), half, 0.05), 0.0)  # cuts the plane

    def test_exact_distance_bounds_the_sampled_one(self):
        from mne.surface import _get_ico_surface

        ico = _get_ico_surface(4)
        surf = dict(rr=ico["rr"] * 0.08, tris=ico["tris"])
        cm, md = anatomy.CubeMeshDistance(surf), anatomy.MeshDistance(surf)
        rng = np.random.default_rng(0)
        half = 0.005
        u = np.linspace(-1, 1, 41) * half
        ga, gb = np.meshgrid(u, u)
        faces = []
        for ax in range(3):
            for sgn in (-1, 1):
                q = np.zeros((ga.size, 3))
                q[:, ax] = sgn * half
                q[:, [i for i in range(3) if i != ax]] = np.column_stack([ga.ravel(), gb.ravel()])
                faces.append(q)
        faces = np.concatenate(faces)
        for _ in range(10):
            n = rng.normal(size=3)
            n /= np.linalg.norm(n)
            frame = np.linalg.qr(rng.normal(size=(3, 3)))[0].T
            centre = n * (0.08 + rng.uniform(0.001, 0.005) + half * np.abs(frame @ n).sum())
            sampled = md.unsigned(centre + faces @ frame).min()
            exact = cm.distance(centre, frame, half, sampled)
            self.assertLessEqual(exact, sampled + 1e-12)
            self.assertGreater(exact, sampled - 0.25 * 2 * half / 40)  # within the sampling resolution


if __name__ == "__main__":
    unittest.main()
