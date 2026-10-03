"""FreeSurfer readers without nibabel and the modelled skull (A-BEM-CHILD)."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from opmsquid import anatomy, fsio, paths

HAVE_SAMPLE = (paths.SUBJECTS_DIR / "sample" / "surf" / "lh.white").exists()


def _sphere(radius, grade=4):
    import mne

    s = mne.surface._get_ico_surface(grade)
    return s["rr"] * radius, s["tris"]


class TestTalairach(unittest.TestCase):
    def test_reads_the_linear_transform(self):
        text = ("MNI Transform File\n% avi2talxfm\n\nTransform_Type = Linear;\nLinear_Transform =\n"
                " 1.0 0.1 0.0 2.0\n 0.0 0.9 0.2 -3.0\n 0.0 0.0 1.1 4.5;\n")
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "talairach.xfm"
            f.write_text(text)
            m = fsio.read_talairach(f)
        np.testing.assert_allclose(m, [[1, 0.1, 0, 2], [0, 0.9, 0.2, -3], [0, 0, 1.1, 4.5], [0, 0, 0, 1]])

    def test_conformed_volume_maps_by_c_ras(self):
        info = dict(volume=np.array([256, 256, 256]), voxelsize=np.array([1.0, 1.0, 1.0]), xras=np.array([-1.0, 0, 0]),
                    yras=np.array([0, 0, -1.0]), zras=np.array([0, 1.0, 0]), cras=np.array([3.0, -2.0, 7.5]))
        m = fsio.surface_ras_to_scanner(info)
        np.testing.assert_allclose(m[:3, :3], np.eye(3), atol=1e-12)
        np.testing.assert_allclose(m[:3, 3], info["cras"], atol=1e-12)


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestAgainstMNE(unittest.TestCase):
    """Our reader against MNE's stored source space of the sample subject (made by MNE-C)."""

    @classmethod
    def setUpClass(cls):
        import mne

        sd = paths.SUBJECTS_DIR / "sample"
        cls.white = [fsio.read_geometry(sd / "surf" / f"{h}.white") for h in ("lh", "rh")]
        cls.sphere = [fsio.read_geometry(sd / "surf" / f"{h}.sphere") for h in ("lh", "rh")]
        cls.ref = mne.read_source_spaces(sd / "bem" / "sample-oct-6-src.fif", verbose=False)

    def test_white_surface_matches_the_stored_one(self):
        for (rr, tris), s in zip(self.white, self.ref):
            self.assertLess(np.abs(rr / 1000.0 - s["rr"]).max(), 1e-7)
            np.testing.assert_array_equal(tris, s["tris"])

    def test_oct6_selection_matches_mne_c_but_for_double_occupation(self):
        src = fsio.surface_source_space(self.white, self.sphere, "sample", "oct6")
        for ours, ref in zip(src, self.ref):
            self.assertEqual(ours["nuse"], ref["nuse"])
            differ = np.setdiff1d(ours["vertno"], ref["vertno"])
            self.assertLess(len(differ), 0.005 * ref["nuse"])  # 7 and 9 of 4,098: displaced to another free neighbour
            near = np.linalg.norm(ours["rr"][differ][:, None, :] - ref["rr"][ref["vertno"]][None], axis=2).min(axis=1)
            self.assertLess(near.max(), 0.003)
            self.assertIsNotNone(ours["pinfo"])

    def test_volume_information_is_read(self):
        _, _, info = fsio.read_geometry(paths.SUBJECTS_DIR / "sample" / "surf" / "lh.white", read_metadata=True)
        self.assertEqual(list(info["volume"]), [256, 256, 256])
        self.assertEqual(len(info["cras"]), 3)


class TestModelSkull(unittest.TestCase):
    """A failed segmentation's inner skull (just below the scalp) moved to a set depth."""

    def test_spheres(self):
        from mne.surface import complete_surface_info

        srr, stris = _sphere(0.090, 5)
        scalp = anatomy._outward(complete_surface_info(dict(rr=srr, tris=stris, np=len(srr), ntri=len(stris)), copy=False))
        irr, itris = _sphere(0.0885, 4)  # 1.5 mm below the scalp
        inner = complete_surface_info(dict(id=1, coord_frame=5, rr=irr, tris=itris, np=len(irr), ntri=len(itris)), copy=False)
        crr, _ = _sphere(0.070, 4)  # cortex 20 mm below the scalp
        new, outer, rep = anatomy.model_skull(inner, scalp, crr, depth=0.008, min_cortex=0.002)
        r_new = np.linalg.norm(new["rr"], axis=1)
        np.testing.assert_allclose(r_new, 0.082, atol=2e-4)  # 8 mm below the 90-mm scalp
        r_out = np.linalg.norm(outer["rr"], axis=1)
        np.testing.assert_allclose(r_out, 0.086, atol=2e-4)  # halfway to the scalp (4 mm of skull, 4 mm of scalp)
        self.assertAlmostEqual(rep["scalp_depth_after_mm"][1], 8.0, delta=0.3)

    def test_cortex_limit(self):
        from mne.surface import complete_surface_info

        srr, stris = _sphere(0.090, 5)
        scalp = anatomy._outward(complete_surface_info(dict(rr=srr, tris=stris, np=len(srr), ntri=len(stris)), copy=False))
        irr, itris = _sphere(0.0885, 4)
        inner = complete_surface_info(dict(id=1, coord_frame=5, rr=irr, tris=itris, np=len(irr), ntri=len(itris)), copy=False)
        crr, _ = _sphere(0.084, 5)  # cortex 6 mm below the scalp: the inner skull stops 2 mm outside it
        new, _, _ = anatomy.model_skull(inner, scalp, crr, depth=0.008, min_cortex=0.002)
        self.assertGreater(np.linalg.norm(new["rr"], axis=1).min(), 0.0855)


SCHOOL = ("sub-Z213", "sub-Z209", "sub-Z226")
HAVE_SCHOOL = all((paths.EXTERNAL / anatomy.SCHOOL_SUBJECTS / s / "bem" / f"{s}-oct-6-src.fif").exists() for s in SCHOOL)


@unittest.skipUnless(HAVE_SCHOOL, "school-aged subjects not prepared (scripts/prepare_school_subjects.py)")
class TestSchoolSubjects(unittest.TestCase):
    """The prepared children load like the templates: head surface on the scalp, nested BEM, full source space."""

    def test_load(self):
        from mne.io.constants import FIFF
        from mne.surface import _CheckInside

        for name in SCHOOL:
            with self.subTest(name=name):
                s = anatomy.load_school(name)
                self.assertTrue(s.head_conform["identity"])  # conformed when prepared
                self.assertEqual([h["nuse"] for h in s.src], [4098, 4098])
                self.assertEqual(set(s.fiducials), {"lpa", "nasion", "rpa"})
                surf = {b["id"]: b for b in s.bem_surfaces}
                self.assertTrue(_CheckInside(surf[FIFF.FIFFV_BEM_SURF_ID_SKULL])(surf[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["rr"]).all())
                self.assertTrue(_CheckInside(surf[FIFF.FIFFV_BEM_SURF_ID_HEAD])(surf[FIFF.FIFFV_BEM_SURF_ID_SKULL]["rr"]).all())
                cortex = anatomy.full_resolution(s)
                self.assertGreater(cortex.usable.mean(), 0.85)
                ear = np.linalg.norm(s.fiducials["rpa"] - s.fiducials["lpa"])
                self.assertTrue(0.10 < ear < 0.16)  # inter-auricular distance of a school-aged head [m]


if __name__ == "__main__":
    unittest.main()
