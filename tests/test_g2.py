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


@unittest.skipUnless(HAVE_SAMPLE, "MNE sample data not available")
class TestArrayComposition(unittest.TestCase):
    """Pins the G2 arrays on the sample head: any change to the placement rules must show up here
    (and invalidates the cached full-resolution lead fields)."""

    @classmethod
    def setUpClass(cls):
        cls.subject = anatomy.load_sample()
        cls.cortex = anatomy.full_resolution(cls.subject)
        cls.dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
        cls.arrays = g2.build_arrays(cls.subject, cls.dig)

    def _geometry(self, a):
        pos = np.array([ch["loc"][:3] for ch in a.info["chs"]])
        axis = np.array([ch["loc"][9:12] for ch in a.info["chs"]])
        return pos, axis

    def test_channel_counts_and_coil_types(self):
        types = {name: sorted({ch["coil_type"] for ch in a.info["chs"]}) for name, a in self.arrays.items()}
        self.assertEqual({name: a.n for name, a in self.arrays.items()},
                         {"squid": 306, "opm_matched": 95, "opm204": 204, "opm_dense": 205})
        self.assertEqual(types["squid"], [3014, 3024])
        self.assertEqual(sum(ch["coil_type"] == 3014 for ch in self.arrays["squid"].info["chs"]), 204)
        for name in ("opm_matched", "opm204", "opm_dense"):
            self.assertEqual(types[name], [9901])

    def test_opm_physical_placement(self):
        from scipy.spatial import cKDTree

        from opmsquid import opm

        skin = g2.skin_surface(self.subject)
        t = self.subject.trans["trans"]
        fids = opm.fiducials_head(self.dig)
        zmin = self.subject.scalp.rr[:, 2].min()
        for name in ("opm_matched", "opm204", "opm_dense"):
            pos, axis = self._geometry(self.arrays[name])
            pos_mri = mne.transforms.apply_trans(t, pos)
            d, i = cKDTree(skin.rr).query(pos_mri)
            angle = np.degrees(np.arccos(np.clip(np.sum((axis @ t[:3, :3].T) * skin.nn[i], axis=1), -1, 1)))
            ear = np.minimum(np.linalg.norm(pos - fids["lpa"], axis=1), np.linalg.norm(pos - fids["rpa"], axis=1))
            with self.subTest(name=name):
                self.assertGreaterEqual(opm.min_spacing(pos).min(), 0.017 - 1e-9 if name != "opm_matched" else 0.020)
                self.assertTrue(np.all(opm.signed_distance(pos_mri, next(
                    s for s in self.subject.bem_surfaces if s["id"] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD)) > 0.004 - 1e-6))
                self.assertLess(d.max(), 0.012)  # nominal 7 mm + at most 5 mm clearance shift
                # axis = BEM head-surface normal averaged within 15 mm (A-OPM-AXIS); against the nearest vertex's normal the
                # matched array stays within 4.4 deg and the dense one within 6.6 deg (a lower occipital site)
                self.assertLess(angle.max(), 6.0 if name == "opm_matched" else 7.5)
                self.assertGreater(ear.min(), 0.019)  # no sensor on the ear
                self.assertGreater(pos_mri[:, 2].min(), zmin)  # nothing at the MRI field-of-view cut
        self.assertLessEqual(self.arrays["opm_dense"].meta["max_extra_shift_mm"], 5.0 + 1e-9)

    def test_opm_clearance_to_the_mri_scalp(self):
        from scipy.spatial import cKDTree

        from opmsquid import opm

        tree = cKDTree(self.subject.scalp.rr)
        t = self.subject.trans["trans"]
        for name in ("opm_matched", "opm204", "opm_dense"):
            pos, _ = self._geometry(self.arrays[name])
            with self.subTest(name=name):  # A-OPM-CLEAR: standoff - 1 mm from every MRI scalp point
                self.assertGreaterEqual(tree.query(mne.transforms.apply_trans(t, pos))[0].min(), opm.STANDOFF - 0.001 - 1e-6)

    def test_cell_integration_points_outside_the_head_surface(self):
        # A-OPM-CLEAR: every integration point of every cell, as the forward model builds it (MNE
        # coil 9901, 'accurate'), >= 1 mm outside the BEM head surface by exact distance
        from mne.forward._make_forward import _create_meg_coils

        from opmsquid import anatomy, opm

        skin = next(s for s in self.subject.bem_surfaces if s["id"] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD)
        md = anatomy.MeshDistance(skin)
        t = self.subject.trans["trans"]
        for name in ("opm_matched", "opm_dense"):
            with mne.use_coil_def(opm.coil_def_file()):
                coils = _create_meg_coils(self.arrays[name].info["chs"], "accurate")
            pts = np.concatenate([c["rmag"] for c in coils])
            self.assertEqual(len(pts), 27 * self.arrays[name].n)
            # the whole cell, not only its quadrature nodes: corners, edges and faces (head frame, as the forward model)
            cells = np.concatenate([opm.cell_volume_points(ch["loc"][:3], ch["loc"][9:12]) for ch in self.arrays[name].info["chs"]])
            with self.subTest(name=name):
                self.assertGreaterEqual(md.signed(mne.transforms.apply_trans(t, pts)).min(), opm.CELL_CLEARANCE - 1e-9)
                self.assertGreaterEqual(md.signed(mne.transforms.apply_trans(t, cells)).min(), opm.CELL_CLEARANCE - 1e-9)
                # the rule samples the cell (it guarantees >= 0.29 mm in general); the exact cube-to-mesh distance
                # of these final arrays is >= 1 mm as well (1.001 mm matched, 1.005 mm dense)
                self.assertGreaterEqual(opm.exact_cell_clearance(self.arrays[name].info, skin, t).min(), opm.CELL_CLEARANCE)

    def test_stored_lead_fields_match_the_arrays(self):
        # every cached full-resolution matrix must have been computed for the arrays built here
        from opmsquid import fullres

        if not (fullres.directory() / "valid_index.npy").exists():
            self.skipTest("no full-resolution cache: run scripts/compute_fullres_forwards.py to check the stored lead fields")
        idx = np.load(fullres.directory() / "valid_index.npy")
        for job, (kind, bem_name) in fullres.JOBS.items():
            if fullres.stored_fingerprint(job) is None:
                continue
            with self.subTest(job=job):
                info = fullres.array_info(kind, self.subject)
                fp = fullres.fingerprint(info, self.subject, self.subject.bem_model(fullres.BEMS[bem_name]), idx, self.cortex)
                self.assertEqual(fullres.stored_fingerprint(job), fp)

    def test_lead_field_fingerprint_tracks_sensors_and_mesh(self):
        from opmsquid import fullres

        a = self.arrays["opm_matched"]
        idx = np.arange(10)
        bem = self.subject.bem_model(g2.BEM_CONDUCTIVITY)
        fp = fullres.fingerprint(a.info, self.subject, bem, idx, self.cortex)
        self.assertEqual(fp, fullres.fingerprint(a.info, self.subject, self.subject.bem_model(g2.BEM_CONDUCTIVITY), idx, self.cortex))
        self.assertNotEqual(fp, fullres.fingerprint(g2.with_scalp_gap(a, 1e-4).info, self.subject, bem, idx, self.cortex))
        self.assertNotEqual(fp, fullres.fingerprint(a.info, self.subject, self.subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=0), idx, self.cortex))
        self.assertEqual(len(bem[0]["tris"]), 20480)  # A-BEM-SKIN: refined head surface by default
        # a changed standard coil definition (one integration point of a built-in coil moved) changes it too
        import tempfile
        from pathlib import Path
        from unittest import mock

        from opmsquid import forward

        lines = forward.standard_coil_file().read_text().splitlines(keepends=True)
        k = next(i for i, ln in enumerate(lines) if ln.strip() and not ln.startswith("#") and len(ln.split()) == 7)
        parts = lines[k].split()
        parts[1] = repr(float(parts[1]) + 1e-4)
        lines[k] = "  ".join(parts) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            changed = Path(tmp) / "coil_def.dat"
            changed.write_text("".join(lines))
            with mock.patch.object(forward, "standard_coil_file", return_value=changed):
                self.assertNotEqual(fp, fullres.fingerprint(a.info, self.subject, bem, idx, self.cortex))

    def test_scalp_gap_variant_keeps_the_sites(self):
        a = self.arrays["opm_dense"]
        b = g2.with_scalp_gap(a, 0.003)
        pa, axis = self._geometry(a)
        pb, axis_b = self._geometry(b)
        self.assertEqual(a.n, b.n)
        np.testing.assert_allclose(pb - pa, 0.003 * axis, atol=1e-12)  # moved outward along the axes, nothing rebuilt
        np.testing.assert_array_equal(axis, axis_b)
        pa2, _ = self._geometry(a)
        np.testing.assert_array_equal(pa, pa2)  # the primary array is untouched

    def test_opm204_is_a_subset_of_the_dense_array(self):
        dense, _ = self._geometry(self.arrays["opm_dense"])
        sub, _ = self._geometry(self.arrays["opm204"])
        self.assertEqual(self.arrays["opm204"].meta["subset_of"], 205)
        self.assertTrue(all(np.any(np.all(np.isclose(dense, p, atol=1e-12), axis=1)) for p in sub))


if __name__ == "__main__":
    unittest.main()
