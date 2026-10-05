"""Report facts of the Methods additions (scripts/report_facts_methods.py): formatting, sources, the implementation
constants of docs/methods.md section 13 against the code lines they cite, key values against the stored results, and the
location-level sign-test values against the study's own test (opmsquid.detection.sign_flip_p)."""
import importlib.util
import json
import math
import re
import tempfile
import tomllib
import unittest
from pathlib import Path

import numpy as np

from opmsquid import detection

ROOT = Path(__file__).resolve().parents[1]
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
NEEDED = ["results/g2/g2_summary.json", "results/g2/head_surface_effect.json", "results/g2/g2_band_sensitivity.json",
          "results/g3b/g3b_summary.json", "results/g4/g4_pediatric_comparison.json", "results/g4/g4_localization_summary.json",
          "configs/g4_confirmatory.toml", "configs/g2_noise_sensitivity.toml", "configs/g3b_pediatric.toml",
          "configs/g4_motion.toml", "docs/methods.md", "docs/provenance_register.md", "docs/literature/jas2026.md"]
NEEDED += [f"results/g4/g4_{lab}_summary.json" for lab in LABELS]
NEEDED += [f"results/g4/g4_localization_{lab}_summary.json" for lab in LABELS[1:]]
HAVE = all((ROOT / p).exists() for p in NEEDED)

# Every implementation constant (docs/methods.md section 13) against the code: the line where the code sets the value (or
# the rule) must lie in a range the row cites, and must hold the token. A moved or changed line fails here.
CODE = {  # numeric rows: (file, line, the value as the code writes it)
    "IC-BG-GRID-MM": ("src/opmsquid/g2.py", 184, "0.007"),
    "IC-ENV-UNIFORM": ("src/opmsquid/environment.py", 41, "unit homogeneous-field components"),
    "IC-ENV-GRADIENT": ("src/opmsquid/environment.py", 24, "symmetric, traceless 3x3 gradient tensors"),
    "IC-ENV-TRIM-S": ("src/opmsquid/environment.py", 67, "trim_s: float = 2.0"),
    "IC-WHITEN-TOL": ("src/opmsquid/metrics.py", 71, "rel_tol: float = 1e-10"),
    "IC-TD-EXPONENT": ("src/opmsquid/background.py", 65, "exponent: float = 1.0"),
    "IC-TD-FMIN-HZ": ("src/opmsquid/background.py", 66, "f_min: float = 0.5"),
    "IC-TD-CSD-NPERSEG": ("src/opmsquid/ied.py", 41, "nperseg: int = 4096"),
    "IC-TD-PAD-S": ("src/opmsquid/ied.py", 93, "pad_s: float = 2.0"),
    "IC-BOOT-G2": ("scripts/g2_adult_comparison.py", 64, "n_boot=1000"),
    "IC-BOOT-G2-SENS": ("scripts/g2_adult_comparison.py", 372, "n_boot=200"),
    "IC-BOOT-G2-PATCH": ("scripts/g2_adult_comparison.py", 606, "rng, 200)"),
    "IC-BOOT-G2-BAND": ("scripts/g2_band_sensitivity.py", 105, "n_boot=200"),
    "IC-BOOT-HSE": ("scripts/study_head_surface_effect.py", 98, "st.rng, 1000)"),
    "IC-BOOT-G3B-MATCHED": ("scripts/g3b_pediatric_helmet.py", 610, 'if o == "opm_dense" else 200'),
    "IC-BOOT-G3B-PLACEMENTS": ("scripts/g3b_pediatric_helmet.py", 635, "rng, 200)"),
    "IC-BOOT-G3B-CHANNELS": ("scripts/g3b_pediatric_helmet.py", 667, "rng, 200)"),
    "IC-BOOT-G3B-USEFUL": ("scripts/g3b_pediatric_helmet.py", 505, "200, cfg"),
    "IC-BOOT-G4": ("scripts/g4_epilepsy_adult.py", 327, "range(1000)"),
    "IC-BOOT-LOC": ("scripts/g4_localization.py", 273, "n_boot=2000"),
    "IC-BOOT-CGAP": ("scripts/study_g3b_constant_gap.py", 861, "secondary=200"),
    "IC-BOOT-COVVAL": ("scripts/study_covariance_validation.py", 439, "default=1000"),
    "IC-SIGNFLIP-EXACT": ("src/opmsquid/detection.py", 169, "len(x) <= 20"),
    "IC-SIGNFLIP-MC": ("src/opmsquid/detection.py", 162, "n_mc=20000"),
    "IC-SEED-G3B": ("scripts/g3b_pediatric_helmet.py", 567, "default_rng(7)"),
    "IC-SEED-G3B-SCHOOL": ("scripts/g3b_pediatric_helmet.py", 568, "default_rng(8)"),
    "IC-SEED-G4-BOOT": ("scripts/g4_epilepsy_adult.py", 287, "default_rng(7)"),
    "IC-SEED-LOC-HELDOUT": ("scripts/g4_localization.py", 157, 'lc["seed"] + 1'),
    "IC-SEED-LOC-SECONDARY": ("scripts/g4_localization.py", 235, 'lc["seed"] + 2'),
    "IC-SEED-CONFIRM": ("scripts/g4_confirmatory.py", 93, 'PURPOSES = ("locations"'),
    "IC-DSPM-MNE-LIMIT": ("scripts/g4_localization.py", 171, "make_inverse_operator"),  # the value: MNE's default, below
}
RULES = {  # rule and formula rows: (file, line, a token of the rule)
    "IC-BG-AREA": ("src/opmsquid/noisemodel.py", 61, "np.bincount(owner, weights=vertex_areas"),
    "IC-BG-VAR": ("src/opmsquid/background.py", 33, "return areas"),
    "IC-COV": ("src/opmsquid/noisemodel.py", 39, "c = p @ c @ p.T"),
    "IC-PROJ": ("src/opmsquid/noisemodel.py", 49, "np.eye(len(w)) - e @ np.linalg.solve(m, e.T * w[None, :])"),
    "IC-TD-DECIMATE": ("src/opmsquid/ied.py", 114, "pad:-pad:self.decimate"),
    "IC-TD-RESPONSE": ("src/opmsquid/ied.py", 113, "y = spec.grid_gain @ q + spec.env_basis @ e + eps"),
    "IC-SEED-G3B-GRID": ("scripts/g3b_pediatric_helmet.py", 138, 'seed = g2cfg["sources"]["seed"]'),
    "IC-SEED-G4-SIM": ("scripts/g4_epilepsy_adult.py", 138, 'default_rng(cfg["simulation"]["seed"])'),
    "IC-DSPM-STUDY": ("src/opmsquid/localization.py", 33, "** (-depth)"),
    "IC-DSPM-MNE-CHS": None,  # MNE-Python's own code: checked against the installed package below
}


def cited(where: str, rel: str, line: int) -> bool:
    """Whether ``where`` ('a.py:1-3, 7; b.py:5') cites ``line`` of ``rel``."""
    for part in where.split(";"):
        part = part.strip()
        if part.startswith(rel + ":"):
            for rng in part[len(rel) + 1:].split(","):
                lo, _, hi = rng.strip().partition("-")
                if int(lo) <= line <= int(hi or lo):
                    return True
    return False


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("report_facts_methods")


class TestFormats(unittest.TestCase):
    def test_numbers(self):
        self.assertEqual(M.sig(0.0353559, 3), "0.0354")
        self.assertEqual(M.sig(0.05 / 9, 2), "0.0056")
        self.assertEqual(M.sci(2.0 ** -17), "7.6 × 10⁻⁶")
        self.assertEqual(M.sci(1 / 20000, 1), "5 × 10⁻⁵")
        self.assertEqual(M.pval(0.0308837890625), "0.031")
        self.assertEqual(M.pval(0.00131225), "0.0013")
        self.assertEqual(M.pval(2.0 ** -35), "2.9 × 10⁻¹¹")
        self.assertEqual(M.db(-0.0591), "−0.06")
        self.assertEqual(M.db(0.004), "0.00")
        self.assertEqual(M.const([0.3, 0.006, 0.3]), "0.3, 0.006 and 0.3")
        self.assertEqual(M.const(20.0), "20")
        self.assertEqual(M.count(20000.0), "20,000")
        self.assertEqual(M.seed(20261005), "20261005")
        self.assertEqual(M.span(441, 716, M.count), "441 to 716")

    def test_quote_must_occur_and_hold_the_number(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "docs").mkdir()
            (Path(d) / "docs" / "x.md").write_text("The skull is 0.25 mm\nthick at its thinnest.")
            self.assertEqual(M.quoted(Path(d), "docs/x.md", "0.25 mm thick", "0.25"), "0.25 mm thick")
            with self.assertRaises(ValueError):
                M.quoted(Path(d), "docs/x.md", "0.3 mm thick")
            with self.assertRaises(ValueError):
                M.quoted(Path(d), "docs/x.md", "0.25 mm thick", "0.3")


@unittest.skipUnless(HAVE, "result files, configurations or documents not available")
class TestMethodsFacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.F = M.facts(ROOT)
        cls.v = {k: f["value"] for k, f in cls.F.items()}
        cls.ic = M.implementation_constants(ROOT)

    def test_every_fact_has_a_value_and_an_existing_source(self):
        self.assertGreater(len(self.F), 150)
        for name, f in self.F.items():
            self.assertRegex(name, r"^meth_[a-z0-9_]+$")
            self.assertEqual(set(f), {"value", "raw", "source"}, name)
            self.assertTrue(str(f["value"]).strip(), name)
            files, sep, how = f["source"].partition(" :: ")
            self.assertTrue(sep and how.strip(), f"{name}: source lacks ' :: <key path>'")
            paths = [p for p in files.replace(",", " ").split() if p.startswith(("results/", "configs/", "docs/"))]
            self.assertTrue(paths, f"{name}: no source file")
            for p in paths:
                self.assertTrue((ROOT / p).exists(), f"{name}: {p} not found")
            self.assertNotRegex(str(f["value"]), r"(^|[\s\[(])-\d", f"{name}: ASCII minus")
            json.dumps(f["raw"])
        self.assertEqual(_load("report_facts").check(self.F, ROOT), [])

    def test_the_combined_facts_include_this_module_without_clashes(self):
        rf = _load("report_facts")
        self.assertIn("report_facts_methods", rf.MODULES)
        combined = rf.build_facts(ROOT)
        self.assertTrue(all(combined[k]["module"] == "report_facts_methods" for k in self.F))

    def test_implementation_constants_match_the_code(self):
        self.assertEqual(set(self.ic), set(CODE) | set(RULES), "every row of section 13 is checked here")
        for ident, spec in {**CODE, **RULES}.items():
            if ident in CODE:
                float(self.ic[ident]["value"])  # numeric
            else:
                self.assertIn(self.ic[ident]["value"], ("formula", "rule"), ident)
            if spec is None:
                continue
            rel, line, token = spec
            text = (ROOT / rel).read_text().splitlines()[line - 1]
            self.assertIn(token, text, f"{ident}: {rel}:{line} reads {text.strip()!r}")
            self.assertTrue(cited(self.ic[ident]["where"], rel, line), f"{ident}: {rel}:{line} not cited in "
                            f"{self.ic[ident]['where']!r}")

    def test_cited_ranges(self):
        self.assertTrue(cited("a.py:1-3, 7; b.py:5", "a.py", 2))
        self.assertTrue(cited("a.py:1-3, 7; b.py:5", "a.py", 7))
        self.assertTrue(cited("a.py:1-3, 7; b.py:5", "b.py", 5))
        self.assertFalse(cited("a.py:1-3, 7; b.py:5", "a.py", 5))
        self.assertFalse(cited("a.py:1-3, 7; b.py:5", "c.py", 1))

    def test_numeric_constants_equal_the_code_values(self):
        """The values themselves, read from the code (not only the cited line's text)."""
        ic = {k: float(r["value"]) for k, r in self.ic.items() if r["value"] not in ("formula", "rule")}
        self.assertEqual(ic["IC-BG-GRID-MM"], 7)
        self.assertEqual(ic["IC-ENV-UNIFORM"] + ic["IC-ENV-GRADIENT"], 8)
        from opmsquid import environment
        self.assertEqual(environment.N_EXT, 8)
        self.assertEqual(len(environment.GRADIENT_BASIS), ic["IC-ENV-GRADIENT"])
        for g in environment.GRADIENT_BASIS:  # symmetric and traceless
            self.assertTrue(np.allclose(g, g.T) and abs(np.trace(g)) < 1e-12)
        g4c = (ROOT / "scripts" / "g4_confirmatory.py").read_text()
        purposes = re.search(r"PURPOSES = \(([^)]*)\)", g4c).group(1)
        self.assertEqual(len(re.findall(r'"[a-z]+"', purposes)), ic["IC-SEED-CONFIRM"])
        # the sign-flip test: exact up to 20 non-zero locations, Monte Carlo (20,000 patterns) beyond
        x21 = np.ones(21)
        self.assertEqual(detection.sign_flip_p(np.ones(20)), 2.0 ** -19)
        self.assertLess(detection.sign_flip_p(x21), 0.001)  # sampled: almost no random pattern reaches the observed sum
        try:
            from mne.defaults import DEFAULTS
        except ImportError:  # pragma: no cover
            return
        self.assertEqual(DEFAULTS["depth_mne"]["limit"], ic["IC-DSPM-MNE-LIMIT"])
        self.assertIs(DEFAULTS["depth_mne"]["limit_depth_chs"], True)
        self.assertEqual(DEFAULTS["depth_mne"]["exp"], float(self.v["meth_dspm_depth"]))

    def test_key_values(self):
        expect = {
            "meth_bg_n_sources": "1,755", "meth_bg_grid_mm": "7", "meth_bg_sd_per_sqrt_mm2_nam": "0.188",
            "meth_bg_var_per_mm2_nam2": "0.0354", "meth_bg_area_cm2": "1,878", "meth_bg_mean_area_mm2": "107",
            "meth_bg_rms_sd_nam": "1.95", "meth_bg_calibration_grad_ft_cm": "37.1", "meth_enbw_hz": "35.1",
            "meth_env_terms": "8", "meth_env_r0_mm": "(0, 0, 40)", "meth_whiten_tol": "10⁻¹⁰",
            "meth_plugin60s_samples": "4,680", "meth_rank_squid_proj": "298", "meth_rank_squid_mag_proj": "99",
            "meth_rank_dense_proj": "200", "meth_rank_matched_proj": "90", "meth_opm_axis_radius_mm": "15",
            "meth_hse_radial_dense_proj_ratio_3dp": "1.073", "meth_hse_final_dense_proj_ratio_3dp": "1.115",
            "meth_hse_radial_dense_proj_ci_3dp": "[1.017, 1.115]", "meth_hse_radial_dense_proj_45_70_ratio": "0.82",
            "meth_hse_final_dense_proj_45_70_ratio": "0.98", "meth_squid_mag_asd": "3.5", "meth_squid_grad_asd": "3.6",
            "meth_squid_meas_mag_rms_ft": "25.9", "meth_coil_grad_type": "3014", "meth_coil_mag_type": "3024",
            "meth_td_fs_sim_hz": "600.6", "meth_td_fs_out_hz": "150.2", "meth_td_decimate": "4",
            "meth_td_candidates_range": "441 to 716", "meth_td_candidates_adult": "716", "meth_td_candidates_infant12mo": "441",
            "meth_null_evaluation_min": "20", "meth_loc_cov_min": "5", "meth_dspm_lambda2": "1/9", "meth_mne_depth_limit": "10",
            "meth_boot_g2": "1,000", "meth_boot_g2_sens": "200", "meth_boot_g2_band": "200", "meth_boot_g3b_dense": "1,000",
            "meth_boot_g3b_matched": "200", "meth_boot_loc": "2,000", "meth_boot_confirm": "1,000", "meth_boot_g2_parcels": "70",
            "meth_seed_g2": "2026", "meth_seed_g4_detection": "44", "meth_seed_g4_localization": "45",
            "meth_seed_loc_heldout": "46", "meth_seed_loc_secondary": "47", "meth_seed_g4_boot": "7",
            "meth_seed_confirm_root": "20261005", "meth_confirm_purposes": "10", "meth_infant12mo_skull_min_mm": "0.25",
            "meth_infant12mo_d_db": "+2.04", "meth_infant12mo_d_bem1_db": "+2.03", "meth_infant12mo_d_minus_adult_db": "+1.05",
            "meth_infant12mo_d_minus_adult_bem1_db": "+1.06", "meth_bem1_change_db_range": "−0.06 to +0.01",
            "meth_locs_per_band": "18", "meth_confirm_locs": "36", "meth_signflip_min_p_18": "7.6 × 10⁻⁶",
            "meth_signflip_min_n": "6", "meth_signtest_18_k": "14", "meth_signtest_18_k_p": "0.031",
            "meth_signtest_18_kminus1_p": "0.096", "meth_signtest_18_k_holm": "16", "meth_signtest_36_k": "25",
            "meth_signtest_36_k_holm": "27", "meth_adult_p": "0.0039", "meth_adult_equal_p": "0.0063",
            "meth_prov_fig9_asd_range": "28 to 35", "meth_prov_pack_mm": "17",
        }
        for name, value in expect.items():
            self.assertEqual(self.v[name], value, name)

    def test_raw_values_are_the_stored_ones(self):
        g2 = json.loads((ROOT / "results/g2/g2_summary.json").read_text())
        g3b = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())
        s2 = g2["noise_validation"]["brain_scale"] * 1e12  # nAm^2 per mm^2
        self.assertAlmostEqual(self.F["meth_bg_var_per_mm2_nam2"]["raw"], s2, places=12)
        self.assertAlmostEqual(math.sqrt(s2), g2["noise_validation"]["brain_source_density_nAm_per_sqrt_mm2"], places=9)
        a = g3b["anatomies"]["adult"]
        self.assertAlmostEqual(self.F["meth_bg_mean_area_mm2"]["raw"], a["cortical_area_cm2"] * 100 / a["n_background_grid"])
        self.assertEqual(a["n_background_grid"], g2["n_background_grid"])
        hse = json.loads((ROOT / "results/g2/head_surface_effect.json").read_text())["variants"]
        self.assertEqual(self.F["meth_hse_radial_dense_proj_ratio_3dp"]["raw"], hse["v4_sites_v3_axes"]["opm_dense/combined/projected"]["ratio"])
        cfg = tomllib.loads((ROOT / "configs/g4_confirmatory.toml").read_text())
        self.assertEqual(self.F["meth_seed_confirm_root"]["raw"], cfg["design"]["root_seed"])
        for lab in LABELS:
            s = json.loads((ROOT / f"results/g4/g4_{lab}_summary.json").read_text())
            self.assertEqual(s["config"]["simulation"]["seed"], self.F["meth_seed_g4_detection"]["raw"], lab)

    def test_sign_test_values_are_those_of_the_study_test(self):
        """With equal per-location differences the sign-flip p is the exact sign test's: check each stored value."""
        for n, k, name in ((18, 14, "meth_signtest_18_k_p"), (18, 13, "meth_signtest_18_kminus1_p"),
                           (18, 16, "meth_signtest_18_k_holm_p")):
            p = detection.sign_flip_p(np.array([1.0] * k + [-1.0] * (n - k)))
            self.assertAlmostEqual(self.F[name]["raw"], p, places=15, msg=name)
        self.assertEqual(self.F["meth_signflip_min_p_18"]["raw"], detection.sign_flip_p(np.full(18, 3.0)))
        self.assertEqual(self.F["meth_adult_equal_p"]["raw"], detection.sign_flip_p(np.array([1.0] * 11 + [-1.0])))
        # 36 locations: beyond the exact limit the implementation samples sign patterns; the exact values are the sign test's
        for k, name in ((25, "meth_signtest_36_k_p"), (24, "meth_signtest_36_kminus1_p"), (27, "meth_signtest_36_k_holm_p")):
            exact = 2 * sum(math.comb(36, j) for j in range(k, 37)) / 2**36
            self.assertAlmostEqual(self.F[name]["raw"], exact, places=15, msg=name)
            mc = detection.sign_flip_p(np.array([1.0] * k + [-1.0] * (36 - k)))
            self.assertLess(abs(mc - exact), 4 * math.sqrt(exact * (1 - exact) / 20000) + 1e-12, name)
        # the critical counts are the smallest that reach the level
        self.assertLess(self.F["meth_signtest_18_k_p"]["raw"], 0.05)
        self.assertGreaterEqual(self.F["meth_signtest_18_kminus1_p"]["raw"], 0.05)
        self.assertLess(self.F["meth_signtest_18_k_holm_p"]["raw"], 0.05 / 9)


if __name__ == "__main__":
    unittest.main()
