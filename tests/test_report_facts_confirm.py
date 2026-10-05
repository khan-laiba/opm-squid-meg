"""Report facts of the confirmatory spike run (scripts/report_facts_confirm.py).

Always: the formatting of censored ratios, unreached S50 and flags, and a clear FileNotFoundError when results/g4_confirm
(or its combined summary) is missing. When results/g4_confirm/g4_confirm_summary.json exists: every fact has a value and
an existing source, the combiner's check() is clean, names and formats follow the conventions of the other fact modules,
the Holm adjustment of the endpoint follows from the stored p values, and key values equal the formatted value at their
place in the result files, given as (fact name, file and key path, formatting) triples: no number is written into this
test. The facts are read from the repository root, or from the root named by REPORT_FACTS_CONFIRM_ROOT (e.g. a pilot run
on test seeds in a temporary tree whose results/g4_confirm holds its files)."""
import importlib.util
import json
import os
import re
import tempfile
import tomllib
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("REPORT_FACTS_CONFIRM_ROOT") or REPO)
DIR = ROOT / "results" / "g4_confirm"
HAVE = (DIR / "g4_confirm_summary.json").exists()


def _load(name):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("report_facts_confirm")
FMT = {"pval": M.pval, "count": M.count, "ratio": M.ratio, "integer": M.integer, "share": M.share, "mm": M.mm,
       "const": M.const, "prob": M.prob, "yes": M.yes, "str": str, "len": lambda x: M.count(len(x)),
       "ratio_ci": lambda c: M.iv(c[0], c[1], M.ratio), "integer_ci": lambda c: M.iv(c[0], c[1], M.integer),
       "share_ci": lambda c: M.iv(c[0], c[1], M.share), "prob_ci": lambda c: M.iv(c[0], c[1], M.prob)}
CENSORED = re.compile(r"^(> |< )\d+\.\d\d$|^not estimable$")

# (fact name, file, key path, formatting); placeholders: {lab} anatomy key, {t} its name token, {ep} the endpoint's stored
# comparison, {dense} / {matched} the stored pair names, {rt} the endpoint's rate tag, {base} the endpoint without scope
PER_ANATOMY = [
    ("cf_{t}_dense_vs_combined_practical_p", "g4c_{lab}_summary.json", ("comparisons", "{ep}", "location_sign_flip_p"), "pval"),
    ("cf_{t}_dense_vs_combined_practical_locs_opm", "g4c_{lab}_summary.json", ("comparisons", "{ep}", "locations_favouring_opm"),
     "count"),
    ("cf_{t}_dense_vs_combined_practical_locs_squid", "g4c_{lab}_summary.json",
     ("comparisons", "{ep}", "locations_favouring_squid"), "count"),
    ("cf_{t}_dense_vs_combined_practical_events_opm_only", "g4c_{lab}_summary.json", ("comparisons", "{ep}", "detected_only_opm"),
     "count"),
    ("cf_{t}_dense_vs_combined_practical_ratio", "g4c_{lab}_summary.json",
     ("comparisons", "{ep}", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_dense_vs_combined_practical_ratio_ci", "g4c_{lab}_summary.json",
     ("comparisons", "{ep}", "s50_ratio_squid_over_opm", "ci95"), "ratio_ci"),
    ("cf_{t}_dense_vs_combined_practical_p_holm", "g4_confirm_summary.json", ("families", "endpoint", "holm_p", "{lab}"), "pval"),
    ("cf_{t}_dense_vs_combined_practical_p_holm", "g4_confirm_summary.json", ("anatomy", "{lab}", "endpoint", "holm_p"), "pval"),
    ("cf_{t}_dense_vs_combined_practical_ratio", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "endpoint", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_s50_combined_practical", "g4c_{lab}_summary.json", ("s50", "squid/combined/practical@{rt}/replicate0", "value"),
     "integer"),
    ("cf_{t}_s50_dense_practical", "g4c_{lab}_summary.json", ("s50", "opm_dense/opm/practical@{rt}/replicate0", "value"),
     "integer"),
    ("cf_{t}_s50_dense_practical", "g4_confirm_summary.json", ("anatomy", "{lab}", "s50_dense_nAm", "value"), "integer"),
    ("cf_{t}_s50_combined_practical_ci", "g4c_{lab}_summary.json", ("s50", "squid/combined/practical@{rt}/replicate0", "ci95"),
     "integer_ci"),
    ("cf_{t}_dense_vs_combined_oracle_ratio", "g4c_{lab}_summary.json",
     ("comparisons", "{dense}/oracle/replicate0", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_dense_vs_combined_oracle_ratio_ci", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "oracle", "s50_ratio_squid_over_opm", "ci95"), "ratio_ci"),
    ("cf_{t}_dense_vs_combined_oracle_p_holm", "g4_confirm_summary.json",
     ("families", "secondary/{dense}/oracle/replicate0", "holm_p", "{lab}"), "pval"),
    ("cf_{t}_dense_vs_combined_mismatch_ratio", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "mismatch", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_dense_vs_combined_mismatch_p", "g4c_{lab}_summary.json",
     ("comparisons", "{dense}/mismatch@{rt}/replicate0", "location_sign_flip_p"), "pval"),
    ("cf_{t}_dense_vs_combined_matchedrate_ratio", "g4c_{lab}_summary.json",
     ("comparisons", "{dense}/practical@{rt}_matched/replicate0", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_matched_vs_combined_practical_ratio", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "matched_array", "s50_ratio_squid_over_opm", "value"), "ratio"),
    ("cf_{t}_dense_vs_combined_practical_pooled_ratio", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "pooled", "value"), "ratio"),
    ("cf_{t}_dense_vs_combined_practical_pooled_p_holm", "g4_confirm_summary.json",
     ("families", "secondary/{base}/pooled", "holm_p", "{lab}"), "pval"),
    ("cf_{t}_mismatch_cost_dense", "g4c_{lab}_summary.json", ("mismatch", "mismatch/cost/opm_dense/opm/mismatch@{rt}/replicate0",
                                                              "value"), "ratio"),
    ("cf_{t}_mismatch_cost_combined_ci", "g4_confirm_summary.json", ("anatomy", "{lab}", "mismatch_cost", "squid/combined", "ci95"),
     "ratio_ci"),
    ("cf_{t}_mismatch_ratio_change_dense", "g4_confirm_summary.json", ("anatomy", "{lab}", "mismatch_ratio_change", "value"),
     "ratio"),
    ("cf_{t}_mc_n_p05", "g4c_{lab}_summary.json", ("monte_carlo", "{base}", "n_p_below_alpha"), "count"),
    ("cf_{t}_rate_dense_practical_evaluation_matched", "g4c_{lab}_summary.json",
     ("false_events", "opm_dense/opm|primary", "evaluation_matched", "rate_per_min"), "share"),
    ("cf_{t}_rate_dense_practical_evaluation_matched", "g4_confirm_summary.json",
     ("anatomy", "{lab}", "false_events", "opm_dense/opm|primary", "evaluation_matched"), "share"),
    ("cf_{t}_rate_combined_practical_heldout_frozen_ci", "g4c_{lab}_summary.json",
     ("false_events", "squid/combined|primary", "heldout_frozen", "ci95"), "share_ci"),
    ("cf_{t}_rate_combined_mismatch_evaluation_frozen", "g4c_{lab}_summary.json",
     ("false_events", "squid/combined|mismatch", "evaluation_frozen", "rate_per_min"), "share"),
    ("cf_{t}_oracle_fpp_dense", "g4c_{lab}_summary.json", ("oracle", "opm_dense/opm", "heldout_false_positive_probability"),
     "prob"),
    ("cf_{t}_oracle_fpp_dense_ci", "g4c_{lab}_summary.json", ("oracle", "opm_dense/opm", "ci95"), "prob_ci"),
    ("cf_{t}_median_depth_mm", "g4c_{lab}_summary.json", ("location_checks", "depth_mm_median"), "mm"),
    ("cf_{t}_n_dictionary", "g4c_{lab}_summary.json", ("n_dictionary",), "count"),
]
GLOBAL = [
    ("cf_confirmatory", "g4_confirm_summary.json", ("confirmatory",), "yes"),
    ("cf_complete", "g4_confirm_summary.json", ("complete",), "yes"),
    ("cf_n_anatomies", "g4_confirm_summary.json", ("anatomies",), "len"),
    ("cf_n_commits", "g4_confirm_summary.json", ("simulated_at_commits",), "len"),
    ("cf_alpha", "g4_confirm_summary.json", ("endpoint", "definition", "alpha"), "const"),
    ("cf_operating_point_per_minute", "g4_confirm_summary.json", ("endpoint", "definition", "false_events_per_min"), "const"),
    ("cf_band_lo", "g4_confirm_summary.json", ("endpoint", "definition", "depth_band_mm", 0), "const"),
    ("cf_n_anatomies_declared", "g4_confirm_summary.json", ("endpoint", "definition", "anatomies"), "len"),
    ("cf_dense_vs_combined_practical_n_holm_pass", "g4_confirm_summary.json", ("families", "endpoint", "n_pass"), "count"),
    ("cf_dense_vs_combined_oracle_n_holm_pass", "g4_confirm_summary.json",
     ("families", "secondary/{dense}/oracle/replicate0", "n_pass"), "count"),
    ("cf_dense_vs_combined_mismatch_n_holm_pass", "g4_confirm_summary.json",
     ("families", "secondary/{dense}/mismatch@{rt}/replicate0", "n_pass"), "count"),
    ("cf_matched_vs_combined_practical_n_holm_pass", "g4_confirm_summary.json",
     ("families", "secondary/{matched}/practical@{rt}/replicate0", "n_pass"), "count"),
    ("cf_mc_n_replicates", "g4_confirm_summary.json", ("monte_carlo_summary", "n_replicates"), "count"),
    ("cf_mc_n_replicates_all_pass", "g4_confirm_summary.json", ("monte_carlo_summary", "n_replicates_all_pass"), "count"),
    ("cf_endpoint_code_check_exact", "g4c_endpoint_code_check.json", ("all_exact",), "yes"),
    ("cf_root_seed", "configs/g4_confirmatory.toml", ("design", "root_seed"), "str"),
    ("cf_n_replicates", "configs/g4_confirmatory.toml", ("design", "noise_replicates"), "count"),
    ("cf_bootstrap_resamples", "configs/g4_confirmatory.toml", ("design", "bootstrap_resamples"), "count"),
    ("cf_null_evaluation_minutes_declared", "configs/g4_confirmatory.toml", ("null", "evaluation_min"), "const"),
]


def _doc(name: str, cache: dict):
    if name not in cache:
        path = ROOT / name if "/" in name else DIR / name
        cache[name] = tomllib.loads(path.read_text()) if name.endswith(".toml") else json.loads(path.read_text())
    return cache[name]


def _at(doc, path):
    for k in path:
        doc = doc[k]
    return doc


class TestFormats(unittest.TestCase):
    def test_censored_and_uncensored_ratios(self):
        F = M.Facts()
        F.ratio("cf_a_ratio", {"value": None, "value_bounds": [1.2345, None], "value_censored": "Neuromag does not reach 50 %",
                               "ci95": [1.1, None]}, "results/x.json :: k")
        F.ratio("cf_b_ratio", {"value": None, "value_bounds": [None, 0.8], "ci95": [None, 0.95]}, "results/x.json :: k")
        F.ratio("cf_c_ratio", {"value": None, "value_bounds": None, "ci95": [None, None]}, "results/x.json :: k")
        F.ratio("cf_d_ratio", {"value": 1.0234, "ci95": [0.95, 1.11]}, "results/x.json :: k")
        F.ratio("cf_e_ratio", {"value": 1.3643, "ci95": [1.1016, 1.5756]}, "results/x.json :: k")
        v = {k: x["value"] for k, x in F.items()}
        self.assertEqual((v["cf_a_ratio"], v["cf_a_ratio_ci"]), ("> 1.23", "[1.10, open]"))
        self.assertEqual((v["cf_b_ratio"], v["cf_b_ratio_ci"]), ("< 0.80", "[open, 0.95]"))
        self.assertEqual((v["cf_c_ratio"], v["cf_c_ratio_ci"]), ("not estimable", "[open, open]"))
        self.assertEqual((v["cf_d_ratio"], v["cf_d_ratio_3dp"], v["cf_d_ratio_ci_3dp"]), ("1.02", "1.023", "[0.950, 1.110]"))
        self.assertEqual((v["cf_e_ratio"], v["cf_e_ratio_ci"]), ("1.36", "[1.10, 1.58]"))
        self.assertNotIn("cf_e_ratio_3dp", v)
        self.assertIn("Neuromag does not reach 50 %", F["cf_a_ratio"]["source"])

    def test_s50_flags_and_names(self):
        F = M.Facts(" [run not confirmatory: test]")
        F.s50("cf_a_s50", {"value": None, "ci95": [241.46, None]}, "results/x.json :: s50")
        F.s50("cf_b_s50", {"value": 46.505, "ci95": [36.2, 63.1]}, "results/x.json :: s50")
        self.assertEqual((F["cf_a_s50"]["value"], F["cf_a_s50_ci"]["value"]), ("not reached", "[241, open]"))
        self.assertEqual((F["cf_b_s50"]["value"], F["cf_b_s50_ci"]["value"]), ("47", "[36, 63]"))
        self.assertTrue(F["cf_b_s50"]["source"].endswith("[run not confirmatory: test]"))
        self.assertEqual((M.yes(True), M.yes(False), M.listing([]), M.names(["adult", "childB"])), ("yes", "no", "none",
                                                                                                   "adult and child B"))
        self.assertEqual((M.pval(0.0), M.pval(4.5e-5), M.pval(0.0103), M.pval(1.0)), ("<0.0001", "<0.0001", "0.010", "1.0"))
        with self.assertRaises(ValueError):
            F.add("g4_wrong_prefix", "1", 1, "results/x.json :: k")
        with self.assertRaises(ValueError):
            F.add("cf_no_key_path", "1", 1, "results/x.json")
        with self.assertRaises(ValueError):
            F.add("cf_b_s50", "1", 1, "results/x.json :: k")  # defined twice

    def test_missing_results_raise_a_clear_error(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError) as cm:
                M.facts(Path(d))
            self.assertIn("results/g4_confirm", str(cm.exception))
            (Path(d) / "results" / "g4_confirm").mkdir(parents=True)
            with self.assertRaises(FileNotFoundError) as cm:
                M.facts(Path(d))
            self.assertIn("g4_confirm_summary.json", str(cm.exception))


@unittest.skipUnless(HAVE, "results/g4_confirm/g4_confirm_summary.json not available (the confirmatory run is not in this tree)")
class TestConfirmFacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f = M.facts(ROOT)
        cls.docs = {}
        cf = _doc("g4_confirm_summary.json", cls.docs)
        ep = cf["endpoint"]["comparison"]
        base = ep.rsplit("/", 1)[0]
        rt = base.rsplit("@", 1)[1]
        cls.subst = dict(ep=ep, base=base, rt=rt, dense="opm_dense/opm_vs_squid/combined",
                         matched="opm_matched/opm_vs_squid/combined")
        cls.labs = list(cf["anatomies"])

    def _expand(self, x, **kw):
        if isinstance(x, str):
            return x.format(**self.subst, **kw)
        if isinstance(x, tuple):
            return tuple(self._expand(y, **kw) for y in x)
        return x

    def test_every_fact_has_a_value_and_an_existing_source(self):
        self.assertGreater(len(self.f), 500)
        for name, fact in self.f.items():
            self.assertEqual(set(fact), {"value", "raw", "source"}, name)
            self.assertNotEqual(str(fact["value"]).strip(), "", name)
            files, sep, how = fact["source"].partition(" :: ")
            self.assertTrue(sep and how.strip(), f"{name}: source lacks ' :: <key path>'")
            paths = [p for p in files.replace(",", " ").split() if p.startswith(("results/", "configs/", "docs/"))]
            self.assertTrue(paths, f"{name}: no source file")
            for p in paths:
                self.assertTrue((ROOT / p).exists(), f"{name}: {p} not found")
        self.assertEqual(_load("report_facts").check(self.f, ROOT), [])  # the combiner's own check

    def test_names_and_formats(self):
        conf = _doc("g4_confirm_summary.json", self.docs)["confirmatory"]
        for name, fact in self.f.items():
            v = str(fact["value"])
            self.assertRegex(name, r"^cf_[a-z0-9_]+$")
            self.assertNotRegex(v, r"(^|[\s\[(])-\d", f"{name}: ASCII minus in {v!r}")
            if name.endswith(("_ci", "_ci_3dp")):
                self.assertRegex(v, r"^\[\S+, \S+\]$", name)
            if name.endswith("_range"):
                self.assertTrue(" to " in v or v == self.f[name[:-len("_range")] + "_min"]["value"], name)
            json.dumps(fact["raw"])  # plain JSON types only
            self.assertEqual("[run not confirmatory" in fact["source"], not conf, name)

    def test_values_at_their_json_paths(self):
        checked = 0
        rows = [(n, f, p, fmt, {}) for n, f, p, fmt in GLOBAL]
        rows += [(n, f, p, fmt, dict(lab=lab, t=lab.lower())) for lab in self.labs for n, f, p, fmt in PER_ANATOMY]
        for name, file, path, fmt, kw in rows:
            name, file, path = self._expand(name, **kw), self._expand(file, **kw), self._expand(path, **kw)
            want = _at(_doc(file, self.docs), path)
            got = self.f[name]["value"]
            if want is None:  # a censored point estimate or an S50 beyond the tested strengths
                self.assertTrue(CENSORED.match(got) or got == "not reached", f"{name}: {got!r} for a null value")
            else:
                self.assertEqual(got, FMT[fmt](want), f"{name} vs {file} :: {path}")
            checked += 1
        self.assertEqual(checked, len(GLOBAL) + len(PER_ANATOMY) * len(self.labs))

    def test_endpoint_holm_follows_from_the_stored_p(self):
        p = {lab: _doc(f"g4c_{lab}_summary.json", self.docs)["comparisons"][self.subst["ep"]]["location_sign_flip_p"]
             for lab in self.labs}
        adj = M.holm(p)
        alpha = _doc("g4_confirm_summary.json", self.docs)["endpoint"]["definition"]["alpha"]
        for lab in self.labs:
            fact = self.f[f"cf_{lab.lower()}_dense_vs_combined_practical_p_holm"]
            self.assertAlmostEqual(fact["raw"], adj[lab], places=12)
            self.assertEqual(self.f[f"cf_{lab.lower()}_dense_vs_combined_practical_holm_pass"]["raw"], adj[lab] < alpha)
        n = sum(v < alpha for v in adj.values())
        self.assertEqual(self.f["cf_dense_vs_combined_practical_n_holm_pass"]["raw"], n)

    def test_ranges_come_from_unrounded_values(self):
        vals = [self.f[f"cf_{lab.lower()}_dense_vs_combined_practical_ratio"]["raw"] for lab in self.labs]
        vals = [v for v in vals if isinstance(v, float)]
        if vals:
            self.assertEqual(self.f["cf_all_dense_vs_combined_practical_ratio_min"]["raw"], min(vals))
            self.assertEqual(self.f["cf_all_dense_vs_combined_practical_ratio_max"]["raw"], max(vals))


if __name__ == "__main__":
    unittest.main()
