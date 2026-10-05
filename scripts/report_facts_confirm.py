#!/usr/bin/env python3
"""Report facts for the confirmatory spike run (prefix cf_): scripts/g4_confirmatory.py, whose endpoint, seeds, sample
size, replicates and mismatch variant were declared in configs/g4_confirmatory.toml before the run. Read from
results/g4_confirm/ (g4_confirm_summary.json of the combine step, g4c_<anatomy>_summary.json, g4c_endpoint_code_check.json),
checked against results/g4/g4_pediatric_comparison.json (the exploratory values the run is set beside) and the declared
configuration; nothing is simulated or re-estimated. Not yet merged by scripts/report_facts.py: add it to MODULES there
once results/g4_confirm is committed.

facts(root) -> {name: {"value": text as printed, "raw": unrounded number(s) or text, "source": "results/<file> :: <key
path>" or "results/<file> :: derived: <how>"}}. Raises FileNotFoundError when results/g4_confirm or one of its files is
missing, and ValueError when the files disagree with one another (a stale combined summary, Holm values that do not follow
from the stored p values) or when a run marked confirmatory does not carry the declared configuration. A run that is not
confirmatory (incomplete, test settings, several or uncommitted commits, several configurations) still yields its facts,
with "[run not confirmatory: <reasons>]" appended to every source.

Formats: those of scripts/report_facts_g4.py, whose helpers are imported. S50 in nAm as integers ("not reached" when the
location-pooled curve stays below 50 %); ratios 2 decimals with _3dp twins when 0.95-1.05; a censored ratio as its bound
("> 1.13", "< 0.90"; "not estimable" when neither system reaches 50 %); intervals "[lo, hi]" in _ci facts, "open" for an
open end (bootstrap resamples outside the tested strengths); p-values 2 significant digits down to 0.0001, "<0.0001" below;
false-event rates per minute 2 decimals; probabilities 2 significant digits; mm 1 decimal; percentages integer (1 decimal
below 10 %); counts with thousands separators; configuration constants as declared; flags "yes"/"no"; lists "a, b and c";
ranges "lo to hi" from unrounded values (one value when both ends print alike) with _min and _max twins; U+2212 for
negatives. "raw" holds the unrounded value (a fraction for _pct facts).

Names. Anatomies adult, school, size2yr, infant2yr, infant18mo, infant12mo, childa, childb, childc; groups over anatomies
all, smaller_heads, templates_controls (scaled adults and infant templates), childrenabc. Arrays dense and matched (the
dense and the site-matched OPM arrays), always against Neuromag's 306 channels (combined); S50 also for combined.
Detectors: practical (the endpoint's scanning detector, thresholds frozen at 1 false event per minute on the calibration
null), matchedrate (the same detector, thresholds matched to 1 false event per minute on the held-out null; its realized
rate is evaluated on the independent evaluation null), oracle (known topography, waveform and time), mismatch and
mismatch_matchedrate (the declared detector-mismatch variant with frozen and with matched thresholds). Scope: no suffix =
noise replicate 0 (the confirmatory realization), _pooled = all noise replicates pooled, _rep<r> = replicate r (Monte Carlo).
  cf_<design>                              declared design, status, commit, seeds, sizes, null durations, code check
  cf_<anatomy>_<array>_vs_combined_<detector>[_pooled|_rep<r>]_<what>   per anatomy: ratio (S50 Neuromag / OPM), ratio_ci,
        reduction_pct (1 - 1/ratio), locs_opm, locs_squid, locs_tied, events_opm_only, events_squid_only, p (sign-flip
        over locations, uncorrected), p_holm (Holm over the anatomies within the family), holm_pass
        (the endpoint: cf_<anatomy>_dense_vs_combined_practical_*)
  cf_<array>_vs_combined_<detector>[_pooled]_<count>   per family: n_holm_pass, holm_pass_anatomies, holm_fail_anatomies,
        n_p05, n_ratio_censored, n_ratio_above1, n_ratio_below1, n_ci_above1, n_ci_below1, n_ci_includes1
  cf_<group>_<array>_vs_combined_<detector>[_pooled]_<what>_min|_max|_range   ratio, reduction_pct, p, p_holm, locs_*
  cf_<anatomy>_s50_<combined|matched|dense>_<detector>[_pooled] (+ _ci), cf_<group>_s50_..._range
  cf_mc_*, cf_<anatomy>_mc_*               Monte Carlo variability of the endpoint over the noise replicates
  cf_<anatomy>_mismatch_cost_<array>..., cf_<anatomy>_mismatch_ratio_change_<array>...   the mismatch variant: each
        array's S50 with the mismatched over the primary detector, and the change of the S50 ratio (Neuromag's cost over
        the OPM's), frozen or matched thresholds, replicate 0 or pooled; cf_mismatch_* the declared variant
  cf_<anatomy>_rate_<array>_<practical|mismatch>_<heldout_frozen|evaluation_frozen|evaluation_matched> (+ _ci)
        realized false events per minute; cf_rate_* ranges; cf_*_rate_equality_* equal-rate tests; cf_*_oracle_fpp_*
        the oracle's per-trial false-positive probability on the held-out null

Usage: PYTHONPATH=src .venv/bin/python scripts/report_facts_confirm.py [--root DIR] [--grep TEXT]
       (prints name, value, source)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIR = "results/g4_confirm"
CF = f"{DIR}/g4_confirm_summary.json"
PER = DIR + "/g4c_{}_summary.json"
CHECK = f"{DIR}/g4c_endpoint_code_check.json"
EXACT = f"{DIR}/g4_confirm_exact_p.json"  # exact sign-flip p of every comparison (scripts/g4_confirm_exact_p.py)
EXACT_NOTE = (f" [p: the exact sign-flip p and its Holm adjustment, {EXACT} (from the stored location_differences; "
              "the run's Monte Carlo values in the field named here agree and give the same Holm decisions)]")
SECTIONS = {"oracle": "opm_dense/opm_vs_squid/combined/oracle/replicate0",  # the combined summary's per-anatomy copies
            "mismatch": "opm_dense/opm_vs_squid/combined/mismatch@1/replicate0",
            "matched_array": "opm_matched/opm_vs_squid/combined/practical@1/replicate0"}
EXPLO = "results/g4/g4_pediatric_comparison.json"
CFG = "configs/g4_confirmatory.toml"
BASE_CFG = "configs/g4_epilepsy.toml"


def _load_g4():
    spec = importlib.util.spec_from_file_location("report_facts_g4", Path(__file__).with_name("report_facts_g4.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G4 = _load_g4()  # the formatting helpers (and anatomy groups) of the exploratory spike facts
MINUS, ratio, ratio3, integer, count, share, pct, sig, pval, const, iv, span, mm = (
    G4.MINUS, G4.ratio, G4.ratio3, G4.integer, G4.count, G4.share, G4.pct, G4.sig, G4.pval, G4.const, G4.iv, G4.span, G4.mm)
holm = G4._holm

ANATOMIES = G4.LABELS
GROUPS = {"all": ANATOMIES, "smaller_heads": ANATOMIES[1:], "templates_controls": G4.TEMPLATES_CONTROLS,
          "childrenabc": G4.CHILDREN}
NAMES = {"adult": "adult", "school": "school-age size", "size2yr": "2-year size", "infant2yr": "24-month template",
         "infant18mo": "18-month template", "infant12mo": "12-month template", "childA": "child A", "childB": "child B",
         "childC": "child C"}  # prose names (Table 1 of the report)
COMPARATOR = "squid/combined"
ARRAYS = {"opm_dense/opm": "dense", "opm_matched/opm": "matched"}
DETSETS = {"squid/combined": "combined", "opm_matched/opm": "matched", "opm_dense/opm": "dense"}
SCOPES = {"replicate0": "", "pooled": "_pooled"}
RATE_KINDS = {"heldout_frozen": "frozen thresholds (calibration null) on the held-out null",
              "evaluation_frozen": "frozen thresholds on the independent evaluation null",
              "evaluation_matched": "thresholds matched on the held-out null, on the independent evaluation null"}
LOBES = ("frontal", "temporal", "parietal", "occipital", "cingulate", "insula")
CI_NOTE = "95 % location bootstrap, every resample kept; open = the resamples' bound lies outside the tested strengths"


# ----------------------------------------------------------------------------------------------
# formatting beyond the imported helpers
def yes(b) -> str:
    return "yes" if b else "no"


def listing(items) -> str:
    """'a, b and c'; 'none' for an empty list."""
    items = [str(x) for x in items]
    return "none" if not items else G4._join(items)


def names(labs) -> str:
    return listing([NAMES[lab] for lab in labs])


def prob(x) -> str:
    return sig(x, 2)


def tag(lab: str) -> str:
    return lab.lower()


def _files(labs) -> str:
    return ", ".join(PER.format(lab) for lab in labs)


class Facts(dict):
    """name -> {"value", "raw", "source"}: names cf_<lower_snake>, a ' :: <key path>' in every source, ``note`` (the
    run's status when it is not confirmatory) appended to every source."""

    def __init__(self, note: str = ""):
        super().__init__()
        self.note = note

    def add(self, name: str, value: str, raw, source: str) -> None:
        if name in self:
            raise ValueError(f"fact {name!r} defined twice")
        if not re.fullmatch(r"cf_[a-z0-9_]+", name):
            raise ValueError(f"bad fact name {name!r}")
        if str(value).strip() == "":
            raise ValueError(f"fact {name!r} has an empty value")
        if " :: " not in source:
            raise ValueError(f"fact {name!r}: source lacks ' :: <key path>'")
        self[name] = dict(value=value, raw=raw, source=source + self.note)

    def spread(self, name: str, values, f, source: str) -> None:
        """name_min, name_max and name_range (from unrounded values)."""
        values = list(values)
        lo, hi = min(values), max(values)
        self.add(f"{name}_min", f(lo), lo, f"{source} (minimum)")
        self.add(f"{name}_max", f(hi), hi, f"{source} (maximum)")
        self.add(f"{name}_range", span(values, f), [lo, hi], f"{source} (minimum to maximum)")

    def ratio(self, name: str, sr: dict, source: str) -> float | None:
        """A paired S50 ratio (or a ratio of them) with its interval as <name>_ci, and _3dp twins near 1; a censored
        point estimate is printed as its bound. Returns the point value (None when censored)."""
        v, ci = sr.get("value"), list(sr.get("ci95") or [None, None])
        ci_src = f"{source}.ci95 ({CI_NOTE})"
        if v is not None:
            self.add(name, ratio(v), v, f"{source}.value")
            self.add(f"{name}_ci", iv(ci[0], ci[1], ratio), ci, ci_src)
            if 0.95 <= v <= 1.05:
                self.add(f"{name}_3dp", ratio3(v), v, f"{source}.value")
                self.add(f"{name}_ci_3dp", iv(ci[0], ci[1], ratio3), ci, ci_src)
            return v
        lo, hi = sr.get("value_bounds") or [None, None]
        why = sr.get("value_censored") or "an S50 outside the tested strengths"
        if lo is not None and hi is None:
            txt, raw = "> " + ratio(lo), [lo, None]
        elif hi is not None and lo is None:
            txt, raw = "< " + ratio(hi), [None, hi]
        else:
            txt, raw = "not estimable", "not estimable"
        self.add(name, txt, raw, f"{source}.value_bounds (point estimate censored: {why})")
        self.add(f"{name}_ci", iv(ci[0], ci[1], ratio), ci, ci_src)
        return None

    def s50(self, name: str, v: dict, source: str) -> None:
        if v["value"] is None:
            self.add(name, "not reached", "not reached", f"{source}.value (null: 50 % detection not reached at the strongest "
                     "tested strength)")
        else:
            self.add(name, integer(v["value"]), v["value"], f"{source}.value (nAm, location-pooled detection curve, log "
                     "interpolation)")
        ci = list(v["ci95"])
        self.add(f"{name}_ci", iv(ci[0], ci[1], integer), ci, f"{source}.ci95 (95 % location bootstrap; open = beyond the "
                 "tested strengths)")


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def _above1(sr: dict) -> bool:
    """Point estimate above 1, or a censored estimate whose lower bound is at least 1."""
    v, b = sr.get("value"), sr.get("value_bounds") or [None, None]
    return v > 1 if v is not None else (b[0] is not None and b[0] >= 1)


def _below1(sr: dict) -> bool:
    v, b = sr.get("value"), sr.get("value_bounds") or [None, None]
    return v < 1 if v is not None else (b[1] is not None and b[1] <= 1)


def _ci_side(sr: dict) -> str:
    """'above' (lower end above 1), 'below' (upper end below 1) or 'includes' (open ends unbounded)."""
    lo, hi = sr.get("ci95") or [None, None]
    return "above" if lo is not None and lo > 1 else ("below" if hi is not None and hi < 1 else "includes")


# ----------------------------------------------------------------------------------------------
# inputs and status
def _read(root: Path) -> dict:
    if not (root / DIR).is_dir():
        raise FileNotFoundError(f"{DIR} not found in {root}: the confirmatory spike run (scripts/g4_confirmatory.py) has no "
                                "results in this tree yet")
    if not (root / CF).is_file():
        raise FileNotFoundError(f"{CF} not found in {root}: the combine step (scripts/g4_confirmatory.py --combine-only) "
                                "has not run")
    cf = _json(root / CF)
    labs = list(cf["anatomies"])
    unknown = [lab for lab in labs if lab not in ANATOMIES]
    if unknown:
        raise ValueError(f"{CF}: unknown anatomies {unknown}")
    missing = [PER.format(lab) for lab in labs] + [CHECK, EXPLO, CFG, BASE_CFG]
    missing = [p for p in missing if not (root / p).is_file()]
    if missing:
        raise FileNotFoundError(f"confirmatory spike facts: {', '.join(missing)} not found in {root}")
    labs = [lab for lab in ANATOMIES if lab in labs]
    return dict(cf=cf, labs=labs, per={lab: _older_pilot(_json(root / PER.format(lab))) for lab in labs},
                check=_json(root / CHECK), explo=_json(root / EXPLO)["comparison"],
                cfg=tomllib.loads((root / CFG).read_text()), base=tomllib.loads((root / BASE_CFG).read_text()))


def _older_pilot(s: dict) -> dict:
    """Descriptive keys that summaries of early test runs (scripts/g4_confirmatory.py before the declared design was final)
    lack, filled with the equivalent values the current code writes. A summary marked confirmatory is never filled: a key
    missing there fails."""
    if s["confirmatory"]:
        return s
    dc, mc = s["declared_choices"], next(iter(s["monte_carlo"].values()))
    s["seeds"].setdefault("declared_root_seed", s["config"]["confirmatory"]["design"]["root_seed"])
    dc.setdefault("noise_replicates_run", mc["n_replicates"])
    dc.setdefault("locations_per_stratum_run", max(p["drawn"] for p in s["location_pools"]))
    dc.setdefault("sign_flip_test", "exact (all sign patterns)" if s["n_locations"] <= 20 else
                  "Monte Carlo, 20,000 sign patterns (detection.sign_flip_p, more than 20 locations)")
    for v in dc["variants"].values():
        v.setdefault("n_templates", len(v["stretches"]))
    s.setdefault("mismatch_geometry", {})
    return s


def _status(D: dict) -> list[str]:
    """Why the run is not confirmatory (empty when it is), by the combine step's rule; must agree with its flag."""
    cf, per, labs = D["cf"], D["per"], D["labs"]
    reasons = []
    if not cf["complete"] or len(labs) != len(ANATOMIES):
        reasons.append(f"{len(labs)} of {len(ANATOMIES)} anatomies")
    tests = [lab for lab in labs if not per[lab]["confirmatory"]]
    if tests:
        reasons.append("test settings (" + ", ".join(tests) + ")")
    commits = cf["simulated_at_commits"]
    if len(commits) > 1:
        reasons.append("several commits (" + ", ".join(commits) + ")")
    if any("+dirty" in c for c in commits):
        reasons.append("uncommitted changes")
    if len(cf["config_digests"]) > 1:
        reasons.append("several configurations")
    if bool(cf["confirmatory"]) != (not reasons):
        raise ValueError(f"{CF}: confirmatory = {cf['confirmatory']}, but the combine rule gives {not reasons} ({reasons})")
    if sorted({per[lab]["simulated_at_commit"] for lab in labs}) != sorted(commits):
        raise ValueError(f"{CF}: simulated_at_commits differ from the per-anatomy summaries (stale combined summary?)")
    return reasons


def _consistency(D: dict) -> None:
    """The combined summary is that of these per-anatomy summaries (not stale) and states their endpoint."""
    cf, per = D["cf"], D["per"]
    ep = cf["endpoint"]["comparison"]
    for lab in D["labs"]:
        a, r = cf["anatomy"][lab], per[lab]["comparisons"][ep]
        if (a["endpoint"]["location_sign_flip_p"] != r["location_sign_flip_p"]
                or a["endpoint"]["s50_ratio_squid_over_opm"] != r["s50_ratio_squid_over_opm"]
                or (a["endpoint"]["locations_favouring_opm"], a["endpoint"]["locations_favouring_squid"])
                != (r["locations_favouring_opm"], r["locations_favouring_squid"])
                or a["simulated_at_commit"] != per[lab]["simulated_at_commit"]):
            raise ValueError(f"{CF}: anatomy['{lab}'] differs from {PER.format(lab)} (rerun the combine step)")
        if per[lab]["endpoint"]["comparison"] != ep or per[lab]["endpoint"]["definition"] != cf["endpoint"]["definition"]:
            raise ValueError(f"{PER.format(lab)}: endpoint differs from the combined summary's")


def _as_declared(D: dict, lab: str) -> bool:
    c = D["per"][lab]["config"]
    return c["confirmatory"] == D["cfg"] and all(c["inherited"][k] == D["base"][k] for k in c["inherited"])


def _common(F: Facts, name: str, values: dict, f, source: str, confirmatory: bool) -> None:
    """One value shared by every anatomy (or anatomy and detector); a range when they differ, an error for a confirmatory
    run."""
    vals = sorted(set(values.values()), key=lambda x: (str(type(x)), x))
    if len(vals) == 1:
        F.add(name, f(vals[0]), vals[0], f"{source} (the same in every anatomy)")
    elif confirmatory:
        raise ValueError(f"{source}: differs between the anatomies of a confirmatory run: {values}")
    else:
        F.add(name, span(vals, f), vals, f"{source} (differs between anatomies: {values})")


# ----------------------------------------------------------------------------------------------
# design, status, provenance
def _design(F: Facts, D: dict, reasons: list[str]) -> None:
    cf, per, labs = D["cf"], D["per"], D["labs"]
    conf = not reasons
    first = labs[0]
    P0, ALL = PER.format(first), _files(labs)
    ep = cf["endpoint"]["definition"]
    F.add("cf_confirmatory", yes(conf), conf, f"{CF} :: confirmatory (all nine anatomies, no test settings, one commit "
          "without uncommitted changes, one configuration)")
    F.add("cf_status", "confirmatory" if conf else "not confirmatory", conf, f"{CF} :: confirmatory")
    F.add("cf_not_confirmatory_reasons", listing(reasons), reasons, f"{CF} :: derived: complete, anatomy[*].confirmatory, "
          "simulated_at_commits, config_digests (the combine step's rule)")
    F.add("cf_complete", yes(cf["complete"]), cf["complete"], f"{CF} :: complete")
    F.add("cf_n_anatomies", count(len(labs)), len(labs), f"{CF} :: anatomies (count)")
    F.add("cf_n_anatomies_declared", count(len(ep["anatomies"])), len(ep["anatomies"]),
          f"{CF} :: endpoint.definition.anatomies (count; the Holm family declared in {CFG})")
    commits = cf["simulated_at_commits"]
    F.add("cf_commit", listing(commits), commits, f"{CF} :: simulated_at_commits (code commit of every anatomy's simulation)")
    F.add("cf_n_commits", count(len(commits)), len(commits), f"{CF} :: simulated_at_commits (count)")
    dg = cf["config_digests"]
    F.add("cf_config_digest", listing(dg), dg, f"{CF} :: config_digests (sha1 of {CFG} and {BASE_CFG}, first 12 hex digits)")
    decl = {lab: _as_declared(D, lab) for lab in labs}
    if conf and not all(decl.values()):
        raise ValueError(f"run marked confirmatory, but its stored configuration differs from {CFG} / {BASE_CFG} for "
                         f"{[lab for lab, ok in decl.items() if not ok]}")
    F.add("cf_config_as_declared", yes(all(decl.values())), decl, f"{ALL} :: config.confirmatory, config.inherited; derived: "
          f"equal to {CFG} and to the [simulation], [events], [detector] and [null] tables of {BASE_CFG}")

    # the endpoint as declared
    src = f"{CF} :: endpoint.definition"
    F.add("cf_alpha", const(ep["alpha"]), ep["alpha"], f"{src}.alpha (declared in {CFG} [endpoint])")
    F.add("cf_operating_point_per_minute", const(ep["false_events_per_min"]), ep["false_events_per_min"],
          f"{src}.false_events_per_min (thresholds frozen on the calibration null)")
    lo, hi = ep["depth_band_mm"]
    F.add("cf_band", f"{const(lo)}–{const(hi)}", [lo, hi], f"{src}.depth_band_mm (mm below the scalp)")
    F.add("cf_band_lo", const(lo), lo, f"{src}.depth_band_mm[0]")
    F.add("cf_band_hi", const(hi), hi, f"{src}.depth_band_mm[1]")
    F.add("cf_endpoint_replicate", str(ep["replicate"]), ep["replicate"], f"{src}.replicate (the noise replicate of the "
          "endpoint: the first realization of every event)")
    F.add("cf_endpoint_family", ep["family"], ep["family"], f"{src}.family")

    # seeds (declared values from the configuration file; a confirmatory run carries exactly these, checked above)
    design = D["cfg"]["design"]
    F.add("cf_root_seed", str(design["root_seed"]), design["root_seed"], f"{CFG} :: design.root_seed (declared before the run)")
    used = {lab: per[lab]["seeds"]["root_seed"] for lab in labs}
    ok = all(v == per[lab]["seeds"]["declared_root_seed"] == design["root_seed"] for lab, v in used.items())
    F.add("cf_root_seed_as_declared", yes(ok), used, f"{ALL}, {CFG} :: seeds.root_seed of every anatomy against "
          "seeds.declared_root_seed and design.root_seed")
    purposes = per[first]["seeds"]["purposes"]
    F.add("cf_n_seed_purposes", count(len(purposes)), len(purposes), f"{P0} :: seeds.purposes (count: one random stream per "
          "anatomy and purpose)")
    F.add("cf_seed_purposes", listing(purposes), purposes, f"{P0} :: seeds.purposes")
    F.add("cf_seed_rule", per[first]["seeds"]["rule"], per[first]["seeds"]["rule"], f"{P0} :: seeds.rule")

    # sizes
    dc = {lab: per[lab]["declared_choices"] for lab in labs}
    F.add("cf_n_replicates", count(design["noise_replicates"]), design["noise_replicates"], f"{CFG} :: design.noise_replicates "
          "(declared; replicate 0 is the confirmatory realization)")
    _common(F, "cf_n_replicates_run", {lab: dc[lab]["noise_replicates_run"] for lab in labs}, count,
            f"{ALL} :: declared_choices.noise_replicates_run", conf)
    per_stratum = design.get("locations_per_stratum") or D["base"]["events"]["n_locations"] // 12
    F.add("cf_locations_per_stratum", count(per_stratum), per_stratum, f"{CFG} :: design.locations_per_stratum (declared "
          f"before the run; without it, {BASE_CFG} events.n_locations / 12)")
    _common(F, "cf_locations_per_stratum_run", {lab: dc[lab]["locations_per_stratum_run"] for lab in labs}, count,
            f"{ALL} :: declared_choices.locations_per_stratum_run", conf)
    _common(F, "cf_n_orientation_strata", {lab: len(per[lab]["location_pools"]) for lab in labs}, count,
            f"{ALL} :: location_pools (count: orientation strata of the band, 0-30, 30-60, 60-90 deg)", conf)
    _common(F, "cf_n_locations", {lab: per[lab]["n_locations"] for lab in labs}, count, f"{ALL} :: n_locations", conf)
    ev = per[first]["config"]["inherited"]["events"]
    per_loc = len(ev["strengths_nAm"]) * len(ev["stretches"])
    F.add("cf_n_events_per_location", count(per_loc), per_loc, f"{P0} :: config.inherited.events; derived: strengths_nAm x "
          "stretches (one focal event each, per noise replicate)")
    _common(F, "cf_n_events_per_replicate", {lab: per[lab]["n_events_per_replicate"] for lab in labs}, count,
            f"{ALL} :: n_events_per_replicate", conf)
    n_ev = {lab: per[lab]["n_events_per_replicate"] * dc[lab]["noise_replicates_run"] for lab in labs}
    _common(F, "cf_n_events_per_anatomy", n_ev, count, f"{ALL} :: derived: n_events_per_replicate x "
            "declared_choices.noise_replicates_run", conf)
    F.add("cf_n_events_total", count(sum(n_ev.values())), sum(n_ev.values()), f"{ALL} :: derived: n_events_per_replicate x "
          "declared_choices.noise_replicates_run, summed over the anatomies")
    F.add("cf_strength_range", f"{const(ev['strengths_nAm'][0])} to {const(ev['strengths_nAm'][-1])}",
          [ev["strengths_nAm"][0], ev["strengths_nAm"][-1]], f"{P0} :: config.inherited.events.strengths_nAm (first to last)")
    F.add("cf_n_strengths", count(len(ev["strengths_nAm"])), len(ev["strengths_nAm"]),
          f"{P0} :: config.inherited.events.strengths_nAm (count)")
    for k in ("baseline", "calibration", "heldout", "evaluation"):
        _common(F, f"cf_null_{k}_minutes", {lab: per[lab]["null_minutes"][k] for lab in labs}, const,
                f"{ALL} :: null_minutes.{k}", conf)
    tot = {lab: sum(per[lab]["null_minutes"].values()) for lab in labs}
    _common(F, "cf_null_total_minutes", tot, const, f"{ALL} :: null_minutes; derived: baseline + calibration + held-out + "
            "evaluation", conf)
    F.add("cf_bootstrap_resamples", count(design["bootstrap_resamples"]), design["bootstrap_resamples"],
          f"{CFG} :: design.bootstrap_resamples (location bootstrap of every S50 and ratio interval)")
    F.add("cf_null_evaluation_minutes_declared", const(D["cfg"]["null"]["evaluation_min"]), D["cfg"]["null"]["evaluation_min"],
          f"{CFG} :: null.evaluation_min (declared: the fourth, independent null set)")
    tests = {lab: dc[lab]["sign_flip_test"] for lab in labs}
    _common(F, "cf_sign_flip_test", tests, str, f"{ALL} :: declared_choices.sign_flip_test", conf)
    m = re.match(r"Monte Carlo, ([\d,]+) sign patterns", tests[first])
    F.add("cf_sign_flip_method", "Monte Carlo" if m else "exact", tests[first], f"{P0} :: declared_choices.sign_flip_test "
          "(exact up to 20 locations, Monte Carlo beyond: opmsquid.detection.sign_flip_p)")
    if m:
        n = int(m.group(1).replace(",", ""))
        F.add("cf_sign_flip_n_patterns", count(n), n, f"{P0} :: declared_choices.sign_flip_test (random sign patterns)")
    cal = {f"{lab}|{k}": v["calibration_events_above_frozen"] for lab in labs for k, v in per[lab]["false_events"].items()}
    _common(F, "cf_calibration_events_above_threshold", cal, count, f"{ALL} :: false_events[*].calibration_events_above_frozen "
            "(calibration null events above each frozen threshold: rate x calibration minutes)", conf)
    alphas = {f"{lab}|{k}": o["alpha"] for lab in labs for k, o in per[lab]["oracle"].items()}
    _common(F, "cf_oracle_alpha", alphas, const, f"{ALL} :: oracle[*].alpha (the oracle's per-trial false-positive "
            "probability on the calibration null)", conf)

    # the endpoint code check (this run's statistics on the exploratory states)
    ck = D["check"]
    avail = {lab: v for lab, v in ck["anatomies"].items() if v.get("available", True) and "error" not in v}
    cmp_ = [c for v in avail.values() for c in v.values() if isinstance(c, dict) and "exact_match" in c]
    F.add("cf_endpoint_code_check_exact", yes(ck["all_exact"]), ck["all_exact"], f"{CHECK} :: all_exact (this run's code, "
          "applied to the exploratory run's stored states, reproduces its location counts, sign-flip p and S50-ratio point "
          "estimates exactly)")
    F.add("cf_endpoint_code_check_n_anatomies", count(len(avail)), len(avail), f"{CHECK} :: anatomies[*] (available: "
          f"{names([k for k in ANATOMIES if k in avail])})")
    F.add("cf_endpoint_code_check_n_comparisons", count(len(cmp_)), len(cmp_), f"{CHECK} :: anatomies[*][*].exact_match "
          "(comparisons checked)")
    F.add("cf_endpoint_code_check_n_exact", count(sum(c["exact_match"] for c in cmp_)), sum(c["exact_match"] for c in cmp_),
          f"{CHECK} :: anatomies[*][*].exact_match (true)")


# ----------------------------------------------------------------------------------------------
# endpoint and secondary families, per anatomy and over anatomies
def _modes(D: dict) -> dict:
    """Stored detector labels -> name tokens (the endpoint's rate, then the declared variants)."""
    rt = f"{D['cf']['endpoint']['definition']['false_events_per_min']:g}"
    out = {f"practical@{rt}": "practical", f"practical@{rt}_matched": "matchedrate", "oracle": "oracle"}
    for v in D["per"][D["labs"][0]]["config"]["confirmatory"].get("variant", []):
        name = v["name"]
        if not re.fullmatch(r"[a-z][a-z0-9]*", name):
            raise ValueError(f"variant name {name!r} is not usable in fact names")
        out[f"{name}@{rt}"] = name
        out[f"{name}@{rt}_matched"] = f"{name}_matchedrate"
    return out


def _family_key(D: dict, name: str) -> str:
    return "endpoint" if name == D["cf"]["endpoint"]["comparison"] else f"secondary/{name}"


def _holm_checked(D: dict, fam_key: str, ps: dict) -> dict:
    """The stored Holm-adjusted p of a family, checked against the stored p values and a recomputation."""
    fam = D["cf"]["families"].get(fam_key) if not fam_key.startswith("replicate") else D["cf"]["monte_carlo"].get(fam_key)
    where = f"families['{fam_key}']" if not fam_key.startswith("replicate") else f"monte_carlo['{fam_key}']"
    if fam is None:
        raise ValueError(f"{CF}: {where} missing")
    if fam["p"] != ps:
        raise ValueError(f"{CF}: {where}.p differs from the per-anatomy summaries (stale combined summary?)")
    re_ = holm(ps)
    if any(abs(re_[k] - fam["holm_p"][k]) > 1e-12 for k in ps):
        raise ValueError(f"{CF}: {where}.holm_p is not the Holm adjustment of its p values")
    alpha = D["cf"]["endpoint"]["definition"]["alpha"]
    if fam["n_pass"] != sum(v < alpha for v in fam["holm_p"].values()):
        raise ValueError(f"{CF}: {where}.n_pass does not count the Holm-adjusted p below alpha")
    return fam["holm_p"]


def _comparison_facts(F: Facts, D: dict, lab: str, name: str, base: str, hp: float, hp_src: str) -> float | None:
    r = D["per"][lab]["comparisons"][name]
    key = f"{PER.format(lab)} :: comparisons['{name}']"
    n_loc, o, s = r["n_locations"], r["locations_favouring_opm"], r["locations_favouring_squid"]
    alpha = D["cf"]["endpoint"]["definition"]["alpha"]
    F.add(f"{base}_locs_opm", count(o), o, f"{key}.locations_favouring_opm (of {n_loc} locations)")
    F.add(f"{base}_locs_squid", count(s), s, f"{key}.locations_favouring_squid (of {n_loc} locations)")
    F.add(f"{base}_locs_tied", count(n_loc - o - s), n_loc - o - s, f"{key}; derived: n_locations - locations_favouring_opm - "
          "locations_favouring_squid (equal detection counts)")
    F.add(f"{base}_events_opm_only", count(r["detected_only_opm"]), r["detected_only_opm"],
          f"{key}.detected_only_opm (of {r['n']} focal events)")
    F.add(f"{base}_events_squid_only", count(r["detected_only_squid"]), r["detected_only_squid"],
          f"{key}.detected_only_squid (of {r['n']} focal events)")
    F.add(f"{base}_p", pval(r["location_sign_flip_p"]), r["location_sign_flip_p"], f"{key}.location_sign_flip_p (two-sided "
          "sign-flip test on the per-location differences in detection counts, uncorrected)")
    F.add(f"{base}_p_holm", pval(hp), hp, hp_src)
    F.add(f"{base}_holm_pass", yes(hp < alpha), hp < alpha, f"{hp_src}; derived: below endpoint.definition.alpha "
          f"({const(alpha)})")
    rv = F.ratio(f"{base}_ratio", r["s50_ratio_squid_over_opm"], f"{key}.s50_ratio_squid_over_opm")
    if rv is not None:
        F.add(f"{base}_reduction_pct", pct(1 - 1 / rv), 1 - 1 / rv, f"{key}.s50_ratio_squid_over_opm.value; derived: "
              "1 - 1/ratio (the OPM's S50 is lower by this share)")
    return rv


def _family_summary(F: Facts, D: dict, rows: dict, holm_p: dict, fb: str, gb: str, src: str, hsrc: str) -> None:
    """Counts over the anatomies of one family (fb) and ranges over every anatomy group (gb with {group})."""
    labs = D["labs"]
    alpha = D["cf"]["endpoint"]["definition"]["alpha"]
    sr = {lab: rows[lab]["s50_ratio_squid_over_opm"] for lab in labs}
    passed = [lab for lab in labs if holm_p[lab] < alpha]
    F.add(f"{fb}_n_holm_pass", count(len(passed)), len(passed), f"{hsrc}; derived: anatomies with Holm-adjusted p below "
          f"{const(alpha)} (of {len(labs)})")
    F.add(f"{fb}_holm_pass_anatomies", names(passed), passed, f"{hsrc}; derived: anatomies with Holm-adjusted p below "
          f"{const(alpha)}")
    F.add(f"{fb}_holm_fail_anatomies", names([lab for lab in labs if lab not in passed]), [lab for lab in labs if lab not in passed],
          f"{hsrc}; derived: anatomies with Holm-adjusted p at or above {const(alpha)}")
    n05 = [lab for lab in labs if rows[lab]["location_sign_flip_p"] < 0.05]
    F.add(f"{fb}_n_p05", count(len(n05)), len(n05), f"{src}.location_sign_flip_p; derived: anatomies with uncorrected p < 0.05: "
          + names(n05))
    cens = [lab for lab in labs if sr[lab]["value"] is None]
    F.add(f"{fb}_n_ratio_censored", count(len(cens)), len(cens), f"{src}.s50_ratio_squid_over_opm.value; derived: anatomies whose "
          "point estimate is censored (null value): " + names(cens))
    for nm, test, what in (("n_ratio_above1", _above1, "above 1 (a censored estimate: its lower bound at least 1)"),
                           ("n_ratio_below1", _below1, "below 1 (a censored estimate: its upper bound at most 1)")):
        sel = [lab for lab in labs if test(sr[lab])]
        F.add(f"{fb}_{nm}", count(len(sel)), len(sel), f"{src}.s50_ratio_squid_over_opm; derived: anatomies with the ratio {what}: "
              + names(sel))
    for side in ("above", "below", "includes"):
        sel = [lab for lab in labs if _ci_side(sr[lab]) == side]
        what = {"above": "lower end above 1", "below": "upper end below 1", "includes": "including 1 (open ends unbounded)"}[side]
        F.add(f"{fb}_n_ci_{side}1", count(len(sel)), len(sel), f"{src}.s50_ratio_squid_over_opm.ci95; derived: anatomies with the "
              f"interval's {what}: " + names(sel))
    for g, members in GROUPS.items():
        sel = [lab for lab in members if lab in labs]
        if not sel:
            continue
        b = gb.format(group=g)
        gsrc = src.replace(_files(labs), _files(sel))
        vals = [sr[lab]["value"] for lab in sel if sr[lab]["value"] is not None]
        if vals:
            F.spread(f"{b}_ratio", vals, ratio, f"{gsrc}.s50_ratio_squid_over_opm.value over {names(sel)} "
                     f"({len(vals)} uncensored)")
            F.spread(f"{b}_reduction_pct", [1 - 1 / x for x in vals], pct, f"{gsrc}.s50_ratio_squid_over_opm.value; "
                     f"derived: 1 - 1/ratio over {names(sel)} ({len(vals)} uncensored)")
        F.spread(f"{b}_p", [rows[lab]["location_sign_flip_p"] for lab in sel], pval, f"{gsrc}.location_sign_flip_p over "
                 f"{names(sel)} (uncorrected)")
        F.spread(f"{b}_p_holm", [holm_p[lab] for lab in sel], pval, f"{hsrc} over {names(sel)}")
        F.spread(f"{b}_locs_opm", [rows[lab]["locations_favouring_opm"] for lab in sel], count,
                 f"{gsrc}.locations_favouring_opm over {names(sel)}")
        F.spread(f"{b}_locs_squid", [rows[lab]["locations_favouring_squid"] for lab in sel], count,
                 f"{gsrc}.locations_favouring_squid over {names(sel)}")
        if g != "all":
            n = sum(holm_p[lab] < alpha for lab in sel)
            F.add(f"{b}_n_holm_pass", count(n), n, f"{hsrc}; derived: anatomies of {names(sel)} with Holm-adjusted p below "
                  f"{const(alpha)}")


def _families(F: Facts, D: dict) -> None:
    labs, per = D["labs"], D["per"]
    ALL = _files(labs)
    for a, at in ARRAYS.items():
        for mode, mt in _modes(D).items():
            for scope, st in SCOPES.items():
                name = f"{a}_vs_{COMPARATOR}/{mode}/{scope}"
                rows = {lab: per[lab]["comparisons"][name] for lab in labs}
                fam_key = _family_key(D, name)
                hp = _holm_checked(D, fam_key, {lab: rows[lab]["location_sign_flip_p"] for lab in labs})
                hsrc = (f"{CF} :: families['{fam_key}'].holm_p (Holm over the {len(labs)} anatomies"
                        + (", the declared endpoint)" if fam_key == "endpoint" else " within this secondary family)"))
                for lab in labs:
                    _comparison_facts(F, D, lab, name, f"cf_{tag(lab)}_{at}_vs_combined_{mt}{st}", hp[lab],
                                      hsrc.replace(".holm_p (", f".holm_p['{lab}'] ("))
                _family_summary(F, D, rows, hp, f"cf_{at}_vs_combined_{mt}{st}", f"cf_{{group}}_{at}_vs_combined_{mt}{st}",
                                f"{ALL} :: comparisons['{name}']", hsrc)
    # S50 of every array, detector and scope
    for key, dt in DETSETS.items():
        for mode, mt in _modes(D).items():
            for scope, st in SCOPES.items():
                k = f"{key}/{mode}/{scope}"
                for lab in labs:
                    F.s50(f"cf_{tag(lab)}_s50_{dt}_{mt}{st}", per[lab]["s50"][k], f"{PER.format(lab)} :: s50['{k}']")
                for g, members in GROUPS.items():
                    sel = [lab for lab in members if lab in labs]
                    vals = [per[lab]["s50"][k]["value"] for lab in sel if per[lab]["s50"][k]["value"] is not None]
                    if vals:
                        F.spread(f"cf_{g}_s50_{dt}_{mt}{st}", vals, integer, f"{_files(sel)} :: s50['{k}'].value over "
                                 f"{names(sel)} ({len(vals)} reached 50 %)")


# ----------------------------------------------------------------------------------------------
# Monte Carlo variability of the endpoint (independent noise realizations of the same events)
def _monte_carlo(F: Facts, D: dict) -> None:
    cf, per, labs = D["cf"], D["per"], D["labs"]
    alpha = cf["endpoint"]["definition"]["alpha"]
    ep = cf["endpoint"]["comparison"]
    base_name = ep.rsplit("/", 1)[0]
    mcs = cf["monte_carlo_summary"]
    n_rep = mcs["n_replicates"]
    src = f"{CF} :: monte_carlo_summary"
    F.add("cf_mc_n_replicates", count(n_rep), n_rep, f"{src}.n_replicates (noise replicates present in every anatomy)")
    npr = mcs["n_pass_per_replicate"]
    F.add("cf_mc_n_pass_per_replicate", listing([count(x) for x in npr]), npr, f"{src}.n_pass_per_replicate (anatomies "
          "passing Holm over the anatomies in each noise replicate, replicate 0 first)")
    F.spread("cf_mc_n_pass", npr, count, f"{src}.n_pass_per_replicate")
    F.add("cf_mc_n_replicates_all_pass", count(mcs["n_replicates_all_pass"]), mcs["n_replicates_all_pass"],
          f"{src}.n_replicates_all_pass (replicates in which every anatomy passes)")
    holm_r = {}
    for r in range(n_rep):
        ps = {lab: per[lab]["comparisons"][f"{base_name}/replicate{r}"]["location_sign_flip_p"] for lab in labs}
        holm_r[r] = _holm_checked(D, f"replicate{r}", ps)
        if cf["monte_carlo"][f"replicate{r}"]["n_pass"] != npr[r]:
            raise ValueError(f"{CF}: monte_carlo_summary.n_pass_per_replicate[{r}] differs from monte_carlo['replicate{r}']")
    if holm_r[0] != cf["families"]["endpoint"]["holm_p"]:
        raise ValueError(f"{CF}: monte_carlo['replicate0'] is not the endpoint family")
    # every replicate of the endpoint, per anatomy (replicate 0 is the endpoint itself: cf_<anatomy>_dense_vs_combined_practical_*)
    at, mt = ARRAYS[base_name.split("_vs_")[0]], _modes(D)[base_name.rsplit("/", 1)[1]]
    for r in range(1, n_rep):
        hsrc = f"{CF} :: monte_carlo['replicate{r}'].holm_p (Holm over the {len(labs)} anatomies within noise replicate {r})"
        for lab in labs:
            _comparison_facts(F, D, lab, f"{base_name}/replicate{r}", f"cf_{tag(lab)}_{at}_vs_combined_{mt}_rep{r}",
                              holm_r[r][lab], hsrc.replace(".holm_p (", f".holm_p['{lab}'] ("))
        F.add(f"cf_{at}_vs_combined_{mt}_rep{r}_n_holm_pass", count(npr[r]), npr[r], f"{CF} :: monte_carlo['replicate{r}'].n_pass")
    # per anatomy, over its replicates
    allv, sds, n05, cells, pass_cells, every, never = [], [], 0, 0, 0, [], []
    for lab in labs:
        t = tag(lab)
        m = per[lab]["monte_carlo"][base_name]
        key = f"{PER.format(lab)} :: monte_carlo['{base_name}']"
        vals = [v for v in m["s50_ratio"] if v is not None]
        if vals:
            F.spread(f"cf_{t}_mc_ratio", vals, ratio, f"{key}.s50_ratio over its {m['n_replicates']} replicates "
                     f"({len(vals)} uncensored)")
        F.add(f"cf_{t}_mc_n_censored", count(m["n_censored"]), m["n_censored"], f"{key}.n_censored")
        F.spread(f"cf_{t}_mc_p", m["location_sign_flip_p"], pval, f"{key}.location_sign_flip_p over its {m['n_replicates']} "
                 "replicates (uncorrected)")
        F.add(f"cf_{t}_mc_n_p05", count(m["n_p_below_alpha"]), m["n_p_below_alpha"], f"{key}.n_p_below_alpha (replicates with "
              f"uncorrected p < {const(alpha)}, of {m['n_replicates']})")
        k = sum(holm_r[r][lab] < alpha for r in range(n_rep))
        F.add(f"cf_{t}_mc_n_holm_pass", count(k), k, f"{CF} :: monte_carlo['replicate<r>'].holm_p['{lab}']; derived: "
              f"replicates in which the anatomy passes Holm over the anatomies (of {n_rep})")
        if m["log2_s50_ratio_sd"] is not None:
            sd = m["log2_s50_ratio_sd"]
            F.add(f"cf_{t}_mc_log2_sd", sig(sd, 2), sd, f"{key}.log2_s50_ratio_sd (SD of log2 S50 ratio over the replicates)")
            F.add(f"cf_{t}_mc_sd_pct", pct(2 ** sd - 1), 2 ** sd - 1, f"{key}.log2_s50_ratio_sd; derived: 2**SD - 1 (the "
                  "replicate-to-replicate spread of the ratio as a factor)")
            sds.append(sd)
        allv += vals
        n05 += m["n_p_below_alpha"]
        cells += n_rep
        pass_cells += k
        (every if k == n_rep else never if k == 0 else []).append(lab)
    ALLF = _files(labs)
    if allv:
        F.spread("cf_mc_ratio", allv, ratio, f"{ALLF} :: monte_carlo['{base_name}'].s50_ratio over every anatomy and "
                 f"replicate ({len(allv)} uncensored)")
    if sds:
        F.spread("cf_mc_log2_sd", sds, lambda x: sig(x, 2), f"{ALLF} :: monte_carlo['{base_name}'].log2_s50_ratio_sd over "
                 "the anatomies")
        F.spread("cf_mc_sd_pct", [2 ** x - 1 for x in sds], pct, f"{ALLF} :: monte_carlo['{base_name}'].log2_s50_ratio_sd; "
                 "derived: 2**SD - 1 over the anatomies")
    F.add("cf_mc_n_cells", count(cells), cells, f"{CF} :: monte_carlo_summary.n_replicates; derived: anatomies x replicates")
    F.add("cf_mc_n_p05", count(n05), n05, f"{ALLF} :: monte_carlo['{base_name}'].n_p_below_alpha, summed over the anatomies "
          "(uncorrected)")
    F.add("cf_mc_n_holm_pass_cells", count(pass_cells), pass_cells, f"{CF} :: monte_carlo['replicate<r>'].holm_p; derived: "
          "anatomy-replicate pairs passing Holm over the anatomies")
    F.add("cf_mc_n_anatomies_pass_every_replicate", count(len(every)), len(every), f"{CF} :: monte_carlo['replicate<r>'].holm_p; "
          "derived: anatomies passing in every replicate: " + names(every))
    F.add("cf_mc_anatomies_pass_every_replicate", names(every), every, f"{CF} :: monte_carlo['replicate<r>'].holm_p; derived")
    F.add("cf_mc_n_anatomies_pass_no_replicate", count(len(never)), len(never), f"{CF} :: monte_carlo['replicate<r>'].holm_p; "
          "derived: anatomies passing in no replicate: " + names(never))
    F.add("cf_mc_anatomies_pass_no_replicate", names(never), never, f"{CF} :: monte_carlo['replicate<r>'].holm_p; derived")


# ----------------------------------------------------------------------------------------------
# the declared detector-mismatch variant
def _mismatch(F: Facts, D: dict) -> None:
    per, labs = D["per"], D["labs"]
    first = labs[0]
    P0, ALL = PER.format(first), _files(labs)
    conf_cfg = per[first]["config"]["confirmatory"]
    rt = f"{D['cf']['endpoint']['definition']['false_events_per_min']:g}"
    dv = per[first]["declared_choices"]["variants"]
    inj = dv["primary"]["stretches"]
    F.add("cf_injected_stretches", const(inj), inj, f"{P0} :: declared_choices.variants.primary.stretches (the injected "
          "morphologies' duration scales, the primary detector's templates)")
    for x in conf_cfg.get("variant", []):
        v = x["name"]
        st = dv[v]["stretches"]
        F.add(f"cf_{v}_templates", listing([ratio(s) for s in st]), st, f"{P0} :: declared_choices.variants['{v}'].stretches "
              f"({x['template_rule'].replace('_', ' ')} of the injected stretches)")
        F.add(f"cf_{v}_n_templates", count(dv[v]["n_templates"]), dv[v]["n_templates"],
              f"{P0} :: declared_choices.variants['{v}'].n_templates")
        off = dv[v]["nearest_template_offset"]
        F.spread(f"cf_{v}_template_offset_pct", list(off.values()), pct, f"{P0} :: declared_choices.variants['{v}']."
                 "nearest_template_offset (relative distance of each injected stretch to the nearest template)")
        F.add(f"cf_{v}_bem_layers", count(len(x["dictionary_conductivity"])), len(x["dictionary_conductivity"]),
              f"{P0} :: config.confirmatory.variant['{v}'].dictionary_conductivity (layers of the candidate fields' BEM; the "
              "truth: 3)")
        F.add(f"cf_{v}_coreg_shift_mm", const(x["coreg_shift_mm"]), x["coreg_shift_mm"],
              f"{P0} :: config.confirmatory.variant['{v}'].coreg_shift_mm (one draw per anatomy, shared by all arrays)")
        F.add(f"cf_{v}_coreg_angle_deg", const(x["coreg_angle_deg"]), x["coreg_angle_deg"],
              f"{P0} :: config.confirmatory.variant['{v}'].coreg_angle_deg")
        geo = {lab: per[lab]["mismatch_geometry"][v] for lab in labs if v in per[lab]["mismatch_geometry"]}  # early tests: none
        if D["confirmatory"] and len(geo) != len(labs):
            raise ValueError(f"mismatch_geometry['{v}'] missing in a confirmatory run")
        for lab, g in geo.items():
            F.add(f"cf_{tag(lab)}_{v}_displacement_mm", mm(g["displacement_at_locations_mm_median"]),
                  g["displacement_at_locations_mm_median"], f"{PER.format(lab)} :: mismatch_geometry['{v}']."
                  "displacement_at_locations_mm_median (the analyst's head-to-MRI error at the event locations)")
            F.add(f"cf_{tag(lab)}_{v}_error_rotation_deg", sig(g["error_rotation_deg"], 2), g["error_rotation_deg"],
                  f"{PER.format(lab)} :: mismatch_geometry['{v}'].error_rotation_deg (the drawn coregistration rotation)")
            F.add(f"cf_{tag(lab)}_{v}_error_translation_mm", mm(g["error_translation_mm"]), g["error_translation_mm"],
                  f"{PER.format(lab)} :: mismatch_geometry['{v}'].error_translation_mm (the drawn coregistration shift)")
        if geo:
            GF = _files(list(geo))
            F.spread(f"cf_{v}_displacement_mm", [g["displacement_at_locations_mm_median"] for g in geo.values()], mm,
                     f"{GF} :: mismatch_geometry['{v}'].displacement_at_locations_mm_median over the anatomies")
            F.spread(f"cf_{v}_displacement_locations_mm", [y for g in geo.values() for y in g["displacement_at_locations_mm_range"]],
                     mm, f"{GF} :: mismatch_geometry['{v}'].displacement_at_locations_mm_range over every anatomy's locations")
        for kind, kt in ((f"{v}@{rt}", ""), (f"{v}@{rt}_matched", "_matchedrate")):
            for scope, st in SCOPES.items():
                for key, dt in DETSETS.items():
                    name = f"{v}/cost/{key}/{kind}/{scope}"
                    vals, sides = [], []
                    for lab in labs:
                        sr = per[lab]["mismatch"][name]
                        x_ = F.ratio(f"cf_{tag(lab)}_{v}_cost_{dt}{kt}{st}", sr, f"{PER.format(lab)} :: mismatch['{name}'] "
                                     "(S50 with the mismatched detector over S50 with the primary detector; above 1: the "
                                     "mismatch raises the S50)")
                        vals += [] if x_ is None else [x_]
                        sides.append((lab, _ci_side(sr)))
                    _spread_counts(F, f"cf_{v}_cost_{dt}{kt}{st}", vals, sides, f"{ALL} :: mismatch['{name}']")
                for a, at in ARRAYS.items():
                    name = f"{v}/ratio_change/{a}_vs_{COMPARATOR}/{kind}/{scope}"
                    vals, sides = [], []
                    for lab in labs:
                        sr = per[lab]["mismatch"][name]
                        x_ = F.ratio(f"cf_{tag(lab)}_{v}_ratio_change_{at}{kt}{st}", sr, f"{PER.format(lab)} :: "
                                     f"mismatch['{name}'] (S50 ratio Neuromag / OPM with the mismatched detector over the "
                                     "same with the primary detector = Neuromag's mismatch cost over the OPM's; above 1: the "
                                     "mismatch costs Neuromag more)")
                        vals += [] if x_ is None else [x_]
                        sides.append((lab, _ci_side(sr)))
                    _spread_counts(F, f"cf_{v}_ratio_change_{at}{kt}{st}", vals, sides, f"{ALL} :: mismatch['{name}']")


def _spread_counts(F: Facts, base: str, vals: list, sides: list, src: str) -> None:
    if vals:
        F.spread(base, vals, ratio, f"{src}.value over the anatomies ({len(vals)} uncensored)")
    for side in ("above", "below", "includes"):
        sel = [lab for lab, s in sides if s == side]
        F.add(f"{base}_n_ci_{side}1", count(len(sel)), len(sel), f"{src}.ci95; derived: anatomies with the interval "
              + {"above": "above 1", "below": "below 1", "includes": "including 1"}[side] + ": " + names(sel))


# ----------------------------------------------------------------------------------------------
# realized false-event rates, equal-rate tests, the oracle's false-positive probability
def _false_events(F: Facts, D: dict) -> None:
    per, labs = D["per"], D["labs"]
    ALL = _files(labs)
    target = D["cf"]["endpoint"]["definition"]["false_events_per_min"]
    variants = {"primary": "practical"}
    variants.update({x["name"]: x["name"] for x in per[labs[0]]["config"]["confirmatory"].get("variant", [])})
    groups: dict = {}
    for lab in labs:
        fe = per[lab]["false_events"]
        for key, dt in DETSETS.items():
            for v, vt in variants.items():
                e = fe[f"{key}|{v}"]
                for kind, what in RATE_KINDS.items():
                    x = e[kind]
                    nm = f"cf_{tag(lab)}_rate_{dt}_{vt}_{kind}"
                    src = f"{PER.format(lab)} :: false_events['{key}|{v}'].{kind}"
                    F.add(nm, share(x["rate_per_min"]), x["rate_per_min"], f"{src}.rate_per_min (false events per minute, "
                          f"{what}: {x['count']} in {x['minutes']:g} min)")
                    F.add(f"{nm}_ci", iv(x["ci95"][0], x["ci95"][1], share), x["ci95"], f"{src}.ci95 (exact Poisson)")
                    groups.setdefault((vt, kind, dt), []).append((lab, x))
    for (vt, kind, dt), xs in groups.items():
        F.spread(f"cf_rate_{dt}_{vt}_{kind}", [x["rate_per_min"] for _, x in xs], share,
                 f"{ALL} :: false_events['{_key_of(dt)}|{_variant_of(vt)}'].{kind}.rate_per_min over the anatomies")
    for vt in variants.values():
        for kind in RATE_KINDS:
            xs = [x for dt in DETSETS.values() for _, x in groups[(vt, kind, dt)]]
            src = f"{ALL} :: false_events['<array>|{_variant_of(vt)}'].{kind}"
            F.spread(f"cf_rate_{vt}_{kind}", [x["rate_per_min"] for x in xs], share, f"{src}.rate_per_min over the "
                     "anatomies and the three arrays")
            out = sum(not (x["ci95"][0] <= target <= x["ci95"][1]) for x in xs)
            F.add(f"cf_rate_{vt}_{kind}_n_ci_excludes_target", count(out), out, f"{src}.ci95; derived: anatomy-array pairs "
                  f"whose exact Poisson interval excludes the target {const(target)} per minute (of {len(xs)})")
            ep = [x for dt in ("combined", "dense") for _, x in groups[(vt, kind, dt)]]
            F.spread(f"cf_rate_endpoint_arrays_{vt}_{kind}", [x["rate_per_min"] for x in ep], share,
                     f"{ALL} :: false_events['squid/combined|{_variant_of(vt)}' and 'opm_dense/opm|{_variant_of(vt)}']."
                     f"{kind}.rate_per_min over the anatomies (the endpoint's two arrays)")
    # equal-rate tests between each OPM array and Neuromag on the evaluation null (approximate: shared noise)
    ps = []
    for lab in labs:
        for name, t in per[lab]["false_event_rate_equality_on_evaluation_null"].items():
            a, v, kind = name.split("/")[0] + "/opm", name.split("/")[-2], name.split("/")[-1]
            pkey = "conditional_binomial_p" if "conditional_binomial_p" in t else "exact_conditional_p"
            p = t[pkey]
            F.add(f"cf_{tag(lab)}_rate_equality_{ARRAYS[a]}_{variants[v]}_{kind}_p", pval(p), p,
                  f"{PER.format(lab)} :: false_event_rate_equality_on_evaluation_null['{name}'].{pkey} (equal false-event "
                  f"rates of the OPM array and Neuromag: {t['count_opm']} vs {t['count_squid']} events; approximate, the "
                  "arrays share the background and room noise)")
            ps.append((lab, name, p))
    src = f"{ALL} :: false_event_rate_equality_on_evaluation_null[*].conditional_binomial_p"
    F.add("cf_rate_equality_n", count(len(ps)), len(ps), f"{src} (tests: anatomies x OPM arrays x detectors x thresholds)")
    low = [f"{lab} {n}" for lab, n, p in ps if p < 0.05]
    F.add("cf_rate_equality_n_p05", count(len(low)), len(low), f"{src}; derived: tests with p < 0.05 (uncorrected): "
          + listing(low))
    F.spread("cf_rate_equality_p", [p for *_, p in ps], pval, src)
    # the oracle on the held-out null
    fpp = []
    for lab in labs:
        for key, dt in DETSETS.items():
            o = per[lab]["oracle"][key]
            if o["heldout_false_positive_probability"] is None:
                continue
            src = f"{PER.format(lab)} :: oracle['{key}']"
            F.add(f"cf_{tag(lab)}_oracle_fpp_{dt}", prob(o["heldout_false_positive_probability"]),
                  o["heldout_false_positive_probability"], f"{src}.heldout_false_positive_probability (at the frozen z_crit, "
                  f"{o['heldout_exceed']} of {o['heldout_samples']:,} held-out samples)")
            F.add(f"cf_{tag(lab)}_oracle_fpp_{dt}_ci", iv(o["ci95"][0], o["ci95"][1], prob), o["ci95"],
                  f"{src}.ci95 (exact binomial)")
            fpp.append(o["heldout_false_positive_probability"])
    if fpp:
        F.spread("cf_oracle_fpp", fpp, prob, f"{ALL} :: oracle[*].heldout_false_positive_probability over the anatomies and "
                 "arrays")


def _key_of(dt: str) -> str:
    return next(k for k, v in DETSETS.items() if v == dt)


def _variant_of(vt: str) -> str:
    return "primary" if vt == "practical" else vt


# ----------------------------------------------------------------------------------------------
# locations, dictionaries; the exploratory run beside the confirmatory one
def _locations(F: Facts, D: dict, conf: bool) -> None:
    per, labs = D["per"], D["labs"]
    ALL = _files(labs)
    for lab in labs:
        t, P = tag(lab), PER.format(lab)
        lc = per[lab]["location_checks"]
        F.add(f"cf_{t}_median_depth_mm", mm(lc["depth_mm_median"]), lc["depth_mm_median"],
              f"{P} :: location_checks.depth_mm_median (mm below the scalp)")
        F.add(f"cf_{t}_nearest_exploratory_mm", mm(lc["nearest_exploratory_location_mm_min"]),
              lc["nearest_exploratory_location_mm_min"], f"{P} :: location_checks.nearest_exploratory_location_mm_min (the "
              "smallest distance of a confirmatory location to an exploratory one)")
        F.add(f"cf_{t}_nearest_exploratory_median_mm", mm(lc["nearest_exploratory_location_mm_median"]),
              lc["nearest_exploratory_location_mm_median"], f"{P} :: location_checks.nearest_exploratory_location_mm_median "
              "(median over the confirmatory locations of the distance to the nearest exploratory one)")
        F.add(f"cf_{t}_n_dictionary", count(per[lab]["n_dictionary"]), per[lab]["n_dictionary"],
              f"{P} :: n_dictionary (candidate sources of the scanning detector)")
        lobes = {}
        for L in per[lab]["locations"]:
            lobes[L["lobe"]] = lobes.get(L["lobe"], 0) + 1
        for lobe in LOBES:
            F.add(f"cf_{t}_locations_{lobe}", count(lobes.get(lobe, 0)), lobes.get(lobe, 0),
                  f"{P} :: derived: count of locations[].lobe == '{lobe}'")
    lc = {lab: per[lab]["location_checks"] for lab in labs}
    F.spread("cf_median_depth_mm", [x["depth_mm_median"] for x in lc.values()], mm,
             f"{ALL} :: location_checks.depth_mm_median over the anatomies")
    near = min(x["nearest_exploratory_location_mm_min"] for x in lc.values())
    F.add("cf_nearest_exploratory_mm", mm(near), near, f"{ALL} :: location_checks.nearest_exploratory_location_mm_min, "
          "smallest over the anatomies (no confirmatory location lies closer to an exploratory one)")
    F.spread("cf_nearest_exploratory_median_mm", [x["nearest_exploratory_location_mm_median"] for x in lc.values()], mm,
             f"{ALL} :: location_checks.nearest_exploratory_location_mm_median over the anatomies")
    F.spread("cf_n_dictionary", [per[lab]["n_dictionary"] for lab in labs], count, f"{ALL} :: n_dictionary over the anatomies")
    disjoint = all(x["disjoint_from_exploratory"] for x in lc.values())
    F.add("cf_locations_disjoint", yes(disjoint), disjoint, f"{ALL} :: location_checks.disjoint_from_exploratory (no "
          "confirmatory location on an exploratory vertex, every anatomy)")
    pools = [(lab, p) for lab in labs for p in per[lab]["location_pools"]]
    nd = [f"{lab} {p['orientation_deg'][0]:g}-{p['orientation_deg'][1]:g} deg" for lab, p in pools if not p["disjoint_from_exploratory"]]
    F.add("cf_n_strata_not_disjoint", count(len(nd)), len(nd), f"{ALL} :: location_pools[*].disjoint_from_exploratory (strata "
          "too small to leave the exploratory vertices out: " + listing(nd) + ")")
    F.spread("cf_stratum_pool", [p["pool_without_exploratory"] for _, p in pools], count,
             f"{ALL} :: location_pools[*].pool_without_exploratory (candidate cortical vertices per orientation stratum of the "
             "band, exploratory vertices left out)")
    _common(F, "cf_n_exploratory_vertices_excluded", {lab: x["n_exploratory_vertices_excluded"] for lab, x in lc.items()},
            count, f"{ALL} :: location_checks.n_exploratory_vertices_excluded", conf)


def _exploratory(F: Facts, D: dict) -> None:
    cf, per, labs = D["cf"], D["per"], D["labs"]
    ep = cf["endpoint"]["comparison"]
    xkey = "{}/paired/" + ep.rsplit("/", 2)[0] + "/" + ep.rsplit("/", 2)[1] + "/depth0"
    conf_sr, explo_sr = {}, {}
    for lab in labs:
        e, c = cf["anatomy"][lab]["exploratory"], D["explo"][xkey.format(lab)]
        for k in ("s50_ratio_squid_over_opm", "location_sign_flip_p", "locations_favouring_opm", "locations_favouring_squid"):
            if e[k] != c[k]:
                raise ValueError(f"{CF} :: anatomy['{lab}'].exploratory.{k} differs from {EXPLO} :: comparison"
                                 f"['{xkey.format(lab)}']")
        conf_sr[lab], explo_sr[lab] = per[lab]["comparisons"][ep]["s50_ratio_squid_over_opm"], c["s50_ratio_squid_over_opm"]
    src = (f"{_files(labs)}, {EXPLO} :: comparisons['{ep}'].s50_ratio_squid_over_opm (confirmatory) against "
           f"comparison['<anatomy>/paired/.../depth0'].s50_ratio_squid_over_opm (exploratory, endpoint chosen after the analyses)")
    both = [lab for lab in labs if conf_sr[lab]["value"] is not None and explo_sr[lab]["value"] is not None]
    rr = {lab: conf_sr[lab]["value"] / explo_sr[lab]["value"] for lab in both}
    at, mt = ARRAYS[ep.split("_vs_")[0]], _modes(D)[ep.rsplit("/", 2)[1]]
    for lab in both:
        F.add(f"cf_{tag(lab)}_{at}_vs_combined_{mt}_over_exploratory", ratio(rr[lab]), rr[lab], f"{src}; derived: "
              f"confirmatory / exploratory ratio, {NAMES[lab]}")
    if rr:
        F.spread(f"cf_all_{at}_vs_combined_{mt}_over_exploratory", list(rr.values()), ratio, f"{src}; derived: confirmatory / "
                 "exploratory ratio over the anatomies with both uncensored")
    below = [lab for lab in both if rr[lab] < 1]
    F.add(f"cf_{at}_vs_combined_{mt}_n_below_exploratory", count(len(below)), len(below), f"{src}; derived: anatomies whose "
          f"confirmatory ratio is below the exploratory one (of {len(both)} with both uncensored): " + names(below))

    def inside(v, ci):
        lo, hi = ci or [None, None]
        return (lo is None or lo <= v) and (hi is None or v <= hi)

    ins = [lab for lab in labs if explo_sr[lab]["value"] is not None and inside(explo_sr[lab]["value"], conf_sr[lab]["ci95"])]
    F.add(f"cf_{at}_vs_combined_{mt}_n_exploratory_within_interval", count(len(ins)), len(ins), f"{src}; derived: anatomies whose "
          "exploratory point estimate lies inside the confirmatory 95 % interval (open ends unbounded): " + names(ins))
    explo_n = {D["explo"][f"{lab}/n_locations_per_depth_band"]["depth0"] for lab in labs}
    if len(explo_n) == 1:
        n = explo_n.pop()
        F.add("cf_exploratory_n_locations", count(n), n, f"{EXPLO} :: comparison['<anatomy>/n_locations_per_depth_band']"
              ".depth0 (the exploratory run's locations in the band, the same in every anatomy)")


# ----------------------------------------------------------------------------------------------
def _exact(D: dict, root: Path) -> dict:
    """Replace every stored sign-flip p, and the Holm values and pass counts derived from it, by the exact p of {EXACT}
    (scripts/g4_confirm_exact_p.py: the declared exact test, computed from the stored location differences), after
    checking that the file belongs to these summaries (its Monte Carlo values are the stored ones, cell by cell) and
    that each copy is matched to its comparison by its differences or ratios. Returns the file's agreement block."""
    if not (root / EXACT).is_file():
        raise FileNotFoundError(f"{EXACT} not found in {root}: run scripts/g4_confirm_exact_p.py")
    ex = _json(root / EXACT)
    comp, cf, per, labs = ex["comparisons"], D["cf"], D["per"], D["labs"]
    alpha = cf["endpoint"]["definition"]["alpha"]

    def exact(path: str, lab: str, stored: float) -> float:
        e = comp[path][lab]
        if e["p_monte_carlo"] != stored:
            raise ValueError(f"{EXACT}: {path} {lab} was computed from p {e['p_monte_carlo']}, the summary stores {stored} "
                             "(rerun scripts/g4_confirm_exact_p.py)")
        return e["p_exact"]

    ep = cf["endpoint"]["comparison"]
    for lab in labs:
        rows = per[lab]["comparisons"]
        for path, c in rows.items():
            if "location_differences" in c:
                c["location_sign_flip_p"] = exact(path, lab, c["location_sign_flip_p"])
        for base, m in per[lab]["monte_carlo"].items():
            m["location_sign_flip_p"] = [exact(f"{base}/replicate{r}", lab, v) for r, v in enumerate(m["location_sign_flip_p"])]
            m["n_p_below_alpha"] = sum(v < alpha for v in m["location_sign_flip_p"])
        a = cf["anatomy"][lab]
        a["endpoint"]["location_sign_flip_p"] = exact(ep, lab, a["endpoint"]["location_sign_flip_p"])
        a["endpoint"]["holm_p"] = ex["families"]["endpoint"]["holm_p"][lab]
        for sect, path in SECTIONS.items():
            if (a[sect]["location_differences"] != rows[path]["location_differences"]
                    or a[sect]["s50_ratio_squid_over_opm"] != rows[path]["s50_ratio_squid_over_opm"]):
                raise ValueError(f"{CF}: anatomy['{lab}'].{sect} is not {PER.format(lab)} comparisons['{path}']")
            a[sect]["location_sign_flip_p"] = exact(path, lab, a[sect]["location_sign_flip_p"])
        base = ep.rsplit("/", 1)[0]
        if a["monte_carlo"]["s50_ratio"] != per[lab]["monte_carlo"][base]["s50_ratio"]:
            raise ValueError(f"{CF}: anatomy['{lab}'].monte_carlo is not {PER.format(lab)} monte_carlo['{base}']")
        a["monte_carlo"]["location_sign_flip_p"] = [exact(f"{base}/replicate{r}", lab, v)
                                                    for r, v in enumerate(a["monte_carlo"]["location_sign_flip_p"])]
    for group, key in [("families", k) for k in cf["families"]] + [("monte_carlo", k) for k in ex["monte_carlo"]]:
        fam, e = cf[group][key], ex[group][key]
        if e["holm_p_monte_carlo"] != fam["holm_p"] or e["n_pass_monte_carlo"] != fam["n_pass"]:
            raise ValueError(f"{EXACT}: {group}['{key}'] does not belong to {CF} (rerun scripts/g4_confirm_exact_p.py)")
        fam.update(p=dict(e["p"]), holm_p=dict(e["holm_p"]), n_pass=e["n_pass"])
    mcs = cf["monte_carlo_summary"]
    mcs["n_pass_per_replicate"] = [ex["monte_carlo"][f"replicate{r}"]["n_pass"] for r in range(mcs["n_replicates"])]
    mcs["n_replicates_all_pass"] = sum(n == len(labs) for n in mcs["n_pass_per_replicate"])
    return ex["agreement"]


def _exact_facts(F: Facts, ag: dict) -> None:
    """How the exact sign-flip p compares with the run's Monte Carlo estimates."""
    s = f"{EXACT} :: agreement"
    F.add("cf_exact_n_cells", count(ag["n_cells"]), ag["n_cells"], f"{s}.n_cells (anatomy x comparison)")
    F.add("cf_exact_n_cells_monte_carlo", count(ag["n_cells_monte_carlo"]), ag["n_cells_monte_carlo"],
          f"{s}.n_cells_monte_carlo (more than 20 non-zero location differences: the run sampled 20,000 sign patterns)")
    F.add("cf_exact_n_cells_enumerated", count(ag["n_cells_enumerated_by_the_run"]), ag["n_cells_enumerated_by_the_run"],
          f"{s}.n_cells_enumerated_by_the_run (the run enumerated every pattern itself; reproduced exactly)")
    F.add("cf_exact_max_abs_p_diff", sig(ag["max_abs_p_difference"], 2), ag["max_abs_p_difference"],
          f"{s}.max_abs_p_difference (largest |exact p - Monte Carlo p| over the cells)")
    F.add("cf_exact_n_families", count(ag["n_families"]), ag["n_families"], f"{s}.n_families (the combined summary's "
          "families and noise replicates of the endpoint)")
    F.add("cf_exact_n_families_same_decisions", count(ag["n_families_same_decisions"]), ag["n_families_same_decisions"],
          f"{s}.n_families_same_decisions (every anatomy's Holm decision at alpha the same with exact and Monte Carlo p)")
    F.add("cf_exact_max_p_where_mc_zero", pval(ag["max_exact_p_where_monte_carlo_zero"]),
          ag["max_exact_p_where_monte_carlo_zero"], f"{s}.max_exact_p_where_monte_carlo_zero")
    F.add("cf_exact_max_holm_p_where_mc_zero", sig(ag["max_exact_holm_p_where_monte_carlo_holm_zero"], 2),
          ag["max_exact_holm_p_where_monte_carlo_holm_zero"], f"{s}.max_exact_holm_p_where_monte_carlo_holm_zero")


# ----------------------------------------------------------------------------------------------
def facts(root: Path = ROOT) -> dict:
    """Every fact of the confirmatory spike run: name -> {"value", "raw", "source"}."""
    root = Path(root)
    D = _read(root)
    reasons = _status(D)
    _consistency(D)
    agreement = _exact(D, root)  # from here on every sign-flip p is the exact one
    D["confirmatory"] = not reasons
    F = Facts("" if not reasons else " [run not confirmatory: " + "; ".join(reasons) + "]")
    _design(F, D, reasons)
    _families(F, D)
    _monte_carlo(F, D)
    _mismatch(F, D)
    _false_events(F, D)
    _locations(F, D, not reasons)
    _exploratory(F, D)
    _exact_facts(F, agreement)
    for name, f in F.items():  # values derived from the confirmatory run's sign-flip p are now the exact ones
        if name.startswith("cf_") and not name.startswith("cf_exact_") and "exploratory" not in name and any(
                k in f["source"] for k in ("sign_flip_p", "holm_p", "n_pass", "n_p_below_alpha", "n_replicates_all_pass")):
            f["source"] += EXACT_NOTE
    return dict(F)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=ROOT, help="repository root holding results/g4_confirm (default: this one)")
    ap.add_argument("--grep", default="", help="print only facts whose name contains this text")
    args = ap.parse_args(argv)
    f = facts(args.root)
    for name, x in f.items():
        if args.grep in name:
            print(f"{name}\t{x['value']}\t{x['source']}")
    print(f"{len(f)} facts")


if __name__ == "__main__":
    main()
