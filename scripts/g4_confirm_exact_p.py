#!/usr/bin/env python3
"""Exact p values of the confirmatory spike run's location sign-flip tests (derived from stored results).

The declared test of the confirmatory run (configs/g4_confirmatory.toml, [endpoint].test) is the two-sided exact
sign-flip test on the per-location differences in detection counts. Its implementation, opmsquid.detection.sign_flip_p,
enumerates every sign pattern up to 20 non-zero differences and draws 20,000 random patterns beyond; with 36 locations
per anatomy most cells were therefore Monte Carlo estimates, printed as zero where no sampled pattern was as extreme.

The differences are integers, so the exact p follows from the distribution of the signed sum S = sum_i s_i d_i over all
2^n equiprobable sign patterns, which is a convolution of n two-point distributions (exact integer counts):
p = P(|S| >= |sum_i d_i|), the same statistic and two-sided rule as sign_flip_p. This script computes that p for every
stored comparison of every anatomy (results/g4_confirm/g4c_<anatomy>_summary.json, comparisons[*].location_differences),
the Holm adjustment over the anatomies of every family of the combined summary (results/g4_confirm/
g4_confirm_summary.json: families and monte_carlo), and the agreement with the stored Monte Carlo values (p and the
decisions at the declared alpha). Nothing is simulated again.

Output: results/g4_confirm/g4_confirm_exact_p.json
Usage: .venv/bin/python scripts/g4_confirm_exact_p.py
"""
from __future__ import annotations

import json
import sys
import tomllib
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from opmsquid import detection, io  # noqa: E402

DIR = ROOT / "results" / "g4_confirm"
CF = DIR / "g4_confirm_summary.json"
OUT = DIR / "g4_confirm_exact_p.json"
CFG = ROOT / "configs" / "g4_confirmatory.toml"
ENDPOINT_PATH = "opm_dense/opm_vs_squid/combined/practical@1/replicate0"


def exact_sign_flip(d) -> tuple[Fraction, int]:
    """Exact two-sided sign-flip p of integer differences ``d`` (zeros dropped, as sign_flip_p does), as a Fraction,
    and the number of non-zero differences."""
    x = [int(v) for v in d if v != 0]
    if any(float(v) != float(w) for v, w in zip([v for v in d if v != 0], x)):
        raise ValueError("location differences must be integers")
    if not x:
        return Fraction(1), 0
    tot = sum(abs(v) for v in x)
    counts = [0] * (2 * tot + 1)  # counts[tot + s]: number of sign patterns with signed sum s
    counts[tot] = 1
    for v in x:
        a = abs(v)
        new = [0] * len(counts)
        for i, c in enumerate(counts):
            if c:
                new[i + a] += c
                new[i - a] += c
        counts = new
    obs = abs(sum(x))
    hits = sum(c for i, c in enumerate(counts) if abs(i - tot) >= obs)
    return Fraction(hits, 2 ** len(x)), len(x)


def holm(ps: dict) -> dict:
    """Holm step-down adjustment (the same rule as the combine step's and the report facts')."""
    run, out = 0.0, {}
    items = sorted(ps.items(), key=lambda kv: kv[1])
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (len(items) - i) * p))
        out[k] = run
    return out


def main():
    cf = json.loads(CF.read_text())
    cfg = tomllib.loads(CFG.read_text())
    alpha = float(cfg["endpoint"]["alpha"])
    labs = list(cf["anatomies"])
    per = {lab: json.loads((DIR / f"g4c_{lab}_summary.json").read_text()) for lab in labs}

    comparisons, checks = {}, dict(enumerated_cells_reproduced=0, max_abs_p_difference=0.0)
    for lab in labs:
        for path, c in per[lab]["comparisons"].items():
            if "location_differences" not in c:
                continue
            p, n = exact_sign_flip(c["location_differences"])
            mc = float(c["location_sign_flip_p"])
            if n <= 20:  # sign_flip_p enumerated these patterns itself: its value must be reproduced exactly
                ref = detection.sign_flip_p(c["location_differences"])
                if abs(float(p) - ref) > 1e-12 or abs(mc - ref) > 1e-12:
                    raise ValueError(f"{lab} {path}: exact p {float(p)} differs from the enumerated {ref} (stored {mc})")
                checks["enumerated_cells_reproduced"] += 1
            checks["max_abs_p_difference"] = max(checks["max_abs_p_difference"], abs(float(p) - mc))
            comparisons.setdefault(path, {})[lab] = dict(p_exact=float(p), p_exact_fraction=f"{p.numerator}/{p.denominator}",
                                                         n_nonzero=n, p_monte_carlo=mc,
                                                         monte_carlo=n > 20)

    def family(path: str, stored: dict) -> dict:
        ps = {lab: comparisons[path][lab]["p_exact"] for lab in stored["p"]}
        for lab, mc in stored["p"].items():  # the family's stored p must be the comparison's stored p
            if abs(mc - comparisons[path][lab]["p_monte_carlo"]) > 1e-15:
                raise ValueError(f"{path} {lab}: the family's stored p {mc} is not the comparison's")
        hp = holm(ps)
        n_pass = sum(v < alpha for v in hp.values())
        same = all((hp[lab] < alpha) == (stored["holm_p"][lab] < alpha) for lab in hp)
        return dict(comparison=path, p=ps, holm_p=hp, n_pass=n_pass, n=len(hp),
                    holm_p_monte_carlo=stored["holm_p"], n_pass_monte_carlo=stored["n_pass"], same_decisions=same)

    families = {}
    for key, fam in cf["families"].items():
        path = ENDPOINT_PATH if key == "endpoint" else key.removeprefix("secondary/")
        families[key] = family(path, fam)
    monte_carlo = {}
    for key, fam in cf["monte_carlo"].items():
        if not key.startswith("replicate"):
            continue
        monte_carlo[key] = family(ENDPOINT_PATH.replace("replicate0", key), fam)

    every = list(families.values()) + list(monte_carlo.values())
    cells = [v for path in comparisons.values() for v in path.values()]
    mc_zero = [v for v in cells if v["p_monte_carlo"] == 0.0]
    adj_mc_zero = [f["holm_p"][lab] for f in every for lab in f["holm_p"] if f["holm_p_monte_carlo"][lab] == 0.0]
    out = dict(
        status="DERIVED (exact p of the confirmatory run's sign-flip tests, from the stored location differences; "
               "nothing simulated again)",
        method="two-sided sign-flip test on the per-location differences in detection counts (zeros dropped): "
               "p = P(|sum s_i d_i| >= |sum d_i|) over all 2^n equiprobable sign patterns, the distribution of the "
               "signed sum computed exactly by convolution (integer counts); the statistic and two-sided rule of "
               "opmsquid.detection.sign_flip_p, which enumerates up to 20 non-zero differences and draws 20,000 random "
               "patterns beyond. Holm over the anatomies of each family, as the combine step.",
        alpha=alpha,
        comparisons=comparisons,
        families=families,
        monte_carlo=monte_carlo,
        agreement=dict(
            n_cells=len(cells),
            n_cells_monte_carlo=sum(v["monte_carlo"] for v in cells),
            n_cells_enumerated_by_the_run=checks["enumerated_cells_reproduced"],
            max_abs_p_difference=checks["max_abs_p_difference"],
            n_cells_monte_carlo_zero=len(mc_zero),
            max_exact_p_where_monte_carlo_zero=max((v["p_exact"] for v in mc_zero), default=None),
            max_exact_holm_p_where_monte_carlo_holm_zero=max(adj_mc_zero, default=None),
            n_families=len(every),
            n_families_same_decisions=sum(f["same_decisions"] for f in every),
            endpoint_n_pass_exact=families["endpoint"]["n_pass"],
            endpoint_max_holm_p_exact=max(families["endpoint"]["holm_p"].values()),
        ),
    )
    if out["agreement"]["n_families_same_decisions"] != out["agreement"]["n_families"]:
        print("note: the exact and Monte Carlo Holm decisions differ in at least one family")
    io.write_json(out, OUT)
    a = out["agreement"]
    print(f"{a['n_cells']} cells ({a['n_cells_monte_carlo']} Monte Carlo, {a['n_cells_enumerated_by_the_run']} enumerated by "
          f"the run and reproduced); largest |exact - Monte Carlo| {a['max_abs_p_difference']:.4g}; families with the same "
          f"decisions: {a['n_families_same_decisions']} of {a['n_families']}; endpoint {a['endpoint_n_pass_exact']} of "
          f"{len(families['endpoint']['holm_p'])}, largest Holm p {a['endpoint_max_holm_p_exact']:.3g} -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    np.seterr(all="raise")
    main()
