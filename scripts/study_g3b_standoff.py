#!/usr/bin/env python3
"""G3B sensitivity: the adult-only OPM standoff (end-to-end goal review, 2026-10-02).

The v3 clearance rule (A-OPM-CLEAR) keeps the whole OPM cell outside each anatomy's BEM head
surface. The adult's BEM head surface lies about 1 mm outside its MRI scalp, the templates' slightly
inside theirs, so the rule moved 164 of the adult's 205 dense sites outward (median sensing-centre
height above the scalp 7.78 mm) and almost none of the children's (7.00 mm). D_adult therefore
carries a larger standoff than D_child, and every Delta = D_child - D_adult inherits the difference.

This check bounds the effect without new assumptions: each child's dense and matched OPM arrays are
moved outward along their axes by the adult's median excess standoff (same sites, as the G2
scalp-gap variants), D_child is recomputed with the stored Neuromag results (primary placement and
the counterfactual helmets; the stored G3B state, nothing else re-simulated) and Delta is re-estimated
with the same parcel/vertex bootstrap as G3B. The shift is uniform, the adult's is site-specific
(0.5-5 mm at the moved sites), so this is a first-order bound.

Output: results/g3b/g3b_standoff_sensitivity.json, results/g3b/G3B_standoff_report.md.
"""
from __future__ import annotations

import pickle
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import mne  # noqa: E402
import numpy as np  # noqa: E402

import g3b_pediatric_helmet as G3  # noqa: E402
from opmsquid import g2, io  # noqa: E402

OUT = ROOT / "results" / "g3b"
COMPARATORS = ("squid:top", "squid:counterfactual", "squid:counterfactual_x-centred")
N_BOOT = 500


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    g2cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    with open(G3.STATE, "rb") as fh:
        state = pickle.load(fh)
    runs = state["runs"]
    anats = G3.load_anatomies(cfg, g2cfg)
    com = G3.Common(cfg, g2cfg)
    com.brain_scale = state["brain_scale"]  # the adult calibration of the stored run
    cond = cfg["condition"] if "condition" in cfg else "intrinsic+brain"
    adult_h = {n: runs["adult"]["geometry"][n]["median_mm"] for n in ("opm_dense", "opm_matched") if n in runs["adult"]["geometry"]}
    child_h = {c: {n: runs[c]["geometry"][n]["median_mm"] for n in ("opm_dense", "opm_matched") if n in runs[c]["geometry"]}
               for c in G3.CHILDREN}
    shift = {n: float(adult_h[n] - np.median([child_h[c][n] for c in G3.CHILDREN])) for n in adult_h}
    G3.log(f"adult median sensing-centre height {adult_h}, children {child_h}; shift applied {shift} mm")
    bem3 = {}
    res = dict(status="NEW (G3B sensitivity: the children's OPM arrays moved outward by the adult's median excess standoff; "
                      "stored Neuromag results; same bootstrap as G3B)", shift_mm=shift, adult_height_mm=adult_h,
               child_height_mm=child_h, comparators=list(COMPARATORS), n_boot=N_BOOT, anatomies={})
    rng = np.random.default_rng(g2cfg["sources"]["seed"] + 7)
    d_adult = {}
    for sq in COMPARATORS:
        for n in ("opm_dense", "opm_matched"):
            d_adult[(sq, n)] = np.where(anats["adult"].cortical, G3.d_db(runs["adult"]["res"], n, sq, "combined", cond), np.nan)
    for c in G3.CHILDREN:
        an = anats[c]
        arrays, _ = G3.build_arrays(an, com, cfg)
        bem3[c] = an.subject.bem_model(com.bem)
        out = {}
        for n in ("opm_dense", "opm_matched"):
            moved = g2.with_scalp_gap(arrays[n], shift[n] * 1e-3)
            g = G3.gains(an, moved, bem3[c])
            r_moved = G3.evaluate(an, moved, g[:, :an.nt], g[:, an.nt:], com, [cond])
            for sq in COMPARATORS:
                both = {n: r_moved, sq: runs[c]["res"][sq]}
                x_shift = np.where(an.cortical, G3.d_db(both, n, sq, "combined", cond), np.nan)
                x_orig = np.where(an.cortical, G3.d_db(runs[c]["res"], n, sq, "combined", cond), np.nan)
                cmp_o = G3.compare(an, anats["adult"], x_orig, d_adult[(sq, n)], cfg, rng, N_BOOT)
                cmp_s = G3.compare(an, anats["adult"], x_shift, d_adult[(sq, n)], cfg, rng, N_BOOT)
                out[f"{n}/{sq}"] = dict(delta_original=cmp_o["delta"], delta_shifted=cmp_s["delta"],
                                        d_child_original=cmp_o["d_child"], d_child_shifted=cmp_s["d_child"], d_adult=cmp_o["d_adult"])
            G3.log(f"{c}: {n} done ({time.time() - t0:.0f} s)")
        res["anatomies"][c] = out
    io.write_json(res, OUT / "g3b_standoff_sensitivity.json")
    (OUT / "G3B_standoff_report.md").write_text(report(res))
    G3.log(f"done in {time.time() - t0:.0f} s")


def med(x):
    return x["median"] if isinstance(x, dict) else x


def report(res: dict) -> str:
    L = ["# G3B sensitivity: the adult-only OPM standoff (review check)", "",
         f"The children's OPM arrays moved outward along their axes by the adult's median excess standoff "
         f"({res['shift_mm']['opm_dense']:.2f} mm dense, {res['shift_mm']['opm_matched']:.2f} mm matched); Neuromag results as "
         f"stored; Delta re-estimated with the G3B bootstrap ({res['n_boot']} resamples). D in dB (OPM minus Neuromag combined, "
         "intrinsic + brain noise); Delta = D_child - D_adult.", "",
         "| anatomy | array | comparator | D_child original | D_child shifted | D_adult | Delta original | Delta shifted | change |",
         "|---|---|---|---|---|---|---|---|---|"]
    for c, out in res["anatomies"].items():
        for key, v in out.items():
            n, sq = key.split("/", 1)
            do, ds = med(v["delta_original"]), med(v["delta_shifted"])
            L.append(f"| {c} | {n} | {sq.split(':')[1]} | {med(v['d_child_original']):+.2f} | {med(v['d_child_shifted']):+.2f} | "
                     f"{med(v['d_adult']):+.2f} | {do:+.2f} | {ds:+.2f} | {ds - do:+.2f} |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
