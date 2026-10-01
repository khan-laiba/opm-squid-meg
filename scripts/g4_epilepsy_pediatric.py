#!/usr/bin/env python3
"""G4 (pediatric): IED detection and bounded localization on a child anatomy in the fixed adult
helmet, with the adult framework unchanged (NEW).

The anatomy, arrays and conventions are those of G3B (scripts/g3b_pediatric_helmet.py): the
Neuromag helmet at the primary source-blind placement (top contact), the refitted dense OPM array
and the matched-site OPM array; the adult's background moment variance per unit cortical area, room
field and intrinsic noise; the 3-layer BEM for the truth. The detection study
(``g4_epilepsy_adult.simulate``) and the localization study (``g4_localization.localize``) run
unchanged with the configuration of the adult (configs/g4_epilepsy.toml: same seeds, strengths,
morphologies, detectors, operating points, null durations, inverse settings and coregistration
error); thresholds are calibrated on this anatomy's own null data. Event locations are stratified
by the same depth and orientation bands; strata the anatomy cannot fill are reported as such.

Usage: g4_epilepsy_pediatric.py [infant2yr|school] [--detection] [--localization] [--compare]
(default: infant2yr, all three). Outputs: results/g4/g4_<label>_*, g4_localization_<label>_*,
g4_pediatric_comparison.json and G4_pediatric_report.md.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g3b_pediatric_helmet as G3  # noqa: E402
import g4_epilepsy_adult as G4  # noqa: E402
import g4_localization as G4L  # noqa: E402
from opmsquid import background, forward, goldenholz, io, opm  # noqa: E402

OUT = ROOT / "results" / "g4"
DETECTORS = ("squid/combined", "squid/grad", "squid/mag", "opm_matched/opm", "opm_dense/opm")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def child_context(key: str, names=("squid", "opm_matched", "opm_dense")) -> G4.Context:
    """A G3B anatomy with its primary-placement Neuromag helmet and refitted OPM arrays."""
    cfg3 = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    cfg2 = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    anats = G3.load_anatomies(cfg3, cfg2)
    com = G3.Common(cfg3, cfg2)
    adult = anats["adult"]
    bem_a = adult.subject.bem_model(com.bem)
    squid_a = G3.g2.Array("squid", com.squid_info, com.kinds, None, {})
    g = G3.gains(adult, squid_a, bem_a)  # the one calibration: adult, measured head position (as G2 and G3B)
    unit = background.sensor_covariance(g[:, adult.nt:], background.moment_covariance(adult.src.grid_area))
    com.brain_scale = background.calibrate(unit, com.grads, com.target_var)
    an = anats[key]
    arrays_all, _ = G3.build_arrays(an, com, cfg3)
    primary = f"squid:{cfg3['placement']['primary']}"
    arrays = {"squid": arrays_all[primary], "opm_matched": arrays_all["opm_matched"], "opm_dense": arrays_all["opm_dense"]}
    arrays = {k: v for k, v in arrays.items() if k in names}
    bem3 = an.subject.bem_model(com.bem)
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        gg = G3.gains(an, a, bem3)
        G_t[name], G_g[name] = gg[:, :an.nt], gg[:, an.nt:]

    def gain(name, vertices):
        return forward.chunked_discrete_gain(arrays[name].info, an.subject.trans, an.cortex.rr[vertices], an.cortex.nn[vertices], bem3,
                                             coil_def=opm.coil_def_file()).astype(np.float64)

    def patch_topographies(name, members):
        union = np.unique(np.concatenate(members))
        col = np.full(an.cortex.n, -1)
        col[union] = np.arange(len(union))
        return goldenholz.patch_topographies(gain(name, union), members, col, an.cortex.area)

    notes = (f"{an.subject.description}; Neuromag helmet at the '{cfg3['placement']['primary']}' placement "
             f"(magnetometer-to-scalp median {arrays['squid'].meta['median_dist_mm']:.1f} mm), OPM dense {arrays_all['opm_dense'].n} "
             f"and matched {arrays_all['opm_matched'].n} sites (G3B rules); adult background scale per unit area")
    return G4.Context(key, an.subject, an.cortex, an.src, com.filt, arrays, G_t, G_g, com.brain_scale, com.env, gain,
                      patch_topographies, notes)


# ----------------------------------------------------------------------------------------------
def compare(labels) -> dict:
    """Adult vs child: strength for 50 % detection per depth band and detector, the paired
    location-level OPM-vs-Neuromag results, and localization medians, side by side."""
    out = {}
    det = {lab: json.loads((OUT / f"g4_{lab}_summary.json").read_text()) for lab in ("adult",) + tuple(labels)
           if (OUT / f"g4_{lab}_summary.json").exists()}
    loc = {lab: json.loads((OUT / f"g4_localization{'' if lab == 'adult' else '_' + lab}_summary.json").read_text())
           for lab in ("adult",) + tuple(labels)
           if (OUT / f"g4_localization{'' if lab == 'adult' else '_' + lab}_summary.json").exists()}
    for lab, s in det.items():
        n_band = {}
        for L in s["locations"]:
            n_band[L["stratum"][0]] = n_band.get(L["stratum"][0], 0) + 1
        out[f"{lab}/n_locations_per_depth_band"] = {f"depth{b}": n_band.get(b, 0) for b in range(len(G4.DEPTH_BANDS))}
        out[f"{lab}/simulated_at_commit"] = s.get("simulated_at_commit")
        for d in DETECTORS:
            for mode in ("oracle", "practical@1"):
                for b in range(len(G4.DEPTH_BANDS)):
                    v = s["detectors"][d]["strength_for_50pct_nAm"].get(f"{mode}/depth{b}")
                    out[f"{lab}/{d}/{mode}/depth{b}/s50"] = v
            out[f"{lab}/{d}/heldout_false_per_min"] = s["detectors"][d]["heldout_false_per_min"]
            roc = s["detectors"][d]["roc"]["superficial_10-30mm_40nAm"]  # held-out thresholds: a matched false-event rate
            fp, se = np.array(roc["false_per_min"]), np.array(roc["sensitivity"])
            k = int(np.searchsorted(fp, 1.0, side="right")) - 1
            out[f"{lab}/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min"] = float(se[max(k, 0)])
        for a in ("opm_dense/opm", "opm_matched/opm"):
            for ref in ("squid/combined", "squid/grad", "squid/mag"):
                for mode in ("oracle", "practical@1"):
                    p = s["paired"].get(f"{a}_vs_{ref}/{mode}", {})
                    for b, r in p.items():
                        out[f"{lab}/paired/{a}_vs_{ref}/{mode}/{b}"] = dict(
                            locations_favouring_opm=r["locations_favouring_opm"], locations_favouring_squid=r["locations_favouring_squid"],
                            location_sign_flip_p=r["location_sign_flip_p"], s50_ratio_squid_over_opm=r["s50_ratio_squid_over_opm"],
                            n_locations=r["n_locations"], p_values="uncorrected")
    for lab, s in loc.items():
        for k, v in s["results"].items():
            out[f"{lab}/localization/{k}"] = {kk: v[kk] for kk in ("n", "detected", "dspm_error_mm_median_all", "ecd_error_mm_median_detected",
                                                                     "joint_detect_and_dspm_within_10mm", "joint_detect_and_ecd_within_10mm")}
        for k, v in s["paired"].items():
            out[f"{lab}/localization_paired/{k}"] = {m: v[m] for m in ("dspm_error_mm", "ecd_error_mm")}
    return out


def report(cmp: dict, labels) -> str:
    L = ["# G4 pediatric: IED detection and bounded localization in the fixed adult helmet (NEW)", "",
         "Same framework, configuration and seeds as the adult (configs/g4_epilepsy.toml); the child anatomy, arrays and "
         "placement come from G3B (Neuromag at top contact, refitted OPM arrays). The adult rows are the frozen adult study "
         "(Neuromag at its measured head position; top contact would raise it by 5.5 mm). Thresholds are calibrated on each "
         "anatomy's own null data. p-values are uncorrected (24 paired detection comparisons per OPM array and anatomy); the "
         "location is the statistical unit.", "",
         "## Strength for 50 % detection [nAm] (focal; practical detector at 1 false event/min)", "",
         "| anatomy | detector | " + " | ".join(f"{a:g}-{b:g} mm" for a, b in G4.DEPTH_BANDS) + " |",
         "|---|---|" + "---|" * len(G4.DEPTH_BANDS)]
    for lab in ("adult",) + tuple(labels):
        for d in DETECTORS:
            cells = []
            for b in range(len(G4.DEPTH_BANDS)):
                v = cmp.get(f"{lab}/{d}/practical@1/depth{b}/s50")
                if v is None:
                    cells.append("no locations")
                else:
                    ci = v["ci95"]
                    cells.append(("none" if v["value"] is None else f"{v['value']:.0f}")
                                 + (f" [{ci[0]:.0f}-{ci[1]:.0f}]" if ci and None not in ci else " [open]"))
            L.append(f"| {lab} | {d} | " + " | ".join(cells) + " |")
    L += ["", "## Held-out false events per minute at the 1-per-minute thresholds, and sensitivity at a matched held-out rate", "",
          "Sensitivity for 40-nAm spikes at 10-30 mm with each detector's threshold set on the held-out null to 1 false event per "
          "minute (the frozen thresholds give 0.5-1.7 per minute, unequal between arrays).", "",
          "| anatomy | detector | held-out rate at the frozen threshold | sensitivity at a matched 1 per minute |", "|---|---|---|---|"]
    for lab in ("adult",) + tuple(labels):
        for d in DETECTORS:
            hr = cmp.get(f"{lab}/{d}/heldout_false_per_min", {}).get("1")
            sv = cmp.get(f"{lab}/{d}/sensitivity_superficial_40nAm_at_heldout_1_per_min")
            if hr is not None and sv is not None:
                L.append(f"| {lab} | {d} | {hr:.2f} | {sv:.2f} |")
    L += ["", "## Paired OPM dense vs Neuromag combined (practical detector, 1 false event/min; locations favouring OPM / SQUID, "
          "sign-flip p, S50 ratio SQUID/OPM [95 % CI])", "",
          "| anatomy | " + " | ".join(f"{a:g}-{b:g} mm" for a, b in G4.DEPTH_BANDS) + " |", "|---|" + "---|" * len(G4.DEPTH_BANDS)]
    for lab in ("adult",) + tuple(labels):
        cells = []
        for b in range(len(G4.DEPTH_BANDS)):
            r = cmp.get(f"{lab}/paired/opm_dense/opm_vs_squid/combined/practical@1/depth{b}")
            if r is None:
                cells.append("no locations")
                continue
            sr = r["s50_ratio_squid_over_opm"]
            ratio = "n/a" if sr["value"] is None else f"{sr['value']:.2f}"
            ci = "" if not sr["ci95"] else f" [{sr['ci95'][0]:.2f}-{sr['ci95'][1]:.2f}]"
            cells.append(f"{r['locations_favouring_opm']}/{r['locations_favouring_squid']}, p {r['location_sign_flip_p']:.3g}, "
                         f"{ratio}{ci}")
        L.append(f"| {lab} | " + " | ".join(cells) + " |")
    L += ["", "## Localization (median error [mm]; joint detection + localization within 10 mm)", "",
          "| anatomy | array | condition | dSPM error (all) | ECD error (detected) | detected | joint dSPM | joint ECD |",
          "|---|---|---|---|---|---|---|---|"]
    for lab in ("adult",) + tuple(labels):
        for k, v in cmp.items():
            if k.startswith(f"{lab}/localization/"):
                arr, fam, s = k.split("/")[2:5]
                fmt = (lambda x: "-" if x is None else f"{x:.1f}")
                L.append(f"| {lab} | {arr} | {fam} {s} | {fmt(v['dspm_error_mm_median_all'])} | {fmt(v['ecd_error_mm_median_detected'])} | "
                         f"{v['detected']:.2f} | {v['joint_detect_and_dspm_within_10mm']:.2f} | {v['joint_detect_and_ecd_within_10mm']:.2f} |")
    L += ["", "Simulated IED-source recovery does not identify an epileptogenic zone or establish surgical benefit. One template "
          "is not a population; the scaled adult is a size-only control."]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("anatomy", nargs="?", default="infant2yr", choices=("infant2yr", "school", "size2yr"))
    ap.add_argument("--detection", action="store_true")
    ap.add_argument("--localization", action="store_true")
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()
    run_all = not (args.detection or args.localization or args.compare)
    t0 = time.time()
    mne.set_log_level("WARNING")
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads((ROOT / "configs" / "g4_epilepsy.toml").read_text())
    if run_all or args.detection or args.localization:
        ctx = child_context(args.anatomy)
        log(f"{args.anatomy}: context ready ({time.time() - t0:.0f} s): {ctx.notes}")
        if run_all or args.detection:
            state = G4.simulate(ctx, cfg, np.random.default_rng(cfg["simulation"]["seed"]))
            G4.save_state(state, ROOT / "cache" / "g4" / f"{args.anatomy}_state.pkl")
            G4.summarise(state)
            log(f"detection done ({time.time() - t0:.0f} s)")
        if run_all or args.localization:
            G4L.localize(ctx, cfg, np.random.default_rng(cfg["localization"]["seed"]))
            log(f"localization done ({time.time() - t0:.0f} s)")
    if run_all or args.compare:
        labels = [lab for lab in ("school", "size2yr", "infant2yr") if (OUT / f"g4_{lab}_summary.json").exists()]
        cmp = compare(labels)
        io.write_json(dict(status="NEW (G4 pediatric vs adult, same framework)", labels=["adult"] + labels, comparison=cmp),
                      OUT / "g4_pediatric_comparison.json")
        (OUT / "G4_pediatric_report.md").write_text(report(cmp, labels))
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    if "--resummarise" in sys.argv:
        lab = sys.argv[1]
        with open(ROOT / "cache" / "g4" / f"{lab}_state.pkl", "rb") as fh:
            G4.summarise(pickle.load(fh))
    else:
        main()
