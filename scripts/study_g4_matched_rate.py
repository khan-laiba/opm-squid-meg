#!/usr/bin/env python3
"""G4 detection at matched held-out false-event rates (docs/methods.md sections 9 and 11).

The practical detector's thresholds are frozen on 20 min of calibration null data; on the 20 min of
held-out null data they give 0.4-1.7 false events per minute, unequal between arrays. This check
re-evaluates the stored simulations (cache/g4/<label>_state.pkl, written by the G4 runs; nothing is
simulated again) with each detector's threshold set on the held-out null to 1 false event per minute,
and reports the paired location-level comparisons and the strengths for 50 % detection at that
matched rate next to the frozen-threshold values. The matched thresholds are fitted on the held-out
data, so this is an in-sample check of the frozen-threshold results, not a replacement for them.
Before the matched evaluation it re-summarises every unmodified state and requires the committed
paired results to be reproduced exactly.

Outputs: results/g4/g4_matched_rate.json, results/g4/G4_matched_rate_report.md.
"""
from __future__ import annotations

import copy
import json
import pickle
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import g4_epilepsy_adult as G4  # noqa: E402
from opmsquid import detection, io, paths  # noqa: E402

OUT = ROOT / "results" / "g4"
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo")
PAIRS = (("opm_dense/opm", "squid/combined"), ("opm_matched/opm", "squid/combined"))
DETECTORS = ("squid/combined", "opm_matched/opm", "opm_dense/opm")


def resummarise(state: dict, tmp: Path) -> dict:
    G4.OUT = tmp
    G4.summarise(state)
    return json.loads((tmp / f"g4_{state.get('label', 'adult')}_summary.json").read_text())


def matched(state: dict) -> dict:
    st = copy.deepcopy(state)
    for key in st["thr"]:
        st["thr"][key] = dict(st["thr"][key], **{"1": detection.threshold_for_rate(st["held_heights"][key], st["minutes_held"], 1.0)})
    return st


def pick(summary: dict) -> dict:
    out = {}
    for a, b in PAIRS:
        for band, r in summary["paired"][f"{a}_vs_{b}/practical@1"].items():
            out[f"{a}_vs_{b}/{band}"] = dict(locations_favouring_opm=r["locations_favouring_opm"],
                                             locations_favouring_squid=r["locations_favouring_squid"],
                                             location_sign_flip_p=r["location_sign_flip_p"],
                                             s50_ratio_squid_over_opm=r["s50_ratio_squid_over_opm"])
    for d in DETECTORS:
        for band in range(len(G4.DEPTH_BANDS)):
            out[f"s50/{d}/depth{band}"] = summary["detectors"][d]["strength_for_50pct_nAm"][f"practical@1/depth{band}"]
        out[f"heldout_false_per_min/{d}"] = summary["detectors"][d]["heldout_false_per_min"]["1"]
    return out


def main():
    result = dict(status="NEW check (review): G4 practical detector at thresholds matched on the held-out null (1 false event/min); "
                         "in-sample for the held-out data", labels=list(LABELS), anatomies={})
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for lab in LABELS:
            with open(paths.CACHE / "g4" / f"{lab}_state.pkl", "rb") as fh:
                state = pickle.load(fh)
            committed = json.loads((OUT / f"g4_{lab}_summary.json").read_text())
            if state.get("simulated_at_commit") != committed.get("simulated_at_commit"):
                raise SystemExit(f"{lab}: cached state ({state.get('simulated_at_commit')}) is not the committed simulation")
            frozen = pick(resummarise(state, tmp))
            if frozen != pick(committed):
                raise SystemExit(f"{lab}: re-summarising the cached state does not reproduce the committed summary")
            m = pick(resummarise(matched(state), tmp))
            result["anatomies"][lab] = dict(simulated_at_commit=state["simulated_at_commit"], frozen=frozen, matched=m,
                                            thresholds=dict(frozen={d: state["thr"][d]["1"] for d in DETECTORS},
                                                            matched={d: matched(state)["thr"][d]["1"] for d in DETECTORS}))
            G4.log(f"{lab}: done")
    io.write_json(result, OUT / "g4_matched_rate.json")
    (OUT / "G4_matched_rate_report.md").write_text(report(result))


def cell(r: dict) -> str:
    return (f"{r['locations_favouring_opm']}/{r['locations_favouring_squid']}, p {r['location_sign_flip_p']:.2g}, "
            + detection.format_s50_ratio(r["s50_ratio_squid_over_opm"]))


def report(res: dict) -> str:
    L = ["# G4: practical detector at matched held-out false-event rates (check)", "",
         "The frozen thresholds (set on 20 min of calibration null data) give unequal false-event rates on the 20 min of "
         "held-out null data. Here every detector's threshold is set on the held-out null to 1 false event per minute and the "
         "stored simulations are re-evaluated (no new simulation; in-sample for the held-out data). Paired dense OPM vs "
         "Neuromag combined and matched OPM vs Neuromag combined: locations favouring OPM / Neuromag, exact sign-flip p "
         "(uncorrected), strength ratio Neuromag / OPM [95 % location bootstrap].", "",
         "| anatomy | pair | thresholds | " + " | ".join(f"{a:g}-{b:g} mm" for a, b in G4.DEPTH_BANDS) + " |",
         "|---|---|---|" + "---|" * len(G4.DEPTH_BANDS)]
    for lab, r in res["anatomies"].items():
        for a, b in PAIRS:
            for kind in ("frozen", "matched"):
                L.append(f"| {lab} | {a.split('/')[0]} vs {b} | {kind} | "
                         + " | ".join(cell(r[kind][f"{a}_vs_{b}/depth{i}"]) for i in range(len(G4.DEPTH_BANDS))) + " |")
    L += ["", "Held-out false events per minute at the frozen thresholds (matched: 1.00 by construction):", "",
          "| anatomy | " + " | ".join(DETECTORS) + " |", "|---|" + "---|" * len(DETECTORS)]
    for lab, r in res["anatomies"].items():
        L.append(f"| {lab} | " + " | ".join(f"{r['frozen'][f'heldout_false_per_min/{d}']:.2f}" for d in DETECTORS) + " |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
