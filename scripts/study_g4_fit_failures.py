#!/usr/bin/env python3
"""G4 localization: failed fits under declared criteria (docs/methods.md section 9; end-to-end goal review).

`mne.fit_dipole` returned a dipole for every event (a failure would have stopped the run), so a
"failed fit" has to be defined. This check reads the stored per-event tables of the six anatomies
(`results/g4/g4_localization*_events.csv`, nothing is recomputed) and counts, per anatomy, array
and condition, (i) gross errors: a dipole (ECD) or a dSPM peak more than 30 mm from the true source,
among all events and among the detected ones, and (ii) the dipole confidence volume above 10 cm^3
(an unconstrained fit). The 30-mm and 10-cm^3 limits are operational choices, not standards.

Output: results/g4/g4_fit_failures.json, results/g4/G4_fit_failures_report.md.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid import io  # noqa: E402

OUT = ROOT / "results" / "g4"
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo")
ARRAYS = ("squid", "opm_matched", "opm_dense")
GROSS_MM = 30.0
VOLUME_CM3 = 10.0


def read(lab: str) -> list[dict]:
    return io.read_csv(OUT / f"g4_localization{'' if lab == 'adult' else '_' + lab}_events.csv")


def summarise(rows: list[dict]) -> dict:
    groups = defaultdict(list)
    for r in rows:
        groups[(r["array"], r["family"], float(r["strength_nAm"]))].append(r)
    out = {}
    for (arr, fam, s), rs in sorted(groups.items()):
        det = [r for r in rs if r["detected"] == "True"]
        ecd = [float(r["ecd_error_mm"]) for r in rs]
        dspm = [float(r["dspm_error_mm"]) for r in rs]
        vol = [float(r["ecd_conf_vol_mm3"]) / 1e3 for r in rs]
        ecd_d = [float(r["ecd_error_mm"]) for r in det]
        dspm_d = [float(r["dspm_error_mm"]) for r in det]
        out[f"{arr}/{fam}/{s:g}nAm"] = dict(
            n=len(rs), n_detected=len(det),
            ecd_gross=sum(e > GROSS_MM for e in ecd), ecd_gross_detected=sum(e > GROSS_MM for e in ecd_d),
            dspm_gross=sum(e > GROSS_MM for e in dspm), dspm_gross_detected=sum(e > GROSS_MM for e in dspm_d),
            ecd_unconstrained=sum(v > VOLUME_CM3 for v in vol), ecd_unconstrained_detected=sum(float(r["ecd_conf_vol_mm3"]) / 1e3 > VOLUME_CM3 for r in det))
    return out


def report(res: dict) -> str:
    L = ["# G4 localization: failed fits under declared criteria", "",
         f"Every event has a dipole (`mne.fit_dipole` never failed). Counted here, from the stored per-event tables: dipole (ECD) "
         f"or dSPM-peak errors above {GROSS_MM:g} mm (gross), among all events and among the detected ones, and dipole confidence "
         f"volumes above {VOLUME_CM3:g} cm^3 (unconstrained). The limits are operational choices. Conditions: 24 events each "
         "(focal and 10-mm patches at 80 and 320 nAm).", "",
         "| anatomy | array | condition | detected | ECD > 30 mm (all / detected) | dSPM > 30 mm (all / detected) | ECD volume > 10 cm^3 (all / detected) |",
         "|---|---|---|---|---|---|---|"]
    for lab, r in res["anatomies"].items():
        for key, v in r.items():
            arr, fam, s = key.split("/")
            L.append(f"| {lab} | {arr} | {fam} {s} | {v['n_detected']}/{v['n']} | {v['ecd_gross']} / {v['ecd_gross_detected']} | "
                     f"{v['dspm_gross']} / {v['dspm_gross_detected']} | {v['ecd_unconstrained']} / {v['ecd_unconstrained_detected']} |")
    t = res["totals"]
    L += ["", f"Totals over the six anatomies ({t['n']} events per array): gross ECD errors {t['ecd_gross']} (Neuromag), "
          f"{t['ecd_gross_matched']} (matched OPM), {t['ecd_gross_dense']} (dense OPM), of which among detected events "
          f"{t['ecd_gross_detected']}, {t['ecd_gross_detected_matched']}, {t['ecd_gross_detected_dense']}; gross dSPM errors "
          f"{t['dspm_gross']}, {t['dspm_gross_matched']}, {t['dspm_gross_dense']} (detected: {t['dspm_gross_detected']}, "
          f"{t['dspm_gross_detected_matched']}, {t['dspm_gross_detected_dense']}). Of all gross errors (both estimators, three "
          f"arrays), {100 * t['gross_undetected_share']:.0f} % belong to undetected events and {100 * t['gross_80nAm_share']:.0f} % "
          "to 80-nAm sources, whose estimates are largely those of noise."]
    return "\n".join(L) + "\n"


def main():
    res = dict(status="NEW (G4 localization: failed fits under declared criteria, from the stored per-event tables)",
               criteria=dict(gross_error_mm=GROSS_MM, unconstrained_volume_cm3=VOLUME_CM3), anatomies={})
    tot = defaultdict(int)
    for lab in LABELS:
        rows = read(lab)
        res["anatomies"][lab] = summarise(rows)
        for arr, suffix in (("squid", ""), ("opm_matched", "_matched"), ("opm_dense", "_dense")):
            rs = [r for r in rows if r["array"] == arr]
            det = [r for r in rs if r["detected"] == "True"]
            tot["n" if arr == "squid" else f"n{suffix}"] += len(rs)
            tot[f"ecd_gross{suffix}"] += sum(float(r["ecd_error_mm"]) > GROSS_MM for r in rs)
            tot[f"ecd_gross_detected{suffix}"] += sum(float(r["ecd_error_mm"]) > GROSS_MM for r in det)
            tot[f"dspm_gross{suffix}"] += sum(float(r["dspm_error_mm"]) > GROSS_MM for r in rs)
            tot[f"dspm_gross_detected{suffix}"] += sum(float(r["dspm_error_mm"]) > GROSS_MM for r in det)
    gross = [(r["detected"] == "True", float(r["strength_nAm"]) == 80.0)
             for lab in LABELS for r in read(lab) for err in ("ecd_error_mm", "dspm_error_mm") if float(r[err]) > GROSS_MM]
    tot["gross_undetected_share"] = sum(not d for d, _ in gross) / len(gross)
    tot["gross_80nAm_share"] = sum(w for _, w in gross) / len(gross)
    res["totals"] = dict(tot)
    io.write_json(res, OUT / "g4_fit_failures.json")
    (OUT / "G4_fit_failures_report.md").write_text(report(res))
    print(report(res))


if __name__ == "__main__":
    main()
