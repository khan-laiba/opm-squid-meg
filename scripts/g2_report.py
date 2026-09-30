#!/usr/bin/env python3
"""G2 report: one section per sensor configuration, generated from the G2 outputs.

Reads results/g2/g2_summary.json (scripts/g2_adult_comparison.py) and, if present,
results/g2/g2_band_sensitivity.json (scripts/g2_band_sensitivity.py); writes results/g2/G2_report.md.
Every number is taken from those files; nothing is recomputed here.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "g2"
CONFIGS = [("squid", "mag", "Neuromag magnetometers (102)"), ("squid", "grad", "Neuromag planar gradiometers (204)"),
           ("squid", "combined", "Neuromag combined (306)"), ("opm_matched", "opm", "OPM matched sites"),
           ("opm204", "opm", "OPM channel-budget control (204)"), ("opm_dense", "opm", "OPM dense array (full system)")]
OPMS = ("opm_matched", "opm204", "opm_dense")
REFS = ("combined", "grad", "mag")


def ratio(v):
    share = f", {100 * v['share_parcels_opm_better']:.0f} % of parcels" if "share_parcels_opm_better" in v else ""
    return (f"{2 ** v['median_log2']:.2f}x [{2 ** v['ci95'][0]:.2f}-{2 ** v['ci95'][1]:.2f}], OPM higher for "
            f"{100 * v['share_opm_better']:.0f} % of targets{share}")


def depth_table(d, key):
    rows = [x for x in d["detectability_vs_depth"][key] if x["median"] is not None]
    head = "| depth [mm] | " + " | ".join(f"{int(x['lo'])}-{int(x['hi'])}" for x in rows) + " |"
    sep = "|---|" + "---|" * len(rows)
    val = "| median detectability (10 nAm) | " + " | ".join(f"{x['median']:.2f}" for x in rows) + " |"
    return "\n".join([head, sep, val])


def main():
    d = json.loads((OUT / "g2_summary.json").read_text())
    bands_file = OUT / "g2_band_sensitivity.json"
    bands = json.loads(bands_file.read_text()) if bands_file.exists() else None
    conds = d["config"]["conditions"]["all"]
    headline = d["config"]["conditions"]["headline"]
    v = d["noise_validation"]
    L = []
    L.append("# G2 report: realistic adult OPM vs Neuromag comparison (NEW)\n")
    L.append(f"Generated from `results/g2/g2_summary.json` (code commit {d['provenance']['commit']}, MNE {d['provenance']['mne_version']}). "
             "Methods: `docs/methods.md` section 8; assumptions in `docs/provenance_register.md`. This is a proposed study; "
             "no author of the reproduced papers has reviewed it.\n")
    L.append("## Common setup\n")
    L.append(f"- Anatomy: MNE sample subject, measured head position; {d['n_targets']} target dipoles (10 nAm, cortical normal, "
             f"usable oct-6 vertices); background grid of {d['n_background_grid']} area-weighted sources.")
    L.append(f"- Band 1-40 Hz (ENBW {d['enbw_hz']:.1f} Hz). Intrinsic noise: SQUID magnetometers "
             f"{v['model']['intrinsic_rms_mag_fT']:.1f} fT, gradiometers {v['model']['intrinsic_rms_grad_fT_cm']:.1f} fT/cm; OPM at "
             "15 fT/sqrt(Hz) " + f"{v['model']['intrinsic_rms_opm_fT']['15']:.0f} fT (RMS in band).")
    L.append(f"- Brain noise calibrated on good gradiometers to {v['measured']['brain_rms_grad_fT_cm']:.1f} fT/cm (task baseline "
             f"minus empty room); predicted magnetometer level {v['model']['brain_rms_mag_fT']:.0f} fT vs {v['measured']['brain_rms_mag_fT']:.0f} fT "
             "measured.")
    L.append(f"- Room field: explains {100 * v['environment_explained_fraction']['mag']:.0f} % (magnetometers) and "
             f"{100 * v['environment_explained_fraction']['grad']:.1f} % (gradiometers) of the empty-room variance; model vs measured "
             f"empty room {v['model']['empty_room_rms_mag_fT']:.0f} vs {v['measured']['empty_room_rms_mag_fT']:.0f} fT, "
             f"{v['model']['empty_room_rms_grad_fT_cm']:.1f} vs {v['measured']['empty_room_rms_grad_fT_cm']:.1f} fT/cm.")
    L.append("- Conditions: " + ", ".join(conds) + ". Detectability = known-topography matched-filter SNR sqrt(s^T C^+ s); "
             "not an event detection rate or a localization accuracy.\n")
    for name, cs, title in CONFIGS:
        L.append(f"## {title}\n")
        geo = d["arrays"][name]
        br = d["bridge_to_sphere"]["sensor_distance_mm"][name]
        n_ch = {"mag": 102, "grad": 204}.get(cs, geo["channels"]) if name == "squid" else geo["channels"]
        L.append(f"- Geometry: {n_ch} channels" + (f", {geo.get('n_sites')} sites" if geo.get("n_sites") else "")
                 + f", {geo.get('axes', '102 sites with 1 magnetometer + 2 planar gradiometers')}; scalp-to-sensor distance median "
                 f"{br['median']:.1f} mm (5-95 %: {br['p5']:.1f}-{br['p95']:.1f} mm)"
                 + (f"; role: {geo['role']}" if geo.get("role") else "") + ".")
        rk = f"squid_{cs}" if name == "squid" and cs != "combined" else name
        L.append("- Retained rank: " + ", ".join(f"{c} {d['retained_rank'].get(f'{rk}/{c}', d['retained_rank'][f'{name}/{c}'])}" for c in conds)
                 + (" (channel subset after the full-array projection)." if rk.startswith("squid_") else "."))
        comp_key = name if name != "squid" else f"squid_{cs}" if cs != "combined" else None
        if comp_key and comp_key in d["noise_composition"]:
            c = d["noise_composition"][comp_key]
            unit = 1e13 if cs == "grad" else 1e15
            u = "fT/cm" if cs == "grad" else "fT"
            L.append(f"- Noise composition (median RMS): intrinsic {c['intrinsic_rms'] * unit:.1f} {u}, brain {c['brain_rms'] * unit:.1f} {u}, "
                     f"room {c['env_rms'] * unit:.1f} {u}.")
        for cond in headline:
            L.append(f"\nDetectability vs depth, {cond}:\n")
            L.append(depth_table(d, f"{name}/{cs}/{cond}"))
        if name != "squid":
            L.append("\nOPM / Neuromag, median detectability ratio (95 % CI from a bootstrap over cortical parcels):\n")
            L.append("| condition | vs combined | vs gradiometers | vs magnetometers |\n|---|---|---|---|")
            for cond in conds:
                L.append(f"| {cond} | " + " | ".join(ratio(d["primary"]["oracle"][f"{name}/{ref}/{cond}"]) for ref in REFS) + " |")
            pl = d["plugin_over_oracle_median"]
            L.append("\n- Estimated (plug-in) covariance, detectability relative to oracle: "
                     + ", ".join(f"{c}: {pl[f'{name}/opm/{c}/T10']:.2f} (10 s), {pl[f'{name}/opm/{c}/T60']:.2f} (60 s)" for c in headline) + ".")
            pa = d["patches"]["comparisons"]
            L.append("- Patches vs Neuromag combined (intrinsic+brain): " + ", ".join(
                f"{r:g} mm {2 ** pa[f'{name}/combined/intrinsic+brain/{r:g}mm']['median_log2']:.2f}x" for r in d["patches"]["radii_mm"]) + ".")
            s = d["sensitivity"]
            if name in ("opm_matched", "opm_dense") and "sensitivity_joint_asd_gap" in d:
                L.append("- Joint OPM noise x scalp gap (vs combined, intrinsic+brain): " + "; ".join(
                    f"{k.replace('/', ', ')} {v[name]['ratio']:.2f}x" for k, v in d["sensitivity_joint_asd_gap"].items()) + ".")
            if "break_even_opm_asd_fT" in d:
                L.append("- Intrinsic noise only, break-even OPM noise (ratio = 1): " + ", ".join(
                    f"vs {ref} {d['break_even_opm_asd_fT'][f'{name}/{ref}']:.1f} fT/sqrt(Hz)" for ref in REFS) + ".")
            L.append("- Sensitivity, one factor at a time (vs combined, intrinsic+brain): " + "; ".join(
                f"{k.replace('_', ' ')} {2 ** s[k][f'{name}/combined/intrinsic+brain']['median_log2']:.2f}x" for k in s
                if f"{name}/combined/intrinsic+brain" in s[k]) + ".")
            lobes = d["by_lobe"].get(f"{name}/combined/intrinsic+brain")
            if lobes:
                L.append("- By lobe (vs combined, intrinsic+brain): " + ", ".join(f"{lb} {2 ** x['median_log2']:.2f}x" for lb, x in lobes.items()) + ".")
            mw = d.get("medial_wall")
            if mw:
                L.append(f"- Without the {mw['n_targets']} medial-wall targets ({100 * mw['share']:.1f} %; FreeSurfer 'unknown', not cortex), "
                         "vs combined: " + ", ".join(f"{c} {mw['ratios_without'][f'{name}/combined/{c}']:.2f}x" for c in headline) + ".")
            if bands:
                L.append("- Frequency bands (vs combined, intrinsic+brain; noise recalibrated per band): " + ", ".join(
                    f"{k} {b['ideal'][f'{name}/combined/intrinsic+brain']['ratio']:.2f}x" for k, b in bands["bands"].items()) + ".")
        else:
            pl = d["plugin_over_oracle_median"]
            L.append("\n- Estimated (plug-in) covariance, detectability relative to oracle: "
                     + ", ".join(f"{c}: {pl[f'squid/{cs}/{c}/T10']:.2f} (10 s), {pl[f'squid/{cs}/{c}/T60']:.2f} (60 s)" for c in headline) + ".")
        L.append("")
    L.append("## Convergence\n")
    c = d["convergence"]
    L.append(f"- Background grid vs every usable vertex: median log2 ratios change by <= {c['background_grid_vs_fullres']['max_abs_change_log2']:.3f}.")
    L.append(f"- BEM 5,120 vs 20,480 triangles (1 layer): <= {c['bem']['refinement_max_abs_change_log2']:.3f}; 3 vs 1 layer: "
             f"<= {c['bem']['bem1_5120']['max_abs_change_log2']:.3f}.")
    L.append(f"- Oct-6 vs random full-resolution targets: <= {c['target_sampling']['max_abs_change_log2']:.3f}; Neuromag 4-point vs "
             f"accurate integration: {100 * c['coil_integration']['squid_4pt_vs_accurate_rel_diff_median']:.1f} % in the gains.")
    L.append("\n## Limitations\n")
    for n in d["notes"]:
        L.append(f"- {n}")
    L.append("- One adult anatomy and one measured head position; between-subject variability is not represented.")
    L.append("- OPM intrinsic noise is a declared sweep, not a device specification; OPM movement artefacts, cross-talk and "
             "calibration errors are not modelled.")
    (OUT / "G2_report.md").write_text("\n".join(L) + "\n")
    print(f"wrote {OUT / 'G2_report.md'} ({len(L)} lines)")


if __name__ == "__main__":
    sys.exit(main())
