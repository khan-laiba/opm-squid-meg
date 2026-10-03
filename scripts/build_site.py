#!/usr/bin/env python3
"""G5: build the local static report from the committed result files (never deployed here).

Order as GOAL.md G5 asks: adult reference benchmarks (G1), realistic adult results (G2), pediatric
extension (G3), epilepsy (G4), then methods, uncertainty, limitations and reproduction. Every
number is read from results/*/ (nothing is recomputed); figures are the committed PNGs; JSON, CSV
and the G2 report are copied for download with their size, SHA-256 and the code commit recorded
in the result file. The G1A Fig. 3 PDF and SVG (licensed font subsets and outlines) are not used.

Usage: .venv/bin/python scripts/build_site.py [--out DIR]   (default site/_build, git-ignored)
The build fails if any local link or anchor is broken.
"""
from __future__ import annotations

import argparse
import html
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import jinja2  # noqa: E402
from markupsafe import Markup  # noqa: E402

from opmsquid import detection, io, sitebuild as sb  # noqa: E402

RES = ROOT / "results"
NAV = [("index.html", "Overview"), ("benchmarks.html", "Adult benchmarks (G1)"), ("adult.html", "Realistic adult (G2)"),
       ("pediatric.html", "Pediatric (G3)"), ("epilepsy.html", "Epilepsy (G4)"), ("methods.html", "Methods and limitations"),
       ("register.html", "Parameters and provenance"), ("reproduce.html", "Reproduce and download")]
FROZEN_TAG = "adult-baseline-v2"  # the adult baseline frozen before any pediatric outcome was compared
CURRENT_TAG = "adult-baseline-v4"  # every head-model- and array-dependent result recomputed with equal OPM standoff (adult and pediatric alike)
MARKER = ".opmsquid_site_build"  # marks an output directory as the builder's own (safe to replace)
LABEL = {"squid": "Neuromag", "opm_matched": "OPM matched", "opm204": "OPM 204 (channel budget)",
         "opm_dense": "OPM dense", "combined": "combined", "grad": "gradiometers", "mag": "magnetometers",
         "squid_mag": "Neuromag magnetometers", "squid_grad": "Neuromag gradiometers"}
DETECTORS = {"squid/combined": "Neuromag combined", "squid/grad": "Neuromag gradiometers", "squid/mag": "Neuromag magnetometers",
             "opm_matched/opm": "OPM matched", "opm_dense/opm": "OPM dense"}
DEPTH_BANDS = ("10-20 mm", "20-30 mm", "30-45 mm", "45-70 mm")
ANAT = {"adult": "adult (sample subject)", "school": "school-age size (scaled adult)", "size2yr": "2-year size (scaled adult)",
        "infant2yr": "2-year template", "infant18mo": "18-month template", "infant12mo": "12-month template",
        "childA": "child A (7.8 y)", "childB": "child B (8.3 y)", "childC": "child C (8.7 y)"}
CHILDREN = ("school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
TEMPLATES = ("infant2yr", "infant18mo", "infant12mo")
SCHOOL = ("childA", "childB", "childC")  # individual school-aged children (OpenNeuro ds005234)


def load(rel):
    p = RES / rel
    return json.loads(p.read_text()) if p.exists() else None


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def x(v, nd=2):
    return f"{v:.{nd}f}x"


def cmp_str(c):
    return f"{2 ** c['median_log2']:.2f}x [{2 ** c['ci95'][0]:.2f}-{2 ** c['ci95'][1]:.2f}]"


def fig(out, rel_png, caption_html, alt):
    """Copy a committed PNG into the build and return its <figure>."""
    src = RES / rel_png
    dst = out / "figures" / rel_png
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return sb.figure(f"figures/{rel_png}", caption_html, alt)


def label(kind):
    return f'<span class="label">{kind}</span>'


# ------------------------------------------------------------------------------------------ pages
def page_index(d):
    g2, g4, loc, g1a, g1b, g1c = d["g2"], d["g4"], d["loc"], d["g1a"], d["g1b"], d["g1c"]
    P = g2["primary"]["oracle"]
    dense, matched = P["opm_dense/combined/intrinsic+brain"], P["opm_matched/combined/intrinsic+brain"]
    sens = g2["sensitivity"]
    joint = g2["sensitivity_joint_asd_gap"]
    bem1 = 2 ** g2["convergence"]["bem"]["bem1_5120"]["median_log2"]["opm_dense/combined/intrinsic+brain"]
    depth = {f"{int(r['lo'])}-{int(r['hi'])}": 2 ** r["median"] for r in g2["log2_ratio_vs_depth"]["opm_dense/combined/intrinsic+brain"]
             if r["median"] is not None}
    s50 = {k: g4["detectors"][k]["strength_for_50pct_nAm"] for k in ("squid/combined", "opm_dense/opm", "opm_matched/opm")}
    pr = g4["paired"]["opm_dense/opm_vs_squid/combined/practical@1"]
    pm = g4["paired"]["opm_matched/opm_vs_squid/combined/practical@1"]
    L = loc["results"]
    status = [("G0 audit, provenance, plan", "-", "done"),
              ("G1A Jas et al. 2026 analytical benchmark", "REPRO", "done, internally reviewed"),
              ("G1B Hunold et al. 2016 depth-orientation spike SNR (MEG)", "ADAPT (+ NEW OPM column)", "done, internally reviewed"),
              ("G1C Goldenholz et al. 2009 cortical SNR maps (MEG)", "ADAPT (+ NEW OPM extension)", "done, internally reviewed"),
              ("G2 realistic adult OPM vs Neuromag", "NEW", f"done, internally reviewed; frozen as {FROZEN_TAG} before the pediatric "
               f"outcomes; recomputed with the whole-cell OPM clearance after the final reviews (adult-baseline-v3) and with every "
               f"anatomy's head surface on its MRI scalp, so that every OPM array has the same standoff, as {CURRENT_TAG}"),
              ("G3A Jas head-size benchmark", "REPRO (+ NEW fixed shell)", "done" if d["g3a"] else "not run"),
              ("G3B pediatric fixed helmet vs head-adaptive OPM", "NEW",
               "done (12-, 18- and 24-month infant templates, three school-aged children and scaled-adult size controls)"
               if d.get("g3b") else "in progress"),
              ("G4 epilepsy detection and bounded localization", "NEW",
               ("adult done; pediatric done (" + ", ".join(ANAT[x] for x in d["g4p"]["labels"] if x != "adult") + ")"
                if d.get("g4p") else "adult done; pediatric in progress")
               + ("; head motion and OPM slippage: bounded extension done" if d.get("motion") else "")),
              ("G5 software, reproduction, report", "-", "local report (not deployed); 148 tests; a clean clone of the result commit "
               "passes them (one skip without the lead-field cache); a fresh-environment smoke test passed on 2026-09-30 at an earlier "
               "commit")]
    h = ["<p>Simulation study comparing on-scalp optically pumped magnetometers (OPM) with the Neuromag SQUID system: "
         "an analytical benchmark and adaptations of two published adult studies, a realistic adult comparison, a pediatric "
         "fixed-helmet versus head-adaptive extension, and interictal-spike detection and localization examples in both. The "
         "OPM advantage is tested, not assumed. Labels: "
         f"{label('REPRO')} reproduction with the paper's definitions, {label('ADAPT')} adaptation where data or details are "
         f"unavailable, {label('NEW')} new experiment.</p>",
         "<h2 id=\"status\">Status</h2>", sb.table(["Milestone", "Label", "Status"], status),
         "<h2 id=\"findings\">Adult findings, with their conditions</h2>"]
    n_dense, n_matched = g2["arrays"]["opm_dense"]["channels"], g2["arrays"]["opm_matched"]["channels"]

    def bands(p):
        return "; ".join(f"{DEPTH_BANDS[b]} {p[f'depth{b}']['locations_favouring_opm']}/{p[f'depth{b}']['locations_favouring_squid']} "
                         f"(p = {p[f'depth{b}']['location_sign_flip_p']:.2g})" for b in range(4))

    r0 = pr["depth0"].get("s50_ratio_squid_over_opm")
    ratio_txt = (f" (Neuromag / dense strength ratio and paired location-bootstrap 95 % interval "
                 f"{detection.format_s50_ratio(r0)})") if r0 else ""
    pl = loc["paired"]
    dspm = {a: pl[f"{a}_vs_squid/patch/320nAm"]["dspm_error_mm"] for a in ("opm_matched", "opm_dense")}
    items = [
        f"With modelled brain noise, the dense {n_dense}-site OPM array has <strong>{cmp_str(dense)}</strong> the known-topography "
        f"detectability of the Neuromag system (all 306 channels; median over {g2['n_targets']:,} cortical targets, 95 % CI "
        f"from a bootstrap over cortical parcels; higher for {100 * dense['share_opm_better']:.0f} % of targets). An OPM array at "
        f"Neuromag's own sites ({n_matched} of 102 fit) ties it: {cmp_str(matched)}. Conditions: OPM white noise 15 fT/&radic;Hz, "
        "3-layer BEM, no extra scalp gap, the sample subject's measured head position; the modelled brain noise predicts "
        f"{g2['noise_validation']['model']['brain_mag_model_over_measured']:.2f}x the measured magnetometer level.",
        f"The advantage depends on the assumptions: OPM noise 7 or 30 fT/&radic;Hz gives "
        f"{x(2 ** sens['opm_asd_7fT']['opm_dense/combined/intrinsic+brain']['median_log2'])} or "
        f"{x(2 ** sens['opm_asd_30fT']['opm_dense/combined/intrinsic+brain']['median_log2'])}; with 30 fT/&radic;Hz and a 6-mm "
        f"scalp gap together, {x(joint['gap6mm/asd30fT']['opm_dense']['ratio'])}; a 1-layer head model gives {x(bem1)}. Without "
        f"brain noise (intrinsic sensor noise only) Neuromag wins: {cmp_str(P['opm_dense/combined/intrinsic'])}.",
        f"By depth (dense vs Neuromag): {x(depth['10-15'])} at 10-15 mm below the scalp, {x(depth['30-35'])} at 30-35 mm and "
        f"{x(depth['50-55'])} at 50-55 mm; the matched array goes from "
        f"{x(2 ** g2['log2_ratio_vs_depth']['opm_matched/combined/intrinsic+brain'][0]['median'])} near the scalp to "
        f"{x(2 ** min(r['median'] for r in g2['log2_ratio_vs_depth']['opm_matched/combined/intrinsic+brain'] if r['median'] is not None))} "
        "at its lowest depth bin.",
        f"Simulated interictal spikes, practical detector at 1 false event per minute: 50 % detection needs "
        f"{s50['opm_dense/opm']['practical@1/depth0']['value']:.0f} nAm with the dense OPM array vs "
        f"{s50['squid/combined']['practical@1/depth0']['value']:.0f} nAm with Neuromag at 10-20 mm depth{ratio_txt}. Locations "
        f"with more events detected only by the dense array / only by Neuromag (exact sign-flip p, uncorrected): {bands(pr)}. "
        f"Matched array: {bands(pm)}.",
        f"Bounded localization (24 locations, one event each): for detected 320-nAm focal spikes the dipole error on the MRI is "
        f"{L['squid/focal/320nAm']['ecd_error_mm_median_detected']:.1f}, {L['opm_matched/focal/320nAm']['ecd_error_mm_median_detected']:.1f} "
        f"and {L['opm_dense/focal/320nAm']['ecd_error_mm_median_detected']:.1f} mm (Neuromag, matched, dense), limited by a "
        "2-mm/2-deg coregistration error in every inverse. For 320-nAm patches the dSPM error changes by "
        f"{dspm['opm_matched']['median_difference']:+.1f} mm (matched, p = {dspm['opm_matched']['wilcoxon_p']:.2g}) and "
        f"{dspm['opm_dense']['median_difference']:+.1f} mm (dense, p = {dspm['opm_dense']['wilcoxon_p']:.2g}) relative to "
        f"Neuromag (uncorrected; {sum(1 for k in pl if k.startswith('opm_dense_vs_squid/')) * sum(1 for m in ('dspm_error_mm', 'ecd_error_mm', 'dspm_mne_error_mm', 'joint_dspm_10mm', 'joint_ecd_10mm') if m in pl['opm_dense_vs_squid/patch/320nAm'])} "
        "comparisons per array against Neuromag combined).",
        f"Benchmarks: the analytical sphere model reproduces Jas et al. (equal-SNR depth "
        f"{g1a['d_eq_mm']['3']:.3f} mm at &eta; = 3, printed about 28 mm); the Hunold adaptation reproduces the depth-orientation "
        f"pattern of the published maps (r = {min(v['pearson_r'] for f in ('dipole/p2p', 'patch/p2p') for v in g1b['variants']['fig6_calibrated']['comparison_with_paper'][f].values()):.2f}-"
        f"{max(v['pearson_r'] for f in ('dipole/p2p', 'patch/p2p') for v in g1b['variants']['fig6_calibrated']['comparison_with_paper'][f].values()):.2f}); "
        f"the Goldenholz adaptation gives a noise-source strength of "
        f"{g1c['comparison_with_paper']['source_sd_nAm_at_4000_sources']['probable_default']:.2f} nAm (paper 1.6-1.9).",
    ]
    h.append("".join(f'<div class="finding">{t}</div>' for t in items))
    h.append("<h2 id=\"not-shown\">What these results do not show</h2><ul>"
             "<li>One adult anatomy and one measured head position: no between-subject variability.</li>"
             "<li>OPM intrinsic noise is a declared sweep (7-30 fT/&radic;Hz), not a device specification; cross-talk is not "
             "modelled, and head motion, cap slippage and calibration errors only in a bounded extension "
             "(<a href=\"epilepsy.html#motion\">motion</a>) that does not establish motion robustness.</li>"
             "<li>Detectability is a known-topography matched-filter SNR, not a clinical detection rate; the spike study "
             "uses simulated events in simulated noise.</li>"
             "<li>Confidence intervals resample cortical parcels or locations of one anatomy; they do not include model "
             "uncertainty, which the sensitivity analyses show instead (one factor at a time, plus a joint noise x gap grid).</li>"
             "<li>Pediatric results rest on three average templates of one database (12, 18 and 24 months), three individual "
             "school-aged children of one dataset (7.8-8.7 years; their skull modelled, their fiducials transferred from the "
             "adult) and two scaled copies of the adult: variability between three children is shown, not estimated; no "
             "age-specific background physiology (only a bounded sensitivity), adult conductivities.</li></ul>")
    if d.get("g3b"):
        k = next(i for i, x in enumerate(h) if x.startswith('<h2 id="not-shown">'))
        h.insert(k, pediatric_findings(d))
    return "\n".join(h)


def sens_sign(sens: dict, kids) -> str:
    """'positive' if every sensitivity variant leaves the child-minus-adult difference of the median D above 0, else the
    exceptions (the text must stay true for whatever the results are)."""
    neg = []
    for key, v in sens.items():
        lab, variant = key.split("/")[0], key.split("/")[1]
        if lab in kids and key.endswith("/opm_dense/combined/intrinsic+brain"):
            a = sens.get(key.replace(lab + "/", "adult/", 1))
            if a is not None and v - a <= 0:
                neg.append(f"{ANAT[lab]}, {variant}")
    return "positive" if not neg else "positive except " + "; ".join(neg)


def pediatric_findings(d):
    g3b, g4p = d["g3b"], d.get("g4p")
    C, dec, sens = g3b["comparisons"], g3b["delta_other_placements"], g3b["sensitivity_median_D_dB"]

    def ci(s):
        c = s.get("ci95")
        return f"{s['median']:+.2f} dB" + (f" [{c[0]:+.2f}, {c[1]:+.2f}]" if c else "")

    kids = [c for c in CHILDREN if f"{c}/opm_dense/combined/intrinsic+brain/detect" in C]
    r = {c: C[f"{c}/opm_dense/combined/intrinsic+brain/detect"] for c in kids}
    items = [
        f"In the fixed Neuromag helmet, raised to 20-mm contact with its top, the dense OPM array's known-topography "
        f"detectability relative to Neuromag combined (D, in dB) grows from <strong>{r[kids[0]]['d_adult']['median']:+.2f} dB</strong> "
        "in the adult to " + "; ".join(f"{r[c]['d_child']['median']:+.2f} dB ({ANAT[c]})" for c in kids) + ". Delta = "
        + "; ".join(f"{ci(r[c]['delta'])} ({ANAT[c]})" for c in kids)
        + " (vertex-wise for the scaled adults; by parcel for the templates and the school-aged children). OPM arrays refitted with "
        "the adult rules "
        "(nothing shrunk); background, room field and sensor noise unchanged.",
        "The gain comes mainly from the fixed helmet's fit: left at the adult's ear-line position, Delta is "
        + ", ".join(f"{dec[f'{c}/centred_vs_adult_centred/combined']['delta']['median']:+.2f}" for c in kids)
        + " dB; in a counterfactual helmet scaled with the head it is "
        + ", ".join(f"{dec[f'{c}/counterfactual_vs_adult_counterfactual/combined']['delta']['median']:+.2f}" for c in kids)
        + " dB, and "
        + ", ".join(f"{dec[f'{c}/counterfactual_x-centred_vs_adult_counterfactual_x-centred/combined']['delta']['median']:+.2f}"
                    for c in kids)
        + " dB about the laterally centred head (same order). With the background fixed per unit cortical area, both systems' detectability "
        "rises in the smaller heads and the on-scalp OPM's rises more; in a helmet scaled with the head the SQUID, its "
        "gradiometers most, gains as much or slightly more.",
        f"The child-minus-adult difference of the median D is {sens_sign(sens, kids)} for OPM noise 7-30 fT/&radic;Hz, background "
        f"variance x0.5 or x2 and a 1-layer head model; at 30 fT/&radic;Hz the adult's D is "
        f"{sens['adult/opm_asd_30fT/opm_dense/combined/intrinsic+brain']:+.2f} dB and the children's "
        + ", ".join(f"{sens[f'{k}/opm_asd_30fT/opm_dense/combined/intrinsic+brain']:+.2f}" for k in kids
                    if f"{k}/opm_asd_30fT/opm_dense/combined/intrinsic+brain" in sens)
        + " dB (same order). A positive Delta is a relative gain for the head-adaptive array, not by itself a clinical advantage.",
    ]
    if g4p:
        cmp_ = g4p["comparison"]

        def s50(lab, key):
            v = cmp_.get(f"{lab}/{key}/practical@1/depth0/s50")
            return "-" if not v or v["value"] is None else f"{v['value']:.0f}"

        labs = [x for x in g4p["labels"] if x != "adult"]
        items.append("Simulated spikes in the same framework (practical detector, 1 false event/min, 10-20 mm depth): strength for "
                     "50 % detection, dense OPM vs Neuromag combined: " + "; ".join(
                         f"{ANAT[lab]} {s50(lab, 'opm_dense/opm')} vs {s50(lab, 'squid/combined')} nAm" for lab in ["adult"] + labs)
                     + ". Details on the <a href=\"epilepsy.html#pediatric\">epilepsy page</a>.")
    return ("<h2 id=\"pediatric-findings\">Pediatric findings (G3B), with their conditions</h2>"
            + "".join(f'<div class="finding">{x}</div>' for x in items))


def page_benchmarks(d, out):
    g1a, g1b, g1c = d["g1a"], d["g1b"], d["g1c"]
    h = [f"<p>{label('REPRO')} G1A reproduces the analytical sphere benchmark of Jas et al. (2026). {label('ADAPT')} G1B and G1C "
         "adapt Hunold et al. (2016) and Goldenholz et al. (2009) to the MNE sample subject: the original anatomy, recordings and "
         "some details are unavailable. The OPM columns are new. Methods: sections 1, 6 and 7 of the "
         "<a href=\"methods.html\">methods</a>.</p>",
         "<h2 id=\"g1a\">G1A: Jas et al. 2026 analytical benchmark</h2>",
         fig(out, "g1a/fig3/Figure3_replicated.png",
             "Fig. 3 replica (&eta; = 3, adult sphere h = 95 mm, b = 80 mm, 30-nAm tangential dipole; OPM on the scalp, SQUID "
             "18 mm above it). The dotted line is drawn at 27.53 mm by a grid rule fitted to the published raster; the exact "
             f"Eq. 3 root is {g1a['d_eq_mm']['3']:.3f} mm.", "Replica of Jas et al. Figure 3"),
         sb.table(["Item", "Printed [mm]", "Drawn [mm]", "Exact [mm]"],
                  [[r["item"], r["printed_mm"], f"{r['drawn_mm']:.2f}", f"{r['exact_mm']:.3f}"] for r in g1a["printed_vs_exact"]
                   if "printed_mm" in r],
                  "Printed, drawn and exact equal-SNR depths"),
         f"<p>Checks: Eq. 1 agrees with an independent Sarvas-field maximisation to "
         f"{g1a['numerical_checks']['OPM']['max_rel_err_sarvas_2d']:.1e} and with MNE's sphere model to "
         f"{g1a['numerical_checks']['OPM']['max_rel_err_mne_sphere']:.1e}; a crossing exists only for "
         f"{g1a['eta_range_with_crossing'][0]:.4f} &lt; &eta; &lt; {g1a['eta_range_with_crossing'][1]:.4f}; "
         "the equal-SNR depth does not depend on the absolute noise level.</p>",
         fig(out, "g1a/Figure_G1A_fig4.png", "Fig. 4 reproduction: SNR vs depth for several noise ratios &eta; and the equal-SNR "
             "depth vs &eta; (exact Eq. 3 roots as markers).", "Reproduction of Jas et al. Figure 4"),
         fig(out, "g1a/Figure_G1A_toy_fig6.png", "Fig. 6 toy experiment (explanatory only): a target dipole between two noise "
             "dipoles; SNR as the ratio of peak fields vs sensor standoff.", "Toy experiment"),
         "<h2 id=\"g1b\">G1B: Hunold et al. 2016 depth-orientation spike SNR</h2>"]
    c = g1b["fig6_calibration"]
    h.append(f"<p>3,783 dipoles (600 nAm) and 20-mm&sup2; patches stratified to the paper's per-bin counts; background of "
             f"{g1b['n_background_dipoles']:,} random dipoles. Two background levels: as specified in the text, and calibrated with one "
             f"scalar ({c['scale']:.2f}) to the paper's Fig. 6 baselines, using our expected baseline over "
             f"{c['mag']['realizations']['n']} background realizations (other channel choices and single realizations give "
             f"{c['scale_range'][0]:.2f}-{c['scale_range'][1]:.2f}). The paper's maps are never used for calibration.</p>")
    rows = []
    for var in ("fig6_calibrated", "as_specified"):
        for fam in ("dipole", "patch"):
            for num in ("p2p", "noisy_p2p"):
                for k in ("mag", "grad"):
                    v = g1b["variants"][var]["comparison_with_paper"][f"{fam}/{num}"][k]
                    rows.append([var.replace("_", " "), fam, num.replace("_", " "), "MM" if k == "mag" else "GM", f"{v['pearson_r']:.2f}",
                                 f"{v['mean_ratio_ours_to_paper']:.2f}", f"{v['mean_ratio_strong_bins']:.2f}", f"{v['mean_ratio_weak_bins']:.2f}",
                                 f"{100 * v['threshold_2p5_agreement']:.0f} %"])
    h.append(sb.table(["Background", "Sources", "Numerator", "Sensor", "r", "Mean ratio ours/paper", "Strong bins", "Weak bins",
                       "2.5-threshold agreement"], rows, "Per-bin comparison with the digitised maps (Figs 4a and 5a)"))
    h.append(fig(out, "g1b/Figure_G1B_vs_paper.png", "Our bin means against the paper's digitised colour classes (noise-free "
                 "peak-to-peak numerator, calibrated background).", "G1B versus paper"))
    h.append(fig(out, "g1b/Figure_G1B_bins_fig6_calibrated_dipole_p2p.png", "Mean SNR per depth x orientation bin for "
                 "magnetometers (MM), gradiometers (GM) and the matched OPM array, with unpaired per-bin tests; dipoles, "
                 "calibrated background.", "G1B dipole bins"))
    h.append(fig(out, "g1b/Figure_G1B_bins_fig6_calibrated_patch_p2p.png", "The same for 20-mm&sup2; patches.", "G1B patch bins"))
    cp = g1c["comparison_with_paper"]
    pooled = g1c["distributions"]["focal/model/pooled"]
    inside = cp["probable_default/focal/model/pooled"]["share_m29_to_m19"]
    ext = g1c["extension_opm"]
    h.append("<h2 id=\"g1c\">G1C: Goldenholz et al. 2009 cortical SNR maps</h2>")
    h.append(f"<p>Eq. 1 SNR of 10-nAm dipoles at {g1c['n_usable_vertices']:,} usable vertices and of 10/16-mm patches at "
             f"{g1c['n_centroids']:,} centroids; modelled noise from {g1c['n_noise_sources']:,} cortical sources calibrated with the "
             f"paper's rule (s_s = {cp['source_sd_nAm']['probable_default']:.2f} nAm on our grid, "
             f"{cp['source_sd_nAm_at_4000_sources']['probable_default']:.2f} nAm at 4,000 sources; paper 1.6-1.9 nAm). Focal median "
             f"{pooled['median_db']:.1f} dB (5-95 %: {pooled['p5_db']:.1f} to {pooled['p95_db']:.1f} dB); {100 * inside:.0f} % of "
             "vertices fall inside the paper's -29 to -19 dB display range.</p>")
    h.append(fig(out, "g1c/Figure_G1C_maps_model.png", "Eq. 1 SNR with modelled brain noise (paper's colour limits).", "G1C modelled maps"))
    h.append(fig(out, "g1c/Figure_G1C_maps_recorded.png", "Eq. 1 SNR with recorded noise (task baselines of the sample "
                 "recording, SSP applied) - ADAPT.", "G1C recorded-noise maps"))
    h.append(fig(out, "g1c/Figure_G1C_maps_opm_extension.png", "NEW: matched OPM array vs Neuromag with brain noise (recorded "
                 "minus empty room) plus intrinsic sensor noise.", "G1C OPM extension maps"))
    rows = []
    for k, v in ext.items():  # e.g. "focal/opm@7_minus_mag", "brain_only/focal/opm_minus_grad"
        parts = k.split("/")
        if "_minus_" not in parts[-1] or "focal" not in parts or "share_opm_better" not in v:
            continue
        opm, ref = parts[-1].split("_minus_")
        noise = "brain noise only" if parts[0] == "brain_only" else f"OPM {opm.split('@')[1]} fT/sqrt(Hz)"
        rows.append([noise, LABEL.get(ref, ref), f"{v['median_db']:+.2f} dB", f"{100 * v['share_opm_better']:.0f} %"])
    if rows:
        h.append(sb.table(["Noise", "Matched OPM minus", "Median Eq. 1 difference", "Vertices where OPM is higher"], rows,
                          "OPM extension, focal dipoles"))
    return "\n".join(h)


def page_adult(d, out):
    g2, bands = d["g2"], d["bands"]
    P = g2["primary"]["oracle"]
    conds = g2["config"]["conditions"]["all"]
    rows = []
    for cond in conds:
        for a in ("opm_matched", "opm204", "opm_dense"):
            c = P[f"{a}/combined/{cond}"]
            rows.append([cond, LABEL[a], cmp_str(c), f"{100 * c['share_opm_better']:.0f} %", f"{100 * c.get('share_parcels_opm_better', 0):.0f} %",
                         cmp_str(P[f"{a}/grad/{cond}"]), cmp_str(P[f"{a}/mag/{cond}"])])
    h = [f"<p>{label('NEW')} Realistic adult comparison on the MNE sample subject: Neuromag at the measured head position vs "
         "single-axis OPM arrays (10-mm cell, 7-mm standoff). Detectability is the known-topography matched-filter SNR "
         "sqrt(s<sup>T</sup> C<sup>+</sup> s). The full generated report is <a href=\"g2-report.html\">here</a>; methods in "
         "section 8 of the <a href=\"methods.html\">methods</a>.</p>",
         fig(out, "g2/Figure_G2_arrays.png", "The four arrays on the sample head: Neuromag (306 channels at 102 sites), the matched "
             f"OPM array ({g2['arrays']['opm_matched']['channels']} of the 102 sites fit the OPM placement rules), the 204-channel "
             f"budget control and the dense {g2['arrays']['opm_dense']['channels']}-site array.", "Sensor arrays"),
         "<h2 id=\"headline\">Headline comparison</h2>",
         sb.table(["Noise condition", "OPM array", "vs Neuromag combined", "Targets OPM higher", "Parcels OPM higher", "vs gradiometers",
                   "vs magnetometers"], rows,
                  f"Median detectability ratio over {g2['n_targets']:,} targets, 95 % CI from a bootstrap over cortical parcels; OPM 15 fT/sqrt(Hz)"),
         fig(out, "g2/Figure_G2_depth.png", "Signal amplitude and detectability vs source depth (median and interquartile range over "
             "targets).", "Detectability versus depth"),
         fig(out, "g2/Figure_G2_heatmaps.png", "Median log2 detectability ratio per depth x orientation bin (at least 10 targets).",
             "Depth-orientation heatmaps"),
         fig(out, "g2/Figure_G2_maps_intrinsic_brain.png", "Cortical maps of the log2 detectability ratio with brain noise (red: OPM "
             "higher; black contour: ratio 1; grey: medial wall and vertices within 4 mm of the inner skull).", "Cortical ratio maps"),
         fig(out, "g2/Figure_G2_maps_projected.png", "The same after removing the 8-dimensional external-field subspace from every "
             "array.", "Cortical ratio maps, projected"),
         fig(out, "g2/Figure_G2_maps_absolute.png", "Absolute detectability of a 10-nAm dipole (log10) per array.", "Absolute maps"),
         fig(out, "g2/Figure_G2_heatmaps_comparators.png", "Dense OPM vs each Neuromag channel set, and extended sources.",
             "Comparator heatmaps"),
         fig(out, "g2/Figure_G2_patches.png", "Extended sources: geodesic patches of 5, 10 and 20 mm radius (fixed total moment and "
             "fixed density).", "Patches"),
         "<h2 id=\"sensitivity\">Sensitivity and model dependence</h2>",
         fig(out, "g2/Figure_G2_sensitivity.png", "One-factor-at-a-time sensitivity (OPM noise, correlated background, background "
             "calibration, head position, scalp gap) with parcel-bootstrap CIs; these show dependence, they do not bound it.",
             "Sensitivity analyses")]
    joint = g2["sensitivity_joint_asd_gap"]
    rows = [[k.replace("/", ", ").replace("gap", "gap ").replace("asd", "OPM noise "), x(v["opm_matched"]["ratio"]), x(v["opm_dense"]["ratio"])]
            for k, v in joint.items()]
    h.append(sb.table(["Scalp gap, OPM noise [fT/sqrt(Hz)]", "OPM matched", "OPM dense"], rows,
                      "Joint OPM noise x scalp-gap grid, vs Neuromag combined, brain noise"))
    bem = g2["convergence"]["bem"]
    h.append(f"<p>Head model: with a 1-layer instead of the 3-layer BEM the dense/Neuromag ratio is "
             f"{x(2 ** bem['bem1_5120']['median_log2']['opm_dense/combined/intrinsic+brain'])} "
             f"(3-layer on the same targets: {x(2 ** bem['subset_reference_median_log2']['opm_dense/combined/intrinsic+brain'])}); "
             "refining the 1-layer mesh changes it by at most "
             f"{bem['refinement_max_abs_change_log2']:.3f} in log2.</p>")
    if bands:
        h.append(fig(out, "g2/Figure_G2_bands.png", "Frequency bands (brain scale and room field recalibrated per band) and a 100-Hz "
                     "first-order OPM response.", "Frequency bands"))
    h.append(fig(out, "g2/Figure_G2_bridge.png", "Bridge to the analytical benchmark: the realistic OPM/magnetometer peak-field ratio "
                 "and the equal-SNR depth vs the sphere model with the real standoffs.", "Bridge to the sphere model"))
    return "\n".join(h)


def page_pediatric(d, out):
    g3a = d["g3a"]
    h = []
    if g3a:
        f = g3a["size_following"]
        rows = [[name, f"{v['h_mm']:.0f}/{v['b_mm']:.0f}", f"{v['eta3']['d_eq_exact_mm']:.1f}", f"{v['eta3']['normalized_exact_pct']:.1f} %",
                 f"{g3a['printed_checks'].get(name, {}).get('printed_pct', '-')}", f"{v['eta3']['volume_fraction_pct']:.1f} %"]
                for name, v in f.items()]
        h += [f"<h2 id=\"g3a\">G3A: Jas et al. head-size benchmark {label('REPRO')}</h2>",
              sb.table(["Head", "h/b [mm]", "d_eq at eta 3 [mm]", "Normalized d_eq", "Printed [%]", "Brain volume with OPM ahead"], rows,
                       "Size-following SQUID shell (h + 18 mm), exact Eq. 3 roots"),
              fig(out, "g3a/Figure_G3A_fig5.png", "Fig. 5 reproduction (A, B), the brain-volume fraction (C) and a NEW idealised "
                  "contrast with every head concentric in one adult shell (D).", "Head-size benchmark")]
    else:
        h.append("<p>G3A (the Jas et al. head-size benchmark) has not been run yet.</p>")
    g3b = d.get("g3b")
    if not g3b:
        h.append("<h2 id=\"g3b\">G3B: fixed adult helmet vs head-adaptive OPM</h2><p>G3B has not been run yet.</p>")
        return "\n".join(h)
    A = g3b["anatomies"]
    rows = [[ANAT[k], html.escape(a["scale_note"]), f"{a['head_size']['ofc_mm']:.0f}",
             f"{a['head_size']['breadth_mm']:.0f} x {a['head_size']['length_mm']:.0f}", f"{a['n_targets']:,}",
             f"{a['cortical_area_cm2']:.0f}", f"{g3b['arrays'][k]['opm_dense']['n']} / {g3b['arrays'][k]['opm_matched']['n']}",
             f"{g3b['placements'][k]['centred']['median_dist_mm']:.1f} / {g3b['placements'][k]['top']['median_dist_mm']:.1f}"]
            for k, a in A.items()]
    h += [f"<h2 id=\"g3b\">G3B: fixed adult helmet vs head-adaptive OPM {label('NEW')}</h2>",
          "<p>The same Neuromag helmet (sensors, coil types and intrinsic noise unchanged) holds the adult, two size-only controls "
          "(the adult scaled to school-age and to 2-year head size; every vertex homologous to the adult's), the 24-, 18- and "
          "12-month infant templates of O'Reilly et al. (2021) in their native dimensions and three school-aged children "
          "(OpenNeuro ds005234, Fadeev et al. 2024; their own cortex and scalp, a modelled skull, fiducials transferred from the adult), at placements chosen from the scalp and helmet geometry only (primary: raised to "
          "20-mm contact with the top of the helmet). The OPM arrays are refitted to each head with the adult rules (10-mm cell, "
          "17-mm packing, nothing shrunk). Brain-background variance per unit cortical area, room field and intrinsic noise are "
          "the adult's. D = OPM minus Neuromag known-topography detectability in dB; Delta = D<sub>child</sub> - "
          "D<sub>adult</sub> on homologous sources (vertex-wise for the scaled controls; parcels and depth strata for the "
          "templates and the school-aged children). A positive Delta is a relative gain for OPM, not by itself an OPM advantage in the child. Methods: section "
          "10 of the <a href=\"methods.html\">methods</a>; full tables in the <a href=\"g3b-report.html\">G3B report</a>.</p>",
          sb.table(["Anatomy", "Construction", "Head circumference [mm]", "Breadth x length [mm]", "Targets", "Usable cortex [cm2]",
                    "OPM dense / matched sites", "Neuromag gap centred / top [mm]"], rows,
                   "Anatomies and arrays (gap: median magnetometer coil centre to scalp)"),
          fig(out, "g3b/Figure_G3B_geometry.png", "Sagittal and coronal scalp sections at the centred (grey) and top-contact (black) "
              "placements in the fixed helmet, the counterfactual helmet scaled with the head (green) and the refitted dense OPM "
              "sites (blue).", "Pediatric geometry")]
    rows = []
    for c in [c for c in CHILDREN if f"{c}/opm_dense/combined/intrinsic+brain/detect" in g3b["comparisons"]]:
        for ref in ("combined", "grad", "mag"):
            r = g3b["comparisons"][f"{c}/opm_dense/{ref}/intrinsic+brain/detect"]

            def f(s):
                ci = s.get("ci95")
                return f"{s['median']:+.2f}" + (f" [{ci[0]:+.2f}, {ci[1]:+.2f}]" if ci else "")

            rows.append([ANAT[c], "Neuromag " + LABEL[ref], f(r["d_child"]), f(r["d_adult"]), f(r["delta"])])
    h += [sb.table(["Child anatomy", "Comparator", "D child [dB]", "D adult [dB]", "Delta [dB]"], rows,
                   "Dense OPM vs Neuromag, intrinsic + brain noise, primary placement: area-weighted medians with parcel-bootstrap "
                   "95 % intervals"),
          fig(out, "g3b/Figure_G3B_depth.png", "D vs depth below the scalp in each anatomy (solid: dense OPM, dashed: matched).",
              "D vs depth"),
          fig(out, "g3b/Figure_G3B_delta.png", "Delta per depth stratum (child minus adult) for each comparator.", "Delta vs depth"),
          fig(out, "g3b/Figure_G3B_placements.png", "D for every source-blind placement and the counterfactual helmet (left) and "
              "the regional magnetometer-to-scalp gaps (right).", "Placements"),
          fig(out, "g3b/Figure_G3B_maps_scaled.png", "D on the adult's inflated cortex (adult and scaled controls).", "Maps, scaled"),
          fig(out, "g3b/Figure_G3B_maps_delta.png", "Delta at every vertex (scaled controls).", "Maps, Delta"),
          *[fig(out, f"g3b/Figure_G3B_maps_{k}.png", f"D on the {ANAT[k]}'s inflated cortex (derived from O'Reilly et al. "
                "2021 / Richards et al. 2016).", f"Maps, {ANAT[k]}") for k in TEMPLATES if (RES / "g3b" / f"Figure_G3B_maps_{k}.png").exists()],
          fig(out, "g3b/Figure_G3B_usefulness_adult.png", "Where the dense OPM array, Neuromag, both or neither reach "
              "detectability 5 for a 100-nAm source (operational threshold), adult.", "Usefulness, adult"),
          *[fig(out, f"g3b/Figure_G3B_usefulness_{k}.png", f"The same for the {ANAT[k]}.", f"Usefulness, {ANAT[k]}")
            for k in TEMPLATES if (RES / "g3b" / f"Figure_G3B_usefulness_{k}.png").exists()]]
    return "\n".join(h)


def page_epilepsy(d, out):
    g4, loc = d["g4"], d["loc"]
    rows = []
    for mode in ("practical@1", "oracle"):
        for key, name in DETECTORS.items():
            s = g4["detectors"][key]["strength_for_50pct_nAm"]
            cells = []
            for b in range(4):
                v = s[f"{mode}/depth{b}"]
                if v["value"] is None:
                    cells.append("not reached")
                    continue
                lo, hi = v["ci95"]
                cells.append(f"{v['value']:.0f} [{'-' if lo is None else f'{lo:.0f}'}-{'>320' if hi is None else f'{hi:.0f}'}]")
            rows.append([mode.replace("practical@1", "practical, 1 false event/min"), name, *cells])
    h = [f"<p>{label('NEW')} Interictal-spike-like events (Hunold's spike-wave waveform, three morphologies, 10-320 nAm focal "
         "dipoles and 10-mm patches) in the G2 noise model in the time domain, identical in every array. An oracle knows source "
         "and time; the practical detector scans time and a cortical dictionary with thresholds frozen on independent null data. "
         "Methods: section 9 of the <a href=\"methods.html\">methods</a>.</p>",
         "<h2 id=\"detection\">Detection</h2>",
         fig(out, "g4/Figure_G4_detection.png", "Detection probability vs strength per depth band (focal, three morphologies pooled; "
             "Wilson bands are event-level and descriptive).", "Detection curves"),
         sb.table(["Detector", "Array", *DEPTH_BANDS], rows,
                  "Strength for 50 % detection [nAm], 95 % interval from a bootstrap over locations (>320: beyond the tested range)"),
         fig(out, "g4/Figure_G4_roc.png", "Practical detector: sensitivity vs false events per minute.", "Sensitivity vs false events")]
    rows = []
    for a in ("opm_dense", "opm_matched"):
        for mode in ("practical@1", "oracle"):
            p = g4["paired"][f"{a}/opm_vs_squid/combined/{mode}"]
            rows.append([LABEL[a], mode.replace("practical@1", "practical")] + [
                f"{p[f'depth{b}']['locations_favouring_opm']}/{p[f'depth{b}']['locations_favouring_squid']} (p {p[f'depth{b}']['location_sign_flip_p']:.2g})"
                for b in range(4)])
    h.append(sb.table(["OPM array", "Detector", *DEPTH_BANDS], rows, "Paired with Neuromag combined on identical events: locations "
                      "with more events detected only by OPM / only by Neuromag, exact sign-flip p (18 locations per band)"))
    rows = []
    for k, v in loc["results"].items():
        name, fam, s = k.split("/")

        def f(val):
            return "-" if val is None else f"{val:.1f}"

        rows.append(["Neuromag combined" if name == "squid" else LABEL[name], fam, s.replace("nAm", " nAm"), f"{100 * v['detected']:.0f} %",
                     f(v["ecd_error_mm_median_detected"]),
                     f(v["dspm_error_mm_median_detected"]), f(v.get("dspm_mne_error_mm_median_detected")),
                     f"{100 * v['joint_detect_and_ecd_within_10mm']:.0f} %",
                     f"{100 * v['joint_detect_and_dspm_within_10mm']:.0f} %"])
    held = (loc.get("thresholds_heldout") or {}).get("false_events") or {}
    held_txt = ("; ".join(f"{'Neuromag combined' if k == 'squid' else LABEL.get(k, k)} {v['rate_per_min']:.2f} "
                          f"[{v['ci95'][0]:.2f}-{v['ci95'][1]:.2f}]" for k, v in held.items()))
    h += ["<h2 id=\"localization\">Bounded localization</h2>",
          f"<p>24 locations; inverse with a 1-layer BEM and a 2-mm/2-deg coregistration error ({loc['config']['coreg_draws']} draws "
          f"shared by all arrays; median displacement at the source {loc['coreg_displacement_mm_median']:.1f} mm). Errors are "
          "measured on the MRI through the analyst's transform. No goodness-of-fit cut: whitened GOF rises as the channel count "
          "falls. Neuromag is also localized with one sensor type (its own covariance, detector and inverse; secondary), and "
          "dSPM is also computed with MNE's own <code>mne.minimum_norm</code> (it differs from the study's implementation in the "
          "depth weighting).</p>"
          + (f"<p>Detector thresholds (1 false event per minute, 10 min of null data) on {loc['thresholds_heldout']['minutes']:g} "
             f"min of independent null data: false events per minute [exact 95 % interval]: {held_txt}.</p>" if held else ""),
          fig(out, "g4/Figure_G4_localization.png", "Localization error on the MRI (all events) for dSPM and the equivalent current "
              "dipole, and the share of events both detected and localized within 10 mm.", "Localization errors"),
          sb.table(["Array", "Source", "Strength", "Detected", "ECD error, detected [mm]", "dSPM error, detected [mm]",
                    "dSPM (MNE), detected [mm]", "Detected + ECD <= 10 mm", "Detected + dSPM <= 10 mm"], rows,
                   "Medians per condition (one event per location)")]
    h.append(page_epilepsy_pediatric(d, out))
    h.append(page_motion(d, out))
    return "\n".join(h)


def page_epilepsy_pediatric(d, out):
    g4p = d.get("g4p")
    if not g4p:
        return "<h2 id=\"pediatric\">Pediatric</h2><p>The pediatric epilepsy runs have not been done yet.</p>"
    C, labs = g4p["comparison"], g4p["labels"]

    def s50(lab, key, b):
        v = C.get(f"{lab}/{key}/practical@1/depth{b}/s50")
        if v is None:
            return "no locations"
        if v["value"] is None:
            return "not reached"
        lo, hi = v["ci95"]
        return f"{v['value']:.0f} [{'-' if lo is None else f'{lo:.0f}'}-{'>320' if hi is None else f'{hi:.0f}'}]"

    rows = [[ANAT[lab], name, *[s50(lab, key, b) for b in range(4)]]
            for lab in labs for key, name in (("opm_dense/opm", "OPM dense"), ("squid/combined", "Neuromag combined"),
                                               ("opm_matched/opm", "OPM matched"))]
    prow = []
    for lab in labs:
        cells = []
        for b in range(4):
            r = C.get(f"{lab}/paired/opm_dense/opm_vs_squid/combined/practical@1/depth{b}")
            if r is None:
                cells.append("no locations")
                continue
            cells.append(f"{r['locations_favouring_opm']}/{r['locations_favouring_squid']} (p {r['location_sign_flip_p']:.2g}); "
                         f"ratio {detection.format_s50_ratio(r['s50_ratio_squid_over_opm'])}")
        prow.append([ANAT[lab], *cells])
    h = [f"<h2 id=\"pediatric\">Pediatric: the same framework on smaller heads in the fixed helmet {label('NEW')}</h2>",
         "<p>The adult detection and localization studies rerun unchanged (configuration, seeds, detectors, operating points) on "
         "the G3B anatomies with Neuromag at the primary placement (top contact) and the refitted OPM arrays; thresholds are "
         "calibrated on each anatomy's own null data. The templates are averages, the school-age and 2-year-size heads scaled "
         "adults and children A-C three individuals of one dataset; none is a population. Full tables: <code>results/g4/G4_pediatric_report.md</code>.</p>",
         sb.table(["Anatomy", "Array", *DEPTH_BANDS], rows, "Strength for 50 % detection [nAm], practical detector at 1 false "
                  "event/min, 95 % interval from a bootstrap over locations"),
         sb.table(["Anatomy", *DEPTH_BANDS], prow, "Dense OPM vs Neuromag combined on identical events: locations favouring OPM / "
                  "Neuromag (exact sign-flip p, uncorrected) and the paired strength ratio Neuromag/OPM [95 % CI]")]
    rates = [C[k]["1"] for k in C if k.endswith("/heldout_false_per_min") and "1" in C[k]]
    h.append(f"<p>Operating points: the frozen 1-per-minute thresholds give {min(rates):.2f}-{max(rates):.2f} false events per "
             "minute on each anatomy's held-out null data (event-level rates, distinct from the oracle's per-trial false-positive "
             "probability). With every detector set to 1 per minute on the held-out null the superficial result is unchanged "
             "(<code>results/g4/G4_matched_rate_report.md</code>, in the downloads). The 50 % detection level and the 1-per-minute "
             "operating point are operational study choices, not clinical standards. Failed dipole fits under declared criteria: "
             "<code>results/g4/G4_fit_failures_report.md</code>. Simulated IED-source recovery does not identify an epileptogenic "
             "zone or establish surgical benefit.</p>")
    for lab in labs:
        if lab == "adult":
            continue
        for kind, cap in (("detection", "detection probability vs strength"), ("roc", "sensitivity vs false events per minute"),
                          ("localization", "localization errors")):
            png = RES / "g4" / f"Figure_G4_{kind}_{lab}.png"
            if png.exists():
                h.append(fig(out, f"g4/Figure_G4_{kind}_{lab}.png", f"{ANAT[lab]}: {cap}.", f"{ANAT[lab]} {kind}"))
    return "\n".join(h)


def page_motion(d, out):
    m = d.get("motion")
    if not m:
        return "<h2 id=\"motion\">Head motion and OPM slippage</h2><p>The bounded motion extension has not been run yet.</p>"
    keys = m["anatomies"]
    pick_sq = ("down 5 mm", "down 10 mm", "x+5 mm", "y-5 mm", "pitch +10 deg", "pitch -10 deg", "yaw +10 deg")
    pick_op = ("slip x +3 deg", "slip y +3 deg", "slip z +3 deg")
    rows = []
    for k in keys:
        for sysname, names in (("squid", pick_sq), ("opm", pick_op)):
            for n in names:
                r = m["geometry"][k][sysname].get(n)
                if r is None:
                    continue
                if not r["feasible"]:
                    rows.append([ANAT[k], "Neuromag combined", n, "infeasible", "", f"nearest magnetometer {r['min_dist_mm']:.1f} mm"])
                    continue
                note = (f"nearest magnetometer {r['min_dist_mm']:.1f} mm" if sysname == "squid"
                        else f"sensors moved {r['median_sensor_shift_mm']:.1f} mm (median)")
                rows.append([ANAT[k], "Neuromag combined" if sysname == "squid" else "OPM dense", n, f"{r['known']['median_db']:+.2f}",
                             f"{r['mismatched']['median_db']:+.2f}", note])
    brows = []
    for k in keys:
        for corr in ("none", "homogeneous", "homogeneous+gradient"):
            for cal in m["calibration_labels"]:
                for field in ("uniform", "gradient"):
                    c = m["coupling"][k]["cases"].get(f"{corr}/{cal}/{field}/neck")
                    if c is None:
                        continue
                    be = c["unmodelled"]["thresholds_deg_unit_field"]

                    def f(r):
                        v, q = r["median_curve"], r.get("per_draw_p10_p90")
                        txt = "not reached" if v is None else f"{v:.3g}"
                        if q:
                            lo, hi = (f"beyond {m['config']['coupling']['rotation_rms_deg'][-1]:g}" if x is None else f"{x:.3g}"
                                      for x in q)
                            txt += f" [{lo}-{hi}]"
                        return txt

                    brows.append([ANAT[k], corr, cal.replace("tilt", "").replace("deg_gain", " deg, ").replace("pct", " %"),
                                  "1 nT" if field == "uniform" else "1 nT/m", f(be["loss_1dB"]), f(be["D_0dB"]),
                                  f"{c['oracle']['opm_change_db_median'][-1]:+.2f}"])
    tc = m["timecourse"]
    h = [f"<h2 id=\"motion\">Head motion and OPM slippage (bounded extension) {label('NEW')}</h2>",
         "<p>Static fit is studied first; this extension bounds the two motion mechanisms the static maps cannot show, on the G3B "
         "arrays and noise conventions (intrinsic + brain noise). <strong>A.</strong> A sustained displacement moves the head "
         "inside the fixed Neuromag helmet, while a head-mounted OPM array moves with the head and changes geometry only if the "
         "cap slips. 'Known': the displaced geometry is used (the ideal limit of movement compensation); 'mismatched': the "
         "template of the reference geometry is applied to the displaced data. <strong>B.</strong> In-band head rotation moves "
         "the OPM array through the static residual field of the room; the SQUIDs are fixed and see no such term. Results are "
         "per unit field and scale linearly with rotation x field. It does not establish motion robustness; sensor dynamic "
         "range, gain changes and real head-motion statistics are not modelled. Methods: section 12 of the "
         "<a href=\"methods.html\">methods</a>; all tables in the <a href=\"motion-report.html\">motion report</a>.</p>",
         fig(out, "g4/Figure_G4_motion.png", "A: median change in detectability for sustained displacements (dots: template of the "
             "reference geometry; bars: geometry known). B: OPM detectability vs in-band rotation in a 1-nT uniform field and a "
             "1-nT/m gradient, when the artefact is not part of the noise model (red: no correction; blue: homogeneous-field "
             "projection; green: 8-term projection; lighter: larger calibration errors). B': exact rigid motion over a recording.",
             "Motion and slippage"),
         sb.table(["Anatomy", "System", "Displacement", "Known [dB]", "Mismatched [dB]", "Note"], rows,
                  "A. Sustained displacement: median change in detectability (dense OPM or Neuromag combined), cortical targets"),
         sb.table(["Anatomy", "Correction", "Calibration error (axis tilt, gain)", "Unit field", "1-dB loss at [deg RMS]",
                   "OPM no longer ahead at [deg RMS]", "Oracle loss at 5 deg [dB]"], brows,
                  "B. In-band rotation (per axis, RMS, about a pivot 60 mm below the head origin) in a unit field at which the "
                  "dense OPM array loses 1 dB, or falls to Neuromag's static detectability, with the artefact outside the noise "
                  "model (median curve over draws; 10th-90th percentiles of the per-draw thresholds); divide by the field in nT "
                  "(or nT/m) for another room. Last column: artefact inside the noise model."),
         f"<p>B'. Exact rigid motion over {tc['config']['duration_s']:g} s ({ANAT[tc['config']['anatomy']]}; drift up to "
         f"{tc['config']['drift_deg']:g} deg, in-band jitter {tc['config']['inband_rotation_rms_deg']:g} deg RMS per axis, "
         f"{tc['config']['b0_nT']:g} nT and {tc['config']['gradient_nT_per_m']:g} nT/m): peak field change at a sensor "
         f"{tc['peak_field_change_pT']['median']:.0f} pT (median), in-band residual after the 8-term projection "
         f"{tc['corrections']['homogeneous+gradient']['inband_rms_fT_median']:.0f} fT RMS; per channel the exact in-band artefact "
         "is within " + f"{max(max(abs(r['exact_over_linear_percentiles']['p5'] - 1), abs(r['exact_over_linear_percentiles']['p95'] - 1)) for r in tc['corrections'].values()) * 100:.1f}"
         + " % of the linear prediction for 90 % of the channels (largest deviation "
         + f"{max(r['exact_over_linear_max_abs_deviation'] for r in tc['corrections'].values()) * 100:.1f} %).</p>"]
    return "\n".join(h)


def page_reproduce(out, manifest):
    readme = (ROOT / "README.md").read_text()
    start = readme.find("## Reproduce")
    end = readme.find("\n## ", start + 5)
    h = [sb.md_to_html(readme[start:end], heading_offset=0) if start >= 0 else "",
         "<h3 id=\"run-all\">scripts/run_all.sh</h3>",
         "<pre><code>" + html.escape((ROOT / "scripts" / "run_all.sh").read_text(), quote=False) + "</code></pre>",
         "<h2 id=\"downloads\">Downloads</h2>",
         "<p>Result files as committed. The commit is the one recorded in the result file; tables and reports carry the commit "
         "of the run that wrote them. The full SHA-256 of every file is in <a href=\"data/MANIFEST.json\">data/MANIFEST.json</a>.</p>",
         sb.table(["File", "Size", "SHA-256 (first 16)", "Code commit"],
                  [[f'<a href="{html.escape(m["href"])}">{html.escape(m["path"])}</a>', m["size"], m["sha256"][:16],
                    m["commit"] + (" (historic diagnosis, not reproducible: see the methods)" if m["commit"].endswith("+dirty") else "")]
                   for m in manifest], cls="downloads", html_cols=(0,))]
    return "\n".join(h)


# tables and reports written by the same run as a summary JSON (same directory)
WRITTEN_BY = {"g1a_curves.csv": "g1a_benchmark.json", "g1b_sources.csv": "g1b_summary.json", "g1c_oct6_values.csv": "g1c_summary.json",
              "g2_targets.csv": "g2_summary.json", "g2_patch_targets.csv": "g2_summary.json", "G2_report.md": "g2_summary.json",
              "g3a_deq.csv": "g3a_size_benchmark.json", "g4_adult_events.csv": "g4_adult_summary.json",
              "g4_localization_events.csv": "g4_localization_summary.json", "G3B_report.md": "g3b_summary.json",
              "G4_pediatric_report.md": "g4_pediatric_comparison.json", "G4_motion_report.md": "g4_motion_summary.json",
              "G4_matched_rate_report.md": "g4_matched_rate.json", "g4_motion_timecourse_example.csv": "g4_motion_summary.json",
              "G4_fit_failures_report.md": "g4_fit_failures.json"}
# result folders whose files carry no provenance of their own: the run that wrote them
PARENT_RUN = {"g1a/fig3": "g1a/g1a_benchmark.json"}
for _k in ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo"):
    WRITTEN_BY[f"g3b_targets_{_k}.csv"] = "g3b_summary.json"
    WRITTEN_BY[f"g4_{_k}_events.csv"] = f"g4_{_k}_summary.json"
    WRITTEN_BY[f"g4_localization_{_k}_events.csv"] = f"g4_localization_{_k}_summary.json"


def build_manifest(out):
    """Copy JSON, CSV and Markdown results for download; record size, SHA-256 and commit."""
    manifest = []
    for p in sorted(RES.rglob("*")):
        if p.suffix not in (".json", ".csv", ".md") or p.name == "README.md":
            continue
        rel = p.relative_to(RES)
        dst = out / "data" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dst)
        src = p if p.suffix == ".json" else p.parent / WRITTEN_BY.get(p.name, "")
        commit = "-"
        if src.suffix == ".json" and src.exists():
            commit = (json.loads(src.read_text()).get("provenance") or {}).get("commit", "-")
        parent = PARENT_RUN.get(str(rel.parent))
        if commit == "-" and parent:
            commit = (json.loads((RES / parent).read_text()).get("provenance") or {}).get("commit", "-")
        manifest.append(dict(path=str(rel), href=f"data/{rel}", size=f"{p.stat().st_size / 1024:.0f} kB", sha256=sb.sha256(p),
                             commit=commit))
    (out / "data" / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(ROOT / "site" / "_build"))
    args = ap.parse_args(argv)
    out = Path(args.out).resolve()
    if out.exists():
        # delete only a directory this builder made (its marker file) or an empty one, never a source tree
        if not (out / MARKER).is_file() and any(out.iterdir()):
            raise SystemExit(f"{out} exists and was not made by build_site.py; choose an empty or new directory")
        shutil.rmtree(out)
    (out / "static").mkdir(parents=True)
    (out / MARKER).write_text("generated by scripts/build_site.py; safe to delete\n")
    shutil.copy2(ROOT / "site" / "static" / "style.css", out / "static" / "style.css")
    d = dict(g1a=load("g1a/g1a_benchmark.json"), g1b=load("g1b/g1b_summary.json"), g1c=load("g1c/g1c_summary.json"),
             g2=load("g2/g2_summary.json"), bands=load("g2/g2_band_sensitivity.json"), g3a=load("g3a/g3a_size_benchmark.json"),
             g4=load("g4/g4_adult_summary.json"), loc=load("g4/g4_localization_summary.json"), g3b=load("g3b/g3b_summary.json"),
             g4p=load("g4/g4_pediatric_comparison.json"), motion=load("g4/g4_motion_summary.json"))
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(ROOT / "site" / "templates"), autoescape=True)
    tpl = env.get_template("base.html")
    head = git("rev-parse", "--short", "HEAD") + ("+dirty" if git("status", "--porcelain", "--", *io.CODE_PATHS, "site") else "")
    head_date = git("log", "-1", "--format=%cs")
    nav = [dict(file=f, title=t) for f, t in NAV]
    manifest = build_manifest(out)
    methods = sb.md_to_html((ROOT / "docs" / "methods.md").read_text(), heading_offset=1)
    pages = {
        "index.html": ("OPM vs SQUID MEG: adult baseline and pediatric extension", page_index(d)),
        "benchmarks.html": ("Adult reference benchmarks (G1)", page_benchmarks(d, out)),
        "adult.html": ("Realistic adult comparison (G2)", page_adult(d, out)),
        "g2-report.html": ("G2 report", sb.md_to_html((RES / "g2" / "G2_report.md").read_text(), heading_offset=1)),
        "pediatric.html": ("Pediatric extension (G3)", page_pediatric(d, out)),
        "g3b-report.html": ("G3B report", sb.md_to_html((RES / "g3b" / "G3B_report.md").read_text(), heading_offset=1)
                            if (RES / "g3b" / "G3B_report.md").exists() else "<p>G3B has not been run yet.</p>"),
        "epilepsy.html": ("Epilepsy relevance (G4)", page_epilepsy(d, out)),
        "motion-report.html": ("G4 motion and slippage report", sb.md_to_html((RES / "g4" / "G4_motion_report.md").read_text(), heading_offset=1)
                               if (RES / "g4" / "G4_motion_report.md").exists() else "<p>The motion extension has not been run yet.</p>"),
        "methods.html": ("Methods, uncertainty and limitations", methods),
        "register.html": ("Parameters, provenance and assumptions",
                          sb.md_to_html((ROOT / "docs" / "provenance_register.md").read_text(), heading_offset=1)),
        "reproduce.html": ("Reproduce and download", page_reproduce(out, manifest)),
    }
    for fname, (title, content) in pages.items():
        current = {"g2-report.html": "adult.html", "g3b-report.html": "pediatric.html", "motion-report.html": "epilepsy.html"}.get(fname, fname)
        (out / fname).write_text(tpl.render(title=title, content=Markup(content), nav=nav, current=current,
                                            head=head, head_date=head_date))
    problems = sb.check_links(out)
    if problems:
        print("\n".join(problems))
        raise SystemExit(f"{len(problems)} broken links")
    print(f"built {len(pages)} pages, {len(manifest)} downloads, {len(list((out / 'figures').rglob('*.png')))} figures -> {out}")


if __name__ == "__main__":
    main()
