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

from opmsquid import io, sitebuild as sb  # noqa: E402

RES = ROOT / "results"
NAV = [("index.html", "Overview"), ("benchmarks.html", "Adult benchmarks (G1)"), ("adult.html", "Realistic adult (G2)"),
       ("pediatric.html", "Pediatric (G3)"), ("epilepsy.html", "Epilepsy (G4)"), ("methods.html", "Methods and limitations"),
       ("register.html", "Parameters and provenance"), ("reproduce.html", "Reproduce and download")]
FROZEN_TAG = "adult-baseline-v2"
LABEL = {"squid": "Neuromag", "opm_matched": "OPM matched", "opm204": "OPM 204 (channel budget)",
         "opm_dense": "OPM dense", "combined": "combined", "grad": "gradiometers", "mag": "magnetometers"}
DETECTORS = {"squid/combined": "Neuromag combined", "squid/grad": "Neuromag gradiometers", "squid/mag": "Neuromag magnetometers",
             "opm_matched/opm": "OPM matched", "opm_dense/opm": "OPM dense"}
DEPTH_BANDS = ("10-20 mm", "20-30 mm", "30-45 mm", "45-70 mm")


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
              ("G1A Jas et al. 2026 analytical benchmark", "REPRO", "done, independently reviewed"),
              ("G1B Hunold et al. 2016 depth-orientation spike SNR (MEG)", "ADAPT (+ NEW OPM column)", "done, independently reviewed"),
              ("G1C Goldenholz et al. 2009 cortical SNR maps (MEG)", "ADAPT (+ NEW OPM extension)", "done, independently reviewed"),
              ("G2 realistic adult OPM vs Neuromag", "NEW", f"done, independently reviewed; frozen as {FROZEN_TAG}"),
              ("G3A Jas head-size benchmark", "REPRO (+ NEW fixed shell)", "done" if d["g3a"] else "not run"),
              ("G3B pediatric fixed helmet vs head-adaptive OPM", "NEW",
               "in progress (2-year infant template and scaled-adult size controls; owner decision 2026-09-30)"),
              ("G4 epilepsy detection and bounded localization", "NEW", "adult done; pediatric waits for G3B"),
              ("G5 software, reproduction, report", "-", "in progress (this local report; clean-environment smoke test passed)")]
    h = ["<p>Simulation study comparing on-scalp optically pumped magnetometers (OPM) with the Neuromag SQUID system: "
         "an analytical benchmark and adaptations of two published adult studies, a realistic adult comparison, and "
         "(planned) pediatric and epilepsy extensions. The OPM advantage is tested, not assumed. Labels: "
         f"{label('REPRO')} reproduction with the paper's definitions, {label('ADAPT')} adaptation where data or details are "
         f"unavailable, {label('NEW')} new experiment.</p>",
         "<h2 id=\"status\">Status</h2>", sb.table(["Milestone", "Label", "Status"], status),
         "<h2 id=\"findings\">Adult findings, with their conditions</h2>"]
    n_dense, n_matched = g2["arrays"]["opm_dense"]["channels"], g2["arrays"]["opm_matched"]["channels"]

    def bands(p):
        return "; ".join(f"{DEPTH_BANDS[b]} {p[f'depth{b}']['locations_favouring_opm']}/{p[f'depth{b}']['locations_favouring_squid']} "
                         f"(p = {p[f'depth{b}']['location_sign_flip_p']:.2g})" for b in range(4))

    r0 = pr["depth0"].get("s50_ratio_squid_over_opm")
    ratio_txt = (f" (Neuromag / dense strength ratio {r0['value']:.2f}, paired location-bootstrap 95 % interval "
                 f"{r0['ci95'][0]:.2f}-{r0['ci95'][1]:.2f})") if r0 and r0.get("value") and r0.get("ci95") else ""
    pl = loc["paired"]
    dspm = {a: pl[f"{a}_vs_squid/patch/320nAm"]["dspm_error_mm"] for a in ("opm_matched", "opm_dense")}
    items = [
        f"With modelled brain noise, the dense {n_dense}-site OPM array has <strong>{cmp_str(dense)}</strong> the known-topography "
        f"detectability of the Neuromag system (all 306 channels; median over {g2['n_targets']:,} cortical targets, 95 % CI "
        f"from a bootstrap over cortical parcels; higher for {100 * dense['share_opm_better']:.0f} % of targets). An OPM array at "
        f"Neuromag's own sites ({n_matched} of 102 fit) ties it: {cmp_str(matched)}. Conditions: OPM white noise 15 fT/&radic;Hz, "
        "3-layer BEM, no extra scalp gap, the sample subject's measured head position; the modelled brain noise predicts 0.73x "
        "the measured magnetometer level.",
        f"The advantage depends on the assumptions: OPM noise 7 or 30 fT/&radic;Hz gives "
        f"{x(2 ** sens['opm_asd_7fT']['opm_dense/combined/intrinsic+brain']['median_log2'])} or "
        f"{x(2 ** sens['opm_asd_30fT']['opm_dense/combined/intrinsic+brain']['median_log2'])}; with 30 fT/&radic;Hz and a 6-mm "
        f"scalp gap together, {x(joint['gap6mm/asd30fT']['opm_dense']['ratio'])}; a 1-layer head model gives {x(bem1)}. Without "
        f"brain noise (intrinsic sensor noise only) Neuromag wins: {cmp_str(P['opm_dense/combined/intrinsic'])}.",
        f"By depth (dense vs Neuromag): {x(depth['10-15'])} at 10-15 mm below the scalp, {x(depth['30-35'])} at 30-35 mm and "
        f"{x(depth['50-55'])} at 50-55 mm; the matched array falls from "
        f"{x(2 ** g2['log2_ratio_vs_depth']['opm_matched/combined/intrinsic+brain'][0]['median'])} to below 1 at depth.",
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
        "Neuromag (uncorrected; 16 comparisons per array).",
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
             "<li>OPM intrinsic noise is a declared sweep (7-30 fT/&radic;Hz), not a device specification; OPM movement "
             "artefacts, cross-talk and calibration errors are not modelled.</li>"
             "<li>Detectability is a known-topography matched-filter SNR, not a clinical detection rate; the spike study "
             "uses simulated events in simulated noise.</li>"
             "<li>Confidence intervals resample cortical parcels or locations of one anatomy; they do not include model "
             "uncertainty, which the sensitivity analyses show instead (one factor at a time, plus a joint noise x gap grid).</li>"
             "<li>No pediatric result in this adult baseline (G3B follows it).</li></ul>")
    return "\n".join(h)


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
    h.append("<h2 id=\"g3b\">G3B: fixed adult helmet vs head-adaptive OPM</h2>"
             "<p>In progress, after the adult baseline: the 2-year infant template (O'Reilly et al. 2021) in its native "
             "dimensions and the adult scaled to school-age and 2-year head size as size-only controls (owner decision and "
             "download approval 2026-09-30).</p>")
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
                     f(v["dspm_error_mm_median_detected"]), f"{100 * v['joint_detect_and_ecd_within_10mm']:.0f} %",
                     f"{100 * v['joint_detect_and_dspm_within_10mm']:.0f} %"])
    h += ["<h2 id=\"localization\">Bounded localization</h2>",
          f"<p>24 locations; inverse with a 1-layer BEM and a 2-mm/2-deg coregistration error ({loc['config']['coreg_draws']} draws "
          f"shared by all arrays; median displacement at the source {loc['coreg_displacement_mm_median']:.1f} mm). Errors are "
          "measured on the MRI through the analyst's transform. No goodness-of-fit cut: whitened GOF rises as the channel count "
          "falls.</p>",
          fig(out, "g4/Figure_G4_localization.png", "Localization error on the MRI (all events) for dSPM and the equivalent current "
              "dipole, and the share of events both detected and localized within 10 mm.", "Localization errors"),
          sb.table(["Array", "Source", "Strength", "Detected", "ECD error, detected [mm]", "dSPM error, detected [mm]",
                    "Detected + ECD <= 10 mm", "Detected + dSPM <= 10 mm"], rows, "Medians per condition (one event per location)")]
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
                  [[f'<a href="{html.escape(m["href"])}">{html.escape(m["path"])}</a>', m["size"], m["sha256"][:16], m["commit"]]
                   for m in manifest], cls="downloads", html_cols=(0,))]
    return "\n".join(h)


# tables and reports written by the same run as a summary JSON (same directory)
WRITTEN_BY = {"g1a_curves.csv": "g1a_benchmark.json", "g1b_sources.csv": "g1b_summary.json", "g1c_oct6_values.csv": "g1c_summary.json",
              "g2_targets.csv": "g2_summary.json", "g2_patch_targets.csv": "g2_summary.json", "G2_report.md": "g2_summary.json",
              "g3a_deq.csv": "g3a_size_benchmark.json", "g4_adult_events.csv": "g4_adult_summary.json",
              "g4_localization_events.csv": "g4_localization_summary.json"}


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
        manifest.append(dict(path=str(rel), href=f"data/{rel}", size=f"{p.stat().st_size / 1024:.0f} kB", sha256=sb.sha256(p),
                             commit=commit))
    (out / "data" / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(ROOT / "site" / "_build"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    (out / "static").mkdir(parents=True)
    shutil.copy2(ROOT / "site" / "static" / "style.css", out / "static" / "style.css")
    d = dict(g1a=load("g1a/g1a_benchmark.json"), g1b=load("g1b/g1b_summary.json"), g1c=load("g1c/g1c_summary.json"),
             g2=load("g2/g2_summary.json"), bands=load("g2/g2_band_sensitivity.json"), g3a=load("g3a/g3a_size_benchmark.json"),
             g4=load("g4/g4_adult_summary.json"), loc=load("g4/g4_localization_summary.json"))
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(ROOT / "site" / "templates"), autoescape=True)
    tpl = env.get_template("base.html")
    head = git("rev-parse", "--short", "HEAD") + ("+dirty" if git("status", "--porcelain", "--", *io.CODE_PATHS, "site") else "")
    head_date = git("log", "-1", "--format=%cs")
    nav = [dict(file=f, title=t) for f, t in NAV]
    manifest = build_manifest(out)
    methods = sb.md_to_html((ROOT / "docs" / "methods.md").read_text(), heading_offset=1)
    pages = {
        "index.html": ("OPM vs SQUID MEG: adult baseline", page_index(d)),
        "benchmarks.html": ("Adult reference benchmarks (G1)", page_benchmarks(d, out)),
        "adult.html": ("Realistic adult comparison (G2)", page_adult(d, out)),
        "g2-report.html": ("G2 report", sb.md_to_html((RES / "g2" / "G2_report.md").read_text(), heading_offset=1)),
        "pediatric.html": ("Pediatric extension (G3)", page_pediatric(d, out)),
        "epilepsy.html": ("Epilepsy relevance (G4)", page_epilepsy(d, out)),
        "methods.html": ("Methods, uncertainty and limitations", methods),
        "register.html": ("Parameters, provenance and assumptions",
                          sb.md_to_html((ROOT / "docs" / "provenance_register.md").read_text(), heading_offset=1)),
        "reproduce.html": ("Reproduce and download", page_reproduce(out, manifest)),
    }
    for fname, (title, content) in pages.items():
        current = "adult.html" if fname == "g2-report.html" else fname
        (out / fname).write_text(tpl.render(title=title, content=Markup(content), nav=nav, current=current,
                                            head=head, head_date=head_date))
    problems = sb.check_links(out)
    if problems:
        print("\n".join(problems))
        raise SystemExit(f"{len(problems)} broken links")
    print(f"built {len(pages)} pages, {len(manifest)} downloads, {len(list((out / 'figures').rglob('*.png')))} figures -> {out}")


if __name__ == "__main__":
    main()
