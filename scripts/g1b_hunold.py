#!/usr/bin/env python3
"""G1B: Hunold et al. (2016) depth-orientation spike simulations, MEG part.

Status: ADAPT (MNE sample subject instead of the paper's two subjects; stratified source sampling
instead of dipole traces; unstated details chosen in configs/hunold_reference.toml). The OPM array
is a NEW extension, first under the paper's brain-noise-only convention, then with intrinsic sensor
noise through a common filter.

Two background levels are reported (U-HU-bglevel):
  as_specified     the text's +/-10 nAm per background dipole (stationary, peak-normalised)
  fig6_calibrated  one scalar on all background moments so that the rendered magnetometer
                   baselines at the Fig. 6 channels (0631, 0711, 0741) match the paper's scale bars;
                   the gradiometer baselines (0412, 0413, 0423) are an independent check.
The primary numerator is peak-to-peak (Appendix E of docs/literature/hunold2016.md); 'peak' and
'noisy_peak' are reported as variants.

Inputs: cache/fullres/{neuromag4pt_hunold,opm_hunold}.npy (scripts/compute_fullres_forwards.py).
Outputs: results/g1b/ (figures, g1b_sources.csv, g1b_summary.json).
"""
from __future__ import annotations

import csv
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from scipy import stats  # noqa: E402

from opmsquid import anatomy, fullres, g2, hunold, io, neuromag, noise, paths, plotting  # noqa: E402

OUT = ROOT / "results" / "g1b"
ARRAYS = ("mag", "grad", "opm")
NUMERATORS = ("p2p", "peak", "noisy_peak", "noisy_p2p")  # p2p primary
VARIANTS = ("fig6_calibrated", "as_specified")
SENSOR_ASD = {"mag": 3.5e-15, "grad": 3.6e-13}  # brochure typical white noise (HW-mag/grad-noise)
OPM_ASD = (7e-15, 15e-15, 30e-15)  # declared sweep (A-OPM-NOISE)
SHAPE = (8, 9)


def load_inputs():
    cfg = tomllib.loads((ROOT / "configs" / "hunold_reference.toml").read_text())
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    valid_idx = np.load(fullres.directory() / "valid_index.npy")
    dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    # fingerprint- and column-checked against the arrays used here (a stale matrix is refused)
    g_nm, _ = fullres.load("neuromag4pt_hunold", neuromag.load_info("T3-4pt"), subject, cortex)
    g_opm, _ = fullres.load("opm_hunold", g2.matched_opm(subject, dig).info, subject, cortex)
    info = neuromag.load_info("T3")
    kinds = neuromag.channel_kinds(info)
    gains = {"mag": np.asarray(g_nm[kinds == "mag"]), "grad": np.asarray(g_nm[kinds == "grad"]), "opm": np.asarray(g_opm)}
    names = {k: [n for n, kk in zip(info["ch_names"], kinds) if kk == k] for k in ("mag", "grad")}
    return cfg, subject, cortex, valid_idx, gains, names


def descriptors(subject, cortex, valid_idx):
    scalp = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    skull = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    depth = np.full(cortex.n, np.nan)
    orient = np.full(cortex.n, np.nan)
    d, o = hunold.bem_node_descriptors(cortex.rr[valid_idx], cortex.nn[valid_idx], scalp["rr"], skull["rr"], skull["nn"])
    depth[valid_idx], orient[valid_idx] = d, o
    return depth, orient


def simulate_background(gains, n_bg_sources, bg_cols, n_times, fs, rng, pad_s, chunk=2000):
    traces = {k: np.zeros((g.shape[0], n_times)) for k, g in gains.items()}
    for start in range(0, n_bg_sources, chunk):
        cols = bg_cols[start:start + chunk]
        q = hunold.background_timecourses(len(cols), n_times, fs, rng, pad_s=pad_s)
        for k, g in gains.items():
            traces[k] += g[:, cols].astype(np.float64) @ q
    return traces


def fig6_realizations(gains, names, usable_pos, n_bg, n_t, fs, pad_s, n_real, seed):
    """Rendered baseline SD [px] at the Fig. 6 channels for independent background realizations
    (a new node subset and new time courses each, from their own random streams, so the primary
    realization and the sources are unchanged). Returns {kind: (n_real, n_channels)}."""
    rows = {k: [names[k].index("MEG " + ch) for ch in hunold.FIG6_BASELINE_SD_PX[k]] for k in ("mag", "grad")}
    sub = {k: np.asarray(gains[k][rows[k]]) for k in rows}
    out = {k: np.empty((n_real, len(rows[k]))) for k in rows}
    for r in range(n_real):
        rng = np.random.default_rng([seed, 1, r])
        cols = np.sort(rng.choice(usable_pos, n_bg, replace=False))
        tr = simulate_background(sub, n_bg, cols, n_t, fs, rng, pad_s)
        for k in rows:
            ppu = hunold.FIG6_SCALE_BAR_PX[k] / hunold.FIG6_SCALE_BAR[k]
            out[k][r] = [hunold.rendered_centroid_sd(tr[k][i, 300:-300], fs, ppu) for i in range(len(rows[k]))]
    return out


def fig6_expected(calib, real):
    """Calibration scalar from our expected rendered SD over background realizations (the paper's
    Fig. 6 baseline is a single draw; so is ours in one realization). The spread of the
    single-realization estimator shows how uncertain one draw makes the scalar."""
    for k in ("mag", "grad"):
        paper = np.array(list(hunold.FIG6_BASELINE_SD_PX[k].values()))
        per_draw = np.mean(paper[None, :] / real[k], axis=1)
        calib[k]["realizations"] = dict(n=int(len(real[k])), ours_px_mean=real[k].mean(axis=0).tolist(),
                                        ours_px_sd=real[k].std(axis=0, ddof=1).tolist(),
                                        expected_ratio=float(np.mean(paper / real[k].mean(axis=0))),
                                        single_draw_ratio_p5_p50_p95=np.percentile(per_draw, [5, 50, 95]).tolist())
    calib["scale_single_draw"] = calib["scale"]
    calib["scale"] = calib["mag"]["realizations"]["expected_ratio"]
    calib["alternatives"] = {"mag_by_channel_name": calib["scale"],
                             "grad_by_channel_name": calib["grad"]["realizations"]["expected_ratio"],
                             "mag_all_channel_median": calib["mag"]["ratio_to_all_channel_median"],
                             "grad_all_channel_median": calib["grad"]["ratio_to_all_channel_median"],
                             "mag_single_draw_p5": calib["mag"]["realizations"]["single_draw_ratio_p5_p50_p95"][0],
                             "mag_single_draw_p95": calib["mag"]["realizations"]["single_draw_ratio_p5_p50_p95"][2]}
    calib["scale_range"] = [min(calib["alternatives"].values()), max(calib["alternatives"].values())]
    return calib


def fig6_calibration(bg, names, fs):
    """Rendered baseline SD of our background at the Fig. 6 channels vs the paper's digitised SD.
    Returns per-kind details; the calibration scalar is the magnetometer mean ratio paper/ours."""
    out = {}
    for k in ("mag", "grad"):
        ppu = hunold.FIG6_SCALE_BAR_PX[k] / hunold.FIG6_SCALE_BAR[k]
        chans = {}
        for ch, sd_paper in hunold.FIG6_BASELINE_SD_PX[k].items():
            i = names[k].index("MEG " + ch)
            sd_ours = hunold.rendered_centroid_sd(bg[k][i, 300:-300], fs, ppu)
            chans[ch] = dict(paper_px=sd_paper, ours_px=sd_ours, ratio_paper_to_ours=sd_paper / sd_ours,
                             paper_si=sd_paper / ppu, ours_true_sd_si=float(bg[k][i, 300:-300].std()))
        every = [hunold.rendered_centroid_sd(bg[k][i, 300:-300], fs, ppu) for i in range(len(names[k]))]
        out[k] = dict(channels=chans, mean_ratio=float(np.mean([c["ratio_paper_to_ours"] for c in chans.values()])),
                      ratio_to_all_channel_median=float(np.mean(list(hunold.FIG6_BASELINE_SD_PX[k].values())) / np.median(every)))
    out["scale"] = out["mag"]["mean_ratio"]
    out["alternatives"] = {"mag_by_channel_name": out["mag"]["mean_ratio"], "grad_by_channel_name": out["grad"]["mean_ratio"],
                           "mag_all_channel_median": out["mag"]["ratio_to_all_channel_median"],
                           "grad_all_channel_median": out["grad"]["ratio_to_all_channel_median"]}
    out["scale_range"] = [min(out["alternatives"].values()), max(out["alternatives"].values())]
    return out


def fig6_spike_comparison(gains, cortex, depth, orient, region, w_unit):
    """Noise-free spike peak-to-peak on the best (largest-amplitude) channel for 600-nAm dipoles
    resembling the Table 1 tangential examples (depth +/- 1.5 mm, orientation +/- 5 deg, left
    posterior frontal cortex: DK precentral, caudal middle frontal, pars opercularis), against the
    paper's Fig. 6 traces (centroid and drawn-line p2p, scale bars)."""
    p2p_unit = w_unit.max() - w_unit.min()
    out = {}
    for src, (d0, o0) in ((k, v["tangential"]) for k, v in hunold.FIG6_EXAMPLES.items()):
        m = (cortex.usable & (cortex.hemi == 0) & np.isin(region, hunold.POSTERIOR_FRONTAL)
             & (np.abs(depth - d0) <= 1.5) & (np.abs(orient - o0) <= 5.0))
        cols = np.flatnonzero(m[cortex.valid])  # positions in the valid-vertex gain columns
        row = dict(n_sources=int(len(cols)))
        for k, unit, lab in (("mag", 1e12, "pT"), ("grad", 1e12, "pT/m")):
            amp = np.abs(gains[k][:, cols]).max(axis=0) * 600e-9 * p2p_unit * unit
            cen, line = hunold.FIG6_SPIKE_P2P_PX[k][(src, "tangential")]
            to_si = hunold.FIG6_SCALE_BAR[k] / hunold.FIG6_SCALE_BAR_PX[k] * unit
            row[k] = dict(unit=lab, ours_median=float(np.median(amp)) if len(amp) else None,
                          ours_iqr=[float(x) for x in np.percentile(amp, [25, 75])] if len(amp) else None,
                          paper_centroid=cen * to_si, paper_drawn_line=line * to_si)
        out[src] = row
    return out


def per_bin_test(a, b, rows, cols, shape):
    """Unpaired test per bin as in the paper: Student's t if both groups pass Shapiro-Wilk
    (p > 0.05), else the Wilcoxon rank-sum test. Returns p-values."""
    p = np.full(shape, np.nan)
    for i in range(shape[0]):
        for j in range(shape[1]):
            m = (rows == i) & (cols == j)
            if m.sum() < 3:
                continue
            x, y = a[m], b[m]
            normal = stats.shapiro(x).pvalue > 0.05 and stats.shapiro(y).pvalue > 0.05
            p[i, j] = stats.ttest_ind(x, y).pvalue if normal else stats.ranksums(x, y).pvalue
    return p


def paper_map(kind, family):
    """The paper's digitised bin classes (midpoints) on our 8 x 9 grid (NaN where not shown)."""
    m = np.full(SHAPE, np.nan)
    p = hunold.PAPER_BIN_CLASS[(family, kind)] + 0.25
    if family == "dipole":
        m[:] = p
    else:
        m[:7, 1:] = p
    return m


def compare_with_paper(ours, kind, family):
    """Bin means vs the paper's digitised classes (midpoint = lower edge + 0.25): correlation,
    mean ratio (overall, for bins whose paper class is >= 2.5 ('strong') and below ('weak')),
    agreement on SNR >= 2.5 (paper class lower edge), and the share inside the paper's class."""
    paper = paper_map(kind, family)
    ok = np.isfinite(ours) & np.isfinite(paper)
    o, p = ours[ok], paper[ok]
    lower = p - 0.25
    strong = lower >= 2.5
    return dict(pearson_r=float(np.corrcoef(o, p)[0, 1]), mean_ratio_ours_to_paper=float(np.mean(o / p)),
                median_ratio_ours_to_paper=float(np.median(o / p)), mean_abs_diff=float(np.mean(np.abs(o - p))),
                mean_ratio_strong_bins=float(np.mean(o[strong] / p[strong])), mean_ratio_weak_bins=float(np.mean(o[~strong] / p[~strong])),
                threshold_2p5_agreement=float(np.mean((o >= 2.5) == strong)), share_ge_2p5_ours=float(np.mean(o >= 2.5)),
                share_ge_2p5_paper=float(np.mean(strong)), within_paper_class=float(np.mean((o >= lower) & (o < lower + 0.5))),
                n_bins=int(ok.sum()))


def gm_minus_mm_sign_agreement(maps, family):
    """Share of bins where our GM - MM difference has the paper's sign (bins with a paper
    difference of at least one colour class, 0.5)."""
    d_paper = paper_map("grad", family) - paper_map("mag", family)
    d_ours = maps["grad"] - maps["mag"]
    ok = np.isfinite(d_paper) & np.isfinite(d_ours) & (np.abs(d_paper) >= 0.5)
    return dict(n_bins=int(ok.sum()), agree=float(np.mean(np.sign(d_ours[ok]) == np.sign(d_paper[ok]))) if ok.any() else None)


def plot_bins(maps, pvals, family, variant, numerator, fname):
    names = [("mag", "SQUID magnetometers"), ("grad", "SQUID gradiometers"), ("opm", "OPM (matched)")]
    diffs = [("grad", "mag", "GM - MM"), ("opm", "mag", "OPM - MM"), ("opm", "grad", "OPM - GM")]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.6))
    fig.subplots_adjust(left=0.06, right=0.97, top=0.9, bottom=0.07, wspace=0.3, hspace=0.35)
    ext = [0, 90, 60, 20]
    for ax, (k, title) in zip(axs[0], names):
        im = ax.imshow(maps[k], origin="upper", extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=10)
        if np.nanmax(maps[k]) > 2.5 > np.nanmin(maps[k]):
            ax.contour(np.linspace(5, 85, SHAPE[1]), np.linspace(22.5, 57.5, SHAPE[0]), maps[k], levels=[2.5], colors="w", linewidths=1.2)
        ax.set_title(f"{title}: mean SNR", fontsize=9)
        ax.set_xlabel("Orientation [deg] (0 radial, 90 tangential)")
        ax.set_ylabel("Depth below scalp [mm]")
        fig.colorbar(im, ax=ax, fraction=0.046)
    for ax, (a, b, title) in zip(axs[1], diffs):
        d = maps[a] - maps[b]
        lim = np.nanmax(np.abs(d)) if np.isfinite(d).any() else 1
        im = ax.imshow(d, origin="upper", extent=ext, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
        p = pvals[(a, b)]
        for i in range(SHAPE[0]):
            for j in range(SHAPE[1]):
                if np.isfinite(p[i, j]) and p[i, j] < 0.05:
                    ax.text(5 + 10 * j, 22.5 + 5 * i, "**" if p[i, j] < 0.01 else "*", ha="center", va="center", fontsize=7)
        ax.set_title(f"{title} (* p<0.05, ** p<0.01, unpaired)", fontsize=9)
        ax.set_xlabel("Orientation [deg]")
        ax.set_ylabel("Depth below scalp [mm]")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(f"G1B Hunold adaptation, {family} sources, background {variant.replace('_', ' ')}, numerator '{numerator}' "
                 "(brain noise only; white line: 2.5)", fontsize=10)
    fig.savefig(fname, dpi=160)
    plt.close(fig)


def plot_vs_paper(all_maps, calib, fname):
    """Paper vs ours (Fig. 6-calibrated background, p2p numerator) for MM and GM, and bin-mean
    scatter for both background levels."""
    fig = plt.figure(figsize=(15, 7.4))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 1.15], left=0.05, right=0.98, top=0.88, bottom=0.08,
                          wspace=0.35, hspace=0.4)
    ext = [0, 90, 60, 20]
    for r, fam in enumerate(("dipole", "patch")):
        maps = all_maps[("fig6_calibrated", fam, "p2p")]
        for c, (k, src) in enumerate((("mag", "paper"), ("mag", "ours"), ("grad", "paper"), ("grad", "ours"))):
            ax = fig.add_subplot(gs[r, c])
            m = paper_map(k, fam) if src == "paper" else np.where(np.isfinite(paper_map(k, fam)), maps[k], np.nan)
            im = ax.imshow(m, origin="upper", extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=10)
            ax.set_title(f"{'MM' if k == 'mag' else 'GM'} {fam}: {src}", fontsize=9)
            ax.set_xlabel("Orientation [deg]", fontsize=8)
            if c == 0:
                ax.set_ylabel("Depth below scalp [mm]", fontsize=8)
            ax.tick_params(labelsize=7)
        ax = fig.add_subplot(gs[r, 4])
        for var, mk in (("fig6_calibrated", "o"), ("as_specified", "x")):
            for k, col in (("mag", "tab:blue"), ("grad", "tab:red")):
                o = all_maps[(var, fam, "p2p")][k]
                p = paper_map(k, fam)
                ok = np.isfinite(o) & np.isfinite(p)
                ax.scatter(p[ok] + np.random.default_rng(0).uniform(-0.1, 0.1, ok.sum()), o[ok], s=10, marker=mk, color=col,
                           alpha=0.7, label=f"{'MM' if k == 'mag' else 'GM'}, {var.replace('_', ' ')}")
        ax.plot([0, 10], [0, 10], "k--", lw=0.8)
        ax.set_xlim(0, 10.5)
        ax.set_ylim(0, 10.5)
        ax.set_xlabel("paper bin class midpoint", fontsize=8)
        ax.set_ylabel("ours (p2p)", fontsize=8)
        ax.set_title(f"{fam}: bin means", fontsize=9)
        ax.tick_params(labelsize=7)
        if r == 0:
            ax.legend(fontsize=6.5, loc="upper left")
    fig.colorbar(im, ax=fig.axes[:4], fraction=0.02, pad=0.01)
    fig.suptitle("G1B vs Hunold et al. (2016) Figs 4-5 (digitised classes): p2p numerator, background calibrated with one scalar "
                 f"({calib['scale']:.2f}; {calib['scale_range'][0]:.2f}-{calib['scale_range'][1]:.2f} over channel choices and single draws) "
                 "to the Fig. 6 baselines", fontsize=10)
    fig.savefig(fname, dpi=150)
    plt.close(fig)


def main():
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("Figure_G1B_bins_*.png"):
        old.unlink()
    mne.set_log_level("WARNING")
    cfg, subject, cortex, valid_idx, gains, names = load_inputs()
    col_of = np.full(cortex.n, -1)
    col_of[valid_idx] = np.arange(len(valid_idx))
    depth, orient = descriptors(subject, cortex, valid_idx)
    rng = np.random.default_rng(cfg["sources"]["seed"])
    fs, n_t = float(cfg["waveform"]["sfreq_hz"]), int(cfg["background"]["duration_s"] * cfg["waveform"]["sfreq_hz"])
    onset = int(cfg["snr"]["onset_s"] * fs)
    base = slice(onset - int(cfg["snr"]["baseline_s"] * fs), onset)

    # background: random 10 % of valid nodes, one fixed realization shared by all arrays and sources
    usable_pos = np.flatnonzero(cortex.usable[valid_idx])  # A-BEM-DIST: >= 4 mm from the inner-skull mesh
    usable_idx = valid_idx[usable_pos]
    n_bg = int(round(cfg["background"]["fraction_of_nodes"] * len(usable_idx)))
    bg_cols = np.sort(rng.choice(usable_pos, n_bg, replace=False))
    t0 = time.time()
    bg_spec = simulate_background(gains, n_bg, bg_cols, n_t, fs, rng, cfg["background"]["pad_s"])
    print(f"background: {n_bg} dipoles ({cortex.area[valid_idx][bg_cols].mean() * 1e6 * 10:.2f} mm2 per dipole incl. the "
          f"other 90 %) in {time.time() - t0:.0f} s")
    calib = fig6_calibration(bg_spec, names, fs)
    t0 = time.time()
    real = fig6_realizations(gains, names, usable_pos, n_bg, n_t, fs, cfg["background"]["pad_s"],
                             cfg["background"]["calibration_realizations"], cfg["sources"]["seed"])
    calib = fig6_expected(calib, real)
    print(f"Fig. 6 calibration ({calib['mag']['realizations']['n']} realizations, {time.time() - t0:.0f} s): paper/ours baseline SD "
          f"mag {calib['scale']:.3f} (single draws 5-95 %: {calib['alternatives']['mag_single_draw_p5']:.3f}-"
          f"{calib['alternatives']['mag_single_draw_p95']:.3f}; primary realization {calib['scale_single_draw']:.3f}), "
          f"grad {calib['alternatives']['grad_by_channel_name']:.3f} (all-channel medians "
          f"{calib['mag']['ratio_to_all_channel_median']:.3f}, {calib['grad']['ratio_to_all_channel_median']:.3f})")
    bgs = {"as_specified": bg_spec, "fig6_calibrated": {k: v * calib["scale"] for k, v in bg_spec.items()}}

    # sources: stratified to the paper's per-bin counts
    dip, achieved = hunold.sample_by_bins(depth, orient, usable_idx, hunold.PAPER_DIPOLE_COUNTS, rng)
    w_unit = hunold.spike_waveform(fs, peak=1.0)
    topo = {"dipole": {k: g[:, col_of[dip]].astype(np.float64) * 600e-9 for k, g in gains.items()}}
    # patches: grown from every sampled dipole; fixed density giving a median total of 622 nAm
    patches, keep = [], []
    for s in dip:
        m = hunold.grow_patch(int(s), cortex.adjacency, orient, cortex.area, cortex.usable,
                              target_area=cfg["sources"]["patch_target_area_mm2"] * 1e-6,
                              window=cfg["sources"]["patch_orientation_window_deg"])
        if m is not None:
            patches.append(m)
            keep.append(s)
    p_area = np.array([cortex.area[m].sum() for m in patches])
    density = cfg["sources"]["patch_median_total_nAm"] * 1e-9 / np.median(p_area)
    p_total = density * p_area
    topo["patch"] = {k: np.stack([g[:, col_of[m]].astype(np.float64) @ (density * cortex.area[m]) for m in patches], axis=1)
                     for k, g in gains.items()}
    p_depth = np.array([depth[m].mean() for m in patches])
    p_orient = np.array([orient[m].mean() for m in patches])
    fam_desc = {"dipole": (depth[dip], orient[dip]), "patch": (p_depth, p_orient)}
    bins_of = {fam: hunold.bin_index(*fam_desc[fam]) for fam in fam_desc}

    results, all_maps, variants_out = {}, {}, {}
    for var in VARIANTS:
        bg = bgs[var]
        a_bg = {k: hunold.background_amplitude(v, base) for k, v in bg.items()}
        summary_bins, comparisons, significance, gm_mm = {}, {}, {}, {}
        for fam in ("dipole", "patch"):
            rows, cols = bins_of[fam]
            results[(var, fam)] = {k: hunold.spike_snr(topo[fam][k], w_unit, bg[k], onset, a_bg[k]) for k in ARRAYS}
            for num in NUMERATORS:
                maps = {k: hunold.bin_means(results[(var, fam)][k][num], rows, cols, SHAPE)[0] for k in ARRAYS}
                counts = hunold.bin_means(results[(var, fam)]["mag"][num], rows, cols, SHAPE)[1]
                all_maps[(var, fam, num)] = maps
                pv = {(a, b): per_bin_test(results[(var, fam)][a][num], results[(var, fam)][b][num], rows, cols, SHAPE)
                      for a, b in (("grad", "mag"), ("opm", "mag"), ("opm", "grad"))}
                key = f"{fam}/{num}"
                summary_bins[key] = dict(counts=counts, mean_snr=maps,
                                         share_of_bins_ge_2p5={k: float(np.mean(maps[k][np.isfinite(maps[k])] >= 2.5)) for k in ARRAYS})
                comparisons[key] = {k: compare_with_paper(maps[k], k, fam) for k in ("mag", "grad")}
                gm_mm[key] = gm_minus_mm_sign_agreement(maps, fam)
                significance[key] = {f"{a}-{b}": pv[(a, b)] for a, b in pv}
                plot_bins(maps, pv, fam, var, num, OUT / f"Figure_G1B_bins_{var}_{fam}_{num}.png")
        variants_out[var] = dict(background_amplitude_median={k: float(np.median(v)) for k, v in a_bg.items()},
                                 bins=summary_bins, comparison_with_paper=comparisons, gm_minus_mm_sign_agreement=gm_mm,
                                 significance=significance)
        for key, comp in comparisons.items():
            print(var, key, {k: (round(v["pearson_r"], 3), round(v["mean_ratio_ours_to_paper"], 2)) for k, v in comp.items()},
                  "GM-MM sign agreement", gm_mm[key]["agree"])
    plot_vs_paper(all_maps, calib, OUT / "Figure_G1B_vs_paper.png")
    # the calibration scalar is uncertain (channel choice, single-draw spread): noise-free numerators scale as 1/k
    calib_sens = {}
    for alt, k_alt in calib["alternatives"].items():
        f = calib["scale"] / k_alt
        calib_sens[alt] = dict(scalar=k_alt, factor=f, **{
            f"{fam}/{num}/{k}": float(variants_out["fig6_calibrated"]["comparison_with_paper"][f"{fam}/{num}"][k]["mean_ratio_ours_to_paper"] * f)
            for fam in ("dipole", "patch") for num in ("p2p", "peak") for k in ("mag", "grad")})
    regions = np.concatenate([plotting.read_freesurfer_annot(paths.SUBJECTS_DIR / "sample" / "label" / f"{h}.aparc.annot")
                              for h in ("lh", "rh")])
    spike_cmp = fig6_spike_comparison(gains, cortex, depth, orient, regions, w_unit)

    # extension: intrinsic sensor noise, common 0.5-70 Hz filter on spike, background and noise
    filt = noise.AnalysisFilter(fs=fs, l_freq=0.5, h_freq=70.0, order=4)
    w_f = filt.apply(np.r_[np.zeros(onset), w_unit, np.zeros(n_t - onset - len(w_unit))])[onset:onset + len(w_unit) + 100]
    ext = {}
    for var in VARIANTS:
        ext[var] = {}
        for fam in ("dipole", "patch"):
            ext[var][fam] = {}
            rows, cols = bins_of[fam]
            for k in ARRAYS:
                for asd in (OPM_ASD if k == "opm" else (SENSOR_ASD[k],)):
                    n = noise.white_noise(asd, fs, bgs[var][k].shape, np.random.default_rng(7))
                    tr = filt.apply(bgs[var][k] + n)
                    a = hunold.background_amplitude(tr, base)
                    snr = hunold.spike_snr(topo[fam][k], w_f, tr, onset, a)["p2p"]
                    ext[var][fam][f"{k}@{asd * 1e15:g}"] = hunold.bin_means(snr, rows, cols, SHAPE)[0]

    # per-source table
    with open(OUT / "g1b_sources.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["family", "vertex", "hemi", "depth_mm", "orientation_deg", "area_mm2", "total_nAm"]
                    + [f"SNR_{var}_{k}_{v}" for var in VARIANTS for k in ARRAYS for v in NUMERATORS])
        for fam, verts, areas, totals in (("dipole", dip, np.zeros(len(dip)), np.full(len(dip), 600.0)),
                                          ("patch", np.array(keep), p_area * 1e6, p_total * 1e9)):
            d, o = fam_desc[fam]
            for i, v in enumerate(verts):
                wr.writerow([fam, int(cortex.vertno[v]), int(cortex.hemi[v]), f"{d[i]:.4f}", f"{o[i]:.4f}", f"{areas[i]:.3f}",
                             f"{totals[i]:.1f}"]
                            + [f"{results[(var, fam)][k][vv][i]:.4f}" for var in VARIANTS for k in ARRAYS for vv in NUMERATORS])

    summary = dict(
        status="ADAPT (MEG part of Hunold et al. 2016 on the MNE sample subject); OPM columns are a NEW extension",
        config=cfg, n_valid_vertices=int(len(valid_idx)), n_usable_vertices=int(len(usable_idx)), n_background_dipoles=int(n_bg),
        dipoles=dict(requested=int(hunold.PAPER_DIPOLE_COUNTS.sum()), achieved=int(achieved.sum()),
                     achieved_per_bin=achieved, shortfall_bins=int(np.sum(achieved < hunold.PAPER_DIPOLE_COUNTS))),
        patches=dict(seeds=int(len(dip)), grown=int(len(patches)), density_nAm_per_mm2=density * 1e9 / 1e6,
                     area_mm2=dict(min=float(p_area.min() * 1e6), median=float(np.median(p_area) * 1e6), max=float(p_area.max() * 1e6)),
                     total_nAm=dict(min=float(p_total.min() * 1e9), median=float(np.median(p_total) * 1e9), max=float(p_total.max() * 1e9)),
                     n_dipoles=dict(min=int(min(map(len, patches))), median=float(np.median(list(map(len, patches)))),
                                    max=int(max(map(len, patches))))),
        numerators=dict(primary="p2p", variants=list(NUMERATORS[1:]),
                        basis=("Fig. 6 digitisation (scripts/digitise_hunold_fig6.py, 8 MEG traces with printed SNR > 2.4): "
                               "p2p / (2 mean|hilbert|) = 0.87-1.10 x printed depending on digitisation handling; literal "
                               "noisy peak 0.61-0.73 x")),
        fig6_calibration=calib, calibration_sensitivity=calib_sens, fig6_spike_comparison=spike_cmp, variants=variants_out,
        extension_intrinsic_noise=dict(description="0.5-70 Hz zero-phase Butterworth (order 4) on spike, background and white "
                                                   "sensor noise; SQUID brochure noise; OPM 7/15/30 fT/sqrt(Hz); p2p numerator",
                                       mean_snr=ext),
        runtime_s=time.time() - t_start)
    io.write_json(summary, OUT / "g1b_summary.json")
    print(f"patches: {len(patches)}/{len(dip)} grown; density {density * 1e3:.1f} nAm/mm2; totals "
          f"{p_total.min() * 1e9:.0f}-{p_total.max() * 1e9:.0f} nAm (median {np.median(p_total) * 1e9:.0f})")
    print(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
