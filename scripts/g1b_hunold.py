#!/usr/bin/env python3
"""G1B: Hunold et al. (2016) depth-orientation spike simulations, MEG part.

Status: ADAPT (MNE sample subject instead of the paper's two subjects; stratified source sampling
instead of dipole traces; unstated details chosen in configs/hunold_reference.toml). The OPM array
is a NEW extension, first under the paper's brain-noise-only convention, then with intrinsic sensor
noise through a common filter.

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

from opmsquid import anatomy, hunold, io, neuromag, noise, paths  # noqa: E402

OUT = ROOT / "results" / "g1b"
ARRAYS = ("mag", "grad", "opm")
LABEL = {"mag": "SQUID magnetometers (102)", "grad": "SQUID planar gradiometers (204)", "opm": "OPM (99, matched sites)"}
SENSOR_ASD = {"mag": 3.5e-15, "grad": 3.6e-13}  # brochure typical white noise (HW-mag/grad-noise)
OPM_ASD = (7e-15, 15e-15, 30e-15)  # declared sweep (A-OPM-NOISE)


def load_inputs():
    cfg = tomllib.loads((ROOT / "configs" / "hunold_reference.toml").read_text())
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    fr = paths.CACHE / "fullres"
    valid_idx = np.load(fr / "valid_index.npy")
    g_nm = np.load(fr / "neuromag4pt_hunold.npy", mmap_mode="r")
    g_opm = np.load(fr / "opm_hunold.npy", mmap_mode="r")
    kinds = neuromag.channel_kinds(neuromag.load_info("T3"))
    gains = {"mag": np.asarray(g_nm[kinds == "mag"]), "grad": np.asarray(g_nm[kinds == "grad"]), "opm": np.asarray(g_opm)}
    return cfg, subject, cortex, valid_idx, gains


def descriptors(subject, cortex, valid_idx):
    scalp = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    skull = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    depth = np.full(cortex.n, np.nan)
    orient = np.full(cortex.n, np.nan)
    d, o = hunold.bem_node_descriptors(cortex.rr[valid_idx], cortex.nn[valid_idx], scalp["rr"], skull["rr"], skull["nn"])
    depth[valid_idx], orient[valid_idx] = d, o
    return depth, orient


def simulate_background(gains, n_bg_sources, bg_cols, n_times, fs, rng, chunk=2000):
    traces = {k: np.zeros((g.shape[0], n_times)) for k, g in gains.items()}
    for start in range(0, n_bg_sources, chunk):
        cols = bg_cols[start:start + chunk]
        q = hunold.background_timecourses(len(cols), n_times, fs, rng)
        for k, g in gains.items():
            traces[k] += g[:, cols].astype(np.float64) @ q
    return traces


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
            normal = stats.shapiro(x).pvalue > 0.05 and stats.shapiro(y).pvalue > 0.05 if m.sum() >= 3 else False
            p[i, j] = stats.ttest_ind(x, y).pvalue if normal else stats.ranksums(x, y).pvalue
    return p


def compare_with_paper(ours, kind, family):
    """Our bin means vs the paper's digitised classes (class midpoints x + 0.25)."""
    paper = hunold.PAPER_BIN_CLASS[(family, kind)] + 0.25
    o = ours if family == "dipole" else ours[:7, 1:]
    ok = np.isfinite(o)
    return dict(pearson_r=float(np.corrcoef(o[ok], paper[ok])[0, 1]), mean_ratio_ours_to_paper=float(np.mean(o[ok] / paper[ok])),
                mean_abs_diff=float(np.mean(np.abs(o[ok] - paper[ok]))), within_half_class=float(np.mean(np.abs(o[ok] - paper[ok]) <= 0.5)),
                n_bins=int(ok.sum()))


def plot_bins(maps, pvals, counts, family, variant, fname):
    names = [("mag", "SQUID magnetometers"), ("grad", "SQUID gradiometers"), ("opm", "OPM (matched)")]
    diffs = [("grad", "mag", "GM - MM"), ("opm", "mag", "OPM - MM"), ("opm", "grad", "OPM - GM")]
    rows, cols = (8, 9)
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.6))
    fig.subplots_adjust(left=0.06, right=0.97, top=0.9, bottom=0.07, wspace=0.3, hspace=0.35)
    ext = [0, 90, 60, 20]
    for ax, (k, title) in zip(axs[0], names):
        im = ax.imshow(maps[k], origin="upper", extent=ext, aspect="auto", cmap="viridis", vmin=0, vmax=10)
        ax.contour(np.linspace(5, 85, cols), np.linspace(22.5, 57.5, rows), maps[k], levels=[2.5], colors="w", linewidths=1.2)
        ax.set_title(f"{title}: mean SNR", fontsize=9)
        ax.set_xlabel("Orientation [deg] (0 radial, 90 tangential)")
        ax.set_ylabel("Depth below scalp [mm]")
        fig.colorbar(im, ax=ax, fraction=0.046)
    for ax, (a, b, title) in zip(axs[1], diffs):
        d = maps[a] - maps[b]
        lim = np.nanmax(np.abs(d)) if np.isfinite(d).any() else 1
        im = ax.imshow(d, origin="upper", extent=ext, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
        p = pvals[(a, b)]
        for i in range(rows):
            for j in range(cols):
                if np.isfinite(p[i, j]) and p[i, j] < 0.05:
                    ax.text(5 + 10 * j, 22.5 + 5 * i, "**" if p[i, j] < 0.01 else "*", ha="center", va="center", fontsize=7)
        ax.set_title(f"{title} (* p<0.05, ** p<0.01, unpaired)", fontsize=9)
        ax.set_xlabel("Orientation [deg]")
        ax.set_ylabel("Depth below scalp [mm]")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(f"G1B Hunold adaptation, {family} sources, SNR numerator '{variant}' (brain noise only; white line: 2.5)",
                 fontsize=10)
    fig.savefig(fname, dpi=160)
    plt.close(fig)


def main():
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    mne.set_log_level("WARNING")
    cfg, subject, cortex, valid_idx, gains = load_inputs()
    col_of = np.full(cortex.n, -1)
    col_of[valid_idx] = np.arange(len(valid_idx))
    depth, orient = descriptors(subject, cortex, valid_idx)
    rng = np.random.default_rng(cfg["sources"]["seed"])
    fs, n_t = float(cfg["waveform"]["sfreq_hz"]), int(cfg["background"]["duration_s"] * cfg["waveform"]["sfreq_hz"])
    onset = int(cfg["snr"]["onset_s"] * fs)
    base = slice(onset - int(cfg["snr"]["baseline_s"] * fs), onset)

    # background: random 10 % of valid nodes, one fixed realization shared by all arrays and sources
    n_bg = int(round(cfg["background"]["fraction_of_nodes"] * len(valid_idx)))
    bg_cols = np.sort(rng.choice(len(valid_idx), n_bg, replace=False))
    t0 = time.time()
    bg = simulate_background(gains, n_bg, bg_cols, n_t, fs, rng)
    print(f"background: {n_bg} dipoles ({cortex.area[valid_idx][bg_cols].mean() * 1e6 * 10:.2f} mm2 per dipole incl. the "
          f"other 90 %) in {time.time() - t0:.0f} s")
    a_bg = {k: hunold.background_amplitude(v, base) for k, v in bg.items()}

    # sources: stratified to the paper's per-bin counts
    dip, achieved = hunold.sample_by_bins(depth, orient, valid_idx, hunold.PAPER_DIPOLE_COUNTS, rng)
    w_unit = hunold.spike_waveform(fs, peak=1.0)
    topo = {"dipole": {k: g[:, col_of[dip]].astype(np.float64) * 600e-9 for k, g in gains.items()}}
    # patches: grown from every sampled dipole; fixed density giving a median total of 622 nAm
    patches, keep = [], []
    for s in dip:
        m = hunold.grow_patch(int(s), cortex.adjacency, orient, cortex.area, cortex.valid,
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
    results, summary_bins, comparisons, significance = {}, {}, {}, {}
    for fam in ("dipole", "patch"):
        d, o = fam_desc[fam]
        rows, cols = hunold.bin_index(d, o)
        results[fam] = {k: hunold.spike_snr(topo[fam][k], w_unit, bg[k], onset, a_bg[k]) for k in ARRAYS}
        for variant in ("peak", "p2p", "noisy_peak"):
            maps, counts = {}, None
            for k in ARRAYS:
                maps[k], counts = hunold.bin_means(results[fam][k][variant], rows, cols, (8, 9))
            pv = {(a, b): per_bin_test(results[fam][a][variant], results[fam][b][variant], rows, cols, (8, 9))
                  for a, b in (("grad", "mag"), ("opm", "mag"), ("opm", "grad"))}
            summary_bins[f"{fam}/{variant}"] = dict(counts=counts, mean_snr={k: maps[k] for k in ARRAYS},
                                                    share_of_bins_ge_2p5={k: float(np.nanmean(maps[k] >= 2.5)) for k in ARRAYS})
            comparisons[f"{fam}/{variant}"] = {k: compare_with_paper(maps[k], k, fam) for k in ("mag", "grad")}
            significance[f"{fam}/{variant}"] = {f"{a}-{b}": pv[(a, b)] for a, b in pv}
            plot_bins(maps, pv, counts, fam, variant, OUT / f"Figure_G1B_bins_{fam}_{variant}.png")

    # extension: intrinsic sensor noise, common 0.5-70 Hz filter on spike, background and noise
    filt = noise.AnalysisFilter(fs=fs, l_freq=0.5, h_freq=70.0, order=4)
    w_f = filt.apply(np.r_[np.zeros(onset), w_unit, np.zeros(n_t - onset - len(w_unit))])[onset:onset + len(w_unit) + 100]
    ext = {}
    for fam in ("dipole", "patch"):
        ext[fam] = {}
        for k in ARRAYS:
            asds = OPM_ASD if k == "opm" else (SENSOR_ASD[k],)
            for asd in asds:
                n = noise.white_noise(asd, fs, bg[k].shape, np.random.default_rng(7))
                tr = filt.apply(bg[k] + n)
                a = hunold.background_amplitude(tr, base)
                snr = hunold.spike_snr(topo[fam][k], w_f, tr, onset, a)["peak"]
                d, o = fam_desc[fam]
                rows, cols = hunold.bin_index(d, o)
                ext[fam][f"{k}@{asd * 1e15:g}"] = hunold.bin_means(snr, rows, cols, (8, 9))[0]

    # per-source table
    with open(OUT / "g1b_sources.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["family", "vertex", "hemi", "depth_mm", "orientation_deg", "area_mm2", "total_nAm"]
                    + [f"SNR_{k}_{v}" for k in ARRAYS for v in ("peak", "p2p", "noisy_peak")])
        for fam, verts, areas, totals in (("dipole", dip, np.zeros(len(dip)), np.full(len(dip), 600.0)),
                                          ("patch", np.array(keep), p_area * 1e6, p_total * 1e9)):
            d, o = fam_desc[fam]
            for i, v in enumerate(verts):
                wr.writerow([fam, int(cortex.vertno[v]), int(cortex.hemi[v]), f"{d[i]:.2f}", f"{o[i]:.2f}", f"{areas[i]:.2f}",
                             f"{totals[i]:.1f}"] + [f"{results[fam][k][vv][i]:.4f}" for k in ARRAYS for vv in ("peak", "p2p", "noisy_peak")])

    summary = dict(
        status="ADAPT (MEG part of Hunold et al. 2016 on the MNE sample subject); OPM columns are a NEW extension",
        config=cfg, n_valid_vertices=int(len(valid_idx)), n_background_dipoles=int(n_bg),
        dipoles=dict(requested=int(hunold.PAPER_DIPOLE_COUNTS.sum()), achieved=int(achieved.sum()),
                     achieved_per_bin=achieved, shortfall_bins=int(np.sum(achieved < hunold.PAPER_DIPOLE_COUNTS))),
        patches=dict(seeds=int(len(dip)), grown=int(len(patches)), density_nAm_per_mm2=density * 1e9 / 1e6,
                     area_mm2=dict(min=float(p_area.min() * 1e6), median=float(np.median(p_area) * 1e6), max=float(p_area.max() * 1e6)),
                     total_nAm=dict(min=float(p_total.min() * 1e9), median=float(np.median(p_total) * 1e9), max=float(p_total.max() * 1e9)),
                     n_dipoles=dict(min=int(min(map(len, patches))), median=float(np.median(list(map(len, patches)))),
                                    max=int(max(map(len, patches))))),
        background_amplitude_median={k: float(np.median(v)) for k, v in a_bg.items()},
        bins=summary_bins, comparison_with_paper=comparisons, significance=significance,
        extension_intrinsic_noise=dict(description="0.5-70 Hz zero-phase Butterworth (order 4) on spike, background and white "
                                                   "sensor noise; SQUID brochure noise; OPM 7/15/30 fT/sqrt(Hz); 'peak' numerator",
                                       mean_snr=ext),
        runtime_s=time.time() - t_start)
    io.write_json(summary, OUT / "g1b_summary.json")
    for key, comp in comparisons.items():
        print(key, {k: (round(v["pearson_r"], 3), round(v["mean_ratio_ours_to_paper"], 2)) for k, v in comp.items()})
    print(f"patches: {len(patches)}/{len(dip)} grown; density {density * 1e3:.1f} nAm/mm2; totals "
          f"{p_total.min() * 1e9:.0f}-{p_total.max() * 1e9:.0f} nAm (median {np.median(p_total) * 1e9:.0f})")
    print(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
