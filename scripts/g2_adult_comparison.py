#!/usr/bin/env python3
"""G2: realistic adult OPM vs Neuromag comparison (NEW).

MNE sample subject in its measured Vectorview head position. Arrays (src/opmsquid/g2.py):
Neuromag T3 (102 magnetometers + 204 planar gradiometers; channel sets mag, grad, combined), the
matched-site OPM array (opm_matched, coverage control), 204 OPM sites spread over the densest feasible
array (opm204, channel-budget control) and the densest feasible single-axis array (opm_dense, full
system). Common target sources (10-nAm cortical-normal dipoles at the valid oct-6 vertices and
geodesic patches) and common noise for every array (src/opmsquid/noisemodel.py): intrinsic white
noise, a cortical background calibrated once to the measured gradiometer brain noise, and the room
field fitted to the empty-room recording, all in one analysis band. Conditions: intrinsic,
intrinsic+brain, intrinsic+brain+env, projected (8-term external subspace removed identically).
Metrics: peak-channel SNR, mean-power SNR_dB, known-topography detectability with the oracle and
with plug-in (estimated) covariances. Sensitivity: OPM noise, scalp gap, head position, background
correlation, covariance estimation. Convergence: background discretisation, BEM layers and
refinement, coil integration, OPM cell, target sampling, whitening tolerance.

Configuration: configs/g2_adult.toml. Outputs: results/g2/.
"""
from __future__ import annotations

import csv
import json
import pickle
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

from opmsquid import (anatomy, background, forward, g2, goldenholz, io, metrics, neuromag, noise,  # noqa: E402
                      noisemodel, opm, paths, plotting)

OUT = ROOT / "results" / "g2"
STATE = ROOT / "cache" / "g2" / "state.pkl"
OPMS = ("opm_matched", "opm204", "opm_dense")
REFS = ("combined", "grad", "mag")
LABEL = {"squid": "Neuromag", "opm_matched": "OPM matched", "opm204": "OPM 204", "opm_dense": "OPM dense",
         "combined": "Neuromag combined", "grad": "Neuromag grad", "mag": "Neuromag mag"}
DEPTH_EDGES = np.arange(10.0, 85.0, 5.0)
ORIENT_EDGES = np.arange(0.0, 90.1, 10.0)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------------------------
# statistics helpers
GROUPS = None  # parcel label per target, set in main(): the bootstrap resamples parcels (targets are correlated)


def compare(d_a, d_b, rng, n_boot=1000, groups=None):
    """Median log2(d_a / d_b) over targets with a 95 % CI from a bootstrap over cortical parcels
    (Desikan-Killiany regions resampled with replacement; neighbouring targets are not
    independent), the share of targets and of parcels (parcel median) with d_a > d_b."""
    lr = np.log2(np.asarray(d_a) / np.asarray(d_b))
    g = GROUPS if groups is None else groups
    ok = np.isfinite(lr)
    lr = lr[ok]
    out = dict(median_log2=float(np.median(lr)), share_opm_better=float(np.mean(lr > 0)), n=int(len(lr)))
    if g is None or len(g) != len(ok):
        idx = rng.integers(0, len(lr), (n_boot, len(lr)))
        boot = np.median(lr[idx], axis=1)
        out["ci_method"] = "targets"
    else:
        g = np.asarray(g)[ok]
        labels = np.unique(g)
        per = [lr[g == lab] for lab in labels]
        boot = np.array([np.median(np.concatenate([per[j] for j in rng.integers(0, len(per), len(per))])) for _ in range(n_boot)])
        out["ci_method"] = f"parcels ({len(labels)})"
        out["share_parcels_opm_better"] = float(np.mean([np.median(x) > 0 for x in per]))
    out["ci95"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    return out


def binned(values, x, edges, min_n=10, stat=np.median):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x < b) & np.isfinite(values)
        out.append(dict(lo=float(a), hi=float(b), n=int(m.sum()),
                        median=float(stat(values[m])) if m.sum() >= min_n else None,
                        q25=float(np.percentile(values[m], 25)) if m.sum() >= min_n else None,
                        q75=float(np.percentile(values[m], 75)) if m.sum() >= min_n else None))
    return out


def heat(values, depth, orient, min_n=10):
    h = np.full((len(DEPTH_EDGES) - 1, len(ORIENT_EDGES) - 1), np.nan)
    r = np.searchsorted(DEPTH_EDGES, depth, side="right") - 1
    c = np.clip(np.searchsorted(ORIENT_EDGES, orient, side="right") - 1, 0, len(ORIENT_EDGES) - 2)
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            m = (r == i) & (c == j) & np.isfinite(values)
            if m.sum() >= min_n:
                h[i, j] = np.median(values[m])
    return h


# ----------------------------------------------------------------------------------------------
class Study:
    """Everything shared by the primary analysis and the variants."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg["sources"]["seed"])
        self.subject = anatomy.load_sample()
        self.cortex = anatomy.full_resolution(self.subject)
        self.dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
        b = cfg["band"]
        raw_sf = self.dig["sfreq"]
        self.filt = noise.AnalysisFilter(fs=raw_sf, l_freq=b["l_freq_hz"], h_freq=b["h_freq_hz"], order=b["order"])
        self.enbw = self.filt.enbw()
        self.bandwidth = b["h_freq_hz"] - b["l_freq_hz"]
        self.q = cfg["sources"]["focal_nAm"] * 1e-9
        self.opm_asd = cfg["sensors"]["opm_asd_primary_fT_per_rtHz"] * 1e-15
        self.src = g2.make_sources(self.subject, self.cortex, self.rng)
        self.points = np.concatenate([self.src.target, self.src.grid])
        self.nt = len(self.src.target)

    # gains for targets and the background grid
    def gains(self, array, fullres=True, **kw):
        job = g2.FULLRES_JOBS.get(array.name) if fullres else None
        g = g2.gains(array, self.subject, self.cortex, self.points, job, **kw)
        return g[:, :self.nt], g[:, self.nt:]

    def noise(self, array, g_grid, brain_scale, env, asd=None, corr=None):
        corr_arg = None if not corr else (self.cortex.rr[self.src.grid], corr)
        return g2.array_noise(array, g_grid, self.src.grid_area, brain_scale, env, self.enbw,
                              self.opm_asd if asd is None else asd, corr_arg)

    def unit_brain(self, g_grid, corr=None):
        mc = (background.moment_covariance(self.src.grid_area) if not corr else
              background.moment_covariance(self.src.grid_area, self.cortex.rr[self.src.grid], corr))
        return background.sensor_covariance(g_grid, mc)


def calibrate(study, squid, g_grid_squid, target_var, grads, corr=None):
    return background.calibrate(study.unit_brain(g_grid_squid, corr), grads, target_var)


def evaluate_all(topo, arr, nz, conditions, cov_est=None):
    """{(chset, condition): metrics} for one array."""
    out = {}
    for cond in conditions:
        ce = None if cov_est is None else cov_est.get(cond)
        for cs, m in g2.channel_sets(arr).items():
            out[(cs, cond)] = g2.evaluate(topo, nz, m, cond, ce)
    return out


def detect_of(res, name, ref_or_cs, cond, key="detect"):
    cs = "opm" if name != "squid" else ref_or_cs
    return res[name][(cs, cond)][key]


def comparisons(study, res, conditions, key="detect", n_boot=1000):
    out = {}
    for cond in conditions:
        for a in OPMS:
            if a not in res:
                continue
            for ref in REFS:
                out[f"{a}/{ref}/{cond}"] = compare(detect_of(res, a, None, cond, key), detect_of(res, "squid", ref, cond, key),
                                                   study.rng, n_boot)
    return out


# ----------------------------------------------------------------------------------------------
def main():
    t_start = time.time()
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    mne.set_log_level("WARNING")
    st = Study(cfg)
    global GROUPS
    GROUPS = st.src.region
    conds, headline = cfg["conditions"]["all"], cfg["conditions"]["headline"]
    bads = cfg["sensors"]["bads"]
    log(f"targets {st.nt}, background grid {len(st.src.grid)} sources, ENBW {st.enbw:.2f} Hz")

    # --- arrays (primary: measured head position, no extra scalp gap) ----------------------------
    arrays = g2.build_arrays(st.subject, st.dig)
    squid = arrays["squid"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads, mags = good & (squid.kinds == "grad"), good & (squid.kinds == "mag")
    geometry = {name: {**{k: v for k, v in a.meta.items() if not isinstance(v, (list, np.ndarray)) or k == "excluded_sites"},
                       "channels": a.n} for name, a in arrays.items()}

    # --- measured noise, environment, background calibration ------------------------------------
    meas = g2.measured_noise(squid.info, st.filt, bads)
    env = meas["environment"]
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        G_t[name], G_g[name] = st.gains(a)
        log(f"gains {name}: {G_t[name].shape[0]} channels")
    target_var = float(np.nanmedian(meas["brain"][grads]))
    brain_scale = calibrate(st, squid, G_g["squid"], target_var, grads)
    unit_sq = st.unit_brain(G_g["squid"])
    model_brain = brain_scale * np.diag(unit_sq)
    noise_nom = {name: st.noise(a, G_g[name], brain_scale, env) for name, a in arrays.items()}
    er_model = noise_nom["squid"].intrinsic_var + np.diag(noise_nom["squid"].env_cov)
    validation = dict(
        brain_scale=brain_scale,
        brain_source_density_nAm_per_sqrt_mm2=float(np.sqrt(brain_scale) * 1e9 * 1e-3),
        measured=dict(empty_room_rms_mag_fT=float(np.sqrt(np.nanmedian(meas["empty_room"][mags])) * 1e15),
                      empty_room_rms_grad_fT_cm=float(np.sqrt(np.nanmedian(meas["empty_room"][grads])) * 1e13),
                      baseline_rms_mag_fT=float(np.sqrt(np.nanmedian(meas["baseline"][mags])) * 1e15),
                      baseline_rms_grad_fT_cm=float(np.sqrt(np.nanmedian(meas["baseline"][grads])) * 1e13),
                      brain_rms_mag_fT=float(np.sqrt(np.nanmedian(meas["brain"][mags])) * 1e15),
                      brain_rms_grad_fT_cm=float(np.sqrt(target_var) * 1e13),
                      negative_brain_variance_channels=int(np.sum(meas["brain"][good] < 0)),
                      n_baseline_windows=int(meas["n_windows"]), n_baseline_samples=int(meas["baseline_n_samples"]),
                      n_empty_room_samples=int(meas["empty_room_n_samples"])),
        model=dict(brain_rms_mag_fT=float(np.sqrt(np.median(model_brain[mags])) * 1e15),
                   brain_rms_grad_fT_cm=float(np.sqrt(np.median(model_brain[grads])) * 1e13),
                   brain_mag_model_over_measured=float(np.sqrt(np.median(model_brain[mags]) / np.nanmedian(meas["brain"][mags]))),
                   empty_room_rms_mag_fT=float(np.sqrt(np.median(er_model[mags])) * 1e15),
                   empty_room_rms_grad_fT_cm=float(np.sqrt(np.median(er_model[grads])) * 1e13),
                   intrinsic_rms_mag_fT=float(np.sqrt(3.5e-15**2 * st.enbw) * 1e15),
                   intrinsic_rms_grad_fT_cm=float(np.sqrt(3.6e-13**2 * st.enbw) * 1e13),
                   intrinsic_rms_opm_fT={f"{a:g}": float(a * np.sqrt(st.enbw)) for a in cfg["sensors"]["opm_asd_fT_per_rtHz"]}),
        environment_explained_fraction=env.explained_fraction,
        per_channel_brain_ratio_model_over_measured=dict(
            grad=dict(zip(("p5", "median", "p95"), np.nanpercentile(model_brain[grads] / meas["brain"][grads], [5, 50, 95]).tolist())),
            mag=dict(zip(("p5", "median", "p95"), np.nanpercentile(model_brain[mags] / np.where(meas["brain"][mags] > 0, meas["brain"][mags], np.nan), [5, 50, 95]).tolist()))),
    )
    log(f"brain noise: measured grad {validation['measured']['brain_rms_grad_fT_cm']:.1f} fT/cm (calibrated), mag measured "
        f"{validation['measured']['brain_rms_mag_fT']:.0f} fT vs model {validation['model']['brain_rms_mag_fT']:.0f} fT")

    # per-array rank and noise composition per condition
    ranks = {}
    for name, nz in noise_nom.items():
        for cond in conds:
            c = nz.covariance(cond)
            ranks[f"{name}/{cond}"] = int(metrics.whitener(c).rank)
            if name == "squid":  # channel subsets after the full-array projection
                for k in ("mag", "grad"):
                    m = squid.kinds == k
                    ranks[f"squid_{k}/{cond}"] = int(metrics.whitener(c[np.ix_(m, m)]).rank)
    comp = {name: dict(intrinsic_rms=float(np.sqrt(np.median(nz.intrinsic_var))), brain_rms=float(np.sqrt(np.median(np.diag(nz.brain_cov)))),
                       env_rms=float(np.sqrt(np.median(np.diag(nz.env_cov)))))
            for name, nz in noise_nom.items() if name != "squid"}
    for k, m in (("mag", squid.kinds == "mag"), ("grad", squid.kinds == "grad")):
        nz = noise_nom["squid"]
        comp[f"squid_{k}"] = dict(intrinsic_rms=float(np.sqrt(np.median(nz.intrinsic_var[m]))),
                                  brain_rms=float(np.sqrt(np.median(np.diag(nz.brain_cov)[m]))),
                                  env_rms=float(np.sqrt(np.median(np.diag(nz.env_cov)[m]))))

    # --- primary evaluation: oracle and plug-in covariances --------------------------------------
    n_est = {f"T{int(t)}": int(2 * st.bandwidth * t) for t in cfg["covariance"]["estimate_seconds"]}
    res = {}
    for name, a in arrays.items():
        cov_est = {cond: {lab: g2.sample_covariance(noise_nom[name].covariance(cond), n, st.rng) for lab, n in n_est.items()}
                   for cond in conds}
        res[name] = evaluate_all(G_t[name] * st.q, a, noise_nom[name], conds, cov_est)
        log(f"evaluated {name}")
    primary = dict(oracle=comparisons(st, res, conds))
    for lab in n_est:
        primary[f"plugin_{lab}"] = comparisons(st, res, conds, key=f"detect_{lab}")
    for key in ("peak", "meanpow_db"):
        if key == "meanpow_db":  # convert to a linear amplitude ratio for the same log2 summary
            lin = {n: {k: dict(v, lin=10 ** (v["meanpow_db"] / 20)) for k, v in r.items()} for n, r in res.items()}
            primary[key] = comparisons(st, lin, conds, key="lin")
        else:
            primary[key] = comparisons(st, res, conds, key="peak")
    plugin_loss = {f"{name}/{cs}/{cond}/{lab}": float(np.median(r[f"detect_{lab}"] / r["detect"]))
                   for name, rr in res.items() for (cs, cond), r in rr.items() for lab in n_est}

    # raw amplitudes (T, T/m) and absolute detectability vs depth
    depth, lobe = st.src.depth_mm, st.src.lobe
    amp = {"squid_mag": np.abs(G_t["squid"][squid.kinds == "mag"] * st.q).max(axis=0),
           "squid_grad": np.abs(G_t["squid"][squid.kinds == "grad"] * st.q).max(axis=0)}
    for a in OPMS:
        amp[a] = np.abs(G_t[a] * st.q).max(axis=0)
    amp_vs_depth = {k: binned(v, depth, DEPTH_EDGES) for k, v in amp.items()}
    snr_vs_depth = {f"{n}/{cs}/{cond}": binned(r[(cs, cond)]["detect"], depth, DEPTH_EDGES)
                    for n, r in res.items() for (cs, cond) in r}
    ratio_vs_depth = {f"{a}/{ref}/{cond}": binned(np.log2(detect_of(res, a, None, cond) / detect_of(res, "squid", ref, cond)), depth, DEPTH_EDGES)
                      for a in OPMS for ref in REFS for cond in conds}
    # the projection's cost: noise-normalised signal norm retained, detectability projected / with the room field
    projection = {}
    for name, a in arrays.items():
        nz = noise_nom[name]
        s_ = G_t[name] * st.q
        d_ = 1.0 / np.sqrt(nz.intrinsic_var)[:, None]
        kept = np.linalg.norm(d_ * (nz.projector() @ s_), axis=0) / np.linalg.norm(d_ * s_, axis=0)
        cs = "combined" if name == "squid" else "opm"
        projection[name] = dict(retained_signal_norm_median=float(np.median(kept)),
                                detect_projected_over_env_median=float(np.median(res[name][(cs, "projected")]["detect"]
                                                                                 / res[name][(cs, "intrinsic+brain+env")]["detect"])))
    # intrinsic noise only: the ratio scales exactly as 1 / OPM ASD, so the break-even noise level follows
    break_even = {f"{a}/{ref}": float(st.opm_asd * 1e15 * 2 ** primary["oracle"][f"{a}/{ref}/intrinsic"]["median_log2"])
                  for a in OPMS for ref in REFS}
    strata = {}
    for a in OPMS:
        for cond in headline:
            lr = np.log2(detect_of(res, a, None, cond) / detect_of(res, "squid", "combined", cond))
            for lab, x, edges in (("dist_inner_skull_mm", st.src.dist_inner_skull_mm, (4.0, 5.0, 6.0, 10.0, 100.0)),
                                  ("orientation_deg", st.src.orientation_deg, (0.0, 30.0, 60.0, 90.1))):
                strata[f"{a}/combined/{cond}/{lab}"] = [dict(lo=lo, hi=hi, n=int(np.sum((x >= lo) & (x < hi))),
                                                             ratio=float(2 ** np.median(lr[(x >= lo) & (x < hi)])))
                                                        for lo, hi in zip(edges[:-1], edges[1:])]
    # the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain, not cortex) holds some
    # usable targets; they enter every median above, so the headline ratios are also given without them
    medial = np.char.endswith(st.src.region.astype(str), "unknown")
    medial_wall = dict(n_targets=int(medial.sum()), share=float(medial.mean()), ratios_without={
        f"{a}/combined/{cond}": float(2 ** np.median(np.log2(detect_of(res, a, None, cond)[~medial]
                                                             / detect_of(res, "squid", "combined", cond)[~medial])))
        for a in OPMS for cond in headline})
    bridge = bridge_to_sphere(st, amp, depth, geometry_distances(st, arrays))
    by_lobe = {}
    for a in OPMS:
        for ref in REFS:
            for cond in headline:
                lr = np.log2(detect_of(res, a, None, cond) / detect_of(res, "squid", ref, cond))
                by_lobe[f"{a}/{ref}/{cond}"] = {lb: compare(2 ** lr[lobe == lb], np.ones((lobe == lb).sum()), st.rng,
                                                        groups=st.src.region[lobe == lb]) for lb in plotting.DK_LOBES}
    log("primary comparisons done")

    # --- extended sources (full-resolution lead fields) -----------------------------------------
    patches = patch_analysis(st, arrays, noise_nom, headline, cfg)
    log("patches done")

    # --- sensitivity analyses -------------------------------------------------------------------
    sens = {}
    # OPM intrinsic noise
    for asd in cfg["sensors"]["opm_asd_fT_per_rtHz"]:
        r2 = {"squid": {k: v for k, v in res["squid"].items() if k[1] in headline}}
        for a in OPMS:
            nz = st.noise(arrays[a], G_g[a], brain_scale, env, asd=asd * 1e-15)
            r2[a] = evaluate_all(G_t[a] * st.q, arrays[a], nz, headline)
        sens[f"opm_asd_{asd:g}fT"] = comparisons(st, r2, headline, n_boot=200)
    log("sensitivity: OPM noise done")
    # correlated background
    for lam in cfg["background"]["correlation_lengths_mm"]:
        if not lam:
            continue
        bs = calibrate(st, squid, G_g["squid"], target_var, grads, corr=lam * 1e-3)
        r2 = {}
        for name, a in arrays.items():
            nz = st.noise(a, G_g[name], bs, env, corr=lam * 1e-3)
            r2[name] = evaluate_all(G_t[name] * st.q, a, nz, headline)
        sens[f"background_corr_{lam:g}mm"] = comparisons(st, r2, headline, n_boot=200)
    log("sensitivity: correlated background done")
    # background calibrated on the magnetometers instead (the gradiometer-calibrated model predicts
    # ~0.7x the measured magnetometer brain-noise amplitude)
    bs_mag = background.calibrate(unit_sq, mags, float(np.nanmedian(meas["brain"][mags])))
    r2 = {}
    for name, a in arrays.items():
        nz = st.noise(a, G_g[name], bs_mag, env)
        r2[name] = evaluate_all(G_t[name] * st.q, a, nz, headline)
    sens["background_mag_calibrated"] = comparisons(st, r2, headline, n_boot=200)
    sens["background_mag_calibrated"]["scale_over_primary"] = float(bs_mag / brain_scale)
    log("sensitivity: magnetometer-calibrated background done")
    # head position (SQUID only; OPMs are head-mounted)
    hp = g2.head_position_variants(squid.info, st.subject, cfg["head_position"]["translation_mm"] * 1e-3,
                                   cfg["head_position"]["pitch_deg"], cfg["head_position"]["fit_clearance_mm"] * 1e-3,
                                   cfg["head_position"]["dewar_spacing_mm"] * 1e-3)
    head_positions = {}
    for vname, v in hp.items():
        head_positions[vname] = dict(min_dist_mm=v["min_dist"] * 1e3, median_dist_mm=v["median_dist"] * 1e3, feasible=v["feasible"])
        if vname == "measured":
            continue
        info_v = squid.info.copy()
        with info_v._unlock():
            info_v["dev_head_t"] = mne.transforms.Transform("meg", "head", v["trans"])
        arr_v = g2.Array("squid", info_v, squid.kinds, None, dict(squid.meta))
        gt, gg = st.gains(arr_v, fullres=False)
        nz = st.noise(arr_v, gg, brain_scale, env)
        r2 = {"squid": evaluate_all(gt * st.q, arr_v, nz, headline)}
        for a in OPMS:
            r2[a] = {k: v2 for k, v2 in res[a].items() if k[1] in headline}
        sens[f"head_{vname}"] = comparisons(st, r2, headline, n_boot=200)
        log(f"sensitivity: head position {vname} done")
    # scalp gap (OPM only), and the joint OPM noise x scalp gap grid (intrinsic+brain, not one factor at a time)
    joint = {}
    sq_ref = {"squid": {k: v2 for k, v2 in res["squid"].items() if k[1] == "intrinsic+brain"}}
    for gap in cfg["sensors"]["opm_scalp_gap_mm"]:
        r2 = {"squid": {k: v2 for k, v2 in res["squid"].items() if k[1] in headline}}
        gap_geo, gap_gains = {}, {}
        for a in OPMS:
            if gap:
                arr = g2.with_scalp_gap(arrays[a], gap * 1e-3)  # same sites, moved outward (not rebuilt)
                gt, gg = st.gains(arr, fullres=False)
                nz = st.noise(arr, gg, brain_scale, env)
                r2[a] = evaluate_all(gt * st.q, arr, nz, headline)
                gap_geo[a] = arr.n
            else:
                arr, gt, gg = arrays[a], G_t[a], G_g[a]
            gap_gains[a] = (arr, gt, gg)
        if gap:
            sens[f"gap_{gap:g}mm"] = comparisons(st, r2, headline, n_boot=200)
            sens[f"gap_{gap:g}mm"]["channels"] = gap_geo
        for asd in (15.0, 20.0, 30.0):
            r3 = dict(sq_ref)
            for a in ("opm_matched", "opm_dense"):
                arr, gt, gg = gap_gains[a]
                r3[a] = evaluate_all(gt * st.q, arr, st.noise(arr, gg, brain_scale, env, asd=asd * 1e-15), ["intrinsic+brain"])
            c3 = comparisons(st, r3, ["intrinsic+brain"], n_boot=200)
            joint[f"gap{gap:g}mm/asd{asd:g}fT"] = {a: dict(ratio=float(2 ** c3[f"{a}/combined/intrinsic+brain"]["median_log2"]),
                                                          ci95=[float(2 ** x) for x in c3[f"{a}/combined/intrinsic+brain"]["ci95"]])
                                                   for a in ("opm_matched", "opm_dense")}
        log(f"sensitivity: gap {gap:g} mm done")
    sens_joint = joint

    # --- convergence -----------------------------------------------------------------------------
    conv = convergence(st, arrays, G_t, G_g, noise_nom, res, brain_scale, target_var, grads, env, headline, cfg)
    log("convergence done")

    # --- outputs (checkpointed: `--replot` redraws figures and rewrites outputs without recomputing)
    patch_det = patches.pop("_det", None)
    summary = dict(
        status="NEW (realistic adult OPM vs Neuromag comparison, MNE sample subject)",
        config=cfg, arrays=geometry, head_positions=head_positions, n_targets=st.nt, n_background_grid=int(len(st.src.grid)),
        enbw_hz=st.enbw, n_estimate_samples=n_est, noise_validation=validation, noise_composition=comp, retained_rank=ranks,
        primary=primary, plugin_over_oracle_median=plugin_loss, amplitude_vs_depth=amp_vs_depth, detectability_vs_depth=snr_vs_depth,
        log2_ratio_vs_depth=ratio_vs_depth, by_lobe=by_lobe, strata=strata, medial_wall=medial_wall, projection=projection,
        break_even_opm_asd_fT=break_even,
        bridge_to_sphere=bridge, patches=patches, sensitivity=sens, sensitivity_joint_asd_gap=sens_joint, convergence=conv,
        runtime_s=time.time() - t_start,
        notes=["Bootstrap CIs resample Desikan-Killiany parcels of one anatomy (targets within a parcel are correlated); they do "
               "not include model or between-subject uncertainty. Sensitivity analyses vary one factor at a time except the "
               "joint OPM noise x scalp gap grid; they show the dependence, they do not bound it.",
               "The 'projected' condition removes the simulated room field exactly (it lies in the removed 8-dim subspace); it "
               "measures the projection's cost (rank, signal attenuation), not residual interference.",
               "Head-position variants move the head in the fixed Neuromag helmet; the room field is kept in head coordinates "
               "(8-term model; the change over a 5-mm move is second order).",
               "Detectability is a known-topography matched-filter SNR, not an event detection rate or localization accuracy."])
    summary["provenance"] = dict(commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__)
    state = dict(summary=summary, arrays=arrays, res=res, amp=amp, patch_det=patch_det)
    STATE.parent.mkdir(parents=True, exist_ok=True)
    with open(STATE, "wb") as fh:
        pickle.dump(state, fh)
    outputs(st, state)
    log(f"done in {time.time() - t_start:.0f} s")


def outputs(st, state):
    s = state["summary"]
    cfg = s["config"]
    figures(st, state["arrays"], state["res"], state["amp"], st.src.depth_mm, st.src.orientation_deg, st.src.lobe, s["patches"],
            s["sensitivity"], s["primary"], cfg["conditions"]["headline"], cfg["conditions"]["all"], s["bridge_to_sphere"],
            state.get("patch_det"))
    write_targets_csv(st, state["res"], state["amp"], cfg["conditions"]["all"])
    if state.get("patch_det"):
        keys = sorted(k for k in state["patch_det"] if k[4] == "fixed_total")
        with open(OUT / "g2_patch_targets.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["hemi", "vertno", "depth_mm"] + [f"detect_{n}_{cs}_{cond}_{r:g}mm_fixed_total" for n, cs, cond, r, _ in keys])
            for i, v in enumerate(st.src.target):
                wr.writerow([int(st.cortex.hemi[v]), int(st.cortex.vertno[v]), f"{st.src.depth_mm[i]:.2f}"]
                            + [f"{state['patch_det'][k][i]:.4f}" for k in keys])
    io.write_json(s, OUT / "g2_summary.json")


# ----------------------------------------------------------------------------------------------
def geometry_distances(st, arrays):
    """Scalp-to-sensor distances [m]: Neuromag magnetometer coils (Euclidean to the nearest scalp
    point) and OPM sensing centres, head frame."""
    from scipy.spatial import cKDTree

    mri_head = np.linalg.inv(st.subject.trans["trans"])
    tree = cKDTree(st.subject.scalp.rr @ mri_head[:3, :3].T + mri_head[:3, 3])
    out = {}
    for name, a in arrays.items():
        pos = np.array([c["loc"][:3] for c in a.info["chs"]])
        if name == "squid":
            t = a.info["dev_head_t"]["trans"]
            pos = (pos @ t[:3, :3].T + t[:3, 3])[a.kinds == "mag"]
        out[name] = tree.query(pos)[0]
    return out


def bridge_to_sphere(st, amp, depth, dist):
    """Realistic vs idealized benchmark. In the sphere (Jas Eq. 1-3) the OPM/SQUID peak-field ratio
    falls with depth and the equal-SNR depth solves ratio = eta. The realistic analogue is the
    median over targets of max|B_OPM| / max|B_mag| (10-nAm cortical-normal dipoles, 3-layer BEM,
    real sensor positions and orientations), and d_eq(eta) is where that median crosses eta
    (intrinsic white noise only, peak-channel SNR, as in the benchmark)."""
    from opmsquid import sphere

    centers = 0.5 * (DEPTH_EDGES[:-1] + DEPTH_EDGES[1:])
    out = dict(sensor_distance_mm={k: dict(median=float(np.median(v) * 1e3), p5=float(np.percentile(v, 5) * 1e3),
                                           p95=float(np.percentile(v, 95) * 1e3)) for k, v in dist.items()})
    xi_sq = float(np.median(dist["squid"]))
    etas = np.round(np.arange(1.0, 6.01, 0.25), 2)
    sph = {"jas_xi0_18": (0.0, 0.018), "realistic_standoffs": (float(np.median(dist["opm_matched"])), xi_sq)}
    out["sphere_ratio_vs_depth"] = {k: [float(x) for x in sphere.signal_ratio(centers * 1e-3, 0.095, *xi)] for k, xi in sph.items()}
    out["sphere_d_eq_mm"] = {k: {f"{e:g}": float(sphere.equal_snr_depth(e, 0.095, 0.080, *xi) * 1e3) for e in etas} for k, xi in sph.items()}
    out["depth_centers_mm"] = centers.tolist()
    for a in OPMS:
        ratio = amp[a] / amp["squid_mag"]
        b = binned(ratio, depth, DEPTH_EDGES)
        med = np.array([x["median"] if x["median"] is not None else np.nan for x in b], float)
        out[f"{a}_ratio_vs_depth"] = b
        d_eq, opm_all, squid_all = {}, [], []
        for e in etas:  # first depth where the median ratio falls below eta (linear interpolation)
            ok = np.isfinite(med)
            c, m = centers[ok], med[ok]
            below = np.flatnonzero(m < e)
            d_eq[f"{e:g}"] = None
            if len(below) == 0:
                opm_all.append(float(e))  # OPM ahead in every depth bin: no crossing
            elif below[0] == 0:
                squid_all.append(float(e))  # SQUID ahead already in the shallowest bin: no crossing
            else:
                i = below[0]
                d_eq[f"{e:g}"] = float(np.interp(e, [m[i], m[i - 1]], [c[i], c[i - 1]]))
        out[f"{a}_d_eq_mm"] = d_eq
        out[f"{a}_eta_opm_ahead_at_all_depths"] = opm_all
        out[f"{a}_eta_squid_ahead_at_all_depths"] = squid_all
    return out


# ----------------------------------------------------------------------------------------------
def patch_analysis(st, arrays, noise_nom, headline, cfg):
    """Geodesic patches around every target: detectability under the fixed-total and fixed-density
    conventions (OPM/SQUID ratios do not depend on the convention; absolute values do)."""
    radii = cfg["sources"]["patch_radii_mm"]
    members = {r: goldenholz.geodesic_patches(st.cortex.adjacency, st.src.target, r * 1e-3, st.cortex.usable) for r in radii}
    area = {r: np.array([st.cortex.area[m].sum() for m in members[r]]) for r in radii}
    dens = cfg["sources"]["patch_density_nAm_per_mm2"] * 1e-9 / 1e-6
    total = cfg["sources"]["patch_total_nAm"] * 1e-9
    out = dict(radii_mm=radii, area_cm2={f"{r:g}": dict(median=float(np.median(area[r]) * 1e4), p5=float(np.percentile(area[r], 5) * 1e4),
                                                        p95=float(np.percentile(area[r], 95) * 1e4)) for r in radii})
    det = {}
    for name, a in arrays.items():
        full, col = g2.fullres_matrix(a, st.subject, st.cortex, g2.FULLRES_JOBS[name])
        for r in radii:
            unit = goldenholz.patch_topographies(full, members[r], col, st.cortex.area)  # per unit density (A m / m^2)
            for conv_name, scale in (("fixed_density", np.full(len(area[r]), dens)), ("fixed_total", total / area[r])):
                ev = evaluate_all(unit * scale[None, :], a, noise_nom[name], headline)
                for (cs, cond), v in ev.items():
                    det[(name, cs, cond, r, conv_name)] = v["detect"]
            # net moment relative to the scalar moment (cancellation), per unit density
            if name == "squid":
                out.setdefault("net_over_scalar_moment", {})[f"{r:g}"] = None
        del full
    # cancellation measure independent of arrays: |sum of oriented moments| / sum |moments|
    for r in radii:
        ratio = np.array([np.linalg.norm((st.cortex.nn[m] * st.cortex.area[m][:, None]).sum(axis=0)) / st.cortex.area[m].sum()
                          for m in members[r]])
        out["net_over_scalar_moment"][f"{r:g}"] = dict(median=float(np.median(ratio)), p5=float(np.percentile(ratio, 5)),
                                                       p95=float(np.percentile(ratio, 95)))
    comp, absd = {}, {}
    rng = np.random.default_rng(1)
    for r in radii:
        for cond in headline:
            for a in OPMS:
                for ref in REFS:
                    comp[f"{a}/{ref}/{cond}/{r:g}mm"] = compare(det[(a, "opm", cond, r, "fixed_total")],
                                                                 det[("squid", ref, cond, r, "fixed_total")], rng, 200)
            for conv_name in ("fixed_total", "fixed_density"):
                for name, cs in [("squid", c) for c in REFS] + [(a, "opm") for a in OPMS]:
                    absd[f"{name}/{cs}/{cond}/{r:g}mm/{conv_name}"] = float(np.median(det[(name, cs, cond, r, conv_name)]))
    out["comparisons"] = comp
    out["median_detectability"] = absd
    out["_det"] = det  # removed before writing
    return out


def convergence(st, arrays, G_t, G_g, noise_nom, res, brain_scale, target_var, grads, env, headline, cfg):
    out = {}
    rng = np.random.default_rng(2)
    squid = arrays["squid"]

    def ratio_medians(r2):
        return {f"{a}/{ref}/{cond}": float(np.median(np.log2(detect_of(r2, a, None, cond) / detect_of(r2, "squid", ref, cond))))
                for a in OPMS if a in r2 for ref in REFS for cond in headline}

    base = ratio_medians(res)
    out["primary_median_log2"] = base

    # (a) background discretisation: 7-mm grid vs every full-resolution vertex (area-weighted)
    full_cov = {}
    for name, a in arrays.items():
        full, col = g2.fullres_matrix(a, st.subject, st.cortex, g2.FULLRES_JOBS[name])
        use = np.flatnonzero((col >= 0) & st.cortex.usable)  # A-BEM-DIST
        c = np.zeros((full.shape[0], full.shape[0]))
        for s in range(0, len(use), 20000):
            idx = use[s:s + 20000]
            g = np.asarray(full[:, col[idx]], dtype=np.float64)
            c += (g * st.cortex.area[idx][None, :]) @ g.T
        full_cov[name] = c
        del full
    bs_full = background.calibrate(full_cov["squid"], grads, target_var)
    r2, var_ratio = {}, {}
    for name, a in arrays.items():
        nz0 = noise_nom[name]
        nz = noisemodel.ArrayNoise(nz0.intrinsic_var, bs_full * full_cov[name], nz0.env_cov, nz0.ext_basis)
        r2[name] = evaluate_all(G_t[name] * st.q, a, nz, headline)
        var_ratio[name] = float(np.median(np.diag(nz.brain_cov) / np.diag(nz0.brain_cov)))
    d = np.sqrt(np.diag(full_cov["squid"]))
    d0 = np.sqrt(np.diag(noise_nom["squid"].brain_cov))
    corr_diff = np.linalg.norm(full_cov["squid"] / np.outer(d, d) - noise_nom["squid"].brain_cov / np.outer(d0, d0)) / \
        np.linalg.norm(noise_nom["squid"].brain_cov / np.outer(d0, d0))
    out["background_grid_vs_fullres"] = dict(scale_ratio_full_over_grid=float(bs_full / brain_scale), median_variance_ratio=var_ratio,
                                             squid_correlation_rel_frobenius_diff=float(corr_diff), median_log2=ratio_medians(r2),
                                             max_abs_change_log2=float(max(abs(v - base[k]) for k, v in ratio_medians(r2).items())))

    # subset of targets for the forward-model checks
    sub = np.sort(rng.choice(st.nt, cfg["convergence"]["subset_targets"], replace=False))
    pts = np.concatenate([st.src.target[sub], st.src.grid])
    ns = len(sub)

    def run_forward_variant(label, gain_fn):
        gt, gg = {}, {}
        for name in ("squid", "opm_matched", "opm_dense"):
            g = gain_fn(arrays[name])
            gt[name], gg[name] = g[:, :ns], g[:, ns:]
        bs = background.calibrate(st.unit_brain(gg["squid"]), grads, target_var)
        r = {}
        for name in gt:
            nz = st.noise(arrays[name], gg[name], bs, env)
            r[name] = evaluate_all(gt[name] * st.q, arrays[name], nz, headline)
        return gt, r

    def direct(bem_surfs, coil_def=None, info=None):
        return lambda arr: forward.chunked_discrete_gain(info or arr.info, st.subject.trans, st.cortex.rr[pts], st.cortex.nn[pts],
                                                         bem_surfs, coil_def=coil_def or opm.coil_def_file()).astype(np.float64)

    ref_gt, ref_r = run_forward_variant("bem3_primary", direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY)))
    ref_med = ratio_medians(ref_r)
    variants = {"bem3_head5120": direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=0)),  # the v1 head surface
                "bem1_5120": direct(st.subject.bem_model((0.3,))),
                "bem1_20480": direct(anatomy.refined_inner_skull(st.subject, 1))}
    fwd_checks = dict(subset_reference_median_log2=ref_med)
    for label, fn in variants.items():
        gt, r = run_forward_variant(label, fn)
        med = ratio_medians(r)
        fwd_checks[label] = dict(median_log2=med, max_abs_change_log2=float(max(abs(med[k] - ref_med[k]) for k in med)),
                                 gain_rel_diff_median={n: float(np.median(np.linalg.norm(gt[n] - ref_gt[n], axis=0) / np.linalg.norm(ref_gt[n], axis=0)))
                                                       for n in gt})
        log(f"convergence: {label} done")
    out["bem"] = fwd_checks
    # 1-layer 5120 vs 20480 isolates the discretisation (same shape)
    g1 = fwd_checks["bem1_5120"]["median_log2"]
    g2_ = fwd_checks["bem1_20480"]["median_log2"]
    out["bem"]["refinement_max_abs_change_log2"] = float(max(abs(g1[k] - g2_[k]) for k in g1))

    # coil integration: Neuromag 4-point ('normal') vs accurate; OPM point vs 10-mm cell
    coil = {}
    info4 = neuromag.load_info("T3-4pt")
    g_acc = direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY))(squid)
    g_4pt = direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY), info=info4)(squid)
    coil["squid_4pt_vs_accurate_rel_diff_median"] = float(np.median(np.linalg.norm(g_4pt - g_acc, axis=0) / np.linalg.norm(g_acc, axis=0)))
    for name in ("opm_matched", "opm_dense"):
        gp = direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY), coil_def=opm.coil_def_file(cell_size=1e-6))(arrays[name])
        gc = direct(st.subject.bem_model(g2.BEM_CONDUCTIVITY))(arrays[name])
        rel = np.abs(np.abs(gc[:, :ns]).max(axis=0) / np.abs(gp[:, :ns]).max(axis=0) - 1)
        coil[f"{name}_cell_vs_point_peak_rel_diff"] = dict(median=float(np.median(rel)), p95=float(np.percentile(rel, 95)), max=float(rel.max()))
    out["coil_integration"] = coil
    log("convergence: coil checks done")

    # target sampling: oct-6 targets vs the same number of random full-resolution vertices
    valid = np.flatnonzero(st.cortex.usable)
    rand = np.sort(rng.choice(valid, st.nt, replace=False))
    r3 = {}
    for name, a in arrays.items():
        full, col = g2.fullres_matrix(a, st.subject, st.cortex, g2.FULLRES_JOBS[name])
        gt = np.asarray(full[:, col[rand]], dtype=np.float64)
        del full
        r3[name] = evaluate_all(gt * st.q, a, noise_nom[name], headline)
    med3 = ratio_medians(r3)
    out["target_sampling"] = dict(random_fullres_median_log2=med3, max_abs_change_log2=float(max(abs(med3[k] - base[k]) for k in med3)))

    # whitening tolerance
    tol = {}
    for name in ("squid", "opm_dense"):
        for cond in headline:
            nz = noise_nom[name]
            c = nz.covariance(cond)
            s = nz.signal(G_t[name][:, :500] * st.q, cond)
            a_ = metrics.detectability(s, w=metrics.whitener(c, rel_tol=1e-8))
            b_ = metrics.detectability(s, w=metrics.whitener(c, rel_tol=1e-12))
            tol[f"{name}/{cond}"] = float(np.max(np.abs(a_ / b_ - 1)))
    out["whitening_tolerance_max_rel_change"] = tol
    return out


# ----------------------------------------------------------------------------------------------
def write_targets_csv(st, res, amp, conds):
    cols = [(n, cs, cond) for n, r in res.items() for (cs, cond) in r]
    with open(OUT / "g2_targets.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["hemi", "vertno", "depth_mm", "orientation_deg", "region", "lobe"] + [f"amp_{k}" for k in amp]
                    + [f"detect_{n}_{cs}_{cond}" for n, cs, cond in cols])
        for i, v in enumerate(st.src.target):
            wr.writerow([int(st.cortex.hemi[v]), int(st.cortex.vertno[v]), f"{st.src.depth_mm[i]:.2f}", f"{st.src.orientation_deg[i]:.2f}",
                         st.src.region[i], st.src.lobe[i]] + [f"{amp[k][i]:.4e}" for k in amp] + [f"{res[n][(cs, cond)]['detect'][i]:.4f}" for n, cs, cond in cols])


def figures(st, arrays, res, amp, depth, orient, lobe, patches, sens, primary, headline, conds, bridge, patch_det=None):
    # 1. amplitude and detectability vs depth
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    centers = 0.5 * (DEPTH_EDGES[:-1] + DEPTH_EDGES[1:])

    def curve(ax, v, label, **kw):
        b = binned(v, depth, DEPTH_EDGES)
        med = np.array([x["median"] if x["median"] is not None else np.nan for x in b], float)
        lo = np.array([x["q25"] if x["q25"] is not None else np.nan for x in b], float)
        hi = np.array([x["q75"] if x["q75"] is not None else np.nan for x in b], float)
        ax.plot(centers, med, label=label, **kw)
        ax.fill_between(centers, lo, hi, alpha=0.15, color=kw.get("color"))

    colors = {"squid_mag": "tab:blue", "squid_grad": "tab:red", "opm_matched": "tab:green", "opm204": "tab:olive", "opm_dense": "tab:purple",
              "combined": "k", "grad": "tab:red", "mag": "tab:blue"}
    for k in ("squid_mag", "opm_matched", "opm_dense"):
        curve(axs[0], amp[k] * 1e15, {"squid_mag": "Neuromag mag", "opm_matched": LABEL["opm_matched"], "opm_dense": LABEL["opm_dense"]}[k], color=colors[k])
    axs[0].set_yscale("log")
    axs[0].set_ylabel("peak |B| of a 10-nAm dipole [fT]")
    ax2 = axs[0].twinx()
    curve(ax2, amp["squid_grad"] * 1e13, "Neuromag grad [fT/cm]", color=colors["squid_grad"], ls="--")
    ax2.set_yscale("log")
    ax2.set_ylabel("peak |dB/dx| [fT/cm] (dashed)", color="tab:red")
    axs[0].legend(fontsize=7, loc="upper right")
    for ax, cond in zip(axs[1:], headline):
        for n, cs in [("squid", c) for c in REFS] + [(a, "opm") for a in OPMS]:
            key = cs if n == "squid" else n
            curve(ax, res[n][(cs, cond)]["detect"], LABEL[key], color=colors[key])
        ax.set_yscale("log")
        ax.set_title(f"detectability, {cond}", fontsize=9)
        ax.set_ylabel("sqrt(s^T C^+ s) for 10 nAm")
        ax.axhline(1, color="0.6", lw=0.6)
        ax.legend(fontsize=7)
    for ax in axs:
        ax.set_xlabel("depth below scalp [mm]")
    fig.suptitle("G2 (NEW): signal and detectability vs depth (median and IQR over cortical targets; OPM 15 fT/sqrt(Hz))", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_depth.png", dpi=150)
    plt.close(fig)

    # 1b. bridge to the idealized sphere benchmark (G1A)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    c = np.array(bridge["depth_centers_mm"])
    for k, lab, ls in (("jas_xi0_18", "sphere, standoffs 0 / 18 mm (Jas)", "--"),
                       ("realistic_standoffs", "sphere, median real standoffs", ":")):
        axs[0].plot(c, bridge["sphere_ratio_vs_depth"][k], ls, color="k", label=lab)
    for a in ("opm_matched", "opm_dense"):
        b = bridge[f"{a}_ratio_vs_depth"]
        med = np.array([x["median"] if x["median"] is not None else np.nan for x in b], float)
        lo = np.array([x["q25"] if x["q25"] is not None else np.nan for x in b], float)
        hi = np.array([x["q75"] if x["q75"] is not None else np.nan for x in b], float)
        axs[0].plot(c, med, "-", color=colors[a], label=f"{LABEL[a]} / Neuromag mag (median, IQR)")
        axs[0].fill_between(c, lo, hi, color=colors[a], alpha=0.15)
    axs[0].set_yscale("log")
    axs[0].set_xlabel("depth below scalp [mm]")
    axs[0].set_ylabel("peak |B| ratio OPM / SQUID magnetometer")
    axs[0].legend(fontsize=7)
    axs[0].set_title("signal ratio vs depth (equal-SNR where it equals eta)", fontsize=9)
    etas = [float(e) for e in bridge["sphere_d_eq_mm"]["jas_xi0_18"]]
    for k, ls in (("jas_xi0_18", "--"), ("realistic_standoffs", ":")):
        axs[1].plot(etas, [bridge["sphere_d_eq_mm"][k][f"{e:g}"] for e in etas], ls, color="k", label=k.replace("_", " "))
    for a in ("opm_matched", "opm_dense"):
        y = [bridge[f"{a}_d_eq_mm"][f"{e:g}"] for e in etas]
        axs[1].plot(etas, [np.nan if v is None else v for v in y], "o-", ms=3, color=colors[a], label=LABEL[a])
    axs[1].set_xlabel("eta = sigma_OPM / sigma_SQUID-mag (intrinsic noise only)")
    axs[1].set_ylabel("equal-SNR depth [mm]")
    axs[1].legend(fontsize=7)
    axs[1].set_title("equal peak-channel SNR depth: sphere vs realistic head", fontsize=9)
    fig.suptitle("G2 bridge to the idealized benchmark (G1A): same metric, realistic geometry", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_bridge.png", dpi=150)
    plt.close(fig)

    # 2. depth-orientation heatmaps of log2 detectability ratios
    fig, axs = plt.subplots(len(headline), 3, figsize=(13, 4.2 * len(headline)), squeeze=False)
    for i, cond in enumerate(headline):
        for j, a in enumerate(OPMS):
            h = heat(np.log2(res[a][("opm", cond)]["detect"] / res["squid"][("combined", cond)]["detect"]), depth, orient)
            im = axs[i, j].imshow(h, origin="upper", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
                                  extent=[0, 90, DEPTH_EDGES[-1], DEPTH_EDGES[0]])
            axs[i, j].set_ylim(65, 10)
            axs[i, j].set_title(f"log2 {LABEL[a]} / Neuromag combined, {cond}", fontsize=8)
            axs[i, j].set_xlabel("orientation [deg] (0 radial)")
            axs[i, j].set_ylabel("depth [mm]")
            fig.colorbar(im, ax=axs[i, j], fraction=0.046)
    fig.suptitle("G2: median log2 detectability ratio per depth x orientation bin (>= 10 targets)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_heatmaps.png", dpi=150)
    plt.close(fig)
    # the dense array against each Neuromag channel set, and 10-mm patches against the combined system
    panels = [(f"{ref}", lambda cond, ref=ref: np.log2(res["opm_dense"][("opm", cond)]["detect"] / res["squid"][(ref, cond)]["detect"]))
              for ref in ("grad", "mag")]
    if patch_det:
        panels.append(("combined, 10-mm patches", lambda cond: np.log2(patch_det[("opm_dense", "opm", cond, 10.0, "fixed_total")]
                                                                          / patch_det[("squid", "combined", cond, 10.0, "fixed_total")])))
    fig, axs = plt.subplots(len(headline), len(panels), figsize=(4.4 * len(panels), 4.2 * len(headline)), squeeze=False)
    for i, cond in enumerate(headline):
        for j, (lab, fn) in enumerate(panels):
            h = heat(fn(cond), depth, orient)
            im = axs[i, j].imshow(h, origin="upper", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
                                  extent=[0, 90, DEPTH_EDGES[-1], DEPTH_EDGES[0]])
            axs[i, j].set_ylim(65, 10)
            axs[i, j].set_title(f"log2 OPM dense / Neuromag {lab}, {cond}", fontsize=8)
            axs[i, j].set_xlabel("orientation [deg] (0 radial)")
            axs[i, j].set_ylabel("depth [mm]")
            fig.colorbar(im, ax=axs[i, j], fraction=0.046)
    fig.suptitle("G2: dense OPM vs each Neuromag comparator, and extended sources", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_heatmaps_comparators.png", dpi=150)
    plt.close(fig)

    # 3. cortical maps
    hemis = plotting.inflated_views(st.subject.subjects_dir, "sample", st.subject.src)
    n_lh_full = int(np.sum(st.cortex.hemi == 0))
    oct_global = np.concatenate([st.subject.src[0]["vertno"], st.subject.src[1]["vertno"] + n_lh_full])
    ok = st.cortex.usable[oct_global]
    n_lh = st.subject.src[0]["nuse"]
    # the medial wall (FreeSurfer 'unknown') is not cortex: kept in the statistics, shown grey like excluded vertices
    medial_at = np.flatnonzero(ok)[np.char.endswith(st.src.region.astype(str), "unknown")]

    def on_map(values):
        v = np.full(len(oct_global), np.nan)
        v[ok] = values
        v[medial_at] = np.nan
        return v

    grey = "grey: medial wall and vertices < 4 mm from the inner skull"
    for cond in headline:
        rows = []
        for a, ref in (("opm_matched", "combined"), ("opm_matched", "grad"), ("opm_dense", "combined"), ("opm_dense", "grad"),
                       ("opm_dense", "mag")):
            rows.append((f"{LABEL[a]}\nvs {LABEL[ref]}", on_map(np.log2(res[a][("opm", cond)]["detect"] / res["squid"][(ref, cond)]["detect"]))))
        if patch_det:
            for a in ("opm_matched", "opm_dense"):
                rows.append((f"{LABEL[a]} vs\nNeuromag combined,\n10-mm patches", on_map(np.log2(
                    patch_det[(a, "opm", cond, 10.0, "fixed_total")] / patch_det[("squid", "combined", cond, 10.0, "fixed_total")]))))
        plotting.cortex_map_figure(hemis, rows, n_lh, plt.get_cmap("RdBu_r"), plt.Normalize(-1, 1),
                                   f"G2: log2 detectability ratio OPM / Neuromag, {cond} (red: OPM higher; {grey})", "log2 ratio",
                                   OUT / f"Figure_G2_maps_{cond.replace('+', '_')}.png", contour_level=0.0)

    # 3b. absolute detectability maps (10 nAm, first headline condition)
    rows = []
    for n, cs in (("squid", "combined"), ("squid", "grad"), ("squid", "mag"), ("opm_matched", "opm"), ("opm_dense", "opm")):
        rows.append((LABEL[cs if n == "squid" else n], on_map(np.log10(res[n][(cs, headline[0])]["detect"]))))
    plotting.cortex_map_figure(hemis, rows, n_lh, plt.get_cmap("viridis"), plt.Normalize(-1.5, 0.5),
                               f"G2: log10 detectability of a 10-nAm dipole, {headline[0]} ({grey})", "log10 detectability",
                               OUT / "Figure_G2_maps_absolute.png")

    # 4. patches
    radii = patches["radii_mm"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, conv_name in zip(axs[:2], ("fixed_total", "fixed_density")):
        for n, cs in [("squid", c) for c in REFS] + [(a, "opm") for a in OPMS]:
            key = cs if n == "squid" else n
            y = [patches["median_detectability"][f"{n}/{cs}/{headline[0]}/{r:g}mm/{conv_name}"] for r in radii]
            ax.plot(radii, y, "o-", label=LABEL[key], color=colors[key])
        ax.set_yscale("log")
        ax.set_xlabel("patch geodesic radius [mm]")
        ax.set_ylabel("median detectability")
        ax.set_title(f"{conv_name.replace('_', ' ')} ({headline[0]})", fontsize=9)
        ax.legend(fontsize=7)
    for a in OPMS:
        y = [patches["comparisons"][f"{a}/combined/{headline[0]}/{r:g}mm"]["median_log2"] for r in radii]
        axs[2].plot(radii, y, "o-", label=LABEL[a], color=colors[a])
    axs[2].axhline(0, color="0.5", lw=0.8)
    axs[2].set_xlabel("patch geodesic radius [mm]")
    axs[2].set_ylabel("median log2 ratio vs Neuromag combined")
    axs[2].set_title("OPM / Neuromag (convention-independent)", fontsize=9)
    axs[2].legend(fontsize=7)
    fig.suptitle("G2: extended sources (geodesic patches around each target)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_patches.png", dpi=150)
    plt.close(fig)

    # 5. sensitivity forest plot (opm_dense and opm_matched vs Neuromag combined, headline conditions)
    labels = []
    vals = {a: {c: [] for c in headline} for a in ("opm_matched", "opm_dense")}
    cis = {a: {c: [] for c in headline} for a in ("opm_matched", "opm_dense")}
    entries = [("primary (oracle)", primary["oracle"])] + [(f"plug-in {k.split('_')[1]}", primary[k]) for k in primary if k.startswith("plugin_")]
    entries += [(k.replace("_", " "), v) for k, v in sens.items()]
    for lab, comp in entries:
        labels.append(lab)
        for a in vals:
            for c in headline:
                x = comp.get(f"{a}/combined/{c}")
                vals[a][c].append(x["median_log2"] if x else np.nan)
                cis[a][c].append(x["ci95"] if x else [np.nan, np.nan])
    fig, axs = plt.subplots(1, 2, figsize=(12, 0.28 * len(labels) + 1.5), sharey=True)
    y = np.arange(len(labels))[::-1]
    for ax, c in zip(axs, headline):
        for a, mk in (("opm_matched", "o"), ("opm_dense", "s")):
            v, ci = np.array(vals[a][c]), np.array(cis[a][c])
            ax.errorbar(v, y, xerr=[v - ci[:, 0], ci[:, 1] - v], fmt=mk, color=colors[a], label=LABEL[a], ms=4, capsize=2)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_xlabel("median log2(d_OPM / d_Neuromag combined); bars: bootstrap 95 % CI over targets", fontsize=7)
        ax.set_title(f"median log2 detectability ratio vs Neuromag combined, {c}", fontsize=8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_sensitivity.png", dpi=150)
    plt.close(fig)

    # 6. arrays
    fig = plt.figure(figsize=(14, 4.2))
    for i, (name, a) in enumerate(arrays.items()):
        ax = fig.add_subplot(1, 4, i + 1, projection="3d")
        pos = np.array([c["loc"][:3] for c in a.info["chs"]])
        if name == "squid":
            t = a.info["dev_head_t"]["trans"]
            pos = pos @ t[:3, :3].T + t[:3, 3]
            pos = pos[a.kinds == "mag"]
        mri_head = np.linalg.inv(st.subject.trans["trans"])
        sc = st.subject.scalp.rr[::40] @ mri_head[:3, :3].T + mri_head[:3, 3]
        ax.scatter(*sc.T * 1e3, s=0.3, c="0.8")
        ax.scatter(*pos.T * 1e3, s=6, c=colors.get(name, "k"))
        ax.set_title(f"{LABEL[name]}: {a.n} channels", fontsize=8)
        ax.view_init(20, -60)
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_arrays.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    if "--replot" in sys.argv:
        with open(STATE, "rb") as fh:
            saved = pickle.load(fh)
        # keep the provenance of the computation (older checkpoints: from the summary it wrote)
        if "provenance" not in saved["summary"] and (OUT / "g2_summary.json").exists():
            saved["summary"]["provenance"] = json.loads((OUT / "g2_summary.json").read_text())["provenance"]
        saved["summary"]["replotted_at_commit"] = io.RUN_COMMIT
        mne.set_log_level("WARNING")
        OUT.mkdir(parents=True, exist_ok=True)
        outputs(Study(saved["summary"]["config"]), saved)
    else:
        main()
