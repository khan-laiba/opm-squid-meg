#!/usr/bin/env python3
"""G4 (adult): bounded IED localization study (NEW).

Truth: events simulated as in scripts/g4_epilepsy_adult.py (3-layer BEM, exact positions, shared
noise across arrays). Inverse model with bounded mismatch: 1-layer (inner skull) BEM and a
coregistration error of 2 mm and 2 deg. The error is drawn K times (random translation direction
and rotation axis); every array uses the same K draws and event e uses draw e mod K, so the arrays
are compared under identical errors averaged over K draws rather than one draw per array. Noise
covariance from 5 min of independent null data. Sources off the inverse grid: MNE/dSPM on a 5-mm
Poisson-disk grid of usable vertices that excludes the true source vertices. Equivalent current
dipole with MNE's fit_dipole (same BEM and transform).

Errors are measured where the analyst reads them: the dSPM peak and the dipole are placed on the
MRI with the analyst's (perturbed) head->MRI transform and compared with the true source in MRI
coordinates. The dipole error in the sensor frame (true transform) is kept as a decomposition.
No goodness-of-fit cut is applied: MNE computes GOF on whitened data, where the noise adds about
one unit per channel, so a fixed GOF threshold favours arrays with fewer channels. GOF and the
95 % confidence volume (MNE's linearised estimate, in physical units) are reported descriptively.

Each event is also passed through a practical detector (1 false event per minute, thresholds from
independent null data), so that localization among detected events and joint detection +
localization are reported separately. Paired OPM-minus-Neuromag differences on identical events:
median error difference with a bootstrap CI over events (one event per location and condition),
Wilcoxon signed-rank p, and McNemar exact p for joint detection + localization within 10 mm.

Configuration: configs/g4_epilepsy.toml ([localization]). Outputs: results/g4/g4_localization_*.
"""
from __future__ import annotations

import csv
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402

import g2_adult_comparison as G  # noqa: E402
import g4_epilepsy_adult as G4  # noqa: E402
from opmsquid import background, detection, environment, forward, g2, goldenholz, ied, io, localization, metrics, opm  # noqa: E402

OUT = ROOT / "results" / "g4"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def evoked_for(array, data, fs):
    info = array.info.copy()
    with info._unlock():
        info["sfreq"] = fs
        info["lowpass"], info["highpass"] = 40.0, 1.0
    return mne.EvokedArray(data, info, tmin=0.0, nave=1, verbose=False)


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads((ROOT / "configs" / "g4_epilepsy.toml").read_text())
    cfg2 = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    sim, lc = cfg["simulation"], cfg["localization"]
    st = G.Study(cfg2)
    rng = np.random.default_rng(lc["seed"])
    arrays = {k: v for k, v in g2.build_arrays(st.subject, st.dig).items() if k in lc["arrays"]}
    squid = arrays["squid"]
    bads = cfg2["sensors"]["bads"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads = good & (squid.kinds == "grad")
    meas = g2.measured_noise(squid.info, st.filt, bads)
    env = meas["environment"]
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        G_t[name], G_g[name] = st.gains(a)
    brain_scale = background.calibrate(st.unit_brain(G_g["squid"]), grads, float(np.nanmedian(meas["brain"][grads])))
    specs = {name: ied.ArraySpec(G_g[name], environment.external_basis(a.info, env.r0, a.coil_def),
                                 np.array([g2.SQUID_ASD[k] for k in a.kinds]) if name == "squid"
                                 else np.full(a.n, sim["opm_asd_fT_per_rtHz"] * 1e-15)) for name, a in arrays.items()}
    gen = ied.NoiseGenerator(specs, brain_scale * st.src.grid_area, env.coef_timecourses, env.coef_cov, st.filt,
                             decimate=sim["decimate"], env_fs=env.sfreq)
    fs = gen.fs_out
    tpl, pk = ied.filtered_template(st.filt, 1.0, sim["decimate"])

    # events: 2 locations per depth x orientation stratum, focal and 10-mm patch, two strengths
    loc, strata = G4.stratified_locations(st, lc["n_locations"], rng)
    tgt = st.src.target[loc]
    members = goldenholz.geodesic_patches(st.cortex.adjacency, tgt, cfg["events"]["patch_radius_mm"] * 1e-3, st.cortex.usable)
    topo, centre = {}, {}
    for name, a in arrays.items():
        full, col = g2.fullres_matrix(a, st.subject, st.cortex, g2.FULLRES_JOBS[name])
        unit = goldenholz.patch_topographies(full, members, col, st.cortex.area)
        for i in range(len(loc)):
            topo[(name, i, "focal")] = G_t[name][:, loc[i]]
            topo[(name, i, "patch")] = unit[:, i] / st.cortex.area[members[i]].sum()
        del full
    for i in range(len(loc)):
        centre[i] = st.cortex.rr[tgt[i]]

    # inverse model: 5-mm off-grid source grid, 1-layer BEM, K coregistration-error draws shared by all arrays
    valid = np.flatnonzero(st.cortex.usable)
    grid = np.setdiff1d(valid[goldenholz.poisson_disk(st.cortex.rr[valid], lc["grid_spacing_mm"] * 1e-3, rng)], tgt)
    bem1 = st.subject.bem_model((0.3,))
    bem1_sol = mne.make_bem_solution(bem1, verbose=False)
    n_draws = lc["coreg_draws"]
    trans_used = [mne.transforms.Transform("head", "mri", localization.perturb_trans(
        st.subject.trans["trans"], rng, lc["coreg_shift_mm"] * 1e-3, lc["coreg_angle_deg"])) for _ in range(n_draws)]
    mri_to_head = mne.transforms.invert_transform(st.subject.trans)
    inv = {}
    for name, a in arrays.items():
        for k in range(n_draws):
            g_inv, _ = forward.discrete_gain(a.info, trans_used[k], st.cortex.rr[grid], st.cortex.nn[grid], bem1, opm.coil_def_file())
            inv[name, k] = g_inv.astype(np.float64)
    log(f"{len(loc)} locations, inverse grid {len(grid)} sources (5 mm, off-grid), {n_draws} coregistration draws")

    # noise covariance (5 min) and a practical detector (thresholds from 10 min) per array, independent null data
    def null_data(minutes):
        acc = {name: [] for name in arrays}
        for _ in range(int(round(minutes * 60 / sim["segment_s"]))):
            for name, y in gen.segment(sim["segment_s"], rng).items():
                acc[name].append(y)
        return {name: np.concatenate(v, axis=1) for name, v in acc.items()}

    base = null_data(lc["baseline_min"])
    covs = {name: metrics.ledoit_wolf_covariance(base[name])[0] for name in arrays}
    whit = {name: metrics.whitener(covs[name]) for name in arrays}
    del base
    cal = null_data(lc["calibration_min"])
    valid_c = np.flatnonzero(st.cortex.usable)
    cand = np.setdiff1d(valid_c[goldenholz.poisson_disk(st.cortex.rr[valid_c], cfg["detector"]["dictionary_spacing_mm"] * 1e-3, rng)], tgt)
    templates = {s: ied.filtered_template(st.filt, s, sim["decimate"]) for s in cfg["events"]["stretches"]}
    dets, thr = {}, {}
    refr = int(round(cfg["detector"]["refractory_s"] * fs))
    for name, a in arrays.items():
        cg = g2.gains(a, st.subject, st.cortex, cand, g2.FULLRES_JOBS[name])
        dets[name] = detection.ScanDetector.build(cal[name], whit[name], templates, cg, refr)
        _, h = dets[name].events(dets[name].statistic(cal[name])[0])
        thr[name] = detection.threshold_for_rate(h, lc["calibration_min"], 1.0)
    del cal
    minv = {key: localization.MNEInverse.make(g_, covs[key[0]], snr=lc["snr"], depth=lc["depth"]) for key, g_ in inv.items()}
    mne_cov = {name: mne.Covariance(covs[name], a.info.ch_names, [], [], nfree=int(lc["baseline_min"] * 60 * fs)) for name, a in arrays.items()}
    log("inverse operators and detectors ready")

    events = [(i, fam, s) for i in range(len(loc)) for fam in ("focal", "patch") for s in lc["strengths_nAm"]]
    rows = []
    tol = int(round(cfg["detector"]["hit_tolerance_s"] * fs))
    for e_idx, (i, fam, s) in enumerate(events):
        k = e_idx % n_draws  # coregistration draw, shared by all arrays for this event
        seg = gen.segment(4.0, rng)
        t_peak = int(round(2.0 * fs))
        # displacement the coregistration error alone produces at the true source (MRI frame)
        coreg_mm = 1e3 * float(np.linalg.norm(mne.transforms.apply_trans(
            trans_used[k], mne.transforms.apply_trans(mri_to_head, centre[i])) - centre[i]))
        for name, a in arrays.items():
            y = seg[name]
            ied.inject(y, topo[(name, i, fam)], tpl, pk, s * 1e-9, t_peak)
            stat, _ = dets[name].statistic(y)
            detected = bool(stat[t_peak - tol:t_peak + tol + 1].max() > thr[name])
            d = minv[name, k].apply(y[:, [t_peak]], "dSPM")[:, 0]
            err, j = localization.peak_error(d, st.cortex.rr[grid], centre[i])
            sup = localization.support_recovery(d, st.cortex.rr[grid], centre[i], cfg["events"]["patch_radius_mm"] * 1e-3) if fam == "patch" else np.nan
            ev_ = evoked_for(a, y[:, t_peak:t_peak + 1], fs)
            with mne.use_coil_def(opm.coil_def_file()):
                dip, _ = mne.fit_dipole(ev_, mne_cov[name], bem1_sol, trans_used[k], min_dist=lc["ecd_min_dist_mm"], verbose=False)
            # where the analyst reads the dipole (MRI via the perturbed transform), and in the sensor frame (true transform)
            ecd_mri = mne.transforms.apply_trans(trans_used[k], dip.pos[0])
            ecd_sensor = mne.transforms.apply_trans(st.subject.trans, dip.pos[0])
            rows.append(dict(event=e_idx, location=i, family=fam, strength_nAm=s, array=name, coreg_draw=k, detected=detected,
                             dspm_error_mm=err * 1e3, dspm_peak=float(np.abs(d).max()), support=sup,
                             ecd_error_mm=1e3 * float(np.linalg.norm(ecd_mri - centre[i])),
                             ecd_error_sensor_frame_mm=1e3 * float(np.linalg.norm(ecd_sensor - centre[i])),
                             coreg_displacement_mm=coreg_mm, ecd_gof=float(dip.gof[0]), ecd_conf_vol_mm3=float(dip.conf["vol"][0]) * 1e9,
                             ecd_khi2_per_dof=float(dip.khi2[0] / dip.nfree[0]),
                             depth_mm=float(st.src.depth_mm[loc[i]]), stratum=[int(x) for x in strata[i]]))
        if e_idx % 20 == 0:
            log(f"event {e_idx + 1}/{len(events)}")
    with open(OUT / "g4_localization_events.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "stratum"])
        wr.writeheader()
        for r in rows:
            wr.writerow({k: v for k, v in r.items() if k != "stratum"})
    summary = dict(status="NEW (G4 adult: bounded localization; 1-layer BEM and 2-mm/2-deg coregistration error in the inverse, "
                          f"{n_draws} draws shared by all arrays)",
                   config=cfg["localization"], n_events=len(events), inverse_grid=int(len(grid)), thresholds_1_per_min=thr,
                   coreg_displacement_mm_median=float(np.median([r["coreg_displacement_mm"] for r in rows])),
                   results=summarise(rows, arrays, lc), paired=paired(rows, arrays, lc, rng))
    io.write_json(summary, OUT / "g4_localization_summary.json")
    figure(rows, arrays, lc)
    log(f"done in {time.time() - t0:.0f} s")


def _med(v):
    v = np.asarray(v, float)
    return float(np.median(v)) if v.size else None


def summarise(rows, arrays, lc):
    out = {}
    for name in arrays:
        for fam in ("focal", "patch"):
            for s in lc["strengths_nAm"]:
                sel = [r for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
                det = np.array([r["detected"] for r in sel])
                col = {c: np.array([r[c] for r in sel], float) for c in ("dspm_error_mm", "ecd_error_mm", "ecd_error_sensor_frame_mm",
                                                                        "ecd_gof", "ecd_conf_vol_mm3", "ecd_khi2_per_dof", "support")}
                out[f"{name}/{fam}/{s:g}nAm"] = dict(
                    n=len(sel), detected=float(det.mean()),
                    dspm_error_mm_median_all=_med(col["dspm_error_mm"]), dspm_error_mm_median_detected=_med(col["dspm_error_mm"][det]),
                    ecd_error_mm_median_all=_med(col["ecd_error_mm"]), ecd_error_mm_median_detected=_med(col["ecd_error_mm"][det]),
                    ecd_error_sensor_frame_mm_median_detected=_med(col["ecd_error_sensor_frame_mm"][det]),
                    ecd_gof_median_detected=_med(col["ecd_gof"][det]), ecd_conf_vol_mm3_median_detected=_med(col["ecd_conf_vol_mm3"][det]),
                    ecd_khi2_per_dof_median_detected=_med(col["ecd_khi2_per_dof"][det]),
                    joint_detect_and_dspm_within_10mm=float(np.mean(det & (col["dspm_error_mm"] <= 10))),
                    joint_detect_and_ecd_within_10mm=float(np.mean(det & (col["ecd_error_mm"] <= 10))),
                    support_recovery_median=_med(col["support"]) if fam == "patch" else None)
    return out


def paired(rows, arrays, lc, rng, n_boot=2000):
    """OPM minus Neuromag on identical events (same location, noise and coregistration draw)."""
    from scipy.stats import binomtest, wilcoxon

    by = {(r["array"], r["event"]): r for r in rows}
    out = {}
    for a in [n for n in arrays if n != "squid"]:
        for fam in ("focal", "patch"):
            for s in lc["strengths_nAm"]:
                ev = sorted({r["event"] for r in rows if r["family"] == fam and r["strength_nAm"] == s})
                res = dict(n=len(ev))
                for metric in ("dspm_error_mm", "ecd_error_mm"):
                    diff = np.array([by[a, e][metric] - by["squid", e][metric] for e in ev])
                    boot = np.median(diff[rng.integers(0, len(diff), (n_boot, len(diff)))], axis=1)
                    p = float(wilcoxon(diff).pvalue) if np.any(diff != 0) else 1.0
                    res[metric] = dict(median_difference=float(np.median(diff)), ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                                       wilcoxon_p=p, share_opm_smaller=float(np.mean(diff < 0)))
                for metric in ("dspm", "ecd"):
                    ok = {n: np.array([by[n, e]["detected"] and by[n, e][f"{metric}_error_mm"] <= 10 for e in ev]) for n in (a, "squid")}
                    only_a, only_b = int(np.sum(ok[a] & ~ok["squid"])), int(np.sum(~ok[a] & ok["squid"]))
                    res[f"joint_{metric}_10mm"] = dict(opm=float(ok[a].mean()), squid=float(ok["squid"].mean()), only_opm=only_a, only_squid=only_b,
                                                       mcnemar_exact_p=float(binomtest(only_a, only_a + only_b, 0.5).pvalue) if only_a + only_b else 1.0)
                out[f"{a}_vs_squid/{fam}/{s:g}nAm"] = res
    return out


def figure(rows, arrays, lc):
    colors = {"squid": "k", "opm_matched": "tab:green", "opm204": "tab:blue", "opm_dense": "tab:purple"}
    conds = [(f, s) for f in ("focal", "patch") for s in lc["strengths_nAm"]]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, meth in zip(axs[:2], ("dspm", "ecd")):
        for k_, name in enumerate(arrays):
            for m_, (fam, s) in enumerate(conds):
                v = [r[f"{meth}_error_mm"] for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
                ax.boxplot(v, positions=[m_ * 4 + k_], widths=0.7, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=colors.get(name, "0.5"), alpha=0.4))
        ax.set_xticks([m_ * 4 + 1 for m_ in range(len(conds))])
        ax.set_xticklabels([f"{f}\n{s:g} nAm" for f, s in conds], fontsize=8)
        ax.set_ylabel("localization error on the MRI [mm] (all events)")
        ax.set_title({"dspm": "dSPM peak (5-mm grid; true sources off the grid)", "ecd": "equivalent current dipole"}[meth], fontsize=9)
    ax = axs[2]
    w = 0.8 / len(arrays)
    for k_, name in enumerate(arrays):
        for m_, (fam, s) in enumerate(conds):
            sel = [r for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
            for off, meth, hatch in ((0.0, "dspm", None), (0.5, "ecd", "//")):
                v = np.mean([r["detected"] and r[f"{meth}_error_mm"] <= 10 for r in sel])
                ax.bar(m_ * 2 + off + (k_ - (len(arrays) - 1) / 2) * w / 2, v, width=w / 2, color=colors.get(name, "0.5"),
                       alpha=0.6, hatch=hatch, edgecolor="k", linewidth=0.3)
    ax.set_xticks([m_ * 2 + 0.25 for m_ in range(len(conds))])
    ax.set_xticklabels([f"{f}\n{s:g} nAm" for f, s in conds], fontsize=8)
    ax.set_ylabel("detected and within 10 mm (share of events)")
    ax.set_title("joint detection + localization (plain: dSPM, hatched: ECD)", fontsize=9)
    for name in arrays:
        axs[0].plot([], [], "s", color=colors.get(name, "0.5"), alpha=0.5, label=name)
    axs[0].legend(fontsize=8)
    fig.suptitle("G4 adult: localization with a 1-layer BEM and 2-mm/2-deg coregistration error (same draws for all arrays)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G4_localization.png", dpi=150)
    plt.close(fig)

if __name__ == "__main__":
    main()
