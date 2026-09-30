#!/usr/bin/env python3
"""G4 (adult): bounded IED localization study (NEW).

Truth: events simulated as in scripts/g4_epilepsy_adult.py (3-layer BEM, exact positions, shared
noise across arrays). Inverse model with bounded mismatch: 1-layer (inner skull) BEM and a
coregistration error of 2 mm and 2 deg (random direction and axis per array), noise covariance
estimated from 5 min of independent null data. Sources off the inverse grid: MNE/dSPM on a 5-mm
Poisson-disk grid of usable vertices that excludes the true source vertices. Equivalent current
dipole with MNE's fit_dipole (same BEM and transform). Each event is also passed through a
practical detector (1 false event per minute, thresholds from independent null data), so that
localization among detected events and joint detection + localization are reported separately.

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

    # inverse model: 5-mm off-grid source grid, 1-layer BEM, perturbed coregistration
    valid = np.flatnonzero(st.cortex.usable)
    grid = np.setdiff1d(valid[goldenholz.poisson_disk(st.cortex.rr[valid], lc["grid_spacing_mm"] * 1e-3, rng)], tgt)
    bem1 = st.subject.bem_model((0.3,))
    bem1_sol = mne.make_bem_solution(bem1, verbose=False)
    inv, trans_used = {}, {}
    for name, a in arrays.items():
        t_err = localization.perturb_trans(st.subject.trans["trans"], rng, lc["coreg_shift_mm"] * 1e-3, lc["coreg_angle_deg"])
        trans_used[name] = mne.transforms.Transform("head", "mri", t_err)
        g_inv, _ = forward.discrete_gain(a.info, trans_used[name], st.cortex.rr[grid], st.cortex.nn[grid], bem1, opm.coil_def_file())
        inv[name] = g_inv.astype(np.float64)
    log(f"{len(loc)} locations, inverse grid {len(grid)} sources (5 mm, off-grid)")

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
    minv = {name: localization.MNEInverse.make(inv[name], covs[name], snr=lc["snr"], depth=lc["depth"]) for name in arrays}
    log("inverse operators and detectors ready")

    events = [(i, fam, s) for i in range(len(loc)) for fam in ("focal", "patch") for s in lc["strengths_nAm"]]
    rows = []
    tol = int(round(cfg["detector"]["hit_tolerance_s"] * fs))
    for e_idx, (i, fam, s) in enumerate(events):
        seg = gen.segment(4.0, rng)
        t_peak = int(round(2.0 * fs))
        for name, a in arrays.items():
            y = seg[name]
            ied.inject(y, topo[(name, i, fam)], tpl, pk, s * 1e-9, t_peak)
            stat, _ = dets[name].statistic(y)
            detected = bool(stat[t_peak - tol:t_peak + tol + 1].max() > thr[name])
            d = minv[name].apply(y[:, [t_peak]], "dSPM")[:, 0]
            err, j = localization.peak_error(d, st.cortex.rr[grid], centre[i])
            sup = localization.support_recovery(d, st.cortex.rr[grid], centre[i], cfg["events"]["patch_radius_mm"] * 1e-3) if fam == "patch" else np.nan
            ev_ = evoked_for(a, y[:, t_peak:t_peak + 1], fs)
            cov = mne.Covariance(covs[name], a.info.ch_names, [], [], nfree=int(lc["baseline_min"] * 60 * fs))
            with mne.use_coil_def(opm.coil_def_file()):
                dip, _ = mne.fit_dipole(ev_, cov, bem1_sol, trans_used[name], min_dist=5.0, verbose=False)
            pos_mri = mne.transforms.apply_trans(st.subject.trans, dip.pos[0])  # true head->MRI transform
            ecd_err = float(np.linalg.norm(pos_mri - centre[i]))
            rows.append(dict(event=e_idx, location=i, family=fam, strength_nAm=s, array=name, detected=detected,
                             dspm_error_mm=err * 1e3, dspm_peak=float(np.abs(d).max()), support=sup, ecd_error_mm=ecd_err * 1e3,
                             ecd_gof=float(dip.gof[0]), depth_mm=float(st.src.depth_mm[loc[i]]), stratum=[int(x) for x in strata[i]]))
        if e_idx % 20 == 0:
            log(f"event {e_idx + 1}/{len(events)}")
    with open(OUT / "g4_localization_events.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "stratum"])
        wr.writeheader()
        for r in rows:
            wr.writerow({k: v for k, v in r.items() if k != "stratum"})
    summary = dict(status="NEW (G4 adult: bounded localization; 1-layer BEM and 2-mm/2-deg coregistration error in the inverse)",
                   config=cfg["localization"], n_events=len(events), inverse_grid=int(len(grid)), thresholds_1_per_min=thr, results={})
    for name in arrays:
        for fam in ("focal", "patch"):
            for s in lc["strengths_nAm"]:
                sel = [r for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
                det = np.array([r["detected"] for r in sel])
                e_d = np.array([r["dspm_error_mm"] for r in sel])
                e_e = np.array([r["ecd_error_mm"] for r in sel])
                gof = np.array([r["ecd_gof"] for r in sel])
                failed = gof < lc["ecd_min_gof_pct"]
                key = f"{name}/{fam}/{s:g}nAm"
                summary["results"][key] = dict(
                    n=len(sel), detected=float(det.mean()),
                    dspm_error_mm_median_all=float(np.median(e_d)), dspm_error_mm_median_detected=float(np.median(e_d[det])) if det.any() else None,
                    ecd_error_mm_median_all=float(np.median(e_e)),
                    ecd_error_mm_median_detected_not_failed=float(np.median(e_e[det & ~failed])) if np.any(det & ~failed) else None,
                    ecd_failed_fits=float(failed.mean()), ecd_gof_median=float(np.median(gof)),
                    joint_detect_and_dspm_within_10mm=float(np.mean(det & (e_d <= 10))),
                    joint_detect_and_ecd_within_10mm=float(np.mean(det & ~failed & (e_e <= 10))),
                    support_recovery_median=float(np.nanmedian([r["support"] for r in sel])) if fam == "patch" else None)
    io.write_json(summary, OUT / "g4_localization_summary.json")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    colors = {"squid": "k", "opm_matched": "tab:green", "opm_dense": "tab:purple"}
    for ax, meth in zip(axs, ("dspm", "ecd")):
        for k_, name in enumerate(arrays):
            for m_, (fam, s) in enumerate([(f, s) for f in ("focal", "patch") for s in lc["strengths_nAm"]]):
                v = [r[f"{meth}_error_mm"] for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
                ax.boxplot(v, positions=[m_ * 4 + k_], widths=0.7, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=colors[name], alpha=0.4))
        ax.set_xticks([m_ * 4 + 1 for m_ in range(4)])
        ax.set_xticklabels([f"{f}\n{s:g} nAm" for f in ("focal", "patch") for s in lc["strengths_nAm"]], fontsize=8)
        ax.set_ylabel("localization error [mm] (all events)")
        ax.set_title({"dspm": "dSPM peak (5-mm off-grid)", "ecd": "equivalent current dipole"}[meth], fontsize=9)
    for name, c in colors.items():
        axs[0].plot([], [], "s", color=c, alpha=0.5, label=name)
    axs[0].legend(fontsize=8)
    fig.suptitle("G4 adult: localization with 1-layer BEM and 2-mm/2-deg coregistration error (Neuromag combined vs OPM)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G4_localization.png", dpi=150)
    plt.close(fig)
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
