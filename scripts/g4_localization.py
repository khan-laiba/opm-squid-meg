#!/usr/bin/env python3
"""G4: bounded IED localization study (NEW), adult (this script) and, through
``g4_epilepsy_adult.Context``, any anatomy with its arrays (scripts/g4_epilepsy_pediatric.py).

Truth: events simulated as in scripts/g4_epilepsy_adult.py (3-layer BEM, exact positions, shared
noise across arrays). Inverse model with bounded mismatch: 1-layer (inner skull) BEM and a
coregistration error of 2 mm and 2 deg. The error is drawn K times (random translation direction
and rotation axis); location i uses draw i mod K for all its events, in every array, so the arrays
are compared under identical errors and every condition sees all K draws (3 locations each for
24 locations and K = 8). Noise
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
localization are reported separately; the thresholds are checked on further held-out null data.
Secondary (review, 2026-10-02): Neuromag with magnetometers or gradiometers alone (their own
covariance, detector and inverse), and dSPM through MNE's own ``mne.minimum_norm`` next to the
study's implementation (opmsquid.localization), which differs from it in the depth weighting. Paired OPM-minus-Neuromag differences on identical events:
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

import g4_epilepsy_adult as G4  # noqa: E402
from opmsquid import detection, environment, forward, g2, goldenholz, ied, io, localization, metrics, opm  # noqa: E402

OUT = ROOT / "results" / "g4"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def evoked_for(info_or_array, data, fs):
    info = getattr(info_or_array, "info", info_or_array).copy()
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
    lc = cfg["localization"]
    ctx = G4.adult_context(cfg2, tuple(lc["arrays"]))
    localize(ctx, cfg, np.random.default_rng(lc["seed"]))
    log(f"done in {time.time() - t0:.0f} s")


def localize(ctx, cfg, rng):
    """The bounded localization study for one anatomy (see the module docstring); writes
    results/g4/g4_localization[_<label>]_{events.csv,summary.json} and the figure."""
    sim, lc = cfg["simulation"], cfg["localization"]
    arrays, G_t, G_g, env = ctx.arrays, ctx.G_t, ctx.G_g, ctx.env
    specs = {name: ied.ArraySpec(G_g[name], environment.external_basis(a.info, env.r0, a.coil_def),
                                 np.array([g2.SQUID_ASD[k] for k in a.kinds]) if name == "squid"
                                 else np.full(a.n, sim["opm_asd_fT_per_rtHz"] * 1e-15)) for name, a in arrays.items()}
    gen = ied.NoiseGenerator(specs, ctx.brain_scale * ctx.src.grid_area, env.coef_timecourses, env.coef_cov, ctx.filt,
                             decimate=sim["decimate"], env_fs=env.sfreq)
    fs = gen.fs_out
    tpl, pk = ied.filtered_template(ctx.filt, 1.0, sim["decimate"])

    # events: 2 locations per depth x orientation stratum, focal and 10-mm patch, two strengths
    loc, strata = G4.stratified_locations(ctx, lc["n_locations"], rng)
    tgt = ctx.src.target[loc]
    members = goldenholz.geodesic_patches(ctx.cortex.adjacency, tgt, cfg["events"]["patch_radius_mm"] * 1e-3, ctx.cortex.usable)
    topo, centre = {}, {}
    for name in arrays:
        unit = ctx.patch_topographies(name, members)
        for i in range(len(loc)):
            topo[(name, i, "focal")] = G_t[name][:, loc[i]]
            topo[(name, i, "patch")] = unit[:, i] / ctx.cortex.area[members[i]].sum()
    for i in range(len(loc)):
        centre[i] = ctx.cortex.rr[tgt[i]]

    # inverse model: 5-mm off-grid source grid, 1-layer BEM, K coregistration-error draws shared by all arrays
    valid = np.flatnonzero(ctx.cortex.usable)
    grid = np.setdiff1d(valid[goldenholz.poisson_disk(ctx.cortex.rr[valid], lc["grid_spacing_mm"] * 1e-3, rng)], tgt)
    bem1 = ctx.subject.bem_model((0.3,))
    bem1_sol = mne.make_bem_solution(bem1, verbose=False)
    n_draws = lc["coreg_draws"]
    trans_used = [mne.transforms.Transform("head", "mri", localization.perturb_trans(
        ctx.subject.trans["trans"], rng, lc["coreg_shift_mm"] * 1e-3, lc["coreg_angle_deg"])) for _ in range(n_draws)]
    mri_to_head = mne.transforms.invert_transform(ctx.subject.trans)
    views = localization.channel_views(arrays, lc.get("neuromag_subsets", []))  # name -> (physical array, channel indices)
    inv_phys, fwd_phys = {}, {}
    for name, a in arrays.items():
        for k in range(n_draws):
            g_inv, _ = forward.discrete_gain(a.info, trans_used[k], ctx.cortex.rr[grid], ctx.cortex.nn[grid], bem1, opm.coil_def_file())
            inv_phys[name, k] = g_inv.astype(np.float64)
            if lc.get("mne_dspm"):
                fwd_phys[name, k] = forward.discrete_forward(a.info, trans_used[k], ctx.cortex.rr[grid], ctx.cortex.nn[grid], bem1_sol,
                                                             opm.coil_def_file())
    inv = {(v, k): inv_phys[phys, k][idx] for v, (phys, idx) in views.items() for k in range(n_draws)}
    vinfo = {v: mne.pick_info(arrays[phys].info, idx) for v, (phys, idx) in views.items()}
    log(f"{ctx.label}: {len(loc)} locations, inverse grid {len(grid)} sources (5 mm, off-grid), {n_draws} coregistration draws; "
        f"views {list(views)}")

    # noise covariance (5 min) and a practical detector (thresholds from 10 min) per array, independent null data
    def null_data(minutes, stream=None):
        acc = {name: [] for name in arrays}
        for _ in range(int(round(minutes * 60 / sim["segment_s"]))):
            for name, y in gen.segment(sim["segment_s"], rng if stream is None else stream).items():
                acc[name].append(y)
        return {name: np.concatenate(v, axis=1) for name, v in acc.items()}

    base = null_data(lc["baseline_min"])
    covs = {v: metrics.ledoit_wolf_covariance(base[phys][idx])[0] for v, (phys, idx) in views.items()}
    whit = {v: metrics.whitener(c) for v, c in covs.items()}
    del base
    cal = null_data(lc["calibration_min"])
    valid_c = np.flatnonzero(ctx.cortex.usable)
    cand = np.setdiff1d(valid_c[goldenholz.poisson_disk(ctx.cortex.rr[valid_c], cfg["detector"]["dictionary_spacing_mm"] * 1e-3, rng)], tgt)
    templates = {s: ied.filtered_template(ctx.filt, s, sim["decimate"]) for s in cfg["events"]["stretches"]}
    dets, thr = {}, {}
    refr = int(round(cfg["detector"]["refractory_s"] * fs))
    cand_gain = {name: ctx.gain(name, cand) for name in arrays}
    for v, (phys, idx) in views.items():
        dets[v] = detection.ScanDetector.build(cal[phys][idx], whit[v], templates, cand_gain[phys][idx], refr)
        _, h = dets[v].events(dets[v].statistic(cal[phys][idx])[0])
        thr[v] = detection.threshold_for_rate(h, lc["calibration_min"], 1.0)
    del cal, cand_gain
    heldout = {}
    if lc.get("holdout_min"):  # the thresholds on independent null data (own random stream: the events stay those of the primary run)
        held = null_data(lc["holdout_min"], np.random.default_rng(lc["seed"] + 1))
        for v, (phys, idx) in views.items():
            _, h = dets[v].events(dets[v].statistic(held[phys][idx])[0])
            heldout[v] = detection.rate_with_ci(int(np.sum(h > thr[v])), lc["holdout_min"])
        del held
    minv = {key: localization.MNEInverse.make(g_, covs[key[0]], snr=lc["snr"], depth=lc["depth"]) for key, g_ in inv.items()}
    nfree = int(lc["baseline_min"] * 60 * fs)
    mne_cov = {v: mne.Covariance(covs[v], vinfo[v].ch_names, [], [], nfree=nfree) for v in views}
    mne_inv = {}
    if lc.get("mne_dspm"):
        lambda2 = 1.0 / lc["snr"] ** 2
        for v, (phys, idx) in views.items():
            for k in range(n_draws):
                fwd_v = mne.pick_channels_forward(fwd_phys[phys, k], include=vinfo[v].ch_names, ordered=True, verbose=False)
                op = mne.minimum_norm.make_inverse_operator(vinfo[v], fwd_v, mne_cov[v], loose=0.0, fixed=True, depth=lc["depth"],
                                                            verbose=False)
                mne_inv[v, k] = mne.minimum_norm.prepare_inverse_operator(op, nave=1, lambda2=lambda2, method="dSPM", verbose=False)
        del fwd_phys
    log("inverse operators and detectors ready")

    events = [(i, fam, s) for i in range(len(loc)) for fam in ("focal", "patch") for s in lc["strengths_nAm"]]
    rows = []
    tol = int(round(cfg["detector"]["hit_tolerance_s"] * fs))
    for e_idx, (i, fam, s) in enumerate(events):
        k = i % n_draws  # coregistration draw of this location: shared by all arrays and conditions; every draw in every condition
        seg = gen.segment(4.0, rng)
        t_peak = int(round(2.0 * fs))
        # displacement the coregistration error alone produces at the true source (MRI frame)
        coreg_mm = 1e3 * float(np.linalg.norm(mne.transforms.apply_trans(
            trans_used[k], mne.transforms.apply_trans(mri_to_head, centre[i])) - centre[i]))
        for name in arrays:
            ied.inject(seg[name], topo[(name, i, fam)], tpl, pk, s * 1e-9, t_peak)
        for name, (phys, idx) in views.items():
            y = seg[phys][idx]
            stat, _ = dets[name].statistic(y)
            detected = detection.event_height(*dets[name].events(stat), t_peak, tol) > thr[name]
            d = minv[name, k].apply(y[:, [t_peak]], "dSPM")[:, 0]
            err, j = localization.peak_error(d, ctx.cortex.rr[grid], centre[i])
            sup = localization.support_recovery(d, np.isin(grid, members[i])) if fam == "patch" else np.nan
            ev_ = evoked_for(vinfo[name], y[:, t_peak:t_peak + 1], fs)
            if mne_inv:
                d_mne = mne.minimum_norm.apply_inverse(ev_, mne_inv[name, k], 1.0 / lc["snr"] ** 2, "dSPM", prepared=True,
                                                       verbose=False).data[:, 0]
                err_mne = localization.peak_error(d_mne, ctx.cortex.rr[grid], centre[i])[0]
            else:
                d_mne, err_mne = np.full(1, np.nan), np.nan
            with mne.use_coil_def(opm.coil_def_file()):
                dip, _ = mne.fit_dipole(ev_, mne_cov[name], bem1_sol, trans_used[k], min_dist=lc["ecd_min_dist_mm"], verbose=False)
            # where the analyst reads the dipole (MRI via the perturbed transform), and in the sensor frame (true transform)
            ecd_mri = mne.transforms.apply_trans(trans_used[k], dip.pos[0])
            ecd_sensor = mne.transforms.apply_trans(ctx.subject.trans, dip.pos[0])
            rows.append(dict(event=e_idx, location=i, family=fam, strength_nAm=s, array=name, coreg_draw=k, detected=detected,
                             dspm_error_mm=err * 1e3, dspm_peak=float(np.abs(d).max()), support=sup,
                             dspm_mne_error_mm=err_mne * 1e3, dspm_mne_peak=float(np.abs(d_mne).max()),
                             ecd_error_mm=1e3 * float(np.linalg.norm(ecd_mri - centre[i])),
                             ecd_error_sensor_frame_mm=1e3 * float(np.linalg.norm(ecd_sensor - centre[i])),
                             coreg_displacement_mm=coreg_mm, ecd_gof=float(dip.gof[0]), ecd_conf_vol_mm3=float(dip.conf["vol"][0]) * 1e9,
                             ecd_khi2_per_dof=float(dip.khi2[0] / dip.nfree[0]),
                             depth_mm=float(ctx.src.depth_mm[loc[i]]), stratum=[int(x) for x in strata[i]]))
        if e_idx % 20 == 0:
            log(f"event {e_idx + 1}/{len(events)}")
    tag = "" if ctx.label == "adult" else f"_{ctx.label}"
    status = (f"NEW (G4 {ctx.label}: bounded localization; 1-layer BEM and 2-mm/2-deg coregistration error in the "
              f"inverse, {n_draws} draws shared by all arrays)")
    with open(OUT / f"g4_localization{tag}_events.csv", "w", newline="") as fh:
        io.csv_status(fh, status)
        wr = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "stratum"])
        wr.writeheader()
        for r in rows:
            wr.writerow({k: v for k, v in r.items() if k != "stratum"})
    summary = dict(status=status, anatomy=ctx.notes,
                   config=cfg["localization"], n_events=len(events), inverse_grid=int(len(grid)), thresholds_1_per_min=thr,
                   views={v: dict(array=phys, channels=int(len(idx))) for v, (phys, idx) in views.items()},
                   thresholds_heldout=dict(minutes=lc.get("holdout_min"), false_events=heldout,
                                           note="false events per minute at the calibrated thresholds on independent null data "
                                                "(exact Poisson 95 % interval); target 1 per minute"),
                   coreg_displacement_mm_median=float(np.median([r["coreg_displacement_mm"] for r in rows])),
                   results=summarise(rows, views, lc), paired=paired(rows, views, lc, rng),
                   paired_secondary=paired(rows, views, lc, np.random.default_rng(lc["seed"] + 2),
                                           comparators=[v for v in views if v.startswith("squid_")]))
    io.write_json(summary, OUT / f"g4_localization{tag}_summary.json")
    figure(rows, views, lc, ctx.label)
    return rows, summary


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
                                                                        "ecd_gof", "ecd_conf_vol_mm3", "ecd_khi2_per_dof", "support",
                                                                        "dspm_mne_error_mm")}
                out[f"{name}/{fam}/{s:g}nAm"] = dict(
                    n=len(sel), detected=float(det.mean()),
                    dspm_error_mm_median_all=_med(col["dspm_error_mm"]), dspm_error_mm_median_detected=_med(col["dspm_error_mm"][det]),
                    ecd_error_mm_median_all=_med(col["ecd_error_mm"]), ecd_error_mm_median_detected=_med(col["ecd_error_mm"][det]),
                    ecd_error_sensor_frame_mm_median_detected=_med(col["ecd_error_sensor_frame_mm"][det]),
                    ecd_gof_median_detected=_med(col["ecd_gof"][det]), ecd_conf_vol_mm3_median_detected=_med(col["ecd_conf_vol_mm3"][det]),
                    ecd_khi2_per_dof_median_detected=_med(col["ecd_khi2_per_dof"][det]),
                    dspm_mne_error_mm_median_all=_med(col["dspm_mne_error_mm"]),
                    dspm_mne_error_mm_median_detected=_med(col["dspm_mne_error_mm"][det]),
                    joint_detect_and_dspm_within_10mm=float(np.mean(det & (col["dspm_error_mm"] <= 10))),
                    joint_detect_and_dspm_mne_within_10mm=float(np.mean(det & (col["dspm_mne_error_mm"] <= 10))),
                    joint_detect_and_ecd_within_10mm=float(np.mean(det & (col["ecd_error_mm"] <= 10))),
                    support_recovery_median=_med(col["support"]) if fam == "patch" else None)
    return out


def paired(rows, arrays, lc, rng, n_boot=2000, comparators=("squid",)):
    """OPM minus Neuromag on identical events (same location, noise and coregistration draw);
    ``comparators``: the Neuromag views compared with (primary: combined)."""
    from scipy.stats import binomtest, wilcoxon

    by = {(r["array"], r["event"]): r for r in rows}
    out = {}
    pairs = [(a, ref) for ref in comparators for a in arrays if not a.startswith("squid")]
    for a, ref in pairs:
        for fam in ("focal", "patch"):
            for s in lc["strengths_nAm"]:
                ev = sorted({r["event"] for r in rows if r["family"] == fam and r["strength_nAm"] == s})
                res = dict(n=len(ev))
                for metric in ("dspm_error_mm", "ecd_error_mm", "dspm_mne_error_mm"):
                    if not np.all(np.isfinite([by[a, e][metric] for e in ev])):
                        continue
                    diff = np.array([by[a, e][metric] - by[ref, e][metric] for e in ev])
                    boot = np.median(diff[rng.integers(0, len(diff), (n_boot, len(diff)))], axis=1)
                    p = float(wilcoxon(diff).pvalue) if np.any(diff != 0) else 1.0
                    res[metric] = dict(median_difference=float(np.median(diff)), ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                                       wilcoxon_p=p, share_opm_smaller=float(np.mean(diff < 0)))
                for metric in ("dspm", "ecd"):
                    ok = {n: np.array([by[n, e]["detected"] and by[n, e][f"{metric}_error_mm"] <= 10 for e in ev]) for n in (a, ref)}
                    only_a, only_b = int(np.sum(ok[a] & ~ok[ref])), int(np.sum(~ok[a] & ok[ref]))
                    res[f"joint_{metric}_10mm"] = dict(opm=float(ok[a].mean()), squid=float(ok[ref].mean()), only_opm=only_a, only_squid=only_b,
                                                       mcnemar_exact_p=float(binomtest(only_a, only_a + only_b, 0.5).pvalue) if only_a + only_b else 1.0)
                out[f"{a}_vs_{ref}/{fam}/{s:g}nAm"] = res
    return out


def figure(rows, arrays, lc, label="adult"):
    colors = {"squid": "k", "squid_mag": "tab:blue", "squid_grad": "tab:red", "opm_matched": "tab:green", "opm204": "tab:cyan",
              "opm_dense": "tab:purple"}
    conds = [(f, s) for f in ("focal", "patch") for s in lc["strengths_nAm"]]
    step = len(arrays) + 1
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, meth in zip(axs[:2], ("dspm", "ecd")):
        for k_, name in enumerate(arrays):
            for m_, (fam, s) in enumerate(conds):
                v = [r[f"{meth}_error_mm"] for r in rows if r["array"] == name and r["family"] == fam and r["strength_nAm"] == s]
                ax.boxplot(v, positions=[m_ * step + k_], widths=0.7, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=colors.get(name, "0.5"), alpha=0.4))
        ax.set_xticks([m_ * step + (len(arrays) - 1) / 2 for m_ in range(len(conds))])
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
    fig.suptitle(f"G4 {label}: localization with a 1-layer BEM and 2-mm/2-deg coregistration error (same draws for all arrays)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"Figure_G4_localization{'' if label == 'adult' else '_' + label}.png", dpi=150)
    plt.close(fig)

if __name__ == "__main__":
    main()
