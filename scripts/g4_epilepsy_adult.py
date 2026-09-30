#!/usr/bin/env python3
"""G4 (adult): IED-like event detection with the G2 framework (NEW).

Arrays, gains, noise model and calibration are those of G2 (configs/g2_adult.toml); the recordings
are simulated in the time domain (opmsquid.ied) with one noise realization per segment shared by
every array, and identical events (source, strength, morphology, time) injected into every array.

Detectors (opmsquid.detection), per array/channel set:
* oracle: knows the event topography, waveform and peak time; threshold for a per-trial
  false-positive probability (null samples at random times);
* practical: scans all times, a template bank (three spike-wave stretches) and a dictionary of
  candidate cortical sources that does not contain the true sources; thresholds for 1 and 0.2
  false events per minute set on 20 min of calibration null data, frozen, then checked on 20 min of
  held-out null data.
Events: 72 target locations stratified by depth x orientation; focal dipoles at 6 strengths x 3
morphologies, and 10-mm patches (fixed total) at 6 strengths. Detection probability with Wilson
intervals, strength for 50 % detection per stratum, and sensitivity vs false events per minute.

Configuration: configs/g4_epilepsy.toml. Outputs: results/g4/.
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
from opmsquid import background, detection, environment, g2, goldenholz, ied, io  # noqa: E402

OUT = ROOT / "results" / "g4"
DEPTH_BANDS = ((10.0, 20.0), (20.0, 30.0), (30.0, 45.0), (45.0, 70.0))
ORIENT_BANDS = ((0.0, 30.0), (30.0, 60.0), (60.0, 90.1))


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def stratified_locations(st, n, rng):
    """Up to n/12 targets per depth x orientation stratum (fewer if the stratum is small)."""
    per = n // (len(DEPTH_BANDS) * len(ORIENT_BANDS))
    out, strata = [], []
    for i, (d0, d1) in enumerate(DEPTH_BANDS):
        for j, (o0, o1) in enumerate(ORIENT_BANDS):
            pool = np.flatnonzero((st.src.depth_mm >= d0) & (st.src.depth_mm < d1) & (st.src.orientation_deg >= o0)
                                  & (st.src.orientation_deg < o1))
            pick = rng.choice(pool, min(per, len(pool)), replace=False) if len(pool) else []
            out += list(pick)
            strata += [(i, j)] * len(pick)
    return np.array(out), np.array(strata)


def main(overrides: dict | None = None):
    t_start = time.time()
    mne.set_log_level("WARNING")
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads((ROOT / "configs" / "g4_epilepsy.toml").read_text())
    for section, values in (overrides or {}).items():
        cfg[section].update(values)
    cfg2 = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    sim, ev, dc = cfg["simulation"], cfg["events"], cfg["detector"]
    st = G.Study(cfg2)
    rng = np.random.default_rng(sim["seed"])
    arrays = {k: v for k, v in g2.build_arrays(st.subject, st.dig).items() if k in ("squid", "opm99", "opm_dense")}
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
    specs = {}
    for name, a in arrays.items():
        asd = (np.array([g2.SQUID_ASD[k] for k in a.kinds]) if name == "squid"
               else np.full(a.n, sim["opm_asd_fT_per_rtHz"] * 1e-15))
        specs[name] = ied.ArraySpec(G_g[name], environment.external_basis(a.info, env.r0, a.coil_def), asd)
    gen = ied.NoiseGenerator(specs, brain_scale * st.src.grid_area, env.coef_timecourses, env.coef_cov, st.filt,
                             decimate=sim["decimate"], env_fs=env.sfreq)
    fs = gen.fs_out
    log(f"generator ready: output rate {fs:.1f} Hz, brain scale {brain_scale:.3e}")

    # detectors: array/channel set -> row mask
    det_sets = {}
    for key in dc["arrays"]:
        name, cs = key.split("/")
        det_sets[key] = (name, g2.channel_sets(arrays[name])[cs])
    templates = {s: ied.filtered_template(st.filt, s, sim["decimate"]) for s in ev["stretches"]}

    # event sources and the candidate dictionary
    loc, strata = stratified_locations(st, ev["n_locations"], rng)
    tgt = st.src.target[loc]
    members = goldenholz.geodesic_patches(st.cortex.adjacency, tgt, ev["patch_radius_mm"] * 1e-3, st.cortex.usable)
    topo = {}
    for name, a in arrays.items():
        full, col = g2.fullres_matrix(a, st.subject, st.cortex, g2.FULLRES_JOBS[name])
        unit = goldenholz.patch_topographies(full, members, col, st.cortex.area)
        area = np.array([st.cortex.area[m].sum() for m in members])
        for i in range(len(loc)):
            topo[(name, i, "focal")] = G_t[name][:, loc[i]]
            topo[(name, i, "patch")] = unit[:, i] / area[i]  # per unit total scalar moment
        del full
    valid = np.flatnonzero(st.cortex.usable)
    cand = valid[goldenholz.poisson_disk(st.cortex.rr[valid], dc["dictionary_spacing_mm"] * 1e-3, rng)]
    cand = np.setdiff1d(cand, tgt)
    cand_gain = {name: g2.gains(a, st.subject, st.cortex, cand, g2.FULLRES_JOBS[name]) for name, a in arrays.items()}
    log(f"{len(loc)} event locations, {len(cand)} dictionary candidates")

    # baseline null -> whiteners; independent calibration null -> detector normalisation, oracle
    # null distribution and thresholds (no in-sample optimism); both then frozen
    seg_s = sim["segment_s"]

    def null_data(minutes):
        acc = {name: [] for name in arrays}
        for _ in range(int(round(minutes * 60 / seg_s))):
            for name, y in gen.segment(seg_s, rng).items():
                acc[name].append(y)
        return {name: np.concatenate(v, axis=1) for name, v in acc.items()}

    base = null_data(cfg["null"]["baseline_min"])
    whiteners = {key: detection.whitener_from_null(base[name][m]) for key, (name, m) in det_sets.items()}
    del base
    cal = null_data(cfg["null"]["calibration_min"])
    refr, tol = int(round(dc["refractory_s"] * fs)), int(round(dc["hit_tolerance_s"] * fs))
    dets, oracles, thr, zcrit = {}, {}, {}, {}
    minutes_cal = cfg["null"]["calibration_min"]
    for key, (name, m) in det_sets.items():
        null = cal[name][m]
        w = whiteners[key]
        dets[key] = detection.ScanDetector.build(null, w, templates, cand_gain[name][m], refr)
        tops = {(i, fam): topo[(name, i, fam)][m] for i in range(len(loc)) for fam in ("focal", "patch")}
        oracles[key], samples = detection.Oracle.calibrate(null, w, templates, tops, rng=rng)
        zcrit[key] = float(np.quantile(np.concatenate(list(samples.values())), 1 - dc["oracle_alpha"]))
        stat, _ = dets[key].statistic(null)
        _, h = dets[key].events(stat)
        thr[key] = {f"{r:g}": detection.threshold_for_rate(h, minutes_cal, r) for r in dc["operating_points_per_min"]}
        log(f"calibrated {key}: rank {w.rank}, thresholds {thr[key]}, oracle z_crit {zcrit[key]:.2f}")
    del cal

    # held-out null: false events per minute over a threshold sweep
    n_held = int(round(cfg["null"]["heldout_min"] * 60 / seg_s))
    held_heights = {key: [] for key in det_sets}
    for _ in range(n_held):
        seg = gen.segment(seg_s, rng)
        for key, (name, m) in det_sets.items():
            stat, _ = dets[key].statistic(seg[name][m])
            held_heights[key].append(dets[key].events(stat)[1])
    held_heights = {k: np.concatenate(v) for k, v in held_heights.items()}
    minutes_held = cfg["null"]["heldout_min"]
    log("held-out null done")

    # events: all combinations, shuffled, 14 per 30-s segment, identical across arrays
    events = [(i, "focal", s, x) for i in range(len(loc)) for s in ev["strengths_nAm"] for x in ev["stretches"]]
    events += [(i, "patch", s, 1.0) for i in range(len(loc)) for s in ev["strengths_nAm"]]
    order = rng.permutation(len(events))
    per_seg = int((seg_s - 1.5) // sim["event_spacing_s"])  # events at 1.5, 3.5, ... s, the last >= 2 s before the end
    times = (1.5 + sim["event_spacing_s"] * np.arange(per_seg)) * fs
    rec = {key: dict(oracle_z=np.zeros(len(events)), near_height=np.zeros(len(events))) for key in det_sets}
    for start in range(0, len(order), per_seg):
        batch = order[start:start + per_seg]
        seg = gen.segment(seg_s, rng)
        for name in arrays:
            for e_idx, t0 in zip(batch, times):
                i, fam, s, x = events[e_idx]
                tpl, pk = templates[x]
                ied.inject(seg[name], topo[(name, i, fam)], tpl, pk, s * 1e-9, int(round(t0)))
        for key, (name, m) in det_sets.items():
            y = seg[name][m]
            stat, _ = dets[key].statistic(y)
            for e_idx, t0 in zip(batch, times):
                i, fam, s, x = events[e_idx]
                t0 = int(round(t0))
                rec[key]["near_height"][e_idx] = float(stat[max(t0 - tol, 0):t0 + tol + 1].max())
                rec[key]["oracle_z"][e_idx] = oracles[key].statistic(y, (i, fam), x, topo[(name, i, fam)][m], t0)
        if (start // per_seg) % 20 == 0:
            log(f"events {start + len(batch)}/{len(events)}")

    # summaries
    ev_arr = np.array([(i, s, x) for i, _, s, x in events], float)
    fam_arr = np.array([f for _, f, _, _ in events])
    depth_band = strata[ev_arr[:, 0].astype(int), 0]
    summary = dict(status="NEW (G4 adult: IED detection, G2 noise model in the time domain)", config=cfg, fs_out=fs,
                   n_events=len(events), n_locations=int(len(loc)), n_dictionary=int(len(cand)),
                   locations=[dict(vertex=int(st.cortex.vertno[v]), hemi=int(st.cortex.hemi[v]), depth_mm=float(st.src.depth_mm[li]),
                                   orientation_deg=float(st.src.orientation_deg[li]), lobe=str(st.src.lobe[li]), stratum=[int(a) for a in strata[k]])
                              for k, (v, li) in enumerate(zip(tgt, loc))],
                   detectors={})
    for key in det_sets:
        r = rec[key]
        out = dict(thresholds=thr[key], oracle_z_crit=zcrit[key],
                   heldout_false_per_min={op: float(np.sum(held_heights[key] > v) / minutes_held) for op, v in thr[key].items()},
                   curves={})
        for fam in ("focal", "patch"):
            for x in (ev["stretches"] if fam == "focal" else [1.0]):
                for db in range(len(DEPTH_BANDS)):
                    for s in ev["strengths_nAm"]:
                        sel = (fam_arr == fam) & (ev_arr[:, 2] == x) & (ev_arr[:, 1] == s) & (depth_band == db)
                        n = int(sel.sum())
                        row = dict(n=n)
                        k = int(np.sum(r["oracle_z"][sel] > zcrit[key]))
                        row["oracle"] = dict(p=k / n if n else None, ci=wilson(k, n))
                        for op, v in thr[key].items():
                            k = int(np.sum(r["near_height"][sel] > v))
                            row[f"practical@{op}"] = dict(p=k / n if n else None, ci=wilson(k, n))
                        out["curves"][f"{fam}/stretch{x:g}/depth{db}/{s:g}nAm"] = row
        # sensitivity vs false events per minute (focal, stretch 1, 40 and 80 nAm pooled over locations)
        roc = {}
        for s in (40.0, 80.0):
            sel = (fam_arr == "focal") & (ev_arr[:, 2] == 1.0) & (ev_arr[:, 1] == s)
            ths = np.quantile(held_heights[key], np.linspace(0.5, 0.9999, 60))
            roc[f"{s:g}nAm"] = dict(false_per_min=[float(np.sum(held_heights[key] > t) / minutes_held) for t in ths],
                                    sensitivity=[float(np.mean(r["near_height"][sel] > t)) for t in ths])
        out["roc"] = roc
        # strength for 50 % detection per depth band (focal, all stretches pooled), log interpolation
        s50 = {}
        for mode in ["oracle"] + [f"practical@{op}" for op in thr[key]]:
            for db in range(len(DEPTH_BANDS)):
                p = []
                for s in ev["strengths_nAm"]:
                    sel = (fam_arr == "focal") & (ev_arr[:, 1] == s) & (depth_band == db)
                    crit = zcrit[key] if mode == "oracle" else thr[key][mode.split("@")[1]]
                    val = r["oracle_z"][sel] if mode == "oracle" else r["near_height"][sel]
                    p.append(float(np.mean(val > crit)) if sel.any() else np.nan)
                p = np.array(p)
                ls = np.log(ev["strengths_nAm"])
                above = np.flatnonzero(p >= 0.5)
                if len(above) == 0:
                    v50 = None
                elif above[0] == 0:
                    v50 = float(ev["strengths_nAm"][0])  # at or below the weakest strength
                else:
                    k = above[0]
                    v50 = float(np.exp(np.interp(0.5, [p[k - 1], p[k]], [ls[k - 1], ls[k]])))
                s50[f"{mode}/depth{db}"] = v50
        out["strength_for_50pct_nAm"] = s50
        summary["detectors"][key] = out
    io.write_json(summary, OUT / "g4_adult_summary.json")
    with open(OUT / "g4_adult_events.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["location", "family", "strength_nAm", "stretch", "depth_band"] + [f"{k}_{v}" for k in det_sets for v in ("oracle_z", "near_height")])
        for n, (i, fam, s, x) in enumerate(events):
            wr.writerow([i, fam, s, x, int(depth_band[n])] + [f"{rec[k][v][n]:.3f}" for k in det_sets for v in ("oracle_z", "near_height")])
    figures(summary, ev)
    log(f"done in {time.time() - t_start:.0f} s")


def figures(summary, ev):
    colors = {"squid/combined": "k", "squid/grad": "tab:red", "squid/mag": "tab:blue", "opm99/opm": "tab:green", "opm_dense/opm": "tab:purple"}
    fig, axs = plt.subplots(2, len(DEPTH_BANDS), figsize=(16, 7.5), sharey=True)
    for row, mode in enumerate(("oracle", "practical@1")):
        for db, (d0, d1) in enumerate(DEPTH_BANDS):
            ax = axs[row, db]
            for key, det in summary["detectors"].items():
                p, lo, hi = [], [], []
                for s in ev["strengths_nAm"]:
                    c = det["curves"][f"focal/stretch1/depth{db}/{s:g}nAm"][mode]
                    p.append(np.nan if c["p"] is None else c["p"])
                    lo.append(np.nan if c["ci"][0] is None else c["ci"][0])
                    hi.append(np.nan if c["ci"][1] is None else c["ci"][1])
                ax.plot(ev["strengths_nAm"], p, "o-", color=colors[key], label=key, ms=3)
                ax.fill_between(ev["strengths_nAm"], lo, hi, color=colors[key], alpha=0.1)
            ax.set_xscale("log")
            ax.set_title(f"{mode.replace('@1', ' (1 false event/min)')}, depth {d0:g}-{d1:g} mm", fontsize=8)
            ax.set_xlabel("focal source strength [nAm]")
            if db == 0:
                ax.set_ylabel("detection probability (stretch 1.0)")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("G4 adult: IED detection vs strength (Wilson 95 % bands; oracle knows source and time)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G4_detection.png", dpi=150)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, s in zip(axs, ("40", "80")):
        for key, det in summary["detectors"].items():
            r = det["roc"][f"{s}nAm"]
            ax.plot(r["false_per_min"], r["sensitivity"], "-", color=colors[key], label=key)
        ax.set_xscale("log")
        ax.set_xlabel("false events per minute (held-out null)")
        ax.set_ylabel("sensitivity")
        ax.set_title(f"focal {s} nAm, stretch 1.0, all locations", fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("G4 adult: practical detector, sensitivity vs false events per minute", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G4_roc.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
