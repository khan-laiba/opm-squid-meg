#!/usr/bin/env python3
"""G4: IED-like event detection with the G2 framework (NEW), adult (this script) and, through
``Context``, any anatomy with its arrays (scripts/g4_epilepsy_pediatric.py).

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
import pickle
import sys
import time
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

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
STATE = ROOT / "cache" / "g4" / "adult_state.pkl"
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


@dataclass
class Context:
    """One anatomy with its arrays ('squid', 'opm_matched', 'opm_dense'): what the detection and
    localization studies need. ``gain(name, vertices)`` gives lead fields at global vertices and
    ``patch_topographies(name, members)`` the summed lead fields of patches (per unit density)."""

    label: str
    subject: object
    cortex: object
    src: object
    filt: object
    arrays: dict
    G_t: dict
    G_g: dict
    brain_scale: float
    env: object
    gain: Callable
    patch_topographies: Callable
    notes: str = ""


def adult_context(cfg2, names=("squid", "opm_matched", "opm_dense")) -> Context:
    """The G2 adult: sample subject, measured head position, full-resolution lead fields."""
    st = G.Study(cfg2)
    arrays = {k: v for k, v in g2.build_arrays(st.subject, st.dig).items() if k in names}
    squid = arrays["squid"]
    bads = cfg2["sensors"]["bads"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads = good & (squid.kinds == "grad")
    meas = g2.measured_noise(squid.info, st.filt, bads)
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        G_t[name], G_g[name] = st.gains(a)
    brain_scale = background.calibrate(st.unit_brain(G_g["squid"]), grads, float(np.nanmedian(meas["brain"][grads])))

    def gain(name, vertices):
        return g2.gains(arrays[name], st.subject, st.cortex, vertices, g2.FULLRES_JOBS[name])

    def patch_topographies(name, members):
        full, col = g2.fullres_matrix(arrays[name], st.subject, st.cortex, g2.FULLRES_JOBS[name])
        return goldenholz.patch_topographies(full, members, col, st.cortex.area)

    return Context("adult", st.subject, st.cortex, st.src, st.filt, arrays, G_t, G_g, brain_scale, meas["environment"], gain,
                   patch_topographies, "MNE sample subject, measured head position (G2)")


def stratified_locations(st, n, rng):
    """Up to n/12 targets per depth x orientation stratum (fewer if the stratum is small), never on
    the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum, not cortex)."""
    per = n // (len(DEPTH_BANDS) * len(ORIENT_BANDS))
    cortical = ~np.char.endswith(st.src.region.astype(str), "unknown")
    out, strata = [], []
    for i, (d0, d1) in enumerate(DEPTH_BANDS):
        for j, (o0, o1) in enumerate(ORIENT_BANDS):
            pool = np.flatnonzero(cortical & (st.src.depth_mm >= d0) & (st.src.depth_mm < d1) & (st.src.orientation_deg >= o0)
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
    ctx = adult_context(cfg2)
    state = simulate(ctx, cfg, np.random.default_rng(cfg["simulation"]["seed"]))
    save_state(state, STATE)
    summarise(state)
    log(f"done in {time.time() - t_start:.0f} s")


def save_state(state, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(state, fh)


def simulate(ctx: Context, cfg: dict, rng: np.random.Generator) -> dict:
    """Null data, frozen detectors and the event simulation for one anatomy (see the module
    docstring); returns the checkpoint state that ``summarise`` reads."""
    sim, ev, dc = cfg["simulation"], cfg["events"], cfg["detector"]
    arrays, G_t, G_g, env = ctx.arrays, ctx.G_t, ctx.G_g, ctx.env
    specs = {}
    for name, a in arrays.items():
        asd = (np.array([g2.SQUID_ASD[k] for k in a.kinds]) if name == "squid"
               else np.full(a.n, sim["opm_asd_fT_per_rtHz"] * 1e-15))
        specs[name] = ied.ArraySpec(G_g[name], environment.external_basis(a.info, env.r0, a.coil_def), asd)
    gen = ied.NoiseGenerator(specs, ctx.brain_scale * ctx.src.grid_area, env.coef_timecourses, env.coef_cov, ctx.filt,
                             decimate=sim["decimate"], env_fs=env.sfreq)
    fs = gen.fs_out
    log(f"{ctx.label}: generator ready: output rate {fs:.1f} Hz, brain scale {ctx.brain_scale:.3e}")

    # detectors: array/channel set -> row mask
    det_sets = {}
    for key in dc["arrays"]:
        name, cs = key.split("/")
        det_sets[key] = (name, g2.channel_sets(arrays[name])[cs])
    templates = {s: ied.filtered_template(ctx.filt, s, sim["decimate"]) for s in ev["stretches"]}

    # event sources and the candidate dictionary
    loc, strata = stratified_locations(ctx, ev["n_locations"], rng)
    tgt = ctx.src.target[loc]
    members = goldenholz.geodesic_patches(ctx.cortex.adjacency, tgt, ev["patch_radius_mm"] * 1e-3, ctx.cortex.usable)
    topo = {}
    area = np.array([ctx.cortex.area[m].sum() for m in members])
    for name in arrays:
        unit = ctx.patch_topographies(name, members)
        for i in range(len(loc)):
            topo[(name, i, "focal")] = G_t[name][:, loc[i]]
            topo[(name, i, "patch")] = unit[:, i] / area[i]  # per unit total scalar moment
    valid = np.flatnonzero(ctx.cortex.usable)
    cand = valid[goldenholz.poisson_disk(ctx.cortex.rr[valid], dc["dictionary_spacing_mm"] * 1e-3, rng)]
    cand = np.setdiff1d(cand, tgt)
    cand_gain = {name: ctx.gain(name, cand) for name in arrays}
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
    log("held-out null done")

    # events: all combinations, shuffled, 14 per 30-s segment, identical across arrays
    events = [(i, "focal", s, x) for i in range(len(loc)) for s in ev["strengths_nAm"] for x in ev["stretches"]]
    events += [(i, "patch", s, 1.0) for i in range(len(loc)) for s in ev["strengths_nAm"]]
    order = rng.permutation(len(events))
    per_seg = int((seg_s - 1.5) // sim["event_spacing_s"])  # events at 1.5, 3.5, ... s, the last >= 2 s before the end
    times = (1.5 + sim["event_spacing_s"] * np.arange(per_seg)) * fs
    # injected events are scored with the detector's own emitted events (local maxima >= refractory apart,
    # as on the null data): the height of the highest emitted event within +/- tol of the true peak (0 if none)
    rec = {key: dict(oracle_z=np.zeros(len(events)), event_height=np.zeros(len(events))) for key in det_sets}
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
            emitted = dets[key].events(stat)
            for e_idx, t0 in zip(batch, times):
                i, fam, s, x = events[e_idx]
                t0 = int(round(t0))
                rec[key]["event_height"][e_idx] = detection.event_height(*emitted, t0, tol)
                rec[key]["oracle_z"][e_idx] = oracles[key].statistic(y, (i, fam), x, topo[(name, i, fam)][m], t0)
        if (start // per_seg) % 20 == 0:
            log(f"events {start + len(batch)}/{len(events)}")

    # checkpoint, then summaries (``--resummarise`` redoes them from the checkpoint)
    return dict(cfg=cfg, label=ctx.label, context=ctx.notes, fs=fs, events=events, strata=strata, rec=rec, held_heights=held_heights,
                thr=thr, zcrit=zcrit, minutes_held=cfg["null"]["heldout_min"], n_dictionary=int(len(cand)),
                simulated_at_commit=io.RUN_COMMIT,
                locations=[dict(vertex=int(ctx.cortex.vertno[v]), hemi=int(ctx.cortex.hemi[v]), depth_mm=float(ctx.src.depth_mm[li]),
                                orientation_deg=float(ctx.src.orientation_deg[li]), lobe=str(ctx.src.lobe[li]),
                                region=str(ctx.src.region[li]), stratum=[int(a) for a in strata[k]])
                           for k, (v, li) in enumerate(zip(tgt, loc))])


def summarise(state):
    cfg, events, strata, rec = state["cfg"], state["events"], state["strata"], state["rec"]
    label = state.get("label", "adult")
    ev = cfg["events"]
    strengths = np.array(ev["strengths_nAm"], float)
    thr, zcrit, held, minutes = state["thr"], state["zcrit"], state["held_heights"], state["minutes_held"]
    loc_i = np.array([e[0] for e in events])
    fam = np.array([e[1] for e in events])
    stren = np.array([e[2] for e in events], float)
    stretch = np.array([e[3] for e in events], float)
    band = strata[loc_i, 0]
    keys = list(rec)

    def detected(key, mode):
        return rec[key]["oracle_z"] > zcrit[key] if mode == "oracle" else rec[key]["event_height"] > thr[key][mode.split("@")[1]]

    modes = ["oracle"] + [f"practical@{op}" for op in thr[keys[0]]]
    rng = np.random.default_rng(7)
    tables = {}  # (detector, mode) -> detection rate per location and strength (focal)
    summary = dict(status=f"NEW (G4 {label}: IED detection, G2 noise model in the time domain)", config=cfg, fs_out=state["fs"],
                   anatomy=state.get("context", ""),
                   simulated_at_commit=state.get("simulated_at_commit"),
                   n_events=len(events), n_locations=len(state["locations"]), n_dictionary=state["n_dictionary"],
                   locations=state["locations"], detectors={}, paired={})
    for key in keys:
        out = dict(thresholds=thr[key], oracle_z_crit=zcrit[key],
                   heldout_false_per_min={op: float(np.sum(held[key] > v) / minutes) for op, v in thr[key].items()}, curves={},
                   strength_for_50pct_nAm={})
        for mode in modes:
            det = detected(key, mode)
            for f_ in ("focal", "patch"):
                for db in range(len(DEPTH_BANDS)):
                    for s_ in strengths:
                        sel = (fam == f_) & (band == db) & (stren == s_)
                        n, k = int(sel.sum()), int(det[sel].sum())
                        out["curves"][f"{mode}/{f_}/depth{db}/{s_:g}nAm"] = dict(n=n, p=k / n if n else None, ci=wilson(k, n))
                        if f_ == "focal":
                            for x in ev["stretches"]:
                                sx = sel & (stretch == x)
                                nx, kx = int(sx.sum()), int(det[sx].sum())
                                out["curves"][f"{mode}/focal_stretch{x:g}/depth{db}/{s_:g}nAm"] = dict(n=nx, p=kx / nx if nx else None,
                                                                                                    ci=wilson(kx, nx))
            # strength for 50 % detection (focal, stretches pooled) with a bootstrap over locations
            n_loc = int(loc_i.max()) + 1
            table = np.full((n_loc, len(strengths)), np.nan)  # detection rate per location and strength
            for L in range(n_loc):
                for j, s_ in enumerate(strengths):
                    m = (fam == "focal") & (loc_i == L) & (stren == s_)
                    if m.any():
                        table[L, j] = det[m].mean()
            tables[(key, mode)] = table
            for db in range(len(DEPTH_BANDS)):
                locs = np.unique(loc_i[band == db])
                if len(locs) == 0:
                    continue
                sub = table[locs]
                bnd = np.array([detection.s50_bounds(sub[rng.integers(0, len(locs), len(locs))].mean(axis=0), strengths)
                                for _ in range(1000)])
                boot = np.where(np.isfinite(bnd[:, 1]), np.maximum(bnd[:, 0], strengths[0]), np.inf)
                # resamples that never reach 50 % sort above every tested strength, those already at 50 % at
                # the weakest strength below it; a limit that falls among them is open (None)
                q = np.percentile(np.where(np.isfinite(boot), boot, 1e9), [2.5, 97.5])
                at_weakest = float(np.mean(bnd[:, 0] == 0.0))
                out["strength_for_50pct_nAm"][f"{mode}/depth{db}"] = dict(
                    value=detection.s50_from(sub.mean(axis=0), strengths),
                    ci95=[None if at_weakest >= 0.025 else float(q[0]), float(q[1]) if q[1] <= strengths[-1] else None],
                    share_resamples_not_reached=float(np.mean(~np.isfinite(boot))), share_resamples_at_weakest=at_weakest)
        # sensitivity vs false events per minute: thresholds at every held-out event height
        h = np.sort(held[key])[::-1]
        ths = h[:min(len(h), int(10 * minutes) + 1)]
        roc = {}
        for lab, sel in (("superficial_10-30mm_40nAm", (fam == "focal") & (band <= 1) & (stren == 40.0)),
                         ("deep_30-70mm_160nAm", (fam == "focal") & (band >= 2) & (stren == 160.0))):
            roc[lab] = dict(false_per_min=[float(np.sum(held[key] > t) / minutes) for t in ths],
                            sensitivity=[float(np.mean(rec[key]["event_height"][sel] > t)) for t in ths], n_events=int(sel.sum()))
        out["roc"] = roc
        summary["detectors"][key] = out
    # paired comparisons on identical events: OPM vs each Neuromag channel set
    for a in [k for k in keys if not k.startswith("squid")]:
        for b in [k for k in keys if k.startswith("squid")]:
            for mode in ("oracle", "practical@1"):
                da, db_ = detected(a, mode), detected(b, mode)
                res = {}
                for band_i in range(len(DEPTH_BANDS)):
                    sel = (fam == "focal") & (band == band_i)
                    only_a, only_b = int(np.sum(da[sel] & ~db_[sel])), int(np.sum(~da[sel] & db_[sel]))
                    from scipy.stats import binomtest

                    pval = binomtest(only_a, only_a + only_b, 0.5).pvalue if only_a + only_b else 1.0
                    stat_a = rec[a]["oracle_z" if mode == "oracle" else "event_height"][sel]
                    stat_b = rec[b]["oracle_z" if mode == "oracle" else "event_height"][sel]
                    diff = stat_a - stat_b
                    # events share locations (6 strengths x 3 morphologies each): the location is the unit
                    loc_s = loc_i[sel]
                    locs = np.unique(loc_s)
                    d_loc = np.array([np.sum((da & ~db_)[sel][loc_s == L]) - np.sum((~da & db_)[sel][loc_s == L]) for L in locs])
                    per = [diff[loc_s == L] for L in locs]
                    bs = [np.median(np.concatenate([per[j] for j in rng.integers(0, len(per), len(per))])) for _ in range(1000)]
                    # S50 ratio (Neuromag / OPM; > 1: the OPM array needs less strength), paired over the same
                    # location resamples; None when a resample does not reach 50 % within the tested range
                    # an S50 outside the tested strengths is only bounded, so each resample gives a ratio interval
                    # and every resample is kept (detection.censored_interval)
                    ta, tb = tables[(a, mode)][locs], tables[(b, mode)][locs]
                    pairs = []
                    for _ in range(1000):
                        r_ = rng.integers(0, len(locs), len(locs))
                        pairs.append((detection.s50_bounds(ta[r_].mean(axis=0), strengths), detection.s50_bounds(tb[r_].mean(axis=0), strengths)))
                    rb = np.array([detection.ratio_bounds(bo, bs_) for bo, bs_ in pairs])
                    above_o = np.array([not np.isfinite(bo[1]) for bo, _ in pairs])
                    above_s = np.array([not np.isfinite(bs_[1]) for _, bs_ in pairs])
                    weakest = np.array([bo[0] == 0.0 or bs_[0] == 0.0 for bo, bs_ in pairs])
                    b0 = (detection.s50_bounds(ta.mean(axis=0), strengths), detection.s50_bounds(tb.mean(axis=0), strengths))
                    v_lo, v_hi = detection.ratio_bounds(*b0)
                    s50_ratio = dict(value=v_lo if v_lo == v_hi else None, value_censored=detection.censoring_label(*b0),
                                     value_bounds=None if v_lo == v_hi else [v_lo if v_lo > 0 else None, v_hi if np.isfinite(v_hi) else None],
                                     ci95=detection.censored_interval(rb[:, 0], rb[:, 1]),
                                     share_resamples_neither_reached=float(np.mean(above_o & above_s)),
                                     share_resamples_opm_not_reached=float(np.mean(above_o & ~above_s)),
                                     share_resamples_squid_not_reached=float(np.mean(above_s & ~above_o)),
                                     share_resamples_reached_at_weakest=float(np.mean(weakest)))
                    res[f"depth{band_i}"] = dict(n=int(sel.sum()), n_locations=int(len(locs)), detected_only_opm=only_a,
                                                 detected_only_squid=only_b, locations_favouring_opm=int(np.sum(d_loc > 0)),
                                                 locations_favouring_squid=int(np.sum(d_loc < 0)), location_sign_flip_p=detection.sign_flip_p(d_loc),
                                                 mcnemar_exact_p_event_level=float(pval), median_stat_difference=float(np.median(diff)),
                                                 ci95_location_bootstrap=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                                                 s50_ratio_squid_over_opm=s50_ratio, p_values="uncorrected")
                summary["paired"][f"{a}_vs_{b}/{mode}"] = res
    io.write_json(summary, OUT / f"g4_{label}_summary.json")
    with open(OUT / f"g4_{label}_events.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["location", "family", "strength_nAm", "stretch", "depth_band"] + [f"{k}_{v}" for k in keys for v in ("oracle_z", "event_height")])
        for n, (i, f_, s_, x) in enumerate(events):
            wr.writerow([i, f_, s_, x, int(band[n])] + [f"{rec[k][v][n]:.3f}" for k in keys for v in ("oracle_z", "event_height")])
    figures(summary, ev, label)


def figures(summary, ev, label="adult"):
    suffix = "" if label == "adult" else f"_{label}"
    colors = {"squid/combined": "k", "squid/grad": "tab:red", "squid/mag": "tab:blue", "opm_matched/opm": "tab:green", "opm_dense/opm": "tab:purple"}
    fig, axs = plt.subplots(2, len(DEPTH_BANDS), figsize=(16, 7.5), sharey=True)
    for row, mode in enumerate(("oracle", "practical@1")):
        for db, (d0, d1) in enumerate(DEPTH_BANDS):
            ax = axs[row, db]
            for key, det in summary["detectors"].items():
                c = [det["curves"][f"{mode}/focal/depth{db}/{s:g}nAm"] for s in ev["strengths_nAm"]]
                p = [np.nan if x["p"] is None else x["p"] for x in c]
                lo = [np.nan if x["ci"][0] is None else x["ci"][0] for x in c]
                hi = [np.nan if x["ci"][1] is None else x["ci"][1] for x in c]
                ax.plot(ev["strengths_nAm"], p, "o-", color=colors[key], label=key, ms=3)
                ax.fill_between(ev["strengths_nAm"], lo, hi, color=colors[key], alpha=0.06 if key == "squid/combined" else 0.1, lw=0)
            ax.set_xscale("log")
            ax.set_title(f"{mode.replace('@1', ' (1 false event/min)')}, depth {d0:g}-{d1:g} mm (n = {c[0]['n']}/point)", fontsize=8)
            ax.set_xlabel("focal source strength [nAm]")
            if db == 0:
                ax.set_ylabel("detection probability (3 morphologies pooled)")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(f"G4 {label}: IED detection vs strength (Wilson 95 % bands; the oracle knows source and time)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"Figure_G4_detection{suffix}.png", dpi=150)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, lab in zip(axs, ("superficial_10-30mm_40nAm", "deep_30-70mm_160nAm")):
        for key, det in summary["detectors"].items():
            r = det["roc"][lab]
            ax.step(r["false_per_min"], r["sensitivity"], where="post", color=colors[key], label=key)
        ax.set_xscale("log")
        ax.set_xlim(0.05, 10)
        ax.axvline(1.0, color="0.6", lw=0.8, ls=":")
        ax.set_xlabel("false events per minute (held-out null)")
        ax.set_ylabel("sensitivity")
        ax.set_title(f"{lab.replace('_', ' ')} (focal, 3 morphologies; n = {r['n_events']})", fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle(f"G4 {label}: practical detector, sensitivity vs false events per minute", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"Figure_G4_roc{suffix}.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    if "--resummarise" in sys.argv:
        with open(STATE, "rb") as fh:
            summarise(pickle.load(fh))
    else:
        main()
