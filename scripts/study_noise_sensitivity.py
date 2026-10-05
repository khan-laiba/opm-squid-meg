#!/usr/bin/env python3
"""Adult noise-model sensitivity analyses (revision; NEW).

The analyses: (a) the near-skull cortex, which the primary model omits, added to the background, to
the targets or both; (b) a coloured OPM noise model (1/f plus white, the corner frequency swept)
beside a uniformly worse white sensor; (c) bounded cardiac and ocular far-field terms for both
systems; (d) the primary adult result over the OPM white-noise sweep; and the joint runs. They test
structured and coloured alternatives to the white-noise model and report the resulting range of the
comparison (supplementary text, sections C.6 and C.7).

Everything is the adult comparison of scripts/g2_adult_comparison.py (MNE sample subject, its
measured head position, the same arrays, 7,661 cortical targets of 10 nAm, 7-mm background grid
calibrated on the measured Neuromag gradiometer brain noise, 8-term room field from the empty-room
recording, 1-40 Hz) with one element changed at a time, plus joint runs with the alternatives
together:

  near_skull  (a) the cortex within 4 mm of the inner skull (A-BEM-DIST), at full resolution, in
              the background, as targets, and both; primary 3-layer lead fields from the cached
              full-resolution matrices above a 2-mm floor, checked with a 1-layer BEM whose inner
              skull is subdivided once (the code base's refinement check), down to 0 mm
  coloured    (b) OPM noise PSD w^2 (1 + f_c/f), corners 1, 3, 10 Hz, at every white level w of the
              sweep; frequency-resolved whitening over sub-bands of the 1-40 Hz band (brain scale
              and room field refitted per sub-band from the recordings, as scripts/g2_band_sensitivity.py
              does per band), the band-variance alternative (spatial whitening only) and, as the
              Neuromag counterpart, its sensor noise from the empty-room spectrum per sub-band
  far_field   (c) cardiac (chest) and ocular (eye) sources: field energy left by the 8-term
              projection in each array, and detectability with them added at two declared levels,
              on top of the published background ('added') or with it refitted ('recalibrated')
  sweep       (d) the primary result over the OPM white-noise sweep 7-30 fT/sqrt(Hz), with and
              without projection, and the white level at which each median ratio equals 1
  joint       the alternatives at once, frequency-resolved, every white level: near-skull background
              + 10-Hz corner (the most adverse combination), and with cardiac and ocular sources
              filling the magnetometer brain-noise shortfall or at the room-field level

Ratios follow results/g2/g2_summary.json: median over targets of log2(d_OPM / d_Neuromag) with a
95 % CI from 1,000 parcel-bootstrap resamples (scripts/g2_adult_comparison.compare), d the
known-topography detectability sqrt(s^T C^+ s) with the oracle covariance.

Configuration: configs/g2_adult.toml (shared) and configs/g2_noise_sensitivity.toml (the declared
alternatives). Outputs (default results/g2_noise_sensitivity/, --out for tests):
noise_sensitivity_summary.json, Figure_noise_sensitivity.png, Figure_noise_sensitivity_depth.png.
"""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
import resource
import sys
import time
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402

import g2_adult_comparison as G  # noqa: E402
from opmsquid import (anatomy, background, environment, forward, g2, ied, io, metrics, neuromag, noisemodel, opm,  # noqa: E402
                      paths, plotting)

OUT = ROOT / "results" / "g2_noise_sensitivity"
STATUS = "NEW (revision: adult noise-model sensitivity; near-skull cortex, coloured OPM noise, far-field sources, noise sweep)"
PURPOSE = ("(a) the near-skull cortex added to the background, to the targets or both; (b) a coloured OPM noise model (1/f "
           "plus white, the corner frequency swept) beside a uniformly worse white sensor; (c) bounded cardiac and ocular "
           "far-field terms for both systems; (d) the primary adult result over the OPM white-noise sweep; the range of "
           "the comparison under these structured and coloured alternatives")
PARTS = ("sweep", "near_skull", "coloured", "far_field", "joint")
MU0_4PI = 1e-7  # mu_0 / (4 pi) [T m / A]
T0 = time.time()


def log(msg):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1e9 if sys.platform == "darwin" else 1e6)
    print(f"[{time.strftime('%H:%M:%S')} +{time.time() - T0:5.0f} s, peak {rss:.2f} GB] {msg}", flush=True)


def peak_rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1e9 if sys.platform == "darwin" else 1e6)


# ----------------------------------------------------------------------------------------------
# shared context: the adult comparison's study, arrays, gains and measured noise
@dataclass
class Ctx:
    cfg: dict
    scfg: dict
    args: argparse.Namespace
    st: G.Study
    arrays: dict
    tsel: np.ndarray  # indices of the analysed targets within the 7,661 (all of them in a full run)
    G_t: dict  # name -> (n_ch, n_selected_targets)
    G_g: dict  # name -> (n_ch, n_grid)
    good: np.ndarray
    grads: np.ndarray
    mags: np.ndarray
    meas: dict
    env: environment.EnvironmentModel
    target_var: float
    unit: dict  # name -> unit-scale background covariance (G2 grid)
    brain_scale: float
    noise_nom: dict
    basis: dict  # name -> (n_ch, 8) room-field basis
    rng: np.random.Generator
    conds: list
    opms: list
    refs: list
    n_boot: int
    extra: dict = field(default_factory=dict)

    @property
    def groups(self):
        return self.st.src.region[self.tsel]

    @property
    def depth(self):
        return self.st.src.depth_mm[self.tsel]

    def topo(self, name):
        return self.G_t[name] * self.st.q


def build_arrays(st):
    """The squid, opm_matched and opm_dense arrays exactly as g2.build_arrays makes them (opm204,
    the channel-budget control, is not needed here)."""
    squid_info = neuromag.load_info("T3")
    arrays = {"squid": g2.Array("squid", squid_info, neuromag.channel_kinds(squid_info), None,
                                dict(sites=102, channels=306, axes="1 mag + 2 planar grad per site"))}
    arrays["opm_matched"] = g2.matched_opm(st.subject, st.dig)
    arrays["opm_dense"] = g2.dense_opm(st.subject, st.dig, "opm_dense")
    return arrays


def setup(cfg, scfg, args) -> Ctx:
    st = G.Study(cfg)
    log(f"study: {st.nt} targets, {len(st.src.grid)} background sources, ENBW {st.enbw:.2f} Hz")
    if args.max_targets and args.max_targets < st.nt:
        tsel = np.sort(np.random.default_rng(args.subset_seed).choice(st.nt, args.max_targets, replace=False))
    else:
        tsel = np.arange(st.nt)
    if args.arrays_cache and Path(args.arrays_cache).exists():
        with open(args.arrays_cache, "rb") as fh:
            arrays = pickle.load(fh)  # development aid only; the lead-field fingerprint check below still applies
        log(f"arrays from {args.arrays_cache} (development cache)")
    else:
        arrays = build_arrays(st)
        if args.arrays_cache:
            with open(args.arrays_cache, "wb") as fh:
                pickle.dump(arrays, fh)
        log("arrays built: " + ", ".join(f"{k} {a.n}" for k, a in arrays.items()))
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        gt, gg = st.gains(a)  # full-resolution 3-layer lead fields, fingerprint- and column-checked
        G_t[name], G_g[name] = np.ascontiguousarray(gt[:, tsel]), gg
    squid = arrays["squid"]
    bads = cfg["sensors"]["bads"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads, mags = good & (squid.kinds == "grad"), good & (squid.kinds == "mag")
    meas = g2.measured_noise(squid.info, st.filt, bads)
    env = meas["environment"]
    target_var = float(np.nanmedian(meas["brain"][grads]))
    unit = {name: st.unit_brain(G_g[name]) for name in arrays}
    brain_scale = background.calibrate(unit["squid"], grads, target_var)
    noise_nom = {name: st.noise(a, G_g[name], brain_scale, env) for name, a in arrays.items()}
    basis = {name: noise_nom[name].ext_basis for name in arrays}
    gs = scfg["general"]
    ctx = Ctx(cfg, scfg, args, st, arrays, tsel, G_t, G_g, good, grads, mags, meas, env, target_var, unit, brain_scale, noise_nom,
              basis, np.random.default_rng(gs["bootstrap_seed"]), list(gs["conditions"]), list(gs["arrays"]), list(gs["references"]),
              int(args.n_boot or gs["n_boot"]))
    log(f"noise model: brain scale {brain_scale:.6e} (calibrated on {int(grads.sum())} gradiometers)")
    return ctx


# ----------------------------------------------------------------------------------------------
# detectability and comparisons (the adult conventions)
def detect(topo, nz: noisemodel.ArrayNoise, array, conds):
    """{(channel set, condition): known-topography detectability} as g2.evaluate computes it."""
    out = {}
    for cond in conds:
        s, c = nz.signal(topo, cond), nz.covariance(cond)
        for cs, m in g2.channel_sets(array).items():
            out[(cs, cond)] = metrics.detectability(s[m], c[np.ix_(m, m)])
    return out


def summarize(c):
    out = dict(c)
    out["ratio"] = float(2 ** c["median_log2"])
    out["ci95_ratio"] = [float(2 ** x) for x in c["ci95"]]
    return out


def compare_all(ctx, d, groups=None, conds=None, opms=None, refs=None, n_boot=None):
    """{'array/ref/cond': median paired ratio with its parcel-bootstrap 95 % CI}."""
    groups = ctx.groups if groups is None else groups
    out = {}
    for cond in conds or ctx.conds:
        for a in opms or ctx.opms:
            for ref in refs or ctx.refs:
                c = G.compare(d[a][("opm", cond)], d["squid"][(ref, cond)], ctx.rng, n_boot or ctx.n_boot, groups=np.asarray(groups))
                if not c["ci_method"].startswith("parcels"):
                    raise RuntimeError("parcel bootstrap expected")
                out[f"{a}/{ref}/{cond}"] = summarize(c)
    return out


def depth_profile(ctx, lr, depth, groups, n_boot=None):
    """Median log2 ratio per 5-mm depth bin (scripts/g2_adult_comparison.DEPTH_EDGES, >= 10 targets),
    with a parcel-bootstrap 95 % CI per bin."""
    out = []
    n_boot = n_boot or ctx.n_boot
    groups = np.asarray(groups)
    for lo, hi in zip(G.DEPTH_EDGES[:-1], G.DEPTH_EDGES[1:]):
        m = (depth >= lo) & (depth < hi) & np.isfinite(lr)
        row = dict(lo=float(lo), hi=float(hi), n=int(m.sum()), median_log2=None, ratio=None, ci95_ratio=None)
        if m.sum() >= 10:
            v, g = lr[m], groups[m]
            per = [v[g == lab] for lab in np.unique(g)]
            boot = np.array([np.median(np.concatenate([per[j] for j in ctx.rng.integers(0, len(per), len(per))])) for _ in range(n_boot)])
            row.update(median_log2=float(np.median(v)), ratio=float(2 ** np.median(v)),
                       ci95_ratio=[float(2 ** np.percentile(boot, 2.5)), float(2 ** np.percentile(boot, 97.5))], n_parcels=len(per))
        out.append(row)
    return out


def profiles(ctx, d, depth=None, groups=None, refs=("combined",), conds=None):
    depth = ctx.depth if depth is None else depth
    groups = ctx.groups if groups is None else groups
    return {f"{a}/{ref}/{cond}": depth_profile(ctx, np.log2(d[a][("opm", cond)] / d["squid"][(ref, cond)]), depth, groups)
            for a in ctx.opms for ref in refs for cond in (conds or ctx.conds)}


def ratio_only(comps):
    return {k: dict(ratio=v["ratio"], ci95_ratio=v["ci95_ratio"]) for k, v in comps.items()}


# ----------------------------------------------------------------------------------------------
# (d) the primary result over the OPM white-noise sweep
def crossing(grid, med):
    """White level at which a decreasing median log2 ratio crosses 0 (log-linear interpolation);
    a bound when the crossing lies outside the grid."""
    below = np.flatnonzero(med < 0)
    if len(below) == 0:
        return dict(value=None, bound=f">{grid[-1]:g}")
    i = below[0]
    if i == 0:
        return dict(value=None, bound=f"<{grid[0]:g}")
    x = np.interp(0.0, [med[i], med[i - 1]], [np.log(grid[i]), np.log(grid[i - 1])])
    return dict(value=float(np.exp(x)), bound=None)


def part_sweep(ctx):
    st, cfg = ctx.st, ctx.cfg
    levels = [float(x) for x in cfg["sensors"]["opm_asd_fT_per_rtHz"]]
    stored = json.loads((ROOT / "results" / "g2" / "g2_summary.json").read_text())
    squid_d = detect(ctx.topo("squid"), ctx.noise_nom["squid"], ctx.arrays["squid"], ctx.conds)
    out = dict(levels_fT=levels, entries={}, stored_g2_commit=stored["provenance"]["commit"])
    max_diff = 0.0
    for asd in levels:
        d = {"squid": squid_d}
        for a in ctx.opms:
            d[a] = detect(ctx.topo(a), st.noise(ctx.arrays[a], ctx.G_g[a], ctx.brain_scale, ctx.env, asd=asd * 1e-15), ctx.arrays[a], ctx.conds)
        comps = compare_all(ctx, d)
        ref = stored["sensitivity"][f"opm_asd_{asd:g}fT"]
        for k, v in comps.items():
            v["stored_g2"] = dict(ratio=float(2 ** ref[k]["median_log2"]), ci95_ratio=[float(2 ** x) for x in ref[k]["ci95"]],
                                  n_boot=200)
            if len(ctx.tsel) == st.nt:
                max_diff = max(max_diff, abs(ref[k]["median_log2"] - v["median_log2"]))
        out["entries"][f"{asd:g}"] = dict(comparisons=comps, depth=profiles(ctx, d))
        log(f"sweep {asd:g} fT: dense {comps['opm_dense/combined/intrinsic+brain']['ratio']:.3f} | "
            f"{comps['opm_dense/combined/projected']['ratio']:.3f}, matched {comps['opm_matched/combined/intrinsic+brain']['ratio']:.3f} | "
            f"{comps['opm_matched/combined/projected']['ratio']:.3f}")
    out["reproduces_stored_g2_max_abs_diff_median_log2"] = max_diff if len(ctx.tsel) == st.nt else None
    # range over the sweep: ratios at the noisiest and quietest level, and the CI envelope
    lo_k, hi_k = f"{max(levels):g}", f"{min(levels):g}"
    out["range"] = {k: dict(ratio_at_max_noise=out["entries"][lo_k]["comparisons"][k]["ratio"],
                            ratio_at_min_noise=out["entries"][hi_k]["comparisons"][k]["ratio"],
                            ratio_at_primary=out["entries"]["15"]["comparisons"][k]["ratio"] if "15" in out["entries"] else None,
                            ci95_envelope=[min(out["entries"][x]["comparisons"][k]["ci95_ratio"][0] for x in out["entries"]),
                                           max(out["entries"][x]["comparisons"][k]["ci95_ratio"][1] for x in out["entries"])],
                            levels_with_ci_above_1=[float(x) for x in out["entries"] if out["entries"][x]["comparisons"][k]["ci95_ratio"][0] > 1],
                            levels_with_ci_below_1=[float(x) for x in out["entries"] if out["entries"][x]["comparisons"][k]["ci95_ratio"][1] < 1])
                    for k in out["entries"][lo_k]["comparisons"]}
    # break-even white level (median ratio = 1) on a fine grid, with its parcel-bootstrap CI
    sc = ctx.scfg["sweep"]
    grid = np.geomspace(sc["grid_fT_per_rtHz"][0], sc["grid_fT_per_rtHz"][1], int(sc["grid_points"]))
    groups = np.asarray(ctx.groups)
    labels = np.unique(groups)
    members = [np.flatnonzero(groups == lab) for lab in labels]
    be = {}
    for a in ctx.opms:
        dd = [detect(ctx.topo(a), st.noise(ctx.arrays[a], ctx.G_g[a], ctx.brain_scale, ctx.env, asd=g * 1e-15), ctx.arrays[a], ctx.conds)
              for g in grid]
        for cond in ctx.conds:
            D = np.array([x[("opm", cond)] for x in dd])
            for ref in ctx.refs:
                lr = np.log2(D / squid_d[(ref, cond)][None, :])
                med = np.median(lr, axis=1)
                point = crossing(grid, med)
                vals, n_out = [], 0
                for _ in range(ctx.n_boot):
                    idx = np.concatenate([members[j] for j in ctx.rng.integers(0, len(members), len(members))])
                    c = crossing(grid, np.median(lr[:, idx], axis=1))
                    if c["value"] is None:
                        n_out += 1
                        vals.append(grid[0] if c["bound"].startswith("<") else grid[-1])
                    else:
                        vals.append(c["value"])
                be[f"{a}/{ref}/{cond}"] = dict(break_even_fT=point["value"], bound=point["bound"],
                                               ci95_fT=[float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))],
                                               resamples_outside_grid=n_out, median_ratio_on_grid=[float(2 ** x) for x in med])
    out["break_even"] = dict(grid_fT=[float(x) for x in grid], entries=be,
                             method="median log2 ratio on a log-spaced grid of OPM white levels, crossing by log-linear "
                                    "interpolation; CI from the same parcel resamples applied to the whole curve (a resample "
                                    "whose curve does not cross within the grid counts at the grid end)")
    return out


# ----------------------------------------------------------------------------------------------
# (a) near-skull cortex
def oct6_targets(st):
    """Global indices of the valid oct-6 vertices, and which of them are the 7,661 targets."""
    n_lh = int(np.sum(st.cortex.hemi == 0))
    oct_global = np.concatenate([st.subject.src[0]["vertno"], st.subject.src[1]["vertno"] + n_lh])
    return oct_global[st.cortex.valid[oct_global]]


def region_of(st, vertices):
    names = np.concatenate([plotting.read_freesurfer_annot(st.subject.labels / f"{h}.aparc.annot") for h in ("lh", "rh")])
    return np.array([f"{'lh' if st.cortex.hemi[v] == 0 else 'rh'}.{names[v]}" for v in vertices], dtype=object)


def near_skull_vertices(ctx):
    """Every valid vertex closer than 4 mm to the inner-skull mesh (A-BEM-DIST) and the vertices used
    with their areas: all of them, or in test runs a random fraction with the areas scaled up so
    that the represented area stays unbiased."""
    cx = ctx.st.cortex
    near_all = np.flatnonzero(cx.valid & (cx.dist_inner_skull < anatomy.MIN_BEM_DISTANCE))
    frac = ctx.args.near_fraction
    if frac < 1.0:
        near_v = np.sort(np.random.default_rng(ctx.args.subset_seed).choice(near_all, int(round(frac * len(near_all))), replace=False))
    else:
        near_v = near_all
    return near_all, near_v, cx.area[near_v] / (len(near_v) / len(near_all))


def fullres_columns(ctx, name, vertices, chunk=20000):
    full, col = g2.fullres_matrix(ctx.arrays[name], ctx.st.subject, ctx.st.cortex, g2.FULLRES_JOBS[name])
    if np.any(col[vertices] < 0):
        raise ValueError("vertices outside the valid full-resolution set")
    out = np.empty((full.shape[0], len(vertices)))
    for s in range(0, len(vertices), chunk):
        out[:, s:s + chunk] = np.asarray(full[:, col[vertices[s:s + chunk]]], np.float64)
    return out


def bem_gains(ctx, name, vertices, bem, chunk=20000):
    """Lead fields of cortical-normal dipoles at ``vertices`` for one BEM (cached in cache/fwd in
    a full run only)."""
    a = ctx.arrays[name]
    cx = ctx.st.cortex
    out = np.empty((a.n, len(vertices)))
    for s in range(0, len(vertices), chunk):
        v = vertices[s:s + chunk]
        out[:, s:s + chunk] = forward.discrete_gain(a.info, ctx.st.subject.trans, cx.rr[v], cx.nn[v], bem, coil_def=opm.coil_def_file(),
                                                    use_cache=ctx.args.fwd_cache)[0]
    return out


def neighbour_energy(ctx, energy, vertices, bands):
    """Lead-field energy of each vertex over the median of its mesh neighbours', for vertices whose
    valid neighbours all have known gains, summarised by distance band from the inner skull."""
    cx = ctx.st.cortex
    pos = {int(v): i for i, v in enumerate(vertices)}
    adj = cx.adjacency.tocsr()
    d = cx.dist_inner_skull[vertices] * 1e3
    ratio = np.full(len(vertices), np.nan)
    for i, v in enumerate(vertices):
        nbv = [int(u) for u in adj.indices[adj.indptr[v]:adj.indptr[v + 1]] if cx.valid[u]]
        if len(nbv) >= 2 and all(u in pos for u in nbv):
            ratio[i] = energy[i] / np.median(energy[[pos[u] for u in nbv]])
    out = []
    flag = ctx.scfg["near_skull"]["neighbour_energy_flag"]
    for lo, hi in zip(bands[:-1], bands[1:]):
        m = (d >= lo) & (d < hi) & np.isfinite(ratio)
        r = ratio[m]
        out.append(dict(lo_mm=lo, hi_mm=hi, n=int(m.sum()), median=float(np.median(r)) if m.any() else None,
                        p99=float(np.percentile(r, 99)) if m.any() else None, max=float(r.max()) if m.any() else None,
                        n_above_flag=int(np.sum(r > flag))))
    return out


def noise_normalised_energy(ctx, name, gains):
    w = 1.0 / np.sqrt(ctx.noise_nom[name].intrinsic_var)
    return np.sum((gains * w[:, None]) ** 2, axis=0)


def near_skull_variants(ctx, gains, unit_grid, near_v, near_area, near_t, floors, label):
    """Background / target / both variants for one forward model. ``gains[name]`` holds the lead
    fields of [G2 targets (selected), near-skull vertices (near_v), near-skull targets (near_t)] in
    that order; ``unit_grid[name]`` the unit G2 grid background of the same model."""
    st, cx = ctx.st, ctx.st.cortex
    nt, nv = len(ctx.tsel), len(near_v)
    d_near = cx.dist_inner_skull[near_v] * 1e3
    d_tgt = cx.dist_inner_skull[near_t] * 1e3
    regions_t = region_of(st, near_t)
    depth_t = anatomy.depth_to_surface(cx.rr[near_t], st.subject.scalp) * 1e3
    out = dict(model=label)

    def evaluate(unit, keep_t):
        bs = background.calibrate(unit["squid"], ctx.grads, ctx.target_var)
        d = {}
        for name, a in ctx.arrays.items():
            nz0 = ctx.noise_nom[name]  # intrinsic noise and room field of the adult model; background replaced
            nz = noisemodel.ArrayNoise(nz0.intrinsic_var, bs * unit[name], nz0.env_cov, nz0.ext_basis)
            topo = np.concatenate([gains[name][:, :nt], gains[name][:, nt + nv:][:, keep_t]], axis=1) * st.q
            d[name] = detect(topo, nz, a, ctx.conds)
        return bs, d

    base_bs, base = evaluate(unit_grid, np.zeros(len(near_t), bool))
    out["baseline"] = dict(brain_scale=base_bs, comparisons=compare_all(ctx, base))
    for fl in floors:
        keep_v = d_near >= fl
        keep_t = d_tgt >= fl
        unit_bg = {name: unit_grid[name] + (gains[name][:, nt:nt + nv][:, keep_v] * near_area[keep_v][None, :]) @ gains[name][:, nt:nt + nv][:, keep_v].T
                   for name in ctx.arrays}
        groups_all = np.concatenate([ctx.groups, regions_t[keep_t]])
        depth_all = np.concatenate([ctx.depth, depth_t[keep_t]])
        is_near = np.r_[np.zeros(nt, bool), np.ones(int(keep_t.sum()), bool)]
        row = dict(floor_mm=fl, near_vertices=int(keep_v.sum()), near_area_cm2=float(near_area[keep_v].sum() * 1e4),
                   near_targets=int(keep_t.sum()))
        bs_bg, d_bg = evaluate(unit_bg, np.zeros(len(near_t), bool))
        row["background"] = dict(brain_scale_over_baseline=bs_bg / base_bs, comparisons=compare_all(ctx, d_bg),
                                 brain_rms_fT={name: float(np.sqrt(np.median(np.diag(bs_bg * unit_bg[name])[ctx.arrays[name].kinds == "mag"])) * 1e15)
                                               for name in ctx.arrays},
                                 brain_rms_baseline_fT={name: float(np.sqrt(np.median(np.diag(base_bs * unit_grid[name])[ctx.arrays[name].kinds == "mag"])) * 1e15)
                                                        for name in ctx.arrays})
        _, d_tg = evaluate(unit_grid, keep_t)
        row["targets"] = dict(all=compare_all(ctx, d_tg, groups=groups_all),
                              near_only=compare_all(ctx, {n: {k: v[is_near] for k, v in x.items()} for n, x in d_tg.items()},
                                                    groups=groups_all[is_near]) if is_near.sum() >= 20 else None)
        _, d_both = evaluate(unit_bg, keep_t)
        row["both"] = dict(all=compare_all(ctx, d_both, groups=groups_all),
                           near_only=compare_all(ctx, {n: {k: v[is_near] for k, v in x.items()} for n, x in d_both.items()},
                                                 groups=groups_all[is_near]) if is_near.sum() >= 20 else None)
        if label == "primary 3-layer BEM":
            row["depth"] = dict(baseline=profiles(ctx, base), background=profiles(ctx, d_bg),
                                both=profiles(ctx, d_both, depth=depth_all, groups=groups_all))
        out[f"floor_{fl:g}mm"] = row
        k = "opm_dense/combined/intrinsic+brain"
        log(f"near-skull, {label}, floor {fl:g} mm: dense {out['baseline']['comparisons'][k]['ratio']:.3f} -> background "
            f"{row['background']['comparisons'][k]['ratio']:.3f}, targets {row['targets']['all'][k]['ratio']:.3f}, both {row['both']['all'][k]['ratio']:.3f}")
    return out


def part_near_skull(ctx):
    st, cx, ns = ctx.st, ctx.st.cortex, ctx.scfg["near_skull"]
    rule = ns["rule_mm"] * 1e-3
    if abs(rule - anatomy.MIN_BEM_DISTANCE) > 1e-12:
        raise ValueError("the near-skull rule must be the code base's A-BEM-DIST")
    near_all, near_v, near_area = near_skull_vertices(ctx)
    oct_valid = oct6_targets(st)
    near_t = oct_valid[cx.dist_inner_skull[oct_valid] < rule]
    if ctx.args.max_targets and len(ctx.tsel) < st.nt:
        near_t = np.sort(np.random.default_rng(ctx.args.subset_seed + 1).choice(near_t, max(20, int(len(near_t) * len(ctx.tsel) / st.nt)),
                                                                                replace=False))
    tg = st.src.target[ctx.tsel]
    out = dict(rule="valid white-surface vertices closer than 4 mm to the 5,120-triangle inner-skull mesh (A-BEM-DIST)",
               n_near_vertices=int(len(near_all)), share_of_valid_vertices=float(len(near_all) / cx.valid.sum()),
               near_area_cm2=float(cx.area[near_all].sum() * 1e4), valid_area_cm2=float(cx.area[cx.valid].sum() * 1e4),
               n_near_vertices_used=int(len(near_v)), n_near_oct6_targets=int(np.sum(cx.dist_inner_skull[oct_valid] < rule)),
               n_near_targets_used=int(len(near_t)),
               vertices_by_band=[dict(lo_mm=lo, hi_mm=hi, n=int(np.sum((cx.dist_inner_skull[near_all] * 1e3 >= lo) & (cx.dist_inner_skull[near_all] * 1e3 < hi))),
                                      area_cm2=float(cx.area[near_all][(cx.dist_inner_skull[near_all] * 1e3 >= lo) & (cx.dist_inner_skull[near_all] * 1e3 < hi)].sum() * 1e4))
                                 for lo, hi in zip(ns["diagnostic_bands_mm"][:-1], ns["diagnostic_bands_mm"][1:])],
               background_discretisation="the G2 7-mm grid over the usable cortex (Voronoi areas) plus every near-skull vertex with its "
                                         "own area: an exact partition of the cortex used, moment variance proportional to area")
    bands = list(ns["diagnostic_bands_mm"])
    # primary 3-layer model: the cached full-resolution lead fields (every valid vertex)
    gains = {name: np.concatenate([ctx.G_t[name], fullres_columns(ctx, name, near_v), fullres_columns(ctx, name, near_t)], axis=1)
             for name in ctx.arrays}
    # artefact screen on the cached 3-layer fields: every valid vertex below 6 mm, so that each has its neighbours
    screen_v = np.flatnonzero(cx.valid & (cx.dist_inner_skull < 0.006))
    out["primary_neighbour_energy"] = {name: neighbour_energy(ctx, noise_normalised_energy(ctx, name, fullres_columns(ctx, name, screen_v)),
                                                              screen_v, bands + [5.0, 6.0]) for name in ctx.arrays}
    log("near-skull: 3-layer lead fields and artefact screen")
    out["primary"] = near_skull_variants(ctx, gains, ctx.unit, near_v, near_area, near_t, ns["floor_primary_mm"], "primary 3-layer BEM")
    del gains
    # accuracy check: 1-layer BEM, inner skull 5,120 vs 20,480 triangles, everything recomputed in each
    pts = np.concatenate([tg, near_v, near_t, st.src.grid])
    uniq, inv = np.unique(pts, return_inverse=True)  # the near-skull targets are also near-skull vertices
    models = {"coarse": (st.subject.bem_model((0.3,)), ns["floor_coarse_check_mm"]),
              "refined": (anatomy.refined_inner_skull(st.subject, 1), ns["floor_refined_check_mm"])}
    one = {}
    for label, (bem, floors) in models.items():
        g = {name: bem_gains(ctx, name, uniq, bem)[:, inv] for name in ctx.arrays}
        log(f"near-skull: 1-layer {label} lead fields")
        n1 = len(tg) + len(near_v) + len(near_t)
        unit_grid = {name: background.sensor_covariance(g[name][:, n1:], background.moment_covariance(st.src.grid_area)) for name in ctx.arrays}
        one[label] = dict(gains={name: g[name][:, :n1] for name in ctx.arrays})
        out[f"one_layer_{label}"] = near_skull_variants(ctx, one[label]["gains"], unit_grid, near_v, near_area, near_t, floors,
                                                        f"1-layer BEM, inner skull {'20,480' if label == 'refined' else '5,120'} triangles")
        sl = slice(len(tg), len(tg) + len(near_v))
        out[f"one_layer_{label}"]["neighbour_energy"] = {name: neighbour_energy(ctx, noise_normalised_energy(ctx, name, g[name][:, sl]), near_v, bands)
                                                         for name in ctx.arrays}
        del g
    # per-vertex gain change coarse -> refined, by distance from the inner skull (near-skull vertices and the targets)
    acc = {}
    dist = np.concatenate([cx.dist_inner_skull[tg], cx.dist_inner_skull[near_v]]) * 1e3
    for name in ctx.arrays:
        gc = np.concatenate([one["coarse"]["gains"][name][:, :len(tg)], one["coarse"]["gains"][name][:, len(tg):len(tg) + len(near_v)]], axis=1)
        gr = np.concatenate([one["refined"]["gains"][name][:, :len(tg)], one["refined"]["gains"][name][:, len(tg):len(tg) + len(near_v)]], axis=1)
        w = 1.0 / np.sqrt(ctx.noise_nom[name].intrinsic_var)[:, None]
        rel = np.linalg.norm(w * (gc - gr), axis=0) / np.linalg.norm(w * gr, axis=0)
        rows = []
        for lo, hi in zip(bands[:-1] + [4.0, 5.0], bands[1:] + [5.0, 6.0]):
            m = (dist >= lo) & (dist < hi)
            rows.append(dict(lo_mm=lo, hi_mm=hi, n=int(m.sum()), median=float(np.median(rel[m])) if m.any() else None,
                             p90=float(np.percentile(rel[m], 90)) if m.any() else None, max=float(rel[m].max()) if m.any() else None))
        acc[name] = rows
    out["gain_change_coarse_to_refined"] = acc
    # the inclusion effect (log2 change of the median ratio) in each model, the check of the primary model's result
    eff = {}
    for model_key in ("primary", "one_layer_coarse", "one_layer_refined"):
        m = out[model_key]
        for fk in [k for k in m if k.startswith("floor_")]:
            for var in ("background", "targets", "both"):
                comps = m[fk][var]["comparisons"] if var == "background" else m[fk][var]["all"]
                for k, v in comps.items():
                    eff.setdefault(k, {})[f"{model_key}/{fk}/{var}"] = float(v["median_log2"] - m["baseline"]["comparisons"][k]["median_log2"])
    out["inclusion_effect_log2"] = eff
    return out


# ----------------------------------------------------------------------------------------------
# (b) coloured OPM noise: sub-band partition of the analysis filter
class SubbandFilter:
    """Zero-phase FFT filter whose composite power response is the primary filter's, |H(f)|^4,
    restricted to [lo, hi): the sub-bands partition the primary band exactly (their ENBWs and the
    variances of stationary noise add up to the primary's). Same interface as noise.AnalysisFilter
    (apply, power_response, enbw, fs) so that the code base's measurement and room-field fits use it."""

    def __init__(self, base, lo: float, hi: float, pad_s: float = 30.0):
        self.base, self.lo, self.hi, self.fs, self.pad_s = base, float(lo), float(hi), base.fs, pad_s

    def power_response(self, f):
        f = np.asarray(f, float)
        return self.base.power_response(f) * ((f >= self.lo) & (f < self.hi))

    def enbw(self, n_freqs: int = 400_001) -> float:
        f = np.linspace(0.0, self.fs / 2.0, n_freqs)
        return float(np.trapezoid(self.power_response(f), f))

    def amplitude(self, nfft):
        f = np.fft.rfftfreq(nfft, 1.0 / self.fs)
        return np.sqrt(self.power_response(f))

    def apply(self, x, axis: int = -1):
        from scipy import fft as sfft

        x = np.moveaxis(np.asarray(x, float), axis, -1)
        n = x.shape[-1]
        nfft = sfft.next_fast_len(n + 2 * int(self.pad_s * self.fs), real=True)
        y = sfft.irfft(sfft.rfft(x, n=nfft, axis=-1) * self.amplitude(nfft), n=nfft, axis=-1)[..., :n]
        return np.moveaxis(y, -1, axis)


@dataclass
class FixedProjectorNoise(noisemodel.ArrayNoise):
    """ArrayNoise whose external-field projector keeps given weights (1 / proj_var) whatever the
    intrinsic variance of the sub-band: one fixed operator applied to the data, as in practice. The
    projector formula is noisemodel.ArrayNoise.projector's with these weights."""

    proj_var: np.ndarray | None = None

    def projector(self) -> np.ndarray:
        e = self.ext_basis
        w = 1.0 / np.asarray(self.proj_var, float)
        m = e.T @ (w[:, None] * e)
        return np.eye(len(w)) - e @ np.linalg.solve(m, e.T * w[None, :])


class Recordings:
    """The sample task recording and its empty-room recording, read once, as g2.measured_noise reads
    them (MEG channels in the SQUID order, mean removed, filtered as continuous records, 2-s trims,
    the -200-0 ms pre-stimulus windows of the task recording; bad channels NaN). Copied from
    opmsquid.g2.measured_noise so that many sub-band filters reuse one FFT of each record; the
    room-field fit is opmsquid.environment.fit_empty_room itself."""

    def __init__(self, squid_info, bads, window=(-0.2, 0.0), trim_s=2.0, pad_s=30.0):
        from scipy import fft as sfft

        names = squid_info.ch_names
        self.squid_info, self.bads, self.trim_s = squid_info, list(bads), trim_s
        self.bad = np.isin(names, list(bads))
        self.spec, self.n, self.nfft = {}, {}, {}
        for key, fname in (("empty_room", g2.ER_FILE), ("baseline", neuromag.RAW_FILE)):
            raw = mne.io.read_raw_fif(paths.SAMPLE_MEG / fname, preload=True, verbose=False)
            data = raw.get_data(picks=[raw.ch_names.index(n) for n in names])
            data = data - data.mean(axis=1, keepdims=True)
            sf = raw.info["sfreq"]
            self.n[key] = data.shape[1]
            self.nfft[key] = sfft.next_fast_len(data.shape[1] + 2 * int(pad_s * sf), real=True)
            self.spec[key] = sfft.rfft(data, n=self.nfft[key], axis=1)
            self.sfreq = sf
            if key == "baseline":
                ev = mne.find_events(raw, stim_channel="STI 014", verbose=False)[:, 0] - raw.first_samp
                self.win = (int(round(window[0] * sf)), int(round(window[1] * sf)))
                self.events = ev
            else:
                self.er_raw = raw  # environment.fit_empty_room picks the MEG channels itself, as in g2.measured_noise
            del data

    def filtered(self, key, filt):
        from scipy import fft as sfft

        return sfft.irfft(self.spec[key] * filt.amplitude(self.nfft[key])[None, :], n=self.nfft[key], axis=1)[:, :self.n[key]]

    def measured(self, filt) -> dict:
        out = {}
        n_trim = int(self.trim_s * self.sfreq)
        for key in ("empty_room", "baseline"):
            data = self.filtered(key, filt)
            if key == "baseline":
                a, b = self.win
                segs = [data[:, e + a:e + b] for e in self.events if e + a > n_trim and e + b < data.shape[1] - n_trim]
                seg = np.concatenate(segs, axis=1)
                out["n_windows"] = len(segs)
            else:
                seg = data[:, n_trim:-n_trim]
            out[key] = np.mean(seg ** 2, axis=1)
            out[key][self.bad] = np.nan
        out["brain"] = out["baseline"] - out["empty_room"]
        out["environment"] = environment.fit_empty_room(self.er_raw, self.squid_info, filt, bads=self.bads)
        return out


def opm_shape(filt, corner, exponent, n_freqs=400_001):
    """Integral of (1 + (f_c/f)^exponent) over the filter's power response, divided by its ENBW
    (the factor by which the 1/f part raises the white variance in that (sub-)band)."""
    f = np.linspace(1e-4, filt.fs / 2.0, n_freqs)
    p = filt.power_response(f)
    return float(np.trapezoid(p * (1.0 + (corner / f) ** exponent), f) / np.trapezoid(p, f)) if np.trapezoid(p, f) > 0 else 1.0


def signal_weights(subs, base, kind):
    """Share of the signal energy (after the analysis filter) in each sub-band: 'flat' (a flat
    amplitude spectrum: the adult convention generalised) or 'spike-wave' (opmsquid.ied.ied_waveform)."""
    enbw = np.array([s["enbw"] for s in subs])
    if kind == "flat":
        return enbw / enbw.sum()
    w = ied.ied_waveform(base.fs)
    nfft = 1 << 18
    spec = np.abs(np.fft.rfft(w, n=nfft)) ** 2
    f = np.fft.rfftfreq(nfft, 1.0 / base.fs)
    e = np.array([np.sum(spec * s["filt"].power_response(f)) for s in subs])
    return e / e.sum()


def make_subbands(ctx, rec, edges):
    fs = ctx.st.filt.fs
    edges = list(edges) + [fs / 2.0 + 1e-9]
    subs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        f = SubbandFilter(ctx.st.filt, lo, hi)
        m = rec.measured(f)
        subs.append(dict(lo=lo, hi=min(hi, fs / 2.0), filt=f, enbw=f.enbw(), target_grad=float(np.nanmedian(m["brain"][ctx.grads])),
                         meas_mag=float(np.nanmedian(m["brain"][ctx.mags])), env=m["environment"],
                         brain_rms_grad_fT_cm=float(np.sqrt(max(np.nanmedian(m["brain"][ctx.grads]), 0.0)) * 1e13),
                         brain_rms_mag_fT=float(np.sqrt(max(np.nanmedian(m["brain"][ctx.mags]), 0.0)) * 1e15)))
    return subs


def subband_noises(ctx, subs, unit, ff=None, ff_level=None):
    """Per sub-band: brain scale calibrated on the measured gradiometer brain noise of that
    sub-band (and the far-field covariance, if any, filling the magnetometer shortfall or a fixed
    level), room field refitted in that sub-band. Returns per sub-band {name: ArrayNoise without
    the OPM intrinsic term} and the scales."""
    rows = []
    for s in subs:
        cal = None
        if ff is None:
            scale = background.calibrate(unit["squid"], ctx.grads, s["target_grad"]) if s["target_grad"] > 0 else 0.0
            level2 = 0.0
        else:
            lev = ff_level
            if ff_level == "room_field":  # the room field's RMS on the Neuromag magnetometers in this sub-band
                lev = float(np.sqrt(np.median(np.diag(s["env"].covariance(ctx.basis["squid"]))[ctx.mags])))
            scale, level2, cal = calibrate_with_far_field(ctx, unit["squid"], ff["squid"], s["target_grad"], s["meas_mag"], lev, "added")
        nzs = {}
        for name, a in ctx.arrays.items():
            brain = scale * unit[name] + (level2 * ff[name] if ff is not None else 0.0)
            ivar = np.array([g2.SQUID_ASD[k] ** 2 for k in a.kinds]) * s["enbw"] if name == "squid" else np.full(a.n, np.nan)
            nzs[name] = noisemodel.ArrayNoise(ivar, brain, s["env"].covariance(ctx.basis[name]), ctx.basis[name])
        rows.append(dict(scale=scale, far_field_level2=level2, noises=nzs, calibration=cal))
    return rows


def fr_detect(ctx, subs, rows, weights, white, corner, exponent, names, conds, topo=None):
    """Frequency-resolved detectability d^2 = sum_k u_k (ENBW_k / ENBW) s^T C_k^+ s (u_k: the signal's
    energy share in sub-band k), equal to the band detectability when every noise term is white."""
    enbw = sum(s["enbw"] for s in subs)
    out = {}
    for name in names:
        a = ctx.arrays[name]
        acc = {}
        for s, r, u in zip(subs, rows, weights):
            nz = r["noises"][name]
            if name != "squid":
                nz = noisemodel.ArrayNoise(np.full(a.n, (white * 1e-15) ** 2 * s["enbw"] * opm_shape_cached(s, corner, exponent)),
                                           nz.brain_cov, nz.env_cov, nz.ext_basis)
            dk = detect(ctx.topo(name) if topo is None else topo[name], nz, a, conds)
            for k, v in dk.items():
                acc[k] = acc.get(k, 0.0) + u * (s["enbw"] / enbw) * v ** 2
        out[name] = {k: np.sqrt(v) for k, v in acc.items()}
    return out


_SHAPE: dict = {}


def opm_shape_cached(s, corner, exponent):
    key = (s["lo"], s["hi"], corner, exponent)
    if key not in _SHAPE:
        _SHAPE[key] = opm_shape(s["filt"], corner, exponent) if corner > 0 else 1.0
    return _SHAPE[key]


def part_coloured(ctx, rec):
    st, cn = ctx.st, ctx.scfg["coloured_noise"]
    whites = [float(x) for x in ctx.cfg["sensors"]["opm_asd_fT_per_rtHz"]]
    primary_w = float(ctx.cfg["sensors"]["opm_asd_primary_fT_per_rtHz"])
    corners = [0.0] + [float(c) for c in cn["corner_hz"]]
    expo = float(cn["exponent"])
    subs = make_subbands(ctx, rec, cn["subband_edges_hz"])
    ctx.extra["subbands"] = subs
    rows = subband_noises(ctx, subs, ctx.unit)
    ctx.extra["subband_rows"] = rows
    weights = {kind: signal_weights(subs, st.filt, kind) for kind in cn["signal_spectra"]}
    out = dict(model=f"OPM noise PSD = w^2 (1 + (f_c / f)^{expo:g}); Neuromag white (brochure); brain and room field from the recordings per sub-band",
               corners_hz=corners[1:], white_levels_fT=whites,
               subbands=[dict(lo_hz=s["lo"], hi_hz=s["hi"], enbw_hz=s["enbw"], brain_rms_grad_fT_cm=s["brain_rms_grad_fT_cm"],
                              brain_rms_mag_fT=s["brain_rms_mag_fT"], brain_scale=r["scale"],
                              env_explained=s["env"].explained_fraction, signal_share={k: float(w[i]) for k, w in weights.items()},
                              opm_variance_factor={f"{c:g}": opm_shape_cached(s, c, expo) for c in corners})
                         for i, (s, r) in enumerate(zip(subs, rows))],
               checks=dict(enbw_sum_over_primary=sum(s["enbw"] for s in subs) / st.enbw,
                           brain_grad_variance_sum_over_primary=sum(s["target_grad"] for s in subs) / ctx.target_var,
                           primary_enbw_hz=st.enbw),
               band_variance_factor={f"{c:g}": opm_shape(st.filt, c, expo) if c > 0 else 1.0 for c in corners},
               method=("Frequency-resolved whitening: the 1-40 Hz composite filter response is partitioned into the listed sub-bands; "
                       "in each, every noise term has its own level (OPM: w^2 times the 1/f + white shape integrated over the "
                       "sub-band; Neuromag: white; brain: the model's spatial covariance scaled to the measured gradiometer "
                       "brain-noise variance of that sub-band; room field: the 8-term fit to the empty-room recording in that "
                       "sub-band), and the matched-filter SNR adds over sub-bands: d^2 = sum_k u_k (ENBW_k/ENBW) s^T C_k^+ s with u_k "
                       "the signal's share of energy in sub-band k. Noise in disjoint frequency bands is uncorrelated, so this is the "
                       "optimal (spatio-spectral) detector's SNR when the noise spectra are flat within each sub-band; it equals the "
                       "band detectability of the adult comparison when every noise term is white. The band-variance alternative "
                       "(spatial whitening only) keeps the adult band covariance and raises the OPM variance by the band-averaged "
                       "1/f factor: the penalty a detector that ignores the spectral shapes would pay (an upper bound on the "
                       "penalty). A uniformly worse white sensor (the white-level sweep) is the contrast it is compared with."))
    res = {}
    # Neuromag does not depend on the OPM noise: one evaluation per signal spectrum
    for kind, w in weights.items():
        sq = fr_detect(ctx, subs, rows, w, 0.0, 0.0, expo, ["squid"], ctx.conds)["squid"]
        for white in whites:
            for c in corners:
                d = {"squid": sq}
                d.update(fr_detect(ctx, subs, rows, w, white, c, expo, ctx.opms, ctx.conds))
                full = white == primary_w
                comps = compare_all(ctx, d, refs=None if full else ["combined"])
                res[f"frequency_resolved/{kind}/{white:g}fT/{c:g}Hz"] = dict(comparisons=comps, depth=profiles(ctx, d) if (full and kind == "flat") else None)
        log(f"coloured, frequency-resolved ({kind}): done")
    # band-variance alternative (spatial whitening only): the adult band model, OPM variance x band 1/f factor
    sq = detect(ctx.topo("squid"), ctx.noise_nom["squid"], ctx.arrays["squid"], ctx.conds)
    for white in whites:
        for c in corners:
            asd = white * 1e-15 * np.sqrt(out["band_variance_factor"][f"{c:g}"])
            d = {"squid": sq}
            for a in ctx.opms:
                d[a] = detect(ctx.topo(a), st.noise(ctx.arrays[a], ctx.G_g[a], ctx.brain_scale, ctx.env, asd=asd), ctx.arrays[a], ctx.conds)
            full = white == primary_w
            res[f"band_variance/{white:g}fT/{c:g}Hz"] = dict(comparisons=compare_all(ctx, d, refs=None if full else ["combined"]),
                                                            depth=profiles(ctx, d) if full else None, equivalent_white_fT=float(asd * 1e15))
    log("coloured, band-variance alternative: done")
    # counterpart on the Neuromag side: its sensor noise from the empty-room recording in every sub-band (what the 8-term
    # fit leaves: the adult sensitivity 'squid_measured_spectrum' resolved in frequency, an upper bound on its sensor
    # noise), the projector keeping the weights of the band-level measured variances as in that sensitivity
    kinds = ctx.arrays["squid"].kinds
    proj_var = g2.measured_squid_variance(ctx.env, kinds)
    rows_m = []
    for s, r in zip(subs, rows):
        nz = r["noises"]["squid"]
        sq_nz = FixedProjectorNoise(g2.measured_squid_variance(s["env"], kinds), nz.brain_cov, nz.env_cov, nz.ext_basis, proj_var)
        rows_m.append(dict(r, noises=dict(r["noises"], squid=sq_nz)))
    sq_m = fr_detect(ctx, subs, rows_m, weights["flat"], 0.0, 0.0, expo, ["squid"], ctx.conds)["squid"]
    for c in corners:
        d = {"squid": sq_m}
        d.update(fr_detect(ctx, subs, rows, weights["flat"], primary_w, c, expo, ctx.opms, ctx.conds))
        res[f"frequency_resolved_squid_measured/flat/{primary_w:g}fT/{c:g}Hz"] = dict(comparisons=compare_all(ctx, d), depth=None)
    out["squid_measured_rms"] = [dict(lo_hz=s["lo"], mag_fT=float(np.sqrt(np.median(g2.measured_squid_variance(s["env"], kinds)[kinds == "mag"])) * 1e15),
                                      brochure_mag_fT=float(np.sqrt(g2.SQUID_ASD["mag"] ** 2 * s["enbw"]) * 1e15)) for s in subs]
    log("coloured, measured Neuromag sensor noise counterpart: done")
    out["results"] = res
    # partition convergence: every sub-band split, corner 10 Hz and white, primary white level, flat spectrum
    split = int(cn.get("partition_check_split", 0))
    if split > 1:
        fine_edges = []
        for lo, hi in zip(list(cn["subband_edges_hz"]), list(cn["subband_edges_hz"][1:]) + [st.filt.fs / 2.0]):
            fine_edges += list(np.linspace(lo, hi, split + 1)[:-1])
        fsubs = make_subbands(ctx, rec, fine_edges)
        frows = subband_noises(ctx, fsubs, ctx.unit)
        fw = signal_weights(fsubs, st.filt, "flat")
        chk = {}
        sq_f = fr_detect(ctx, fsubs, frows, fw, 0.0, 0.0, expo, ["squid"], ctx.conds)["squid"]
        for c in (0.0, max(corners)):
            d = {"squid": sq_f}
            d.update(fr_detect(ctx, fsubs, frows, fw, primary_w, c, expo, ctx.opms, ctx.conds))
            fine = compare_all(ctx, d, refs=["combined"])
            coarse = res[f"frequency_resolved/flat/{primary_w:g}fT/{c:g}Hz"]["comparisons"]
            chk[f"{c:g}Hz"] = {k: dict(fine=v["ratio"], primary_partition=coarse[k]["ratio"],
                                       abs_diff_log2=abs(v["median_log2"] - coarse[k]["median_log2"])) for k, v in fine.items()}
        out["partition_check"] = dict(n_subbands=len(fsubs), split=split, results=chk)
        log("coloured, partition check: done")
    return out


# ----------------------------------------------------------------------------------------------
# (c) far-field physiological sources
def read_mgz(fname):
    """Minimal reader of a FreeSurfer .mgz volume (gzip-compressed MGH; the study has no nibabel):
    data as (i, j, k) and the header's voxel sizes."""
    with gzip.open(fname, "rb") as f:
        raw = f.read()
    version, w, h, d, nframes, dtype, _dof = (int(x) for x in np.frombuffer(raw[:28], ">i4"))
    dt = {0: ">u1", 1: ">i4", 3: ">f4", 4: ">i2"}[dtype]
    n = w * h * d * nframes
    data = np.frombuffer(raw[284:284 + n * np.dtype(dt).itemsize], dt).reshape((nframes, d, h, w))[0]
    spacing = np.frombuffer(raw[30:42], ">f4")
    return np.transpose(data, (2, 1, 0)), dict(version=version, shape=(w, h, d), spacing=spacing.tolist())


def locate_eyes(st, ff):
    """Centres of the two eyes of the sample subject: the darkest-centred spheres (vitreous, dark on
    T1) with a bright surround (orbital fat) in the declared head-frame search box, one per side.
    The conformed T1 (256^3, 1 mm) maps voxels to FreeSurfer surface RAS, the MRI frame of MNE."""
    from scipy.signal import fftconvolve

    t1, hdr = read_mgz(paths.SUBJECTS_DIR / "sample" / "mri" / "T1.mgz")
    if tuple(hdr["shape"]) != (256, 256, 256) or not np.allclose(hdr["spacing"], 1.0):
        raise ValueError("expected a conformed 256^3 1-mm T1")
    vox2tkr = np.array([[-1, 0, 0, 128], [0, 0, 1, -128], [0, -1, 0, 128], [0, 0, 0, 1]], float)
    tkr2vox = np.linalg.inv(vox2tkr)
    r_in, (r1, r2) = ff["eye_inner_radius_mm"], ff["eye_shell_mm"]
    rad = int(np.ceil(r2))
    z, y, x = np.mgrid[-rad:rad + 1, -rad:rad + 1, -rad:rad + 1]
    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    k_in = (dist <= r_in).astype(float)
    k_sh = ((dist >= r1) & (dist <= r2)).astype(float)
    head_mri = st.subject.trans["trans"]
    cand = {}
    for side, sgn in (("left", -1.0), ("right", 1.0)):
        xs = np.arange(ff["eye_search_x_mm"][0], ff["eye_search_x_mm"][1] + 0.5, 1.0) * sgn
        ys = np.arange(ff["eye_search_y_mm"][0], ff["eye_search_y_mm"][1] + 0.5, 1.0)
        zs = np.arange(ff["eye_search_z_mm"][0], ff["eye_search_z_mm"][1] + 0.5, 1.0)
        P = np.stack(np.meshgrid(xs, ys, zs, indexing="ij"), axis=-1).reshape(-1, 3) * 1e-3
        mri = (P @ head_mri[:3, :3].T + head_mri[:3, 3]) * 1e3
        cand[side] = (P, mri, np.round(mri @ tkr2vox[:3, :3].T + tkr2vox[:3, 3]).astype(int))
    allvox = np.concatenate([c[2] for c in cand.values()])
    lo, hi = allvox.min(axis=0) - rad - 1, allvox.max(axis=0) + rad + 2
    if np.any(lo < 0) or np.any(hi > 256):
        raise ValueError("eye search box outside the MRI volume")
    sub = t1[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]].astype(float)
    mean_in = fftconvolve(sub, k_in / k_in.sum(), mode="same")  # symmetric kernels: convolution = local mean
    mean_sh = fftconvolve(sub, k_sh / k_sh.sum(), mode="same")
    out = {}
    for side, (P, mri, vox) in cand.items():
        v = tuple((vox - lo).T)
        score = mean_sh[v] - mean_in[v]
        i = int(np.argmax(score))
        out[side] = dict(head_mm=(P[i] * 1e3).tolist(), mri_mm=mri[i].tolist(), inner_mean=float(mean_in[v][i]),
                         shell_mean=float(mean_sh[v][i]), score=float(score[i]))
    out["interocular_distance_mm"] = float(np.linalg.norm(np.subtract(out["left"]["head_mm"], out["right"]["head_mm"])))
    return out


def dipole_field_gain(array, r0, kind):
    """(n_ch, 3) response of every channel (its own coil integration points, as the room-field
    basis) to unit moments along head x, y, z at r0 [m, head frame]: 'current_dipole' = the
    primary-current field of a current dipole in an unbounded homogeneous conductor (1 A m),
    'magnetic_dipole' = a magnetic dipole (1 A m^2)."""
    out = np.zeros((array.n, 3))
    for i, c in enumerate(environment.coil_integration(array.info, array.coil_def)):
        R = c["rmag"] - np.asarray(r0, float)
        Rn = np.linalg.norm(R, axis=1)[:, None]
        for j, q in enumerate(np.eye(3)):
            if kind == "current_dipole":
                B = MU0_4PI * np.cross(q, R) / Rn ** 3
            elif kind == "magnetic_dipole":
                Rh = R / Rn
                B = MU0_4PI * (3.0 * (Rh @ q)[:, None] * Rh - q) / Rn ** 3
            else:
                raise ValueError(kind)
            out[i, j] = c["w"] @ np.sum(B * c["cosmag"], axis=1)
    return out


def eye_bem_gain(ctx, array, eyes_head):
    """(n_ch, 3) response to unit current dipoles along head x, y, z at both eye centres together
    (conjugate eye movements), in a homogeneous head: a 1-layer BEM on the head surface (0.3 S/m)."""
    from mne.surface import complete_surface_info

    st = ctx.st
    head = next(s for s in st.subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    surf = complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_BRAIN, sigma=0.3, coord_frame=head["coord_frame"], rr=head["rr"].copy(),
                                      tris=head["tris"].copy(), np=len(head["rr"]), ntri=len(head["tris"])), copy=False, verbose=False)
    head_mri = st.subject.trans["trans"]
    rr = np.array(eyes_head) @ head_mri[:3, :3].T + head_mri[:3, 3]
    out = np.zeros((array.n, 3))
    for j, q in enumerate(np.eye(3)):
        nn = np.repeat((q @ head_mri[:3, :3].T)[None], len(rr), axis=0)
        g, _ = forward.discrete_gain(array.info, st.subject.trans, rr, nn, [surf], coil_def=opm.coil_def_file(), use_cache=ctx.args.fwd_cache)
        out[:, j] = np.asarray(g, float).sum(axis=1)
    return out


def calibrate_with_far_field(ctx, unit_sq, ff_sq, target_grad, meas_mag, level, mode="added"):
    """Brain scale s and far-field variance factor a^2 for a far field at its level on the Neuromag
    magnetometers: a fixed RMS (in T) or 'magnetometer_shortfall' (median_mag(measured brain noise)
    - median_mag(s unit)). mode 'added': s stays the gradiometer calibration without far field (the
    published background) and the far field comes on top (adverse to the OPM; the model then
    over-predicts the measured gradiometer brain noise by the far field's share). mode
    'recalibrated': s is refitted so that median_grad(s unit + a^2 F) is the measured gradiometer
    brain noise (consistent with the data; infeasible when the far field alone exceeds it).
    Returns (s, a^2, dict(feasible, model/measured RMS ratios on gradiometers and magnetometers))."""
    from scipy.optimize import brentq

    du, df = np.diag(unit_sq), np.diag(ff_sq)
    f_mag = np.median(df[ctx.mags])
    if target_grad <= 0:
        return 0.0, 0.0, dict(feasible=False, reason="no measured gradiometer brain noise in this band")
    s0 = target_grad / np.median(du[ctx.grads])

    def fit(s, a2):
        return dict(grad_rms_model_over_measured=float(np.sqrt(np.median(s * du[ctx.grads] + a2 * df[ctx.grads]) / target_grad)),
                    mag_rms_model_over_measured=float(np.sqrt(np.median(s * du[ctx.mags] + a2 * df[ctx.mags]) / meas_mag)) if meas_mag > 0 else None,
                    scale_over_gradiometer_calibration=float(s / s0))

    if mode == "added":
        if level == "magnetometer_shortfall":
            a2 = max(meas_mag - np.median(s0 * du[ctx.mags]), 0.0) / f_mag
        else:
            a2 = float(level) ** 2 / f_mag
        return float(s0), float(a2), dict(feasible=True, **fit(s0, a2))
    if mode != "recalibrated":
        raise ValueError(mode)

    def scale_for(a2):
        # solved for x = s / s0: brentq's absolute tolerance (2e-12 by default) is far larger than the scale itself
        # (~3.5e-14), so the unnormalised root would come back as ~0
        g = lambda x: np.median(x * s0 * du[ctx.grads] + a2 * df[ctx.grads]) / target_grad - 1.0  # noqa: E731
        return (float(brentq(g, 0.0, 2.0, xtol=1e-13, rtol=1e-13) * s0), True) if g(0.0) < 0 else (0.0, False)

    if level == "magnetometer_shortfall":
        s, a2, ok = s0, 0.0, True
        for _ in range(100):
            a2_new = max(meas_mag - np.median(s * du[ctx.mags]), 0.0) / f_mag
            s_new, ok = scale_for(a2_new)
            done = abs(s_new - s) <= 1e-10 * max(s, 1e-300) and abs(a2_new - a2) <= 1e-10 * max(a2_new, 1e-300)
            s, a2 = s_new, a2_new
            if done or not ok:
                break
        else:
            raise RuntimeError("far-field calibration did not converge")
        expected = max(meas_mag - np.median(s * du[ctx.mags]), 0.0)
    else:
        a2 = float(level) ** 2 / f_mag
        s, ok = scale_for(a2)
        expected = float(level) ** 2
    info = dict(feasible=bool(ok and s > 0), **fit(s, a2))
    if not info["feasible"]:
        info["reason"] = "the far field alone exceeds the measured gradiometer brain noise (no cortical background left)"
        return float(s), float(a2), info
    # self-check: the gradiometer brain noise is matched and the far field has its level on the magnetometers
    if abs(info["grad_rms_model_over_measured"] - 1.0) > 1e-6 or abs(a2 * f_mag - expected) > 1e-6 * max(expected, meas_mag, 1e-300):
        raise RuntimeError(f"far-field calibration failed ({info}, level {a2 * f_mag:.3e} vs {expected:.3e})")
    return float(s), float(a2), info


def part_far_field(ctx, eyes=None):
    st, ff = ctx.st, ctx.scfg["far_field"]
    eyes = eyes or locate_eyes(st, ff)
    eyes_head = [np.array(eyes[s]["head_mm"]) * 1e-3 for s in ("left", "right")]
    log(f"eyes (head frame, mm): left {np.round(eyes['left']['head_mm'], 1)}, right {np.round(eyes['right']['head_mm'], 1)}, "
        f"interocular {eyes['interocular_distance_mm']:.1f} mm")
    sources = {}
    xy = np.array(ff["heart_head_mm"], float) * 1e-3
    for depth in ff["heart_depth_mm"]:
        r0 = np.array([xy[0], xy[1], -depth * 1e-3])
        for model in ff["heart_models"]:
            if model == "magnetic_dipole" and depth != ff["heart_depth_mm"][len(ff["heart_depth_mm"]) // 2]:
                continue
            sources[f"cardiac/{model}/{depth:g}mm"] = {name: dipole_field_gain(a, r0, model) for name, a in ctx.arrays.items()}
    for model in ff["eye_models"]:
        if model == "head_bem":
            sources["ocular/head_bem"] = {name: eye_bem_gain(ctx, a, eyes_head) for name, a in ctx.arrays.items()}
        else:
            sources[f"ocular/{model}"] = {name: sum(dipole_field_gain(a, e, model) for e in eyes_head) for name, a in ctx.arrays.items()}
    log("far-field source fields computed")
    # check: the coil integration reproduces the room-field basis for a uniform field
    chk = {}
    for name, a in ctx.arrays.items():
        b = np.zeros(a.n)
        for i, c in enumerate(environment.coil_integration(a.info, a.coil_def)):
            b[i] = c["w"] @ (c["cosmag"] @ np.array([1.0, 2.0, 3.0]))
        chk[name] = float(np.max(np.abs(b - ctx.basis[name][:, :3] @ np.array([1.0, 2.0, 3.0]))) / np.max(np.abs(b)))
    # covariance per unit variance per moment (isotropic: x, y, z equal)
    cov = {k: {name: g @ g.T for name, g in v.items()} for k, v in sources.items()}
    mid = ff["heart_depth_mm"][len(ff["heart_depth_mm"]) // 2]
    combo_key = f"cardiac/current_dipole/{mid:g}mm"
    share = float(ff["combined_share"])

    def normalised(c):  # unit median variance on the Neuromag magnetometers
        return {name: x / np.median(np.diag(c["squid"])[ctx.mags]) for name, x in c.items()}

    covn = {k: normalised(c) for k, c in cov.items()}
    covn["cardiac+ocular"] = {name: share * covn[combo_key][name] + share * covn["ocular/head_bem"][name] for name in ctx.arrays}
    # projection: noise-normalised energy left by the 8-term projection (primary noise weights)
    proj = {}
    for k, gains in sources.items():
        row = {}
        for name, a in ctx.arrays.items():
            nz = ctx.noise_nom[name]
            w = 1.0 / np.sqrt(nz.intrinsic_var)[:, None]
            g = gains[name]
            pg = nz.projector() @ g
            sets = {"all": np.ones(a.n, bool)} if name != "squid" else {"combined": np.ones(a.n, bool), "mag": a.kinds == "mag", "grad": a.kinds == "grad"}
            for cs, m in sets.items():
                e0 = np.sum((w[m] * g[m]) ** 2, axis=0)
                e1 = np.sum((w[m] * pg[m]) ** 2, axis=0)
                row[f"{name}/{cs}"] = dict(total=float(e1.sum() / e0.sum()), per_axis={ax: float(e1[j] / e0[j]) for j, ax in enumerate("xyz")})
        proj[k] = row
    # context: the cortical targets and the cortical background under the same projection
    context = {}
    for name, a in ctx.arrays.items():
        nz = ctx.noise_nom[name]
        w = 1.0 / np.sqrt(nz.intrinsic_var)[:, None]
        s = ctx.G_t[name]
        frac = np.sum((w * (nz.projector() @ s)) ** 2, axis=0) / np.sum((w * s) ** 2, axis=0)
        bc = nz.brain_cov
        p = nz.projector()
        context[name] = dict(cortical_targets_median=float(np.median(frac)), cortical_targets_p5_p95=[float(np.percentile(frac, 5)), float(np.percentile(frac, 95))],
                             cortical_background=float(np.trace((w * p) @ bc @ (w * p).T) / np.trace(w * bc * w.T)),
                             room_field=float(np.sum((w * (p @ nz.ext_basis)) ** 2) / np.sum((w * nz.ext_basis) ** 2)))
    # levels and detectability with the far-field term added
    env_mag = float(np.sqrt(np.median(np.diag(ctx.noise_nom["squid"].env_cov)[ctx.mags])))
    model_mag = float(np.median(ctx.brain_scale * np.diag(ctx.unit["squid"])[ctx.mags]))
    meas_mag = float(np.nanmedian(ctx.meas["brain"][ctx.mags]))
    levels = {"room_field": env_mag, "magnetometer_shortfall": "magnetometer_shortfall"}
    out = dict(eyes=eyes, heart_head_mm=[[ff["heart_head_mm"][0], ff["heart_head_mm"][1], -d] for d in ff["heart_depth_mm"]],
               coil_integration_check_max_rel=chk, projection_energy_left=proj, projection_context=context,
               level_values_fT=dict(room_field=env_mag * 1e15, magnetometer_shortfall=float(np.sqrt(max(meas_mag - model_mag, 0.0)) * 1e15),
                                    measured_brain_mag=float(np.sqrt(meas_mag) * 1e15), modelled_brain_mag=float(np.sqrt(model_mag) * 1e15)),
               model_notes=dict(
                   cardiac="current dipole: primary-current (Biot-Savart) field of a current dipole in an unbounded homogeneous "
                           "conductor (volume currents of an unbounded medium give no magnetic field; the torso's boundaries are not "
                           "modelled); magnetic dipole: the far field of a bounded current loop. Three orthogonal moments of equal "
                           "variance (the cardiac vector rotates over the beat).",
                   ocular="current dipoles at both eye centres located in the sample T1, the same moment in both eyes (conjugate "
                          "movements), three orthogonal moments of equal variance (horizontal and vertical rotation of the "
                          "corneo-retinal dipole, lid movement); head BEM: homogeneous head inside the BEM head surface, 0.3 S/m."),
               results={})
    modes = ("added", "recalibrated")
    out["calibration_modes"] = dict(
        added="cortical background as published (calibrated on the gradiometers without far field), far field on top: the "
              "adverse bound; the model then over-predicts the measured gradiometer brain noise by the far field's share "
              "(calibration.grad_rms_model_over_measured)",
        recalibrated="cortical scale refitted so that background + far field match the measured gradiometer brain noise: "
                     "consistent with the data; a far field that alone exceeds it is infeasible (no comparisons)")
    for k, c in covn.items():
        for lev_name, lev in levels.items():
            for mode in modes:
                scale, a2, info = calibrate_with_far_field(ctx, ctx.unit["squid"], c["squid"], ctx.target_var, meas_mag, lev, mode)
                key = f"{k}/{lev_name}/{mode}"
                if not info["feasible"]:
                    out["results"][key] = dict(calibration=info, comparisons=None, depth=None)
                    continue
                d, rms = {}, {}
                for name, a in ctx.arrays.items():
                    nz0 = ctx.noise_nom[name]
                    nz = noisemodel.ArrayNoise(nz0.intrinsic_var, scale * ctx.unit[name] + a2 * c[name], nz0.env_cov, nz0.ext_basis)
                    d[name] = detect(ctx.topo(name), nz, a, ctx.conds)
                    dv = np.diag(a2 * c[name])
                    if name == "squid":
                        rms["squid_mag_fT"] = float(np.sqrt(np.median(dv[a.kinds == "mag"])) * 1e15)
                        rms["squid_grad_fT_cm"] = float(np.sqrt(np.median(dv[a.kinds == "grad"])) * 1e13)
                    else:
                        rms[f"{name}_fT"] = float(np.sqrt(np.median(dv)) * 1e15)
                full = k in (combo_key, "ocular/head_bem", "cardiac+ocular")
                out["results"][key] = dict(brain_scale_over_primary=scale / ctx.brain_scale, calibration=info, rms=rms,
                                           comparisons=compare_all(ctx, d, refs=None if full else ["combined"]),
                                           depth=profiles(ctx, d) if (k == "cardiac+ocular" and lev_name == "magnetometer_shortfall"
                                                                      and mode == "added") else None)
        msg = []
        for mode in modes:
            r = out["results"][f"{k}/magnetometer_shortfall/{mode}"]
            msg.append(f"{mode} " + (f"{r['comparisons']['opm_dense/combined/intrinsic+brain']['ratio']:.3f} | "
                                     f"{r['comparisons']['opm_dense/combined/projected']['ratio']:.3f}" if r["comparisons"] else "infeasible"))
        log(f"far field {k}, shortfall level, dense: " + ", ".join(msg))
    ctx.extra["far_field_cov"] = covn
    ctx.extra["eyes"] = eyes
    return out


# ----------------------------------------------------------------------------------------------
# joint: every adverse alternative at once (frequency-resolved, flat signal spectrum)
def part_joint(ctx, rec):
    st, jc, cn = ctx.st, ctx.scfg["joint"], ctx.scfg["coloured_noise"]
    whites = [float(x) for x in ctx.cfg["sensors"]["opm_asd_fT_per_rtHz"]]
    expo = float(cn["exponent"])
    subs = ctx.extra.get("subbands") or make_subbands(ctx, rec, cn["subband_edges_hz"])
    rows0 = ctx.extra.get("subband_rows") or subband_noises(ctx, subs, ctx.unit)
    cx = st.cortex
    _, near_v, area = near_skull_vertices(ctx)
    keep = cx.dist_inner_skull[near_v] * 1e3 >= jc["near_skull_floor_mm"]
    near_v, area = near_v[keep], area[keep]
    unit_ns = {}
    for name in ctx.arrays:
        g = fullres_columns(ctx, name, near_v)
        unit_ns[name] = ctx.unit[name] + (g * area[None, :]) @ g.T
    covn = ctx.extra.get("far_field_cov")
    if covn is None:
        raise RuntimeError("the joint run needs the far-field part")
    ffc = {name: covn["cardiac+ocular"][name] for name in ctx.arrays}
    w = signal_weights(subs, st.filt, "flat")
    corner, lev = float(jc["corner_hz"]), jc["level"]
    variants = {"baseline_white": (ctx.unit, None, 0.0, None), "joint": (unit_ns, ffc, corner, lev),
                "joint_room_field_level": (unit_ns, ffc, corner, "room_field"),
                "joint_without_near_skull": (ctx.unit, ffc, corner, lev),
                "joint_without_far_field": (unit_ns, None, corner, None),
                "joint_without_coloured": (unit_ns, ffc, 0.0, lev)}
    out = dict(definition="All alternatives together, frequency-resolved (flat signal spectrum); near-skull cortex >= {0:g} mm in the "
                          "background; OPM 1/f corner "
                          "{1:g} Hz; cardiac (current dipole, 250 mm) and ocular (head BEM) sources, each half the far-field variance, "
                          "which in every sub-band fills the magnetometer brain-noise shortfall (measured minus modelled median "
                          "magnetometer variance; 'joint') or equals the room field's RMS on the Neuromag magnetometers "
                          "('joint_room_field_level'), added on top of the cortical background calibrated on the gradiometers of each "
                          "sub-band (the far-field mode 'added': adverse to the OPM); G2 targets (the near-skull targets, favourable to "
                          "the OPM, are left out). The shortfall level attributes to heart and eyes everything the cortical model misses on "
                          "the magnetometers, including any environmental or muscular difference between the task and empty-room "
                          "recordings, and then over-predicts the gradiometers (gradiometer_rms_model_over_measured): an upper bound.".format(jc["near_skull_floor_mm"], jc["corner_hz"]), results={})
    cache = {}
    for vname, (unit, ffcov, corner_v, level_v) in variants.items():
        key = (id(unit), ffcov is not None, level_v)
        if key not in cache:
            cache[key] = rows0 if (unit is ctx.unit and ffcov is None) else subband_noises(ctx, subs, unit, ff=ffcov, ff_level=level_v)
        rows = cache[key]
        sq = fr_detect(ctx, subs, rows, w, 0.0, 0.0, expo, ["squid"], ctx.conds)["squid"]
        full = vname in ("baseline_white", "joint", "joint_room_field_level", "joint_without_far_field")
        levels = whites if full else [float(ctx.cfg["sensors"]["opm_asd_primary_fT_per_rtHz"])]
        for white in levels:
            d = {"squid": sq}
            d.update(fr_detect(ctx, subs, rows, w, white, corner_v, expo, ctx.opms, ctx.conds))
            out["results"][f"{vname}/{white:g}fT"] = dict(comparisons=compare_all(ctx, d, refs=None if full else ["combined"]),
                                                          depth=profiles(ctx, d) if (vname in ("joint", "joint_without_far_field")
                                                                                     and white in (15.0, 30.0)) else None,
                                                          far_field_rms_mag_fT=[float(np.sqrt(r["far_field_level2"] * np.median(np.diag(ffcov["squid"])[ctx.mags])) * 1e15)
                                                                                for r in rows] if ffcov is not None else None)
        if ffcov is not None:
            out.setdefault("gradiometer_rms_model_over_measured", {})[vname] = [
                r["calibration"]["grad_rms_model_over_measured"] if r["calibration"] and r["calibration"].get("feasible") else None for r in rows]
        log(f"joint: {vname} done")
    out["most_adverse"] = ("joint_without_far_field: the far-field terms do not act against the OPM in this model (compare joint and "
                           "joint_without_far_field), so the near-skull background with the 1/f corner is the most adverse combination")
    k = "opm_dense/combined/intrinsic+brain"
    log("joint at 15 fT: dense " + ", ".join(f"{v} {out['results'][f'{v}/15fT']['comparisons'][k]['ratio']:.3f}"
                                             for v in ("baseline_white", "joint", "joint_without_far_field")))
    return out


# ----------------------------------------------------------------------------------------------
# figures (manuscript terms only)
COL = {"opm_dense": "#0072B2", "opm_matched": "#009E73"}
LAB = {"opm_dense": "OPM dense (208 sites)", "opm_matched": "OPM matched (98 sites)"}
COND_LAB = {"intrinsic+brain": "sensor + brain noise", "projected": "room interference, after 8-term projection"}
OFFSET = {("opm_dense", "intrinsic+brain"): -0.27, ("opm_dense", "projected"): -0.09,
          ("opm_matched", "intrinsic+brain"): 0.09, ("opm_matched", "projected"): 0.27}


def _ratio_axis(ax, axis="y"):
    """Plain tick labels (ratios near 1 on a linear axis, or a log axis spanning decades)."""
    from matplotlib.ticker import FuncFormatter, NullFormatter

    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(fmt)
    (ax.yaxis if axis == "y" else ax.xaxis).set_minor_formatter(NullFormatter())


def _forest(ax, rows, title):
    """rows: [(label, {(array, cond): comparison} or None for a section heading)]; ratio axis on x (log)."""
    ticks, labels, bold, y = [], [], [], 0.0
    for text, vals in rows:
        if vals is None:
            ticks.append(-y)
            labels.append(text)
            bold.append(True)
            y += 1.0
            continue
        for (a, cond), v in vals.items():
            r, (lo, hi) = v["ratio"], v["ci95_ratio"]
            ax.errorbar(r, -y + OFFSET[(a, cond)], xerr=[[r - lo], [hi - r]], fmt="o" if cond == "intrinsic+brain" else "s",
                        color=COL[a], mfc=COL[a] if cond == "intrinsic+brain" else "white", ms=3.2, capsize=1.3, lw=0.8, mew=0.8)
        ticks.append(-y)
        labels.append(text)
        bold.append(False)
        y += 1.0
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=6.4)
    for t, b in zip(ax.get_yticklabels(), bold):
        if b:
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-y + 0.4, 0.6)
    ax.axvline(1, color="0.5", lw=0.7, zorder=0)
    _ratio_axis(ax, "x")
    ax.set_xlabel("detectability ratio OPM / Neuromag (306 channels)")
    ax.set_title(title, loc="left")


def _vals(get, keys=(("opm_dense", "intrinsic+brain"), ("opm_dense", "projected"), ("opm_matched", "intrinsic+brain"),
                     ("opm_matched", "projected"))):
    return {(a, c): get(f"{a}/combined/{c}") for a, c in keys}


def figures(summary, out_dir):
    import report_style as RS

    RS.apply()
    fig = plt.figure(figsize=(RS.FULL_W, 10.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.15], hspace=0.5, wspace=0.7, left=0.08, right=0.98, top=0.96, bottom=0.1)
    axs = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(3)]
    sw, ns, cn, fr, jt = (summary.get(k) for k in ("sweep", "near_skull", "coloured", "far_field", "joint"))

    # (a) OPM white-noise sweep, and every adverse alternative together
    ax = axs[0][0]
    if sw:
        lv = sw["levels_fT"]
        for a in ("opm_dense", "opm_matched"):
            for cond, ls in (("intrinsic+brain", "-"), ("projected", "--")):
                r = [sw["entries"][f"{x:g}"]["comparisons"][f"{a}/combined/{cond}"] for x in lv]
                ax.plot(lv, [x["ratio"] for x in r], ls, color=COL[a], marker="o", ms=2.5, lw=1.1)
                if cond == "intrinsic+brain":
                    ax.fill_between(lv, [x["ci95_ratio"][0] for x in r], [x["ci95_ratio"][1] for x in r], color=COL[a], alpha=0.13, lw=0)
            if jt:
                for jk, ls, mk in (("joint", "-.", "^"), ("joint_without_far_field", ":", "v")):
                    r = [jt["results"][f"{jk}/{x:g}fT"]["comparisons"][f"{a}/combined/intrinsic+brain"] for x in lv]
                    ax.plot(lv, [x["ratio"] for x in r], ls, color=COL[a], marker=mk, ms=3, lw=1.0)
        ax.axhline(1, color="0.5", lw=0.7, zorder=0)
        ax.set_xscale("log")
        ax.set_xticks(lv)
        _ratio_axis(ax, "x")
        ax.set_xlabel("OPM white noise level [fT/√Hz]")
        ax.set_ylabel("detectability ratio OPM / Neuromag (306)")
        ax.set_title("a  OPM sensor noise level", loc="left")
        h = [plt.Line2D([], [], color="0.3", ls="-", label="sensor + brain noise (band: 95 % CI)"),
             plt.Line2D([], [], color="0.3", ls="--", label="after 8-term projection"),
             plt.Line2D([], [], color="0.3", ls=":", marker="v", ms=3, label="near-skull cortex + 10-Hz 1/f corner (frequency-resolved)"),
             plt.Line2D([], [], color="0.3", ls="-.", marker="^", ms=3, label="the same + heart and eyes (magnetometer shortfall)")]
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - 0.32 * (hi - lo), hi)  # empty space at the bottom for the legend
        ax.legend(handles=h, loc="lower left", fontsize=5.8)

    # (b) near-skull cortex
    if ns:
        rows = []
        p = ns["primary"]
        rows.append(("3-layer BEM (published)", None))
        rows.append(("published model", _vals(lambda k: p["baseline"]["comparisons"][k])))
        for fk in [k for k in p if k.startswith("floor_")]:
            f = p[fk]["floor_mm"]
            rows += [(f"+ near-skull cortex as noise (≥{f:g} mm)", _vals(lambda k, fk=fk: p[fk]["background"]["comparisons"][k])),
                     (f"+ near-skull cortex as sources (≥{f:g} mm)", _vals(lambda k, fk=fk: p[fk]["targets"]["all"][k])),
                     (f"+ both (≥{f:g} mm)", _vals(lambda k, fk=fk: p[fk]["both"]["all"][k]))]
        r1 = ns.get("one_layer_refined")
        if r1:
            rows.append(("check: refined 1-layer BEM", None))
            rows.append(("published model", _vals(lambda k: r1["baseline"]["comparisons"][k])))
            for fk in [k for k in r1 if k.startswith("floor_")]:
                rows.append((f"+ both (≥{r1[fk]['floor_mm']:g} mm)", _vals(lambda k, fk=fk: r1[fk]["both"]["all"][k])))
        _forest(axs[0][1], rows, "b  cortex within 4 mm of the inner skull")

    # (c, d) coloured OPM noise at the primary white level
    if cn:
        corners = [0.0] + cn["corners_hz"]
        xs = np.arange(len(corners))
        methods = (("frequency_resolved/flat", "-", "o"), ("frequency_resolved/spike-wave", "-.", "^"), ("band_variance", ":", "s"))
        for ax, cond, letter in ((axs[1][0], "intrinsic+brain", "c"), (axs[1][1], "projected", "d")):
            for a in ("opm_dense", "opm_matched"):
                for method, ls, mk in methods:
                    vals = [cn["results"][f"{method}/15fT/{c:g}Hz"]["comparisons"][f"{a}/combined/{cond}"] for c in corners]
                    ax.plot(xs, [v["ratio"] for v in vals], ls, marker=mk, ms=2.8, color=COL[a], lw=1.0)
                    if method == "frequency_resolved/flat":
                        ax.fill_between(xs, [v["ci95_ratio"][0] for v in vals], [v["ci95_ratio"][1] for v in vals], color=COL[a], alpha=0.13, lw=0)
                key = "frequency_resolved_squid_measured/flat/15fT/{}Hz"
                if key.format("0") in cn["results"]:
                    vals = [cn["results"][key.format(f"{c:g}")]["comparisons"][f"{a}/combined/{cond}"] for c in corners]
                    ax.plot(xs, [v["ratio"] for v in vals], "none", marker="D", ms=3, mec=COL[a], mfc="white", mew=0.8)
            ax.axhline(1, color="0.5", lw=0.7, zorder=0)
            ax.set_xticks(xs)
            ax.set_xticklabels(["white"] + [f"{c:g} Hz" for c in corners[1:]])
            ax.set_xlabel("OPM 1/f corner frequency (white level 15 fT/√Hz)")
            ax.set_ylabel("detectability ratio OPM / Neuromag (306)")
            ax.set_title(f"{letter}  coloured OPM noise, {'sensor + brain noise' if cond == 'intrinsic+brain' else 'after projection'}", loc="left")
        h = [plt.Line2D([], [], color="0.3", ls="-", marker="o", ms=2.8, label="frequency-resolved, flat signal spectrum (band: 95 % CI)"),
             plt.Line2D([], [], color="0.3", ls="-.", marker="^", ms=2.8, label="frequency-resolved, spike-wave spectrum"),
             plt.Line2D([], [], color="0.3", ls=":", marker="s", ms=2.8, label="band variance (spatial whitening only)"),
             plt.Line2D([], [], color="0.3", ls="none", marker="D", ms=3, mfc="white", label="frequency-resolved, Neuromag noise from the empty room")]
        axs[1][0].legend(handles=h, loc="upper center", bbox_to_anchor=(1.12, -0.2), ncol=2, fontsize=6)

    # (e) far field: energy left by the projection; (f) detectability with the far field added
    if fr:
        ax = axs[2][0]
        keys = list(fr["projection_energy_left"])

        def nice(k):
            p = k.split("/")
            if p[0] == "cardiac":
                return f"heart\n{p[2]}" + ("\n(magn.)" if p[1] == "magnetic_dipole" else "")
            return "eyes" + ("\n(unbounded)" if p[1] != "head_bem" else "")

        sets = [("opm_dense/all", "OPM dense", COL["opm_dense"]), ("opm_matched/all", "OPM matched", COL["opm_matched"]),
                ("squid/combined", "Neuromag (306)", "#000000"), ("squid/mag", "Neuromag magnetometers", "#999999"),
                ("squid/grad", "Neuromag gradiometers", "#555555")]
        xs = np.arange(len(keys))
        wbar = 0.16
        for j, (sk, sl, c) in enumerate(sets):
            ax.bar(xs + (j - 2) * wbar, [fr["projection_energy_left"][k][sk]["total"] for k in keys], wbar, color=c, label=sl)
        ax.set_yscale("log")
        ax.set_ylim(top=12.0)  # room for the legend above the bars
        _ratio_axis(ax, "y")
        ax.set_xticks(xs)
        ax.set_xticklabels([nice(k) for k in keys], fontsize=6)
        ax.set_ylabel("field energy left by the projection")
        ax.set_title("e  heart and eye fields after the 8-term projection", loc="left")
        ax.legend(fontsize=5.6, loc="upper right", ncol=2, columnspacing=0.8, handlelength=1.2)
        base = sw["entries"]["15"]["comparisons"] if sw else None
        rows = []
        if base:
            rows.append(("published model", _vals(lambda k: base[k])))
        for lev, lev_lab in (("room_field", "room-field level"), ("magnetometer_shortfall", "magnetometer-shortfall level")):
            rows.append((lev_lab, None))
            mid = summary["config"]["sensitivity"]["far_field"]["heart_depth_mm"]
            mid = mid[len(mid) // 2]
            shown = (f"cardiac/current_dipole/{mid:g}mm", f"cardiac/magnetic_dipole/{mid:g}mm", "ocular/head_bem", "cardiac+ocular")
            for k in fr["results"]:
                src, lv_, mode = k.rsplit("/", 2)
                if lv_ != lev or src not in shown or mode != "added" or not fr["results"][k]["comparisons"]:
                    continue
                p = src.split("/")
                lab = ("heart + eyes" if src == "cardiac+ocular" else
                       ("heart, magnetic dipole" if p[1] == "magnetic_dipole" else f"heart, {p[2].replace('mm', ' mm')} below")
                       if p[0] == "cardiac" else "eyes" + (" (unbounded)" if p[1] != "head_bem" else ""))
                rows.append((lab, _vals(lambda kk, k=k: fr["results"][k]["comparisons"][kk])))
        _forest(axs[2][1], rows, "f  heart and eye sources added to the published model")
    handles = [plt.Line2D([], [], color=COL[a], marker="o", ls="-", ms=3, label=LAB[a]) for a in ("opm_dense", "opm_matched")]
    handles += [plt.Line2D([], [], color="0.3", marker="o", ls="none", ms=3, label="sensor + brain noise"),
                plt.Line2D([], [], color="0.3", marker="s", mfc="white", ls="none", ms=3, label="with room interference, after 8-term projection")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=6.6, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_dir / "Figure_noise_sensitivity.png")
    plt.close(fig)

    # depth profiles
    fig, axs = plt.subplots(2, 2, figsize=(RS.FULL_W, 5.8), sharex=True)
    curves = []
    if sw:
        for lvl, ls, lw in (("15", "-", 1.3), ("7", (0, (1, 1)), 0.9), ("30", (0, (4, 1.5)), 0.9)):
            curves.append((f"published model, OPM {lvl} fT/√Hz", sw["entries"].get(lvl, {}).get("depth"), "k", ls, lw))
    if ns:
        fk = [k for k in ns["primary"] if k.startswith("floor_")][0]
        curves.append(("near-skull cortex as noise and as sources", ns["primary"][fk]["depth"]["both"], "#E69F00", "-", 1.1))
    if cn:
        cmax = max(cn["corners_hz"])
        curves.append((f"OPM 1/f corner {cmax:g} Hz (frequency-resolved)", cn["results"][f"frequency_resolved/flat/15fT/{cmax:g}Hz"]["depth"],
                       "#CC79A7", "-", 1.1))
    if fr:
        curves.append(("heart + eyes filling the magnetometer shortfall", fr["results"]["cardiac+ocular/magnetometer_shortfall/added"]["depth"],
                       "#D55E00", "-", 1.1))
    if jt:
        curves.append(("near-skull cortex + 10-Hz 1/f corner, OPM 15 fT/√Hz", jt["results"]["joint_without_far_field/15fT"]["depth"],
                       "#56B4E9", "-", 1.1))
    for i, a in enumerate(("opm_dense", "opm_matched")):
        for j, cond in enumerate(("intrinsic+brain", "projected")):
            ax = axs[i, j]
            for name, prof, c, ls, lw in curves:
                p = (prof or {}).get(f"{a}/combined/{cond}")
                if not p:
                    continue
                ok = [b for b in p if b["ratio"] is not None]
                x = [0.5 * (b["lo"] + b["hi"]) for b in ok]
                ax.plot(x, [b["ratio"] for b in ok], color=c, ls=ls, lw=lw, label=name)
                if name.startswith("published model, OPM 15"):
                    ax.fill_between(x, [b["ci95_ratio"][0] for b in ok], [b["ci95_ratio"][1] for b in ok], color=c, alpha=0.12, lw=0)
            ax.axhline(1, color="0.5", lw=0.7, zorder=0)
            ax.set_title(f"{LAB[a]}, {COND_LAB[cond]}", loc="left", fontsize=7.2)
            if i == 1:
                ax.set_xlabel("source depth below the scalp [mm]")
            if j == 0:
                ax.set_ylabel("ratio OPM / Neuromag (306)")
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=6.4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(out_dir / "Figure_noise_sensitivity_depth.png")
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
def headline(summary):
    """The key ratios of every part side by side (OPM / Neuromag 306 channels, median and 95 % CI)."""
    keys = [f"{a}/combined/{c}" for a in ("opm_dense", "opm_matched") for c in ("intrinsic+brain", "projected")]
    rows = {}

    def put(name, comps):
        if comps:
            rows[name] = {k: dict(ratio=comps[k]["ratio"], ci95_ratio=comps[k]["ci95_ratio"]) for k in keys if k in comps}

    sw, ns, cn, fr, jt = (summary.get(k) for k in ("sweep", "near_skull", "coloured", "far_field", "joint"))
    if sw:
        for lv in sw["levels_fT"]:
            put(f"published model, OPM white {lv:g} fT/sqrt(Hz)", sw["entries"][f"{lv:g}"]["comparisons"])
    if ns:
        p = ns["primary"]
        for fk in [k for k in p if k.startswith("floor_")]:
            put(f"near-skull cortex (>= {p[fk]['floor_mm']:g} mm) in the background", p[fk]["background"]["comparisons"])
            put(f"near-skull cortex (>= {p[fk]['floor_mm']:g} mm) as additional targets", p[fk]["targets"]["all"])
            put(f"near-skull cortex (>= {p[fk]['floor_mm']:g} mm) in the background and as targets", p[fk]["both"]["all"])
    if cn:
        for c in [0.0] + cn["corners_hz"]:
            put(f"OPM 1/f corner {c:g} Hz, frequency-resolved (flat signal)", cn["results"][f"frequency_resolved/flat/15fT/{c:g}Hz"]["comparisons"])
            put(f"OPM 1/f corner {c:g} Hz, band variance (spatial whitening only)", cn["results"][f"band_variance/15fT/{c:g}Hz"]["comparisons"])
    if fr:
        for lev in ("room_field", "magnetometer_shortfall"):
            for mode in ("added", "recalibrated"):
                r = fr["results"][f"cardiac+ocular/{lev}/{mode}"]
                if r["comparisons"]:
                    put(f"heart + eyes, {lev.replace('_', ' ')} level, {mode}", r["comparisons"])
    if jt:
        for lv in summary["config"]["adult"]["sensors"]["opm_asd_fT_per_rtHz"]:
            put(f"near-skull cortex + 1/f corner (most adverse combination), OPM white {lv:g} fT/sqrt(Hz)",
                jt["results"][f"joint_without_far_field/{lv:g}fT"]["comparisons"])
            put(f"all alternatives together (far field filling the magnetometer shortfall), OPM white {lv:g} fT/sqrt(Hz)",
                jt["results"][f"joint/{lv:g}fT"]["comparisons"])
            put(f"all alternatives together (far field at the room-field level), OPM white {lv:g} fT/sqrt(Hz)",
                jt["results"][f"joint_room_field_level/{lv:g}fT"]["comparisons"])
            put(f"frequency-resolved reference (white OPM noise) at {lv:g} fT/sqrt(Hz)", jt["results"][f"baseline_white/{lv:g}fT"]["comparisons"])
    return rows


LIMITATIONS = [
    "One adult anatomy; the variants explore model alternatives, the parcel-bootstrap CIs do not include them.",
    "Near-skull cortex: lead fields closer than ~2 mm to the 5,120-triangle inner skull are numerical artefacts in the primary "
    "3-layer model (excluded there); the 0-2 mm cortex enters only the refined 1-layer check.",
    "Coloured noise: only the OPM sensor noise is given a 1/f rise (exponent 1, declared); the Neuromag sensor noise is the "
    "brochure white noise, with the empty-room spectrum as the counterpart. OPM cross-talk and correlated sensor noise are "
    "not modelled.",
    "Far-field sources: idealised conductors (unbounded medium, homogeneous head), declared heart position and isotropic "
    "moments, amplitudes set as stationary band RMS at the Neuromag magnetometers although cardiac and ocular artefacts are "
    "intermittent; Neuromag is charged no SSS/tSSS beyond the common 8-term projection. Two calibrations bracket the term: "
    "'added' (the published cortical background with the far field on top; it over-predicts the measured gradiometer brain "
    "noise by the far field's share) and 'recalibrated' (the cortical scale refitted to the measured gradiometer noise; "
    "infeasible when the far field alone exceeds it). Eye fields with strong gradiometer content lower the Neuromag "
    "detectability more than the OPM's in this model; the magnetometer shortfall they are scaled to also contains whatever "
    "else the cortical model misses (environmental or muscular differences between the task and empty-room recordings).",
    "Detectability only (known-topography matched filter, oracle covariance); the simulated spike detection of the study was "
    "not repeated under these alternatives.",
]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=OUT, help="output directory (default results/g2_noise_sensitivity)")
    p.add_argument("--parts", default=",".join(PARTS), help=f"comma-separated subset of {PARTS}")
    p.add_argument("--max-targets", type=int, default=0, help="test runs: a random subset of the 7,661 targets")
    p.add_argument("--near-fraction", type=float, default=1.0, help="test runs: a random fraction of the near-skull vertices (areas rescaled)")
    p.add_argument("--n-boot", type=int, default=0, help="bootstrap resamples (default from the configuration: 1000)")
    p.add_argument("--subset-seed", type=int, default=1, help="seed of the test subsets")
    p.add_argument("--arrays-cache", default="", help="development aid: pickle of the built arrays (never for results)")
    p.add_argument("--replot", action="store_true", help="redraw the figures from the summary in --out")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.replot:
        figures(json.loads((args.out / "noise_sensitivity_summary.json").read_text()), args.out)
        return
    parts = [x.strip() for x in args.parts.split(",") if x.strip()]
    if unknown := set(parts) - set(PARTS):
        raise SystemExit(f"unknown parts {unknown}")
    if "joint" in parts and "far_field" not in parts:
        raise SystemExit("the joint run needs the far_field part")
    test = bool(args.max_targets) or args.near_fraction < 1.0 or bool(args.n_boot)
    args.fwd_cache = not test  # subsets never write the shared forward cache
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    scfg = tomllib.loads((ROOT / "configs" / "g2_noise_sensitivity.toml").read_text())
    ctx = setup(cfg, scfg, args)
    summary = dict(status=STATUS, test_run=test, config=dict(adult=cfg, sensitivity=scfg),
                   purpose=PURPOSE,
                   conventions=("Median over targets of the paired ratio d_OPM / d_Neuromag (stored as median_log2 and as ratio), 95 % "
                                "CI from a bootstrap over the Desikan-Killiany parcels and the medial-wall labels (scripts/"
                                "g2_adult_comparison.compare), share of targets and of parcels with the OPM higher; d = sqrt(s^T C^+ s) "
                                "with the oracle covariance, 10-nAm cortical-normal dipoles; conditions 'intrinsic+brain' (sensor + "
                                "cortical background) and 'projected' (sensor + background + room field, 8-term external subspace "
                                "removed by the same noise-weighted projection for every array); comparators 'combined' (Neuromag 306 "
                                "channels), 'mag', 'grad' (subsets of the jointly projected data after projection). Depth profiles: "
                                "median ratio per 5-mm depth bin with a parcel-bootstrap 95 % CI per bin."),
                   setup=dict(n_targets=int(len(ctx.tsel)), n_targets_total=int(ctx.st.nt), n_background_grid=int(len(ctx.st.src.grid)),
                              enbw_hz=ctx.st.enbw, brain_scale=ctx.brain_scale, n_boot=ctx.n_boot,
                              arrays={k: a.n for k, a in ctx.arrays.items()}),
                   runtime_s={}, peak_rss_gb={})
    rec = None
    for part in PARTS:
        if part not in parts:
            continue
        t1 = time.time()
        log(f"--- {part}")
        if part in ("coloured", "joint") and rec is None:  # read here, after the near-skull BEM work, to keep the peaks apart
            rec = Recordings(ctx.arrays["squid"].info, cfg["sensors"]["bads"])
            # the copied measurement with the primary filter's response by FFT reproduces g2.measured_noise (Butterworth)
            m = rec.measured(SubbandFilter(ctx.st.filt, 0.0, ctx.st.filt.fs + 1.0))
            summary["setup"]["fft_primary_band_over_butterworth"] = dict(
                brain_grad=float(np.nanmedian(m["brain"][ctx.grads]) / ctx.target_var),
                brain_mag=float(np.nanmedian(m["brain"][ctx.mags]) / np.nanmedian(ctx.meas["brain"][ctx.mags])),
                empty_room_mag=float(np.nanmedian(m["empty_room"][ctx.mags]) / np.nanmedian(ctx.meas["empty_room"][ctx.mags])),
                env_coef_cov_rel_frobenius=float(np.linalg.norm(m["environment"].coef_cov - ctx.env.coef_cov) / np.linalg.norm(ctx.env.coef_cov)))
            log("recordings read (sub-band filtering by FFT): primary band by FFT / Butterworth "
                + ", ".join(f"{k} {v:.4f}" for k, v in summary["setup"]["fft_primary_band_over_butterworth"].items()))
        if part == "sweep":
            summary["sweep"] = part_sweep(ctx)
        elif part == "near_skull":
            summary["near_skull"] = part_near_skull(ctx)
        elif part == "coloured":
            summary["coloured"] = part_coloured(ctx, rec)
        elif part == "far_field":
            summary["far_field"] = part_far_field(ctx)
        elif part == "joint":
            summary["joint"] = part_joint(ctx, rec)
        summary["runtime_s"][part] = time.time() - t1
        summary["peak_rss_gb"][part] = peak_rss_gb()
        io.write_json(summary, args.out / "noise_sensitivity_summary.json")  # checkpoint after every part
    summary["runtime_s"]["total"] = time.time() - T0
    summary["peak_rss_gb"]["total"] = peak_rss_gb()
    summary["headline"] = headline(summary)
    summary["limitations"] = LIMITATIONS
    summary["notes"] = [
        "Every variant changes one element of the adult noise model and recalibrates the cortical background on the measured Neuromag "
        "gradiometer brain noise, as the adult comparison does; the joint run changes them together. Bootstrap CIs resample parcels of "
        "one anatomy; they do not include model uncertainty, which these variants explore.",
        "Near-skull cortex: lead fields of sources closer than ~2 mm to the 5,120-triangle inner-skull mesh are not converged "
        "(see neighbour_energy and gain_change_coarse_to_refined); the primary model includes them only from 2 mm, the refined "
        "1-layer model down to 0 mm, and the inclusion effect is compared between models.",
        "Far-field sources are idealised (an unbounded or homogeneous-head conductor, isotropic moments, declared positions and "
        "levels); their purpose is the projection's differential effect on the arrays and a bounded detectability change, in two "
        "calibrations ('added', 'recalibrated'; see far_field.calibration_modes).",
        "Joint run: every alternative at once, frequency-resolved; its white-noise reference (baseline_white) is the "
        "frequency-resolved detector without them, about 1.5 % below the published band ratio for the dense array.",
    ]
    io.write_json(summary, args.out / "noise_sensitivity_summary.json")
    figures(json.loads((args.out / "noise_sensitivity_summary.json").read_text()), args.out)
    log(f"done in {time.time() - T0:.0f} s, peak RSS {peak_rss_gb():.2f} GB -> {args.out}")


if __name__ == "__main__":
    main()
