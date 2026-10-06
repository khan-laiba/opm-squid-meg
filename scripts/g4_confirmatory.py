#!/usr/bin/env python3
"""G4 confirmatory spike-detection run (pre-specified; Sections 2.5 and S7.5 of the paper).

The exploratory spike study (scripts/g4_epilepsy_adult.py, g4_epilepsy_pediatric.py; results/g4)
chose its endpoint after the analyses (report Section 2.8). Here that endpoint is fixed in advance
(configs/g4_confirmatory.toml [endpoint]) and tested once more on new data:

* Endpoint (primary, confirmatory): dense OPM vs Neuromag's 306 channels; practical scanning
  detector with thresholds frozen at 1 false event per minute on calibration null data; focal
  spikes 10-20 mm below the scalp; the two-sided exact sign-flip test on per-location differences in
  detection counts, Holm-corrected over the nine anatomies; effect size the paired S50 ratio
  (Neuromag / OPM) with its location bootstrap. The statistics are those of
  ``g4_epilepsy_adult.summarise`` (``--check-endpoint-code`` recomputes the committed exploratory
  values from the cached exploratory states with this script's code).
* New data: every random stream is new (spawned from one declared root seed per anatomy and
  purpose); the event locations (the declared stratification rule, restricted to the 10-20 mm band:
  6 per orientation stratum, never the medial wall) are drawn anew, vertex-disjoint from the
  exploratory locations; the candidate dictionary (same construction rule) is drawn anew.
* Null data, each from its own stream: baseline (whiteners), calibration (detector normalisation,
  frozen thresholds, the oracle's z_crit), held-out (realized false-event rates at the frozen
  thresholds; thresholds matched to 1 false event per minute are fitted here, as in
  scripts/study_g4_matched_rate.py) and evaluation (independent of both: realized rates of the frozen
  and of the matched thresholds, so the matched-rate check is out of sample).
* Monte Carlo variability: every event is simulated in ``noise_replicates`` independent noise
  realizations (new event order and noise, same locations and detectors). Replicate 0 is the
  confirmatory run (one realization per event, as declared); the endpoint in every replicate and
  pooled over the replicates is secondary.
* Detector mismatch (configs/g4_confirmatory.toml [[variant]], declared): the same events and null
  data scanned with templates whose stretches are the geometric midpoints of the injected ones and a
  candidate dictionary computed with a 1-layer BEM and a 2-mm / 2-deg coregistration error (the
  localization study's inverse-model mismatch; one draw per anatomy shared by all arrays), with its
  own calibrated thresholds; the mismatch cost of every array and the change of the S50 ratio are
  reported with paired location bootstraps.
* Secondary: the matched-site OPM array against the same comparator, the oracle detector, and the
  matched thresholds.

Anatomies, arrays, placements and the noise model come unchanged from the exploratory contexts
(``g4_epilepsy_adult.adult_context``: the adult at its measured head position;
``g4_epilepsy_pediatric.child_context``: the other heads at top contact). Fragments copied rather
than imported (they live inside larger functions): the noise-generator set-up of
``g4_epilepsy_adult.simulate``, the S50 and paired-ratio statistics of
``g4_epilepsy_adult.summarise`` and ``_holm`` of scripts/report_facts_g4.py.

Usage:
  g4_confirmatory.py [anatomy ...] [--out DIR] [--no-combine]    simulate + summarise (default: all nine)
  g4_confirmatory.py --combine-only [--out DIR]                   Holm over the anatomies, report, figure
  g4_confirmatory.py <anatomy> --resummarise [--out DIR]          summaries again from the saved state
  g4_confirmatory.py --check-endpoint-code [--out DIR]            this code on the exploratory states
Test settings (the outputs are then marked non-confirmatory): --per-stratum N, --replicates R,
--null-scale F, and always --seed S (S other than the declared root seed: no test run draws from the
confirmatory random streams). Outputs (default results/g4_confirm/): g4c_<anatomy>_summary.json and
g4c_<anatomy>_events.csv; g4_confirm_summary.json, G4_confirm_report.md and Figure_G4_confirm.png;
g4c_endpoint_code_check.json. States: cache/g4_confirm/<anatomy>_state.pkl (with --out: <out>/state/).
"""
from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import pickle
import resource
import sys
import time
import tomllib
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import mne  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import g4_epilepsy_adult as G4  # noqa: E402
from opmsquid import detection, environment, forward, g2, goldenholz, ied, io, localization, opm, paths  # noqa: E402

OUT = paths.RESULTS / "g4_confirm"
STATE_DIR = paths.CACHE / "g4_confirm"
CFG_FILE = paths.CONFIGS / "g4_confirmatory.toml"
BASE_CFG_FILE = paths.CONFIGS / "g4_epilepsy.toml"
ANATOMIES = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
DETECTOR_SETS = ("squid/combined", "opm_matched/opm", "opm_dense/opm")
PAIRS = (("opm_dense/opm", "squid/combined"), ("opm_matched/opm", "squid/combined"))
PURPOSES = ("locations", "dictionary", "coregistration", "baseline", "calibration", "oracle", "heldout", "evaluation", "events",
            "bootstrap")
STATUS = ("NEW (G4 confirmatory: the declared spike endpoint re-tested with new seeds, newly drawn locations, independent null "
          "data, noise replicates and a detector-mismatch variant)")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def peak_rss_gb() -> float:
    """Peak resident memory of this process (macOS reports bytes, Linux kilobytes)."""
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r / 1e6


# ----------------------------------------------------------------------------------------------
# configuration, random streams
def load_config() -> tuple[dict, dict, str]:
    text = CFG_FILE.read_text()
    cfg, base = tomllib.loads(text), tomllib.loads(BASE_CFG_FILE.read_text())
    ep = cfg["endpoint"]
    if tuple(ep["anatomies"]) != ANATOMIES:
        raise ValueError("[endpoint] anatomies must be the nine anatomies of the exploratory study")
    if ep["false_events_per_min"] not in base["detector"]["operating_points_per_min"]:
        raise ValueError("the endpoint's operating point is not one of the exploratory study's")
    if tuple(ep["depth_band_mm"]) != G4.DEPTH_BANDS[0]:
        raise ValueError("the endpoint's band must be the exploratory study's first depth band")
    return cfg, base, hashlib.sha1((text + BASE_CFG_FILE.read_text()).encode()).hexdigest()[:12]


def stream(root_seed: int, label: str, purpose: str, *sub: int) -> np.random.Generator:
    """Independent random stream of one anatomy and purpose (SeedSequence spawn key): adding, removing
    or reordering a computation leaves every other stream unchanged."""
    return np.random.default_rng(np.random.SeedSequence(root_seed, spawn_key=(ANATOMIES.index(label), PURPOSES.index(purpose), *sub)))


def item_rng(root_seed: int, label: str, name: str) -> np.random.Generator:
    """Bootstrap stream of one reported statistic (keyed by its name, independent of evaluation order)."""
    return stream(root_seed, label, "bootstrap", zlib.crc32(name.encode()))


def per_stratum(cfg: dict, base: dict) -> int:
    """Locations per orientation stratum of the band: [design] locations_per_stratum if declared, else
    the exploratory design's per depth x orientation stratum (72 / 12 = 6, 18 in the band)."""
    return int(cfg["design"].get("locations_per_stratum") or base["events"]["n_locations"] // (len(G4.DEPTH_BANDS) * len(G4.ORIENT_BANDS)))


def variant_templates(var: dict, injected) -> list[float]:
    if var["template_rule"] != "geometric_midpoints":
        raise ValueError(f"unknown template rule {var['template_rule']!r}")
    s = sorted(float(x) for x in injected)
    return [float(np.sqrt(a * b)) for a, b in zip(s[:-1], s[1:])]


def mode_label(variant: str, rate: float, matched: bool = False) -> str:
    return f"{'practical' if variant == 'primary' else variant}@{rate:g}" + ("_matched" if matched else "")


# ----------------------------------------------------------------------------------------------
# anatomy, locations
def context(label: str) -> G4.Context:
    """The exploratory study's anatomy context: the adult at its measured head position (G2), every
    other head at the G3B primary placement (top contact) with its refitted OPM arrays."""
    if label == "adult":
        return G4.adult_context(tomllib.loads((paths.CONFIGS / "g2_adult.toml").read_text()))
    import g4_epilepsy_pediatric as G4P  # imports G3B; only needed for the smaller heads

    return G4P.child_context(label)


def exploratory_summary(label: str) -> dict:
    return json.loads((paths.RESULTS / "g4" / f"g4_{label}_summary.json").read_text())


def global_index(cortex, hemi: int, vertno: int) -> int:
    n_lh = int(np.sum(cortex.hemi == 0))
    g = int(vertno) + int(hemi) * n_lh
    if cortex.hemi[g] != hemi or cortex.vertno[g] != vertno:
        raise ValueError(f"vertex {hemi}/{vertno} not found in the full-resolution cortex")
    return g


def draw_locations(ctx, band_mm, per: int, exclude: set, rng: np.random.Generator):
    """The exploratory stratification rule (``g4_epilepsy_adult.stratified_locations``) restricted to
    one depth band: per orientation stratum up to ``per`` cortical targets (never the medial wall,
    FreeSurfer 'unknown'), drawn without replacement; targets in ``exclude`` ((hemi, vertno) of the
    exploratory locations) are left out of a stratum's pool whenever ``per`` others remain."""
    d0, d1 = band_mm
    src, cortex = ctx.src, ctx.cortex
    cortical = ~np.char.endswith(src.region.astype(str), "unknown")
    keys = list(zip(cortex.hemi[src.target].tolist(), cortex.vertno[src.target].tolist()))
    original = np.array([k in exclude for k in keys])
    loc, strata, pools = [], [], []
    for j, (o0, o1) in enumerate(G4.ORIENT_BANDS):
        pool = np.flatnonzero(cortical & (src.depth_mm >= d0) & (src.depth_mm < d1) & (src.orientation_deg >= o0)
                              & (src.orientation_deg < o1))
        fresh = pool[~original[pool]]
        use = fresh if len(fresh) >= per else pool
        pick = rng.choice(use, min(per, len(use)), replace=False) if len(use) else np.array([], int)
        loc += [int(x) for x in pick]
        strata += [(G4.DEPTH_BANDS.index(tuple(band_mm)), j)] * len(pick)
        pools.append(dict(orientation_deg=[o0, o1], pool=int(len(pool)), pool_without_exploratory=int(len(fresh)), drawn=int(len(pick)),
                          disjoint_from_exploratory=bool(len(fresh) >= per)))
    return np.array(loc, int), np.array(strata, int), pools


# ----------------------------------------------------------------------------------------------
# simulation
def null_concat(gen, minutes: float, seg_s: float, rng) -> dict:
    acc = {}
    for _ in range(int(round(minutes * 60 / seg_s))):
        for name, y in gen.segment(seg_s, rng).items():
            acc.setdefault(name, []).append(y)
    return {name: np.concatenate(v, axis=1) for name, v in acc.items()}


def simulate(ctx: G4.Context, label: str, cfg: dict, base: dict, settings: dict) -> dict:
    """Null data, frozen detectors (primary and mismatch variants, oracle), independent rate checks and
    the event replicates for one anatomy; returns the state that ``summarise`` reads."""
    t_sim = time.time()
    timing = {}
    sim, ev, dc = base["simulation"], base["events"], base["detector"]
    des, ep = cfg["design"], cfg["endpoint"]
    root = settings["root_seed"]  # the declared root seed, or a test seed (``--seed``)
    arrays, G_t, G_g, env = ctx.arrays, ctx.G_t, ctx.G_g, ctx.env
    minutes = {k: base["null"][f"{k}_min"] * settings["null_scale"] for k in ("baseline", "calibration", "heldout")}
    minutes["evaluation"] = cfg["null"]["evaluation_min"] * settings["null_scale"]
    rate = ep["false_events_per_min"]

    # noise generator (as g4_epilepsy_adult.simulate): one realization shared by all arrays per segment
    specs = {}
    for name, a in arrays.items():
        asd = (np.array([g2.SQUID_ASD[k] for k in a.kinds]) if name == "squid"
               else np.full(a.n, sim["opm_asd_fT_per_rtHz"] * 1e-15))
        specs[name] = ied.ArraySpec(G_g[name], environment.external_basis(a.info, env.r0, a.coil_def), asd)
    gen = ied.NoiseGenerator(specs, ctx.brain_scale * ctx.src.grid_area, env.coef_timecourses, env.coef_cov, ctx.filt,
                             decimate=sim["decimate"], env_fs=env.sfreq)
    fs = gen.fs_out
    det_sets = {}
    for key in DETECTOR_SETS:
        name, cs = key.split("/")
        det_sets[key] = (name, g2.channel_sets(arrays[name])[cs])

    # event locations: new draw, 10-20 mm band, vertex-disjoint from the exploratory locations
    explo = exploratory_summary(label)["locations"]
    exclude = {(int(L["hemi"]), int(L["vertex"])) for L in explo} if des["exclude_original_locations"] else set()
    loc, strata, pools = draw_locations(ctx, ep["depth_band_mm"], settings["per_stratum"], exclude, stream(root, label, "locations"))
    tgt = ctx.src.target[loc]
    topo = {(name, i): G_t[name][:, loc[i]] for name in arrays for i in range(len(loc))}
    explo_g = np.array([global_index(ctx.cortex, L["hemi"], L["vertex"]) for L in explo])
    explo_band = np.array([L["stratum"][0] == G4.DEPTH_BANDS.index(tuple(ep["depth_band_mm"])) for L in explo])
    dist = np.linalg.norm(ctx.cortex.rr[tgt][:, None, :] - ctx.cortex.rr[explo_g][None, :, :], axis=-1) * 1e3

    # detectors: primary (the exploratory construction, new dictionary draw) and the declared mismatch variants
    t0 = time.time()
    valid = np.flatnonzero(ctx.cortex.usable)
    cand = valid[goldenholz.poisson_disk(ctx.cortex.rr[valid], dc["dictionary_spacing_mm"] * 1e-3, stream(root, label, "dictionary"))]
    cand = np.setdiff1d(cand, tgt)
    variants = {"primary": dict(stretches=[float(x) for x in ev["stretches"]], gains={name: ctx.gain(name, cand) for name in arrays},
                                dictionary="truth forward model (3-layer BEM, exact coregistration)")}
    mri_to_head = mne.transforms.invert_transform(ctx.subject.trans)
    for var in cfg.get("variant", []):
        # the analyst's head->MRI transform: the true one with a rigid error (random rotation axis through the head
        # origin and translation direction; the same draw for every variant of this anatomy, shared by all arrays)
        err = localization.perturb_trans(ctx.subject.trans["trans"], stream(root, label, "coregistration"), var["coreg_shift_mm"] * 1e-3,
                                         var["coreg_angle_deg"])
        delta = np.linalg.inv(ctx.subject.trans["trans"]) @ err  # the error alone (head frame)
        trans_err = mne.transforms.Transform("head", "mri", err)
        bem = ctx.subject.bem_model(tuple(var["dictionary_conductivity"]))
        gains = {name: forward.chunked_discrete_gain(a.info, trans_err, ctx.cortex.rr[cand], ctx.cortex.nn[cand], bem,
                                                     coil_def=opm.coil_def_file()).astype(np.float64) for name, a in arrays.items()}
        shift = [1e3 * float(np.linalg.norm(mne.transforms.apply_trans(trans_err, mne.transforms.apply_trans(mri_to_head, ctx.cortex.rr[v]))
                                            - ctx.cortex.rr[v])) for v in tgt]
        stretches = variant_templates(var, ev["stretches"])
        nearest = {f"{x:g}": float(min(abs(t / x - 1.0) for t in stretches)) for x in ev["stretches"]}
        variants[var["name"]] = dict(stretches=stretches, gains=gains, analyst_head_to_mri=err, displacement_at_locations_mm=shift,
                                     error_rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(delta[:3, :3]) - 1) / 2, -1, 1)))),
                                     error_translation_mm=float(1e3 * np.linalg.norm(delta[:3, 3])), nearest_template_offset=nearest,
                                     n_templates=len(stretches),
                                     dictionary=(f"{len(var['dictionary_conductivity'])}-layer BEM, coregistration error "
                                                 f"{var['coreg_shift_mm']:g} mm / {var['coreg_angle_deg']:g} deg"))
    for v in variants.values():
        v["templates"] = {s: ied.filtered_template(ctx.filt, s, sim["decimate"]) for s in v["stretches"]}
    templates = variants["primary"]["templates"]
    timing["dictionary_s"] = time.time() - t0
    memory = dict(after_dictionaries_gb=peak_rss_gb())
    forward._SOLUTIONS.clear()  # no forward computation follows; a refined 3-layer BEM solution holds ~2 GB
    gc.collect()
    log(f"{label}: {len(loc)} locations (pools {[p['pool_without_exploratory'] for p in pools]}), {len(cand)} dictionary candidates, "
        f"variants {list(variants)} ({timing['dictionary_s']:.0f} s)")

    seg_s = sim["segment_s"]
    refr, tol = int(round(dc["refractory_s"] * fs)), int(round(dc["hit_tolerance_s"] * fs))
    t0 = time.time()
    base_null = null_concat(gen, minutes["baseline"], seg_s, stream(root, label, "baseline"))
    whiteners = {key: detection.whitener_from_null(base_null[name][m]) for key, (name, m) in det_sets.items()}
    del base_null
    cal = null_concat(gen, minutes["calibration"], seg_s, stream(root, label, "calibration"))
    dets, thr, oracles, zcrit, ncal = {}, {}, {}, {}, {}
    tops = {key: {i: topo[(name, i)][m] for i in range(len(loc))} for key, (name, m) in det_sets.items()}
    for key, (name, m) in det_sets.items():
        null = cal[name][m]
        for vname, v in variants.items():
            d = detection.ScanDetector.build(null, whiteners[key], v["templates"], v["gains"][name][m], refr)
            _, h = d.events(d.statistic(null)[0])
            dets[key, vname] = d
            thr[f"{key}|{vname}"] = detection.threshold_for_rate(h, minutes["calibration"], rate)
            ncal[f"{key}|{vname}"] = int(np.sum(h > thr[f"{key}|{vname}"]))
        oracles[key], samples = detection.Oracle.calibrate(null, whiteners[key], templates, tops[key], rng=stream(root, label, "oracle", 0))
        zcrit[key] = float(np.quantile(np.concatenate(list(samples.values())), 1 - dc["oracle_alpha"]))
        log(f"calibrated {key}: rank {whiteners[key].rank}, thresholds "
            f"{ {v: round(thr[f'{key}|{v}'], 3) for v in variants} }, oracle z_crit {zcrit[key]:.2f}")
    del cal
    timing["baseline_calibration_s"] = time.time() - t0
    memory["after_calibration_gb"] = peak_rss_gb()

    # held-out null (frozen-threshold rates; the matched thresholds are fitted on it) and the oracle's
    # per-trial false-positive probability at its frozen z_crit on the same independent data: in every
    # segment, its statistic at random times (as Oracle.calibrate draws them, as many in all as on the
    # calibration null), normalised by the calibration's null SD
    t0 = time.time()
    rng_h, rng_o = stream(root, label, "heldout"), stream(root, label, "oracle", 1)
    n_held = int(round(minutes["heldout"] * 60 / seg_s))
    n_oracle = int(np.ceil(20000 / n_held))  # Oracle.calibrate's default of 20,000 samples per topography and template, over the segments
    held = {k: [] for k in thr}
    oracle_held = {key: dict(samples=0, exceed=0) for key in det_sets}
    for _ in range(n_held):
        seg = gen.segment(seg_s, rng_h)
        for key, (name, m) in det_sets.items():
            y = seg[name][m]
            for vname in variants:
                held[f"{key}|{vname}"].append(dets[key, vname].events(dets[key, vname].statistic(y)[0])[1])
            o_h, s_h = detection.Oracle.calibrate(y, whiteners[key], templates, tops[key], n_samples=n_oracle, rng=rng_o)
            z = np.concatenate([s_h[k] * o_h.sigma[k] / oracles[key].sigma[k] for k in s_h])
            oracle_held[key]["samples"] += int(z.size)
            oracle_held[key]["exceed"] += int(np.sum(z > zcrit[key]))
    held = {k: np.concatenate(v) for k, v in held.items()}
    timing["heldout_s"] = time.time() - t0
    memory["after_heldout_gb"] = peak_rss_gb()
    t0 = time.time()
    rng_e = stream(root, label, "evaluation")
    evaluation = {k: [] for k in thr}
    for _ in range(int(round(minutes["evaluation"] * 60 / seg_s))):
        seg = gen.segment(seg_s, rng_e)
        for key, (name, m) in det_sets.items():
            for vname in variants:
                evaluation[f"{key}|{vname}"].append(dets[key, vname].events(dets[key, vname].statistic(seg[name][m])[0])[1])
    evaluation = {k: np.concatenate(v) for k, v in evaluation.items()}
    timing["evaluation_s"] = time.time() - t0
    memory["after_evaluation_gb"] = peak_rss_gb()
    log(f"held-out and evaluation null done ({timing['heldout_s']:.0f} + {timing['evaluation_s']:.0f} s)")

    # events: focal, every location x strength x morphology, in every noise replicate; identical across
    # arrays and detectors; 14 per 30-s segment as in the exploratory study
    events = [(i, float(s), float(x)) for i in range(len(loc)) for s in ev["strengths_nAm"] for x in ev["stretches"]]
    n_rep = settings["replicates"]
    per_seg = int((seg_s - 1.5) // sim["event_spacing_s"])
    times = (1.5 + sim["event_spacing_s"] * np.arange(per_seg)) * fs
    heights = {v: {key: np.zeros((n_rep, len(events))) for key in det_sets} for v in variants}
    oracle_z = {key: np.zeros((n_rep, len(events))) for key in det_sets}
    orders, timing["events_s"] = [], []
    for r in range(n_rep):
        t0 = time.time()
        rng = stream(root, label, "events", r)
        order = rng.permutation(len(events))
        orders.append(order)
        for start in range(0, len(order), per_seg):
            batch = order[start:start + per_seg]
            seg = gen.segment(seg_s, rng)
            for name in arrays:
                for e_idx, t_ in zip(batch, times):
                    i, s, x = events[e_idx]
                    tpl, pk = templates[x]
                    ied.inject(seg[name], topo[(name, i)], tpl, pk, s * 1e-9, int(round(t_)))
            for key, (name, m) in det_sets.items():
                y = seg[name][m]
                for vname in variants:
                    emitted = dets[key, vname].events(dets[key, vname].statistic(y)[0])
                    for e_idx, t_ in zip(batch, times):
                        heights[vname][key][r, e_idx] = detection.event_height(*emitted, int(round(t_)), tol)
                for e_idx, t_ in zip(batch, times):
                    i, s, x = events[e_idx]
                    oracle_z[key][r, e_idx] = oracles[key].statistic(y, i, x, topo[(name, i)][m], int(round(t_)))
        timing["events_s"].append(time.time() - t0)
        log(f"replicate {r + 1}/{n_rep}: {len(events)} events ({timing['events_s'][-1]:.0f} s)")
    timing["simulate_s"] = time.time() - t_sim
    memory["after_events_gb"] = peak_rss_gb()
    timing["peak_rss_so_far_gb"] = memory

    locations = [dict(vertex=int(ctx.cortex.vertno[v]), hemi=int(ctx.cortex.hemi[v]), depth_mm=float(ctx.src.depth_mm[li]),
                      orientation_deg=float(ctx.src.orientation_deg[li]), lobe=str(ctx.src.lobe[li]), region=str(ctx.src.region[li]),
                      stratum=[int(a) for a in strata[k]], nearest_exploratory_location_mm=float(dist[k].min()),
                      nearest_exploratory_location_same_band_mm=float(dist[k][explo_band].min()) if explo_band.any() else None)
                 for k, (v, li) in enumerate(zip(tgt, loc))]
    return dict(label=label, context=ctx.notes, cfg=cfg, base=base, settings=settings, fs=fs, minutes=minutes, events=events,
                orders=[o.tolist() for o in orders], heights=heights, oracle_z=oracle_z, thresholds=thr, calibration_events_above=ncal,
                oracle_z_crit=zcrit, oracle_heldout=oracle_held, null_heights=dict(heldout=held, evaluation=evaluation),
                locations=locations, location_pools=pools, n_exploratory_excluded=len(exclude), n_dictionary=int(len(cand)),
                variants={k: {kk: vv for kk, vv in v.items() if kk not in ("gains", "templates")} for k, v in variants.items()},
                whitener_rank={key: int(w.rank) for key, w in whiteners.items()}, simulated_at_commit=io.RUN_COMMIT, timing=timing)


# ----------------------------------------------------------------------------------------------
# statistics (as g4_epilepsy_adult.summarise; the location is the unit)
def location_table(det, loc_i, stren, strengths, n_loc) -> np.ndarray:
    """Detection rate per location and strength (the location's events at that strength pooled)."""
    table = np.full((n_loc, len(strengths)), np.nan)
    for L in range(n_loc):
        for j, s in enumerate(strengths):
            m = (loc_i == L) & (stren == s)
            if m.any():
                table[L, j] = det[m].mean()
    return table


def s50_summary(table, strengths, rng, n_boot) -> dict:
    """Strength for 50 % detection of the location-pooled curve with its location bootstrap
    (g4_epilepsy_adult.summarise, 'strength_for_50pct_nAm')."""
    n = len(table)
    bnd = np.array([detection.s50_bounds(table[rng.integers(0, n, n)].mean(axis=0), strengths) for _ in range(n_boot)])
    boot = np.where(np.isfinite(bnd[:, 1]), np.maximum(bnd[:, 0], strengths[0]), np.inf)
    q = np.percentile(np.where(np.isfinite(boot), boot, 1e9), [2.5, 97.5])
    at_weakest = float(np.mean(bnd[:, 0] == 0.0))
    return dict(value=detection.s50_from(table.mean(axis=0), strengths),
                ci95=[None if at_weakest >= 0.025 or q[0] > strengths[-1] else float(q[0]), float(q[1]) if q[1] <= strengths[-1] else None],
                share_resamples_not_reached=float(np.mean(~np.isfinite(boot))), share_resamples_at_weakest=at_weakest)


def _censoring(b_den, b_num, names) -> str | None:
    def state(b):
        return None if b[0] == b[1] else ("does not reach 50 %" if not np.isfinite(b[1]) else "reaches 50 % at the weakest strength")

    parts = [f"{n} {s}" for n, s in ((names[1], state(b_num)), (names[0], state(b_den))) if s is not None]
    return "; ".join(parts) if parts else None


def paired_s50_ratio(t_den, t_num, strengths, rng, n_boot, names=("the OPM", "Neuromag")) -> dict:
    """S50(num) / S50(den) of two detection tables on the same locations, with a paired location
    bootstrap; an S50 outside the tested strengths is only bounded, so every resample gives a ratio
    interval and every resample is kept (detection.censored_interval). With den = OPM and
    num = Neuromag this is g4_epilepsy_adult.summarise's 's50_ratio_squid_over_opm'."""
    n = len(t_den)
    pairs = []
    for _ in range(n_boot):
        r_ = rng.integers(0, n, n)
        pairs.append((detection.s50_bounds(t_den[r_].mean(axis=0), strengths), detection.s50_bounds(t_num[r_].mean(axis=0), strengths)))
    rb = np.array([detection.ratio_bounds(bd, bn) for bd, bn in pairs])
    above_d = np.array([not np.isfinite(bd[1]) for bd, _ in pairs])
    above_n = np.array([not np.isfinite(bn[1]) for _, bn in pairs])
    weakest = np.array([bd[0] == 0.0 or bn[0] == 0.0 for bd, bn in pairs])
    b0 = (detection.s50_bounds(t_den.mean(axis=0), strengths), detection.s50_bounds(t_num.mean(axis=0), strengths))
    v_lo, v_hi = detection.ratio_bounds(*b0)
    endpoint_names = names == ("the OPM", "Neuromag")
    return dict(value=v_lo if v_lo == v_hi else None,
                value_censored=detection.censoring_label(*b0) if endpoint_names else _censoring(*b0, names),
                value_bounds=None if v_lo == v_hi else [v_lo if v_lo > 0 else None, v_hi if np.isfinite(v_hi) else None],
                ci95=detection.censored_interval(rb[:, 0], rb[:, 1]),
                share_resamples_neither_reached=float(np.mean(above_d & above_n)),
                **{("share_resamples_opm_not_reached" if endpoint_names else "share_resamples_denominator_not_reached"):
                   float(np.mean(above_d & ~above_n)),
                   ("share_resamples_squid_not_reached" if endpoint_names else "share_resamples_numerator_not_reached"):
                   float(np.mean(above_n & ~above_d))},
                share_resamples_reached_at_weakest=float(np.mean(weakest)))


def paired(det_opm, det_squid, loc_i, stren, strengths, n_loc, rng, n_boot) -> dict:
    """Paired OPM vs Neuromag on identical events: per-location differences in detection counts
    (exact sign-flip test) and the paired S50 ratio Neuromag / OPM."""
    only_o, only_s = det_opm & ~det_squid, ~det_opm & det_squid
    d_loc = np.array([int(only_o[loc_i == L].sum()) - int(only_s[loc_i == L].sum()) for L in range(n_loc)])
    t_o = location_table(det_opm, loc_i, stren, strengths, n_loc)
    t_s = location_table(det_squid, loc_i, stren, strengths, n_loc)
    return dict(n=int(len(det_opm)), n_locations=int(n_loc), detected_only_opm=int(only_o.sum()), detected_only_squid=int(only_s.sum()),
                locations_favouring_opm=int(np.sum(d_loc > 0)), locations_favouring_squid=int(np.sum(d_loc < 0)),
                location_differences=d_loc.tolist(), location_sign_flip_p=detection.sign_flip_p(d_loc),
                s50_ratio_squid_over_opm=paired_s50_ratio(t_o, t_s, strengths, rng, n_boot), p_values="uncorrected")


def ratio_change(tables_ref, tables_var, strengths, rng, n_boot) -> dict:
    """Change of the S50 ratio Neuromag / OPM from the primary detector (ref) to a mismatch variant
    (var): [S50_squid,var / S50_opm,var] / [S50_squid,ref / S50_opm,ref] = (Neuromag's mismatch cost) /
    (the OPM's), with bounds when an S50 is censored and a paired location bootstrap (the same
    resamples for all four curves). Above 1: the mismatch costs Neuromag more than the OPM."""
    (to_r, ts_r), (to_v, ts_v) = tables_ref, tables_var
    n = len(to_r)

    def bounds(ix):
        b = [detection.s50_bounds(t[ix].mean(axis=0), strengths) for t in (to_r, ts_r, to_v, ts_v)]
        lo_r, hi_r = detection.ratio_bounds(b[0], b[1])
        lo_v, hi_v = detection.ratio_bounds(b[2], b[3])
        lo = 0.0 if lo_v == 0.0 or not np.isfinite(hi_r) else lo_v / hi_r
        hi = float("inf") if not np.isfinite(hi_v) or lo_r == 0.0 else hi_v / lo_r
        return lo, hi

    lo0, hi0 = bounds(np.arange(n))
    rb = np.array([bounds(rng.integers(0, n, n)) for _ in range(n_boot)])
    return dict(value=lo0 if lo0 == hi0 else None, value_bounds=None if lo0 == hi0 else [lo0 if lo0 > 0 else None, hi0 if np.isfinite(hi0) else None],
                ci95=detection.censored_interval(rb[:, 0], rb[:, 1]))


def holm(ps: dict) -> dict:
    """Holm-adjusted p-values (copied from scripts/report_facts_g4.py, ``_holm``)."""
    run, out = 0.0, {}
    items = sorted(ps.items(), key=lambda kv: kv[1])
    for i, (k, p) in enumerate(items):
        run = max(run, min(1.0, (len(items) - i) * p))
        out[k] = run
    return out


def rate_equality(k_a: int, k_b: int) -> float:
    """Conditional test of equal false-event rates of two detectors on the same null data (equal
    durations): binomial p of k_a among k_a + k_b at 0.5. Approximate: it treats the counts as
    independent Poisson counts, while the arrays share the background and room noise."""
    return float(binomtest(k_a, k_a + k_b, 0.5).pvalue) if k_a + k_b else 1.0


# ----------------------------------------------------------------------------------------------
def summarise(state: dict, out_dir: Path) -> dict:
    t0 = time.time()
    cfg, base, label = state["cfg"], state["base"], state["label"]
    root, n_boot = state["settings"].get("root_seed", cfg["design"]["root_seed"]), cfg["design"]["bootstrap_resamples"]
    rate = cfg["endpoint"]["false_events_per_min"]
    strengths = np.array(base["events"]["strengths_nAm"], float)
    events = state["events"]
    loc_i = np.array([e[0] for e in events])
    stren = np.array([e[1] for e in events], float)
    n_loc = len(state["locations"])
    variants = list(state["heights"])
    n_rep = next(iter(state["oracle_z"].values())).shape[0]
    minutes, thr = state["minutes"], state["thresholds"]
    held, evl = state["null_heights"]["heldout"], state["null_heights"]["evaluation"]
    thr_m = {k: detection.threshold_for_rate(held[k], minutes["heldout"], rate) for k in thr}
    modes = ["oracle"] + [mode_label(v, rate, m) for v in variants for m in (False, True)]

    def detected(key, mode):  # (n_rep, n_events)
        if mode == "oracle":
            return state["oracle_z"][key] > state["oracle_z_crit"][key]
        matched = mode.endswith("_matched")
        v = next(x for x in variants if mode_label(x, rate, matched) == mode)
        t = thr_m[f"{key}|{v}"] if matched else thr[f"{key}|{v}"]
        return state["heights"][v][key] > t

    scopes = [f"replicate{r}" for r in range(n_rep)] + ["pooled"]

    def scoped(det, scope):
        if scope == "pooled":
            return det.reshape(-1), np.tile(loc_i, n_rep), np.tile(stren, n_rep)
        return det[int(scope[len("replicate"):])], loc_i, stren

    comparisons, s50, tables = {}, {}, {}
    for key in DETECTOR_SETS:
        for mode in modes:
            det = detected(key, mode)
            for scope in scopes:
                d, li, st = scoped(det, scope)
                tables[key, mode, scope] = location_table(d, li, st, strengths, n_loc)
                if scope in ("replicate0", "pooled"):
                    name = f"{key}/{mode}/{scope}"
                    s50[name] = s50_summary(tables[key, mode, scope], strengths, item_rng(root, label, f"s50/{name}"), n_boot)
    for a, b in PAIRS:
        for mode in modes:
            da, db = detected(a, mode), detected(b, mode)
            for scope in scopes:
                (xa, li, st), (xb, _, _) = scoped(da, scope), scoped(db, scope)
                name = f"{a}_vs_{b}/{mode}/{scope}"
                comparisons[name] = paired(xa, xb, li, st, strengths, n_loc, item_rng(root, label, name), n_boot)

    # Monte Carlo variability over the noise replicates (same locations and detectors)
    monte_carlo = {}
    for a, b in PAIRS:
        for mode in modes:
            reps = [comparisons[f"{a}_vs_{b}/{mode}/replicate{r}"] for r in range(n_rep)]
            vals = np.array([np.nan if c["s50_ratio_squid_over_opm"]["value"] is None else c["s50_ratio_squid_over_opm"]["value"] for c in reps])
            ps = np.array([c["location_sign_flip_p"] for c in reps])
            fin = vals[np.isfinite(vals)]
            monte_carlo[f"{a}_vs_{b}/{mode}"] = dict(
                n_replicates=n_rep, s50_ratio=vals.tolist(), n_censored=int(np.sum(~np.isfinite(vals))),
                s50_ratio_median=float(np.median(fin)) if fin.size else None,
                s50_ratio_range=[float(fin.min()), float(fin.max())] if fin.size else None,
                log2_s50_ratio_sd=float(np.std(np.log2(fin), ddof=1)) if fin.size > 1 else None,
                locations_favouring=[[c["locations_favouring_opm"], c["locations_favouring_squid"]] for c in reps],
                location_sign_flip_p=ps.tolist(), n_p_below_alpha=int(np.sum(ps < cfg["endpoint"]["alpha"])),
                note="spread over independent noise realizations of the same events (locations, dictionaries and thresholds fixed)")

    # detector mismatch: every array's cost (S50 mismatch / S50 primary) and the change of the S50 ratio
    mismatch = {}
    for v in variants[1:]:
        for matched in (False, True):
            mp, mv = mode_label("primary", rate, matched), mode_label(v, rate, matched)
            for scope in ("replicate0", "pooled"):
                for key in DETECTOR_SETS:
                    name = f"{v}/cost/{key}/{mv}/{scope}"
                    mismatch[name] = paired_s50_ratio(tables[key, mp, scope], tables[key, mv, scope], strengths,
                                                      item_rng(root, label, name), n_boot, names=("primary", v))
                for a, b in PAIRS:
                    name = f"{v}/ratio_change/{a}_vs_{b}/{mv}/{scope}"
                    mismatch[name] = ratio_change((tables[a, mp, scope], tables[b, mp, scope]), (tables[a, mv, scope], tables[b, mv, scope]),
                                                  strengths, item_rng(root, label, name), n_boot)

    # false events: frozen thresholds on the held-out and on the evaluation null; matched thresholds (fitted on
    # the held-out null) on the evaluation null; equality of rates between arrays on the evaluation null
    false_events = {}
    for k in thr:
        false_events[k] = dict(threshold_frozen=thr[k], threshold_matched=thr_m[k],
                               calibration_events_above_frozen=state["calibration_events_above"][k],
                               heldout_frozen=detection.rate_with_ci(int(np.sum(held[k] > thr[k])), minutes["heldout"]),
                               evaluation_frozen=detection.rate_with_ci(int(np.sum(evl[k] > thr[k])), minutes["evaluation"]),
                               evaluation_matched=detection.rate_with_ci(int(np.sum(evl[k] > thr_m[k])), minutes["evaluation"]))
    rate_tests = {}
    for v in variants:
        for a, b in PAIRS:
            for kind in ("frozen", "matched"):
                ka = false_events[f"{a}|{v}"][f"evaluation_{kind}"]["count"]
                kb = false_events[f"{b}|{v}"][f"evaluation_{kind}"]["count"]
                rate_tests[f"{a}_vs_{b}/{v}/{kind}"] = dict(count_opm=ka, count_squid=kb, conditional_binomial_p=rate_equality(ka, kb),
                                                            note="approximate: the arrays share the background and room noise")
    oracle = {}
    for key in DETECTOR_SETS:
        oh = state["oracle_heldout"][key]
        ci = binomtest(oh["exceed"], oh["samples"]).proportion_ci() if oh["samples"] else None
        oracle[key] = dict(z_crit=state["oracle_z_crit"][key], alpha=base["detector"]["oracle_alpha"], heldout_samples=oh["samples"],
                           heldout_exceed=oh["exceed"], heldout_false_positive_probability=oh["exceed"] / oh["samples"] if oh["samples"] else None,
                           ci95=[ci.low, ci.high] if ci else None)

    ep_name = f"{cfg['endpoint']['opm']}_vs_{cfg['endpoint']['comparator']}/{mode_label('primary', rate)}/replicate{cfg['endpoint']['replicate']}"
    confirmatory = not state["settings"]["overrides"]
    summary = dict(
        status=STATUS, anatomy=label, context=state["context"], confirmatory=confirmatory,
        test_overrides=state["settings"]["overrides"] or None,
        endpoint=dict(definition=cfg["endpoint"], comparison=ep_name, result=comparisons[ep_name]),
        config=dict(confirmatory=cfg, inherited=dict(simulation=base["simulation"], events=base["events"], detector=base["detector"],
                                                    null=base["null"]), config_digest=state["settings"]["config_digest"]),
        declared_choices=dict(root_seed=cfg["design"]["root_seed"], noise_replicates=cfg["design"]["noise_replicates"],
                              noise_replicates_run=n_rep, evaluation_min=cfg["null"]["evaluation_min"],
                              locations_per_stratum=cfg["design"].get("locations_per_stratum", "exploratory design (n_locations / 12)"),
                              locations_per_stratum_run=state["settings"]["per_stratum"],
                              sign_flip_test=("exact (all sign patterns)" if n_loc <= 20 else
                                              "Monte Carlo, 20,000 sign patterns (detection.sign_flip_p, more than 20 locations)"),
                              variants={k: {kk: vv for kk, vv in v.items() if kk in ("stretches", "n_templates", "dictionary",
                                                                                     "nearest_template_offset")}
                                        for k, v in state["variants"].items()},
                              localization_study_reference={k: base["localization"][k] for k in ("coreg_shift_mm", "coreg_angle_deg",
                                                                                                 "inverse_bem")},
                              note=("a variant is a complete detector: its thresholds are calibrated on the same null data to the same "
                                    "false-event rate, so every detector is compared at the same operating point; a bank with fewer "
                                    "templates searches less and gets a lower threshold, a shift common to the arrays")),
        seeds=dict(root_seed=root, declared_root_seed=cfg["design"]["root_seed"], anatomy_index=ANATOMIES.index(label), purposes=list(PURPOSES),
                   rule="np.random.SeedSequence(root_seed, spawn_key=(anatomy_index, purpose_index[, replicate or statistic]))"),
        simulated_at_commit=state["simulated_at_commit"], fs_out=state["fs"], null_minutes=minutes, n_events_per_replicate=len(events),
        n_locations=n_loc, n_dictionary=state["n_dictionary"], whitener_rank=state["whitener_rank"],
        locations=state["locations"], location_pools=state["location_pools"],
        location_checks=dict(n_exploratory_vertices_excluded=state["n_exploratory_excluded"],
                             disjoint_from_exploratory=bool(all(L["nearest_exploratory_location_mm"] > 0 for L in state["locations"])),
                             nearest_exploratory_location_mm_min=float(min(L["nearest_exploratory_location_mm"] for L in state["locations"])),
                             nearest_exploratory_location_mm_median=float(np.median([L["nearest_exploratory_location_mm"] for L in state["locations"]])),
                             depth_mm_median=float(np.median([L["depth_mm"] for L in state["locations"]]))),
        mismatch_geometry={k: dict(analyst_head_to_mri=v["analyst_head_to_mri"], error_rotation_deg=v["error_rotation_deg"],
                                   error_translation_mm=v["error_translation_mm"],
                                   displacement_at_locations_mm_median=float(np.median(v["displacement_at_locations_mm"])),
                                   displacement_at_locations_mm_range=[float(np.min(v["displacement_at_locations_mm"])),
                                                                       float(np.max(v["displacement_at_locations_mm"]))])
                           for k, v in state["variants"].items() if "analyst_head_to_mri" in v},
        false_events=false_events, false_event_rate_equality_on_evaluation_null=rate_tests, oracle=oracle,
        comparisons=comparisons, monte_carlo=monte_carlo, mismatch=mismatch, s50=s50,
        timing_s=dict(state["timing"], summarise_s=time.time() - t0), peak_rss_gb=state.get("peak_rss_gb"))
    out_dir.mkdir(parents=True, exist_ok=True)
    io.write_json(summary, out_dir / f"g4c_{label}_summary.json")
    with open(out_dir / f"g4c_{label}_events.csv", "w", newline="") as fh:
        io.csv_status(fh, f"{STATUS}; {label}" + ("" if confirmatory else " [test settings]"))
        wr = csv.writer(fh)
        wr.writerow(["replicate", "event", "location", "strength_nAm", "stretch"]
                    + [f"{k}_oracle_z" for k in DETECTOR_SETS] + [f"{k}_height_{v}" for k in DETECTOR_SETS for v in variants])
        for r in range(n_rep):
            for e, (i, s, x) in enumerate(events):
                wr.writerow([r, e, i, f"{s:g}", f"{x:g}"] + [f"{state['oracle_z'][k][r, e]:.4f}" for k in DETECTOR_SETS]
                            + [f"{state['heights'][v][k][r, e]:.4f}" for k in DETECTOR_SETS for v in variants])
    res = comparisons[ep_name]
    log(f"{label} endpoint: {res['locations_favouring_opm']}/{res['locations_favouring_squid']} locations, p {res['location_sign_flip_p']:.3g}, "
        f"S50 ratio {detection.format_s50_ratio(res['s50_ratio_squid_over_opm'])} (summarised in {time.time() - t0:.0f} s)")
    return summary


# ----------------------------------------------------------------------------------------------
def combine(out_dir: Path) -> dict | None:
    """Holm over the anatomies for the confirmatory endpoint (and, family by family, for the secondary
    analyses), Monte Carlo variability of the family's conclusion, side by side with the exploratory
    results; writes g4_confirm_summary.json, G4_confirm_report.md and Figure_G4_confirm.png."""
    per = {lab: json.loads((out_dir / f"g4c_{lab}_summary.json").read_text()) for lab in ANATOMIES
           if (out_dir / f"g4c_{lab}_summary.json").exists()}
    if not per:
        log(f"combine: no per-anatomy summaries in {out_dir}")
        return None
    labs = list(per)
    first = per[labs[0]]
    cfg = first["config"]["confirmatory"]
    alpha = cfg["endpoint"]["alpha"]
    commits = sorted({s["simulated_at_commit"] for s in per.values()})
    digests = sorted({s["config"]["config_digest"] for s in per.values()})
    complete = len(per) == len(ANATOMIES)
    confirmatory = (complete and all(s["confirmatory"] for s in per.values()) and len(commits) == 1 and "+dirty" not in commits[0]
                    and len(digests) == 1)
    ep = first["endpoint"]["comparison"]
    n_rep = min(s["monte_carlo"][ep.rsplit("/", 1)[0]]["n_replicates"] for s in per.values())

    def family(name):
        ps = {lab: per[lab]["comparisons"][name]["location_sign_flip_p"] for lab in labs}
        adj = holm(ps)
        return dict(p=ps, holm_p=adj, n_pass=int(sum(v < alpha for v in adj.values())), n=len(ps))

    out = dict(status=STATUS, confirmatory=confirmatory, complete=complete, anatomies=labs, simulated_at_commits=commits,
               config_digests=digests, endpoint=dict(definition=cfg["endpoint"], comparison=ep), families={}, anatomy={})
    out["families"]["endpoint"] = family(ep)
    for a, b in PAIRS:
        for mode in ("practical@1", "practical@1_matched", "oracle", "mismatch@1", "mismatch@1_matched"):
            for scope in ("replicate0", "pooled"):
                name = f"{a}_vs_{b}/{mode}/{scope}"
                if name != ep and all(name in s["comparisons"] for s in per.values()):
                    out["families"][f"secondary/{name}"] = family(name)
    out["monte_carlo"] = {f"replicate{r}": family(ep.replace("/replicate0", f"/replicate{r}")) for r in range(n_rep)}
    out["monte_carlo_summary"] = dict(n_replicates=n_rep, n_pass_per_replicate=[out["monte_carlo"][f"replicate{r}"]["n_pass"] for r in range(n_rep)],
                                      n_replicates_all_pass=int(sum(out["monte_carlo"][f"replicate{r}"]["n_pass"] == len(labs) for r in range(n_rep))),
                                      note="Holm over the anatomies within each noise replicate (same locations and detectors)")
    explo_p = {}
    for lab in labs:
        s = per[lab]
        res = s["comparisons"][ep]
        mc = s["monte_carlo"][ep.rsplit("/", 1)[0]]
        e = exploratory_summary(lab)["paired"]["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]
        explo_p[lab] = e["location_sign_flip_p"]
        fe = s["false_events"]
        out["anatomy"][lab] = dict(
            confirmatory=s["confirmatory"], simulated_at_commit=s["simulated_at_commit"],
            s50_squid_nAm=s["s50"][f"squid/combined/practical@1/replicate0"], s50_dense_nAm=s["s50"][f"opm_dense/opm/practical@1/replicate0"],
            endpoint=dict(locations_favouring_opm=res["locations_favouring_opm"], locations_favouring_squid=res["locations_favouring_squid"],
                          location_sign_flip_p=res["location_sign_flip_p"], holm_p=out["families"]["endpoint"]["holm_p"][lab],
                          s50_ratio_squid_over_opm=res["s50_ratio_squid_over_opm"]),
            pooled=s["comparisons"][ep.replace("/replicate0", "/pooled")]["s50_ratio_squid_over_opm"],
            monte_carlo=dict(s50_ratio=mc["s50_ratio"], location_sign_flip_p=mc["location_sign_flip_p"], log2_s50_ratio_sd=mc["log2_s50_ratio_sd"]),
            matched_array=s["comparisons"]["opm_matched/opm_vs_squid/combined/practical@1/replicate0"],
            oracle=s["comparisons"]["opm_dense/opm_vs_squid/combined/oracle/replicate0"],
            mismatch=s["comparisons"].get("opm_dense/opm_vs_squid/combined/mismatch@1/replicate0"),
            mismatch_ratio_change=s["mismatch"].get("mismatch/ratio_change/opm_dense/opm_vs_squid/combined/mismatch@1/replicate0"),
            mismatch_cost={k: s["mismatch"].get(f"mismatch/cost/{k}/mismatch@1/replicate0") for k in DETECTOR_SETS},
            false_events={k: dict(heldout_frozen=v["heldout_frozen"]["rate_per_min"], evaluation_frozen=v["evaluation_frozen"]["rate_per_min"],
                                  evaluation_matched=v["evaluation_matched"]["rate_per_min"]) for k, v in fe.items()},
            rate_equality=s["false_event_rate_equality_on_evaluation_null"],
            exploratory=dict(locations_favouring_opm=e["locations_favouring_opm"], locations_favouring_squid=e["locations_favouring_squid"],
                             location_sign_flip_p=e["location_sign_flip_p"], s50_ratio_squid_over_opm=e["s50_ratio_squid_over_opm"],
                             source=f"results/g4/g4_{lab}_summary.json :: paired['opm_dense/opm_vs_squid/combined/practical@1'].depth0"),
            location_checks=s["location_checks"], timing_s=s["timing_s"], peak_rss_gb=s["peak_rss_gb"])
    explo_adj = holm(explo_p)
    for lab in labs:
        out["anatomy"][lab]["exploratory"]["holm_p_over_present_anatomies"] = explo_adj[lab]
    io.write_json(out, out_dir / "g4_confirm_summary.json")
    (out_dir / "G4_confirm_report.md").write_text(report(out))
    figure(out, out_dir)
    log(f"combine: {len(labs)} anatomies, endpoint Holm passes {out['families']['endpoint']['n_pass']}/{len(labs)}, confirmatory={confirmatory}")
    return out


def _ratio_txt(sr):
    return detection.format_s50_ratio(sr) if sr else "-"


def report(out: dict) -> str:
    fam = out["families"]["endpoint"]
    ep = out["endpoint"]["definition"]
    L = ["# G4 confirmatory spike-detection run", "",
         f"Status: {'confirmatory (all nine anatomies, one clean commit, declared settings)' if out['confirmatory'] else 'NOT confirmatory (incomplete, test settings, mixed or dirty commits)'}; "
         f"commits {', '.join(out['simulated_at_commits'])}.", "",
         f"Endpoint (fixed before the run, configs/g4_confirmatory.toml): {ep['opm']} vs {ep['comparator']}, {ep['detector']} detector at "
         f"{ep['false_events_per_min']:g} false event per minute (thresholds frozen on the calibration null), focal spikes "
         f"{ep['depth_band_mm'][0]:g}-{ep['depth_band_mm'][1]:g} mm, noise replicate {ep['replicate']}; {ep['test']}; {ep['multiplicity']} "
         f"(alpha {ep['alpha']:g}); effect: {ep['effect']}. New seeds, newly drawn locations (vertex-disjoint from the exploratory "
         "ones), independent null data for fitting thresholds and for evaluating false-event rates.", "",
         "## Endpoint", "",
         "| anatomy | S50 Neuromag / dense [nAm] | S50 ratio [95 % CI] | locations OPM / Neuromag | p | Holm p | exploratory: ratio [CI], locations, p |",
         "|---|---|---|---|---|---|---|"]
    for lab, a in out["anatomy"].items():
        e, x = a["endpoint"], a["exploratory"]
        s_sq, s_de = (("-" if v is None else f"{v:.0f}") for v in (a["s50_squid_nAm"]["value"], a["s50_dense_nAm"]["value"]))
        L.append(f"| {lab} | {s_sq} / {s_de} | "
                 f"{_ratio_txt(e['s50_ratio_squid_over_opm'])} | {e['locations_favouring_opm']} / {e['locations_favouring_squid']} | "
                 f"{e['location_sign_flip_p']:.3g} | {e['holm_p']:.3g} | {_ratio_txt(x['s50_ratio_squid_over_opm'])}, "
                 f"{x['locations_favouring_opm']}/{x['locations_favouring_squid']}, p {x['location_sign_flip_p']:.3g} |")
    L += ["", f"Holm over the anatomies: {fam['n_pass']} of {fam['n']} below {ep['alpha']:g}.", "",
          "## Monte Carlo variability (independent noise realizations of the same events)", "",
          f"Anatomies passing Holm per replicate: {out['monte_carlo_summary']['n_pass_per_replicate']} "
          f"(all pass in {out['monte_carlo_summary']['n_replicates_all_pass']} of {out['monte_carlo_summary']['n_replicates']} replicates).", "",
          "| anatomy | S50 ratio per replicate | p per replicate | SD log2 ratio | pooled ratio [CI] |", "|---|---|---|---|---|"]
    for lab, a in out["anatomy"].items():
        mc = a["monte_carlo"]
        ratios = ", ".join("cens." if v is None else f"{v:.2f}" for v in mc["s50_ratio"])
        ps = ", ".join(f"{p:.2g}" for p in mc["location_sign_flip_p"])
        sd = "-" if mc["log2_s50_ratio_sd"] is None else f"{mc['log2_s50_ratio_sd']:.3f}"
        L.append(f"| {lab} | {ratios} | {ps} | {sd} | {_ratio_txt(a['pooled'])} |")
    L += ["", "## Secondary families (Holm over the anatomies within each family; not confirmatory)", "",
          "| family | anatomies passing | largest Holm p |", "|---|---|---|"]
    for name, f in out["families"].items():
        if name.startswith("secondary/"):
            L.append(f"| {name[len('secondary/'):]} | {f['n_pass']} / {f['n']} | {max(f['holm_p'].values()):.3g} |")
    L += ["", "## Detector mismatch (templates at the geometric midpoints of the injected stretches; dictionary with a 1-layer BEM "
          "and a 2-mm / 2-deg coregistration error)", "",
          "| anatomy | ratio, mismatch [CI] | locations, p | ratio change (Neuromag cost / OPM cost) [CI] | cost Neuromag | cost dense |",
          "|---|---|---|---|---|---|"]
    for lab, a in out["anatomy"].items():
        m = a["mismatch"]
        if m is None:
            continue
        rc = a["mismatch_ratio_change"]
        L.append(f"| {lab} | {_ratio_txt(m['s50_ratio_squid_over_opm'])} | {m['locations_favouring_opm']}/{m['locations_favouring_squid']}, "
                 f"p {m['location_sign_flip_p']:.3g} | {_ratio_txt(rc)} | {_ratio_txt(a['mismatch_cost']['squid/combined'])} | "
                 f"{_ratio_txt(a['mismatch_cost']['opm_dense/opm'])} |")
    L += ["", "## False events per minute (thresholds for 1 per minute)", "",
          "Frozen thresholds (calibration null) on the held-out and on the evaluation null; matched thresholds (fitted on the held-out "
          "null) on the evaluation null, which neither threshold saw.", "",
          "| anatomy | array | detector | held-out, frozen | evaluation, frozen | evaluation, matched |", "|---|---|---|---|---|---|"]
    for lab, a in out["anatomy"].items():
        for k, v in a["false_events"].items():
            arr, var = k.split("|")
            L.append(f"| {lab} | {arr} | {var} | {v['heldout_frozen']:.2f} | {v['evaluation_frozen']:.2f} | {v['evaluation_matched']:.2f} |")
    L += ["", "Simulated spikes in a model; the detector knows the morphology family and (primary variant) the forward model; "
          "no claim about clinical detection is made from this run alone."]
    return "\n".join(L) + "\n"


def figure(out: dict, out_dir: Path):
    labs = list(out["anatomy"])
    rows = [("exploratory (results/g4)", lambda a: a["exploratory"]["s50_ratio_squid_over_opm"], "0.55", "s"),
            ("confirmatory (replicate 0)", lambda a: a["endpoint"]["s50_ratio_squid_over_opm"], "k", "o"),
            ("pooled noise replicates", lambda a: a["pooled"], "tab:blue", "D"),
            ("detector mismatch", lambda a: (a["mismatch"] or {}).get("s50_ratio_squid_over_opm"), "tab:red", "v"),
            ("matched-site array", lambda a: a["matched_array"]["s50_ratio_squid_over_opm"], "tab:green", "^")]
    fig, ax = plt.subplots(figsize=(7.5, 0.55 * len(labs) * len(rows) / 2.2 + 1.5))
    lo_ax, hi_ax = 0.6, 2.5
    for k, lab in enumerate(labs):
        for j, (name, get, col, mk) in enumerate(rows):
            sr = get(out["anatomy"][lab])
            if not sr:
                continue
            y = k + (j - (len(rows) - 1) / 2) * 0.15
            ci = sr.get("ci95") or [None, None]
            a_, b_ = (lo_ax if ci[0] is None else ci[0]), (hi_ax if ci[1] is None else ci[1])
            ax.plot([a_, b_], [y, y], color=col, lw=1.2)
            v = sr.get("value")
            if v is not None:
                ax.plot(v, y, mk, color=col, ms=5, label=name if k == 0 else None)
                continue
            vb = sr.get("value_bounds") or [None, None]  # censored point estimate: its bound, open marker
            v = vb[0] if vb[0] is not None else vb[1]
            if v is not None:
                ax.plot(v, y, mk, color=col, mfc="white", ms=5)
    ax.axvline(1.0, color="0.5", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(lo_ax, hi_ax)
    ticks = [0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels([f"{lab} (Holm p {out['anatomy'][lab]['endpoint']['holm_p']:.2g})" for lab in labs], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("paired S50 ratio Neuromag / OPM (above 1: the OPM detects at lower strength)")
    ax.set_title("G4 confirmatory run: focal spikes 10-20 mm, practical detector at 1 false event/min\n"
                 "(dense OPM unless marked; 95 % location-bootstrap intervals, open ends drawn to the axis; open markers: censored)",
                 fontsize=8)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(out_dir / "Figure_G4_confirm.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
def _explo_detected(st: dict, key: str, mode: str) -> np.ndarray:
    """Detections in an exploratory state (g4_epilepsy_adult.simulate): oracle or practical at 1 per minute."""
    if mode == "oracle":
        return st["rec"][key]["oracle_z"] > st["zcrit"][key]
    return st["rec"][key]["event_height"] > st["thr"][key]["1"]


def check_endpoint_code(out_dir: Path) -> dict:
    """This script's statistics on the exploratory states (cache/g4/<anatomy>_state.pkl): the location
    counts, the sign-flip p and the S50 ratio point estimate must equal the committed exploratory
    values exactly; the bootstrap intervals (another random stream) are compared."""
    cfg, base, _ = load_config()
    strengths = np.array(base["events"]["strengths_nAm"], float)
    res, ok_all = {}, True
    for lab in ANATOMIES:
        path = paths.CACHE / "g4" / f"{lab}_state.pkl"
        if not path.exists():
            res[lab] = dict(available=False)
            continue
        with open(path, "rb") as fh:
            st = pickle.load(fh)
        committed = exploratory_summary(lab)
        if st.get("simulated_at_commit") != committed.get("simulated_at_commit"):
            res[lab] = dict(available=True, error="cached state is not the committed simulation")
            ok_all = False
            continue
        ev = st["events"]
        fam = np.array([e[1] for e in ev])
        li = np.array([e[0] for e in ev])
        band = st["strata"][li, 0]
        sel = (fam == "focal") & (band == 0)
        uniq = np.unique(li[sel])
        loc_i = np.searchsorted(uniq, li[sel])
        stren = np.array([e[2] for e in ev], float)[sel]
        res[lab] = {}
        for a, b in PAIRS:
            for mode in ("practical@1", "oracle"):
                mine = paired(_explo_detected(st, a, mode)[sel], _explo_detected(st, b, mode)[sel], loc_i, stren, strengths, len(uniq),
                              np.random.default_rng(7), 1000)
                theirs = committed["paired"][f"{a}_vs_{b}/{mode}"]["depth0"]
                exact = (mine["locations_favouring_opm"] == theirs["locations_favouring_opm"]
                         and mine["locations_favouring_squid"] == theirs["locations_favouring_squid"]
                         and mine["location_sign_flip_p"] == theirs["location_sign_flip_p"]
                         and mine["s50_ratio_squid_over_opm"]["value"] == theirs["s50_ratio_squid_over_opm"]["value"]
                         and mine["s50_ratio_squid_over_opm"]["value_bounds"] == theirs["s50_ratio_squid_over_opm"]["value_bounds"])
                ci_m, ci_t = mine["s50_ratio_squid_over_opm"]["ci95"], theirs["s50_ratio_squid_over_opm"]["ci95"]
                res[lab][f"{a}_vs_{b}/{mode}"] = dict(exact_match=bool(exact), mine=dict(
                    locations=[mine["locations_favouring_opm"], mine["locations_favouring_squid"]], p=mine["location_sign_flip_p"],
                    ratio=mine["s50_ratio_squid_over_opm"]["value"], ci95=ci_m), committed=dict(
                    locations=[theirs["locations_favouring_opm"], theirs["locations_favouring_squid"]], p=theirs["location_sign_flip_p"],
                    ratio=theirs["s50_ratio_squid_over_opm"]["value"], ci95=ci_t))
                ok_all &= bool(exact)
        log(f"{lab}: endpoint code check {'exact' if all(v['exact_match'] for v in res[lab].values()) else 'MISMATCH'}")
    out = dict(status="check: this script's endpoint statistics on the exploratory states reproduce the committed exploratory values",
               all_exact=bool(ok_all), anatomies=res,
               note="the 95 % intervals use another bootstrap stream than the exploratory summaries and agree only approximately")
    io.write_json(out, out_dir / "g4c_endpoint_code_check.json")
    return out


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("anatomies", nargs="*", help=f"any of {', '.join(ANATOMIES)} (default: all nine)")
    ap.add_argument("--out", type=Path, default=None, help=f"output folder (default {OUT})")
    ap.add_argument("--combine-only", action="store_true")
    ap.add_argument("--no-combine", action="store_true")
    ap.add_argument("--resummarise", action="store_true")
    ap.add_argument("--check-endpoint-code", action="store_true")
    ap.add_argument("--per-stratum", type=int, default=None, help="test setting: locations per orientation stratum")
    ap.add_argument("--replicates", type=int, default=None, help="test setting: noise replicates")
    ap.add_argument("--null-scale", type=float, default=None, help="test setting: factor on every null duration")
    ap.add_argument("--seed", type=int, default=None,
                    help="test setting, required with any other: a root seed other than the declared one, so that no test run "
                         "draws from the confirmatory random streams")
    args = ap.parse_args()
    bad = [a for a in args.anatomies if a not in ANATOMIES]
    if bad:
        ap.error(f"unknown anatomies {bad}; choose from {ANATOMIES}")
    mne.set_log_level("WARNING")
    out_dir = OUT if args.out is None else args.out.resolve()
    state_dir = STATE_DIR if args.out is None else out_dir / "state"
    cfg, base, digest = load_config()
    overrides = {k: v for k, v in (("per_stratum", args.per_stratum), ("replicates", args.replicates), ("null_scale", args.null_scale),
                                   ("seed", args.seed)) if v is not None}
    if overrides and (args.seed is None or args.seed == cfg["design"]["root_seed"]):
        ap.error("test settings need --seed with a value other than the declared root seed (configs/g4_confirmatory.toml)")
    if args.check_endpoint_code:
        r = check_endpoint_code(out_dir)
        log(f"endpoint code check: all exact = {r['all_exact']}")
    labels = list(args.anatomies) or ([] if (args.combine_only or args.check_endpoint_code) else list(ANATOMIES))
    settings = dict(per_stratum=args.per_stratum or per_stratum(cfg, base), replicates=args.replicates or cfg["design"]["noise_replicates"],
                    null_scale=args.null_scale or 1.0, root_seed=cfg["design"]["root_seed"] if args.seed is None else args.seed,
                    overrides=overrides, config_digest=digest)
    for lab in labels if not args.combine_only else []:
        t0 = time.time()
        spath = state_dir / f"{lab}_state.pkl"
        if args.resummarise:
            with open(spath, "rb") as fh:
                state = pickle.load(fh)
        else:
            ctx = context(lab)
            gc.collect()  # the anatomies the context builder loaded besides this one
            t_ctx = time.time() - t0
            rss_ctx = peak_rss_gb()
            log(f"{lab}: context ready ({t_ctx:.0f} s, peak RSS so far {rss_ctx:.2f} GB): {ctx.notes}")
            state = simulate(ctx, lab, cfg, base, settings)
            state["timing"]["context_s"] = t_ctx
            state["timing"]["peak_rss_so_far_gb"] = dict(after_context_gb=rss_ctx, **state["timing"]["peak_rss_so_far_gb"])
            state["peak_rss_gb"] = peak_rss_gb()
            spath.parent.mkdir(parents=True, exist_ok=True)
            with open(spath, "wb") as fh:
                pickle.dump(state, fh)
            del ctx
        summarise(state, out_dir)
        forward._SOLUTIONS.clear()  # refined 3-layer BEM solutions are ~2 GB each
        gc.collect()
        log(f"{lab}: done in {time.time() - t0:.0f} s, peak RSS {peak_rss_gb():.2f} GB")
    if args.combine_only or not args.no_combine:
        combine(out_dir)


if __name__ == "__main__":
    main()
