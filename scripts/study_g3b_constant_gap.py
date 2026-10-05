#!/usr/bin/env python3
"""Helmet fit at a constant gap: the fixed adult helmet against a helmet fitted to each head at the
adult's gap (G3B control; revision of the pediatric comparison).

The G3B counterfactual helmet is scaled with the head (head-circumference ratio), which also brings
the Neuromag sensors closer to most smaller heads than to the adult. This study holds the absolute
gap, the placement convention and each head's OPM arrays fixed:

* Gap-matched helmet ("fitted at the adult's gap"): the Neuromag helmet scaled (``pediatric.
  scaled_helmet``) about the head origin of the laterally centred head (G3B's 'x-centred' pose before
  its top contact, the centre of G3B's 'counterfactual_x-centred' helmet) by the factor k at which the
  median magnetometer-coil-centre-to-scalp distance equals the adult's, subject to the feasibility
  rule of every G3B helmet (no coil centre within the 18-mm Dewar spacing of the scalp and none
  inside the head; where the target gap violates it, k is the smallest feasible factor above it and
  the gap actually reached is larger). Two adult targets, each with the adult under the same rule:
    gap_matched      the adult's own helmet about its laterally centred head (its
                     'counterfactual_x-centred', k = 1): every head at the adult's gap under the
                     construction of the scaled helmet of the headline result (primary);
    gap_matched_top  the adult at top contact in the fixed helmet (the primary placement): every head
                     at the gap of the adult's primary reference (secondary). The adult under this rule
                     is its own helmet shrunk about its laterally centred head to that gap; Delta is also
                     given against the adult at top contact itself.
* Arrays: the five helmets per head (fixed helmet at top contact; G3B's two counterfactual helmets
  scaled with the head; the two gap-matched helmets) against the same two OPM arrays per head (G3B's
  dense and site-matched arrays, unchanged), so the OPM channel count is fixed within each head; every
  convention (band, noise, room field, background scale calibrated once on the adult, 3-layer BEM)
  is G3B's, and the lead fields of G3B's own arrays come from the forward cache.
* Estimands (G3B definitions): D = 20 log10(d_OPM / d_Neuromag) of a 10-nAm cortical-normal dipole
  (known-topography detectability, oracle covariance), area-weighted medians without the medial wall,
  95 % intervals from a bootstrap over parcels (one anatomy; no between-subject variability).
  - Delta under the same rule: D_child - D_adult with each head in its own gap-matched helmet (vertex-
    wise for the scaled adults; parcel-matched for the templates and the children), and at an equal
    OPM site count (the adult's dense array subsampled to the child's site count, as G3B's channel-
    count control).
  - Within-head contrast: D(fixed helmet, top contact) - D(fitted helmet) at every target of one head
    (the OPM array cancels: it equals the change of Neuromag's own detectability between the helmets),
    for each fitted helmet (the two scaled with the head and the two gap-matched).
  - Interaction of helmet fit with head size: the child's within-head contrast minus the adult's,
    which equals Delta(fixed helmet) - Delta(fitted helmet) at every vertex of a scaled adult; the
    medians are not additive, and the residual is reported.
  - Placement band (no new forward model): from the stored G3B placements, Delta over the 12 source-
    blind placements of the fixed helmet, each head against the adult in the same placement
    (difference of the area-weighted medians; placements infeasible for either head left out).
* Regions: area-weighted median D per lobe and per parcel group (the regions of the report's
  region-by-head figure) for every helmet, and per-target D for every helmet in
  targets_<anatomy>.csv.

Configuration: configs/g3b_pediatric.toml and configs/g2_adult.toml (nothing new). Inputs also read:
results/g3b/g3b_summary.json (stored placements, the placement band and the reproduction checks).
Outputs (default results/g3b_constant_gap/; --out for tests): g3b_constant_gap_summary.json,
targets_<anatomy>.csv, Figure_constant_gap.png. State: cache/g3b_constant_gap/state.pkl (with --out:
<out>/state.pkl); --replot redoes summaries, files and figure from it.

Usage: PYTHONPATH=src .venv/bin/python scripts/study_g3b_constant_gap.py [--anatomies adult school ...]
       [--out DIR] [--n-boot N] [--replot]
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import pickle
import resource
import sys
import time
import tomllib
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g3b_pediatric_helmet as G3B  # noqa: E402  (selects the Agg backend)
import report_style as style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from opmsquid import anatomy, background, forward, g2, io, opm, paths, pediatric as P, plotting  # noqa: E402

OUT = ROOT / "results" / "g3b_constant_gap"
STATE = ROOT / "cache" / "g3b_constant_gap" / "state.pkl"
G3B_SUMMARY = ROOT / "results" / "g3b" / "g3b_summary.json"
STATUS = ("NEW (G3B control: the fixed adult helmet against a helmet fitted to each head at the adult's gap; within-head "
          "fixed-versus-fitted contrasts, their interaction with head size, and the placement band of Delta)")

ANATOMIES, CHILDREN, SCALED, SCHOOL, TEMPLATES = G3B.ANATOMIES, G3B.CHILDREN, G3B.SCALED, G3B.SCHOOL, G3B.TEMPLATES
REFS, OPMS = G3B.REFS, G3B.OPMS
# the gap-matched helmets: name -> the adult condition whose median magnetometer-to-scalp gap is the target
RULES = {"gap_matched": "counterfactual_x-centred", "gap_matched_top": "top"}
PRIMARY_RULE = "gap_matched"
HELMETS = ("top", "counterfactual", "counterfactual_x-centred") + tuple(RULES)
FITTED = HELMETS[1:]  # within-head contrasts: the fixed helmet at top contact minus each of these
HELMET_RULE = {
    "top": "fixed adult helmet, head raised (device +z) to 20-mm contact (G3B primary placement)",
    "counterfactual": "helmet scaled by the head-circumference ratio about the head origin of the centred head (G3B counterfactual)",
    "counterfactual_x-centred": "helmet scaled by the head-circumference ratio about the head origin of the laterally centred head "
                                "(G3B counterfactual_x-centred)",
    "gap_matched": "helmet scaled about the head origin of the laterally centred head to the adult's median gap in its own (unscaled) "
                   "helmet at its measured height, laterally centred (the adult's counterfactual_x-centred, k = 1)",
    "gap_matched_top": "helmet scaled about the head origin of the laterally centred head to the adult's median gap at top contact "
                       "in the fixed helmet",
}
# the 12 source-blind placements of the fixed helmet (docs/methods.md G3B 'Placements'; report_facts_g3.FAMILY)
FAMILY = ("top", "back", "x+5mm", "x-5mm", "y+5mm", "y-5mm", "pitch+10deg", "pitch-10deg", "roll+5deg", "roll-5deg",
          "yaw+10deg", "yaw-10deg")
# regions: the six lobes and the parcels and parcel groups of the report's regions by head (report_facts_g3.REGIONS)
REGIONS = {**{lb: ("lobe", lb) for lb in plotting.DK_LOBES},
           "precentral": ("parcels", ("precentral",)), "superiortemporal": ("parcels", ("superiortemporal",)),
           "parahippocampal": ("parcels", ("parahippocampal",)), "entorhinal": ("parcels", ("entorhinal",)),
           "mesial_temporal": ("parcels", ("parahippocampal", "entorhinal")),
           "lateral_temporal": ("parcels", ("superiortemporal", "middletemporal", "inferiortemporal", "bankssts",
                                            "transversetemporal"))}
COND_TOKEN = {"intrinsic+brain": "ib", "projected": "proj"}
OPM_TOKEN = {"opm_dense": "dense", "opm_matched": "matched"}
FAM_CODE = {"D": 1, "within": 2, "same_rule": 3, "vs_adult_top": 4, "equal_channels": 5, "interaction": 6}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r / 1e6  # bytes on macOS, kB on Linux


# ----------------------------------------------------------------------------------------------
# anatomies: G3B's, loaded one at a time (G3B.load_anatomies holds all nine at once)
def load_anatomy(key: str, cfg: dict, g2cfg: dict, adult: G3B.Anatomy | None = None) -> G3B.Anatomy:
    """The anatomy ``key`` exactly as ``g3b_pediatric_helmet.load_anatomies`` builds it (same
    sources, grid and seeds, so G3B's lead fields are found in the forward cache)."""
    seed = g2cfg["sources"]["seed"]
    if list(cfg["anatomy"]["templates"]) != list(TEMPLATES.values()):
        raise ValueError("configs/g3b_pediatric.toml [anatomy] templates must match G3B.TEMPLATES")
    if {c["key"]: c["subject"] for c in cfg["anatomy"]["school"]} != SCHOOL:
        raise ValueError("configs/g3b_pediatric.toml [[anatomy.school]] must match G3B.SCHOOL")
    if key == "adult":
        sub = anatomy.load_sample()
        cor = anatomy.full_resolution(sub)
        an = G3B.Anatomy("adult", sub, cor, g2.make_sources(sub, cor, np.random.default_rng(seed)), P.head_size(sub),
                         "MNE sample subject (as G2)")
    elif key in TEMPLATES:
        sub = anatomy.load_template(TEMPLATES[key])
        cor = anatomy.full_resolution(sub)
        an = G3B.Anatomy(key, sub, cor, g2.make_sources(sub, cor, np.random.default_rng(seed)), P.head_size(sub),
                         f"{G3B.LABEL[key]} ({TEMPLATES[key]}), native dimensions")
    elif key in SCALED:
        if key == "school":
            f, note = cfg["anatomy"]["school_age_scale"], "adult x 85/95 (Jas Table 1 child/adult head radius)"
        else:  # the 2-year template's occipitofrontal circumference over the adult's (G3B 'size_matched_rule')
            f = P.head_size(anatomy.load_template(TEMPLATES["infant2yr"]))["ofc_mm"] / adult.size["ofc_mm"]
            note = "adult x template/adult occipitofrontal circumference"
        sub = anatomy.scaled(adult.subject, f, f"sample_x{f:.4f}")
        cor = anatomy.full_resolution(sub)
        an = G3B.Anatomy(key, sub, cor, G3B.scaled_sources(adult, sub, cor), P.head_size(sub), f"{note}: {f:.4f}")
    elif key in SCHOOL:
        sub = anatomy.load_school(SCHOOL[key])
        cor = anatomy.full_resolution(sub)
        an = G3B.Anatomy(key, sub, cor, g2.make_sources(sub, cor, np.random.default_rng(seed)), P.head_size(sub),
                         f"{G3B.LABEL[key]} ({SCHOOL[key]}), individual MRI, modelled skull")
    else:
        raise ValueError(f"unknown anatomy {key}")
    an.ofc_ratio = an.size["ofc_mm"] / (an.size if adult is None else adult.size)["ofc_mm"]
    return an


# ----------------------------------------------------------------------------------------------
# helmets
def gap_matched_helmet(info: mne.Info, subject, pose: np.ndarray, target: float, bracket=(0.5, 1.5), tol: float = 1e-7,
                       step: float = 1e-3, scan: float = 0.01, grid: float = 1e-6) -> dict:
    """The helmet of ``info`` scaled by k about the head origin of ``pose`` (device-to-head) so that
    the median magnetometer-coil-centre-to-scalp distance is ``target`` [m], then, if a coil centre is
    within the 18-mm Dewar spacing of the scalp or inside the head (``HelmetFit.describe``), the
    smallest feasible k above it (upward scan in ``step``, bisection, then rounded up to a multiple of
    ``grid``, so that two targets below the same limit give the same helmet). The target is bracketed by
    scanning k from 1 in steps of ``scan`` (towards smaller k if the unscaled helmet's gap is larger):
    the gap rises with k while the coils are outside the head (inside it, the unsigned distance rises
    again, so the first crossing from k = 1 is the one sought), then bisection to ``tol``. A target the
    unscaled helmet meets exactly (the adult under its own rule) gives k = 1."""
    centre = np.linalg.inv(pose)[:3, 3]

    def at(k):
        inf = P.scaled_helmet(info, k, centre)
        return inf, P.HelmetFit(inf, subject).describe(pose)

    def gap(k):
        return at(k)[1]["median_dist_mm"] * 1e-3

    g1 = gap(1.0)
    if g1 == target:
        k_t = 1.0
    else:
        d = -scan if g1 > target else scan
        a = 1.0
        while True:
            b = a + d
            if not bracket[0] <= b <= bracket[1]:
                raise RuntimeError(f"target gap {target * 1e3:.2f} mm not reached within the scale range {bracket}")
            if (gap(b) < target) == (d < 0):
                break
            a = b
        lo, hi = sorted((a, b))
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if gap(mid) < target else (lo, mid)
        k_t = hi  # median gap at or just above the target
    inf, desc = at(k_t)
    k = k_t
    if not desc["feasible"]:
        a = b = k_t
        while not at(b)[1]["feasible"]:
            a, b = b, b + step
            if b > bracket[1]:
                raise RuntimeError("no feasible scale factor within the range")
        while b - a > tol:
            mid = 0.5 * (a + b)
            a, b = (a, mid) if at(mid)[1]["feasible"] else (mid, b)
        k = float(np.ceil(b / grid) * grid)  # the boundary on a fixed grid: the same helmet whichever target led to it
        inf, desc = at(k)
        if not desc["feasible"]:
            k = b
            inf, desc = at(k)
    return dict(k=float(k), k_at_target=float(k_t), target_gap_mm=float(target * 1e3), clearance_binding=bool(k != k_t),
                info=inf, trans=pose, **desc)


def helmet_geometry(info: mne.Info, subject, trans: np.ndarray, regions: dict) -> dict:
    """Magnetometer coil-centre-to-scalp distances [mm]: median, percentiles, extremes, the helmet
    region (MNE Vectorview selections) of the closest coil and the median per region."""
    d = P.HelmetFit(info, subject).distances(trans) * 1e3
    i = int(np.argmin(d))
    return dict(median_mm=float(np.median(d)), p10_mm=float(np.percentile(d, 10)), p90_mm=float(np.percentile(d, 90)),
                min_mm=float(d.min()), max_mm=float(d.max()), closest_region=next((k for k, v in regions.items() if i in v), None),
                regions={k: float(np.median(d[v])) for k, v in regions.items()})


def forward_cached(an: G3B.Anatomy, array: g2.Array, bem: list) -> bool:
    """Whether the lead fields of ``array`` at this anatomy's points are in the forward cache (one
    chunk of ``forward.chunked_discrete_gain``)."""
    pts = an.points
    if len(pts) > 20000:
        return False
    key = forward.cache_key(array.info, an.subject.trans, bem, (np.asarray(an.cortex.rr[pts], float), np.asarray(an.cortex.nn[pts], float)),
                            opm.coil_def_file(), "discrete")
    return (paths.CACHE / "fwd" / f"{key}.npz").exists()


# ----------------------------------------------------------------------------------------------
def run_anatomy(an: G3B.Anatomy, com: G3B.Common, cfg: dict, targets: dict | None) -> tuple[dict, dict | None]:
    """Helmets, arrays, lead fields and detectabilities of one anatomy. ``targets`` (rule -> gap [m])
    is None for the adult, whose gaps define it. Returns the state record (and, for the adult, its
    dense OPM array and lead fields for the channel-count control)."""
    t0 = time.time()
    conds = cfg["conditions"]["headline"]
    pc = cfg["placement"]
    pl = P.placements(com.squid_info, an.subject, com.base, pc["translation_mm"] * 1e-3, pc["pitch_deg"], pc["roll_deg"],
                      pc["clearance_mm"] * 1e-3, yaw_deg=pc["yaw_deg"])
    helmets = {"top": dict(info=com.squid_info, trans=pl["top"]["trans"], rule=pl["top"]["rule"], k=1.0, feasible=pl["top"]["feasible"])}
    helmets["counterfactual"] = P.counterfactual_helmet(com.squid_info, an.subject, com.base, an.ofc_ratio)
    cfx = P.counterfactual_helmet(com.squid_info, an.subject, pl["x-centred"]["pose"], an.ofc_ratio)
    cfx["rule"] = cfx["rule"].replace("centred placement", "laterally centred placement (x-centred, before the top contact)")
    helmets["counterfactual_x-centred"] = cfx
    if targets is None:  # the adult: its own gaps are the targets
        targets = {"gap_matched": cfx["median_dist_mm"] * 1e-3, "gap_matched_top": pl["top"]["median_dist_mm"] * 1e-3}
    for rule, target in targets.items():
        helmets[rule] = gap_matched_helmet(com.squid_info, an.subject, pl["x-centred"]["pose"], target)
        helmets[rule]["rule"] = HELMET_RULE[rule] + f"; k = {helmets[rule]['k']:.4f}"
    arrays = {}
    if an.key == "adult":  # the background calibration (G2's rule, at the adult's measured position) comes first
        arrays["squid:centred"] = g2.Array("squid", P.with_dev_head(com.squid_info, pl["centred"]["trans"]), com.kinds, None, {})
    for h, v in helmets.items():
        arrays[f"squid:{h}"] = g2.Array("squid", P.with_dev_head(v["info"], v["trans"]), com.kinds, None, {})
    dig = an.subject.digitisation()
    arrays["opm_dense"] = g2.dense_opm(an.subject, dig, "opm_dense")
    arrays["opm_matched"] = g2.matched_opm(an.subject, dig, squid_info=arrays["squid:top"].info)  # sites at the primary placement
    log(f"{an.key}: {an.nt} targets, {len(an.src.grid)} grid sources; OPM dense {arrays['opm_dense'].n}, matched "
        f"{arrays['opm_matched'].n} sites; helmets k = " + ", ".join(f"{h} {v['k']:.4f}" for h, v in helmets.items())
        + f" ({time.time() - t0:.0f} s)")
    regions = P.helmet_regions(com.squid_info)
    geometry = {h: helmet_geometry(v["info"], an.subject, v["trans"], regions) for h, v in helmets.items()}
    bem3 = an.subject.bem_model(com.bem)
    det, cached, keep = {}, {}, None
    for name, a in arrays.items():
        cached[name] = forward_cached(an, a, bem3)
        t1 = time.time()
        g = G3B.gains(an, a, bem3)
        G_t, G_g = g[:, :an.nt], g[:, an.nt:]
        if name == "squid:centred":
            unit = background.sensor_covariance(G_g, background.moment_covariance(an.src.grid_area))
            com.brain_scale = background.calibrate(unit, com.grads, com.target_var)
            continue
        for (cs, cond), v in G3B.evaluate(an, a, G_t, G_g, com, conds).items():
            det[(name, cs, cond)] = v["detect"]
        if an.key == "adult" and name == "opm_dense":
            keep = dict(array=a, G_t=G_t, G_g=G_g)
        log(f"{an.key}: {name} lead fields {'from the cache' if cached[name] else 'computed'} and evaluated ({time.time() - t1:.0f} s)")
        del g, G_t, G_g
    rec = dict(key=an.key, scale_note=an.scale_note, size=an.size, ofc_ratio=float(an.ofc_ratio), nt=an.nt,
               target=an.src.target, region=an.src.region.astype(str), lobe=an.src.lobe.astype(str), depth_mm=an.src.depth_mm,
               orientation_deg=an.src.orientation_deg, weights=an.weights, cortical=an.cortical,
               hemi=an.cortex.hemi[an.src.target], vertno=an.cortex.vertno[an.src.target], det=det,
               helmets={h: {k: v for k, v in d.items() if k not in ("info", "trans", "pose")} for h, d in helmets.items()},
               geometry=geometry, targets_mm={r: float(t * 1e3) for r, t in targets.items()},
               opm={n: {**{k: v for k, v in arrays[n].meta.items() if isinstance(v, (int, float, str))}, "n_sites": arrays[n].n}
                    for n in OPMS},
               forward_cached=cached, runtime_s=time.time() - t0)
    log(f"{an.key}: done ({rec['runtime_s']:.0f} s; peak RSS so far {peak_rss_gb():.1f} GB)")
    return rec, keep


def compute(keys, cfg, g2cfg) -> dict:
    t_start = time.time()
    com = G3B.Common(cfg, g2cfg)
    adult = load_anatomy("adult", cfg, g2cfg)
    runs, keep, targets = {}, None, None
    for key in keys:  # the adult first: it fixes the background scale and the target gaps
        t0 = time.time()
        an = adult if key == "adult" else load_anatomy(key, cfg, g2cfg, adult)
        log(f"{key}: {an.scale_note}; OFC {an.size['ofc_mm']:.0f} mm; loaded ({time.time() - t0:.0f} s)")
        runs[key], k = run_anatomy(an, com, cfg, targets)
        if key == "adult":
            keep, targets = k, {r: v * 1e-3 for r, v in runs[key]["targets_mm"].items()}
        else:
            del an
        forward._SOLUTIONS.clear()  # the 3-layer BEM solution (~2 GB) of the anatomy just done
        anatomy._REFINED.clear()
        gc.collect()
    # channel count: the adult's dense array subsampled to each smaller head's site count (G3B's control, same lead fields)
    counts = sorted({r["opm"]["opm_dense"]["n_sites"] for k, r in runs.items() if k != "adult"})
    if counts:
        cc = G3B.channel_count_control(adult, dict(_arrays={"opm_dense": keep["array"]}, _G=({"opm_dense": keep["G_t"]},
                                                                                                {"opm_dense": keep["G_g"]})),
                                       com, cfg, counts)
        runs["adult"]["channel_count"] = {n: {k: v["detect"] for k, v in ev.items()} for n, ev in cc.items()}
    return dict(runs=runs, brain_scale=com.brain_scale, enbw=com.enbw, anatomies=list(keys),
                provenance=dict(commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__),
                runtime_s=time.time() - t_start, peak_rss_gb=peak_rss_gb())


# ----------------------------------------------------------------------------------------------
# summaries
def view(rec: dict) -> SimpleNamespace:
    """What G3B's comparison and strata functions read from an anatomy."""
    src = SimpleNamespace(target=rec["target"], region=rec["region"], lobe=rec["lobe"], depth_mm=rec["depth_mm"],
                          orientation_deg=rec["orientation_deg"])
    return SimpleNamespace(key=rec["key"], src=src, weights=rec["weights"], cortical=rec["cortical"], nt=rec["nt"])


def d_target(rec: dict, o: str, helmet: str, ref: str, cond: str) -> np.ndarray:
    """D [dB] per target (OPM array ``o`` minus Neuromag comparator ``ref`` in ``helmet``); NaN on the medial wall."""
    det = rec["det"]
    return np.where(rec["cortical"], 20.0 * np.log10(det[(o, "opm", cond)] / det[(f"squid:{helmet}", ref, cond)]), np.nan)


def fit_target(rec: dict, helmet: str, ref: str, cond: str) -> np.ndarray:
    """Within-head contrast per target: D(fixed helmet, top contact) - D(helmet) = 20 log10(d_helmet / d_top) of
    Neuromag (the same for every OPM array); NaN on the medial wall."""
    det = rec["det"]
    return np.where(rec["cortical"], 20.0 * np.log10(det[(f"squid:{helmet}", ref, cond)] / det[("squid:top", ref, cond)]), np.nan)


def rng_for(seed: int, *codes: int) -> np.random.Generator:
    """One bootstrap stream per summary, fixed by what it summarises (not by the order of computation, so a run on a
    subset of the anatomies gives the same intervals as the full run)."""
    return np.random.default_rng(np.random.SeedSequence([int(seed)] + [int(c) for c in codes]))


def lean(cfg: dict) -> dict:
    """The G3B configuration without depth and orientation strata: ``G3B.compare`` then returns D_child, D_adult and
    Delta (and the parcel table) only."""
    return dict(cfg, strata=dict(cfg["strata"], depth_edges_mm=[], orientation_edges_deg=[]))


def trim(res: dict, keep_parcels: bool) -> dict:
    out = {k: v for k, v in res.items() if not k.startswith("_") and (keep_parcels or k != "parcels")}
    return {k: v for k, v in out.items() if not (isinstance(v, list) and not v)}


def region_medians(rec: dict, x: np.ndarray) -> dict:
    parcel = np.array([r.split(".", 1)[-1] for r in rec["region"]])
    out = {}
    for name, (kind, what) in REGIONS.items():
        m = rec["cortical"] & ((rec["lobe"] == what) if kind == "lobe" else np.isin(parcel, what))
        out[name] = P.weighted_median(x[m], rec["weights"][m]) if m.any() else None
    return out


def region_counts(rec: dict) -> dict:
    parcel = np.array([r.split(".", 1)[-1] for r in rec["region"]])
    return {name: int(np.sum(rec["cortical"] & ((rec["lobe"] == what) if kind == "lobe" else np.isin(parcel, what))))
            for name, (kind, what) in REGIONS.items()}


def placement_band(g3b: dict) -> dict:
    """Delta over the fixed helmet's placement family, each smaller head against the adult in the same placement
    (difference of the stored area-weighted medians of D, dense OPM), from results/g3b/g3b_summary.json."""
    plc, bad, comp = g3b["placement_D"], g3b["infeasible_placements"], g3b["comparisons"]
    out = {}
    for c in CHILDREN:
        for ref in REFS:
            for cond in g3b["config"]["conditions"]["headline"]:
                dc = {p: plc[f"{c}/{p}/{ref}/{cond}"]["median"] for p in FAMILY if p not in bad[c]}
                da = {p: plc[f"adult/{p}/{ref}/{cond}"]["median"] for p in FAMILY if p not in bad["adult"]}
                both = [p for p in FAMILY if p in dc and p in da]
                dlt = {p: dc[p] - da[p] for p in both}
                v = np.array(list(dlt.values()))
                prim = comp[f"{c}/opm_dense/{ref}/{cond}/detect"]["delta"]
                out[f"{c}/{ref}/{cond}"] = dict(
                    n_placements=len(both), excluded_infeasible=sorted(set(FAMILY) - set(both)), delta_by_placement=dlt,
                    min=float(v.min()), max=float(v.max()), span=float(np.ptp(v)), median=float(np.median(v)),
                    top=float(dlt["top"]), rank_of_top_from_lowest=int(sorted(dlt, key=dlt.get).index("top") + 1),
                    crossed=dict(lo=float(min(dc.values()) - max(da.values())), hi=float(max(dc.values()) - min(da.values()))),
                    child_family=dict(min=float(min(dc.values())), max=float(max(dc.values())), n=len(dc)),
                    adult_family=dict(min=float(min(da.values())), max=float(max(da.values())), n=len(da)),
                    primary_delta=dict(median=prim["median"], ci95=prim["ci95"]),
                    primary_minus_top_difference_of_medians=float(prim["median"] - dlt["top"]))
    return out


def summarise(state: dict, cfg: dict, g2cfg: dict, g3b: dict | None, n_boot: dict) -> dict:
    t0 = time.time()
    runs = state["runs"]
    keys = [k for k in ANATOMIES if k in runs]
    kids = [k for k in keys if k != "adult"]
    conds = cfg["conditions"]["headline"]
    seed = g2cfg["sources"]["seed"]
    a, va = runs["adult"], view(runs["adult"])
    ia = {k: ANATOMIES.index(k) for k in keys}
    nb_p, nb_s = n_boot["primary"], n_boot["secondary"]
    out = dict(status=STATUS, anatomies_computed=keys, brain_scale=state["brain_scale"], enbw_hz=state["enbw"])
    out["target_gaps"] = {r: dict(gap_mm=a["targets_mm"][r], adult_condition=RULES[r], rule=HELMET_RULE[r]) for r in RULES}
    out["helmets"] = {k: {h: dict(rule=runs[k]["helmets"][h].get("rule"), k=runs[k]["helmets"][h].get("k"),
                                  k_nominal=runs[k]["helmets"][h].get("k_nominal"), k_at_target=runs[k]["helmets"][h].get("k_at_target"),
                                  target_gap_mm=runs[k]["helmets"][h].get("target_gap_mm"),
                                  clearance_binding=runs[k]["helmets"][h].get("clearance_binding"),
                                  feasible=runs[k]["helmets"][h].get("feasible", True), **runs[k]["geometry"][h])
                          for h in HELMETS} for k in keys}
    out["anatomies"] = {k: dict(scale_note=runs[k]["scale_note"], head_size=runs[k]["size"], ofc_ratio_to_adult=runs[k]["ofc_ratio"],
                                n_targets=runs[k]["nt"], n_cortical_targets=int(runs[k]["cortical"].sum()),
                                opm_arrays=runs[k]["opm"], forward_cached=runs[k]["forward_cached"], runtime_s=runs[k]["runtime_s"],
                                region_targets=region_counts(runs[k])) for k in keys}

    # D per head, helmet, OPM array, comparator and condition
    D = {}
    for k in keys:
        r = runs[k]
        for ih, h in enumerate(HELMETS):
            for io_, o in enumerate(OPMS):
                for ir, ref in enumerate(REFS):
                    for ic, cond in enumerate(conds):
                        x = d_target(r, o, h, ref, cond)
                        nb = nb_p if (h in ("top", PRIMARY_RULE) and o == "opm_dense") else nb_s
                        D[f"{k}/{h}/{o}/{ref}/{cond}"] = G3B.summary_stats(x, r["weights"], r["region"],
                                                                           rng_for(seed, FAM_CODE["D"], ia[k], ih, io_, ir, ic), nb)
    out["D"] = D
    log(f"summaries: D ({time.time() - t0:.0f} s)")

    # within-head contrasts: D(top) - D(fitted helmet), the same for both OPM arrays
    W = {}
    edges = np.array(cfg["strata"]["depth_edges_mm"])
    for k in keys:
        r = runs[k]
        for ih, h in enumerate(FITTED):
            for ir, ref in enumerate(REFS):
                for ic, cond in enumerate(conds):
                    x = fit_target(r, h, ref, cond)
                    nb = nb_p if h == PRIMARY_RULE else nb_s
                    rng = rng_for(seed, FAM_CODE["within"], ia[k], ih, ir, ic)
                    e = G3B.summary_stats(x, r["weights"], r["region"], rng, nb)
                    e["by_depth"] = G3B.strata_rows(x, r["weights"], r["region"], r["depth_mm"], edges, rng,
                                                    nb if (h == PRIMARY_RULE and ref == "combined") else 0, cfg["strata"]["min_n"])
                    e["by_lobe"] = {lb: P.weighted_median(x[r["lobe"] == lb], r["weights"][r["lobe"] == lb]) for lb in plotting.DK_LOBES}
                    W[f"{k}/{h}/{ref}/{cond}"] = e
    out["within_head"] = W
    log(f"summaries: within-head contrasts ({time.time() - t0:.0f} s)")

    # Delta under the same rule (each head in its own gap-matched helmet), against the adult at top contact, and at equal
    # OPM site count; the interaction of the fit with head size
    S, T, E, X = {}, {}, {}, {}
    cc = a.get("channel_count", {})
    for c in kids:
        r, vc = runs[c], view(runs[c])
        for ih, h in enumerate(RULES):
            for io_, o in enumerate(OPMS):
                for ir, ref in enumerate(REFS):
                    for ic, cond in enumerate(conds):
                        primary = h == PRIMARY_RULE and o == "opm_dense"
                        res = G3B.compare(vc, va, d_target(r, o, h, ref, cond), d_target(a, o, h, ref, cond), cfg if primary else lean(cfg),
                                          rng_for(seed, FAM_CODE["same_rule"], ia[c], ih, io_, ir, ic), nb_p if primary else nb_s)
                        S[f"{c}/{h}/{o}/{ref}/{cond}"] = trim(res, primary and ref == "combined" and cond == "intrinsic+brain")
                        if h == "gap_matched_top":
                            res = G3B.compare(vc, va, d_target(r, o, h, ref, cond), d_target(a, o, "top", ref, cond), lean(cfg),
                                              rng_for(seed, FAM_CODE["vs_adult_top"], ia[c], ih, io_, ir, ic), nb_s)
                            T[f"{c}/{h}/{o}/{ref}/{cond}"] = trim(res, False)
        n = r["opm"]["opm_dense"]["n_sites"]
        for ih, h in enumerate(("top",) + tuple(RULES)):
            for ir, ref in enumerate(REFS):
                for ic, cond in enumerate(conds):
                    if n not in cc:
                        continue
                    xa = np.where(a["cortical"], 20 * np.log10(cc[n][("opm", cond)] / a["det"][(f"squid:{h}", ref, cond)]), np.nan)
                    nb = nb_p if h == PRIMARY_RULE else nb_s
                    res = G3B.compare(vc, va, d_target(r, "opm_dense", h, ref, cond), xa, lean(cfg),
                                      rng_for(seed, FAM_CODE["equal_channels"], ia[c], ih, ir, ic), nb)
                    E[f"{c}/{h}/{ref}/{cond}"] = dict(n_sites=int(n), d_child=res["d_child"], d_adult_subsampled=res["d_adult"],
                                                      delta=trim(res["delta"], False) if isinstance(res["delta"], dict) else res["delta"])
        for ih, h in enumerate(FITTED):
            for ir, ref in enumerate(REFS):
                for ic, cond in enumerate(conds):
                    primary = h == PRIMARY_RULE and ref == "combined"
                    res = G3B.compare(vc, va, fit_target(r, h, ref, cond), fit_target(a, h, ref, cond), cfg if primary else lean(cfg),
                                      rng_for(seed, FAM_CODE["interaction"], ia[c], ih, ir, ic), nb_p if h == PRIMARY_RULE else nb_s)
                    # additivity of the medians: Delta(top) - Delta(h) against the interaction (dense OPM; medians only)
                    d_top = G3B.compare(vc, va, d_target(r, "opm_dense", "top", ref, cond), d_target(a, "opm_dense", "top", ref, cond),
                                        lean(cfg), None, 0)["delta"]["median"]
                    d_h = G3B.compare(vc, va, d_target(r, "opm_dense", h, ref, cond), d_target(a, "opm_dense", h, ref, cond),
                                      lean(cfg), None, 0)["delta"]["median"]
                    e = trim(res, primary and cond == "intrinsic+brain")
                    X[f"{c}/{h}/{ref}/{cond}"] = dict(
                        fit_child=e["d_child"], fit_adult=e["d_adult"], interaction=e["delta"], homology=e["homology"],
                        **{k_.replace("delta_", "interaction_"): v for k_, v in e.items() if k_ not in ("d_child", "d_adult", "delta", "homology")},
                        delta_fixed_median=d_top, delta_fitted_median=d_h, additivity_residual=float(d_top - d_h - e["delta"]["median"]))
        log(f"summaries: {c} Delta, equal site count and interaction ({time.time() - t0:.0f} s)")
    out["delta_same_rule"], out["delta_vs_adult_top"], out["delta_equal_channels"], out["interaction"] = S, T, E, X

    # regions: area-weighted median D per lobe and parcel group, every helmet
    out["regions"] = {f"{k}/{h}/{o}/{ref}/{cond}": region_medians(runs[k], d_target(runs[k], o, h, ref, cond))
                      for k in keys for h in HELMETS for o in OPMS for ref in REFS for cond in conds}
    if g3b is not None:
        out["placement_band"] = placement_band(g3b)
        out["checks"] = checks(runs, keys, conds, cfg, g3b)
    out["headline"] = headline(out, keys, g3b)
    out["n_boot"] = dict(n_boot, rule=f"{nb_p} resamples for the dense OPM array in the fixed helmet and the primary gap-matched "
                                      f"helmet ('{PRIMARY_RULE}') and for every summary of that helmet; {nb_s} otherwise (G3B: 1000 for "
                                      "the primary comparisons, 200 for the site-matched array and the other placements)",
                         streams=f"numpy default_rng(SeedSequence([{seed} (configs/g2_adult.toml sources.seed), family, anatomy, "
                                 "helmet, array, comparator, condition indices]))")
    out["notes"] = notes(out, keys)
    out["runtime"] = dict(compute_s=state.get("runtime_s"), compute_peak_rss_gb=state.get("peak_rss_gb"), summaries_s=time.time() - t0)
    return out


def headline(s: dict, keys, g3b: dict | None) -> dict:
    """Dense OPM against Neuromag (306), sensor + brain noise, per head: the gaps reached, D in each helmet, Delta in the
    fixed helmet and in the helmet scaled with the head (stored G3B values) and in the gap-matched helmets (this study), the
    within-head contrasts, the interaction and the placement band. Medians with 95 % intervals."""
    ref, cond, o = "combined", "intrinsic+brain", "opm_dense"

    def mc(e):
        return dict(median=e["median"], ci95=e.get("ci95"))

    out = {}
    for k in keys:
        e = dict(gap_mm={h: s["helmets"][k][h]["median_mm"] for h in HELMETS}, k={h: s["helmets"][k][h]["k"] for h in HELMETS},
                 clearance_binding={r: s["helmets"][k][r]["clearance_binding"] for r in RULES},
                 D={h: mc(s["D"][f"{k}/{h}/{o}/{ref}/{cond}"]) for h in HELMETS},
                 within_head={h: mc(s["within_head"][f"{k}/{h}/{ref}/{cond}"]) for h in FITTED})
        if k != "adult":
            e["delta"] = {r: mc(s["delta_same_rule"][f"{k}/{r}/{o}/{ref}/{cond}"]["delta"]) for r in RULES}
            e["delta"]["gap_matched_top_vs_adult_top"] = mc(s["delta_vs_adult_top"][f"{k}/gap_matched_top/{o}/{ref}/{cond}"]["delta"])
            eq = s["delta_equal_channels"].get(f"{k}/{PRIMARY_RULE}/{ref}/{cond}")
            if eq:
                e["delta"][f"{PRIMARY_RULE}_equal_sites"] = dict(mc(eq["delta"]), n_sites=eq["n_sites"])
            e["interaction"] = {h: mc(s["interaction"][f"{k}/{h}/{ref}/{cond}"]["interaction"]) for h in FITTED}
            if g3b is not None:
                e["delta"]["top_stored_g3b"] = mc(g3b["comparisons"][f"{k}/{o}/{ref}/{cond}/detect"]["delta"])
                for h in ("counterfactual", "counterfactual_x-centred"):
                    e["delta"][f"{h}_stored_g3b"] = mc(g3b["delta_other_placements"][f"{k}/{h}_vs_adult_{h}/{ref}"]["delta"])
                b = s["placement_band"][f"{k}/{ref}/{cond}"]
                e["placement_band"] = dict(min=b["min"], max=b["max"], span=b["span"], n_placements=b["n_placements"])
        out[k] = e
    return dict(array=o, comparator=ref, condition=cond, by_anatomy=out)


def checks(runs: dict, keys, conds, cfg, g3b: dict) -> dict:
    """The recomputed fixed and scaled helmets against the stored G3B results (same code, cache and inputs: the
    medians, which use no random numbers, must agree to rounding)."""
    plc, comp, dop, ccs = g3b["placement_D"], g3b["comparisons"], g3b["delta_other_placements"], g3b["channel_count_control"]
    diff = dict(placement_D_median=0.0, placement_D_by_lobe=0.0, matched_D_median=0.0, primary_delta_median=0.0,
                counterfactual_delta_median=0.0, channel_count_delta_median=0.0, helmet_k=0.0, helmet_median_gap_mm=0.0)
    n = dict.fromkeys(diff, 0)

    def upd(name, v):
        diff[name] = max(diff[name], abs(float(v)))
        n[name] += 1

    a, va = runs["adult"], view(runs["adult"])
    for k in keys:
        r = runs[k]
        for h in ("top", "counterfactual", "counterfactual_x-centred"):
            stored = g3b["placements"][k][h]
            upd("helmet_median_gap_mm", r["geometry"][h]["median_mm"] - stored["median_dist_mm"])
            if h != "top":
                upd("helmet_k", r["helmets"][h]["k"] - stored["k"])
            for ref in REFS:
                for cond in conds:
                    x = d_target(r, "opm_dense", h, ref, cond)
                    e = plc[f"{k}/{h}/{ref}/{cond}"]
                    upd("placement_D_median", P.weighted_median(x, r["weights"]) - e["median"])
                    for lb, v in e["by_lobe"].items():
                        m = r["lobe"] == lb
                        upd("placement_D_by_lobe", P.weighted_median(x[m], r["weights"][m]) - v)
        for ref in REFS:
            for cond in conds:
                upd("matched_D_median", P.weighted_median(d_target(r, "opm_matched", "top", ref, cond), r["weights"])
                    - g3b["D_median_dB"][f"{k}/opm_matched/{ref}/{cond}/detect"])
        if k == "adult":
            continue
        vc = view(r)
        for o in OPMS:
            for ref in REFS:
                for cond in conds:
                    res = G3B.compare(vc, va, d_target(r, o, "top", ref, cond), d_target(a, o, "top", ref, cond), lean(cfg), None, 0)
                    upd("primary_delta_median", res["delta"]["median"] - comp[f"{k}/{o}/{ref}/{cond}/detect"]["delta"]["median"])
        for h in ("counterfactual", "counterfactual_x-centred"):
            for ref in REFS:
                res = G3B.compare(vc, va, d_target(r, "opm_dense", h, ref, "intrinsic+brain"),
                                  d_target(a, "opm_dense", h, ref, "intrinsic+brain"), lean(cfg), None, 0)
                upd("counterfactual_delta_median", res["delta"]["median"] - dop[f"{k}/{h}_vs_adult_{h}/{ref}"]["delta"]["median"])
        cc = a.get("channel_count", {}).get(r["opm"]["opm_dense"]["n_sites"])
        if cc is not None:
            for ref in REFS:
                xa = np.where(a["cortical"], 20 * np.log10(cc[("opm", "intrinsic+brain")] / a["det"][("squid:top", ref, "intrinsic+brain")]),
                              np.nan)
                res = G3B.compare(vc, va, d_target(r, "opm_dense", "top", ref, "intrinsic+brain"), xa, lean(cfg), None, 0)
                upd("channel_count_delta_median", res["delta"]["median"] - ccs[f"{k}/{ref}"]["delta"]["median"])
    return dict(max_abs_difference=diff, n_compared=n,
                what="recomputed medians (dB; gaps mm; k) minus results/g3b/g3b_summary.json for the fixed helmet at top contact, "
                     "the two helmets scaled with the head, the site-matched array, the primary Delta, the counterfactual Delta and "
                     "the channel-count control")


def notes(s: dict, keys) -> list[str]:
    tg = s["target_gaps"]
    bind = [k for k in keys if any(s["helmets"][k][r]["clearance_binding"] for r in RULES)]
    return [
        "D = 20 log10(d_OPM / d_Neuromag) of a 10-nAm cortical-normal dipole (known-topography detectability with the oracle noise "
        "covariance), area-weighted median over cortical targets (medial wall left out). Delta = D_child - D_adult (G3B estimators: "
        "vertex-wise for the scaled adults, area-weighted median of parcel differences for the templates and children). Intervals: "
        "95 %, bootstrap over parcels of one anatomy; no between-subject variability.",
        f"Gap-matched helmets: the Neuromag helmet scaled about the head origin of each head's laterally centred pose (G3B 'x-centred' "
        f"before the top contact) to the adult's median magnetometer-coil-centre-to-scalp gap: {tg['gap_matched']['gap_mm']:.2f} mm "
        f"('{PRIMARY_RULE}': the adult in its own helmet about its laterally centred head, i.e. the adult's counterfactual_x-centred, "
        f"k = 1) and {tg['gap_matched_top']['gap_mm']:.2f} mm ('gap_matched_top': the adult at top contact). No coil centre within "
        "the 18-mm Dewar spacing of the scalp or inside the head (the G3B feasibility rule)"
        + (f"; it binds before the target in {', '.join(bind)}: the gap reached is given under 'helmets'." if bind else "."),
        "The OPM arrays are G3B's (dense, refitted to each head; site-matched, the Neuromag sites at top contact projected to the "
        "scalp) and are not changed with the helmet, so a within-head contrast D(top) - D(fitted) is a property of Neuromag alone "
        "and is the same for both arrays.",
        "Equal site count: the adult's dense array subsampled by farthest-point sampling to the child's dense site count (G3B "
        "channel-count control), with the adult in the same helmet rule.",
        "Interaction = within-head contrast of the child minus that of the adult (same estimator as Delta; fit_child, fit_adult "
        "and interaction take the places of D_child, D_adult and Delta, interaction_by_* its strata, 'parcels' its parcel table). "
        "At every vertex of a scaled adult it equals Delta(fixed helmet) - Delta(fitted helmet); the area-weighted medians are "
        "not additive (additivity_residual = Delta(fixed) - Delta(fitted) - interaction, medians).",
        "Placement band: Delta over the 12 source-blind placements of the fixed helmet (top and back contact, +-5 mm x/y, pitch "
        "+-10 deg, roll +-5 deg, yaw +-10 deg), each head against the adult in the same placement, as the difference of the stored "
        "area-weighted medians of D (dense OPM); placements infeasible for either head left out. 'crossed' pairs any feasible "
        "placement of the child with any of the adult.",
        "A mechanistic control, not a pediatric SQUID system: a scaled helmet keeps the adult coils' sizes, orientations, "
        "integration points and noise."]


# ----------------------------------------------------------------------------------------------
# outputs
def write_targets(state: dict, out_dir: Path, conds) -> list[str]:
    """Per-target D [dB] for every helmet, both OPM arrays, the three Neuromag comparators and both noise conditions."""
    files = []
    for k, r in state["runs"].items():
        cols = [(o, h, ref, cond) for o in OPMS for h in HELMETS for ref in REFS for cond in conds]
        vals = [d_target(dict(r, cortical=np.ones(r["nt"], bool)), o, h, ref, cond) for o, h, ref, cond in cols]
        path = out_dir / f"targets_{k}.csv"
        with open(path, "w", newline="") as fh:
            io.csv_status(fh, STATUS)
            fh.write("# D_<array>_<helmet>_<comparator>_<condition> = 20 log10(d_OPM / d_Neuromag) [dB], 10-nAm cortical-normal "
                     "dipole; helmets: top (fixed adult helmet, top contact), counterfactual and counterfactual_x-centred (scaled "
                     "with the head about the centred and the laterally centred head), gap_matched and gap_matched_top (scaled about "
                     "the laterally centred head to the adult's gap in its own helmet, laterally centred, and at top contact); "
                     "cortical = 0: medial wall (left out of every summary)\n")
            wr = csv.writer(fh)
            wr.writerow(["hemi", "vertno", "cortical", "depth_mm", "orientation_deg", "region", "lobe", "area_mm2"]
                        + [f"D_{OPM_TOKEN[o]}_{h}_{ref}_{COND_TOKEN[cond]}" for o, h, ref, cond in cols])
            for i in range(r["nt"]):
                wr.writerow([int(r["hemi"][i]), int(r["vertno"][i]), int(r["cortical"][i]), f"{r['depth_mm'][i]:.6f}",
                             f"{r['orientation_deg'][i]:.3f}", r["region"][i], r["lobe"][i], f"{r['weights'][i] * 1e6:.4f}"]
                            + [f"{v[i]:.4f}" for v in vals])
        files.append(path.name)
    return files


HSTYLE = {"top": dict(marker="s", color="#000000", mfc="#000000", label="fixed adult helmet, top contact"),
          "counterfactual_x-centred": dict(marker="o", color="#777777", mfc="white", label="helmet scaled with the head (laterally centred)"),
          "gap_matched": dict(marker="o", color="#0072B2", mfc="#0072B2", label="helmet fitted at the adult's gap"),
          "gap_matched_top": dict(marker="^", color="#0072B2", mfc="white", label="helmet fitted at the adult's top-contact gap")}


def figure(s: dict, g3b: dict | None, path: Path) -> dict:
    """Four panels in the manuscript's terms: gaps, D per helmet, Delta, and the within-head contrast."""
    style.apply()
    keys = [k for k in style.ANAT_ORDER if k in s["helmets"]]
    kids = [k for k in keys if k != "adult"]
    xs = {k: i for i, k in enumerate(keys)}
    cond, ref, o = "intrinsic+brain", "combined", "opm_dense"
    tg = s["target_gaps"]
    fig = plt.figure(figsize=(style.FULL_W, 8.2))
    gs = fig.add_gridspec(2, 2, left=0.11, right=0.985, top=0.965, bottom=0.31, hspace=0.62, wspace=0.34)
    axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)])

    def xaxis(ax, ks, zero=True):
        ax.set_xticks([xs[k] for k in ks], [style.ANAT_SHORT[k] for k in ks], rotation=40, ha="right")
        ax.set_xlim(min(xs[k] for k in ks) - 0.6, max(xs[k] for k in ks) + 0.6)
        if zero:
            ax.axhline(0, color="0.75", lw=0.6, zorder=0)

    def point(ax, x, med, ci, st, dx=0.0, label=None):
        lo, hi = (ci if ci else (med, med))
        ax.errorbar([x + dx], [med], yerr=[[med - lo], [hi - med]], fmt=st["marker"], color=st["color"], mfc=st["mfc"], ms=4.2,
                    lw=0.9, capsize=1.6, label=label)

    # (a) gaps
    ax = axs[0, 0]
    shown = ("top", "counterfactual_x-centred", "gap_matched", "gap_matched_top")
    for j, h in enumerate(shown):
        st = HSTYLE[h]
        dx = (j - 1.5) * 0.13
        ys = [s["helmets"][k][h]["median_mm"] for k in keys]
        ax.plot([xs[k] + dx for k in keys], ys, st["marker"], color=st["color"], mfc=st["mfc"], ms=4.2, ls="none", label=st["label"])
        for k in keys:
            if h in RULES and s["helmets"][k][h]["clearance_binding"]:
                ax.annotate("*", (xs[k] + dx, s["helmets"][k][h]["median_mm"]), xytext=(3, 1), textcoords="offset points", fontsize=9,
                            color=st["color"])
    for r, ls in (("gap_matched", "--"), ("gap_matched_top", ":")):
        ax.axhline(tg[r]["gap_mm"], color="#0072B2", lw=0.7, ls=ls, zorder=0)
    ax.set_ylabel("median magnetometer-to-scalp\ngap (mm)")
    ax.set_title("(a)  Neuromag gap in each helmet", loc="left")
    xaxis(ax, keys, zero=False)

    # (b) D per helmet
    ax = axs[0, 1]
    for j, h in enumerate(("top", "counterfactual_x-centred", "gap_matched")):
        st = HSTYLE[h]
        for k in keys:
            e = s["D"][f"{k}/{h}/{o}/{ref}/{cond}"]
            point(ax, xs[k], e["median"], e["ci95"], st, (j - 1) * 0.18)
    ax.set_ylabel("D, dense OPM minus\nNeuromag (dB)")
    ax.set_title("(b)  OPM advantage D in each helmet", loc="left")
    xaxis(ax, keys)

    # (c) Delta
    ax = axs[1, 0]
    band = s.get("placement_band", {})
    for k in kids:
        b = band.get(f"{k}/{ref}/{cond}")
        if b:
            ax.add_patch(plt.Rectangle((xs[k] - 0.42, b["min"]), 0.84, b["max"] - b["min"], fc="0.88", ec="none", zorder=0))
    series = []
    if g3b is not None:
        series.append(("top", "fixed adult helmet, top contact", lambda k: g3b["comparisons"][f"{k}/{o}/{ref}/{cond}/detect"]["delta"]))
        series.append(("counterfactual_x-centred", "helmet scaled with the head",
                       lambda k: g3b["delta_other_placements"][f"{k}/counterfactual_x-centred_vs_adult_counterfactual_x-centred/{ref}"]["delta"]))
    series.append(("gap_matched", "helmet fitted at the adult's gap", lambda k: s["delta_same_rule"][f"{k}/gap_matched/{o}/{ref}/{cond}"]["delta"]))
    eq = dict(marker="D", color="#009E73", mfc="white", label="fitted at the adult's gap, equal OPM site count")
    series.append(("equal", eq["label"], lambda k: s["delta_equal_channels"][f"{k}/gap_matched/{ref}/{cond}"]["delta"]))
    for j, (h, lab, get) in enumerate(series):
        st = eq if h == "equal" else HSTYLE[h]
        dx = (j - (len(series) - 1) / 2) * 0.17
        for i, k in enumerate(kids):
            try:
                e = get(k)
            except KeyError:
                continue
            point(ax, xs[k], e["median"], e.get("ci95"), st, dx, lab if i == 0 else None)
    ax.set_ylabel("Δ = D(smaller head) −\nD(adult, same helmet rule) (dB)")
    ax.set_title("(c)  Change from the adult", loc="left")
    if kids:
        xaxis(ax, kids)

    # (d) within-head contrast
    ax = axs[1, 1]
    for j, h in enumerate(("counterfactual_x-centred", "gap_matched")):
        st = HSTYLE[h]
        for k in keys:
            e = s["within_head"][f"{k}/{h}/{ref}/{cond}"]
            point(ax, xs[k], e["median"], e["ci95"], st, (j - 0.5) * 0.22)
    ax.set_ylabel("D(fixed helmet, top contact) −\nD(fitted helmet), same head (dB)")
    ax.set_title("(d)  Fixed minus fitted helmet in each head", loc="left")
    xaxis(ax, keys)

    handles = [Line2D([], [], marker=HSTYLE[h]["marker"], color=HSTYLE[h]["color"], mfc=HSTYLE[h]["mfc"], ls="none", ms=4.5,
                      label=HSTYLE[h]["label"]) for h in shown]
    handles += [Line2D([], [], marker=eq["marker"], color=eq["color"], mfc=eq["mfc"], ls="none", ms=4.5, label=eq["label"]),
                Patch(fc="0.88", ec="none", label="(c) fixed helmet: range over its source-blind placements"),
                Line2D([], [], color="#0072B2", ls="--", lw=0.8, label=f"(a) adult's gap, laterally centred ({tg['gap_matched']['gap_mm']:.1f} mm)"),
                Line2D([], [], color="#0072B2", ls=":", lw=0.8, label=f"(a) adult's gap at top contact ({tg['gap_matched_top']['gap_mm']:.1f} mm)")]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.205), fontsize=7)
    bind = [style.ANAT_SHORT[k] for k in keys if any(s["helmets"][k][r]["clearance_binding"] for r in RULES)]
    foot = ("D: known-topography detectability of a 10-nAm cortical dipole, dense OPM array against Neuromag (306 channels), sensor "
            "plus brain noise; area-weighted median over the cortex; bars: 95 % intervals, parcels resampled within each head. "
            "Fitted helmet: the adult helmet scaled about the laterally centred head until the median magnetometer-to-scalp gap "
            "equals the adult's, with no magnetometer within 18 mm of the scalp"
            + (f" (* the 18-mm limit is reached first: {', '.join(bind)})" if bind else "") + ". Each head keeps its own OPM arrays. "
            "Grey range in (c): \u0394 at each of the 12 source-blind placements of the fixed helmet (top and back contact, shifts and "
            "rotations; infeasible ones left out), the smaller head and the adult in the same placement, difference of the medians.")
    fig.text(0.5, 0.005, "\n".join(_wrap(foot, 140)), ha="center", va="bottom", fontsize=6.5, color="0.3")
    fig.savefig(path, dpi=style.DPI)
    plt.close(fig)
    return dict(file=path.name, anatomies=keys, comparator=ref, condition=cond, array=o)


def _wrap(text: str, width: int) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width)


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--anatomies", nargs="+", default=list(ANATOMIES), choices=list(ANATOMIES),
                    help="anatomies to compute (the adult is always computed first); default: all nine")
    ap.add_argument("--out", type=Path, default=None, help=f"output directory (default {OUT.relative_to(ROOT)}; state next to it)")
    ap.add_argument("--n-boot", type=int, default=None,
                    help="bootstrap resamples for every interval (tests only; default: configs/g3b_pediatric.toml strata.n_boot "
                         "for the primary summaries, 200 for the others, as G3B)")
    ap.add_argument("--replot", action="store_true", help="redo summaries, files and figure from the stored state")
    args = ap.parse_args()
    t_start = time.time()
    mne.set_log_level("WARNING")
    out_dir = OUT if args.out is None else args.out.resolve()
    state_file = STATE if args.out is None else out_dir / "state.pkl"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    g2cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    keys = ["adult"] + [k for k in ANATOMIES if k in args.anatomies and k != "adult"]
    if args.replot:
        with open(state_file, "rb") as fh:
            state = pickle.load(fh)
    else:
        state = compute(keys, cfg, g2cfg)
        state_file.parent.mkdir(parents=True, exist_ok=True)
        with open(state_file, "wb") as fh:
            pickle.dump(state, fh)
        log(f"state written to {state_file} ({time.time() - t_start:.0f} s)")
    n_boot = dict(primary=cfg["strata"]["n_boot"], secondary=200) if args.n_boot is None else dict(primary=args.n_boot, secondary=args.n_boot)
    g3b = json.loads(G3B_SUMMARY.read_text()) if G3B_SUMMARY.exists() else None
    summary = summarise(state, cfg, g2cfg, g3b, n_boot)
    summary["files"] = dict(targets=write_targets(state, out_dir, cfg["conditions"]["headline"]))
    summary["files"]["figure"] = figure(summary, g3b, out_dir / "Figure_constant_gap.png")
    summary["config"] = dict(g3b=cfg, g2_conventions=dict(band=g2cfg["band"], sensors=g2cfg["sensors"], sources=g2cfg["sources"]),
                             inputs=[str(G3B_SUMMARY.relative_to(ROOT))])
    if args.replot:
        summary["replotted_at_commit"] = io.RUN_COMMIT
    io.write_json(dict(summary, provenance=state["provenance"]), out_dir / "g3b_constant_gap_summary.json")
    log(f"done in {time.time() - t_start:.0f} s (peak RSS {peak_rss_gb():.1f} GB)")


if __name__ == "__main__":
    main()
