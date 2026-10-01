#!/usr/bin/env python3
"""G3B: the fixed adult Neuromag helmet versus a head-adaptive OPM array on smaller heads (NEW).

Anatomies (src/opmsquid/anatomy.py):
  adult      the MNE sample subject (the G2 anatomy and targets)
  school     size-only control: the adult scaled by 85/95 (Jas et al. 2026 Table 1, child 8 y vs adult
             head radius); every vertex homologous to the adult's
  size2yr    size-only control: the adult scaled by the 2-year template's occipitofrontal
             circumference over the adult's
  infant2yr  the 2-year infant template (O'Reilly et al. 2021), native dimensions; an average
             template, not an individual child
  infant18mo, infant12mo  the 18- and 12-month templates of the same series (variability across
             templates; still averages from one database)
Arrays: the same Neuromag helmet (coils 3014/3024, 'accurate' integration, intrinsic noise) at
source-blind placements (src/opmsquid/pediatric.py: centred, top and back contact, bounded
translations and rotations) and a counterfactual helmet scaled with the head (mechanistic
control); OPM arrays refitted to each head under the G2 rules (10-mm cell, 7-mm standoff, 17-mm
packing, clearance and coverage rules, nothing shrunk): the dense array (head-adaptive) and the
matched-site array (the Neuromag sites at the primary placement, projected to the scalp).
Common conventions (G2): analysis band, intrinsic noise, room field, background moment variance per
unit cortical area (calibrated once, on the adult's Neuromag gradiometers), 3-layer BEM
conductivities. Nothing is changed with age except the geometry.
Metric: known-topography detectability d of a 10-nAm cortical-normal dipole, in dB (20 log10 d);
D = dB_OPM - dB_SQUID for each SQUID comparator (mag, grad, combined); Delta = D_child - D_adult on
homologous sources: vertex-wise for the scaled controls; Desikan-Killiany parcels and declared
depth/orientation strata for the template. Summaries are area-weighted medians with intervals from
a bootstrap over parcels. Peak-channel and mean-power SNR (dB) are secondary metrics.
Configuration: configs/g3b_pediatric.toml (and configs/g2_adult.toml). Outputs: results/g3b/.
"""
from __future__ import annotations

import argparse
import csv
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

from opmsquid import (anatomy, background, forward, g2, goldenholz, io, neuromag, noise, opm, paths,  # noqa: E402
                      pediatric as P, plotting)

OUT = ROOT / "results" / "g3b"
STATE = ROOT / "cache" / "g3b" / "state.pkl"
TEMPLATES = {"infant2yr": "ANTS2-0Years3T", "infant18mo": "ANTS18-0Months3T", "infant12mo": "ANTS12-0Months3T"}
ANATOMIES = ("adult", "school", "size2yr") + tuple(TEMPLATES)
CHILDREN = ANATOMIES[1:]
SCALED = ("school", "size2yr")
REFS = ("combined", "grad", "mag")
OTHER_PLACEMENTS = ("centred", "back", "x-centred", "top-18mm", "counterfactual", "counterfactual_x-centred")
PLACEMENT_ORDER = ("centred", "top", "back", "x+5mm", "x-5mm", "y+5mm", "y-5mm", "pitch+10deg", "pitch-10deg", "roll+5deg",
                   "roll-5deg", "x-centred", "top-18mm", "counterfactual", "counterfactual_x-centred")
OPMS = ("opm_dense", "opm_matched")
LABEL = {"adult": "adult", "school": "school-age size (scaled adult)", "size2yr": "2-year size (scaled adult)",
         "infant2yr": "2-year template", "infant18mo": "18-month template", "infant12mo": "12-month template",
         "opm_dense": "OPM dense (refitted)", "opm_matched": "OPM matched",
         "combined": "Neuromag combined", "grad": "Neuromag grad", "mag": "Neuromag mag"}
COLORS = {"adult": "k", "school": "tab:blue", "size2yr": "tab:green", "infant2yr": "tab:red", "infant18mo": "tab:orange",
          "infant12mo": "tab:purple"}
METRICS = ("detect", "peak", "meanpow_db")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------------------------
class Common:
    """Conventions shared by every anatomy (G2's band, noise levels, room field and background scale)."""

    def __init__(self, cfg, g2cfg):
        self.cfg, self.g2cfg = cfg, g2cfg
        b = g2cfg["band"]
        dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
        self.filt = noise.AnalysisFilter(fs=dig["sfreq"], l_freq=b["l_freq_hz"], h_freq=b["h_freq_hz"], order=b["order"])
        self.enbw = self.filt.enbw()
        self.q = g2cfg["sources"]["focal_nAm"] * 1e-9
        self.opm_asd = g2cfg["sensors"]["opm_asd_primary_fT_per_rtHz"] * 1e-15
        self.squid_info = neuromag.load_info("T3")
        self.kinds = neuromag.channel_kinds(self.squid_info)
        self.base = np.array(self.squid_info["dev_head_t"]["trans"])
        bads = g2cfg["sensors"]["bads"]
        good = ~np.isin(self.squid_info.ch_names, bads)
        self.grads = good & (self.kinds == "grad")
        meas = g2.measured_noise(self.squid_info, self.filt, bads)
        self.env = meas["environment"]
        self.target_var = float(np.nanmedian(meas["brain"][self.grads]))
        self.brain_scale = None  # set on the adult (moment variance per unit area, then fixed)
        self.bem = tuple(cfg["anatomy"]["bem_conductivity"])


class Anatomy:
    def __init__(self, key, subject, cortex, src, size, scale_note):
        self.key, self.subject, self.cortex, self.src, self.size, self.scale_note = key, subject, cortex, src, size, scale_note
        self.weights = P.target_areas(cortex, src.target)
        self.nt = len(src.target)
        # the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain) is not cortex:
        # its targets are computed but left out of every G3B summary (they dominate the deepest strata)
        self.cortical = ~np.char.endswith(src.region.astype(str), "unknown")
        self.points = np.concatenate([src.target, src.grid])


def scaled_sources(adult: Anatomy, subject, cortex) -> g2.Sources:
    """The adult's targets and background grid on a scaled copy (those still usable there)."""
    keep_t = cortex.usable[adult.src.target]
    target = adult.src.target[keep_t]
    grid = adult.src.grid[cortex.usable[adult.src.grid]]
    valid = np.flatnonzero(cortex.usable)
    from opmsquid import noisemodel

    grid_area = noisemodel.grid_areas(cortex.rr[grid], cortex.rr[valid], cortex.area[valid])
    depth = anatomy.depth_to_surface(cortex.rr[target], subject.scalp) * 1e3
    orient = anatomy.orientation_angle(cortex.rr[target], cortex.nn[target], subject.inner_skull)
    return g2.Sources(target, grid, grid_area, depth, orient, adult.src.lobe[keep_t], adult.src.region[keep_t],
                      cortex.dist_inner_skull[target] * 1e3)


def load_anatomies(cfg, g2cfg) -> dict:
    seed = g2cfg["sources"]["seed"]
    a_sub = anatomy.load_sample()
    a_cor = anatomy.full_resolution(a_sub)
    adult = Anatomy("adult", a_sub, a_cor, g2.make_sources(a_sub, a_cor, np.random.default_rng(seed)), P.head_size(a_sub),
                    "MNE sample subject (as G2)")
    if list(cfg["anatomy"]["templates"]) != list(TEMPLATES.values()):
        raise ValueError("configs/g3b_pediatric.toml [anatomy] templates must match TEMPLATES")
    templates = {}
    for key, name in TEMPLATES.items():
        t_sub = anatomy.load_template(name)
        t_cor = anatomy.full_resolution(t_sub)
        templates[key] = Anatomy(key, t_sub, t_cor, g2.make_sources(t_sub, t_cor, np.random.default_rng(seed)), P.head_size(t_sub),
                                 f"{LABEL[key]} ({name}), native dimensions")
    infant = templates["infant2yr"]
    out = {"adult": adult}
    ofc = infant.size["ofc_mm"] / adult.size["ofc_mm"]
    for key, f, note in (("school", cfg["anatomy"]["school_age_scale"], "adult x 85/95 (Jas Table 1 child/adult head radius)"),
                         ("size2yr", ofc, "adult x template/adult occipitofrontal circumference")):
        sub = anatomy.scaled(a_sub, f, f"sample_x{f:.4f}")
        cor = anatomy.full_resolution(sub)
        out[key] = Anatomy(key, sub, cor, scaled_sources(adult, sub, cor), P.head_size(sub), f"{note}: {f:.4f}")
    out.update(templates)
    for an in out.values():
        an.ofc_ratio = an.size["ofc_mm"] / adult.size["ofc_mm"]
    return out


# ----------------------------------------------------------------------------------------------
def build_arrays(an: Anatomy, com: Common, cfg) -> tuple[dict, dict]:
    """SQUID arrays (one per placement, plus the counterfactual helmet) and the refitted OPM arrays."""
    pc = cfg["placement"]
    pl = P.placements(com.squid_info, an.subject, com.base, pc["translation_mm"] * 1e-3, pc["pitch_deg"], pc["roll_deg"],
                      pc["clearance_mm"] * 1e-3)
    arrays = {}
    for name, v in pl.items():
        arrays[f"squid:{name}"] = g2.Array("squid", P.with_dev_head(com.squid_info, v["trans"]), com.kinds, None,
                                           {k: x for k, x in v.items() if k not in ("trans", "pose")})
    cf = P.counterfactual_helmet(com.squid_info, an.subject, com.base, an.ofc_ratio)
    arrays["squid:counterfactual"] = g2.Array("squid", P.with_dev_head(cf["info"], cf["trans"]), com.kinds, None,
                                              {k: x for k, x in cf.items() if k not in ("info", "trans")})
    pl["counterfactual"] = {k: x for k, x in cf.items() if k != "info"}
    cfx = P.counterfactual_helmet(com.squid_info, an.subject, pl["x-centred"]["pose"], an.ofc_ratio)
    cfx["rule"] = cfx["rule"].replace("centred placement", "laterally centred placement (x-centred, before the top contact)")
    arrays["squid:counterfactual_x-centred"] = g2.Array("squid", P.with_dev_head(cfx["info"], cfx["trans"]), com.kinds, None,
                                                        {k: x for k, x in cfx.items() if k not in ("info", "trans")})
    pl["counterfactual_x-centred"] = {k: x for k, x in cfx.items() if k != "info"}
    dig = an.subject.digitisation()
    arrays["opm_dense"] = g2.dense_opm(an.subject, dig, "opm_dense")
    arrays["opm_matched"] = g2.matched_opm(an.subject, dig, squid_info=arrays[f"squid:{pc['primary']}"].info)
    return arrays, pl


def gains(an: Anatomy, array: g2.Array, bem: list, points: np.ndarray | None = None, label: str = "") -> np.ndarray:
    pts = an.points if points is None else points
    return forward.chunked_discrete_gain(array.info, an.subject.trans, an.cortex.rr[pts], an.cortex.nn[pts], bem,
                                         coil_def=opm.coil_def_file(), label=label).astype(np.float64)


def evaluate(an: Anatomy, array: g2.Array, g_t, g_g, com: Common, conds, asd=None, bg=1.0, q=None) -> dict:
    nz = g2.array_noise(array, g_g, an.src.grid_area, com.brain_scale * bg, com.env, com.enbw, com.opm_asd if asd is None else asd)
    topo = g_t * (com.q if q is None else q)
    return {(cs, cond): g2.evaluate(topo, nz, m, cond) for cond in conds for cs, m in g2.channel_sets(array).items()}


def d_db(res: dict, opm_name: str, squid_name: str, ref: str, cond: str, metric: str = "detect") -> np.ndarray:
    """D [dB] per source: OPM minus SQUID comparator."""
    o, s = res[opm_name][("opm", cond)][metric], res[squid_name][(ref, cond)][metric]
    return o - s if metric == "meanpow_db" else 20.0 * np.log10(o / s)


# ----------------------------------------------------------------------------------------------
def run_anatomy(an: Anatomy, com: Common, cfg) -> dict:
    t0 = time.time()
    conds = cfg["conditions"]["all"]
    headline = cfg["conditions"]["headline"]
    primary = f"squid:{cfg['placement']['primary']}"
    arrays, pl = build_arrays(an, com, cfg)
    log(f"{an.key}: {an.nt} targets, {len(an.src.grid)} grid sources; OPM dense {arrays['opm_dense'].n}, matched "
        f"{arrays['opm_matched'].n} sites; arrays built ({time.time() - t0:.0f} s)")
    bem3 = an.subject.bem_model(com.bem)
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        g = gains(an, a, bem3)
        G_t[name], G_g[name] = g[:, :an.nt], g[:, an.nt:]
    log(f"{an.key}: 3-layer gains done ({time.time() - t0:.0f} s)")
    if an.key == "adult":  # the one background calibration (G2's rule, adult measured head position)
        unit = background.sensor_covariance(G_g["squid:centred"], background.moment_covariance(an.src.grid_area))
        com.brain_scale = background.calibrate(unit, com.grads, com.target_var)
    res = {name: evaluate(an, a, G_t[name], G_g[name], com, conds) for name, a in arrays.items()}
    # sensitivity: OPM noise, background level (primary placement and OPM arrays only)
    sens = {}
    for asd in com.g2cfg["sensors"]["opm_asd_fT_per_rtHz"]:
        sens[f"opm_asd_{asd:g}fT"] = {n: evaluate(an, arrays[n], G_t[n], G_g[n], com, headline, asd=asd * 1e-15) for n in OPMS}
    for f in cfg["background"]["variance_factors"]:
        sens[f"background_x{f:g}"] = {n: evaluate(an, arrays[n], G_t[n], G_g[n], com, headline, bg=f) for n in OPMS + (primary,)}
    # 1-layer BEM (conductivity-free check), primary arrays
    bem1 = an.subject.bem_model((com.bem[-1],))
    res_bem1 = {}
    for name in OPMS + (primary, "squid:centred"):
        g = gains(an, arrays[name], bem1)
        res_bem1[name] = evaluate(an, arrays[name], g[:, :an.nt], g[:, an.nt:], com, headline)
    log(f"{an.key}: 1-layer BEM done ({time.time() - t0:.0f} s)")
    amp = {n: np.abs(G_t[n][arrays[n].kinds == "mag"] * com.q).max(axis=0) for n in (primary, "opm_dense", "opm_matched")}
    noise_rms = {}  # median per-channel RMS of the brain background and of the intrinsic noise (fT; gradiometers fT/cm)
    for n in (primary, "squid:centred", "opm_dense", "opm_matched"):
        nz = g2.array_noise(arrays[n], G_g[n], an.src.grid_area, com.brain_scale, com.env, com.enbw, com.opm_asd)
        for kind in np.unique(arrays[n].kinds):
            m = arrays[n].kinds == kind
            unit = 1e13 if kind == "grad" else 1e15
            noise_rms[f"{n}/{kind}"] = dict(brain=float(np.median(np.sqrt(np.diag(nz.brain_cov)[m])) * unit),
                                           intrinsic=float(np.median(np.sqrt(nz.intrinsic_var[m])) * unit))
    geometry = sensor_geometry(an, arrays, pl)
    opm_pos = {n: np.array([c["loc"][:3] for c in arrays[n].info["chs"]]) for n in OPMS}  # head frame
    geometry["source_to_sensor_mm"] = source_sensor_distances(an, arrays, (primary, "squid:centred") + OPMS)
    geometry["opm_coverage"] = opm_coverage(an, arrays)
    return dict(arrays={n: dict(n=a.n, meta=a.meta) for n, a in arrays.items()}, placements=pl, res=res, sens=sens, res_bem1=res_bem1,
                opm_pos=opm_pos, noise_rms=noise_rms,
                amp=amp, geometry=geometry, runtime_s=time.time() - t0, _arrays=arrays, _G=(G_t, G_g))


def channel_count_control(an: Anatomy, run: dict, com: Common, cfg, counts) -> dict:
    """The adult's dense array subsampled by farthest-point sampling to each child's dense site
    count (same lead fields, fewer channels): does the smaller channel count alone change D?"""
    a = run["_arrays"]["opm_dense"]
    G_t, G_g = run["_G"]
    pos = np.array([c["loc"][:3] for c in a.info["chs"]])
    out = {}
    for n in counts:
        idx = opm.farthest_point_subset(pos, int(n))
        sub = g2.Array("opm_dense", mne.pick_info(a.info, idx), a.kinds[idx], a.coil_def, dict(n_sites=int(n), subset_of=a.n))
        out[int(n)] = evaluate(an, sub, G_t["opm_dense"][idx], G_g["opm_dense"][idx], com, cfg["conditions"]["headline"])
    return out


def run_patches(an: Anatomy, run: dict, com: Common, cfg, centres: np.ndarray) -> dict:
    """Extended sources on a bounded subset of targets: fixed-total geodesic patches (signed sums of
    cortical-normal dipoles) for the primary SQUID placement and the dense OPM array."""
    t0 = time.time()
    headline = cfg["conditions"]["headline"]
    primary = f"squid:{cfg['placement']['primary']}"
    radii = cfg["sources"]["patch_radii_mm"]
    total = cfg["sources"]["patch_total_nAm"] * 1e-9
    members = {r: goldenholz.geodesic_patches(an.cortex.adjacency, an.src.target[centres], r * 1e-3, an.cortex.usable) for r in radii}
    union = np.unique(np.concatenate([np.concatenate(m) for m in members.values()]))
    col = np.full(an.cortex.n, -1)
    col[union] = np.arange(len(union))
    area = {r: np.array([an.cortex.area[m].sum() for m in members[r]]) for r in radii}
    bem3 = an.subject.bem_model(com.bem)
    G_t, G_g = run["_G"]
    out = dict(centres=centres, area_cm2={f"{r:g}": area[r] * 1e4 for r in radii}, det={})
    for name in (primary, "opm_dense"):
        a = run["_arrays"][name]
        g = gains(an, a, bem3, union)
        for r in radii:
            topo = goldenholz.patch_topographies(g, members[r], col, an.cortex.area) * (total / area[r])[None, :]
            ev = evaluate(an, a, topo, G_g[name], com, headline, q=1.0)
            for (cs, cond), v in ev.items():
                out["det"][(name, cs, cond, r)] = v["detect"]
    log(f"{an.key}: patches ({len(union)} vertices) done ({time.time() - t0:.0f} s)")
    return out


def sensor_geometry(an: Anatomy, arrays: dict, pl: dict) -> dict:
    """Scalp-to-sensor distances: magnetometer coil centres (each placement) and OPM sensing
    centres, plus per helmet region (MNE Vectorview selections) for every placement."""
    from scipy.spatial import cKDTree

    tree = cKDTree(P.scalp_head_frame(an.subject))
    regions = P.helmet_regions(arrays["squid:centred"].info)
    out = {}
    for name, a in arrays.items():
        if name.startswith("squid:"):
            d = P.HelmetFit(a.info, an.subject).distances(a.info["dev_head_t"]["trans"])
            out[name] = dict(median_mm=float(np.median(d) * 1e3), min_mm=float(d.min() * 1e3),
                             regions={k: float(np.median(d[idx]) * 1e3) for k, idx in regions.items()})
        else:
            d = tree.query(np.array([c["loc"][:3] for c in a.info["chs"]]))[0]
            out[name] = dict(median_mm=float(np.median(d) * 1e3), min_mm=float(d.min() * 1e3))
    return out


def source_sensor_distances(an: Anatomy, arrays: dict, names) -> dict:
    """Distance [mm] from every target to the nearest magnetometer coil centre (SQUID) or sensing
    centre (OPM): a geometric quantity kept separate from depth below the scalp and orientation."""
    from scipy.spatial import cKDTree

    hm = an.subject.trans["trans"]  # head -> MRI
    out = {}
    for n in names:
        a = arrays[n]
        pos = np.array([c["loc"][:3] for c, k in zip(a.info["chs"], a.kinds) if k == "mag"])
        if n.startswith("squid:"):
            t = a.info["dev_head_t"]["trans"]
            pos = pos @ t[:3, :3].T + t[:3, 3]
        pos = pos @ hm[:3, :3].T + hm[:3, 3]
        out[n] = cKDTree(pos).query(an.cortex.rr[an.src.target])[0] * 1e3
    return out


def opm_coverage(an: Anatomy, arrays: dict) -> dict:
    """Sites per 100 cm^2 of the scalp above the brow plane (the coverage region of A-OPM-COVER)."""
    rr = P.scalp_head_frame(an.subject)
    tris = an.subject.scalp.tris
    above = opm.above_brow_plane(rr, opm.fiducials_head(an.subject.digitisation()))
    area = 0.5 * np.linalg.norm(np.cross(rr[tris[:, 1]] - rr[tris[:, 0]], rr[tris[:, 2]] - rr[tris[:, 0]]), axis=1)
    cover_cm2 = float(area[above[tris].all(axis=1)].sum() * 1e4)
    return dict(scalp_above_brow_plane_cm2=cover_cm2, **{f"{n}_sites_per_100cm2": 100.0 * arrays[n].n / cover_cm2 for n in OPMS})


# ----------------------------------------------------------------------------------------------
# summaries
def summary_stats(x, w, groups, rng, n_boot):
    ok = np.isfinite(x)
    ci = P.grouped_bootstrap(x[ok], groups[ok], w[ok], rng, n_boot) if n_boot else None
    return dict(median=P.weighted_median(x[ok], w[ok]), ci95=ci, share_positive=float(np.sum(w[ok] * (x[ok] > 0)) / np.sum(w[ok])),
                n=int(ok.sum()))


def strata_rows(x, w, groups, values, edges, rng, n_boot, min_n):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (values >= lo) & (values < hi) & np.isfinite(x)
        row = dict(lo=float(lo), hi=float(hi), n=int(m.sum()))
        if m.sum() >= min_n:
            row.update(median=P.weighted_median(x[m], w[m]),
                       ci95=P.grouped_bootstrap(x[m], groups[m], w[m], rng, n_boot) if n_boot else None)
        rows.append(row)
    return rows


def delta_strata(xc, wc, gc, vc, xa, wa, ga, va, edges, rng, n_boot, min_n):
    """Delta per stratum between two anatomies without vertex correspondence: difference of the
    area-weighted medians, with an interval from resampling parcels within each anatomy."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mc = (vc >= lo) & (vc < hi) & np.isfinite(xc)
        ma = (va >= lo) & (va < hi) & np.isfinite(xa)
        row = dict(lo=float(lo), hi=float(hi), n_child=int(mc.sum()), n_adult=int(ma.sum()))
        if mc.sum() >= min_n and ma.sum() >= min_n:
            row["d_child"] = P.weighted_median(xc[mc], wc[mc])
            row["d_adult"] = P.weighted_median(xa[ma], wa[ma])
            row["delta"] = row["d_child"] - row["d_adult"]
            row["ci95"] = delta_bootstrap(xc[mc], wc[mc], gc[mc], xa[ma], wa[ma], ga[ma], rng, n_boot) if n_boot else None
        rows.append(row)
    return rows


def delta_bootstrap(xc, wc, gc, xa, wa, ga, rng, n_boot):
    def groups(g):
        return [np.flatnonzero(g == lab) for lab in np.unique(g)]

    ic, ia = groups(gc), groups(ga)
    boot = []
    for _ in range(n_boot):
        tc = np.concatenate([ic[j] for j in rng.integers(0, len(ic), len(ic))])
        ta = np.concatenate([ia[j] for j in rng.integers(0, len(ia), len(ia))])
        boot.append(P.weighted_median(xc[tc], wc[tc]) - P.weighted_median(xa[ta], wa[ta]))
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]


def parcel_table(xc, wc, gc, xa, wa, ga, min_n):
    """Delta per Desikan-Killiany parcel (homologous region): area-weighted medians in each anatomy."""
    rows = []
    for lab in sorted(set(gc) | set(ga)):
        if str(lab).endswith("unknown"):
            continue
        mc, ma = (gc == lab) & np.isfinite(xc), (ga == lab) & np.isfinite(xa)
        row = dict(parcel=str(lab), n_child=int(mc.sum()), n_adult=int(ma.sum()), area_adult_cm2=float(wa[ma].sum() * 1e4),
                   area_child_cm2=float(wc[mc].sum() * 1e4))
        if mc.sum() >= min_n and ma.sum() >= min_n:
            row.update(d_child=P.weighted_median(xc[mc], wc[mc]), d_adult=P.weighted_median(xa[ma], wa[ma]))
            row["delta"] = row["d_child"] - row["d_adult"]
        rows.append(row)
    return rows


def compare(an_c: Anatomy, an_a: Anatomy, xc: np.ndarray, xa: np.ndarray, cfg, rng, n_boot: int | None = None) -> dict:
    """D_child, D_adult and Delta for one OPM array / SQUID comparator / condition (intervals from
    ``n_boot`` parcel resamples; 0: medians only)."""
    st = cfg["strata"]
    n_boot = st["n_boot"] if n_boot is None else n_boot
    min_n = st["min_n"]
    dedges, oedges = np.array(st["depth_edges_mm"]), np.array(st["orientation_edges_deg"])
    gc, ga = an_c.src.region.astype(str), an_a.src.region.astype(str)
    out = dict(d_child=summary_stats(xc, an_c.weights, gc, rng, n_boot), d_adult=summary_stats(xa, an_a.weights, ga, rng, n_boot))
    if an_c.key in SCALED:  # vertex-wise homology
        idx = np.searchsorted(an_a.src.target, an_c.src.target)
        assert np.array_equal(an_a.src.target[idx], an_c.src.target)
        dv = xc - xa[idx]
        out["homology"] = "vertex-wise (same cortical vertex)"
        out["delta"] = summary_stats(dv, an_c.weights, gc, rng, n_boot)
        out["delta_by_adult_depth"] = strata_rows(dv, an_c.weights, gc, an_a.src.depth_mm[idx], dedges, rng, n_boot, min_n)
        out["delta_by_orientation"] = strata_rows(dv, an_c.weights, gc, an_c.src.orientation_deg, oedges, rng, n_boot, min_n)
        out["_delta_vertex"] = dv
    else:
        out["homology"] = "Desikan-Killiany parcels and declared depth/orientation strata (no vertex correspondence)"
        parcels = parcel_table(xc, an_c.weights, gc, xa, an_a.weights, ga, min_n)
        ok = [p for p in parcels if "delta" in p]
        dp = np.array([p["delta"] for p in ok])
        wp = np.array([p["area_adult_cm2"] for p in ok])
        boot = [P.weighted_median(dp[j], wp[j]) for j in (rng.integers(0, len(dp), len(dp)) for _ in range(n_boot))]
        ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if n_boot else None
        out["delta"] = dict(median=P.weighted_median(dp, wp), ci95=ci,
                            share_positive=float(np.sum(wp * (dp > 0)) / np.sum(wp)), n_parcels=len(ok),
                            sparse_parcels=[p["parcel"] for p in parcels if "delta" not in p],
                            method="area-weighted median over parcels of (D_child - D_adult), parcel bootstrap")
        out["parcels"] = parcels
    out["delta_by_depth"] = delta_strata(xc, an_c.weights, gc, an_c.src.depth_mm, xa, an_a.weights, ga, an_a.src.depth_mm,
                                         dedges, rng, n_boot, min_n)
    out["delta_by_orientation_strata"] = delta_strata(xc, an_c.weights, gc, an_c.src.orientation_deg, xa, an_a.weights, ga,
                                                      an_a.src.orientation_deg, oedges, rng, n_boot, min_n)
    return out


def usefulness(an: Anatomy, run: dict, cfg, primary: str) -> dict:
    """Where OPM, SQUID, both or neither reach the detectability threshold for a reference moment
    (shares of the usable cortical area represented by the targets, area-weighted)."""
    u = cfg["usefulness"]
    thr = u["detectability_threshold"]
    q0 = 10.0  # nAm, the moment of the evaluated dipoles
    out = {}
    for cond in ("intrinsic+brain", "projected"):
        for ref in REFS:
            d_o = run["res"]["opm_dense"][("opm", cond)]["detect"]
            d_s = run["res"][primary][(ref, cond)]["detect"]
            for qr in u["reference_moments_nAm"]:
                o, s = d_o * qr / q0 >= thr, d_s * qr / q0 >= thr
                w = an.weights * an.cortical / np.sum(an.weights * an.cortical)
                out[f"{ref}/{cond}/{qr:g}nAm"] = dict(both=float(w[o & s].sum()), opm_only=float(w[o & ~s].sum()),
                                                     squid_only=float(w[~o & s].sum()), neither=float(w[~o & ~s].sum()))
            for name, d in (("opm_dense", d_o), (ref, d_s)):
                q5 = np.where(an.cortical, q0 * thr / d, np.nan)  # moment [nAm] at which d reaches the threshold
                out[f"q_threshold_vs_depth/{name}/{cond}"] = strata_rows(q5, an.weights, an.src.region.astype(str), an.src.depth_mm,
                                                                         np.array(cfg["strata"]["depth_edges_mm"]), np.random.default_rng(0),
                                                                         200, cfg["strata"]["min_n"])
    return out


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replot", action="store_true", help="redo summaries, figures and the report from cache/g3b/state.pkl")
    args = ap.parse_args()
    t_start = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    g2cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    anats = load_anatomies(cfg, g2cfg)
    for k, an in anats.items():
        log(f"{k}: {an.scale_note}; OFC {an.size['ofc_mm']:.0f} mm; {an.nt} targets")
    if args.replot:
        with open(STATE, "rb") as fh:
            state = pickle.load(fh)
    else:
        com = Common(cfg, g2cfg)
        runs = {}
        for key in ANATOMIES:  # the adult first: it fixes the background scale
            runs[key] = run_anatomy(anats[key], com, cfg)
        # patch centres: homologous for the adult and the scaled controls, own for the template
        rng = np.random.default_rng(g2cfg["sources"]["seed"] + 1)
        a = anats["adult"]
        common_t = np.intersect1d(np.intersect1d(a.src.target, anats["school"].src.target), anats["size2yr"].src.target)
        chosen = np.sort(rng.choice(common_t, cfg["sources"]["n_patch_centres"], replace=False))
        patches = {}
        for key in ANATOMIES:
            an = anats[key]
            idx = (np.searchsorted(an.src.target, chosen) if key not in TEMPLATES
                   else np.sort(rng.choice(an.nt, cfg["sources"]["n_patch_centres"], replace=False)))
            patches[key] = run_patches(an, runs[key], com, cfg, idx)
        counts = sorted({runs[c]["arrays"]["opm_dense"]["n"] for c in CHILDREN})
        runs["adult"]["channel_count"] = channel_count_control(anats["adult"], runs["adult"], com, cfg, counts)
        for r in runs.values():
            r.pop("_G", None)
            r.pop("_arrays", None)
        state = dict(runs=runs, patches=patches, brain_scale=com.brain_scale, enbw=com.enbw, provenance=dict(
            commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__), runtime_s=time.time() - t_start)
        STATE.parent.mkdir(parents=True, exist_ok=True)
        with open(STATE, "wb") as fh:
            pickle.dump(state, fh)
    summary = summarise(anats, state, cfg)
    if args.replot:  # the computation's commit stays in provenance; record the one that redrew summaries and figures
        summary["replotted_at_commit"] = io.RUN_COMMIT
    figures(anats, state, summary, cfg)
    write_targets_csv(anats, state, cfg)
    write_report(anats, summary, cfg)
    io.write_json(dict(summary, provenance=state["provenance"]), OUT / "g3b_summary.json")
    log(f"done in {time.time() - t_start:.0f} s")


def summarise(anats, state, cfg) -> dict:
    runs, patches = state["runs"], state["patches"]
    primary = f"squid:{cfg['placement']['primary']}"
    headline = cfg["conditions"]["headline"]
    conds = cfg["conditions"]["all"]
    rng = np.random.default_rng(7)
    out = dict(status="NEW (G3B: fixed adult Neuromag helmet vs head-adaptive OPM on smaller heads; size-only controls and "
                      "the 12-, 18- and 24-month templates)", config=cfg, brain_scale=state["brain_scale"], enbw_hz=state["enbw"])
    out["anatomies"] = {k: dict(description=an.subject.description, scale_note=an.scale_note, head_size=an.size,
                                ofc_ratio_to_adult=an.ofc_ratio, n_targets=an.nt, n_background_grid=int(len(an.src.grid)),
                                cortical_area_cm2=float(an.cortex.area[an.cortex.usable].sum() * 1e4),
                                target_depth_mm=dict(zip(("p5", "median", "p95"),
                                                         np.percentile(an.src.depth_mm[an.cortical], [5, 50, 95]).tolist())))
                        for k, an in anats.items()}
    out["arrays"] = {k: {n: dict(n=v["n"], **{kk: vv for kk, vv in v["meta"].items() if not isinstance(vv, (list, np.ndarray, dict))})
                         for n, v in r["arrays"].items()} for k, r in runs.items()}
    out["placements"] = {k: {n: {kk: vv for kk, vv in v.items() if kk not in ("trans", "info", "motion")} for n, v in r["placements"].items()}
                         for k, r in runs.items()}
    out["sensor_distances"] = {k: {n: v for n, v in r["geometry"].items() if n != "source_to_sensor_mm"} for k, r in runs.items()}
    edges = np.array(cfg["strata"]["depth_edges_mm"])
    out["source_to_sensor_mm_by_depth"] = {}
    for k, r in runs.items():
        an = anats[k]
        for n, dist in r["geometry"].get("source_to_sensor_mm", {}).items():
            out["source_to_sensor_mm_by_depth"][f"{k}/{n}"] = [
                dict(lo=float(lo), hi=float(hi), n=int(m.sum()), median=float(np.median(dist[m])) if m.sum() else None)
                for lo, hi in zip(edges[:-1], edges[1:]) for m in [(an.src.depth_mm >= lo) & (an.src.depth_mm < hi) & an.cortical]]
    # D per anatomy (primary placement), all comparators/conditions/metrics
    D = {}
    for k, r in runs.items():
        for o in OPMS:
            for ref in REFS:
                for cond in conds:
                    for metric in METRICS:
                        D[(k, o, primary, ref, cond, metric)] = np.where(anats[k].cortical, d_db(r["res"], o, primary, ref, cond, metric),
                                                                         np.nan)
    out["D_median_dB"] = {f"{k}/{o}/{ref}/{cond}/{metric}": P.weighted_median(v, anats[k].weights)
                          for (k, o, s, ref, cond, metric), v in D.items()}
    # primary comparisons: D_child, D_adult, Delta (dense and matched OPM; every comparator; headline conditions)
    comp = {}
    for c in CHILDREN:
        for o in OPMS:
            for ref in REFS:
                for cond in headline:
                    for metric in METRICS:
                        nb = 0 if metric != "detect" else (cfg["strata"]["n_boot"] if o == "opm_dense" else 200)
                        res = compare(anats[c], anats["adult"], D[(c, o, primary, ref, cond, metric)],
                                      D[("adult", o, primary, ref, cond, metric)], cfg, rng, nb)
                        res.pop("_delta_vertex", None)
                        comp[f"{c}/{o}/{ref}/{cond}/{metric}"] = res
    out["comparisons"] = comp
    # placements: D (dense OPM vs each comparator) per placement and per lobe; Delta at matched rules
    plc = {}
    for k, r in runs.items():
        an = anats[k]
        for name in r["placements"]:
            sq = f"squid:{name}"
            for ref in REFS:
                for cond in headline:
                    x = np.where(an.cortical, d_db(r["res"], "opm_dense", sq, ref, cond), np.nan)
                    plc[f"{k}/{name}/{ref}/{cond}"] = dict(
                        median=P.weighted_median(x, an.weights),
                        by_lobe={lb: P.weighted_median(x[an.src.lobe == lb], an.weights[an.src.lobe == lb]) for lb in plotting.DK_LOBES})
    out["placement_D"] = plc
    for c in CHILDREN:  # the adult at the same rule (its counterfactual helmet has factor 1: the adult itself)
        for name in OTHER_PLACEMENTS:
            for ref in REFS:
                xc = np.where(anats[c].cortical, d_db(runs[c]["res"], "opm_dense", f"squid:{name}", ref, "intrinsic+brain"), np.nan)
                xa = np.where(anats["adult"].cortical, d_db(runs["adult"]["res"], "opm_dense", f"squid:{name}", ref, "intrinsic+brain"),
                              np.nan)
                res = compare(anats[c], anats["adult"], xc, xa, cfg, rng, 200)
                out.setdefault("delta_other_placements", {})[f"{c}/{name}_vs_adult_{name}/{ref}"] = dict(
                    d_child=res["d_child"], d_adult=res["d_adult"], delta=res["delta"])
    # absolute detectability of each system (dB of d for 10 nAm), and its vertex-wise change in the scaled controls
    absd, vw = {}, {}
    for k, r in runs.items():
        an = anats[k]
        for cond in headline:
            for name, cs in (("opm_dense", "opm"), ("opm_matched", "opm")) + tuple((primary, ref) for ref in REFS):
                x = np.where(an.cortical, 20 * np.log10(r["res"][name][(cs, cond)]["detect"]), np.nan)
                absd[f"{k}/{name}/{cs}/{cond}"] = P.weighted_median(x, an.weights)
    for c in SCALED:
        an, a = anats[c], anats["adult"]
        idx = np.searchsorted(a.src.target, an.src.target)
        for cond in headline:
            for name, cs in (("opm_dense", "opm"),) + tuple((primary, ref) for ref in REFS):
                x = 20 * np.log10(runs[c]["res"][name][(cs, cond)]["detect"] / runs["adult"]["res"][name][(cs, cond)]["detect"][idx])
                vw[f"{c}/{name}/{cs}/{cond}"] = P.weighted_median(np.where(an.cortical, x, np.nan), an.weights)
    out["absolute_detectability_dB"] = absd
    out["vertexwise_change_dB"] = vw
    out["template_depth_checks"] = template_depth_checks(anats, D, primary, cfg)
    # channel count: the adult's dense array subsampled to each child's site count
    cc = {}
    ra = runs["adult"]
    for c in CHILDREN:
        n = runs[c]["arrays"]["opm_dense"]["n"]
        res_n = ra.get("channel_count", {}).get(n)
        if res_n is None:
            continue
        for ref in REFS:
            xa_n = np.where(anats["adult"].cortical, 20 * np.log10(res_n[("opm", "intrinsic+brain")]["detect"]
                                                                   / ra["res"][primary][(ref, "intrinsic+brain")]["detect"]), np.nan)
            res = compare(anats[c], anats["adult"], D[(c, "opm_dense", primary, ref, "intrinsic+brain", "detect")], xa_n, cfg, rng, 200)
            cc[f"{c}/{ref}"] = dict(n_sites=int(n), d_adult_subsampled=res["d_adult"], d_child=res["d_child"], delta=res["delta"])
    out["channel_count_control"] = cc
    # sensitivity: OPM noise, background, BEM
    sens = {}
    for k, r in runs.items():
        for key, s in r["sens"].items():
            for o in OPMS:
                for cond in headline:
                    for ref in REFS:
                        sq = s.get(primary, r["res"][primary])
                        x = np.where(anats[k].cortical, 20 * np.log10(s[o][("opm", cond)]["detect"] / sq[(ref, cond)]["detect"]), np.nan)
                        sens[f"{k}/{key}/{o}/{ref}/{cond}"] = P.weighted_median(x, anats[k].weights)
        for o in OPMS:
            for cond in headline:
                for ref in REFS:
                    x = np.where(anats[k].cortical,
                                 20 * np.log10(r["res_bem1"][o][("opm", cond)]["detect"] / r["res_bem1"][primary][(ref, cond)]["detect"]),
                                 np.nan)
                    sens[f"{k}/bem1/{o}/{ref}/{cond}"] = P.weighted_median(x, anats[k].weights)
    for c in CHILDREN:
        for key in [kk for kk in runs[c]["sens"]] + ["bem1"]:
            for o in OPMS:
                for cond in headline:
                    for ref in REFS:
                        sens[f"delta/{c}/{key}/{o}/{ref}/{cond}"] = sens[f"{c}/{key}/{o}/{ref}/{cond}"] - sens[f"adult/{key}/{o}/{ref}/{cond}"]
    out["sensitivity_median_D_dB"] = sens
    # patches
    pt = {}
    for k, p in patches.items():
        for r_ in cfg["sources"]["patch_radii_mm"]:
            pt[f"{k}/area_cm2/{r_:g}mm"] = float(np.median(p["area_cm2"][f"{r_:g}"]))
            for ref in REFS:
                for cond in headline:
                    x = np.where(anats[k].cortical[p["centres"]],
                                 20 * np.log10(p["det"][("opm_dense", "opm", cond, r_)] / p["det"][(primary, ref, cond, r_)]), np.nan)
                    w = anats[k].weights[p["centres"]]
                    pt[f"{k}/{r_:g}mm/{ref}/{cond}"] = P.weighted_median(x, w)
                    focal = D[(k, "opm_dense", primary, ref, cond, "detect")][p["centres"]]
                    pt[f"{k}/focal_same_centres/{ref}/{cond}"] = P.weighted_median(focal, w)
    for c in CHILDREN:
        for r_ in cfg["sources"]["patch_radii_mm"]:
            for ref in REFS:
                for cond in headline:
                    pt[f"delta/{c}/{r_:g}mm/{ref}/{cond}"] = pt[f"{c}/{r_:g}mm/{ref}/{cond}"] - pt[f"adult/{r_:g}mm/{ref}/{cond}"]
    out["patches_median_D_dB"] = pt
    out["usefulness"] = {k: usefulness(anats[k], runs[k], cfg, primary) for k in ANATOMIES}
    out["amplitude_median_fT"] = {k: {n: float(np.median(v[anats[k].cortical]) * 1e15) for n, v in r["amp"].items()}
                                  for k, r in runs.items()}
    out["noise_rms_median"] = {k: r.get("noise_rms") for k, r in runs.items()}
    a = anats["adult"]
    ratio_g2 = runs["adult"]["res"]["opm_dense"][("opm", "intrinsic+brain")]["detect"] / \
        runs["adult"]["res"]["squid:centred"][("combined", "intrinsic+brain")]["detect"]
    out["link_to_g2"] = dict(
        adult_centred_unweighted_all_targets=float(np.median(ratio_g2)),
        adult_centred_unweighted_cortical=float(np.median(ratio_g2[a.cortical])),
        adult_centred_area_weighted_cortical_dB=P.weighted_median(20 * np.log10(ratio_g2[a.cortical]), a.weights[a.cortical]),
        note="G2's headline (1.13x) is the unweighted median over all targets at the measured (= centred) position; G3B "
             "summaries are area-weighted and leave out the medial wall, and its primary placement is top contact")
    out["medial_wall_targets"] = {k: int(np.sum(~an.cortical)) for k, an in anats.items()}
    out["notes"] = [
        "D = 20 log10(d_OPM / d_SQUID) of a 10-nAm cortical-normal dipole (known-topography detectability with the oracle noise "
        "covariance; independent of the moment). Delta = D_child - D_adult. A positive Delta is an increase in relative OPM "
        "performance under these matching assumptions; it does not by itself mean that OPM beats SQUID in the child.",
        "Scaled controls: the adult's vertices, so Delta is vertex-wise. The absolute 4-mm usable-source rule drops "
        f"{anats['adult'].nt - anats['school'].nt} and {anats['adult'].nt - anats['size2yr'].nt} superficial adult targets in the "
        "scaled copies, so D_child and D_adult are medians over slightly different target sets while Delta uses the common vertices. "
        "The templates: no vertex correspondence; Delta is computed per Desikan-Killiany parcel and per declared depth/orientation "
        "stratum from area-weighted medians.",
        "Intervals: bootstrap over parcels of one anatomy (or of each anatomy, for between-anatomy strata); they do not include "
        f"between-subject variability. {('One', 'Two', 'Three', 'Four')[len(TEMPLATES) - 1]} average templates of one database ("
        + ", ".join(LABEL[k] for k in TEMPLATES)
        + ") are not a population: template results are conditional simulations.",
        "Every child array uses the adult's conventions: background moment variance per unit cortical area, room field, "
        "intrinsic noise, sensor sizes and the 3-layer BEM conductivities; only geometry changes. Both systems' detectability "
        "rises in the smaller heads, the OPM's more (absolute detectability table), by different routes: the on-scalp OPM sees "
        "more signal from a cortex that is closer in absolute terms at about the same brain noise, while the SQUIDs' brain noise "
        "falls (the cortex is farther from the fixed helmet and, with the background fixed per unit area, smaller) more than "
        "their signal. The templates' averaged white surfaces are smoother than an individual cortex (usable area "
        + ", ".join(f"{out['anatomies'][k]['cortical_area_cm2']:,.0f}" for k in TEMPLATES) + " cm^2 for the "
        + ", ".join(LABEL[k] for k in TEMPLATES) + f" vs {out['anatomies']['adult']['cortical_area_cm2']:,.0f} cm^2 for the adult), "
        "which lowers their background power and their patch cancellation further; scaling the background variance x0.5 or x2 "
        "leaves D_child almost unchanged.",
        "Placements are chosen from the scalp and helmet geometry only. Under the adult's measured pose a head with other "
        "fiducials need not be centred laterally; 'x-centred' shifts each head along device x to equal left/right median gaps "
        "before the top contact (shift: " + ", ".join(f"{LABEL[k]} {out['placements'][k]['x-centred']['shift_x_mm']:+.1f} mm"
                                                       for k in ANATOMIES)
        + "; negative = to the left), and 'counterfactual_x-centred' scales the helmet about that laterally centred head. The "
        "counterfactual helmet (scaled with the head) is a mechanistic control, not a pediatric SQUID system.",
        "Targets on the medial wall (FreeSurfer 'unknown': the cut through the corpus callosum and midbrain, not cortex) are left "
        "out of every summary; they would otherwise dominate the deepest strata."]
    return out


# ----------------------------------------------------------------------------------------------
def template_depth_checks(anats, D, primary, cfg) -> dict:
    """Templates (no vertex correspondence): how much of the pooled difference of the medians
    (dense OPM vs Neuromag combined, intrinsic + brain) reflects the template's shallower cortex
    (its targets reweighted to the adult's area share per depth stratum), and radial (0-30 deg)
    vs tangential (60-90 deg) sources at matched depth. Deterministic, no intervals."""
    edges = np.array(cfg["strata"]["depth_edges_mm"])
    min_n = cfg["strata"]["min_n"]
    a = anats["adult"]
    key = ("opm_dense", primary, "combined", "intrinsic+brain", "detect")
    xa = D[("adult",) + key]
    da, oa, wa, oka = a.src.depth_mm, a.src.orientation_deg, a.weights, np.isfinite(xa)

    def share(d, w, ok, lo, hi):
        return float(w[ok & (d >= lo) & (d < hi)].sum() / w[ok].sum())

    out = {}
    for k in TEMPLATES:
        an = anats[k]
        xc = D[(k,) + key]
        dc, oc, wc, okc = an.src.depth_mm, an.src.orientation_deg, an.weights, np.isfinite(xc)
        wr = np.zeros_like(wc)
        for lo, hi in zip(edges[:-1], edges[1:]):
            mc, ma = okc & (dc >= lo) & (dc < hi), oka & (da >= lo) & (da < hi)
            if mc.any() and ma.any():
                wr[mc] = wc[mc] * (wa[ma].sum() / wa[oka].sum()) / (wc[mc].sum() / wc[okc].sum())
        rows = []
        for lo, hi in ((0, 15), (15, 25), (25, 40), (40, 90)):
            row = dict(lo=lo, hi=hi)
            for lab, (olo, ohi) in (("radial", (0, 30)), ("tangential", (60, 90.1))):
                mc = okc & (dc >= lo) & (dc < hi) & (oc >= olo) & (oc < ohi)
                ma = oka & (da >= lo) & (da < hi) & (oa >= olo) & (oa < ohi)
                row[lab] = (P.weighted_median(xc[mc], wc[mc]) - P.weighted_median(xa[ma], wa[ma])
                            if mc.sum() >= min_n and ma.sum() >= min_n else None)
            rows.append(row)
        out[k] = dict(pooled_difference_db=P.weighted_median(xc, wc) - P.weighted_median(xa, wa),
                      depth_reweighted_difference_db=P.weighted_median(np.where(wr > 0, xc, np.nan), wr) - P.weighted_median(xa, wa),
                      area_share_10_20mm=share(dc, wc, okc, 10, 20), adult_area_share_10_20mm=share(da, wa, oka, 10, 20),
                      area_share_deeper_50mm=share(dc, wc, okc, 50, np.inf),
                      median_depth_mm=P.weighted_median(dc[okc], wc[okc]), adult_median_depth_mm=P.weighted_median(da[oka], wa[oka]),
                      radial_share=share(oc, wc, okc, 0, 30), orientation_at_matched_depth=rows)
    return out


def write_targets_csv(anats, state, cfg):
    primary = f"squid:{cfg['placement']['primary']}"
    for k, an in anats.items():
        r = state["runs"][k]["res"]
        cols = [(n, cs, cond) for n in (primary, "squid:centred", "squid:counterfactual") + OPMS for (cs, cond) in r[n]]
        with open(OUT / f"g3b_targets_{k}.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["hemi", "vertno", "depth_mm", "orientation_deg", "region", "lobe", "area_mm2"]
                        + [f"detect_{n.replace('squid:', 'squid_')}_{cs}_{cond}" for n, cs, cond in cols])
            for i, v in enumerate(an.src.target):
                wr.writerow([int(an.cortex.hemi[v]), int(an.cortex.vertno[v]), f"{an.src.depth_mm[i]:.2f}", f"{an.src.orientation_deg[i]:.2f}",
                             an.src.region[i], an.src.lobe[i], f"{an.weights[i] * 1e6:.3f}"]
                            + [f"{r[n][(cs, cond)]['detect'][i]:.5f}" for n, cs, cond in cols])


def fmt_ci(s):
    if not s.get("ci95"):
        return f"{s['median']:+.2f} dB"
    return f"{s['median']:+.2f} dB [{s['ci95'][0]:+.2f}, {s['ci95'][1]:+.2f}]"


def fmt_opt(v):
    return "-" if v is None else f"{v:+.2f}"


def write_report(anats, s, cfg):
    primary = cfg["placement"]["primary"]
    L = ["# G3B: fixed adult Neuromag helmet versus head-adaptive OPM on smaller heads (NEW)", "",
         f"Code commit: see `g3b_summary.json` (provenance). Configuration: `configs/g3b_pediatric.toml`. Primary placement: "
         f"`{primary}` (head raised to 20-mm contact). Metric: known-topography detectability in dB; D = OPM - SQUID; "
         "Delta = D_child - D_adult. Intervals: parcel bootstrap (one anatomy each; no between-subject variability).", "",
         "## Anatomies", "", "| anatomy | description | OFC [mm] | breadth x length [mm] | targets | usable cortex [cm2] | OPM dense / matched sites |",
         "|---|---|---|---|---|---|---|"]
    for k, a in s["anatomies"].items():
        arr = s["arrays"][k]
        L.append(f"| {LABEL[k]} | {a['scale_note']} | {a['head_size']['ofc_mm']:.0f} | {a['head_size']['breadth_mm']:.0f} x "
                 f"{a['head_size']['length_mm']:.0f} | {a['n_targets']} | {a['cortical_area_cm2']:.0f} | {arr['opm_dense']['n']} / "
                 f"{arr['opm_matched']['n']} |")
    L += ["", "## Placements in the fixed helmet (magnetometer coil centre to scalp)", "",
          "| anatomy | placement | moved [mm] | min [mm] | median [mm] | feasible |", "|---|---|---|---|---|---|"]
    # x-centred: lateral shift then top contact; counterfactual_x-centred: helmet scaled about the laterally centred head
    for k, pl in s["placements"].items():
        for n, v in pl.items():
            if n in ("centred", "top", "back", "x-centred", "top-18mm", "counterfactual", "counterfactual_x-centred"):
                L.append(f"| {LABEL[k]} | {n} | {v.get('moved_mm', 0.0):.1f} | {v['min_dist_mm']:.1f} | {v['median_dist_mm']:.1f} | "
                         f"{v['feasible']} |")
    g = s["link_to_g2"]
    L += ["", f"Link to G2: at the adult's measured (= centred) position the dense/combined detectability ratio is "
          f"{g['adult_centred_unweighted_all_targets']:.3f}x as an unweighted median over all targets (G2's headline), "
          f"{g['adult_centred_unweighted_cortical']:.3f}x without the medial wall and {g['adult_centred_area_weighted_cortical_dB']:+.2f} dB "
          "area-weighted without it (the G3B convention).", "",
          "## D_child, D_adult and Delta (dense OPM; intrinsic + brain noise; detectability dB)", "",
          "| child anatomy | comparator | D_child | D_adult | Delta | homology |", "|---|---|---|---|---|---|"]
    for c in CHILDREN:
        for ref in REFS:
            r = s["comparisons"][f"{c}/opm_dense/{ref}/intrinsic+brain/detect"]
            L.append(f"| {LABEL[c]} | {LABEL[ref]} | {fmt_ci(r['d_child'])} | {fmt_ci(r['d_adult'])} | {fmt_ci(r['delta'])} | "
                     f"{'vertex' if c in SCALED else 'parcel'} |")
    L += ["", "Projected condition (room-field subspace removed):", "", "| child anatomy | comparator | D_child | D_adult | Delta |",
          "|---|---|---|---|---|"]
    for c in CHILDREN:
        for ref in REFS:
            r = s["comparisons"][f"{c}/opm_dense/{ref}/projected/detect"]
            L.append(f"| {LABEL[c]} | {LABEL[ref]} | {fmt_ci(r['d_child'])} | {fmt_ci(r['d_adult'])} | {fmt_ci(r['delta'])} |")
    L += ["", "Matched-site OPM (coverage control), intrinsic + brain:", "", "| child anatomy | comparator | D_child | D_adult | Delta |",
          "|---|---|---|---|---|"]
    for c in CHILDREN:
        for ref in REFS:
            r = s["comparisons"][f"{c}/opm_matched/{ref}/intrinsic+brain/detect"]
            L.append(f"| {LABEL[c]} | {LABEL[ref]} | {fmt_ci(r['d_child'])} | {fmt_ci(r['d_adult'])} | {fmt_ci(r['delta'])} |")
    L += ["", "## What drives Delta: placement and helmet fit (dense OPM vs Neuromag combined, intrinsic + brain)", "",
          "Each child placement is compared with the adult at the same rule (the adult's counterfactual helmet has factor 1).", "",
          "| child anatomy | top (primary) | centred | x-centred | top-18mm | back | counterfactual | counterfactual, x-centred |",
          "|---|---|---|---|---|---|---|---|"]
    for c in CHILDREN:
        cells = [fmt_ci(s["comparisons"][f"{c}/opm_dense/combined/intrinsic+brain/detect"]["delta"])]
        for name in ("centred", "x-centred", "top-18mm", "back", "counterfactual", "counterfactual_x-centred"):
            cells.append(fmt_ci(s["delta_other_placements"][f"{c}/{name}_vs_adult_{name}/combined"]["delta"]))
        L.append(f"| {LABEL[c]} | " + " | ".join(cells) + " |")
    L += ["", "Counterfactual Delta by comparator (helmet scaled with the head; the dependence on the comparator points to the SQUID "
          "side of the change):", "", "| child anatomy | comparator | counterfactual | counterfactual, x-centred |", "|---|---|---|---|"]
    for c in CHILDREN:
        for ref in REFS:
            L.append(f"| {LABEL[c]} | {LABEL[ref]} | {fmt_ci(s['delta_other_placements'][f'{c}/counterfactual_vs_adult_counterfactual/{ref}']['delta'])} | "
                     f"{fmt_ci(s['delta_other_placements'][f'{c}/counterfactual_x-centred_vs_adult_counterfactual_x-centred/{ref}']['delta'])} |")
    ab, vw = s["absolute_detectability_dB"], s["vertexwise_change_dB"]
    L += ["", "## Absolute detectability (median 20 log10 d of a 10-nAm dipole, intrinsic + brain, primary placement)", "",
          "| anatomy | OPM dense | OPM matched | Neuromag combined | Neuromag grad | Neuromag mag |", "|---|---|---|---|---|---|"]
    for k in ANATOMIES:
        L.append(f"| {LABEL[k]} | {ab[f'{k}/opm_dense/opm/intrinsic+brain']:+.2f} | {ab[f'{k}/opm_matched/opm/intrinsic+brain']:+.2f} | "
                 + " | ".join(f"{ab[f'{k}/squid:{primary}/{ref}/intrinsic+brain']:+.2f}" for ref in REFS) + " |")
    L += ["", "Vertex-wise change from the adult (scaled controls; same vertex): both systems gain, the OPM more.", "",
          "| child anatomy | OPM dense | Neuromag combined | Neuromag grad | Neuromag mag |", "|---|---|---|---|---|"]
    for c in SCALED:
        L.append(f"| {LABEL[c]} | {vw[f'{c}/opm_dense/opm/intrinsic+brain']:+.2f} | "
                 + " | ".join(f"{vw[f'{c}/squid:{primary}/{ref}/intrinsic+brain']:+.2f}" for ref in REFS) + " |")
    cc = s.get("channel_count_control", {})
    if cc:
        L += ["", "## Channel count: the adult's dense array subsampled to each child's site count", "",
              "| child anatomy | sites | D_child | D_adult, subsampled | Delta at equal channel count |", "|---|---|---|---|---|"]
        for c in CHILDREN:
            for ref in ("combined",):
                r = cc.get(f"{c}/{ref}")
                if r:
                    L.append(f"| {LABEL[c]} | {r['n_sites']} | {fmt_ci(r['d_child'])} | {fmt_ci(r['d_adult_subsampled'])} | {fmt_ci(r['delta'])} |")
    L += [""]
    L += ["## Delta by depth stratum (dense OPM vs Neuromag combined, intrinsic + brain)", "",
          "| child anatomy | depth [mm] | n child / adult | D_child | D_adult | Delta [95 % CI] |", "|---|---|---|---|---|---|"]
    for c in CHILDREN:
        for row in s["comparisons"][f"{c}/opm_dense/combined/intrinsic+brain/detect"]["delta_by_depth"]:
            if "delta" in row:
                L.append(f"| {LABEL[c]} | {row['lo']:g}-{row['hi']:g} | {row['n_child']} / {row['n_adult']} | {row['d_child']:+.2f} | "
                         f"{row['d_adult']:+.2f} | {row['delta']:+.2f} [{row['ci95'][0]:+.2f}, {row['ci95'][1]:+.2f}] |")
            else:
                L.append(f"| {LABEL[c]} | {row['lo']:g}-{row['hi']:g} | {row['n_child']} / {row['n_adult']} | sparse | | |")
    L += ["", "Scaled controls, vertex-wise (homologous) Delta by the adult's depth:", "",
          "| child anatomy | adult depth [mm] | n | Delta [95 % CI] |", "|---|---|---|---|"]
    for c in SCALED:
        for row in s["comparisons"][f"{c}/opm_dense/combined/intrinsic+brain/detect"].get("delta_by_adult_depth", []):
            cell = f"{row['median']:+.2f} [{row['ci95'][0]:+.2f}, {row['ci95'][1]:+.2f}]" if row.get("ci95") else "sparse"
            L.append(f"| {LABEL[c]} | {row['lo']:g}-{row['hi']:g} | {row['n']} | {cell} |")
    tdc = s.get("template_depth_checks", {})
    if tdc:
        L += ["", "Templates: the pooled difference of the medians (D_child - D_adult over all targets) and the same with the "
              "template's targets reweighted to the adult's area share per depth stratum; then radial (0-30 deg) and tangential "
              "(60-90 deg) sources at matched depth (difference of the medians; '-': fewer than "
              f"{cfg['strata']['min_n']} targets):", "",
              "| template | pooled | depth-reweighted | area at 10-20 mm (adult) | median depth [mm] (adult) | radial 0-15 / 15-25 / "
              "25-40 / 40-90 mm | tangential 0-15 / 15-25 / 25-40 / 40-90 mm |", "|---|---|---|---|---|---|---|"]
        for k, r in tdc.items():
            rad = " / ".join(fmt_opt(row["radial"]) for row in r["orientation_at_matched_depth"])
            tan = " / ".join(fmt_opt(row["tangential"]) for row in r["orientation_at_matched_depth"])
            L.append(f"| {LABEL[k]} | {r['pooled_difference_db']:+.2f} | {r['depth_reweighted_difference_db']:+.2f} | "
                     f"{r['area_share_10_20mm']:.0%} ({r['adult_area_share_10_20mm']:.0%}) | {r['median_depth_mm']:.1f} "
                     f"({r['adult_median_depth_mm']:.1f}) | {rad} | {tan} |")
    L += ["", "## Delta by orientation stratum (0 deg = radial to the inner skull; dense OPM vs Neuromag combined)", "",
          "| child anatomy | orientation [deg] | n child / adult | D_child | D_adult | Delta [95 % CI] |", "|---|---|---|---|---|---|"]
    for c in CHILDREN:
        for row in s["comparisons"][f"{c}/opm_dense/combined/intrinsic+brain/detect"]["delta_by_orientation_strata"]:
            if "delta" in row:
                L.append(f"| {LABEL[c]} | {row['lo']:g}-{row['hi']:g} | {row['n_child']} / {row['n_adult']} | {row['d_child']:+.2f} | "
                         f"{row['d_adult']:+.2f} | {row['delta']:+.2f} [{row['ci95'][0]:+.2f}, {row['ci95'][1]:+.2f}] |")
            else:
                L.append(f"| {LABEL[c]} | {row['lo']:g}-{row['hi']:g} | {row['n_child']} / {row['n_adult']} | sparse | | |")
    L += ["", "## Placement, counterfactual helmet and sensitivity (median D, dense OPM vs Neuromag combined, intrinsic + brain)", "",
          "| anatomy | " + " | ".join(PLACEMENT_ORDER) + " |", "|---" * (len(PLACEMENT_ORDER) + 1) + "|"]
    for k in ANATOMIES:
        L.append(f"| {LABEL[k]} | " + " | ".join(f"{s['placement_D'][f'{k}/{n}/combined/intrinsic+brain']['median']:+.2f}"
                                                for n in PLACEMENT_ORDER) + " |")
    sens = s["sensitivity_median_D_dB"]
    keys = [f"opm_asd_{a:g}fT" for a in (7, 10, 15, 20, 30)] + ["background_x0.5", "background_x2", "bem1"]
    L += ["", "| anatomy | " + " | ".join(keys) + " |", "|---" * (len(keys) + 1) + "|"]
    for k in ANATOMIES:
        L.append(f"| {LABEL[k]} | " + " | ".join(f"{sens[f'{k}/{kk}/opm_dense/combined/intrinsic+brain']:+.2f}" for kk in keys) + " |")
    L += ["", "Difference of these medians, child minus adult, with the same variant applied to both (a sensitivity of the "
          "medians, not the paired Delta estimator; the background variants scale the adult too):", "",
          "| child anatomy | " + " | ".join(keys) + " |", "|---" * (len(keys) + 1) + "|"]
    for c in CHILDREN:
        L.append(f"| {LABEL[c]} | " + " | ".join(f"{sens[f'delta/{c}/{kk}/opm_dense/combined/intrinsic+brain']:+.2f}" for kk in keys) + " |")
    L += ["", "## Regions that gain or lose with the placement (median D by lobe, dense OPM vs Neuromag combined, intrinsic + brain)", "",
          "| anatomy | placement | " + " | ".join(plotting.DK_LOBES) + " |", "|---|---|" + "---|" * len(plotting.DK_LOBES)]
    for k in ANATOMIES:
        for n in ("centred", "top", "x-centred", "back", "counterfactual"):
            by = s["placement_D"][f"{k}/{n}/combined/intrinsic+brain"]["by_lobe"]
            L.append(f"| {LABEL[k]} | {n} | " + " | ".join(f"{by[lb]:+.2f}" for lb in plotting.DK_LOBES) + " |")
    L += ["", "## Source-to-sensor distance (median, mm) by depth below the scalp", "",
          "| anatomy | sensors | " + " | ".join(f"{a:g}-{b:g}" for a, b in zip(cfg["strata"]["depth_edges_mm"][:-1],
                                                                              cfg["strata"]["depth_edges_mm"][1:])) + " |",
          "|---|---|" + "---|" * (len(cfg["strata"]["depth_edges_mm"]) - 1)]
    for key, rows in s.get("source_to_sensor_mm_by_depth", {}).items():
        k, n = key.split("/", 1)
        if n in (f"squid:{primary}", "opm_dense"):
            L.append(f"| {LABEL[k]} | {n} | " + " | ".join("-" if r["median"] is None else f"{r['median']:.0f}" for r in rows) + " |")
    if any(s.get("noise_rms_median", {}).values()):
        L += ["", "## Noise composition (median per-channel RMS in the band; magnetometers and OPM in fT, gradiometers in fT/cm)", "",
              "| anatomy | channels | brain background | intrinsic |", "|---|---|---|---|"]
        for k, nr in s["noise_rms_median"].items():
            for key, v in (nr or {}).items():
                L.append(f"| {LABEL[k]} | {key} | {v['brain']:.1f} | {v['intrinsic']:.1f} |")
    L += ["", "## Usefulness (dense OPM vs Neuromag combined, intrinsic + brain): share of usable cortical area", "",
          f"A source counts as usable when its detectability reaches {cfg['usefulness']['detectability_threshold']:g} at the reference "
          "moment (an operational choice, not a clinical standard).", "",
          "| anatomy | moment [nAm] | both | OPM only | SQUID only | neither |", "|---|---|---|---|---|---|"]
    for k in ANATOMIES:
        for qr in cfg["usefulness"]["reference_moments_nAm"]:
            u = s["usefulness"][k][f"combined/intrinsic+brain/{qr:g}nAm"]
            L.append(f"| {LABEL[k]} | {qr:g} | {u['both']:.2f} | {u['opm_only']:.2f} | {u['squid_only']:.2f} | {u['neither']:.2f} |")
    L += ["", "## Notes", ""] + [f"- {n}" for n in s["notes"]]
    (OUT / "G3B_report.md").write_text("\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------------
def figures(anats, state, s, cfg):
    runs = state["runs"]
    primary = f"squid:{cfg['placement']['primary']}"
    edges = np.array(cfg["strata"]["depth_edges_mm"])
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 1. D vs depth per anatomy (dense OPM vs each comparator; both headline conditions)
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.2), sharex=True)
    for i, cond in enumerate(cfg["conditions"]["headline"]):
        for j, ref in enumerate(REFS):
            ax = axs[i, j]
            for k in ANATOMIES:
                an = anats[k]
                x = np.where(an.cortical, d_db(runs[k]["res"], "opm_dense", primary, ref, cond), np.nan)
                rows = strata_rows(x, an.weights, an.src.region.astype(str), an.src.depth_mm, edges, np.random.default_rng(0), 200,
                                   cfg["strata"]["min_n"])
                med = np.array([r.get("median", np.nan) for r in rows])
                lo = np.array([r["ci95"][0] if "ci95" in r else np.nan for r in rows])
                hi = np.array([r["ci95"][1] if "ci95" in r else np.nan for r in rows])
                ax.plot(centers, med, "o-", color=COLORS[k], label=LABEL[k], ms=4)
                ax.fill_between(centers, lo, hi, color=COLORS[k], alpha=0.12)
                xm = np.where(an.cortical, d_db(runs[k]["res"], "opm_matched", primary, ref, cond), np.nan)
                rows_m = strata_rows(xm, an.weights, an.src.region.astype(str), an.src.depth_mm, edges, np.random.default_rng(0), 50,
                                     cfg["strata"]["min_n"])
                ax.plot(centers, [r.get("median", np.nan) for r in rows_m], "--", color=COLORS[k], lw=1, alpha=0.7)
            ax.axhline(0, color="0.5", lw=0.8)
            ax.set_title(f"OPM vs {LABEL[ref]}, {cond}", fontsize=9)
            if i == 1:
                ax.set_xlabel("source depth below the scalp [mm] (native)")
            if j == 0:
                ax.set_ylabel("D = OPM - SQUID [dB]")
    axs[0, 0].legend(fontsize=7, title="solid: dense OPM; dashed: matched", title_fontsize=7)
    fig.suptitle("G3B: relative OPM performance vs depth in each anatomy (fixed helmet, top contact; area-weighted medians "
                 "without the medial wall, parcel-bootstrap 95 % bands)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G3B_depth.png", dpi=150)
    plt.close(fig)

    # 2. Delta by depth stratum (child - adult), each comparator, intrinsic+brain
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
    for j, ref in enumerate(REFS):
        ax = axs[j]
        for c in CHILDREN:
            rows = s["comparisons"][f"{c}/opm_dense/{ref}/intrinsic+brain/detect"]["delta_by_depth"]
            x = np.array([0.5 * (r["lo"] + r["hi"]) for r in rows])
            y = np.array([r.get("delta", np.nan) for r in rows])
            lo = np.array([r["ci95"][0] if "ci95" in r else np.nan for r in rows])
            hi = np.array([r["ci95"][1] if "ci95" in r else np.nan for r in rows])
            off = (CHILDREN.index(c) - (len(CHILDREN) - 1) / 2) * 0.5
            ax.errorbar(x + off, y, yerr=[y - lo, hi - y], fmt="o-", color=COLORS[c], label=LABEL[c], ms=4, capsize=2)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlabel("depth stratum [mm] (native)")
        ax.set_title(f"Delta = D_child - D_adult, dense OPM vs {LABEL[ref]}", fontsize=9)
    axs[0].set_ylabel("Delta [dB]")
    axs[0].legend(fontsize=7)
    fig.suptitle("G3B: change of the OPM-SQUID difference from adult to child, per depth stratum (intrinsic + brain noise)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G3B_delta.png", dpi=150)
    plt.close(fig)

    # 3. placements: D per placement and helmet-region gaps
    names = list(PLACEMENT_ORDER)
    fig, axs = plt.subplots(1, 2, figsize=(15, 4.8))
    for i, k in enumerate(ANATOMIES):
        y = [s["placement_D"][f"{k}/{n}/combined/intrinsic+brain"]["median"] for n in names]
        axs[0].plot(np.arange(len(names)) + (i - (len(ANATOMIES) - 1) / 2) * 0.1, y, "o", color=COLORS[k], label=LABEL[k], ms=5)
    for xv in np.arange(len(names) - 1) + 0.5:
        axs[0].axvline(xv, color="0.9", lw=0.6, zorder=0)
    axs[0].set_xticks(range(len(names)))
    axs[0].set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    axs[0].axhline(0, color="0.5", lw=0.8)
    axs[0].set_ylabel("median D, dense OPM vs Neuromag combined [dB]")
    axs[0].legend(fontsize=7)
    axs[0].set_title("helmet placement (source-blind) and the counterfactual helmet", fontsize=9)
    regions = list(s["sensor_distances"]["adult"]["squid:centred"]["regions"])
    w = 0.2
    for i, k in enumerate(ANATOMIES):
        for jn, (n, mk) in enumerate((("squid:centred", "o"), ("squid:top", "s"))):
            v = [s["sensor_distances"][k][n]["regions"][r] for r in regions]
            axs[1].plot(np.arange(len(regions)) + (i - (len(ANATOMIES) - 1) / 2) * w, v, mk, color=COLORS[k],
                        mfc="none" if jn == 0 else COLORS[k],
                        label=f"{LABEL[k]}, {n[6:]}")
    axs[1].set_xticks(range(len(regions)))
    axs[1].set_xticklabels(regions, rotation=30, ha="right", fontsize=8)
    axs[1].set_ylabel("median magnetometer-to-scalp distance [mm]")
    axs[1].legend(fontsize=6, ncol=2)
    axs[1].set_title("regional gaps in the fixed helmet (open: centred, filled: top contact)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G3B_placements.png", dpi=150)
    plt.close(fig)

    # 4. cortex maps: D (dense vs combined, intrinsic+brain, primary) for each anatomy; Delta for the scaled controls
    norm = plt.Normalize(-6, 6)
    a = anats["adult"]
    hemis_a = plotting.inflated_views(a.subject.subjects_dir, "sample", a.subject.src)
    hemis_t = {k: plotting.inflated_views(anats[k].subject.subjects_dir, anats[k].subject.name, anats[k].subject.src)
               for k in TEMPLATES}

    def on_map(an, values, target=None):
        n_lh = int(np.sum(an.cortex.hemi == 0))
        oct_global = np.concatenate([an.subject.src[0]["vertno"], an.subject.src[1]["vertno"] + n_lh])
        v = np.full(len(oct_global), np.nan)
        pos = {g: i for i, g in enumerate(oct_global)}
        for g, x in zip(an.src.target if target is None else target, values):
            v[pos[g]] = x
        v[np.char.endswith(np.array([str(r) for r in region_of(an, oct_global)]), "unknown")] = np.nan
        return v, len(an.subject.src[0]["vertno"])

    rows = []
    for k in ("adult", "school", "size2yr"):
        v, n_lh = on_map(anats[k], d_db(runs[k]["res"], "opm_dense", primary, "combined", "intrinsic+brain"))
        rows.append((f"{LABEL[k]}\nD", v))
    plotting.cortex_map_figure(hemis_a, rows, n_lh, plt.get_cmap("RdBu_r"), norm,
                               "G3B: dense OPM vs Neuromag combined (top contact, intrinsic + brain), D [dB] on the adult's inflated "
                               "cortex (scaled controls are vertex-homologous)", "D [dB]", OUT / "Figure_G3B_maps_scaled.png")
    rows = []
    for c in SCALED:
        an = anats[c]
        idx = np.searchsorted(a.src.target, an.src.target)
        dv = d_db(runs[c]["res"], "opm_dense", primary, "combined", "intrinsic+brain") - \
            d_db(runs["adult"]["res"], "opm_dense", primary, "combined", "intrinsic+brain")[idx]
        v, n_lh = on_map(an, dv)
        rows.append((f"{LABEL[c]}\nDelta", v))
    plotting.cortex_map_figure(hemis_a, rows, n_lh, plt.get_cmap("RdBu_r"), plt.Normalize(-2, 2),
                               "G3B: Delta = D_child - D_adult [dB] at every vertex (dense OPM vs Neuromag combined, top contact, "
                               "intrinsic + brain)", "Delta [dB]", OUT / "Figure_G3B_maps_delta.png")
    for k in TEMPLATES:
        v, n_lh_t = on_map(anats[k], d_db(runs[k]["res"], "opm_dense", primary, "combined", "intrinsic+brain"))
        vm, _ = on_map(anats[k], d_db(runs[k]["res"], "opm_matched", primary, "combined", "intrinsic+brain"))
        plotting.cortex_map_figure(hemis_t[k], [(f"{LABEL[k]}\nD", v), (f"{LABEL[k]}\nD, matched OPM", vm)], n_lh_t,
                                   plt.get_cmap("RdBu_r"), norm,
                                   f"G3B: {LABEL[k]}, OPM vs Neuromag combined (top contact, intrinsic + brain), D [dB]", "dB",
                                   OUT / f"Figure_G3B_maps_{k}.png")

    # 5. usefulness maps at the primary reference moment
    qr = cfg["usefulness"]["primary_reference_nAm"]
    thr = cfg["usefulness"]["detectability_threshold"]
    cmap = matplotlib.colors.ListedColormap(["#efe3c8", "tab:orange", "tab:blue", "tab:purple"])  # neither, SQUID only, OPM only, both
    for k, hem in [("adult", hemis_a)] + [(k, hemis_t[k]) for k in TEMPLATES]:
        an = anats[k]
        rows = []
        for cond in ("intrinsic+brain",):
            for ref in REFS:
                o = runs[k]["res"]["opm_dense"][("opm", cond)]["detect"] * qr / 10.0 >= thr
                sq = runs[k]["res"][primary][(ref, cond)]["detect"] * qr / 10.0 >= thr
                v, n_lh = on_map(an, (o.astype(float) * 2 + sq.astype(float)))
                rows.append((f"{LABEL[k]}\nOPM dense vs\n{LABEL[ref]}", v))
        plotting.cortex_map_figure(hem, rows, n_lh, cmap, matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4),
                                   f"G3B usefulness at {qr:g} nAm (detectability >= {thr:g}; beige neither, orange SQUID only, blue OPM "
                                   f"only, purple both; grey: medial wall)", "0 neither, 1 SQUID only, 2 OPM only, 3 both",
                                   OUT / f"Figure_G3B_usefulness_{k}.png", categorical=True)

    # 6. geometry: sagittal and coronal sections with both helmets and the dense OPM sites
    info = neuromag.load_info("T3")
    mags = P.magnetometer_positions(info)
    fig, axs = plt.subplots(2, len(ANATOMIES), figsize=(4 * len(ANATOMIES), 8))
    for j, k in enumerate(ANATOMIES):
        an = anats[k]
        rr = P.scalp_head_frame(an.subject)
        pls = runs[k]["placements"]
        centre = np.linalg.inv(pls["centred"]["trans"])[:3, 3]
        mags_cf = centre + pls["counterfactual"]["k"] * (mags - centre)
        top = np.linalg.inv(pls["top"]["trans"])  # head -> device at top contact
        opm_dev = runs[k]["opm_pos"]["opm_dense"] @ top[:3, :3].T + top[:3, 3]
        for i, (ax_idx, plane) in enumerate(((0, "sagittal"), (1, "coronal"))):
            ax = axs[i, j]
            for name, col in (("centred", "0.6"), ("top", "k")):
                hd = np.linalg.inv(pls[name]["trans"])
                pts = rr[np.abs(rr[:, ax_idx]) < 0.004] @ hd[:3, :3].T + hd[:3, 3]
                ax.plot(pts[:, 1 - ax_idx], pts[:, 2], ".", ms=0.5, color=col, alpha=0.5)
            for m, mk, col, lab in ((mags, "s", "tab:orange", "Neuromag magnetometers"),
                                    (mags_cf, "x", "tab:green", f"counterfactual helmet (x{pls['counterfactual']['k']:.3f}), "
                                                                 "around the grey head"),
                                    (opm_dev, "o", "tab:blue", "dense OPM sites (top contact)")):
                sel = np.abs(m[:, ax_idx] - centre[ax_idx]) < 0.02
                ax.plot(m[sel, 1 - ax_idx], m[sel, 2], mk, color=col, ms=3.5, mfc="none" if mk == "o" else col, label=lab)
            ax.set_aspect("equal")
            ax.set_title(f"{LABEL[k]}, {plane}", fontsize=8)
            ax.set_xlabel(("device y [m] (anterior +)" if plane == "sagittal" else "device x [m] (right +)"), fontsize=7)
            ax.set_ylabel("device z [m] (up +)", fontsize=7)
            ax.tick_params(labelsize=7)
        axs[0, j].legend(fontsize=6, loc="lower left")
    fig.suptitle("G3B geometry: scalp sections at the centred (grey) and top-contact (black) placements in the fixed helmet; the "
                 "counterfactual helmet is drawn around the centred head; sensors within 20 mm of the section", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G3B_geometry.png", dpi=150)
    plt.close(fig)


def region_of(an, oct_global):
    names = np.concatenate([plotting.read_freesurfer_annot(an.subject.labels / f"{h}.aparc.annot") for h in ("lh", "rh")])
    return names[oct_global]


if __name__ == "__main__":
    main()
