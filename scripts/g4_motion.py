#!/usr/bin/env python3
"""G4 bounded secondary extension: head motion and OPM slippage (NEW; docs/methods.md section 12).

Static fit is studied first (G2, G3B, G4). This script adds two motion mechanisms that the static
maps cannot show (src/opmsquid/motion.py), on the G3B arrays and noise conventions:

A. Sustained displacements (geometry). Neuromag at the primary G3B placement (top contact): the
   head displaced in the fixed helmet. Dense OPM array: the head moves with the array (no change)
   unless the cap slips (rigid rotation about the head origin). Detectability with the displaced
   geometry known (ideal movement compensation, or a known slip) and with the template of the
   reference geometry (mismatched; noise covariance from the displaced data).
B. In-band motion (room-field coupling). The dense OPM array moving rigidly with the head in a static
   residual field: in-band artefact covariance J Sigma J^T per unit field (1 nT uniform, 1 nT/m
   gradient), after no correction, the homogeneous-field projection (3 terms) or the 8-term
   projection, with sensor calibration errors. Its effect on detectability, and the in-band motion
   x field at which the OPM loses 1 or 3 dB or falls to the Neuromag's static detectability. The
   SQUIDs are fixed in the room and have no such term. A time-domain simulation with exact rigid
   motion (slow drift plus in-band jitter) checks the linear model and gives the peak field change.

Not modelled: sensor dynamic range, gain changes with the operating field and cross-axis
projection beyond the declared calibration errors; head-motion statistics of real children;
movement-compensation algorithms (the 'known geometry' rows are their ideal limit).
Configuration: configs/g4_motion.toml. Outputs: results/g4/g4_motion_summary.json,
g4_motion_timecourse_example.csv, G4_motion_report.md, Figure_G4_motion.png.
Usage: g4_motion.py [--quick] [--replot]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g3b_pediatric_helmet as G3  # noqa: E402
from opmsquid import background, g2, io, metrics, motion, noise, opm  # noqa: E402
from opmsquid import pediatric as P  # noqa: E402

OUT = ROOT / "results" / "g4"
CORRECTIONS = ("none", "homogeneous", "homogeneous+gradient")
log = G3.log


# ----------------------------------------------------------------------------------------------
def setup(keys):
    cfg3 = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    cfg2 = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    anats = G3.load_anatomies(cfg3, cfg2)
    com = G3.Common(cfg3, cfg2)
    adult = anats["adult"]
    g = G3.gains(adult, G3.g2.Array("squid", com.squid_info, com.kinds, None, {}), adult.subject.bem_model(com.bem))
    unit = background.sensor_covariance(g[:, adult.nt:], background.moment_covariance(adult.src.grid_area))
    com.brain_scale = background.calibrate(unit, com.grads, com.target_var)  # G2's rule, as G3B and pediatric G4
    return {k: anats[k] for k in keys}, com, cfg3


def reference(an, com, cfg3) -> dict:
    """Primary G3B Neuromag placement and the refitted dense OPM array, with their gains."""
    arrays, pl = G3.build_arrays(an, com, cfg3)
    primary = cfg3["placement"]["primary"]
    bem3 = an.subject.bem_model(com.bem)
    out = dict(squid=arrays[f"squid:{primary}"], opm=arrays["opm_dense"], dev_head=np.array(pl[primary]["trans"]), bem3=bem3, G={})
    for name in ("squid", "opm"):
        g = G3.gains(an, out[name], bem3)
        out["G"][name] = (g[:, :an.nt], g[:, an.nt:])
    return out


def noise_cov(an, array, g_grid, com, cond) -> tuple:
    nz = g2.array_noise(array, g_grid, an.src.grid_area, com.brain_scale, com.env, com.enbw, com.opm_asd)
    return nz, nz.covariance(cond)


def change_db(d, d0, floor):
    r = np.divide(d, d0, out=np.zeros_like(d), where=d0 > 0)
    return np.where(r > 0, 20.0 * np.log10(np.maximum(r, 1e-30)), floor).clip(min=floor)


def summarise_change(an, x, cfg3) -> dict:
    """Area-weighted median over cortical targets, share losing > 3 dB, and medians by depth."""
    m, w = an.cortical, an.weights
    edges = cfg3["strata"]["depth_edges_mm"]
    by_depth = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = m & (an.src.depth_mm >= lo) & (an.src.depth_mm < hi)
        by_depth.append(dict(lo=lo, hi=hi, n=int(k.sum()), median=P.weighted_median(x[k], w[k]) if k.sum() >= cfg3["strata"]["min_n"] else None))
    return dict(median_db=P.weighted_median(x[m], w[m]), share_loss_gt_3db=float(np.sum(w[m] * (x[m] < -3.0)) / np.sum(w[m])),
                by_depth=by_depth)


# ----------------------------------------------------------------------------------------------
# A. sustained displacements
def squid_displacements(cfg, origin) -> list:
    g = cfg["geometry"]
    out = []
    for d in g["translations_mm"]:
        out.append((f"down {d:g} mm", P.translate([0.0, 0.0, -d * 1e-3])))
        for ax, lab in ((0, "x"), (1, "y")):
            for sgn in (1, -1):
                out.append((f"{lab}{'+' if sgn > 0 else '-'}{d:g} mm", P.translate(sgn * d * 1e-3 * np.eye(3)[ax])))
    for deg in g["rotations_deg"]:
        for ax, lab in ((0, "pitch"), (1, "roll"), (2, "yaw")):
            for sgn in (1, -1):
                out.append((f"{lab} {sgn * deg:+g} deg", P.rotate_about(ax, sgn * deg, origin)))
    return out


def slipped_opm(an, array, rot: np.ndarray) -> tuple:
    """The dense array rotated rigidly about the head origin (head frame); sensors that would enter
    the scalp moved outward along their axes to the A-OPM-CLEAR clearances."""
    pos = np.array([c["loc"][:3] for c in array.info["chs"]]) @ rot.T
    ax = np.array([c["loc"][9:12] for c in array.info["chs"]]) @ rot.T
    hm = an.subject.trans["trans"]  # head -> MRI
    skin = next(s for s in an.subject.bem_surfaces if s["id"] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD)
    pos_mri, ax_mri = pos @ hm[:3, :3].T + hm[:3, 3], ax @ hm[:3, :3].T
    pos_mri, extra, _ = opm.resolve_clearance(pos_mri, ax_mri, an.subject.scalp, opm.STANDOFF - 0.001, max_extra=0.03,
                                              model_surface=skin, model_clearance=opm.MODEL_CLEARANCE,
                                              cell_clearance=opm.CELL_CLEARANCE, to_head=np.linalg.inv(hm)[:3, :3])
    mh = np.linalg.inv(hm)
    from scipy.spatial import cKDTree

    scalp = P.scalp_head_frame(an.subject)
    final = pos_mri @ mh[:3, :3].T + mh[:3, 3]
    arr = opm.OPMArray(final, ax, scalp[cKDTree(scalp).query(final)[1]], np.arange(len(pos)), opm.STANDOFF, 0.0, "slipped dense OPM")
    info = opm.make_info(arr)
    return g2.Array("opm_dense", info, array.kinds.copy(), array.coil_def, dict(array.meta, slipped=True)), extra


def part_a(an, ref, com, cfg, cfg3) -> dict:
    cond, floor = cfg["condition"], cfg["geometry"]["floor_db"]
    out = {"squid": {}, "opm": {}}
    fit = P.HelmetFit(com.squid_info, an.subject)
    base = ref["dev_head"]
    origin = np.linalg.inv(base)[:3, 3]
    for sysname in ("squid", "opm"):
        arr0 = ref[sysname]
        gt0, gg0 = ref["G"][sysname]
        _, c0 = noise_cov(an, arr0, gg0, com, cond)
        s0 = gt0 * com.q
        d0 = np.linalg.norm(metrics.whitener(c0).apply(s0), axis=0)
        if sysname == "squid":
            cases = [(name, m, None) for name, m in squid_displacements(cfg, origin)]
        else:
            cases = [(f"slip {lab} {sgn * deg:+g} deg", None, motion.rotation(np.radians(sgn * deg) * np.eye(3)[ax]))
                     for deg in cfg["geometry"]["slip_deg"] for ax, lab in ((0, "x"), (1, "y"), (2, "z")) for sgn in (1, -1)]
        for name, m, rot in cases:
            row = {}
            if sysname == "squid":
                t = P.moved(base, m)
                dist = fit.distances(t)
                row.update(min_dist_mm=float(dist.min() * 1e3), median_dist_mm=float(np.median(dist) * 1e3))
                if dist.min() < P.DEWAR:
                    out[sysname][name] = dict(row, feasible=False)
                    continue
                arr = g2.Array("squid", P.with_dev_head(com.squid_info, t), com.kinds, None, {})
            else:
                arr, extra = slipped_opm(an, arr0, rot)
                p0 = np.array([c["loc"][:3] for c in arr0.info["chs"]])
                p1 = np.array([c["loc"][:3] for c in arr.info["chs"]])
                row.update(median_sensor_shift_mm=float(np.median(np.linalg.norm(p1 - p0, axis=1)) * 1e3),
                           n_lifted=int(np.sum(extra > 0)), max_lift_mm=float(extra.max() * 1e3))
            g = G3.gains(an, arr, ref["bem3"])
            gt, gg = g[:, :an.nt], g[:, an.nt:]
            _, ck = noise_cov(an, arr, gg, com, cond)
            wk = metrics.whitener(ck)
            sk = wk.apply(gt * com.q)
            d_known = np.linalg.norm(sk, axis=0)
            d_mis = motion.mismatched_detectability(wk.apply(s0), sk)
            row.update(feasible=True, known=summarise_change(an, change_db(d_known, d0, floor), cfg3),
                       mismatched=summarise_change(an, change_db(d_mis, d0, floor), cfg3),
                       share_wrong_polarity=float(np.sum(an.weights[an.cortical] * (d_mis[an.cortical] <= 0)) / an.weights[an.cortical].sum()))
            out[sysname][name] = row
        log(f"{an.key}: part A {sysname} done ({len(cases)} cases)")
    return out


# ----------------------------------------------------------------------------------------------
# B. in-band motion in a static residual field
def part_b(an, ref, com, cfg) -> dict:
    """Common random numbers: field draw i is the same for every correction, pivot, calibration level
    and anatomy (seeded by (seed, field, i)); calibration draw i by (seed, level, i), shared by every
    correction and pivot. Thresholds come from the median curve over draws; the 10th-90th
    percentiles of the per-draw thresholds show the draw-to-draw spread."""
    cb = cfg["coupling"]
    cond = cfg["condition"]
    op = ref["opm"]
    gt, gg = ref["G"]["opm"]
    nz, c_ib = noise_cov(an, op, gg, com, cond)
    s_nom = motion.sensors(op.info, op.coil_def)
    basis = motion.external_basis(s_nom, com.env.r0)
    proj = {"none": np.eye(op.n), "homogeneous": motion.projector(basis[:, :3]), "homogeneous+gradient": motion.projector(basis)}
    s = gt * com.q
    m, w = an.cortical, an.weights
    # Neuromag's static detectability at the same condition (no motion term: fixed sensors; no projection charged to it)
    gq, gqg = ref["G"]["squid"]
    _, cq = noise_cov(an, ref["squid"], gqg, com, cond)
    d_sq = np.linalg.norm(metrics.whitener(cq).apply(gq * com.q), axis=0)
    grid = cb["rotation_rms_deg"]
    thetas = np.radians(grid)
    out = dict(rotation_rms_deg=grid, n_draws=cb["n_draws"], cases={}, coupling_fT={})
    d_ib = np.linalg.norm(metrics.whitener(c_ib).apply(s), axis=0)
    fields = {f: [motion.random_field(np.random.default_rng([cb["seed"], fi, i]), 1e-9 if f == "uniform" else 0.0,
                                      1e-9 if f == "gradient" else 0.0) for i in range(cb["n_draws"])]
              for fi, f in enumerate(("uniform", "gradient"))}
    cals = {tuple(c): [motion.with_calibration_errors(s_nom, np.random.default_rng([cb["seed"], 100 + ci, i]), *c)
                       for i in range(cb["n_draws"])] for ci, c in enumerate(cb["calibration"])}

    def med_db(num, den):
        return P.weighted_median(20 * np.log10(num[m] / den[m]), w[m])

    def thresholds(curves_l, curves_d):
        res = {}
        for name, cv, level in [(f"loss_{x:g}dB", curves_l, -x) for x in cb["loss_db"]] + [("D_0dB", curves_d, 0.0)]:
            # draws that never reach the level within the tested rotations are censored (beyond the grid), not dropped
            per = np.array([np.inf if v is None else v for v in (motion.crossing(grid, c, level) for c in cv)])
            q = np.percentile(per, [10, 90], method="inverted_cdf")
            res[name] = dict(median_curve=motion.crossing(grid, np.median(cv, axis=0), level),
                             per_draw_p10_p90=[None if not np.isfinite(x) else float(x) for x in q],
                             draws_reached=int(np.isfinite(per).sum()))
        return res

    for corr, pm in proj.items():
        cs, ss = pm @ c_ib @ pm.T, pm @ s
        ws = metrics.whitener(cs)
        d_static = np.linalg.norm(ws.apply(ss), axis=0)
        out[f"static/{corr}"] = dict(opm_change_db=med_db(d_static, d_ib), D_db=med_db(d_static, d_sq))
        for pivot_name, pivot in (("neck", cb["pivot_head_m"]), ("origin", cb["pivot_sensitivity_m"])):
            for ax_err, gain_err in cb["calibration"]:
                cal = f"tilt{ax_err:g}deg_gain{gain_err * 100:g}pct"
                for field in ("uniform", "gradient"):
                    if pivot_name == "origin" and (field == "uniform" or corr == "none"):
                        continue  # the pivot does not matter for a uniform field; one sensitivity row per correction
                    curves = {"unmodelled": ([], []), "oracle": ([], [])}
                    rms = []
                    for i in range(cb["n_draws"]):
                        b0, grad = fields[field][i]
                        act, gains = cals[(ax_err, gain_err)][i]
                        jac = motion.jacobian(act, b0, grad, com.env.r0, np.asarray(pivot, float), gains)[:, :3]
                        cu = pm @ jac @ jac.T @ pm.T  # per rad^2 of per-axis in-band rotation
                        rms.append(np.median(np.sqrt(np.maximum(np.diag(cu), 0))) * np.radians(1.0) * 1e15)  # fT per deg per unit field
                        for metric, ds in (("unmodelled", motion.unmodelled_detectability(ws, ss, cs, cu, thetas)),
                                           ("oracle", motion.oracle_detectability(ss, cs, cu, thetas))):
                            curves[metric][0].append([med_db(d, d_static) for d in ds])
                            curves[metric][1].append([med_db(d, d_sq) for d in ds])
                    key = f"{corr}/{cal}/{field}/{pivot_name}"
                    out["cases"][key] = {}
                    for metric, (loss, dd) in curves.items():
                        loss, dd = np.array(loss), np.array(dd)
                        out["cases"][key][metric] = dict(
                            opm_change_db_median=np.median(loss, axis=0), opm_change_db_range=[loss.min(axis=0), loss.max(axis=0)],
                            D_db_median=np.median(dd, axis=0), D_db_range=[dd.min(axis=0), dd.max(axis=0)],
                            thresholds_deg_unit_field=thresholds(loss, dd))
                    out["coupling_fT"][key] = dict(median=float(np.median(rms)), range=[float(np.min(rms)), float(np.max(rms))])
        log(f"{an.key}: part B {corr} done")
    return out


# ----------------------------------------------------------------------------------------------
# B'. exact rigid motion over a recording
def trajectory(rng, n, fs, filt, drift_deg, inband_deg) -> np.ndarray:
    """(3, n) rotation vector [rad]: a slow drift (below 0.1 Hz, peak ``drift_deg`` per axis) plus
    jitter whose in-band (``filt``) RMS is ``inband_deg`` per axis."""
    from scipy import signal

    lp = signal.butter(2, 0.1, btype="lowpass", fs=fs, output="sos")
    out = np.empty((3, n))
    for k in range(3):
        d = signal.sosfiltfilt(lp, rng.standard_normal(n))
        d *= drift_deg / np.abs(d).max()
        j = rng.standard_normal(n)
        jb = filt.apply(j)
        j *= inband_deg / np.sqrt(np.mean(jb[n // 10:-n // 10] ** 2))
        out[k] = np.radians(d + j)
    return out


def timecourse(an, ref, com, cfg, rng) -> dict:
    tc, cb = cfg["timecourse"], cfg["coupling"]
    g2cfg = com.g2cfg["band"]
    fs = tc["fs_hz"]
    n = int(tc["duration_s"] * fs)
    filt = noise.AnalysisFilter(fs=fs, l_freq=g2cfg["l_freq_hz"], h_freq=g2cfg["h_freq_hz"], order=g2cfg["order"])
    rot = trajectory(rng, n, fs, filt, tc["drift_deg"], tc["inband_rotation_rms_deg"])
    op = ref["opm"]
    s_nom = motion.sensors(op.info, op.coil_def)
    basis = motion.external_basis(s_nom, com.env.r0)
    proj = {"none": np.eye(op.n), "homogeneous": motion.projector(basis[:, :3]), "homogeneous+gradient": motion.projector(basis)}
    b0, grad = motion.random_field(rng, tc["b0_nT"] * 1e-9, tc["gradient_nT_per_m"] * 1e-9)
    act, gains = motion.with_calibration_errors(s_nom, rng, *tc["calibration"])
    pivot = np.asarray(cb["pivot_head_m"], float)
    y = np.empty((op.n, n))
    for i in range(n):
        y[:, i] = motion.readings(act, b0, grad, com.env.r0, motion.rotation(rot[:, i]), np.zeros(3), pivot, gains)
    raw = y - y[:, :1]
    yf = filt.apply(y - y.mean(axis=1, keepdims=True))
    trim = slice(int(2 * fs), n - int(2 * fs))
    rf = filt.apply(rot)[:, trim]
    sigma = np.cov(rf)  # in-band motion covariance actually present
    jac = motion.jacobian(act, b0, grad, com.env.r0, pivot, gains)[:, :3]
    out = dict(peak_field_change_pT=dict(median=float(np.median(np.abs(raw).max(axis=1)) * 1e12),
                                         max=float(np.abs(raw).max() * 1e12)),
               inband_rotation_rms_deg=np.degrees(np.sqrt(np.diag(sigma))).tolist(), corrections={})
    example = int(np.argmax(np.abs(raw).max(axis=1)))
    t_show = slice(int(10 * fs), int(20 * fs))
    out["_example"] = dict(t=np.arange(n)[t_show] / fs, raw=raw[example, t_show], filtered={})
    for corr, pm in proj.items():
        exact = np.sqrt(np.mean((pm @ yf)[:, trim] ** 2, axis=1))
        lin = np.sqrt(np.maximum(np.diag(pm @ jac @ sigma @ jac.T @ pm.T), 0))
        ratio = exact / np.maximum(lin, 1e-30)
        out["corrections"][corr] = dict(inband_rms_fT_median=float(np.median(exact) * 1e15),
                                        linear_prediction_fT_median=float(np.median(lin) * 1e15),
                                        exact_over_linear_median=float(np.median(ratio)),
                                        exact_over_linear_percentiles={f"p{q:g}": float(np.percentile(ratio, q)) for q in (1, 5, 50, 95, 99)},
                                        exact_over_linear_max_abs_deviation=float(np.max(np.abs(ratio - 1.0))))
        out["_example"]["filtered"][corr] = (pm @ yf)[example, t_show]
    return out


# ----------------------------------------------------------------------------------------------
FIG_SQUID = ("down 2 mm", "down 5 mm", "down 10 mm", "x+5 mm", "x+10 mm", "y-5 mm", "y-10 mm", "pitch +10 deg", "pitch -10 deg",
             "roll +10 deg", "yaw +10 deg")
FIG_OPM = ("slip x +1 deg", "slip x +3 deg", "slip y +1 deg", "slip y +3 deg", "slip z +1 deg", "slip z +3 deg")
CORR_COLOR = {"none": "tab:red", "homogeneous": "tab:blue", "homogeneous+gradient": "tab:green"}
EXAMPLE_COLUMNS = ("t_s", "raw_pT", "none_pT", "homogeneous_pT", "homogeneous+gradient_pT")


def write_example(examples, path, status):
    cols = [examples["t"], examples["raw"] * 1e12] + [examples["filtered"][c] * 1e12 for c in CORRECTIONS]
    with open(path, "w", newline="") as fh:
        io.csv_status(fh, status)
        wr = csv.writer(fh)
        wr.writerow(EXAMPLE_COLUMNS)
        for row in zip(*cols):
            wr.writerow([f"{v:.6g}" for v in row])


def read_example(path) -> dict:
    rows = io.read_csv(path)
    col = {k: np.array([float(r[k]) for r in rows]) for k in EXAMPLE_COLUMNS}
    return dict(t=col["t_s"], raw=col["raw_pT"] * 1e-12, filtered={c: col[f"{c}_pT"] * 1e-12 for c in CORRECTIONS})


def figure(summary, examples, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    keys = summary["anatomies"]
    fig, axs = plt.subplots(2, 3, figsize=(18, 10.5))
    # A: sustained displacement, one panel per anatomy (a subset of the cases; all in the report)
    lo = min([0.0] + [min(r["known"]["median_db"], r["mismatched"]["median_db"]) for k in keys for s in ("squid", "opm")
                      for n, r in summary["geometry"][k][s].items() if r["feasible"]])
    for i, k in enumerate(keys[:3]):
        ax = axs[0, i]
        rows = [("Neuromag", n, summary["geometry"][k]["squid"].get(n)) for n in FIG_SQUID]
        rows += [("OPM", n, summary["geometry"][k]["opm"].get(n)) for n in FIG_OPM]
        rows = [(s, n, r) for s, n, r in rows if r is not None]
        y = np.arange(len(rows))[::-1]
        for yy, (s, n, r) in zip(y, rows):
            col = "tab:orange" if s == "Neuromag" else "tab:blue"
            if not r["feasible"]:
                ax.text(0.25, yy, "infeasible (no room in the helmet)", va="center", ha="right", fontsize=7, color="0.4")
                continue
            ax.plot([r["mismatched"]["median_db"], r["known"]["median_db"]], [yy, yy], color=col, lw=1, alpha=0.6)
            ax.plot(r["mismatched"]["median_db"], yy, "o", color=col, ms=5)
            ax.plot(r["known"]["median_db"], yy, "|", color=col, ms=10, mew=2)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{s} {n}" for s, n, _ in rows], fontsize=7)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set_xlim(lo - 0.3, 0.3)
        ax.set_xlabel("median change in detectability [dB]")
        ax.set_title(f"A. Sustained displacement, {G3.LABEL[k]}", fontsize=10)
        if i == 0:
            ax.legend(handles=[Line2D([], [], marker="o", ls="", color="0.3", label="template of the reference geometry"),
                               Line2D([], [], marker="|", ls="", mew=2, ms=10, color="0.3", label="geometry known"),
                               Line2D([], [], color="tab:orange", label="head moved in the Neuromag helmet"),
                               Line2D([], [], color="tab:blue", label="OPM cap slipped")], fontsize=7, loc="lower left")
    # B: in-band rotation in a unit field, artefact outside the noise model
    th = np.array(summary["coupling"][keys[0]]["rotation_rms_deg"])
    styles = dict(zip(keys, ("-", "--", ":")))
    cal_show = summary["calibration_labels"]
    for j, field in enumerate(("uniform", "gradient")):
        ax = axs[1, j]
        for k in keys:
            for corr in CORRECTIONS:
                for ci, cal in enumerate(cal_show):
                    if (corr == "none") != (ci == 0):  # 'none': calibration does not matter; corrections: with errors only
                        continue
                    c = summary["coupling"][k]["cases"].get(f"{corr}/{cal}/{field}/neck")
                    if c is None:
                        continue
                    ax.plot(th, c["unmodelled"]["opm_change_db_median"], styles[k], color=CORR_COLOR[corr],
                            alpha=1.0 if ci <= 1 else 0.5, lw=1.4)
        ax.set_xscale("log")
        ax.set_ylim(-12, 0.5)
        ax.axhline(-1, color="0.6", lw=0.7)
        ax.axhline(-3, color="0.6", lw=0.7, ls="--")
        unit = "1 nT uniform field" if field == "uniform" else "1 nT/m gradient"
        ax.set_xlabel(f"in-band head rotation [deg RMS per axis] in a {unit}")
        ax.set_ylabel("median OPM detectability change [dB]")
        ax.set_title(f"B. In-band rotation, {field} field (artefact outside the noise model)", fontsize=10)
        if j == 0:
            hs = [Line2D([], [], color=CORR_COLOR[c], label=c + (" (calibration irrelevant: perfect shown)" if c == "none" else ""))
                  for c in CORRECTIONS]
            hs += [Line2D([], [], color="0.3", ls=styles[k], label=G3.LABEL[k]) for k in keys]
            hs += [Line2D([], [], color="0.3", alpha=a, label=lab) for a, lab in
                   ((1.0, "projections: calibration 1 deg / 1 %"), (0.5, "projections: calibration 3 deg / 3 %"))]
            ax.legend(handles=hs, fontsize=7, loc="lower left")
    # B': exact rigid motion
    ax = axs[1, 2]
    ex = examples
    ax.plot(ex["t"], (ex["raw"] - ex["raw"].mean()) * 1e12, color="0.65", lw=0.7, label="raw field change (DC removed)")
    for corr in CORRECTIONS:
        ax.plot(ex["t"], ex["filtered"][corr] * 1e12, color=CORR_COLOR[corr], lw=0.7, label=f"in band, {corr}")
    tc = summary["timecourse"]["config"]
    ax.set_xlabel("time [s]")
    ax.set_ylabel("field at one sensor [pT]")
    ax.set_title(f"B'. Exact rigid motion ({G3.LABEL[tc['anatomy']]}): drift {tc['drift_deg']:g} deg, in band "
                 f"{tc['inband_rotation_rms_deg']:g} deg RMS,\n{tc['b0_nT']:g} nT, {tc['gradient_nT_per_m']:g} nT/m, calibration "
                 f"{tc['calibration'][0]:g} deg / {tc['calibration'][1] * 100:g} %", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("G4 extension: head motion and OPM slippage (bounded; detectability = known-topography matched filter, "
                 "intrinsic + brain noise)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fmt(v, spec="+.1f"):
    return "-" if v is None else format(v, spec)


def report(s) -> str:
    L = ["# G4 extension: head motion and OPM slippage (NEW, bounded secondary analysis)", "",
         "Static fit was studied first (G2, G3B, G4). This analysis adds what the static maps cannot show, on the G3B arrays "
         "and noise conventions (Neuromag at top contact; the refitted dense OPM array; intrinsic + brain noise, the G3B "
         "primary condition; 10-nAm cortical-normal dipoles; area-weighted medians over cortical targets). It bounds two "
         "mechanisms; it does not establish motion robustness. Configuration: `configs/g4_motion.toml`; code: "
         "`scripts/g4_motion.py`, `src/opmsquid/motion.py`.", "",
         "## A. Sustained displacement: the head in the fixed helmet (Neuromag) vs a slipped cap (dense OPM)", "",
         "'known': the displaced geometry is known (template and noise of the displaced geometry: ideal movement "
         "compensation or a known slip). 'mismatched': the template of the reference geometry applied to the displaced data "
         "(noise covariance from those data). A head-mounted array moving with the head has no geometry change; only slips "
         f"appear for it. A wrong-polarity template counts as {s['config']['geometry']['floor_db']:g} dB.", "",
         "| anatomy | system | displacement | known [dB] | mismatched [dB] | share > 3 dB loss (mismatched) | note |",
         "|---|---|---|---|---|---|---|"]
    for k in s["anatomies"]:
        for sysname in ("squid", "opm"):
            for n, r in s["geometry"][k][sysname].items():
                if not r["feasible"]:
                    L.append(f"| {G3.LABEL[k]} | Neuromag | {n} | infeasible | | | nearest magnetometer {r['min_dist_mm']:.1f} mm |")
                    continue
                note = (f"nearest magnetometer {r['min_dist_mm']:.1f} mm" if sysname == "squid"
                        else f"sensors moved {r['median_sensor_shift_mm']:.1f} mm (median); {r['n_lifted']} lifted (max {r['max_lift_mm']:.1f} mm)")
                L.append(f"| {G3.LABEL[k]} | {'Neuromag combined' if sysname == 'squid' else 'OPM dense'} | {n} | "
                         f"{r['known']['median_db']:+.2f} | {r['mismatched']['median_db']:+.2f} | {r['mismatched']['share_loss_gt_3db']:.2f} | {note} |")
    cb = s["config"]["coupling"]
    L += ["", "## B. In-band motion of the head-mounted array in a static residual field", "",
          "The dense OPM array moves rigidly with the head (rotation about a pivot "
          f"{-cb['pivot_head_m'][2] * 1e3:g} mm below the head origin) in a static field: 1 nT uniform, or a 1 nT/m "
          "symmetric traceless gradient (Frobenius norm), random isotropic draws. The artefact covariance is J Sigma J^T for "
          "in-band rotation with the given RMS per axis; it is added to the noise after the correction ('none'; "
          "'homogeneous': 3-term homogeneous-field projection; 'homogeneous+gradient': the 8-term projection of G2), applied "
          "to signal and noise alike. Calibration errors: RMS tilt of each sensitive axis and RMS gain error, unknown to the "
          "analyst. Everything is linear in rotation x field: for a field of B nT the rotation thresholds below divide by B. "
          "Neuromag is fixed in the room and has no such term; its static detectability is the comparator for D.", "",
          "Static reference (no motion): OPM detectability change from the correction alone and D (dense OPM vs Neuromag "
          "combined, both intrinsic + brain):", "", "| anatomy | correction | OPM change [dB] | D [dB] |", "|---|---|---|---|"]
    for k in s["anatomies"]:
        for corr in CORRECTIONS:
            st = s["coupling"][k][f"static/{corr}"]
            L.append(f"| {G3.LABEL[k]} | {corr} | {st['opm_change_db']:+.2f} | {st['D_db']:+.2f} |")
    grid = cb["rotation_rms_deg"]

    def th(r):
        v = r["median_curve"]
        txt = "-" if v is None else (f"<= {grid[0]:g}" if v == grid[0] else f"{v:.3g}")
        q = r.get("per_draw_p10_p90")
        if q:
            lo, hi = (f"> {grid[-1]:g}" if x is None else f"{x:.3g}" for x in q)
            txt += f" [{lo}-{hi}]"
        return txt

    L += ["", f"In-band artefact per channel (median over channels, then over {s['coupling'][s['anatomies'][0]].get('n_draws', '?')} "
          "draws) for 1 deg RMS rotation per axis in the unit field, after each correction [fT]. Then the in-band rotation (deg "
          "RMS per axis, in the unit field) at which the median OPM detectability falls by 1 or 3 dB, or D falls to 0 dB, when "
          "the artefact is not part of the analyst's noise model (matched filter of the static covariance applied to data that "
          f"contain it), from the median curve over draws, with the 10th-90th percentiles of the per-draw thresholds (draws that "
          f"do not reach the level within the tested rotations count as beyond them; '-': not reached up to {max(grid):g} deg). The draws are common random numbers: the same fields and calibration errors in "
          f"every correction, pivot and anatomy. Last column: the loss at {max(grid):g} deg when the artefact is part of the "
          "known noise covariance (the optimal filter nulls its at most 3 spatial patterns: the bound for data-driven nulling or "
          "regression on measured head motion).", "",
          "| anatomy | field | pivot | correction | calibration (tilt, gain) | artefact [fT per deg] | 1 dB | 3 dB | D = 0 | "
          f"oracle loss at {max(grid):g} deg [dB] |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k in s["anatomies"]:
        for key, c in s["coupling"][k]["cases"].items():
            corr, cal, field, pivot = key.split("/")
            be = c["unmodelled"]["thresholds_deg_unit_field"]
            cf = s["coupling"][k]["coupling_fT"][key]
            L.append(f"| {G3.LABEL[k]} | {field} | {pivot} | {corr} | {cal} | {cf['median']:.3g} | {th(be['loss_1dB'])} | "
                     f"{th(be['loss_3dB'])} | {th(be['D_0dB'])} | {c['oracle']['opm_change_db_median'][-1]:+.2f} |")
    tc = s["timecourse"]
    L += ["", f"## B'. Exact rigid motion over {tc['config']['duration_s']:g} s ({G3.LABEL[tc['config']['anatomy']]})", "",
          f"Slow drift up to {tc['config']['drift_deg']:g} deg per axis plus in-band jitter of "
          f"{tc['config']['inband_rotation_rms_deg']:g} deg RMS per axis (measured: "
          + ", ".join(f"{v:.3f}" for v in tc["inband_rotation_rms_deg"]) + f" deg), in B0 = {tc['config']['b0_nT']:g} nT and "
          f"G = {tc['config']['gradient_nT_per_m']:g} nT/m, calibration errors {tc['config']['calibration']} (tilt deg, gain). "
          f"Peak field change at a sensor (drift included): median {tc['peak_field_change_pT']['median']:.0f} pT, maximum "
          f"{tc['peak_field_change_pT']['max']:.0f} pT; this offset moves the sensors' operating point (dynamic range, gain), "
          "which is not modelled beyond the calibration errors.", "",
          "| correction | in-band RMS, exact [fT] | linear prediction [fT] | exact / linear, per channel: median [5th-95th percentile] | largest deviation |",
          "|---|---|---|---|---|"]
    for corr, r in tc["corrections"].items():
        q = r["exact_over_linear_percentiles"]
        L.append(f"| {corr} | {r['inband_rms_fT_median']:.1f} | {r['linear_prediction_fT_median']:.1f} | {q['p50']:.4f} "
                 f"[{q['p5']:.4f}-{q['p95']:.4f}] | {100 * r['exact_over_linear_max_abs_deviation']:.1f} % |")
    L += ["", "## Notes", ""] + [f"- {n}" for n in s["notes"]]
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="adult only, two field draws (smoke test; writes to cache/)")
    ap.add_argument("--replot", action="store_true", help="redo the report and figure from the summary and the example time course")
    args = ap.parse_args()
    t0 = time.time()
    mne.set_log_level("WARNING")
    out_dir = ROOT / "cache" / "g4_motion_quick" if args.quick else OUT
    if args.replot:
        summary = json.loads((out_dir / "g4_motion_summary.json").read_text())
        summary["replotted_at_commit"] = io.RUN_COMMIT
        (out_dir / "G4_motion_report.md").write_text(report(summary))
        figure(summary, read_example(out_dir / "g4_motion_timecourse_example.csv"), out_dir / "Figure_G4_motion.png")
        io.write_json(summary, out_dir / "g4_motion_summary.json")
        log(f"replotted in {time.time() - t0:.0f} s -> {out_dir}")
        return
    cfg = tomllib.loads((ROOT / "configs" / "g4_motion.toml").read_text())
    if args.quick:
        cfg["anatomies"] = ["adult"]
        cfg["coupling"]["n_draws"] = 2
        cfg["geometry"]["translations_mm"] = [5.0]
        cfg["geometry"]["rotations_deg"] = [5.0]
        cfg["geometry"]["slip_deg"] = [3.0]
        cfg["timecourse"]["duration_s"] = 20.0
    keys = cfg["anatomies"]
    anats, com, cfg3 = setup(keys)
    refs, geometry, coupling = {}, {}, {}
    for k in keys:
        refs[k] = reference(anats[k], com, cfg3)
        log(f"{k}: reference arrays ready (Neuromag at {cfg3['placement']['primary']}, OPM dense {refs[k]['opm'].n} sites)")
        geometry[k] = part_a(anats[k], refs[k], com, cfg, cfg3)
        coupling[k] = part_b(anats[k], refs[k], com, cfg)
    tk = cfg["timecourse"]["anatomy"] if cfg["timecourse"]["anatomy"] in keys else keys[0]
    tc = timecourse(anats[tk], refs[tk], com, cfg, np.random.default_rng(cfg["timecourse"]["seed"]))  # its own stream
    examples = tc.pop("_example")
    tc["config"] = dict(cfg["timecourse"], anatomy=tk)
    summary = dict(
        status="NEW (G4 bounded secondary extension: head motion and OPM slippage; not a motion-robustness result)",
        config=cfg, anatomies=keys, calibration_labels=[f"tilt{a:g}deg_gain{g * 100:g}pct" for a, g in cfg["coupling"]["calibration"]],
        geometry=geometry, coupling=coupling, timecourse=tc,
        notes=[
            "A head-mounted array moving rigidly with the head keeps its sensor-to-head geometry, so sustained head displacement "
            "changes only the SQUID geometry; the OPM's geometry changes only if the cap slips. The 'known' rows are the ideal "
            "limit of movement compensation (continuous head-position tracking for the SQUID, a measured slip for the OPM).",
            "In a perfectly calibrated array any rigid motion changes the readings by a uniform field plus a symmetric traceless "
            "gradient in the head frame, which the 8-term projection removes exactly; in a uniform field the change is uniform and "
            "the homogeneous projection removes it exactly; translation in a gradient is uniform, rotation in a gradient is not "
            "(tests/test_motion.py, finite rotations included). Calibration errors leave residuals proportional to the field change.",
            "D compares the OPM after its correction with Neuromag without any projection or SSS (none is charged to it), which is "
            "conservative for the OPM. The gradient is defined about x_ref = r0 of G2 (0, 0, 40 mm, head frame), which sets the "
            "uniform part of a translated gradient.",
            "A slipped rigid cap must lift where the head is in its way: sensors that would enter the scalp are moved out along "
            "their axes without the 5-mm limit of A-OPM-CLEAR (a flexible cap or sliding holders), so the slip rows include "
            "the lifts (count and largest lift per row in the table).",
            "The static room field is treated as in G2 for both systems (removed exactly by the projection in the projected "
            "condition); the same calibration errors would also leave part of it, for both systems, which is not modelled here.",
            "Not modelled: sensor dynamic range, gain change with the operating field and cross-axis projection beyond the "
            "declared calibration errors; head-motion statistics of real children; specific movement-compensation algorithms; "
            "field changes from moving magnetic material.",
            "Field strengths, motion amplitudes, pivot and calibration errors are declared sweeps (A-MOT-*), not measurements of "
            "a particular room or device; results scale linearly with rotation x field."])
    out_dir.mkdir(parents=True, exist_ok=True)
    io.write_json(summary, out_dir / "g4_motion_summary.json")
    write_example(examples, out_dir / "g4_motion_timecourse_example.csv", summary["status"])
    (out_dir / "G4_motion_report.md").write_text(report(summary))
    figure(summary, examples, out_dir / "Figure_G4_motion.png")
    log(f"done in {time.time() - t0:.0f} s -> {out_dir}")


if __name__ == "__main__":
    main()
