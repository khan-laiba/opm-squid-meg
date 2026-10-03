#!/usr/bin/env python3
"""A-BEM-CHILD checked on the adult: how much does the modelled skull change the G2 headline?

The school-aged children's watershed segmentations (OpenNeuro ds005234) put the inner skull just
below the scalp (a median 0.8-2.5 mm over the upper head, against 9.7 mm in the adult and 5.6-8.5
mm in the infant templates: scripts/study_school_anatomy.py), so their skull is modelled
(``anatomy.model_skull``): the inner skull is moved inward to a set depth below the MRI scalp where
it is shallower, no vertex closer than 2 mm to a white-surface vertex, and
the outer skull is put halfway to the scalp. Here the same procedure is applied to the adult, whose
segmented skull is known: its inner skull is first degraded the way the children's failed (the
part above its centroid moved to 1.5 mm below the scalp, blended in over 15 mm), then modelled at
depths of 6, 8 (the children's value), 9.7 (the adult's own median) and 10 mm. For each variant:
the error of the modelled inner skull against the segmented one, and the G2 headline (oracle
detectability of the dense and matched OPM arrays relative to Neuromag combined, intrinsic + brain
and projected; G2's arrays, measured-noise calibration refitted per head model, room field and
parcel bootstrap; point gains) on the targets and background grid that are at least 4 mm inside
every variant's inner skull.

Output: results/g3b/child_bem_validation.json.
"""
from __future__ import annotations

import dataclasses
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import mne  # noqa: E402
import numpy as np  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from mne.surface import complete_surface_info  # noqa: E402

import g2_adult_comparison as G2  # noqa: E402
from opmsquid import anatomy, g2, io  # noqa: E402

OUT = ROOT / "results" / "g3b"
CONDS = ("intrinsic+brain", "projected")
DEPTHS = (0.006, 0.008, 0.0097, 0.010)
DEPTH_EDGES = (10.0, 20.0, 30.0, 45.0, 70.0)
FAILED_DEPTH = 0.0015  # the degraded inner skull's depth below the scalp (children: a median 0.8-2.5 mm)
BLEND = 0.015  # [m] above the inner skull's centroid over which the degradation is blended in


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def surf_of(subject, sid):
    return next(s for s in subject.bem_surfaces if s["id"] == sid)


def degrade(inner: dict, scalp: anatomy.Surface) -> tuple[dict, np.ndarray]:
    """The inner skull's upper part moved outward to FAILED_DEPTH below the scalp (blended in)."""
    rr, tris = np.asarray(inner["rr"], float), np.asarray(inner["tris"])
    nn = anatomy._outward(inner).nn
    d_s = anatomy.MeshDistance(dict(rr=scalp.rr, tris=scalp.tris)).unsigned(rr)
    w = np.clip((rr[:, 2] - rr[:, 2].mean()) / BLEND, 0.0, 1.0)
    new = rr + (w * np.maximum(d_s - FAILED_DEPTH, 0.0))[:, None] * nn
    if anatomy.crossing_edges(new, tris):
        raise RuntimeError("the degraded inner skull folds")
    out = complete_surface_info(dict(id=FIFF.FIFFV_BEM_SURF_ID_BRAIN, coord_frame=inner["coord_frame"], rr=new, tris=tris.copy(),
                                     np=len(new), ntri=len(tris)), copy=False, verbose=False)
    return out, w > 0


def restrict(st, keep_t: np.ndarray, keep_g: np.ndarray) -> None:
    s = st.src
    st.src = g2.Sources(s.target[keep_t], s.grid[keep_g], s.grid_area[keep_g], s.depth_mm[keep_t], s.orientation_deg[keep_t],
                        s.lobe[keep_t], s.region[keep_t], s.dist_inner_skull_mm[keep_t])
    st.points = np.concatenate([st.src.target, st.src.grid])
    st.nt = len(st.src.target)


def headline(st, subject, arrays: dict, meas, cfg) -> dict:
    """Oracle comparisons with this head model; the background is calibrated on its Neuromag gradiometers (G2's rule)."""
    st.subject = subject
    squid = arrays["squid"]
    good = ~np.isin(squid.info.ch_names, cfg["sensors"]["bads"])
    grads = good & (squid.kinds == "grad")
    G = {n: st.gains(a, fullres=False) for n, a in arrays.items()}
    bs = G2.calibrate(st, squid, G["squid"][1], float(np.nanmedian(meas["brain"][grads])), grads)
    res = {n: G2.evaluate_all(G[n][0] * st.q, a, st.noise(a, G[n][1], bs, meas["environment"]), CONDS) for n, a in arrays.items()}
    out = dict(brain_scale=bs)
    depth = st.src.depth_mm
    for n in arrays:
        if n == "squid":
            continue
        for cond in CONDS:
            o, s = G2.detect_of(res, n, None, cond), G2.detect_of(res, "squid", "combined", cond)
            c = G2.compare(o, s, st.rng, 1000)
            lr = np.log2(o / s)
            out[f"{n}/combined/{cond}"] = dict(ratio=float(2 ** c["median_log2"]), ci95=[float(2 ** x) for x in c["ci95"]],
                                               by_depth=[dict(lo=lo, hi=hi, n=int(np.sum((depth >= lo) & (depth < hi))),
                                                              ratio=float(2 ** np.median(lr[(depth >= lo) & (depth < hi)])))
                                                         for lo, hi in zip(DEPTH_EDGES[:-1], DEPTH_EDGES[1:])])
    out["squid_combined_detectability_median"] = float(np.median(G2.detect_of(res, "squid", "combined", "intrinsic+brain")))
    out["opm_dense_detectability_median"] = float(np.median(G2.detect_of(res, "opm_dense", None, "intrinsic+brain")))
    return out


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    st = G2.Study(cfg)
    G2.GROUPS = st.src.region
    sub = st.subject
    cortex = st.cortex
    head = surf_of(sub, FIFF.FIFFV_BEM_SURF_ID_HEAD)
    inner = surf_of(sub, FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    skull = surf_of(sub, FIFF.FIFFV_BEM_SURF_ID_SKULL)
    failed, degraded = degrade(inner, sub.scalp)
    to_true = anatomy.MeshDistance(inner)
    to_scalp = anatomy.MeshDistance(dict(rr=sub.scalp.rr, tris=sub.scalp.tris))
    variants = {"segmented": (inner, skull)}
    report = dict(status="NEW check: A-BEM-CHILD (the modelled skull of the school-aged children) applied to the adult",
                  failed_depth_mm=FAILED_DEPTH * 1e3, degraded_vertices=int(degraded.sum()), n_vertices=len(degraded),
                  segmented_scalp_depth_mm=[float(x) for x in np.percentile(to_scalp.unsigned(inner["rr"][degraded]) * 1e3, [10, 50, 90])],
                  failed_scalp_depth_mm=[float(x) for x in np.percentile(to_scalp.unsigned(failed["rr"][degraded]) * 1e3, [10, 50, 90])],
                  models={}, variants={})
    for depth in DEPTHS:
        new, outer, rep = anatomy.model_skull(failed, sub.scalp, cortex.rr[cortex.valid], depth=depth)
        name = f"modelled_{depth * 1e3:g}mm"
        err = to_true.signed(new["rr"][degraded]) * 1e3  # positive: outside the segmented inner skull
        report["models"][name] = dict(rep, inner_skull_error_mm_p10_p50_p90=[float(x) for x in np.percentile(err, [10, 50, 90])],
                                      inner_skull_abs_error_median_mm=float(np.median(np.abs(err))))
        variants[name] = (new, outer)
        log(f"{name}: inner-skull error p10/50/90 {np.round(np.percentile(err, [10, 50, 90]), 1)} mm")
    # targets and grid points at least 4 mm inside every variant's inner skull
    keep = np.ones(len(st.points), bool)
    for name, (inn, _) in variants.items():
        keep &= anatomy.MeshDistance(inn).signed(cortex.rr[st.points]) <= -anatomy.MIN_BEM_DISTANCE
    nt = st.nt
    report["targets_kept"], report["grid_kept"] = int(keep[:nt].sum()), int(keep[nt:].sum())
    restrict(st, keep[:nt], keep[nt:])
    G2.GROUPS = st.src.region
    squid_info = g2.neuromag.load_info("T3")
    meas = g2.measured_noise(squid_info, st.filt, cfg["sensors"]["bads"])
    arrays = g2.build_arrays(sub, st.dig)
    arrays = {n: arrays[n] for n in ("squid", "opm_matched", "opm_dense")}
    for name, (inn, outer) in variants.items():
        subject = dataclasses.replace(sub, bem_surfaces=[head, outer, inn])
        report["variants"][name] = headline(st, subject, arrays, meas, cfg)
        v = report["variants"][name]
        log(f"{name}: " + "; ".join(f"{k} {x['ratio']:.3f}" for k, x in v.items() if isinstance(x, dict) and "ratio" in x))
    report["runtime_s"] = time.time() - t0
    io.write_json(report, OUT / "child_bem_validation.json")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
