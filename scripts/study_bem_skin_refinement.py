#!/usr/bin/env python3
"""Does the G2 headline depend on the head-surface discretization? (pre-freeze review F2 follow-up)

scripts/study_opm_near_mesh.py shows that, on the sample subject's 3-layer BEM, the field at OPM
cell integration points changes by several percent (95th percentile 8-12 % of the array field
scale at 3-4 mm) when only the head surface is refined from 5,120 to 20,480 triangles, while an
exact sphere test (scripts/study_bem_sphere_accuracy.py) shows the same BEM code is accurate near
a regular sphere. This script measures what that means for the G2 comparison: the median
detectability ratios OPM / Neuromag (intrinsic+brain and projected) on the G2 convergence subset
(the same 1,000 targets, seed 2) and the G2 background grid, recomputed with
  bem3_5120        the G2 primary model (3 layers, 5,120 triangles per surface)
  bem3_skin20480   the same with the head surface subdivided once (flat midpoints)
  bem1_5120        inner skull only (the standard MEG model; no boundary near the sensors)
Each model gets its own background calibration on the Neuromag gradiometers, as in G2.
Output: results/g2/bem_skin_refinement.json.
"""
from __future__ import annotations

import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g2_adult_comparison as G  # noqa: E402
from opmsquid import background, forward, g2, io, opm  # noqa: E402

OUT = ROOT / "results" / "g2"


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    st = G.Study(cfg)
    headline = cfg["conditions"]["headline"]
    arrays = {k: v for k, v in g2.build_arrays(st.subject, st.dig).items() if k in ("squid", "opm_matched", "opm_dense")}
    squid = arrays["squid"]
    bads = cfg["sensors"]["bads"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads = good & (squid.kinds == "grad")
    meas = g2.measured_noise(squid.info, st.filt, bads)
    env = meas["environment"]
    target_var = float(np.nanmedian(meas["brain"][grads]))
    rng = np.random.default_rng(2)  # the G2 convergence subset
    sub = np.sort(rng.choice(st.nt, cfg["convergence"]["subset_targets"], replace=False))
    pts = np.concatenate([st.src.target[sub], st.src.grid])
    ns = len(sub)
    models = {"bem3_5120": st.subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=0),
              "bem3_skin20480": st.subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=1), "bem1_5120": st.subject.bem_model((0.3,))}
    out = dict(status="NEW check (review F2 follow-up)", n_targets=int(ns), n_grid=int(len(st.src.grid)), models={})
    gains = {}
    for label, bem in models.items():
        t1 = time.time()
        gt, gg = {}, {}
        for name, a in arrays.items():
            g = forward.chunked_discrete_gain(a.info, st.subject.trans, st.cortex.rr[pts], st.cortex.nn[pts], bem,
                                              coil_def=opm.coil_def_file()).astype(np.float64)
            gt[name], gg[name] = g[:, :ns], g[:, ns:]
        gains[label] = gt
        bs = background.calibrate(st.unit_brain(gg["squid"]), grads, target_var)
        res = {name: G.evaluate_all(gt[name] * st.q, a, st.noise(a, gg[name], bs, env), headline) for name, a in arrays.items()}
        med = {f"{a}/{ref}/{cond}": float(np.median(np.log2(G.detect_of(res, a, None, cond) / G.detect_of(res, "squid", ref, cond))))
               for a in ("opm_matched", "opm_dense") for ref in ("combined", "grad", "mag") for cond in headline}
        out["models"][label] = dict(median_log2=med, ratio={k: float(2 ** v) for k, v in med.items()}, brain_scale=float(bs))
        print(f"{label} ({time.time() - t1:.0f} s): " + ", ".join(f"{k} {2 ** v:.3f}" for k, v in med.items() if "/combined/" in k),
              flush=True)
    ref = gains["bem3_5120"]
    for label in ("bem3_skin20480", "bem1_5120"):
        out["models"][label]["gain_rel_diff_vs_bem3_5120"] = {
            n: dict(median=float(np.median(np.linalg.norm(gains[label][n] - ref[n], axis=0) / np.linalg.norm(ref[n], axis=0))),
                    p95=float(np.percentile(np.linalg.norm(gains[label][n] - ref[n], axis=0) / np.linalg.norm(ref[n], axis=0), 95)))
            for n in ref}
        print(label, "gain difference vs bem3_5120:", {n: round(v["median"], 4) for n, v in out["models"][label]["gain_rel_diff_vs_bem3_5120"].items()})
    out["runtime_s"] = time.time() - t0
    io.write_json(out, OUT / "bem_skin_refinement.json")
    print(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
