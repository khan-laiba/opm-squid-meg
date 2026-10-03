#!/usr/bin/env python3
"""G2 headline under the v3 and v4 head models: what the equal-standoff fix changed (A-BEM-CONFORM).

v4 puts the sample subject's BEM head surface on its MRI scalp. Against v3 this changes three
things at once: the forward model's outer boundary, the OPM sites the whole-cell clearance rule
admits (the dense array reaches the lower occipital scalp again) and the OPM sensitive axes (the
conformed surface's own triangle normals instead of the stored outer skin's, which lie close to
radial). This check separates them on the G2 headline: oracle known-topography detectability of
the dense and matched OPM arrays relative to Neuromag combined, intrinsic + brain and projected,
with G2's targets, background grid, measured-noise calibration (refitted for each head model),
room field and parcel bootstrap; point gains (no full-resolution matrices).

  v3                     stored outer skin; arrays built on it (205 dense, 95 matched sites)
  v3_arrays_v4_bem       the same v3 arrays, the conformed surface in the forward model
  v4_sites_v3_axes       conformed surface, v4 sites, axes from the stored outer skin's normals
                         (averaged within 15 mm at the site's nearest scalp point, as in v3; the
                         cells turn with the axis, the sensing centres stay)
  v4                     conformed surface, v4 arrays (208 dense, 98 matched sites)
  v4_dense_above_v3_low  the v4 dense array without the sites below v3's lowest dense site
                         (the lower occipital scalp that v3's outer skin excluded)

Output: results/g2/head_surface_effect.json.
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
from scipy.spatial import cKDTree  # noqa: E402

import g2_adult_comparison as G2  # noqa: E402
from opmsquid import anatomy, g2, io, neuromag, opm  # noqa: E402

OUT = ROOT / "results" / "g2"
CONDS = ("intrinsic+brain", "projected")
DEPTH_EDGES = (10.0, 20.0, 30.0, 45.0, 70.0)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stored_subject(sub: anatomy.Subject) -> anatomy.Subject:
    """The sample subject with its stored BEM surfaces (v1-v3: the outer skin as segmented)."""
    surfs = mne.read_bem_surfaces(sub.subjects_dir / "sample" / "bem" / "sample-5120-5120-5120-bem.fif", verbose=False)
    return dataclasses.replace(sub, bem_surfaces=surfs, head_conform=None)


def with_axes(array: g2.Array, axes_head: np.ndarray) -> g2.Array:
    """The same sensing centres with other sensitive axes (cell frames as ``opm.make_info`` derives them)."""
    info = array.info.copy()
    with info._unlock():
        for ch, n in zip(info["chs"], axes_head):
            n = n / np.linalg.norm(n)
            ex, ey = opm.cell_frame(n)
            ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = ex, ey, n
    return g2.Array(array.name, info, array.kinds.copy(), array.coil_def, dict(array.meta))


def subset(array: g2.Array, keep: np.ndarray) -> g2.Array:
    idx = np.flatnonzero(keep)
    return g2.Array(array.name, mne.pick_info(array.info, idx), array.kinds[idx], array.coil_def, dict(array.meta, n_sites=int(len(idx))))


def site_feet(array: g2.Array, subject) -> np.ndarray:
    """Nearest MRI-scalp point of every sensing centre (MRI frame)."""
    pos = mne.transforms.apply_trans(subject.trans["trans"], np.array([c["loc"][:3] for c in array.info["chs"]]))
    tree = cKDTree(subject.scalp.rr)
    return subject.scalp.rr[tree.query(pos)[1]]


def headline(st, subject, arrays: dict, meas, cfg) -> dict:
    """Oracle comparisons with this head model; the background is calibrated on its Neuromag gradiometers (G2's rule)."""
    st.subject = subject
    squid = arrays["squid"]
    good = ~np.isin(squid.info.ch_names, cfg["sensors"]["bads"])
    grads = good & (squid.kinds == "grad")
    G = {n: st.gains(a, fullres=False) for n, a in arrays.items()}
    bs = G2.calibrate(st, squid, G["squid"][1], float(np.nanmedian(meas["brain"][grads])), grads)
    res = {n: G2.evaluate_all(G[n][0] * st.q, a, st.noise(a, G[n][1], bs, meas["environment"]), CONDS) for n, a in arrays.items()}
    out = dict(brain_scale=bs, sites={n: a.n for n, a in arrays.items() if n != "squid"})
    depth = st.src.depth_mm
    for n in arrays:
        if n == "squid":
            continue
        for cond in CONDS:
            o, s = G2.detect_of(res, n, None, cond), G2.detect_of(res, "squid", "combined", cond)
            c = G2.compare(o, s, st.rng, 1000)
            lr = np.log2(o / s)
            out[f"{n}/combined/{cond}"] = dict(ratio=float(2 ** c["median_log2"]), ci95=[float(2 ** x) for x in c["ci95"]],
                                               share_opm_better=c["share_opm_better"],
                                               by_depth=[dict(lo=lo, hi=hi, n=int(np.sum((depth >= lo) & (depth < hi))),
                                                              ratio=float(2 ** np.median(lr[(depth >= lo) & (depth < hi)])))
                                                         for lo, hi in zip(DEPTH_EDGES[:-1], DEPTH_EDGES[1:])])
    return out


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    st = G2.Study(cfg)
    G2.GROUPS = st.src.region
    sub4 = st.subject
    sub3 = stored_subject(sub4)
    squid_info = neuromag.load_info("T3")
    meas = g2.measured_noise(squid_info, st.filt, cfg["sensors"]["bads"])
    a3 = g2.build_arrays(sub3, st.dig)
    a4 = g2.build_arrays(sub4, st.dig)
    log(f"arrays: v3 dense {a3['opm_dense'].n}, matched {a3['opm_matched'].n}; v4 dense {a4['opm_dense'].n}, matched {a4['opm_matched'].n}")
    # v3-style axes at the v4 sites: the stored outer skin's normals averaged within 15 mm (A-OPM-AXIS before v4)
    stored = anatomy._outward(next(s for s in sub3.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD))
    rot = np.linalg.inv(sub4.trans["trans"])[:3, :3]
    old_axes = {n: opm.smoothed_normals(stored.rr, stored.nn, site_feet(a4[n], sub4), g2.AXIS_RADIUS) @ rot.T
                for n in ("opm_matched", "opm_dense")}
    low3 = site_feet(a3["opm_dense"], sub3)[:, 2].min()
    keep = site_feet(a4["opm_dense"], sub4)[:, 2] >= low3
    pick = ("squid", "opm_matched", "opm_dense")
    variants = {
        "v3": (sub3, {n: a3[n] for n in pick}),
        "v3_arrays_v4_bem": (sub4, {n: a3[n] for n in pick}),
        "v4_sites_v3_axes": (sub4, {"squid": a4["squid"], **{n: with_axes(a4[n], old_axes[n]) for n in ("opm_matched", "opm_dense")}}),
        "v4": (sub4, {n: a4[n] for n in pick}),
        "v4_dense_above_v3_low": (sub4, {"squid": a4["squid"], "opm_dense": subset(a4["opm_dense"], keep)}),
    }
    out = dict(status="NEW check: the G2 headline under the v3 and v4 head models (A-BEM-CONFORM), oracle, point gains",
               variants={}, v3_lowest_dense_site_mri_z_mm=float(low3 * 1e3), v4_dense_sites_below_it=int(np.sum(~keep)))
    for name, (subject, arrays) in variants.items():
        out["variants"][name] = headline(st, subject, arrays, meas, cfg)
        v = out["variants"][name]
        log(f"{name}: " + "; ".join(f"{k} {x['ratio']:.3f}" for k, x in v.items() if isinstance(x, dict) and "ratio" in x))
    out["runtime_s"] = time.time() - t0
    io.write_json(out, OUT / "head_surface_effect.json")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
