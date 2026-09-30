#!/usr/bin/env python3
"""Accuracy of the 3-layer BEM field near the head surface at the OPM cell's integration points
(pre-freeze review F2; PLAN v2 item 1; docs/methods.md section 8).

Convergence test within the same physical model: the 3-layer BEM of G2 (5,120-triangle head
surface) against the same model with the head surface subdivided once (20,480 triangles, the
same flat geometry, finer potential discretization and integration). Both are evaluated at point
magnetometers placed at every integration point of every OPM cell (27 per site, along the
channel's sensitive axis) and at every sensing centre, and for the cell channels themselves, for
the matched and dense arrays as built by opmsquid.g2. Sources: random usable cortical dipoles.

Errors are normalized by the source's field scale over the array (RMS of the fine-model centre
fields), so that zero crossings of a single field do not inflate them; channel errors are also
given relative to the channel's own field where that field is at least the array RMS. MNE refuses
point sensors inside the head surface, so only integration points outside it are probed.

Output: results/g2/near_mesh_check.json and Figure_G2_near_mesh.png.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402

from opmsquid import anatomy, forward, g2, io, neuromag, opm, paths  # noqa: E402

OUT = ROOT / "results" / "g2"
N_SOURCES = 300
SEED = 11
PROBE_MIN_HEIGHT_MM = 0.2  # probes must be outside the head surface
HEIGHT_EDGES_MM = (0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 20.0)


def cell_points(info):
    """Integration points (head frame), quadrature weights and axis of every OPM channel."""
    from mne.forward._make_forward import _create_meg_coils

    with mne.use_coil_def(opm.coil_def_file()):
        coils = _create_meg_coils(info["chs"], "accurate")
    return [(c["rmag"].copy(), c["w"].copy(), c["cosmag"][0] / np.linalg.norm(c["cosmag"][0])) for c in coils]


def probe_info(points, axes):
    """Point magnetometers (MNE coil 2000) at ``points`` (head frame) along ``axes``."""
    info = mne.create_info([f"P{i:05d}" for i in range(len(points))], 1000.0, "mag")
    with info._unlock():
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")
        for ch, p, n in zip(info["chs"], points, axes):
            a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ex = np.cross(a, n)
            ex /= np.linalg.norm(ex)
            ch["loc"][:3], ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = p, ex, np.cross(n, ex), n
            ch["coil_type"] = FIFF.FIFFV_COIL_POINT_MAGNETOMETER
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
    return info


def summarise(values, heights):
    out = []
    for lo, hi in zip(HEIGHT_EDGES_MM[:-1], HEIGHT_EDGES_MM[1:]):
        m = (heights >= lo) & (heights < hi)
        v = values[m].ravel()
        out.append(dict(lo_mm=lo, hi_mm=hi, n_points=int(m.sum()), median=float(np.median(v)) if v.size else None,
                        p95=float(np.percentile(v, 95)) if v.size else None, max=float(v.max()) if v.size else None))
    return out


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    subject = anatomy.load_sample()
    dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    cortex = anatomy.full_resolution(subject)
    rng = np.random.default_rng(SEED)
    src = rng.choice(np.flatnonzero(cortex.usable), N_SOURCES, replace=False)
    rr, nn = cortex.rr[src], cortex.nn[src]
    skin = next(s for s in subject.bem_surfaces if s["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    coarse = subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=0)  # 5,120-triangle head surface (v1)
    fine = subject.bem_model(g2.BEM_CONDUCTIVITY, head_refine=1)  # 20,480 (A-BEM-SKIN, v2 primary)
    arrays = {k: v for k, v in g2.build_arrays(subject, dig).items() if k in ("opm_matched", "opm_dense")}
    t = subject.trans["trans"]
    out = dict(status="NEW check (pre-freeze review F2): 3-layer BEM, head surface 5,120 vs 20,480 triangles",
               n_sources=N_SOURCES, seed=SEED, arrays={})
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for name, a in arrays.items():
        cells = cell_points(a.info)
        n_ch, n_pt = len(cells), len(cells[0][0])
        ipts = np.concatenate([c[0] for c in cells])
        iaxes = np.concatenate([np.repeat(c[2][None], n_pt, axis=0) for c in cells])
        height = opm.signed_distance(mne.transforms.apply_trans(t, ipts), skin).reshape(n_ch, n_pt) * 1e3
        probed = height.ravel() > PROBE_MIN_HEIGHT_MM
        centres = np.array([ch["loc"][:3] for ch in a.info["chs"]])
        info = probe_info(np.concatenate([ipts[probed], centres]), np.concatenate([iaxes[probed], [c[2] for c in cells]]))
        pc = forward.discrete_gain(info, subject.trans, rr, nn, coarse, None, use_cache=False)[0].astype(float)
        pf = forward.discrete_gain(info, subject.trans, rr, nn, fine, None, use_cache=False)[0].astype(float)
        cc = forward.discrete_gain(a.info, subject.trans, rr, nn, coarse, opm.coil_def_file(), use_cache=False)[0].astype(float)
        cf = forward.discrete_gain(a.info, subject.trans, rr, nn, fine, opm.coil_def_file(), use_cache=False)[0].astype(float)
        n_probe = int(probed.sum())
        scale = np.sqrt(np.mean(pf[n_probe:] ** 2, axis=0))  # per source: RMS of the fine centre fields over the array
        point_err = np.abs(pc[:n_probe] - pf[:n_probe]) / scale  # (n_probe, n_src)
        centre_err = np.abs(pc[n_probe:] - pf[n_probe:]) / scale
        chan_err = np.abs(cc - cf) / scale  # (n_ch, n_src)
        strong = np.abs(cf) >= scale
        chan_rel = np.where(strong, np.abs(cc / np.where(strong, cf, 1.0) - 1.0), np.nan)
        chan_rel_max = np.nanmax(np.where(np.isnan(chan_rel), -np.inf, chan_rel), axis=1)
        chan_rel_max[~np.isfinite(chan_rel_max)] = np.nan
        low = height.min(axis=1)
        order = np.argsort(-np.nan_to_num(chan_rel_max, nan=-1))
        out["arrays"][name] = dict(
            n_sites=n_ch, n_probed_points=n_probe, lowest_point_height_mm=dict(min=float(low.min()), median=float(np.median(low))),
            sites_lowest_point_below_mm={f"{d:g}": int(np.sum(low < d)) for d in (0.0, 1.0, 2.0, 3.0)},
            point_error_by_height=summarise(point_err, height.ravel()[probed]),
            centre_error=dict(median=float(np.median(centre_err)), p95=float(np.percentile(centre_err, 95)), max=float(centre_err.max())),
            channel_error_scaled=dict(median=float(np.median(chan_err)), p95=float(np.percentile(chan_err, 95)), max=float(chan_err.max())),
            channel_error_relative=dict(n_channels_above_1pct=int(np.nansum(chan_rel_max > 0.01)),
                                        n_channels_above_2pct=int(np.nansum(chan_rel_max > 0.02)),
                                        max=float(np.nanmax(chan_rel_max))),
            worst_channels=[dict(channel=int(c), lowest_point_height_mm=float(low[c]), max_relative_error=float(chan_rel_max[c]))
                            for c in order[:10]])
        rows = out["arrays"][name]["point_error_by_height"]
        mids = [(r["lo_mm"] + r["hi_mm"]) / 2 for r in rows]
        axs[0].semilogy(mids, [r["p95"] or np.nan for r in rows], "o-", label=f"{name}: 95th percentile")
        axs[0].semilogy(mids, [r["median"] or np.nan for r in rows], "s--", label=f"{name}: median")
        axs[1].semilogy(low, np.nanmax(chan_err, axis=1), "o", ms=3, alpha=0.6, label=name)
        print(f"{name}: {n_ch} sites, {n_probe} probed points; centre error p95 {100 * np.percentile(centre_err, 95):.3f} %; "
              f"channels relative error > 1 %: {np.nansum(chan_rel_max > 0.01)}, > 2 %: {np.nansum(chan_rel_max > 0.02)}, "
              f"max {100 * np.nanmax(chan_rel_max):.1f} %", flush=True)
        for r in rows:
            if r["n_points"]:
                print(f"   height {r['lo_mm']:5.1f}..{r['hi_mm']:5.1f} mm: n {r['n_points']:5d}  median {100 * r['median']:.3f} %"
                      f"  p95 {100 * r['p95']:.3f} %  max {100 * r['max']:.2f} %  (of the array field RMS)", flush=True)
        for w in out["arrays"][name]["worst_channels"][:5]:
            print(f"   worst channel {w['channel']}: lowest point {w['lowest_point_height_mm']:.2f} mm, relative error "
                  f"{100 * w['max_relative_error']:.1f} %", flush=True)
    axs[0].set_xlabel("integration-point height above the 5,120-triangle head surface [mm]")
    axs[0].set_ylabel("|coarse - refined| / array field RMS")
    axs[0].legend(fontsize=7)
    axs[0].set_title("point level: head surface 5,120 vs 20,480 triangles", fontsize=9)
    axs[1].set_xlabel("lowest integration point of the cell [mm above the head surface]")
    axs[1].set_ylabel("channel |coarse - refined| / array field RMS (max over sources)")
    axs[1].axvline(0, color="k", lw=0.6)
    axs[1].legend(fontsize=7)
    axs[1].set_title("channel level (27-point cell)", fontsize=9)
    fig.suptitle("G2 check: convergence of the 3-layer BEM field at OPM cell integration points", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_near_mesh.png", dpi=150)
    plt.close(fig)
    out["runtime_s"] = time.time() - t0
    io.write_json(out, OUT / "near_mesh_check.json")
    print(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
