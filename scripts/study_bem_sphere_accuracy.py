#!/usr/bin/env python3
"""Exact test of MNE's 3-layer BEM field close to the outer surface (supports the v2 fix of
pre-freeze review F2; docs/methods.md section 8).

Three concentric spheres (radii 80, 85, 90 mm; 0.3, 0.006, 0.3 S/m) meshed like the sample
subject's BEM (ico 4, 5,120 triangles) and with the outer sphere subdivided once (20,480). Outside
a spherically symmetric conductor the radial magnetic field does not depend on the volume
currents, so the analytic Sarvas field is exact for radial sensors. Point magnetometers point
radially at heights h above the outer sphere; random dipoles at 30-75 mm radius.

Error per sensor and source: |B_bem - B_exact|, normalized by the RMS of the exact field over all
sensors at that height (so zero crossings of single fields do not inflate it). Output:
results/g2/bem_sphere_check.json and Figure_G2_bem_sphere.png.
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

from opmsquid import anatomy, io  # noqa: E402

OUT = ROOT / "results" / "g2"
RADII = (0.080, 0.085, 0.090)  # inner skull, outer skull, scalp [m]
SIGMA = (0.3, 0.006, 0.3)  # brain, skull, scalp
HEIGHTS_MM = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0)
N_SENSORS, N_DIPOLES, SEED = 150, 120, 5


def sphere_surface(radius, grade, sid, sigma, subdivide=0):
    s = mne.surface._get_ico_surface(grade)
    rr, tris = s["rr"] / np.linalg.norm(s["rr"], axis=1, keepdims=True), s["tris"]
    if subdivide:
        fine = anatomy.Surface(rr, tris, rr.copy()).subdivided(subdivide)
        rr, tris = fine.rr / np.linalg.norm(fine.rr, axis=1, keepdims=True), fine.tris  # new vertices on the sphere
    surf = dict(id=sid, sigma=sigma, coord_frame=FIFF.FIFFV_COORD_MRI, rr=rr * radius, tris=tris.copy(), np=len(rr), ntri=len(tris))
    return mne.surface.complete_surface_info(surf, copy=False, verbose=False)


def radial_info(points):
    info = mne.create_info([f"S{i:04d}" for i in range(len(points))], 1000.0, "mag")
    with info._unlock():
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")
        for ch, p in zip(info["chs"], points):
            n = p / np.linalg.norm(p)
            a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ex = np.cross(a, n)
            ex /= np.linalg.norm(ex)
            ch["loc"][:3], ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = p, ex, np.cross(n, ex), n
            ch["coil_type"] = FIFF.FIFFV_COIL_POINT_MAGNETOMETER
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
    return info


def forward_field(info, bem_or_sphere, pos, ori):
    src = mne.setup_volume_source_space(pos=dict(rr=pos, nn=ori), verbose=False)
    ident = mne.transforms.Transform("head", "mri")
    fwd = mne.make_forward_solution(info, ident, src, bem_or_sphere, meg=True, eeg=False, mindist=0.0, verbose=False)
    g = fwd["sol"]["data"].reshape(len(info.ch_names), -1, 3)
    return np.einsum("cik,ik->ci", g, ori)


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    rng = np.random.default_rng(SEED)
    dirs = rng.normal(size=(N_SENSORS, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    d = rng.normal(size=(N_DIPOLES, 3))
    pos = d / np.linalg.norm(d, axis=1, keepdims=True) * rng.uniform(0.030, 0.075, (N_DIPOLES, 1))
    ori = rng.normal(size=(N_DIPOLES, 3))
    ori -= np.sum(ori * pos, axis=1, keepdims=True) * pos / np.sum(pos**2, axis=1, keepdims=True)  # tangential (radial is silent)
    ori /= np.linalg.norm(ori, axis=1, keepdims=True)
    sphere = mne.make_sphere_model(r0=(0.0, 0.0, 0.0), head_radius=None, verbose=False)
    ids = (FIFF.FIFFV_BEM_SURF_ID_BRAIN, FIFF.FIFFV_BEM_SURF_ID_SKULL, FIFF.FIFFV_BEM_SURF_ID_HEAD)
    models = {}
    for label, sub in (("scalp 5,120 triangles", 0), ("scalp 20,480 triangles", 1)):
        surfs = [sphere_surface(RADII[2], 4, ids[2], SIGMA[2], sub), sphere_surface(RADII[1], 4, ids[1], SIGMA[1]),
                 sphere_surface(RADII[0], 4, ids[0], SIGMA[0])]
        models[label] = mne.make_bem_solution(surfs, verbose=False)
    out = dict(status="NEW check: 3-layer sphere BEM vs the exact radial field", radii_mm=[r * 1e3 for r in RADII], sigma=SIGMA,
               n_sensors=N_SENSORS, n_dipoles=N_DIPOLES, heights_mm=list(HEIGHTS_MM), models={})
    for label, bem in models.items():
        rows = []
        for h in HEIGHTS_MM:
            info = radial_info(dirs * (RADII[2] + h * 1e-3))
            exact = forward_field(info, sphere, pos, ori)
            b = forward_field(info, bem, pos, ori)
            err = np.abs(b - exact) / np.sqrt(np.mean(exact**2, axis=0))
            rows.append(dict(height_mm=h, median=float(np.median(err)), p95=float(np.percentile(err, 95)), max=float(err.max())))
            print(f"{label}: h {h:4.1f} mm  median {100 * rows[-1]['median']:.3f} %  p95 {100 * rows[-1]['p95']:.3f} %  "
                  f"max {100 * rows[-1]['max']:.2f} %", flush=True)
        out["models"][label] = rows
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for label, rows in out["models"].items():
        ax.semilogy([r["height_mm"] for r in rows], [r["p95"] for r in rows], "o-", label=f"{label}: 95th percentile")
        ax.semilogy([r["height_mm"] for r in rows], [r["median"] for r in rows], "s--", label=f"{label}: median")
    ax.set_xlabel("sensor height above the outer sphere [mm]")
    ax.set_ylabel("|BEM - exact| / RMS of the exact field")
    ax.legend(fontsize=7)
    ax.set_title("3-layer sphere BEM (MNE linear collocation) vs the exact radial field", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_bem_sphere.png", dpi=150)
    plt.close(fig)
    out["runtime_s"] = time.time() - t0
    io.write_json(out, OUT / "bem_sphere_check.json")
    print(f"done in {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
