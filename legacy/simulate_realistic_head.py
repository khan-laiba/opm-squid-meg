#!/usr/bin/env python3
"""
Figure 3 of Jas et al. (2026, bioRxiv 2026.08.17.744953) with realistic anatomy: a real
cortical surface and a boundary-element (BEM) head model instead of the spherical head, for
single dipoles and extended cortical patches of different sizes.

Anatomy and forward model (MNE-Python)
--------------------------------------
* FreeSurfer ``fsaverage`` (the MNE template brain): cortical source space on the white-matter
  surface (ico-5, 20484 sources), MRI-derived scalp (``fsaverage-head.fif``) and inner skull.
  ``--subjects-dir`` / ``--subject`` select any FreeSurfer subject that has the same files
  (``bem/<subject>-ico-5-src.fif``, ``-5120-5120-5120-bem.fif``, ``-head.fif``,
  ``-fiducials.fif``; ``surf/?h.inflated`` and ``?h.sulc`` for the maps).
* Single-compartment BEM (inner skull, 5120 triangles, 0.3 S/m), the standard MEG head model
  (``mne.make_bem_solution``); fields from ``mne.make_forward_solution``.
* Sources point along the cortical normal (fixed orientation from cortical patch statistics,
  ``mne.convert_forward_solution(..., surf_ori=True, force_fixed=True, use_cps=True)``), i.e.
  realistic orientations: tangential in sulcal walls, radial (nearly MEG-silent) on gyral
  crowns.  A "dipole" is one source with Q = 30 nAm.  A patch of radius rho is the set of
  sources within the geodesic distance rho of a centre source, measured along the cortical mesh
  as MNE does (``mne.add_source_space_distances``: shortest paths along mesh edges, ~5 % longer
  than true geodesics, so a patch covers ~15 % less cortex than a disc of radius rho).  It
  carries a uniform current density (weights = cortical area of each source, ``pinfo``) with
  total moment 30 nAm; its field is the superposition of the source fields, so the opposite
  walls of a sulcus cancel.  Feeding the patch source estimates to
  ``mne.simulation.simulate_evoked`` reproduces the superposition (a consistency check of the
  MNE simulation path, not an independent test of the fields).

Sensors (the paper's idealisation, on a real head)
--------------------------------------------------
Point magnetometers measuring the field component normal to the scalp, covering the scalp
above the plane through the pre-auricular points and 3 cm above the nasion (helmet
coverage), on the MRI scalp mesh subdivided once (1.87 mm spacing, 28112 sensors): OPM on the
scalp (xi = 0) and SQUID 18 mm off the scalp along the scalp normal (a conformal off-scalp
array, as in the sphere model).

Signal, noise, SNR, depth
-------------------------
Signal = peak |B| over the array (as in the paper).  Noise as in Fig. 3: sigma_SQUID =
0.3546 pT, sigma_OPM = eta * sigma_SQUID (eta = 3).  Depth d = distance from the source (patch
centre) to the sensor-covered scalp (= the scalp for sources under the helmet; the distance to
the whole scalp is also written to the CSV).  With real anatomy the SNR ratio is not a function
of depth alone (depth explains 95 % of the variance of its logarithm for dipoles, 73 % for
20-mm patches), so the equal-SNR depth d_eq of the sphere model becomes a statistic: d_50, the
depth at which the cortex at that depth is split half and half between OPM and SQUID (area-
weighted fraction with SNR_OPM > SNR_SQUID in 2-mm depth bins; not defined -- reported as
censored -- if even the shallowest bin favours SQUID or no bin does), next to the fraction of
the whole cortex where SNR_OPM > SNR_SQUID.  The sphere model (``simulate_extended_sources.py``)
is shown for comparison with caps of the same area as the realistic patches; realistic d_50
and sphere d_eq agree to ~1 mm at eta = 3 but not in general (see the summary figure).

Accuracy
--------
* Sensor sampling: peaks on the 1.87-mm arrays are within 0.3 % (median 0.06 %) of those on
  0.93-mm arrays; the 3.7-mm arrays change d_50 by <= 0.2 mm.
* BEM discretisation: all sources are >= 5 mm from the inner skull; with the inner skull
  refined to 20480 triangles (same shape) the OPM/SQUID ratio changes by 0.04 % (median; 99th
  percentile 0.6 %), the OPM-favoured area by <= 0.02 percentage points and d_50 by <= 0.03 mm.

Outputs (``realistic_head_output/``)
------------------------------------
Figure_realistic_fig3.{png,pdf}     Fig. 3 for all cortical dipoles: signal, SNR, SNR ratio
Figure_realistic_fig3_patch_rho_XXmm.{png,pdf}   the same for the patches of radius XX mm
Figure_realistic_maps.{png,pdf}     SNR_OPM / SNR_SQUID on the inflated cortex, per patch size
Figure_realistic_summary.{png,pdf}  patch-size and eta dependence, comparison with the sphere
realistic_head_sources.csv          per-source depth, area, signals and SNR ratio per size
realistic_head_summary.json         model, checks and summary statistics
realistic_head_results.npz          simulated peak fields (``--replot`` re-analyses them)

Usage
-----
    python simulate_realistic_head.py                       # ~10 min; fsaverage from SUBJECTS_DIR
    python simulate_realistic_head.py --replot --eta 2.5    # re-analyse the saved fields
    python simulate_realistic_head.py --radii 0 10 20 --subjects-dir /path/to/subjects
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.spatial import cKDTree

import replicate_figure3 as R  # noqa: E402  (Agg backend, fonts, paper parameters)
import simulate_extended_sources as E  # noqa: E402  (sphere-model patches, for comparison)
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402

PATCH_RADII_MM = (0.0, 5.0, 10.0, 15.0, 20.0)
MAP_RADII_MM = (0.0, 10.0, 20.0)
ARRAYS = {"OPM": R.XI_OPM, "SQUID": R.XI_SQUID}  # scalp-to-sensor distance [m]
SENSOR_SUBDIVISIONS = 1  # scalp mesh (3.7 mm) subdivided once -> 1.87 mm sensor spacing
SENSOR_CHUNK = 2000  # sensors per forward computation (memory)
HELMET_BROW_OFFSET = 0.030  # helmet edge: plane through LPA, RPA and nasion + 3 cm
TRANS = mne.transforms.Transform("head", "mri")  # sensors are given in MRI coordinates
DEPTH_BIN_MM = 2.0
MIN_BIN_COUNT = 30  # sources per depth bin for binned statistics
ETA_GRID = np.round(np.linspace(1.0, 6.0, 101), 6)
SPHERE_DEPTHS_MM = np.arange(15.0, 85.5, 1.0)  # sphere-model reference curves
COLORS = {"OPM": "b", "SQUID": "r"}
RESULTS_FILE = "realistic_head_results.npz"


# =========================================================================================
# Anatomy
# =========================================================================================
def read_freesurfer_surface(fname: str) -> tuple[np.ndarray, np.ndarray]:
    """Vertices [mm] and triangles of a FreeSurfer binary triangle surface (?h.inflated)."""
    with open(fname, "rb") as fh:
        if fh.read(3) != b"\xff\xff\xfe":
            raise ValueError(f"not a FreeSurfer triangle surface: {fname}")
        fh.readline()  # "created by ..." line and the empty line after it
        fh.readline()
        n_vert, n_face = np.fromfile(fh, ">i4", 2)
        verts = np.fromfile(fh, ">f4", 3 * n_vert).reshape(-1, 3).astype(float)
        faces = np.fromfile(fh, ">i4", 3 * n_face).reshape(-1, 3)
    return verts, faces


def read_freesurfer_curv(fname: str) -> np.ndarray:
    """Per-vertex values of a FreeSurfer curvature file (?h.sulc; > 0 in sulci)."""
    with open(fname, "rb") as fh:
        if fh.read(3) != b"\xff\xff\xff":
            raise ValueError(f"not a FreeSurfer curvature file: {fname}")
        n_vert, _, _ = np.fromfile(fh, ">i4", 3)
        return np.fromfile(fh, ">f4", n_vert).astype(float)


@dataclass
class Mesh:
    rr: np.ndarray  # vertices [m]
    nn: np.ndarray  # outward unit normals
    tris: np.ndarray


@dataclass
class Anatomy:
    subjects_dir: str
    subject: str
    src: mne.SourceSpaces
    bem: mne.bem.ConductorModel
    scalp: Mesh
    fiducials: dict


def resolve_subjects_dir(subjects_dir: str | None, subject: str) -> str:
    if subjects_dir:
        return subjects_dir
    configured = mne.get_config("SUBJECTS_DIR") or os.environ.get("SUBJECTS_DIR")
    if configured and Path(configured, subject, "bem").is_dir():
        return configured
    if subject != "fsaverage":
        raise SystemExit("pass --subjects-dir for a subject other than fsaverage")
    return str(Path(mne.datasets.fetch_fsaverage(verbose=True)).parent)  # MNE download if missing


def load_source_space(subjects_dir: str, subject: str) -> mne.SourceSpaces:
    src = mne.read_source_spaces(Path(subjects_dir, subject, "bem", f"{subject}-ico-5-src.fif"), verbose=False)
    if src.kind != "surface" or len(src) != 2:
        raise ValueError("expected a two-hemisphere cortical surface source space")
    return src


def load_anatomy(subjects_dir: str, subject: str, max_radius: float) -> Anatomy:
    bem_dir = Path(subjects_dir, subject, "bem")
    src = load_source_space(subjects_dir, subject)
    # patch statistics: cortical area of each source and its normal (cortical patch statistics)
    mne.add_source_space_distances(src, dist_limit=0.0, verbose=False)
    if max_radius > 0:  # geodesic distances along the cortex up to the largest patch radius
        mne.add_source_space_distances(src, dist_limit=max_radius + 1e-3, verbose=False)
    surfs = mne.read_bem_surfaces(bem_dir / f"{subject}-5120-5120-5120-bem.fif", verbose=False)
    inner = next(s for s in surfs if s["id"] == FIFF.FIFFV_BEM_SURF_ID_BRAIN)
    inner["sigma"] = 0.3
    bem = mne.make_bem_solution([inner], verbose=False)  # single-compartment MEG model
    head = mne.read_bem_surfaces(bem_dir / f"{subject}-head.fif", verbose=False)[0]
    # MNE's vertex normals follow the (consistent) triangle orientation: make them point outwards
    outward = np.median(np.sum((head["rr"] - inner["rr"].mean(axis=0)) * head["nn"], axis=1)) > 0
    nn = head["nn"] if outward else -head["nn"]
    fids, _ = mne.io.read_fiducials(bem_dir / f"{subject}-fiducials.fif")
    names = {FIFF.FIFFV_POINT_LPA: "lpa", FIFF.FIFFV_POINT_NASION: "nasion", FIFF.FIFFV_POINT_RPA: "rpa"}
    fiducials = {names[f["ident"]]: np.asarray(f["r"], float) for f in fids if f["ident"] in names}
    return Anatomy(subjects_dir, subject, src, bem, Mesh(head["rr"], nn, head["tris"]), fiducials)


def subdivide(mesh: Mesh, times: int = 1) -> Mesh:
    """Split every triangle into four (edge midpoints; normals averaged), `times` times."""
    rr, nn, tris = mesh.rr, mesh.nn, mesh.tris
    for _ in range(times):
        edges = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
        edges, inv = np.unique(edges, axis=0, return_inverse=True)
        mid_nn = nn[edges[:, 0]] + nn[edges[:, 1]]
        mid_nn /= np.linalg.norm(mid_nn, axis=1, keepdims=True)
        m = inv.reshape(3, -1).T + len(rr)  # midpoints of edges (01, 12, 20) of each triangle
        tris = np.concatenate([np.c_[tris[:, 0], m[:, 0], m[:, 2]], np.c_[tris[:, 1], m[:, 1], m[:, 0]],
                               np.c_[tris[:, 2], m[:, 2], m[:, 1]], m])
        rr = np.vstack([rr, 0.5 * (rr[edges[:, 0]] + rr[edges[:, 1]])])
        nn = np.vstack([nn, mid_nn])
    return Mesh(rr, nn, tris)


def helmet_mesh(anat: Anatomy) -> Mesh:
    """Scalp above the plane through LPA, RPA and a point 3 cm above the nasion."""
    f = anat.fiducials
    brow = f["nasion"] + np.array([0.0, 0.0, HELMET_BROW_OFFSET])
    normal = np.cross(f["rpa"] - f["lpa"], brow - f["lpa"])
    normal *= np.sign(normal[2]) / np.linalg.norm(normal)
    above = (anat.scalp.rr - f["lpa"]) @ normal > 0
    return Mesh(anat.scalp.rr, anat.scalp.nn, anat.scalp.tris[above[anat.scalp.tris].all(axis=1)])


def mesh_points(mesh: Mesh) -> tuple[np.ndarray, np.ndarray]:
    used = np.unique(mesh.tris)
    return mesh.rr[used], mesh.nn[used]


def make_sensor_info(pos: np.ndarray, normals: np.ndarray) -> mne.Info:
    """Point magnetometers at `pos` measuring the field component along `normals`."""
    info = mne.create_info([f"MAG{i:05d}" for i in range(len(pos))], 1000.0, "mag")
    with info._unlock():  # sensor geometry set explicitly (device = head = MRI coordinates)
        info["dev_head_t"] = mne.transforms.Transform("meg", "head")
        for ch, p, n in zip(info["chs"], pos, normals):
            a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ex = np.cross(a, n)
            ex /= np.linalg.norm(ex)
            ch["loc"][:3], ch["loc"][3:6], ch["loc"][6:9], ch["loc"][9:12] = p, ex, np.cross(n, ex), n
            ch["coil_type"] = FIFF.FIFFV_COIL_POINT_MAGNETOMETER
            ch["coord_frame"] = FIFF.FIFFV_COORD_DEVICE
    return info


def fixed_forward(anat: Anatomy, pos: np.ndarray, normals: np.ndarray):
    """BEM forward solution with sources fixed along the cortical normal."""
    info = make_sensor_info(pos, normals)
    fwd = mne.make_forward_solution(info, TRANS, anat.src, anat.bem, meg=True, eeg=False, mindist=0.0,
                                    verbose=False)
    fwd = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True, use_cps=True, verbose=False)
    if any(not np.array_equal(a["vertno"], b["vertno"]) for a, b in zip(fwd["src"], anat.src)):
        raise RuntimeError("the forward model dropped cortical sources (outside the inner skull?)")
    return fwd, info


# =========================================================================================
# Sources and simulation
# =========================================================================================
def vertex_areas(src: mne.SourceSpaces) -> list[np.ndarray]:
    """Cortical area [m^2] represented by each source (its patch of the full-resolution surface)."""
    out = []
    for s in src:
        rr, tris = s["rr"], s["tris"]
        tri_area = 0.5 * np.linalg.norm(np.cross(rr[tris[:, 1]] - rr[tris[:, 0]], rr[tris[:, 2]] - rr[tris[:, 0]]),
                                        axis=1)
        v_area = np.bincount(tris.ravel(), np.repeat(tri_area / 3.0, 3), minlength=len(rr))
        out.append(np.array([v_area[s["pinfo"][k]].sum() for k in s["patch_inds"]]))
    return out


def patch_operator(src: mne.SourceSpaces, radius: float, areas: list[np.ndarray]) -> sp.csr_matrix:
    """(n_sources, n_sources) weights: row i = uniform current density over the sources within
    the geodesic distance `radius` of source i, i.e. weights proportional to cortical area,
    summing to 1 (the patch carries the total moment)."""
    blocks = []
    for s, area in zip(src, areas):
        n = s["nuse"]
        near = sp.identity(n, format="csr")
        if radius > 0:
            d = sp.csr_matrix(s["dist"])[s["vertno"]][:, s["vertno"]]  # stored: 0 < d <= limit
            d.data = (d.data <= radius).astype(float)
            d.eliminate_zeros()
            near = near + d
        w = sp.csr_matrix(near) @ sp.diags(area)
        blocks.append(sp.diags(1.0 / np.asarray(w.sum(axis=1)).ravel()) @ w)
    return sp.block_diag(blocks, format="csr")


def peak_fields(anat: Anatomy, pos, normals, operators: dict, label: str) -> dict:
    """Peak |B| [T] over the array for every source/patch of every operator (running maximum
    over sensor chunks of the BEM forward model)."""
    peaks = {k: np.zeros(op.shape[0]) for k, op in operators.items()}
    for start in range(0, len(pos), SENSOR_CHUNK):
        sl = slice(start, min(start + SENSOR_CHUNK, len(pos)))
        t0 = time.time()
        fwd, _ = fixed_forward(anat, pos[sl], normals[sl])
        gain_t = fwd["sol"]["data"].T  # (n_sources, n_sensors) [T / (A m)]
        for k, op in operators.items():
            peaks[k] = np.maximum(peaks[k], np.abs(op @ gain_t).max(axis=1) * R.Q_DIPOLE)
        print(f"  {label}: sensors {sl.start}-{sl.stop} of {len(pos)} ({time.time() - t0:.1f} s)", flush=True)
    return peaks


def check_with_simulate_evoked(anat: Anatomy, pos, normals, operator: sp.csr_matrix, n_check: int = 5,
                               n_sensors: int = 500, seed: int = 0) -> float:
    """Simulate a few patches as source estimates with mne.simulation.simulate_evoked; largest
    deviation from the superposition used in peak_fields (relative to the peak)."""
    rng = np.random.default_rng(seed)
    sel = rng.choice(len(pos), n_sensors, replace=False)
    fwd, info = fixed_forward(anat, pos[sel], normals[sel])
    vertices = [s["vertno"] for s in fwd["src"]]
    err = 0.0
    for i in rng.choice(operator.shape[0], n_check, replace=False):
        amp = operator[i].toarray().ravel() * R.Q_DIPOLE
        stc = mne.SourceEstimate(amp[:, None], vertices, tmin=0.0, tstep=1e-3)
        sim = mne.simulation.simulate_evoked(fwd, stc, info, cov=None, nave=np.inf, verbose=False).data[:, 0]
        ref = operator[i] @ fwd["sol"]["data"].T * R.Q_DIPOLE
        err = max(err, float(np.max(np.abs(sim - ref.ravel())) / np.max(np.abs(ref))))
    return err


def simulate(anat: Anatomy, radii_mm, subdivisions: int = SENSOR_SUBDIVISIONS) -> dict:
    """Peak fields of all dipoles and patches on both arrays, source depths and areas (the
    expensive part; saved to RESULTS_FILE)."""
    t_start = time.time()
    src = anat.src
    areas = vertex_areas(src)
    area = np.concatenate(areas)
    operators = {r: patch_operator(src, r * 1e-3, areas) for r in radii_mm}
    helmet = helmet_mesh(anat)
    sensor_mesh = subdivide(helmet, subdivisions)
    pos, normals = mesh_points(sensor_mesh)
    spacing = np.linalg.norm(sensor_mesh.rr[sensor_mesh.tris[:, 0]] - sensor_mesh.rr[sensor_mesh.tris[:, 1]], axis=1)
    print(f"{anat.subject}: {len(area)} cortical sources ({area.sum() * 1e4:.0f} cm2 of cortex), "
          f"{len(pos)} sensors per array", flush=True)
    peaks = {name: peak_fields(anat, pos + xi * normals, normals, operators, name) for name, xi in ARRAYS.items()}
    check = check_with_simulate_evoked(anat, pos, normals, operators[max(radii_mm)])
    src_rr = np.concatenate([s["rr"][s["vertno"]] for s in src])
    raw = dict(subject=anat.subject, subjects_dir=anat.subjects_dir, radii_mm=np.array(radii_mm, float),
               depth_mm=cKDTree(mesh_points(subdivide(helmet, 3))[0]).query(src_rr)[0] * 1e3,
               depth_scalp_mm=cKDTree(subdivide(anat.scalp, 2).rr).query(src_rr)[0] * 1e3, area_m2=area,
               hemi=np.concatenate([np.full(s["nuse"], h) for h, s in enumerate(src)]),
               vertno=np.concatenate([s["vertno"] for s in src]), n_sensors=len(pos),
               sensor_spacing_mm=float(np.median(spacing) * 1e3), simulate_evoked_rel_err=check)
    for r in radii_mm:
        members = (operators[r] > 0).astype(float)
        raw[f"b_opm_{r:g}"], raw[f"b_squid_{r:g}"] = peaks["OPM"][r], peaks["SQUID"][r]
        raw[f"patch_area_cm2_{r:g}"] = np.asarray(members @ area).ravel() * 1e4
        raw[f"n_vertices_{r:g}"] = np.asarray(members.sum(axis=1)).ravel()
    raw["runtime_s"] = time.time() - t_start
    return raw


def save_raw(raw: dict, path: Path) -> None:
    np.savez_compressed(path, **{k: np.asarray(v) for k, v in raw.items()})


def load_raw(path: Path) -> dict:
    with np.load(path) as f:
        raw = {k: f[k] for k in f.files}
    for key, cast in (("subject", str), ("subjects_dir", str), ("n_sensors", int), ("sensor_spacing_mm", float),
                      ("simulate_evoked_rel_err", float), ("runtime_s", float)):
        raw[key] = cast(raw[key])
    return raw


# =========================================================================================
# Analysis
# =========================================================================================
def weighted_quantiles(values, weights, qs):
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cdf = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(qs, cdf, v)


def _bins(depth):
    return np.floor(depth / DEPTH_BIN_MM).astype(int)


def depth_profile(depth, values, weights, favoured=None):
    """Area-weighted quartiles of `values` (and the area fraction where `favoured`) in 2-mm
    depth bins that hold at least MIN_BIN_COUNT sources with finite values."""
    idx, ok = _bins(depth), np.isfinite(values)
    rows = []
    for b in range(idx.max() + 1):
        m = (idx == b) & ok
        if m.sum() >= MIN_BIN_COUNT:
            q = weighted_quantiles(values[m], weights[m], [0.25, 0.5, 0.75])
            frac = np.sum(weights[m] * favoured[m]) / np.sum(weights[m]) if favoured is not None else np.nan
            rows.append([(b + 0.5) * DEPTH_BIN_MM, *q, frac])
    c, q25, q50, q75, frac = np.array(rows).T
    return dict(center=c, q25=q25, median=q50, q75=q75, fraction=frac)


def fraction_profile(depth, favoured, weights):
    """Bin centres and area fraction where `favoured` (bins as in depth_profile)."""
    idx = _bins(depth)
    count, w, wf = np.bincount(idx), np.bincount(idx, weights), np.bincount(idx, weights * favoured)
    keep = count >= MIN_BIN_COUNT
    return ((np.arange(count.size) + 0.5) * DEPTH_BIN_MM)[keep], wf[keep] / w[keep]


def half_depth(center, fraction) -> float:
    """d_50: depth at which the area fraction favouring OPM falls below 1/2 (linear
    interpolation between bin centres); NaN (censored) if it is below 1/2 already in the
    shallowest bin or never falls below 1/2."""
    below = np.flatnonzero(fraction < 0.5)
    if below.size == 0 or below[0] == 0:
        return np.nan
    j = below[0]
    f0, f1 = fraction[j - 1] - 0.5, fraction[j] - 0.5
    return float(center[j - 1] + (center[j] - center[j - 1]) * f0 / (f0 - f1))


def d50_text(d50: float, center, fraction) -> str:
    if np.isfinite(d50):
        return f"{d50:.1f} mm"
    return f"< {center[0]:g} mm" if fraction[0] < 0.5 else f"> {center[-1]:g} mm"


def depth_r2(depth, y, weights) -> float:
    """Area-weighted fraction of the variance of y explained by its 2-mm depth-bin means."""
    idx = _bins(depth)
    sw, swy = np.bincount(idx, weights), np.bincount(idx, weights * y)
    fitted = (swy / np.where(sw > 0, sw, 1.0))[idx]
    mean = np.sum(weights * y) / np.sum(weights)
    return float(1.0 - np.sum(weights * (y - fitted) ** 2) / np.sum(weights * (y - mean) ** 2))


def equal_area_radius_mm(area_cm2: float, r_q: float = R.BRAIN_RADIUS) -> float:
    """Geodesic radius [mm] of the spherical cap of the given area on the brain surface (r = b)."""
    return float(r_q * np.arccos(1.0 - area_cm2 * 1e-4 / (2.0 * np.pi * r_q**2)) * 1e3)


def sphere_peak(d_mm: float, radius_mm: float, xi: float) -> float:
    """Peak radial field [T] of the spherical model (h = 95 mm) for a tangential 30 nAm dipole
    (paper Eq. 1) or spherical-cap patch (simulate_extended_sources.py) at depth d_mm; NaN
    where the cap would exceed a hemisphere."""
    r_q, s = R.HEAD_RADIUS - d_mm * 1e-3, R.HEAD_RADIUS + xi
    if radius_mm == 0:
        return float(R.bmax_analytic(r_q, s)[0])
    if not E.patch_defined(r_q, radius_mm * 1e-3):
        return np.nan
    return E.peak_field_analytic(r_q, radius_mm * 1e-3, s)


def sphere_curve(radius_mm: float, eta: float) -> dict:
    """Sphere-model signals, signal and SNR ratio and equal-SNR depth d_eq (paper Eq. 3) of a
    dipole or cap of the given radius, on SPHERE_DEPTHS_MM."""
    d = SPHERE_DEPTHS_MM
    b = {k: np.array([sphere_peak(x, radius_mm, xi) for x in d]) for k, xi in ARRAYS.items()}
    ratio = b["OPM"] / b["SQUID"]
    g = ratio - eta
    k = np.flatnonzero(np.isfinite(g[:-1]) & np.isfinite(g[1:]) & (np.sign(g[:-1]) != np.sign(g[1:])))
    d_eq = np.nan
    if k.size:
        d_eq = brentq(lambda x: sphere_peak(x, radius_mm, ARRAYS["OPM"]) / sphere_peak(x, radius_mm, ARRAYS["SQUID"])
                      - eta, d[k[0]], d[k[0] + 1], xtol=1e-6)
    ok = np.isfinite(ratio)  # d_eq(eta) by inverting the (decreasing) ratio; NaN outside its range
    d_eq_vs_eta = np.interp(ETA_GRID, ratio[ok][::-1], d[ok][::-1], left=np.nan, right=np.nan)
    return dict(radius_mm=radius_mm, depth=d, b_opm=b["OPM"], b_squid=b["SQUID"], signal_ratio=ratio,
                snr_ratio=ratio / eta, d_eq=float(d_eq), d_eq_vs_eta=d_eq_vs_eta)


def analyse(raw: dict, eta: float) -> dict:
    """SNRs, depth profiles, d_50 (also vs eta) and sphere-model comparisons from the peaks."""
    depth, area = raw["depth_mm"], raw["area_m2"]
    sigma_squid = float(R.bmax_analytic(R.NOISE_REF_RADIUS, R.HEAD_RADIUS + R.XI_SQUID)[0])  # as Fig. 3
    sigma_opm = eta * sigma_squid
    sizes, curves = {}, {}
    for r in (float(x) for x in raw["radii_mm"]):
        b_opm, b_squid = raw[f"b_opm_{r:g}"], raw[f"b_squid_{r:g}"]
        signal_ratio = b_opm / b_squid
        ratio = signal_ratio / eta
        patch_area = raw[f"patch_area_cm2_{r:g}"]
        rho_eq = 0.0 if r == 0 else round(equal_area_radius_mm(float(np.median(patch_area))), 2)
        for rr in (r, rho_eq):
            if rr not in curves:
                curves[rr] = sphere_curve(rr, eta)
        sph = curves[rho_eq]
        prof = depth_profile(depth, ratio, area, favoured=ratio > 1)
        d50 = half_depth(prof["center"], prof["fraction"])
        ok = np.isfinite(sph["b_opm"])
        loss = b_opm / np.interp(depth, sph["depth"][ok], sph["b_opm"][ok], left=np.nan, right=np.nan)
        sizes[r] = dict(
            b_opm=b_opm, b_squid=b_squid, snr_opm=b_opm / sigma_opm, snr_squid=b_squid / sigma_squid,
            signal_ratio=signal_ratio, snr_ratio=ratio, patch_area_cm2=patch_area,
            n_vertices=raw[f"n_vertices_{r:g}"], profile=prof, d50=d50,
            d50_text=d50_text(d50, prof["center"], prof["fraction"]),
            d50_vs_eta=np.array([half_depth(*fraction_profile(depth, signal_ratio > e, area)) for e in ETA_GRID]),
            opm_area_fraction=float(np.sum(area * (ratio > 1)) / area.sum()),
            opm_area_fraction_vs_eta=np.array([np.sum(area * (signal_ratio > e)) / area.sum() for e in ETA_GRID]),
            depth_r2=depth_r2(depth, np.log(ratio), area), loss_profile=depth_profile(depth, loss, area),
            sphere=curves[r], sphere_eq=sph)
    return dict(raw=raw, eta=eta, subject=raw["subject"], depth_mm=depth, area_m2=area,
                n_sensors=raw["n_sensors"], sigma_squid_pT=sigma_squid * 1e12, sigma_opm_pT=sigma_opm * 1e12,
                sizes=sizes)


# =========================================================================================
# Figures
# =========================================================================================
def _rc(fonts):
    family, serif = fonts.regular.get_name(), fonts.italic_serif.get_name()
    return dict(R.FIGURE_RC, **{"mathtext.fontset": "custom", "mathtext.rm": family,
                                "mathtext.bf": f"{family}:bold", "mathtext.it": f"{serif}:italic",
                                "mathtext.default": "it", "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                                "ytick.major.width": 0.8, "xtick.labelsize": 8, "ytick.labelsize": 8})


def _small(fonts, size=7.0):
    p = fonts.regular.copy()
    p.set_size(size)
    return p


def _panel(ax, letter, title, fonts):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontproperties=fonts.bold, va="bottom")
    ax.text(0.0, 1.05, title, transform=ax.transAxes, fontproperties=fonts.regular, va="bottom")


def _size_label(r):
    return "dipole" if r == 0 else f"ρ = {r:g} mm"


def _xmax(res):
    return 10 * np.ceil(np.percentile(res["depth_mm"], 99.9) / 10)


def _save(fig, outdir: Path, stem: str) -> list[str]:
    paths = []
    for ext in ("png", "pdf"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300, facecolor="white")
        paths.append(str(p))
    plt.close(fig)
    return paths


def figure_fig3(res: dict, radius: float, fonts, outdir: Path, stem: str) -> list[str]:
    """Paper Fig. 3 (A: signal, B: SNR) plus C: SNR ratio, for the dipoles (radius 0) or the
    patches of one size centred on every cortical source."""
    s = res["sizes"][radius]
    depth, area = res["depth_mm"], res["area_m2"]
    sph = s["sphere_eq"]
    sigma = {"OPM": res["sigma_opm_pT"], "SQUID": res["sigma_squid_pT"]}
    small = _small(fonts)
    with plt.rc_context(_rc(fonts)):
        fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6))
        fig.subplots_adjust(left=0.055, right=0.99, top=0.84, bottom=0.15, wspace=0.27)
        for ax, key in ((axs[0], "b"), (axs[1], "snr")):
            scale = 1e12 if key == "b" else 1.0
            for arr in ("OPM", "SQUID"):
                vals = s[f"{key}_{arr.lower()}"] * scale
                ax.scatter(depth, vals, s=1.0, color=COLORS[arr], alpha=0.12, lw=0, rasterized=True)
                prof = depth_profile(depth, vals, area)
                ax.plot(prof["center"], prof["median"], color=COLORS[arr], lw=1.5,
                        ls="-" if arr == "OPM" else (0, (3.7, 1.6)))
                ref = sph[f"b_{arr.lower()}"] * 1e12
                ax.plot(sph["depth"], ref if key == "b" else ref / sigma[arr], color=COLORS[arr], lw=0.9, ls=":")
            top = max(np.percentile(s[f"{key}_opm"], 99.7), np.percentile(s[f"{key}_squid"], 99.7)) * scale
            ax.set_ylim(0, top * 1.08)
        ax = axs[2]
        ax.scatter(depth, s["snr_ratio"], s=1.0, color="0.35", alpha=0.12, lw=0, rasterized=True)
        prof = s["profile"]
        ax.fill_between(prof["center"], prof["q25"], prof["q75"], color="0.55", alpha=0.35, lw=0)
        ax.plot(prof["center"], prof["median"], color="k", lw=1.5)
        ax.plot(sph["depth"], sph["snr_ratio"], color="k", lw=0.9, ls=":")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_yscale("log")
        lo, hi = np.percentile(s["snr_ratio"], [0.2, 99.8])
        ax.set_ylim(min(lo, 0.5) / 1.1, max(hi, 2.0) * 1.1)
        ticks = [t for t in (0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10) if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()
        d50 = s["d50"]
        for a in axs:
            a.set_xlim(0, _xmax(res))
            a.set_xlabel("Source depth below scalp ($d$) [mm]", fontproperties=fonts.regular)
            if np.isfinite(d50):
                a.axvline(d50, color="k", ls=":", lw=1.2)
        _panel(axs[0], "A", "Signal [pT]", fonts)
        _panel(axs[1], "B", "Signal-to-Noise Ratio", fonts)
        _panel(axs[2], "C", r"SNR$_\mathrm{OPM}$ / SNR$_\mathrm{SQUID}$", fonts)
        handles = [Line2D([], [], color="b", lw=1.5), Line2D([], [], color="r", lw=1.5, ls=(0, (3.7, 1.6))),
                   Line2D([], [], color="0.3", lw=0.9, ls=":")]
        sphere_label = ("sphere model (tangential dipole)" if radius == 0 else
                        "sphere model (tangential cap, same area:\n"
                        f"ρ = {sph['radius_mm']:.1f} mm on the brain surface)")
        axs[0].legend(handles, ["OPM (on-scalp)", "SQUID (off-scalp)", sphere_label],
                      prop=small, frameon=False, loc="upper right")
        if np.isfinite(d50):
            axs[2].text(d50 + 1.0, axs[2].get_ylim()[0] * 1.06, f"$d_{{50}}$ = {d50:.1f} mm", fontproperties=small)
        d50_line = "" if np.isfinite(d50) else f"\n$d_{{50}}$ {s['d50_text']} (censored)"
        d_eq = f"= {sph['d_eq']:.1f} mm" if np.isfinite(sph["d_eq"]) else "none (no crossing at 15-85 mm)"
        axs[2].text(0.98, 0.97, f"OPM better for {100 * s['opm_area_fraction']:.0f} % of the cortex\n"
                    f"sphere model: $d_\\mathrm{{eq}}$ {d_eq}{d50_line}", transform=axs[2].transAxes,
                    ha="right", va="top", fontproperties=small, linespacing=1.4)
        what = (f"30 nAm dipoles normal to the cortex at each of {len(depth)} sources" if radius == 0 else
                f"30 nAm patches of geodesic radius {radius:g} mm (median area {np.median(s['patch_area_cm2']):.1f} "
                f"cm²) centred on each of {len(depth)} sources")
        fig.text(0.5, 0.985, f"Realistic head: {res['subject']} cortex, {what}; single-layer BEM; {res['n_sensors']} "
                 f"sensors per array; η = {res['eta']:g}", ha="center", va="top", fontproperties=small,
                 color="0.3")
        fig.text(0.5, 0.945, "Lines: area-weighted median per 2-mm depth bin (C: with interquartile range); "
                 "dotted curves: sphere model; dotted line: $d_{50}$, where OPM is better for half of the cortex "
                 "at that depth", ha="center", va="top", fontproperties=small, color="0.3")
        return _save(fig, outdir, stem)


def _render_view(ax, verts, tris, values, gyral, camera_x: float, cmap, norm):
    """Orthographic view of a surface from x = camera_x * infinity (z up), back faces culled
    and triangles painted back to front; values drawn per triangle, equal-SNR line on top."""
    d = np.array([camera_x, 0.0, 0.0])
    screen = np.column_stack([verts[:, 1] * camera_x, verts[:, 2]])  # camera at -x: anterior left
    tri_n = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]], verts[tris[:, 2]] - verts[tris[:, 0]])
    cos_view = tri_n @ d / np.linalg.norm(tri_n, axis=1)
    t = tris[cos_view > 0]
    t = t[np.argsort(verts[t].mean(axis=1) @ d)]
    rgba = cmap(norm(values[t].mean(axis=1)))
    rgba[:, :3] *= (0.80 + 0.20 * gyral[t].mean(axis=1))[:, None]  # sulci shaded darker
    ax.add_collection(PolyCollection(screen[t], facecolors=rgba, edgecolors="none", antialiased=False,
                                     rasterized=True))
    # equal-SNR line, away from the silhouette (edge-on triangles give broken contours)
    ax.tricontour(mtri.Triangulation(screen[:, 0], screen[:, 1], tris[cos_view > 0.35]), values, levels=[0.0],
                  colors="k", linewidths=0.5)
    ax.set_xlim(screen[:, 0].min() - 2, screen[:, 0].max() + 2)
    ax.set_ylim(screen[:, 1].min() - 2, screen[:, 1].max() + 2)
    ax.set_aspect("equal")
    ax.axis("off")


def figure_maps(res: dict, src: mne.SourceSpaces, subjects_dir: str, fonts, outdir: Path,
                radii=MAP_RADII_MM) -> list[str]:
    radii = [r for r in radii if r in res["sizes"]] or sorted(res["sizes"])[:3]
    surf_dir = Path(subjects_dir, res["subject"], "surf")
    hemis = []
    for h, name in enumerate(("lh", "rh")):
        verts, _ = read_freesurfer_surface(str(surf_dir / f"{name}.inflated"))
        sulc = read_freesurfer_curv(str(surf_dir / f"{name}.sulc"))
        s = src[h]
        remap = np.full(s["np"], -1)
        remap[s["vertno"]] = np.arange(s["nuse"])
        tris = remap[s["use_tris"]]
        hemis.append((verts[s["vertno"]], tris[(tris >= 0).all(axis=1)], (sulc[s["vertno"]] < 0).astype(float)))
    logs = np.concatenate([np.log10(res["sizes"][r]["snr_ratio"]) for r in radii])
    lim = max(np.log10(2.0), np.ceil(np.percentile(np.abs(logs), 99) / np.log10(2.0) * 2) / 2 * np.log10(2.0))
    norm = plt.Normalize(-lim, lim)
    cmap = plt.get_cmap("RdBu")  # blue = OPM better (the OPM colour of Fig. 3), red = SQUID better
    small = _small(fonts)
    n_lh = src[0]["nuse"]
    aspect = max(np.ptp(v[:, 1]) / np.ptp(v[:, 2]) for v, _, _ in hemis)  # width / height of a view
    left, right, wspace, hspace, top_in, bottom_in = 0.13, 0.9, 0.02, 0.08, 0.5, 0.06
    row = 9.0 * (right - left) / (4 + 3 * wspace) / aspect * 1.04
    height = row * (len(radii) + (len(radii) - 1) * hspace) + top_in + bottom_in
    with plt.rc_context(_rc(fonts)):
        fig, axs = plt.subplots(len(radii), 4, figsize=(9.0, height), squeeze=False)
        fig.subplots_adjust(left=left, right=right, top=1 - top_in / height, bottom=bottom_in / height,
                            wspace=wspace, hspace=hspace)
        for i, r in enumerate(radii):
            logr = np.log10(res["sizes"][r]["snr_ratio"])
            for j, (h, cam) in enumerate(((0, -1.0), (0, 1.0), (1, -1.0), (1, 1.0))):
                verts, tris, gyral = hemis[h]
                vals = logr[:n_lh] if h == 0 else logr[n_lh:]
                _render_view(axs[i, j], verts, tris, vals, gyral, cam, cmap, norm)
            st = res["sizes"][r]
            area_txt = "" if r == 0 else f" ({np.median(st['patch_area_cm2']):.1f} cm²)"
            axs[i, 0].text(-0.06, 0.5, f"{_size_label(r)}{area_txt}\nOPM better for\n"
                           f"{100 * st['opm_area_fraction']:.0f} % of the cortex", transform=axs[i, 0].transAxes,
                           ha="right", va="center", fontproperties=small, linespacing=1.35)
        for j, title in enumerate(("left hemisphere, lateral", "left hemisphere, medial",
                                   "right hemisphere, medial", "right hemisphere, lateral")):
            axs[0, j].set_title(title, fontproperties=small, pad=2)
        y0, y1 = bottom_in / height, 1 - top_in / height
        cax = fig.add_axes([0.915, y0 + 0.12 * (y1 - y0), 0.012, 0.76 * (y1 - y0)])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        ticks = [t for t in (0.1, 0.2, 0.25, 1 / 3, 0.5, 1, 2, 3, 4, 5, 10) if abs(np.log10(t)) <= lim + 1e-9]
        cb.set_ticks(np.log10(ticks))
        cb.set_ticklabels(["1/3" if abs(t - 1 / 3) < 1e-9 else f"{t:g}" for t in ticks])
        cb.ax.tick_params(labelsize=7, width=0.6)
        cb.outline.set_linewidth(0.6)
        cb.set_label(r"SNR$_\mathrm{OPM}$ / SNR$_\mathrm{SQUID}$" + f"  (η = {res['eta']:g})",
                     fontproperties=small)
        fig.text(0.515, 0.995, f"Inflated {res['subject']} cortex (sulci shaded darker); each source coloured by the "
                 "SNR ratio of the dipole or patch centred on it.  Blue: OPM better, red: SQUID better, black line: "
                 "equal SNR.", ha="center", va="top", fontproperties=small, color="0.3")
        return _save(fig, outdir, "Figure_realistic_maps")


def figure_summary(res: dict, fonts, outdir: Path) -> list[str]:
    radii = sorted(res["sizes"])
    greys = plt.cm.Greys(np.linspace(0.95, 0.4, len(radii)))
    small = _small(fonts)
    eta = res["eta"]
    with plt.rc_context(_rc(fonts)):
        fig, axs = plt.subplots(2, 3, figsize=(12.6, 7.0))
        fig.subplots_adjust(left=0.05, right=0.985, top=0.9, bottom=0.08, wspace=0.27, hspace=0.45)
        (ax_a, ax_b, ax_c), (ax_d, ax_e, ax_f) = axs
        handles, labels = [], []
        for k, r in enumerate(radii):
            st, sph = res["sizes"][r], res["sizes"][r]["sphere_eq"]
            prof = st["profile"]
            (line,) = ax_a.plot(prof["center"], prof["median"], color=greys[k], lw=1.5)
            ax_a.plot(sph["depth"], sph["snr_ratio"], color=greys[k], lw=0.9, ls=":")
            ax_b.plot(prof["center"], 100 * prof["fraction"], color=greys[k], lw=1.5)
            ax_c.plot(ETA_GRID, 100 * st["opm_area_fraction_vs_eta"], color=greys[k], lw=1.5)
            ax_d.plot(ETA_GRID, st["d50_vs_eta"], color=greys[k], lw=1.5)
            ax_d.plot(ETA_GRID, sph["d_eq_vs_eta"], color=greys[k], lw=0.9, ls=":")
            ax_f.plot(st["loss_profile"]["center"], st["loss_profile"]["median"], color=greys[k], lw=1.5)
            handles.append(line)
            area_txt = "" if r == 0 else f" ({np.median(st['patch_area_cm2']):.1f} cm²)"
            labels.append(_size_label(r) + area_txt)
        handles.append(Line2D([], [], color="0.3", lw=0.9, ls=":"))
        labels.append("sphere model (caps of equal area)")
        ax_a.axhline(1.0, color="k", lw=0.6)
        ax_a.set_yscale("log")
        ax_a.set_ylim(0.25, 4)
        ax_a.set_yticks([0.25, 0.5, 1, 2, 4])
        ax_a.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
        ax_a.minorticks_off()
        ax_a.legend(handles, labels, prop=small, frameon=False, loc="upper right")
        ax_b.axhline(50, color="k", lw=0.6, ls=":")
        ax_b.set_ylim(0, 100)
        for ax in (ax_a, ax_b, ax_f):
            ax.set_xlim(0, _xmax(res))
            ax.set_xlabel("Source depth below scalp ($d$) [mm]", fontproperties=fonts.regular)
        for ax in (ax_c, ax_d):
            ax.axvline(eta, color="k", lw=0.6, ls=":")
            ax.set_xlim(ETA_GRID[0], ETA_GRID[-1])
            ax.set_xlabel(r"Relative noise level ($\eta$)", fontproperties=fonts.regular)
        ax_c.set_ylim(0, 100)
        finite = np.concatenate([np.r_[res["sizes"][r]["d50_vs_eta"], res["sizes"][r]["sphere_eq"]["d_eq_vs_eta"]]
                                 for r in radii])
        finite = finite[np.isfinite(finite)]
        ax_d.set_ylim(10 * np.ceil(finite.max() / 10) if finite.size else 90, 10)  # depth increases downwards
        ax_d.legend([Line2D([], [], color="k", lw=1.5), Line2D([], [], color="k", lw=0.9, ls=":")],
                    ["realistic head: $d_{50}$", "sphere model: $d_\\mathrm{eq}$ (caps of equal area)"],
                    prop=small, frameon=False, loc="lower right")
        rr = np.array(radii)
        d50 = np.array([res["sizes"][r]["d50"] for r in radii])
        eq = np.array([res["sizes"][r]["sphere_eq"]["d_eq"] for r in radii])
        nominal = np.array([res["sizes"][r]["sphere"]["d_eq"] for r in radii])
        ax_e.plot(rr, d50, "o-", color="k", lw=1.5, ms=4, label="realistic head: $d_{50}$")
        ax_e.plot(rr, eq, "s:", color="0.35", lw=1.2, ms=4, mfc="white",
                  label="sphere model: $d_\\mathrm{eq}$, cap of equal area")
        ax_e.plot(rr, nominal, "^:", color="0.6", lw=1.0, ms=4, mfc="white",
                  label="sphere model: $d_\\mathrm{eq}$, cap of radius ρ")
        for x, v, text in zip(rr, d50, (res["sizes"][r]["d50_text"] for r in radii)):
            if not np.isfinite(v):
                ax_e.annotate(f"$d_{{50}}$ {text}", (x, np.nanmin(np.r_[eq, nominal])), fontproperties=small,
                              ha="center", va="top", xytext=(0, -6), textcoords="offset points")
        ax_e.set_xlim(-1, max(rr) + 1)
        ax_e.set_xticks(rr)
        ax_e.set_xlabel("Patch radius (ρ) [mm]  (0 = dipole)", fontproperties=fonts.regular)
        vals = np.r_[d50, eq, nominal]
        if np.isfinite(vals).any():
            ax_e.set_ylim(5 * np.floor(np.nanmin(vals) / 5 - 0.5), 5 * np.ceil(np.nanmax(vals) / 5 + 0.5))
        ax_e.yaxis.set_major_locator(MultipleLocator(5))
        ax_e.legend(prop=small, frameon=False, loc="lower left")
        ax_f.axhline(1.0, color="k", lw=0.6)
        ax_f.set_yscale("log")
        top = max(2.0, 1.2 * max(np.nanmax(res["sizes"][r]["loss_profile"]["median"]) for r in radii))
        ax_f.set_ylim(0.1, top)
        f_ticks = [t for t in (0.1, 0.2, 0.5, 1, 2, 5, 10) if t <= top]
        ax_f.set_yticks(f_ticks)
        ax_f.set_yticklabels([f"{t:g}" for t in f_ticks])
        ax_f.minorticks_off()
        ax_f.text(0.98, 0.04, "> 1 for deep sources: deep sphere-model sources lie\nnear the centre, where "
                  "tangential dipoles are nearly silent", transform=ax_f.transAxes, ha="right", va="bottom",
                  fontproperties=small, color="0.35")
        _panel(ax_a, "A", r"Median SNR$_\mathrm{OPM}$ / SNR$_\mathrm{SQUID}$", fonts)
        _panel(ax_b, "B", "Cortex where OPM is better [% per depth bin]", fonts)
        _panel(ax_c, "C", "Cortex where OPM is better [% of area]", fonts)
        _panel(ax_d, "D", "Equal-SNR depth vs noise level [mm]", fonts)
        _panel(ax_e, "E", f"Equal-SNR depth at η = {eta:g} [mm]", fonts)
        _panel(ax_f, "F", "OPM signal, realistic / sphere model (median)", fonts)
        fig.text(0.5, 0.985, f"Realistic head ({res['subject']}, single-layer BEM) vs sphere model; "
                 f"η = {eta:g} unless varied.  $d_{{50}}$: depth at which OPM is better for half of the cortex "
                 "at that depth.  Sphere-model caps have the median area of the realistic patches (panel E also "
                 "shows the nominal radius).", ha="center", va="top", fontproperties=small, color="0.3")
        return _save(fig, outdir, "Figure_realistic_summary")


# =========================================================================================
# Driver
# =========================================================================================
def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_tables(res: dict, outdir: Path, files) -> dict:
    raw, radii = res["raw"], sorted(res["sizes"])
    with open(outdir / "realistic_head_sources.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        header = ["hemi", "vertno", "depth_mm", "depth_whole_scalp_mm", "area_mm2"]
        for r in radii:
            header += [f"B_OPM_pT_rho{r:g}", f"B_SQUID_pT_rho{r:g}", f"SNR_ratio_rho{r:g}"]
        wr.writerow(header)
        for i in np.argsort(raw["depth_mm"], kind="stable"):
            row = ["lh" if raw["hemi"][i] == 0 else "rh", int(raw["vertno"][i]), f"{raw['depth_mm'][i]:.3f}",
                   f"{raw['depth_scalp_mm'][i]:.3f}", f"{raw['area_m2'][i] * 1e6:.4f}"]
            for r in radii:
                s = res["sizes"][r]
                row += [f"{s['b_opm'][i] * 1e12:.6g}", f"{s['b_squid'][i] * 1e12:.6g}", f"{s['snr_ratio'][i]:.6g}"]
            wr.writerow(row)
    sizes = {}
    for r in radii:
        s = res["sizes"][r]
        sizes[f"{r:g}"] = dict(
            opm_better_cortex_percent=100 * s["opm_area_fraction"], d50_mm=s["d50"], d50=s["d50_text"],
            sphere_model_d_eq_mm={"cap_of_equal_area": s["sphere_eq"]["d_eq"],
                                  "cap_of_radius_rho": s["sphere"]["d_eq"]},
            equal_area_cap_radius_mm=s["sphere_eq"]["radius_mm"],
            log_snr_ratio_variance_explained_by_depth=s["depth_r2"],
            median_patch_area_cm2=float(np.median(s["patch_area_cm2"])),
            median_sources_per_patch=float(np.median(s["n_vertices"])),
            median_signal_pT={"OPM": float(np.median(s["b_opm"]) * 1e12),
                              "SQUID": float(np.median(s["b_squid"]) * 1e12)},
            max_signal_pT={"OPM": float(s["b_opm"].max() * 1e12), "SQUID": float(s["b_squid"].max() * 1e12)},
            vs_eta={"eta": ETA_GRID, "opm_better_cortex_percent": 100 * s["opm_area_fraction_vs_eta"],
                    "d50_mm": s["d50_vs_eta"], "sphere_d_eq_mm_equal_area": s["sphere_eq"]["d_eq_vs_eta"]})
    summary = dict(
        subject=res["subject"], subjects_dir="SUBJECTS_DIR (local; not recorded)", mne_version=mne.__version__,
        head_model="single-compartment BEM (inner skull, 5120 triangles, 0.3 S/m)",
        source_space=f"ico-5 white-matter surface, {len(raw['depth_mm'])} sources, fixed cortical-normal orientation",
        cortex_area_cm2=float(raw["area_m2"].sum() * 1e4),
        sensors=dict(n_per_array=raw["n_sensors"], spacing_mm=raw["sensor_spacing_mm"],
                     xi_mm={k: v * 1e3 for k, v in ARRAYS.items()},
                     coverage="scalp above the plane through LPA, RPA and nasion + 3 cm",
                     measured_component="normal to the scalp"),
        moment_nAm=R.Q_DIPOLE * 1e9,
        patches=("geodesic radius along the cortical mesh (MNE: shortest paths along mesh edges), "
                 "uniform (area-weighted) current density"),
        eta=res["eta"], sigma_squid_pT=res["sigma_squid_pT"], sigma_opm_pT=res["sigma_opm_pT"],
        depth="distance to the sensor-covered scalp",
        depth_range_mm=[float(raw["depth_mm"].min()), float(raw["depth_mm"].max())],
        d50=("depth at which the area fraction of the cortex at that depth with SNR_OPM > SNR_SQUID drops below "
             "50 % (area-weighted, 2-mm bins with >= 30 sources); null if censored (see 'd50' per size)"),
        sphere_model="tangential dipole / spherical cap of the same area on the brain surface (h = 95 mm, b = 80 mm)",
        simulate_evoked_consistency_max_rel_err=raw["simulate_evoked_rel_err"], simulation_runtime_s=raw["runtime_s"],
        sizes=sizes, outputs=[Path(f).name for f in files])  # relative to the summary
    with open(outdir / "realistic_head_summary.json", "w") as fh:
        json.dump(_json_safe(summary), fh, indent=2, allow_nan=False)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description="Figure 3 of Jas et al. (2026) with a realistic head: FreeSurfer "
                                             "cortex (default fsaverage), BEM forward model (MNE-Python), dipoles "
                                             "and cortical patches.")
    ap.add_argument("--subjects-dir", default=None, help="FreeSurfer SUBJECTS_DIR (default: MNE config / fsaverage)")
    ap.add_argument("--subject", default="fsaverage")
    ap.add_argument("--radii", type=float, nargs="+", default=list(PATCH_RADII_MM),
                    help="patch radii [mm] along the cortex; 0 = single dipole (default: 0 5 10 15 20)")
    ap.add_argument("--eta", type=R._positive_float, default=R.ETA, help="sigma_OPM / sigma_SQUID (default 3)")
    ap.add_argument("--sensor-subdivisions", type=int, default=SENSOR_SUBDIVISIONS, choices=(0, 1, 2),
                    help="scalp-mesh subdivisions for the sensor arrays (3.7 / 1.87 / 0.93 mm; default 1)")
    ap.add_argument("--replot", action="store_true",
                    help=f"re-analyse {RESULTS_FILE} in --outdir (e.g. another --eta) instead of simulating")
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent / "realistic_head_output"))
    ap.add_argument("--artwork-fonts", action="store_true",
                    help="licensed artwork fonts (Myriad Pro; see replicate_figure3.py); default: open fonts")
    ap.add_argument("--font-dir", default=None, help="folder containing MyriadPro-Regular.otf / -Bold.otf (implies --artwork-fonts)")
    args = ap.parse_args(argv)
    if any(not 0 <= r <= 40 for r in args.radii):
        ap.error("patch radii must be between 0 and 40 mm")

    mne.set_log_level("WARNING")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.replot:
        if not (outdir / RESULTS_FILE).is_file():
            ap.error(f"--replot needs {outdir / RESULTS_FILE} from a previous run")
        raw = load_raw(outdir / RESULTS_FILE)
        subjects_dir = args.subjects_dir or raw["subjects_dir"]
        src = load_source_space(subjects_dir, raw["subject"])
        print(f"re-analysing {outdir / RESULTS_FILE} ({raw['subject']}, radii "
              f"{', '.join(f'{r:g}' for r in raw['radii_mm'])} mm; --radii and --sensor-subdivisions are ignored)")
    else:
        subjects_dir = resolve_subjects_dir(args.subjects_dir, args.subject)
        radii = sorted(set(float(r) for r in args.radii))
        anat = load_anatomy(subjects_dir, args.subject, max(radii) * 1e-3)
        raw = simulate(anat, radii, args.sensor_subdivisions)
        save_raw(raw, outdir / RESULTS_FILE)
        src = anat.src
    res = analyse(raw, args.eta)
    fonts = R.setup_fonts(args.font_dir, args.artwork_fonts)
    radii = sorted(res["sizes"])
    files = [str(outdir / RESULTS_FILE)]
    for r in radii:
        stem = "Figure_realistic_fig3" if r == 0 else f"Figure_realistic_fig3_patch_rho_{r:02g}mm"
        files += figure_fig3(res, r, fonts, outdir, stem)
    files += figure_maps(res, src, subjects_dir, fonts, outdir) + figure_summary(res, fonts, outdir)
    summary = write_tables(res, outdir, files)

    print(f"MNE-Python {mne.__version__}, {res['subject']}: {len(res['depth_mm'])} sources "
          f"(depth {res['depth_mm'].min():.1f}-{res['depth_mm'].max():.1f} mm), single-layer BEM, "
          f"{res['n_sensors']} sensors per array ({raw['sensor_spacing_mm']:.2f} mm); simulate_evoked consistency "
          f"{raw['simulate_evoked_rel_err']:.1e}; eta = {args.eta:g}; simulation {raw['runtime_s']:.0f} s")
    print(f"{'source':>14s} | {'median area':>11s} | {'OPM better':>10s} | {'d50':>9s} | "
          f"{'sphere d_eq, equal area / radius rho':>36s}")
    for r in radii:
        st = summary["sizes"][f"{r:g}"]
        sph = st["sphere_model_d_eq_mm"]
        print(f"{_size_label(r):>14s} | {st['median_patch_area_cm2']:7.2f} cm2 | "
              f"{st['opm_better_cortex_percent']:8.1f} % | {st['d50']:>9s} | "
              f"{_fmt_mm(sph['cap_of_equal_area']):>17s} / {_fmt_mm(sph['cap_of_radius_rho']):>16s}")
    print("Wrote:", outdir)
    return res, summary


def _fmt_mm(value) -> str:
    return "none" if value is None or not np.isfinite(value) else f"{value:.2f} mm"


if __name__ == "__main__":
    main()
