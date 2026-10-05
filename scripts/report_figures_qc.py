#!/usr/bin/env python3
"""Supplementary figure of the MRI quality check of the school-aged children (referee round 1), drawn
in the report's style from the check's own geometry.

  Figure_S_children_qc  (a) For the adult, children A-C and the 2-year template, a 50 x 50 mm section
      through the white-surface vertex closest to the scalp used (the size of the check's zoom panels;
      sagittal or coronal, whichever holds more of the scalp's outward normal at that point, so that the
      section crosses the layers): the T1, the dense MRI scalp used, the MRI head boundary, the outer and
      inner skull of the BEM used and the white surface. (b) Cumulative share of white-surface vertices
      against their distance to the scalp used and to the MRI head boundary, lowest quarter (as the
      check's summary figure). (c) Offset of the MRI head boundary from the scalp used over the head above
      the fiducial plane (the check's cap): median and quartiles.

Nothing is fitted or chosen here beyond the section plane: every surface and distance comes from the
functions of scripts/study_children_qc.py, imported and called as its analyse() calls them: the
anatomies as G3B and G4 load them (load), the MRI head boundary along the scalp's smoothed outward
normals (the steps of scalp_vs_mri: sample_along, otsu_threshold, head_boundary in two passes with
consensus; the Otsu crossing is the check's primary edge; see boundary_offsets) and as a surface
(corrected_scalp: the scalp moved onto it, the surface the check measures the cortex against), the exact
point-to-mesh distance (SignedDistance), the distance records (near_surface), the offset summaries
(summarise_offsets) and the plane sections (plane_segments). The script stops unless the records it
recomputes equal those stored in results/g3b_children_qc/children_qc.json (relative tolerance 1e-9): for
every anatomy the scalp's vertex counts, the closest white vertex (distance, position, parcel), the full
distance records to the scalp used and to the MRI head boundary (percentiles, counts, areas and parcels
below 8 and 10 mm, sector percentiles), the offset summaries over the cap and the whole head (Otsu,
half-maximum and steepest edges), the Otsu and air levels, and the MRI-boundary surface's construction
report; and the curves of (b) must pass through the stored counts below 8 and 10 mm. Titles, labels and
caption numbers are read from the stored file after that check.

Inputs: results/g3b_children_qc/children_qc.json, configs/g3b_pediatric.toml, the five anatomies under
$OPMSQUID_DATA/external (MNE sample subject; infant template ANTS2-0Years3T; OpenNeuro ds005234
sub-Z213, sub-Z209, sub-Z226 = children A-C) with their T1 volumes, and their full-resolution cortices
($OPMSQUID_CACHE/anatomy; built there by opmsquid.anatomy if missing).
Outputs (results/report/): Figure_S_children_qc.png and figures_qc.json (inputs, description, alt text,
caption draft, plotted values and the checks).
About 5 min and 1.3 GB of memory (mostly the exact distances of about 300,000 white-surface vertices to
two scalp meshes per anatomy, as the check computes them).

  PYTHONPATH=src .venv/bin/python scripts/report_figures_qc.py
"""
from __future__ import annotations

import json
import math
import re
import resource
import sys
import textwrap
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import report_style as style  # noqa: E402  (selects the Agg backend before pyplot is used)
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from scipy import ndimage  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import study_children_qc as qc  # noqa: E402  (the check: loading, boundary, distances, sections)
from opmsquid import anatomy, io  # noqa: E402

QC_JSON = "results/g3b_children_qc/children_qc.json"
CONFIG = "configs/g3b_pediatric.toml"
NAME = "Figure_S_children_qc"
OUT_JSON = "figures_qc.json"
KEYS = ("adult", "childA", "childB", "childC", "infant2yr")  # the check's anatomies: the children between the references
LABEL = {**style.ANAT_LABEL, "infant2yr": "2-year template"}  # the manuscript's name of the 24-month template
HALF_MM = 25.0  # half-width of the sections: the check's zoom panels (study_children_qc.overlay_figure)
PIXEL_MM = 0.2  # T1 sampling of the sections, as the check's zoom panels
CDF_STEP_MM = 0.02  # (b) drawn on this grid (exact counts at every grid point)
X_MAX_MM, SHARE_MAX = 25.0, 0.25  # (b): the lowest quarter of the distribution, as the check's summary figure
TOL = 1e-9  # relative: same code, data and libraries, so the recomputation must reproduce the stored values
BLOCK = 40000  # scalp profiles per call of the check's head_boundary (memory only; see boundary_offsets)
# (a): surface, legend label, colour (Okabe-Ito, none of the head colours of (b) and (c)), line width [pt]
SURFACES = (("scalp", "scalp used (dense MRI scalp)", "#F0E442", 1.1),
            ("mri", "MRI head boundary (T1 edge)", "#CC79A7", 1.1),
            ("outer_skull", "outer skull (BEM used)", "#56B4E9", 0.9),
            ("inner_skull", "inner skull (BEM used)", "#E69F00", 0.9),
            ("white", "white surface", "#0072B2", 0.8))
PLANES = {0: ("sagittal", (1, 2)), 1: ("coronal", (0, 2))}  # section axis -> name, in-plane axes (MRI frame)
DK = dict(  # Desikan-Killiany parcels (FreeSurfer aparc) in words
    bankssts="banks of the superior temporal sulcus", caudalanteriorcingulate="caudal anterior cingulate",
    caudalmiddlefrontal="caudal middle frontal", cuneus="cuneus", entorhinal="entorhinal", fusiform="fusiform",
    inferiorparietal="inferior parietal", inferiortemporal="inferior temporal", isthmuscingulate="isthmus cingulate",
    lateraloccipital="lateral occipital", lateralorbitofrontal="lateral orbitofrontal", lingual="lingual",
    medialorbitofrontal="medial orbitofrontal", middletemporal="middle temporal", parahippocampal="parahippocampal",
    paracentral="paracentral", parsopercularis="pars opercularis", parsorbitalis="pars orbitalis",
    parstriangularis="pars triangularis", pericalcarine="pericalcarine", postcentral="postcentral",
    posteriorcingulate="posterior cingulate", precentral="precentral", precuneus="precuneus",
    rostralanteriorcingulate="rostral anterior cingulate", rostralmiddlefrontal="rostral middle frontal",
    superiorfrontal="superior frontal", superiorparietal="superior parietal", superiortemporal="superior temporal",
    supramarginal="supramarginal", frontalpole="frontal pole", temporalpole="temporal pole",
    transversetemporal="transverse temporal", insula="insula")
TEXT_NAME = dict(adult="adult", childA="child A", childB="child B", childC="child C", infant2yr="2-year template")
MINUS, EN, TIMES, PM = "−", "–", "×", "±"


# ----------------------------------------------------------------------------------------------
# helpers
def need(ok: bool, what: str) -> None:
    if not ok:
        raise ValueError(what)


def close(a, b, tol: float = TOL) -> bool:
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a)), abs(float(b)))


def compare(new, old, path: str) -> int:
    """Recomputed record ``new`` against the stored ``old`` (as io.write_json stored it: NaN as null),
    value by value; the number of values compared. Stops at the first difference."""
    if isinstance(old, dict):
        need(isinstance(new, dict) and set(new) == set(old),
             f"{path}: keys differ ({sorted(set(new) ^ set(old)) if isinstance(new, dict) else type(new).__name__})")
        return sum(compare(new[k], old[k], f"{path}.{k}") for k in old)
    if isinstance(old, list):
        need(isinstance(new, (list, tuple, np.ndarray)) and len(new) == len(old), f"{path}: length differs")
        return sum(compare(a, b, f"{path}[{i}]") for i, (a, b) in enumerate(zip(new, old)))
    if old is None:
        need(new is None or (isinstance(new, float) and math.isnan(new)), f"{path}: recomputed {new!r}, stored null")
        return 1
    if isinstance(old, (str, bool)):
        need(new == old and type(new) is type(old), f"{path}: recomputed {new!r}, stored {old!r}")
        return 1
    need(isinstance(new, (int, float, np.integer, np.floating)) and not isinstance(new, bool) and close(new, old),
         f"{path}: recomputed {new!r}, stored {old!r}")
    return 1


def signed(x: float, nd: int = 2) -> str:
    s = f"{x:+.{nd}f}"
    return s.replace("-", MINUS)


def region_words(region: str) -> str:
    """'rh.lateraloccipital' -> 'right lateral occipital'."""
    hemi, name = region.split(".", 1)
    return f"{dict(lh='left', rh='right')[hemi]} {DK[name]}"


# ----------------------------------------------------------------------------------------------
# the check's geometry, recomputed
def boundary_offsets(vol: dict, rr_m: np.ndarray, nn: np.ndarray, cap: np.ndarray, tris: np.ndarray) -> dict:
    """study_children_qc.scalp_vs_mri's edges, Otsu and air levels, with its head_boundary applied to
    blocks of BLOCK profiles: head_boundary treats every profile on its own, so the offsets are the same
    to the bit, and the temporaries of the 267,000-profile adult scalp stay well below 1 GB (all
    profiles at once take the process above 3 GB). The check below compares them with the stored ones."""
    t = np.arange(qc.PROFILE_MM[0], qc.PROFILE_MM[1] + qc.PROFILE_MM[2] / 2, qc.PROFILE_MM[2])
    otsu = qc.otsu_threshold(vol["data"].ravel())
    prof = qc.sample_along(vol, rr_m, nn, t)
    air = float(np.median(prof[np.ix_(cap, t >= qc.AIR_FROM_MM)]))  # the cap's samples far outside, as scalp_vs_mri
    adj, deg = qc.adjacency(np.asarray(tris), len(rr_m))
    blocks = [slice(i, i + BLOCK) for i in range(0, len(prof), BLOCK)]
    first = np.concatenate([qc.head_boundary(prof[b], t, otsu, air)[qc.PRIMARY] for b in blocks])
    reference = qc.consensus(first, adj, deg)  # second pass: the fall nearest to the neighbourhood's first-pass offset
    parts = [qc.head_boundary(prof[b], t, otsu, air, reference=reference[b]) for b in blocks]
    edge = {k: np.concatenate([p[k] for p in parts]) for k in ("half_max", "otsu", "steepest")}
    return dict(edge=edge, otsu=otsu, air=air)


def measure(key: str, cfg: dict) -> dict:
    """One anatomy as the check measures it (the calls of study_children_qc.analyse, scalp_vs_mri's over
    blocks of profiles), the records it stores, and the section through the closest white vertex."""
    t0 = time.time()
    subject, t1_path, _, _, _ = qc.load(key, cfg)
    cortex = anatomy.full_resolution(subject)
    regions = qc.region_names(subject, cortex)
    scalp = subject.scalp
    tris = np.asarray(scalp.tris)
    adj, deg = qc.adjacency(tris, len(scalp.rr))
    nn_s = qc.smooth_vectors(scalp.nn, adj, deg)
    height = qc.head_frame(subject, scalp.rr)[:, 2]
    cap = height >= 0
    d_used = qc.SignedDistance(scalp.rr, scalp.tris)(cortex.rr) * 1e3
    hf_cortex = qc.head_frame(subject, cortex.rr)
    closest = int(np.argmin(d_used))
    t1 = qc.read_mgh(t1_path)
    res = boundary_offsets(t1, scalp.rr, nn_s, cap, scalp.tris)
    edge = res["edge"]
    mri, mri_rep = qc.corrected_scalp(scalp, edge[qc.PRIMARY], nn_s, cap, qc.taper_weight(height))
    d_mri = qc.SignedDistance(mri.rr, mri.tris)(cortex.rr) * 1e3
    rec = dict(scalp=dict(n_vertices=int(len(scalp.rr)), cap_vertices=int(cap.sum())),
               closest_white_vertex=dict(region=str(regions[closest]), distance_mm=float(d_used[closest]),
                                         mri_frame_mm=(cortex.rr[closest] * 1e3).tolist()),
               white_to_scalp_used=qc.near_surface(d_used, cortex, regions, hf_cortex),
               white_to_mri_boundary=qc.near_surface(d_mri, cortex, regions, hf_cortex),
               scalp_vs_mri=dict(otsu_threshold=res["otsu"], air_level=res["air"], offsets=qc.summarise_offsets(edge, cap)),
               mri_boundary_surface=mri_rep)

    # the section: through the closest white vertex, normal to the MRI axis (x: sagittal, y: coronal) along which
    # the scalp's smoothed outward normal at its vertex nearest to the white vertex has the smaller component
    j = int(cKDTree(scalp.rr).query(cortex.rr[closest])[1])
    normal = nn_s[j]
    axis = 0 if abs(normal[0]) < abs(normal[1]) else 1
    plane, (u, v) = PLANES[axis]
    centre = cortex.rr[closest] * 1e3
    value = float(centre[axis])
    gu = centre[u] + np.arange(-HALF_MM, HALF_MM + PIXEL_MM / 2, PIXEL_MM)
    gv = centre[v] + np.arange(-HALF_MM, HALF_MM + PIXEL_MM / 2, PIXEL_MM)
    uu, vv = np.meshgrid(gu, gv)
    p = np.zeros(uu.shape + (3,))
    p[..., u], p[..., v], p[..., axis] = uu, vv, value
    inv = np.linalg.inv(t1["vox2ras_tkr"])  # MRI frame [mm] -> voxel, as the check's overlays
    image = ndimage.map_coordinates(t1["data"], (p.reshape(-1, 3) @ inv[:3, :3].T + inv[:3, 3]).T, order=1,
                                    mode="constant").reshape(uu.shape)
    bem = {b["id"]: b for b in subject.bem_surfaces}
    meshes = dict(scalp=(scalp.rr, scalp.tris), mri=(mri.rr, mri.tris),
                  outer_skull=(bem[FIFF.FIFFV_BEM_SURF_ID_SKULL]["rr"], bem[FIFF.FIFFV_BEM_SURF_ID_SKULL]["tris"]),
                  inner_skull=(bem[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["rr"], bem[FIFF.FIFFV_BEM_SURF_ID_BRAIN]["tris"]),
                  white=(cortex.rr, cortex.tris))
    lo, hi = np.array([gu[0], gv[0]]) - 1.0, np.array([gu[-1], gv[-1]]) + 1.0
    segments = {}
    for name, (rr, t) in meshes.items():
        s = qc.plane_segments(np.asarray(rr, float) * 1e3, np.asarray(t), axis, value)[:, :, [u, v]]
        segments[name] = s[np.all((s >= lo) & (s <= hi), axis=2).any(axis=1)]  # those that reach the section
    section = dict(plane=plane, axis="xyz"[axis], coordinate_mm=value, in_plane_axes=["xyz"[u], "xyz"[v]],
                   centre=(float(centre[u]), float(centre[v])), image=image.astype(np.float32),
                   extent=(gu[0] - PIXEL_MM / 2, gu[-1] + PIXEL_MM / 2, gv[0] - PIXEL_MM / 2, gv[-1] + PIXEL_MM / 2),
                   vmax=float(np.percentile(t1["data"][t1["data"] > 0], 99.5)),  # the check's grey scale
                   segments=segments, scalp_normal=normal.tolist(),
                   normal_in_plane_share=float(np.sqrt(1.0 - normal[axis] ** 2)))
    box = np.all(np.abs(scalp.rr * 1e3 - centre) <= HALF_MM + 3.0, axis=1)  # the scalp in and near the section's cube
    section["scalp_height_in_section_min_mm"] = float(height[box].min() * 1e3)
    # (b): share of the white vertices closer than x (strictly, as the check's counts below 8 and 10 mm)
    grid = np.arange(0.0, X_MAX_MM + CDF_STEP_MM / 2, CDF_STEP_MM)
    cdf = {w: np.searchsorted(np.sort(d), grid, side="left") / len(d) for w, d in (("used", d_used), ("mri", d_mri))}
    out = dict(rec=rec, section=section, grid=grid, cdf=cdf, n_white=int(cortex.n), subject=subject.name,
               runtime_s=time.time() - t0)
    print(f"{key}: measured in {out['runtime_s']:.0f} s", flush=True)
    return out


def check(M: dict, J: dict, cfg: dict) -> dict:
    """Every recomputed record against the stored one; stops at the first difference."""
    for c in cfg["anatomy"]["school"]:  # the children's labels carry their configured ages
        need(LABEL[c["key"]].endswith(f"({c['age_y']:.1f} y)"), f"{c['key']}: label {LABEL[c['key']]!r}, age {c['age_y']}")
    counts = {}
    for key, m in M.items():
        st = J["anatomies"][key]
        need(st["t1_available"], f"{key}: the check ran without its T1")
        sv = st["scalp_vs_mri"]
        need(sv["primary"] == qc.PRIMARY, f"{key}: stored primary edge {sv['primary']}, the check's is {qc.PRIMARY}")
        stored = dict(scalp=st["scalp"], closest_white_vertex=st["closest_white_vertex"],
                      white_to_scalp_used=st["white_to_scalp_used"], white_to_mri_boundary=st["white_to_mri_boundary"],
                      scalp_vs_mri={k: sv[k] for k in ("otsu_threshold", "air_level", "offsets")},
                      mri_boundary_surface={k: v for k, v in st["mri_boundary_surface"].items()
                                            if k != "remeasured_cap_offset_mm"})  # the closure: not recomputed
        counts[key] = {name: compare(m["rec"][name], stored[name], f"{key}.{name}") for name in stored}
        # the curves of (b) pass through the stored counts below 8 and 10 mm
        for w, rk in (("used", "white_to_scalp_used"), ("mri", "white_to_mri_boundary")):
            for thr in J["parameters"]["near_mm"]:
                i = int(round(thr / CDF_STEP_MM))
                need(close(m["grid"][i], thr), f"{thr} mm is not on the grid of (b)")
                need(int(round(m["cdf"][w][i] * m["n_white"])) == st[rk][f"below_{thr:g}mm"]["n_vertices"],
                     f"{key}: the curve of (b) ({w}) at {thr:g} mm does not give the stored count")
                counts[key][f"curve_{w}_at_{thr:g}mm"] = 1
        # the description says the MRI-boundary surface follows the measured offsets in every section
        need(m["section"]["scalp_height_in_section_min_mm"] >= 0.0,
             f"{key}: the section reaches below the fiducial plane, where the MRI-boundary surface is extrapolated")
    return counts


# ----------------------------------------------------------------------------------------------
# the figure
def draw(M: dict, J: dict) -> dict:
    style.apply()
    A = J["anatomies"]
    near = [float(x) for x in J["parameters"]["near_mm"]]
    tol = float(J["parameters"]["tolerance_mm"])
    W, H = style.FULL_W, 6.05
    fig = plt.figure(figsize=(W, H))

    def rect(x, top, w, h):  # inches from the left edge and from the top edge -> figure fraction
        return [x / W, (H - top - h) / H, w / W, h / H]

    # (a) sections
    n, s, x0 = len(KEYS), 1.36, 0.02
    gap = (W - 2 * x0 - n * s) / (n - 1)
    top_a = 0.70
    fig.text(x0 / W, 1 - 0.02 / H, f"(a)  Sections of {2 * HALF_MM:g} {TIMES} {2 * HALF_MM:g} mm through the white-surface "
             "vertex closest to the scalp used", ha="left", va="top", fontsize=8.5)
    drawn = {}
    for i, key in enumerate(KEYS):
        sec = M[key]["section"]
        ax = fig.add_axes(rect(x0 + i * (s + gap), top_a, s, s))
        ax.imshow(sec["image"], cmap="gray", origin="lower", extent=sec["extent"], vmin=0, vmax=sec["vmax"],
                  interpolation="antialiased")
        for name, _, color, lw in SURFACES:
            ax.add_collection(LineCollection(sec["segments"][name], colors=color, linewidths=lw, capstyle="round"))
        cu, cv = sec["centre"]
        ax.plot(cu, cv, "o", ms=5.4, mfc="none", mec="black", mew=1.7)  # a white ring outlined in black: visible on
        ax.plot(cu, cv, "o", ms=5.4, mfc="none", mec="white", mew=0.8)  # the bright white matter and on the dark CSF
        ax.set_xlim(cu - HALF_MM, cu + HALF_MM)
        ax.set_ylim(cv - HALF_MM, cv + HALF_MM)
        ax.set_aspect("equal")
        ticks = np.arange(-20.0, 20.5, 10.0)
        ax.set_xticks(cu + ticks, [])
        ax.set_yticks(cv + ticks, [])
        ax.tick_params(length=2.2, width=0.6, direction="out")
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_linewidth(0.6)
        ax.text(0.035, 0.035, sec["plane"], transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, color="white")
        d = A[key]["closest_white_vertex"]["distance_mm"]
        ax.set_title(f"{LABEL[key]}\nclosest white vertex\n{d:.1f} mm below the scalp", fontsize=7.8, pad=3.5,
                     linespacing=1.2)
        drawn[key] = dict(plane=sec["plane"], section_axis=sec["axis"], section_coordinate_mm=sec["coordinate_mm"],
                          centre_mm=dict(zip(sec["in_plane_axes"], sec["centre"])),
                          scalp_normal_at_nearest_vertex=sec["scalp_normal"],
                          normal_share_in_section=sec["normal_in_plane_share"],
                          segments_drawn={k: len(v) for k, v in sec["segments"].items()}, title_distance_mm=d,
                          region=A[key]["closest_white_vertex"]["region"])
    handles = [Line2D([], [], color=c, lw=1.8, label=lab) for _, lab, c, _ in SURFACES]
    handles.append(Line2D([], [], ls="none", marker="o", ms=5.4, mfc="none", mec="black", mew=1.2, label="closest white vertex"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1 - (top_a + s + 0.06) / H), ncol=3,
               fontsize=7.3, handlelength=1.9, columnspacing=2.0, borderaxespad=0.0)

    # (b) and (c)
    top_bc = top_a + s + 0.66  # their headings
    top_ax = top_bc + 0.52
    h_ax = 1.85
    xb, wb = 0.52, 3.2
    xc, wc = 4.78, 2.37
    fig.text(0.02 / W, 1 - top_bc / H, "(b)  White-surface vertices near the head surface", ha="left", va="top", fontsize=8.5)
    fig.text((xc - 1.08) / W, 1 - top_bc / H, "(c)  MRI head boundary against the scalp used", ha="left", va="top",
             fontsize=8.5)
    heads = [Line2D([], [], color=style.ANAT_COLOR[k], lw=1.3, marker=style.CLASS_MARKER[style.ANAT_CLASS[k]], ms=4.2,
                    label=LABEL[k]) for k in KEYS]
    fig.legend(handles=heads, loc="upper center", bbox_to_anchor=(0.5, 1 - (top_bc + 0.2) / H), ncol=len(KEYS),
               fontsize=7.3, handlelength=2.0, columnspacing=1.5, borderaxespad=0.0)

    axb = fig.add_axes(rect(xb, top_ax, wb, h_ax))
    curves = {}
    for key in KEYS:
        m, c = M[key], style.ANAT_COLOR[key]
        axb.plot(m["grid"], 100 * m["cdf"]["used"], color=c, lw=1.25, solid_capstyle="butt")
        axb.plot(m["grid"], 100 * m["cdf"]["mri"], color=c, lw=1.25, ls=(0, (3.2, 1.6)))
        step = int(round(0.5 / CDF_STEP_MM))
        curves[key] = {w: [float(x) for x in m["cdf"][w][::step]] for w in ("used", "mri")}
    for thr in near:
        axb.axvline(thr, color="0.55", lw=0.7, ls=":", zorder=0)
        axb.text(thr, 100 * SHARE_MAX * 1.01, f"{thr:g}", ha="center", va="bottom", fontsize=6.8, color="0.35")
    axb.set_xlim(0, X_MAX_MM)
    axb.set_ylim(0, 100 * SHARE_MAX)
    axb.set_xticks(np.arange(0, X_MAX_MM + 0.1, 5))
    axb.set_xlabel("Distance, white-surface vertex to head surface (mm)")
    axb.set_ylabel("White-surface vertices\n(cumulative %)")
    style_handles = [Line2D([], [], color="0.2", lw=1.25, label="scalp used"),
                     Line2D([], [], color="0.2", lw=1.25, ls=(0, (3.2, 1.6)), label="MRI head boundary")]
    axb.legend(handles=style_handles, loc="upper left", fontsize=7.2, handlelength=2.2, borderaxespad=0.2, frameon=True,
               facecolor="white", edgecolor="none", framealpha=1.0, borderpad=0.25)  # over the dotted 8- and 10-mm lines

    axc = fig.add_axes(rect(xc, top_ax, wc, h_ax))
    ref = A["adult"]["scalp_vs_mri"]["offsets"]["cap"][qc.PRIMARY]["median_mm"]
    axc.axvspan(ref - tol, ref + tol, color="0.9", lw=0, zorder=0)
    axc.axvline(0, color="0.45", lw=0.7, zorder=1)
    offsets = {}
    for i, key in enumerate(KEYS):
        o = A[key]["scalp_vs_mri"]["offsets"]["cap"][qc.PRIMARY]
        c = style.ANAT_COLOR[key]
        axc.plot([o["p25_mm"], o["p75_mm"]], [-i, -i], color=c, lw=6.0, alpha=0.5, solid_capstyle="butt", zorder=2)
        axc.plot(o["median_mm"], -i, ls="none", marker=style.CLASS_MARKER[style.ANAT_CLASS[key]], ms=5.0, color=c,
                 mec="white", mew=0.6, zorder=3)
        axc.text(o["p75_mm"] + 0.15, -i, signed(o["median_mm"]), ha="left", va="center", fontsize=7.0, color="0.15")
        offsets[key] = dict(p25_mm=o["p25_mm"], median_mm=o["median_mm"], p75_mm=o["p75_mm"])
    axc.text(ref, 0.62, f"adult {PM} {tol:g} mm", ha="center", va="bottom", fontsize=6.8, color="0.35")
    axc.set_yticks(-np.arange(len(KEYS)), [LABEL[k] for k in KEYS])
    axc.tick_params(axis="y", length=0)
    axc.set_ylim(-len(KEYS) + 0.45, 0.95)
    xlim = (-2.6, 2.9)
    need(all(xlim[0] < o["p25_mm"] and o["p75_mm"] + 0.75 < xlim[1] for o in offsets.values()), "(c): an offset outside the axis")
    axc.set_xlim(*xlim)
    axc.set_xticks([-2, -1, 0, 1, 2])
    axc.spines["left"].set_visible(False)
    axc.set_xlabel("MRI head boundary outside\nthe scalp used (mm)")

    foot = (f"Sections: superior up; sagittal with anterior to the right, coronal with the subject's right to the right; "
            f"ticks every 10 mm. MRI head boundary: where the T1 falls through the volume's air/tissue (Otsu) level along "
            f"the scalp's outward normals; in (a) and (b) the scalp moved onto it. "
            f"(c) Over the head above the fiducial plane; bars: quartiles; shaded: the adult's median {PM} {tol:g} mm. "
            f"Data: OpenNeuro ds005234 (Fadeev et al., 2024) for children A{EN}C; 2-year template: O'Reilly et al. (2021); "
            f"adult: MNE sample data.")
    fig.text(0.02 / W, 1 - (top_ax + h_ax + 0.62) / H, textwrap.fill(foot, 138), ha="left", va="top", fontsize=6.6,
             color="0.25")
    # what a reader sees: no internal keys, milestone codes or file names; the data acknowledgement
    shown = [t.get_text() for t in fig.findobj(lambda a: hasattr(a, "get_text") and a.get_visible())]
    internal = [t for t in shown if re.search(r"child[ABC]\b|infant\d|sub-Z|\bG\d[A-Z]?\b|squid|opm_|_mm\b|\.mgz|\.fif", t)]
    need(not internal, f"internal labels in the figure: {internal}")
    need(any("Data: OpenNeuro ds005234 (Fadeev et al., 2024)" in " ".join(t.split()) for t in shown),
         "footnote: the data acknowledgement is missing")
    path = style.save(fig, NAME)
    print(path.relative_to(ROOT))
    return dict(sections=drawn, cumulative_share_at_0p5mm=curves, offsets_cap_otsu=offsets, adult_reference_mm=ref,
                tolerance_mm=tol, footnote=foot)


# ----------------------------------------------------------------------------------------------
# provenance and caption
def entry(M: dict, J: dict, cfg: dict, drawn: dict, counts: dict) -> dict:
    A, P = J["anatomies"], J["parameters"]
    skull = float(cfg["anatomy"]["child_skull_depth_mm"])
    near = [float(x) for x in P["near_mm"]]
    tol = float(P["tolerance_mm"])
    kids = ("childA", "childB", "childC")
    cw = {k: A[k]["closest_white_vertex"] for k in KEYS}
    off = {k: A[k]["scalp_vs_mri"]["offsets"]["cap"] for k in KEYS}
    used8 = {k: A[k]["white_to_scalp_used"][f"below_{near[0]:g}mm"]["area_cm2"] for k in KEYS}
    mri8 = {k: A[k]["white_to_mri_boundary"][f"below_{near[0]:g}mm"]["area_cm2"] for k in KEYS}
    verdict = {k: A[k]["verdict"]["class"] for k in KEYS}
    need(all(verdict[k] == "usable" for k in ("adult",) + kids) and verdict["infant2yr"] == "outside",
         f"verdicts {verdict}: the caption's wording assumes children and adult usable, the template outside")
    need(used8["adult"] == 0 and used8["infant2yr"] == 0 and mri8["adult"] == 0,
         "caption: the adult or the template has white surface within 8 mm")
    need(near[0] == skull, f"caption: the check's first threshold ({near[0]:g} mm) is not the modelled skull depth "
                           f"({skull:g} mm)")
    sections = "; ".join(f"{TEXT_NAME[k]}: {region_words(cw[k]['region'])}, {cw[k]['distance_mm']:.1f} mm, "
                         f"{drawn['sections'][k]['plane']}" for k in KEYS)
    med = {k: off[k][qc.PRIMARY]["median_mm"] for k in KEYS}
    half = {k: off[k]["half_max"]["median_mm"] for k in KEYS}
    caption = (
        "MRI quality check of the anatomy of the three school-aged children, with the adult and the 2-year (24-month) infant "
        "template as references; the scalp, skull and white surfaces as the forward models use them. "
        f"(a) Sections of {2 * HALF_MM:g} {TIMES} {2 * HALF_MM:g} mm through the white-surface vertex closest to the scalp used "
        f"(open circle; {sections}), sagittal or coronal, whichever contains more of the scalp's outward normal there, so that "
        "the section crosses the layers. Grey: T1 intensity. Yellow: the dense MRI scalp used for sensor placement and source "
        "depths. Pink: the MRI head boundary, where the T1 falls through the volume's air/tissue (Otsu) level along the "
        "scalp's outward normals, drawn as the scalp moved onto it. Light blue and orange: the outer and inner skull of the "
        f"boundary-element model used (for the children modelled: the inner skull {skull:g} mm below the scalp where the "
        "dataset's lay shallower, the outer skull about halfway to the scalp). Dark blue: the white surface. "
        "(b) Cumulative share of white-surface vertices against their distance to the scalp used (solid) and to the MRI head "
        f"boundary (dashed), lowest quarter; dotted lines at the check's thresholds, {near[0]:g} mm (also the depth of the "
        f"children's modelled inner skull) and {near[1]:g} mm. Within {near[0]:g} mm of the "
        f"scalp used lie {used8['childA']:.1f}, {used8['childB']:.1f} and {used8['childC']:.1f} cm² of white surface in "
        f"children A{EN}C (adult and template: none), within {near[0]:g} mm of the MRI head boundary {mri8['childA']:.1f}, "
        f"{mri8['childB']:.1f} and {mri8['childC']:.1f} cm² (template {mri8['infant2yr']:.1f} cm², adult none). "
        "(c) Offset of the MRI head boundary from the scalp used over the head above the fiducial plane, along the scalp's "
        "outward normals: median (marker and value) and quartiles (bar); positive where the boundary lies outside the scalp "
        f"used. Shaded: the adult's median {PM} {tol:g} mm, the tolerance of the check's verdict (one voxel). The "
        f"children's scalps lie a median {med['childA']:.2f}, {med['childB']:.2f} and {med['childC']:.2f} mm inside their MRI "
        f"head boundary, as the adult's does ({med['adult']:.2f} mm; at the half-maximum edge {half['childA']:.2f}, "
        f"{half['childB']:.2f}, {half['childC']:.2f} and {half['adult']:.2f} mm), and the check rates their surfaces usable "
        f"(frames, white-surface registration and scalp within its tolerances); the template's scalp lies "
        f"{-med['infant2yr']:.2f} mm outside its boundary (verdict: outside). The children's shallow cortex is therefore not "
        "the product of a misplaced scalp. Children A–C: OpenNeuro ds005234 (Fadeev et al., 2024; de-identified, data "
        "licence CC0); "
        "2-year template: O'Reilly et al. (2021), built from the Neurodevelopmental MRI Database (Richards et al., 2016); "
        "adult: MNE sample subject.")
    alt = ("Three-part figure. (a) Five small grey-scale MRI sections, each 50 mm wide, of the adult, children A to C and the "
           "2-year template, each centred on the white-surface vertex closest to the scalp, with coloured outlines of the scalp "
           "used, the MRI head boundary just outside it (inside it for the template), the outer and inner skull and the white "
           f"surface; the children's closest vertices lie {cw['childA']['distance_mm']:.1f} to "
           f"{max(cw[k]['distance_mm'] for k in kids):.1f} mm below the scalp, the adult's "
           f"{cw['adult']['distance_mm']:.1f} mm. (b) Line chart of the cumulative share of white-surface vertices (0 to 25 %) "
           "against distance (0 to 25 mm): the children's curves rise first, child B's earliest; the dashed curves (to the MRI "
           "head boundary) lie to the right of the solid ones (deeper) for the adult and children and to the left for the "
           "template. (c) Dot-and-bar chart of the offset of the MRI head boundary from the scalp used: adult and children at "
           f"{signed(min(med[k] for k in ('adult',) + kids))} to {signed(max(med[k] for k in ('adult',) + kids))} mm, the "
           f"template {signed(med['infant2yr'])} mm.")
    description = (
        "Clean redraw of the MRI quality check of the school-aged children (results/g3b_children_qc/Figure_QC_overlay_*.png "
        "and Figure_QC_summary.png, which carry internal labels and whole MRI slices), recomputed with the check's own "
        "functions (scripts/study_children_qc.py: load; sample_along, otsu_threshold, head_boundary and consensus, the steps "
        "of scalp_vs_mri, over blocks of profiles; corrected_scalp, SignedDistance, near_surface, summarise_offsets, "
        "plane_segments) and checked equal to its stored records (see checks). (a) The section planes "
        "are this figure's only choice: through the closest white vertex, normal to the MRI-frame axis (x: sagittal, y: "
        "coronal) along which the scalp's smoothed outward normal at the scalp vertex nearest to that white vertex has the "
        "smaller component; zoom and T1 grey scale as the check's zoom panels (50 mm, 0.2-mm trilinear sampling, 0 to the "
        "99.5th percentile of the non-zero voxels). The MRI head boundary is drawn as the check's surface fitted to it "
        "(corrected_scalp), the surface (b) measures against; every section lies above the fiducial plane, where that "
        f"surface follows the measured offsets (cleaned and smoothed over {qc.SMOOTH_RINGS} neighbour rings). (b) The "
        "empirical cumulative "
        f"share (vertices closer than the distance, as the check's counts below {near[0]:g} and {near[1]:g} mm, through "
        f"which the curves pass) on a {CDF_STEP_MM:g}-mm grid; the y range is the check's summary "
        f"figure's (0 to {SHARE_MAX:g}). (c) The stored quartiles and median of the primary (Otsu) edge over the cap.")
    inputs = [
        f"{QC_JSON} :: anatomies.<adult|childA|childB|childC|infant2yr>.(scalp, closest_white_vertex, white_to_scalp_used, "
        "white_to_mri_boundary, scalp_vs_mri.(primary, otsu_threshold, air_level, offsets), mri_boundary_surface) "
        "(recomputed and compared), .verdict.class, .t1_available; parameters.(near_mm, tolerance_mm)",
        f"{CONFIG} :: anatomy.school (the children's subjects and BEM folders, via study_children_qc.load)",
        "$OPMSQUID_DATA/external: MNE-sample-data/subjects/sample (T1.mgz, surf/*.white, bem: sample-head.fif, "
        "sample-5120-5120-5120-bem.fif, sample-oct-6-src.fif, label/*.aparc.annot), MNE-sample-data/MEG/sample/"
        "sample_audvis_raw-trans.fif; infant_subjects/ANTS2-0Years3T; school_subjects/sub-Z213, sub-Z209, sub-Z226 "
        "(children A-C, OpenNeuro ds005234: mri/T1.mgz and the prepared bem/ and label/ files)",
        "$OPMSQUID_CACHE/anatomy/<subject>_fullres.npz (full-resolution white surfaces, opmsquid.anatomy.full_resolution)"]
    values = dict(
        order=list(KEYS), labels={k: LABEL[k] for k in KEYS}, subjects={k: M[k]["subject"] for k in KEYS},
        white_vertices={k: M[k]["n_white"] for k in KEYS},
        closest_white_vertex={k: dict(region=cw[k]["region"], region_words=region_words(cw[k]["region"]),
                                      distance_mm=cw[k]["distance_mm"], mri_frame_mm=cw[k]["mri_frame_mm"]) for k in KEYS},
        sections=drawn["sections"],
        cumulative_share_grid_mm=dict(start=0.0, stop=X_MAX_MM, step=0.5),
        cumulative_share=drawn["cumulative_share_at_0p5mm"],
        white_within_mm={k: {f"{w}/{thr:g}": dict(n_vertices=A[k][rk][f"below_{thr:g}mm"]["n_vertices"],
                                                   area_cm2=A[k][rk][f"below_{thr:g}mm"]["area_cm2"])
                             for w, rk in (("scalp_used", "white_to_scalp_used"), ("mri_boundary", "white_to_mri_boundary"))
                             for thr in near} for k in KEYS},
        offsets_cap_otsu=drawn["offsets_cap_otsu"],
        offsets_cap_half_max_median_mm=half,
        adult_reference_median_mm=drawn["adult_reference_mm"], tolerance_mm=tol, near_mm=near,
        verdicts=verdict, footnote=drawn["footnote"])
    checks = dict(
        rule=f"recomputed with the check's functions and compared value by value with {QC_JSON} (relative tolerance "
             f"{TOL:g}); the script stops at the first difference. Not recomputed: the verdicts (read), the closure "
             "(remeasured_cap_offset_mm) and the registration checks, which the figure does not show",
        values_equal={k: v for k, v in counts.items()},
        qc_results_commit=J["provenance"]["commit"],
        versions=dict(mne=mne.__version__, numpy=np.__version__, stored=dict(mne=J["provenance"]["mne_version"],
                                                                             numpy=J["provenance"]["numpy_version"])),
        runtime_s={k: round(M[k]["runtime_s"], 1) for k in KEYS},
        peak_memory_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1e9 if sys.platform == "darwin" else 1e6), 2))
    return dict(inputs=inputs, description=description, alt=alt, caption_draft=caption, values=values, checks=checks)


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    J = json.loads((ROOT / QC_JSON).read_text())
    cfg = tomllib.loads((ROOT / CONFIG).read_text())
    M = {key: measure(key, cfg) for key in KEYS}
    counts = check(M, J, cfg)
    print(f"stored records reproduced: {sum(sum(c.values()) for c in counts.values()):,} values", flush=True)
    drawn = draw(M, J)
    style.write_provenance(OUT_JSON, io.json_safe({NAME: entry(M, J, cfg, drawn, counts)}))
    print(f"results/report/{OUT_JSON} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
