#!/usr/bin/env python3
"""The G3B helmet geometry as a committed result: the scalp and white-surface sections, the Neuromag
magnetometer coil centres and the dense OPM sites that the report's geometry figure draws
(scripts/report_figures_clean.py, Figure_R12_geometry.png), exported from the stored state of the G3B
run so that the figure is drawn from results/ only. Nothing is simulated or fitted again: the
transforms, scale factors and OPM sites are copied from the state, and the sections and coil centres
are computed from them and the MRI surfaces exactly as the figure computed them before.

Inputs (all read only):
  cache/g3b/state.pkl (OPMSQUID_CACHE, or --state): per anatomy the device-to-head transforms of the
      placements 'top' (fixed adult helmet, head raised to top contact) and 'counterfactual_x-centred'
      (the helmet scaled by k about the head origin of the laterally centred pose, before any top
      contact), k, and the dense OPM sites (head frame). The script stops unless the state's commit is
      that of results/g3b/g3b_summary.json and its moves, scale factors, gaps and site counts are the
      stored ones.
  results/g3b/g3b_summary.json: those checks, the anatomies' descriptions and the two scale factors
      of the size-only controls (config.anatomy.school_age_scale; the 24-month template's head
      circumference over the adult's, as scripts/g3b_pediatric_helmet.py computes it).
  OPMSQUID_DATA: the magnetometer coil centres of the sample recording (opmsquid.neuromag.load_info('T3'),
      device frame) and each anatomy's MRI scalp and white surfaces, loaded as the G3B run loads them
      (opmsquid.anatomy: load_sample, load_template, load_school; the scaled controls by anatomy.scaled).

Output: results/g3b/g3b_geometry_sections.json (or --out). For each of the nine G3B anatomies, in its
head frame (Neuromag convention of its own fiducials: x right, y anterior, z up), coordinates in mm
rounded to 0.01 mm:
  sections.sagittal.scalp, sections.coronal.scalp: polylines where the MRI scalp mesh crosses the
      planes x = 0 (points: y, z) and y = 0 (points: x, z);
  sections.coronal.white.lh / .rh: the same for the full white surface of each hemisphere (the surface
      held by the oct-6 source space), which the figure draws in its coronal sections;
  coils_mm.fixed_top, coils_mm.scaled_x_centred: the 102 magnetometer coil centres (x, y, z) of the
      fixed helmet at top contact and of the helmet scaled with the head, laterally centred;
  opm_dense_mm: the dense OPM sites (sensing centres, x, y, z);
  transforms: the device-to-head transforms of both placements and the head-to-MRI transform (4 x 4,
      metres, as stored), k and the scaling centre;
  gaps_mm: the coil-to-scalp gaps (nearest MRI scalp vertex) recomputed here from the unrounded
      geometry; they must equal the stored ones within 1e-6 mm.
A polyline's consecutive points are the segments of the section (one per crossed triangle, between the
two points where its edges cross the plane; a vertex on the plane counts as just above it), joined
where neighbouring triangles share the crossed edge; a closed loop repeats its first point at the end.
The polylines hold exactly the segments that report_figures_clean.py drew before (checked here).
Provenance: this export's commit (opmsquid.io) and the state's and the summary's commits.

Usage: OPMSQUID_DATA=<data dir> OPMSQUID_CACHE=<cache dir> PYTHONPATH=src .venv/bin/python \\
           scripts/export_g3b_geometry.py [--state <state.pkl>] [--out <json>]
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import mne  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from opmsquid import anatomy, io, neuromag, paths, pediatric as P  # noqa: E402

SUMMARY = "results/g3b/g3b_summary.json"
OUT = ROOT / "results" / "g3b" / "g3b_geometry_sections.json"
STATE_LABEL = "cache/g3b/state.pkl"
# as scripts/g3b_pediatric_helmet.py (ANATOMIES, TEMPLATES, SCHOOL)
ANATOMIES = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
TEMPLATES = {"infant2yr": "ANTS2-0Years3T", "infant18mo": "ANTS18-0Months3T", "infant12mo": "ANTS12-0Months3T"}
SCHOOL = {"childA": "sub-Z213", "childB": "sub-Z209", "childC": "sub-Z226"}
PLACEMENTS = ("top", "counterfactual_x-centred")
PLANES = (("sagittal", 0), ("coronal", 1))  # plane name, axis set to 0
DECIMALS = 2  # mm -> 0.01 mm
GAP_TOL_MM = 1e-6
STATE_TOL = 1e-9


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mm(x: np.ndarray) -> list:
    """Coordinates [mm] rounded to 0.01 mm (no negative zeros)."""
    return (np.round(np.asarray(x, float), DECIMALS) + 0.0).tolist()


# ----------------------------------------------------------------------------------------------
# sections
def trails(seg: np.ndarray, n_nodes: int) -> list[list[int]]:
    """Node sequences that cover every segment (node pair) exactly once: maximal trails, started first
    at nodes of odd degree (open chains); a closed trail repeats its first node at the end."""
    adj = [[] for _ in range(n_nodes)]
    for i, (a, b) in enumerate(seg):
        adj[a].append(i)
        adj[b].append(i)
    used = np.zeros(len(seg), bool)
    deg = np.array([len(x) for x in adj])
    out = []
    for start in list(np.flatnonzero(deg % 2 == 1)) + list(range(n_nodes)):
        while any(not used[i] for i in adj[start]):
            path, cur = [int(start)], int(start)
            while True:
                i = next((i for i in adj[cur] if not used[i]), None)
                if i is None:
                    break
                used[i] = True
                a, b = seg[i]
                cur = int(b if a == cur else a)
                path.append(cur)
            out.append(path)
    if not used.all():
        raise SystemExit("section: a segment was left out of the polylines")
    return out


def plane_section_reference(rr: np.ndarray, tris: np.ndarray, axis: int) -> np.ndarray:
    """Segments (n, 2, 3) where the mesh crosses the plane rr[:, axis] = 0, computed as the geometry figure
    computed them before this export (report_figures_clean.plane_section): the reference of ``section``."""
    s = rr[:, axis].copy()
    s[s == 0.0] = 1e-12
    st = s[tris]
    t = tris[(st.min(axis=1) < 0) & (st.max(axis=1) > 0)]
    pts = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a, b = t[:, i], t[:, j]
        cross = np.sign(s[a]) != np.sign(s[b])
        with np.errstate(divide="ignore", invalid="ignore"):  # edges that do not cross are discarded below
            f = s[a] / (s[a] - s[b])
            p = rr[a] + f[:, None] * (rr[b] - rr[a])
        pts.append(np.where(cross[:, None], p, np.nan))
    pts = np.stack(pts, axis=1)
    two = np.argsort(np.isnan(pts[:, :, 0]), axis=1, kind="stable")[:, :2]
    return pts[np.arange(len(pts))[:, None], two]


def section(rr: np.ndarray, tris: np.ndarray, axis: int, what: str) -> tuple[list[np.ndarray], int, float]:
    """Polylines ((n, 3) arrays) where the triangle mesh (rr, tris) crosses the plane rr[:, axis] = 0, the
    number of segments, and the largest difference [mm] of the segments from ``plane_section_reference``
    (same crossing rule, a vertex on the plane counting as just above it; same points, each crossed edge
    computed once). Stops unless the polylines hold every reference segment exactly once."""
    s = rr[:, axis].copy()
    s[s == 0.0] = 1e-12
    st = s[tris]
    t = tris[(st.min(axis=1) < 0) & (st.max(axis=1) > 0)]
    if not len(t):
        return [], 0, 0.0
    crossed = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a, b = t[:, i], t[:, j]
        crossed.append(np.where((np.sign(s[a]) != np.sign(s[b]))[:, None], np.sort(np.column_stack([a, b]), axis=1), -1))
    crossed = np.stack(crossed, axis=1)  # (n, 3 edges, 2 vertices); -1 where the edge does not cross
    keep = crossed[:, :, 0] >= 0
    if not np.all(keep.sum(axis=1) == 2):
        raise SystemExit(f"{what}: a crossed triangle without exactly two crossing edges")
    edges, node = np.unique(crossed[keep], axis=0, return_inverse=True)
    node = np.asarray(node).reshape(-1, 2)  # per crossed triangle, its two crossed edges (edge order) as point indices
    a, b = edges[:, 0], edges[:, 1]
    f = s[a] / (s[a] - s[b])
    pts = rr[a] + f[:, None] * (rr[b] - rr[a])
    lines = trails(node, len(edges))

    def rows(x):  # unordered index pairs, sorted
        x = np.sort(np.asarray(x), axis=1)
        return x[np.lexsort(x.T[::-1])]

    ref = plane_section_reference(rr, tris, axis)
    walked = np.concatenate([np.column_stack([p[:-1], p[1:]]) for p in lines])
    err = float(np.abs(pts[node] - ref).max())  # same triangle order and edge order
    if len(ref) != len(node) or err > 1e-9 or np.abs(pts[:, axis]).max() > 1e-9 or not np.array_equal(rows(walked), rows(node)):
        raise SystemExit(f"{what}: the polylines do not hold the reference segments (max point difference {err:.2e} mm)")
    return [pts[p] for p in lines], len(node), err


# ----------------------------------------------------------------------------------------------
def load_state(path: Path, g: dict) -> dict:
    with open(paths.require(path, f"stored G3B state ({STATE_LABEL}; set OPMSQUID_CACHE or --state)"), "rb") as fh:
        st = pickle.load(fh)
    if st["provenance"]["commit"] != g["provenance"]["commit"]:
        raise SystemExit(f"{path}: run commit {st['provenance']['commit']} is not that of {SUMMARY} ({g['provenance']['commit']})")
    for k in ANATOMIES:
        run, stored = st["runs"][k], g["placements"][k]
        for name in PLACEMENTS:
            keys = ("min_dist_mm", "median_dist_mm", "max_dist_mm") + (("moved_mm",) if name == "top" else ("k", "k_nominal"))
            for key in keys:
                if abs(run["placements"][name][key] - stored[name][key]) > STATE_TOL:
                    raise SystemExit(f"{path}: {k} {name} {key} differs from {SUMMARY}")
        n = len(run["opm_pos"]["opm_dense"])
        if n != g["arrays"][k]["opm_dense"]["n_sites"]:
            raise SystemExit(f"{path}: {k} has {n} dense OPM sites, {SUMMARY} {g['arrays'][k]['opm_dense']['n_sites']}")
    return st


def scale_factors(g: dict) -> dict:
    """The size-only controls' factors as scripts/g3b_pediatric_helmet.py load_anatomies sets them (after
    checking that the run's anatomies are the ones listed here)."""
    a, cfg = g["anatomies"], g["config"]["anatomy"]
    if (list(a) != list(ANATOMIES) or list(cfg["templates"]) != list(TEMPLATES.values())
            or {c["key"]: c["subject"] for c in cfg["school"]} != SCHOOL):
        raise SystemExit(f"{SUMMARY}: the run's anatomies are not {ANATOMIES} with templates {TEMPLATES} and children {SCHOOL}")
    f = {"school": float(cfg["school_age_scale"]),
         "size2yr": a["infant2yr"]["head_size"]["ofc_mm"] / a["adult"]["head_size"]["ofc_mm"]}
    for k, v in f.items():
        if not a[k]["scale_note"].endswith(f": {v:.4f}"):
            raise SystemExit(f"{SUMMARY}: {k} scale note '{a[k]['scale_note']}' does not give the factor {v:.4f}")
    return f


def load_subjects(factors: dict) -> dict:
    adult = anatomy.load_sample()
    out = {"adult": (adult, "MNE-sample-data/subjects/sample: bem/sample-head.fif (MRI scalp), bem/sample-oct-6-src.fif "
                            "(white surfaces); MEG/sample/sample_audvis_raw-trans.fif (head to MRI)")}
    for k, f in factors.items():
        out[k] = (anatomy.scaled(adult, f, f"sample_x{f:.4f}"),
                  f"the adult's surfaces and head-to-MRI translation scaled by {f!r} about the MRI origin (opmsquid.anatomy.scaled)")
    for k, name in TEMPLATES.items():
        out[k] = (anatomy.load_template(name), f"infant_subjects/{name}/bem: {name}-head.fif (MRI scalp), {name}-oct-6-src.fif "
                                                f"(white surfaces), {name}-fiducials.fif (head frame)")
    for k, name in SCHOOL.items():
        out[k] = (anatomy.load_school(name), f"school_subjects/{name}/bem: {name}-head.fif (MRI scalp), {name}-oct-6-src.fif "
                                              f"(white surfaces), {name}-fiducials.fif (head frame)")
    return {k: out[k] for k in ANATOMIES}


def export_anatomy(k: str, sub, source: str, run: dict, stored: dict, mags: np.ndarray, g: dict) -> dict:
    """One anatomy in its head frame, computed as report_figures_clean.py computed it before this export."""
    mh = np.linalg.inv(sub.trans["trans"])  # MRI -> head
    scalp = P.scalp_head_frame(sub)  # m
    top = np.asarray(run["placements"]["top"]["trans"], float)
    cfx = np.asarray(run["placements"]["counterfactual_x-centred"]["trans"], float)
    kk = float(run["placements"]["counterfactual_x-centred"]["k"])
    fixed = mags @ top[:3, :3].T + top[:3, 3]  # fixed helmet, head at top contact (device -> head)
    centre = np.linalg.inv(cfx)[:3, 3]  # head origin, device frame (as opmsquid.pediatric.counterfactual_helmet)
    scaled = (centre + kk * (mags - centre)) @ cfx[:3, :3].T + cfx[:3, 3]
    tree = cKDTree(scalp)
    gaps = {}
    for name, c in (("top", fixed), ("counterfactual_x-centred", scaled)):
        dist = tree.query(c)[0] * 1e3
        gaps[name] = dict(min_dist_mm=float(dist.min()), median_dist_mm=float(np.median(dist)), max_dist_mm=float(dist.max()))
        for key, val in gaps[name].items():
            if abs(val - stored[name][key]) > GAP_TOL_MM:
                raise SystemExit(f"{k}: recomputed {name} {key} {val:.6f} differs from the stored {stored[name][key]:.6f}")
    scalp_mm, tris = scalp * 1e3, sub.scalp.tris
    white = {h: ((s["rr"] @ mh[:3, :3].T + mh[:3, 3]) * 1e3, s["tris"]) for h, s in zip(("lh", "rh"), sub.src)}
    sections, n_seg, err = {}, {}, 0.0
    for plane, axis in PLANES:
        cols = [1 - axis, 2]  # in-plane coordinates: (y, z) in the sagittal, (x, z) in the coronal section
        lines, n, e = section(scalp_mm, tris, axis, f"{k} {plane} scalp")
        sections[plane] = {"scalp": [mm(p[:, cols]) for p in lines]}
        n_seg[f"{plane}_scalp"], err = n, max(err, e)
        if plane == "coronal":  # the midsagittal plane passes between the hemispheres: white sections in the coronal plane only
            sections[plane]["white"] = {}
            for h, (rr, tr) in white.items():
                lines, n, e = section(rr, tr, axis, f"{k} {plane} {h} white")
                sections[plane]["white"][h] = [mm(p[:, cols]) for p in lines]
                n_seg[f"{plane}_white_{h}"], err = n, max(err, e)
    opm = np.asarray(run["opm_pos"]["opm_dense"], float)
    cf = run["placements"]["counterfactual_x-centred"]
    return dict(
        description=g["anatomies"][k]["description"], scale_note=g["anatomies"][k]["scale_note"], surfaces=source,
        transforms={"head_to_mri": np.asarray(sub.trans["trans"], float),
                    "top": dict(rule=run["placements"]["top"]["rule"], device_to_head=top,
                                moved_mm=float(run["placements"]["top"]["moved_mm"])),
                    "counterfactual_x-centred": dict(rule=cf["rule"], device_to_head=cfx, k=kk, k_nominal=float(cf["k_nominal"]),
                                                     scaling_centre_device_m=centre)},
        coils_mm=dict(fixed_top=mm(fixed * 1e3), scaled_x_centred=mm(scaled * 1e3)),
        opm_dense_mm=mm(opm * 1e3),
        sections=sections,
        gaps_mm=gaps,
        counts=dict(scalp_vertices=len(scalp_mm), scalp_triangles=len(tris),
                    white_vertices={h: len(rr) for h, (rr, _) in white.items()},
                    white_triangles={h: len(tr) for h, (_, tr) in white.items()},
                    section_segments=n_seg, polylines={f"{p}_{w}": (len(v) if w == "scalp" else {h: len(x) for h, x in v.items()})
                                                       for p, d in sections.items() for w, v in d.items()},
                    coils=len(fixed), opm_dense_sites=len(opm)),
        section_check_max_point_difference_mm=err)


# ----------------------------------------------------------------------------------------------
def _flat(x: list) -> bool:
    """A list of scalars, or of lists of scalars (a point list, a matrix): written on one line."""
    return all(not isinstance(v, (dict, list)) or (isinstance(v, list) and not any(isinstance(w, (dict, list)) for w in v))
               for v in x)


def dumps(obj, level: int = 0) -> str:
    """JSON with dicts and lists of containers indented and each point list or matrix on one line."""
    pad = " " * (level + 1)
    if isinstance(obj, dict) and obj:
        return "{\n" + ",\n".join(f"{pad}{json.dumps(str(k))}: {dumps(v, level + 1)}" for k, v in obj.items()) + "\n" + " " * level + "}"
    if isinstance(obj, list) and obj and not _flat(obj):
        return "[\n" + ",\n".join(pad + dumps(v, level + 1) for v in obj) + "\n" + " " * level + "]"
    return json.dumps(obj, separators=(",", ":"), allow_nan=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--state", type=Path, default=paths.CACHE / "g3b" / "state.pkl",
                    help="stored G3B state (default: $OPMSQUID_CACHE/g3b/state.pkl; only read)")
    ap.add_argument("--out", type=Path, default=OUT, help=f"output JSON (default: {OUT.relative_to(ROOT)})")
    args = ap.parse_args()
    mne.set_log_level("WARNING")
    g = json.loads((ROOT / SUMMARY).read_text())
    log(f"reading {args.state}")
    st = load_state(args.state, g)
    info = neuromag.load_info("T3")
    mags = P.magnetometer_positions(info)  # device frame [m]
    names = [c["ch_name"] for c, kd in zip(info["chs"], neuromag.channel_kinds(info)) if kd == "mag"]
    factors = scale_factors(g)
    subjects = load_subjects(factors)
    heads = {}
    for k, (sub, source) in subjects.items():
        heads[k] = export_anatomy(k, sub, source, st["runs"][k], g["placements"][k], mags, g)
        c = heads[k]["counts"]
        log(f"{k}: gaps reproduced; sections {c['section_segments']} segments; {c['opm_dense_sites']} OPM sites")
    result = dict(
        status="NEW export (revision): the G3B helmet geometry drawn in Figure_R12_geometry.png, from the stored state of the G3B run",
        description=(
            "Per anatomy, in its head frame (Neuromag convention of its own fiducials: x right, y anterior, z up), mm rounded to "
            "0.01 mm: sections.<plane>.scalp, polylines where the MRI scalp mesh crosses the plane through the head origin "
            "(sagittal x = 0, points (y, z); coronal y = 0, points (x, z)); sections.coronal.white.<hemi>, the same for the full "
            "white surface of each hemisphere (the surface of the oct-6 source space; the midsagittal plane passes between the "
            "hemispheres); a polyline's consecutive points are the section's segments (one per crossed triangle, a vertex on the "
            "plane counting as just above it), joined where neighbouring triangles share the crossed edge, and a closed loop "
            "repeats its first point at the end. coils_mm.fixed_top: the 102 Neuromag magnetometer coil centres of the fixed "
            "adult helmet with the head at top contact (device_to_head of transforms.top applied to neuromag.coil_centres_device_mm); "
            "coils_mm.scaled_x_centred: the counterfactual helmet scaled with the head and laterally centred (coil centres scaled "
            "by k about scaling_centre_device_m, the head origin of the laterally centred pose in the device frame, then "
            "device_to_head of transforms.counterfactual_x-centred). opm_dense_mm: the dense OPM sites (sensing centres) refitted "
            "to the head. Transforms: 4 x 4, metres, unrounded, as stored in the G3B state (device_to_head) or as the anatomy "
            "is loaded (head_to_mri). gaps_mm: magnetometer coil centre to the nearest MRI scalp vertex, recomputed here from the "
            "unrounded geometry, equal to results/g3b/g3b_summary.json placements within 1e-6 mm (checked). The polylines hold "
            "exactly the segments the geometry figure drew from the same surfaces before this export (checked)."),
        inputs=[f"{STATE_LABEL} (OPMSQUID_CACHE; read only) :: runs[<anatomy>].placements['top'|'counterfactual_x-centred'] "
                "(trans, k, k_nominal, moved_mm, rule, min/median/max_dist_mm), runs[<anatomy>].opm_pos['opm_dense'], provenance",
                f"{SUMMARY} :: provenance.commit, placements[<anatomy>]['top'|'counterfactual_x-centred'], "
                "arrays[<anatomy>].opm_dense.n_sites, anatomies[<anatomy>] (description, scale_note, head_size.ofc_mm of the adult "
                "and the 24-month template), config.anatomy.school_age_scale",
                "MNE-sample-data/MEG/sample/sample_audvis_raw.fif (opmsquid.neuromag.load_info('T3'): magnetometer coil centres, "
                "device frame)",
                "MRI scalps and white surfaces via opmsquid.anatomy (OPMSQUID_DATA): see anatomies[<anatomy>].surfaces"],
        units=dict(coordinates="mm, head frame of each anatomy, rounded to 0.01 mm",
                   transforms="4 x 4 homogeneous, metres, unrounded", gaps="mm, unrounded"),
        planes={"sagittal": "x = 0; polyline points (y, z)", "coronal": "y = 0; polyline points (x, z)"},
        neuromag=dict(source="MNE-sample-data/MEG/sample/sample_audvis_raw.fif, opmsquid.neuromag.load_info('T3')",
                      magnetometers=len(mags), names=names, coil_centres_device_mm=mm(mags * 1e3)),
        scale_factors=factors,
        anatomies=heads,
        provenance=dict(commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__,
                        source_state=dict(file=STATE_LABEL, **st["provenance"]),
                        summary=dict(file=SUMMARY, commit=g["provenance"]["commit"], replotted_at_commit=g.get("replotted_at_commit"))))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(io.json_safe(result)) + "\n")
    log(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
