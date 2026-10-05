#!/usr/bin/env python3
"""The adult's sensor arrays as a committed result: the Neuromag magnetometer coil centres, the
site-matched and the dense OPM arrays and the head surface that the report's arrays figure draws
(scripts/report_figures_clean.py, Figure_R15_arrays.png), so that the figure is drawn from results/
only. The G2 run kept no copy of its arrays (cache/g2/state.pkl holds an earlier run's arrays, with 95
and 205 sites), so they are rebuilt here from the anatomy by the G2 run's own code and checked against
every descriptor of them that the run stored. Nothing is simulated.

Arrays, as scripts/g2_adult_comparison.py builds them (g2.build_arrays on the MNE sample subject,
anatomy.load_sample, with the digitisation of the sample recording):
  Neuromag      opmsquid.neuromag.load_info('T3'): the 102 magnetometer coil centres of the sample
                recording (one per sensor site; each site also holds two planar gradiometers, 306
                channels), at the adult's measured head position (the recording's device-to-head
                transform: G2's primary position);
  site-matched  opmsquid.g2.matched_opm: the Neuromag sites that fit the OPM placement rules, projected
                onto the scalp (98 of 102);
  dense         opmsquid.g2.dense_opm(..., 'opm_dense'): the densest feasible single-axis array (208 sites).
G2's channel-budget control (204 sites, a subset of the dense array) is not exported: the main text does
not analyse it. The ray casting inside g2.matched_opm (opmsquid.neuromag.ray_mesh_distance) runs here 8
rays at a time instead of its default 64: each ray's intersections are computed independently of the
others, so the result is the same, with about an eighth of the peak memory (each temporary holds rays x
539,108 scalp triangles x 3 values).

Checks (the script stops unless every one holds):
  results/g2/g2_summary.json :: arrays.<array> (every stored descriptor: site counts, excluded sites,
      per-site outward shifts, spacings, standoff, scalp gap, the exact cell-to-head-surface clearance
      computed as the G2 run computes it; numbers within 1e-9), arrays.head_surface_conform (what
      conforming the BEM head surface to the MRI scalp did, within 1e-9 mm), and
      bridge_to_sphere.sensor_distance_mm (median, 5th and 95th percentile of the sensor-to-scalp
      distances of the Neuromag magnetometers and of both OPM arrays, as the G2 run computes them:
      nearest MRI scalp vertex; within 1e-9 mm);
  results/g3b/g3b_geometry_sections.json :: anatomies.adult.opm_dense_mm: the G3B run built the adult's
      dense array with the same code; its stored sites (rounded to 0.01 mm) equal these within the
      rounding.

Output: results/g2/g2_arrays.json (or --out). Head frame of the sample subject (Neuromag convention from
the digitised fiducials: x towards the right preauricular point, y towards the nasion, z up), mm rounded
to 0.01 mm; unit vectors rounded to 1e-6:
  fiducials_mm: LPA, nasion, RPA of the sample recording;
  neuromag: magnetometer names, coil centres and coil normals, and the device-to-head transform (4 x 4,
      metres, unrounded);
  opm_matched, opm_dense: sensing centres and sensitive axes, the Neuromag site (magnetometer index) of
      each site-matched site, and the stored descriptors;
  head_surface: the boundary-element head surface on the MRI scalp (anatomy.load_sample: the 3-layer
      BEM's outer surface, conformed to the MRI scalp by anatomy.head_on_scalp; the forward models
      subdivide it once, the same flat geometry), vertices and triangles (outward-oriented: checked),
      which the figure renders as the head.
Provenance: this export's commit (opmsquid.io) and the commits of the G2 summary and the G3B export.

Usage: OPMSQUID_DATA=<data dir> PYTHONPATH=src .venv/bin/python scripts/export_g2_arrays.py [--out <json>]
"""
from __future__ import annotations

import argparse
import contextlib
import functools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import mne  # noqa: E402
from mne.io.constants import FIFF  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from opmsquid import anatomy, g2, io, neuromag, opm, paths  # noqa: E402

SUMMARY = "results/g2/g2_summary.json"
GEOMETRY = "results/g3b/g3b_geometry_sections.json"  # scripts/export_g3b_geometry.py
OUT = ROOT / "results" / "g2" / "g2_arrays.json"
OPMS = ("opm_matched", "opm_dense")
RAY_CHUNK = 8  # rays per step of neuromag.ray_mesh_distance (its default: 64)
TOL = 1e-9  # stored descriptors (mm, or dimensionless)
DECIMALS = 2  # mm -> 0.01 mm
AXIS_DECIMALS = 6


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mm(x: np.ndarray, decimals: int = DECIMALS) -> list:
    """Values rounded to ``decimals`` (no negative zeros)."""
    return (np.round(np.asarray(x, float), decimals) + 0.0).tolist()


@contextlib.contextmanager
def ray_chunk(n: int):
    """opmsquid.neuromag.ray_mesh_distance with ``chunk=n`` (module docstring)."""
    orig = neuromag.ray_mesh_distance
    neuromag.ray_mesh_distance = functools.partial(orig, chunk=n)
    try:
        yield
    finally:
        neuromag.ray_mesh_distance = orig


def plain(x):
    """JSON round trip (string keys, lists), for comparing recomputed and stored descriptors."""
    return json.loads(json.dumps(io.json_safe(x)))


def same(have, want, path: str, bad: list) -> None:
    if isinstance(want, dict):
        if not isinstance(have, dict) or set(have) != set(want):
            bad.append(f"{path}: keys {sorted(have) if isinstance(have, dict) else have} vs {sorted(want)}")
            return
        for k in want:
            same(have[k], want[k], f"{path}.{k}", bad)
    elif isinstance(want, list):
        if not isinstance(have, list) or len(have) != len(want):
            bad.append(f"{path}: {have} vs {want}")
            return
        for i, (a, b) in enumerate(zip(have, want)):
            same(a, b, f"{path}[{i}]", bad)
    elif isinstance(want, bool) or isinstance(want, str) or want is None:
        if have != want:
            bad.append(f"{path}: {have!r} vs {want!r}")
    elif isinstance(want, (int, float)):
        if not isinstance(have, (int, float)) or isinstance(have, bool) or abs(have - want) > TOL:
            bad.append(f"{path}: {have!r} vs {want!r}")
    else:
        bad.append(f"{path}: cannot compare {type(want).__name__}")


def descriptors(a: g2.Array, head: dict, head_mri: np.ndarray) -> dict:
    """The descriptors the G2 run stores for an array (scripts/g2_adult_comparison.py main: the array's metadata
    without its per-site lists except the excluded sites, its channel count and, for an OPM array, the exact distance
    of its cells from the BEM head surface)."""
    out = {k: v for k, v in a.meta.items() if not isinstance(v, (list, np.ndarray)) or k == "excluded_sites"}
    out["channels"] = a.n
    if a.name != "squid":
        out["exact_cell_clearance_min_mm"] = float(opm.exact_cell_clearance(a.info, head, head_mri).min() * 1e3)
    return out


def sensor_distances(scalp_head: np.ndarray, pos: dict) -> dict:
    """Sensor-to-scalp distances (nearest MRI scalp vertex, head frame) as the G2 run summarises them
    (scripts/g2_adult_comparison.py geometry_distances and bridge_to_sphere)."""
    tree = cKDTree(scalp_head)
    out = {}
    for name, p in pos.items():
        v = tree.query(p)[0]
        out[name] = dict(median=float(np.median(v) * 1e3), p5=float(np.percentile(v, 5) * 1e3), p95=float(np.percentile(v, 95) * 1e3))
    return out


def outward_volume(rr: np.ndarray, tris: np.ndarray) -> float:
    """Signed volume enclosed by a closed triangle mesh: positive when its triangles are outward-oriented."""
    a, b, c = rr[tris[:, 0]], rr[tris[:, 1]], rr[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# ----------------------------------------------------------------------------------------------
def _flat(x: list) -> bool:
    """A list of scalars, or of lists of scalars (a point list, a matrix): written on one line."""
    return all(not isinstance(v, (dict, list)) or (isinstance(v, list) and not any(isinstance(w, (dict, list)) for w in v))
               for v in x)


def dumps(obj, level: int = 0) -> str:
    """JSON with dicts and lists of containers indented and each point list or matrix on one line (the layout of
    scripts/export_g3b_geometry.py)."""
    pad = " " * (level + 1)
    if isinstance(obj, dict) and obj:
        return "{\n" + ",\n".join(f"{pad}{json.dumps(str(k))}: {dumps(v, level + 1)}" for k, v in obj.items()) + "\n" + " " * level + "}"
    if isinstance(obj, list) and obj and not _flat(obj):
        return "[\n" + ",\n".join(pad + dumps(v, level + 1) for v in obj) + "\n" + " " * level + "]"
    return json.dumps(obj, separators=(",", ":"), allow_nan=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=OUT, help=f"output JSON (default: {OUT.relative_to(ROOT)})")
    args = ap.parse_args()
    mne.set_log_level("WARNING")
    s = json.loads((ROOT / SUMMARY).read_text())
    geo = json.loads((ROOT / GEOMETRY).read_text())
    t0 = time.time()

    # the G2 run's anatomy and arrays (scripts/g2_adult_comparison.py: Study and g2.build_arrays)
    sub = anatomy.load_sample()
    dig = mne.io.read_info(paths.require(paths.SAMPLE_MEG / neuromag.RAW_FILE, "MNE sample recording"), verbose=False)
    squid_info = neuromag.load_info("T3")
    arrays = {"squid": g2.Array("squid", squid_info, neuromag.channel_kinds(squid_info), None,
                                dict(sites=102, channels=306, axes="1 mag + 2 planar grad per site"))}  # as g2.build_arrays
    with ray_chunk(RAY_CHUNK):
        arrays["opm_matched"] = g2.matched_opm(sub, dig)
    log(f"site-matched OPM array rebuilt: {arrays['opm_matched'].n} sites ({time.time() - t0:.0f} s)")
    arrays["opm_dense"] = g2.dense_opm(sub, dig, "opm_dense")
    log(f"dense OPM array rebuilt: {arrays['opm_dense'].n} sites ({time.time() - t0:.0f} s)")

    head_mri = sub.trans["trans"]
    mri_head = np.linalg.inv(head_mri)
    head = next(x for x in sub.bem_surfaces if x["id"] == FIFF.FIFFV_BEM_SURF_ID_HEAD)
    dev_head = np.asarray(squid_info["dev_head_t"]["trans"], float)
    kinds = arrays["squid"].kinds
    mag_names = [c["ch_name"] for c, kd in zip(squid_info["chs"], kinds) if kd == "mag"]
    mag_loc = np.array([c["loc"] for c, kd in zip(squid_info["chs"], kinds) if kd == "mag"])
    mags = mag_loc[:, :3] @ dev_head[:3, :3].T + dev_head[:3, 3]  # coil centres, head frame [m]
    mag_nn = mag_loc[:, 9:12] @ dev_head[:3, :3].T  # coil normals (the coil's z axis), head frame
    pos = {n: np.array([c["loc"][:3] for c in arrays[n].info["chs"]]) for n in OPMS}  # head frame (device = head)
    ax = {n: np.array([c["loc"][9:12] for c in arrays[n].info["chs"]]) for n in OPMS}

    # checks against the G2 run's stored descriptors
    bad = []
    for name, a in arrays.items():
        same(plain(descriptors(a, head, head_mri)), s["arrays"][name], f"arrays.{name}", bad)
    conform = plain(sub.head_conform)  # head_on_scalp's report; 'reassigned' and 'repaired' were added after the G2 run
    stored_conform = s["arrays"]["head_surface_conform"]
    same({k: conform.get(k) for k in stored_conform}, stored_conform, "arrays.head_surface_conform", bad)
    scalp_head = sub.scalp.rr @ mri_head[:3, :3].T + mri_head[:3, 3]
    dist = sensor_distances(scalp_head, {"squid": mags, **pos})
    same(dist, {k: s["bridge_to_sphere"]["sensor_distance_mm"][k] for k in dist}, "bridge_to_sphere.sensor_distance_mm", bad)
    g3b_dense = np.asarray(geo["anatomies"]["adult"]["opm_dense_mm"], float)
    g3b_err = float(np.abs(g3b_dense - pos["opm_dense"] * 1e3).max()) if g3b_dense.shape == pos["opm_dense"].shape else np.inf
    if g3b_err > 0.5 * 10.0 ** -DECIMALS + 1e-6:  # the export rounds to 0.01 mm
        bad.append(f"{GEOMETRY} adult opm_dense_mm: differs by {g3b_err} mm")
    head_rr = head["rr"] @ mri_head[:3, :3].T + mri_head[:3, 3]
    vol = outward_volume(head_rr, head["tris"])
    if vol <= 0:
        bad.append("the BEM head surface's triangles are not outward-oriented")
    if bad:
        raise SystemExit("the rebuilt arrays differ from the stored G2 run:\n  " + "\n  ".join(bad))
    log(f"checks passed: descriptors, sensor distances, head surface; G3B adult dense sites within {g3b_err:.4f} mm")

    # the Neuromag site of each site-matched OPM site: opm.matched_to_neuromag keeps the magnetometers in their order
    # (site_id = magnetometer index) less the excluded ones; checked geometrically: each sensing centre lies nearest to
    # the line through its own magnetometer's coil centre along the coil normal (the projection's ray)
    excluded = set(s["arrays"]["opm_matched"]["excluded_sites"])
    neuromag_site = np.array([i for i in range(len(mags)) if i not in excluded])
    rel = pos["opm_matched"][:, None] - mags[None]  # (OPM site, magnetometer, 3)
    off_ray = np.linalg.norm(rel - np.einsum("ijk,jk->ij", rel, mag_nn)[..., None] * mag_nn[None], axis=2)
    nearest = np.argmin(off_ray, axis=1)
    if len(neuromag_site) != len(pos["opm_matched"]) or not np.array_equal(nearest, neuromag_site):
        raise SystemExit("the site-matched OPM sites do not follow the magnetometers they are projected from")
    own = off_ray[np.arange(len(nearest)), nearest]
    log(f"site-matched sites: each nearest to its own magnetometer's normal line (at most {own.max() * 1e3:.2f} mm off it)")
    fids = opm.fiducials_head(dig)
    result = dict(
        status=("NEW export (revision): the adult's arrays drawn in Figure_R15_arrays.png, rebuilt with the G2 run's code and "
                "checked against its stored descriptors"),
        description=(
            "Head frame of the MNE sample subject (Neuromag convention from the digitised fiducials: x towards the right "
            "preauricular point, y towards the nasion, z up), mm rounded to 0.01 mm, unit vectors to 1e-6. neuromag: the 102 "
            "magnetometer coil centres and coil normals of the sample recording (opmsquid.neuromag.load_info('T3'); each of "
            "the 102 sensor sites also holds two planar gradiometers: 306 channels) at the adult's measured head position "
            "(dev_head_t, the recording's device-to-head transform: G2's primary position). opm_matched: the site-matched "
            "OPM array (opmsquid.g2.matched_opm: the Neuromag sites that fit the OPM placement rules, each magnetometer "
            "projected along its inward normal onto the MRI scalp, the sensing centre placed at the standoff along the "
            "smoothed head-surface normal and moved out where needed; neuromag_site = the magnetometer index, order of "
            "neuromag.names). opm_dense: the dense OPM array (opmsquid.g2.dense_opm 'opm_dense'). Sensing centres and "
            "sensitive axes (the BEM head-surface normal averaged within 15 mm). descriptors: as stored in "
            f"{SUMMARY} arrays (equal, checked). head_surface: the boundary-element head surface on the MRI scalp "
            "(opmsquid.anatomy.load_sample: the outer surface of sample-5120-5120-5120-bem.fif conformed to the MRI scalp "
            "sample-head.fif by anatomy.head_on_scalp, conform = what that did, equal to the stored "
            "arrays.head_surface_conform in every stored field ('reassigned' and 'repaired' were added to that report after "
            "the G2 run); the forward models subdivide it once, the same flat geometry), triangles "
            "outward-oriented. The arrays are rebuilt here (the G2 run kept no copy of them); every descriptor of them that "
            f"{SUMMARY} stores, and the median, 5th and 95th percentile of each array's sensor-to-scalp distance "
            "(bridge_to_sphere.sensor_distance_mm), are reproduced within 1e-9, and the dense sites equal the G3B run's "
            f"adult dense array ({GEOMETRY}) within the 0.01-mm rounding (checks)."),
        inputs=[f"{SUMMARY} :: arrays.squid, arrays.opm_matched, arrays.opm_dense, arrays.head_surface_conform, "
                "bridge_to_sphere.sensor_distance_mm (squid, opm_matched, opm_dense), provenance.commit (checks)",
                f"{GEOMETRY} :: anatomies.adult.opm_dense_mm, provenance (check)",
                "MNE-sample-data/MEG/sample/sample_audvis_raw.fif (opmsquid.neuromag.load_info('T3'): coil centres, normals, "
                "device-to-head transform; digitised fiducials)",
                "MNE-sample-data/subjects/sample/bem: sample-head.fif (MRI scalp), sample-5120-5120-5120-bem.fif (BEM surfaces); "
                "MEG/sample/sample_audvis_raw-trans.fif (head to MRI): opmsquid.anatomy.load_sample (OPMSQUID_DATA)",
                "code: opmsquid.g2.matched_opm and dense_opm (as g2.build_arrays), opmsquid.opm.exact_cell_clearance"],
        units=dict(coordinates="mm, head frame, rounded to 0.01 mm", axes="unit vectors, head frame, rounded to 1e-6",
                   transforms="4 x 4 homogeneous, metres, unrounded"),
        fiducials_mm={k: mm(v * 1e3) for k, v in fids.items()},
        neuromag=dict(source="MNE-sample-data/MEG/sample/sample_audvis_raw.fif, opmsquid.neuromag.load_info('T3')",
                      sites=int(s["arrays"]["squid"]["sites"]), channels=int(s["arrays"]["squid"]["channels"]),
                      head_position="the adult's measured head position (device-to-head transform of the sample recording)",
                      dev_head_t=dev_head, names=mag_names, coil_centres_mm=mm(mags * 1e3),
                      coil_normals=mm(mag_nn, AXIS_DECIMALS)),
        opm_matched=dict(descriptors=s["arrays"]["opm_matched"], neuromag_site=neuromag_site,
                         sensing_centres_mm=mm(pos["opm_matched"] * 1e3), sensitive_axes=mm(ax["opm_matched"], AXIS_DECIMALS)),
        opm_dense=dict(descriptors=s["arrays"]["opm_dense"], sensing_centres_mm=mm(pos["opm_dense"] * 1e3),
                       sensitive_axes=mm(ax["opm_dense"], AXIS_DECIMALS)),
        head_surface=dict(source="MNE-sample-data/subjects/sample/bem/sample-5120-5120-5120-bem.fif (head surface), conformed to "
                                 "bem/sample-head.fif by opmsquid.anatomy.load_sample",
                          n_vertices=len(head_rr), n_triangles=len(head["tris"]), conform=conform,
                          enclosed_volume_cm3=vol * 1e6, vertices_mm=mm(head_rr * 1e3), triangles=np.asarray(head["tris"], int)),
        checks=dict(stored_descriptors=f"{SUMMARY} arrays.squid/opm_matched/opm_dense/head_surface_conform: equal within {TOL:g}",
                    sensor_distance_mm=dist, g3b_adult_dense_max_difference_mm=g3b_err,
                    matched_site_max_distance_from_its_magnetometer_ray_mm=float(own.max() * 1e3), ray_chunk=RAY_CHUNK),
        provenance=dict(commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__,
                        summary=dict(file=SUMMARY, commit=s["provenance"]["commit"]),
                        g3b_geometry=dict(file=GEOMETRY, commit=geo["provenance"]["commit"],
                                          source_state_commit=geo["provenance"]["source_state"]["commit"])))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(dumps(io.json_safe(result)) + "\n")
    log(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.2f} MB; {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
