#!/usr/bin/env python3
"""Full-precision depth and explicit bin membership for the downloadable per-target tables.

The per-target tables print each target's depth below the scalp rounded (results/g2/g2_targets.csv and
results/g3b/g3b_targets_<anatomy>.csv to 0.01 mm, results/g3b_constant_gap/targets_<anatomy>.csv to
1e-6 mm), while every summary bins the unrounded depth into half-open bins, lo <= depth < hi. A target
9.996 mm below the scalp is counted below 10 mm but printed as 10.00, so a count taken from a rounded
table can differ from the one stated (child B: 247 cortical targets below 10 mm in g3b_summary.json,
245 below 10.00 mm in g3b_targets_childB.csv). Next to each table this script writes a compact
companion, <table>_depth.csv: the table's keys (hemi, vertno) in its row order, the depth as computed
(a double, in the shortest decimal form that reads back as the same double) and the bins that the
summaries of that table use (empty: outside every bin):

  results/g2/g2_targets_depth.csv
      depth_bin_5mm        the 5-mm bins of g2_summary.json (amplitude, detectability and log2 ratio vs
                           depth, 10-80 mm), also those of g2_covariance_validation/covariance_validation.json
      spike_depth_band     the depth bands of the simulated spikes (g4_epilepsy_adult.DEPTH_BANDS: 10-20,
                           20-30, 30-45, 45-70 mm), also those of head_surface_effect.json
      orientation_stratum  the orientation strata (g4_epilepsy_adult.ORIENT_BANDS: 0-30, 30-60, 60-90.1 deg
                           from the inner-skull normal), also those of g2_summary.json; the spike locations
                           were drawn among the cortical targets per depth band and orientation stratum
  results/g3b/g3b_targets_<anatomy>_depth.csv
      cortical             0 on the medial wall, which every summary of these heads leaves out
      depth_stratum        the depth strata of g3b_summary.json (configs/g3b_pediatric.toml
                           strata.depth_edges_mm: 0-10 = below 10 mm, 10-15, ..., 60-90 mm)
      spike_depth_band, orientation_stratum  as above (the orientation strata are also those of
                           g3b_summary.json; the spike locations of every head were drawn from these targets)
  results/g3b_constant_gap/targets_<anatomy>_depth.csv
      cortical, depth_stratum (the strata of g3b_constant_gap_summary.json; its table prints the
                           orientation to 0.001 deg, which places every target in its stratum)

g2_patch_targets.csv has the keys, row order and printed depth of g2_targets.csv (checked), so
g2_targets_depth.csv serves it too.

Source of the depth. Neither cached state holds it: cache/g2/state.pkl keeps the summary, arrays,
detectabilities, amplitudes and patch results, cache/g3b/state.pkl the per-anatomy runs (arrays,
placements, detectabilities, geometry) and patches, without target descriptors. Both analysis scripts
recompute the targets, their depth, orientation and region from the anatomy at every run and --replot
(g2_adult_comparison.Study and g3b_pediatric_helmet.load_anatomies, through opmsquid.g2.make_sources, or
g3b_pediatric_helmet.scaled_sources for the scaled adults; depth: the distance from a target to the
nearest vertex of the dense scalp surface). This script does the same, one anatomy at a time, with
study_g3b_constant_gap.load_anatomy (each anatomy as load_anatomies builds it; the adult as both
comparisons build it: g2.make_sources on the sample subject with the seed of configs/g2_adult.toml).

Checks (nothing is written if one fails):
  * the keys and their order equal those of each stored table, and the cortical flag its regions;
  * the depth and the orientation printed as each table prints them equal every stored value;
  * every depth-bin and orientation-stratum count stored in g2_summary.json, head_surface_effect.json and
    covariance_validation.json (all adult targets), g3b_summary.json and g3b_constant_gap_summary.json
    (cortical targets: 'n' and 'n_child' by the head's own bins, 'n_adult' by the adult's, a list named
    '*adult_depth' by the adult's depth at the same vertex of a scaled adult) follows from the companions;
  * the depth percentiles, area shares and area-weighted median depths of g3b_summary.json follow from
    the full-precision depth;
  * where the constant-gap state is in the cache (OPMSQUID_CACHE/g3b_constant_gap/state.pkl), the depth
    and orientation equal those stored there, bit for bit.
It also reports (without failing) the targets printed on the other side of a bin edge, and the
below-8/10-mm counts that results/g3b_children_qc/children_qc.json took from the rounded tables, against
the full-precision counts.

Inputs: the anatomies (OPMSQUID_DATA) and their full-resolution cortex (OPMSQUID_CACHE/anatomy; read,
never written when present). About 6 minutes on one core (the two scaled adults take most), 1 GB.
Usage: PYTHONPATH=src .venv/bin/python scripts/export_target_precision.py [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import pickle
import re
import sys
import time
import tomllib
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import mne  # noqa: E402

import g2_adult_comparison as G2  # noqa: E402  (selects the Agg backend)
import g3b_pediatric_helmet as G3B  # noqa: E402
import g4_epilepsy_adult as G4  # noqa: E402
import study_g3b_constant_gap as CG  # noqa: E402
from opmsquid import io, paths, pediatric as P  # noqa: E402

RES = ROOT / "results"
ANATOMIES = tuple(G3B.ANATOMIES)  # the adult first: the scaled adults are built from it
SCALED = tuple(G3B.SCALED)
NATIVE = tuple(G3B.NATIVE)
SUMMARIES = {"g2_summary": RES / "g2" / "g2_summary.json",
             "head_surface_effect": RES / "g2" / "head_surface_effect.json",
             "covariance_validation": RES / "g2_covariance_validation" / "covariance_validation.json",
             "g3b_summary": RES / "g3b" / "g3b_summary.json",
             "constant_gap_summary": RES / "g3b_constant_gap" / "g3b_constant_gap_summary.json"}
ADULT_SUMMARIES = ("g2_summary", "head_surface_effect", "covariance_validation")
PATCH_TABLE = RES / "g2" / "g2_patch_targets.csv"
CHILDREN_QC = RES / "g3b_children_qc" / "children_qc.json"
CG_STATE = paths.CACHE / "g3b_constant_gap" / "state.pkl"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------------------------
# bins
def edges_of(bands) -> tuple:
    """Edges of contiguous half-open bands ((lo, hi), ...)."""
    for (_, hi), (lo, _) in zip(bands[:-1], bands[1:]):
        if hi != lo:
            raise ValueError(f"bands not contiguous: {bands}")
    return tuple(float(b[0]) for b in bands) + (float(bands[-1][1]),)


ADULT_BINS = tuple(float(x) for x in G2.DEPTH_EDGES)  # the 5-mm bins of the adult depth curves, 10-80 mm
SPIKE_EDGES = edges_of(G4.DEPTH_BANDS)  # the depth bands of the simulated spikes, 10-20, 20-30, 30-45, 45-70 mm
ORIENT_EDGES = edges_of(G4.ORIENT_BANDS)  # orientation strata [deg]: 0-30, 30-60, 60-90.1 (90.1 closes the last at 90)


def strata_edges(cfg: dict, what: str = "depth_edges_mm") -> tuple:
    """The depth (or orientation) strata of the smaller-head summaries (configs/g3b_pediatric.toml)."""
    return tuple(float(x) for x in cfg["strata"][what])


def label(lo: float, hi: float) -> str:
    return f"{lo:g}-{hi:g}"


def bin_labels(values, edges) -> np.ndarray:
    """Label 'lo-hi' of the half-open bin lo <= value < hi holding each value; '' outside every bin."""
    e = np.asarray(edges, float)
    i = np.searchsorted(e, np.asarray(values, float), side="right") - 1
    labs = np.array([label(a, b) for a, b in zip(e[:-1], e[1:])] + [""], dtype=object)
    return labs[np.where((i >= 0) & (i < len(e) - 1), i, len(e) - 1)]


def labels_by_rule(values, edges) -> np.ndarray:
    """The same labels by the analysis scripts' own test, (x >= lo) & (x < hi), bin by bin (a check of ``bin_labels``)."""
    x = np.asarray(values, float)
    out = np.full(len(x), "", dtype=object)
    for lo, hi in zip(edges[:-1], edges[1:]):
        out[(x >= lo) & (x < hi)] = label(lo, hi)
    return out


def exact(x: float) -> str:
    """The shortest decimal form that reads back as the same double."""
    s = repr(float(x))
    if float(s) != float(x):
        raise ValueError(f"{x!r} does not round-trip")
    return s


# ----------------------------------------------------------------------------------------------
# tables
COLUMNS = {"adult": ("hemi", "vertno", "depth_mm", "depth_bin_5mm", "spike_depth_band", "orientation_stratum"),
           "pediatric": ("hemi", "vertno", "cortical", "depth_mm", "depth_stratum", "spike_depth_band", "orientation_stratum"),
           "constant_gap": ("hemi", "vertno", "cortical", "depth_mm", "depth_stratum")}
DEPTH_BINS = {"depth_bin_5mm": lambda strata: ADULT_BINS, "depth_stratum": lambda strata: strata,
              "spike_depth_band": lambda strata: SPIKE_EDGES}  # bins of the depth; orientation_stratum bins the orientation
SUMMARY_OF = {"adult": "g2_summary.json", "pediatric": "g3b_summary.json", "constant_gap": "g3b_constant_gap_summary.json"}


class Spec(NamedTuple):
    """One stored per-target table and its companion."""
    family: str  # 'adult' (adult comparison), 'pediatric' (smaller heads, fixed helmet), 'constant_gap' (helmet at the adult's gap)
    anatomy: str
    table: Path
    decimals: int  # decimals of depth_mm in the stored table

    @property
    def companion(self) -> Path:
        return self.table.with_name(f"{self.table.stem}_depth.csv")

    @property
    def columns(self) -> tuple:
        return COLUMNS[self.family]


def specs() -> list[Spec]:
    out = [Spec("adult", "adult", RES / "g2" / "g2_targets.csv", 2)]
    out += [Spec("pediatric", k, RES / "g3b" / f"g3b_targets_{k}.csv", 2) for k in ANATOMIES]
    out += [Spec("constant_gap", k, RES / "g3b_constant_gap" / f"targets_{k}.csv", 6) for k in ANATOMIES]
    return out


def read_table(path: Path, columns=("hemi", "vertno", "cortical", "depth_mm", "orientation_deg", "region")) -> dict:
    """The '#' lines, the commit of the first one and the requested columns (strings) of a result table."""
    with open(path, newline="") as fh:
        lines = fh.read().splitlines()
    comments = [x for x in lines if x.startswith("#")]
    rows = csv.reader(x for x in lines if not x.startswith("#"))
    head = next(rows)
    idx = {c: head.index(c) for c in columns if c in head}
    out = {c: [] for c in idx}
    for r in rows:
        for c, i in idx.items():
            out[c].append(r[i])
    m = re.search(r"\|\s*commit\s+(\S+)\s*$", comments[0]) if comments else None
    return dict(out, _comments=comments, _commit=m.group(1) if m else "unknown", _header=head)


def load_summaries() -> dict:
    return {k: json.loads(p.read_text()) for k, p in SUMMARIES.items()}


class Checks:
    """Counts the checks per source and keeps the failures."""

    def __init__(self):
        self.n, self.failed = {}, []

    def __call__(self, source: str, ok: bool, msg: str = "") -> None:
        self.n[source] = self.n.get(source, 0) + 1
        if not ok:
            self.failed.append(f"{source}: {msg}")

    @property
    def total(self) -> int:
        return sum(self.n.values())


# ----------------------------------------------------------------------------------------------
# targets, as the analysis scripts build them
class Targets(NamedTuple):
    """The targets of one anatomy as the analysis scripts build them."""
    key: str
    hemi: np.ndarray
    vertno: np.ndarray
    target: np.ndarray  # global full-resolution vertex index
    depth_mm: np.ndarray  # double, as computed
    orientation_deg: np.ndarray
    cortical: np.ndarray
    weights: np.ndarray  # area represented [m^2] (the weights of the area-weighted summaries)


def compute_targets(cfg: dict, g2cfg: dict, keys=ANATOMIES) -> dict[str, Targets]:
    out, adult = {}, None
    for key in keys:
        t0 = time.time()
        an = CG.load_anatomy(key, cfg, g2cfg, adult)
        t = an.src.target
        out[key] = Targets(key, np.asarray(an.cortex.hemi[t], int), np.asarray(an.cortex.vertno[t], int), np.asarray(t).copy(),
                           np.asarray(an.src.depth_mm, float).copy(), np.asarray(an.src.orientation_deg, float).copy(),
                           np.asarray(an.cortical, bool).copy(), np.asarray(an.weights, float).copy())
        if key == "adult":
            adult = an
        del an
        gc.collect()
        log(f"{key}: {len(t)} targets, depth {out[key].depth_mm.min():.3f}-{out[key].depth_mm.max():.3f} mm ({time.time() - t0:.0f} s)")
    return out


def companion_columns(spec: Spec, t: Targets, strata: tuple) -> dict:
    """The companion of ``spec`` as columns of strings (as ``read_table`` returns them)."""
    d = t.depth_mm
    col = {"hemi": [str(x) for x in t.hemi], "vertno": [str(x) for x in t.vertno], "cortical": [str(int(x)) for x in t.cortical],
           "depth_mm": [exact(x) for x in d], "orientation_stratum": list(bin_labels(t.orientation_deg, ORIENT_EDGES))}
    for c, e in DEPTH_BINS.items():
        col[c] = list(bin_labels(d, e(strata)))
    return {c: col[c] for c in spec.columns}


# ----------------------------------------------------------------------------------------------
# checks on the stored files and the companions (also run by tests/test_target_precision.py)
def check_against_table(chk: Checks, spec: Spec, comp: dict, table: dict, strata: tuple) -> None:
    """Keys, row order, printed depth, cortical flag and bin labels of a companion against its stored table."""
    src = f"{spec.companion.name} vs {spec.table.name}"
    n = len(table["depth_mm"])
    chk(src, len(comp["depth_mm"]) == n, f"{len(comp['depth_mm'])} rows, table {n}")
    if len(comp["depth_mm"]) != n:
        return
    chk(src, comp["hemi"] == table["hemi"] and comp["vertno"] == table["vertno"], "keys or row order differ")
    d = np.array([float(x) for x in comp["depth_mm"]])
    bad = [i for i in range(n) if f"{d[i]:.{spec.decimals}f}" != table["depth_mm"][i]]
    chk(src, not bad, f"{len(bad)} depths do not print as stored" + (f" (row {bad[0]}: {comp['depth_mm'][bad[0]]}, stored "
                                                                       f"{table['depth_mm'][bad[0]]})" if bad else ""))
    if "cortical" in comp:
        chk(src, comp["cortical"] == [str(int(not r.endswith("unknown"))) for r in table["region"]],
            "cortical differs from the table's regions (medial wall: '<hemi>.unknown')")
        if "cortical" in table:
            chk(src, comp["cortical"] == table["cortical"], "cortical differs from the table's cortical column")
    for c in spec.columns:
        if c in DEPTH_BINS:
            chk(src, comp[c] == list(labels_by_rule(d, DEPTH_BINS[c](strata))), f"{c} does not follow from depth_mm")
    if "orientation_stratum" in comp:  # the stratum must contain the printed orientation, give or take half its last digit
        h = 0.5 * 10.0 ** -len(table["orientation_deg"][0].split(".")[1])
        o = np.array([float(x) for x in table["orientation_deg"]])
        lab = np.array(comp["orientation_stratum"], dtype=object)
        bad = np.flatnonzero((lab != bin_labels(o - h, ORIENT_EDGES)) & (lab != bin_labels(o + h, ORIENT_EDGES)))
        chk(src, len(bad) == 0, f"orientation_stratum excludes the printed orientation for {len(bad)} targets")


def binned_lists(obj, path=()):
    """(path, rows) of every list of bin rows (dicts with 'lo' and 'hi') in a summary."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from binned_lists(v, path + (str(k),))
    elif isinstance(obj, list):
        if obj and all(isinstance(r, dict) and "lo" in r and "hi" in r for r in obj):
            yield path, obj
        else:
            for i, v in enumerate(obj):
                yield from binned_lists(v, path + (str(i),))


def row_edges(rows) -> tuple | None:
    if any(a["hi"] != b["lo"] for a, b in zip(rows[:-1], rows[1:])):
        return None
    return tuple(float(r["lo"]) for r in rows) + (float(rows[-1]["hi"]),)


def kind_of(path: tuple, edges, kinds: dict) -> str | None:
    """The companion column whose bins a stored list uses: ``kinds`` maps a column to (edges, a word of the path)."""
    for col, (e, word) in kinds.items():
        if edges == e and any(word in p for p in path):
            return col
    return None


def check_adult_counts(chk: Checks, name: str, summary: dict, comp: dict) -> None:
    """Every bin count ('n') of an adult-comparison summary, over all its targets, against the companion's bins."""
    kinds = {"depth_bin_5mm": (ADULT_BINS, "depth"), "spike_depth_band": (SPIKE_EDGES, "depth"),
             "orientation_stratum": (ORIENT_EDGES, "orientation")}
    lab = {c: np.array(comp[c], dtype=object) for c in kinds}
    for path, rows in binned_lists(summary):
        col = kind_of(path, row_edges(rows), kinds)
        if col is None or "n" not in rows[0]:
            continue
        for r in rows:
            want = int(np.sum(lab[col] == label(r["lo"], r["hi"])))
            chk(SUMMARIES[name].name, r["n"] == want, f"{'/'.join(path)} [{label(r['lo'], r['hi'])}] n = {r['n']}, companion {want}")


def check_strata_counts(chk: Checks, name: str, summary: dict, comps: dict[str, dict], strata: tuple) -> None:
    """Every depth-stratum (and, where the companions have it, orientation-stratum) count of a smaller-head summary
    against the companions, over cortical targets: 'n' and 'n_child' by the head's own bins, 'n_adult' by the
    adult's, and a list named '*adult_depth' by the adult's depth at the same vertex (the same hemi and vertno) of
    a scaled adult."""
    kinds = {c: k for c, k in {"depth_stratum": (strata, "depth"), "orientation_stratum": (ORIENT_EDGES, "orientation")}.items()
             if all(c in x for x in comps.values())}
    lab = {(k, c): np.array(x[c], dtype=object) for k, x in comps.items() for c in kinds}
    cort = {k: np.array(x["cortical"]) == "1" for k, x in comps.items()}
    a = comps["adult"]
    row_of = {key: i for i, key in enumerate(zip(a["hemi"], a["vertno"], strict=True))}
    homolog = {k: np.array([row_of[key] for key in zip(comps[k]["hemi"], comps[k]["vertno"], strict=True)]) for k in SCALED if k in comps}
    for path, rows in binned_lists(summary):
        col = kind_of(path, row_edges(rows), kinds)
        keys = [x for x in ("n", "n_child", "n_adult") if x in rows[0]]
        if col is None or not keys:
            continue
        src = SUMMARIES[name].name
        anat = next((tok for p in path for tok in p.split("/") if tok in comps), None)
        by_adult = col == "depth_stratum" and "adult_depth" in path[-1]
        chk(src, anat is not None and (not by_adult or anat in homolog), f"{'/'.join(path)}: no anatomy, or not a scaled adult")
        if anat is None or (by_adult and anat not in homolog):
            continue
        for r in rows:
            s = label(r["lo"], r["hi"])
            for k in keys:
                if k == "n_adult":
                    want = np.sum(cort["adult"] & (lab[("adult", col)] == s))
                elif by_adult:
                    want = np.sum(cort[anat] & (lab[("adult", col)][homolog[anat]] == s))
                else:
                    want = np.sum(cort[anat] & (lab[(anat, col)] == s))
                chk(src, r[k] == int(want), f"{'/'.join(path)} [{s}] {k} = {r[k]}, companion {int(want)}")


def verify(chk: Checks, tables: dict, comps: dict, summaries: dict, strata: tuple, patch: dict | None = None) -> None:
    """All checks that need only the stored files and the companions."""
    for s in specs():
        check_against_table(chk, s, comps[s], tables[s], strata)
    adult = next(s for s in specs() if s.family == "adult")
    if patch is not None:  # g2_patch_targets.csv: the same rows and printed depth as g2_targets.csv
        t = tables[adult]
        chk(f"{PATCH_TABLE.name} vs {adult.table.name}",
            (patch["hemi"], patch["vertno"], patch["depth_mm"]) == (t["hemi"], t["vertno"], t["depth_mm"]), "keys, order or depth differ")
    for name in ADULT_SUMMARIES:
        check_adult_counts(chk, name, summaries[name], comps[adult])
    for fam, name in (("pediatric", "g3b_summary"), ("constant_gap", "constant_gap_summary")):
        check_strata_counts(chk, name, summaries[name], {s.anatomy: comps[s] for s in specs() if s.family == fam}, strata)


# ----------------------------------------------------------------------------------------------
# checks that need the computed targets (script only)
def check_printing(chk: Checks, T: dict[str, Targets], tables: dict) -> None:
    """The computed depth and orientation print as every stored table prints them."""
    for s in specs():
        t, tab = T[s.anatomy], tables[s]
        src = f"{s.table.name} (computed)"
        chk(src, [f"{x:.{s.decimals}f}" for x in t.depth_mm] == tab["depth_mm"], "the computed depth does not print as stored")
        dec = len(tab["orientation_deg"][0].split(".")[1])
        chk(src, [f"{x:.{dec}f}" for x in t.orientation_deg] == tab["orientation_deg"], "the computed orientation does not print as stored")
        chk(src, [str(x) for x in t.hemi] == tab["hemi"] and [str(x) for x in t.vertno] == tab["vertno"], "keys or row order differ")


def check_exact(chk: Checks, T: dict[str, Targets], summaries: dict) -> None:
    """Quantities of the summaries that need the full-precision depth and the area weights."""

    def same(src, what, stored, value, rel=1e-12):
        chk(src, stored == value or abs(stored - value) <= rel * max(abs(stored), abs(value)), f"{what}: stored {stored!r}, "
            f"from the depth {value!r}")

    g2s, g3, cg = summaries["g2_summary"], summaries["g3b_summary"], summaries["constant_gap_summary"]
    same("g2_summary.json", "n_targets", g2s["n_targets"], len(T["adult"].depth_mm), 0)
    for k, t in T.items():
        a = g3["anatomies"][k]
        same("g3b_summary.json", f"anatomies.{k}.n_targets", a["n_targets"], len(t.depth_mm), 0)
        for q, v in zip(("p5", "median", "p95"), np.percentile(t.depth_mm[t.cortical], [5, 50, 95])):
            same("g3b_summary.json", f"anatomies.{k}.target_depth_mm.{q}", a["target_depth_mm"][q], float(v))
        same("g3b_constant_gap_summary.json", f"anatomies.{k}.n_targets", cg["anatomies"][k]["n_targets"], len(t.depth_mm), 0)
        same("g3b_constant_gap_summary.json", f"anatomies.{k}.n_cortical_targets", cg["anatomies"][k]["n_cortical_targets"],
             int(t.cortical.sum()), 0)

    def share(t, lo, hi):  # as g3b_pediatric_helmet.template_depth_checks (finite D: the cortical targets)
        m = t.cortical
        return float(t.weights[m & (t.depth_mm >= lo) & (t.depth_mm < hi)].sum() / t.weights[m].sum())

    a = T["adult"]
    for k in NATIVE:
        e, t, p = g3["template_depth_checks"][k], T[k], f"template_depth_checks.{k}"
        same("g3b_summary.json", f"{p}.area_share_10_20mm", e["area_share_10_20mm"], share(t, 10, 20))
        same("g3b_summary.json", f"{p}.adult_area_share_10_20mm", e["adult_area_share_10_20mm"], share(a, 10, 20))
        same("g3b_summary.json", f"{p}.area_share_deeper_50mm", e["area_share_deeper_50mm"], share(t, 50, np.inf))
        same("g3b_summary.json", f"{p}.median_depth_mm", e["median_depth_mm"], P.weighted_median(t.depth_mm[t.cortical], t.weights[t.cortical]))
        same("g3b_summary.json", f"{p}.adult_median_depth_mm", e["adult_median_depth_mm"],
             P.weighted_median(a.depth_mm[a.cortical], a.weights[a.cortical]))


def check_state(chk: Checks, T: dict[str, Targets]) -> None:
    """The targets stored in the constant-gap state (the one cached state that keeps them), bit for bit."""
    if not CG_STATE.exists():
        log(f"{CG_STATE} not found: the bitwise comparison with the stored depth is skipped")
        return
    with open(CG_STATE, "rb") as fh:
        runs = pickle.load(fh)["runs"]
    for k, t in T.items():
        r = runs.get(k)
        if r is None:
            continue
        for what in ("hemi", "vertno", "target", "depth_mm", "orientation_deg", "cortical", "weights"):
            chk(f"{CG_STATE.parent.name}/{CG_STATE.name}", np.array_equal(np.asarray(r[what]), np.asarray(getattr(t, what))),
                f"runs.{k}.{what} differs from the recomputed one")
    del runs
    gc.collect()


# ----------------------------------------------------------------------------------------------
# what rounding does (reported, not checked)
def report_crossings(T: dict[str, Targets], tables: dict, strata: tuple) -> None:
    """Targets whose printed depth (or orientation) lies in another bin than the computed one."""
    for s in specs():
        t, tab = T[s.anatomy], tables[s]
        lines = []
        bins = [(c, t.depth_mm, tab["depth_mm"], e(strata), "mm") for c, e in DEPTH_BINS.items() if c in s.columns]
        if "orientation_stratum" in s.columns:
            bins.append(("orientation_stratum", t.orientation_deg, tab["orientation_deg"], ORIENT_EDGES, "deg"))
        for c, full, printed, e, unit in bins:
            a, b = bin_labels(full, e), bin_labels([float(x) for x in printed], e)
            for i in np.flatnonzero(a != b):
                lines.append(f"    {c}: {t.hemi[i]}/{t.vertno[i]}{'' if t.cortical[i] else ' (medial wall)'} {float(full[i])!r} {unit} "
                             f"printed {printed[i]}: {a[i] or 'none'}, printed {b[i] or 'none'}")
        if lines:
            log(f"{s.table.name}: {len(lines)} bin memberships differ when read from the printed values\n" + "\n".join(lines))


def report_children_qc(T: dict[str, Targets]) -> None:
    """children_qc.json counts its below-8/10-mm targets on the rounded g3b tables (all targets, depth < limit)."""
    if not CHILDREN_QC.exists():
        return
    for k, a in json.loads(CHILDREN_QC.read_text())["anatomies"].items():
        for key, v in (a.get("g3b_targets") or {}).items():
            m = re.fullmatch(r"below_(\d+(?:\.\d+)?)mm", key)
            if m and k in T:
                full = int(np.sum(T[k].depth_mm < float(m.group(1))))
                log(f"{CHILDREN_QC.name} anatomies.{k}.g3b_targets.{key}.n = {v['n']} (from the rounded table), full precision {full}"
                    + ("" if v["n"] == full else "  <- differs"))


# ----------------------------------------------------------------------------------------------
def status_line(spec: Spec, table_commit: str) -> str:
    also = f"; {PATCH_TABLE.name} has the same rows and printed depth" if spec.family == "adult" else ""
    return (f"Depth below the scalp at full precision and its bins for the targets of {spec.table.name} (commit {table_commit}), "
            f"same rows in the same order{also}")


def legend_line(spec: Spec, strata: tuple) -> str:
    def bins(e):
        return ", ".join(label(a, b) for a, b in zip(e[:-1], e[1:]))

    out = ["cortical: 0 on the medial wall, which every summary of these heads leaves out"] if "cortical" in spec.columns else []
    out.append(f"depth_mm: the depth [mm] as computed, in the shortest decimal form that reads back as the same double "
               f"({spec.table.name} prints it to {spec.decimals} decimals)")
    out.append("bins are half-open, lo <= value < hi (empty: outside every bin)")
    if "depth_bin_5mm" in spec.columns:
        out.append(f"depth_bin_5mm: the 5-mm bins of the depth curves of {SUMMARY_OF['adult']} ({ADULT_BINS[0]:g}-{ADULT_BINS[-1]:g} mm)")
    if "depth_stratum" in spec.columns:
        out.append(f"depth_stratum: the depth strata of {SUMMARY_OF[spec.family]} ({bins(strata)} mm; {label(*strata[:2])} = below "
                   f"{strata[1]:g} mm)")
    if "spike_depth_band" in spec.columns:
        out.append(f"spike_depth_band: the depth bands of the simulated spikes ({bins(SPIKE_EDGES)} mm"
                   + ("; also those of head_surface_effect.json)" if spec.family == "adult" else ")"))
    if "orientation_stratum" in spec.columns:
        out.append(f"orientation_stratum: the orientation strata of {SUMMARY_OF[spec.family]} and of the simulated spikes ({bins(ORIENT_EDGES)} "
                   f"deg from the inner-skull normal; {ORIENT_EDGES[-1]:g} closes the last at 90), from the orientation as computed; the "
                   "spike locations were drawn among the cortical targets per depth band and orientation stratum")
    return "; ".join(out)


def write_companion(spec: Spec, comp: dict, table_commit: str, strata: tuple) -> None:
    with open(spec.companion, "w", newline="") as fh:
        io.csv_status(fh, status_line(spec, table_commit))
        fh.write(f"# {legend_line(spec, strata)}\n")
        wr = csv.writer(fh)
        wr.writerow(spec.columns)
        wr.writerows(zip(*(comp[c] for c in spec.columns), strict=True))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true", help="compute and check, write nothing")
    args = ap.parse_args(argv)
    t_start = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g3b_pediatric.toml").read_text())
    g2cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    strata = strata_edges(cfg)
    summaries = load_summaries()
    stored_cfg = {"g3b_summary": summaries["g3b_summary"]["config"], "constant_gap_summary": summaries["constant_gap_summary"]["config"]["g3b"]}
    for name, c in stored_cfg.items():  # the strata of the stored runs are those of the configuration and of the spike study
        if strata_edges(c) != strata or strata_edges(c, "orientation_edges_deg") != ORIENT_EDGES:
            raise SystemExit(f"{SUMMARIES[name]}: strata differ from configs/g3b_pediatric.toml or from the spike study's")
    tables = {s: read_table(s.table) for s in specs()}
    patch = read_table(PATCH_TABLE, ("hemi", "vertno", "depth_mm"))
    T = compute_targets(cfg, g2cfg)
    comps = {s: companion_columns(s, T[s.anatomy], strata) for s in specs()}

    chk = Checks()
    verify(chk, tables, comps, summaries, strata, patch)
    check_printing(chk, T, tables)
    check_exact(chk, T, summaries)
    check_state(chk, T)
    report_crossings(T, tables, strata)
    report_children_qc(T)
    b = next(s for s in specs() if s.family == "pediatric" and s.anatomy == "childB")
    cort = np.array(comps[b]["cortical"]) == "1"
    printed = np.array([float(x) for x in tables[b]["depth_mm"]])
    log(f"child B, cortical targets below {strata[1]:g} mm: {int(np.sum(cort & (np.array(comps[b]['depth_stratum']) == label(*strata[:2]))))} "
        f"(companion), {int(np.sum(cort & (printed < strata[1])))} printed below {strata[1]:.2f} and "
        f"{int(np.sum(cort & (printed == strata[1])))} printed as {strata[1]:.2f} in {b.table.name}")
    log("checks per source: " + "; ".join(f"{k} {v}" for k, v in chk.n.items() if "_depth.csv" not in k)
        + f"; companions against their tables {sum(v for k, v in chk.n.items() if '_depth.csv' in k)}")
    if chk.failed:
        for x in chk.failed[:60]:
            print("  FAIL", x)
        raise SystemExit(f"{len(chk.failed)} of {chk.total} checks failed; nothing written")
    log(f"all {chk.total} checks passed")
    if args.dry_run:
        return
    for s in specs():
        write_companion(s, comps[s], tables[s]["_commit"], strata)
        back = read_table(s.companion, s.columns)
        if any(back[c] != comps[s][c] for c in s.columns):
            raise SystemExit(f"{s.companion}: read back differs from what was written")
        log(f"wrote {s.companion.relative_to(ROOT)} ({s.companion.stat().st_size / 1024:.0f} kB)")
    log(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
