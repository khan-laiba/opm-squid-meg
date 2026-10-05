#!/usr/bin/env python3
"""Report facts for G3A (Jas et al. head-size benchmark), G3B (fixed adult Neuromag helmet versus
head-adaptive OPM arrays on smaller heads) and the regions by head; prefixes g3a_, g3b_, regh_.

facts(root) -> {name: {"value": text as printed, "raw": unrounded number(s) or text, "source":
"<file> :: <key path>" or "<file> :: derived: <method>"}}, merged by scripts/report_facts.py.
Read only: results/g3a/g3a_size_benchmark.json; results/g3b/g3b_summary.json, g3b_targets_<anatomy>.csv,
school_anatomy_checks.json, child_bem_validation.json, school_subjects_preparation.json; the G3B
configuration as stored in the summary ("config" = configs/g3b_pediatric.toml), configs/g2_adult.toml
(conventions G3B shares with G2) and one history number quoted in docs/methods.md.

Name tokens. Anatomies: adult, school (school-age size control), size2yr (2-year size control),
infant2yr, infant18mo, infant12mo (24-, 18-, 12-month templates), child_a, child_b, child_c.
Groups (ranges): all, smaller (the eight smaller heads), templates_scaled (templates and size
controls), templates, scaled, children (A-C). OPM arrays: dense, matched; Neuromag comparators:
combined (306), grad, mag. Conditions: intr, ib (intrinsic + brain), proj (projected). Metrics:
detectability (no token), peak (peak-channel SNR), meanpow (mean-power SNR). Placements: top (the
primary; implied when no placement token is given), centred, back, xp5mm, xm5mm, yp5mm, ym5mm,
pitchp10, pitchm10, rollp5, rollm5, yawp10, yawm10, xcentred, top18mm, cf (counterfactual helmet
scaled about the centred head), cfx (scaled about the laterally centred head). Quantities: d (median
D in dB), delta (the stored Delta estimator), d_minus_adult (difference of the two medians), _ci
(stored parcel-bootstrap 95 % interval), _share_pos (area share with a positive value), _n;
_min/_max/_range over anatomies or conditions, from unrounded values.
Formats: dB signed, 2 decimals, U+2212 ('0.00' when it rounds to zero); ratios and detectabilities 2
decimals (+ _3dp twin near 1); mm and cm 1 decimal (_2dp twins where finer values are quoted);
counts with thousands separators; percentages integer with % (1 decimal below 10 %; _1dp twins where
one decimal is quoted); area shares as fractions, 2 decimals (usefulness 3, with _pct twins);
configuration constants as configured.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
import tomllib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid.pediatric import weighted_median  # noqa: E402  (the area-weighted median of every G3B summary)
from opmsquid.plotting import DK_LOBES  # noqa: E402

G3A = "results/g3a/g3a_size_benchmark.json"
G3B = "results/g3b/g3b_summary.json"
CHECKS = "results/g3b/school_anatomy_checks.json"
BEMVAL = "results/g3b/child_bem_validation.json"
PREP = "results/g3b/school_subjects_preparation.json"
G2CFG = "configs/g2_adult.toml"
METHODS = "docs/methods.md"

ANATS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
SMALLER = ANATS[1:]
SCALED = ("school", "size2yr")
TEMPLATES = ("infant2yr", "infant18mo", "infant12mo")
SCHOOL = ("childA", "childB", "childC")
NATIVE = TEMPLATES + SCHOOL
TOK = {a: a for a in ANATS} | {"childA": "child_a", "childB": "child_b", "childC": "child_c"}
GROUPS = {"all": ANATS, "smaller": SMALLER, "templates_scaled": SCALED + TEMPLATES, "templates": TEMPLATES,
          "scaled": SCALED, "children": SCHOOL}
GROUP_DESC = {"all": "all nine anatomies", "smaller": "the eight smaller heads",
              "templates_scaled": "the three templates and the two size controls", "templates": "the three infant templates",
              "scaled": "the two size controls", "children": "children A-C"}
MAIN_GROUPS = ("smaller", "children", "templates_scaled")
OPM = {"opm_dense": "dense", "opm_matched": "matched"}
REFS = ("combined", "grad", "mag")
COND = {"intrinsic": "intr", "intrinsic+brain": "ib", "projected": "proj"}
HEADLINE = ("intrinsic+brain", "projected")
METRIC = {"detect": "", "peak": "_peak", "meanpow_db": "_meanpow"}
PLACE = {"centred": "centred", "top": "top", "back": "back", "x+5mm": "xp5mm", "x-5mm": "xm5mm", "y+5mm": "yp5mm",
         "y-5mm": "ym5mm", "pitch+10deg": "pitchp10", "pitch-10deg": "pitchm10", "roll+5deg": "rollp5", "roll-5deg": "rollm5",
         "yaw+10deg": "yawp10", "yaw-10deg": "yawm10", "x-centred": "xcentred", "top-18mm": "top18mm",
         "counterfactual": "cf", "counterfactual_x-centred": "cfx"}
OTHER = ("centred", "back", "x-centred", "top-18mm", "counterfactual", "counterfactual_x-centred")
# the 12 source-blind placements of the placement family (docs/methods.md G3B 'Placements': top, back and the ten
# bounded variants; 'centred', 'x-centred', 'top-18mm' and the counterfactual helmets are not members)
FAMILY = ("top", "back", "x+5mm", "x-5mm", "y+5mm", "y-5mm", "pitch+10deg", "pitch-10deg", "roll+5deg", "roll-5deg",
          "yaw+10deg", "yaw-10deg")
SENS = {"opm_asd_7fT": "asd7", "opm_asd_10fT": "asd10", "opm_asd_15fT": "asd15", "opm_asd_20fT": "asd20",
        "opm_asd_30fT": "asd30", "background_x0.5": "bg0p5", "background_x2": "bg2", "bem1": "bem1"}
SENS_DESC = {"asd7": "OPM noise 7 fT/sqrt(Hz)", "asd10": "OPM noise 10 fT/sqrt(Hz)", "asd15": "OPM noise 15 fT/sqrt(Hz) (primary)",
             "asd20": "OPM noise 20 fT/sqrt(Hz)", "asd30": "OPM noise 30 fT/sqrt(Hz)", "bg0p5": "background variance x0.5",
             "bg2": "background variance x2", "bem1": "1-layer BEM"}
REGION_TOK = {"left frontal": "left_frontal", "left temporal": "left_temporal", "left parietal": "left_parietal",
              "left occipital": "left_occipital", "right frontal": "right_frontal", "right temporal": "right_temporal",
              "right parietal": "right_parietal", "right occipital": "right_occipital"}
# regions by head: the six lobes (the CSV 'lobe' column, plotting.DK_LOBES), four parcels and two parcel groups
REGIONS = {**{lb: ("lobe", lb) for lb in DK_LOBES},
           "superiortemporal": ("parcels", ("superiortemporal",)), "parahippocampal": ("parcels", ("parahippocampal",)),
           "precentral": ("parcels", ("precentral",)), "entorhinal": ("parcels", ("entorhinal",)),
           "mesial_temporal": ("parcels", ("parahippocampal", "entorhinal")),
           "lateral_temporal": ("parcels", ("superiortemporal", "middletemporal", "inferiortemporal", "bankssts",
                                            "transversetemporal"))}
REGH_PLACE = {"top": "squid_top", "cf": "squid_counterfactual"}  # placements with per-target columns in the CSV
PLACE_KEY = {"top": "top", "cf": "counterfactual", "cfx": "counterfactual_x-centred"}

MINUS = "\u2212"


# ----------------------------------------------------------------------------------------------
# formatting
def _minus(text: str) -> str:
    return text.replace("-", MINUS)


def fixed(x, nd=2) -> str:
    """Unsigned fixed point (U+2212 for a negative value)."""
    return _minus(f"{x:.{nd}f}")


def signed(x, nd=2) -> str:
    """dB difference: '+0.44', '−0.38'; '0.00' when the value rounds to zero."""
    return f"{0:.{nd}f}" if round(x, nd) == 0 else _minus(f"{x:+.{nd}f}")


def signed1(x) -> str:
    return signed(x, 1)


def ratio(x) -> str:
    return fixed(x, 2)


def ratio3(x) -> str:
    return fixed(x, 3)


def mm(x) -> str:
    return fixed(x, 1)


def mm2(x) -> str:
    return fixed(x, 2)


def count(n) -> str:
    return f"{int(round(n)):,}"


def pct(x) -> str:
    """Percentage (x in %): integer with '%', one decimal below 10 %."""
    return fixed(x, 1) + "%" if round(abs(x), 1) < 10 else fixed(x, 0) + "%"


def pct1(x) -> str:
    return fixed(x, 1) + "%"


def share(x) -> str:
    return fixed(x, 2)


def share3(x) -> str:
    return fixed(x, 3)


def const(x) -> str:
    """A configuration constant as configured (no added decimals)."""
    return f"{round(float(x), 6):g}"


def interval(fmt, lo, hi) -> str:
    return f"[{fmt(lo)}, {fmt(hi)}]"


def span(fmt, lo, hi) -> str:
    return f"{fmt(lo)} to {fmt(hi)}"


def db_ratio(d) -> float:
    """Amplitude (detectability) ratio of a dB difference."""
    return 10.0 ** (d / 20.0)


# ----------------------------------------------------------------------------------------------
NAME = re.compile(r"(g3a|g3b|regh)(_[a-z0-9]+)+")


def _plain(v):
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    return v


class Facts(dict):
    def add(self, name, value, raw, source):
        if not NAME.fullmatch(name):
            raise ValueError(f"fact name {name!r} is not lower_snake_case with a g3a_/g3b_/regh_ prefix")
        if name in self:
            raise ValueError(f"fact {name!r} defined twice")
        self[name] = {"value": value, "raw": _plain(raw), "source": source}

    def ratio(self, name, x, source, quoted3=False):
        """A ratio (2 decimals) and, when it lies within 0.95-1.05 or is quoted to 3 decimals, its _3dp twin."""
        self.add(name, ratio(x), x, source)
        if quoted3 or 0.95 <= x <= 1.05:
            self.add(f"{name}_3dp", ratio3(x), x, source)

    def interval(self, name, fmt, ci, source):
        if ci:
            self.add(name, interval(fmt, *ci), list(ci), source)

    def range(self, name, items, fmt, file, key, what):
        """name_min, name_max and name_range ('lo to hi') over items {label: raw}, from unrounded values."""
        vals = {k: v for k, v in items.items() if v is not None and math.isfinite(v)}
        if len(vals) < 2:
            return
        klo, khi = min(vals, key=vals.get), max(vals, key=vals.get)
        lo, hi = vals[klo], vals[khi]
        src = f"{file} :: derived: {{}} over {what} of {key}"
        self.add(f"{name}_min", fmt(lo), lo, src.format("min") + f" (at {klo})")
        self.add(f"{name}_max", fmt(hi), hi, src.format("max") + f" (at {khi})")
        self.add(f"{name}_range", span(fmt, lo, hi), [lo, hi], src.format("min and max") + f" (at {klo} and {khi})")

    def group_ranges(self, tpl, series, fmt, file, key, groups=MAIN_GROUPS):
        """Ranges of series {anatomy: raw} over each group; tpl holds '{g}' for the group token."""
        for g in groups:
            self.range(tpl.format(g=g), {a: series[a] for a in GROUPS[g] if a in series}, fmt, file, key, GROUP_DESC[g])


# ----------------------------------------------------------------------------------------------
# G3A
HEAD3A = {"Infant (newborn)": "newborn", "Infant (1-yr)": "infant1yr", "Child (8-yr)": "child8yr", "Adult": "adult"}
SHELL = {"size_following": "follow", "fixed_adult_shell": "fixed"}
SHELL_DESC = {"size_following": "SQUID on the size-following shell s = h + 18 mm",
              "fixed_adult_shell": "every head concentric in one fixed adult shell (s = 113 mm)"}
ETA3A = {"eta2.5": "eta2p5", "eta3": "eta3", "eta4": "eta4", "eta6": "eta6"}


def g3a_facts(F: Facts, s: dict) -> None:
    for head, h in HEAD3A.items():
        row = s["size_following"][head]
        F.add(f"g3a_{h}_h", const(row["h_mm"]), row["h_mm"], f"{G3A} :: size_following['{head}'].h_mm (Table 1 scalp radius, mm)")
        F.add(f"g3a_{h}_b", const(row["b_mm"]), row["b_mm"], f"{G3A} :: size_following['{head}'].b_mm (Table 1 brain radius, mm)")
        F.add(f"g3a_{h}_layer", const(row["layer_mm"]), row["layer_mm"],
              f"{G3A} :: size_following['{head}'].layer_mm (extracerebral layer h - b, mm)")
        for shell, sh in SHELL.items():
            r = s[shell][head]
            F.add(f"g3a_{h}_{sh}_xi", const(r["xi_squid_mm"]), r["xi_squid_mm"],
                  f"{G3A} :: {shell}['{head}'].xi_squid_mm (SQUID shell standoff above the scalp, i.e. the gap, mm)")
            for q, desc in (("eta0", "noise ratio below which the OPM is ahead at every depth (centre limit)"),
                            ("eta1", "noise ratio above which the SQUID is ahead at every depth (brain surface)")):
                src = f"{G3A} :: {shell}['{head}'].{q} ({desc}; {SHELL_DESC[shell]})"
                F.add(f"g3a_{h}_{sh}_{q}", ratio(r[q]), r[q], src)
                F.add(f"g3a_{h}_{sh}_{q}_1dp", fixed(r[q], 1), r[q], src)
            for e, et in ETA3A.items():
                v, key, base = r[e], f"{shell}['{head}']['{e}']", f"g3a_{h}_{sh}_{et}"
                vol = f"{G3A} :: {key}.volume_fraction_pct (brain-sphere volume with SNR_OPM > SNR_SQUID, %; {SHELL_DESC[shell]})"
                F.add(f"{base}_volume", pct(v["volume_fraction_pct"]), v["volume_fraction_pct"], vol)
                F.add(f"{base}_volume_1dp", pct1(v["volume_fraction_pct"]), v["volume_fraction_pct"], vol)
                if v["d_eq_exact_mm"] is None:  # no crossing: one system ahead at every depth (see eta0 / eta1)
                    continue
                F.add(f"{base}_deq", mm(v["d_eq_exact_mm"]), v["d_eq_exact_mm"],
                      f"{G3A} :: {key}.d_eq_exact_mm (equal-SNR depth below the scalp, exact root of Eq. 3, mm)")
                F.add(f"{base}_deq_grid", mm(v["d_eq_grid_mm"]), v["d_eq_grid_mm"],
                      f"{G3A} :: {key}.d_eq_grid_mm (the paper's grid rule, mm)")
                norm = f"{G3A} :: {key}.normalized_exact_pct ((d_eq - (h - b)) / b, Fig. 5B 'normalized d_eq', %)"
                F.add(f"{base}_norm", pct(v["normalized_exact_pct"]), v["normalized_exact_pct"], norm)
                F.add(f"{base}_norm_1dp", pct1(v["normalized_exact_pct"]), v["normalized_exact_pct"], norm)
                F.add(f"{base}_norm_grid", pct(v["normalized_grid_pct"]), v["normalized_grid_pct"],
                      f"{G3A} :: {key}.normalized_grid_pct (grid rule, %)")
                F.add(f"{base}_deq_over_b", pct(v["d_eq_over_b_pct"]), v["d_eq_over_b_pct"], f"{G3A} :: {key}.d_eq_over_b_pct (%)")
                below = v["d_eq_exact_mm"] - r["layer_mm"]
                F.add(f"{base}_below_brain", mm(below), below,
                      f"{G3A} :: derived: {key}.d_eq_exact_mm - {shell}['{head}'].layer_mm (equal-SNR depth below the brain surface, mm)")
    f3 = {head: s["size_following"][head]["eta3"] for head in HEAD3A}
    what = "the four Table 1 heads (size-following shell, eta = 3)"
    key = "size_following['<head>']['eta3']"
    F.range("g3a_follow_eta3_deq", {h: f3[h]["d_eq_exact_mm"] for h in HEAD3A}, mm, G3A, f"{key}.d_eq_exact_mm", what)
    F.range("g3a_follow_eta3_norm_1dp", {h: f3[h]["normalized_exact_pct"] for h in HEAD3A}, pct1, G3A,
            f"{key}.normalized_exact_pct", what)
    F.range("g3a_follow_eta3_volume_1dp", {h: f3[h]["volume_fraction_pct"] for h in HEAD3A}, pct1, G3A,
            f"{key}.volume_fraction_pct", what)
    F.range("g3a_follow_eta3_below_brain", {h: f3[h]["d_eq_exact_mm"] - s["size_following"][h]["layer_mm"] for h in HEAD3A},
            mm, G3A, f"{key}.d_eq_exact_mm - size_following['<head>'].layer_mm", what)
    F.range("g3a_layer", {h: s["size_following"][h]["layer_mm"] for h in HEAD3A}, const, G3A,
            "size_following['<head>'].layer_mm", "the four Table 1 heads")
    radii = {round(s["fixed_adult_shell"][h]["h_mm"] + s["fixed_adult_shell"][h]["xi_squid_mm"], 6) for h in HEAD3A}
    if len(radii) != 1:
        raise ValueError(f"{G3A}: the fixed shell radius differs between heads {radii}")
    radius = radii.pop()
    F.add("g3a_fixed_shell_radius", const(radius), radius,
          f"{G3A} :: derived: fixed_adult_shell['<head>'].h_mm + xi_squid_mm (the same for every head, mm)")
    for head, h in (("Infant (newborn)", "newborn"), ("Adult", "adult")):
        v = s["printed_checks"][head]["printed_pct"]
        F.add(f"g3a_printed_{h}_norm", pct(v), v, f"{G3A} :: printed_checks['{head}'].printed_pct (value printed in the paper, %)")
    F.add("g3a_commit", s["provenance"]["commit"], s["provenance"]["commit"], f"{G3A} :: provenance.commit")


# ----------------------------------------------------------------------------------------------
# G3B: anatomies, arrays and placements
def anatomy_facts(F: Facts, s: dict) -> None:
    an, arr, cfg = s["anatomies"], s["arrays"], s["config"]
    adult_area = an["adult"]["cortical_area_cm2"]
    for nt, key in (("scaled", "size-only control"), ("templates", "infant template"), ("children", "school-aged child")):
        ks = [a for a, v in an.items() if key in v["description"]]
        F.add(f"g3b_n_{nt}_heads", count(len(ks)), len(ks),
              f"{G3B} :: derived: number of anatomies whose description contains '{key}' ({', '.join(ks)})")
    for a in ANATS:
        t, x = TOK[a], an[a]
        hs = x["head_size"]
        F.add(f"g3b_{t}_ofc_cm", fixed(hs["ofc_mm"] / 10, 1), hs["ofc_mm"] / 10,
              f"{G3B} :: anatomies['{a}'].head_size.ofc_mm / 10 (occipitofrontal circumference, cm)")
        for q in ("ofc_mm", "breadth_mm", "length_mm", "vertex_height_mm", "inter_auricular_mm"):
            F.add(f"g3b_{t}_{q}", mm(hs[q]), hs[q], f"{G3B} :: anatomies['{a}'].head_size.{q}")
        F.add(f"g3b_{t}_cap_volume_cm3", count(hs["cap_volume_cm3"]), hs["cap_volume_cm3"],
              f"{G3B} :: anatomies['{a}'].head_size.cap_volume_cm3 (cap volume above the fiducial plane, cm^3)")
        F.add(f"g3b_{t}_ofc_ratio", ratio3(x["ofc_ratio_to_adult"]), x["ofc_ratio_to_adult"],
              f"{G3B} :: anatomies['{a}'].ofc_ratio_to_adult (3 decimals, as the scale factors are quoted)")
        F.add(f"g3b_{t}_cortex_area_cm2", count(x["cortical_area_cm2"]), x["cortical_area_cm2"],
              f"{G3B} :: anatomies['{a}'].cortical_area_cm2 (usable cortical area, cm^2)")
        share_area = 100 * x["cortical_area_cm2"] / adult_area
        F.add(f"g3b_{t}_cortex_area_pct_of_adult", pct(share_area), share_area,
              f"{G3B} :: derived: anatomies['{a}'].cortical_area_cm2 / anatomies['adult'].cortical_area_cm2 (%; with the background "
              "fixed per unit area, also the share of the adult's background power)")
        F.add(f"g3b_{t}_n_targets", count(x["n_targets"]), x["n_targets"], f"{G3B} :: anatomies['{a}'].n_targets")
        mw = s["medial_wall_targets"][a]
        F.add(f"g3b_{t}_medial_wall_targets", count(mw), mw, f"{G3B} :: medial_wall_targets['{a}'] (FreeSurfer 'unknown', excluded)")
        F.add(f"g3b_{t}_n_cortical_targets", count(x["n_targets"] - mw), x["n_targets"] - mw,
              f"{G3B} :: derived: anatomies['{a}'].n_targets - medial_wall_targets['{a}'] (targets in every G3B summary)")
        F.add(f"g3b_{t}_n_background_grid", count(x["n_background_grid"]), x["n_background_grid"],
              f"{G3B} :: anatomies['{a}'].n_background_grid")
        for q in ("p5", "median", "p95"):
            v = x["target_depth_mm"][q]
            F.add(f"g3b_{t}_target_depth_{q}", mm(v), v,
                  f"{G3B} :: anatomies['{a}'].target_depth_mm.{q} (unweighted, cortical targets, mm below the scalp)")
        for o, ot in OPM.items():
            meta = arr[a][o]
            F.add(f"g3b_{t}_{ot}_sites", count(meta["n"]), meta["n"], f"{G3B} :: arrays['{a}']['{o}'].n")
            F.add(f"g3b_{t}_{ot}_moved_out", count(meta["n_moved_out"]), meta["n_moved_out"],
                  f"{G3B} :: arrays['{a}']['{o}'].n_moved_out (sites moved outward by the clearance rule)")
            F.add(f"g3b_{t}_{ot}_min_spacing", mm(meta["min_spacing_mm"]), meta["min_spacing_mm"],
                  f"{G3B} :: arrays['{a}']['{o}'].min_spacing_mm")
            geo = s["sensor_distances"][a][o]
            src = f"{G3B} :: sensor_distances['{a}']['{o}'].median_mm (median sensing-centre height above the MRI scalp, mm)"
            F.add(f"g3b_{t}_{ot}_standoff", mm(geo["median_mm"]), geo["median_mm"], src)
            F.add(f"g3b_{t}_{ot}_standoff_2dp", mm2(geo["median_mm"]), geo["median_mm"], src)
            src = f"{G3B} :: sensor_distances['{a}']['{o}'].exact_cell_clearance_min_mm (whole-cell distance to the head surface, mm)"
            F.add(f"g3b_{t}_{ot}_cell_clearance", mm(geo["exact_cell_clearance_min_mm"]), geo["exact_cell_clearance_min_mm"], src)
            F.add(f"g3b_{t}_{ot}_cell_clearance_3dp", fixed(geo["exact_cell_clearance_min_mm"], 3),
                  geo["exact_cell_clearance_min_mm"], src)
        F.add(f"g3b_{t}_dense_median_spacing", mm(arr[a]["opm_dense"]["median_spacing_mm"]), arr[a]["opm_dense"]["median_spacing_mm"],
              f"{G3B} :: arrays['{a}']['opm_dense'].median_spacing_mm")
        cov = s["sensor_distances"][a]["opm_coverage"]
        F.add(f"g3b_{t}_covered_scalp_cm2", count(cov["scalp_above_brow_plane_cm2"]), cov["scalp_above_brow_plane_cm2"],
              f"{G3B} :: sensor_distances['{a}'].opm_coverage.scalp_above_brow_plane_cm2")
        for o, ot in OPM.items():
            v = cov[f"{o}_sites_per_100cm2"]
            F.add(f"g3b_{t}_{ot}_sites_per_100cm2", fixed(v, 1), v, f"{G3B} :: sensor_distances['{a}'].opm_coverage.{o}_sites_per_100cm2")
    for i, c in enumerate(cfg["anatomy"]["school"]):
        k = c["key"]
        F.add(f"g3b_{TOK[k]}_age", fixed(c["age_y"], 1), c["age_y"], f"{G3B} :: config.anatomy.school[{i}].age_y (years)")
        F.add(f"g3b_{TOK[k]}_subject", c["subject"], c["subject"], f"{G3B} :: config.anatomy.school[{i}].subject (OpenNeuro ds005234)")
    for k, name in zip(TEMPLATES, cfg["anatomy"]["templates"]):
        F.add(f"g3b_{TOK[k]}_template", name, name, f"{G3B} :: config.anatomy.templates (O'Reilly et al. 2021)")
    F.range("g3b_children_age", {c["key"]: c["age_y"] for c in cfg["anatomy"]["school"]}, lambda v: fixed(v, 1), G3B,
            "config.anatomy.school[].age_y", GROUP_DESC["children"])
    key = "anatomies['<anat>']"
    for g in ("all", "smaller", "templates", "children", "templates_scaled"):
        F.range(f"g3b_{g}_ofc_cm", {a: an[a]["head_size"]["ofc_mm"] / 10 for a in GROUPS[g]}, lambda v: fixed(v, 1), G3B,
                f"{key}.head_size.ofc_mm / 10", GROUP_DESC[g])
        F.range(f"g3b_{g}_cortex_area_cm2", {a: an[a]["cortical_area_cm2"] for a in GROUPS[g]}, count, G3B,
                f"{key}.cortical_area_cm2", GROUP_DESC[g])
        F.range(f"g3b_{g}_cortex_area_pct_of_adult", {a: 100 * an[a]["cortical_area_cm2"] / adult_area for a in GROUPS[g]}, pct,
                G3B, f"{key}.cortical_area_cm2 / anatomies['adult'].cortical_area_cm2 (%)", GROUP_DESC[g])
        for o, ot in OPM.items():
            F.range(f"g3b_{g}_{ot}_sites", {a: arr[a][o]["n"] for a in GROUPS[g]}, count, G3B, f"arrays['<anat>']['{o}'].n",
                    GROUP_DESC[g])
            F.range(f"g3b_{g}_{ot}_standoff_2dp", {a: s["sensor_distances"][a][o]["median_mm"] for a in GROUPS[g]}, mm2, G3B,
                    f"sensor_distances['<anat>']['{o}'].median_mm", GROUP_DESC[g])
    for a in SMALLER:
        row = s["comparisons"][f"{a}/opm_dense/combined/intrinsic+brain/detect"]["delta_by_depth"][0]
        F.add(f"g3b_{TOK[a]}_shallow_targets_lt10mm", count(row["n_child"]), row["n_child"],
              f"{G3B} :: comparisons['{a}/opm_dense/combined/intrinsic+brain/detect'].delta_by_depth[0].n_child "
              f"(cortical targets {row['lo']:g}-{row['hi']:g} mm below the scalp)")
    row = s["comparisons"]["school/opm_dense/combined/intrinsic+brain/detect"]["delta_by_depth"][0]
    F.add("g3b_adult_shallow_targets_lt10mm", count(row["n_adult"]), row["n_adult"],
          f"{G3B} :: comparisons['school/opm_dense/combined/intrinsic+brain/detect'].delta_by_depth[0].n_adult")
    F.add("g3b_n_anatomies", count(len(ANATS)), len(ANATS), f"{G3B} :: derived: number of keys of anatomies")


def placement_facts(F: Facts, s: dict) -> None:
    pl, sd = s["placements"], s["sensor_distances"]
    n = {len(pl[a]) for a in ANATS}
    if n != {len(PLACE)} or any(set(pl[a]) != set(PLACE) for a in ANATS):
        raise ValueError(f"{G3B}: unexpected placements {n}")
    F.add("g3b_n_placements", count(len(PLACE)), len(PLACE), f"{G3B} :: derived: number of keys of placements['<anat>'] (every anatomy)")
    F.add("g3b_n_family_placements", count(len(FAMILY)), len(FAMILY),
          f"{G3B} :: derived: placement family = top, back and the ten bounded variants of placements['<anat>'] "
          "(docs/methods.md G3B 'Placements')")
    infeasible = [(a, p) for a in ANATS for p in s["infeasible_placements"][a]]
    F.add("g3b_n_infeasible_placements", count(len(infeasible)), len(infeasible),
          f"{G3B} :: derived: number of entries of infeasible_placements (" + ", ".join(f"{a} {p}" for a, p in infeasible) + ")")
    for a in ANATS:
        t = TOK[a]
        for p, pt in PLACE.items():
            g = sd[a][f"squid:{p}"]
            src = f"{G3B} :: sensor_distances['{a}']['squid:{p}']"
            F.add(f"g3b_{t}_{pt}_gap", mm(g["median_mm"]), g["median_mm"],
                  f"{src}.median_mm (median magnetometer coil centre to scalp, mm)")
            F.add(f"g3b_{t}_{pt}_gap_nearest", mm(g["min_mm"]), g["min_mm"],
                  f"{src}.min_mm (nearest magnetometer coil centre to scalp, mm)")
            if "moved_mm" in pl[a][p]:
                F.add(f"g3b_{t}_{pt}_moved", mm(pl[a][p]["moved_mm"]), pl[a][p]["moved_mm"],
                      f"{G3B} :: placements['{a}']['{p}'].moved_mm (head moved from the start pose, mm)")
            if p in ("centred", "top"):
                for reg, rt in REGION_TOK.items():
                    v = g["regions"][reg]
                    F.add(f"g3b_{t}_{pt}_gap_{rt}", mm(v), v, f"{src}.regions['{reg}'] (median gap in the Vectorview selection, mm)")
        for p, pt in (("counterfactual", "cf"), ("counterfactual_x-centred", "cfx")):
            for q in ("k", "k_nominal"):
                v = pl[a][p][q]
                F.add(f"g3b_{t}_{pt}_{q}", ratio3(v), v,
                      f"{G3B} :: placements['{a}']['{p}'].{q} (helmet scale factor; k_nominal = head-circumference ratio)")
        v = pl[a]["x-centred"]["shift_x_mm"]
        F.add(f"g3b_{t}_xcentred_shift", signed1(v), v,
              f"{G3B} :: placements['{a}']['x-centred'].shift_x_mm (lateral shift along device x, mm; negative = to the left)")
    for p, pt in PLACE.items():
        series = {a: sd[a][f"squid:{p}"]["median_mm"] for a in ANATS}
        F.group_ranges(f"g3b_{{g}}_{pt}_gap", series, mm, G3B, f"sensor_distances['<anat>']['squid:{p}'].median_mm")
    F.group_ranges("g3b_{g}_top_moved", {a: pl[a]["top"]["moved_mm"] for a in ANATS}, mm, G3B,
                   "placements['<anat>']['top'].moved_mm")
    for p, pt in (("counterfactual", "cf"), ("counterfactual_x-centred", "cfx")):
        F.group_ranges(f"g3b_{{g}}_{pt}_k", {a: pl[a][p]["k"] for a in ANATS}, ratio3, G3B, f"placements['<anat>']['{p}'].k")
        adult = sd["adult"][f"squid:{p}"]["median_mm"]
        below = {a: sd[a][f"squid:{p}"]["median_mm"] for a in SMALLER if sd[a][f"squid:{p}"]["median_mm"] < adult}
        F.add(f"g3b_{pt}_gap_below_adult_count", count(len(below)), len(below),
              f"{G3B} :: derived: number of the eight smaller heads with sensor_distances['<anat>']['squid:{p}'].median_mm below the "
              f"adult's ({', '.join(below)})")
        F.range(f"g3b_{pt}_gap_below_adult", below, mm, G3B, f"sensor_distances['<anat>']['squid:{p}'].median_mm",
                "the smaller heads whose gap is below the adult's")


# ----------------------------------------------------------------------------------------------
# G3B: D, Delta and their summaries
def _delta_desc(c: str) -> str:
    if c in SCALED:
        return "vertex-wise: area-weighted median over the common cortical vertices of D_child - D_adult"
    return ("adult-area-weighted median over Desikan-Killiany parcels (>= 10 targets in both anatomies) of the parcel "
            "difference D_child - D_adult")


def _delta_share_desc(c: str) -> str:
    if c in SCALED:
        return "area share of the common vertices with Delta > 0"
    return "share of the parcels' adult area with parcel Delta > 0"


def d_median_facts(F: Facts, s: dict) -> None:
    dm = s["D_median_dB"]
    for o, ot in OPM.items():
        for r in REFS:
            for k, kt in COND.items():
                for m, mt in METRIC.items():
                    stem = f"{ot}_vs_{r}_{kt}{mt}"
                    key = f"{o}/{r}/{k}/{m}"
                    series = {a: dm[f"{a}/{key}"] for a in ANATS}
                    for a, v in series.items():
                        src = (f"{G3B} :: D_median_dB['{a}/{key}'] (area-weighted median over cortical targets, medial wall "
                               "excluded; Neuromag at the primary placement 'top')")
                        F.add(f"g3b_{TOK[a]}_{stem}_d", signed(v), v, src)
                        if m != "meanpow_db":
                            F.ratio(f"g3b_{TOK[a]}_{stem}_d_ratio", db_ratio(v), src.replace(" :: ", " :: derived: 10^(D/20) of ", 1))
                        if a != "adult":
                            diff = v - series["adult"]
                            F.add(f"g3b_{TOK[a]}_{stem}_d_minus_adult", signed(diff), diff,
                                  f"{G3B} :: derived: D_median_dB['{a}/{key}'] - D_median_dB['adult/{key}'] (difference of the medians)")
                    if m == "meanpow_db" or k == "intrinsic":  # ranges only for the series the report ranges over
                        continue
                    F.group_ranges(f"g3b_{{g}}_{stem}_d", series, signed, G3B, f"D_median_dB['<anat>/{key}']")
                    F.group_ranges(f"g3b_{{g}}_{stem}_d_minus_adult", {a: series[a] - series["adult"] for a in SMALLER}, signed,
                                   G3B, f"D_median_dB['<anat>/{key}'] - D_median_dB['adult/{key}']")
                    F.group_ranges(f"g3b_{{g}}_{stem}_d_ratio", {a: db_ratio(v) for a, v in series.items()}, ratio, G3B,
                                   f"10^(D_median_dB['<anat>/{key}'] / 20)", groups=("smaller",))
    key = "opm_dense/combined/intrinsic+brain/peak"
    neg = [a for a in SMALLER if dm[f"{a}/{key}"] < 0]
    F.add("g3b_smaller_dense_vs_combined_ib_peak_d_negative_count", count(len(neg)), len(neg),
          f"{G3B} :: derived: number of the eight smaller heads with D_median_dB['<anat>/{key}'] < 0 ({', '.join(neg)})")
    lg = s["link_to_g2"]
    F.ratio("g3b_link_g2_ratio_all_targets", lg["adult_centred_unweighted_all_targets"],
            f"{G3B} :: link_to_g2.adult_centred_unweighted_all_targets (dense/combined, i+b, adult at the measured position, "
            "unweighted median over all targets: G2's headline)", quoted3=True)
    F.ratio("g3b_link_g2_ratio_cortical", lg["adult_centred_unweighted_cortical"],
            f"{G3B} :: link_to_g2.adult_centred_unweighted_cortical (as above without the medial wall)", quoted3=True)
    F.add("g3b_link_g2_db", signed(lg["adult_centred_area_weighted_cortical_dB"]), lg["adult_centred_area_weighted_cortical_dB"],
          f"{G3B} :: link_to_g2.adult_centred_area_weighted_cortical_dB (area-weighted, medial wall excluded: the G3B convention)")


def comparison_facts(F: Facts, s: dict) -> None:
    comp = s["comparisons"]
    for o, ot in OPM.items():
        for r in REFS:
            for k in HEADLINE:
                for m, mt in METRIC.items():
                    key = f"{o}/{r}/{k}/{m}"
                    stem = f"{ot}_vs_{r}_{COND[k]}{mt}"
                    for c in SMALLER:
                        e, t = comp[f"{c}/{key}"], TOK[c]
                        path = f"{G3B} :: comparisons['{c}/{key}']"
                        F.interval(f"g3b_{t}_{stem}_d_ci", signed, e["d_child"]["ci95"],
                                   f"{path}.d_child.ci95 (parcel-bootstrap 95 % interval of the area-weighted median D)")
                        F.add(f"g3b_{t}_{stem}_d_share_pos", share(e["d_child"]["share_positive"]), e["d_child"]["share_positive"],
                              f"{path}.d_child.share_positive (area share of the cortical targets with D > 0: OPM ahead)")
                        F.add(f"g3b_{t}_{stem}_delta", signed(e["delta"]["median"]), e["delta"]["median"],
                              f"{path}.delta.median ({_delta_desc(c)})")
                        F.interval(f"g3b_{t}_{stem}_delta_ci", signed, e["delta"]["ci95"], f"{path}.delta.ci95 (parcel bootstrap)")
                        F.add(f"g3b_{t}_{stem}_delta_share_pos", share(e["delta"]["share_positive"]), e["delta"]["share_positive"],
                              f"{path}.delta.share_positive ({_delta_share_desc(c)})")
                        nk = "n" if c in SCALED else "n_parcels"
                        F.add(f"g3b_{t}_{stem}_delta_n", count(e["delta"][nk]), e["delta"][nk],
                              f"{path}.delta.{nk} ({'common cortical vertices' if c in SCALED else 'parcels'})")
                    e = comp[f"school/{key}"]
                    if len({comp[f"{c}/{key}"]["d_adult"]["share_positive"] for c in SMALLER}) != 1:
                        raise ValueError(f"{G3B}: d_adult.share_positive differs between comparisons of {key}")
                    F.add(f"g3b_adult_{stem}_d_share_pos", share(e["d_adult"]["share_positive"]), e["d_adult"]["share_positive"],
                          f"{G3B} :: comparisons['school/{key}'].d_adult.share_positive (the same in every comparison)")
                    if e["d_adult"]["ci95"]:
                        cis = [comp[f"{c}/{key}"]["d_adult"]["ci95"] for c in SMALLER]
                        F.interval(f"g3b_adult_{stem}_d_ci", signed, e["d_adult"]["ci95"],
                                   f"{G3B} :: comparisons['school/{key}'].d_adult.ci95 (the adult's interval as stored with the "
                                   "school-age comparison, as docs/methods.md quotes it; each comparison resamples the adult afresh: "
                                   f"bounds {min(x[0] for x in cis):.3f}-{max(x[0] for x in cis):.3f} / "
                                   f"{min(x[1] for x in cis):.3f}-{max(x[1] for x in cis):.3f} dB)")
                    if m == "meanpow_db":
                        continue
                    F.group_ranges(f"g3b_{{g}}_{stem}_delta", {c: comp[f"{c}/{key}"]["delta"]["median"] for c in SMALLER}, signed,
                                   G3B, f"comparisons['<anat>/{key}'].delta.median")
                    if m != "detect":
                        continue
                    F.group_ranges(f"g3b_{{g}}_{stem}_d_share_pos", {c: comp[f"{c}/{key}"]["d_child"]["share_positive"] for c in SMALLER},
                                   share, G3B, f"comparisons['<anat>/{key}'].d_child.share_positive")
                    F.range(f"g3b_smaller_{stem}_delta_ci_hi", {c: comp[f"{c}/{key}"]["delta"]["ci95"][1] for c in SMALLER}, signed,
                            G3B, f"comparisons['<anat>/{key}'].delta.ci95[1]", GROUP_DESC["smaller"])
    strata_facts(F, s)


def _tok(lo, hi, unit):
    return f"{lo:g}_{math.floor(hi):g}{unit}"


def strata_facts(F: Facts, s: dict) -> None:
    """Delta by depth and by orientation stratum (dense vs combined, i+b), and the scaled controls' vertex-wise
    Delta by the adult's depth and by orientation."""
    comp = s["comparisons"]
    key, stem = "opm_dense/combined/intrinsic+brain/detect", "dense_vs_combined_ib"
    for c in SMALLER:
        e, t = comp[f"{c}/{key}"], TOK[c]
        for field, kind, unit in (("delta_by_depth", "depth", "mm"), ("delta_by_orientation_strata", "orient", "deg")):
            for i, row in enumerate(e[field]):
                base = f"g3b_{t}_{stem}_{kind}_{_tok(row['lo'], row['hi'], unit)}"
                path = f"{G3B} :: comparisons['{c}/{key}'].{field}[{i}]"
                stratum = f"{row['lo']:g}-{row['hi']:g} {unit} below the scalp" if unit == "mm" else \
                    f"{row['lo']:g}-{row['hi']:g} deg from radial"
                for q in ("n_child", "n_adult"):
                    F.add(f"{base}_{q}", count(row[q]), row[q], f"{path}.{q} (cortical targets, {stratum})")
                if "delta" not in row:
                    continue
                F.add(f"{base}_delta", signed(row["delta"]), row["delta"],
                      f"{path}.delta (difference of the area-weighted medians, {stratum})")
                F.interval(f"{base}_delta_ci", signed, row["ci95"], f"{path}.ci95 (parcels resampled within each anatomy)")
                F.add(f"{base}_d_child", signed(row["d_child"]), row["d_child"], f"{path}.d_child")
                F.add(f"{base}_d_adult", signed(row["d_adult"]), row["d_adult"], f"{path}.d_adult")
        for field, kind, unit in (("delta_by_adult_depth", "adultdepth", "mm"), ("delta_by_orientation", "vorient", "deg")):
            for i, row in enumerate(e.get(field, [])):
                base = f"g3b_{t}_{stem}_{kind}_{_tok(row['lo'], row['hi'], unit)}"
                path = f"{G3B} :: comparisons['{c}/{key}'].{field}[{i}]"
                F.add(f"{base}_n", count(row["n"]), row["n"], f"{path}.n")
                if "median" in row:
                    F.add(f"{base}_vdelta", signed(row["median"]), row["median"],
                          f"{path}.median (vertex-wise Delta, area-weighted median in the stratum)")
                    F.interval(f"{base}_vdelta_ci", signed, row["ci95"], f"{path}.ci95")
    rows = {c: comp[f"{c}/{key}"]["delta_by_depth"] for c in SMALLER}

    def pick(cs, lo, hi):
        return {f"{c} {r['lo']:g}-{r['hi']:g} mm": r["delta"] for c in cs for r in rows[c]
                if "delta" in r and lo <= r["lo"] and r["hi"] <= hi}

    path = f"comparisons['<anat>/{key}'].delta_by_depth[].delta"
    F.range("g3b_children_dense_vs_combined_ib_depth_10_30mm_delta", pick(SCHOOL, 10, 30), signed, G3B, path,
            "children A-C and the strata 10-15 to 25-30 mm")
    F.range("g3b_children_dense_vs_combined_ib_depth_30_60mm_delta", pick(SCHOOL, 30, 60), signed, G3B, path,
            "children A-C and the strata 30-40 to 50-60 mm")
    F.range("g3b_templates_scaled_dense_vs_combined_ib_depth_delta", pick(SCALED + TEMPLATES, 0, 90), signed, G3B, path,
            "the templates and size controls and every populated stratum")
    F.range("g3b_children_dense_vs_combined_ib_depth_delta", pick(SCHOOL, 0, 90), signed, G3B, path,
            "children A-C and every populated stratum")
    below = [c for c in SCHOOL for r in rows[c] if r["lo"] == 20 and r["hi"] == 25 and "ci95" in r and r["ci95"][1] < 0]
    F.add("g3b_children_depth_20_25mm_ci_below_zero_count", count(len(below)), len(below),
          f"{G3B} :: derived: number of children A-C whose comparisons['<anat>/{key}'].delta_by_depth (20-25 mm) ci95 lies "
          f"entirely below 0 ({', '.join(below)})")
    # the children's shallow strata (10-15 to 25-30 mm): how many child-stratum cells, how many negative, and the
    # largest lower interval bound (at or below 0 means no interval lies entirely above 0)
    cells = [(c, r) for c in SCHOOL for r in rows[c] if "delta" in r and 10 <= r["lo"] and r["hi"] <= 30]
    sh = f"{G3B} :: derived: children A-C and the strata 10-15 to 25-30 mm of comparisons['<anat>/{key}'].delta_by_depth[]"
    neg = [f"{c} {r['lo']:g}-{r['hi']:g} mm" for c, r in cells if r["delta"] < 0]
    F.add("g3b_children_depth_10_30mm_cells_n", count(len(cells)), len(cells), f"{sh} (child-stratum cells with a delta)")
    F.add("g3b_children_depth_10_30mm_delta_negative_count", count(len(neg)), len(neg), f"{sh}.delta < 0 ({', '.join(neg)})")
    lo_max = max(cells, key=lambda cr: cr[1]["ci95"][0])
    F.add("g3b_children_depth_10_30mm_ci_lo_max", signed(lo_max[1]["ci95"][0]), lo_max[1]["ci95"][0],
          f"{sh}.ci95[0], largest lower bound (at {lo_max[0]} {lo_max[1]['lo']:g}-{lo_max[1]['hi']:g} mm)")


def other_placement_facts(F: Facts, s: dict) -> None:
    """Delta at the other placements (each against the adult at the same rule) and the derived counterfactual contrasts."""
    dop, plc = s["delta_other_placements"], s["placement_D"]
    for p in OTHER:
        pt = PLACE[p]
        for r in REFS:
            stem = f"{pt}_dense_vs_{r}_ib"
            key = f"{p}_vs_adult_{p}/{r}"
            for c in SMALLER:
                e, t = dop[f"{c}/{key}"], TOK[c]
                path = f"{G3B} :: delta_other_placements['{c}/{key}']"
                F.interval(f"g3b_{t}_{stem}_d_ci", signed, e["d_child"]["ci95"], f"{path}.d_child.ci95 (parcel bootstrap)")
                F.add(f"g3b_{t}_{stem}_d_share_pos", share(e["d_child"]["share_positive"]), e["d_child"]["share_positive"],
                      f"{path}.d_child.share_positive (area share with D > 0)")
                F.add(f"g3b_{t}_{stem}_delta", signed(e["delta"]["median"]), e["delta"]["median"],
                      f"{path}.delta.median ({_delta_desc(c)})")
                F.interval(f"g3b_{t}_{stem}_delta_ci", signed, e["delta"]["ci95"], f"{path}.delta.ci95 (parcel bootstrap)")
                F.add(f"g3b_{t}_{stem}_delta_share_pos", share(e["delta"]["share_positive"]), e["delta"]["share_positive"],
                      f"{path}.delta.share_positive ({_delta_share_desc(c)})")
                if p in ("counterfactual", "counterfactual_x-centred"):
                    d_c = e["d_child"]["median"]
                    for ref_p, ref_t in (("top", "top"), (p, pt)):
                        d_a = plc[f"adult/{ref_p}/{r}/intrinsic+brain"]["median"]
                        F.add(f"g3b_{t}_{stem}_d_minus_adult_{ref_t}", signed(d_c - d_a), d_c - d_a,
                              f"{G3B} :: derived: delta_other_placements['{c}/{key}'].d_child.median - "
                              f"placement_D['adult/{ref_p}/{r}/intrinsic+brain'].median (difference of the medians; "
                              + ("against the primary adult reference, top contact)" if ref_p == "top" else "against the adult's own "
                                 "counterfactual, which has factor 1)"))
            e = dop[f"school/{key}"]
            F.interval(f"g3b_adult_{stem}_d_ci", signed, e["d_adult"]["ci95"],
                       f"{G3B} :: delta_other_placements['school/{key}'].d_adult.ci95 (the adult's interval as stored with the "
                       "school-age comparison)")
            F.add(f"g3b_adult_{stem}_d_share_pos", share(e["d_adult"]["share_positive"]), e["d_adult"]["share_positive"],
                  f"{G3B} :: delta_other_placements['school/{key}'].d_adult.share_positive")
            F.group_ranges(f"g3b_{{g}}_{stem}_delta", {c: dop[f"{c}/{key}"]["delta"]["median"] for c in SMALLER}, signed, G3B,
                           f"delta_other_placements['<anat>/{key}'].delta.median")
            F.range(f"g3b_smaller_{stem}_delta_ci_hi", {c: dop[f"{c}/{key}"]["delta"]["ci95"][1] for c in SMALLER}, signed, G3B,
                    f"delta_other_placements['<anat>/{key}'].delta.ci95[1]", GROUP_DESC["smaller"])
            nonneg = [c for c in SMALLER if dop[f"{c}/{key}"]["delta"]["median"] >= 0]
            F.add(f"g3b_smaller_{stem}_delta_nonneg_count", count(len(nonneg)), len(nonneg),
                  f"{G3B} :: derived: number of the eight smaller heads with delta_other_placements['<anat>/{key}'].delta.median >= 0 "
                  f"({', '.join(nonneg) or 'none'})")
            if p in ("counterfactual", "counterfactual_x-centred"):
                for ref_p, ref_t in (("top", "top"), (p, pt)):
                    d_a = plc[f"adult/{ref_p}/{r}/intrinsic+brain"]["median"]
                    F.group_ranges(f"g3b_{{g}}_{stem}_d_minus_adult_{ref_t}",
                                   {c: dop[f"{c}/{key}"]["d_child"]["median"] - d_a for c in SMALLER}, signed, G3B,
                                   f"delta_other_placements['<anat>/{key}'].d_child.median - "
                                   f"placement_D['adult/{ref_p}/{r}/intrinsic+brain'].median")


def placement_d_facts(F: Facts, s: dict) -> None:
    """Median D per anatomy and placement (dense OPM vs each comparator), by lobe (i+b), and the placement family."""
    plc = s["placement_D"]
    for a in ANATS:
        t = TOK[a]
        for p, pt in PLACE.items():
            for r in REFS:
                for k in HEADLINE:
                    e = plc[f"{a}/{p}/{r}/{k}"]
                    path = f"{G3B} :: placement_D['{a}/{p}/{r}/{k}']"
                    stem = f"g3b_{t}_{pt}_dense_vs_{r}_{COND[k]}"
                    F.add(f"{stem}_d", signed(e["median"]), e["median"],
                          f"{path}.median (area-weighted median D, medial wall excluded)")
                    if k == "intrinsic+brain" and r in ("combined", "mag"):
                        for lb, v in e["by_lobe"].items():
                            F.add(f"{stem}_d_{lb}", signed(v), v, f"{path}.by_lobe['{lb}'] (area-weighted median D over the lobe)")
    for p, pt in PLACE.items():
        for k in HEADLINE:
            F.group_ranges(f"g3b_{{g}}_{pt}_dense_vs_combined_{COND[k]}_d",
                           {a: plc[f"{a}/{p}/combined/{k}"]["median"] for a in ANATS}, signed, G3B,
                           f"placement_D['<anat>/{p}/combined/{k}'].median")
    for r in ("combined", "mag"):  # the counterfactual helmets by lobe, over children A-C and every lobe
        for p in ("counterfactual", "counterfactual_x-centred"):
            pt = PLACE[p]
            F.range(f"g3b_children_{pt}_dense_vs_{r}_ib_d_lobes",
                    {f"{c} {lb}": v for c in SCHOOL for lb, v in plc[f"{c}/{p}/{r}/intrinsic+brain"]["by_lobe"].items()}, signed, G3B,
                    f"placement_D['<anat>/{p}/{r}/intrinsic+brain'].by_lobe", "children A-C and the six lobes")
    # placement family (dense vs combined, i+b): the 12 source-blind placements, feasible ones only
    fam = {}
    for a in ANATS:
        t = TOK[a]
        vals = {p: plc[f"{a}/{p}/combined/intrinsic+brain"]["median"] for p in FAMILY if p not in s["infeasible_placements"][a]}
        fam[a] = float(np.median(list(vals.values())))
        key = f"placement_D['{a}/<placement>/combined/intrinsic+brain'].median"
        what = f"the {len(vals)} feasible family placements of {a} (top, back, +-5 mm x/y, pitch +-10, roll +-5, yaw +-10 deg)"
        F.add(f"g3b_{t}_family_n", count(len(vals)), len(vals), f"{G3B} :: derived: number of {what}")
        F.add(f"g3b_{t}_family_dense_vs_combined_ib_d_median", signed(fam[a]), fam[a],
              f"{G3B} :: derived: median (numpy) over {what} of {key}")
        F.range(f"g3b_{t}_family_dense_vs_combined_ib_d", vals, signed, G3B, key, what)
        F.add(f"g3b_{t}_family_dense_vs_combined_ib_d_span", fixed(max(vals.values()) - min(vals.values())),
              max(vals.values()) - min(vals.values()), f"{G3B} :: derived: max - min over {what} of {key} (dB, unsigned)")
        order = sorted(vals, key=vals.get)
        F.add(f"g3b_{t}_top_rank_in_family", count(order.index("top") + 1), order.index("top") + 1,
              f"{G3B} :: derived: rank of 'top' from the lowest over {what} of {key}")
    ranks = {a: F[f"g3b_{TOK[a]}_top_rank_in_family"]["raw"] for a in ANATS}
    F.group_ranges("g3b_{g}_top_rank_in_family", ranks, count, G3B,
                   "the rank of 'top' among the feasible family placements (placement_D combined i+b)")
    ext = {p: plc[f"adult/{p}/combined/intrinsic+brain"]["median"] for p in FAMILY + ("x-centred", "top-18mm")}
    F.range("g3b_adult_family_ext_dense_vs_combined_ib_d", ext, signed, G3B,
            "placement_D['adult/<placement>/combined/intrinsic+brain'].median", "the adult's 12 family placements, x-centred and top-18mm")
    d_top = {a: plc[f"{a}/top/combined/intrinsic+brain"]["median"] for a in ANATS}
    shrink = {}
    for c in SMALLER:
        t = TOK[c]
        fd, td = fam[c] - fam["adult"], d_top[c] - d_top["adult"]
        shrink[c] = td - fd
        F.add(f"g3b_{t}_family_d_minus_adult", signed(fd), fd,
              f"{G3B} :: derived: family median of placement_D['{c}/<placement>/combined/intrinsic+brain'].median minus the adult's")
        F.add(f"g3b_{t}_family_vs_top_shrink", fixed(td - fd), td - fd,
              f"{G3B} :: derived: (D_top child - D_top adult) - (family median child - family median adult), placement_D "
              "combined i+b (dB; positive = the child-minus-adult difference is smaller against the family medians)")
    for g in ("smaller", "children", "templates_scaled"):
        F.range(f"g3b_{g}_family_vs_top_shrink", {c: shrink[c] for c in GROUPS[g]}, fixed, G3B,
                "(top difference - family-median difference), placement_D combined i+b", GROUP_DESC[g])
        F.range(f"g3b_{g}_family_d_minus_adult", {c: fam[c] - fam["adult"] for c in GROUPS[g]}, signed, G3B,
                "family median child - family median adult, placement_D combined i+b", GROUP_DESC[g])
    spread = max(ext[p] for p in FAMILY) - min(ext[p] for p in FAMILY)
    prim = s["comparisons"]
    small = [c for c in SMALLER if prim[f"{c}/opm_dense/combined/intrinsic+brain/detect"]["delta"]["median"] < spread]
    F.add("g3b_smaller_delta_below_adult_family_span_count", count(len(small)), len(small),
          f"{G3B} :: derived: number of the eight primary Deltas (comparisons['<anat>/opm_dense/combined/intrinsic+brain/detect']"
          f".delta.median) below the adult's family span ({', '.join(small)})")
    cen = d_top["adult"] - plc["adult/centred/combined/intrinsic+brain"]["median"]
    src = (f"{G3B} :: derived: placement_D['adult/centred/combined/intrinsic+brain'].median - "
           "placement_D['adult/top/combined/intrinsic+brain'].median (the adult's measured pose against top contact)")
    F.add("g3b_adult_centred_minus_top_d", signed(-cen), -cen, src)
    F.add("g3b_adult_centred_minus_top_d_abs", fixed(abs(cen)), abs(cen), src + " (magnitude)")


# ----------------------------------------------------------------------------------------------
# G3B: mechanism, depth mix, channel count, sensitivity, patches, usefulness
def mechanism_facts(F: Facts, s: dict) -> None:
    ab, vw, amp, nz = s["absolute_detectability_dB"], s["vertexwise_change_dB"], s["amplitude_median_fT"], s["noise_rms_median"]
    systems = (("opm_dense/opm", "dense"), ("opm_matched/opm", "matched"), ("squid:top/combined", "combined"),
               ("squid:top/grad", "grad"), ("squid:top/mag", "mag"))
    for a in ANATS:
        t = TOK[a]
        for key, st in systems:
            for k in HEADLINE:
                v = ab[f"{a}/{key}/{k}"]
                F.add(f"g3b_{t}_absdet_{st}_{COND[k]}", signed(v), v,
                      f"{G3B} :: absolute_detectability_dB['{a}/{key}/{k}'] (median 20 log10 d of a 10-nAm dipole, dB)")
        for key, st in (("opm_dense", "dense"), ("opm_matched", "matched"), ("squid:top", "mag")):
            v = amp[a][key]
            F.add(f"g3b_{t}_peak_field_{st}", count(v), v,
                  f"{G3B} :: amplitude_median_fT['{a}']['{key}'] (median over targets of the largest magnetometer field of a 10-nAm "
                  "dipole, fT)")
        rr = amp[a]["opm_dense"] / amp[a]["squid:top"]
        F.ratio(f"g3b_{t}_peak_field_ratio_dense_mag", rr,
                f"{G3B} :: derived: amplitude_median_fT['{a}']['opm_dense'] / amplitude_median_fT['{a}']['squid:top']")
        for key, nt in (("squid:top/mag", "top_mag"), ("squid:top/grad", "top_grad"), ("squid:centred/mag", "centred_mag"),
                        ("squid:centred/grad", "centred_grad"), ("opm_dense/mag", "dense"), ("opm_matched/mag", "matched")):
            for comp in ("brain", "intrinsic"):
                v = nz[a][key][comp]
                grad = key.endswith("grad")
                F.add(f"g3b_{t}_noise_{nt}_{comp}", fixed(v, 1) if grad else count(v), v,
                      f"{G3B} :: noise_rms_median['{a}']['{key}'].{comp} (median per-channel RMS in the band, "
                      + ("fT/cm)" if grad else "fT)"))
    for c in SCALED:
        for key, st in systems[:1] + systems[2:]:
            for k in HEADLINE:
                v = vw[f"{c}/{key}/{k}"]
                F.add(f"g3b_{TOK[c]}_vwchange_{st}_{COND[k]}", signed(v), v,
                      f"{G3B} :: vertexwise_change_dB['{c}/{key}/{k}'] (vertex-wise change of the system's own detectability "
                      "from the adult)")
    for key, st in (("opm_dense", "dense"), ("opm_matched", "matched"), ("squid:top", "mag")):
        F.group_ranges(f"g3b_{{g}}_peak_field_{st}", {a: amp[a][key] for a in ANATS}, count, G3B, f"amplitude_median_fT['<anat>']['{key}']")
    for key, nt in (("squid:top/mag", "top_mag"), ("squid:top/grad", "top_grad"), ("opm_dense/mag", "dense")):
        fmt = (lambda v: fixed(v, 1)) if key.endswith("grad") else count
        F.group_ranges(f"g3b_{{g}}_noise_{nt}_brain", {a: nz[a][key]["brain"] for a in ANATS}, fmt, G3B,
                       f"noise_rms_median['<anat>']['{key}'].brain")
    for key, st in systems:
        F.group_ranges(f"g3b_{{g}}_absdet_{st}_ib", {a: ab[f"{a}/{key}/intrinsic+brain"] for a in ANATS}, signed, G3B,
                       f"absolute_detectability_dB['<anat>/{key}/intrinsic+brain']")


def depth_check_facts(F: Facts, s: dict) -> None:
    tdc = s["template_depth_checks"]
    path = f"{G3B} :: template_depth_checks"
    for k in NATIVE:
        e, t = tdc[k], TOK[k]
        F.add(f"g3b_{t}_pooled_diff", signed(e["pooled_difference_db"]), e["pooled_difference_db"],
              f"{path}['{k}'].pooled_difference_db (difference of the pooled area-weighted medians, dense vs combined, i+b)")
        F.add(f"g3b_{t}_reweighted_diff", signed(e["depth_reweighted_difference_db"]), e["depth_reweighted_difference_db"],
              f"{path}['{k}'].depth_reweighted_difference_db (child targets reweighted to the adult's area share per depth stratum)")
        for q, nt in (("area_share_10_20mm", "area_share_10_20mm"), ("area_share_deeper_50mm", "area_share_gt50mm"),
                      ("radial_share", "radial_share")):
            F.add(f"g3b_{t}_{nt}", share(e[q]), e[q], f"{path}['{k}'].{q} (fraction of the cortical area)")
            F.add(f"g3b_{t}_{nt}_pct", pct(100 * e[q]), 100 * e[q], f"{path}['{k}'].{q} x 100 (%)")
        F.add(f"g3b_{t}_depth_wmedian", mm(e["median_depth_mm"]), e["median_depth_mm"],
              f"{path}['{k}'].median_depth_mm (area-weighted median depth below the scalp, mm)")
        for i, row in enumerate(e["orientation_at_matched_depth"]):
            for lab in ("radial", "tangential"):
                if row[lab] is not None:
                    F.add(f"g3b_{t}_matched_depth_{row['lo']:g}_{row['hi']:g}mm_{lab}", signed(row[lab]), row[lab],
                          f"{path}['{k}'].orientation_at_matched_depth[{i}].{lab} (difference of the medians at matched depth, dB)")
    if len({tdc[k]["adult_area_share_10_20mm"] for k in NATIVE}) != 1 or len({tdc[k]["adult_median_depth_mm"] for k in NATIVE}) != 1:
        raise ValueError(f"{G3B}: the adult's entries of template_depth_checks differ between anatomies")
    e = tdc["infant2yr"]
    F.add("g3b_adult_area_share_10_20mm", share(e["adult_area_share_10_20mm"]), e["adult_area_share_10_20mm"],
          f"{path}['infant2yr'].adult_area_share_10_20mm (the same in every entry)")
    F.add("g3b_adult_area_share_10_20mm_pct", pct(100 * e["adult_area_share_10_20mm"]), 100 * e["adult_area_share_10_20mm"],
          f"{path}['infant2yr'].adult_area_share_10_20mm x 100 (%)")
    F.add("g3b_adult_depth_wmedian", mm(e["adult_median_depth_mm"]), e["adult_median_depth_mm"],
          f"{path}['infant2yr'].adult_median_depth_mm (area-weighted, mm)")
    for g in ("templates", "children"):
        for q, nt, fmt in (("pooled_difference_db", "pooled_diff", signed), ("depth_reweighted_difference_db", "reweighted_diff", signed),
                           ("area_share_10_20mm", "area_share_10_20mm", share), ("median_depth_mm", "depth_wmedian", mm)):
            F.range(f"g3b_{g}_{nt}", {k: tdc[k][q] for k in GROUPS[g]}, fmt, G3B, f"template_depth_checks['<anat>'].{q}", GROUP_DESC[g])
        F.range(f"g3b_{g}_area_share_10_20mm_pct", {k: 100 * tdc[k]["area_share_10_20mm"] for k in GROUPS[g]}, pct, G3B,
                "template_depth_checks['<anat>'].area_share_10_20mm x 100", GROUP_DESC[g])


def channel_count_facts(F: Facts, s: dict) -> None:
    """The adult's dense array subsampled to each child's site count; the site-count penalty and the counterfactual Delta plus it."""
    cc, comp, dop = s["channel_count_control"], s["comparisons"], s["delta_other_placements"]
    pen, plus = {}, {}
    for c in SMALLER:
        t = TOK[c]
        F.add(f"g3b_{t}_eqcount_sites", count(cc[f"{c}/combined"]["n_sites"]), cc[f"{c}/combined"]["n_sites"],
              f"{G3B} :: channel_count_control['{c}/combined'].n_sites")
        for r in REFS:
            e, path = cc[f"{c}/{r}"], f"{G3B} :: channel_count_control['{c}/{r}']"
            stem = f"g3b_{t}_eqcount_dense_vs_{r}_ib"
            F.add(f"{stem}_d_adult_sub", signed(e["d_adult_subsampled"]["median"]), e["d_adult_subsampled"]["median"],
                  f"{path}.d_adult_subsampled.median (adult's dense array subsampled by farthest-point sampling to {e['n_sites']} sites)")
            F.interval(f"{stem}_d_adult_sub_ci", signed, e["d_adult_subsampled"]["ci95"], f"{path}.d_adult_subsampled.ci95")
            F.add(f"{stem}_delta", signed(e["delta"]["median"]), e["delta"]["median"],
                  f"{path}.delta.median (Delta at equal channel count; {_delta_desc(c)})")
            F.interval(f"{stem}_delta_ci", signed, e["delta"]["ci95"], f"{path}.delta.ci95")
            F.add(f"{stem}_delta_share_pos", share(e["delta"]["share_positive"]), e["delta"]["share_positive"],
                  f"{path}.delta.share_positive ({_delta_share_desc(c)})")
            full = comp[f"{c}/opm_dense/{r}/intrinsic+brain/detect"]["d_adult"]["median"]
            p = full - e["d_adult_subsampled"]["median"]
            pen[(c, r)] = p
            F.add(f"g3b_{t}_site_penalty_vs_{r}", fixed(p), p,
                  f"{G3B} :: derived: comparisons['{c}/opm_dense/{r}/intrinsic+brain/detect'].d_adult.median - "
                  f"channel_count_control['{c}/{r}'].d_adult_subsampled.median (dB the adult loses with {e['n_sites']} sites; unsigned)")
            for pl, pt in (("counterfactual", "cf"), ("counterfactual_x-centred", "cfx")):
                d = dop[f"{c}/{pl}_vs_adult_{pl}/{r}"]["delta"]["median"]
                plus[(c, r, pt)] = d + p
                F.add(f"g3b_{t}_{pt}_dense_vs_{r}_ib_delta_plus_penalty", signed(d + p), d + p,
                      f"{G3B} :: derived: delta_other_placements['{c}/{pl}_vs_adult_{pl}/{r}'].delta.median + the site-count penalty "
                      f"(g3b_{t}_site_penalty_vs_{r})")
    for r in REFS:
        F.group_ranges(f"g3b_{{g}}_eqcount_dense_vs_{r}_ib_delta", {c: cc[f"{c}/{r}"]["delta"]["median"] for c in SMALLER}, signed,
                       G3B, f"channel_count_control['<anat>/{r}'].delta.median")
        F.group_ranges(f"g3b_{{g}}_eqcount_dense_vs_{r}_ib_d_adult_sub",
                       {c: cc[f"{c}/{r}"]["d_adult_subsampled"]["median"] for c in SMALLER}, signed, G3B,
                       f"channel_count_control['<anat>/{r}'].d_adult_subsampled.median")
        F.group_ranges(f"g3b_{{g}}_site_penalty_vs_{r}", {c: pen[(c, r)] for c in SMALLER}, fixed, G3B,
                       f"comparisons['<anat>/opm_dense/{r}/intrinsic+brain/detect'].d_adult.median - "
                       f"channel_count_control['<anat>/{r}'].d_adult_subsampled.median")
        for pt in ("cf", "cfx"):
            series = {c: plus[(c, r, pt)] for c in SMALLER}
            pl = "counterfactual" if pt == "cf" else "counterfactual_x-centred"
            key = f"delta_other_placements['<anat>/{pl}_vs_adult_{pl}/{r}'].delta.median + site-count penalty"
            F.group_ranges(f"g3b_{{g}}_{pt}_dense_vs_{r}_ib_delta_plus_penalty", series, signed, G3B, key)
            positive = {c: v for c, v in series.items() if v > 0}
            F.add(f"g3b_{pt}_dense_vs_{r}_ib_delta_plus_penalty_positive_count", count(len(positive)), len(positive),
                  f"{G3B} :: derived: number of the eight smaller heads with {key} > 0 ({', '.join(positive)})")
            F.range(f"g3b_{pt}_dense_vs_{r}_ib_delta_plus_penalty_positive", positive, signed, G3B, key,
                    "the smaller heads where it is positive")


def sensitivity_facts(F: Facts, s: dict) -> None:
    sens = s["sensitivity_median_D_dB"]
    for key, v in sens.items():
        parts = key.split("/")
        if parts[0] == "delta":
            _, c, var, o, r, k = parts
            F.add(f"g3b_{TOK[c]}_{SENS[var]}_{OPM[o]}_vs_{r}_{COND[k]}_d_minus_adult", signed(v), v,
                  f"{G3B} :: sensitivity_median_D_dB['{key}'] (child minus adult difference of the medians, {SENS_DESC[SENS[var]]} "
                  "applied to both)")
        else:
            a, var, o, r, k = parts
            F.add(f"g3b_{TOK[a]}_{SENS[var]}_{OPM[o]}_vs_{r}_{COND[k]}_d", signed(v), v,
                  f"{G3B} :: sensitivity_median_D_dB['{key}'] (median D, {SENS_DESC[SENS[var]]}; Neuromag at top contact)")
    asd = [v for v in SENS.values() if v.startswith("asd")]
    var_of = {v: kk for kk, v in SENS.items()}
    for o, ot in OPM.items():
        for r in ("combined",):  # ranges for the comparator the report ranges over
            for k in HEADLINE:
                stem = f"{ot}_vs_{r}_{COND[k]}"
                for var, vt in SENS.items():
                    F.group_ranges(f"g3b_{{g}}_{vt}_{stem}_d", {a: sens[f"{a}/{var}/{o}/{r}/{k}"] for a in ANATS}, signed, G3B,
                                   f"sensitivity_median_D_dB['<anat>/{var}/{o}/{r}/{k}']")
                    F.group_ranges(f"g3b_{{g}}_{vt}_{stem}_d_minus_adult", {c: sens[f"delta/{c}/{var}/{o}/{r}/{k}"] for c in SMALLER},
                                   signed, G3B, f"sensitivity_median_D_dB['delta/<anat>/{var}/{o}/{r}/{k}']")
                for a in ANATS:
                    t = TOK[a]
                    F.range(f"g3b_{t}_asd_{stem}_d", {f"{vt}": sens[f"{a}/{var_of[vt]}/{o}/{r}/{k}"] for vt in asd}, signed, G3B,
                            f"sensitivity_median_D_dB['{a}/opm_asd_<n>fT/{o}/{r}/{k}']", "the OPM noise levels 7-30 fT/sqrt(Hz)")
                    if a != "adult":
                        F.range(f"g3b_{t}_asd_{stem}_d_minus_adult", {vt: sens[f"delta/{a}/{var_of[vt]}/{o}/{r}/{k}"] for vt in asd},
                                signed, G3B, f"sensitivity_median_D_dB['delta/{a}/opm_asd_<n>fT/{o}/{r}/{k}']",
                                "the OPM noise levels 7-30 fT/sqrt(Hz)")
    for a in ANATS:  # how far the background variants move D itself
        prim = s["D_median_dB"][f"{a}/opm_dense/combined/intrinsic+brain/detect"]
        mv = max(abs(sens[f"{a}/{var}/opm_dense/combined/intrinsic+brain"] - prim) for var in ("background_x0.5", "background_x2"))
        F.add(f"g3b_{TOK[a]}_bg_d_largest_shift", fixed(mv), mv,
              f"{G3B} :: derived: max over background_x0.5 and background_x2 of |sensitivity_median_D_dB['{a}/<variant>/opm_dense/"
              f"combined/intrinsic+brain'] - D_median_dB['{a}/opm_dense/combined/intrinsic+brain/detect']| (dB)")


def patch_facts(F: Facts, s: dict) -> None:
    pt, cfg = s["patches_median_D_dB"], s["config"]["sources"]
    radii = [f"{r:g}" for r in cfg["patch_radii_mm"]]
    for a in ANATS:
        t = TOK[a]
        for r in radii:
            F.add(f"g3b_{t}_patch{r}mm_area", fixed(pt[f"{a}/area_cm2/{r}mm"], 2), pt[f"{a}/area_cm2/{r}mm"],
                  f"{G3B} :: patches_median_D_dB['{a}/area_cm2/{r}mm'] (median patch area, cm^2)")
            F.add(f"g3b_{t}_patch{r}mm_n_centres", count(pt[f"{a}/n_centres/{r}mm"]), pt[f"{a}/n_centres/{r}mm"],
                  f"{G3B} :: patches_median_D_dB['{a}/n_centres/{r}mm']")
            for ref in REFS:
                for k in HEADLINE:
                    key = f"{a}/{r}mm/{ref}/{k}"
                    F.add(f"g3b_{t}_patch{r}mm_dense_vs_{ref}_{COND[k]}_d", signed(pt[key]), pt[key],
                          f"{G3B} :: patches_median_D_dB['{key}'] (median D over the patch centres, fixed 10-nAm total)")
                    if a != "adult":
                        dk = f"delta/{a}/{r}mm/{ref}/{k}"
                        F.add(f"g3b_{t}_patch{r}mm_dense_vs_{ref}_{COND[k]}_d_minus_adult", signed(pt[dk]), pt[dk],
                              f"{G3B} :: patches_median_D_dB['{dk}'] (child minus adult difference of the medians)")
            for lab, lt in (("squid", "combined"), ("opm_dense", "dense")):
                for k in HEADLINE:
                    base = f"{a}/fixed_density/{r}mm/{lab}/{k}"
                    v = pt[f"{base}/median_detectability"]
                    F.add(f"g3b_{t}_patch{r}mm_fixdens_{lt}_{COND[k]}_det", ratio(v), v,
                          f"{G3B} :: patches_median_D_dB['{base}/median_detectability'] (median detectability at "
                          f"{s['config']['sources']['fixed_density_nAm_per_mm2']:g} nAm/mm^2)")
                    v = pt[f"{base}/share_usable"]
                    src = f"{G3B} :: patches_median_D_dB['{base}/share_usable'] (share of patch centres with d >= 5)"
                    F.add(f"g3b_{t}_patch{r}mm_fixdens_{lt}_{COND[k]}_usable", share(v), v, src)
                    F.add(f"g3b_{t}_patch{r}mm_fixdens_{lt}_{COND[k]}_usable_pct", pct(100 * v), 100 * v, src + " x 100 (%)")
        for ref in REFS:
            for k in HEADLINE:
                for tag, rt in (("focal_same_centres", "patch"), ("focal_same_centres_20mm", "patch20mm")):
                    key = f"{a}/{tag}/{ref}/{k}"
                    F.add(f"g3b_{t}_{rt}_focal_dense_vs_{ref}_{COND[k]}_d", signed(pt[key]), pt[key],
                          f"{G3B} :: patches_median_D_dB['{key}'] (focal D at the same centres)")


def usefulness_facts(F: Facts, s: dict) -> None:
    us = s["usefulness"]
    for a in ANATS:
        t = TOK[a]
        for k in HEADLINE:
            for q in s["config"]["usefulness"]["reference_moments_nAm"]:
                key = f"combined/{k}/{q:g}nAm"
                for cat, v in us[a][key].items():
                    src = (f"{G3B} :: usefulness['{a}']['{key}'].{cat} (share of usable cortical area, dense OPM vs Neuromag combined, "
                           "d >= 5)")
                    F.add(f"g3b_{t}_useful_{COND[k]}_{q:g}nam_{cat}", share3(v), v, src)
                    F.add(f"g3b_{t}_useful_{COND[k]}_{q:g}nam_{cat}_pct", pct(100 * v), 100 * v, src + " x 100 (%)")
        for name, nt in (("opm_dense", "dense"), ("combined", "combined")):
            for i, row in enumerate(us[a][f"q_threshold_vs_depth/{name}/intrinsic+brain"]):
                if "median" in row:
                    F.add(f"g3b_{t}_qthr_{nt}_ib_{_tok(row['lo'], row['hi'], 'mm')}", count(row["median"]), row["median"],
                          f"{G3B} :: usefulness['{a}']['q_threshold_vs_depth/{name}/intrinsic+brain'][{i}].median (moment for d = 5, "
                          "nAm, area-weighted median in the depth stratum)")
    for cat in ("both", "opm_only", "squid_only", "neither"):
        F.group_ranges(f"g3b_{{g}}_useful_ib_100nam_{cat}", {a: us[a]["combined/intrinsic+brain/100nAm"][cat] for a in ANATS}, share3, G3B,
                       f"usefulness['<anat>']['combined/intrinsic+brain/100nAm'].{cat}")


# ----------------------------------------------------------------------------------------------
# G3B: the school-aged children's anatomy checks, the child BEM validation and the preparation
def school_check_facts(F: Facts, root: Path) -> None:
    ch = json.loads((root / CHECKS).read_text())
    names = {"adult": "adult", "ANTS2-0Years3T": "infant2yr", "ANTS18-0Months3T": "infant18mo", "ANTS12-0Months3T": "infant12mo",
             **{f"{c}/{v}": f"{TOK[c]}_{v}" for c in SCHOOL for v in ("watershed", "modelled")}}
    isd = ch["inner_skull_scalp_depth_upper_head"]
    for key, nt in names.items():
        e = isd[key]
        for q, v in zip(("p10", "p50", "p90"), e["p10_p50_p90_mm"]):
            F.add(f"g3b_{nt}_inner_skull_depth_{q}", mm(v), v,
                  f"{CHECKS} :: inner_skull_scalp_depth_upper_head['{key}'].p10_p50_p90_mm ({q}; inner skull below the scalp over "
                  "the upper head, mm)")
    p50 = {nt: isd[key]["p10_p50_p90_mm"][1] for key, nt in names.items()}
    key = "inner_skull_scalp_depth_upper_head['<entry>'].p10_p50_p90_mm[1]"
    F.range("g3b_children_watershed_inner_skull_depth_p50", {c: p50[f"{TOK[c]}_watershed"] for c in SCHOOL}, mm, CHECKS, key,
            "the children's watershed inner skulls")
    F.range("g3b_children_modelled_inner_skull_depth_p50", {c: p50[f"{TOK[c]}_modelled"] for c in SCHOOL}, mm, CHECKS, key,
            "the children's modelled inner skulls")
    F.range("g3b_templates_inner_skull_depth_p50", {k: p50[k] for k in TEMPLATES}, mm, CHECKS, key, "the three templates")
    thin = ch["children_thinnest_layers"]
    layers = (("white_to_scalp_min_mm", "white_to_scalp_closest"), ("white_to_inner_skull_min_mm", "white_to_inner_skull_closest"),
              ("inner_to_outer_skull_min_mm", "inner_to_outer_skull_thinnest"),
              ("outer_skull_to_head_min_mm", "outer_skull_to_head_thinnest"))
    for c in SCHOOL:
        for q, nt in layers:
            v, src = thin[c][q], f"{CHECKS} :: children_thinnest_layers['{c}'].{q} (mm)"
            F.add(f"g3b_{TOK[c]}_{nt}", mm(v), v, src)
            F.add(f"g3b_{TOK[c]}_{nt}_2dp", mm2(v), v, src)
        for q in ("white_vertices_outside_inner_skull", "n_white_vertices"):
            F.add(f"g3b_{TOK[c]}_{q}", count(thin[c][q]), thin[c][q], f"{CHECKS} :: children_thinnest_layers['{c}'].{q}")
    for q, nt in layers:
        for fmt, sfx in ((mm, ""), (mm2, "_2dp")):
            F.range(f"g3b_children_{nt}{sfx}", {c: thin[c][q] for c in SCHOOL}, fmt, CHECKS, f"children_thinnest_layers['<child>'].{q}",
                    GROUP_DESC["children"])
    fid = ch["fiducial_transfer_on_templates"]
    tnames = dict(zip(("ANTS2-0Years3T", "ANTS18-0Months3T", "ANTS12-0Months3T"), TEMPLATES))
    for key, k in tnames.items():
        e, path = fid[key], f"{CHECKS} :: fiducial_transfer_on_templates['{key}']"
        for f_, v in e["distance_mm"].items():
            F.add(f"g3b_{k}_fid_transfer_{f_}", mm(v), v, f"{path}.distance_mm.{f_} (transferred vs the template's own fiducial, mm)")
            F.add(f"g3b_{k}_fid_scalpfit_{f_}", mm(e["scalp_fit"]["distance_mm"][f_]), e["scalp_fit"]["distance_mm"][f_],
                  f"{path}.scalp_fit.distance_mm.{f_} (the same transfer fitted between the scalps, mm)")
        F.add(f"g3b_{k}_fid_transfer_rotation", fixed(e["head_frame_rotation_deg"], 1), e["head_frame_rotation_deg"],
              f"{path}.head_frame_rotation_deg (deg)")
        F.add(f"g3b_{k}_fid_scalpfit_rotation", fixed(e["scalp_fit"]["head_frame_rotation_deg"], 1),
              e["scalp_fit"]["head_frame_rotation_deg"], f"{path}.scalp_fit.head_frame_rotation_deg (deg)")
        F.add(f"g3b_{k}_fid_cortex_fit_median", mm(e["cortex_fit_median_mm"]), e["cortex_fit_median_mm"], f"{path}.cortex_fit_median_mm")
    path = "fiducial_transfer_on_templates['<template>']"
    F.range("g3b_templates_fid_transfer", {f"{k} {f_}": v for key, k in tnames.items() for f_, v in fid[key]["distance_mm"].items()},
            mm, CHECKS, f"{path}.distance_mm (lpa, nasion, rpa)", "the three templates and three fiducials")
    F.range("g3b_templates_fid_transfer_rotation", {k: fid[key]["head_frame_rotation_deg"] for key, k in tnames.items()},
            lambda v: fixed(v, 1), CHECKS, f"{path}.head_frame_rotation_deg", "the three templates")
    F.range("g3b_templates_fid_scalpfit_nasion", {k: fid[key]["scalp_fit"]["distance_mm"]["nasion"] for key, k in tnames.items()},
            mm, CHECKS, f"{path}.scalp_fit.distance_mm.nasion", "the three templates")
    mni = ch["mni_fiducials_on_adult"]
    for f_, v in mni["distance_mm"].items():
        F.add(f"g3b_adult_mni_fid_{f_}", mm(v), v, f"{CHECKS} :: mni_fiducials_on_adult.distance_mm.{f_} (mm)")
    F.range("g3b_adult_mni_fid", mni["distance_mm"], mm, CHECKS, "mni_fiducials_on_adult.distance_mm", "the three fiducials")
    F.add("g3b_adult_mni_fid_rotation", fixed(mni["head_frame_rotation_deg"], 1), mni["head_frame_rotation_deg"],
          f"{CHECKS} :: mni_fiducials_on_adult.head_frame_rotation_deg (deg)")
    tal = ch["talairach_files"]
    own = tal["childA/sub-Z209/mri/transforms/talairach.xfm"]
    F.add("g3b_child_a_own_talairach_offset", mm(own["centroid_offset_from_adult_mm"]), own["centroid_offset_from_adult_mm"],
          f"{CHECKS} :: talairach_files['childA/sub-Z209/mri/transforms/talairach.xfm'].centroid_offset_from_adult_mm (child A's own file)")
    other = {k: v["centroid_offset_from_adult_mm"] for k, v in tal.items() if not v["own_file"]}
    F.range("g3b_children_folder_talairach_offset", other, mm, CHECKS, "talairach_files['<file>'].centroid_offset_from_adult_mm",
            "the files stored in the children's own folders (other subjects')")
    F.add("g3b_checks_commit", ch["provenance"]["commit"], ch["provenance"]["commit"], f"{CHECKS} :: provenance.commit")

    bv = json.loads((root / BEMVAL).read_text())
    F.add("g3b_bemval_failed_depth", mm(bv["failed_depth_mm"]), bv["failed_depth_mm"],
          f"{BEMVAL} :: failed_depth_mm (the adult's inner skull degraded to this depth below the scalp, mm)")
    F.add("g3b_bemval_degraded_vertices", count(bv["degraded_vertices"]), bv["degraded_vertices"], f"{BEMVAL} :: degraded_vertices")
    for q in ("targets_kept", "grid_kept"):
        F.add(f"g3b_bemval_{q}", count(bv[q]), bv[q], f"{BEMVAL} :: {q}")
    F.add("g3b_bemval_runtime_s", count(bv["runtime_s"]), bv["runtime_s"], f"{BEMVAL} :: runtime_s (s)")
    for m, e in bv["models"].items():
        mt = m.replace("modelled_", "model").replace(".", "p")
        F.add(f"g3b_bemval_{mt}_inner_skull_error", mm(e["inner_skull_abs_error_median_mm"]), e["inner_skull_abs_error_median_mm"],
              f"{BEMVAL} :: models['{m}'].inner_skull_abs_error_median_mm (median distance from the segmented inner skull, mm)")
    seg = bv["variants"]["segmented"]
    changes = {}
    for var, e in bv["variants"].items():
        vt = var.replace("modelled_", "model").replace(".", "p")
        for key in ("opm_dense/combined/intrinsic+brain", "opm_dense/combined/projected", "opm_matched/combined/intrinsic+brain",
                    "opm_matched/combined/projected"):
            o, _, k = key.split("/")
            stem = f"g3b_bemval_{vt}_{OPM[o]}_vs_combined_{COND[k]}"
            v, src = e[key]["ratio"], f"{BEMVAL} :: variants['{var}']['{key}'].ratio (G2 headline ratio on the adult, this inner skull)"
            F.ratio(f"{stem}_ratio", v, src)
            F.add(f"{stem}_ratio_4dp", fixed(v, 4), v, src)
            F.interval(f"{stem}_ratio_ci", ratio, e[key]["ci95"], f"{BEMVAL} :: variants['{var}']['{key}'].ci95")
            if var != "segmented":
                changes[f"{var} {key}"] = abs(v - seg[key]["ratio"])
                for i, row in enumerate(e[key]["by_depth"]):
                    changes[f"{var} {key} {row['lo']:g}-{row['hi']:g} mm"] = abs(row["ratio"] - seg[key]["by_depth"][i]["ratio"])
    head = {k: v for k, v in changes.items() if not k.endswith(" mm")}
    k_max = max(head, key=head.get)
    F.add("g3b_bemval_max_ratio_change", fixed(head[k_max], 3), head[k_max],
          f"{BEMVAL} :: derived: max over the modelled 6-10 mm variants and the four array/condition ratios of |ratio - segmented "
          f"ratio| (at {k_max})")
    k_max = max(changes, key=changes.get)
    F.add("g3b_bemval_max_ratio_change_by_depth", fixed(changes[k_max], 3), changes[k_max],
          f"{BEMVAL} :: derived: as g3b_bemval_max_ratio_change, also over the four depth bands (at {k_max})")
    F.add("g3b_bemval_commit", bv["provenance"]["commit"], bv["provenance"]["commit"], f"{BEMVAL} :: provenance.commit")

    prep = json.loads((root / PREP).read_text())
    for i, e in enumerate(prep["children"]):
        t, path = TOK[e["key"]], f"{PREP} :: children[{i}]"
        sm, fd, hc = e["skull_model"], e["fiducials"], e["head_conform"]
        F.add(f"g3b_{t}_skull_cortex_clearance_closest", mm2(sm["cortex_clearance_min_mm"]), sm["cortex_clearance_min_mm"],
              f"{path}.skull_model.cortex_clearance_min_mm (modelled inner skull to the nearest white-surface vertex, mm)")
        F.add(f"g3b_{t}_skull_thickness_p50", mm(sm["skull_thickness_mm"][1]), sm["skull_thickness_mm"][1],
              f"{path}.skull_model.skull_thickness_mm[1] (median, mm)")
        F.add(f"g3b_{t}_skull_moved_share", pct(100 * sm["moved_share"]), 100 * sm["moved_share"],
              f"{path}.skull_model.moved_share x 100 (inner-skull vertices moved, %)")
        F.add(f"g3b_{t}_skull_moved_median", mm(sm["moved_median_mm"]), sm["moved_median_mm"], f"{path}.skull_model.moved_median_mm")
        F.add(f"g3b_{t}_fid_scale", ratio3(fd["scale"]), fd["scale"], f"{path}.fiducials.scale (similarity fit adult -> child cortex)")
        F.add(f"g3b_{t}_fid_rotation", fixed(fd["rotation_deg"], 1), fd["rotation_deg"], f"{path}.fiducials.rotation_deg (deg)")
        F.add(f"g3b_{t}_fid_cortex_fit_median", mm(fd["cortex_fit_median_mm"]), fd["cortex_fit_median_mm"],
              f"{path}.fiducials.cortex_fit_median_mm")
        F.add(f"g3b_{t}_head_conform_median", mm(hc["moved_median_mm"]), hc["moved_median_mm"], f"{path}.head_conform.moved_median_mm")
        F.add(f"g3b_{t}_cortex_usable_pct", pct(100 * e["cortex_usable_share"]), 100 * e["cortex_usable_share"],
              f"{path}.cortex_usable_share x 100 (%)")
    F.add("g3b_prep_commit", prep["provenance"]["commit"], prep["provenance"]["commit"], f"{PREP} :: provenance.commit")


# ----------------------------------------------------------------------------------------------
# G3B: configuration and provenance
def config_facts(F: Facts, s: dict, root: Path) -> None:
    cfg = s["config"]
    g2 = tomllib.loads((root / G2CFG).read_text())
    path = f"{G3B} :: config"
    F.add("g3b_cfg_primary_placement", cfg["placement"]["primary"], cfg["placement"]["primary"], f"{path}.placement.primary")
    for q, nt in (("clearance_mm", "clearance"), ("dewar_spacing_mm", "dewar_spacing"), ("translation_mm", "translation"),
                  ("pitch_deg", "pitch"), ("roll_deg", "roll"), ("yaw_deg", "yaw")):
        F.add(f"g3b_cfg_{nt}", const(cfg["placement"][q]), cfg["placement"][q], f"{path}.placement.{q}")
    src = cfg["sources"]
    F.add("g3b_cfg_patch_radii", ", ".join(const(r) for r in src["patch_radii_mm"]), src["patch_radii_mm"],
          f"{path}.sources.patch_radii_mm")
    for q, nt in (("patch_total_nAm", "patch_total"), ("fixed_density_nAm_per_mm2", "fixed_density"),
                  ("n_patch_centres", "n_patch_centres")):
        F.add(f"g3b_cfg_{nt}", const(src[q]), src[q], f"{path}.sources.{q}")
    F.add("g3b_cfg_patch20_stride", const(src["patch_centre_stride"][-1]), src["patch_centre_stride"][-1],
          f"{path}.sources.patch_centre_stride[-1] (20-mm patches on every k-th centre)")
    u = cfg["usefulness"]
    F.add("g3b_cfg_useful_threshold", const(u["detectability_threshold"]), u["detectability_threshold"],
          f"{path}.usefulness.detectability_threshold")
    F.add("g3b_cfg_useful_moments", ", ".join(const(q) for q in u["reference_moments_nAm"]), u["reference_moments_nAm"],
          f"{path}.usefulness.reference_moments_nAm (nAm)")
    F.add("g3b_cfg_useful_primary_moment", const(u["primary_reference_nAm"]), u["primary_reference_nAm"],
          f"{path}.usefulness.primary_reference_nAm (nAm)")
    F.add("g3b_cfg_bg_factors", ", ".join(const(f) for f in cfg["background"]["variance_factors"]),
          cfg["background"]["variance_factors"], f"{path}.background.variance_factors")
    vf = cfg["background"]["variance_factors"]
    F.add("g3b_cfg_bg_factor_low", const(min(vf)), min(vf), f"{path}.background.variance_factors (smallest factor)")
    F.add("g3b_cfg_bg_factor_high", const(max(vf)), max(vf), f"{path}.background.variance_factors (largest factor)")
    st = cfg["strata"]
    F.add("g3b_cfg_depth_edges", ", ".join(const(x) for x in st["depth_edges_mm"]), st["depth_edges_mm"], f"{path}.strata.depth_edges_mm")
    edges = st["depth_edges_mm"]
    for lo, hi in zip(edges[:-1], edges[1:]):  # depth strata below each head's scalp, as labels ('20–25')
        F.add(f"g3b_depth_stratum_{const(lo)}_{const(hi)}_mm", f"{const(lo)}–{const(hi)}", [lo, hi],
              f"{path}.strata.depth_edges_mm (stratum {const(lo)}-{const(hi)} mm below the scalp)")
    if edges[1] != 10 or edges[5] != 30:
        raise ValueError("depth strata changed: the span facts below assume edges 10 and 30 mm")
    F.add("g3b_depth_span_10_30_mm", f"{const(edges[1])}–{const(edges[5])}", [edges[1], edges[5]],
          f"{path}.strata.depth_edges_mm (the strata from {const(edges[1])}-{const(edges[2])} to {const(edges[4])}-{const(edges[5])} mm)")
    F.add("g3b_depth_span_10_20_mm", "10–20", [10, 20],
          f"{G3B} :: template_depth_checks['<anat>'].area_share_10_20mm (the area share 10-20 mm below the scalp)")
    F.add("g3b_shallow_target_limit_mm", const(edges[1]), edges[1],
          f"{path}.strata.depth_edges_mm[1] (the shallowest stratum, {const(edges[0])}-{const(edges[1])} mm)")
    F.add("g3b_cfg_orientation_edges", ", ".join(const(x) for x in st["orientation_edges_deg"]), st["orientation_edges_deg"],
          f"{path}.strata.orientation_edges_deg (deg; 90.1 closes the last stratum at 90)")
    F.add("g3b_cfg_min_n", const(st["min_n"]), st["min_n"], f"{path}.strata.min_n (fewer targets: sparse)")
    F.add("g3b_cfg_n_boot", count(st["n_boot"]), st["n_boot"], f"{path}.strata.n_boot (parcel bootstrap resamples)")
    an = cfg["anatomy"]
    F.add("g3b_cfg_bem_conductivity", ", ".join(const(x) for x in an["bem_conductivity"]), an["bem_conductivity"],
          f"{path}.anatomy.bem_conductivity (S/m: scalp, skull, brain)")
    F.add("g3b_cfg_school_scale", ratio3(an["school_age_scale"]), an["school_age_scale"],
          f"{path}.anatomy.school_age_scale (85/95: Jas et al. Table 1 8-year / adult head radius)")
    F.add("g3b_cfg_child_skull_depth", const(an["child_skull_depth_mm"]), an["child_skull_depth_mm"],
          f"{path}.anatomy.child_skull_depth_mm (A-BEM-CHILD)")
    F.add("g3b_cfg_child_skull_min_cortex", const(an["child_skull_min_cortex_mm"]), an["child_skull_min_cortex_mm"],
          f"{path}.anatomy.child_skull_min_cortex_mm (A-BEM-CHILD)")
    F.add("g3b_cfg_conditions_headline", ", ".join(cfg["conditions"]["headline"]), cfg["conditions"]["headline"],
          f"{path}.conditions.headline")
    gs = f"{G2CFG} ::"
    F.add("g3b_opm_noise_primary", const(g2["sensors"]["opm_asd_primary_fT_per_rtHz"]), g2["sensors"]["opm_asd_primary_fT_per_rtHz"],
          f"{gs} sensors.opm_asd_primary_fT_per_rtHz (fT/sqrt(Hz); G3B reads it from the G2 configuration)")
    F.add("g3b_opm_noise_sweep", ", ".join(const(x) for x in g2["sensors"]["opm_asd_fT_per_rtHz"]), g2["sensors"]["opm_asd_fT_per_rtHz"],
          f"{gs} sensors.opm_asd_fT_per_rtHz (fT/sqrt(Hz))")
    F.add("g3b_squid_noise_mag", const(g2["sensors"]["squid_asd"]["mag_fT_per_rtHz"]), g2["sensors"]["squid_asd"]["mag_fT_per_rtHz"],
          f"{gs} sensors.squid_asd.mag_fT_per_rtHz (fT/sqrt(Hz))")
    F.add("g3b_squid_noise_grad", const(g2["sensors"]["squid_asd"]["grad_fT_per_cm_rtHz"]),
          g2["sensors"]["squid_asd"]["grad_fT_per_cm_rtHz"], f"{gs} sensors.squid_asd.grad_fT_per_cm_rtHz (fT/cm/sqrt(Hz))")
    F.add("g3b_band_low", const(g2["band"]["l_freq_hz"]), g2["band"]["l_freq_hz"], f"{gs} band.l_freq_hz (Hz)")
    F.add("g3b_band_high", const(g2["band"]["h_freq_hz"]), g2["band"]["h_freq_hz"], f"{gs} band.h_freq_hz (Hz)")
    F.add("g3b_focal_moment", const(g2["sources"]["focal_nAm"]), g2["sources"]["focal_nAm"], f"{gs} sources.focal_nAm (nAm)")
    F.add("g3b_background_grid_spacing", const(g2["background"]["grid_spacing_mm"]), g2["background"]["grid_spacing_mm"],
          f"{gs} background.grid_spacing_mm (mm)")
    F.add("g3b_enbw", fixed(s["enbw_hz"], 1), s["enbw_hz"], f"{G3B} :: enbw_hz (equivalent noise bandwidth of the analysis filter, Hz)")
    meta = s["arrays"]["adult"]["opm_dense"]
    F.add("g3b_opm_standoff_nominal", const(meta["standoff_mm"]), meta["standoff_mm"], f"{G3B} :: arrays['adult']['opm_dense'].standoff_mm")
    text = (root / METHODS).read_text()
    m = re.search(r"(\d+)-mm cell, (\d+)-mm\s+standoff, (\d+)-mm packing", text)
    if m:
        quote = " ".join(m.group(0).split())
        F.add("g3b_opm_cell", m.group(1), int(m.group(1)), f"{METHODS} :: G3B 'Head-adaptive OPM arrays', quoted ('{quote}', mm)")
        F.add("g3b_opm_packing", m.group(3), int(m.group(3)), f"{METHODS} :: G3B 'Head-adaptive OPM arrays', quoted ('{quote}', mm)")
    m = re.search(r"every Delta fell by ([0-9.]+)-([0-9.]+) dB", text)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        src = f"{METHODS} :: G3B results paragraph, quoted ('{m.group(0)}'; v3 -> v4 change, the v3 results are not committed)"
        F.add("g3b_v3_v4_delta_drop_min", fixed(lo), lo, src)
        F.add("g3b_v3_v4_delta_drop_max", fixed(hi), hi, src)
        F.add("g3b_v3_v4_delta_drop_range", span(fixed, lo, hi), [lo, hi], src)
    F.add("g3b_commit", s["provenance"]["commit"], s["provenance"]["commit"], f"{G3B} :: provenance.commit (computation)")
    F.add("g3b_replot_commit", s["replotted_at_commit"], s["replotted_at_commit"], f"{G3B} :: replotted_at_commit (summaries redrawn)")
    F.add("g3b_mne_version", s["provenance"]["mne_version"], s["provenance"]["mne_version"], f"{G3B} :: provenance.mne_version")


# ----------------------------------------------------------------------------------------------
# regions by head (regh_): dense OPM vs Neuromag (i+b) per lobe, parcel and parcel group. Parcels and parcel groups from
# the per-target CSVs; lobes from the summary's placement_D by_lobe (unrounded inputs; the CSV recomputation, recorded
# in each source, agrees within 0.001 dB); counterfactual_x-centred has no per-target column, so lobes only.
def read_targets(path: Path) -> dict:
    """Columns of a g3b_targets_<anatomy>.csv (comment line skipped): text columns as arrays of str, the rest as floats."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(line for line in fh if not line.startswith("#")))
    head, data = rows[0], rows[1:]
    cols = {}
    for j, name in enumerate(head):
        vals = [r[j] for r in data]
        cols[name] = np.array(vals) if name in ("region", "lobe") else np.array(vals, dtype=float)
    return cols


def region_facts(F: Facts, s: dict, root: Path) -> None:
    plc = s["placement_D"]
    series = {}  # (placement, ref, region) -> {anatomy: D}
    for a in ANATS:
        t, rel = TOK[a], f"results/g3b/g3b_targets_{a}.csv"
        cols = read_targets(root / rel)
        cortical = ~np.char.endswith(cols["region"].astype(str), "unknown")
        parcel = np.array([r.split(".", 1)[-1] for r in cols["region"]])
        w = cols["area_mm2"]
        opm = cols["detect_opm_dense_opm_intrinsic+brain"]
        for reg, (kind, what) in REGIONS.items():
            m = cortical & ((cols["lobe"] == what) if kind == "lobe" else np.isin(parcel, what))
            desc = f"lobe '{what}'" if kind == "lobe" else "parcels " + "+".join(what) + " (both hemispheres)"
            sel = f"targets in the {desc}, medial wall ('unknown') excluded"
            F.add(f"regh_{t}_{reg}_n", count(m.sum()), int(m.sum()), f"{rel} :: derived: number of {sel}")
            F.add(f"regh_{t}_{reg}_area", fixed(w[m].sum() / 100, 1), w[m].sum() / 100,
                  f"{rel} :: derived: sum of area_mm2 / 100 over {sel} (cm^2)")
            dep = weighted_median(cols["depth_mm"][m], w[m])
            F.add(f"regh_{t}_{reg}_depth", mm(dep), dep,
                  f"{rel} :: derived: area-weighted median (weights area_mm2) of depth_mm over {sel} (mm below the scalp)")
            dep = float(np.median(cols["depth_mm"][m]))
            F.add(f"regh_{t}_{reg}_depth_unweighted", mm(dep), dep,
                  f"{rel} :: derived: median (unweighted, numpy) of depth_mm over {sel} (mm below the scalp)")
            for pl, col in REGH_PLACE.items():
                for r in REFS:
                    sq = cols[f"detect_{col}_{r}_intrinsic+brain"]
                    v = weighted_median(20 * np.log10(opm[m] / sq[m]), w[m])
                    how = (f"area-weighted median (weights area_mm2: sort, cumulative weight, first value reaching half the total) of "
                           f"20 log10(detect_opm_dense_opm_intrinsic+brain / detect_{col}_{r}_intrinsic+brain) over {sel}; "
                           f"n = {int(m.sum())}")
                    if kind == "lobe":  # the summary holds the same median from unrounded values (the CSV rounds d to 5 decimals)
                        key = f"{a}/{'top' if pl == 'top' else 'counterfactual'}/{r}/intrinsic+brain"
                        exact = plc[key]["by_lobe"][what]
                        series.setdefault((pl, r, reg), {})[a] = exact
                        F.add(f"regh_{t}_{pl}_vs_{r}_{reg}", signed(exact), exact,
                              f"{G3B} :: placement_D['{key}'].by_lobe['{what}'] (cross-check from {rel}: {how} gives {v:+.4f} dB)")
                        continue
                    series.setdefault((pl, r, reg), {})[a] = v
                    F.add(f"regh_{t}_{pl}_vs_{r}_{reg}", signed(v), v, f"{rel} :: derived: {how}")
        for r in REFS:
            for lb, v in plc[f"{a}/counterfactual_x-centred/{r}/intrinsic+brain"]["by_lobe"].items():
                series.setdefault(("cfx", r, lb), {})[a] = v
                F.add(f"regh_{t}_cfx_vs_{r}_{lb}", signed(v), v,
                      f"{G3B} :: placement_D['{a}/counterfactual_x-centred/{r}/intrinsic+brain'].by_lobe['{lb}'] (the targets CSV has no "
                      "counterfactual_x-centred column, so parcels and parcel groups are not available for it)")
    for (pl, r, reg), vals in series.items():
        if r == "grad":
            continue
        if reg in DK_LOBES:
            key = f"placement_D['<anat>/{PLACE_KEY[pl]}/{r}/intrinsic+brain'].by_lobe['{reg}']"
            file = G3B
        else:
            key = f"area-weighted median D of regh_<anat>_{pl}_vs_{r}_{reg}"
            file = ", ".join(f"results/g3b/g3b_targets_{a}.csv" for a in ANATS)
        F.group_ranges(f"regh_{{g}}_{pl}_vs_{r}_{reg}", vals, signed, file, key)
    for reg in REGIONS:
        F.range(f"regh_all_{reg}_n", {a: F[f"regh_{TOK[a]}_{reg}_n"]["raw"] for a in ANATS}, count,
                ", ".join(f"results/g3b/g3b_targets_{a}.csv" for a in ANATS), f"the number of targets of regh_<anat>_{reg}_n",
                GROUP_DESC["all"])
    for r in REFS:
        F.range(f"regh_children_cfx_vs_{r}_lobes", {f"{c} {lb}": series[("cfx", r, lb)][c] for c in SCHOOL for lb in DK_LOBES},
                signed, G3B, f"placement_D['<anat>/counterfactual_x-centred/{r}/intrinsic+brain'].by_lobe",
                "children A-C and the six lobes")
    # the cells of the regions-by-head table (results/report/Figure_R5_regions_heads.png): the six lobes and four parcel rows
    # at top contact; the six lobes in the scaled, laterally centred counterfactual helmet
    rows = list(DK_LOBES) + ["precentral", "superiortemporal", "parahippocampal", "mesial_temporal"]
    top = {reg: series[("top", "combined", reg)] for reg in rows}
    files = ", ".join(f"results/g3b/g3b_targets_{a}.csv" for a in ANATS)
    key = "regh_<anat>_top_vs_combined_<row> (lobes: placement_D by_lobe; parcels and the mesial temporal group: area-weighted medians)"
    F.range("regh_top_vs_combined_cells", {f"{a} {reg}": top[reg][a] for reg in rows for a in ANATS}, signed, f"{G3B}, {files}",
            key, "all nine heads and the ten rows (six lobes, precentral, superior temporal, parahippocampal, mesial temporal)")
    above = sum(top[reg][a] > top[reg]["adult"] for reg in rows for a in SMALLER)
    src = (f"{G3B}, {files} :: derived: number of cells (smaller head, row) of {key} that exceed the adult's value in the same row, "
           "over the eight smaller heads and the ten rows")
    F.add("regh_top_vs_combined_smaller_cells_above_adult_count", count(above), above, src)
    F.add("regh_top_vs_combined_smaller_cells_n", count(len(rows) * len(SMALLER)), len(rows) * len(SMALLER),
          f"{G3B}, {files} :: derived: eight smaller heads x ten rows")
    cfx = {lb: series[("cfx", "combined", lb)] for lb in DK_LOBES}
    at_or_below = sum(cfx[lb][a] <= cfx[lb]["adult"] for lb in DK_LOBES for a in SMALLER)
    key = "placement_D['<anat>/counterfactual_x-centred/combined/intrinsic+brain'].by_lobe"
    F.add("regh_cfx_vs_combined_smaller_lobe_cells_at_or_below_adult_count", count(at_or_below), at_or_below,
          f"{G3B} :: derived: number of cells (smaller head, lobe) of {key} at or below the adult's value for the same lobe, "
          "over the eight smaller heads and the six lobes")
    F.add("regh_cfx_vs_combined_smaller_lobe_cells_n", count(len(DK_LOBES) * len(SMALLER)), len(DK_LOBES) * len(SMALLER),
          f"{G3B} :: derived: eight smaller heads x six lobes of {key}")


# ----------------------------------------------------------------------------------------------
def extra_facts(F: Facts, s: dict, s3a: dict) -> None:
    """Anatomy labels as numbers (template ages, sphere ages), the feasible placement families, the NumPy version, Delta
    and D ranges for the scaled adults and the templates separately, and the lobe cells of the scaled, laterally centred
    helmet counted against the adult at top contact."""
    an, cfg = s["anatomies"], s["config"]
    for a in TEMPLATES:
        m = re.search(r"ANTS(\d+)-0(Years|Months)3T", an[a]["description"])
        if not m:
            raise ValueError(f"{G3B}: template name not found in anatomies['{a}'].description")
        months = int(m.group(1)) * (12 if m.group(2) == "Years" else 1)
        F.add(f"g3b_{a}_age_months", str(months), months,
              f"{G3B} :: anatomies['{a}'].description and config.anatomy.templates (template {m.group(0)}: age in months)")
        if m.group(2) == "Years":
            F.add(f"g3b_{a}_age_years", m.group(1), int(m.group(1)),
                  f"{G3B} :: anatomies['{a}'].description (template {m.group(0)}: age in years; the 2-year size control is "
                  "the adult scaled to this template's head circumference, anatomies['size2yr'].scale_note)")
    for head, h in HEAD3A.items():
        m = re.search(r"\((\d+)-yr\)", head)
        if m:
            F.add(f"g3a_{h}_age_years", m.group(1), int(m.group(1)), f"{G3A} :: size_following['{head}'] (the head's age label, years)")
    inf = s["infeasible_placements"]
    for a in ANATS:
        n = len([p for p in FAMILY if p not in inf.get(a, [])])
        F.add(f"g3b_{TOK[a]}_n_family_feasible", count(n), n,
              f"{G3B} :: derived: the {len(FAMILY)} family placements minus infeasible_placements['{a}'] ({', '.join(inf.get(a, [])) or 'none'})")
    F.add("g3b_numpy_version", s["provenance"]["numpy_version"], s["provenance"]["numpy_version"], f"{G3B} :: provenance.numpy_version")
    comp, dm = s["comparisons"], s["D_median_dB"]
    for k, kt in (("intrinsic+brain", "ib"), ("projected", "proj")):
        key = f"opm_dense/combined/{k}/detect"
        F.group_ranges(f"g3b_{{g}}_dense_vs_combined_{kt}_delta", {c: comp[f"{c}/{key}"]["delta"]["median"] for c in SMALLER},
                       signed, G3B, f"comparisons['<anat>/{key}'].delta.median", groups=("scaled", "templates"))
        F.group_ranges(f"g3b_{{g}}_dense_vs_combined_{kt}_d", {a: dm[f"{a}/{key}"] for a in SMALLER}, signed, G3B,
                       f"D_median_dB['<anat>/{key}']", groups=("scaled", "templates"))
    plc = s["placement_D"]
    # the counterfactual helmet with the room field projected out: no Delta estimator is stored for it
    # (delta_other_placements holds intrinsic + brain only), so the difference of the pooled medians
    for k, kt in (("projected", "proj"),):
        key = f"placement_D['<anat>/counterfactual_x-centred/combined/{k}'].median - placement_D['adult/counterfactual_x-centred/combined/{k}'].median"
        ref = plc[f"adult/counterfactual_x-centred/combined/{k}"]["median"]
        F.group_ranges(f"g3b_{{g}}_cfx_dense_vs_combined_{kt}_d_minus_adult_cfx",
                       {a: plc[f"{a}/counterfactual_x-centred/combined/{k}"]["median"] - ref for a in SMALLER}, signed, G3B, key)
    top = plc["adult/top/combined/intrinsic+brain"]["by_lobe"]
    val = lambda x: x["median"] if isinstance(x, dict) else x  # noqa: E731
    key = "placement_D['<anat>/counterfactual_x-centred/combined/intrinsic+brain'].by_lobe"
    for g, labs in (("smaller", SMALLER), ("templates_scaled", SCALED + TEMPLATES), ("children", SCHOOL)):
        cells = [(a, lb) for a in labs for lb in DK_LOBES]
        n = sum(val(plc[f"{a}/counterfactual_x-centred/combined/intrinsic+brain"]["by_lobe"][lb]) <= val(top[lb]) for a, lb in cells)
        F.add(f"regh_cfx_vs_combined_{g}_lobe_cells_at_or_below_adult_top_count", count(n), n,
              f"{G3B} :: derived: number of cells ({GROUP_DESC[g]} x six lobes) of {key} at or below the adult's value for the "
              "same lobe at top contact (placement_D['adult/top/combined/intrinsic+brain'].by_lobe)")
        if g != "smaller":
            F.add(f"regh_cfx_vs_combined_{g}_lobe_cells_n", count(len(cells)), len(cells),
                  f"{G3B} :: derived: {GROUP_DESC[g]} x six lobes of {key}")


# ----------------------------------------------------------------------------------------------
def facts(root: Path = ROOT) -> dict:
    """Every G3A, G3B and regions-by-head fact, name -> {"value", "raw", "source"}."""
    root = Path(root)
    F = Facts()
    s3a = json.loads((root / G3A).read_text())
    g3a_facts(F, s3a)
    s = json.loads((root / G3B).read_text())
    anatomy_facts(F, s)
    placement_facts(F, s)
    d_median_facts(F, s)
    comparison_facts(F, s)
    other_placement_facts(F, s)
    placement_d_facts(F, s)
    mechanism_facts(F, s)
    depth_check_facts(F, s)
    channel_count_facts(F, s)
    sensitivity_facts(F, s)
    patch_facts(F, s)
    usefulness_facts(F, s)
    school_check_facts(F, root)
    config_facts(F, s, root)
    region_facts(F, s, root)
    extra_facts(F, s, s3a)
    return dict(F)


if __name__ == "__main__":
    out = facts(ROOT)
    for k, f in out.items():
        print(f"{k}\t{f['value']}\t{f['source']}")
    print(f"{len(out)} facts", file=sys.stderr)
