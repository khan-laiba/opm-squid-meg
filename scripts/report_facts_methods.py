#!/usr/bin/env python3
"""Report facts for the Methods details asked for in referee round 1 (prefix meth_), merged by scripts/report_facts.py.

The referees asked for: the cortical background's moment variance per source and its grid, the oracle covariance, the
OPM sensitive-axis construction and its effect, the external-field projection, the Neuromag system modelled, the
parameters of the time-domain (spike) model, the two dSPM implementations, the number of bootstrap resamples of every
interval family, the random seeds of the spike runs, the 12-month template's thin skull and the single-layer check, what
18 (or 36) locations per band allow a location-level test to show, and the provenance of the declared parameters.
report/methods_additions.md uses these facts.

facts(root) -> {name: {"value": text as printed, "raw": unrounded number(s) or text, "source": "<file> :: <key path>" or
"<file> :: derived: <how>"}}. Sources are the committed result files (with key paths), the configurations (TOML key
paths) and the documentation: values set in code are read from the 'Implementation constants' table of docs/methods.md
(section 13; one row per value, by ID, naming the file and line it comes from); register and literature values are
quoted, and every quote is checked against the document when the facts are built (a missing quote raises), so an edited
document cannot leave a stale fact. Derived values are exact arithmetic on stored values: unit conversions, ratios of
stored totals and the exact two-sided sign test implied by the stored location design. Nothing is simulated or
re-estimated.

Formats: counts and resample numbers with thousands separators; seeds as integers without separators; dB signed with 2
decimals; ratios 3 decimals (the head-surface decomposition) or 2; p-values 2 significant digits (scientific notation
below 0.0001); U+2212 for minus; configuration constants as declared; ranges "lo to hi" from unrounded values.

Usage: .venv/bin/python scripts/report_facts_methods.py [--grep TEXT]   (prints name, value, source)
"""
from __future__ import annotations

import argparse
import json
import math
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MINUS = "−"
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
LABELS = ("adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC")
SMALLER = LABELS[1:]
TOK = {a: a.lower() for a in LABELS}

G2 = "results/g2/g2_summary.json"
HSE = "results/g2/head_surface_effect.json"
BAND = "results/g2/g2_band_sensitivity.json"
G3B = "results/g3b/g3b_summary.json"
PED = "results/g4/g4_pediatric_comparison.json"
CFG_CONFIRM = "configs/g4_confirmatory.toml"
CFG_NOISE = "configs/g2_noise_sensitivity.toml"
CFG_G3B = "configs/g3b_pediatric.toml"
CFG_MOTION = "configs/g4_motion.toml"
METHODS = "docs/methods.md"
REGISTER = "docs/provenance_register.md"
JAS = "docs/literature/jas2026.md"
IC_HEADING = "## 13. Implementation constants"
DENSE_IB, DENSE_PROJ = "opm_dense/combined/intrinsic+brain", "opm_dense/combined/projected"
MATCHED_PROJ = "opm_matched/combined/projected"


def g4_summary(label: str) -> str:
    return f"results/g4/g4_{label}_summary.json"


def loc_summary(label: str) -> str:
    return "results/g4/g4_localization_summary.json" if label == "adult" else f"results/g4/g4_localization_{label}_summary.json"


ALL_G4 = ", ".join(g4_summary(lab) for lab in LABELS)


# ------------------------------------------------------------------------------------------------
# formatting
def _fmt(x, nd: int, sep: bool = False, sign: bool = False) -> str:
    s = format(x, f"{',' if sep else ''}.{nd}f")
    if s.lstrip("-").strip("0.,") == "":  # rounds to zero: no sign
        return s.lstrip("-")
    if sign and not s.startswith("-"):
        s = "+" + s
    return s.replace("-", MINUS)


def num(x, nd: int) -> str:
    return _fmt(x, nd)


def db(x) -> str:
    return _fmt(x, 2, sign=True)


def count(n) -> str:
    return f"{int(n):,d}"


def seed(n) -> str:
    return str(int(n))


def sig(x, n: int = 2) -> str:
    """n significant digits in plain decimals (integers keep all their digits)."""
    if x == 0:
        return "0"
    e = int(format(abs(x), f".{n - 1}e").split("e")[1])
    return _fmt(x, max(n - 1 - e, 0), sep=True)


def sci(x, n: int = 2) -> str:
    """'7.6 × 10⁻⁶'."""
    m, e = format(x, f".{n - 1}e").split("e")
    return f"{m} × 10{str(int(e)).translate(SUP)}"


def pval(p) -> str:
    return sci(p) if p < 1e-4 else sig(p, 2)


def const(x) -> str:
    """A configuration constant as declared (no added decimals); a list as 'a, b and c'."""
    if isinstance(x, (list, tuple)):
        items = [const(v) for v in x]
        return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]
    if isinstance(x, float) and x.is_integer():
        x = int(x)
    return (f"{x:,d}" if isinstance(x, int) else format(x, "g")).replace("-", MINUS)


def span(lo, hi, f) -> str:
    a, b = f(lo), f(hi)
    return a if a == b else f"{a} to {b}"


# ------------------------------------------------------------------------------------------------
class Facts(dict):
    def add(self, name: str, value: str, raw, source: str) -> None:
        if name in self:
            raise ValueError(f"fact {name!r} defined twice")
        if not re.fullmatch(r"meth_[a-z0-9_]+", name):
            raise ValueError(f"bad fact name {name!r}")
        if str(value).strip() == "":
            raise ValueError(f"fact {name!r} has an empty value")
        if " :: " not in source:
            raise ValueError(f"fact {name!r}: source lacks ' :: <key path>'")
        self[name] = dict(value=value, raw=raw, source=source)

    def spread(self, base: str, values: list, f, source: str) -> None:
        lo, hi = float(min(values)), float(max(values))
        self.add(f"{base}_min", f(lo), lo, f"{source} (minimum)")
        self.add(f"{base}_max", f(hi), hi, f"{source} (maximum)")
        self.add(f"{base}_range", span(lo, hi, f), [lo, hi], f"{source} (minimum to maximum)")


def _json(root: Path, rel: str) -> dict:
    return json.loads((root / rel).read_text())


def _toml(root: Path, rel: str) -> dict:
    return tomllib.loads((root / rel).read_text())


def _flat(text: str) -> str:
    return " ".join(text.split())


def quoted(root: Path, rel: str, quote: str, *numbers: str) -> str:
    """``quote`` must occur in the document (whitespace normalised) and contain each printed number: a fact quoted from a
    document fails to build when the document changes. Returns the quote."""
    if _flat(quote) not in _flat((root / rel).read_text()):
        raise ValueError(f"{rel}: quote not found: {quote!r}")
    for n in numbers:
        if n not in quote:
            raise ValueError(f"{rel}: {n!r} is not in the quote {quote!r}")
    return quote


def implementation_constants(root: Path) -> dict:
    """Rows of the 'Implementation constants' table (docs/methods.md section 13): ID -> {constant, value, where}."""
    text = (root / METHODS).read_text()
    if IC_HEADING not in text:
        raise ValueError(f"{METHODS}: section {IC_HEADING!r} not found")
    rows = {}
    for line in text[text.index(IC_HEADING) + len(IC_HEADING):].splitlines():
        if line.startswith("## "):
            break
        if line.startswith("| IC-"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) != 4:
                raise ValueError(f"{METHODS}: malformed implementation-constant row {line[:60]!r}")
            if cells[0] in rows:
                raise ValueError(f"{METHODS}: implementation constant {cells[0]} listed twice")
            rows[cells[0]] = dict(constant=cells[1], value=cells[2], where=cells[3])
    if not rows:
        raise ValueError(f"{METHODS}: no implementation constants found")
    return rows


class IC:
    """Numeric values of the implementation-constant rows with their sources."""

    def __init__(self, root: Path):
        self.rows = implementation_constants(root)

    def value(self, ident: str) -> float:
        v = self.rows[ident]["value"]
        try:
            return float(v)
        except ValueError:
            raise ValueError(f"{METHODS}: {ident} has no numeric value ({v!r})") from None

    def src(self, *idents: str) -> str:
        for i in idents:
            if i not in self.rows:
                raise KeyError(f"{METHODS}: implementation constant {i} not found")
        where = "; ".join(f"{i}, set in {self.rows[i]['where']}" for i in idents)
        return f"{METHODS} :: section 13 'Implementation constants', {where}"


# ------------------------------------------------------------------------------------------------
def background_facts(F: Facts, root: Path, ic: IC) -> None:
    """Cortical background: grid, moment variance per source, calibration; the oracle covariance's terms."""
    g2, g3b = _json(root, G2), _json(root, G3B)
    nv = g2["noise_validation"]
    n = g2["n_background_grid"]
    F.add("meth_bg_n_sources", count(n), n, f"{G2} :: n_background_grid")
    grid = g2["config"]["background"]["grid_spacing_mm"]
    if grid != ic.value("IC-BG-GRID-MM"):
        raise ValueError("background grid spacing: configuration and implementation constant differ")
    F.add("meth_bg_grid_mm", const(grid), grid, f"{G2}, {METHODS} :: config.background.grid_spacing_mm; section 13 "
          f"IC-BG-GRID-MM (greedy Poisson-disk rule, set in {ic.rows['IC-BG-GRID-MM']['where']})")
    sd = nv["brain_source_density_nAm_per_sqrt_mm2"]
    F.add("meth_bg_sd_per_sqrt_mm2_nam", num(sd, 3), sd, f"{G2} :: noise_validation.brain_source_density_nAm_per_sqrt_mm2 "
          "(moment SD per square root of cortical area, nAm per sqrt(mm^2); IC-BG-VAR)")
    s2 = nv["brain_scale"] * 1e18 * 1e-6  # (A m)^2 per m^2 -> nAm^2 per mm^2
    F.add("meth_bg_var_per_mm2_nam2", sig(s2, 3), s2, f"{G2} :: derived: noise_validation.brain_scale (moment variance per "
          "unit cortical area, (A m)^2 per m^2) x 1e18 nAm^2/(A m)^2 x 1e-6 m^2/mm^2: nAm^2 per mm^2 (docs/methods.md IC-BG-VAR)")
    ad = g3b["anatomies"]["adult"]
    if ad["n_background_grid"] != n or abs(g3b["brain_scale"] - nv["brain_scale"]) > 1e-12 * nv["brain_scale"]:
        raise ValueError("G2 and G3B adult background grids or brain scales differ")
    area = ad["cortical_area_cm2"]
    F.add("meth_bg_area_cm2", count(round(area)), area, f"{G3B} :: anatomies.adult.cortical_area_cm2 (the adult's usable "
          "white-surface area, which the Voronoi cells of the background sources partition: docs/methods.md IC-BG-AREA; the G2 "
          f"grid, anatomies.adult.n_background_grid = {n:,})")
    mean_area = area * 100.0 / n
    F.add("meth_bg_mean_area_mm2", count(round(mean_area)), mean_area, f"{G3B} :: derived: anatomies.adult.cortical_area_cm2 "
          "x 100 / anatomies.adult.n_background_grid (mean cell area a_i, mm^2)")
    rms = math.sqrt(s2 * mean_area)
    F.add("meth_bg_rms_sd_nam", num(rms, 2), rms, f"{G2}, {G3B} :: derived: sqrt(brain variance per mm^2 x mean cell area): the "
          "root-mean-square moment SD over the background sources, nAm (the moment variance s^2 a_i is linear in the area)")
    grad = nv["measured"]["brain_rms_grad_fT_cm"]
    F.add("meth_bg_calibration_grad_ft_cm", num(grad, 1), grad, f"{G2} :: noise_validation.measured.brain_rms_grad_fT_cm (square "
          "root of the median over good gradiometers of the measured in-band brain-noise variance, task baseline minus empty "
          "room, fT/cm; the level to which the median modelled gradiometer variance is fitted)")
    F.add("meth_enbw_hz", num(g2["enbw_hz"], 1), g2["enbw_hz"], f"{G2} :: enbw_hz (equivalent noise bandwidth of the "
          "composite 1-40 Hz zero-phase filter)")
    b = g2["config"]["band"]
    F.add("meth_band_lo_hz", const(b["l_freq_hz"]), b["l_freq_hz"], f"{G2} :: config.band.l_freq_hz")
    F.add("meth_band_hi_hz", const(b["h_freq_hz"]), b["h_freq_hz"], f"{G2} :: config.band.h_freq_hz")
    F.add("meth_filter_order", const(b["order"]), b["order"], f"{G2} :: config.band.order (Butterworth, zero phase)")
    u, gr = ic.value("IC-ENV-UNIFORM"), ic.value("IC-ENV-GRADIENT")
    F.add("meth_env_uniform_terms", const(u), u, ic.src("IC-ENV-UNIFORM"))
    F.add("meth_env_gradient_terms", const(gr), gr, ic.src("IC-ENV-GRADIENT"))
    rr = g2["retained_rank"]
    if u + gr != rr["squid/intrinsic+brain+env"] - rr["squid/projected"]:
        raise ValueError("external-field terms and the projected rank loss differ")
    F.add("meth_env_terms", const(u + gr), u + gr, ic.src("IC-ENV-UNIFORM", "IC-ENV-GRADIENT") + "; derived: their sum "
          f"(the rank the projection removes: {G2} retained_rank['squid/intrinsic+brain+env'] - ['squid/projected'])")
    r0 = [1e3 * x for x in g2["config"]["environment"]["r0_head_m"]]
    F.add("meth_env_r0_mm", "(" + ", ".join(const(round(x, 6)) for x in r0) + ")", r0,
          f"{G2} :: config.environment.r0_head_m (expansion point of the external field, head frame, mm)")
    tol = ic.value("IC-WHITEN-TOL")
    F.add("meth_whiten_tol", "10" + str(int(round(math.log10(tol)))).translate(SUP), tol, ic.src("IC-WHITEN-TOL"))
    ne = g2["n_estimate_samples"]
    secs = g2["config"]["covariance"]["estimate_seconds"]
    F.add("meth_plugin_seconds", const(secs), secs, f"{G2} :: config.covariance.estimate_seconds")
    F.add("meth_plugin10s_samples", count(ne["T10"]), ne["T10"], f"{G2} :: n_estimate_samples.T10 (2 x 39 Hz x 10 s)")
    F.add("meth_plugin60s_samples", count(ne["T60"]), ne["T60"], f"{G2} :: n_estimate_samples.T60 (2 x 39 Hz x 60 s)")


def projection_facts(F: Facts, root: Path) -> None:
    g2 = _json(root, G2)
    sq = g2["arrays"]["squid"]
    if sq["axes"] != "1 mag + 2 planar grad per site" or sq["channels"] != 3 * sq["sites"]:
        raise ValueError("unexpected Neuromag array description")
    F.add("meth_squid_channels", count(sq["channels"]), sq["channels"], f"{G2} :: arrays.squid.channels")
    F.add("meth_squid_mags", count(sq["sites"]), sq["sites"], f"{G2} :: arrays.squid.sites (one magnetometer per site)")
    F.add("meth_squid_grads", count(2 * sq["sites"]), 2 * sq["sites"], f"{G2} :: derived: 2 x arrays.squid.sites "
          "(arrays.squid.axes '1 mag + 2 planar grad per site')")
    for key, name in (("opm_dense", "dense"), ("opm_matched", "matched")):
        v = g2["arrays"][key]["n_sites"]
        F.add(f"meth_{name}_sites", count(v), v, f"{G2} :: arrays.{key}.n_sites")
    rr = g2["retained_rank"]
    for key, name in (("squid/projected", "squid"), ("squid_mag/projected", "squid_mag"), ("squid_grad/projected", "squid_grad"),
                      ("opm_dense/projected", "dense"), ("opm_matched/projected", "matched")):
        F.add(f"meth_rank_{name}_proj", count(rr[key]), rr[key], f"{G2} :: retained_rank['{key}'] (rank of the projected "
              "oracle covariance; the Neuromag subsets are channel subsets of the jointly projected data)")


def axis_facts(F: Facts, root: Path) -> None:
    """OPM sensitive axes: construction (register A-OPM-AXIS) and the effect of nearly radial axes (head_surface_effect)."""
    q = quoted(root, REGISTER, "normal of the BEM head surface (5,120 triangles, on the MRI scalp: A-BEM-CONFORM; vertex "
               "normals from its own triangles) averaged within 15 mm", "15")
    F.add("meth_opm_axis_radius_mm", "15", 15, f"{REGISTER} :: A-OPM-AXIS ('{q}')")
    q = quoted(root, REGISTER, "the axes deviate by a median 1.2 deg (95th pct 4.9 deg matched, 7.4 deg dense; the largest, up "
               "to 37 deg, at lower occipital sites and next to the pinna, where no plane fits the scalp)", "1.2", "4.9", "7.4", "37")
    for name, v in (("median", "1.2"), ("p95_matched", "4.9"), ("p95_dense", "7.4"), ("max", "37")):
        F.add(f"meth_opm_axis_dev_{name}_deg", v, float(v), f"{REGISTER} :: A-OPM-AXIS, against the plane fitted to the dense "
              f"MRI scalp within 15 mm of the site ('{q}')")
    F.add("meth_opm_axis_dev_pctl", "95th", 95, f"{REGISTER} :: A-OPM-AXIS ('{q}': the percentile of the p95 values)")
    q = quoted(root, REGISTER, "v1-v3 used the sample outer skin's stored vertex normals, which are not those of its triangles "
               "(within 5 deg of radial from the surface centroid): their axes deviated from that plane by a median 8.4-9.0 deg",
               "5 deg", "8.4-9.0")
    F.add("meth_opm_axis_radial_deg", "5", 5, f"{REGISTER} :: A-OPM-AXIS ('{q}')")
    F.add("meth_opm_axis_radial_dev_range_deg", "8.4 to 9.0", [8.4, 9.0], f"{REGISTER} :: A-OPM-AXIS ('{q}')")
    h = _json(root, HSE)["variants"]
    src = (f"{HSE} :: variants.<v>['<key>'].ratio and ci95 (oracle known-topography detectability, Neuromag 306 channels, "
           "median over targets with a parcel bootstrap; point gains; v4 = the final arrays and axes, v4_sites_v3_axes = the "
           "final sites with the earlier, nearly radial axes)")
    for v, tag in (("v4", "final"), ("v4_sites_v3_axes", "radial")):
        for key, kname in ((DENSE_PROJ, "dense_proj"), (DENSE_IB, "dense_ib"), (MATCHED_PROJ, "matched_proj")):
            e = h[v][key]
            s = src.replace("<v>", v).replace("<key>", key)
            F.add(f"meth_hse_{tag}_{kname}_ratio_3dp", num(e["ratio"], 3), e["ratio"], s)
            if key == DENSE_PROJ:
                F.add(f"meth_hse_{tag}_{kname}_ci_3dp", f"[{num(e['ci95'][0], 3)}, {num(e['ci95'][1], 3)}]", e["ci95"], s)
                deep = max(e["by_depth"], key=lambda r: r["hi"])  # the deepest band
                if (deep["lo"], deep["hi"]) != (45.0, 70.0):
                    raise ValueError("head_surface_effect: the deepest band is not 45-70 mm")
                F.add(f"meth_hse_{tag}_{kname}_45_70_ratio", num(deep["ratio"], 2), deep["ratio"],
                      f"{HSE} :: variants.{v}['{key}'].by_depth[45-70 mm].ratio (median over the band's targets)")
                if tag == "final":
                    for end in ("lo", "hi"):
                        F.add(f"meth_hse_deep_{end}_mm", const(deep[end]), deep[end], f"{HSE} :: variants.{v}['{key}']."
                              f"by_depth[-1].{end} (the deepest band, mm below the scalp)")


def neuromag_facts(F: Facts, root: Path) -> None:
    g2 = _json(root, G2)
    asd = g2["config"]["sensors"]["squid_asd"]
    q = quoted(root, REGISTER, "manufacturer's typical value: TRIUX datasheet, document NM23083B-A")
    F.add("meth_squid_mag_asd", const(asd["mag_fT_per_rtHz"]), asd["mag_fT_per_rtHz"], f"{G2}, {REGISTER} :: config.sensors."
          f"squid_asd.mag_fT_per_rtHz (fT/sqrt(Hz)); HW-mag-noise ('{q}'), U-HW1")
    F.add("meth_squid_grad_asd", const(asd["grad_fT_per_cm_rtHz"]), asd["grad_fT_per_cm_rtHz"], f"{G2}, {REGISTER} :: "
          f"config.sensors.squid_asd.grad_fT_per_cm_rtHz (fT/(cm sqrt(Hz))); HW-grad-noise ('{q}'), U-HW1")
    m = g2["noise_validation"]["model"]
    F.add("meth_squid_mag_rms_ft", num(m["intrinsic_rms_mag_fT"], 1), m["intrinsic_rms_mag_fT"],
          f"{G2} :: noise_validation.model.intrinsic_rms_mag_fT (in-band RMS of the typical white noise)")
    F.add("meth_squid_grad_rms_ft_cm", num(m["intrinsic_rms_grad_fT_cm"], 1), m["intrinsic_rms_grad_fT_cm"],
          f"{G2} :: noise_validation.model.intrinsic_rms_grad_fT_cm")
    ms = g2["sensitivity"]["squid_measured_spectrum"]["intrinsic_rms"]
    F.add("meth_squid_meas_mag_rms_ft", num(ms["mag_fT"], 1), ms["mag_fT"], f"{G2} :: sensitivity.squid_measured_spectrum."
          "intrinsic_rms.mag_fT (median in-band RMS of the empty-room residual after the 8-term fit; A-G2-SQUIDMEAS)")
    F.add("meth_squid_meas_grad_rms_ft_cm", num(ms["grad_fT_cm"], 1), ms["grad_fT_cm"], f"{G2} :: sensitivity."
          "squid_measured_spectrum.intrinsic_rms.grad_fT_cm")
    q = quoted(root, REGISTER, "planar gradiometer 3014, magnetometer 3024 (T3)", "3014", "3024")
    F.add("meth_coil_grad_type", "3014", 3014, f"{REGISTER} :: HW-T3 ('{q}')")
    F.add("meth_coil_mag_type", "3024", 3024, f"{REGISTER} :: HW-T3 ('{q}')")
    q = quoted(root, REGISTER, "MNE sample recording (MGH Vectorview, 306 channels; stored coils 3012 x 204, 3024 x 102)",
               "3012")
    F.add("meth_coil_stored_grad_type", "3012", 3012, f"{REGISTER} :: HW-geometry ('{q}'); D-T3 (3012 replaced by 3014)")
    q = quoted(root, REGISTER, "MNE 1.13.2 `coil_def.dat`, 'accurate' level (3014: 8 points, 26.39-mm coils, 16.80-mm "
               "baseline; 3024: 16 points, 21-mm coil)", "8 points", "16 points")
    F.add("meth_coil_grad_points", "8", 8, f"{REGISTER} :: HW-coildef ('{q}')")
    F.add("meth_coil_mag_points", "16", 16, f"{REGISTER} :: HW-coildef ('{q}')")
    d = g2["config"]["head_position"]["dewar_spacing_mm"]
    q = quoted(root, REGISTER, "manufacturer's specification: the same TRIUX datasheet")
    F.add("meth_dewar_mm", const(d), d, f"{G2}, {REGISTER} :: config.head_position.dewar_spacing_mm; HW-18mm ('{q}')")


def time_domain_facts(F: Facts, root: Path, ic: IC) -> None:
    """The spike study's time-domain model: one value per parameter, identical in the nine anatomies' runs."""
    g2 = _json(root, G2)
    s = {lab: _json(root, g4_summary(lab)) for lab in LABELS}
    a = s["adult"]
    for lab, x in s.items():  # what every anatomy's run recorded
        for part in ("simulation", "null", "detector", "events"):
            if x["config"][part] != a["config"][part]:
                raise ValueError(f"{g4_summary(lab)}: config.{part} differs from the adult's")
        if x["fs_out"] != a["fs_out"]:
            raise ValueError(f"{g4_summary(lab)}: fs_out differs from the adult's")
    c = a["config"]
    A = g4_summary("adult")
    same = " (identical in all nine g4_<anatomy>_summary.json)"
    if c["simulation"]["band_hz"] != [g2["config"]["band"]["l_freq_hz"], g2["config"]["band"]["h_freq_hz"]]:
        raise ValueError("the spike study's band differs from the G2 band")
    F.add("meth_td_exponent", const(ic.value("IC-TD-EXPONENT")), ic.value("IC-TD-EXPONENT"), ic.src("IC-TD-EXPONENT"))
    F.add("meth_td_fmin_hz", const(ic.value("IC-TD-FMIN-HZ")), ic.value("IC-TD-FMIN-HZ"), ic.src("IC-TD-FMIN-HZ"))
    F.add("meth_td_csd_nperseg", count(ic.value("IC-TD-CSD-NPERSEG")), ic.value("IC-TD-CSD-NPERSEG"), ic.src("IC-TD-CSD-NPERSEG"))
    F.add("meth_td_pad_s", const(ic.value("IC-TD-PAD-S")), ic.value("IC-TD-PAD-S"), ic.src("IC-TD-PAD-S"))
    dec = c["simulation"]["decimate"]
    F.add("meth_td_decimate", const(dec), dec, f"{A} :: config.simulation.decimate{same}; rule: docs/methods.md IC-TD-DECIMATE")
    F.add("meth_td_fs_out_hz", num(a["fs_out"], 1), a["fs_out"], f"{A} :: fs_out (after decimation){same}")
    F.add("meth_td_fs_sim_hz", num(a["fs_out"] * dec, 1), a["fs_out"] * dec, f"{A} :: derived: fs_out x config.simulation."
          "decimate (the sample recording's rate, at which the noise is simulated)")
    for key, name in (("segment_s", "segment_s"), ("event_spacing_s", "event_spacing_s"), ("opm_asd_fT_per_rtHz", "opm_asd")):
        v = c["simulation"][key]
        F.add(f"meth_td_{name}", const(v), v, f"{A} :: config.simulation.{key}{same}")
    v = c["detector"]["dictionary_spacing_mm"]
    F.add("meth_td_dictionary_mm", const(v), v, f"{A} :: config.detector.dictionary_spacing_mm{same}")
    for key, name in (("baseline_min", "whitener"), ("calibration_min", "calibration"), ("heldout_min", "heldout")):
        v = c["null"][key]
        F.add(f"meth_null_{name}_min", const(v), v, f"{A} :: config.null.{key}{same}")
    ev = _toml(root, CFG_CONFIRM)["null"]["evaluation_min"]
    F.add("meth_null_evaluation_min", const(ev), ev, f"{CFG_CONFIRM} :: null.evaluation_min (confirmatory run)")
    nd = {lab: x["n_dictionary"] for lab, x in s.items()}
    for lab, v in nd.items():
        F.add(f"meth_td_candidates_{TOK[lab]}", count(v), v, f"{g4_summary(lab)} :: n_dictionary (practical detector's candidate "
              "field patterns: the 10-mm Poisson-disk grid without the true sources)")
    F.spread("meth_td_candidates", list(nd.values()), count, f"{ALL_G4} :: n_dictionary over the nine anatomies")
    bw = _json(root, BAND)["opm_response"]
    F.add("meth_opm_bw_corner_hz", const(bw["corner_hz"]), bw["corner_hz"], f"{BAND} :: opm_response.corner_hz ({bw['model']}; "
          "A-OPM-BW, band sensitivity analysis only: docs/methods.md IC-TD-RESPONSE)")
    lc = {lab: _json(root, loc_summary(lab))["config"] for lab in LABELS}
    for lab, x in lc.items():
        if x != lc["adult"]:
            raise ValueError(f"{loc_summary(lab)}: config differs from the adult's")
    L = loc_summary("adult")
    for key, name in (("baseline_min", "cov"), ("calibration_min", "threshold"), ("holdout_min", "heldout")):
        v = lc["adult"][key]
        F.add(f"meth_loc_{name}_min", const(v), v, f"{L} :: config.{key} (identical in all nine localization summaries)")
    F.add("meth_dspm_depth", const(lc["adult"]["depth"]), lc["adult"]["depth"], f"{L} :: config.depth (depth-weighting "
          "exponent, both dSPM implementations)")
    snr = lc["adult"]["snr"]
    F.add("meth_dspm_snr", const(snr), snr, f"{L} :: config.snr")
    F.add("meth_dspm_lambda2", f"1/{const(snr * snr)}", 1.0 / snr**2, f"{L} :: derived: lambda^2 = 1 / config.snr^2")
    F.add("meth_mne_depth_limit", const(ic.value("IC-DSPM-MNE-LIMIT")), ic.value("IC-DSPM-MNE-LIMIT"),
          ic.src("IC-DSPM-MNE-LIMIT", "IC-DSPM-MNE-CHS"))


def bootstrap_facts(F: Facts, root: Path, ic: IC) -> None:
    """Resamples of every interval family (configurations where they hold the number, else the implementation constant)."""
    for ident, name in (("IC-BOOT-G2", "g2"), ("IC-BOOT-G2-SENS", "g2_sens"), ("IC-BOOT-G2-PATCH", "g2_patch"),
                        ("IC-BOOT-G2-BAND", "g2_band"), ("IC-BOOT-HSE", "hse"), ("IC-BOOT-G3B-MATCHED", "g3b_matched"),
                        ("IC-BOOT-G3B-PLACEMENTS", "g3b_placements"), ("IC-BOOT-G3B-CHANNELS", "g3b_channels"),
                        ("IC-BOOT-G3B-USEFUL", "g3b_useful"), ("IC-BOOT-LOC", "loc"), ("IC-BOOT-CGAP", "cgap_other"),
                        ("IC-BOOT-COVVAL", "covval")):
        v = ic.value(ident)
        F.add(f"meth_boot_{name}", count(v), v, ic.src(ident))
    v = ic.value("IC-BOOT-G4")
    q = quoted(root, METHODS, f"S50 intervals come from a bootstrap over locations ({count(v)} resamples)", count(v))
    F.add("meth_boot_g4", count(v), v, ic.src("IC-BOOT-G4") + f"; section 9 ('{q}')")
    g2 = _json(root, G2)
    methods = {e["ci_method"] for e in g2["primary"]["oracle"].values()}
    if len(methods) != 1 or not (m := re.fullmatch(r"parcels \((\d+)\)", methods.pop())):
        raise ValueError(f"{G2}: unexpected ci_method of the primary comparisons")
    F.add("meth_boot_g2_parcels", m.group(1), int(m.group(1)), f"{G2} :: primary.oracle[*].ci_method ('parcels (70)': the 68 "
          "Desikan-Killiany parcels and the two medial-wall labels)")
    g3b = _json(root, G3B)
    nb = g3b["config"]["strata"]["n_boot"]
    if nb != _toml(root, CFG_G3B)["strata"]["n_boot"]:
        raise ValueError("G3B: stored and configured n_boot differ")
    F.add("meth_boot_g3b_dense", count(nb), nb, f"{G3B} :: config.strata.n_boot (the dense array at top contact: D, Delta and "
          "the depth and orientation strata; docs/methods.md IC-BOOT-G3B-MATCHED)")
    F.add("meth_boot_cgap_primary", count(nb), nb, f"{CFG_G3B} :: strata.n_boot (the gap-matched helmet's primary comparisons: "
          "docs/methods.md IC-BOOT-CGAP)")
    v = _toml(root, CFG_CONFIRM)["design"]["bootstrap_resamples"]
    F.add("meth_boot_confirm", count(v), v, f"{CFG_CONFIRM} :: design.bootstrap_resamples (location bootstrap)")
    v = _toml(root, CFG_NOISE)["general"]["n_boot"]
    F.add("meth_boot_noise_sens", count(v), v, f"{CFG_NOISE} :: general.n_boot (parcel bootstrap of the noise-model sensitivity "
          "analyses)")
    v = _toml(root, CFG_MOTION)["coupling"]["n_draws"]
    F.add("meth_motion_draws", count(v), v, f"{CFG_MOTION} :: coupling.n_draws (field and calibration draws; no resampling)")
    q = quoted(root, METHODS, "with the 10th-90th percentiles of the per-draw thresholds", "10th-90th")
    F.add("meth_motion_pct_lo", "10th", 10, f"{METHODS} :: section 12 ('{q}')")
    F.add("meth_motion_pct_hi", "90th", 90, f"{METHODS} :: section 12 ('{q}')")


def seed_facts(F: Facts, root: Path, ic: IC) -> None:
    g2 = _json(root, G2)
    v = g2["config"]["sources"]["seed"]
    F.add("meth_seed_g2", seed(v), v, f"{G2} :: config.sources.seed (the adult's targets, background grid, plug-in samples and "
          "bootstrap; every smaller head's targets and grid: docs/methods.md IC-SEED-G3B-GRID)")
    F.add("meth_seed_g3b_boot", seed(ic.value("IC-SEED-G3B")), ic.value("IC-SEED-G3B"), ic.src("IC-SEED-G3B"))
    F.add("meth_seed_g3b_school_boot", seed(ic.value("IC-SEED-G3B-SCHOOL")), ic.value("IC-SEED-G3B-SCHOOL"),
          ic.src("IC-SEED-G3B-SCHOOL"))
    sims = {lab: _json(root, g4_summary(lab))["config"]["simulation"]["seed"] for lab in LABELS}
    if len(set(sims.values())) != 1:
        raise ValueError(f"spike detection seeds differ between anatomies: {sims}")
    s = sims["adult"]
    F.add("meth_seed_g4_detection", seed(s), s, f"{ALL_G4} :: config.simulation.seed (the same in all nine anatomies' runs; "
          "docs/methods.md IC-SEED-G4-SIM)")
    locs = {lab: _json(root, loc_summary(lab))["config"]["seed"] for lab in LABELS}
    if len(set(locs.values())) != 1:
        raise ValueError(f"localization seeds differ between anatomies: {locs}")
    ls = locs["adult"]
    F.add("meth_seed_g4_localization", seed(ls), ls, f"{loc_summary('adult')} :: config.seed (the same in all nine localization "
          "summaries)")
    for ident, name in (("IC-SEED-LOC-HELDOUT", "loc_heldout"), ("IC-SEED-LOC-SECONDARY", "loc_secondary")):
        v = ls + ic.value(ident)
        F.add(f"meth_seed_{name}", seed(v), v, f"{loc_summary('adult')}, {METHODS} :: derived: config.seed + the offset of "
              f"section 13 {ident} (set in {ic.rows[ident]['where']})")
    F.add("meth_seed_g4_boot", seed(ic.value("IC-SEED-G4-BOOT")), ic.value("IC-SEED-G4-BOOT"), ic.src("IC-SEED-G4-BOOT"))
    v = _toml(root, CFG_CONFIRM)["design"]["root_seed"]
    F.add("meth_seed_confirm_root", seed(v), v, f"{CFG_CONFIRM} :: design.root_seed ([declared] root of every confirmatory stream)")
    F.add("meth_confirm_purposes", count(ic.value("IC-SEED-CONFIRM")), ic.value("IC-SEED-CONFIRM"), ic.src("IC-SEED-CONFIRM"))


def skull_facts(F: Facts, root: Path) -> None:
    """The 12-month template's thin skull and the conductivity-free (1-layer) check of every smaller head."""
    q = quoted(root, REGISTER, "The 12-month template's skull is 0.25 mm thick at its thinnest (23 of 2,562 inner-skull vertices "
               "within 1 mm of the outer skull), as distributed; the 1-layer BEM check does not depend on it.",
               "0.25 mm", "23 of 2,562", "1 mm")
    g3b = _json(root, G3B)
    tpl = g3b["config"]["anatomy"]["templates"]
    if "ANTS12-0Months3T" not in tpl:
        raise ValueError(f"{G3B}: the 12-month template is not among {tpl}")
    F.add("meth_infant12mo_age_months", "12", 12, f"{G3B} :: config.anatomy.templates ('ANTS12-0Months3T': the 12-month "
          "template, O'Reilly et al. 2021)")
    for name, v, raw in (("meth_infant12mo_skull_min_mm", "0.25", 0.25), ("meth_infant12mo_thin_vertices", "23", 23),
                         ("meth_infant12mo_inner_skull_vertices", "2,562", 2562), ("meth_infant12mo_thin_within_mm", "1", 1)):
        F.add(name, v, raw, f"{REGISTER} :: D-G3-ANAT ('{q}')")
    q = quoted(root, METHODS, "MNE-Python 1.13.2 BEM forward models (single-compartment inner skull, 0.3 S/m, unless a paper "
               "configuration prescribes otherwise", "0.3 S/m")
    F.add("meth_bem1_conductivity", "0.3", 0.3, f"{METHODS} :: section 3 ('{q}')")
    cond = g3b["config"]["anatomy"]["bem_conductivity"]
    F.add("meth_conductivity", const(cond), cond, f"{G3B} :: config.anatomy.bem_conductivity (S/m, every head)")
    D, S = g3b["D_median_dB"], g3b["sensitivity_median_D_dB"]
    key = "opm_dense/combined/intrinsic+brain"
    d12, d12_1 = D[f"infant12mo/{key}/detect"], S[f"infant12mo/bem1/{key}"]
    F.add("meth_infant12mo_d_db", db(d12), d12, f"{G3B} :: D_median_dB['infant12mo/{key}/detect'] (three-layer BEM, top contact)")
    F.add("meth_infant12mo_d_bem1_db", db(d12_1), d12_1, f"{G3B} :: sensitivity_median_D_dB['infant12mo/bem1/{key}'] (1-layer BEM)")
    da = D[f"adult/{key}/detect"]
    F.add("meth_infant12mo_d_minus_adult_db", db(d12 - da), d12 - da, f"{G3B} :: derived: D_median_dB['infant12mo/{key}/detect'] "
          f"- D_median_dB['adult/{key}/detect'] (difference of the medians)")
    dd = S[f"delta/infant12mo/bem1/{key}"]
    F.add("meth_infant12mo_d_minus_adult_bem1_db", db(dd), dd, f"{G3B} :: sensitivity_median_D_dB['delta/infant12mo/bem1/{key}'] "
          "(difference of the medians, 1-layer BEM for both heads)")
    change = [S[f"delta/{lab}/bem1/{key}"] - (D[f"{lab}/{key}/detect"] - da) for lab in SMALLER]
    F.spread("meth_bem1_change_db", change, db, f"{G3B} :: derived: sensitivity_median_D_dB['delta/<head>/bem1/{key}'] - "
             f"(D_median_dB['<head>/{key}/detect'] - D_median_dB['adult/{key}/detect']) over the eight smaller heads (how much "
             "the 1-layer BEM changes each head's difference from the adult)")


def location_facts(F: Facts, root: Path, ic: IC) -> None:
    """What the stored location design allows the exact sign-flip test to show (no power is computed)."""
    ped = _json(root, PED)["comparison"]
    per = {v for lab in LABELS for v in ped[f"{lab}/n_locations_per_depth_band"].values()}
    if len(per) != 1:
        raise ValueError(f"{PED}: locations per depth band differ: {per}")
    n18 = per.pop()
    F.add("meth_locs_per_band", count(n18), n18, f"{PED} :: comparison['<anatomy>/n_locations_per_depth_band'] (every band of "
          "all nine anatomies)")
    adult = _json(root, g4_summary("adult"))
    n_or = len({loc["stratum"][1] for loc in adult["locations"]})
    F.add("meth_orientation_strata", count(n_or), n_or, f"{g4_summary('adult')} :: derived: distinct locations[].stratum[1]")
    cc = _toml(root, CFG_CONFIRM)
    k = cc["design"]["locations_per_stratum"]
    F.add("meth_confirm_locs_per_stratum", count(k), k, f"{CFG_CONFIRM} :: design.locations_per_stratum ([declared])")
    n36 = k * n_or
    F.add("meth_confirm_locs", count(n36), n36, f"{CFG_CONFIRM}, {g4_summary('adult')} :: derived: design.locations_per_stratum x "
          "the orientation strata of the 10-20 mm band")
    alpha = cc["endpoint"]["alpha"]
    n_anat = len(cc["endpoint"]["anatomies"])
    F.add("meth_alpha", const(alpha), alpha, f"{CFG_CONFIRM} :: endpoint.alpha")
    F.add("meth_n_anatomies", count(n_anat), n_anat, f"{CFG_CONFIRM} :: endpoint.anatomies (count)")
    holm1 = alpha / n_anat
    F.add("meth_holm_first_alpha", sig(holm1, 2), holm1, f"{CFG_CONFIRM} :: derived: endpoint.alpha / the number of anatomies "
          "(the first, strictest step of Holm's procedure)")
    exact_max, n_mc = ic.value("IC-SIGNFLIP-EXACT"), ic.value("IC-SIGNFLIP-MC")
    F.add("meth_signflip_exact_max", count(exact_max), exact_max, ic.src("IC-SIGNFLIP-EXACT"))
    F.add("meth_signflip_mc_patterns", count(n_mc), n_mc, ic.src("IC-SIGNFLIP-MC"))
    F.add("meth_signflip_mc_resolution", sci(1.0 / n_mc, 1), 1.0 / n_mc, ic.src("IC-SIGNFLIP-MC") + "; derived: 1 / the number "
          "of sign patterns")
    design = f"{PED}, {CFG_CONFIRM}, {METHODS} :: derived from the design (18 and 36 locations; alpha) and the test of section 13 "
    for n in (n18, n36):  # all n locations non-zero and concordant: the two extreme sign patterns of 2^n
        p = 2.0 ** (1 - n)
        F.add(f"meth_signflip_min_p_{n}", sci(p), p, design + f"IC-SIGNFLIP-EXACT: smallest attainable two-sided p, 2 / 2^{n}")
    n_min = next(n for n in range(1, 64) if 2.0 ** (1 - n) < alpha)
    F.add("meth_signflip_min_n", count(n_min), n_min, design + "IC-SIGNFLIP-EXACT: fewest non-zero locations for which p < alpha "
          "is attainable at all (all concordant)")
    F.add("meth_signflip_min_n_p", pval(2.0 ** (1 - n_min)), 2.0 ** (1 - n_min), design + f"IC-SIGNFLIP-EXACT: 2 / 2^{n_min}")

    def sign_test(n, kk):  # two-sided exact sign test = the sign-flip test when every |difference| is equal
        return min(1.0, 2.0 * sum(math.comb(n, j) for j in range(kk, n + 1)) / 2.0**n)

    for n in (n18, n36):
        for lab, a in (("", alpha), ("_holm", holm1)):
            kk = next(kk for kk in range(n // 2 + 1, n + 1) if sign_test(n, kk) < a)
            how = (design + f"IC-SIGNFLIP-EXACT, with equal per-location differences (the exact two-sided sign test): fewest of {n} "
                   f"non-zero locations favouring one array for p < {sig(a, 2)}")
            F.add(f"meth_signtest_{n}_k{lab}", count(kk), kk, how)
            F.add(f"meth_signtest_{n}_k{lab}_p", pval(sign_test(n, kk)), sign_test(n, kk), how + " (its p)")
            if not lab:
                F.add(f"meth_signtest_{n}_kminus1", count(kk - 1), kk - 1, how + " (one location fewer)")
                F.add(f"meth_signtest_{n}_kminus1_p", pval(sign_test(n, kk - 1)), sign_test(n, kk - 1), how + " (p one location "
                      "fewer)")
    e = adult["paired"]["opm_dense/opm_vs_squid/combined/practical@1"]["depth0"]
    o, q = e["locations_favouring_opm"], e["locations_favouring_squid"]
    src = f"{g4_summary('adult')} :: paired['opm_dense/opm_vs_squid/combined/practical@1'].depth0"
    F.add("meth_adult_locs_opm", count(o), o, f"{src}.locations_favouring_opm")
    F.add("meth_adult_locs_squid", count(q), q, f"{src}.locations_favouring_squid")
    F.add("meth_adult_p", pval(e["location_sign_flip_p"]), e["location_sign_flip_p"], f"{src}.location_sign_flip_p (stored)")
    pe = sign_test(o + q, max(o, q))
    F.add("meth_adult_equal_p", pval(pe), pe, f"{src}; derived: the exact two-sided sign test of the same split (the "
          "sign-flip p if every location's difference had the same size)")


def provenance_facts(F: Facts, root: Path) -> None:
    """Declared parameters and their status (register entries; figure readings in the preprint extraction)."""
    s = _json(root, G2)["config"]["sensors"]
    v = s["opm_asd_primary_fT_per_rtHz"]
    q = quoted(root, REGISTER, "15 fT/sqrt(Hz), the middle of the A-OPM-NOISE sweep, not a device value", "15")
    F.add("meth_prov_opm_asd_primary", const(v), v, f"{G2}, {REGISTER} :: config.sensors.opm_asd_primary_fT_per_rtHz; "
          f"A-G2-OPMNOISE ('{q}')")
    sw = s["opm_asd_fT_per_rtHz"]
    F.add("meth_prov_opm_asd_sweep", span(min(sw), max(sw), const), [min(sw), max(sw)], f"{G2} :: config.sensors."
          "opm_asd_fT_per_rtHz (the declared sweep, fT/sqrt(Hz))")
    q = quoted(root, REGISTER, "~7 mm (2-mm helmet shell + half of a 10-mm cell) | paper p. 10", "7 mm")
    F.add("meth_prov_standoff_mm", "7", 7, f"{REGISTER} :: J-opm-standoff ('{q}')")
    q = quoted(root, REGISTER, "10-mm cube, FieldLine Gen2, single axis, normal component | paper p. 10", "10-mm")
    F.add("meth_prov_cell_mm", "10", 10, f"{REGISTER} :: J-opm-cell ('{q}')")
    q = quoted(root, REGISTER, "Minimum sensing-centre spacing 17 mm (10-mm cell in a ~12-17 mm package)", "17 mm", "12-17")
    F.add("meth_prov_pack_mm", "17", 17, f"{REGISTER} :: U-OPM-PACK ('{q}'; a modelling choice; device footprints in the register row)")
    F.add("meth_prov_package_range_mm", "12 to 17", [12, 17], f"{REGISTER} :: U-OPM-PACK ('{q}')")
    q = quoted(root, REGISTER, "SERF OPMs of this class have ~100-150 Hz bandwidth", "100-150")
    F.add("meth_prov_bw_range_hz", "100 to 150", [100, 150], f"{REGISTER} :: A-OPM-BW ('{q}')")
    q = quoted(root, JAS, "OPM empty room ≈29–31 dB re fT²/Hz, i.e. about 28–35 fT/√Hz, if the axis unit is dB re 1 fT²/Hz.",
               "29–31", "28–35")
    F.add("meth_prov_fig9_db_range", "29 to 31", [29, 31], f"{JAS} :: section 5.6 'Spectra (p. 20; Fig. 9)', Levels [fig] ('{q}')")
    F.add("meth_prov_fig9_asd_range", "28 to 35", [28, 35], f"{JAS} :: section 5.6 'Spectra (p. 20; Fig. 9)', Levels [fig] ('{q}')")
    q = quoted(root, JAS, "I measured these from the 600-dpi raster; each should be within about ±2 %.", "600-dpi", "±2 %")
    F.add("meth_prov_raster_dpi", "600", 600, f"{JAS} :: section 5.6 'Fig. 8 bar values [fig]' ('{q}')")
    F.add("meth_prov_fig8_precision_pct", "2%", 2, f"{JAS} :: section 5.6 'Fig. 8 bar values [fig]' ('{q}')")


def facts(root: Path = ROOT) -> dict:
    root = Path(root)
    F, ic = Facts(), IC(root)
    background_facts(F, root, ic)
    projection_facts(F, root)
    axis_facts(F, root)
    neuromag_facts(F, root)
    time_domain_facts(F, root, ic)
    bootstrap_facts(F, root, ic)
    seed_facts(F, root, ic)
    skull_facts(F, root)
    location_facts(F, root, ic)
    provenance_facts(F, root)
    return dict(F)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grep", default="", help="print only facts whose name contains this text")
    args = ap.parse_args(argv)
    for name, f in facts(ROOT).items():
        if args.grep in name:
            print(f"{name} = {f['value']}    [{f['source']}]")


if __name__ == "__main__":
    main()
