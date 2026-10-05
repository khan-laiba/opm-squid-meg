#!/usr/bin/env python3
"""Report facts for G1A, G1B, G1C, G2 and the adult regions (prefixes g1a_, g1b_, g1c_, g2_, reg_, lit_).

facts(root) -> {name: {"value", "raw", "source"}}, merged by scripts/report_facts.py. Every number is read or
derived from the committed result files (results/g1a, g1b, g1c, g2), configs/ or docs/ (lit_: the preprint values
transcribed or read from its figures in docs/literature/jas2026.md); nothing is re-run.

Formats: ratios 2 decimals, with a _3dp twin near 1 (0.95-1.05) or where 3 decimals are quoted elsewhere; dB
differences signed (+0.44, −0.38), dB levels 2 decimals; intervals "[lo, hi]" in separate _ci facts; mm 1 decimal;
nAm integer; counts with thousands separators; percentages integer with "%" (56%; 1 decimal below 10 %, 2 significant
digits below 0.1 %); ranges "lo to hi" from unrounded values, with _min and _max facts; negatives use U+2212.
"raw" holds the unrounded value in the printed unit (percent for % facts, ratio for log2 summaries).
G2 names: g2_[estimator_]<array>_vs_<comparator>_<condition>_<what>; arrays dense (208 sites), matched (98),
opm204; comparators combined (306 channels), mag (102), grad (204); conditions int (intrinsic), ib (intrinsic +
brain), ibenv (+ room field), proj (projected); estimators peak_, meanpow_, plugin10s_, plugin60s_ (none: oracle
known-topography detectability). Depth bins are named <lo>_<hi> in mm.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np

MINUS = "−"
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
G1A, FIG3 = "results/g1a/g1a_benchmark.json", "results/g1a/fig3/figure3_summary.json"
G1B, G1C = "results/g1b/g1b_summary.json", "results/g1c/g1c_summary.json"
G2, BAND = "results/g2/g2_summary.json", "results/g2/g2_band_sensitivity.json"
TARGETS, PATCHES = "results/g2/g2_targets.csv", "results/g2/g2_patch_targets.csv"
HSE, NEAR = "results/g2/head_surface_effect.json", "results/g2/near_mesh_check.json"
SKIN, SKIN_V1 = "results/g2/bem_skin_refinement.json", "results/g2/bem_skin_refinement_v1_arrays.json"
SPHERE, JAS = "results/g2/bem_sphere_check.json", "docs/literature/jas2026.md"
METHODS, REGISTER, CFG = "docs/methods.md", "docs/provenance_register.md", "configs/g2_adult.toml"
LITJSON, LITMD = "docs/literature/epilepsy_opm_studies.json", "docs/literature/epilepsy_opm_studies.md"

ARRAYS = {"opm_dense": "dense", "opm_matched": "matched", "opm204": "opm204"}
COMPS = ("combined", "mag", "grad")
CONDS = {"intrinsic": "int", "intrinsic+brain": "ib", "intrinsic+brain+env": "ibenv", "projected": "proj"}
ESTIMATORS = {"oracle": "", "peak": "peak_", "meanpow_db": "meanpow_", "plugin_T10": "plugin10s_", "plugin_T60": "plugin60s_"}
REGIONS = {"superiortemporal": ["superiortemporal"], "parahippocampal": ["parahippocampal"], "entorhinal": ["entorhinal"],
           "fusiform": ["fusiform"], "temporalpole": ["temporalpole"], "precentral": ["precentral"],
           "mesial_temporal": ["parahippocampal", "entorhinal"],
           "lateral_temporal": ["superiortemporal", "middletemporal", "inferiortemporal", "bankssts", "transversetemporal"]}


# ------------------------------------------------------------------------------------------------
# formatting
def num(x, nd):
    """Fixed decimals, U+2212 for negatives, no sign on a value that rounds to zero."""
    s = f"{x:.{nd}f}"
    return s.lstrip("-") if float(s) == 0 else s.replace("-", MINUS)


def signed(x, nd=2):
    s = num(x, nd)
    return s if float(s.replace(MINUS, "-")) == 0 or s.startswith(MINUS) else "+" + s


def sig2(x):
    """Two significant digits, trailing zeros kept (0.0070, 0.48, 6.0, 12)."""
    if x == 0:
        return "0"
    nd = max(0, 1 - math.floor(math.log10(abs(x))))
    s = f"{x:.{nd}f}"
    if nd > 0 and abs(float(s)) >= 10 ** (2 - nd):  # rounding carried into a new digit (9.96 -> 10)
        s = f"{x:.{nd - 1}f}"
    return s.replace("-", MINUS)


def pct(p):
    """A percentage (already x 100): integer from 10 %, 1 decimal below, 2 significant digits below 0.1 %."""
    if p == 0:
        return "0%"
    if abs(p) < 0.095:
        return sig2(p) + "%"
    s = f"{p:.1f}"
    return (f"{p:.0f}" if abs(float(s)) >= 10 else s).replace("-", MINUS) + "%"


def count(n):
    return f"{int(round(n)):,}"


def sci(x, sig=2):
    m, e = f"{x:.{sig - 1}e}".split("e")
    return f"{m} × 10{str(int(e)).translate(SUP)}"


def eta_text(x):
    return f"{x:.2f}".rstrip("0").rstrip(".") if x != int(x) else f"{x:.1f}"


def tag(x):
    """A number as a name part: 4.5 -> 4p5, 10.0 -> 10."""
    return (f"{x:g}").replace(".", "p").replace("-", "m")


def r2(x):
    return num(x, 2)


def r3(x):
    return num(x, 3)


def _py(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_py(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return float(x)
    return x


class Facts(dict):
    def add(self, name, value, raw, source):
        if name in self:
            raise ValueError(f"fact {name!r} defined twice")
        if not name.replace("_", "").isalnum() or name != name.lower():
            raise ValueError(f"bad fact name {name!r}")
        self[name] = {"value": value, "raw": _py(raw), "source": source}

    def ratio(self, base, r, src, ci=None, three=False):
        """<base>_ratio (2 dp), _3dp twin near 1 or on request, and the interval as <base>_ci (and _ci_3dp)."""
        three = three or 0.95 <= r <= 1.05
        self.add(f"{base}_ratio", r2(r), r, src)
        if three:
            self.add(f"{base}_ratio_3dp", r3(r), r, src)
        if ci is not None:
            self.add(f"{base}_ci", f"[{r2(ci[0])}, {r2(ci[1])}]", list(ci), src)
            if three:
                self.add(f"{base}_ci_3dp", f"[{r3(ci[0])}, {r3(ci[1])}]", list(ci), src)

    def pct2(self, name, p, src):
        """A two-significant-digit twin of a percentage, added only where it differs from pct()."""
        if sig2(p) + "%" != pct(p):
            self.add(name, sig2(p) + "%", p, src)

    def span(self, base, values, fmt, src, extra=()):
        """<base>_min, _max and _range ('lo to hi') from unrounded values; extra: (suffix, fmt) range twins."""
        lo, hi = float(min(values)), float(max(values))
        self.add(f"{base}_min", fmt(lo), lo, src)
        self.add(f"{base}_max", fmt(hi), hi, src)
        self.add(f"{base}_range", f"{fmt(lo)} to {fmt(hi)}", [lo, hi], src)
        for suffix, f in extra:
            self.add(f"{base}_range_{suffix}", f"{f(lo)} to {f(hi)}", [lo, hi], src)


def load(root, rel):
    return json.loads((Path(root) / rel).read_text())


def read_csv(root, rel):
    """Rows of a result CSV; the first line is a '#' provenance comment."""
    with open(Path(root) / rel) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def cmp_entry(F, base, e, file, key, three=False, shares=True, parcels=False, counts=False):
    """A stored comparison (median_log2, ci95, share_opm_better, ...) as ratio, CI and shares."""
    r = 2 ** e["median_log2"]
    ci = [2 ** e["ci95"][0], 2 ** e["ci95"][1]]
    how = e.get("ci_method", "")
    how = "parcel bootstrap, " + how if how.startswith("parcels") else ("target bootstrap" if how == "targets" else how)
    src = (f"{file} :: {key} (ratio = 2**median_log2, the median over {e['n']:,} targets of d_OPM/d_SQUID; "
           f"95 % CI = 2**ci95, {how})")
    F.ratio(base, r, src, ci, three)
    if shares:
        F.add(f"{base}_share", pct(100 * e["share_opm_better"]), 100 * e["share_opm_better"],
              f"{file} :: {key}.share_opm_better (share of targets with d_OPM > d_SQUID)")
    if parcels and "share_parcels_opm_better" in e:
        F.add(f"{base}_share_parcels", pct(100 * e["share_parcels_opm_better"]), 100 * e["share_parcels_opm_better"],
              f"{file} :: {key}.share_parcels_opm_better (share of parcels whose median favours the OPM)")
    if counts:
        n_opm = e["share_opm_better"] * e["n"]
        F.add(f"{base}_n_opm_better", count(n_opm), round(n_opm),
              f"{file} :: derived: {key}.share_opm_better x n ({e['n']:,} targets)")
        F.add(f"{base}_n_squid_better", count(e["n"] - n_opm), round(e["n"] - n_opm),
              f"{file} :: derived: (1 - {key}.share_opm_better) x n ({e['n']:,} targets)")


# ------------------------------------------------------------------------------------------------
def lit_facts(F):
    """Preprint values (Jas et al. 2026) as transcribed or read from its figures in docs/literature/jas2026.md."""
    s = JAS + " :: "
    for name, v, where in (("lit_deq_eta3_printed_mm", 28, "section 3.4, Fig. 3 d_eq as printed (pp. 13-14)"),
                           ("lit_deq_eta2p5_printed_mm", 34, "section 3.4, Fig. 4B d_eq as printed (p. 14)"),
                           ("lit_deq_eta4_printed_mm", 19, "section 3.4, Fig. 4C d_eq as printed (p. 14)"),
                           ("lit_q_nam", 30, "section 2.1, dipole moment Q (p. 6)"),
                           ("lit_xi_squid_mm", 18, "section 2.1, SQUID standoff xi (Triux neo minimum, p. 8)"),
                           ("lit_xi_opm_mm", 0, "section 2.1, OPM on-scalp standoff (p. 8)"),
                           ("lit_eta_min", 1, "section 2.3, eta swept from 1 (p. 9)"),
                           ("lit_eta_max", 6, "section 2.3, eta swept to 6 (p. 9)"),
                           ("lit_dewar_wall_mm", 18, "section 5.3, Triux neo Dewar wall (p. 11)"),
                           ("lit_opm_standoff_mm", 7, "section 5.2, OPM sensing centre ~7 mm from the helmet inner surface (p. 10)"),
                           ("lit_opm_cell_mm", 10, "section 5.2, cubical vapour cell side length (p. 10)"),
                           ("lit_sef_n_opm", 11, "section 5.2, single-axis FieldLine Gen2 OPMs (p. 10)"),
                           ("lit_sef_n_trials_per_run", 500, "section 5.1, trials per run (p. 10)"),
                           ("lit_sef_age_years", 62, "section 5.1, the one participant's age (p. 10)"),
                           ("lit_sef_highpass_hz", 4, "section 5.5, high-pass filter (p. 12)"),
                           ("lit_sef_antialias_hz", 330, "sections 5.2-5.3, anti-alias low-pass (pp. 10-11)"),
                           ("lit_n20_latency_ms", 24, "section 5.5, N20 latency from the Fig. 7B labels (p. 18)"),
                           ("lit_n20_amplitude_ratio_printed", 3, "section 7 item 9, N20 'about 3 times' larger in OPM (p. 18)")):
        F.add(name, str(v), v, s + where)
    F.add("lit_eta0_printed", "1.7", 1.7, s + "section 3.4, Fig. 4E label eta0 (p. 15)")
    F.add("lit_eta1_printed", "5.3", 5.3, s + "section 3.4, Fig. 4E label eta1 (p. 15)")
    F.add("lit_eta_range", "1 to 6", [1, 6], s + "section 2.3 (p. 9)")
    F.add("lit_eta_reference", "3", 3, s + "section 7 item 1 and section 3.4: the adult model's eta = 3, at which d_eq is printed "
          "as 28 mm (Fig. 3, pp. 13-14) and the normalised d_eq of Fig. 5B is labelled (p. 16)")
    table1 = {"newborn": (55, 48), "1yr": (70, 62), "8yr": (85, 73), "adult": (95, 80)}
    for head, (h, b) in table1.items():
        F.add(f"lit_table1_{head}_h_mm", str(h), h, s + "section 3.1, Table 1 head radius h (p. 7)")
        F.add(f"lit_table1_{head}_b_mm", str(b), b, s + "section 3.1, Table 1 brain radius b (p. 7)")
        F.add(f"lit_table1_{head}_extracerebral_mm", str(h - b), h - b, s + "section 3.1, Table 1 (h - b), CSF + skull + scalp (p. 7)")
    F.add("lit_table1_extracerebral_range_mm", "7 to 15", [7, 15], s + "section 3.1, Table 1 (h - b), newborn to adult")
    F.add("lit_norm_deq_newborn_printed_pct", "50%", 50, s + "section 3.3, Fig. 5B normalised d_eq at eta = 3, 'infant/newborn' (p. 16)")
    F.add("lit_norm_deq_adult_printed_pct", "15%", 15, s + "section 3.3, Fig. 5B normalised d_eq at eta = 3, adult (p. 16)")
    F.add("lit_noise_squid_range", "2 to 5", [2, 5], s + "section 5.6, literature intrinsic noise of SQUIDs, fT/sqrt(Hz) (Brookes 2022; p. 4)")
    F.add("lit_noise_opm_range", "7 to 30", [7, 30], s + "section 5.6, literature intrinsic noise of OPMs, fT/sqrt(Hz) (Brookes 2022; p. 4)")
    F.add("lit_opm_signal_gain_range", "4 to 8", [4, 8], s + "section 5.6, OPM signals 4-8x off-scalp SQUID (Brickwedde 2024; p. 3)")
    # Fig. 8 bar values (figure readings, about +/-2 % each); C3 (OPM) and MEG0431 (SQUID)
    fig8 = {"opm_sigma_n20_xi0_ft": 36.1, "squid_sigma_n20_xi0_ft": 7.9, "opm_sigma_n20_xi0rep_ft": 40.2,
            "squid_sigma_n20_xi0rep_ft": 9.7, "opm_sigma_n20_xi2_ft": 36.3, "opm_sigma_n20_empty_room_ft": 16.3,
            "squid_sigma_n20_empty_room_ft": 1.6, "opm_b_n20_xi0_ft": 347, "squid_b_n20_xi0_ft": 132,
            "opm_snr_xi0": 9.51, "squid_snr_xi0": 16.6}
    for k, v in fig8.items():
        F.add(f"lit_sef_{k}", num(v, 0 if k.startswith(("opm_b", "squid_b")) else 1), v,
              s + f"section 5.6, Fig. 8 bar value (figure reading, about +/-2 %): {k}")
    eta = fig8["opm_sigma_n20_xi0_ft"] / fig8["squid_sigma_n20_xi0_ft"]
    F.add("lit_sef_noise_ratio_closest", num(eta, 1), eta, s + "derived: OPM / SQUID sigma_N20 at the closest runs, 36.1 / 7.9 "
          "(section 5.6 'Measured eta'; register J-eta-meas ~4.6; 1 decimal, figure readings)")
    eta_rep = fig8["opm_sigma_n20_xi0rep_ft"] / fig8["squid_sigma_n20_xi0rep_ft"]
    F.add("lit_sef_noise_ratio_repeat", num(eta_rep, 1), eta_rep, s + "derived: 40.2 / 9.7 at the repeat runs (section 5.6, ~4.1)")
    eta_er = fig8["opm_sigma_n20_empty_room_ft"] / fig8["squid_sigma_n20_empty_room_ft"]
    F.add("lit_sef_noise_ratio_empty_room", num(eta_er, 0), eta_er, s + "derived: 16.3 / 1.6 in the empty room (section 5.6, ~10)")
    for k, v, what in (("opm_noise_xi1", 91, "sigma_N20 OPM xi1/xi0"), ("opm_noise_xi2", 101, "sigma_N20 OPM xi2/xi0"),
                       ("squid_noise_xi1", 79, "sigma_N20 SQUID xi1/xi0"), ("squid_noise_xi2", 59, "sigma_N20 SQUID xi2/xi0"),
                       ("opm_snr_xi1", 68, "SNR OPM xi1/xi0"), ("opm_snr_xi2", 31, "SNR OPM xi2/xi0"),
                       ("squid_snr_xi1", 87, "SNR SQUID xi1/xi0"), ("squid_snr_xi2", 73, "SNR SQUID xi2/xi0")):
        F.add(f"lit_sef_{k}_pct", pct(v), v, s + f"section 5.6, printed ratio {what} (p. 19)")
    F.add("lit_sef_opm_floor_db_range", "29 to 31", [29, 31], s + "section 5.6, Fig. 9 OPM empty-room PSD, dB re fT^2/Hz (figure reading)")
    F.add("lit_sef_opm_floor_asd_range", "28 to 35", [28, 35], s + "section 5.6, Fig. 9 OPM empty-room floor in fT/sqrt(Hz) "
          "(figure reading, if the axis is dB re 1 fT^2/Hz; 'roughly 30')")
    lf = [10 ** (db / 20) for db in (46, 47)]
    F.add("lit_sef_opm_lowfreq_peak_asd_range", f"{lf[0]:.0f} to {lf[1]:.0f}", lf, s + "derived: 10**(dB/20) of the OPM "
          "low-frequency PSD peak, 46-47 dB at ~5 Hz (section 5.6, figure reading), fT/sqrt(Hz)")
    mid = 10 ** ((40 - 34) / 20)
    F.add("lit_sef_midband_noise_ratio", num(mid, 1), mid, s + "derived: 10**((40 - 34)/20), the OPM (40 dB, ~16 Hz) over the "
          "SQUID (34 dB, ~15 Hz) mid-band PSD peaks at xi0 (section 5.6, figure readings; amplitude ratio)")
    F.add("lit_sef_deq_eta4p6_mm", "17.2", 17.2, s + "section 5.6, adult sphere d_eq at eta = 4.6 [derived in the note]")
    F.add("lit_sef_deq_eta4p1_mm", "19.3", 19.3, s + "section 5.6, adult sphere d_eq at eta = 4.1 [derived in the note]")
    F.add("lit_sef_band_range_hz", "4 to 330", [4, 330], s + "section 5.5, 4-Hz high-pass to 330-Hz anti-alias filter")


def lit_study_facts(F, root):
    """Counts of the verified literature table (docs/literature/epilepsy_opm_studies.json, one record per study)."""
    recs = load(root, LITJSON)
    s = LITJSON + " :: derived: "
    n_clin = sum(r["category"] == "clinical" for r in recs)
    n_model = sum(r["category"] == "modelling" for r in recs)
    n_full = sum(r["verification_level"] != "abstract_only" for r in recs)
    F.add("lit_studies_n", count(len(recs)), len(recs), s + "number of records (studies)")
    F.add("lit_studies_clinical_n", count(n_clin), n_clin, s + "records with category 'clinical'")
    F.add("lit_studies_modelling_n", count(n_model), n_model, s + "records with category 'modelling'")
    F.add("lit_studies_full_text_n", count(n_full), n_full, s + "records whose verification_level is not 'abstract_only' (full text "
          "of the version of record, or of the preprint where only the published abstract was accessible)")
    n_vor = sum(r["verification_level"] == "full_text" for r in recs)
    F.add("lit_studies_record_full_text_n", count(n_vor), n_vor, s + "records whose verification_level is 'full_text' (full text of "
          "the version of record)")
    md = (Path(root) / LITMD).read_text()
    m = re.search(r"Searches and reading were done on (\d{4})-(\d{2})-(\d{2})", md)
    if not m:
        raise ValueError(f"{LITMD}: search date not found")
    y, mo, dy = (int(x) for x in m.groups())
    months = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November",
              "December")
    F.add("lit_search_date", f"{dy} {months[mo - 1]} {y}", m.group(0).rsplit(" ", 1)[1],
          f"{LITMD} :: introduction ('Searches and reading were done on {m.group(0).rsplit(' ', 1)[1]}')")


# ------------------------------------------------------------------------------------------------
def g1a_facts(F, root):
    d, f3 = load(root, G1A), load(root, FIG3)
    p = d["parameters"]
    for k, name in (("h_mm", "g1a_h_mm"), ("b_mm", "g1a_b_mm"), ("Q_nAm", "g1a_q_nam"), ("xi_opm_mm", "g1a_xi_opm_mm"),
                    ("xi_squid_mm", "g1a_xi_squid_mm")):
        F.add(name, num(p[k], 0), p[k], f"{G1A} :: parameters.{k}")
    F.add("g1a_sigma_squid_pt", num(p["sigma_squid_pT"], 4), p["sigma_squid_pT"], f"{G1A} :: parameters.sigma_squid_pT (recovered, R-sigma)")
    F.add("g1a_sigma_opm_eta3_pt", num(f3["sigma_opm_pT"], 3), f3["sigma_opm_pT"], f"{FIG3} :: sigma_opm_pT (eta = 3)")
    for eta, v in d["d_eq_mm"].items():
        if v is None:
            continue
        name = f"g1a_deq_eta{tag(float(eta))}_mm"
        F.add(name, num(v, 1), v, f"{G1A} :: d_eq_mm['{eta}'] (exact Eq. 3 root)")
        F.add(name + "_2dp", num(v, 2), v, f"{G1A} :: d_eq_mm['{eta}'] (exact Eq. 3 root)")
    F.add("g1a_deq_eta3_mm_3dp", num(d["d_eq_mm"]["3"], 3), d["d_eq_mm"]["3"], f"{G1A} :: d_eq_mm['3'] (exact Eq. 3 root)")
    none = [float(k) for k, v in d["d_eq_mm"].items() if v is None]
    F.add("g1a_no_crossing_etas", " and ".join(f"{e:g}" for e in none), none, f"{G1A} :: d_eq_mm (null entries: no crossing)")
    lo, hi = d["eta_range_with_crossing"]
    for nd, sfx in ((3, ""), (4, "_4dp")):
        F.add(f"g1a_eta0_exact{sfx}", num(lo, nd), lo, f"{G1A} :: eta_range_with_crossing[0]")
        F.add(f"g1a_eta1_exact{sfx}", num(hi, nd), hi, f"{G1A} :: eta_range_with_crossing[1]")
        F.add(f"g1a_eta_crossing_range{sfx}", f"{num(lo, nd)} to {num(hi, nd)}", [lo, hi], f"{G1A} :: eta_range_with_crossing")
    names = {"Fig. 3 d_eq (eta 3)": "fig3", "Fig. 4B d_eq (eta 2.5)": "fig4b", "Fig. 4C d_eq (eta 4)": "fig4c"}
    for i, item in enumerate(d["printed_vs_exact"]):
        key = f"{G1A} :: printed_vs_exact[{i}]"
        if item["item"] in names:
            n = names[item["item"]]
            F.add(f"g1a_{n}_deq_printed_mm", f"{item['printed_mm']:g}", item["printed_mm"], key + ".printed_mm")
            F.add(f"g1a_{n}_deq_drawn_mm", f"{item['drawn_mm']:g}", item["drawn_mm"], key + ".drawn_mm (marker position in the artwork)")
            F.add(f"g1a_{n}_deq_exact_mm", num(item["exact_mm"], 3 if n == "fig3" else 2), item["exact_mm"], key + ".exact_mm")
        elif item["item"] == "eta0 / eta1":
            F.add("g1a_eta0_printed", f"{item['printed'][0]:g}", item["printed"][0], key + ".printed[0]")
            F.add("g1a_eta1_printed", f"{item['printed'][1]:g}", item["printed"][1], key + ".printed[1]")
        elif item["item"] == "Fig. 6 depths":
            F.add("g1a_fig6_depths_text_mm", ", ".join(f"{x:g}" for x in item["printed_text_mm"]), item["printed_text_mm"],
                  key + ".printed_text_mm (r_Q = 0.4 b, 0.6 b, 0.8 b)")
            F.add("g1a_fig6_depths_caption_mm", ", ".join(f"{x:g}" for x in item["printed_caption_mm"]),
                  item["printed_caption_mm"], key + ".printed_caption_mm (reproduce the curves)")
    nc = d["numerical_checks"]
    for sensor in ("OPM", "SQUID"):
        for k, label in (("max_rel_err_sarvas_2d", "sarvas"), ("max_rel_err_mne_sphere", "mne")):
            F.add(f"g1a_eq1_vs_{label}_rel_err_{sensor.lower()}", sci(nc[sensor][k]), nc[sensor][k],
                  f"{G1A} :: numerical_checks.{sensor}.{k} (maximum relative error of Eq. 1)")
    for xi, x in d["toy_experiment"].items():
        if not xi.startswith("xi_") or not isinstance(x, dict) or "B2/B1" not in x:
            continue
        mm = xi[3:-2]
        F.ratio(f"g1a_toy_b2_b1_xi{mm}mm", x["B2/B1"], f"{G1A} :: toy_experiment.{xi}['B2/B1'] (Fig. 6, deep noise dipole)", three=True)
        F.ratio(f"g1a_toy_b2_b3_xi{mm}mm", x["B2/B3"], f"{G1A} :: toy_experiment.{xi}['B2/B3'] (Fig. 6, superficial noise dipole)", three=True)
    ind = d["d_eq_independent_of_sigma"]
    F.add("g1a_deq_sigma_independent_mm", num(ind["sigma_x1"], 3), [ind[k] for k in ("sigma_x0.5", "sigma_x1", "sigma_x2")],
          f"{G1A} :: d_eq_independent_of_sigma (identical for sigma_SQUID x 0.5, 1, 2)")
    F.add("g1a_fig3_deq_grid_mm", num(f3["d_eq_grid_mm"], 2), f3["d_eq_grid_mm"], f"{FIG3} :: d_eq_grid_mm (raster-matched grid rule, U-J1)")
    F.add("g1a_fig3_deq_paper_caption_mm", f"{f3['d_eq_paper_caption_mm']:g}", f3["d_eq_paper_caption_mm"], f"{FIG3} :: d_eq_paper_caption_mm")
    F.add("g1a_fig3_max_rel_err_forward", sci(f3["max_rel_err_forward"]), f3["max_rel_err_forward"],
          f"{FIG3} :: max_rel_err_forward (MNE sphere forward vs Biot-Savart)")
    F.add("g1a_fig3_max_rel_err_peak", sci(f3["max_rel_err_peak"]), f3["max_rel_err_peak"],
          f"{FIG3} :: max_rel_err_peak (interpolated peak vs Eq. 1)")
    for k in ("OPM", "SQUID"):
        F.add(f"g1a_fig3_signal_brain_surface_{k.lower()}_pt", num(f3["signal_at_brain_surface_pT"][k], 2),
              f3["signal_at_brain_surface_pT"][k], f"{FIG3} :: signal_at_brain_surface_pT.{k} (d = 15 mm)")
        F.add(f"g1a_fig3_snr_brain_surface_{k.lower()}", num(f3["snr_at_brain_surface"][k], 2), f3["snr_at_brain_surface"][k],
              f"{FIG3} :: snr_at_brain_surface.{k} (d = 15 mm, eta = 3)")


# ------------------------------------------------------------------------------------------------
def g1b_facts(F, root):
    d = load(root, G1B)
    cfg = d["config"]
    for name, v, key in (("g1b_n_dipoles", d["dipoles"]["achieved"], "dipoles.achieved"),
                         ("g1b_n_patches", d["patches"]["grown"], "patches.grown"),
                         ("g1b_n_patch_seeds", d["patches"]["seeds"], "patches.seeds"),
                         ("g1b_n_background", d["n_background_dipoles"], "n_background_dipoles"),
                         ("g1b_n_usable_vertices", d["n_usable_vertices"], "n_usable_vertices"),
                         ("g1b_n_valid_vertices", d["n_valid_vertices"], "n_valid_vertices"),
                         ("g1b_dipole_shortfall_bins", d["dipoles"]["shortfall_bins"], "dipoles.shortfall_bins")):
        F.add(name, count(v), v, f"{G1B} :: {key}")
    pt = d["patches"]["total_nAm"]
    F.add("g1b_patch_total_median_nam", num(pt["median"], 0), pt["median"], f"{G1B} :: patches.total_nAm.median")
    F.span("g1b_patch_total_nam", [pt["min"], pt["max"]], lambda x: num(x, 0), f"{G1B} :: patches.total_nAm.min/max")
    F.add("g1b_paper_patch_total_range_nam", "612 to 678", [612, 678], f"{REGISTER} :: HU-patch (totals as printed, pp. 1149-1151)")
    F.add("g1b_patch_area_median_mm2", num(d["patches"]["area_mm2"]["median"], 1), d["patches"]["area_mm2"]["median"],
          f"{G1B} :: patches.area_mm2.median")
    F.add("g1b_patch_target_area_mm2", f"{cfg['sources']['patch_target_area_mm2']:g}", cfg["sources"]["patch_target_area_mm2"],
          f"{G1B} :: config.sources.patch_target_area_mm2")
    F.add("g1b_patch_density_nam_mm2", num(d["patches"]["density_nAm_per_mm2"], 1), d["patches"]["density_nAm_per_mm2"],
          f"{G1B} :: patches.density_nAm_per_mm2")
    F.add("g1b_dipole_peak_nam", num(cfg["sources"]["dipole_peak_nAm"], 0), cfg["sources"]["dipole_peak_nAm"],
          f"{G1B} :: config.sources.dipole_peak_nAm")
    F.add("g1b_background_duration_s", f"{cfg['background']['duration_s']:g}", cfg["background"]["duration_s"],
          f"{G1B} :: config.background.duration_s")
    F.add("g1b_background_share_pct", pct(100 * cfg["background"]["fraction_of_nodes"]), 100 * cfg["background"]["fraction_of_nodes"],
          f"{G1B} :: config.background.fraction_of_nodes")
    F.add("g1b_background_peak_nam", num(cfg["background"]["peak_nAm"], 0), cfg["background"]["peak_nAm"], f"{G1B} :: config.background.peak_nAm")
    F.add("g1b_threshold", f"{cfg['snr']['threshold']:g}", cfg["snr"]["threshold"], f"{G1B} :: config.snr.threshold")
    F.add("g1b_baseline_s", f"{cfg['snr']['baseline_s']:g}", cfg["snr"]["baseline_s"], f"{G1B} :: config.snr.baseline_s")
    F.add("g1b_onset_s", f"{cfg['snr']['onset_s']:g}", cfg["snr"]["onset_s"], f"{G1B} :: config.snr.onset_s")
    F.add("g1b_sfreq_hz", count(cfg["waveform"]["sfreq_hz"]), cfg["waveform"]["sfreq_hz"], f"{G1B} :: config.waveform.sfreq_hz")
    F.add("g1b_bem_triangles", count(cfg["head_model"]["bem_triangles_per_surface"]), cfg["head_model"]["bem_triangles_per_surface"],
          f"{G1B} :: config.head_model.bem_triangles_per_surface (head surface refined to 20,480 triangles, {REGISTER} A-BEM-SKIN)")
    F.add("g1b_extension_band_hz", "0.5 to 70", [0.5, 70], f"{G1B} :: extension_intrinsic_noise.description ('0.5-70 Hz zero-phase Butterworth')")
    F.add("g1b_min_patches_per_bin", str(cfg["descriptors"]["min_patches_per_bin"]), cfg["descriptors"]["min_patches_per_bin"],
          f"{G1B} :: config.descriptors.min_patches_per_bin")
    db = cfg["descriptors"]["depth_bins_mm"]
    F.add("g1b_depth_bins_range_mm", f"{db[0]} to {db[-1]}", [db[0], db[-1]], f"{G1B} :: config.descriptors.depth_bins_mm (5-mm bins)")
    F.add("g1b_conductivity", ", ".join(f"{c:g}" for c in cfg["head_model"]["conductivity_scalp_skull_brain"]),
          cfg["head_model"]["conductivity_scalp_skull_brain"], f"{G1B} :: config.head_model.conductivity_scalp_skull_brain (S/m)")
    cal = d["fig6_calibration"]
    F.add("g1b_fig6_calibration_scale", r2(cal["scale"]), cal["scale"], f"{G1B} :: fig6_calibration.scale (background multiplier)")
    F.add("g1b_fig6_calibration_scale_3dp", r3(cal["scale"]), cal["scale"], f"{G1B} :: fig6_calibration.scale (background multiplier)")
    F.span("g1b_calibration_alternatives", cal["scale_range"], r2, f"{G1B} :: fig6_calibration.scale_range")
    F.add("g1b_calibration_single_draw", r2(cal["scale_single_draw"]), cal["scale_single_draw"], f"{G1B} :: fig6_calibration.scale_single_draw")
    sd = cal["mag"]["realizations"]["single_draw_ratio_p5_p50_p95"]
    F.span("g1b_calibration_single_draw_p5_p95", [sd[0], sd[2]], r2, f"{G1B} :: fig6_calibration.mag.realizations.single_draw_ratio_p5_p50_p95")
    F.add("g1b_calibration_realizations", str(cal["mag"]["realizations"]["n"]), cal["mag"]["realizations"]["n"],
          f"{G1B} :: fig6_calibration.mag.realizations.n")
    # comparison with the digitised paper maps
    for variant, vtag in (("fig6_calibrated", ""), ("as_specified", "asspec_")):
        cw = d["variants"][variant]["comparison_with_paper"]
        for fam in ("dipole", "patch"):
            for kind in ("mag", "grad"):
                c = cw[f"{fam}/p2p"][kind]
                base, key = f"g1b_{vtag}{fam}_{kind}", f"{G1B} :: variants.{variant}.comparison_with_paper['{fam}/p2p'].{kind}"
                F.add(f"{base}_r", r2(c["pearson_r"]), c["pearson_r"], key + ".pearson_r")
                F.add(f"{base}_r_3dp", r3(c["pearson_r"]), c["pearson_r"], key + ".pearson_r")
                F.ratio(f"{base}_level", c["mean_ratio_ours_to_paper"], key + ".mean_ratio_ours_to_paper", three=True)
                F.ratio(f"{base}_level_median", c["median_ratio_ours_to_paper"], key + ".median_ratio_ours_to_paper")
                F.ratio(f"{base}_level_strong", c["mean_ratio_strong_bins"], key + ".mean_ratio_strong_bins (paper SNR >= 2.5)")
                F.ratio(f"{base}_level_weak", c["mean_ratio_weak_bins"], key + ".mean_ratio_weak_bins (paper SNR < 2.5)")
                for k, nm in (("threshold_2p5_agreement", "threshold_agreement"), ("within_paper_class", "within_class"),
                              ("share_ge_2p5_ours", "share_ge2p5_ours"), ("share_ge_2p5_paper", "share_ge2p5_paper")):
                    F.add(f"{base}_{nm}", pct(100 * c[k]), 100 * c[k], f"{key}.{k} (share of bins)")
                    F.add(f"{base}_{nm}_3dp", r3(c[k]), c[k], f"{key}.{k} (fraction of bins)")
                F.add(f"{base}_n_bins", str(c["n_bins"]), c["n_bins"], key + ".n_bins")
        for fam in ("dipole", "patch"):
            vals = [cw[f"{fam}/p2p"][k] for k in ("mag", "grad")]
            F.span(f"g1b_{vtag}{fam}_r", [v["pearson_r"] for v in vals], r2,
                   f"{G1B} :: variants.{variant}.comparison_with_paper['{fam}/p2p'].*.pearson_r")
            F.span(f"g1b_{vtag}{fam}_level", [v["mean_ratio_ours_to_paper"] for v in vals], r2,
                   f"{G1B} :: variants.{variant}.comparison_with_paper['{fam}/p2p'].*.mean_ratio_ours_to_paper")
        every = [cw[f"{fam}/p2p"][k]["mean_ratio_ours_to_paper"] for fam in ("dipole", "patch") for k in ("mag", "grad")]
        F.span(f"g1b_{vtag}level_all", every, r2,
               f"{G1B} :: variants.{variant}.comparison_with_paper['dipole/p2p', 'patch/p2p'].*.mean_ratio_ours_to_paper")
    sens = [v for alt in d["calibration_sensitivity"].values() for k, v in alt.items() if "/p2p/" in k]
    F.span("g1b_calibration_sensitivity_level", sens, r2, f"{G1B} :: calibration_sensitivity.*['dipole|patch/p2p/mag|grad'] "
           "(mean ratio ours/paper under every calibration alternative)")
    for fam, nb in (("dipole", 54), ("patch", 49)):
        s = d["variants"]["fig6_calibrated"]["gm_minus_mm_sign_agreement"][f"{fam}/p2p"]
        F.add(f"g1b_sign_agree_{fam}_n", str(round(s["agree"] * s["n_bins"])), round(s["agree"] * s["n_bins"]),
              f"{G1B} :: derived: variants.fig6_calibrated.gm_minus_mm_sign_agreement['{fam}/p2p'].agree x n_bins")
        F.add(f"g1b_sign_agree_{fam}_bins", str(s["n_bins"]), s["n_bins"],
              f"{G1B} :: variants.fig6_calibrated.gm_minus_mm_sign_agreement['{fam}/p2p'].n_bins")
    sp = d["fig6_spike_comparison"]["dipole superficial"]
    F.add("g1b_fig6_n_superficial", str(sp["n_sources"]), sp["n_sources"], f"{G1B} :: fig6_spike_comparison['dipole superficial'].n_sources")
    for kind, unit_nd in (("mag", 2), ("grad", 0)):
        for k, nm in (("ours_median", "ours"), ("paper_centroid", "paper_centroid"), ("paper_drawn_line", "paper_drawn")):
            F.add(f"g1b_fig6_{kind}_spike_{nm}", num(sp[kind][k], unit_nd), sp[kind][k],
                  f"{G1B} :: fig6_spike_comparison['dipole superficial'].{kind}.{k} ({sp[kind]['unit']})")
    g = sp["grad"]
    F.span("g1b_fig6_grad_spike_over_paper", [g["ours_median"] / g["paper_drawn_line"], g["ours_median"] / g["paper_centroid"]], r2,
           f"{G1B} :: derived: fig6_spike_comparison['dipole superficial'].grad.ours_median / paper_drawn_line and / paper_centroid",
           extra=(("1dp", lambda x: num(x, 1)),))
    hs = d["hilbert_segment_only"]
    for kind in ("mag", "grad"):
        m = hs[f"amplitude_ratio_segment_over_whole/{kind}"]["median"]
        F.add(f"g1b_hilbert_segment_{kind}_pct", pct(100 * (1 - m)), 100 * (1 - m),
              f"{G1B} :: derived: 1 - hilbert_segment_only['amplitude_ratio_segment_over_whole/{kind}'].median")
    bg = d["variants"]["fig6_calibrated"]["background_amplitude_median"]
    for kind, scale, nd in (("mag", 1e12, 2), ("grad", 1e12, 1), ("opm", 1e12, 2)):
        F.add(f"g1b_background_amplitude_{kind}", num(bg[kind] * scale, nd), bg[kind] * scale,
              f"{G1B} :: variants.fig6_calibrated.background_amplitude_median.{kind} (pT; grad pT/m)")
    # OPM columns: count-weighted row means (depth rows) of the OPM over the SQUID SNR maps
    labels = [f"{a}_{b}" for a, b in zip(db[:-1], db[1:])]
    for fam in ("dipole", "patch"):
        b = d["variants"]["fig6_calibrated"]["bins"][f"{fam}/p2p"]
        c = np.array(b["counts"], float)
        m = {k: np.array(v, float) for k, v in b["mean_snr"].items()}
        ok = np.isfinite(m["mag"]) & np.isfinite(m["opm"]) & np.isfinite(m["grad"])
        pre = "g1b_opm" if fam == "dipole" else "g1b_patch_opm"
        for comp in ("mag", "grad"):
            rows = np.where(ok, c * m["opm"], 0).sum(1) / np.where(ok, c * m[comp], 0).sum(1)
            src = (f"{G1B} :: derived: variants.fig6_calibrated.bins['{fam}/p2p'], per depth row sum(counts x mean_snr.opm) / "
                   f"sum(counts x mean_snr.{comp}) (ratio of count-weighted row means; matched 98-site OPM array, {METHODS} section 6)")
            for lab, r in zip(labels, rows):
                F.ratio(f"{pre}_vs_{comp}_row_{lab}", r, src)
            F.span(f"{pre}_vs_{comp}_rows", rows, r2, src)
        sig = d["variants"]["fig6_calibrated"]["significance"][f"{fam}/p2p"]
        for pair in ("opm-mag", "opm-grad", "grad-mag"):
            pv = np.array(sig[pair], float)
            F.add(f"g1b_{fam}_{pair.replace('-', '_')}_bins_p05", str(int(np.nansum(pv < 0.05))), int(np.nansum(pv < 0.05)),
                  f"{G1B} :: derived: count of variants.fig6_calibrated.significance['{fam}/p2p']['{pair}'] < 0.05 (uncorrected per-bin tests)")
            F.add(f"g1b_{fam}_{pair.replace('-', '_')}_bins_tested", str(int(np.sum(np.isfinite(pv)))), int(np.sum(np.isfinite(pv))),
                  f"{G1B} :: derived: finite entries of variants.fig6_calibrated.significance['{fam}/p2p']['{pair}']")
    # extension: white sensor noise (0.5-70 Hz) on every array
    b = d["variants"]["fig6_calibrated"]["bins"]["dipole/p2p"]
    c = np.array(b["counts"], float)
    m0 = {k: np.array(v, float) for k, v in b["mean_snr"].items()}
    base = {comp: (c * m0["opm"]).sum(1) / (c * m0[comp]).sum(1) for comp in ("mag", "grad")}
    ext = d["extension_intrinsic_noise"]["mean_snr"]["fig6_calibrated"]["dipole"]
    moves = {"mag": [], "grad": []}
    for asd in (7, 15, 30):
        o = np.array(ext[f"opm@{asd}"], float)
        for comp, key in (("mag", "mag@3.5"), ("grad", "grad@360")):
            rows = (c * o).sum(1) / (c * np.array(ext[key], float)).sum(1)
            src = (f"{G1B} :: derived: extension_intrinsic_noise.mean_snr.fig6_calibrated.dipole['opm@{asd}'] over ['{key}'], "
                   "ratio of count-weighted row means (white sensor noise, 0.5-70 Hz)")
            for lab, r in zip(labels, rows):
                F.ratio(f"g1b_opm{asd}ft_vs_{comp}_row_{lab}", r, src)
            moves[comp].extend(np.abs(rows - base[comp]))
    for comp in ("mag", "grad"):
        mx = float(max(moves[comp]))
        F.add(f"g1b_opm_vs_{comp}_rows_noise_max_change", num(mx, 2), mx,
              f"{G1B} :: derived: max |row ratio with sensor noise (OPM 7/15/30 fT/sqrt(Hz)) - brain-noise-only row ratio|, dipole/p2p")
        F.add(f"g1b_opm_vs_{comp}_rows_noise_max_change_3dp", num(mx, 3), mx,
              f"{G1B} :: derived: as g1b_opm_vs_{comp}_rows_noise_max_change")


# ------------------------------------------------------------------------------------------------
def g1c_facts(F, root):
    d = load(root, G1C)
    for name, v, key in (("g1c_n_usable_vertices", d["n_usable_vertices"], "n_usable_vertices"),
                         ("g1c_n_valid_vertices", d["n_valid_vertices"], "n_valid_vertices"),
                         ("g1c_n_centroids", d["n_centroids"], "n_centroids"),
                         ("g1c_n_noise_sources", d["n_noise_sources"], "n_noise_sources"),
                         ("g1c_n_mag", d["channels"]["mag"], "channels.mag"), ("g1c_n_grad", d["channels"]["grad"], "channels.grad"),
                         ("g1c_n_pooled", d["channels"]["pooled"], "channels.pooled"),
                         ("g1c_recorded_noise_samples", d["recorded_noise"]["n_samples"], "recorded_noise.n_samples")):
        F.add(name, count(v), v, f"{G1C} :: {key}")
    F.add("g1c_bad_channel", d["channels"]["excluded"][0], d["channels"]["excluded"][0], f"{G1C} :: channels.excluded")
    F.add("g1c_focal_nam", num(d["config"]["sources"]["focal_nAm"], 0), d["config"]["sources"]["focal_nAm"], f"{G1C} :: config.sources.focal_nAm")
    hm = d["config"]["head_model"]["variants"]
    F.add("g1c_skull_as_printed", f"{hm['as_printed'][1]:g}", hm["as_printed"][1], f"{G1C} :: config.head_model.variants.as_printed[1] (S/m)")
    F.add("g1c_skull_probable", f"{hm['probable_default'][1]:g}", hm["probable_default"][1],
          f"{G1C} :: config.head_model.variants.probable_default[1] (S/m)")
    band = d["config"]["recorded_noise"]["band_hz"]
    F.add("g1c_recorded_band_hz", f"{band[0]:g} to {band[1]:g}", band, f"{G1C} :: config.recorded_noise.band_hz")
    F.add("g1c_patch_radii_mm", " and ".join(f"{r:g}" for r in d["config"]["sources"]["patch_radii_mm"]), d["config"]["sources"]["patch_radii_mm"],
          f"{G1C} :: config.sources.patch_radii_mm (geodesic radius)")
    F.add("g1c_noise_grid_mm", "7", 7, f"{G1C} :: config.noise_model.sources (Poisson-disk grid, 7 mm minimum spacing)")
    F.add("g1c_reference_grid_sources", "4,000", 4000, f"{G1C} :: comparison_with_paper.source_sd_nAm_at_4000_sources (key: normalised to the "
          f"~4,000 sources of a 7-mm MNE grid, {METHODS} section 7)")
    F.add("g1c_patch_density_pam_mm2", num(d["config"]["sources"]["patch_density_pAm_per_mm2"], 0),
          d["config"]["sources"]["patch_density_pAm_per_mm2"], f"{G1C} :: config.sources.patch_density_pAm_per_mm2")
    cp = d["comparison_with_paper"]
    F.add("g1c_source_sd_nam", num(cp["source_sd_nAm"]["probable_default"], 2), cp["source_sd_nAm"]["probable_default"],
          f"{G1C} :: comparison_with_paper.source_sd_nAm.probable_default (skull 0.006 S/m, our 7-mm grid)")
    F.add("g1c_source_sd_as_printed_nam", num(cp["source_sd_nAm"]["as_printed"], 2), cp["source_sd_nAm"]["as_printed"],
          f"{G1C} :: comparison_with_paper.source_sd_nAm.as_printed (skull 0.06 S/m)")
    F.add("g1c_source_sd_4000_nam", num(cp["source_sd_nAm_at_4000_sources"]["probable_default"], 2),
          cp["source_sd_nAm_at_4000_sources"]["probable_default"], f"{G1C} :: comparison_with_paper.source_sd_nAm_at_4000_sources.probable_default")
    F.add("g1c_source_sd_4000_as_printed_nam", num(cp["source_sd_nAm_at_4000_sources"]["as_printed"], 2),
          cp["source_sd_nAm_at_4000_sources"]["as_printed"], f"{G1C} :: comparison_with_paper.source_sd_nAm_at_4000_sources.as_printed")
    F.span("g1c_paper_source_sd_nam", cp["paper_source_sd_nAm"], lambda x: num(x, 1), f"{G1C} :: comparison_with_paper.paper_source_sd_nAm")
    F.add("g1c_extension_source_sd_nam", num(cp["extension_source_sd_nAm_brain_only"], 2), cp["extension_source_sd_nAm_brain_only"],
          f"{G1C} :: comparison_with_paper.extension_source_sd_nAm_brain_only")
    # distributions (Eq. 1 SNR, dB) on every usable vertex (focal) or oct-6 centroid (patches)
    for key, x in d["distributions"].items():
        src, noise, chans = key.split("/")
        base = f"g1c_{src}_{noise}_{chans}"
        k = f"{G1C} :: distributions['{key}']"
        if src.startswith("patch16_minus"):
            F.add(f"g1c_patch_gain_{noise}_{chans}_db", signed(x["median_db"]), x["median_db"], k + ".median_db (paired 16 - 10 mm patch SNR)")
            F.add(f"g1c_patch_gain_{noise}_{chans}_area_scaling_db", signed(x["area_scaling_db"]), x["area_scaling_db"],
                  k + ".area_scaling_db (20 log10 of our patch-area ratio)")
            continue
        F.add(f"{base}_median_db", num(x["median_db"], 2), x["median_db"], k + ".median_db")
        F.add(f"{base}_p5_db", num(x["p5_db"], 2), x["p5_db"], k + ".p5_db")
        F.add(f"{base}_p95_db", num(x["p95_db"], 2), x["p95_db"], k + ".p95_db")
        sk = x["skull_0p06_minus_0p006_db"]
        F.add(f"{base}_skull_change_db", signed(sk["median"]), sk["median"], k + ".skull_0p06_minus_0p006_db.median (0.06 minus 0.006 S/m)")
        F.add(f"{base}_skull_change_p95abs_db", num(sk["p95_abs"], 2), sk["p95_abs"], k + ".skull_0p06_minus_0p006_db.p95_abs")
        F.add(f"{base}_skull_change_p95abs_db_1dp", num(sk["p95_abs"], 1), sk["p95_abs"], k + ".skull_0p06_minus_0p006_db.p95_abs")
    fp = d["distributions"]["focal/model/pooled"]
    F.add("g1c_focal_model_pooled_p5_p95_range_1dp", f"{num(fp['p5_db'], 1)} to {num(fp['p95_db'], 1)}", [fp["p5_db"], fp["p95_db"]],
          f"{G1C} :: distributions['focal/model/pooled'].p5_db, p95_db")
    F.add("g1c_focal_model_pooled_centroids_median_db", num(fp["oct6_centroids_median_db"], 2), fp["oct6_centroids_median_db"],
          f"{G1C} :: distributions['focal/model/pooled'].oct6_centroids_median_db")
    for chans in ("mag", "grad", "pooled"):
        c = cp[f"probable_default/focal/model/{chans}"]
        k = f"{G1C} :: comparison_with_paper['probable_default/focal/model/{chans}']"
        for s in ("share_below_m29", "share_m29_to_m19", "share_above_m19", "share_below_m29_deep_medial_regions"):
            F.add(f"g1c_focal_{chans}_{s.replace('share_', '').replace('_regions', '')}_pct", pct(100 * c[s]), 100 * c[s], f"{k}.{s}")
        for lobe, v in c["median_by_lobe"].items():
            F.add(f"g1c_focal_{chans}_{lobe}_median_db", num(v, 2), v, f"{k}.median_by_lobe.{lobe}")
        g = cp[f"probable_default/patch16_minus_patch10/model/{chans}"]
        k = f"{G1C} :: comparison_with_paper['probable_default/patch16_minus_patch10/model/{chans}']"
        F.add(f"g1c_patch_gain_mesial_temporal_{chans}_db", signed(g["mesial_temporal_median_db"]), g["mesial_temporal_median_db"],
              k + ".mesial_temporal_median_db")
        for lobe, v in g["by_lobe"].items():
            F.add(f"g1c_patch_gain_{lobe}_{chans}_db", signed(v), v, f"{k}.by_lobe.{lobe}")
    g = cp["probable_default/patch16_minus_patch10/model/pooled"]
    F.add("g1c_patch_gain_mesial_temporal_n", str(g["n_mesial_temporal"]), g["n_mesial_temporal"],
          f"{G1C} :: comparison_with_paper['probable_default/patch16_minus_patch10/model/pooled'].n_mesial_temporal")
    F.add("g1c_patch_gain_nominal_area_scaling_db", signed(g["nominal_area_scaling_db"]), g["nominal_area_scaling_db"],
          f"{G1C} :: comparison_with_paper['probable_default/patch16_minus_patch10/model/pooled'].nominal_area_scaling_db (3 vs 8 cm^2)")
    F.add("g1c_paper_patch_gain_mesial_temporal_db", signed(g["paper_mesial_temporal_db"], 0), g["paper_mesial_temporal_db"],
          f"{G1C} :: comparison_with_paper['probable_default/patch16_minus_patch10/model/pooled'].paper_mesial_temporal_db")
    shares = [cp[f"probable_default/focal/model/{c}"]["share_m29_to_m19"] * 100 for c in ("mag", "grad", "pooled")]
    F.span("g1c_focal_in_paper_range_pct", shares, pct, f"{G1C} :: comparison_with_paper['probable_default/focal/model/*'].share_m29_to_m19")
    for r, x in d["patch_area_cm2"].items():
        F.add(f"g1c_patch{r}mm_area_cm2", num(x["median"], 2), x["median"], f"{G1C} :: patch_area_cm2['{r}'].median")
        F.add(f"g1c_patch{r}mm_flat_disc_cm2", num(x["flat_disc"], 2), x["flat_disc"], f"{G1C} :: patch_area_cm2['{r}'].flat_disc")
    for r, x in d["patch_truncation"].items():
        k = f"{G1C} :: patch_truncation['{r}']"
        F.add(f"g1c_patch{r}mm_lose5pct_pct", pct(100 * x["share_losing_over_5pct"]), 100 * x["share_losing_over_5pct"],
              k + ".share_losing_over_5pct")
        F.add(f"g1c_patch{r}mm_lose20pct_pct", pct(100 * x["share_losing_over_20pct"]), 100 * x["share_losing_over_20pct"],
              k + ".share_losing_over_20pct")
        F.add(f"g1c_patch{r}mm_truncation_snr_change_db", signed(x["snr_change_db"]["median"]), x["snr_change_db"]["median"],
              k + ".snr_change_db.median")
        F.add(f"g1c_patch{r}mm_truncation_snr_change_p5_db", signed(x["snr_change_db"]["p5"], 1), x["snr_change_db"]["p5"], k + ".snr_change_db.p5")
        F.add(f"g1c_patch{r}mm_truncation_over1db_pct", pct(100 * x["snr_change_db"]["share_abs_over_1db"]),
              100 * x["snr_change_db"]["share_abs_over_1db"], k + ".snr_change_db.share_abs_over_1db")
    pt = d["patch_truncation"]
    F.span("g1c_truncation_lose5pct", [100 * pt[r]["share_losing_over_5pct"] for r in pt], pct,
           f"{G1C} :: patch_truncation['10', '16'].share_losing_over_5pct")
    F.span("g1c_truncation_over1db", [100 * pt[r]["snr_change_db"]["share_abs_over_1db"] for r in pt], pct,
           f"{G1C} :: patch_truncation['10', '16'].snr_change_db.share_abs_over_1db")
    w = d["without_medial_wall"]
    F.add("g1c_wall_share_pct", pct(100 * w["share_of_usable_vertices"]), 100 * w["share_of_usable_vertices"],
          f"{G1C} :: without_medial_wall.share_of_usable_vertices")
    F.add("g1c_wall_centroids", count(w["centroids_on_wall"]), w["centroids_on_wall"], f"{G1C} :: without_medial_wall.centroids_on_wall")
    F.add("g1c_wall_noise_sources", count(w["noise_sources_on_wall"]), w["noise_sources_on_wall"],
          f"{G1C} :: without_medial_wall.noise_sources_on_wall")
    F.add("g1c_nowall_focal_pooled_median_db", num(w["focal_model_pooled"]["median_db"], 2), w["focal_model_pooled"]["median_db"],
          f"{G1C} :: without_medial_wall.focal_model_pooled.median_db")
    F.add("g1c_nowall_focal_pooled_in_paper_range_pct", pct(100 * w["focal_model_pooled"]["share_in_paper_range"]),
          100 * w["focal_model_pooled"]["share_in_paper_range"], f"{G1C} :: without_medial_wall.focal_model_pooled.share_in_paper_range")
    F.add("g1c_wall_focal_pooled_in_paper_range_pct", pct(100 * w["with_wall_focal_model_pooled"]["share_in_paper_range"]),
          100 * w["with_wall_focal_model_pooled"]["share_in_paper_range"], f"{G1C} :: without_medial_wall.with_wall_focal_model_pooled.share_in_paper_range")
    v = d["variant_no_wall_untruncated"]
    F.add("g1c_untruncated_n_centroids", count(v["n_centroids"]), v["n_centroids"], f"{G1C} :: variant_no_wall_untruncated.n_centroids")
    for key, x in v.items():
        if isinstance(x, dict) and "median_db" in x:
            src, noise, chans = key.split("/")
            k = f"{G1C} :: variant_no_wall_untruncated['{key}']"
            F.add(f"g1c_untruncated_{src}_{chans}_median_db", num(x["median_db"], 2), x["median_db"], k + ".median_db")
            F.add(f"g1c_untruncated_{src}_{chans}_in_paper_range_pct", pct(100 * x["share_in_paper_range"]), 100 * x["share_in_paper_range"],
                  k + ".share_in_paper_range")
    rn = d["recorded_noise"]["empty_room_over_recorded_variance"]
    for kind in ("mag", "grad"):
        F.add(f"g1c_empty_room_share_{kind}_pct", pct(100 * rn[kind]), 100 * rn[kind],
              f"{G1C} :: recorded_noise.empty_room_over_recorded_variance.{kind}")
    ld = d["lead_field_diagnostics"]
    F.add("g1c_leadfield_energy_ratio_max", num(ld["usable_columns_energy_ratio_max"], 1), ld["usable_columns_energy_ratio_max"],
          f"{G1C} :: lead_field_diagnostics.usable_columns_energy_ratio_max")
    F.span("g1c_noise_top_source_share", [100 * ld["noise_grid_top_source_power_share"][k] for k in ("mag", "grad")], pct,
           f"{G1C} :: lead_field_diagnostics.noise_grid_top_source_power_share.mag/grad")
    # OPM extension: paired per-vertex medians of OPM minus SQUID Eq. 1 SNR
    e = d["extension_opm"]
    for src in ("focal", "patch10", "patch16"):
        for kind in ("mag", "grad"):
            x = e[f"brain_only/{src}/opm_minus_{kind}"]
            F.add(f"g1c_opm_minus_{kind}_{src}_brain_only_db", signed(x["median_db"]), x["median_db"],
                  f"{G1C} :: extension_opm['brain_only/{src}/opm_minus_{kind}'].median_db")
            F.add(f"g1c_opm_minus_{kind}_{src}_brain_only_share_pct", pct(100 * x["share_opm_better"]), 100 * x["share_opm_better"],
                  f"{G1C} :: extension_opm['brain_only/{src}/opm_minus_{kind}'].share_opm_better")
            vals = []
            for asd in (7, 10, 15, 20, 30):
                x = e[f"{src}/opm@{asd}_minus_{kind}"]
                vals.append(x["median_db"])
                F.add(f"g1c_opm{asd}ft_minus_{kind}_{src}_db", signed(x["median_db"]), x["median_db"],
                      f"{G1C} :: extension_opm['{src}/opm@{asd}_minus_{kind}'].median_db")
                F.add(f"g1c_opm{asd}ft_minus_{kind}_{src}_share_pct", pct(100 * x["share_opm_better"]), 100 * x["share_opm_better"],
                      f"{G1C} :: extension_opm['{src}/opm@{asd}_minus_{kind}'].share_opm_better")
            F.span(f"g1c_opm_minus_{kind}_{src}_noise_sweep_db", vals, signed,
                   f"{G1C} :: extension_opm['{src}/opm@7..30_minus_{kind}'].median_db (OPM 7-30 fT/sqrt(Hz))")
    a, b = e["focal/opm@10_minus_mag"]["median_db"], e["focal/opm@15_minus_mag"]["median_db"]
    cross = 10 + 5 * a / (a - b)
    F.add("g1c_opm_mag_crossing_asd", num(cross, 0), cross, f"{G1C} :: derived: linear interpolation of extension_opm['focal/opm@10_minus_mag'] "
          "and ['focal/opm@15_minus_mag'].median_db to 0 dB (fT/sqrt(Hz))")


# ------------------------------------------------------------------------------------------------
def g2_config_facts(F, d):
    c = d["config"]
    k = f"{G2} :: config"
    F.add("g2_n_targets", count(d["n_targets"]), d["n_targets"], f"{G2} :: n_targets")
    F.add("g2_n_background", count(d["n_background_grid"]), d["n_background_grid"], f"{G2} :: n_background_grid")
    F.add("g2_n_parcels", "70", 70, f"{G2} :: primary.oracle[*].ci_method ('parcels (70)': 68 Desikan-Killiany parcels + 2 medial-wall labels)")
    F.add("g2_n_convergence_subset", count(c["convergence"]["subset_targets"]), c["convergence"]["subset_targets"], k + ".convergence.subset_targets")
    F.add("g2_band_lo_hz", f"{c['band']['l_freq_hz']:g}", c["band"]["l_freq_hz"], k + ".band.l_freq_hz")
    F.add("g2_band_hi_hz", f"{c['band']['h_freq_hz']:g}", c["band"]["h_freq_hz"], k + ".band.h_freq_hz")
    F.add("g2_band_range_hz", f"{c['band']['l_freq_hz']:g} to {c['band']['h_freq_hz']:g}", [c["band"]["l_freq_hz"], c["band"]["h_freq_hz"]],
          k + ".band")
    F.add("g2_filter_order", str(c["band"]["order"]), c["band"]["order"], k + ".band.order (Butterworth, zero phase)")
    F.add("g2_enbw_hz", num(d["enbw_hz"], 1), d["enbw_hz"], f"{G2} :: enbw_hz")
    F.add("g2_squid_mag_asd", f"{c['sensors']['squid_asd']['mag_fT_per_rtHz']:g}", c["sensors"]["squid_asd"]["mag_fT_per_rtHz"],
          k + ".sensors.squid_asd.mag_fT_per_rtHz (fT/sqrt(Hz), HW-mag-noise)")
    F.add("g2_squid_grad_asd", f"{c['sensors']['squid_asd']['grad_fT_per_cm_rtHz']:g}", c["sensors"]["squid_asd"]["grad_fT_per_cm_rtHz"],
          k + ".sensors.squid_asd.grad_fT_per_cm_rtHz (fT/(cm sqrt(Hz)), HW-grad-noise)")
    F.add("g2_opm_asd_primary", f"{c['sensors']['opm_asd_primary_fT_per_rtHz']:g}", c["sensors"]["opm_asd_primary_fT_per_rtHz"],
          k + ".sensors.opm_asd_primary_fT_per_rtHz (fT/sqrt(Hz), A-G2-OPMNOISE)")
    sweep = c["sensors"]["opm_asd_fT_per_rtHz"]
    F.add("g2_opm_asd_sweep", ", ".join(f"{x:g}" for x in sweep), sweep, k + ".sensors.opm_asd_fT_per_rtHz")
    F.span("g2_opm_asd_sweep", sweep, lambda x: f"{x:g}", k + ".sensors.opm_asd_fT_per_rtHz")
    for x in sweep:  # each level of the sweep, as named in the text (fT/sqrt(Hz))
        F.add(f"g2_opm_asd_sweep_{tag(x)}", f"{x:g}", x, k + f".sensors.opm_asd_fT_per_rtHz (the level {x:g})")
    gaps = [g for g in c["sensors"]["opm_scalp_gap_mm"] if g > 0]
    F.add("g2_opm_gap_sweep_mm", " and ".join(f"{g:g}" for g in gaps), gaps, k + ".sensors.opm_scalp_gap_mm (0 = primary)")
    F.add("g2_opm_gap_small_mm", f"{min(gaps):g}", min(gaps), k + ".sensors.opm_scalp_gap_mm (the smaller non-zero scalp gap)")
    F.add("g2_opm_gap_large_mm", f"{max(gaps):g}", max(gaps), k + ".sensors.opm_scalp_gap_mm (the larger scalp gap)")
    F.add("g2_ci_level_pct", pct(95), 95, f"{METHODS} :: section 8 ('95 % CIs from a bootstrap over 70 groups'); "
          f"{G2} :: primary.oracle[*].ci95 (2.5th and 97.5th percentiles of the bootstrap)")
    F.add("g2_bad_channel", c["sensors"]["bads"][0], c["sensors"]["bads"][0], k + ".sensors.bads")
    F.add("g2_focal_nam", num(c["sources"]["focal_nAm"], 0), c["sources"]["focal_nAm"], k + ".sources.focal_nAm")
    F.add("g2_patch_radii_mm", ", ".join(f"{x:g}" for x in c["sources"]["patch_radii_mm"]), c["sources"]["patch_radii_mm"],
          k + ".sources.patch_radii_mm")
    F.add("g2_patch_total_nam", num(c["sources"]["patch_total_nAm"], 0), c["sources"]["patch_total_nAm"], k + ".sources.patch_total_nAm")
    F.add("g2_patch_density_nam_mm2", f"{c['sources']['patch_density_nAm_per_mm2']:g}", c["sources"]["patch_density_nAm_per_mm2"],
          k + ".sources.patch_density_nAm_per_mm2")
    F.add("g2_bg_grid_mm", f"{c['background']['grid_spacing_mm']:g}", c["background"]["grid_spacing_mm"], k + ".background.grid_spacing_mm")
    corr = [x for x in c["background"]["correlation_lengths_mm"] if x > 0]
    F.add("g2_bg_corr_lengths_mm", " and ".join(f"{x:g}" for x in corr), corr, k + ".background.correlation_lengths_mm")
    F.add("g2_bg_calibration_rule", c["background"]["calibration"], c["background"]["calibration"],
          k + ".background.calibration (the one fitted brain-noise level)")
    F.add("g2_room_terms", "8", 8, f"{CFG} :: environment (8-term external expansion: homogeneous + linear gradient; register A-G2-ENV)")
    F.add("g2_room_origin_mm", ", ".join(f"{1e3 * x:g}" for x in c["environment"]["r0_head_m"]), [1e3 * x for x in c["environment"]["r0_head_m"]],
          k + ".environment.r0_head_m (head frame, mm)")
    F.add("g2_cov_seconds", " and ".join(f"{x:g}" for x in c["covariance"]["estimate_seconds"]), c["covariance"]["estimate_seconds"],
          k + ".covariance.estimate_seconds (plug-in Ledoit-Wolf)")
    for lab, s in (("T10", "10s"), ("T60", "60s")):
        n = d["n_estimate_samples"][lab]
        F.add(f"g2_plugin{s}_samples", count(n), n, f"{G2} :: n_estimate_samples.{lab} (independent Gaussian samples, 2 x 39 Hz x T, "
              "the 2BT equivalent of 10/60 s; not covariance estimates from coloured time series)")
    hp = c["head_position"]
    F.add("g2_head_translation_mm", f"{hp['translation_mm']:g}", hp["translation_mm"],
          k + ".head_position.translation_mm (+/- along each device axis)")
    F.add("g2_head_pitch_deg", f"{hp['pitch_deg']:g}", hp["pitch_deg"], k + ".head_position.pitch_deg (+/-)")
    F.add("g2_head_fit_clearance_mm", f"{hp['fit_clearance_mm']:g}", hp["fit_clearance_mm"], k + ".head_position.fit_clearance_mm ('well fitted')")
    F.add("g2_dewar_spacing_mm", f"{hp['dewar_spacing_mm']:g}", hp["dewar_spacing_mm"], k + ".head_position.dewar_spacing_mm (HW-18mm)")
    F.add("g2_conductivity", ", ".join(f"{x:g}" for x in c["anatomy"]["bem_conductivity"]), c["anatomy"]["bem_conductivity"],
          k + ".anatomy.bem_conductivity (S/m, scalp/skull/brain)")
    for name, v, txt, where in (("g2_opm_cell_mm", 10, "10", "A-OPM-CELL (10-mm cubic sensing volume, 27 integration points)"),
                                ("g2_opm_axis_radius_mm", 15, "15", "A-OPM-AXIS (head-surface normal averaged within 15 mm)"),
                                ("g2_cell_clearance_rule_mm", 1, "1", "A-OPM-CLEAR (cell >= 1 mm outside the BEM head surface)"),
                                ("g2_centre_clearance_scalp_mm", 6, "6", "A-OPM-CLEAR (sensing centre >= standoff - 1 mm from the MRI scalp)"),
                                ("g2_centre_clearance_bem_mm", 4, "4", "A-OPM-CLEAR (sensing centre >= 4 mm from the BEM head surface)"),
                                ("g2_preauricular_exclusion_mm", 20, "20", "A-OPM-COVER (no site within 20 mm of the preauricular points)"),
                                ("g2_dense_packing_mm", 17, "17", "U-OPM-PACK (minimum sensing-centre spacing of the dense arrays)"),
                                ("g2_source_skull_clearance_mm", 4, "4", "A-BEM-DIST (usable sources >= 4 mm from the inner skull)"),
                                ("g2_head_triangles_coarse", 5120, "5,120", "A-BEM-SKIN (BEM head surface before subdivision)"),
                                ("g2_head_triangles", 20480, "20,480", "A-BEM-SKIN (head surface subdivided once in every 3-layer model)"),
                                ("g2_plugin_bandwidth_hz", 39, "39", "A-G2-COVEST (2 x 39 Hz x T independent samples)")):
        F.add(name, txt, v, f"{REGISTER} :: {where}")
    F.add("g2_usable_vertices_pct", pct(91.3), 91.3, f"{REGISTER} :: A-BEM-DIST (91.3 % of the valid vertices are usable)")
    F.add("g2_usable_vertices_pct_1dp", "91.3%", 91.3, f"{REGISTER} :: A-BEM-DIST (91.3 % of the valid vertices are usable)")
    F.add("g2_runtime_s", count(d["runtime_s"]), d["runtime_s"], f"{G2} :: runtime_s")
    F.add("g2_runtime_min", num(d["runtime_s"] / 60, 0), d["runtime_s"] / 60, f"{G2} :: derived: runtime_s / 60")
    F.add("g2_commit", d["provenance"]["commit"], d["provenance"]["commit"], f"{G2} :: provenance.commit")
    # arrays, standoffs and the conformed head surface
    a = d["arrays"]
    F.add("g2_squid_sites", str(a["squid"]["sites"]), a["squid"]["sites"], f"{G2} :: arrays.squid.sites")
    F.add("g2_squid_channels", str(a["squid"]["channels"]), a["squid"]["channels"], f"{G2} :: arrays.squid.channels")
    for arr, short in ARRAYS.items():
        x = a[arr]
        F.add(f"g2_{short}_sites", str(x["n_sites"]), x["n_sites"], f"{G2} :: arrays.{arr}.n_sites")
        F.add(f"g2_{short}_standoff_mm", f"{x['standoff_mm']:g}", x["standoff_mm"],
              f"{G2} :: arrays.{arr}.standoff_mm (sensing centre from the helmet inner surface)")
        F.add(f"g2_{short}_cell_clearance_mm", num(x["exact_cell_clearance_min_mm"], 1), x["exact_cell_clearance_min_mm"],
              f"{G2} :: arrays.{arr}.exact_cell_clearance_min_mm (whole 10-mm cell outside the BEM head surface)")
        F.add(f"g2_{short}_cell_clearance_mm_3dp", num(x["exact_cell_clearance_min_mm"], 3), x["exact_cell_clearance_min_mm"],
              f"{G2} :: arrays.{arr}.exact_cell_clearance_min_mm")
        F.add(f"g2_{short}_min_spacing_mm", num(x["min_spacing_mm"], 1), x["min_spacing_mm"], f"{G2} :: arrays.{arr}.min_spacing_mm")
        if "median_spacing_mm" in x:
            F.add(f"g2_{short}_median_spacing_mm", num(x["median_spacing_mm"], 1), x["median_spacing_mm"], f"{G2} :: arrays.{arr}.median_spacing_mm")
        moved = x.get("n_moved_out", x.get("parent_n_moved_out"))
        F.add(f"g2_{short}_moved_out", str(moved), moved, f"{G2} :: arrays.{arr}.n_moved_out (sites moved outward for clearance)")
    F.add("g2_matched_excluded_sites", str(len(a["opm_matched"]["excluded_sites"])), len(a["opm_matched"]["excluded_sites"]),
          f"{G2} :: arrays.opm_matched.excluded_sites (Neuromag sites outside the OPM coverage rule)")
    F.add("g2_dense_max_shift_mm", num(a["opm_dense"]["max_extra_shift_mm"], 1), a["opm_dense"]["max_extra_shift_mm"],
          f"{G2} :: arrays.opm_dense.max_extra_shift_mm")
    h = a["head_surface_conform"]
    for key, nm in (("moved_median_mm", "median"), ("moved_p90_mm", "p90"), ("moved_max_mm", "max"), ("outer_skull_min_mm", "outer_skull_min")):
        F.add(f"g2_conform_{nm}_mm", num(h[key], 1), h[key], f"{G2} :: arrays.head_surface_conform.{key}")
        F.add(f"g2_conform_{nm}_mm_2dp", num(h[key], 2), h[key], f"{G2} :: arrays.head_surface_conform.{key}")
    F.add("g2_conform_vertices", count(h["n_vertices"]), h["n_vertices"], f"{G2} :: arrays.head_surface_conform.n_vertices")
    for arr, x in d["bridge_to_sphere"]["sensor_distance_mm"].items():
        short = ARRAYS.get(arr, arr)
        for stat in ("median", "p5", "p95"):
            F.add(f"g2_{short}_scalp_dist_{stat}_mm", num(x[stat], 1), x[stat],
                  f"{G2} :: bridge_to_sphere.sensor_distance_mm.{arr}.{stat} (nearest MRI-scalp distance; SQUID: magnetometer coils)")
        F.add(f"g2_{short}_scalp_dist_median_mm_2dp", num(x["median"], 2), x["median"], f"{G2} :: bridge_to_sphere.sensor_distance_mm.{arr}.median")
    F.span("g2_squid_scalp_dist_p5_p95_mm", [d["bridge_to_sphere"]["sensor_distance_mm"]["squid"][s] for s in ("p5", "p95")],
           lambda x: num(x, 1), f"{G2} :: bridge_to_sphere.sensor_distance_mm.squid.p5/p95")
    F.add("g2_squid_normal_dist_range_mm", "25.1 to 41.4", [25.1, 41.4],
          f"{METHODS} :: section 2.1 (magnetometer coil to MRI scalp along the inward normal)")
    F.add("g2_squid_normal_dist_median_mm", "30.9", 30.9, f"{METHODS} :: section 2.1 (magnetometer coil to MRI scalp along the inward normal)")
    for pos, x in d["head_positions"].items():
        nm = pos.replace("+", "p").replace("-", "m").lower()
        F.add(f"g2_headpos_{nm}_min_dist_mm", num(x["min_dist_mm"], 1), x["min_dist_mm"],
              f"{G2} :: head_positions['{pos}'].min_dist_mm (scalp to nearest magnetometer)")
        F.add(f"g2_headpos_{nm}_median_dist_mm", num(x["median_dist_mm"], 1), x["median_dist_mm"], f"{G2} :: head_positions['{pos}'].median_dist_mm")
    t = d["triaxial_control"]
    F.add("g2_triax_sites", str(t["sites"]), t["sites"], f"{G2} :: triaxial_control.sites")
    F.add("g2_triax_channels", str(t["channels"]), t["channels"], f"{G2} :: triaxial_control.channels")
    F.add("g2_triax_asd", f"{t['opm_asd_fT']:.0f}", t["opm_asd_fT"], f"{G2} :: triaxial_control.opm_asd_fT (fT/sqrt(Hz), every axis)")
    for k2, r in d["retained_rank"].items():
        arr, cond = k2.split("/")
        if cond in ("projected", "intrinsic"):
            nm = ARRAYS.get(arr, arr)
            F.add(f"g2_rank_{nm}_{CONDS[cond]}", str(r), r, f"{G2} :: retained_rank['{k2}']")


def g2_noise_facts(F, d, band):
    nv = d["noise_validation"]
    m, mo = nv["measured"], nv["model"]
    k = f"{G2} :: noise_validation"
    for key, nm, nd in (("empty_room_rms_mag_fT", "empty_room_mag_ft", 1), ("empty_room_rms_grad_fT_cm", "empty_room_grad_ftcm", 1),
                        ("baseline_rms_mag_fT", "baseline_mag_ft", 0), ("baseline_rms_grad_fT_cm", "baseline_grad_ftcm", 1),
                        ("brain_rms_mag_fT", "brain_mag_ft", 0), ("brain_rms_grad_fT_cm", "brain_grad_ftcm", 1)):
        F.add(f"g2_measured_{nm}", num(m[key], nd), m[key], f"{k}.measured.{key}")
    for key, nm, nd in (("empty_room_rms_mag_fT", "empty_room_mag_ft", 1), ("empty_room_rms_grad_fT_cm", "empty_room_grad_ftcm", 1),
                        ("brain_rms_mag_fT", "brain_mag_ft", 0), ("brain_rms_grad_fT_cm", "brain_grad_ftcm", 1),
                        ("intrinsic_rms_mag_fT", "intrinsic_mag_ft", 1), ("intrinsic_rms_grad_fT_cm", "intrinsic_grad_ftcm", 1)):
        F.add(f"g2_model_{nm}", num(mo[key], nd), mo[key], f"{k}.model.{key}")
    F.add("g2_model_brain_grad_fitted_ftcm", num(mo["brain_rms_grad_fT_cm"], 1), mo["brain_rms_grad_fT_cm"],
          f"{k}.model.brain_rms_grad_fT_cm (fitted: equal to the measured level by the calibration, A-G2-BRAIN)")
    for asd, v in mo["intrinsic_rms_opm_fT"].items():
        F.add(f"g2_model_intrinsic_opm{asd}ft_ft", num(v, 1), v, f"{k}.model.intrinsic_rms_opm_fT['{asd}'] (in-band RMS, fT)")
    F.ratio("g2_brain_mag_model_over_measured", mo["brain_mag_model_over_measured"], f"{k}.model.brain_mag_model_over_measured "
            "(the one independent prediction: magnetometer brain-noise RMS)", three=True)
    under = 100 * (1 - mo["brain_mag_model_over_measured"])
    F.add("g2_brain_mag_underprediction_pct", pct(under), under, f"{G2} :: derived: 1 - noise_validation.model.brain_mag_model_over_measured")
    for kind in ("mag", "grad"):
        F.add(f"g2_room_explained_{kind}_pct", pct(100 * nv["environment_explained_fraction"][kind]), 100 * nv["environment_explained_fraction"][kind],
              f"{k}.environment_explained_fraction.{kind} (fitted room field, empty-room variance)")
        p = nv["per_channel_brain_ratio_model_over_measured"][kind]
        src = f"{k}.per_channel_brain_ratio_model_over_measured.{kind} (per-channel modelled / measured brain-noise variance)"
        for stat in ("p5", "median", "p95"):
            F.add(f"g2_brain_var_ratio_{kind}_{stat}", r2(p[stat]), p[stat], f"{src}.{stat}")
        F.add(f"g2_brain_var_ratio_{kind}_p5_p95_range", f"{r2(p['p5'])} to {r2(p['p95'])}", [p["p5"], p["p95"]], src)
    F.add("g2_negative_brain_variance_channels", str(m["negative_brain_variance_channels"]), m["negative_brain_variance_channels"],
          f"{k}.measured.negative_brain_variance_channels")
    for key, nm in (("n_baseline_windows", "baseline_windows"), ("n_baseline_samples", "baseline_samples"), ("n_empty_room_samples", "empty_room_samples")):
        F.add(f"g2_measured_{nm}", count(m[key]), m[key], f"{k}.measured.{key}")
    F.add("g2_noise_fitted", "the gradiometer brain-noise level (median good-gradiometer variance) and the room field (empty-room fit)",
          "fitted", f"{REGISTER} :: A-G2-BRAIN, A-G2-ENV; {G2} :: config.background.calibration, config.environment.fit")
    F.add("g2_noise_predicted", "the magnetometer brain-noise level; nothing validates the OPM covariance", "predicted",
          f"{REGISTER} :: A-G2-BRAIN ('the magnetometer level is a validation, not a fit'); {G2} :: noise_validation")
    nc = d["noise_composition"]
    for arr, x in nc.items():
        nm = ARRAYS.get(arr, arr)
        grad = arr == "squid_grad"
        scale, nd, unit = (1e13, 1, "fT/cm") if grad else (1e15, 0, "fT")
        for term in ("intrinsic", "brain", "env"):
            v = x[f"{term}_rms"] * scale
            F.add(f"g2_noise_{nm}_{'room' if term == 'env' else term}", num(v, nd), v,
                  f"{G2} :: noise_composition.{arr}.{term}_rms ({unit}, in-band RMS)")
    tot = {a: np.sqrt(sum(nc[a][f"{t}_rms"] ** 2 for t in ("intrinsic", "brain", "env"))) for a in nc}
    totnr = {a: np.sqrt(sum(nc[a][f"{t}_rms"] ** 2 for t in ("intrinsic", "brain"))) for a in nc}
    for arr in ("opm_dense", "opm_matched"):
        nm = ARRAYS[arr]
        r, rn = tot[arr] / tot["squid_mag"], totnr[arr] / totnr["squid_mag"]
        F.ratio(f"g2_noise_rms_{nm}_vs_mag", r, f"{G2} :: derived: sqrt of the summed noise_composition.{arr} variances (intrinsic, brain, "
                "env) over the same for squid_mag (single-channel noise RMS ratio)")
        F.add(f"g2_noise_rms_{nm}_vs_mag_ratio_1dp", num(r, 1), r, f"{G2} :: derived: as g2_noise_rms_{nm}_vs_mag_ratio")
        F.ratio(f"g2_noise_rms_{nm}_vs_mag_noroom", rn, f"{G2} :: derived: as g2_noise_rms_{nm}_vs_mag_ratio without the env term")
        F.add(f"g2_noise_rms_{nm}_vs_mag_noroom_ratio_1dp", num(rn, 1), rn, f"{G2} :: derived: as g2_noise_rms_{nm}_vs_mag_noroom_ratio")
    F.add("g2_noise_rms_dense_vs_mag_range_1dp",
          f"{num(tot['opm_dense'] / tot['squid_mag'], 1)} to {num(totnr['opm_dense'] / totnr['squid_mag'], 1)}",
          [tot["opm_dense"] / tot["squid_mag"], totnr["opm_dense"] / totnr["squid_mag"]], f"{G2} :: derived: with and without the room field (as above)")
    F.add("g2_brain_noise_density", num(nv["brain_source_density_nAm_per_sqrt_mm2"], 2), nv["brain_source_density_nAm_per_sqrt_mm2"],
          f"{k}.brain_source_density_nAm_per_sqrt_mm2 (nAm/sqrt(mm^2))")
    # noise regime per band (g2_band_sensitivity.json)
    for b, x in band["bands"].items():
        bt = b.lower().replace("-", "_")
        kb = f"{BAND} :: bands['{b}']"
        lo, hi = b.lower().replace("hz", "").split("-")
        F.add(f"g2_band_{bt}_label", f"{lo}–{hi}", [float(lo), float(hi)], f"{kb} (the band's key: {lo} to {hi} Hz)")
        for arr, v in x["intrinsic_share_of_variance"].items():
            nm = ARRAYS.get(arr, arr)
            F.add(f"g2_band_{bt}_intrinsic_share_{nm}_pct", pct(100 * v), 100 * v,
                  f"{kb}.intrinsic_share_of_variance.{arr} (share of in-band channel variance)")
            F.add(f"g2_band_{bt}_intrinsic_share_{nm}_pct_1dp", f"{100 * v:.1f}%", 100 * v, f"{kb}.intrinsic_share_of_variance.{arr}")
            F.add(f"g2_band_{bt}_brain_room_share_{nm}_pct", pct(100 * (1 - v)), 100 * (1 - v),
                  f"{BAND} :: derived: 1 - bands['{b}'].intrinsic_share_of_variance.{arr}")
        F.add(f"g2_band_{bt}_enbw_hz", num(x["enbw_hz"], 1), x["enbw_hz"], f"{kb}.enbw_hz")
        F.add(f"g2_band_{bt}_brain_grad_ftcm", num(x["brain_rms_grad_fT_cm"], 1), x["brain_rms_grad_fT_cm"], f"{kb}.brain_rms_grad_fT_cm")
        F.add(f"g2_band_{bt}_brain_mag_measured_ft", num(x["brain_rms_mag_fT_measured"], 0), x["brain_rms_mag_fT_measured"],
              f"{kb}.brain_rms_mag_fT_measured")
        for kind in ("mag", "grad"):
            F.add(f"g2_band_{bt}_room_explained_{kind}_pct", pct(100 * x["env_explained"][kind]), 100 * x["env_explained"][kind],
                  f"{kb}.env_explained.{kind}")
        F.ratio(f"g2_band_{bt}_opm_response_gain", x["opm_response_gain"], f"{kb}.opm_response_gain (100-Hz first-order OPM response)")
        F.ratio(f"g2_band_{bt}_simulated_over_predicted_var", x["simulated_over_predicted_variance"], f"{kb}.simulated_over_predicted_variance")
        for variant in ("ideal", "opm_response"):
            for key, e in x[variant].items():
                arr, comp, cond = key.split("/")
                if comp != "combined":
                    continue
                base = f"g2_band_{bt}_{'' if variant == 'ideal' else 'resp_'}{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
                src = (f"{kb}.{variant}['{key}'] (ratio and ci95 as stored; ci95 is a TARGET bootstrap, not a parcel bootstrap; "
                       f"{'100-Hz first-order OPM response' if variant != 'ideal' else 'ideal OPM response'})")
                F.ratio(base, e["ratio"], src, e["ci95"], three=True)
                F.add(f"{base}_share", pct(100 * e["share_opm_better"]), 100 * e["share_opm_better"], f"{kb}.{variant}['{key}'].share_opm_better")
    F.add("g2_band_opm_response_corner_hz", f"{band['opm_response']['corner_hz']:g}", band["opm_response"]["corner_hz"],
          f"{BAND} :: opm_response.corner_hz")


def g2_primary_facts(F, d, root):
    for est, pre in ESTIMATORS.items():
        for key, e in d["primary"][est].items():
            arr, comp, cond = key.split("/")
            base = f"g2_{pre}{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
            oracle = est == "oracle"
            cmp_entry(F, base, e, G2, f"primary.{est}['{key}']", three=oracle, shares=True,
                      parcels=est in ("oracle", "peak"), counts=oracle)
    p = d["primary"]
    for arr in ("opm_dense", "opm_matched"):
        r = 2 ** p["peak"][f"{arr}/mag/intrinsic+brain"]["median_log2"]
        F.add(f"g2_peak_{ARRAYS[arr]}_vs_mag_ib_gain_pct", pct(100 * (r - 1)), 100 * (r - 1),
              f"{G2} :: derived: 2**primary.peak['{arr}/mag/intrinsic+brain'].median_log2 - 1 (peak-channel SNR gain over 102 magnetometers)")
    vals = [2 ** p["oracle"][f"{a}/grad/intrinsic+brain"]["median_log2"] for a in ("opm_matched", "opm_dense")]
    F.span("g2_vs_grad_ib", vals, r2, f"{G2} :: primary.oracle['opm_matched|opm_dense/grad/intrinsic+brain'] (ratio = 2**median_log2)")
    vals = [2 ** p["oracle"][f"{a}/mag/intrinsic+brain"]["median_log2"] for a in ("opm_matched", "opm_dense")]
    F.span("g2_vs_mag_ib", vals, r2, f"{G2} :: primary.oracle['opm_matched|opm_dense/mag/intrinsic+brain'] (ratio = 2**median_log2)")
    vals = [2 ** p["oracle"][f"{a}/grad/intrinsic"]["median_log2"] for a in ARRAYS]
    F.span("g2_vs_grad_int", vals, r2, f"{G2} :: primary.oracle['*/grad/intrinsic'] (ratio = 2**median_log2)")
    ints = [p["oracle"][f"{a}/combined/intrinsic"]["share_opm_better"] * 100 for a in ARRAYS]
    F.add("g2_vs_combined_int_share_max", pct(max(ints)), max(ints), f"{G2} :: primary.oracle['*/combined/intrinsic'].share_opm_better (largest)")
    # plug-in over oracle (median per array of plug-in / oracle detectability)
    for key, v in d["plugin_over_oracle_median"].items():
        arr, cs, cond, lab = key.split("/")
        if cs not in ("combined", "opm", "mag", "grad"):
            continue
        nm = ARRAYS.get(arr, arr) if arr != "squid" else f"squid_{cs}"
        F.add(f"g2_plugin{lab[1:]}s_loss_{nm}_{CONDS[cond]}_pct", pct(100 * (1 - v)), 100 * (1 - v),
              f"{G2} :: derived: 1 - plugin_over_oracle_median['{key}'] (median loss of detectability with the plug-in covariance)")
    for k2, v in d["break_even_opm_asd_fT"].items():
        arr, comp = k2.split("/")
        F.add(f"g2_breakeven_{ARRAYS[arr]}_vs_{comp}_asd", num(v, 2), v, f"{G2} :: break_even_opm_asd_fT['{k2}'] (fT/sqrt(Hz); intrinsic noise only)")
    be = [d["break_even_opm_asd_fT"][f"{a}/grad"] for a in ARRAYS]
    F.span("g2_breakeven_vs_grad_asd", be, lambda x: num(x, 2), f"{G2} :: break_even_opm_asd_fT['*/grad']", extra=(("int", lambda x: num(x, 0)),))
    for arr, x in d["projection"].items():
        nm = ARRAYS.get(arr, arr)
        F.ratio(f"g2_proj_retained_norm_{nm}", x["retained_signal_norm_median"], f"{G2} :: projection.{arr}.retained_signal_norm_median")
        cost = 100 * (1 - x["detect_projected_over_env_median"])
        F.add(f"g2_proj_cost_{nm}_pct", pct(cost), cost, f"{G2} :: derived: 1 - projection.{arr}.detect_projected_over_env_median")
        F.pct2(f"g2_proj_cost_{nm}_pct_2sf", cost, f"{G2} :: derived: 1 - projection.{arr}.detect_projected_over_env_median")
    mw = d["medial_wall"]
    F.add("g2_medialwall_n", count(mw["n_targets"]), mw["n_targets"], f"{G2} :: medial_wall.n_targets")
    F.add("g2_medialwall_share_pct", pct(100 * mw["share"]), 100 * mw["share"], f"{G2} :: medial_wall.share")
    for key, r in mw["ratios_without"].items():
        arr, comp, cond = key.split("/")
        F.ratio(f"g2_nowall_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}", r,
                f"{G2} :: medial_wall.ratios_without['{key}'] (median ratio without the medial-wall targets)", three=True)
    # history of the dense headline (dense vs combined, intrinsic + brain) over the model versions
    hist = [("v1", load(root, SKIN_V1)["models"]["bem3_5120"]["ratio"]["opm_dense/combined/intrinsic+brain"],
             f"{SKIN_V1} :: models.bem3_5120.ratio['opm_dense/combined/intrinsic+brain'] (v1 arrays, 5,120-triangle head surface, 1,000-target subset)"),
            ("v2", 1.129, f"{REGISTER} :: A-BEM-SKIN (v2 arrays: 1.147x coarse vs 1.129x refined; 1,000-target subset)"),
            ("v3", load(root, HSE)["variants"]["v3"]["opm_dense/combined/intrinsic+brain"]["ratio"],
             f"{HSE} :: variants.v3['opm_dense/combined/intrinsic+brain'].ratio"),
            ("v4", 2 ** d["primary"]["oracle"]["opm_dense/combined/intrinsic+brain"]["median_log2"],
             f"{G2} :: primary.oracle['opm_dense/combined/intrinsic+brain'] (2**median_log2)")]
    for ver, r, src in hist:
        F.ratio(f"g2_history_{ver}_dense_vs_combined_ib", r, src, three=True)


def g2_depth_facts(F, d):
    L = d["log2_ratio_vs_depth"]
    first = L["opm_dense/combined/intrinsic+brain"]
    for b in first:
        if b["n"] > 0:
            F.add(f"g2_depth_{b['lo']:.0f}_{b['hi']:.0f}_n", count(b["n"]), b["n"],
                  f"{G2} :: log2_ratio_vs_depth[*][{first.index(b)}].n (targets with depth in [lo, hi) mm)")
            F.add(f"g2_depth_bin_{b['lo']:.0f}_{b['hi']:.0f}_mm", f"{b['lo']:.0f}–{b['hi']:.0f}", [b["lo"], b["hi"]],
                  f"{G2} :: log2_ratio_vs_depth[*][{first.index(b)}].lo, .hi (depth bin below the scalp, mm)")
    for key, bins in L.items():
        arr, comp, cond = key.split("/")
        base = f"g2_depth_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
        src = f"{G2} :: log2_ratio_vs_depth['{key}'] (ratio = 2**median, the median over the bin's targets of log2(d_OPM/d_SQUID))"
        three = comp == "combined" and cond in ("intrinsic+brain", "projected")
        pop = [(b["lo"], b["hi"], 2 ** b["median"]) for b in bins if b["median"] is not None]
        for lo, hi, r in pop:
            F.ratio(f"{base}_{lo:.0f}_{hi:.0f}", r, src, three=three)
        F.span(f"{base}_bins", [r for *_, r in pop], r2, src + "; over the populated bins (n >= 10)")
        below = [i for i, (*_, r) in enumerate(pop) if r < 1]
        if below and below == list(range(below[0], len(pop))):
            lo, hi, _ = pop[below[0]]
            F.add(f"{base}_below1_from_mm", f"{lo:.0f} to {hi:.0f}", [lo, hi],
                  f"{G2} :: derived: log2_ratio_vs_depth['{key}'], the shallowest bin from which every deeper populated bin has ratio < 1")
            F.add(f"{base}_below1_from_bin", f"{lo:.0f}–{hi:.0f}", [lo, hi],
                  f"{G2} :: derived: log2_ratio_vs_depth['{key}'], the shallowest bin from which every deeper populated bin has "
                  "ratio < 1 (bin label, mm)")
    for arr in ARRAYS:
        for cond in ("intrinsic+brain", "projected"):
            key = f"{arr}/combined/{cond}"
            deep = [2 ** b["median"] for b in L[key] if b["median"] is not None and 35 <= b["lo"] < 60]
            F.span(f"g2_depth_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}_35_60", deep, r2,
                   f"{G2} :: log2_ratio_vs_depth['{key}'] bins 35-40 to 55-60 mm (2**median)", extra=(("3dp", r3),))
    # brain-noise multiplier: (i+b ratio) / (intrinsic ratio) per depth bin
    for arr in ARRAYS:
        for comp in COMPS:
            ib, it = L[f"{arr}/{comp}/intrinsic+brain"], L[f"{arr}/{comp}/intrinsic"]
            vals = []
            for a, b in zip(ib, it):
                if a["median"] is None or b["median"] is None:
                    continue
                m = 2 ** (a["median"] - b["median"])
                vals.append(m)
                F.ratio(f"g2_brainmult_{ARRAYS[arr]}_vs_{comp}_{a['lo']:.0f}_{a['hi']:.0f}", m,
                        f"{G2} :: derived: 2**(log2_ratio_vs_depth['{arr}/{comp}/intrinsic+brain'][bin].median - "
                        f"log2_ratio_vs_depth['{arr}/{comp}/intrinsic'][bin].median) (i+b over intrinsic, ratio of bin medians)")
            F.span(f"g2_brainmult_{ARRAYS[arr]}_vs_{comp}", vals, r2,
                   f"{G2} :: derived: i+b over intrinsic bin-median ratio, log2_ratio_vs_depth['{arr}/{comp}/*'], populated bins")
    # absolute detectability and peak amplitude by depth (dense and Neuromag)
    for key, bins in d["detectability_vs_depth"].items():
        arr, cs, cond = key.split("/")
        if cond not in ("intrinsic+brain", "projected") or cs not in ("combined", "opm"):
            continue
        nm = ARRAYS.get(arr, "squid")
        for b in bins:
            if b["median"] is not None:
                F.add(f"g2_detect_{nm}_{CONDS[cond]}_{b['lo']:.0f}_{b['hi']:.0f}", num(b["median"], 2), b["median"],
                      f"{G2} :: detectability_vs_depth['{key}'] bin {b['lo']:g}-{b['hi']:g} mm .median (d of a 10-nAm dipole)")
    for arr, bins in d["amplitude_vs_depth"].items():
        nm = ARRAYS.get(arr, arr)
        unit, scale = ("pT/m", 1e12) if arr == "squid_grad" else ("fT", 1e15)
        for b in bins:
            if b["median"] is not None:
                F.add(f"g2_amp_{nm}_{b['lo']:.0f}_{b['hi']:.0f}", num(b["median"] * scale, 0), b["median"] * scale,
                      f"{G2} :: amplitude_vs_depth['{arr}'] bin {b['lo']:g}-{b['hi']:g} mm .median ({unit}, peak channel, 10 nAm)")


def g2_bridge_facts(F, d):
    b = d["bridge_to_sphere"]
    centres = b["depth_centers_mm"]
    for arr in ARRAYS:
        bins = b[f"{arr}_ratio_vs_depth"]
        src = f"{G2} :: bridge_to_sphere.{arr}_ratio_vs_depth (median per depth bin of the OPM / Neuromag magnetometer peak-field ratio)"
        vals = []
        for x in bins:
            if x["median"] is None:
                continue
            vals.append(x["median"])
            nm = f"g2_leadfield_{ARRAYS[arr]}_vs_mag_{x['lo']:.0f}_{x['hi']:.0f}"
            F.ratio(nm, x["median"], src)
            F.add(nm + "_ratio_1dp", num(x["median"], 1), x["median"], src)
        F.span(f"g2_leadfield_{ARRAYS[arr]}_vs_mag_bins", vals, r2, src + "; populated bins", extra=(("1dp", lambda v: num(v, 1)),))
    for model, nm in (("realistic_standoffs", "real"), ("jas_xi0_18", "jas")):
        src = f"{G2} :: bridge_to_sphere.sphere_ratio_vs_depth.{model} (sphere peak-field ratio at the bin centre)"
        for i, x in enumerate(b["opm_dense_ratio_vs_depth"]):
            if x["median"] is None:
                continue
            c, r = centres[i], b["sphere_ratio_vs_depth"][model][i]
            base = f"g2_leadfield_sphere_{nm}_{x['lo']:.0f}_{x['hi']:.0f}"
            F.ratio(base, r, src + f", {c:g} mm")
            F.add(base + "_ratio_1dp", num(r, 1), r, src + f", {c:g} mm")
    for model, nm in (("realistic_standoffs", "sphere_real"), ("jas_xi0_18", "sphere_jas")):
        for eta, v in b["sphere_d_eq_mm"][model].items():
            if v is not None:
                F.add(f"g2_bridge_deq_{nm}_eta{tag(float(eta))}_mm", num(v, 1), v, f"{G2} :: bridge_to_sphere.sphere_d_eq_mm.{model}['{eta}']")
    v = b["sphere_d_eq_mm"]["realistic_standoffs"]["3"]
    F.add("g2_bridge_deq_sphere_real_eta3_mm_2dp", num(v, 2), v, f"{G2} :: bridge_to_sphere.sphere_d_eq_mm.realistic_standoffs['3']")
    v = b["sphere_d_eq_mm"]["jas_xi0_18"]["3"]
    F.add("g2_bridge_deq_sphere_jas_eta3_mm_3dp", num(v, 3), v, f"{G2} :: bridge_to_sphere.sphere_d_eq_mm.jas_xi0_18['3']")
    for arr in ARRAYS:
        nm = ARRAYS[arr]
        for eta, v in b[f"{arr}_d_eq_mm"].items():
            if v is not None:
                F.add(f"g2_bridge_deq_{nm}_eta{tag(float(eta))}_mm", num(v, 1), v,
                      f"{G2} :: bridge_to_sphere.{arr}_d_eq_mm['{eta}'] (intrinsic noise, peak-channel SNR, from the depth-binned ratio)")
        ahead, behind = b[f"{arr}_eta_opm_ahead_at_all_depths"], b[f"{arr}_eta_squid_ahead_at_all_depths"]
        F.add(f"g2_bridge_{nm}_eta_opm_ahead_max", eta_text(max(ahead)), max(ahead),
              f"{G2} :: bridge_to_sphere.{arr}_eta_opm_ahead_at_all_depths (largest)")
        F.add(f"g2_bridge_{nm}_eta_squid_ahead_min", eta_text(min(behind)), min(behind),
              f"{G2} :: bridge_to_sphere.{arr}_eta_squid_ahead_at_all_depths (smallest)")
    deqs = [b[f"{a}_d_eq_mm"]["3"] for a in ARRAYS]
    F.span("g2_bridge_deq_eta3_arrays_mm", deqs, lambda x: num(x, 1), f"{G2} :: bridge_to_sphere.*_d_eq_mm['3'] (matched, opm204, dense)")


def g2_strata_lobe_facts(F, d):
    for key, bins in d["strata"].items():
        arr, comp, cond, what = key.split("/")
        for x in bins:
            if what == "dist_inner_skull_mm":
                lab = "gt10" if x["hi"] > 50 else f"{x['lo']:.0f}_{x['hi']:.0f}"
                base = f"g2_skulldist_{lab}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
                nlab = f"g2_skulldist_{lab}_n"
            else:
                lab = f"{x['lo']:.0f}_{min(x['hi'], 90):.0f}"
                base = f"g2_orient_{lab}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
                nlab = f"g2_orient_{lab}_n"
            F.ratio(base, x["ratio"], f"{G2} :: strata['{key}'] bin {x['lo']:g}-{x['hi']:g} .ratio (2**median of log2 ratios)")
            if nlab not in F:
                F.add(nlab, count(x["n"]), x["n"], f"{G2} :: strata['{key}'] bin {x['lo']:g}-{x['hi']:g} .n")
    for key, lobes in d["by_lobe"].items():
        arr, comp, cond = key.split("/")
        for lobe, e in lobes.items():
            base = f"g2_lobe_{lobe}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
            cmp_entry(F, base, e, G2, f"by_lobe['{key}'].{lobe}", three=False, shares=comp == "combined")
            if f"g2_lobe_{lobe}_n" not in F:
                F.add(f"g2_lobe_{lobe}_n", count(e["n"]), e["n"], f"{G2} :: by_lobe['{key}'].{lobe}.n")
                F.add(f"g2_lobe_{lobe}_n_parcels", e["ci_method"].split("(")[1].rstrip(")"), int(e["ci_method"].split("(")[1].rstrip(")")),
                      f"{G2} :: by_lobe['{key}'].{lobe}.ci_method")
        if comp == "combined":
            F.span(f"g2_lobe_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", [2 ** e["median_log2"] for e in lobes.values()], r2,
                   f"{G2} :: by_lobe['{key}'].*.median_log2 (2**, over the six lobes)")


def g2_controls_facts(F, d):
    t = d["triaxial_control"]
    for key, e in t.items():
        if not isinstance(e, dict) or "median_log2" not in e:
            continue
        variant, comp, cond = key.split("/")
        vt = {"equal_noise": "equal", "tangential_noise_x2": "tang2x"}[variant]
        cn = "matched" if comp == "over_matched_normal_only" else comp
        cmp_entry(F, f"g2_triax_{vt}_vs_{cn}_{CONDS[cond]}", e, G2, f"triaxial_control['{key}']", three=True, shares=True)
    p = d["patches"]
    for r in ("5", "10", "20"):
        a, n = p["area_cm2"][r], p["net_over_scalar_moment"][r]
        F.add(f"g2_patch{r}mm_area_cm2", num(a["median"], 2) if a["median"] < 10 else num(a["median"], 1), a["median"],
              f"{G2} :: patches.area_cm2['{r}'].median")
        F.span(f"g2_patch{r}mm_area_p5_p95_cm2", [a["p5"], a["p95"]], lambda x: num(x, 2), f"{G2} :: patches.area_cm2['{r}'].p5/p95")
        F.ratio(f"g2_patch{r}mm_net_moment", n["median"], f"{G2} :: patches.net_over_scalar_moment['{r}'].median (net / scalar moment)")
    for key, e in p["comparisons"].items():
        arr, comp, cond, rad = key.split("/")
        cmp_entry(F, f"g2_patch{rad}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}", e, G2, f"patches.comparisons['{key}']", three=True)
    for arr in ("opm_dense", "opm_matched"):
        for cond in ("intrinsic+brain", "projected"):
            vals = [2 ** p["comparisons"][f"{arr}/combined/{cond}/{r}mm"]["median_log2"] for r in (5, 10, 20)]
            F.span(f"g2_patches_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", vals, r2,
                   f"{G2} :: patches.comparisons['{arr}/combined/{cond}/5|10|20mm'] (2**median_log2)", extra=(("3dp", r3),))
    for key, v in p["median_detectability"].items():
        arr, cs, cond, rad, conv = key.split("/")
        nm = ARRAYS.get(arr, f"squid_{cs}")
        F.add(f"g2_patch{rad}_detect_{nm}_{CONDS[cond]}_{conv}", num(v, 2), v, f"{G2} :: patches.median_detectability['{key}']")


SENS_NAMES = {"opm_asd_7fT": "asd7ft", "opm_asd_10fT": "asd10ft", "opm_asd_15fT": "asd15ft", "opm_asd_20fT": "asd20ft",
              "opm_asd_30fT": "asd30ft", "background_corr_5mm": "bgcorr5mm", "background_corr_10mm": "bgcorr10mm",
              "background_mag_calibrated": "bgmagcal", "squid_measured_spectrum": "squidmeas", "head_x+5mm": "head_xp5mm",
              "head_x-5mm": "head_xm5mm", "head_y+5mm": "head_yp5mm", "head_y-5mm": "head_ym5mm", "head_z+5mm": "head_zp5mm",
              "head_z-5mm": "head_zm5mm", "head_pitch+5deg": "head_pitchp5deg", "head_pitch-5deg": "head_pitchm5deg",
              "head_well_fitted": "head_wellfitted", "gap_3mm": "gap3mm", "gap_6mm": "gap6mm"}


def g2_sensitivity_facts(F, d):
    s = d["sensitivity"]
    for variant, entries in s.items():
        vt = SENS_NAMES[variant]
        for key, e in entries.items():
            if isinstance(e, dict) and "median_log2" in e:
                arr, comp, cond = key.split("/")
                if comp == "combined":
                    cmp_entry(F, f"g2_sens_{vt}_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", e, G2, f"sensitivity.{variant}['{key}']", three=True)
    for arr in ARRAYS:
        for cond in ("intrinsic+brain", "projected"):
            vals = [2 ** s[v][f"{arr}/combined/{cond}"]["median_log2"] for v in s if v.startswith("head_") and v != "head_well_fitted"]
            F.span(f"g2_sens_headpos_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", vals, r2,
                   f"{G2} :: sensitivity.head_[xyz]+-5mm, head_pitch+-5deg['{arr}/combined/{cond}'] (2**median_log2)", extra=(("3dp", r3),))
    sc = s["background_mag_calibrated"]["scale_over_primary"]
    F.ratio("g2_sens_bgmagcal_scale", sc, f"{G2} :: sensitivity.background_mag_calibrated.scale_over_primary (background variance scale)")
    ir = s["squid_measured_spectrum"]["intrinsic_rms"]
    for k, nm in (("mag_fT", "mag_ft"), ("grad_fT_cm", "grad_ftcm"), ("brochure_mag_fT", "brochure_mag_ft"), ("brochure_grad_fT_cm", "brochure_grad_ftcm")):
        F.add(f"g2_sens_squidmeas_intrinsic_{nm}", num(ir[k], 1), ir[k],
              f"{G2} :: sensitivity.squid_measured_spectrum.intrinsic_rms.{k} (median in-band RMS)")
    key = "opm_dense/combined/intrinsic+brain"
    gap6 = 2 ** (s["gap_6mm"][key]["median_log2"] - d["primary"]["oracle"][key]["median_log2"])
    F.add("g2_gap6mm_dense_cost_pct", pct(100 * (1 - gap6)), 100 * (1 - gap6),
          f"{G2} :: derived: 1 - 2**sensitivity.gap_6mm['opm_dense/combined/intrinsic+brain'].median_log2 / 2**primary.oracle[same].median_log2")
    j = d["sensitivity_joint_asd_gap"]
    for key, x in j.items():
        gap, asd = key.split("/")
        for arr, e in x.items():
            F.ratio(f"g2_joint_{gap}_{asd.lower()}_{ARRAYS[arr]}_vs_combined_ib", e["ratio"],
                    f"{G2} :: sensitivity_joint_asd_gap['{key}'].{arr} (ratio and ci95 as stored, intrinsic + brain)", e["ci95"], three=True)
    for arr in ("opm_dense", "opm_matched"):
        F.span(f"g2_joint_{ARRAYS[arr]}_vs_combined_ib", [x[arr]["ratio"] for x in j.values()], r2,
               f"{G2} :: sensitivity_joint_asd_gap[*].{arr}.ratio (3 gaps x 3 noise levels)", extra=(("3dp", r3),))


def g2_convergence_facts(F, d, root):
    c = d["convergence"]
    b = c["background_grid_vs_fullres"]
    F.add("g2_conv_bggrid_max_change_log2", num(b["max_abs_change_log2"], 3), b["max_abs_change_log2"],
          f"{G2} :: convergence.background_grid_vs_fullres.max_abs_change_log2 (largest change of a median log2 ratio)")
    F.add("g2_conv_bggrid_scale_ratio", num(b["scale_ratio_full_over_grid"], 3), b["scale_ratio_full_over_grid"],
          f"{G2} :: convergence.background_grid_vs_fullres.scale_ratio_full_over_grid")
    F.add("g2_conv_bem1_refinement_max_change_log2", num(c["bem"]["refinement_max_abs_change_log2"], 4), c["bem"]["refinement_max_abs_change_log2"],
          f"{G2} :: convergence.bem.refinement_max_abs_change_log2 (1-layer BEM, 5,120 vs 20,480 triangles)")
    for model in ("bem3_head5120", "bem1_5120", "bem1_20480"):
        for key, v in c["bem"][model]["median_log2"].items():
            arr, comp, cond = key.split("/")
            if comp == "combined":
                F.ratio(f"g2_conv_{model}_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", 2 ** v,
                        f"{G2} :: convergence.bem.{model}.median_log2['{key}'] (2**, 1,000-target subset)", three=True)
        g = c["bem"][model]["gain_rel_diff_median"]
        for arr, v in g.items():
            F.add(f"g2_conv_{model}_gain_diff_{ARRAYS.get(arr, arr)}_pct", pct(100 * v), 100 * v,
                  f"{G2} :: convergence.bem.{model}.gain_rel_diff_median.{arr} (median relative gain change vs the primary model)")
    for key, v in c["bem"]["subset_reference_median_log2"].items():
        arr, comp, cond = key.split("/")
        if comp == "combined":
            F.ratio(f"g2_conv_subset_{ARRAYS[arr]}_vs_combined_{CONDS[cond]}", 2 ** v,
                    f"{G2} :: convergence.bem.subset_reference_median_log2['{key}'] (2**, primary model on the 1,000-target subset)", three=True)
    ci = c["coil_integration"]
    F.add("g2_conv_squid_4pt_pct", pct(100 * ci["squid_4pt_vs_accurate_rel_diff_median"]), 100 * ci["squid_4pt_vs_accurate_rel_diff_median"],
          f"{G2} :: convergence.coil_integration.squid_4pt_vs_accurate_rel_diff_median")
    for arr in ("opm_matched", "opm_dense"):
        x = ci[f"{arr}_cell_vs_point_peak_rel_diff"]
        for stat in ("median", "p95", "max"):
            F.add(f"g2_conv_cell_vs_point_{ARRAYS[arr]}_{stat}_pct", pct(100 * x[stat]), 100 * x[stat],
                  f"{G2} :: convergence.coil_integration.{arr}_cell_vs_point_peak_rel_diff.{stat} (peak-channel field)")
    F.add("g2_conv_target_sampling_max_change_log2", num(c["target_sampling"]["max_abs_change_log2"], 3), c["target_sampling"]["max_abs_change_log2"],
          f"{G2} :: convergence.target_sampling.max_abs_change_log2 (oct-6 vs random full-resolution targets)")
    w = max(c["whitening_tolerance_max_rel_change"].values())
    F.add("g2_conv_whitening_max_change", num(w, 0), w, f"{G2} :: convergence.whitening_tolerance_max_rel_change (largest)")
    # head-surface refinement (bem_skin_refinement.json; v4 arrays and v1 arrays)
    for rel, pre in ((SKIN, "g2_skin"), (SKIN_V1, "g2_skin_v1")):
        sk = load(root, rel)
        for model, x in sk["models"].items():
            for key, r in x["ratio"].items():
                arr, comp, cond = key.split("/")
                F.ratio(f"{pre}_{model}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}", r,
                        f"{rel} :: models.{model}.ratio['{key}'] (1,000-target subset)", three=True)
            for arr, g in x.get("gain_rel_diff_vs_bem3_5120", {}).items():
                F.add(f"{pre}_{model}_gain_diff_{ARRAYS.get(arr, arr)}_pct", pct(100 * g["median"]), 100 * g["median"],
                      f"{rel} :: models.{model}.gain_rel_diff_vs_bem3_5120.{arr}.median")
    sph = load(root, SPHERE)
    for model, rows in sph["models"].items():
        nm = "5120" if "5,120" in model else "20480"
        for x in rows:
            h = tag(x["height_mm"])
            for stat in ("median", "p95"):
                v = 100 * x[stat]
                F.add(f"g2_bemsphere_{nm}_h{h}mm_{stat}_pct", pct(v), v,
                      f"{SPHERE} :: models['{model}'][height {x['height_mm']:g} mm].{stat} (relative field error)")
                F.pct2(f"g2_bemsphere_{nm}_h{h}mm_{stat}_pct_2sf", v, f"{SPHERE} :: models['{model}'][height {x['height_mm']:g} mm].{stat}")
    nm_ = load(root, NEAR)
    for arr, x in nm_["arrays"].items():
        nm = ARRAYS[arr]
        lp = x["lowest_point_height_mm"]
        F.add(f"g2_nearmesh_{nm}_lowest_point_mm", num(lp["min"], 1), lp["min"],
              f"{NEAR} :: arrays.{arr}.lowest_point_height_mm.min (integration points)")
        F.add(f"g2_nearmesh_{nm}_lowest_point_mm_2dp", num(lp["min"], 2), lp["min"], f"{NEAR} :: arrays.{arr}.lowest_point_height_mm.min")
        for row in x["point_error_by_height"]:
            if row["n_points"] == 0:
                continue
            lab = f"{tag(row['lo_mm'])}_{tag(row['hi_mm'])}mm"
            for stat in ("median", "p95"):
                v = 100 * row[stat]
                F.add(f"g2_nearmesh_{nm}_{lab}_{stat}_pct", pct(v), v,
                      f"{NEAR} :: arrays.{arr}.point_error_by_height[{row['lo_mm']:g}-{row['hi_mm']:g} mm].{stat} (5,120 vs 20,480-triangle head surface)")
                F.pct2(f"g2_nearmesh_{nm}_{lab}_{stat}_pct_2sf", v,
                       f"{NEAR} :: arrays.{arr}.point_error_by_height[{row['lo_mm']:g}-{row['hi_mm']:g} mm].{stat}")
            F.add(f"g2_nearmesh_{nm}_{lab}_n", count(row["n_points"]), row["n_points"],
                  f"{NEAR} :: arrays.{arr}.point_error_by_height[{row['lo_mm']:g}-{row['hi_mm']:g} mm].n_points")
        ce = x["channel_error_relative"]
        F.add(f"g2_nearmesh_{nm}_channels_over1pct", str(ce["n_channels_above_1pct"]), ce["n_channels_above_1pct"],
              f"{NEAR} :: arrays.{arr}.channel_error_relative.n_channels_above_1pct")
        F.add(f"g2_nearmesh_{nm}_channels_over2pct", str(ce["n_channels_above_2pct"]), ce["n_channels_above_2pct"],
              f"{NEAR} :: arrays.{arr}.channel_error_relative.n_channels_above_2pct")
        F.add(f"g2_nearmesh_{nm}_channel_error_max_pct", pct(100 * ce["max"]), 100 * ce["max"], f"{NEAR} :: arrays.{arr}.channel_error_relative.max")
    lows = [x["lowest_point_height_mm"]["min"] for x in nm_["arrays"].values()]
    F.add("g2_nearmesh_lowest_point_mm_2dp", num(min(lows), 2), min(lows), f"{NEAR} :: arrays.*.lowest_point_height_mm.min (both arrays)")
    # v3 -> v4 head-model change (head_surface_effect.json)
    hse = load(root, HSE)
    vnames = {"v3": "v3", "v3_arrays_v4_bem": "v3arrays_v4bem", "v4_sites_v3_axes": "v4sites_v3axes", "v4": "v4",
              "v4_dense_above_v3_low": "v4_nolowsites"}
    for v, x in hse["variants"].items():
        vt = vnames[v]
        for arr, n in x["sites"].items():
            F.add(f"g2_hse_{vt}_{ARRAYS[arr]}_sites", str(n), n, f"{HSE} :: variants.{v}.sites.{arr}")
        for key, e in x.items():
            if not isinstance(e, dict) or "ratio" not in e:
                continue
            arr, comp, cond = key.split("/")
            base = f"g2_hse_{vt}_{ARRAYS[arr]}_vs_{comp}_{CONDS[cond]}"
            F.ratio(base, e["ratio"], f"{HSE} :: variants.{v}['{key}'].ratio and ci95 (oracle, point gains)", e["ci95"], three=True)
            F.add(f"{base}_share", pct(100 * e["share_opm_better"]), 100 * e["share_opm_better"], f"{HSE} :: variants.{v}['{key}'].share_opm_better")
            for bb in e["by_depth"]:
                F.ratio(f"{base}_{bb['lo']:.0f}_{bb['hi']:.0f}", bb["ratio"],
                        f"{HSE} :: variants.{v}['{key}'].by_depth[{bb['lo']:g}-{bb['hi']:g} mm].ratio")
    F.add("g2_hse_v4_dense_sites_below_v3_low", str(hse["v4_dense_sites_below_it"]), hse["v4_dense_sites_below_it"],
          f"{HSE} :: v4_dense_sites_below_it")


# ------------------------------------------------------------------------------------------------
def region_facts(F, root):
    """Per-parcel and per-group medians over the stored per-target rows (both hemispheres pooled)."""
    rows = read_csv(root, TARGETS)
    parcel = np.array([r["region"].split(".", 1)[1] for r in rows])
    col = lambda c: np.array([float(r[c]) for r in rows])  # noqa: E731
    depth, amag = col("depth_mm"), col("amp_squid_mag")
    amps = {"dense": col("amp_opm_dense"), "matched": col("amp_opm_matched")}
    det = {(a, cond): col(f"detect_opm_{a}_opm_{cond}") for a in ("dense", "matched") for cond in ("intrinsic", "intrinsic+brain", "projected")}
    sq = {(c, cond): col(f"detect_squid_{c}_{cond}") for c in COMPS for cond in ("intrinsic", "intrinsic+brain", "projected")}
    prow = {(r["hemi"], r["vertno"]): r for r in read_csv(root, PATCHES)}
    keyed = [prow[(r["hemi"], r["vertno"])] for r in rows]
    pcol = lambda c: np.array([float(r[c]) for r in keyed])  # noqa: E731
    peak3 = []
    for g, parcels in REGIONS.items():
        m = np.isin(parcel, parcels)
        what = f"targets with region lh./rh.{'|'.join(parcels)} ({int(m.sum())} targets, both hemispheres)"
        F.add(f"reg_{g}_n", count(m.sum()), int(m.sum()), f"{TARGETS} :: derived: count of {what}")
        F.add(f"reg_{g}_depth_mm", num(np.median(depth[m]), 1), float(np.median(depth[m])), f"{TARGETS} :: derived: median depth_mm over {what}")
        for a, amp in amps.items():
            r = float(np.median(amp[m] / amag[m]))
            src = f"{TARGETS} :: derived: median over {what} of amp_opm_{a} / amp_squid_mag (per-target peak-field ratio; CSV values as stored)"
            F.ratio(f"reg_{g}_peakfield_{a}_vs_mag", r, src)
            F.add(f"reg_{g}_peakfield_{a}_vs_mag_ratio_1dp", num(r, 1), r, src)
            if a == "dense" and g in ("superiortemporal", "parahippocampal", "precentral"):
                peak3.append(r)
        for (a, cond), dv in det.items():
            for c in COMPS:
                r = float(np.median(dv[m] / sq[(c, cond)][m]))
                F.ratio(f"reg_{g}_{a}_vs_{c}_{CONDS[cond]}", r, f"{TARGETS} :: derived: median over {what} of the per-target ratio "
                        f"detect_opm_{a}_opm_{cond} / detect_squid_{c}_{cond} (unweighted; CSV values as stored, 4 decimals)")
        for a in ("dense", "matched"):
            for c in COMPS:
                for cond in ("intrinsic+brain", "projected"):
                    num_ = pcol(f"detect_opm_{a}_opm_{cond}_10mm_fixed_total")[m]
                    den = pcol(f"detect_squid_{c}_{cond}_10mm_fixed_total")[m]
                    r = float(np.median(num_ / den))
                    F.ratio(f"reg_{g}_patch10mm_{a}_vs_{c}_{CONDS[cond]}", r, f"{PATCHES}, {TARGETS} :: derived: median over {what} "
                            f"(joined on hemi, vertno) of detect_opm_{a}_opm_{cond}_10mm_fixed_total / detect_squid_{c}_{cond}_10mm_fixed_total")
    src = f"{TARGETS} :: derived: min/max of reg_<parcel>_peakfield_dense_vs_mag over superiortemporal, parahippocampal, precentral"
    F.span("reg_three_parcels_peakfield_dense_vs_mag", peak3, r2, src, extra=(("1dp", lambda x: num(x, 1)),))
    mes = [float(np.median(depth[parcel == p])) for p in ("parahippocampal", "entorhinal")]
    F.span("reg_mesial_parcels_depth", mes, lambda x: num(x, 1),
           f"{TARGETS} :: derived: median depth_mm of parahippocampal and of entorhinal targets")


# ------------------------------------------------------------------------------------------------
def anatomy_noise_extra_facts(F, d, root):
    """Neuromag's gradiometer count; the noise budget per channel (variance shares, the brain-noise inflation over the
    sensor floor, the single-channel noise ratio across the OPM sweep); lobe depths and fields from the per-target CSV;
    regional fields and within-system contrasts; every Desikan-Killiany parcel; the bridge's agreement window; the
    cortex the 4-mm rule leaves out."""
    sq = d["arrays"]["squid"]
    F.add("g2_squid_grad_channels", str(sq["channels"] - sq["sites"]), sq["channels"] - sq["sites"],
          f"{G2} :: derived: arrays.squid.channels - arrays.squid.sites (two planar gradiometers per sensor site)")
    nc = d["noise_composition"]
    terms = (("intrinsic", "sensor"), ("brain", "brain"), ("env", "room"))
    for arr, x in nc.items():
        nm = ARRAYS.get(arr, arr)
        tot = sum(x[f"{t}_rms"] ** 2 for t, _ in terms)
        for t, lab in terms:
            v = 100 * x[f"{t}_rms"] ** 2 / tot
            F.add(f"g2_noise_share_{nm}_{lab}_pct", pct(v), v,
                  f"{G2} :: derived: noise_composition.{arr}.{t}_rms**2 over the sum of the squared intrinsic, brain and env RMS "
                  "(share of the in-band channel variance, each component's median channel variance)")
        infl = math.sqrt(1 + (x["brain_rms"] / x["intrinsic_rms"]) ** 2)
        F.add(f"g2_noise_inflation_{nm}_1dp", num(infl, 1), infl,
              f"{G2} :: derived: sqrt(1 + (noise_composition.{arr}.brain_rms / intrinsic_rms)**2) (in-band noise RMS with brain "
              "noise over that of sensor noise alone)")
    infl = {a: math.sqrt(1 + (nc[a]["brain_rms"] / nc[a]["intrinsic_rms"]) ** 2) for a in ("squid_mag", "opm_dense", "opm_matched")}
    F.ratio("g2_noise_inflation_mag_over_dense", infl["squid_mag"] / infl["opm_dense"],
            f"{G2} :: derived: the brain-noise inflation of the Neuromag magnetometers over that of the dense OPM channels "
            "(sqrt(1 + (brain_rms / intrinsic_rms)**2) of noise_composition.squid_mag over the same of opm_dense)")
    opm_rms = d["noise_validation"]["model"]["intrinsic_rms_opm_fT"]
    mag = nc["squid_mag"]
    for asd, v in opm_rms.items():
        o, vi = nc["opm_dense"], v * 1e-15
        with_room = math.sqrt(vi ** 2 + o["brain_rms"] ** 2 + o["env_rms"] ** 2) / math.sqrt(
            sum(mag[f"{t}_rms"] ** 2 for t, _ in terms))
        no_room = math.sqrt(vi ** 2 + o["brain_rms"] ** 2) / math.sqrt(mag["intrinsic_rms"] ** 2 + mag["brain_rms"] ** 2)
        src = (f"{G2} :: derived: single-channel in-band noise RMS of a dense-array OPM over a Neuromag magnetometer with the OPM "
               f"sensor noise noise_validation.model.intrinsic_rms_opm_fT['{asd}'] and noise_composition.opm_dense / squid_mag")
        F.ratio(f"g2_noise_rms_dense_vs_mag_asd{asd}", with_room, src + " (sensor, brain and room field)")
        F.ratio(f"g2_noise_rms_dense_vs_mag_noroom_asd{asd}", no_room, src + " (sensor and brain noise)")
    rows = read_csv(root, TARGETS)
    parcel = np.array([r["region"].split(".", 1)[1] for r in rows])
    lobe = np.array([r["lobe"] for r in rows])
    col = lambda c: np.array([float(r[c]) for r in rows])  # noqa: E731
    depth, amag, adense = col("depth_mm"), col("amp_squid_mag"), col("amp_opm_dense")
    dd, dm = col("detect_opm_dense_opm_intrinsic+brain"), col("detect_opm_matched_opm_intrinsic+brain")
    dc = col("detect_squid_combined_intrinsic+brain")
    for lb in sorted(set(lobe) - {"other"}):
        m = lobe == lb
        what = f"the {int(m.sum())} targets with lobe '{lb}'"
        F.add(f"g2_lobe_{lb}_depth_mm", num(np.median(depth[m]), 1), float(np.median(depth[m])),
              f"{TARGETS} :: derived: median depth_mm over {what} (unweighted)")
        F.add(f"g2_lobe_{lb}_amp_mag_ft", num(1e15 * np.median(amag[m]), 0), float(1e15 * np.median(amag[m])),
              f"{TARGETS} :: derived: median amp_squid_mag over {what} (peak field of the best magnetometer, 10-nAm dipole, fT)")
        F.add(f"g2_lobe_{lb}_amp_dense_ft", num(1e15 * np.median(adense[m]), 0), float(1e15 * np.median(adense[m])),
              f"{TARGETS} :: derived: median amp_opm_dense over {what} (peak field of the best dense-array OPM, 10-nAm dipole, fT)")
    for g, parcels in REGIONS.items():
        m = np.isin(parcel, parcels)
        what = f"targets with region lh./rh.{'|'.join(parcels)}"
        F.add(f"reg_{g}_amp_mag_ft", num(1e15 * np.median(amag[m]), 0), float(1e15 * np.median(amag[m])),
              f"{TARGETS} :: derived: median amp_squid_mag over {what} (best magnetometer, 10-nAm dipole, fT)")
        F.add(f"reg_{g}_amp_dense_ft", num(1e15 * np.median(adense[m]), 0), float(1e15 * np.median(adense[m])),
              f"{TARGETS} :: derived: median amp_opm_dense over {what} (best dense-array OPM, 10-nAm dipole, fT)")
    # contrasts within one system (absolute detectability, region over region)
    ph, st = parcel == "parahippocampal", parcel == "superiortemporal"
    tl, fl = lobe == "temporal", lobe == "frontal"
    for nm, x, xn in (("dense", dd, "detect_opm_dense_opm_intrinsic+brain"), ("combined", dc, "detect_squid_combined_intrinsic+brain")):
        F.ratio(f"reg_parahippocampal_vs_superiortemporal_{nm}_abs", float(np.median(x[ph]) / np.median(x[st])),
                f"{TARGETS} :: derived: median {xn} over the parahippocampal targets over the same over the superior temporal "
                "targets (one system's own detectability, region against region)")
        F.ratio(f"g2_lobe_temporal_vs_frontal_{nm}_abs", float(np.median(x[tl]) / np.median(x[fl])),
                f"{TARGETS} :: derived: median {xn} over the temporal-lobe targets over the same over the frontal-lobe targets "
                "(one system's own detectability, lobe against lobe)")
    # every Desikan-Killiany parcel (the medial wall, 'unknown', is not a parcel)
    vals = {}
    for p in sorted(set(parcel) - {"unknown"}):
        m = parcel == p
        rd, rm = float(np.median(dd[m] / dc[m])), float(np.median(dm[m] / dc[m]))
        src = (f"{TARGETS} :: derived: median over the {int(m.sum())} targets of parcel {p} (both hemispheres) of the per-target "
               "ratio detect_opm_<array>_opm_intrinsic+brain / detect_squid_combined_intrinsic+brain")
        F.ratio(f"reg_dk_{p}_dense_vs_combined_ib", rd, src)
        F.ratio(f"reg_dk_{p}_matched_vs_combined_ib", rm, src)
        F.add(f"reg_dk_{p}_n", count(m.sum()), int(m.sum()), f"{TARGETS} :: derived: targets of parcel {p} (both hemispheres)")
        vals[p] = (rd, rm)
    F.add("reg_dk_n_parcels", count(len(vals)), len(vals), f"{TARGETS} :: derived: Desikan-Killiany parcels with targets (medial "
          "wall excluded)")
    below = sorted(p for p, (rd, _) in vals.items() if rd < vals["parahippocampal"][0])
    F.add("reg_dk_dense_vs_combined_ib_below_parahippocampal_count", count(len(below)), len(below),
          f"{TARGETS} :: derived: parcels whose median dense ratio (reg_dk_<parcel>_dense_vs_combined_ib) is below the "
          f"parahippocampal parcel's: {', '.join(below)}")
    mes = np.isin(parcel, ("parahippocampal", "entorhinal"))
    rmes = float(np.median(dm[mes] / dc[mes]))
    belowm = sorted(p for p, (_, rm) in vals.items() if rm < rmes and p not in ("parahippocampal", "entorhinal"))
    F.add("reg_dk_matched_vs_combined_ib_below_mesial_temporal_count", count(len(belowm)), len(belowm),
          f"{TARGETS} :: derived: parcels outside the mesial temporal group whose median matched ratio is below the group's "
          f"({rmes:.4f}): {', '.join(belowm)}")
    # the bridge: the eta window in which the sphere (real standoffs) and both real arrays have an equal-SNR depth
    b = d["bridge_to_sphere"]
    sph = b["sphere_d_eq_mm"]["realistic_standoffs"]
    etas = [e for e in sph if sph[e] is not None and all(b[f"{a}_d_eq_mm"].get(e) is not None for a in ("opm_dense", "opm_matched"))]
    diffs = [b[f"{a}_d_eq_mm"][e] - sph[e] for a in ("opm_dense", "opm_matched") for e in etas]
    ev = sorted(float(e) for e in etas)
    src = (f"{G2} :: derived: bridge_to_sphere.opm_dense_d_eq_mm and opm_matched_d_eq_mm minus sphere_d_eq_mm.realistic_standoffs "
           "where all three have a crossing")
    F.add("g2_bridge_common_eta_min", eta_text(ev[0]), ev[0], src + " (smallest eta)")
    F.add("g2_bridge_common_eta_max", eta_text(ev[-1]), ev[-1], src + " (largest eta)")
    F.add("g2_bridge_deq_arrays_minus_sphere_real_mm_range", f"{signed(min(diffs), 1)} to {signed(max(diffs), 1)}",
          [min(diffs), max(diffs)], src + " (mm)")
    F.add("g2_excluded_vertices_pct_1dp", "8.7%", round(100 - 91.3, 1), f"{METHODS}, {REGISTER} :: section 3 ('The rule removes the "
          "cortex nearest the skull (8.7 %)'); A-BEM-DIST (91.3 % of the valid vertices are usable)")


# ------------------------------------------------------------------------------------------------
def facts(root: Path) -> dict:
    """Every G1A, G1B, G1C, G2 and adult-region fact, name -> {"value", "raw", "source"}."""
    root = Path(root)
    F = Facts()
    lit_facts(F)
    lit_study_facts(F, root)
    g1a_facts(F, root)
    g1b_facts(F, root)
    g1c_facts(F, root)
    d, band = load(root, G2), load(root, BAND)
    g2_config_facts(F, d)
    g2_noise_facts(F, d, band)
    g2_primary_facts(F, d, root)
    g2_depth_facts(F, d)
    g2_bridge_facts(F, d)
    g2_strata_lobe_facts(F, d)
    g2_controls_facts(F, d)
    g2_sensitivity_facts(F, d)
    g2_convergence_facts(F, d, root)
    region_facts(F, root)
    anatomy_noise_extra_facts(F, d, root)
    return dict(F)


if __name__ == "__main__":
    for name, f in facts(Path(__file__).resolve().parents[1]).items():
        print(f"{name}\t{f['value']}\t{f['source']}")
