#!/usr/bin/env python3
"""Export the manuscript's figures at print quality from the stored results (no computation is repeated).

Each figure is drawn by the same code that drew the report's figures (scripts/report_figures_*.py,
scripts/study_g3b_constant_gap.py, scripts/study_covariance_validation.py, scripts/study_noise_sensitivity.py), reading
the same result files, with five differences: the output is a vector PDF (embedded raster layers, such as the cortical
maps, at 600 dpi); panel labels written as '(a) Title' become a bold upper-case letter at the top left ('A  Title'), the
journal's convention; British spellings in the figures' text become American, the manuscript's spelling; the
explanatory notes drawn below the report's panels are left out, since the journal keeps that information in the
caption; and nothing is written into results/ (the drawing scripts' provenance records go to
paper/build/figure_export/). Figures that the manuscript shows as one figure but the report drew as two images are
stacked into one PDF (one file per figure, as the journal asks).

Legibility: the figures listed in FIXES are re-laid out at save time (size, type, legends; contents unchanged) so that
no text prints smaller than 7 pt under the documents' \\includegraphics options (PRINT below), or 6.5 pt in the figures
whose layout leaves no more room (LEGIBLE); every text drawn into each PDF is recorded
(build/figure_export/font_audit.json) and the printed sizes are checked (build/figure_export/print_audit.txt). Every
figure is also checked, after its fixes, for texts that overlap or touch (clashes(); recorded in
build/figure_export/clash_audit.json and reported in print_audit.txt), the re-laid-out ones also for keys that cover
data.

The main text shows the adult and the four principal smaller heads (PRINCIPAL) in its pediatric figures (Figs. 8-10):
those figures are drawn a second time with report_style.ANAT_ORDER restricted to them and saved with the suffix
'_principal'; the versions with all nine anatomies are supplementary figures (Figs. S16-S18), in which the 18-month
template, reported with the children in Section S10, is set apart with them (set_apart()), as it is in the other
supplementary figures that group the heads (Figs. S9, S12, S14, S15); in Figs. S7, S8 and S11 it keeps its place among
the templates, labelled misregistered (misregistered()). In every figure the children stand in the order of their
labels, A, B, C, as in the tables (pediatric()). In Fig. 2, child B's column title names Section S10, where the text
reports it.

Outputs: paper/figures/Figure_1.pdf ... Figure_10.pdf (main text) and Figure_S1.pdf ... Figure_S18.pdf (supplementary
material), each with a 300-dpi PNG preview beside it in paper/build/figure_export/; paper/figure_stacks.json gives where
each part of a stacked figure sits, so that a document can show one part alone (build_paper.py's trim()).

Usage: .venv/bin/python paper/export_figures.py [--no-draw]
"""
from __future__ import annotations

import argparse
import contextlib
import functools
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.axes  # noqa: E402
import matplotlib.axis  # noqa: E402
import matplotlib.category  # noqa: E402
import matplotlib.collections  # noqa: E402
import matplotlib.legend  # noqa: E402
import matplotlib.colors  # noqa: E402
import matplotlib.patches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.spines  # noqa: E402
import matplotlib.ticker  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

import report_style as style  # noqa: E402

EXPORT = HERE / "build" / "figure_export"   # drawing scripts' outputs and previews (not committed)
FIGDIR = HERE / "figures"                   # the manuscript's figure files
PDF_DPI, PNG_DPI = 600, 300
PANEL = re.compile(r"^(?:\(([a-l])\)(?:\s+|$)|([a-l])(?:\s{2,}|$))(.*)$", re.S)  # "(a) Title", "a  Title" or "a"
# figures drawn without panel letters whose panels the manuscript letters (continuing the figure they are stacked under):
# letters given to the axes that carry a left title, in reading order
EXTRA_LETTERS = {"Figure_noise_sensitivity_depth": "ghij"}
# single-panel figures of the manuscript that the report drew as one panel of a larger figure: that panel letter goes
DROP_LETTER = {"Figure_R12_geometry": "(d)", "Figure_R11_maps_adult": "(c)"}

# the main text's pediatric figures (Figs. 8-10) show the adult and these four smaller heads: the scaled adults and the
# 24- and 12-month templates; the 18-month template and the three children are shown in the supplement (Figs. S16-S18)
PRINCIPAL = ("adult", "school", "size2yr", "infant2yr", "infant12mo")
VARIANT = ""  # appended to the name of every figure saved while principal_heads() is in force

# manuscript figure -> drawn figures stacked top to bottom (names as the drawing scripts save them, plus VARIANT)
MAIN = {"Figure_1": ["Figure_R15_arrays"], "Figure_2": ["Figure_R12_geometry"], "Figure_3": ["Figure_R3_conditions"],
        "Figure_4": ["Figure_R17_noise_checks"], "Figure_5": ["Figure_R1_adult_depth"], "Figure_6": ["Figure_R11_maps_adult"],
        "Figure_7": ["Figure_R4_regions_adult"], "Figure_8": ["Figure_R6_pediatric_D_principal", "Figure_R13_maps_heads"],
        "Figure_9": ["Figure_helmet_fit_principal"], "Figure_10": ["Figure_R16_confirm_principal"]}
SUPP = {"Figure_S1": ["Figure_R0_sphere"], "Figure_S2": ["Figure_R2_noise_model"],
        "Figure_S3": ["Figure_covariance_structure"], "Figure_S4": ["Figure_covariance_bands"],
        "Figure_S5": ["Figure_noise_sensitivity", "Figure_noise_sensitivity_depth"], "Figure_S6": ["Figure_S_children_qc"],
        "Figure_S7": ["Figure_R8_depth_matched"], "Figure_S8": ["Figure_R7_helmet_fit"], "Figure_S9": ["Figure_R5_regions_heads"],
        "Figure_S10": ["Figure_R14_maps_scaled"], "Figure_S11": ["Figure_R9_noise_floor"], "Figure_S12": ["Figure_R10_spikes"],
        "Figure_S13": ["Figure_S_detection_curves_adult"], "Figure_S14": ["Figure_S_localization_effects"],
        "Figure_S15": ["Figure_S_joint_detection_localization"],
        "Figure_S16": ["Figure_R6_pediatric_D", "Figure_R13_maps_heads"], "Figure_S17": ["Figure_helmet_fit"],
        "Figure_S18": ["Figure_R16_confirm"]}

# print geometry (preamble.tex: A4, 2.4-cm side margins, 2.4/2.6-cm top/bottom), bp
TEXT_W, TEXT_H = 162.0 / 25.4 * 72.0, 247.0 / 25.4 * 72.0
MIN_PRINTED_PT = 7.0
# the documents' \includegraphics options per figure: width=\textwidth with, where given, height=<cap>\textheight,
# keepaspectratio; for a stacked figure shown one part per page (trim()), one cap per part. The caps of Figs. 3, 4, S2,
# S3, S5, S12 and S18 are those this export recommends (the .tex.j2 files must carry the same); the others are the
# documents' own.
PRINT = {"Figure_3": 0.72, "Figure_4": 0.72, "Figure_5": 0.62, "Figure_7": 0.66, "Figure_8": 0.6, "Figure_10": 0.62,
         "Figure_S2": 0.6, "Figure_S3": 0.72, "Figure_S5": (0.72, 0.75), "Figure_S12": 0.62, "Figure_S18": 0.88}
# figures whose legibility this export fixes (or checks), with the smallest printed type each must reach (per part where
# PRINT gives one cap per part): 7 pt, or 6.5 pt where the layout leaves no more room; the others keep the type the
# report drew them with
LEGIBLE = {"Figure_1": 7.0, "Figure_2": 7.0, "Figure_3": 7.0, "Figure_4": 7.0, "Figure_5": 6.5, "Figure_6": 6.5,
           "Figure_7": 7.0, "Figure_8": 7.0, "Figure_9": 7.0, "Figure_10": 7.0, "Figure_S1": 6.5, "Figure_S2": 7.0,
           "Figure_S3": 7.0, "Figure_S4": 7.0, "Figure_S5": (6.5, 7.0), "Figure_S6": 7.0, "Figure_S7": 6.5,
           "Figure_S8": 6.5, "Figure_S9": 6.5, "Figure_S10": 6.5, "Figure_S11": 6.5, "Figure_S12": 7.0,
           "Figure_S13": 6.5, "Figure_S14": 7.0, "Figure_S15": 7.0, "Figure_S16": 7.0, "Figure_S17": 7.0,
           "Figure_S18": 7.0}


def legend_texts(fig: Figure) -> set[int]:
    """ids of the texts of every legend of the figure (subfigures and axes included)."""
    legends = [a for a in fig.findobj(lambda o: isinstance(o, matplotlib.legend.Legend))]
    return {id(t) for lg in legends for t in lg.get_texts()} | {id(lg.get_title()) for lg in legends}


def relabel(fig: Figure) -> int:
    """'(a) Title' -> bold 'A' + title, for every text of the figure that starts with a panel label (legend entries
    keep a plain reference, '(C) ...', set by americanize())."""
    n = 0
    in_legend = legend_texts(fig)
    for t in fig.findobj(Text):
        if id(t) in in_legend:
            continue
        m = PANEL.match(t.get_text() or "")
        if not m:
            continue
        letter, rest = (m.group(1) or m.group(2)).upper(), m.group(3)
        t.set_text(rf"$\mathbf{{{letter}}}$" + (f"   {rest}" if rest else ""))
        n += 1
    return n


def add_letters(fig: Figure, letters: str) -> None:
    """Prefix '(g) ', '(h) ', ... to the left titles of the figure's axes, in reading order (top to bottom, then left
    to right); relabel() then sets them as bold capitals."""
    if not letters:
        return
    axes = [ax for ax in fig.axes if ax.get_visible() and ax.get_title(loc="left")]
    axes.sort(key=lambda ax: (-round(ax.get_position().y1, 2), round(ax.get_position().x0, 2)))
    if len(axes) != len(letters):
        raise SystemExit(f"{len(axes)} titled panels, {len(letters)} letters")
    for ax, letter in zip(axes, letters):
        ax.set_title(f"({letter}) {ax.get_title(loc='left')}", loc="left")


# British spellings in the figures' text -> the manuscript's American spelling (whole words; case of the first letter kept)
AMERICAN = {"centred": "centered", "centre": "center", "centres": "centers", "colour": "color", "colours": "colors",
            "coloured": "colored", "grey": "gray", "greys": "grays", "modelled": "modeled", "modelling": "modeling",
            "favour": "favor", "favours": "favors", "favoured": "favored", "favouring": "favoring",
            "favourable": "favorable", "normalised": "normalized", "normalisation": "normalization",
            "labelled": "labeled", "analysed": "analyzed", "analyse": "analyze", "behaviour": "behavior",
            "neighbour": "neighbor", "neighbours": "neighbors", "neighbouring": "neighboring", "summarised": "summarized",
            "summarise": "summarize", "summarises": "summarizes", "millimetre": "millimeter", "millimetres": "millimeters",
            "centimetre": "centimeter", "centimetres": "centimeters", "vapour": "vapor", "artefact": "artifact",
            "artefacts": "artifacts", "haemorrhage": "hemorrhage", "paediatric": "pediatric", "licence": "license",
            "optimised": "optimized", "minimised": "minimized", "maximised": "maximized", "recognised": "recognized",
            "characterised": "characterized", "visualised": "visualized", "organised": "organized"}
BRITISH = re.compile(r"\b(" + "|".join(AMERICAN) + r")\b", re.I)


# the analysis scripts' internal labels -> the manuscript's terms
TERMS = {"published model": "primary model", "room interference": "room field", "OPM dense": "dense OPM",
         "OPM matched": "site-matched OPM", "BEM (published)": "BEM (primary model)", "not stored": "not computed",
         "standard placements": "source-blind placements", "individual child": "school-aged child",
         "Fixed minus fitted helmet, within each head": "Fixed minus other helmet, within each head",
         "D(fixed) − D(fitted)": "D(fixed) − D(other)",
         ", thresholds frozen for a nominal 1 false event/min;": ";", "confirmatory": "pre-specified"}
# whole labels replaced as they stand (figure rows and legend entries)
EXACT = {"+ room field": "Sensor, brain and room noise",
         "+ room field, after the 8-term projection": "Sensor, brain and room noise, room field projected out",
         "1-layer BEM (1,000-target subset; no interval)": "Single-compartment BEM (1,000-target subset; no interval)",
         "3-layer BEM (primary model)": "Three-layer BEM (primary model)",
         "medial wall (not cortex)": "medial wall (not colored)",
         "within 4 mm of the inner skull (not simulated)": "within 4 mm of the inner skull (neither target nor background)",
         "primary 95 % interval, per array (shaded)": "primary 95 % interval (in each array's color)",
         "panels f and i": "panels F and I",
         # Fig. 3: the gap rows add a gap to the nominal OPM standoff (the primary condition has none beyond it)
         "OPM scalp gap 3 mm": "OPM scalp gap, additional 3 mm", "OPM scalp gap 6 mm": "OPM scalp gap, additional 6 mm",
         "30 fT/√Hz and a 3-mm gap (joint)": "30 fT/√Hz and an additional 3-mm gap (combined)",
         # Fig. S11: 30 fT/√Hz is the end of the noise sweep, a sensitivity setting, not a measured bound
         "≈ OPM floor,\nJas et al. 2026": "upper end\nof the sweep",
         # Fig. S5: the bands are parcel-bootstrap intervals of one anatomy, not confidence intervals
         "sensor + brain noise (band: 95 % CI)": "sensor + brain noise (band: 95 % parcel-bootstrap interval)",
         "frequency-resolved, flat signal spectrum (band: 95 % CI)":
             "frequency-resolved, flat signal spectrum (band: 95 % parcel-bootstrap interval)"}
UNIT_BRACKETS = re.compile(r"\[((?:fT|mm|dB|Hz|nAm|cm|ms|s)\b[^\]]{0,12})\]")
# a panel referred to inside a text, "(c)", "(b, c)", "(a-c)" or "(d and e)" -> upper case, as the panel letters are
PANEL_REF = re.compile(r"\(([a-l](?:(?:,\s*|\s*[-\u2013]\s*|\s+and\s+)[a-l])*)\)")


def americanize_text(s: str) -> str:
    """British spellings -> American; lower-case panel references -> upper case."""
    def swap(m):
        word = AMERICAN[m.group(1).lower()]
        return word[0].upper() + word[1:] if m.group(1)[0].isupper() else word
    s = BRITISH.sub(swap, s)
    for internal, term in TERMS.items():
        s = s.replace(internal, term).replace(internal[0].upper() + internal[1:], term[0].upper() + term[1:])
    s = UNIT_BRACKETS.sub(r"(\1)", s)
    s = PANEL_REF.sub(lambda m: "(" + re.sub(r"\b[a-l]\b", lambda k: k.group(0).upper(), m.group(1)) + ")", s)
    return EXACT.get(s, s)


def americanize(fig: Figure) -> int:
    """americanize_text() on every text of the figure, and on the fixed tick labels (which are redrawn from their
    formatter, so the formatter is rewritten)."""
    n = 0
    for t in fig.findobj(Text):
        s = t.get_text() or ""
        new = americanize_text(s)
        if new != s:
            t.set_text(new)
            n += 1
    for ax in fig.findobj(lambda a: isinstance(a, matplotlib.axes.Axes)):  # subfigures' axes included
        for axis in (ax.xaxis, ax.yaxis):
            for fmt in (axis.get_major_formatter(), axis.get_minor_formatter()):
                if isinstance(fmt, matplotlib.ticker.FixedFormatter):
                    fmt.seq = [americanize_text(str(x)) for x in fmt.seq]
                elif (isinstance(fmt, matplotlib.ticker.FuncFormatter) and isinstance(fmt.func, functools.partial)
                      and fmt.func.args and isinstance(fmt.func.args[0], dict)):  # set_ticklabels on fixed ticks
                    labels = fmt.func.args[0]
                    for k in labels:
                        labels[k] = americanize_text(str(labels[k]))
                elif isinstance(fmt, matplotlib.category.StrCategoryFormatter):  # categorical axis
                    fmt._units = {americanize_text(str(k)): v for k, v in fmt._units.items()}
    return n


def audit_text(fig: Figure, stem: str) -> list[str]:
    """After drawing: any British spelling or lower-case panel reference left in the figure's visible text."""
    left = []
    for t in fig.findobj(Text):
        s = t.get_text() or ""
        if t.get_visible() and s and americanize_text(s) != s:
            left.append(f"{stem}: {s[:120]!r}")
    return left


TEXT_AUDIT: list[str] = []


def strip_notes(fig: Figure) -> int:
    """Hide the explanatory notes that the report's figures carry below their panels (long, small, grey text): the
    journal puts that information in the caption. Axes left empty by it (axes used only to hold a note) are hidden too,
    so that the tight bounding box drops their space."""
    n = 0
    for t in fig.findobj(Text):
        s = (t.get_text() or "").strip()
        if not t.get_visible() or len(s) < 120 or t.get_fontsize() > 7.6:
            continue
        r, g, b, _ = matplotlib.colors.to_rgba(t.get_color())
        if abs(r - g) < 1e-6 and abs(g - b) < 1e-6 and r <= 0.5 and r > 0.0:
            t.set_visible(False)
            n += 1
    for ax in fig.axes:
        if ax.axison:
            continue
        shown = [c for c in ax.get_children() if c.get_visible() and c is not ax.patch
                 and not isinstance(c, (matplotlib.spines.Spine, matplotlib.axis.Axis))
                 and not (isinstance(c, Text) and not (c.get_text() or "").strip())]
        if not shown:
            ax.set_visible(False)
    return n


_orig_savefig = Figure.savefig
_orig_text_draw = Text.draw
FONTS: dict[str, list] = {}  # per saved figure: (size in pt, text) of every text drawn into its PDF
CLASH: dict[str, list] = {}  # per saved figure: the pairs of its texts that overlap or touch (clashes())
_DRAWN: list | None = None   # the list being filled while a PDF is written


_SEEN: list | None = None   # while clashes() draws a figure: every visible, non-empty text drawn, with its extent


def _recording_draw(self, renderer):
    """Text.draw during the export: notes the size of every visible, non-empty text drawn while a PDF is written (and
    the drawn extent of every text while clashes() looks for overlapping texts)."""
    shown = self.get_visible() and (self.get_text() or "").strip()
    if _DRAWN is not None and shown:
        _DRAWN.append((round(float(self.get_fontsize()), 3), self.get_text()))
    result = _orig_text_draw(self, renderer)
    if _SEEN is not None and shown:
        _SEEN.append((self, self.get_window_extent(renderer)))
    return result


def text_ink(t: Text, renderer, box=None) -> list:
    """A text's glyphs, one display box per line (the extent of its glyphs); a text at 90 degrees by its box."""
    from matplotlib.textpath import TextPath
    if round(t.get_rotation()) % 180:
        return [box if box is not None else t.get_window_extent(renderer)]
    k = t.get_figure(root=True).dpi / 72
    px, py = t.get_transform().transform(t.get_unitless_position())
    out = []
    for line, _, (x, y) in t._get_layout(renderer)[1]:
        if line.strip():
            e = TextPath((0, 0), line, prop=t.get_fontproperties()).get_extents()
            out.append(matplotlib.transforms.Bbox.from_extents(px + x + e.x0 * k, py + y + e.y0 * k,
                                                               px + x + e.x1 * k, py + y + e.y1 * k))
    return out


def clashes(fig: Figure, gap_pt: float = 1.0, tol_pt: float = 0.3) -> list[str]:
    """Pairs of texts drawn into the figure whose glyphs overlap, or come closer than gap_pt points side by side or one
    above the other (overlapping by more than tol_pt the other way); each line of a text is taken by the extent of its
    glyphs. The entries of one legend, set at its own spacing, count only when they overlap. Texts set at an angle
    other than 0 or 90 degrees are left out, and texts at 90 degrees are taken by their boxes."""
    global _SEEN
    _SEEN = []
    try:
        fig.canvas.draw()
        seen = [(t, b) for t, b in _SEEN if round(t.get_rotation()) % 90 == 0]
    finally:
        _SEEN = None
    legend_of = {id(t): id(lg) for lg in fig.findobj(lambda o: isinstance(o, matplotlib.legend.Legend))
                 for t in lg.get_texts() + [lg.get_title()]}
    renderer = fig.canvas.get_renderer()
    boxes = [(t, text_ink(t, renderer, b)) for t, b in seen]
    gap, tol = gap_pt * fig.dpi / 72, tol_pt * fig.dpi / 72
    out = []
    for i, (s, sa) in enumerate(boxes):
        for t, tb in boxes[i + 1:]:
            if t is s:  # an artist drawn twice
                continue
            g = 0.0 if legend_of.get(id(s), 0) == legend_of.get(id(t), 1) else gap
            for a in sa:
                for b in tb:
                    dx, dy = min(a.x1, b.x1) - max(a.x0, b.x0), min(a.y1, b.y1) - max(a.y0, b.y0)  # < 0: the gap
                    if (dx > -g and dy > tol) or (dy > -g and dx > tol):
                        out.append(f"{s.get_text()[:30]!r} / {t.get_text()[:30]!r}")
                        break
                else:
                    continue
                break
    return out


def savefig(self, fname, *args, **kwargs):
    """Every figure saved during the export: panel labels converted, spelling made American, notes removed, the
    legibility fixes of FIXES applied, then a vector PDF and a PNG preview in EXPORT (cropped to what is drawn). A figure
    saved while principal_heads() is in force gets the suffix VARIANT."""
    global _DRAWN
    base = Path(str(fname)).stem
    stem = base + VARIANT
    if base in RECOLOR:
        recolor(self, RECOLOR[base])
    if base == "Figure_R16_confirm":
        plain_p_values(self)
    if base in DROP_LETTER:
        for t in self.findobj(Text):
            txt = t.get_text() or ""
            if txt.strip() == DROP_LETTER[base]:
                t.set_visible(False)
            elif txt.startswith(DROP_LETTER[base] + " "):
                t.set_text(txt[len(DROP_LETTER[base]):].lstrip())
    add_letters(self, EXTRA_LETTERS.get(base, ""))
    relabel(self)
    americanize(self)
    strip_notes(self)
    if getattr(self, "_suptitle", None) is not None:  # the caption carries the figure's title, not the figure
        self._suptitle.set_visible(False)
    if base in FIXES:
        FIXES[base](self, principal=bool(VARIANT))
    CLASH[stem] = clashes(self)
    EXPORT.mkdir(parents=True, exist_ok=True)
    kwargs.pop("dpi", None)
    kwargs.update(bbox_inches="tight", pad_inches=0.03)
    _DRAWN = []
    try:
        _orig_savefig(self, EXPORT / f"{stem}.pdf", *args, dpi=PDF_DPI, **kwargs)
        FONTS[stem] = sorted(set(_DRAWN))
    finally:
        _DRAWN = None
    _orig_savefig(self, EXPORT / f"{stem}.png", *args, dpi=PNG_DPI, **kwargs)
    TEXT_AUDIT.extend(audit_text(self, stem))


def save(fig, name: str) -> Path:
    """report_style.save during the export: returns the preview PNG, whose pixel size at style.DPI (= PNG_DPI) some
    scripts check against the page width."""
    fig.savefig(EXPORT / f"{name}.pdf")
    plt.close(fig)
    return EXPORT / f"{name}{VARIANT}.png"


@contextlib.contextmanager
def principal_heads():
    """report_style.ANAT_ORDER restricted to PRINCIPAL (order kept) while the figures are drawn, which are saved with
    the suffix '_principal'."""
    global VARIANT
    order = style.ANAT_ORDER
    style.ANAT_ORDER = [k for k in order if k in PRINCIPAL]
    VARIANT = "_principal"
    try:
        yield
    finally:
        style.ANAT_ORDER, VARIANT = order, ""


# the 18-month template, which the MRI check found misregistered (its white surface fits its T1 image best after a
# 2.3-mm shift), is reported with children A-C in Section S10: in the figures of all nine anatomies (Figs. S9, S12 and
# S14-S18) it is drawn below the 12-month template, with the children, under a divider or class label APART_HEAD where
# the figure labels its classes; in the others (Figs. S7, S8 and S11) it keeps its place and is labelled MISREGISTERED
APART, APART_HEAD = "infant18mo", "Reported separately"
MISREGISTERED = " (misregistered)"


@contextlib.contextmanager
def misregistered():
    """While a figure is drawn: APART's full label (report_style.ANAT_LABEL, the label of keys and row labels) ends in
    MISREGISTERED. Nothing changes when APART is not drawn (principal_heads()) or its label already says so."""
    saved = style.ANAT_LABEL[APART]
    if APART in style.ANAT_ORDER and not saved.endswith(MISREGISTERED):
        style.ANAT_LABEL[APART] = saved + MISREGISTERED
    try:
        yield
    finally:
        style.ANAT_LABEL[APART] = saved


@contextlib.contextmanager
def set_apart():
    """While a figure of all nine anatomies is drawn: APART after the other templates, in the children's class (so that
    the drawing scripts group it with them; its marker stays a triangle, keep_marker()), whose class label becomes
    APART_HEAD, and labelled '(misregistered)' (misregistered()). Nothing changes when APART is not drawn
    (principal_heads())."""
    if APART not in style.ANAT_ORDER:
        yield
        return
    saved = style.ANAT_ORDER, style.ANAT_CLASS[APART], style.CLASS_LABEL["child"]
    order = [k for k in style.ANAT_ORDER if k != APART]
    i = next(i for i, k in enumerate(order) if style.ANAT_CLASS[k] == "child")
    style.ANAT_ORDER = order[:i] + [APART] + order[i:]
    style.ANAT_CLASS[APART], style.CLASS_LABEL["child"] = "child", APART_HEAD
    try:
        with misregistered():
            yield
    finally:
        style.ANAT_ORDER, style.ANAT_CLASS[APART], style.CLASS_LABEL["child"] = saved


def keep_marker(mod) -> None:
    """A drawing script's mstyle() (an anatomy's marker, by class) keeping the template's triangle for APART while
    set_apart() puts it in the children's class."""
    mstyle = mod.mstyle

    def kept(k, *args, **kw):
        cls = style.ANAT_CLASS[k]
        if k != APART or cls == "template":
            return mstyle(k, *args, **kw)
        style.ANAT_CLASS[k] = "template"
        try:
            return mstyle(k, *args, **kw)
        finally:
            style.ANAT_CLASS[k] = cls
    mod.mstyle = kept


@contextlib.contextmanager
def class_gaps(mod, extra: float, apart_extra: float = 0.0):
    """A drawing script's slots() (row positions: one unit per anatomy, `gap` units between classes) with `extra` units
    added to every gap between classes and `apart_extra` more before APART: room for the class headings."""
    slots = mod.slots

    def wider(keys, gap):
        pos, x, prev = {}, 0.0, None
        for k in keys:
            if prev is not None and style.ANAT_CLASS[k] != prev:
                x += gap + extra + (apart_extra if k == APART else 0.0)
            pos[k], x, prev = x, x + 1.0, style.ANAT_CLASS[k]
        return pos
    mod.slots = wider
    try:
        yield
    finally:
        mod.slots = slots


R6_GAP, R6_APART_GAP = 0.2, 0.45  # R6: units added between classes (headings and their rules) and before APART
CONFIRM_APART_GAP = 0.3           # R16: units added before APART (its divider and heading)


def pediatric():
    """report_figures_pediatric for the export: within each class the heads ordered as drawn before, by head
    circumference, except the children, in the order of their labels (A, B, C, as in the tables); R6 drawn with room for
    its class headings and, with all nine anatomies (Fig. S16), with APART set apart, as in R5 (Fig. S9, under the class
    label APART_HEAD) and R10 (Fig. S12, in the children's block of rows). R7-R9 (Figs. S7, S8, S11) keep APART among the
    templates, labelled '(misregistered)' (draw_all())."""
    ped = load("report_figures_pediatric")

    def order(s: dict) -> list[str]:
        ofc = {k: s["anatomies"][k]["head_size"]["ofc_mm"] for k in style.ANAT_ORDER}
        rank = {k: i for i, k in enumerate(style.ANAT_ORDER)}
        return sorted(style.ANAT_ORDER, key=lambda k: (ped.CLASSES.index(style.ANAT_CLASS[k]),
                                                       rank[k] if style.ANAT_CLASS[k] == "child" else -ofc[k]))
    ped.order = order
    keep_marker(ped)
    fig_r5, fig_r6, fig_r10 = ped.fig_r5, ped.fig_r6, ped.fig_r10

    def r5(s: dict, keys, *args) -> dict:
        with set_apart():
            return fig_r5(s, order(s), *args)

    def r6(s: dict, keys) -> dict:
        with set_apart(), class_gaps(ped, R6_GAP, R6_APART_GAP):
            return fig_r6(s, order(s))

    def r10(s: dict, keys) -> dict:
        with set_apart():
            return fig_r10(s, order(s))
    ped.fig_r5, ped.fig_r6, ped.fig_r10 = r5, r6, r10
    return ped


def draw_confirm() -> None:
    """R16 (Fig. S18; Fig. 10 under principal_heads(), its provenance beside the full figure's), APART set apart."""
    conf = load("report_figures_confirm")
    keep_marker(conf)
    if VARIANT:
        conf.OUT_JSON = "figures_confirm_principal.json"
    with set_apart(), class_gaps(conf, 0.0, CONFIRM_APART_GAP):
        conf.main(["--out-dir", str(EXPORT)])


# ----------------------------------------------------------------------------------------------
# legibility: helpers for the per-figure fixes (FIXES), applied at save time to the finished figure
def floor_fonts(fig: Figure, pt: float) -> None:
    """Every text of the figure (tick labels and legends included) at least `pt` points."""
    for t in fig.findobj(Text):
        if t.get_fontsize() < pt:
            t.set_fontsize(pt)
    for ax in fig.findobj(lambda a: isinstance(a, matplotlib.axes.Axes)):
        for name, axis in (("x", ax.xaxis), ("y", ax.yaxis)):
            for which, ticks in (("major", axis.get_major_ticks()), ("minor", axis.get_minor_ticks())):
                if ticks and ticks[0].label1.get_fontsize() < pt:  # ticks created at draw time copy the first one
                    ax.tick_params(axis=name, which=which, labelsize=pt)


def inches(fig: Figure, artist) -> matplotlib.transforms.Bbox:
    """The artist's drawn extent in inches from the figure's lower left corner."""
    return artist.get_window_extent(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())


def place(ax, x0: float, y0: float, w: float, h: float) -> None:
    """Axes position in inches."""
    W, H = ax.get_figure(root=True).get_size_inches()
    ax.set_position([x0 / W, y0 / H, w / W, h / H])


def plotted_boxes(ax, skip=(), guides: bool = False) -> list:
    """Display boxes of what the axes draws: one per marker (its size and edge), per line segment (its width), per
    patch, per segment of a line collection (error bars), and per text other than those in `skip`. Guide lines across
    the axes (axhline, axvline: not in data coordinates) count only with guides=True."""
    import numpy as np
    fig = ax.get_figure(root=True)
    r, k = fig.canvas.get_renderer(), fig.dpi / 72
    Bbox = matplotlib.transforms.Bbox
    out = []

    def segments(xy, w):
        return [Bbox.from_extents(min(a[0], b[0]) - w, min(a[1], b[1]) - w, max(a[0], b[0]) + w, max(a[1], b[1]) + w)
                for a, b in zip(xy[:-1], xy[1:])]
    for ln in ax.get_lines():
        if not ln.get_visible() or (ln.get_transform() is not ax.transData and not guides):
            continue
        xy = ln.get_transform().transform(np.asarray(ln.get_xydata(), float))
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if ln.get_marker() not in ("None", "none", "", " ", None):
            m = (ln.get_markersize() + ln.get_markeredgewidth()) / 2 * k
            out += [Bbox.from_extents(x - m, y - m, x + m, y + m) for x, y in xy]
        if ln.get_linestyle() not in ("None", "none", "", " ") and len(xy) > 1:
            out += segments(xy, ln.get_linewidth() / 2 * k)
    for p in ax.patches:
        if p.get_visible():
            b, w = p.get_window_extent(r), (p.get_linewidth() / 2 * k if p.get_edgecolor()[3] > 0 else 0.0)
            out.append(Bbox.from_extents(b.x0 - w, b.y0 - w, b.x1 + w, b.y1 + w))
    for c in ax.collections:
        if not c.get_visible():
            continue
        if isinstance(c, matplotlib.collections.LineCollection):
            w = max(c.get_linewidths(), default=1.0) / 2 * k
            for seg in c.get_segments():
                out += segments(c.get_transform().transform(np.asarray(seg, float)), w)
        else:
            out.append(c.get_window_extent(r))
    for t in ax.texts:
        if t.get_visible() and (t.get_text() or "").strip() and all(t is not s for s in skip):
            out += text_ink(t, r)
    return out


def overlaps(a, b, pad: float = 0.0) -> bool:
    """Two display boxes closer than pad (display units) or overlapping."""
    return a.x0 - pad < b.x1 and b.x0 - pad < a.x1 and a.y0 - pad < b.y1 and b.y0 - pad < a.y1


def room_above(ax, texts, pad_pt: float = 2.0) -> None:
    """The axes' upper y limit (linear axis) raised until the given texts, set at its top in axes coordinates, clear
    what the axes draws below them by pad_pt points."""
    fig = ax.get_figure(root=True)
    pad = pad_pt * fig.dpi / 72
    for _ in range(12):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        ink = [b for t in texts for b in text_ink(t, r)]
        under = [b for b in plotted_boxes(ax, skip=texts) if any(overlaps(b, i, pad) for i in ink)]
        if not under:
            return
        lo, hi = ax.get_ylim()
        need = max(b.y1 for b in under) + pad - min(i.y0 for i in ink)  # display units the data must come down
        ax.set_ylim(lo, hi + 1.25 * need / ax.get_window_extent(r).height * (hi - lo))
    raise SystemExit(f"room_above: no room for {texts[0].get_text()!r}")


def label_width(fig: Figure, texts) -> float:
    """Widest of the given (visible, non-empty) texts, inches."""
    fig.canvas.draw()
    return max((inches(fig, t).width for t in texts if t.get_visible() and (t.get_text() or "").strip()), default=0.0)


def rebuild_legend(old, handles=None, labels=None, **kw):
    """A legend replaced by one with the same entries (or those given) and the given options; returns the new one."""
    handles = list(old.legend_handles) if handles is None else handles
    labels = [t.get_text() for t in old.get_texts()] if labels is None else labels
    parent = old.axes if old.axes is not None else old.get_figure(root=False)
    title = old.get_title().get_text()
    old.remove()
    if title:
        kw.setdefault("title", title)
    return parent.legend(handles, labels, **kw)


def fit_width(fig: Figure, gridspecs, width_in: float) -> None:
    """Shift the left and right edges of the given gridspecs so that what is drawn spans width_in inches (the saved
    figure, its 0.03-in pad included, is then 0.06 in wider)."""
    for _ in range(3):
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        W = fig.get_size_inches()[0]
        d_left, d_right = 0.0 - bb.x0, width_in - bb.x1  # wanted: x0 at 0, x1 at width_in
        if abs(d_left) < 0.005 and abs(d_right) < 0.005:
            return
        for gs in gridspecs:
            gs.update(left=gs.left + d_left / W, right=gs.right + d_right / W)
        for ax in fig.axes:
            ss = ax.get_subplotspec()
            if ss is not None:
                ax.set_position(ss.get_position(fig))


def fit_constrained(fig: Figure, width_in: float, height_in: float) -> None:
    """A figure laid out by the constrained layout engine, sized so that what is drawn spans width_in inches (the
    saved figure, its 0.03-in pad included, is then the text width when width_in is DRAWN_W) and height_in high."""
    W = width_in + 0.08  # the engine leaves its pad (3 pt) at either side
    for _ in range(4):
        fig.set_size_inches(W, height_in)
        fig.canvas.draw()
        drawn = fig.get_tightbbox(fig.canvas.get_renderer()).width
        if abs(drawn - width_in) < 0.003:
            return
        W += width_in - drawn


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def run_main(name: str, argv: list[str]):
    mod = load(name)
    old = sys.argv
    sys.argv = [name] + argv
    try:
        mod.main()
    finally:
        sys.argv = old


# the scaled adults get colors no array uses (the report drew them in the arrays' blues); labels as in Table 2
HEAD_COLORS = {"school": "#332288", "size2yr": "#882255"}
HEAD_SHORT = {"school": "School-age size", "size2yr": "2-year size"}
# per figure: colors replaced in every artist ({old: new}, lower-case hex); Fig. 7's fitted helmet as in Fig. 1
RECOLOR: dict[str, dict] = {}


def recolor(fig: Figure, mapping: dict) -> None:
    """Replace colors (as hex) in lines, markers and line collections of the figure."""
    to_hex = matplotlib.colors.to_hex
    for a in fig.findobj():
        for get, put in (("get_color", "set_color"), ("get_markerfacecolor", "set_markerfacecolor"),
                         ("get_markeredgecolor", "set_markeredgecolor")):
            if isinstance(a, matplotlib.collections.Collection) or not (hasattr(a, get) and hasattr(a, put)):
                continue
            try:
                c = getattr(a, get)()
                if isinstance(c, str) and c not in ("none", "None") and to_hex(c) in mapping:
                    getattr(a, put)(mapping[to_hex(c)])
            except (ValueError, TypeError):
                pass
        if isinstance(a, matplotlib.collections.Collection):
            for get, put in (("get_edgecolor", "set_edgecolor"), ("get_facecolor", "set_facecolor")):
                cols = getattr(a, get)()
                if len(cols):
                    new = [mapping.get(to_hex(c, keep_alpha=False), None) for c in cols]
                    if any(new):
                        getattr(a, put)([n if n else tuple(c) for n, c in zip(new, cols)])


def setup():
    """The drawing scripts' output, saving and look redirected to the export."""
    Figure.savefig = savefig
    Text.draw = _recording_draw
    style.OUT = EXPORT
    style.save = save
    style.DPI = PNG_DPI  # the preview's resolution; the PDFs' raster layers use PDF_DPI
    style.ANAT_COLOR.update(HEAD_COLORS)
    style.ANAT_SHORT.update(HEAD_SHORT)


def draw_noise_checks():
    """Fig. 4 (R17): report_figures_noise draws it at its page width (175 mm, checked on the saved preview) and
    HEIGHT_IN high, with its type sizes FS_* (rows, group headers, notes). Printed at the text width (0.926 of 175 mm) under
    a height cap of 0.72, it needs type of 7.6 pt (fix_noise_checks raises the rest) and 7.55 in instead of 7.7. (Drawn
    at the text width instead, the forest plots of B and C, whose row labels and group headers need the width the report
    gives them, are squeezed by the constrained layout.) Scenario C's gloss in (B), the report's 'pessimistic', reads
    'least favorable to the OPM', on a second line."""
    noise = load("report_figures_noise")
    noise.HEIGHT_IN = 7.55
    noise.FS_ROW = noise.FS_HEAD = noise.FS_NOTE = 7.6
    # (B): scenario C is the one least favorable to the OPM (the report's gloss: 'pessimistic'), on two lines
    glosses = {letter: gloss for letter, _, gloss in noise.SCENARIOS}
    if glosses.get("C") != "OPM bears the excess (pessimistic)":
        raise SystemExit(f"Figure_R17_noise_checks: scenario C's gloss is no longer the one replaced: {glosses.get('C')!r}")
    noise.SCENARIOS = tuple((letter, key, "OPM bears the excess\n(least favorable to the OPM)" if letter == "C" else gloss)
                            for letter, key, gloss in noise.SCENARIOS)
    noise.main()


REPORT_W = style.FULL_W  # the report's figure width (inches), as report_style sets it


@contextlib.contextmanager
def full_width(width_in: float):
    """report_style.FULL_W set to width_in while a figure is drawn whose script lays it out from it."""
    old = style.FULL_W
    style.FULL_W = width_in
    try:
        yield
    finally:
        style.FULL_W = old


def draw_clean():
    """report_figures_clean, with Figs. 1 (R15) and 2 (R12) drawn at the text width: their script sets their panels'
    width and height from report_style.FULL_W at the type sizes it gives, so the panels shrink to the text width and the
    type keeps its size (7 pt printed instead of 6.2; a little added for the margins their layout leaves out of what it
    draws)."""
    clean = load("report_figures_clean")

    def at_text_width(draw, extra):
        def drawn():
            with full_width(DRAWN_W + extra):
                return draw()
        return drawn

    clean.figure_arrays = at_text_width(clean.figure_arrays, 0.06)
    clean.figure_geometry = at_text_width(clean.figure_geometry, 0.1)
    old = sys.argv
    sys.argv = ["report_figures_clean"]
    try:
        clean.main()
    finally:
        sys.argv = old


def draw_pediatric_principal():
    """Fig. 8 (top): R6 with the principal heads. report_figures_pediatric.main() cannot draw with them: it draws every
    figure, then stops in provenance(), whose caption draft of R6 (and alt text of R8) takes the minimum over the
    individual children, of whom there are none (ValueError: min() iterable argument is empty). main()'s set-up is
    repeated here and fig_r6 alone is called (its own checks run)."""
    ped = pediatric()
    style.apply()
    s = ped.load(ped.G3B)
    keys = ped.order(s)
    ped.CM.update({k: s["anatomies"][k]["head_size"]["ofc_mm"] / 10 for k in keys})
    ped.fig_r6(s, keys)


def draw_helmet_fit() -> None:
    """Fig. S17 (APART set apart) and, under principal_heads(), Fig. 9."""
    summary = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())
    g3b = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())
    with set_apart():
        helmet_fit_figure(summary, g3b, EXPORT / "Figure_helmet_fit.png")


def draw_all():
    setup()
    run_main("report_figures_adult", [])
    draw_clean()
    draw_noise_checks()
    with misregistered():
        pediatric().main()  # R5-R10: Figs. S7-S9, S11, S12 and the top of Fig. S16
    draw_confirm()
    with principal_heads():
        draw_pediatric_principal()
        draw_confirm()
    run_main("report_figures_qc", [])
    with set_apart():  # Figs. S14 and S15 (and S13, the adult's alone)
        run_main("report_figures_supplement", [])
    draw_helmet_fit()
    with principal_heads():
        draw_helmet_fit()
    cv = load("study_covariance_validation")
    cvs = json.loads((ROOT / "results/g2_covariance_validation/covariance_validation.json").read_text())
    cv.figure_structure(cvs, EXPORT / "Figure_covariance_structure.png")
    cv.figure_bands(cvs, EXPORT / "Figure_covariance_bands.png")
    ns = load("study_noise_sensitivity")
    nss = json.loads((ROOT / "results/g2_noise_sensitivity/noise_sensitivity_summary.json").read_text())
    ns.figures(nss, EXPORT)
    (EXPORT / "font_audit.json").write_text(json.dumps({k: v for k, v in sorted(FONTS.items())}, ensure_ascii=False,
                                                       indent=0) + "\n")
    (EXPORT / "clash_audit.json").write_text(json.dumps(dict(sorted(CLASH.items())), ensure_ascii=False,
                                                        indent=0) + "\n")


def helmet_fit_figure(s: dict, g3b: dict, path: Path) -> None:
    """The manuscript's helmet-fit figure, drawn from the stored results of the helmet-fit study (the report drew the same
    quantities in four panels with a third helmet): (a) the Neuromag gap in the fixed adult helmet at top contact and in
    the helmet fitted at the adult's gap; (b) Delta against the adult placed by the same rule in each helmet, with the
    fixed helmet's range over its source-blind placements; (c) the interaction, the primary helmet-fit quantity (a
    head's fixed-minus-fitted contrast minus the adult's), which does not depend on the OPM noise. Dense OPM array
    against Neuromag's 306 channels, sensor plus brain noise; the heads reported separately (Section S10) are set apart. The
    heads are those of report_style.ANAT_ORDER (all nine: Fig. S17; principal_heads(): Fig. 9). Drawn at the text width,
    no text under 7.2 pt; bars and caps 1.2 pt wide (thinner, a bar can print as a pale hairline). Under set_apart()
    the heads right of the dotted line, APART with the children, are labelled APART_HEAD in each panel (the upper limit
    raised for the label where needed) and APART is drawn with open markers, keyed '18-month template
    (misregistered)'."""
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Patch
    style.apply()
    cond, ref, o = "intrinsic+brain", "combined", "opm_dense"
    keys = [k for k in style.ANAT_ORDER if k in s["helmets"]]
    kids = [k for k in keys if k != "adult"]
    xs = {k: i for i, k in enumerate(keys)}
    fixed = dict(marker="s", color="#000000", mfc="#000000")
    fitted = dict(marker="o", color="#D55E00", mfc="#D55E00")
    apart = APART in keys and style.ANAT_CLASS[APART] == "child"  # set_apart() in force

    def face(k, st):  # APART, set apart: open markers
        return "white" if apart and k == APART else st["mfc"]
    fig, axs = plt.subplots(1, 3, figsize=(DRAWN_W, 3.4),
                            gridspec_kw=dict(wspace=0.45, left=0.085, right=0.99, top=0.9, bottom=0.36))
    heads = []  # the APART_HEAD labels, one per panel

    def xaxis(ax, ks, zero=True):
        ax.set_xticks([xs[k] for k in ks], [style.ANAT_SHORT[k] for k in ks], rotation=45, ha="right")
        ax.set_xlim(min(xs[k] for k in ks) - 0.6, max(xs[k] for k in ks) + 0.6)
        if zero:
            ax.axhline(0, color="0.75", lw=0.6, zorder=0)
        first_child = min((xs[k] for k in ks if style.ANAT_CLASS[k] == "child"), default=None)
        if first_child is not None:
            ax.axvline(first_child - 0.5, color="0.55", lw=0.6, ls=":", zorder=0)
            if apart:
                heads.append(ax.text(first_child - 0.5, 1.0, APART_HEAD.replace(" ", "\n"), ha="left", va="top",
                                     fontsize=7.2, color="0.45", style="italic",
                                     transform=offset_copy(ax.get_xaxis_transform(), fig=fig, x=2.5, y=-1.0,
                                                           units="points")))

    def point(ax, x, e, st, dx, k):
        med = e["median"]
        lo, hi = e.get("ci95") or (med, med)
        ax.errorbar([x + dx], [med], yerr=[[med - lo], [hi - med]], fmt=st["marker"], color=st["color"],
                    mfc=face(k, st), ms=4.2, lw=1.2, capsize=1.8, capthick=1.2)

    ax = axs[0]
    for h, st, dx in (("top", fixed, -0.13), ("gap_matched", fitted, 0.13)):
        for k in keys:
            v = s["helmets"][k][h]
            ax.plot([xs[k] + dx], [v["median_mm"]], st["marker"], color=st["color"], mfc=face(k, st), ms=4.2, ls="none")
            if h == "gap_matched" and v.get("clearance_binding"):
                ax.annotate("*", (xs[k] + dx, v["median_mm"]), xytext=(3, 1), textcoords="offset points", fontsize=9,
                            color=st["color"])
    target = s["target_gaps"]["gap_matched"]["gap_mm"]
    ax.axhline(target, color="#D55E00", lw=0.7, ls="--", zorder=0)
    ax.set_ylabel("median magnetometer-to-scalp gap (mm)")
    ax.set_title("(a)  Neuromag gap", loc="left")
    xaxis(ax, keys, zero=False)

    ax = axs[1]
    band = s.get("placement_band", {})
    for k in kids:
        b = band.get(f"{k}/{ref}/{cond}")
        if b:
            ax.add_patch(plt.Rectangle((xs[k] - 0.38, b["min"]), 0.76, b["max"] - b["min"], fc="0.88", ec="none", zorder=0))
    for k in kids:
        point(ax, xs[k], g3b["comparisons"][f"{k}/{o}/{ref}/{cond}/detect"]["delta"], fixed, -0.13, k)
        point(ax, xs[k], s["delta_same_rule"][f"{k}/gap_matched/{o}/{ref}/{cond}"]["delta"], fitted, 0.13, k)
    ax.set_ylabel("Δ, paired change from the adult (dB)")
    ax.set_title("(b)  Change from the adult", loc="left")
    xaxis(ax, kids)

    ax = axs[2]
    for k in kids:
        point(ax, xs[k], s["interaction"][f"{k}/gap_matched/{ref}/{cond}"]["interaction"], fitted, 0.0, k)
    ax.set_ylabel("interaction (dB)")
    ax.set_title("(c)  Interaction", loc="left")
    xaxis(ax, kids)

    handles = [Line2D([], [], ls="none", ms=4.5, label="fixed adult helmet, top contact", **fixed),
               Line2D([], [], ls="none", ms=4.5, label="helmet fitted at the adult's gap", **fitted)]
    if apart:
        handles.append((Line2D([], [], ls="none", ms=4.5, **dict(fixed, mfc="white")),
                        Line2D([], [], ls="none", ms=4.5, **dict(fitted, mfc="white"))))
    handles += [Line2D([], [], color="#D55E00", ls="--", lw=0.8),
                Patch(fc="0.88", ec="none")]
    labels = ["fixed adult helmet, top contact", "helmet fitted at the adult's gap"]
    labels += [f"open: {style.ANAT_LABEL[APART]}"] if apart else []
    labels += [f"(a) the adult's gap ({target:.1f} mm)",
               "(b) fixed helmet over its source-blind placements\n(difference of medians, not paired)"]
    fig.legend(handles, labels, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0), fontsize=7.2,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)})
    if apart:  # each panel's label clear of what it plots: the upper limit raised where needed
        for t in heads:
            room_above(t.axes, [t], pad_pt=2.0)
    fig.savefig(path)
    plt.close(fig)


def plain_p_values(fig: Figure) -> None:
    """Remove the check marks (and the bold) that the report's spike figure puts on significant p values."""
    for t in fig.findobj(Text):
        s = t.get_text() or ""
        if "\u2713" in s:
            t.set_text(s.replace(" \u2713", "").replace("\u2713", "").rstrip())
            t.set_fontweight("normal")


# ----------------------------------------------------------------------------------------------
# legibility fixes, per drawn figure (FIXES, by the name the drawing script saves it under); each is applied to the
# finished figure at save time and changes only its size, layout, type sizes and keys, never what is plotted
DRAWN_W = TEXT_W / 72.0 - 0.06  # inches: what is drawn, so that the saved figure (0.03-in pad each side) is the text width
DOTTED = (0, (1.0, 1.4))         # report_figures_noise's line style of the projected condition


def square_keys(fig: Figure) -> None:
    """Principal versions: the generic keys (filled/open, dark/light) drawn with the children's circle get a square,
    a marker the figure shows."""
    for lg in fig.legends + [a.get_legend() for a in fig.axes if a.get_legend() is not None]:
        for h in lg.legend_handles:
            if isinstance(h, Line2D) and h.get_marker() == "o":
                h.set_marker("s")


def fix_conditions(fig: Figure, principal: bool = False) -> None:
    """Fig. 3 (R3): drawn at its printed size (the text width, 7.0 in high: height cap 0.72), no text under 7.2 pt,
    laid out by hand (constrained layout cannot: the group headers and the axis label are longer than the forest plot is
    wide at this width, which it answers by squeezing the plot). The headers 'dense' and 'site-matched' apart: each value
    column starts a fixed distance into the value axes (2 and 36 pt) and 'site-matched' is set on two lines, broken at
    its hyphen, bottom-aligned with 'dense'. The ratio axis is labelled at 0.5-0.8, 1, 1.2 and 1.4, with 0.9, 1.1 and
    1.3 as unlabelled ticks; its label on two lines. The dB axis is checked again at the new size."""
    adult = sys.modules["report_figures_adult"]
    ax, tw = next((p, t) for p, t, _, _ in adult.DB_AXES if p.get_figure(root=True) is fig)
    vx = next(a for a in fig.axes if not a.axison)
    major = (0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.4)
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(major))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter([f"{t:g}" for t in major]))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0.9, 1.1, 1.3]))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="x", which="minor", length=2.0, width=0.6)
    heads = [t for t in vx.texts if t.get_text() in ("dense", "site-matched")]
    if len(heads) != 2:
        raise SystemExit("Figure_R3_conditions: value column headers not found")
    for t in vx.texts:  # values and headers: dense column at x = 0.06, site-matched at 0.5 (axes fraction) as drawn
        dx = 2.0 if t.get_position()[0] < 0.3 else 36.0
        t.set_x(0.0)
        t.set_transform(offset_copy(vx.get_yaxis_transform(), fig=fig, x=dx, units="points"))
    for t in heads:
        t.set_text(t.get_text().replace("site-", "site-\n"))
        t.set_verticalalignment("bottom")
        t.set_y(0.45)
    ax.set_xlabel(ax.get_xlabel().replace(" (log scale", "\n(log scale"))
    floor_fonts(fig, 7.2)
    fig.set_layout_engine("none")
    W, H = DRAWN_W, 7.0
    fig.set_size_inches(W, H)
    groups = [t for t in ax.texts if t.get_fontweight() == "bold"]  # the group headers, drawn from the plot's left edge
    if len(groups) != 3:
        raise SystemExit(f"Figure_R3_conditions: {len(groups)} group headers, 3 expected")
    for t in groups:  # longer than the plot is wide: set as section heads in the row-label column, from its left edge
        t.set_transform(matplotlib.transforms.blended_transform_factory(fig.transFigure, ax.transData))
        t.set_x(0.01 / W)
    old = fig.legends[0]
    lg = rebuild_legend(old, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=old.get_texts()[0].get_fontsize(),
                        title_fontsize=old.get_title().get_fontsize())
    fig.canvas.draw()
    lab = label_width(fig, ax.get_yticklabels()) + (ax.yaxis.get_tick_padding() + ax.yaxis.majorTicks[0].get_pad()) / 72
    need = max(inches(fig, t).x1 for t in vx.texts) - inches(fig, vx).x0 + 0.03  # the value columns' width
    gap = 0.05
    L, X = lab + 0.01, W - 0.01 - lab - gap - need
    top, bottom = 0.4, 1.0  # inches above and below the plots, then set from what is drawn there
    for _ in range(3):
        place(ax, L, bottom, X, H - top - bottom)  # the dB axis is twinned: it follows
        place(vx, L + X + gap, bottom, need, H - top - bottom)
        fig.canvas.draw()
        bottom += inches(fig, lg).y1 + 0.1 - inches(fig, ax.xaxis.label).y0  # axis label 0.1 in above the key
        top += inches(fig, tw.xaxis.label).y1 - (H - 0.01)
    adult.check_db_axes(fig)
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    if bb.width > DRAWN_W + 0.01 or max(inches(fig, t).x1 for t in vx.texts) > inches(fig, vx).x1 + 0.005:
        raise SystemExit(f"Figure_R3_conditions: {bb.width:.3f} in drawn, or the value columns overrun their axes")


def fix_noise_checks(fig: Figure, principal: bool = False) -> None:
    """Fig. 4 (R17), its type raised to 7.6 pt (draw_noise_checks sets report_figures_noise's FS_* and HEIGHT_IN; the
    rest here): in (D) the break-even intervals after the projection are dotted, the projection's line style, and the
    open diamonds get their own key entry beside the filled ones; the break-even values, which the report wrote beside
    their diamonds among the curves, are given instead in a small key of the diamonds in the empty lower left corner of
    (D), checked to cover no plotted point. The dB axes are checked again."""
    axd = next(a for a in fig.findobj(lambda o: isinstance(o, matplotlib.axes.Axes))
               if a.get_legend() is not None and any("break-even" in t.get_text() for t in a.get_legend().get_texts()))
    n = 0
    for ln in axd.get_lines():  # interval bars: two points, 2.4 pt; the projected ones just below a ratio of 1
        if ln.get_linewidth() == 2.4 and len(ln.get_xdata()) == 2 and ln.get_ydata()[0] < 1:
            ln.set_linestyle(DOTTED)
            n += 1
    if n != 2:
        raise SystemExit(f"Figure_R17_noise_checks: {n} projected break-even intervals found, 2 expected")
    old = axd.get_legend()
    handles, labels = list(old.legend_handles), [t.get_text() for t in old.get_texts()]
    i = next(i for i, lab in enumerate(labels) if lab.startswith("break-even"))

    def key(filled: bool):
        bar = Line2D([], [], color="0.25", lw=2.4, alpha=0.6, solid_capstyle="butt", ls="-" if filled else DOTTED)
        return bar, Line2D([], [], ls="none", marker="D", ms=4.4, mec="0.25", mfc="0.25" if filled else "white", mew=1.1)

    handles[i:i + 1] = [key(True), key(False)]
    labels[i:i + 1] = ["break-even level, 95 % interval:\nsensor plus brain noise (filled)", "after the projection (open, dotted)"]
    rebuild_legend(old, handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=7.6, handlelength=2.4,
                   borderaxespad=0.0, labelspacing=0.45, ncol=1)
    # (D): the break-even values, written at their diamonds where the curves cross a ratio of 1, in a key of the
    # diamonds (one column per array, dense first; filled above open) in the lower left corner, where nothing is plotted
    diamonds = {(float(ln.get_xdata()[0]), float(ln.get_ydata()[0])): ln
                for ln in axd.get_lines() if ln.get_marker() == "D"}
    values = [t for t in axd.texts if re.fullmatch(r"\d+\.\d", t.get_text())]
    if len(values) != 4 or any(tuple(map(float, t.xy)) not in diamonds for t in values):
        raise SystemExit("Figure_R17_noise_checks: four break-even values, each at its diamond, expected in (D)")
    order = [style.ARRAY_COLOR[a] for a in ("opm_dense", "opm_matched")]
    values.sort(key=lambda t: (order.index(t.get_color()), t.xy[1] < 1))
    marks = [diamonds[tuple(map(float, t.xy))] for t in values]
    key = matplotlib.legend.Legend(
        axd, [Line2D([], [], ls="none", marker="D", ms=d.get_markersize(), mew=d.get_markeredgewidth(),
                     mec=d.get_markeredgecolor(), mfc=d.get_markerfacecolor()) for d in marks],
        [t.get_text() for t in values], loc="lower left", ncol=2, fontsize=7.6, title="break-even level",
        title_fontsize=7.6, alignment="left", handlelength=0.9, handletextpad=0.3, columnspacing=1.0, borderpad=0.2,
        labelspacing=0.4, borderaxespad=0.3)
    axd.add_artist(key)
    for t in values:
        t.remove()
    # (A): its key, below the three rows, would cover the last row's label at 7.6 pt: more empty room below the rows
    axa = next(a for a in fig.findobj(lambda o: isinstance(o, matplotlib.axes.Axes))
               if a.get_legend() is not None and any("measured noise" in t.get_text() for t in a.get_legend().get_texts()))
    lo, hi = axa.get_ylim()
    axa.set_ylim(lo + 0.9, hi)
    rebuild_legend(axa.get_legend(), loc="lower right", fontsize=7.6, handlelength=1.7, borderaxespad=0.15,
                   labelspacing=0.45, frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    # group headers of (B) and (C), drawn from the plot's left edge: at 7.6 pt longer than the plots are wide, so set as
    # section heads in the row-label column, from the panel's left edge
    for ax in fig.findobj(lambda o: isinstance(o, matplotlib.axes.Axes)):
        for t in ax.texts:
            if t.get_fontweight() == "bold":
                sf = ax.get_figure(root=False)
                t.set_transform(matplotlib.transforms.blended_transform_factory(sf.transSubfigure, ax.transData))
                t.set_x(0.012)
    floor_fonts(fig, 7.6)
    sys.modules["report_figures_adult"].check_db_axes(fig)
    if legend_hits(axd, key, full=True):
        raise SystemExit(f"Figure_R17_noise_checks: D's key of the break-even values covers "
                         f"{legend_hits(axd, key, full=True)} plotted points")


ROW_IN = 0.24  # R6: inches per row unit (the report's 2.24 in for the 9.4 units of nine anatomies)


def class_headings(fig: Figure, a1, a2) -> float | None:
    """R6's class headings: each class's label (A's label column) and the name of its Δ estimator (B, from its left
    edge, on white so that the zero line does not cross it) set on one line in the gap above the class's first row,
    underlined by the class's rule (light gray, across both panels), which lies between the labels and the markers.
    Spacing in points at ROW_IN inches per row unit: rule 2.5 pt above the highest marker of the class, labels 1.5 pt
    above the rule and at least 2 pt below the class above. The upper limit is raised for the first class's heading.
    Returns the height of the divider above the class headed APART_HEAD (2.5 pt above its labels), or None."""
    import numpy as np
    pt = 1 / (ROW_IN * 72)  # row units per point
    labels = sorted(a1.texts, key=lambda t: -t.get_position()[1])
    names = sorted(a2.texts, key=lambda t: -t.get_position()[1])
    if [round(t.get_position()[1], 6) for t in labels] != [round(t.get_position()[1], 6) for t in names] or not labels:
        raise SystemExit("Figure_R6_pediatric_D: class labels and estimator names do not pair up")
    for ax in (a1, a2):  # the report's rules, half-way between the classes, through the labels: removed
        for ln in [ln for ln in ax.get_lines() if ln.get_transform() is not ax.transData
                   and list(ln.get_xdata()) == [0, 1] and ln.get_ydata()[0] == ln.get_ydata()[1]]:
            ln.remove()
    rows = sorted((float(y) for y in a1.get_yticks()), reverse=True)
    top, bottom = {}, {}  # highest and lowest ink of each row (markers and interval lines of both panels)
    for ax in (a1, a2):
        for ln in ax.get_lines():
            if ln.get_transform() is not ax.transData:
                continue
            marked = ln.get_marker() not in ("None", "none", "", " ", None)
            ext = ((ln.get_markersize() + ln.get_markeredgewidth()) / 2 if marked else ln.get_linewidth() / 2) * pt
            for y in np.asarray(ln.get_ydata(), float):
                row = min(rows, key=lambda r: abs(r - y))
                top[row], bottom[row] = max(top.get(row, -np.inf), y + ext), min(bottom.get(row, np.inf), y - ext)
    fig.canvas.draw()
    divider, lo, hi = None, *a1.get_ylim()
    for lab, name in zip(labels, names):  # the class's first row: the report set its labels half a unit above it
        first = min(rows, key=lambda r: abs(r - (lab.get_position()[1] - 0.5)))
        above = [r for r in rows if r > first + 0.5]
        rule = top[first] + 2.5 * pt
        base = rule + 1.5 * pt  # the labels' lower edge (descenders included)
        height = max(inches(fig, t).height for t in (lab, name)) * 72 * pt
        lab.set_y(base)
        name.set_position((0.0, base))
        name.set_transform(offset_copy(a2.get_yaxis_transform(), fig=fig, x=3.0, units="points"))
        name.set_horizontalalignment("left")
        name.set_bbox(dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
        ceiling = base + height + (2.5 * pt if lab.get_text() == APART_HEAD else 0.0)
        if lab.get_text() == APART_HEAD:
            divider = ceiling
        if above and ceiling + 2.0 * pt > bottom[min(above)]:
            raise SystemExit(f"Figure_R6_pediatric_D: no room for the heading '{lab.get_text()}' (widen R6_GAP)")
        hi = max(hi, ceiling + 1.0 * pt)
        for ax in (a1, a2):
            ax.axhline(rule, color="0.88", lw=0.7, zorder=0)
    a1.set_ylim(lo, hi)
    return divider


def fix_pediatric_d(fig: Figure, principal: bool = False) -> None:
    """R6 (top of Figs. 8 and S16): drawn at the text width, no text under 7.2 pt; the height follows the rows (the
    report's spacing), so that the principal version, with four smaller heads, is shorter. B's axis label names Δ as
    the text does (on two lines: on one it would reach A's); the class headings set in the gaps between the classes
    (class_headings()); with all nine anatomies, the class headed APART_HEAD set off by a divider (0.45 gray) across
    the label column and both panels. The panel titles 0.04 in below the key's box (the report left about 0.16 in), so
    that Fig. 8 keeps its 7-pt type under its height cap."""
    a1, a2 = fig.axes[0], fig.axes[1]
    gs = a1.get_subplotspec().get_gridspec()
    if principal:
        square_keys(fig)
    a2.set_xlabel("Δ, paired change\nfrom the adult (dB)")
    floor_fonts(fig, 7.2)
    divider = class_headings(fig, a1, a2)
    W0, H0 = fig.get_size_inches()
    top_in, bottom_in = (1 - gs.top) * H0, gs.bottom * H0
    lo, hi = a1.get_ylim()
    for _ in range(2):  # the second time with the room above the panels cut to what the key and the titles need
        H = top_in + (hi - lo) * ROW_IN + bottom_in
        fig.set_size_inches(DRAWN_W, H)
        gs.update(top=1 - top_in / H, bottom=bottom_in / H)
        fit_width(fig, [gs], DRAWN_W)
        fig.canvas.draw()
        top_in -= inches(fig, fig.legends[0]).y0 - max(inches(fig, ax._left_title).y1 for ax in (a1, a2)) - 0.04
    if divider is not None:  # from the figure's left edge (fit_width() put what is drawn there) to B's right edge
        tr = matplotlib.transforms.blended_transform_factory(fig.transFigure, a1.transData)
        fig.add_artist(Line2D([0.0, a2.get_position().x1], [divider, divider], transform=tr, color="0.45", lw=0.8))
    fig.canvas.draw()
    marks = plotted_boxes(a2, skip=a2.texts)
    r, pad = fig.canvas.get_renderer(), fig.dpi / 72
    if any(overlaps(b, i, pad) for t in a2.texts for i in text_ink(t, r) for b in marks):
        raise SystemExit("Figure_R6_pediatric_D: a heading touches a marker or an interval")
    no_clashes(fig, "Figure_R6_pediatric_D")


def fix_maps_heads(fig: Figure, principal: bool = False) -> None:
    """R13 (bottom of Figs. 8 and S16): the whole image scaled to the text width (its type, 7.5 pt and up, kept), the
    legend of the greys moved left to stay inside, and the note moved below the color bar, into 0.075 in added at the
    bottom."""
    W0, H0 = fig.get_size_inches()
    dy = 0.075
    axes = [(ax, ax.get_position(original=True).bounds) for ax in fig.axes]
    texts = [(t, t.get_position()) for t in fig.texts]
    note = next(t for t in fig.texts if t.get_text().startswith("Fixed adult helmet"))
    handles, labels = list(fig.legends[0].legend_handles), [t.get_text() for t in fig.legends[0].get_texts()]
    fig.legends[0].remove()

    def scale(f):
        W, H = W0 * f, H0 * f + dy

        def up(y):  # a height of the report's layout (fraction of H0) -> fraction of the new figure, raised by dy
            return (y * H0 * f + dy) / H

        for ax, (x0, y0, w, h) in axes:
            ax.set_position([x0, up(y0), w, h * H0 * f / H])
        for t, (x, y) in texts:
            t.set_position((x, up(y)))
        fig.set_size_inches(W, H)
        note.set_position((0.04 / W0, 0.03 / H))
        for lg in fig.legends:
            lg.remove()
        fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(2.70 / W, up(0.24 / H0)), ncol=1, handlelength=1.2,
                   handleheight=0.9, fontsize=7.5, borderaxespad=0.0)
        fig.canvas.draw()
        return fig.get_tightbbox(fig.canvas.get_renderer()).width

    f = (DRAWN_W + 0.08) / W0  # the maps leave 0.04 in at either side; the color bar's end labels overhang a little
    over = scale(f) - (W0 - 0.08) * f
    width = scale((DRAWN_W - over) / (W0 - 0.08))
    if width > DRAWN_W + 0.005:
        raise SystemExit(f"Figure_R13_maps_heads: {width:.3f} in drawn, wider than the text width")
    space_map_titles(fig, "Figure_R13_maps_heads")


def insert_space(fig: Figure, cuts) -> None:
    """Vertical room added inside a figure laid out in figure fractions: at each (y, dy) of `cuts` (inches from the
    bottom), the axes, texts and legends lying above y move up by dy inches and the figure grows by the total; what lies
    below keeps its place."""
    fig.canvas.draw()
    W, H = fig.get_size_inches()
    H2 = H + sum(dy for _, dy in cuts)

    def lift(y):  # inches added below an element whose lower edge is y inches from the bottom
        return sum(dy for c, dy in cuts if y >= c)
    for ax in fig.axes:
        p = ax.get_position(original=True)
        ax.set_position([p.x0, (p.y0 * H + lift(p.y0 * H)) / H2, p.width, p.height * H / H2])
    for t in fig.texts:
        x, y = t.get_position()
        t.set_position((x, (y * H + lift(inches(fig, t).y0)) / H2))
    for lg in fig.legends:
        b = lg.get_bbox_to_anchor().transformed(fig.transFigure.inverted())
        lg.set_bbox_to_anchor((b.x0, (b.y0 * H + lift(b.y0 * H)) / H2), transform=fig.transFigure)
    fig.set_size_inches(W, H2)


def space_map_titles(fig: Figure, name: str, min_pt: float = 8.0) -> None:
    """Cortical-map figures (Figs. 6, 8C, S10, S16C: rows of four views, each under its title): room added above each
    row title, so that it stands at least min_pt clear of what lies above it (the views of the row above, by the extent
    of their cortex, or the column labels; the report's layout left 4-9 pt of white), everything above it moving up and
    the figure growing by as much (insert_space())."""
    import numpy as np
    fig.canvas.draw()
    r, k = fig.canvas.get_renderer(), fig.dpi / 72
    titles = [t for t in fig.texts if " vs Neuromag" in t.get_text()]
    if not titles:
        raise SystemExit(f"{name}: no row titles found")
    views = [matplotlib.transforms.Bbox.from_extents(*xy.min(axis=0), *xy.max(axis=0)) for ax in fig.axes
             for c in ax.collections if isinstance(c, matplotlib.collections.PolyCollection) and not ax.axison
             for xy in [ax.transData.transform(np.vstack([p.vertices for p in c.get_paths()]))]]
    cuts = []
    for t in titles:
        ink = text_ink(t, r)
        top, x0, x1 = max(b.y1 for b in ink), min(b.x0 for b in ink), max(b.x1 for b in ink)
        above = [b.y0 for b in views if b.y0 > top and b.x1 > x0 and b.x0 < x1]
        above += [b.y0 for s in fig.texts if s is not t for b in text_ink(s, r)
                  if b.y0 > top and b.x1 > x0 and b.x0 < x1]
        gap = (min(above) - top) / k if above else min_pt
        if gap < min_pt:
            cuts.append((top / fig.dpi + 0.001, (min_pt - gap) / 72))
    insert_space(fig, cuts)
    no_clashes(fig, name)


# Fig. S10 (R14): the scaled adults' row titles name them as the text does
SCALED_TITLES = {"School-age size (scaled adult): ": "School-age scaled adult: ",
                 "2-year size (scaled adult): ": "2-year scaled adult: "}


def fix_maps(name: str):
    """Figs. 6 (R11) and S10 (R14): their row titles set clear of the maps above them (space_map_titles()); in Fig. S10
    the scaled adults' titles name them as the text does ('School-age scaled adult', '2-year scaled adult')."""
    def fix(fig: Figure, principal: bool = False) -> None:
        if name == "Figure_R14_maps_scaled":
            for old, new in SCALED_TITLES.items():
                found = [t for t in fig.texts if old in t.get_text()]
                if len(found) != 1:
                    raise SystemExit(f"{name}: one row title with {old!r} expected, {len(found)} found")
                found[0].set_text(found[0].get_text().replace(old, new))
        space_map_titles(fig, name)
    return fix


def fix_confirm(fig: Figure, principal: bool = False) -> None:
    """R16 (Figs. 10 and S18): the Holm p columns left out (their values are in the text and Table S11; the header
    named the anatomies of the figure, not the nine of the adjustment), the location counts kept (7.2 pt, dark and light
    rows 8.8 pt apart), the panels widened into the space freed, the whole drawn at the text width with no text under
    7.2 pt; the principal version keys the runs with a square. With all nine anatomies (set_apart()), the class rule
    above APART is a divider (0.45 gray) across the label column and the panels of each row, labelled APART_HEAD below
    it in the label column (in A's and B's), as the class labels of Fig. S16."""
    ax_a, ax_ta, ax_b, ax_tb, ax_c, ax_tc = fig.axes[:6]
    if [a.axison for a in fig.axes[:6]] != [True, False, True, False, True, False]:
        raise SystemExit("Figure_R16_confirm: unexpected axes")
    for t in ax_ta.texts:
        if t.get_position()[0] > 0.6:  # the p column (x = 0.8) and its header
            t.set_visible(False)
        else:
            t.set_x(0.5)
    ax_tb.set_visible(False)
    ax_tc.set_visible(False)
    if principal:
        square_keys(fig)
    floor_fonts(fig, 7.2)
    H = fig.get_size_inches()[1]
    fig.set_size_inches(DRAWN_W, H)
    ya, yb = ax_a.get_position(), ax_b.get_position()
    pad = (ax_a.yaxis.get_tick_padding() + ax_a.yaxis.majorTicks[0].get_pad()) / 72
    lab = label_width(fig, ax_a.get_yticklabels() + ax_b.get_yticklabels()) + pad
    cnt = label_width(fig, [t for t in ax_ta.texts if t.get_visible()]) + 0.06
    L, R = lab + 0.01, DRAWN_W - 0.035
    place(ax_a, L, ya.y0 * H, R - cnt - 0.08 - L, ya.height * H)
    place(ax_ta, R - cnt, ya.y0 * H, cnt, ya.height * H)
    gap = 0.3
    wb = (R - L - gap) / 2
    place(ax_b, L, yb.y0 * H, wb, yb.height * H)
    place(ax_c, L + wb + gap, yb.y0 * H, wb, yb.height * H)
    apart = [y for y, t in zip(ax_a.get_yticks(), ax_a.get_yticklabels()) if t.get_text().endswith("(misregistered)")]
    for row in ((ax_a,), (ax_b, ax_c)) if apart else ():
        rules = [(ln, ln.get_ydata()[0]) for ax in row for ln in ax.get_lines()
                 if ln.get_transform() is not ax.transData and list(ln.get_xdata()) == [0, 1]
                 and ln.get_ydata()[0] == ln.get_ydata()[1]]
        y = min(v for _, v in rules if v > apart[0])  # the class rule above APART
        for ln, v in rules:
            ln.set_visible(v != y)
        tr = matplotlib.transforms.blended_transform_factory(fig.transFigure, row[0].transData)
        fig.add_artist(Line2D([0.01 / DRAWN_W, R / DRAWN_W], [y, y], transform=tr, color="0.45", lw=0.8))
        row[0].text(0.0, y, APART_HEAD, ha="right", va="top", fontsize=7.2, color="0.45", style="italic",
                    transform=offset_copy(row[0].get_yaxis_transform(), fig=fig, x=-pad * 72, y=-2.0, units="points"))
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    if bb.width > DRAWN_W + 0.01:
        raise SystemExit(f"Figure_R16_confirm: {bb.width:.3f} in drawn, wider than the text width")
    no_clashes(fig, "Figure_R16_confirm")


def legend_hits(ax, lg, full: bool = False) -> int:
    """Plotted points (line vertices, densified along drawn segments, and the vertices of collections) that fall inside
    the legend's box, its border pad excluded (included with full=True, for a legend whose filled box hides what lies
    under it): a legend set at a larger type size must not cover data. Guide lines (two vertices, no marker: reference
    levels, the identity line) are not data."""
    import numpy as np
    r = ax.get_figure(root=True).canvas.get_renderer()
    box = lg.get_window_extent(r)
    p = 0.0 if full else lg.borderpad * lg._fontsize * ax.get_figure(root=True).dpi / 72
    x0, x1, y0, y1 = box.x0 + p, box.x1 - p, box.y0 + p, box.y1 - p
    pts = []
    for ln in ax.get_lines():
        if not ln.get_visible() or (len(ln.get_xydata()) == 2 and ln.get_marker() in ("None", "none", "", None)):
            continue
        xy = np.asarray(ln.get_xydata(), float)
        if ln.get_linestyle() not in ("None", "none", "", " ") and len(xy) > 1:
            t = np.linspace(0, 1, 12)[:, None, None]
            xy = (xy[:-1] * (1 - t) + xy[1:] * t).reshape(-1, 2)
        pts.append(ln.get_transform().transform(xy[np.all(np.isfinite(xy), axis=1)]))
    for c in ax.collections:
        if c.get_visible():
            for path in c.get_paths():
                pts.append(c.get_transform().transform(path.vertices))
    pts = np.concatenate(pts) if pts else np.zeros((0, 2))
    inside = (pts[:, 0] > x0) & (pts[:, 0] < x1) & (pts[:, 1] > y0) & (pts[:, 1] < y1)
    return int(inside.sum())


def bars_under(ax, lg) -> int:
    """Bars (rectangles) of the axes that the legend's box, its border pad excluded, overlaps."""
    r = ax.get_figure(root=True).canvas.get_renderer()
    box = lg.get_window_extent(r)
    p = lg.borderpad * lg._fontsize * ax.get_figure(root=True).dpi / 72
    n = 0
    for patch in ax.patches:
        if isinstance(patch, matplotlib.patches.Rectangle) and patch.get_visible():
            b = patch.get_window_extent(r)
            if min(b.x1, box.x1 - p) > max(b.x0, box.x0 + p) and min(b.y1, box.y1 - p) > max(b.y0, box.y0 + p):
                n += 1
    return n


def tight_at_text_width(fig: Figure, H: float, keys, name: str, slack: float = 0.0, **kw) -> None:
    """tight_layout of a figure drawn at the text width (less `slack` inches: text measured on the screen canvas can come
    out a little wider in the PDF) and H inches high, the panels kept clear of the keys at the bottom: keys() places them
    for the figure's current size and returns their height (inches). Stops when what is drawn is wider than the text
    width (a key or a title too long) or when a panel's legend overruns the panel or covers plotted points."""
    pad = 0.3 * matplotlib.rcParams["font.size"] / 72  # tight_layout's pad, inches, cropped again when saved
    W = DRAWN_W - slack + 2 * pad
    inside = [a.get_legend() for a in fig.axes if a.get_legend() is not None]
    for lg in inside:  # legends inside the panels do not size them (wider than a panel at its old size, they would
        lg.set_in_layout(False)  # shrink it further at every pass); they are checked to lie inside below
    for _ in range(2):  # the second time narrower by what the first drew beyond the text width (labels at the edges)
        fig.set_size_inches(W, H)
        for _ in range(2):
            fig.canvas.draw()
            fig.tight_layout(rect=(0, (keys() + 0.08) / H, 1, 1), pad=0.3, **kw)
        fig.canvas.draw()
        W -= max(0.0, fig.get_tightbbox(fig.canvas.get_renderer()).width - (DRAWN_W - slack) + 0.002)
    for lg in inside:  # the legend's box, less its border pad, inside the panel and clear of the data
        lg.set_in_layout(True)
        a, b = inches(fig, lg), inches(fig, lg.axes)
        tol = lg.borderpad * lg._fontsize / 72 + 0.005
        if a.x0 < b.x0 - tol or a.x1 > b.x1 + tol or a.y0 < b.y0 - tol or a.y1 > b.y1 + tol:
            raise SystemExit(f"{name}: a legend overruns its panel ({a.x0:.2f}-{a.x1:.2f} in, panel {b.x0:.2f}-{b.x1:.2f})")
        hits = legend_hits(lg.axes, lg)
        if hits:
            raise SystemExit(f"{name}: the legend '{lg.get_texts()[0].get_text()[:30]}' covers {hits} plotted points")
    width = fig.get_tightbbox(fig.canvas.get_renderer()).width
    if width > DRAWN_W + 0.005:
        raise SystemExit(f"{name}: {width:.3f} in drawn, wider than the text width")


def fix_covariance_structure(fig: Figure, principal: bool = False) -> None:
    """Fig. S3: drawn at the text width and 7.0 in high (a page with its caption: height cap 0.72), type at least 7 pt
    (tick labels and keys 7, axis labels and panel titles 7.5). At this size the panels are 1.2 in wide: the keys of A-C
    stay inside their panels, compact and placed where they cover no data (C given more room above its rows); the keys
    that D and G carried for D-E and G-H, which no longer fit there, join that of F and I below the panels, each titled
    with its panels; axis labels and C's title longer than the panels are set on two lines. The keys of B and C, which
    the line at a ratio of 1 crosses, on white (no data under them)."""
    axs = fig.axes
    if len(axs) != 9:
        raise SystemExit(f"Figure_covariance_structure: {len(axs)} axes, 9 expected")
    wrap = {"magnetometers (fT)": "magnetometers\n(fT)", "gradiometers (fT/cm)": "gradiometers\n(fT/cm)",
            "left/right, front/back, upper/lower halves": "left/right, front/back,\nupper/lower halves",
            "random halves of the sites (95 % range)": "random halves of the\nsites (95 % range)",
            "exact model, as it would be measured": "exact model, as it\nwould be measured"}
    shared = {}
    for i, ax in enumerate(axs):
        lg = ax.get_legend()
        if lg is None:
            continue
        if i in (3, 6):  # D's key (lines of D and E) and G's (G and H): below the panels
            shared[i] = (list(lg.legend_handles), [t.get_text() for t in lg.get_texts()])
            lg.remove()
            continue
        compact = dict(handlelength=0.8, handletextpad=0.3, borderpad=0.2, borderaxespad=0.1, labelspacing=0.1)
        if i in (1, 2):  # B and C: the line at a ratio of 1 passes behind the key, not through its words
            compact.update(frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
        rebuild_legend(lg, labels=[wrap.get(t.get_text(), t.get_text()) for t in lg.get_texts()],
                       loc="upper left" if i == 0 else "best", fontsize=7.0, **compact)
    lo, hi = axs[2].get_ylim()
    axs[2].set_ylim(lo, hi + 0.55)  # (C): room for its key above the two rows
    old = fig.legends[0]  # F and I
    fi = (list(old.legend_handles), [t.get_text() for t in old.get_texts()], old.get_title().get_text())
    old.remove()
    titles = {"Spread over channels (5-95 %)": "Spread over\nchannels (5-95 %)"}
    labels = {"model / measured variance per channel": "model / measured variance\nper channel",
              "sensor distance (mm; 1st point: same site)": "sensor distance\n(mm; 1st point: same site)",
              "dimension of the dominant subspace": "dimension of the\ndominant subspace",
              "share of brain-noise variance": "share of brain-noise\nvariance", "agreement of correlations (r)": "agreement of\ncorrelations (r)"}
    for ax in axs:
        for lab in (ax.xaxis.label, ax.yaxis.label):  # longer than the panels at this width: on two lines
            lab.set_text(labels.get(lab.get_text(), lab.get_text()))
            lab.set_fontsize(7.5)
        title = ax._left_title
        letter, _, rest = title.get_text().partition("   ")
        title.set_text(f"{letter}   {titles.get(rest, rest)}")
        title.set_fontsize(7.5)
        ax.tick_params(axis="both", which="both", labelsize=7.0)
    floor_fonts(fig, 7.0)
    opts = dict(fontsize=7.0, title_fontsize=7.0, handlelength=1.6, handletextpad=0.4, columnspacing=1.0,
                labelspacing=0.3, borderpad=0.2, borderaxespad=0.0, frameon=False)
    made = {}

    def keys():  # D-E and G-H side by side above F-I, centered; returns their height (inches)
        W, H = fig.get_size_inches()
        for lg in made.values():
            lg.remove()
        made["fi"] = fig.legend(fi[0], fi[1], title=fi[2], ncol=2, loc="lower center", bbox_to_anchor=(0.5, 0.0), **opts)
        made["de"] = fig.legend(*shared[3], title="panels D and E", ncol=3, loc="lower left", **opts)
        made["gh"] = fig.legend(*shared[6], title="panels G and H", ncol=1, loc="lower left", **opts)
        fig.canvas.draw()
        h_fi = inches(fig, made["fi"]).height
        w_de, w_gh = inches(fig, made["de"]).width, inches(fig, made["gh"]).width
        x0 = (W - w_de - 0.35 - w_gh) / 2
        y0 = (h_fi + 0.06) / H
        made["de"].set_bbox_to_anchor((x0 / W, y0), transform=fig.transFigure)
        made["gh"].set_bbox_to_anchor(((x0 + w_de + 0.35) / W, y0), transform=fig.transFigure)
        fig.canvas.draw()
        return max(inches(fig, lg).y1 for lg in made.values())

    tight_at_text_width(fig, 7.0, keys, "Figure_covariance_structure", slack=0.06, h_pad=0.5, w_pad=0.5)
    for ax in axs[1:3]:
        if legend_hits(ax, ax.get_legend(), full=True):
            raise SystemExit("Figure_covariance_structure: the white key of B or C hides plotted points")


def fix_covariance_bands(fig: Figure, principal: bool = False) -> None:
    """Fig. S4: drawn at the text width and 3.6 in high, type at least 7 pt (tick labels and keys 7, axis labels and
    panel titles 7.5), the titles on two lines. The keys, which at this size would cover data in A and C, go below the
    panels side by side, each titled with its panels (A's lines are those of A and B)."""
    axs = fig.axes
    if len(axs) != 3 or any(ax.get_legend() is None for ax in (axs[0], axs[2])):
        raise SystemExit("Figure_covariance_bands: three panels with keys in A and C expected")
    shared = {}
    for i in (0, 2):
        lg = axs[i].get_legend()
        shared[i] = (list(lg.legend_handles), [t.get_text() for t in lg.get_texts()])
        lg.remove()
    floor_fonts(fig, 7.0)
    for ax in axs:
        title = ax._left_title
        title.set_text(title.get_text().replace(" (brain:", "\n(brain:").replace("Magnetometer brain",
                                                                                   "Magnetometer\nbrain"))
        title.set_fontsize(7.5)
        for lab in (ax.xaxis.label, ax.yaxis.label):
            lab.set_fontsize(7.5)
        ax.tick_params(axis="both", which="both", labelsize=7.0)
    opts = dict(fontsize=7.0, title_fontsize=7.0, handlelength=2.0, handletextpad=0.5, columnspacing=1.2,
                labelspacing=0.3, borderpad=0.2, borderaxespad=0.0, frameon=False, alignment="left")
    made = {}

    def keys():  # A and B's key (two columns) and C's beside it, tops level, centered as a pair; returns their height
        W, H = fig.get_size_inches()
        for lg in made.values():
            lg.remove()
        made["ab"] = fig.legend(*shared[0], title="panels A and B", ncol=2, loc="upper left", **opts)
        made["c"] = fig.legend(*shared[2], title="panel C", ncol=1, loc="upper left", **opts)
        fig.canvas.draw()
        (w_ab, h_ab), (w_c, h_c) = [(inches(fig, lg).width, inches(fig, lg).height) for lg in made.values()]
        x0, top = (W - w_ab - 0.4 - w_c) / 2, max(h_ab, h_c)
        made["ab"].set_bbox_to_anchor((x0 / W, top / H), transform=fig.transFigure)
        made["c"].set_bbox_to_anchor(((x0 + w_ab + 0.4) / W, top / H), transform=fig.transFigure)
        fig.canvas.draw()
        return max(inches(fig, lg).y1 for lg in made.values())

    tight_at_text_width(fig, 3.6, keys, "Figure_covariance_bands", slack=0.02, w_pad=1.0)
    no_clashes(fig, "Figure_covariance_bands")


def decorations(fig: Figure, ax) -> tuple[float, float, float]:
    """Extent (inches) of an axes' tick labels and axis labels to its left and below it, and of its title above it."""
    box = inches(fig, ax)
    shown = lambda ts: [inches(fig, t) for t in ts if t.get_visible() and (t.get_text() or "").strip()]  # noqa: E731
    left = shown(ax.get_yticklabels() + [ax.yaxis.label])
    below = shown(ax.get_xticklabels() + [ax.xaxis.label])
    above = shown([ax._left_title])
    return (box.x0 - min((b.x0 for b in left), default=box.x0), box.y0 - min((b.y0 for b in below), default=box.y0),
            max((b.y1 for b in above), default=box.y1) - box.y1)


def data_y(ax) -> tuple[float, float]:
    """Lowest and highest y of what the axes plots (lines, their markers' centers, and filled bands; guide lines with
    two points and no marker left out)."""
    import numpy as np
    ys = []
    for ln in ax.get_lines():
        if ln.get_visible() and not (len(ln.get_xydata()) == 2 and ln.get_marker() in ("None", "none", "", None)):
            ys.append(np.asarray(ln.get_ydata(), float))
    for c in ax.collections:
        for path in c.get_paths():
            ys.append(path.vertices[:, 1])
    ys = np.concatenate(ys)
    ys = ys[np.isfinite(ys)]
    return float(ys.min()), float(ys.max())


def fix_noise_sensitivity(fig: Figure, principal: bool = False) -> None:
    """Fig. S5 (A-F), recomposed at the text width and 6.9 in high (height cap 0.72: its long caption shares the page;
    the report drew it 7.2 x 10.2 in): type at least 6.5 pt (tick labels, keys and the rows of B and F 6.5, axis labels
    7, panel titles 7.2 as in G-J), every panel and key kept. Laid out by hand: A beside B and E beside F, the plots of
    B and F aligned at the right, their rows about 11 pt apart (10 pt in the report's print); C and D side by side, of
    equal width; each title from its column's left edge. Keys: A's inside A, its long entries on two lines and A's axis
    extended downwards to hold it; C and D's below them in one row, and the arrays' and markers' at the bottom in one
    row, entries on two lines (the first of C and D's on three); E's inside E in two columns. The bands of A, C and D are
    named parcel-bootstrap intervals (the report: '95 % CI'; EXACT). Axis labels longer than their axes go on two lines
    (the same words); B and F are labelled every 0.1 with unlabelled ticks between, their bold section heads set from the
    left edge of the row-label column; E's axis ends just beyond its bars."""
    import numpy as np
    panels = fig.axes[:6]
    A, B, C, D, E, F = panels
    if len(fig.axes) != 6 or any(ax.get_legend() is None for ax in (A, C, E)) or len(fig.legends) != 1:
        raise SystemExit("Figure_noise_sensitivity: six panels with keys in A, C and E and one below them expected")
    FS = 6.5
    fig.set_layout_engine("none")
    W, H = DRAWN_W, 6.9
    fig.set_size_inches(W, H)
    floor_fonts(fig, FS)
    two = {s.replace("\n", " "): s for s in (  # axis labels and key entries on two lines (the same words)
        "detectability ratio\nOPM / Neuromag (306)", "detectability ratio OPM /\nNeuromag (306 channels)",
        "field energy left\nby the projection", "near-skull cortex + 10-Hz 1/f corner\n(frequency-resolved)",
        "the same + heart and eyes\n(magnetometer shortfall)",
        "sensor + brain noise (band:\n95 % parcel-bootstrap interval)",
        "frequency-resolved, flat signal\nspectrum (band: 95 %\nparcel-bootstrap interval)",
        "frequency-resolved,\nspike-wave spectrum", "band variance (spatial\nwhitening only)",
        "frequency-resolved, Neuromag\nnoise from the empty room", "with room field, after\n8-term projection")}
    for ax in panels:
        ax._left_title.set_fontsize(7.2)  # as the titles of G-J
        for lab in (ax.xaxis.label, ax.yaxis.label):
            lab.set_fontsize(7.0)
        rewrap([ax.xaxis.label, ax.yaxis.label], two)
        ax.tick_params(axis="both", which="both", labelsize=FS)
    heads = {}
    for ax in (B, F):
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.05))
        ax.tick_params(axis="x", which="minor", length=2.0, width=0.6)
        # the bold section heads (rows without data) from the left edge of the row-label column, so that they do not
        # set its width
        rows = [(t.get_position()[1], t.get_text(), t.get_fontweight() == "bold") for t in ax.get_yticklabels()]
        if [y for y, _, _ in rows] != [float(y) for y in ax.get_yticks()] or sum(b for _, _, b in rows) != 2:
            raise SystemExit("Figure_noise_sensitivity: B and F should each have two bold section heads")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["" if b else s for _, s, b in rows]))
        tr = matplotlib.transforms.blended_transform_factory(fig.transFigure, ax.transData)
        heads[ax] = [ax.text(0.0, y, s, transform=tr, ha="left", va="center", fontsize=FS, fontweight="bold")
                     for y, s, b in rows if b]
    compact = dict(fontsize=FS, borderpad=0.2, labelspacing=0.35, handletextpad=0.5, frameon=False)

    def keep(lg):
        return list(lg.legend_handles), [rewrapped(t.get_text(), two) for t in lg.get_texts()]

    key_a = rebuild_legend(A.get_legend(), *keep(A.get_legend()), loc="lower left", handlelength=2.6,
                           borderaxespad=0.3, **compact)
    key_e = rebuild_legend(E.get_legend(), loc="upper center", ncol=2, handlelength=1.0, columnspacing=0.7,
                           borderaxespad=0.3, **compact)
    hc, lc = keep(C.get_legend())
    C.get_legend().remove()
    key_cd = fig.legend(hc, lc, loc="upper center", ncol=4, handlelength=2.6, columnspacing=1.0, borderaxespad=0.0,
                        **compact)
    hf, lf = keep(fig.legends[0])
    fig.legends[0].remove()
    key_all = fig.legend(hf, lf, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, handlelength=1.6,
                         columnspacing=1.4, borderaxespad=0.0, **compact)
    r, gap, vgap, w_bf, h_cd = 0.08, 0.18, 0.1, 1.45, 1.0  # right margin (the last tick label overhangs), column gap,
    # row gap, plot width of B and F, plot height of C and D
    bars = [(p.get_x(), p.get_x() + p.get_width()) for p in E.patches]
    E.set_xlim(min(a for a, _ in bars) - 0.05, max(b for _, b in bars) + 0.05)  # E's margins: its labels need the width
    units_b, units_f = np.ptp(B.get_ylim()), np.ptp(F.get_ylim())  # rows of B and F (plus their margins)
    for _ in range(4):
        fig.canvas.draw()
        dec = {ax: decorations(fig, ax) for ax in panels}
        for ax in (B, F):  # the row-label column holds the section heads too
            col = max(dec[ax][0], max(inches(fig, t).width for t in heads[ax]) + 0.06)
            dec[ax] = (col,) + dec[ax][1:]
        k_cd, k_all = inches(fig, key_cd).height, inches(fig, key_all).height
        lw = max(dec[ax][0] for ax in (A, C, E))
        x_bf = W - r - w_bf
        w_cd = (W - r - lw - dec[D][0] - gap) / 2
        t1, d1 = max(dec[A][2], dec[B][2]), max(dec[A][1], dec[B][1])
        t2, d2 = max(dec[C][2], dec[D][2]), max(dec[C][1], dec[D][1])
        t3, d3 = max(dec[E][2], dec[F][2]), max(dec[E][1], dec[F][1])
        fixed = t1 + d1 + vgap + t2 + h_cd + d2 + 0.06 + k_cd + vgap + t3 + d3 + vgap + k_all
        u = (H - fixed) / (units_b + units_f)
        h1, h3 = units_b * u, units_f * u
        y = H - t1 - h1
        place(A, lw, y, x_bf - dec[B][0] - gap - lw, h1)
        place(B, x_bf, y, w_bf, h1)
        y -= d1 + vgap + t2 + h_cd
        place(C, lw, y, w_cd, h_cd)
        place(D, W - r - w_cd, y, w_cd, h_cd)
        y -= d2 + 0.06
        key_cd.set_bbox_to_anchor((0.5, y / H), transform=fig.transFigure)
        y -= k_cd + vgap + t3 + h3
        place(E, lw, y, x_bf - dec[F][0] - gap - lw, h3)
        place(F, x_bf, y, w_bf, h3)
        for ax in (B, F):
            for t in heads[ax]:
                t.set_x((x_bf - dec[ax][0]) / W)
    fig.canvas.draw()
    for ax in (B, F):  # the section heads end before the plot
        if max(t.get_window_extent().x1 for t in heads[ax]) > ax.get_window_extent().x0 - 2:
            raise SystemExit("Figure_noise_sensitivity: a section head reaches into its plot")
    for ax in panels:  # each title from its column's left edge
        left = 0.0 if ax in (A, C, E) else inches(fig, ax).x0 - dec[ax][0]
        ax._left_title.set_x((left - inches(fig, ax).x0) / inches(fig, ax).width)
    # A's axis extended downwards so that its key lies below the curves; E's upwards so that its key lies above the bars
    fig.canvas.draw()
    lo, hi = A.get_ylim()
    frac = (inches(fig, key_a).height + 0.08) / inches(fig, A).height
    A.set_ylim((data_y(A)[0] - frac * hi) / (1 - frac), hi)
    top_bar = max(p.get_y() + p.get_height() for p in E.patches)
    lo, hi = np.log10(E.get_ylim())
    frac = (inches(fig, key_e).height + 0.04) / inches(fig, E).height
    if (hi - np.log10(top_bar)) / (hi - lo) < frac:
        E.set_ylim(10 ** lo, 10 ** ((np.log10(top_bar) - frac * lo) / (1 - frac)))
    fig.canvas.draw()
    if legend_hits(A, key_a) or bars_under(E, key_e):
        raise SystemExit("Figure_noise_sensitivity: a key covers data")
    width = fig.get_tightbbox(fig.canvas.get_renderer()).width
    if width > DRAWN_W + 0.005:
        raise SystemExit(f"Figure_noise_sensitivity: {width:.3f} in drawn, wider than the text width")
    no_clashes(fig, "Figure_noise_sensitivity")


def fix_noise_depth(fig: Figure, principal: bool = False) -> None:
    """Fig. S5 (G-J), printed upright on a portrait page instead of rotated: the 2 x 2 panels drawn at the text width and
    7.0 in high, type at least 7 pt, each title on two lines (array; condition), the key below the panels in two
    columns."""
    floor_fonts(fig, 7.0)
    for ax in fig.axes:
        t = ax._left_title
        t.set_text(re.sub(r"(\(\d+ sites\)), ", r"\1,\n", t.get_text()))
    leg = rebuild_legend(fig.legends[0], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=7.0,
                         handlelength=2.2, columnspacing=1.2)
    tight_at_text_width(fig, 7.0, lambda: inches(fig, leg).y1, "Figure_noise_sensitivity_depth", h_pad=1.0, w_pad=1.2)


def thin_tick_labels(fig: Figure, gap_pt: float = 2.0) -> None:
    """On every x axis whose labels (fixed labels) come closer than gap_pt, every other label is left blank, keeping 0's."""
    fig.canvas.draw()
    for ax in fig.axes:
        fmt = ax.xaxis.get_major_formatter()
        fixed = isinstance(fmt, matplotlib.ticker.FixedFormatter)
        bydict = (isinstance(fmt, matplotlib.ticker.FuncFormatter) and isinstance(fmt.func, functools.partial)
                  and fmt.func.args and isinstance(fmt.func.args[0], dict))  # set_ticklabels on fixed ticks
        if not (fixed or bydict):
            continue
        shown = [t for t in ax.get_xticklabels() if t.get_visible() and t.get_text()]
        boxes = sorted(inches(fig, t).x0 * 72 for t in shown), sorted(inches(fig, t).x1 * 72 for t in shown)
        if not any(x0 - x1 < gap_pt for x0, x1 in zip(boxes[0][1:], boxes[1][:-1])):
            continue
        ticks = [float(t) for t in ax.get_xticks()]
        keep = ticks.index(0.0) % 2 if 0.0 in ticks else 0
        if fixed:
            fmt.seq = [s if i % 2 == keep else "" for i, s in enumerate(fmt.seq)]
        else:
            labels = fmt.func.args[0]
            for i, t in enumerate(ticks):
                if i % 2 != keep and t in labels:
                    labels[t] = ""


def fix_supplement_spikes(fig: Figure, principal: bool = False) -> None:
    """Figs. S14 and S15, printed at the text width (0.88-0.90 of their drawn size): type at least 8 pt; where tick
    labels then crowd, every other one is left blank. The p marks (* and **), which the report set for 6.6- and 7-pt
    type, placed again for 8-pt type (put_marks()): S14's beside its interval, left or, where that place is taken, right
    of it (as the report did), its panels' x ranges extended a little to the left where a mark left of its interval
    would reach the axis (panels that share a range keep sharing it); S15's right of its bar. What is drawn keeps the
    report's figure width (the panels narrowed from the left where the row labels need more room: with the 18-month
    template set apart, its label, '(misregistered)', is the longest)."""
    floor_fonts(fig, 8.0)
    fit_width(fig, [fig.axes[0].get_subplotspec().get_gridspec()], fig.get_size_inches()[0] - 0.06)
    thin_tick_labels(fig)
    joint = all(ax.patches for ax in fig.axes)  # S15: bars; S14: intervals
    found = {}
    for ax in fig.axes:
        items = []
        for t in [t for t in ax.texts if t.get_text() in ("*", "**")]:
            if joint:  # S15: written at the row's height, 1.5 % right of the bar of detected events
                y = t.get_position()[1]
                ends = [p.get_x() + p.get_width() for p in ax.patches if abs(p.get_y() + p.get_height() / 2 - y) < 1e-9]
                items.append((t, y, [(max(ends), 1)]))
            else:  # S14: written 0.07 below the row, beside its interval
                y = t.get_position()[1] + 0.07
                spans = [sorted(ln.get_xdata()) for ln in ax.get_lines() if ln.get_transform() is ax.transData
                         and len(ln.get_xdata()) == 2 and ln.get_marker() in ("None", "none", "", None)
                         and abs(ln.get_ydata()[0] - y) < 1e-9 and abs(ln.get_ydata()[1] - y) < 1e-9]
                if len(spans) != 1:
                    raise SystemExit(f"{fig.axes[0].get_title()}: no single interval at the row of a mark")
                items.append((t, y, [(spans[0][0], -1), (spans[0][1], 1)]))
        if items:
            found[ax] = items
    if not joint:  # S14: room on the left for a mark beside each interval
        fig.canvas.draw()
        ranges = {}
        for ax in fig.axes:
            ranges.setdefault(tuple(ax.get_xlim()), []).append(ax)
        for (x0, x1), group in ranges.items():
            need = 0.0
            for ax in group:
                per_pt = (x1 - x0) / (inches(fig, ax).width * 72)
                for t, _, cands in found.get(ax, []):
                    need = max(need, x0 - (cands[0][0] - (1.3 + mark_glyphs(t).width + 1.5) * per_pt))
            for ax in group if need > 0 else []:
                ax.set_xlim(x0 - 1.1 * need, x1)
    for ax, items in found.items():
        right = 0.0
        if joint:  # S15: up to 2 pt before the next panel (its bars start at its edge) or the figure's edge
            fig.canvas.draw()
            x1 = ax.get_window_extent().x1
            nxt = [a.get_window_extent().x0 for a in fig.axes if a.get_window_extent().x0 > x1]
            right = ((min(nxt) if nxt else fig.bbox.x1) - x1) * 72 / fig.dpi - 2.0
        put_marks(ax, items, right_over=right)
    no_clashes(fig, "Figure_S_" + ("joint_detection_localization" if joint else "localization_effects"))


def mark_glyphs(t: Text):
    """A mark's glyphs (TextPath extents, points from its anchor set at ha='left', va='baseline')."""
    from matplotlib.textpath import TextPath
    return TextPath((0, 0), t.get_text(), prop=t.get_fontproperties()).get_extents()


def put_marks(ax, items, gap: float = 1.3, pad_drawn: float = 0.4, pad_marks: float = 1.0,
              right_over: float = 0.0) -> None:
    """Significance marks placed on their rows. items: (text, row's y, candidates), a candidate being (x, side): the
    mark's glyphs `gap` pt left (side -1) or right (+1) of x (data units), then farther out in steps of their width plus
    1.5 pt. Each mark, from the top row down, takes the first place, nearest first and candidates in turn, where its
    glyphs, centered on the row, lie inside the axes (on the right up to right_over pt beyond it) and clear what the
    axes draws (markers, lines with their width, bars) by pad_drawn pt and the marks already placed by pad_marks pt;
    where a place next to its interval is taken only by one mark placed before it, that mark moves one step farther out
    if it can, so that both stay near their intervals. Guide lines across the axes (the zero line, the class rules) are
    drawn over by the marks, which are set above them on a thin white outline and below the data."""
    from matplotlib.patheffects import withStroke
    fig = ax.get_figure(root=True)
    k = fig.dpi / 72
    fig.canvas.draw()
    drawn = plotted_boxes(ax, skip=[t for t, _, _ in items])
    frame = ax.get_window_extent(fig.canvas.get_renderer())
    placed = {}  # text -> (glyph box, x, y, side, j)

    def spot(t, x, y, side, j):  # (dx, dy) of the anchor in points and the glyph box (display units)
        e = mark_glyphs(t)
        off = gap + j * (e.width + 1.5)
        dx, dy = (off - e.x0 if side > 0 else -off - e.x1), -(e.y0 + e.y1) / 2
        X, Y = ax.transData.transform((x, y))
        return dx, dy, matplotlib.transforms.Bbox.from_extents(X + (dx + e.x0) * k, Y + (dy + e.y0) * k,
                                                               X + (dx + e.x1) * k, Y + (dy + e.y1) * k)

    def clear(ink):  # inside the axes and clear of the data
        return (frame.x0 + k <= ink.x0 and ink.x1 <= frame.x1 + right_over * k
                and not any(overlaps(ink, b, pad_drawn * k) for b in drawn))

    def put(t, x, y, side, j):
        dx, dy, ink = spot(t, x, y, side, j)
        t.set_position((x, y))
        t.set_transform(offset_copy(ax.transData, fig=fig, x=dx, y=dy, units="points"))
        t.set_horizontalalignment("left")
        t.set_verticalalignment("baseline")
        t.set_zorder(1.5)  # guide lines 0-1, data 2 and up
        t.set_path_effects([withStroke(linewidth=1.6, foreground="white")])
        placed[t] = (ink, x, y, side, j)

    for t, y, cands in sorted(items, key=lambda it: -it[1]):  # top row first
        for j, (x, side) in ((j, c) for j in range(8) for c in cands):
            ink = spot(t, x, y, side, j)[2]
            if not clear(ink):
                continue
            taken = [s for s, p in placed.items() if overlaps(ink, p[0], pad_marks * k)]
            if not taken:
                put(t, x, y, side, j)
                break
            if j == 0 and len(taken) == 1:  # the one mark in the way one step farther out, if it can
                s = taken[0]
                _, xs, ys, sides, js = placed[s]
                moved = spot(s, xs, ys, sides, js + 1)[2]
                others = [p[0] for u, p in placed.items() if u is not s] + [ink]
                if clear(moved) and not any(overlaps(moved, b, pad_marks * k) for b in others):
                    put(s, xs, ys, sides, js + 1)
                    put(t, x, y, side, j)
                    break
        else:
            raise SystemExit(f"put_marks: no place for the mark at y = {y:g} in '{ax.get_title()}'")


def fix_helmet_fit(fig: Figure, principal: bool = False) -> None:
    """Figs. 9 and S17 (helmet_fit_figure, drawn at the text width): what is drawn spans the text width exactly."""
    fit_width(fig, [fig.axes[0].get_subplotspec().get_gridspec()], DRAWN_W)


def no_clashes(fig: Figure, name: str, **kw) -> None:
    """Stops when texts of the figure overlap or touch (clashes())."""
    found = clashes(fig, **kw)
    if found:
        raise SystemExit(f"{name}: texts overlap or touch: " + "; ".join(found))


def ratio_ticks(ax, labeled, unlabeled) -> None:
    """A logarithmic ratio axis (x) labelled at the given ratios only, with unlabelled ticks at the others."""
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(labeled))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter([f"{t:g}" for t in labeled]))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(unlabeled))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="x", which="minor", length=2.0, width=0.6)


def rewrapped(s: str, breaks: dict) -> str:
    """The text with its line breaks reset ({old text: new text}); stops if the new text changes the words."""
    new = breaks.get(s, s)
    if new.replace("\n", " ").replace("- ", "-") != s.replace("\n", " ").replace("- ", "-"):
        raise SystemExit(f"rewrap: {s!r} -> {new!r} changes the words")
    return new


def rewrap(texts, breaks: dict) -> int:
    """Line breaks of the given texts reset ({old text: new text}, the same words); returns how many were reset."""
    n = 0
    for t in texts:
        new = rewrapped(t.get_text(), breaks)
        if new != t.get_text():
            t.set_text(new)
            n += 1
    return n


def fix_noise_model(fig: Figure, principal: bool = False) -> None:
    """Fig. S2 (R2): drawn at the text width and 5.75 in high (height cap 0.6), no text under 7.0 pt, panel titles 8 pt
    and axis labels 7.5 pt; the constrained layout re-flows the panels (the empty row of the left-out note dropped).
    Texts that the narrower panels would make collide are set on more lines (the same words): the status notes of B
    (B's key moved below them) and the labels under B and C. Under C's bars the word 'Neuromag', shared by its two
    channel types, is set once below 'mag.' and 'grad.', on the line where 'OPM' stands under the OPM bars."""
    panels = [ax for ax in fig.axes if ax.get_visible()]
    if len(panels) != 5:
        raise SystemExit(f"Figure_R2_noise_model: {len(panels)} panels, 5 expected")
    a, g, b1, b2, c = panels
    floor_fonts(fig, 7.0)
    breaks = {"room field:\n94 % of variance": "room field:\n94 % of\nvariance",
              "datasheet sensor noise;\nroom field 1.6 %": "datasheet\nsensor noise;\nroom field 1.6 %",
              "sensor-noise variance over the summed median\nvariances of the three components, as in (A)":
                  "sensor-noise variance over the\nsummed median variances of the\nthree components, as in (A)"}
    for ax in panels:
        ax._left_title.set_fontsize(8.0)
        ax.yaxis.label.set_fontsize(7.5)
        lab = ax.xaxis.label  # under B: '<channels>; per channel, model /\nmeasured brain variance\n<range>, median <m>'
        lab.set_text(lab.get_text().replace("; per channel, model /\nmeasured brain variance\n",
                                            "; per channel,\nmodel / measured brain variance\n"))
        rewrap([lab] + ax.texts, breaks)
    old = b1.get_legend()  # above the empty-room bars, clear of the brain-noise bars at this width
    rebuild_legend(old, loc="upper left", bbox_to_anchor=(0.0, 0.62), fontsize=old.get_texts()[0].get_fontsize(),
                   handlelength=1.4, handletextpad=0.5, borderpad=0.2, borderaxespad=0.2)
    labels = [t.get_text() for t in c.get_xticklabels()]
    if labels[2:] != ["Neuromag\nmag.", "Neuromag\ngrad."]:
        raise SystemExit(f"Figure_R2_noise_model: C's bar labels {labels}")
    c.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(
        [labels[0], labels[1].replace("site-matched", "site-\nmatched"), "mag.", "grad."]))
    fig.canvas.draw()
    first = c.get_xticklabels()[0]  # 'dense\nOPM': 'Neuromag' level with its second line
    dy = (first.get_window_extent().y0 - c.get_window_extent().y0) * 72 / fig.dpi
    c.annotate("Neuromag", xy=(2.5, 0), xycoords=c.get_xaxis_transform(), xytext=(0, dy), textcoords="offset points",
               ha="center", va="bottom", fontsize=first.get_fontsize())
    outer = a.get_subplotspec().get_topmost_subplotspec().get_gridspec()
    outer.set_height_ratios([1.0, 1.15, 0.001])
    fit_constrained(fig, DRAWN_W - 0.02, 5.75)
    no_clashes(fig, "Figure_R2_noise_model")
    for ax in panels:
        if ax.get_legend() is not None and bars_under(ax, ax.get_legend()):
            raise SystemExit(f"Figure_R2_noise_model: a key covers {bars_under(ax, ax.get_legend())} bars")


def fix_children_qc(fig: Figure, principal: bool = False) -> None:
    """Fig. S6: laid out by its script in inches at the report's 7.2-in width (five MRI sections side by side), printed
    at about 0.88 of that: its small type raised to 8 pt (7 pt printed). A's section titles raised to stand 3.5 pt
    clear of their sections (at the report's pad their descenders came within 1.8 pt). B is narrowed by 0.18 in, so
    that its axis label stays inside the figure and the row labels of C (moved 0.05 in to the right) clear it."""
    floor_fonts(fig, 8.0)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ax in [ax for ax in fig.axes if ax.images]:
        gap = (ax.get_window_extent(r).y1 - min(b.y0 for b in text_ink(ax.title, r))) * 72 / fig.dpi
        ax._set_title_offset_trans(ax.titleOffsetTrans.get_matrix()[1, 2] * 72 / fig.dpi + 3.5 - gap)
    b = next(ax for ax in fig.axes if ax.get_xlabel().startswith("Distance"))
    c = next(ax for ax in fig.axes if ax.get_xlabel().startswith("MRI head boundary"))
    pb, pc = inches(fig, b), inches(fig, c)
    place(b, pb.x0 + 0.08, pb.y0, pb.width - 0.18, pb.height)
    place(c, pc.x0 + 0.05, pc.y0, pc.width - 0.05, pc.height)
    old = b.get_legend()  # its white box (over the dotted 8- and 10-mm lines) made compact: it must not hide the curves
    key = rebuild_legend(old, loc="upper left", fontsize=8.0, handlelength=1.5, handletextpad=0.5, borderaxespad=0.1,
                         borderpad=0.2, frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    fig.canvas.draw()
    rows_c = min(inches(fig, t).x0 for t in c.get_yticklabels())
    if inches(fig, b.yaxis.label).x0 < 0.0 or rows_c < inches(fig, b).x1 + 0.05:
        raise SystemExit("Figure_S_children_qc: B's axis label leaves the figure or C's row labels reach B")
    if legend_hits(b, key, full=True):
        raise SystemExit(f"Figure_S_children_qc: B's key covers {legend_hits(b, key, full=True)} plotted points")
    no_clashes(fig, "Figure_S_children_qc")


def fix_helmet_routes(fig: Figure, principal: bool = False) -> None:
    """Fig. S8 (R7): drawn by its script at the report's 7.2-in width, printed at about 0.89 of that: its small type
    raised to 7.4 pt (6.5 pt printed; at 8 pt the key of B would overrun the figure). A's third condition label on four
    lines (the same words), clear of its neighbor. A's note on the gray band names the adult's placements as everywhere
    else, its 12 source-blind placements (the report: 'standard'), on three lines, below the band and clear of the
    curves and of the 'adult' label. C's title says what 'ahead' meant ('where the OPM's detectability is higher'), on two
    lines."""
    floor_fonts(fig, 7.4)
    c = [ax for ax in fig.axes if ax.get_title(loc="left").endswith("Share of cortex where the OPM is ahead")]
    if len(c) != 1:
        raise SystemExit("Figure_R7_helmet_fit: C's title not found")
    c[0]._left_title.set_text(c[0]._left_title.get_text().replace("where the OPM is ahead",
                                                                  "where the OPM's\ndetectability is higher"))
    a = next(ax for ax in fig.axes if ax.get_title(loc="left").endswith("OPM advantage at five helmet conditions"))
    labels = [t.get_text() for t in a.get_xticklabels()]
    a.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(
        [rewrapped(s, {"laterally\ncentered,\ntop contact": "laterally\ncentered,\ntop\ncontact"}) for s in labels]))
    note = [t for t in a.texts if t.get_text().startswith("gray band: adult over its ")]
    m = note and re.fullmatch(r"gray band: adult over its (\d+) standard\nplacements (\(.+\))", note[0].get_text())
    if not m:
        raise SystemExit("Figure_R7_helmet_fit: A's note on the gray band not found")
    note[0].set_text(f"gray band:\nadult over its {m.group(1)} source-blind\nplacements {m.group(2)}")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    if any(overlaps(b, i, fig.dpi / 72) for i in text_ink(note[0], r) for b in plotted_boxes(a, skip=note)):
        raise SystemExit("Figure_R7_helmet_fit: A's note on the gray band reaches what A draws")
    no_clashes(fig, "Figure_R7_helmet_fit")


def fix_regions_heads(fig: Figure, principal: bool = False) -> None:
    """Fig. S9 (R5): drawn by its script at the report's 7.2-in width, printed at about 0.876 of that: its small type
    raised to 7.45 pt (6.5 pt printed; at 8 pt the row labels and the column of target counts would overrun the
    figure). The column headers of the two scaled adults, too wide for their columns, on two lines like the others
    ('School-age' and '2-year' over the head circumference: 'size' is said by their class label, 'scaled adult (size
    only)'), each moved 2 pt off its column's center, away from the other, so that they stand as far apart as the other
    headers; the heat maps narrowed from the left so that the longest row label stays inside the figure, and the note
    across B's hatched row on two lines; the class brackets raised to clear the headers, and each title kept as far
    above the brackets as the report set it."""
    import numpy as np
    from matplotlib.transforms import ScaledTranslation
    heat = [ax for ax in fig.axes if ax.get_title(loc="left")]
    if len(heat) != 2:
        raise SystemExit(f"Figure_R5_regions_heads: {len(heat)} titled panels, 2 expected")
    classes = set(americanize_text(s) for s in style.CLASS_LABEL.values())

    def parts(ax):
        brackets = [ln for ln in ax.get_lines() if not ln.get_clip_on() and len(ln.get_xdata()) == 2
                    and ln.get_marker() in ("None", "none", "", None)]
        names = [t for t in ax.texts if t.get_text() in classes]
        return brackets, names

    fig.canvas.draw()
    gap0 = {ax: inches(fig, ax._left_title).y0 - max(inches(fig, t).y1 for t in parts(ax)[1]) for ax in heat}
    floor_fonts(fig, 7.45)
    short = {f"{s}\n": f"{s.removesuffix(' size')}\n" for s in HEAD_SHORT.values()}  # their class says 'size'
    for ax in heat:
        labels = [t.get_text() for t in ax.get_xticklabels()]
        new = [functools.reduce(lambda x, kv: x.replace(*kv), short.items(), s) for s in labels]
        if sum(n != s for s, n in zip(labels, new)) != len(short):
            raise SystemExit(f"Figure_R5_regions_heads: the scaled adults' column headers not found in {labels}")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(new))
    fig.canvas.draw()
    # the heat maps narrowed from the left so that the longest row label stays inside the figure
    shift = 0.02 - min(inches(fig, t).x0 for ax in heat for t in ax.get_yticklabels() if t.get_text())
    if shift > 0:
        for ax in heat:
            b = inches(fig, ax)
            place(ax, b.x0 + shift, b.y0, b.width - shift, b.height)
    fig.canvas.draw()
    for ax in heat:  # the two scaled adults' headers 2 pt apart from their columns' centers (tick marks have no length)
        ticks = [float(x) for x in ax.get_xticks()]
        i = next(i for i, t in enumerate(ax.get_xticklabels()) if t.get_text().startswith("School-age\n"))
        if not ax.get_xticklabels()[i + 1].get_text().startswith("2-year\n"):
            raise SystemExit("Figure_R5_regions_heads: the 2-year size's header does not follow the school-age size's")
        d = 2.0 / 72 / inches(fig, ax).width * abs(np.diff(ax.get_xlim())[0])
        ticks[i], ticks[i + 1] = ticks[i] - d, ticks[i + 1] + d
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
    fig.canvas.draw()
    for ax in heat:  # the note across the hatched row (B) on two lines when wider than the row
        for t in ax.texts:
            s = t.get_text()
            if s.endswith(": not computed for this helmet") and inches(fig, t).width > inches(fig, ax).width - 0.2:
                t.set_text(rewrapped(s, {s: s.replace(": not computed", ":\nnot computed")}))
    fig.canvas.draw()
    for ax in heat:
        brackets, names = parts(ax)
        if len(brackets) != len(names) or not names:
            raise SystemExit("Figure_R5_regions_heads: class brackets and names not found")
        head = max(inches(fig, t).y1 for t in ax.get_xticklabels() if t.get_visible() and t.get_text())
        low = min(inches(fig, ln).y0 for ln in brackets)
        dy = max(0.0, head + 0.05 - low)  # inches
        for a in brackets + names:
            a.set_transform(a.get_transform() + ScaledTranslation(0.0, dy, fig.dpi_scale_trans))
    fig.canvas.draw()
    for ax in heat:  # the title as far above the class names as before
        t = ax._left_title
        rise = max(inches(fig, n).y1 for n in parts(ax)[1]) + gap0[ax] - inches(fig, t).y0
        pad = ax.titleOffsetTrans.get_matrix()[1, 2] * 72 / fig.dpi + rise * 72
        ax._set_title_offset_trans(pad)
    no_clashes(fig, "Figure_R5_regions_heads")


def spread_rows(ax, centers, f: float) -> None:
    """Everything the axes draws in data coordinates off a row center (lines, arrows, texts) moved f times as far from
    the nearest center: the rows' entries spread apart, the rows themselves (their labels) in place. Lines and texts in
    other coordinates (guide lines across the axes, labels in axes coordinates) are left alone."""
    import numpy as np
    c = np.asarray(centers, float)

    def m(y):
        y = np.asarray(y, float)
        near = c[np.abs(y[..., None] - c).argmin(axis=-1)]
        return near + (y - near) * f

    for ln in ax.get_lines():
        if ln.get_transform() is ax.transData:
            ln.set_ydata(m(ln.get_ydata()))
    for t in ax.texts:
        if isinstance(t, matplotlib.text.Annotation):
            if t.xycoords == "data":
                t.xy = (t.xy[0], float(m(t.xy[1])))
            if t.anncoords == "data":
                t.xyann = (t.xyann[0], float(m(t.xyann[1])))
        elif t.get_transform() is ax.transData:
            x, y = t.get_position()
            t.set_position((x, float(m(y))))


def fix_spikes(fig: Figure, principal: bool = False) -> None:
    """Fig. S12 (R10): drawn at the text width and 6.0 in high (height cap 0.62; the report drew it 7.2 x 8.0 in), no
    text under 7.0 pt (panel titles 8, axis and band titles 7.5), laid out by hand: A with its columns of location counts
    and Holm p values, the three deeper bands of B below, each with its own row-label column. Within each anatomy the
    dense and site-matched entries are set a quarter row from its center instead of a fifth (everything of the entry
    moved together), so that the rows can be closer; 'neither reaches 50 %' keeps the line at a ratio of 1 off its
    letters with a thin white outline instead of a white box, which would cover the next entry at this spacing (as
    does '← favors Neuromag', which at this width reaches across that line)."""
    from matplotlib.patheffects import withStroke
    a = next(ax for ax in fig.axes if ax.get_title(loc="left").endswith("below the scalp"))
    cols = next(ax for ax in fig.axes if not ax.axison)
    bs = [ax for ax in fig.axes if ax.get_title().endswith(" mm")]
    if len(bs) != 3 or len(fig.legends) != 1:
        raise SystemExit("Figure_R10_spikes: panel A, its columns, three deeper bands and one key expected")
    centers = [float(y) for y in a.get_yticks()]
    for ax in [a, cols] + bs:
        spread_rows(ax, centers, 1.25)
    floor_fonts(fig, 7.0)
    a._left_title.set_fontsize(8.0)
    deeper = next(t for t in bs[0].texts if isinstance(t, matplotlib.text.Annotation) and t.get_text())
    deeper.set_fontsize(8.0)
    for ax in [a] + bs:
        ax.title.set_fontsize(7.5)
        ax.xaxis.label.set_fontsize(7.5)
        ax.tick_params(axis="both", which="both", labelsize=7.0)
        for t in ax.texts:
            if t.get_text() == "neither reaches 50 %":
                t.set_bbox(None)
                t.set_path_effects([withStroke(linewidth=2.2, foreground="white")])
            elif "favors" in t.get_text():  # '← favors Neuromag' now reaches across the line at 1
                t.set_path_effects([withStroke(linewidth=2.2, foreground="white")])
    fig.set_layout_engine("none")
    W, H = DRAWN_W, 6.0
    fig.set_size_inches(W, H)
    key = fig.legends[0]
    units_a, units_b = abs(a.get_ylim()[1] - a.get_ylim()[0]), abs(bs[0].get_ylim()[1] - bs[0].get_ylim()[0])
    l, r, gb, gap = 0.04, 0.03, 0.12, 0.12  # margins (the PDF's text runs a little wider), gap between the bands, A-B
    loc_texts = [t for t in cols.texts if t.get_transform() is cols.transData]
    heads = [t for t in cols.texts if t.get_transform() is not cols.transData]
    for _ in range(3):
        fig.canvas.draw()
        da, db = decorations(fig, a), decorations(fig, bs[0])
        col1 = max(inches(fig, t).width for t in loc_texts + heads if t.get_position()[0] < 0.5)
        col2 = max(inches(fig, t).width for t in loc_texts + heads if t.get_position()[0] > 0.5)
        w_c = max(2 * ((col1 + col2) / 2 + 0.1), 2 * col2 + 0.02)  # columns at 1/4 and 3/4 of its width
        x_c = W - r - w_c
        h_key = inches(fig, key).height
        t_a = max(da[2], max(inches(fig, t).y1 for t in heads) - inches(fig, cols).y1)
        t_b = max(inches(fig, deeper).y1, max(inches(fig, ax.title).y1 for ax in bs)) - inches(fig, bs[0]).y1
        d_b = max(decorations(fig, ax)[1] for ax in bs)
        u = (H - h_key - 0.06 - t_a - da[1] - gap - t_b - d_b) / (units_a + units_b)
        y = H - h_key - 0.06 - t_a - units_a * u
        place(a, l + da[0], y, x_c + w_c / 4 - col1 / 2 - 0.08 - l - da[0], units_a * u)  # 0.08 in before column 1
        place(cols, x_c, y, w_c, units_a * u)
        y -= da[1] + gap + t_b + units_b * u
        w_b = (W - r - l - db[0] - 2 * gb) / 3
        for i, ax in enumerate(bs):
            place(ax, l + db[0] + i * (w_b + gb), y, w_b, units_b * u)
        w_key = inches(fig, key).width
        mid = min(max((l + da[0] + W - r) / 2, w_key / 2 + 0.1), W - w_key / 2 - 0.1)
        key.set_bbox_to_anchor((mid / W, 1.0), transform=fig.transFigure)
    fig.canvas.draw()
    width = fig.get_tightbbox(fig.canvas.get_renderer()).width
    if width > DRAWN_W - 0.03:  # with room for the PDF's slightly wider text
        raise SystemExit(f"Figure_R10_spikes: {width:.3f} in drawn, wider than the text width")
    no_clashes(fig, "Figure_R10_spikes")


def fix_arrays(fig: Figure, principal: bool = False) -> None:
    """Fig. 1 (R15), drawn at the text width by draw_clean(): its sensor markers scaled with the panels (the ratio of
    their width to the report's; figure_arrays makes the three panels W - 0.62 - 0.04 - 2 x 0.10 in wide), so that they
    keep their size relative to the head and print as before."""
    f = (fig.get_size_inches()[0] - 0.86) / (REPORT_W - 0.86)
    n = 0
    for ax in fig.axes:
        for ln in ax.get_lines():
            if ln.get_marker() not in ("None", "none", "", None):
                ln.set_markersize(ln.get_markersize() * f)
                ln.set_markeredgewidth(ln.get_markeredgewidth() * f)
                n += 1
    if n != 6:
        raise SystemExit(f"Figure_R15_arrays: {n} sets of sensor markers, 6 expected")
    no_clashes(fig, "Figure_R15_arrays")


def fix_geometry(fig: Figure, principal: bool = False) -> None:
    """Fig. 2 (R12), drawn at the text width by draw_clean(): the sensor markers scaled with the panels (figure_geometry
    makes the three panels W - 0.66 - 0.04 - 2 x 0.10 in wide), so that they print as before; the note under each
    column title (gaps, scale factor, sites), wider than a column at this width, on four lines (the same words), and the
    figure made taller at the top to hold it. Child B's column title names Section S10, where the text reports it."""
    W, H = fig.get_size_inches()
    f = (W - 0.9) / (REPORT_W - 0.9)
    marked = [ln for ax in fig.axes for ln in ax.get_lines() if ln.get_marker() not in ("None", "none", "", None)]
    if len(marked) != 18:
        raise SystemExit(f"Figure_R12_geometry: {len(marked)} sets of sensor markers, 18 expected")
    for ln in marked:
        ln.set_markersize(ln.get_markersize() * f)
        ln.set_markeredgewidth(ln.get_markeredgewidth() * f)
    notes = [t for t in fig.texts if t.get_text().startswith("median gap:")]
    if len(notes) != 3:
        raise SystemExit(f"Figure_R12_geometry: {len(notes)} column notes, 3 expected")
    # child B, one of the heads the text reports separately: its column title says where (Section S10)
    title_b = [t for t in fig.texts if t.get_text() == style.ANAT_LABEL["childB"]]
    if len(title_b) != 1 or not style.ANAT_LABEL["childB"].endswith(")"):
        raise SystemExit(f"Figure_R12_geometry: one column title {style.ANAT_LABEL['childB']!r} expected")
    title_b[0].set_text(style.ANAT_LABEL["childB"][:-1] + "; Section S10)")
    for t in notes:
        s = t.get_text()
        t.set_text(rewrapped(s, {s: s.replace(", fitted ", ",\nfitted ", 1).replace("; ", ";\n", 1)}))
    fig.canvas.draw()
    dh = max(0.0, max(inches(fig, ax).y1 for ax in fig.axes) + 0.08 - min(inches(fig, t).y0 for t in notes))
    if dh:
        H2 = H + dh
        for ax in fig.axes:  # the panels kept where they are, in inches from the bottom
            b = ax.get_position(original=True)
            ax.set_position([b.x0, b.y0 * H / H2, b.width, b.height * H / H2])
        for t in fig.texts:  # the titles and notes kept as far from the top
            x, y = t.get_position()
            t.set_position((x, (H2 - (1 - y) * H) / H2))
        fig.set_size_inches(W, H2)
    fig.canvas.draw()
    if min(inches(fig, t).y0 for t in notes) < max(inches(fig, ax).y1 for ax in fig.axes) + 0.07:
        raise SystemExit("Figure_R12_geometry: the column notes reach the panels")
    no_clashes(fig, "Figure_R12_geometry")


def fix_adult_depth(fig: Figure, principal: bool = False) -> None:
    """Fig. 5 (R1): drawn by its script at the report's 7.2-in width, printed at about 0.89 of that: its 7-pt type (the
    counts per bin and the keys) raised to 7.35 pt (6.5 pt printed); the keys checked to stay clear of the curves (their
    boxes, border pads included, cover no plotted point) and the dB axes checked again. A's key, two columns, set with
    slightly shorter handles and column spacing to keep its footprint, and its longest entry on two lines (the same
    words): narrower, the key stands in the empty upper right of A instead of over the shallowest bins."""
    floor_fonts(fig, 7.35)
    a = next(ax for ax in fig.axes if ax.get_legend() is not None and ax.get_legend()._ncols == 2)
    two = {s.replace("\n", " "): s for s in ("Sensor, brain and room noise,\nroom field projected out",)}
    old = a.get_legend()
    rebuild_legend(old, labels=[rewrapped(t.get_text(), two) for t in old.get_texts()], loc="upper right", ncol=2,
                   fontsize=7.35, handlelength=2.2, columnspacing=1.0, handletextpad=0.6)
    fig.canvas.draw()
    for ax in fig.axes:
        lg = ax.get_legend()
        if lg is not None and legend_hits(ax, lg, full=True):
            raise SystemExit(f"Figure_R1_adult_depth: a key covers {legend_hits(ax, lg, full=True)} plotted points")
    sys.modules["report_figures_adult"].check_db_axes(fig)
    no_clashes(fig, "Figure_R1_adult_depth")


def fix_small_type(name: str, pt: float, top_room: bool = False):
    """A figure drawn by its script at the report's 7.2-in width and printed at about 0.88-0.90 of it, whose type under
    `pt` (its keys and notes at 7 pt) is raised to `pt` (6.5 pt printed); its keys checked to stay clear of the data
    (with top_room, the first axes' upper limit, shared by its panels, raised until they do) and its texts of each
    other."""
    def fix(fig: Figure, principal: bool = False) -> None:
        floor_fonts(fig, pt)
        for _ in range(12):
            fig.canvas.draw()
            hits = sum(legend_hits(ax, ax.get_legend()) for ax in fig.axes if ax.get_legend() is not None)
            if not hits:
                break
            if not top_room:
                raise SystemExit(f"{name}: a key covers {hits} plotted points")
            lo, hi = fig.axes[0].get_ylim()
            fig.axes[0].set_ylim(lo, hi + 0.04 * (hi - lo))
        else:
            raise SystemExit(f"{name}: the keys still cover plotted points")
        no_clashes(fig, name)
    return fix


def fix_sphere(fig: Figure, principal: bool = False) -> None:
    """Fig. S1 (R0): drawn by its script at the report's 7.2-in width, printed at about 0.89 of that: its 7-pt type (the
    notes) raised to 7.4 pt (6.5 pt printed, fix_small_type()); B's notes, larger now, moved in 1-mm steps until they
    clear the curves and the key by 1.5 pt: the note on where Neuromag leads deeper (down), the other shallower (up)."""
    fix_small_type("Figure_R0_sphere", 7.4)(fig)
    bx = next(ax for ax in fig.axes if ax.get_title(loc="left").endswith("Realistic adult head"))
    for start, step in (("Neuromag ahead at every depth", 1.0), ("η ≤", -1.0)):
        note = next(t for t in bx.texts if t.get_text().startswith(start))
        for _ in range(12):
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            drawn = plotted_boxes(bx, skip=[note]) + [bx.get_legend().get_window_extent(r)]
            if not any(overlaps(i, b, 1.5 * fig.dpi / 72) for i in text_ink(note, r) for b in drawn):
                break
            note.set_y(note.get_position()[1] + step)  # depth axis, deeper downwards
        else:
            raise SystemExit(f"Figure_R0_sphere: no place for B's note {start!r}")
    no_clashes(fig, "Figure_R0_sphere")


def fix_regions_adult(fig: Figure, principal: bool = False) -> None:
    """Fig. 7 (R4): drawn at the text width, 6.2 in high as the report drew it (under its height cap of 0.66 the width
    binds), no text under 7.0 pt, panel titles 8 pt and axis labels 7.5 pt; the constrained layout re-flows the four
    panels. The ratio axes of A and B are labelled at 0.8, 1, 1.2 and 1.4, with 0.9, 1.1 and 1.3 as unlabelled ticks
    (as in Fig. 3: at this width the labels would touch). The dB axes are checked again."""
    panels = [ax for ax in fig.axes if ax.get_title(loc="left")]
    if len(panels) != 4:
        raise SystemExit(f"Figure_R4_regions_adult: {len(panels)} titled panels, 4 expected")
    floor_fonts(fig, 7.0)
    for ax in panels[:2]:
        ratio_ticks(ax, (0.8, 1, 1.2, 1.4), (0.9, 1.1, 1.3))
    for ax in fig.axes:
        ax._left_title.set_fontsize(8.0)
        for lab in (ax.xaxis.label, ax.yaxis.label):
            lab.set_fontsize(min(lab.get_fontsize(), 7.5))
    panels[3].xaxis.label.set_fontsize(7.0)  # (D): its label, longer than the panel is wide, given more room
    panels[0].get_subplotspec().get_gridspec().set_width_ratios([1.6, 1.6, 1.0, 1.15])
    fit_constrained(fig, DRAWN_W - 0.02, 6.2)  # text measured on the screen canvas comes out a little wider in the PDF
    sys.modules["report_figures_adult"].check_db_axes(fig)
    no_clashes(fig, "Figure_R4_regions_adult")


FIXES = {"Figure_R15_arrays": fix_arrays, "Figure_R12_geometry": fix_geometry,                    # Figs. 1, 2
         "Figure_R3_conditions": fix_conditions, "Figure_R17_noise_checks": fix_noise_checks,         # Figs. 3, 4
         "Figure_R1_adult_depth": fix_adult_depth, "Figure_R4_regions_adult": fix_regions_adult,      # Figs. 5, 7
         "Figure_R6_pediatric_D": fix_pediatric_d, "Figure_R13_maps_heads": fix_maps_heads,           # Figs. 8, S16
         "Figure_helmet_fit": fix_helmet_fit, "Figure_R16_confirm": fix_confirm,                      # Figs. 9, 10, S17, S18
         "Figure_R2_noise_model": fix_noise_model, "Figure_covariance_structure": fix_covariance_structure,  # S2, S3
         "Figure_covariance_bands": fix_covariance_bands,                                             # Fig. S4
         "Figure_noise_sensitivity": fix_noise_sensitivity, "Figure_noise_sensitivity_depth": fix_noise_depth,  # S5
         "Figure_S_children_qc": fix_children_qc, "Figure_R7_helmet_fit": fix_helmet_routes,          # Figs. S6, S8
         "Figure_R5_regions_heads": fix_regions_heads, "Figure_R10_spikes": fix_spikes,               # Figs. S9, S12
         "Figure_S_localization_effects": fix_supplement_spikes,                                      # Fig. S14
         "Figure_S_joint_detection_localization": fix_supplement_spikes,                              # Fig. S15
         "Figure_R11_maps_adult": fix_maps("Figure_R11_maps_adult"),                                  # Fig. 6
         "Figure_R14_maps_scaled": fix_maps("Figure_R14_maps_scaled"),                                # Fig. S10
         "Figure_R0_sphere": fix_sphere,                                                              # Fig. S1
         "Figure_R8_depth_matched": fix_small_type("Figure_R8_depth_matched", 7.5, top_room=True),    # Fig. S7
         # Fig. S11: its key, with the 18-month template labelled misregistered, is wider than the panels; the 6.5-pt
         # labels above the panels raised to 6.8 pt (6.5 pt printed)
         "Figure_R9_noise_floor": fix_small_type("Figure_R9_noise_floor", 6.8),
         "Figure_S_detection_curves_adult": fix_small_type("Figure_S_detection_curves_adult", 7.35)}  # Fig. S13


def stack(parts: list[Path], out: Path):
    """One PDF holding the parts top to bottom at their own sizes (vector content kept)."""
    if len(parts) == 1:
        shutil.copy2(parts[0], out)
        return
    tex = EXPORT / f"{out.stem}_stack.tex"
    body = "\\\\[3mm]\n".join(f"\\includegraphics{{{p.name}}}" for p in parts)
    tex.write_text("\\documentclass[varwidth=\\maxdimen,border=0pt]{standalone}\n\\usepackage{graphicx}\n"
                   "\\begin{document}\n\\centering\n" + body + "\n\\end{document}\n")
    r = subprocess.run(["tectonic", "-X", "compile", tex.name], cwd=EXPORT, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"stacking {out.name} failed:\n{r.stderr[-2000:]}")
    shutil.copy2(EXPORT / f"{tex.stem}.pdf", out)


def page_size(pdf: Path) -> tuple[float, float]:
    """Width and height of a one-page PDF in bp (pdfinfo)."""
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
    w, h = re.search(r"Page size:\s+([\d.]+) x ([\d.]+) pts", out).groups()
    return float(w), float(h)


def stack_layout(name: str, parts: list[Path]) -> dict:
    """Where each part sits in the stacked figure (bp, from the top; parts are centered), so that a document can show
    one part alone by trimming the others away."""
    W, H = page_size(FIGDIR / f"{name}.pdf")
    sizes = [page_size(p) for p in parts]
    gap = (H - sum(h for _, h in sizes)) / max(len(parts) - 1, 1)
    out, y = [], 0.0
    for p, (w, h) in zip(parts, sizes):
        out.append({"part": p.stem, "width": round(w, 3), "height": round(h, 3), "top": round(y, 3)})
        y += h + gap
    return {"width": round(W, 3), "height": round(H, 3), "parts": out}


def print_audit(fonts: dict, clash: dict) -> tuple[list[str], list[str], list[str]]:
    """Per manuscript figure (per part when its parts are printed on separate pages): page size, print scale under the
    \\includegraphics options of PRINT (width = the text width; height capped where a cap is given; aspect kept), the
    smallest text drawn into it, drawn and printed, and its texts that overlap or touch (clashes()). Returns the
    report's lines, the figures (parts) of LEGIBLE that print text under the size LEGIBLE wants for them, and the
    figures (parts) with texts that overlap or touch."""
    lines, short = [f"text width {TEXT_W:.1f} bp, text height {TEXT_H:.1f} bp; smallest printed text wanted: "
                    f"{MIN_PRINTED_PT:g} pt, or the size given where the layout leaves no more room (sub- and superscripts "
                    "of math labels not counted); texts that overlap or touch: clashes()"], []
    overlapping = []
    for name, parts in {**MAIN, **SUPP}.items():
        cap = PRINT.get(name)
        if isinstance(cap, tuple):  # one part per page (build_paper.py's trim())
            groups = [(f"{name} part {i + 1} ({p})", [p], page_size(EXPORT / f"{p}.pdf"), c, i)
                      for i, (p, c) in enumerate(zip(parts, cap))]
        else:
            groups = [(name, parts, page_size(FIGDIR / f"{name}.pdf"), cap, None)]
        for label, ps, (w, h), c, part in groups:
            scale = min(TEXT_W / w, TEXT_H * c / h) if c else TEXT_W / w
            drawn = [(s, t) for p in ps for s, t in fonts.get(p, [])]
            if not drawn:
                lines.append(f"{label}: no font record (drawn with --no-draw before?)")
                continue
            size, text = min(drawn)
            want = LEGIBLE.get(name)
            want = want[part] if isinstance(want, tuple) else want
            under = size * scale < (want or MIN_PRINTED_PT) - 1e-6
            flag = ("  << UNDER" if want else "  (not revised)") if under else ""
            found = [c for p in ps for c in clash.get(p, ["not checked"])]
            lines.append(f"{label}: {w:.1f} x {h:.1f} bp, " + (f"height cap {c:g}" if c else "no height cap")
                         + f", print scale {scale:.3f}; smallest text {size:.2f} pt drawn, {size * scale:.2f} pt printed "
                         f"({text[:50]!r})" + (f" [wanted {want:g} pt]" if want and want != MIN_PRINTED_PT else "") + flag
                         + "; texts overlapping: " + ("; ".join(found) + "  << OVERLAP" if found else "none"))
            if want and under:
                short.append(label)
            if found:
                overlapping.append(label)
    return lines, short, overlapping


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-draw", action="store_true", help="only stack the already drawn PDFs")
    args = ap.parse_args(argv)
    if not args.no_draw:
        draw_all()
        (EXPORT / "text_audit.txt").write_text("\n".join(TEXT_AUDIT) + "\n")
        print(f"text left to fix in the figures: {len(TEXT_AUDIT)} (build/figure_export/text_audit.txt)")
    FIGDIR.mkdir(exist_ok=True)
    layouts = {}
    for name, parts in {**MAIN, **SUPP}.items():
        paths = [EXPORT / f"{p}.pdf" for p in parts]
        missing = [p.name for p in paths if not p.exists()]
        if missing:
            raise SystemExit(f"{name}: not drawn: {missing}")
        stack(paths, FIGDIR / f"{name}.pdf")
        print(f"{name}.pdf <- {', '.join(parts)}")
        if len(paths) > 1:
            layouts[name] = stack_layout(name, paths)
    (HERE / "figure_stacks.json").write_text(json.dumps(layouts, indent=1) + "\n")
    fonts, clash = ({}, {}) if not (EXPORT / "font_audit.json").exists() else (
        json.loads((EXPORT / "font_audit.json").read_text()), json.loads((EXPORT / "clash_audit.json").read_text()))
    lines, short, overlapping = print_audit(fonts, clash)
    (EXPORT / "print_audit.txt").write_text("\n".join(lines) + "\n")
    print(f"figures revised for legibility printing text under the size wanted: {len(short)} "
          f"(build/figure_export/print_audit.txt){': ' + ', '.join(short) if short else ''}")
    print(f"figures with texts that overlap or touch: {len(overlapping)}"
          f"{': ' + ', '.join(overlapping) if overlapping else ''} (build/figure_export/print_audit.txt)")


if __name__ == "__main__":
    main()
