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

Outputs: paper/figures/Figure_1.pdf ... Figure_8.pdf (main text) and Figure_S1.pdf ... (supplementary material), each
with a 300-dpi PNG preview beside it in paper/build/figure_export/; paper/figure_stacks.json gives where each part of a
stacked figure sits, so that a document can show one part alone (build_paper.py's trim()).

Usage: .venv/bin/python paper/export_figures.py [--only NAME ...]
"""
from __future__ import annotations

import argparse
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
import matplotlib.colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.spines  # noqa: E402
import matplotlib.ticker  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.text import Text  # noqa: E402

import report_style as style  # noqa: E402

EXPORT = HERE / "build" / "figure_export"   # drawing scripts' outputs and previews (not committed)
FIGDIR = HERE / "figures"                   # the manuscript's figure files
PDF_DPI, PNG_DPI = 600, 300
PANEL = re.compile(r"^(?:\(([a-i])\)(?:\s+|$)|([a-i])(?:\s{2,}|$))(.*)$", re.S)  # "(a) Title", "a  Title" or "a"
# figures drawn without panel letters whose panels the manuscript letters (continuing the figure they are stacked under):
# letters given to the axes that carry a left title, in reading order
EXTRA_LETTERS = {"Figure_noise_sensitivity_depth": "ghij"}

# manuscript figure -> drawn figures stacked top to bottom (names as the drawing scripts save them)
MAIN = {"Figure_1": ["Figure_R15_arrays", "Figure_R12_geometry"], "Figure_2": ["Figure_R3_conditions"],
        "Figure_3": ["Figure_R17_noise_checks"], "Figure_4": ["Figure_R1_adult_depth", "Figure_R11_maps_adult"],
        "Figure_5": ["Figure_R4_regions_adult"], "Figure_6": ["Figure_R6_pediatric_D", "Figure_R13_maps_heads"],
        "Figure_7": ["Figure_constant_gap"], "Figure_8": ["Figure_R16_confirm"]}
SUPP = {"Figure_S1": ["Figure_R0_sphere"], "Figure_S2": ["Figure_R2_noise_model"],
        "Figure_S3": ["Figure_covariance_structure"], "Figure_S4": ["Figure_covariance_bands"],
        "Figure_S5": ["Figure_noise_sensitivity", "Figure_noise_sensitivity_depth"], "Figure_S6": ["Figure_S_children_qc"],
        "Figure_S7": ["Figure_R8_depth_matched"], "Figure_S8": ["Figure_R7_helmet_fit"], "Figure_S9": ["Figure_R5_regions_heads"],
        "Figure_S10": ["Figure_R14_maps_scaled"], "Figure_S11": ["Figure_R9_noise_floor"], "Figure_S12": ["Figure_R10_spikes"],
        "Figure_S13": ["Figure_S_detection_curves_adult"], "Figure_S14": ["Figure_S_localization_effects"],
        "Figure_S15": ["Figure_S_joint_detection_localization"]}


def relabel(fig: Figure) -> int:
    """'(a) Title' -> bold 'A' + title, for every text of the figure that starts with a panel label."""
    n = 0
    for t in fig.findobj(Text):
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
         "OPM matched": "site-matched OPM"}
# a panel referred to inside a text, "(c)", "(b, c)", "(a-c)" or "(d and e)" -> upper case, as the panel letters are
PANEL_REF = re.compile(r"\(([a-i](?:(?:,\s*|\s*[-\u2013]\s*|\s+and\s+)[a-i])*)\)")


def americanize_text(s: str) -> str:
    """British spellings -> American; lower-case panel references -> upper case."""
    def swap(m):
        word = AMERICAN[m.group(1).lower()]
        return word[0].upper() + word[1:] if m.group(1)[0].isupper() else word
    s = BRITISH.sub(swap, s)
    for internal, term in TERMS.items():
        s = s.replace(internal, term).replace(internal[0].upper() + internal[1:], term[0].upper() + term[1:])
    return PANEL_REF.sub(lambda m: "(" + re.sub(r"\b[a-i]\b", lambda k: k.group(0).upper(), m.group(1)) + ")", s)


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


def savefig(self, fname, *args, **kwargs):
    """Every figure saved during the export: panel labels converted, spelling made American, notes removed, then a
    vector PDF and a PNG preview in EXPORT (cropped to what is drawn)."""
    stem = Path(str(fname)).stem
    add_letters(self, EXTRA_LETTERS.get(stem, ""))
    relabel(self)
    americanize(self)
    strip_notes(self)
    if getattr(self, "_suptitle", None) is not None:  # the caption carries the figure's title, not the figure
        self._suptitle.set_visible(False)
    EXPORT.mkdir(parents=True, exist_ok=True)
    kwargs.pop("dpi", None)
    kwargs.update(bbox_inches="tight", pad_inches=0.03)
    _orig_savefig(self, EXPORT / f"{stem}.pdf", *args, dpi=PDF_DPI, **kwargs)
    _orig_savefig(self, EXPORT / f"{stem}.png", *args, dpi=PNG_DPI, **kwargs)
    TEXT_AUDIT.extend(audit_text(self, stem))


def save(fig, name: str) -> Path:
    """report_style.save during the export: returns the preview PNG, whose pixel size at style.DPI (= PNG_DPI) some
    scripts check against the page width."""
    fig.savefig(EXPORT / f"{name}.pdf")
    plt.close(fig)
    return EXPORT / f"{name}.png"


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


def draw_all():
    Figure.savefig = savefig
    style.OUT = EXPORT
    style.save = save
    style.DPI = PNG_DPI  # the preview's resolution; the PDFs' raster layers use PDF_DPI
    run_main("report_figures_adult", [])
    run_main("report_figures_clean", [])
    run_main("report_figures_noise", [])
    run_main("report_figures_pediatric", [])
    load("report_figures_confirm").main(["--out-dir", str(EXPORT)])
    run_main("report_figures_qc", [])
    run_main("report_figures_supplement", [])
    cg = load("study_g3b_constant_gap")
    summary = json.loads((ROOT / "results/g3b_constant_gap/g3b_constant_gap_summary.json").read_text())
    g3b = json.loads((ROOT / "results/g3b/g3b_summary.json").read_text())
    cg.figure(summary, g3b, EXPORT / "Figure_constant_gap.png")
    cv = load("study_covariance_validation")
    cvs = json.loads((ROOT / "results/g2_covariance_validation/covariance_validation.json").read_text())
    cv.figure_structure(cvs, EXPORT / "Figure_covariance_structure.png")
    cv.figure_bands(cvs, EXPORT / "Figure_covariance_bands.png")
    ns = load("study_noise_sensitivity")
    nss = json.loads((ROOT / "results/g2_noise_sensitivity/noise_sensitivity_summary.json").read_text())
    ns.figures(nss, EXPORT)


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


if __name__ == "__main__":
    main()
