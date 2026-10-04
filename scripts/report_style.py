"""Shared look of the report's new figures (scripts/report_figures_*.py): open fonts only, one
colour and marker per array and per anatomy (Okabe-Ito, colour-blind safe), the same sizes.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "report"
DPI = 200
FULL_W, HALF_W = 7.2, 3.5  # inches: full width and one column

RC = {
    "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans", "font.size": 8.5,
    "axes.titlesize": 9, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "savefig.dpi": DPI, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    "figure.dpi": DPI,
}

# arrays (colour) and Neuromag channel sets
ARRAY_COLOR = {"opm_dense": "#0072B2", "opm_matched": "#56B4E9", "opm204": "#56B4E9",
               "squid": "#000000", "combined": "#000000", "mag": "#999999", "grad": "#555555"}
ARRAY_LABEL = {"opm_dense": "OPM dense", "opm_matched": "OPM matched", "squid": "Neuromag",
               "combined": "Neuromag (306)", "mag": "Neuromag magnetometers (102)", "grad": "Neuromag gradiometers (204)"}

# anatomies, in the order of the pediatric results; class decides the marker
ANAT_ORDER = ["adult", "school", "size2yr", "infant2yr", "infant18mo", "infant12mo", "childA", "childB", "childC"]
ANAT_LABEL = {"adult": "Adult", "school": "School-age size (scaled adult)", "size2yr": "2-year size (scaled adult)",
              "infant2yr": "24-month template", "infant18mo": "18-month template", "infant12mo": "12-month template",
              "childA": "Child A (7.8 y)", "childB": "Child B (8.3 y)", "childC": "Child C (8.7 y)"}
ANAT_SHORT = {"adult": "Adult", "school": "School size", "size2yr": "2-y size", "infant2yr": "24 mo", "infant18mo": "18 mo",
              "infant12mo": "12 mo", "childA": "Child A", "childB": "Child B", "childC": "Child C"}
ANAT_CLASS = {"adult": "adult", "school": "scaled", "size2yr": "scaled", "infant2yr": "template", "infant18mo": "template",
              "infant12mo": "template", "childA": "child", "childB": "child", "childC": "child"}
CLASS_LABEL = {"adult": "adult", "scaled": "scaled adult (size only)", "template": "infant template", "child": "individual child"}
CLASS_MARKER = {"adult": "D", "scaled": "s", "template": "^", "child": "o"}
ANAT_COLOR = {"adult": "#000000", "school": "#56B4E9", "size2yr": "#0072B2", "infant2yr": "#D55E00", "infant18mo": "#E69F00",
              "infant12mo": "#CC79A7", "childA": "#009E73", "childB": "#8C564B", "childC": "#7F7F7F"}


def apply() -> None:
    plt.rcParams.update(RC)


def commit() -> str:
    """Short HEAD commit, '+dirty' when the figure code or its inputs have uncommitted changes."""
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "scripts", "results"], cwd=ROOT, capture_output=True,
                           text=True).stdout.strip()
    return head + ("+dirty" if dirty else "")


def save(fig, name: str) -> Path:
    """Write results/report/<name>.png and close the figure."""
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def write_provenance(filename: str, entries: dict) -> None:
    """results/report/<filename>: {figure: {inputs, description}} plus the commit that drew them."""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / filename).write_text(json.dumps({"provenance": {"commit": commit()}, "figures": entries}, indent=1))
