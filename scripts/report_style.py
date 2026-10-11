"""Shared look of the report's new figures (scripts/report_figures_*.py): open fonts only, one
colour and marker per array and per anatomy (colour-blind safe: the Okabe-Ito colours, with brown and grey for
children B and C), the same sizes.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.legend_handler import HandlerPatch, HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

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

# arrays (colour) and Neuromag channel sets; the labels are the manuscript's terms
ARRAY_COLOR = {"opm_dense": "#0072B2", "opm_matched": "#56B4E9", "opm204": "#56B4E9",
               "squid": "#000000", "combined": "#000000", "mag": "#999999", "grad": "#555555"}
ARRAY_LABEL = {"opm_dense": "dense OPM array", "opm_matched": "site-matched OPM array", "squid": "Neuromag",
               "combined": "Neuromag (306 channels)", "mag": "Neuromag magnetometers (102)", "grad": "Neuromag gradiometers (204)"}
ARRAY_SHORT = {"opm_dense": "dense OPM", "opm_matched": "site-matched OPM"}  # the manuscript's short forms

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


class KeyPair(tuple):
    """A legend key of two handles set side by side, one per OPM array (drawn by KEY_HANDLERS)."""


class BandLine(Patch):
    """The legend handle of a reference drawn as a shaded band (an interval) with a line at its middle (the median)."""


class _HandlerBandLine(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        band = Rectangle((-xdescent, -ydescent), width, height, facecolor=orig_handle.get_facecolor(),
                         edgecolor="none", transform=trans)
        x = -xdescent + width / 2
        line = Line2D([x, x], [-ydescent, -ydescent + height], color=orig_handle.get_edgecolor(),
                      lw=orig_handle.get_linewidth(), transform=trans)
        return [band, line]


# legend handler map of the composite keys above (pass it to every legend that holds one)
KEY_HANDLERS = {KeyPair: HandlerTuple(ndivide=None, pad=0.3), BandLine: _HandlerBandLine()}


def band_key(array: str, alpha: float, label: str | None = None) -> Patch:
    """The key of a band shaded in an array's colour at `alpha`, as drawn, with a thin outline so that it shows."""
    c = ARRAY_COLOR[array]
    return Patch(facecolor=to_rgba(c, alpha), edgecolor=to_rgba(c, 0.6), lw=0.6, label=label)


def band_pair(alpha: float, arrays=("opm_dense", "opm_matched")) -> KeyPair:
    """One key for the bands of both OPM arrays: one swatch in each array's colour."""
    return KeyPair(band_key(a, alpha) for a in arrays)


def apply() -> None:
    plt.rcParams.update(RC)


def commit() -> str:
    """Short HEAD commit, '+dirty' when the figure code or its inputs have uncommitted changes (the figures' own
    outputs in results/report/ do not count)."""
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "scripts", "results", ":!results/report"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
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
