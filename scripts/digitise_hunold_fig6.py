#!/usr/bin/env python3
"""Digitise Hunold et al. (2016) Fig. 6 (p. 1157) for an absolute calibration of the background.

Reads the locally held reference PDF (never committed), extracts the embedded 200-ppi raster of
page 13 with poppler's pdfimages, measures the three scale-bar brackets and, for each of the 24
traces, the per-column centroid of the drawn line. Prints the baseline standard deviation (pixels,
display time < 0.90 s) and the spike-window extreme, per trace and per channel, and the values in
physical units. The constants in src/opmsquid/hunold.py (FIG6_*) were copied from this output.

Usage: .venv/bin/python scripts/digitise_hunold_fig6.py Hunold_2016_Physiol._Meas._37_1146.pdf
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from opmsquid import hunold  # noqa: E402

BLOCKS = ("eeg", "mag", "grad")
ROWS = ("dipole superficial", "dipole deep", "patch superficial", "patch deep")
# panel x extents [px] and t = 0 column (the time axis has ticks at 0/1/2 s, 130.7 px/s)
SIDES = {"radial": (74, 333, 71.5), "tangential": (404, 663, 402.0)}
BRACKET_COLS = {"eeg": (500, 512), "mag": (500, 512), "grad": (500, 512)}
BRACKET_ROWS = {"eeg": (425, 475), "mag": (880, 935), "grad": (1340, 1384)}


def extract(pdf: Path) -> np.ndarray:
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdfimages", "-f", "13", "-l", "13", "-png", str(pdf), f"{tmp}/fig"], check=True)
        return np.asarray(Image.open(f"{tmp}/fig-000.png").convert("L")).astype(float)


def axis_rows(g: np.ndarray, xa: int, xb: int) -> list[float]:
    """Rows of the 12 horizontal trace axes (darkness-weighted centre of two-row runs)."""
    sub = g[:, xa:xb]
    d = (255.0 - sub).mean(axis=1) / 255.0
    cover = (np.minimum(sub[:-1], sub[1:]) < 180).mean(axis=1)  # continuous line within a 2-row band
    lines = np.flatnonzero(cover > 0.9)
    groups = np.split(lines, np.flatnonzero(np.diff(lines) > 1) + 1)
    out = []
    for r in groups:
        if len(r):
            rr = np.arange(r[0], r[-1] + 2)
            out.append(float(np.average(rr, weights=d[rr])))
    return out


def bracket_height(g: np.ndarray, block: str) -> float:
    """Distance between the top and bottom serifs of the '[' scale bar [px]."""
    (x0, x1), (y0, y1) = BRACKET_COLS[block], BRACKET_ROWS[block]
    sub = 255.0 - g[y0:y1, x0:x1]
    serif = sub[:, sub.sum(axis=0).argmax() + 1:].sum(axis=1)  # the serifs extend right of the stem
    rows = np.flatnonzero(serif > 0.3 * serif.max())
    groups = np.split(rows, np.flatnonzero(np.diff(rows) > 2) + 1)
    top, bottom = groups[0], groups[-1]
    centre = [float(np.average(r, weights=serif[r])) for r in (top, bottom)]
    return centre[1] - centre[0]


def trace_centroids(g: np.ndarray, axis_row: float, side: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per pixel column: time [s], darkness-weighted centroid, and the upper and lower edge of the
    drawn line (up = positive, px)."""
    xa, xb, x0 = SIDES[side]
    top, bot = int(axis_row - 46), int(axis_row - 2)
    t, c, hi, lo = [], [], [], []
    for x in range(xa, xb):
        col = g[top:bot, x]
        ys = np.flatnonzero(col < 200)
        if len(ys):
            w = 255.0 - col[ys]
            t.append((x - x0) / hunold.FIG6_PX_PER_S)
            c.append(-float(np.sum((ys + top) * w) / w.sum()))
            hi.append(-float(ys.min() + top))
            lo.append(-float(ys.max() + top))
    return np.array(t), np.array(c), np.array(hi), np.array(lo)


def centroid_sd_bias() -> float:
    """Rendered-centroid SD / true SD for a Hunold-spectrum background drawn at ~3 px SD (the
    paper's baselines): the per-column centroid averages 7.65 ms and reads the SD low."""
    x = hunold.background_timecourses(8, 6000, 1000.0, np.random.default_rng(0))
    ratios = []
    for tr in x:
        tr = tr / tr.std()
        per = 900
        true = np.sqrt(np.mean([np.var(tr[i:i + per]) for i in range(0, 6000 - per + 1, per)]))
        ratios.append(hunold.rendered_centroid_sd(tr, 1000.0, px_per_unit=3.0) / (3.0 * true))
    return float(np.mean(ratios))


def main(pdf: Path):
    g = extract(pdf)
    rows = axis_rows(g, 74, 330)
    assert len(rows) == 12, rows
    bars = {b: round(bracket_height(g, b), 2) for b in BLOCKS}
    print("scale-bar brackets [px]:", bars)
    per_channel: dict[tuple[str, str], list[float]] = {}
    bias = centroid_sd_bias()
    cand = {"noisy peak (centroid) / 2<e>": [], "p2p (centroid) / 2<e>": [], "p2p (centroid) / 2<e>, SD bias-corrected": [],
            "p2p (drawn-line extremes) / 2<e>": [], "noisy peak (drawn-line extremes) / 2<e>": []}
    print(f"{'block':5s} {'source':20s} {'side':10s} {'channel':7s} {'SNR':>5s} {'sd_px':>6s} {'peak_px':>8s} {'p2p_px':>7s} "
          f"{'p2p_line':>8s}")
    for b, block in enumerate(BLOCKS):
        for r, src in enumerate(ROWS):
            for side in SIDES:
                t, c, hi, lo = trace_centroids(g, rows[4 * b + r], side)
                base = c[t < 0.90]
                w = (t >= 0.90) & (t < 1.15)
                mu, sd = base.mean(), base.std()
                peak, p2p = np.abs(c[w] - mu).max(), c[w].max() - c[w].min()
                p2p_line = hi[w].max() - lo[w].min()
                peak_line = max(hi[w].max() - mu, mu - lo[w].min())
                ch, snr = hunold.FIG6_TRACES[block][(src, side)]
                per_channel.setdefault((block, ch), []).append(sd)
                print(f"{block:5s} {src:20s} {side:10s} {ch:7s} {snr:5.2f} {sd:6.2f} {peak:8.1f} {p2p:7.1f} {p2p_line:8.1f}")
                if block != "eeg" and snr > 2.4:
                    e2 = 2 * sd * np.sqrt(np.pi / 2)  # 2 <|hilbert|> for a Gaussian baseline
                    cand["noisy peak (centroid) / 2<e>"].append(peak / e2 / snr)
                    cand["p2p (centroid) / 2<e>"].append(p2p / e2 / snr)
                    cand["p2p (centroid) / 2<e>, SD bias-corrected"].append(p2p / (e2 / bias) / snr)
                    cand["p2p (drawn-line extremes) / 2<e>"].append(p2p_line / e2 / snr)
                    cand["noisy peak (drawn-line extremes) / 2<e>"].append(peak_line / e2 / snr)
    print(f"candidate SNR / printed SNR, MEG traces with printed SNR > 2.4 (n = {len(cand['p2p (centroid) / 2<e>'])}); "
          f"centroid SD bias factor {bias:.3f}:")
    for k, v in cand.items():
        print(f"  {k:42s} {np.mean(v):.2f} +/- {np.std(v):.2f}")
    print("per channel: mean baseline sd [px] over the traces sharing that channel's background")
    for (block, ch), v in per_channel.items():
        unit = hunold.FIG6_SCALE_BAR[block] / bars[block]
        print(f"  {block:5s} {ch:6s} n={len(v)} sd={np.mean(v):.2f} px = {np.mean(v) * unit:.3e} (SI)")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
