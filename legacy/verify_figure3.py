#!/usr/bin/env python3
"""
Quantitative comparison of the replicated Figure 3 with the figure published in
bioRxiv 2026.08.17.744953 (page 14 of the preprint PDF).

The preprint stores Fig. 3 as lossless-grey and JPEG raster strips, placed at 1097 ppi on
the page; stroke widths and dash periods show that the artwork was rasterized at 1200 px
per figure inch.  This script

1. extracts the raster from the PDF (poppler ``pdfimages``) and stitches the strips,
2. loads the replica written by ``replicate_figure3.py`` (same 1200-ppi canvas, so both
   images share one pixel grid; it is rendered first if missing) and passes it through the
   preprint's image storage (same strips, JPEG quality 90, 4:2:0) so that both images carry
   the same compression artefacts,
3. measures both images with the *same* code:
   - axis calibration (tick marks, spines) and y-axis limits,
   - digitized curves (OPM blue, SQUID red), vertical and curve-normal deviations, also for
     a replica rendered without the fitted sub-pixel raster registration (physics only),
   - position of the dotted equal-SNR line (searched over the whole panel),
   - per-colour ink masks (blue / red / black / legend frame grey): IoU, 2-px tolerance,
   - connected components of black ink (text, arrows), matched one-to-one,
   - sub-pixel registration (normalised cross-correlation) of a representative set of
     elements: text runs, arrowheads, some tick labels, legend parts,
4. writes ``verification_report.json``, an overlay, a difference image and a side-by-side
   image, and prints PASS/FAIL checks (the exit status is 0 only if all pass).

``--self-test`` runs negative controls (spurious blob, missing curve, shifted frame, wrong
d_eq marker, wrong eta) and confirms that each one makes the checks fail.

Usage
-----
    python verify_figure3.py
    python verify_figure3.py --pdf 2026.08.17.744953.full.pdf \
        --replica figure3_output/Figure3_replicated.png --outdir figure3_output
    python verify_figure3.py --self-test
"""
from __future__ import annotations

import argparse
import glob
import io
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
PAGE = 14
# the preprint stores Fig. 3 as horizontal strips: (height, kind) with kind "L" = lossless
# grey, "J" = JPEG quality 90 with 4:2:0 chroma subsampling
PUBLISHED_STRIPS = [(214, "L"), (214, "J"), (214, "J"), (214, "J"), (214, "J"), (214, "J"),
                    (428, "J"), (428, "J"), (214, "J"), (428, "J"), (425, "L")]
LEGEND_BOX = (1000, 290, 2200, 880)  # x0, y0, x1, y1 of the legend in panel A [px]
LEGEND_FRAME = (1028, 316, 2165, 855)  # frame stroke centre-line [px]

# regions (y0, y1, x0, x1) [px] of individual elements for sub-pixel registration
ELEMENTS = {
    "panel letter A": (20, 145, 5, 120), "title A": (15, 180, 295, 1030),
    "panel letter B": (20, 145, 2458, 2555), "title B": (15, 180, 2770, 4230),
    "x-label A 'Dipole depth ('": (3030, 3195, 585, 1578), "x-label A 'd'": (3035, 3160, 1578, 1662),
    "x-label A ') [mm]'": (3030, 3195, 1700, 2110), "x-label B 'Dipole depth ('": (3030, 3195, 3010, 4002),
    "x-label B 'd'": (3035, 3160, 4002, 4086), "x-label B ') [mm]'": (3030, 3195, 4124, 4534),
    "legend 'OPM'": (375, 505, 1535, 1870), "legend 'SQUID'": (618, 760, 1535, 1995),
    "legend frame corner": (300, 345, 1015, 1060),
    "annotation A 'S_OPM'": (1030, 1222, 1045, 1325), "annotation A '> S_SQUID'": (1030, 1222, 1360, 1845),
    "annotation B1 'SNR_OPM'": (960, 1150, 3352, 3700), "annotation B1 '> SNR_SQUID'": (960, 1150, 3870, 4550),
    "annotation B2 'SNR_SQUID'": (1570, 1760, 3382, 3930), "annotation B2 '> SNR_OPM'": (1570, 1760, 3970, 4580),
    "arrow A left head": (1200, 1290, 605, 730), "arrow A vertex": (1170, 1225, 950, 1005),
    "arrow A down head": (1990, 2115, 1040, 1115), "arrow B1 head": (1215, 1295, 3060, 3180),
    "arrow B1 tail": (1150, 1180, 3300, 3345), "arrow B2 head": (2200, 2325, 3810, 3890),
    "arrow B2 tail": (1865, 1900, 3955, 3985),
    "y tick label A '4'": (645, 755, 0, 80), "y tick label B '4'": (518, 628, 2432, 2515),
    "x tick label A '20'": (2795, 2908, 548, 708), "x tick label B '80'": (2795, 2908, 4282, 4443),
    "spine top A": (215, 260, 180, 210),
}
# acceptance thresholds (1 px of the 1200-ppi raster = 0.021 mm of figure, 0.046 mm of depth)
MAX_CURVE_MEDIAN_PX = 0.5
MAX_DEQ_LINE_PX = 0.5
MIN_MASK_AGREEMENT = 0.99  # recall and precision of each ink mask at 2-px tolerance
MAX_ELEMENT_SHIFT_PX = 1.5
MIN_ELEMENT_NCC = 0.85  # lowest: the italic d (the published font version differs slightly)


# -----------------------------------------------------------------------------------------
# image loading
# -----------------------------------------------------------------------------------------
def extract_published_figure(pdf: Path, page: int = PAGE) -> np.ndarray:
    exe = shutil.which("pdfimages")
    if exe is None:
        raise SystemExit("poppler 'pdfimages' is required to extract the published figure")
    pdf = Path(pdf).resolve()
    if not pdf.is_file():
        raise SystemExit(f"preprint PDF not found: {pdf}")
    with tempfile.TemporaryDirectory() as tmp:
        try:
            subprocess.run([exe, "-f", str(page), "-l", str(page), "-png", str(pdf), f"{tmp}/s"],
                           check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as err:
            raise SystemExit(f"pdfimages failed on {pdf}: {err.stderr.strip()}") from err
        strips = [np.asarray(Image.open(f).convert("RGB")) for f in sorted(glob.glob(f"{tmp}/s-*.png"))]
    heights = [s.shape[0] for s in strips]
    if heights != [h for h, _ in PUBLISHED_STRIPS]:
        raise SystemExit(f"unexpected image strips on page {page}: heights {heights}")
    width = min(s.shape[1] for s in strips)
    if max(s.shape[1] for s in strips) - width > 1:
        raise SystemExit("image strips differ in width by more than one pixel")
    for s in strips:  # the grey strips carry one extra (blank) column on the right
        if s.shape[1] > width and s[:, width:].min() < 250:
            raise SystemExit("cropped strip column is not blank")
    return np.concatenate([s[:, :width] for s in strips], axis=0)


def publication_pipeline(rgb: np.ndarray) -> np.ndarray:
    """Apply the preprint's storage of Fig. 3 (grey / JPEG q90 4:2:0 strips) to an image, so
    that the replica and the published raster carry the same compression artefacts."""
    out, y = [], 0
    for height, kind in PUBLISHED_STRIPS:
        strip = Image.fromarray(np.ascontiguousarray(rgb[y:y + height]))
        if kind == "L":
            out.append(np.asarray(strip.convert("L").convert("RGB")))
        else:
            buf = io.BytesIO()
            strip.save(buf, format="JPEG", quality=90, subsampling=2)
            buf.seek(0)
            out.append(np.asarray(Image.open(buf).convert("RGB")))
        y += height
    return np.concatenate(out, axis=0)


def load_rgb(path: Path) -> np.ndarray:
    im = Image.open(path)
    if im.mode in ("RGBA", "LA"):
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    return np.asarray(im.convert("RGB"))


# -----------------------------------------------------------------------------------------
# measurements (identical for both images)
# -----------------------------------------------------------------------------------------
def masks(rgb: np.ndarray) -> dict:
    r, g, b = (rgb[..., i].astype(float) for i in range(3))
    frame = np.zeros(rgb.shape[:2], bool)  # 30-px band around the legend frame stroke
    x0, y0, x1, y1 = LEGEND_FRAME
    frame[y0 - 15:y1 + 16, x0 - 15:x1 + 16] = True
    frame[y0 + 15:y1 - 14, x0 + 15:x1 - 14] = False
    return dict(
        blue=(b > 150) & (r < 120) & (g < 120),
        red=(r > 150) & (g < 120) & (b < 120),
        black=np.maximum(np.maximum(r, g), b) < 128,
        grey=frame & (np.abs(r - g) < 12) & (np.abs(g - b) < 12) & (r > 180) & (r < 235),
        any_ink=np.minimum(np.minimum(r, g), b) < 200,  # catch-all: any non-white pixel
    )


def _runs(idx: np.ndarray) -> list[tuple[int, int]]:
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > 1)
    starts = np.r_[idx[0], idx[cut + 1]]
    ends = np.r_[idx[cut], idx[-1]]
    return list(zip(starts.tolist(), ends.tolist()))


def calibrate(rgb: np.ndarray) -> dict:
    """Locate the two y-spines, the x-axis, the tick marks and the axis limits."""
    ink = (255.0 - rgb.mean(axis=2)) / 255.0
    black = masks(rgb)["black"]
    h, w = black.shape
    spines = [(a + b) / 2 for a, b in _runs(np.flatnonzero(black.sum(axis=0) > 0.45 * h)) if b - a < 30]
    if len(spines) < 2:
        raise RuntimeError("could not locate the two y-axes")
    y0 = None
    for a, b in _runs(np.flatnonzero(black.sum(axis=1) > 0.35 * w)):  # the x-axis line
        if 5 < b - a < 25:
            rows = np.arange(a - 3, b + 4)
            prof = ink[a - 3:b + 4, int(spines[0]) + 40:int(spines[0]) + 400].mean(axis=1)
            y0 = float((prof * rows).sum() / prof.sum())
            break
    if y0 is None:
        raise RuntimeError("could not locate the x-axis")
    out = {"x_axis_y": y0, "panels": []}
    for sx in spines[:2]:
        row = black[int(y0) + 30]  # x ticks: short vertical marks below the axis
        xt = [(a + b) / 2 for a, b in _runs(np.flatnonzero(row)) if b - a < 25 and a > sx - 10]
        xt = [x for x in xt if x < sx + 2300][:5]
        col = black[: int(y0) + 10, int(sx) - 25]  # y ticks: short marks left of the spine
        yt = [(a + b) / 2 for a, b in _runs(np.flatnonzero(col)) if b - a < 25]
        if len(xt) < 3 or len(yt) < 3:
            raise RuntimeError("could not locate the tick marks")
        p = np.polyfit(np.arange(len(xt)) * 20.0, xt, 1)
        yv = np.arange(len(yt))[::-1].astype(float)
        q = np.polyfit(yv[:-1], yt[:-1], 1)  # exclude the tick at 0 (overlaps the axis)
        top = float(np.flatnonzero(black[:, int(sx)])[0]) - 0.5
        width = float(ink[int(y0) - 600, int(sx) - 20:int(sx) + 21].sum())
        out["panels"].append(dict(
            spine_x=float(sx), x0_px=float(p[1]), px_per_mm=float(p[0]),
            y_px_per_unit=float(-q[0]), y_ticks=len(yt), spine_width_px=width,
            ylim_top=float((y0 - (top + width / 2)) / -q[0]),  # projecting cap removed
        ))
    return out


def digitize(rgb: np.ndarray, colour: str, cal: dict, panel: int, depths_mm: np.ndarray,
             exclude_box=None) -> np.ndarray:
    """Vertical-centroid digitization of a thick anti-aliased curve (data units)."""
    r, g, b = (rgb[..., i].astype(float) for i in range(3))
    if colour == "blue":
        wgt = np.clip((b - np.maximum(r, g)) / 255.0, 0, 1)
    else:
        wgt = np.clip((r - np.maximum(g, b)) / 255.0, 0, 1)
    if exclude_box is not None:
        x0, y0_, x1, y1 = exclude_box
        wgt[y0_:y1, x0:x1] = 0
    n_rows, n_cols = wgt.shape
    pc = cal["panels"][panel]
    out = np.full(depths_mm.size, np.nan)
    for k, d in enumerate(depths_mm):
        x = pc["x0_px"] + d * pc["px_per_mm"]
        i = int(np.floor(x))
        if not 0 <= i < n_cols - 1:
            continue
        f = x - i
        prof = (1 - f) * wgt[:, i] + f * wgt[:, i + 1]
        idx = np.flatnonzero(prof > 0.5)
        if idx.size < 10:
            continue
        a, bb = max(_runs(idx), key=lambda t: t[1] - t[0])
        if bb - a < 18:  # dash gaps / curve ends
            continue
        lo, hi = max(a - 3, 0), min(bb + 4, n_rows)
        seg = prof[lo:hi]
        out[k] = (cal["x_axis_y"] - (seg * np.arange(lo, hi)).sum() / seg.sum()) / pc["y_px_per_unit"]
    return out


def curve_statistics(published: dict, replica: dict, cal: dict, depths_mm: np.ndarray) -> dict:
    """|published - replica| of digitized curves: vertical and curve-normal distances [px]."""
    out = {}
    for key in published:
        pc = cal["panels"][0 if key.startswith("A") else 1]
        a, b = published[key], replica[key]
        ok = np.isfinite(a) & np.isfinite(b)
        n_pub = int(np.isfinite(a).sum())
        if ok.sum() < 3:
            out[key] = dict(n=int(ok.sum()), n_published=n_pub, median_abs_diff_px=np.nan,
                            max_abs_diff_px=np.nan, median_normal_diff_px=np.nan,
                            max_normal_diff_px=np.nan, median_rel_diff=np.nan)
            continue
        dy = np.abs(a[ok] - b[ok]) * pc["y_px_per_unit"]
        slope = np.gradient(a[ok] * pc["y_px_per_unit"], depths_mm[ok] * pc["px_per_mm"])
        dn = dy / np.sqrt(1.0 + slope**2)
        out[key] = dict(n=int(ok.sum()), n_published=n_pub,
                        median_abs_diff_px=float(np.median(dy)), max_abs_diff_px=float(dy.max()),
                        median_normal_diff_px=float(np.median(dn)), max_normal_diff_px=float(dn.max()),
                        median_rel_diff=float(np.median(np.abs(a[ok] - b[ok]) / np.abs(a[ok]))))
    return out


def dotted_line_x(rgb: np.ndarray, cal: dict, panel: int) -> float:
    """Horizontal position [mm] of the dotted equal-SNR line of a panel, NaN if absent.

    The line is a column of ~25-px black dots repeating every ~66 px: it is located as the
    column band with the most black ink outside the coloured curves, anywhere in the panel,
    and refined by the ink centroid of the dots (rows fully covered across the search window,
    e.g. crossing arrow shafts, are ignored)."""
    pc = cal["panels"][panel]
    r0, r1 = 250, int(cal["x_axis_y"]) - 20
    x_lo = int(pc["spine_x"]) + 30
    x_hi = int(pc["x0_px"] + 95.0 * pc["px_per_mm"]) - 10
    sub = rgb[r0:r1, x_lo:x_hi].astype(float)
    colour = np.abs(sub[..., 0] - sub[..., 2]) > 60
    black = (sub.max(axis=2) < 128) & ~colour
    band = ndi.uniform_filter1d(black.sum(axis=0).astype(float), 25)
    c = int(np.argmax(band))
    if band[c] < 400:  # less ink than ~16 dots
        return float("nan")
    w0, w1 = max(c - 20, 0), min(c + 21, black.shape[1])
    win = black[:, w0:w1]
    rows = win.any(axis=1) & (win.sum(axis=1) < 0.9 * (w1 - w0))
    dots = [(a, b) for a, b in _runs(np.flatnonzero(rows)) if 15 <= b - a + 1 <= 35]
    if len(dots) < 15:
        return float("nan")
    dot_rows = np.concatenate([np.arange(a, b + 1) for a, b in dots])
    ink = (255.0 - sub[dot_rows, w0:w1].mean(axis=2)) / 255.0
    ink[colour[dot_rows, w0:w1]] = 0
    prof = ink.sum(axis=0)
    x = x_lo + w0 + (prof * np.arange(w1 - w0)).sum() / prof.sum()
    return float((x - pc["x0_px"]) / pc["px_per_mm"])


def components(rgb: np.ndarray, cal: dict, d_eq_mm: list[float]) -> list[tuple]:
    """Bounding boxes of black ink blobs (text, arrows) with axes/ticks/dotted line removed."""
    m = masks(rgb)["black"].copy()
    y0 = int(round(cal["x_axis_y"]))
    m[y0 - 11:y0 + 12, :] = False
    for pc, deq in zip(cal["panels"], d_eq_mm):
        sx = int(round(pc["spine_x"]))
        m[: y0 + 12, sx - 11:sx + 12] = False
        for k in range(5):  # x ticks
            x = int(round(pc["x0_px"] + 20 * k * pc["px_per_mm"]))
            m[y0:y0 + 70, x - 9:x + 10] = False
        for v in range(pc["y_ticks"]):  # y ticks
            y = int(round(y0 - v * pc["y_px_per_unit"]))
            m[max(y - 9, 0):y + 10, sx - 70:sx] = False
        if np.isfinite(deq):  # dotted line
            xd = int(round(pc["x0_px"] + deq * pc["px_per_mm"]))
            m[250:y0, xd - 16:xd + 17] = False
    lab, _ = ndi.label(m, structure=np.ones((3, 3)))
    boxes = []
    for i, sl in enumerate(ndi.find_objects(lab)):
        area = int((lab[sl] == i + 1).sum())
        if area >= 30:
            boxes.append((sl[0].start, sl[0].stop, sl[1].start, sl[1].stop, area))
    return boxes


def match_components(a: list, b: list, max_dist: float = 40.0) -> tuple[list, list, list]:
    """One-to-one matching of blobs (Hungarian assignment on centre distance; blobs with an
    area ratio outside 0.5-2 cannot match).  Returns (matches, unmatched a, unmatched b)."""
    if not a or not b:
        return [], list(a), list(b)
    ca = np.array([((y0 + y1) / 2, (x0 + x1) / 2, ar) for y0, y1, x0, x1, ar in a])
    cb = np.array([((y0 + y1) / 2, (x0 + x1) / 2, ar) for y0, y1, x0, x1, ar in b])
    dist = np.hypot(ca[:, None, 0] - cb[None, :, 0], ca[:, None, 1] - cb[None, :, 1])
    ratio = cb[None, :, 2] / ca[:, None, 2]
    cost = np.where((dist < max_dist) & (ratio > 0.5) & (ratio < 2.0), dist, 1e6)
    rows, cols = linear_sum_assignment(cost)
    matches, used_a, used_b = [], set(), set()
    for i, j in zip(rows, cols):
        if cost[i, j] >= 1e6:
            continue
        used_a.add(i)
        used_b.add(j)
        ya0, ya1, xa0, xa1, area_a = a[i]
        yb0, yb1, xb0, xb1, area_b = b[j]
        matches.append(dict(orig=[ya0, ya1, xa0, xa1], repl=[yb0, yb1, xb0, xb1],
                            d_left=xb0 - xa0, d_right=xb1 - xa1, d_top=yb0 - ya0, d_bottom=yb1 - ya1,
                            area_ratio=area_b / area_a))
    return (matches, [c for i, c in enumerate(a) if i not in used_a],
            [c for j, c in enumerate(b) if j not in used_b])


def iou(a: np.ndarray, b: np.ndarray, tol: int = 0) -> dict:
    if tol:
        st = np.ones((2 * tol + 1, 2 * tol + 1), bool)
        a_d, b_d = ndi.binary_dilation(a, st), ndi.binary_dilation(b, st)
        return dict(recall=float((a & b_d).sum() / max(a.sum(), 1)),
                    precision=float((b & a_d).sum() / max(b.sum(), 1)))
    inter, union = (a & b).sum(), (a | b).sum()
    return dict(iou=float(inter / max(union, 1)), dice=float(2 * inter / max(a.sum() + b.sum(), 1)))


def register(a: np.ndarray, b: np.ndarray, max_shift: int = 5) -> tuple[float, float, float, float]:
    """Sub-pixel shift (dx, dy) that best aligns patch b to patch a (normalised
    cross-correlation, integer search by FFT + parabolic refinement).  Returns
    (ncc at zero shift, best ncc, dx, dy); a positive dx means b must move right.
    Shifts are NaN if a patch is blank or the best match lies on the search border."""
    a = a - a.mean()
    b = b - b.mean()
    norm = np.sqrt((a**2).sum() * (b**2).sum())
    if norm == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    cc = fftconvolve(a, b[::-1, ::-1], mode="full") / norm
    cy, cx = b.shape[0] - 1, b.shape[1] - 1  # zero-shift index
    win = cc[cy - max_shift:cy + max_shift + 1, cx - max_shift:cx + max_shift + 1]
    iy, ix = np.unravel_index(np.argmax(win), win.shape)
    if iy in (0, win.shape[0] - 1) or ix in (0, win.shape[1] - 1):
        return float(cc[cy, cx]), float(win[iy, ix]), float("nan"), float("nan")

    def refine(v_minus, v0, v_plus):
        den = v_minus - 2 * v0 + v_plus
        return 0.5 * (v_minus - v_plus) / den if den < 0 else 0.0

    dy = iy - max_shift + refine(*win[iy - 1:iy + 2, ix])
    dx = ix - max_shift + refine(*win[iy, ix - 1:ix + 2])
    return float(cc[cy, cx]), float(win[iy, ix]), float(dx), float(dy)


def _curves(im: np.ndarray, cal: dict, depths: np.ndarray) -> dict:
    return {"A_OPM": digitize(im, "blue", cal, 0, depths, LEGEND_BOX),
            "A_SQUID": digitize(im, "red", cal, 0, depths, LEGEND_BOX),
            "B_OPM": digitize(im, "blue", cal, 1, depths),
            "B_SQUID": digitize(im, "red", cal, 1, depths)}


def _json_safe(value):
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


# -----------------------------------------------------------------------------------------
# evaluation
# -----------------------------------------------------------------------------------------
def evaluate(orig: np.ndarray, repl_raw: np.ndarray, unregistered_raw: np.ndarray | None = None,
             cal_pub: dict | None = None) -> dict:
    """All measurements and PASS/FAIL checks for one replica (raw 1200-dpi render)."""
    if repl_raw.shape != orig.shape:
        raise SystemExit(f"shape mismatch: published {orig.shape} vs replica {repl_raw.shape}")
    repl = publication_pipeline(repl_raw)  # compare like with like
    report = {"image_shape": list(orig.shape)}
    try:
        cal = {"published": cal_pub or calibrate(orig), "replica": calibrate(repl)}
    except RuntimeError as err:
        return {"error": str(err), "checks": {"replica axes could be calibrated": False}}
    report["calibration"] = cal

    # curves: the replica is digitized with the published calibration (same pixel grid)
    depths = np.arange(15.5, 93.0, 0.5)
    c = cal["published"]
    curves_pub = _curves(orig, c, depths)
    report["curves"] = curve_statistics(curves_pub, _curves(repl, c, depths), c, depths)
    if unregistered_raw is not None:  # physics only: no fitted sub-pixel registration
        report["curves_without_registration"] = curve_statistics(
            curves_pub, _curves(publication_pipeline(unregistered_raw), c, depths), c, depths)
    deq = {n: [dotted_line_x(im, cal[n], p) for p in (0, 1)] for n, im in (("published", orig), ("replica", repl))}
    report["d_eq_line_mm"] = deq

    mo, mr, mraw = masks(orig), masks(repl), masks(repl_raw)
    report["masks"] = {k: {**iou(mo[k], mr[k]), **iou(mo[k], mr[k], tol=2)} for k in mo}
    report["masks_replica_raw_vs_pipeline"] = {k: iou(mraw[k], mr[k])["iou"] for k in mo}

    comp_o = components(orig, cal["published"], deq["published"])
    comp_r = components(repl, cal["replica"], deq["replica"])
    matches, miss_o, miss_r = match_components(comp_o, comp_r)
    offs = np.array([[m["d_left"], m["d_right"], m["d_top"], m["d_bottom"]] for m in matches]).reshape(-1, 4)
    report["components"] = dict(
        n_published=len(comp_o), n_replica=len(comp_r), n_matched=len(matches),
        median_abs_offset_px=np.median(np.abs(offs), axis=0).tolist() if len(offs) else None,
        max_abs_offset_px=np.max(np.abs(offs), axis=0).tolist() if len(offs) else None,
        worst=sorted(matches, key=lambda m: -max(abs(m["d_left"]), abs(m["d_right"]),
                                                  abs(m["d_top"]), abs(m["d_bottom"])))[:8],
        unmatched_published=miss_o[:20], unmatched_replica=miss_r[:20],
    )

    ink_o = (255.0 - orig.mean(axis=2)) / 255.0
    ink_r = (255.0 - repl.mean(axis=2)) / 255.0
    elements = {}
    for name, (y0, y1, x0, x1) in ELEMENTS.items():
        n0, nb, dx, dy = register(ink_o[y0:y1, x0:x1], ink_r[y0:y1, x0:x1])
        elements[name] = dict(ncc=n0, ncc_best=nb, shift_px=[dx, dy])
    report["elements"] = elements

    gray_o, gray_r = orig.mean(axis=2), repl.mean(axis=2)
    report["pixels"] = dict(
        mean_abs_diff=float(np.mean(np.abs(gray_o - gray_r))),
        frac_pixels_diff_gt_128=float(np.mean(np.abs(gray_o - gray_r) > 128)),
        corr_blurred=float(np.corrcoef(ndi.gaussian_filter(gray_o, 3).ravel(),
                                       ndi.gaussian_filter(gray_r, 3).ravel())[0, 1]),
    )

    pp, pr = cal["published"]["panels"], cal["replica"]["panels"]
    shifts = np.array([np.hypot(*e["shift_px"]) for e in elements.values()])
    nccs = np.array([e["ncc"] for e in elements.values()])
    cs = report["curves"].values()
    deq_px = [abs(a - b) * pp[k]["px_per_mm"] for k, (a, b) in enumerate(zip(deq["published"], deq["replica"]))]
    report["checks"] = {
        "axes frame identical (x0 and px/mm within 0.1 px)": all(
            abs(pp[k]["x0_px"] - pr[k]["x0_px"]) < 0.1 and abs(pp[k]["px_per_mm"] - pr[k]["px_per_mm"]) * 95 < 0.1
            for k in range(2)),
        "y-axis limits within 0.1 %": all(abs(pp[k]["ylim_top"] / pr[k]["ylim_top"] - 1) < 1e-3 for k in range(2)),
        f"curves digitized and median deviation < {MAX_CURVE_MEDIAN_PX} px": all(
            v["n"] >= 0.9 * v["n_published"] and v["median_abs_diff_px"] < MAX_CURVE_MEDIAN_PX for v in cs),
        f"equal-SNR line found and within {MAX_DEQ_LINE_PX} px": all(
            np.isfinite(d) and d < MAX_DEQ_LINE_PX for d in deq_px),
        f"each ink colour within 2 px (recall and precision > {MIN_MASK_AGREEMENT})": all(
            v["recall"] > MIN_MASK_AGREEMENT and v["precision"] > MIN_MASK_AGREEMENT for v in report["masks"].values()),
        "text/arrow blobs matched one-to-one (none missing, none extra)": not miss_o and not miss_r,
        f"elements registered within {MAX_ELEMENT_SHIFT_PX} px with NCC > {MIN_ELEMENT_NCC}": bool(
            np.all(np.isfinite(shifts)) and np.all(shifts < MAX_ELEMENT_SHIFT_PX) and np.all(nccs > MIN_ELEMENT_NCC)),
    }
    report["_images"] = (orig, repl, ink_o, ink_r)
    return report


def print_summary(report: dict):
    if "error" in report:
        print("Evaluation error:", report["error"])
        return
    cal = report["calibration"]
    pp, pr = cal["published"]["panels"], cal["replica"]["panels"]
    print("Axis calibration (published | replica):")
    for k in range(2):
        print(f"  panel {'AB'[k]}: x0 {pp[k]['x0_px']:.2f} | {pr[k]['x0_px']:.2f} px, "
              f"{pp[k]['px_per_mm']:.4f} | {pr[k]['px_per_mm']:.4f} px/mm, "
              f"ylim top {pp[k]['ylim_top']:.4f} | {pr[k]['ylim_top']:.4f}")
    print("Curves |published - replica| [px], median / max (vertical), median (curve-normal):")
    for k, v in report["curves"].items():
        print(f"  {k:8s} n={v['n']:3d}  {v['median_abs_diff_px']:.2f} / {v['max_abs_diff_px']:.2f}   "
              f"{v['median_normal_diff_px']:.2f}")
    if "curves_without_registration" in report:
        print("  without the fitted raster registration (physics only), median / max:",
              {k: (round(v["median_abs_diff_px"], 2), round(v["max_abs_diff_px"], 2))
               for k, v in report["curves_without_registration"].items()})
    deq = report["d_eq_line_mm"]
    print("Dotted d_eq line [mm]: published", np.round(deq["published"], 3), "replica", np.round(deq["replica"], 3))
    print("Ink masks IoU:", {k: round(v["iou"], 3) for k, v in report["masks"].items()},
          " within 2 px (recall/precision):",
          {k: (round(v["recall"], 3), round(v["precision"], 3)) for k, v in report["masks"].items()})
    c = report["components"]
    print(f"Text/arrow blobs: {c['n_matched']} matched of {c['n_published']} published / {c['n_replica']} replica; "
          f"median |offset| (l,r,t,b) {c['median_abs_offset_px']} px, max {c['max_abs_offset_px']} px")
    el = report["elements"]
    shifts = {n: np.hypot(*e["shift_px"]) for n, e in el.items()}
    worst = max(shifts, key=lambda n: shifts[n] if np.isfinite(shifts[n]) else np.inf)
    lowest = min(el, key=lambda n: el[n]["ncc"])
    print(f"Elements ({len(el)}): median |shift| {np.nanmedian(list(shifts.values())):.2f} px, max "
          f"{shifts[worst]:.2f} px ({worst}); median NCC {np.median([e['ncc'] for e in el.values()]):.4f}, "
          f"min {el[lowest]['ncc']:.3f} ({lowest})")
    print("Pixels:", {k: round(v, 5) for k, v in report["pixels"].items()})
    print("Checks:")
    for k, v in report["checks"].items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")


def write_images(report: dict, outdir: Path):
    orig, repl, ink_o, ink_r = report["_images"]
    overlay = np.stack([1 - ink_r, 1 - ink_o, 1 - 0.5 * (ink_o + ink_r)], axis=2)  # published magenta, replica green
    Image.fromarray((np.clip(overlay, 0, 1) * 255).astype(np.uint8)).save(outdir / "Figure3_overlay.png")
    diff = np.abs(orig.astype(int) - repl.astype(int)).max(axis=2)
    Image.fromarray((255 - np.clip(diff, 0, 255)).astype(np.uint8)).save(outdir / "Figure3_difference.png")
    side = np.concatenate([orig, np.full((orig.shape[0], 40, 3), 255, np.uint8), repl], axis=1)
    Image.fromarray(side).resize((side.shape[1] // 3, side.shape[0] // 3),
                                 Image.Resampling.LANCZOS).save(outdir / "Figure3_side_by_side.png")


def _render(path: Path, **options) -> np.ndarray:
    import replicate_figure3
    replicate_figure3.render_png(path, **options)
    return load_rgb(path)


def self_test(orig: np.ndarray, repl_raw: np.ndarray, options: dict, base_passed: bool) -> bool:
    """Negative controls: every perturbed replica must fail at least one check (meaningful
    only if the unperturbed replica passes all checks)."""
    if not base_passed:
        print("Negative controls skipped: the unperturbed replica does not pass all checks")
        return False
    cal_pub = calibrate(orig)
    controls = {}
    img = repl_raw.copy()
    img[1500:1545, 1600:1645] = 0  # spurious 45 x 45 px blob in an empty part of panel A
    controls["spurious blob"] = img
    img = repl_raw.copy()
    r, g, b = (img[..., i].astype(int) for i in range(3))
    img[(b > 150) & (r < 120) & (g < 120)] = 255  # OPM curves removed
    controls["missing OPM curves"] = img
    controls["whole figure shifted by 2 px"] = np.roll(repl_raw, 2, axis=1)
    with tempfile.TemporaryDirectory() as tmp:
        base = dict(options, raster_registration=True)
        controls["d_eq marker at the root of Eq. (3)"] = _render(Path(tmp) / "exact.png", **dict(base, deq_marker="exact"))
        controls["eta = 3.1"] = _render(Path(tmp) / "eta.png", **dict(base, eta=3.1))
    ok = True
    print("Negative controls (each must FAIL):")
    for name, img in controls.items():
        checks = evaluate(orig, img, cal_pub=cal_pub)["checks"]
        failed = [k for k, v in checks.items() if not v]
        detected = bool(failed)
        ok &= detected
        print(f"  [{'OK' if detected else 'NOT DETECTED'}] {name}: failed -> {failed}")
    return ok


# -----------------------------------------------------------------------------------------
def main(argv=None) -> bool:
    ap = argparse.ArgumentParser(description="Compare the replicated Figure 3 with the preprint raster.")
    ap.add_argument("--pdf", default=str(HERE / "2026.08.17.744953.full.pdf"))
    ap.add_argument("--replica", default=str(HERE / "figure3_output" / "Figure3_replicated.png"))
    ap.add_argument("--outdir", default=str(HERE / "figure3_output"))
    ap.add_argument("--skip-unregistered", action="store_true",
                    help="skip the comparison with a replica rendered without raster registration")
    ap.add_argument("--self-test", action="store_true", help="run the negative controls as well")
    ap.add_argument("--artwork-fonts", action="store_true",
                    help="render replicas with the artwork fonts when the replica's summary does not say (see replicate_figure3.py)")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    replica = Path(args.replica)

    orig = extract_published_figure(Path(args.pdf))
    Image.fromarray(orig).save(outdir / "Figure3_published_extracted.png")

    # options of the replica under test (written by replicate_figure3.py next to the PNG)
    summary_file = replica.parent / "figure3_summary.json"
    options = dict(eta=3.0, deq_marker="published", dtheta_deg=0.05, artwork_fonts=args.artwork_fonts)
    if summary_file.is_file():
        s = json.loads(summary_file.read_text())
        png = s.get("outputs", {}).get("png")
        if png and Path(png).resolve() == replica.resolve():  # summary written with this PNG
            options = dict(eta=s.get("eta", 3.0), deq_marker=s.get("deq_marker", "published"),
                           dtheta_deg=s.get("dtheta_deg", 0.05),
                           artwork_fonts=s.get("artwork_fonts", "Myriad" in s.get("font", "")))
    if not replica.is_file():  # end-to-end: render the replica first
        replica.parent.mkdir(parents=True, exist_ok=True)
        _render(replica, **options)
    repl_raw = load_rgb(replica)

    unregistered = None
    if not args.skip_unregistered:
        with tempfile.TemporaryDirectory() as tmp:
            unregistered = _render(Path(tmp) / "plain.png", raster_registration=False, **options)

    report = evaluate(orig, repl_raw, unregistered)
    report["replica_options"] = options
    ok = all(report["checks"].values())
    if "_images" in report:
        write_images(report, outdir)
    print_summary(report)
    if args.self_test:
        report["self_test_passed"] = self_test(orig, repl_raw, options, base_passed=ok)
        print("Self-test:", "PASS" if report["self_test_passed"] else "FAIL")
        ok &= report["self_test_passed"]
    with open(outdir / "verification_report.json", "w") as fh:
        json.dump(_json_safe({k: v for k, v in report.items() if k != "_images"}), fh, indent=2, allow_nan=False)
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
