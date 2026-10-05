#!/usr/bin/env python3
"""Noise-model covariance validation against the measured Neuromag noise (revision, referee round 1).

Requests
* Codex, major issue 1: "Extend the checks beyond median channel RMS to held-out channelwise
  variances, spatial correlations, covariance eigenstructure and spectral behavior where the
  supplied measurement framework permits. Distinguish quantities fitted to data from independent
  predictions."
* Fable, major issue 1: the cortical background "is calibrated on one scalar (median Neuromag
  gradiometer variance), under-predicts the measured magnetometer brain noise by 27 % (0.73x),
  spans 0.13-0.78 of the measured per-channel variance on the magnetometers".

Measured noise, as G2 builds it (``g2.measured_noise``): the MNE sample recording's 320 pre-stimulus
windows (-200 to 0 ms) and its empty-room recording; each channel's record mean removed, filtered
as one continuous record (configs/g2_adult.toml band, zero-phase Butterworth), no SSP, MEG 2443
(bad) left out. Here the full channel x channel covariances, not only their diagonals (the
diagonals reproduce G2's variances; checked). Measured brain covariance = baseline - empty room.

Model, as G2 builds it for Neuromag (``g2.array_noise``): brochure white sensor noise + cortical
background (1,755 independent, area-weighted cortical sources; ONE scale fitted so that the median
good-gradiometer variance equals the measured one) + room field (8-term external expansion, its
coefficient covariance fitted to the empty-room recording). Variants (the G2 sensitivity values):
background calibrated on the magnetometers; correlated background (configs/g2_adult.toml lengths).

Analyses, per channel type (102 magnetometers, 203 good gradiometers)
1. channelwise variances: ratio distribution and correlation of the spatial pattern; held out: the
   scale fitted on half the gradiometer sites predicts the other half (random and spatial splits)
   and the magnetometers;
2. spatial correlation vs inter-sensor distance;
3. eigenvalue spectra (normalised) and principal angles between the dominant subspaces;
4. sub-bands of the analysis band (declared edges): band-averaged spectra, the magnetometer
   prediction and the channel pattern per band (model rescaled per band on the gradiometers, as in
   scripts/g2_band_sensitivity.py);
5. Neuromag known-topography detectability of the 7,661 G2 targets with the measured covariance
   (task baseline; empty room + model brain) vs the model covariance: median ratio, parcel CI;
6. the implied OPM/Neuromag ratio under stated assumptions (S1-S5 in the JSON).
Diagnostics of the measured noise that the model omits: the heart's field (heartbeats detected from
the magnetometers; the cardiac-locked average removed from every baseline sample at its cardiac
phase), how concentrated the measured excess over the model is and how much of it lies in the
8-dimensional far-field subspace, and any time-locked residue of the preceding stimulus.
References for every statistic: split halves of the measured data; stationary surrogates of the
recordings (multivariate phase randomisation: the measured cross-spectra are kept, non-stationarity
is removed), which give the finite-sample bias of detectability computed with an estimated
covariance; and the model "as it would be measured": the surrogates recoloured to the model
covariance, i.e. the null distribution of each statistic if the model were exact.

Copied rather than refactored (shared modules unchanged): the reading and filtering of
g2.measured_noise (``Recording``, ``baseline_windows``; it returns only variances) and
noisemodel.ArrayNoise.projector on a channel subset (``projector``); both checked against the
originals at run time or in the unit checks.

Outputs (results/g2_covariance_validation/ unless --out): covariance_validation.json,
covariance_validation_targets.csv, Figure_covariance_structure.png, Figure_covariance_bands.png,
Figure_covariance_detectability.png.
"""
from __future__ import annotations

import argparse
import csv
import inspect
import resource
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402
from scipy import stats  # noqa: E402

import g2_adult_comparison as G  # noqa: E402
import report_style as RS  # noqa: E402
from opmsquid import background, environment, g2, io, metrics, neuromag, noise, paths  # noqa: E402

OUT = paths.RESULTS / "g2_covariance_validation"
STATUS = "REVISION (noise-model covariance validation against the measured Neuromag noise, MNE sample subject)"

# Declared analysis choices (not model parameters; all written to the JSON). Everything else comes
# from configs/g2_adult.toml or from G2's code (window and trim: g2.measured_noise's defaults).
SUB_BANDS = ((1.0, 4.0), (4.0, 8.0), (8.0, 13.0), (13.0, 20.0), (20.0, 30.0), (30.0, 40.0))
SUB_BAND_NOTE = "conventional delta, theta, alpha, low beta, high beta and low gamma edges inside the 1-40 Hz analysis band"
DIST_EDGES_MM = np.arange(25.0, 276.0, 25.0)  # pair-distance bins; nearest Neuromag sites are 33-39 mm apart
K_SUBSPACE = (1, 2, 3, 5, 10, 20, 40)
N_EIG_STORED = 60
WINDOW = inspect.signature(g2.measured_noise).parameters["window"].default
TRIM_S = inspect.signature(g2.measured_noise).parameters["trim_s"].default

REQUESTS = {
    "codex_major_1": "Extend the checks beyond median channel RMS to held-out channelwise variances, spatial correlations, "
                     "covariance eigenstructure and spectral behavior where the supplied measurement framework permits. "
                     "Distinguish quantities fitted to data from independent predictions.",
    "fable_major_1": "That model is calibrated on one scalar (median Neuromag gradiometer variance), under-predicts the measured "
                     "magnetometer brain noise by 27 % (0.73x), spans 0.13-0.78 of the measured per-channel variance on the "
                     "magnetometers",
}
FITTED_VS_PREDICTED = [
    dict(quantity="cortical background scale (one scalar; one per sub-band in the band analysis)", role="fitted",
         data="median over the good gradiometers of the measured brain variance (task baseline - empty room) in the band"),
    dict(quantity="room-field coefficient covariance (8 x 8)", role="fitted",
         data="empty-room recording; acts only inside the 8-dimensional external subspace of the model's empty-room and total covariance"),
    dict(quantity="sensor white noise (3.5 fT/sqrt(Hz), 3.6 fT/(cm sqrt(Hz)))", role="not fitted (brochure values HW-mag-noise, HW-grad-noise)",
         data="compared with the empty-room variances outside the room-field subspace"),
    dict(quantity="median good-gradiometer brain variance", role="fitted (equal by construction)", data="-"),
    dict(quantity="gradiometer brain variances at held-out sites (scale fitted on the other half)", role="independent prediction", data="-"),
    dict(quantity="gradiometer channel pattern (relative variances)", role="independent prediction (only the overall scale is fitted)", data="-"),
    dict(quantity="magnetometer brain variances (scale and pattern)", role="independent prediction", data="-"),
    dict(quantity="spatial correlations vs inter-sensor distance", role="independent prediction", data="-"),
    dict(quantity="eigenvalue spectra and dominant subspaces", role="independent prediction", data="-"),
    dict(quantity="magnetometer variance and channel patterns per sub-band", role="independent prediction (scale refitted per band on the gradiometers)",
         data="-"),
    dict(quantity="Neuromag detectability with the measured covariance", role="independent prediction", data="-"),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r * 1024 / 1e9


# ----------------------------------------------------------------------------------------------
# measured data (the reading and filtering of g2.measured_noise, keeping the full covariances)
class Recording:
    """One recording of the sample dataset read as ``g2.measured_noise`` reads it: the given SQUID
    channels in G2's order, each channel's mean over the whole record removed, filtered as one
    continuous record. Copied from g2.measured_noise (which returns only the per-channel variances)."""

    def __init__(self, fname: str, names: list[str]):
        self.raw = mne.io.read_raw_fif(paths.require(paths.SAMPLE_MEG / fname, fname), preload=False, verbose=False)
        self.sfreq = float(self.raw.info["sfreq"])
        data = self.raw.get_data(picks=[self.raw.ch_names.index(n) for n in names])
        self.data = data - data.mean(axis=1, keepdims=True)
        self.n_times = self.data.shape[1]

    def filtered(self, filt: noise.AnalysisFilter, chunk: int = 64) -> np.ndarray:
        out = np.empty_like(self.data)
        for i in range(0, len(out), chunk):  # row-wise filter: identical to filtering all rows at once
            out[i:i + chunk] = filt.apply(self.data[i:i + chunk])
        return out


def baseline_windows(rec: Recording) -> tuple[list[np.ndarray], np.ndarray]:
    """Sample indices of the pre-stimulus windows used by g2.measured_noise (same rule) and the event
    code of each window."""
    sf, n_trim = rec.sfreq, int(TRIM_S * rec.sfreq)
    events = mne.find_events(rec.raw, stim_channel="STI 014", verbose=False)
    ev = events[:, 0] - rec.raw.first_samp
    a, b = int(round(WINDOW[0] * sf)), int(round(WINDOW[1] * sf))
    keep = [i for i, e in enumerate(ev) if e + a > n_trim and e + b < rec.n_times - n_trim]
    return [np.arange(ev[i] + a, ev[i] + b) for i in keep], events[keep, 2]


def time_locked_removed_cov(x: np.ndarray, wins: list[np.ndarray], codes: np.ndarray) -> np.ndarray:
    """Covariance of the windows after removing, per event type, their time-locked average (any residual
    evoked response from the preceding stimulus), normalised by the remaining degrees of freedom."""
    xw = np.stack([x[:, w] for w in wins], axis=1)  # (channels, windows, samples)
    acc, dof = np.zeros((x.shape[0], x.shape[0])), 0
    for c in np.unique(codes):
        xc = xw[:, codes == c, :]
        xc = (xc - xc.mean(axis=1, keepdims=True)).reshape(x.shape[0], -1)
        acc += xc @ xc.T
        dof += (int(np.sum(codes == c)) - 1) * xw.shape[2]
    return acc / dof


def cardiac_field(rec: Recording, x: np.ndarray) -> dict:
    """The heart's field in the task recording: heartbeats detected from the magnetometers (MNE's
    find_ecg_events with its defaults; the sample data have no ECG channel) and the cardiac-locked
    average of the band-passed data over one median R-R cycle centred on the R peak. Its covariance
    (the average waveform's second moment over the cycle) is what a heartbeat at random phase relative
    to the baseline windows adds to their covariance. ``x``: the recording band-passed (``Recording.filtered``)."""
    ev = mne.preprocessing.find_ecg_events(rec.raw, verbose=False)[0]
    r_all = np.sort(ev[:, 0] - rec.raw.first_samp)
    rr = int(np.median(np.diff(r_all)))
    hw, n_trim = rr // 2, int(TRIM_S * rec.sfreq)
    r = r_all[(r_all > hw + n_trim) & (r_all < rec.n_times - hw - n_trim)]
    acc = np.zeros((x.shape[0], 2 * hw))
    for i in r:
        acc += x[:, i - hw:i + hw]
    avg = acc / len(r)
    return dict(average=avg, cov=cov(avg), n_beats=int(len(r)), rr_s=float(rr / rec.sfreq), r_all=r_all, hw=hw)


def heart_template(card: dict, idx: np.ndarray) -> np.ndarray:
    """The cardiac-locked average at the samples ``idx``, placed by each sample's offset from its nearest
    R peak (zero beyond half a median R-R interval): subtracting it removes the average heartbeat field
    from those samples at their actual cardiac phase."""
    r_all, hw, avg = card["r_all"], card["hw"], card["average"]
    pos = np.searchsorted(r_all, idx)
    left, right = r_all[np.clip(pos - 1, 0, len(r_all) - 1)], r_all[np.clip(pos, 0, len(r_all) - 1)]
    off = idx - np.where(np.abs(idx - left) <= np.abs(right - idx), left, right)
    ok = (off >= -hw) & (off < hw)
    out = np.zeros((avg.shape[0], len(idx)))
    out[:, ok] = avg[:, hw + off[ok]]
    return out


def cov(x: np.ndarray) -> np.ndarray:
    """Second-moment matrix of band-passed data (n_channels, n_samples), as g2.measured_noise's variances."""
    return x @ x.T / x.shape[1]


class Surrogates:
    """Stationary surrogates of a filtered multichannel record: every channel's Fourier phases rotated
    by the same random phase per frequency (multivariate phase randomisation). This keeps every
    auto- and cross-spectrum, so the full-record covariance is exactly that of the data and the
    temporal correlation is the measured one; non-stationarity and non-Gaussianity are removed."""

    def __init__(self, x: np.ndarray):
        self.n = x.shape[1]
        self.cov_full = cov(x)
        self.spec = np.fft.rfft(x, axis=1)

    def sample(self, idx: np.ndarray, rng: np.random.Generator, chunk: int = 64) -> np.ndarray:
        """Covariance of one surrogate over the samples ``idx`` (as for the data)."""
        ph = np.exp(2j * np.pi * rng.random(self.spec.shape[1]))
        ph[0] = 1.0
        if self.n % 2 == 0:
            ph[-1] = 1.0
        x = np.empty((self.spec.shape[0], len(idx)))
        for i in range(0, len(x), chunk):
            x[i:i + chunk] = np.fft.irfft(self.spec[i:i + chunk] * ph[None, :], n=self.n, axis=1)[:, idx]
        return cov(x)


def recolour(c_target: np.ndarray, c_source: np.ndarray) -> np.ndarray:
    """A with A c_source A^T = c_target (symmetric square roots, on unit-diagonal matrices): data with
    covariance c_source become data with covariance c_target, keeping their temporal structure."""
    ds, dt = 1.0 / np.sqrt(np.diag(c_source)), np.sqrt(np.diag(c_target))
    ws, us = np.linalg.eigh(c_source * np.outer(ds, ds))
    wt, ut = np.linalg.eigh(c_target / np.outer(dt, dt))
    inv_sqrt = (us / np.sqrt(ws)) @ us.T
    sqrt_t = (ut * np.sqrt(np.clip(wt, 0.0, None))) @ ut.T
    return (dt[:, None] * sqrt_t) @ (inv_sqrt * ds[None, :])


# ----------------------------------------------------------------------------------------------
# structure statistics
def eig_desc(c):
    w, u = np.linalg.eigh(0.5 * (c + c.T))
    o = np.argsort(w)[::-1]
    return w[o], u[:, o]


def spectrum_summary(w: np.ndarray) -> dict:
    frac = w / w.sum()
    cum = np.cumsum(frac)
    pos = w[w > 0]
    return dict(top1=float(frac[0]), top5=float(cum[4]), top10=float(cum[9]),
                n_for_50pct=int(np.searchsorted(cum, 0.5) + 1), n_for_90pct=int(np.searchsorted(cum, 0.9) + 1),
                n_for_95pct=int(np.searchsorted(cum, 0.95) + 1),
                participation_ratio=float(pos.sum() ** 2 / np.sum(pos**2)),
                negative_eigenvalue_share=float(-w[w < 0].sum() / pos.sum()))


def subspace_overlap(ua, ub, k):
    """Mean cos^2 of the principal angles between the top-k subspaces (1 = identical, k/n = random),
    the largest principal angle [deg]."""
    s = np.clip(np.linalg.svd(ua[:, :k].T @ ub[:, :k], compute_uv=False), 0.0, 1.0)
    return float(np.mean(s**2)), float(np.degrees(np.arccos(s.min())))


def pct(x, q=(5, 25, 50, 75, 95)):
    return {f"p{p}": float(v) for p, v in zip(q, np.percentile(x, q))}


def pair_stats(a, b, idx, dist_mm, full=True):
    """Agreement of covariance ``b`` (reference, e.g. the model) with ``a`` (e.g. measured) on the
    channels ``idx``: variances, spatial correlations vs distance, eigenstructure. Ratios are b/a."""
    a_, b_ = a[np.ix_(idx, idx)], b[np.ix_(idx, idx)]
    va, vb = np.diag(a_), np.diag(b_)
    pos = (va > 0) & (vb > 0)
    out = dict(n_channels=int(len(idx)), n_nonpositive=int(np.sum(va <= 0)),
               amp_ratio_of_medians=float(np.sqrt(np.median(vb) / np.median(va))),
               var_ratio=pct(vb[pos] / va[pos]),
               log_var_pearson=float(np.corrcoef(np.log(va[pos]), np.log(vb[pos]))[0, 1]),
               log_var_spearman=float(stats.spearmanr(va[pos], vb[pos])[0]))
    ia = np.flatnonzero(pos)
    ra = a_[np.ix_(ia, ia)] / np.sqrt(np.outer(va[ia], va[ia]))
    rb = b_[np.ix_(ia, ia)] / np.sqrt(np.outer(vb[ia], vb[ia]))
    iu = np.triu_indices(len(ia), 1)
    xa, xb, d = ra[iu], rb[iu], dist_mm[np.ix_(idx[ia], idx[ia])][iu]
    out["corr_pattern_pearson"] = float(np.corrcoef(xa, xb)[0, 1])
    bins = [("same site", d < 1.0)] + [(f"{lo:g}-{hi:g}", (d >= lo) & (d < hi)) for lo, hi in zip(DIST_EDGES_MM[:-1], DIST_EDGES_MM[1:])]
    curve = []
    for lab, m in bins:
        n = int(m.sum())
        row = dict(bin=lab, n_pairs=n)
        if n >= 10:
            row.update(median_r_a=float(np.median(xa[m])), median_r_b=float(np.median(xb[m])),
                       median_abs_r_a=float(np.median(np.abs(xa[m]))), median_abs_r_b=float(np.median(np.abs(xb[m]))),
                       agreement=float(np.corrcoef(xa[m], xb[m])[0, 1]))
        curve.append(row)
    out["corr_vs_distance"] = curve
    if full:
        wa, ua = eig_desc(a_)
        wb, ub = eig_desc(b_)
        out["eig_a"], out["eig_b"] = spectrum_summary(wa), spectrum_summary(wb)
        out["eig_a_normalised"] = (wa / wa.sum())[:N_EIG_STORED].tolist()
        out["eig_b_normalised"] = (wb / wb.sum())[:N_EIG_STORED].tolist()
        ov = {}
        for k in K_SUBSPACE:
            if k < len(idx):
                o, ang = subspace_overlap(ua, ub, k)
                captured = float(np.trace(ub[:, :k].T @ a_ @ ub[:, :k]) / wa[:k].sum())
                ov[str(k)] = dict(overlap=o, largest_angle_deg=ang, random_overlap=k / len(idx), captured_variance=captured)
        out["subspace"] = ov
    return out


def eig_stats(a, b, idx):
    """Eigenstructure part of ``pair_stats`` only (combined channel set)."""
    a_, b_ = a[np.ix_(idx, idx)], b[np.ix_(idx, idx)]
    wa, ua = eig_desc(a_)
    wb, ub = eig_desc(b_)
    ov = {}
    for k in K_SUBSPACE:
        o, ang = subspace_overlap(ua, ub, k)
        ov[str(k)] = dict(overlap=o, largest_angle_deg=ang, random_overlap=k / len(idx),
                          captured_variance=float(np.trace(ub[:, :k].T @ a_ @ ub[:, :k]) / wa[:k].sum()))
    return dict(n_channels=int(len(idx)), eig_a=spectrum_summary(wa), eig_b=spectrum_summary(wb),
                eig_a_normalised=(wa / wa.sum())[:N_EIG_STORED].tolist(), eig_b_normalised=(wb / wb.sum())[:N_EIG_STORED].tolist(),
                subspace=ov)


def excess_stats(a, b, idx, ext_n, ks=(1, 3, 5, 10)):
    """Eigen-decomposition of a - b (e.g. measured - model, normalised units) on the channels idx: how
    concentrated the excess is (share of its positive part in the top-k components, relative to b's
    trace) and the share of its leading components inside the 8-dimensional external-field subspace
    (homogeneous field + linear gradients, ``ext_n`` = normalised basis; a random direction: 8 / n)."""
    e = (a - b)[np.ix_(idx, idx)]
    w, u = eig_desc(e)
    pos = w[w > 0]
    q = np.linalg.qr(ext_n[idx])[0]
    tr_b = float(np.trace(b[np.ix_(idx, idx)]))
    return dict(positive_part_over_model_trace=float(pos.sum() / tr_b), negative_part_over_model_trace=float(-w[w < 0].sum() / tr_b),
                top_k_share_of_positive_part={str(k): float(pos[:k].sum() / pos.sum()) for k in ks},
                top_k_over_model_trace={str(k): float(pos[:k].sum() / tr_b) for k in ks},
                leading_component_external_share={str(j + 1): float(np.sum((q.T @ u[:, j]) ** 2)) for j in range(3)},
                external_share_random=float(ext_n.shape[1] / len(idx)))


def flat_scalars(d, prefix=""):
    """{dotted key: float} of every scalar in a nested pair_stats result (lists of bins by bin label)."""
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat_scalars(v, f"{prefix}{k}."))
    elif isinstance(d, list) and d and isinstance(d[0], dict) and "bin" in d[0]:
        for row in d:
            out.update(flat_scalars({k: v for k, v in row.items() if k != "bin"}, f"{prefix}{row['bin']}."))
    elif isinstance(d, (int, float)) and not isinstance(d, bool):
        out[prefix[:-1]] = float(d)
    return out


def summarise_null(samples: list[dict]) -> dict:
    """2.5 / 50 / 97.5 percentiles of every scalar over the null (or surrogate) samples."""
    flats = [flat_scalars(s) for s in samples]
    keys = [k for k in flats[0] if all(k in f for f in flats)]
    return {k: [float(np.percentile([f[k] for f in flats], q)) for q in (2.5, 50, 97.5)] for k in keys}


# ----------------------------------------------------------------------------------------------
def held_out(meas_var, unit_var, grads, mags, site, xyz, n_splits, rng):
    """Calibrate the background scale on half the gradiometer sites (median rule of G2) and predict
    the other half and the magnetometers. ``meas_var``/``unit_var``: per channel (good channels)."""

    def evaluate(cal, ho):
        scale = np.median(meas_var[cal]) / np.median(unit_var[cal])
        pred = scale * unit_var
        pos = ho[meas_var[ho] > 0]
        lr = np.log10(pred[pos] / meas_var[pos])
        return dict(heldout_amp_ratio_of_medians=float(np.sqrt(np.median(pred[ho]) / np.median(meas_var[ho]))),
                    heldout_median_abs_log10_ratio=float(np.median(np.abs(lr))),
                    heldout_log_var_pearson=float(np.corrcoef(np.log(pred[pos]), np.log(meas_var[pos]))[0, 1]),
                    mag_amp_ratio_of_medians=float(np.sqrt(np.median(pred[mags]) / np.median(meas_var[mags]))),
                    n_calibration=int(len(cal)), n_heldout=int(len(ho)))

    sites = np.unique(site[grads])
    rand = []
    for _ in range(n_splits):
        cal_sites = rng.permutation(sites)[:len(sites) // 2]
        in_cal = np.isin(site[grads], cal_sites)
        rand.append(evaluate(grads[in_cal], grads[~in_cal]))
    out = dict(random_site_splits=dict(n_splits=n_splits, n_sites=int(len(sites)),
                                       **{k: dict(zip(("p2.5", "p50", "p97.5"), np.percentile([r[k] for r in rand], [2.5, 50, 97.5]).tolist()))
                                          for k in rand[0] if not k.startswith("n_")}))
    spatial = {}
    for ax, name, cut in ((0, "left_right", 0.0), (1, "posterior_anterior", None), (2, "inferior_superior", None)):
        c = np.median(xyz[np.unique(site[grads]), ax]) if cut is None else cut
        side = xyz[site[grads], ax] >= c
        lo_name, hi_name = name.split("_")
        spatial[f"fit_{lo_name}_predict_{hi_name}"] = evaluate(grads[~side], grads[side])
        spatial[f"fit_{hi_name}_predict_{lo_name}"] = evaluate(grads[side], grads[~side])
    out["spatial_splits"] = spatial
    out["spatial_split_rule"] = "device frame: left/right at x = 0; posterior/anterior and inferior/superior at the median site coordinate"
    out["all_gradiometers"] = evaluate(grads, grads)
    return out


# ----------------------------------------------------------------------------------------------
def detect(s, c):
    return metrics.detectability(s, c)


def ratio_summary(d_a, d_b, groups, rng, n_boot):
    """Median ratio d_a/d_b over targets with the G2 parcel-bootstrap 95 % CI."""
    c = G.compare(d_a, d_b, rng, n_boot, groups=groups)
    return dict(ratio=float(2 ** c["median_log2"]), ci95=[float(2 ** x) for x in c["ci95"]], share_above_1=c["share_opm_better"],
                share_parcels_above_1=c.get("share_parcels_opm_better"), n=c["n"], ci_method=c["ci_method"])


def projector(basis, ivar):
    """G2's noise-weighted external-subspace projector (noisemodel.ArrayNoise.projector) on a channel subset."""
    w = 1.0 / ivar
    m = basis.T @ (w[:, None] * basis)
    return np.eye(len(w)) - basis @ np.linalg.solve(m, basis.T * w[None, :])


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Noise-model covariance validation against the measured Neuromag noise")
    ap.add_argument("--out", type=Path, default=OUT, help="output directory (default results/g2_covariance_validation)")
    ap.add_argument("--surrogates", type=int, default=40, help="stationary surrogates in the analysis band")
    ap.add_argument("--band-surrogates", type=int, default=10, help="stationary surrogates per sub-band")
    ap.add_argument("--splits", type=int, default=1000, help="random half splits of the gradiometer sites")
    ap.add_argument("--boot", type=int, default=1000, help="parcel-bootstrap resamples")
    ap.add_argument("--target-stride", type=int, default=1, help="use every n-th target (tests only)")
    ap.add_argument("--skip-opm", action="store_true", help="skip the OPM arrays and the implied ratios (tests only)")
    ap.add_argument("--replot", action="store_true", help="redraw the figures from the JSON in --out, no computation")
    args = ap.parse_args()
    t0 = time.time()
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.replot:
        import json

        figures(json.loads((out_dir / "covariance_validation.json").read_text()), out_dir)
        return
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((paths.CONFIGS / "g2_adult.toml").read_text())
    seeds = np.random.SeedSequence(cfg["sources"]["seed"]).spawn(5)
    rng_split, rng_surr, rng_band, rng_boot, rng_boot2 = (np.random.default_rng(s) for s in seeds)

    st = G.Study(cfg)  # targets, background grid, analysis filter: exactly G2's
    tsel = np.arange(0, st.nt, args.target_stride)
    groups = st.src.region[tsel]
    bads = cfg["sensors"]["bads"]
    info = neuromag.load_info("T3")
    squid = g2.Array("squid", info, neuromag.channel_kinds(info), None, dict(sites=102, channels=306))
    kinds = squid.kinds
    good = ~np.isin(squid.info.ch_names, bads)
    gidx = np.flatnonzero(good)
    names_good = [squid.info.ch_names[i] for i in gidx]
    k_good = kinds[good]
    mags, grads = np.flatnonzero(k_good == "mag"), np.flatnonzero(k_good == "grad")
    all_good = np.arange(len(gidx))
    sets = {"combined": all_good, "mag": mags, "grad": grads}
    pos_dev = np.array([c["loc"][:3] for c in squid.info["chs"]])[good]  # device frame [m]
    dist_mm = np.linalg.norm(pos_dev[:, None] - pos_dev[None], axis=-1) * 1e3
    site = neuromag.sensor_geometry(squid.info).site[good]  # index of the sensor location (0..101, magnetometer order)
    site_xyz = np.array([c["loc"][:3] for c in squid.info["chs"]])[kinds == "mag"]  # device frame, indexed by site
    log(f"targets {len(tsel)} of {st.nt}, good channels {len(gidx)} ({len(mags)} mag, {len(grads)} grad)")

    gt_all, gg_all = st.gains(squid)  # G2's verified full-resolution lead fields
    s_all = gt_all[:, tsel] * st.q
    s_good = s_all[good]

    # --- G2's measured noise and model (main band) ---------------------------------------------
    meas_g2 = g2.measured_noise(squid.info, st.filt, bads)
    env = meas_g2["environment"]
    target_var = float(np.nanmedian(meas_g2["brain"][good & (kinds == "grad")]))
    unit = {"independent": st.unit_brain(gg_all)}
    brain_scale = background.calibrate(unit["independent"], good & (kinds == "grad"), target_var)
    nz = st.noise(squid, gg_all, brain_scale, env)
    ivar = nz.intrinsic_var[good]
    basis = nz.ext_basis[good]
    M_br = nz.brain_cov[np.ix_(gidx, gidx)]
    M_er = np.diag(ivar) + nz.env_cov[np.ix_(gidx, gidx)]
    M_tot = M_er + M_br
    unit_good = unit["independent"][np.ix_(gidx, gidx)]
    log(f"brain scale {brain_scale:.6e} (G2 summary 3.535585942806047e-14)")

    # --- model variants (G2 sensitivity values) --------------------------------------------------
    variants = {"independent": dict(corr=None, scale=brain_scale, label="independent sources (G2 primary)"),
                "independent_mag_calibrated": dict(corr=None, scale=background.calibrate(unit["independent"], good & (kinds == "mag"),
                                                                                         float(np.nanmedian(meas_g2["brain"][good & (kinds == "mag")]))),
                                                   label="independent sources, calibrated on the magnetometers")}
    for lam in cfg["background"]["correlation_lengths_mm"]:
        if lam:
            u = st.unit_brain(gg_all, corr=lam * 1e-3)
            unit[f"correlated_{lam:g}mm"] = u
            variants[f"correlated_{lam:g}mm"] = dict(corr=lam * 1e-3, scale=background.calibrate(u, good & (kinds == "grad"), target_var),
                                                     label=f"correlated sources ({lam:g} mm), calibrated on the gradiometers")
    for name, v in variants.items():
        u = unit["independent"] if v["corr"] is None else unit[name]
        v["M_br"] = v["scale"] * u[np.ix_(gidx, gidx)]
        v["M_tot"] = M_er + v["M_br"]
    k_mag = variants["independent_mag_calibrated"]["scale"] / brain_scale
    log(f"variants: {list(variants)}; magnetometer/gradiometer calibration scale ratio {k_mag:.4f}")
    g2_rows = io.read_csv(paths.RESULTS / "g2" / "g2_targets.csv")

    # --- OPM detectability (built first: the dense array's construction is this script's memory peak)
    d_opm, opm_repro = {}, {}
    if not args.skip_opm:
        d_opm, opm_repro = opm_detectability(st, env, variants, tsel, g2_rows)
        log("OPM detectability done; G2 per-target values reproduced")

    # --- measured covariances (full) -------------------------------------------------------------
    task = Recording(neuromag.RAW_FILE, names_good)
    er = Recording(neuromag.ER_FILE, names_good)
    wins, codes = baseline_windows(task)
    widx = np.concatenate(wins)
    half_w = len(wins) // 2
    widx1, widx2 = np.concatenate(wins[:half_w]), np.concatenate(wins[half_w:])
    n_trim_er = int(TRIM_S * er.sfreq)
    eidx = np.arange(n_trim_er, er.n_times - n_trim_er)
    eidx1, eidx2 = eidx[:len(eidx) // 2], eidx[len(eidx) // 2:]

    def measure(filt):
        xt = task.filtered(filt)
        xe = er.filtered(filt)
        m = dict(base=cov(xt[:, widx]), er=cov(xe[:, eidx]), base1=cov(xt[:, widx1]), base2=cov(xt[:, widx2]),
                 er1=cov(xe[:, eidx1]), er2=cov(xe[:, eidx2]))
        m["brain"], m["brain1"], m["brain2"] = m["base"] - m["er"], m["base1"] - m["er1"], m["base2"] - m["er2"]
        m["base_tl"] = time_locked_removed_cov(xt, wins, codes)
        return m, xt, xe

    meas, xt, xe = measure(st.filt)
    checks = dict(
        diag_baseline_vs_g2_max_rel_dev=float(np.max(np.abs(np.diag(meas["base"]) / meas_g2["baseline"][good] - 1))),
        diag_empty_room_vs_g2_max_rel_dev=float(np.max(np.abs(np.diag(meas["er"]) / meas_g2["empty_room"][good] - 1))),
        n_windows=len(wins), n_baseline_samples=int(len(widx)), n_empty_room_samples=int(len(eidx)),
        window_event_codes={str(int(c)): int(np.sum(codes == c)) for c in np.unique(codes)},
        brain_scale=brain_scale, brain_scale_g2_summary=3.535585942806047e-14,
        brain_scale_rel_dev_from_g2=float(brain_scale / 3.535585942806047e-14 - 1))
    if checks["diag_baseline_vs_g2_max_rel_dev"] > 1e-9 or checks["diag_empty_room_vs_g2_max_rel_dev"] > 1e-9:
        raise RuntimeError(f"measured covariances do not reproduce G2's variances: {checks}")
    log(f"measured covariances: {len(wins)} windows ({len(widx)} samples), empty room {len(eidx)} samples; diagonals = G2's")

    card = cardiac_field(task, xt)
    meas["base_noheart"] = cov(xt[:, widx] - heart_template(card, widx))  # baseline windows with the average heartbeat removed
    log(f"cardiac field: {card['n_beats']} heartbeats, R-R {card['rr_s']:.3f} s")

    # noise normalisation for structure statistics: channels in units of the brochure sensor-noise SD
    dn = 1.0 / np.sqrt(ivar)

    def nrm(c):
        return c * np.outer(dn, dn)

    # --- 1-3: variances, spatial correlation, eigenstructure (main band) -------------------------
    structure = {}
    for comp, (a_key, b_mat) in dict(brain=("brain", M_br), total=("base", M_tot), empty_room=("er", M_er)).items():
        a = nrm(meas[a_key])
        b = nrm(b_mat)
        structure[comp] = {t: pair_stats(a, b, sets[t], dist_mm) for t in ("mag", "grad")}
        structure[comp]["combined_eigen"] = eig_stats(a, b, all_good)
        h1, h2 = nrm(meas[a_key + "1"]), nrm(meas[a_key + "2"])
        structure[comp]["split_half"] = {t: pair_stats(h1, h2, sets[t], dist_mm) for t in ("mag", "grad")}
        structure[comp]["split_half"]["combined_eigen"] = eig_stats(h1, h2, all_good)
    for t in ("mag", "grad"):  # reliability of the measured channel pattern (Spearman-Brown) and the disattenuated model correlation
        r_half = structure["brain"]["split_half"][t]["log_var_pearson"]
        rel = 2 * r_half / (1 + r_half) if r_half > 0 else np.nan
        structure["brain"][t]["measured_pattern_reliability"] = float(rel)
        structure["brain"][t]["log_var_pearson_disattenuated"] = float(structure["brain"][t]["log_var_pearson"] / np.sqrt(rel))
    ext_n = basis * dn[:, None]
    for comp, (a_key, b_mat) in dict(brain=("brain", M_br), total=("base", M_tot)).items():
        structure[comp]["excess"] = {t: excess_stats(nrm(meas[a_key]), nrm(b_mat), sets[t], ext_n) for t in ("mag", "grad", "combined")}
    # the heart's share of the measured brain noise and of the model's shortfall: the cardiac-locked average's own
    # second moment (vc_avg) and the variance its removal at the actual cardiac phase takes out of the baseline windows (vc)
    vb, vm = np.diag(meas["brain"]), np.diag(M_br)
    vc_avg, vc = np.diag(card["cov"]), np.diag(meas["base"]) - np.diag(meas["base_noheart"])
    cardiac = dict(n_beats=card["n_beats"], rr_s=card["rr_s"], heart_rate_bpm=60.0 / card["rr_s"],
                   method="mne.preprocessing.find_ecg_events (defaults; synthetic ECG from the magnetometers); cardiac-locked average "
                          "of the band-passed task recording over one median R-R cycle, subtracted from every baseline sample at its "
                          "offset from the nearest R peak; 'heart' variance = baseline variance before - after the subtraction")
    for t in ("mag", "grad"):
        idx = sets[t]
        unit_t = 1e15 if t == "mag" else 1e13
        q = np.linalg.qr(ext_n[idx])[0]
        uc, sc, _ = np.linalg.svd(nrm(card["cov"])[np.ix_(idx, idx)], hermitian=True)
        ue = eig_desc((nrm(meas["brain"]) - nrm(M_br))[np.ix_(idx, idx)])[1]
        pos = vb[idx] > 0
        shortfall = np.median(vb[idx]) - np.median(vm[idx])
        cardiac[t] = dict(median_rms_in_baseline=float(np.sqrt(max(np.median(vc[idx]), 0.0)) * unit_t),
                          max_rms_in_baseline=float(np.sqrt(max(vc[idx].max(), 0.0)) * unit_t),
                          median_rms_average_waveform=float(np.sqrt(np.median(vc_avg[idx])) * unit_t),
                          median_share_of_measured_brain_variance=float(np.median(vc[idx][pos] / vb[idx][pos])),
                          share_of_model_shortfall_at_the_median=(float(np.median(vc[idx]) / shortfall) if shortfall > 0.05 * np.median(vb[idx])
                                                                  else None),
                          brain_amp_ratio_model_over_measured_without_heart=float(np.sqrt(np.median(vm[idx]) / np.median(vb[idx] - vc[idx]))),
                          leading_component_power_share=float(sc[0] / sc.sum()),
                          leading_component_external_share=float(np.sum((q.T @ uc[:, 0]) ** 2)),
                          leading_component_cos2_with_leading_excess=float((uc[:, 0] @ ue[:, 0]) ** 2))
    nh_brain = nrm(meas["base_noheart"] - meas["er"])
    structure["brain_without_heart"] = {t: pair_stats(nh_brain, nrm(M_br), sets[t], dist_mm) for t in ("mag", "grad")}
    structure["brain_without_heart"]["excess"] = {t: excess_stats(nh_brain, nrm(M_br), sets[t], ext_n) for t in ("mag", "grad", "combined")}
    variant_structure = {name: {t: pair_stats(nrm(meas["brain"]), nrm(v["M_br"]), sets[t], dist_mm) for t in ("mag", "grad")}
                         for name, v in variants.items() if name != "independent"}
    heldout = held_out(np.diag(meas["brain"]), np.diag(unit_good), grads, mags, site, site_xyz, args.splits, rng_split)
    tl_brain = np.diag(meas["base_tl"]) - np.diag(meas["er"])
    time_locked = {t: dict(brain_amp_ratio_time_locked_removed_over_primary=float(np.sqrt(np.median(tl_brain[sets[t]])
                                                                                          / np.median(np.diag(meas["brain"])[sets[t]]))),
                           model_over_measured_time_locked_removed=float(np.sqrt(np.median(np.diag(M_br)[sets[t]]) / np.median(tl_brain[sets[t]]))))
                   for t in ("mag", "grad")}
    log("structure statistics done")

    # --- 5: detectability (main band) ------------------------------------------------------------
    p_good = projector(basis, ivar)
    d = {}
    for cs, idx in sets.items():
        sub = np.ix_(idx, idx)
        si = s_good[idx]
        d[("model", cs)] = detect(si, M_tot[sub])
        d[("model_ib", cs)] = detect(si, (np.diag(ivar) + M_br)[sub])
        d[("measured", cs)] = detect(si, meas["base"][sub])
        d[("hybrid", cs)] = detect(si, (meas["er"] + M_br)[sub])
        d[("measured_half1", cs)] = detect(si, meas["base1"][sub])
        d[("measured_half2", cs)] = detect(si, meas["base2"][sub])
        d[("measured_time_locked_removed", cs)] = detect(si, meas["base_tl"][sub])
        d[("measured_without_heart", cs)] = detect(si, meas["base_noheart"][sub])
        for name, v in variants.items():
            if name != "independent":
                d[(f"model_{name}", cs)] = detect(si, v["M_tot"][sub])
    # the projected condition (8-term external subspace removed, G2's projector on the 305 good channels)
    ps = p_good @ s_good
    d[("model_projected", "combined")] = detect(ps, p_good @ M_tot @ p_good.T)
    d[("measured_projected", "combined")] = detect(ps, p_good @ meas["base"] @ p_good.T)
    d[("hybrid_projected", "combined")] = detect(ps, p_good @ (meas["er"] + M_br) @ p_good.T)
    # G2 reproduction: the model on all 306 channels against the stored per-target values
    repro = {}
    for cs, m306 in (("combined", np.ones(306, bool)), ("mag", kinds == "mag"), ("grad", kinds == "grad")):
        for cond in ("intrinsic+brain", "intrinsic+brain+env"):
            idx = np.flatnonzero(m306)
            dd = detect(s_all[idx], nz.covariance(cond)[np.ix_(idx, idx)])
            stored = np.array([float(g2_rows[i][f"detect_squid_{cs}_{cond}"]) for i in tsel])
            repro[f"{cs}/{cond}"] = float(np.max(np.abs(dd - stored)))
            if cs == "combined" and cond == "intrinsic+brain+env":
                d[("model306", cs)] = dd
            if cs == "combined" and cond == "intrinsic+brain":
                d[("model306_ib", cs)] = dd
    checks["g2_detectability_max_abs_dev_vs_stored_csv"] = repro
    if max(repro.values()) > 6e-5:  # the CSV keeps 4 decimals
        raise RuntimeError(f"the model does not reproduce G2's stored Neuromag detectability: {repro}")
    log("detectability (data) done; G2 per-target values reproduced")

    # --- surrogates: finite-sample bias of the measured detectability, and the model null --------
    sur_t, sur_e = Surrogates(xt), Surrogates(xe)
    del xt, xe
    a_tot = recolour(M_tot, sur_t.cov_full)
    a_er = recolour(M_er, sur_e.cov_full)
    full_d = {}
    for cs, idx in sets.items():
        sub = np.ix_(idx, idx)
        full_d[("data", cs)] = detect(s_good[idx], sur_t.cov_full[sub])
        full_d[("data_hybrid", cs)] = detect(s_good[idx], (sur_e.cov_full + M_br)[sub])
    acc = {k: np.zeros(len(tsel)) for k in full_d}
    sur_rows, null_struct = [], []
    null_d = {k: [] for k in [("null", cs) for cs in sets] + [("null_hybrid", cs) for cs in sets]}
    for i in range(args.surrogates):
        cx = sur_t.sample(widx, rng_surr)
        ce = sur_e.sample(eidx, rng_surr)
        null_tot = a_tot @ cx @ a_tot.T
        null_er = a_er @ ce @ a_er.T
        row = {}
        for cs, idx in sets.items():
            sub = np.ix_(idx, idx)
            si = s_good[idx]
            dd = detect(si, cx[sub])
            acc[("data", cs)] += (dd / full_d[("data", cs)]) ** 2
            row[f"data/{cs}"] = float(np.median(dd / full_d[("data", cs)]))
            dh = detect(si, (ce + M_br)[sub])
            acc[("data_hybrid", cs)] += (dh / full_d[("data_hybrid", cs)]) ** 2
            row[f"data_hybrid/{cs}"] = float(np.median(dh / full_d[("data_hybrid", cs)]))
            null_d[("null", cs)].append(detect(si, null_tot[sub]))
            null_d[("null_hybrid", cs)].append(detect(si, (null_er + M_br)[sub]))
            row[f"null/{cs}"] = float(np.median(null_d[("null", cs)][-1] / d[("model", cs)]))
            row[f"null_hybrid/{cs}"] = float(np.median(null_d[("null_hybrid", cs)][-1] / d[("model", cs)]))
        row["null_projected/combined"] = float(np.median(detect(ps, p_good @ null_tot @ p_good.T) / d[("model_projected", "combined")]))
        sur_rows.append(row)
        nb = nrm(null_tot - null_er)
        null_struct.append(dict(brain={t: pair_stats(nb, nrm(M_br), sets[t], dist_mm) for t in ("mag", "grad")},
                                total={t: pair_stats(nrm(null_tot), nrm(M_tot), sets[t], dist_mm) for t in ("mag", "grad")},
                                empty_room={t: pair_stats(nrm(null_er), nrm(M_er), sets[t], dist_mm) for t in ("mag", "grad")}))
        for comp, mod in (("brain", M_br), ("total", M_tot)):
            nc = nb if comp == "brain" else nrm(null_tot)
            null_struct[-1][comp]["combined_eigen"] = eig_stats(nc, nrm(mod), all_good)
            null_struct[-1][comp]["excess"] = {t: excess_stats(nc, nrm(mod), sets[t], ext_n) for t in ("mag", "grad", "combined")}
        if (i + 1) % 10 == 0 or i + 1 == args.surrogates:
            log(f"surrogate {i + 1}/{args.surrogates}")
    bias = {k: np.sqrt(v / args.surrogates) for k, v in acc.items()}  # per-target finite-sample bias of d (estimated covariance)
    for cs in sets:
        d[("measured_corrected", cs)] = d[("measured", cs)] / bias[("data", cs)]
        d[("measured_time_locked_removed_corrected", cs)] = d[("measured_time_locked_removed", cs)] / bias[("data", cs)]
        d[("measured_without_heart_corrected", cs)] = d[("measured_without_heart", cs)] / bias[("data", cs)]
        d[("hybrid_corrected", cs)] = d[("hybrid", cs)] / bias[("data_hybrid", cs)]
        for i, row in enumerate(sur_rows):  # the null corrected as the data are
            row[f"null_corrected/{cs}"] = float(np.median(null_d[("null", cs)][i] / bias[("data", cs)] / d[("model", cs)]))
            row[f"null_hybrid_corrected/{cs}"] = float(np.median(null_d[("null_hybrid", cs)][i] / bias[("data_hybrid", cs)] / d[("model", cs)]))
    del null_d
    sur_summary = {k: [float(np.percentile([r[k] for r in sur_rows], q)) for q in (2.5, 50, 97.5)] for k in sur_rows[0]}
    null_summary = summarise_null(null_struct)
    null_spectra = {f"{comp}.{t}": np.median([ns[comp][t]["eig_a_normalised"] for ns in null_struct], axis=0).tolist()
                    for comp in ("brain", "total", "empty_room") for t in ("mag", "grad")}
    del sur_t, sur_e
    log("surrogates done")

    # --- detectability summaries ----------------------------------------------------------------
    det = dict(channels=dict(combined=int(len(all_good)), mag=int(len(mags)), grad=int(len(grads))), comparisons={})
    for cs in sets:
        for lab, num, den in (("measured_over_model", "measured", "model"), ("measured_corrected_over_model", "measured_corrected", "model"),
                              ("hybrid_over_model", "hybrid", "model"), ("hybrid_corrected_over_model", "hybrid_corrected", "model"),
                              ("measured_corrected_over_model_intrinsic_brain", "measured_corrected", "model_ib"),
                              ("half1_over_model", "measured_half1", "model"), ("half2_over_model", "measured_half2", "model"),
                              ("time_locked_removed_corrected_over_model", "measured_time_locked_removed_corrected", "model")):
            det["comparisons"][f"{lab}/{cs}"] = ratio_summary(d[(num, cs)], d[(den, cs)], groups, rng_boot, args.boot)
        if ("measured_without_heart_corrected", cs) in d:
            det["comparisons"][f"without_heart_corrected_over_model/{cs}"] = ratio_summary(
                d[("measured_without_heart_corrected", cs)], d[("model", cs)], groups, rng_boot, args.boot)
        det["comparisons"][f"finite_sample_bias/{cs}"] = dict(median=float(np.median(bias[("data", cs)])),
                                                             p5=float(np.percentile(bias[("data", cs)], 5)),
                                                             p95=float(np.percentile(bias[("data", cs)], 95)),
                                                             hybrid_median=float(np.median(bias[("data_hybrid", cs)])))
        for name in variants:
            if name != "independent":
                det["comparisons"][f"measured_corrected_over_model_{name}/{cs}"] = ratio_summary(
                    d[("measured_corrected", cs)], d[(f"model_{name}", cs)], groups, rng_boot, args.boot)
    det["comparisons"]["measured_over_model/combined/projected"] = ratio_summary(
        d[("measured_projected", "combined")], d[("model_projected", "combined")], groups, rng_boot, args.boot)
    det["comparisons"]["hybrid_over_model/combined/projected"] = ratio_summary(
        d[("hybrid_projected", "combined")], d[("model_projected", "combined")], groups, rng_boot, args.boot)
    det["comparisons"]["model305_over_model306/combined"] = ratio_summary(d[("model", "combined")], d[("model306", "combined")], groups,
                                                                         rng_boot, args.boot)
    det["surrogate_medians"] = sur_summary
    det["note"] = ("measured = sample covariance of the task baselines; corrected = divided per target by the finite-sample bias of "
                   "detectability with a covariance estimated from these samples (stationary surrogates with the measured cross-spectra); "
                   "null = the model covariance 'as it would be measured' (surrogates recoloured to the model), so measured/model "
                   "inside the null range is what an exact model would give; hybrid = measured empty room + model brain. "
                   "Parcel CIs cover the choice of targets only, not the finite recording (see surrogate_medians).")
    depth = st.src.depth_mm[tsel]
    lr = np.log2(d[("measured_corrected", "combined")] / d[("model", "combined")])
    det["by_depth_combined"] = G.binned(2 ** lr, depth, G.DEPTH_EDGES)
    lobes = st.src.lobe[tsel]
    det["by_lobe_combined"] = {}
    for lb in np.unique(lobes):
        m = lobes == lb
        row = dict(n_targets=int(m.sum()), measured_corrected=ratio_summary(d[("measured_corrected", "combined")][m], d[("model", "combined")][m],
                                                                             groups[m], rng_boot, args.boot))
        if ("measured_without_heart_corrected", "combined") in d:
            row["without_heart_corrected"] = ratio_summary(d[("measured_without_heart_corrected", "combined")][m], d[("model", "combined")][m],
                                                           groups[m], rng_boot, args.boot)
        det["by_lobe_combined"][lb] = row
    log("detectability summaries done")

    # --- 6: implied OPM / Neuromag ratios --------------------------------------------------------
    implication = None
    if not args.skip_opm:
        implication = implied_ratios(d_opm, opm_repro, variants, d, groups, rng_boot2, args.boot, k_mag, sets)
        log("implied OPM/Neuromag ratios done")

    # --- 4: sub-bands ----------------------------------------------------------------------------
    bands = sub_band_analysis(task, er, squid, good, gidx, ivar, unit_good, basis, sets, dist_mm, cfg, st, bads,
                              args.band_surrogates, rng_band, widx, eidx, meas)
    log("sub-bands done")

    # --- outputs -----------------------------------------------------------------------------------
    var_tab = {}
    for comp in ("brain", "total", "empty_room"):
        for t in ("mag", "grad"):
            r = structure[comp][t]
            var_tab[f"{comp}/{t}"] = dict(amp_ratio_of_medians_model_over_measured=r["amp_ratio_of_medians"],
                                          null95=null_summary.get(f"{comp}.{t}.amp_ratio_of_medians"),
                                          var_ratio_model_over_measured=r["var_ratio"],
                                          null_var_ratio_p5=null_summary.get(f"{comp}.{t}.var_ratio.p5"),
                                          null_var_ratio_p95=null_summary.get(f"{comp}.{t}.var_ratio.p95"),
                                          log_var_pearson=r["log_var_pearson"], log_var_spearman=r["log_var_spearman"],
                                          null_log_var_pearson=null_summary.get(f"{comp}.{t}.log_var_pearson"),
                                          split_half_log_var_pearson=structure[comp]["split_half"][t]["log_var_pearson"],
                                          n_nonpositive_measured=r["n_nonpositive"])
    summary = dict(
        status=STATUS,
        requests=REQUESTS,
        fitted_vs_predicted=FITTED_VS_PREDICTED,
        declared_choices=dict(sub_bands_hz=[list(b) for b in SUB_BANDS], sub_band_note=SUB_BAND_NOTE,
                              distance_bins_mm=DIST_EDGES_MM.tolist(), subspace_dimensions=list(K_SUBSPACE),
                              n_surrogates=args.surrogates, n_band_surrogates=args.band_surrogates, n_random_splits=args.splits,
                              n_parcel_bootstrap=args.boot, seed=cfg["sources"]["seed"], target_stride=args.target_stride,
                              structure_units="each channel divided by its brochure sensor-noise SD (unit-free; only the combined "
                                              "eigenstructure depends on it)",
                              surrogate_method="multivariate phase randomisation of the filtered continuous records (same random phase "
                                               "for every channel at each frequency)"),
        data=dict(task_file=neuromag.RAW_FILE, empty_room_file=neuromag.ER_FILE, band_hz=[cfg["band"]["l_freq_hz"], cfg["band"]["h_freq_hz"]],
                  window_s=list(WINDOW), trim_s=TRIM_S, bad_channels=bads, ssp="not applied (as G2)", **{k: checks[k] for k in
                                                                                                        ("n_windows", "n_baseline_samples",
                                                                                                         "n_empty_room_samples")}),
        model=dict(brain_scale=brain_scale, magnetometer_calibration_scale_ratio=k_mag, enbw_hz=st.enbw,
                   variants={k: dict(label=v["label"], scale=v["scale"], correlation_length_mm=None if v["corr"] is None else v["corr"] * 1e3)
                             for k, v in variants.items()},
                   environment_explained_fraction=env.explained_fraction),
        checks=checks,
        variances=var_tab,
        held_out=heldout,
        cardiac_check=cardiac,
        time_locked_check=dict(note="baseline windows after removing their time-locked average per event type (any residual evoked "
                                    "response from the preceding stimulus); the model is not refitted", **time_locked),
        structure=structure,
        structure_variants=variant_structure,
        structure_null=null_summary,
        structure_null_spectra=null_spectra,
        detectability=det,
        opm_implication=implication,
        sub_bands=bands,
        channels=dict(note="good channels, native units (T^2 for magnetometers, (T/m)^2 for gradiometers), variances in the analysis band",
                      names=names_good, kind=k_good.tolist(), site=site.tolist(), measured_baseline_var=np.diag(meas["base"]).tolist(),
                      measured_empty_room_var=np.diag(meas["er"]).tolist(), measured_brain_var=np.diag(meas["brain"]).tolist(),
                      model_brain_var=np.diag(M_br).tolist(), model_empty_room_var=np.diag(M_er).tolist(),
                      cardiac_var=np.diag(card["cov"]).tolist()),
    )
    summary["plain_summary"] = plain_summary(summary)
    summary["runtime_s"] = time.time() - t0
    summary["peak_rss_gb"] = peak_rss_gb()
    summary["provenance"] = dict(commit=io.RUN_COMMIT, mne_version=mne.__version__, numpy_version=np.__version__)
    io.write_json(summary, out_dir / "covariance_validation.json")
    write_targets_csv(out_dir / "covariance_validation_targets.csv", st, tsel, d, d_opm)
    figures(summary, out_dir)
    log(f"done in {time.time() - t0:.0f} s, peak RSS {peak_rss_gb():.2f} GB -> {out_dir}")


def figures(summary, out_dir):
    RS.apply()
    figure_structure(summary, out_dir / "Figure_covariance_structure.png")
    figure_bands(summary, out_dir / "Figure_covariance_bands.png")
    figure_detectability(summary, out_dir / "Figure_covariance_detectability.png")


# ----------------------------------------------------------------------------------------------
OPM_ARRAYS = ("opm_dense", "opm_matched")


def opm_detectability(st, env, variants, tsel, g2_rows):
    """Model detectability of the G2 targets for the dense and the site-matched OPM arrays (G2's
    construction, lead fields and noise), for every background variant; checked against G2's
    stored per-target values."""
    d_opm, repro = {}, {}
    for name in OPM_ARRAYS:
        arr = g2.dense_opm(st.subject, st.dig, name) if name in g2.DENSE else g2.matched_opm(st.subject, st.dig)
        gt, gg = st.gains(arr)
        s = gt[:, tsel] * st.q
        for vname, v in variants.items():
            nz_o = st.noise(arr, gg, v["scale"], env, corr=v["corr"])
            for cond in ("intrinsic+brain", "intrinsic+brain+env", "projected"):
                if vname != "independent" and cond == "projected":
                    continue
                d_opm[(name, vname, cond)] = detect(nz_o.signal(s, cond), nz_o.covariance(cond))
        for cond in ("intrinsic+brain", "intrinsic+brain+env"):
            stored = np.array([float(g2_rows[i][f"detect_{name}_opm_{cond}"]) for i in tsel])
            repro[f"{name}/{cond}"] = float(np.max(np.abs(d_opm[(name, "independent", cond)] - stored)))
        del arr, gt, gg
    if max(repro.values()) > 6e-5:  # the CSV keeps 4 decimals
        raise RuntimeError(f"the OPM model does not reproduce G2's stored detectability: {repro}")
    return d_opm, repro


def implied_ratios(d_opm, repro, variants, d, groups, rng, n_boot, k_mag, sets):
    """OPM/Neuromag detectability ratios with the measured Neuromag covariance, under stated assumptions."""
    arrays = OPM_ARRAYS
    out = dict(checks=dict(g2_opm_detectability_max_abs_dev_vs_stored_csv=repro, magnetometer_variance_excess_k=k_mag),
               assumptions={
                   "model_306": "G2 as published: model covariance for both systems, Neuromag with all 306 channels",
                   "model_305": "model covariance for both systems, Neuromag without the bad channel (the reference for the rows below)",
                   "S1_same_relative_change": "the OPM's real covariance changes its detectability by the same factor, target by target, "
                                              "as the measured covariance changes Neuromag's (the model's error is common to both systems): "
                                              "the ratio is that of the model (identical to model_305 by construction)",
                   "S2_same_variance_excess_both_modelled": "both systems' cortical background scaled by the measured magnetometer variance "
                                                            "excess k (model structure kept; = G2's magnetometer-calibrated sensitivity)",
                   "S3_neuromag_measured_opm_variance_excess": "Neuromag with its measured covariance (finite-sample corrected); the OPM with "
                                                               "its cortical background scaled by k: the OPM bears the full magnetometer "
                                                               "variance excess as cortical-like noise and none of the measured structure "
                                                               "(a pessimistic bound for the OPM)",
                   "S4_neuromag_measured_opm_as_modelled": "Neuromag with its measured covariance (finite-sample corrected); the OPM exactly "
                                                           "as modelled (the model's error is Neuromag-specific)",
                   "S5_neuromag_measured_empty_room": "Neuromag with the measured empty room + model brain (finite-sample corrected); the "
                                                      "OPM as modelled (isolates the sensor and room part of Neuromag's noise)",
                   "variants": "the OPM and Neuromag background both from the named variant (model rows), or Neuromag measured and the OPM "
                               "from the variant (measured rows)"},
               ratios={})
    cond = "intrinsic+brain+env"  # like for like with the measured covariance, which holds the room field
    for name in arrays:
        for cs in sets:
            r = out["ratios"]
            dn = d_opm[(name, "independent", cond)]
            dk = d_opm[(name, "independent_mag_calibrated", cond)]
            key = f"{name}/{cs}"
            if cs == "combined":
                r[f"{key}/model_306"] = ratio_summary(dn, d[("model306", cs)], groups, rng, n_boot)
                r[f"{key}/model_306/intrinsic+brain"] = ratio_summary(d_opm[(name, "independent", "intrinsic+brain")],
                                                                      d[("model306_ib", cs)], groups, rng, n_boot)
            r[f"{key}/model_305"] = ratio_summary(dn, d[("model", cs)], groups, rng, n_boot)
            r[f"{key}/S1_same_relative_change"] = r[f"{key}/model_305"]
            r[f"{key}/S2_same_variance_excess_both_modelled"] = ratio_summary(dk, d[("model_independent_mag_calibrated", cs)], groups, rng, n_boot)
            r[f"{key}/S3_neuromag_measured_opm_variance_excess"] = ratio_summary(dk, d[("measured_corrected", cs)], groups, rng, n_boot)
            r[f"{key}/S4_neuromag_measured_opm_as_modelled"] = ratio_summary(dn, d[("measured_corrected", cs)], groups, rng, n_boot)
            r[f"{key}/S4_uncorrected"] = ratio_summary(dn, d[("measured", cs)], groups, rng, n_boot)
            r[f"{key}/S5_neuromag_measured_empty_room"] = ratio_summary(dn, d[("hybrid_corrected", cs)], groups, rng, n_boot)
            for vname in variants:
                if vname.startswith("correlated"):
                    dv = d_opm[(name, vname, cond)]
                    r[f"{key}/variant_{vname}/model"] = ratio_summary(dv, d[(f"model_{vname}", cs)], groups, rng, n_boot)
                    r[f"{key}/variant_{vname}/neuromag_measured"] = ratio_summary(dv, d[("measured_corrected", cs)], groups, rng, n_boot)
        r = out["ratios"]
        r[f"{name}/combined/projected/model_305"] = ratio_summary(d_opm[(name, "independent", "projected")],
                                                                  d[("model_projected", "combined")], groups, rng, n_boot)
        r[f"{name}/combined/projected/S4_uncorrected"] = ratio_summary(d_opm[(name, "independent", "projected")],
                                                                       d[("measured_projected", "combined")], groups, rng, n_boot)
    out["condition"] = ("OPM: sensor + brain + room (G2 'intrinsic+brain+env'), like for like with the measured Neuromag covariance; "
                        "model_306/intrinsic+brain is the G2 headline condition; 'projected' rows remove the 8-term external subspace "
                        "(Neuromag measured rows there are not finite-sample corrected)")
    return out


# ----------------------------------------------------------------------------------------------
def sub_band_analysis(task, er, squid, good, gidx, ivar_main, unit_good, basis, sets, dist_mm, cfg, st, bads, n_surr, rng,
                      widx, eidx, meas_main):
    """Per sub-band: measured empty room, baseline and brain covariances; the model with its scale
    refitted on the gradiometers and the room field refitted to the empty room (as
    scripts/g2_band_sensitivity.py); magnetometer prediction, channel patterns, dominant subspaces."""
    kinds_good = squid.kinds[good]
    mags, grads = sets["mag"], sets["grad"]
    dn = 1.0 / np.sqrt(ivar_main)
    asd2 = np.array([g2.SQUID_ASD[k] ** 2 for k in kinds_good])
    out = {}
    main_brain = meas_main["brain"] * np.outer(dn, dn)
    main_u = {t: eig_desc(main_brain[np.ix_(sets[t], sets[t])])[1] for t in ("mag", "grad")}
    unit_diag = np.diag(unit_good)
    full_key = f"{cfg['band']['l_freq_hz']:g}-{cfg['band']['h_freq_hz']:g}Hz"
    diag_base, diag_er = {}, {}
    for lo, hi in ((cfg["band"]["l_freq_hz"], cfg["band"]["h_freq_hz"]),) + SUB_BANDS:
        key = f"{lo:g}-{hi:g}Hz"
        filt = noise.AnalysisFilter(fs=task.sfreq, l_freq=lo, h_freq=hi, order=cfg["band"]["order"])
        enbw = filt.enbw()
        xt, xe = task.filtered(filt), er.filtered(filt)
        cb, ce = cov(xt[:, widx]), cov(xe[:, eidx])
        half = len(widx) // 2
        cb1, cb2 = cov(xt[:, widx[:half]]), cov(xt[:, widx[half:]])
        ce1, ce2 = cov(xe[:, eidx[:len(eidx) // 2]]), cov(xe[:, eidx[len(eidx) // 2:]])
        env_b = environment.fit_empty_room(er.raw, squid.info, filt, bads=list(bads))
        m_er = np.diag(asd2 * enbw) + (basis @ env_b.coef_cov @ basis.T)
        brain = cb - ce
        diag_base[key], diag_er[key] = np.diag(cb).copy(), np.diag(ce).copy()
        scale = float(np.median(np.diag(brain)[grads]) / np.median(unit_diag[grads]))  # G2's rule (nanmedian over good gradiometers)
        m_br = scale * unit_good
        row = dict(enbw_hz=enbw, brain_scale=scale)
        for t, idx in (("mag", mags), ("grad", grads)):
            unit_t = 1e15 if t == "mag" else 1e13  # fT, fT/cm
            v = {lab: np.diag(c)[idx] for lab, c in (("empty_room_measured", ce), ("empty_room_model", m_er), ("baseline_measured", cb),
                                                      ("brain_measured", brain), ("brain_model", m_br))}
            row[t] = {lab: dict(median_amplitude=float(np.sqrt(np.median(x)) * unit_t),
                                median_asd=float(np.sqrt(np.median(x) / enbw) * unit_t)) for lab, x in v.items()}
            row[t]["brain_amp_ratio_model_over_measured"] = float(np.sqrt(np.median(v["brain_model"]) / np.median(v["brain_measured"])))
            row[t]["empty_room_amp_ratio_model_over_measured"] = float(np.sqrt(np.median(v["empty_room_model"]) / np.median(v["empty_room_measured"])))
            ps = pair_stats((brain) * np.outer(dn, dn), m_br * np.outer(dn, dn), idx, dist_mm)
            row[t]["brain_log_var_pearson"] = ps["log_var_pearson"]
            row[t]["brain_corr_pattern_pearson"] = ps["corr_pattern_pearson"]
            row[t]["brain_subspace_overlap_k5"] = ps["subspace"]["5"]["overlap"]
            row[t]["brain_subspace_overlap_k10"] = ps["subspace"]["10"]["overlap"]
            hs = pair_stats((cb1 - ce1) * np.outer(dn, dn), (cb2 - ce2) * np.outer(dn, dn), idx, dist_mm)
            row[t]["split_half_log_var_pearson"] = hs["log_var_pearson"]
            row[t]["split_half_subspace_overlap_k5"] = hs["subspace"]["5"]["overlap"]
            ub = eig_desc((brain * np.outer(dn, dn))[np.ix_(idx, idx)])[1]
            row[t]["overlap_k5_with_full_band_measured"] = subspace_overlap(ub, main_u[t], 5)[0]
            row[t]["n_nonpositive_brain"] = int(np.sum(v["brain_measured"] <= 0))
        row["mag_over_grad_brain_variance_measured_over_model"] = float(
            (np.median(np.diag(brain)[mags]) / np.median(np.diag(brain)[grads])) / (np.median(unit_diag[mags]) / np.median(unit_diag[grads])))
        # null: the band's model as it would be measured
        if n_surr:
            sur_t, sur_e = Surrogates(xt), Surrogates(xe)
            m_tot = m_er + m_br
            a_tot, a_er = recolour(m_tot, sur_t.cov_full), recolour(m_er, sur_e.cov_full)
            nulls = []
            for _ in range(n_surr):
                nb = a_tot @ sur_t.sample(widx, rng) @ a_tot.T - a_er @ sur_e.sample(eidx, rng) @ a_er.T
                scale_n = np.median(np.diag(nb)[grads]) / np.median(unit_diag[grads])  # refit as for the data
                r_ = {}
                for t, idx in (("mag", mags), ("grad", grads)):
                    ps = pair_stats(nb * np.outer(dn, dn), scale_n * unit_good * np.outer(dn, dn), idx, dist_mm, full=False)
                    r_[f"{t}.brain_amp_ratio_model_over_measured"] = float(np.sqrt(np.median(scale_n * unit_diag[idx]) / np.median(np.diag(nb)[idx])))
                    r_[f"{t}.brain_log_var_pearson"] = ps["log_var_pearson"]
                nulls.append(r_)
            row["null95"] = {k: [float(np.percentile([n_[k] for n_ in nulls], q)) for q in (2.5, 50, 97.5)] for k in nulls[0]}
            del sur_t, sur_e
        del xt, xe
        out[key] = row
        log(f"band {key}: mag brain model/measured {row['mag']['brain_amp_ratio_model_over_measured']:.3f}, "
            f"pattern r mag {row['mag']['brain_log_var_pearson']:.2f} grad {row['grad']['brain_log_var_pearson']:.2f}")
    sub = [k for k in out if k != full_key]
    out["closure_median"] = {f"{lab}/{t}": float(np.median(sum(dg[k][sets[t]] for k in sub) / dg[full_key][sets[t]]))
                             for lab, dg in (("baseline", diag_base), ("empty_room", diag_er)) for t in ("mag", "grad")}
    out["note"] = ("the model has no brain spectrum: its cortical background has one spatial covariance at every frequency, rescaled per "
                   "band on the gradiometers; it therefore predicts a band-independent magnetometer/gradiometer brain-variance ratio and "
                   "channel pattern. 'closure_median' = per channel, sum of the sub-band variances / full-band variance (median over "
                   "channels; filter-bank check: the sub-bands' skirts overlap).")
    return out


# ----------------------------------------------------------------------------------------------
def plain_summary(s) -> list[str]:
    """Data-driven sentences (every number read from the summary)."""
    v, det = s["variances"], s["detectability"]["comparisons"]
    ho = s["held_out"]["random_site_splits"]
    sp = s["held_out"]["spatial_splits"]
    ex = s["structure"]["brain"]["excess"]["mag"]
    card = s["cardiac_check"]
    b = s["sub_bands"]
    bands = [k for k in b if k.endswith("Hz")][1:]
    lines = [
        f"Fitted: one scale of the cortical background (the median good-gradiometer brain-noise variance) and the room field's 8 x 8 "
        f"coefficient covariance (empty room). Everything below except the gradiometer median is an independent prediction.",
        f"Magnetometer brain-noise amplitude (median): model/measured {v['brain/mag']['amp_ratio_of_medians_model_over_measured']:.2f}; an exact "
        f"model, measured the same way, would give {v['brain/mag']['null95'][0]:.2f}-{v['brain/mag']['null95'][2]:.2f}.",
        f"Gradiometers: the scale fitted on half the sites predicts the other half's median amplitude to "
        f"{ho['heldout_amp_ratio_of_medians']['p50']:.2f} [{ho['heldout_amp_ratio_of_medians']['p2.5']:.2f}, "
        f"{ho['heldout_amp_ratio_of_medians']['p97.5']:.2f}] over random halves, but "
        f"{sp['fit_inferior_predict_superior']['heldout_amp_ratio_of_medians']:.2f} (fitted below, predicting above) and "
        f"{sp['fit_superior_predict_inferior']['heldout_amp_ratio_of_medians']:.2f} (fitted above, predicting below).",
        f"Per-channel pattern (Pearson r of log variances, model vs measured): magnetometers {v['brain/mag']['log_var_pearson']:.2f}, "
        f"gradiometers {v['brain/grad']['log_var_pearson']:.2f}, against a split-half reliability of the measurement of "
        f"{v['brain/mag']['split_half_log_var_pearson']:.2f} and {v['brain/grad']['split_half_log_var_pearson']:.2f}; per-channel "
        f"model/measured variance 5th-95th percentile {v['brain/mag']['var_ratio_model_over_measured']['p5']:.2f}-"
        f"{v['brain/mag']['var_ratio_model_over_measured']['p95']:.2f} (magnetometers) and "
        f"{v['brain/grad']['var_ratio_model_over_measured']['p5']:.2f}-{v['brain/grad']['var_ratio_model_over_measured']['p95']:.2f} "
        f"(gradiometers), where an exact model would give about {v['brain/mag']['null_var_ratio_p5'][1]:.2f}-"
        f"{v['brain/mag']['null_var_ratio_p95'][1]:.2f}.",
        f"The magnetometers' excess over the model is concentrated: its leading 3 components hold "
        f"{ex['top_k_share_of_positive_part']['3']:.0%} of it, and {ex['leading_component_external_share']['1']:.0%} of the leading "
        f"component lies in the 8-dimensional far-field (external) subspace (a random direction: {ex['external_share_random']:.0%}).",
        f"The heart (cardiac-locked average, {card['n_beats']} beats) adds {card['mag']['median_rms_in_baseline']:.0f} fT RMS (median "
        f"magnetometer) to the baseline windows: {card['mag']['median_share_of_measured_brain_variance']:.0%} of the measured magnetometer "
        f"brain variance and {card['mag']['share_of_model_shortfall_at_the_median']:.0%} of the model's shortfall; without it the "
        f"magnetometer prediction is {card['mag']['brain_amp_ratio_model_over_measured_without_heart']:.2f} of the measured amplitude and the "
        f"channel-pattern correlation rises to {s['structure']['brain_without_heart']['mag']['log_var_pearson']:.2f} (magnetometers) and "
        f"{s['structure']['brain_without_heart']['grad']['log_var_pearson']:.2f} (gradiometers).",
        "Magnetometer prediction per band (model/measured amplitude, model rescaled per band on the gradiometers): "
        + ", ".join(f"{k[:-2]} Hz {b[k]['mag']['brain_amp_ratio_model_over_measured']:.2f}" for k in bands) + ".",
    ]
    for cs in ("combined", "mag", "grad"):
        c = det[f"measured_corrected_over_model/{cs}"]
        h = det[f"hybrid_corrected_over_model/{cs}"]
        nh = det.get(f"without_heart_corrected_over_model/{cs}")
        nl = s["detectability"]["surrogate_medians"][f"null_corrected/{cs}"]
        raw = det[f"measured_over_model/{cs}"]["ratio"]
        lines.append(f"Neuromag {cs}: detectability with the measured covariance / model {c['ratio']:.3f} [{c['ci95'][0]:.3f}, "
                     f"{c['ci95'][1]:.3f}] after the finite-sample correction (uncorrected {raw:.3f}; an exact model, measured the same "
                     f"way, would give {nl[0]:.3f}-{nl[2]:.3f}); measured empty room + model brain {h['ratio']:.3f}"
                     + (f"; heart removed {nh['ratio']:.3f}" if nh else "") + ".")
    bl = s["detectability"]["by_lobe_combined"]
    lines.append("By lobe (all channels, measured/model): " + ", ".join(f"{lb} {r['measured_corrected']['ratio']:.2f}" for lb, r in bl.items()) + ".")
    if s["opm_implication"]:
        r = s["opm_implication"]["ratios"]
        for a in ("opm_dense", "opm_matched"):
            parts = [f"{lab} {r[f'{a}/combined/{k}']['ratio']:.3f} [{r[f'{a}/combined/{k}']['ci95'][0]:.3f}, {r[f'{a}/combined/{k}']['ci95'][1]:.3f}]"
                     for k, lab in (("model_306", "both modelled"), ("S1_same_relative_change", "S1 same relative change in both"),
                                    ("S2_same_variance_excess_both_modelled", "S2 both backgrounds x k"),
                                    ("S5_neuromag_measured_empty_room", "S5 Neuromag measured empty room"),
                                    ("S4_neuromag_measured_opm_as_modelled", "S4 Neuromag measured, OPM as modelled"),
                                    ("S3_neuromag_measured_opm_variance_excess", "S3 Neuromag measured, OPM background x k"))]
            lines.append(f"{a} / Neuromag (all channels, sensor + brain + room): " + "; ".join(parts) + ".")
    return lines


def write_targets_csv(path, st, tsel, d, d_opm):
    cols = [(f"detect_neuromag_{cs}_{lab}", d[(lab, cs)]) for cs in ("combined", "mag", "grad")
            for lab in ("model", "measured", "measured_corrected", "hybrid", "hybrid_corrected", "measured_without_heart_corrected")]
    cols += [(f"detect_{a}_{v}_{c}".replace("+", "_"), x) for (a, v, c), x in d_opm.items() if v in ("independent", "independent_mag_calibrated")]
    with open(path, "w", newline="") as fh:
        io.csv_status(fh, STATUS)
        wr = csv.writer(fh)
        wr.writerow(["hemi", "vertno", "depth_mm", "region", "lobe"] + [c for c, _ in cols])
        for j, i in enumerate(tsel):
            v = st.src.target[i]
            wr.writerow([int(st.cortex.hemi[v]), int(st.cortex.vertno[v]), f"{st.src.depth_mm[i]:.4f}", st.src.region[i], st.src.lobe[i]]
                        + [f"{x[j]:.7g}" for _, x in cols])


# ----------------------------------------------------------------------------------------------
# figures (plain manuscript terms)
C_MEAS, C_MODEL, C_NULL, C_HALF = "#000000", "#0072B2", "#999999", "#E69F00"
C_VAR = {"correlated_5mm": "#009E73", "correlated_10mm": "#CC79A7"}
TYPE_LABEL = {"mag": "magnetometers", "grad": "gradiometers"}
LOBE_LABEL = {"frontal": "frontal", "parietal": "parietal", "temporal": "temporal", "occipital": "occipital", "insula": "insula",
              "cingulate": "cingulate", "other": "medial wall"}


def _bins(curve):
    return np.array([12.5 if r["bin"] == "same site" else 0.5 * sum(map(float, r["bin"].split("-"))) for r in curve])


def _log_ticks(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def figure_structure(s, path):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    st, nl, heldout = s["structure"], s["structure_null"], s["held_out"]
    ch = s["channels"]
    meas_var, model_var = np.array(ch["measured_brain_var"]), np.array(ch["model_brain_var"])
    kind = np.array(ch["kind"])
    mags, grads = np.flatnonzero(kind == "mag"), np.flatnonzero(kind == "grad")
    fig, axs = plt.subplots(3, 3, figsize=(RS.FULL_W, 8.0))
    # (a) channel variances
    ax = axs[0, 0]
    for idx, t, mk, unit_t, col in ((mags, "mag", "o", 1e15, C_MEAS), (grads, "grad", "s", 1e13, C_MODEL)):
        ok = meas_var[idx] > 0
        ax.loglog(np.sqrt(meas_var[idx][ok]) * unit_t, np.sqrt(model_var[idx][ok]) * unit_t, mk, ms=2.2, mfc="none", color=col,
                  label=f"{TYPE_LABEL[t]} ({'fT' if t == 'mag' else 'fT/cm'})")
    lim = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lim, lim, color="0.6", lw=0.7)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("measured brain noise (RMS)")
    ax.set_ylabel("model brain noise (RMS)")
    ax.set_title("a  Brain noise per channel", loc="left")
    ax.legend(loc="upper left", fontsize=6, handletextpad=0.2)
    # (b) held-out prediction
    ax = axs[0, 1]
    rs = heldout["random_site_splits"]
    items = [("held-out\ngradiometers", rs["heldout_amp_ratio_of_medians"]), ("magnetometers", rs["mag_amp_ratio_of_medians"])]
    for j, (lab, q) in enumerate(items):
        ax.errorbar([q["p50"]], [j], xerr=[[q["p50"] - q["p2.5"]], [q["p97.5"] - q["p50"]]], fmt="o", color=C_MODEL, ms=4, capsize=2,
                    label="random halves of the sites (95 % range)" if j == 0 else None)
    for k, r in enumerate(heldout["spatial_splits"].values()):
        ax.plot(r["heldout_amp_ratio_of_medians"], 0.3, "|", color=C_HALF, ms=8, label="left/right, front/back, upper/lower halves" if k == 0 else None)
        ax.plot(r["mag_amp_ratio_of_medians"], 1.3, "|", color=C_HALF, ms=8)
    ax.axvline(1.0, color="0.6", lw=0.7)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([lab for lab, _ in items])
    ax.set_ylim(-0.6, 3.0)
    ax.set_xlabel("model / measured (median RMS)")
    ax.set_title("b  Prediction from half the sites", loc="left")
    ax.legend(loc="upper center", fontsize=5.5, handletextpad=0.3)
    # (c) per-channel spread vs an exact model
    ax = axs[0, 2]
    for j, t in enumerate(("mag", "grad")):
        r = st["brain"][t]["var_ratio"]
        ax.plot([r["p5"], r["p95"]], [j, j], color=C_MEAS, lw=2, label="measured" if j == 0 else None)
        ax.plot([r["p50"]], [j], "o", color=C_MEAS, ms=4)
        n5, n50, n95 = (nl[f"brain.{t}.var_ratio.{q}"][1] for q in ("p5", "p50", "p95"))
        ax.plot([n5, n95], [j + 0.3, j + 0.3], color=C_NULL, lw=2, label="exact model, as it would be measured" if j == 0 else None)
        ax.plot([n50], [j + 0.3], "o", color=C_NULL, ms=4)
    _log_ticks(ax, [0.1, 0.2, 0.5, 1, 2])
    ax.axvline(1.0, color="0.6", lw=0.7)
    ax.set_yticks([0.15, 1.15])
    ax.set_yticklabels(["magnetometers", "gradiometers"])
    ax.set_ylim(-0.5, 2.2)
    ax.set_xlabel("model / measured variance per channel")
    ax.set_title("c  Spread over channels (5-95 %)", loc="left")
    ax.legend(loc="upper left", fontsize=5.5)
    # (d, e) spatial correlation vs distance
    for j, t in enumerate(("mag", "grad")):
        ax = axs[1, j]
        cur = st["brain"][t]["corr_vs_distance"]
        x = _bins(cur)
        key = "median_r" if t == "mag" else "median_abs_r"
        ok = np.array([key + "_a" in r for r in cur])
        ax.plot(x[ok], [r[key + "_a"] for r, o in zip(cur, ok) if o], "o-", color=C_MEAS, ms=3, label="measured")
        hc = st["brain"]["split_half"][t]["corr_vs_distance"]
        ax.plot(x[ok], [r[key + "_a"] for r, o in zip(hc, ok) if o], ":", color=C_HALF, lw=1.2, label="measured, 2nd half")
        ax.plot(x[ok], [r[key + "_b"] for r, o in zip(cur, ok) if o], "s-", color=C_MODEL, ms=3, label="model, independent")
        for vname, col in C_VAR.items():
            if vname in s["structure_variants"]:
                cv = s["structure_variants"][vname][t]["corr_vs_distance"]
                ax.plot(x[ok], [r[key + "_b"] for r, o in zip(cv, ok) if o], "-", color=col, lw=1,
                        label=f"correlated {vname.split('_')[1].replace('mm', ' mm')}")
        ax.axhline(0, color="0.8", lw=0.6)
        ax.set_xlabel("sensor distance (mm)" if t == "mag" else "sensor distance (mm; 1st point: same site)")
        ax.set_ylabel("median correlation" if t == "mag" else "median |correlation|")
        ax.set_title(f"{'de'[j]}  Correlation, {TYPE_LABEL[t]}", loc="left")
        if j == 0:
            ax.set_ylim(-1.0, 1.0)
            ax.legend(fontsize=5, loc="lower left", ncol=2, columnspacing=0.8, handlelength=1.5)
    # (f) agreement of the pairwise correlations vs distance
    ax = axs[1, 2]
    for t, mk in (("mag", "o"), ("grad", "s")):
        cur = st["brain"][t]["corr_vs_distance"]
        x = _bins(cur)
        ok = np.array(["agreement" in r for r in cur])
        nn = [nl.get(f"brain.{t}.corr_vs_distance.{r['bin']}.agreement") for r, o in zip(cur, ok) if o]
        if all(v is not None for v in nn):
            ax.fill_between(x[ok], [v[0] for v in nn], [v[2] for v in nn], color=C_NULL, alpha=0.35, lw=0)
        hc = st["brain"]["split_half"][t]["corr_vs_distance"]
        ax.plot(x[ok], [r["agreement"] for r, o in zip(hc, ok) if o], mk + ":", color=C_HALF, ms=2.5, mfc="none")
        ax.plot(x[ok], [r["agreement"] for r, o in zip(cur, ok) if o], mk + "-", color=C_MEAS, ms=3)
    ax.set_ylim(-0.2, 1.05)
    ax.set_xlabel("sensor distance (mm)")
    ax.set_ylabel("agreement of correlations (r)")
    ax.set_title("f  Correlation pattern", loc="left")
    # (g, h) eigenvalue spectra
    for j, t in enumerate(("mag", "grad")):
        ax = axs[2, j]
        ea, eb = np.array(st["brain"][t]["eig_a_normalised"]), np.array(st["brain"][t]["eig_b_normalised"])
        k = np.arange(1, len(ea) + 1)
        ax.semilogy(k[ea > 0], ea[ea > 0], "o", color=C_MEAS, ms=2.5, label="measured")
        ax.semilogy(k, eb, "-", color=C_MODEL, lw=1.2, label="model")
        nsp = np.array(s["structure_null_spectra"][f"brain.{t}"])
        ax.semilogy(k[nsp > 0], nsp[nsp > 0], "--", color=C_MODEL, lw=1, label="model, as it would be measured")
        ax.set_xlabel("component")
        ax.set_ylabel("share of brain-noise variance")
        ax.set_title(f"{'gh'[j]}  Eigenvalues, {TYPE_LABEL[t]}", loc="left")
        if j == 0:
            ax.legend(fontsize=5.5, loc="lower left")
    # (i) subspace overlap
    ax = axs[2, 2]
    for t, mk in (("mag", "o"), ("grad", "s")):
        sub = st["brain"][t]["subspace"]
        ks = [int(k) for k in sub]
        nn = [nl.get(f"brain.{t}.subspace.{k}.overlap") for k in ks]
        if all(v is not None for v in nn):
            ax.fill_between(ks, [v[0] for v in nn], [v[2] for v in nn], color=C_NULL, alpha=0.35, lw=0)
        hs = st["brain"]["split_half"][t]["subspace"]
        ax.plot(ks, [hs[str(k)]["overlap"] for k in ks], mk + ":", color=C_HALF, ms=2.5, mfc="none")
        ax.plot(ks, [sub[str(k)]["overlap"] for k in ks], mk + "-", color=C_MEAS, ms=3)
        ax.plot(ks, [sub[str(k)]["random_overlap"] for k in ks], "-", color="0.75", lw=0.8)
    _log_ticks(ax, [1, 2, 5, 10, 20, 40])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("dimension of the dominant subspace")
    ax.set_ylabel("overlap (mean cos$^2$)")
    ax.set_title("i  Dominant subspaces", loc="left")
    handles = [Line2D([], [], color=C_MEAS, marker="o", ms=3, label="model vs measured: magnetometers (o), gradiometers (s)"),
               Line2D([], [], color=C_HALF, ls=":", marker="o", mfc="none", ms=2.5, label="first vs second half of the measured data"),
               Patch(color=C_NULL, alpha=0.35, label="exact model, as it would be measured (95 %)"),
               Line2D([], [], color="0.75", lw=0.8, label="random subspaces")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=6.5, title="panels f and i", title_fontsize=6.5)
    fig.suptitle("Measured and modelled Neuromag brain noise (task baseline minus empty room), MNE sample subject, 1-40 Hz", fontsize=9)
    fig.tight_layout(rect=(0, 0.065, 1, 0.98))
    fig.savefig(path)
    plt.close(fig)


def _nice_log_y(ax):
    lo, hi = ax.get_ylim()
    ticks = [t for t in (0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500) if lo <= t <= hi]
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def figure_bands(s, path):
    b = s["sub_bands"]
    full = [k for k in b if k.endswith("Hz")][0]
    keys = [k for k in b if k.endswith("Hz") and k != full]
    x = np.arange(len(keys))
    xl = [k[:-2] for k in keys]
    fig, axs = plt.subplots(1, 3, figsize=(RS.FULL_W, 2.9))
    for j, (t, unit) in enumerate((("mag", "fT/$\\sqrt{Hz}$"), ("grad", "fT/(cm $\\sqrt{Hz}$)"))):
        ax = axs[j]
        for lab, col, ls, name in (("empty_room_measured", C_MEAS, "o-", "empty room, measured"),
                                   ("empty_room_model", C_MEAS, "--", "empty room, model"),
                                   ("brain_measured", C_MODEL, "s-", "brain noise, measured"),
                                   ("brain_model", C_MODEL, "--", "brain noise, model")):
            ax.semilogy(x, [b[k][t][lab]["median_asd"] for k in keys], ls, color=col, ms=3, label=name)
        ax.set_xticks(x)
        ax.set_xticklabels(xl, fontsize=6.5, rotation=40)
        ax.set_xlabel("band (Hz)")
        ax.set_ylabel(f"median amplitude ({unit})")
        _nice_log_y(ax)
        ax.set_title(f"{'ab'[j]}  {TYPE_LABEL[t].capitalize()}" + (" (brain: predicted)" if t == "mag" else " (brain: fitted per band)"),
                     loc="left", fontsize=8)
        if j == 0:
            ax.legend(fontsize=5.5, loc="lower left")
    ax = axs[2]
    y = [b[k]["mag"]["brain_amp_ratio_model_over_measured"] for k in keys]
    if "null95" in b[keys[0]]:
        ax.fill_between(x, [b[k]["null95"]["mag.brain_amp_ratio_model_over_measured"][0] for k in keys],
                        [b[k]["null95"]["mag.brain_amp_ratio_model_over_measured"][2] for k in keys], color=C_NULL, alpha=0.35, lw=0,
                        label="exact model, as it would be measured")
    ax.plot(x, y, "o-", color=C_MEAS, ms=3, label="measured")
    ax.axhline(b[full]["mag"]["brain_amp_ratio_model_over_measured"], color=C_MODEL, lw=0.8, ls="--", label=f"whole band ({full[:-2]} Hz)")
    ax.axhline(1.0, color="0.6", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(xl, fontsize=6.5, rotation=40)
    ax.set_xlabel("band (Hz)")
    ax.set_ylabel("model / measured (median RMS)")
    ax.set_title("c  Magnetometer brain noise", loc="left", fontsize=8)
    ax.legend(fontsize=5.5, loc="lower left")
    fig.suptitle("Noise per frequency band; the model's cortical background is rescaled in each band on the gradiometers", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def figure_detectability(s, path):
    det = s["detectability"]["comparisons"]
    sm = s["detectability"]["surrogate_medians"]
    imp = s["opm_implication"]
    n_t = det["measured_corrected_over_model/combined"]["n"]
    fig = plt.figure(figsize=(RS.FULL_W, 6.4 if imp else 3.4))
    top = fig.add_gridspec(1, 2, left=0.17, right=0.97, top=0.9, bottom=0.6 if imp else 0.3, wspace=0.6)
    axs = [fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1])]
    if imp:
        axs.append(fig.add_subplot(fig.add_gridspec(1, 1, left=0.36, right=0.97, top=0.43, bottom=0.08)[0, 0]))
    ax = axs[0]
    series = [("measured covariance", "measured_corrected_over_model", C_MEAS, "o"),
              ("measured, heart's field removed", "without_heart_corrected_over_model", "#D55E00", "D"),
              ("measured empty room + model brain", "hybrid_corrected_over_model", C_MODEL, "s")]
    labels = []
    for j, cs in enumerate(("combined", "mag", "grad")):
        for k, (lab, key, col, mk) in enumerate(series):
            c = det.get(f"{key}/{cs}")
            if c is None:
                continue
            y = j + 0.2 * (k - 1)
            ax.errorbar([c["ratio"]], [y], xerr=[[c["ratio"] - c["ci95"][0]], [c["ci95"][1] - c["ratio"]]], fmt=mk, color=col, ms=3.5, capsize=2,
                        label=lab if j == 0 else None)
        nl = sm[f"null_corrected/{cs}"]
        ax.plot([nl[0], nl[2]], [j + 0.38, j + 0.38], color=C_NULL, lw=2.5, label="exact model, as it would be measured" if j == 0 else None)
        labels.append({"combined": "all 305 channels", "mag": "102 magnetometers", "grad": "203 gradiometers"}[cs])
    ax.axvline(1.0, color="0.6", lw=0.7)
    ax.set_yticks(range(3))
    ax.set_yticklabels(labels)
    ax.set_ylim(2.65, -0.45)
    ax.set_xlabel(f"measured / model (median over {n_t:,} targets)")
    ax.set_title("a  Neuromag detectability with its measured noise", loc="left", fontsize=8)
    ax.legend(fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    # (b) by lobe
    ax = axs[1]
    bl = s["detectability"]["by_lobe_combined"]
    lobes = [lb for lb in ("frontal", "parietal", "temporal", "occipital", "insula", "cingulate", "other") if lb in bl]
    for j, lb in enumerate(lobes):
        for k, (key, col, mk) in enumerate((("measured_corrected", C_MEAS, "o"), ("without_heart_corrected", "#D55E00", "D"))):
            c = bl[lb].get(key)
            if c is None:
                continue
            ax.errorbar([c["ratio"]], [j + 0.18 * (k - 0.5)], xerr=[[c["ratio"] - c["ci95"][0]], [c["ci95"][1] - c["ratio"]]], fmt=mk,
                        color=col, ms=3.2, capsize=1.5)
    ax.axvline(1.0, color="0.6", lw=0.7)
    ax.set_yticks(range(len(lobes)))
    ax.set_yticklabels([f"{LOBE_LABEL.get(lb, lb)} ({bl[lb]['n_targets']})" for lb in lobes], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("measured / model, all 305 channels")
    ax.set_title("b  By lobe (targets)", loc="left", fontsize=8)
    if imp:
        ax = axs[2]
        r = imp["ratios"]
        rows = [("model_306", "both modelled (as published)"),
                ("S1_same_relative_change", "same relative change in both (S1)"),
                ("S2_same_variance_excess_both_modelled", "both modelled, cortical background x k (S2)"),
                ("S5_neuromag_measured_empty_room", "Neuromag: measured empty room + model brain (S5)"),
                ("S4_neuromag_measured_opm_as_modelled", "Neuromag measured, OPM as modelled (S4)"),
                ("S3_neuromag_measured_opm_variance_excess", "Neuromag measured, OPM background x k (S3)")]
        for j, (key, lab) in enumerate(rows):
            for k, (a, col, mk) in enumerate((("opm_dense", "#0072B2", "o"), ("opm_matched", "#56B4E9", "s"))):
                c = r[f"{a}/combined/{key}"]
                ax.errorbar([c["ratio"]], [j + 0.18 * (k - 0.5)], xerr=[[c["ratio"] - c["ci95"][0]], [c["ci95"][1] - c["ratio"]]], fmt=mk,
                            color=col, ms=3.5, capsize=2,
                            label=({"opm_dense": "dense OPM array (208 sites)", "opm_matched": "site-matched OPM array (98 sites)"}[a]
                                   if j == 0 else None))
        ax.axvline(1.0, color="0.6", lw=0.7)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for _, lab in rows], fontsize=6.3)
        ax.invert_yaxis()
        ax.set_xlabel("OPM / Neuromag detectability (median over targets, 95 % parcel interval; all 305 Neuromag channels)")
        ax.set_title(f"c  Implied OPM / Neuromag ratio under stated assumptions (k = {imp['checks']['magnetometer_variance_excess_k']:.2f}, "
                     "the measured / modelled magnetometer brain-noise variance)", loc="left", fontsize=8)
        ax.legend(fontsize=6, loc="lower right")
    fig.suptitle("Known-topography detectability with the measured Neuromag noise (sensor + brain + room noise, 1-40 Hz)", fontsize=9)
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
