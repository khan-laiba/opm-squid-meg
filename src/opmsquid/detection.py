"""IED detection (G4, NEW): a known-source oracle and a practical scanning detector.

Both work on spatially whitened data (whitener from a covariance *estimated on independent null
data*) and temporal matched filtering with the analysis-filtered IED template.

Oracle: knows the event topography s, its waveform (template) and its peak time t0. Statistic
    z = u^T (W y * h)(t0) / sigma_null,   u = W s / ||W s||,
one test per event, one-sided (known polarity). Its null distribution is sampled at random times
of null data, so the threshold is a per-trial false-positive probability, not an event rate.

Practical: knows neither time nor source. For every sample, the statistic is the maximum over a
template bank (waveform stretches) and a dictionary of candidate cortical sources (fixed normal
orientation, a grid distinct from the true sources) of |d_j^T (W y * h_k)(t)| / sigma_jk. Local
maxima at least ``refractory`` apart are events; the threshold for a target false-event rate
(events per minute) is set on calibration null data and frozen, and verified on held-out null data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from . import metrics


def matched_filter(y: np.ndarray, template: np.ndarray, peak_index: int) -> np.ndarray:
    """out[..., t] = sum_k y[..., t - peak_index + k] template[k]: the template's peak aligned with
    sample t (zero outside the data)."""
    full = signal.fftconvolve(y, template[::-1][None, :] if y.ndim == 2 else template[::-1], mode="full", axes=-1)
    start = len(template) - 1 - peak_index
    return full[..., start:start + y.shape[-1]]


def whitener_from_null(null: np.ndarray) -> metrics.Whitener:
    """Rank-aware whitener from a Ledoit-Wolf covariance estimate of null data (n_ch, n_t)."""
    cov, _ = metrics.ledoit_wolf_covariance(null - null.mean(axis=1, keepdims=True))
    return metrics.whitener(cov)


@dataclass
class Oracle:
    w: metrics.Whitener
    templates: dict  # stretch -> (template, peak index)
    sigma: dict  # (source key, stretch) -> null SD

    @staticmethod
    def calibrate(null: np.ndarray, w: metrics.Whitener, templates: dict, topographies: dict, n_samples: int = 20000,
                  rng: np.random.Generator | None = None) -> tuple["Oracle", dict]:
        """Null SD of each (topography, template) statistic, and its null samples (normalised)."""
        rng = rng or np.random.default_rng(0)
        yw = w.apply(null)
        idx = rng.choice(np.arange(200, null.shape[1] - 200), min(n_samples, null.shape[1] - 400), replace=False)
        sigma, samples = {}, {}
        for k, (tpl, pk) in templates.items():
            yf = matched_filter(yw, tpl, pk)[:, idx]
            for key, s in topographies.items():
                u = w.apply(s)
                u = u / np.linalg.norm(u)
                z = u @ yf
                sigma[(key, k)] = float(z.std())
                samples[(key, k)] = z / sigma[(key, k)]
        return Oracle(w, templates, sigma), samples

    def statistic(self, data: np.ndarray, key, stretch, topography: np.ndarray, t0: int) -> float:
        tpl, pk = self.templates[stretch]
        u = self.w.apply(topography)
        u = u / np.linalg.norm(u)
        lo, hi = max(t0 - len(tpl), 0), min(t0 + len(tpl), data.shape[1])
        seg = self.w.apply(data[:, lo:hi])
        return float((u @ matched_filter(seg, tpl, pk))[t0 - lo] / self.sigma[(key, stretch)])


@dataclass
class ScanDetector:
    w: metrics.Whitener
    templates: dict  # stretch -> (template, peak index)
    dictionary: np.ndarray  # (r, J) unit-norm whitened candidate topographies
    sigma: np.ndarray  # (n_templates, J) null SD of each output
    refractory: int  # samples

    @staticmethod
    def build(null_cal: np.ndarray, w: metrics.Whitener, templates: dict, candidates: np.ndarray, refractory: int,
              chunk: int = 4000) -> "ScanDetector":
        d = w.apply(candidates)
        d = d / np.linalg.norm(d, axis=0, keepdims=True)
        det = ScanDetector(w, templates, d, np.ones((len(templates), d.shape[1])), refractory)
        acc = np.zeros((len(templates), d.shape[1]))
        n = 0
        for part in det._outputs(null_cal, chunk):
            acc += np.sum(part**2, axis=2)
            n += part.shape[2]
        det.sigma = np.sqrt(acc / n)
        return det

    def _outputs(self, data: np.ndarray, chunk: int):
        """Yield normalised outputs (n_templates, J, n) chunk by chunk (overlapping margins so the
        convolution is exact inside each chunk)."""
        yw = self.w.apply(data)
        margin = max(len(t) for t, _ in self.templates.values())
        for start in range(0, yw.shape[1], chunk):
            a, b = max(start - margin, 0), min(start + chunk + margin, yw.shape[1])
            parts = []
            for i, (tpl, pk) in enumerate(self.templates.values()):
                yf = matched_filter(yw[:, a:b], tpl, pk)[:, start - a:start - a + min(chunk, yw.shape[1] - start)]
                parts.append((self.dictionary.T @ yf) / self.sigma[i][:, None])
            yield np.stack(parts)

    def statistic(self, data: np.ndarray, chunk: int = 4000) -> tuple[np.ndarray, np.ndarray]:
        """Per-sample max |normalised output| over templates and candidates, and the argmax
        candidate index."""
        stat, arg = [], []
        for part in self._outputs(data, chunk):
            a = np.abs(part).max(axis=0)  # over templates
            stat.append(a.max(axis=0))
            arg.append(a.argmax(axis=0))
        return np.concatenate(stat), np.concatenate(arg)

    def events(self, stat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Local maxima at least ``refractory`` samples apart: (sample indices, heights)."""
        peaks, props = signal.find_peaks(stat, height=0.0, distance=self.refractory)
        return peaks, props["peak_heights"]


def threshold_for_rate(heights: np.ndarray, minutes: float, rate_per_min: float) -> float:
    """Smallest threshold with at most rate x minutes null events above it."""
    h = np.sort(np.asarray(heights))[::-1]
    k = int(np.floor(rate_per_min * minutes))
    return float(h[k]) if k < len(h) else float(h[-1])


def match_events(peaks: np.ndarray, heights: np.ndarray, truth: np.ndarray, tol: int, threshold: float) -> np.ndarray:
    """For each true peak sample: detected if an event above ``threshold`` lies within +/- tol."""
    sel = peaks[heights > threshold]
    if len(sel) == 0:
        return np.zeros(len(truth), bool)
    pos = np.searchsorted(sel, truth)
    near = np.full(len(truth), np.inf)
    for off in (-1, 0):
        j = np.clip(pos + off, 0, len(sel) - 1)
        near = np.minimum(near, np.abs(sel[j] - truth))
    return near <= tol
