"""Sensor-level SNR metrics.

Notation: ``s`` is a declared noise-free target topography, shape (n_channels,) or
(n_channels, n_sources), in each channel's native unit (T for magnetometers, T/m for planar
gradiometers). ``cov`` is the matching noise covariance (native units squared), or a vector
of per-channel variances where only the diagonal is used.

Every metric here is dimensionless and invariant to a consistent per-channel change of units
(s_i -> a_i s_i, C -> A C A); tests/test_metrics.py checks this. Unscaled T and T/m values are
never averaged together: each metric works on noise-normalised quantities.

General metrics (GOAL G2)
    peak_channel_snr      max_i |s_i| / sqrt(C_ii)
    mean_power_snr_db     10 log10[(1/N) sum_i s_i^2 / C_ii]
    detectability         sqrt(s^T C^+ s) with a rank-aware whitener (known-topography,
                          oracle-covariance matched-filter SNR; not an event-detection rate or
                          a localization accuracy)

Paper-specific metrics live in the modules implementing each paper (e.g. the Goldenholz
Eq. 1 and Hunold Hilbert-envelope SNR), so that their definitions stay separate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _topographies(s):
    s = np.asarray(s, dtype=float)
    return (s[:, None], True) if s.ndim == 1 else (s, False)


def _variances(cov) -> np.ndarray:
    cov = np.asarray(cov, dtype=float)
    var = np.diag(cov) if cov.ndim == 2 else cov
    if np.any(var <= 0):
        raise ValueError("noise variances must be positive")
    return var


def peak_channel_snr(s, cov) -> np.ndarray | float:
    """Peak-channel RMS SNR, max_i |s_i| / sqrt(C_ii), per topography."""
    s2, single = _topographies(s)
    snr = np.max(np.abs(s2) / np.sqrt(_variances(cov))[:, None], axis=0)
    return float(snr[0]) if single else snr


def mean_power_snr_db(s, cov) -> np.ndarray | float:
    """Mean-power SNR in dB, 10 log10[(1/N) sum_i s_i^2 / C_ii], per topography."""
    s2, single = _topographies(s)
    ratio = np.mean(s2**2 / _variances(cov)[:, None], axis=0)
    with np.errstate(divide="ignore"):
        db = 10.0 * np.log10(ratio)
    return float(db[0]) if single else db


@dataclass
class Whitener:
    """Rank-aware whitening operator ``matrix`` (rank x n_channels): W C W^T = I_rank."""

    matrix: np.ndarray
    rank: int
    eigenvalues: np.ndarray  # eigenvalues of the noise-normalised covariance, descending
    scale: np.ndarray  # 1 / sqrt(C_ii) applied before the eigen-decomposition

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x


def whitener(cov, rel_tol: float = 1e-10, rank: int | None = None) -> Whitener:
    """Rank-aware whitener without explicit inversion of an ill-conditioned matrix.

    The channels are first scaled by 1/sqrt(C_ii), which removes units and the large
    magnetometer/gradiometer scale differences; the resulting correlation-like matrix
    R = D C D is eigen-decomposed, eigenvalues below ``rel_tol`` times the largest (or beyond
    ``rank``) are discarded, and W = Lambda_r^(-1/2) U_r^T D. Signal components in the
    discarded (noise-free) subspace are dropped, as in MNE's whitening, rather than producing
    unbounded SNR.
    """
    cov = np.asarray(cov, dtype=float)
    cov = 0.5 * (cov + cov.T)
    d = 1.0 / np.sqrt(_variances(cov))
    rmat = cov * np.outer(d, d)
    evals, evecs = np.linalg.eigh(rmat)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    keep = int(np.sum(evals > rel_tol * evals[0])) if rank is None else int(rank)
    if not 0 < keep <= len(evals) or evals[keep - 1] <= 0:
        raise ValueError(f"invalid rank {keep} for a covariance with eigenvalues {evals[:3]}...")
    matrix = (evecs[:, :keep] / np.sqrt(evals[:keep])).T * d[None, :]
    return Whitener(matrix=matrix, rank=keep, eigenvalues=evals, scale=d)


def detectability(s, cov=None, w: Whitener | None = None) -> np.ndarray | float:
    """Known-topography detectability sqrt(s^T C^+ s) = ||W s||, per topography.

    This is the SNR of an ideal matched filter that knows the topography and the (oracle or
    estimated) noise covariance; it is not a detection rate and not a localization accuracy.
    """
    if w is None:
        if cov is None:
            raise ValueError("pass cov or a precomputed whitener")
        w = whitener(cov)
    s2, single = _topographies(s)
    d = np.linalg.norm(w.apply(s2), axis=0)
    return float(d[0]) if single else d


def plugin_detectability(s, cov_true, cov_est, rel_tol: float = 1e-10) -> np.ndarray | float:
    """Output SNR of the matched filter h = C_est^+ s built from an estimated covariance and
    applied to data whose noise has the true covariance: (h^T s) / sqrt(h^T C_true h). Equals
    ``detectability(s, cov_true)`` when C_est = C_true and is never larger (Cauchy-Schwarz)."""
    s2, single = _topographies(s)
    w = whitener(cov_est, rel_tol=rel_tol)
    h = w.matrix.T @ (w.matrix @ s2)
    c = 0.5 * (np.asarray(cov_true, float) + np.asarray(cov_true, float).T)
    num = np.sum(h * s2, axis=0)
    den = np.sqrt(np.maximum(np.sum(h * (c @ h), axis=0), 0.0))
    out = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return float(out[0]) if single else out


def empirical_covariance(x: np.ndarray) -> np.ndarray:
    """Covariance of zero-mean data x (n_channels, n_times), normalised by n_times."""
    x = np.asarray(x, dtype=float)
    return x @ x.T / x.shape[1]


def ledoit_wolf_covariance(x: np.ndarray) -> tuple[np.ndarray, float]:
    """Ledoit-Wolf shrinkage of the covariance of zero-mean data, computed on channels scaled
    to unit variance (so the shrinkage target is unit-free) and scaled back. Returns
    (covariance, shrinkage coefficient)."""
    x = np.asarray(x, dtype=float)
    std = np.sqrt(np.mean(x**2, axis=1))
    z = x / std[:, None]
    n_ch, n_t = z.shape
    s = z @ z.T / n_t
    mu = np.trace(s) / n_ch
    delta2 = np.sum((s - mu * np.eye(n_ch)) ** 2)
    # sum_t ||z_t z_t^T - S||_F^2 = sum_t ||z_t||^4 - n_t ||S||_F^2 (no n_ch^2 x n_t array needed)
    beta2 = min((np.sum(np.sum(z**2, axis=0) ** 2) - n_t * np.sum(s**2)) / n_t**2, delta2)
    shrink = beta2 / delta2 if delta2 > 0 else 1.0
    s_shrunk = shrink * mu * np.eye(n_ch) + (1.0 - shrink) * s
    return s_shrunk * np.outer(std, std), float(shrink)
