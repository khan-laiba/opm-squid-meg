"""Cortical background (brain-noise) sources, common to all arrays (G2; background for G1B/G4).

Model: a random cortical current density J (A m per m^2 of cortex), normal to the cortex, with
spatial covariance s^2 rho(|x - x'|). The moment of a source representing area a_i is
m_i = integral of J over its patch, so Cov(m_i, m_j) = s^2 a_i a_j rho_ij / A_c, with A_c the
correlation area of rho:

* independent baseline: rho -> delta, Cov = s^2 a_i delta_ij. Moment variance is proportional
  to area, so refining the source mesh does not add background power.
* correlated extension (bounded): rho = exp(-d / lambda) with Euclidean 3-D distance d
  (positive definite), A_c = 2 pi lambda^2 (planar value), which tends to the independent model
  as lambda -> 0. Opposite walls of a sulcus closer than lambda are correlated (documented).

The source scale s is not a known constant; ``calibrate`` fixes it so that the background
variance of a chosen channel set (default: Neuromag gradiometers, which are least affected by
cardiac, ocular and environmental fields) matches a measured brain-noise level. The same s is
then used for every array.
"""
from __future__ import annotations

import numpy as np


def moment_covariance(areas: np.ndarray, rr: np.ndarray | None = None, length: float | None = None):
    """Unit-scale (s = 1) moment covariance: a vector of variances (independent model) or a
    dense matrix (correlated model with correlation length ``length`` [m])."""
    areas = np.asarray(areas, float)
    if not length:
        return areas
    d = np.linalg.norm(rr[:, None, :] - rr[None, :, :], axis=-1)
    return np.outer(areas, areas) * np.exp(-d / length) / (2.0 * np.pi * length**2)


def sensor_covariance(gain: np.ndarray, moment_cov) -> np.ndarray:
    """Background sensor covariance L Sigma L^T for gain (n_channels, n_sources)."""
    mc = np.asarray(moment_cov)
    if mc.ndim == 1:
        return (gain * mc[None, :]) @ gain.T
    return gain @ mc @ gain.T


def joint_factor(gains: list[np.ndarray], moment_cov, rel_tol: float = 1e-12) -> tuple[np.ndarray, list[slice]]:
    """Factor F of the joint background covariance of several arrays, [L_1; L_2; ...] Sigma
    [...]^T = F F^T, so that y = F w with white w reproduces *one* common source realization seen
    by every array (identical background across arrays, exact cross-array covariance)."""
    stacked = np.vstack(gains)
    cov = sensor_covariance(stacked, moment_cov)
    d = np.sqrt(np.diag(cov))
    ev, u = np.linalg.eigh(cov / np.outer(d, d))
    keep = ev > rel_tol * ev.max()
    f = (u[:, keep] * np.sqrt(ev[keep])) * d[:, None]
    bounds = np.cumsum([0] + [g.shape[0] for g in gains])
    return f, [slice(a, b) for a, b in zip(bounds[:-1], bounds[1:])]


def calibrate(unit_cov: np.ndarray, channels, target_variance: float, statistic=np.median) -> float:
    """Scale s^2 such that statistic(diag(s^2 unit_cov)[channels]) equals target_variance."""
    return float(target_variance / statistic(np.diag(unit_cov)[channels]))


def pink_noise(n_series: int, n_times: int, fs: float, rng: np.random.Generator, exponent: float = 1.0,
               f_min: float = 0.5) -> np.ndarray:
    """Unit-variance Gaussian noise with a 1/f^exponent power spectrum above f_min (flat below)."""
    f = np.fft.rfftfreq(n_times, 1.0 / fs)
    shape = np.maximum(f, f_min) ** (-exponent / 2.0)
    shape[0] = 0.0
    spec = (rng.standard_normal((n_series, f.size)) + 1j * rng.standard_normal((n_series, f.size))) * shape
    x = np.fft.irfft(spec, n=n_times, axis=1)
    return x / x.std(axis=1, keepdims=True)
