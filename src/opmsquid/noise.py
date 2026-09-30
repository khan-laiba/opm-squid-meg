"""Noise spectra, analysis filters and PSD-to-variance conversion.

Conventions
-----------
* Spectra are one-sided: a white-noise amplitude spectral density ``asd`` (e.g. T/sqrt(Hz) or
  T/(m sqrt(Hz))) has one-sided PSD asd^2 on [0, fs/2], so unfiltered sampled white noise has
  variance asd^2 * fs / 2.
* The variance after the analysis filter is the integral of the one-sided PSD times the
  squared magnitude of the *composite* response actually applied. A zero-phase
  (forward-backward) IIR filter applies |H|^2 in amplitude, i.e. |H|^4 in power; this is
  handled by ``AnalysisFilter.power_response``. ``tests/test_noise.py`` checks the analytic
  variance against simulated, filtered noise.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class AnalysisFilter:
    """Butterworth band-pass (or low-pass if l_freq is None) used identically on signal and
    noise. ``zero_phase`` applies it forward and backward (sosfiltfilt)."""

    fs: float
    l_freq: float | None
    h_freq: float
    order: int = 4
    zero_phase: bool = True

    @property
    def sos(self) -> np.ndarray:
        if self.l_freq is None:
            return signal.butter(self.order, self.h_freq, btype="lowpass", fs=self.fs, output="sos")
        return signal.butter(self.order, [self.l_freq, self.h_freq], btype="bandpass", fs=self.fs, output="sos")

    def apply(self, x: np.ndarray, axis: int = -1) -> np.ndarray:
        if self.zero_phase:
            return signal.sosfiltfilt(self.sos, x, axis=axis)
        return signal.sosfilt(self.sos, x, axis=axis)

    def power_response(self, freqs: np.ndarray) -> np.ndarray:
        """|H_composite(f)|^2 of the operation performed by ``apply``."""
        _, h = signal.sosfreqz(self.sos, worN=np.asarray(freqs, dtype=float), fs=self.fs)
        p = np.abs(h) ** 2
        return p**2 if self.zero_phase else p

    def enbw(self, n_freqs: int = 200_001) -> float:
        """Equivalent noise bandwidth [Hz]: integral of the power response over [0, fs/2]."""
        f = np.linspace(0.0, self.fs / 2.0, n_freqs)
        return float(np.trapezoid(self.power_response(f), f))

    def variance(self, psd, n_freqs: int = 200_001) -> float:
        """Variance of noise with one-sided PSD ``psd`` (callable of frequency, or a constant)
        after this filter."""
        f = np.linspace(0.0, self.fs / 2.0, n_freqs)
        p = psd(f) if callable(psd) else np.full_like(f, float(psd))
        return float(np.trapezoid(p * self.power_response(f), f))


def white_variance(asd: float, filt: AnalysisFilter) -> float:
    """Variance of white noise with amplitude spectral density ``asd`` after ``filt``."""
    return asd**2 * filt.enbw()


def white_noise(asd: float, fs: float, shape, rng: np.random.Generator) -> np.ndarray:
    """Sampled white noise with one-sided amplitude spectral density ``asd``."""
    return rng.standard_normal(shape) * asd * np.sqrt(fs / 2.0)
