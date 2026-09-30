"""IED-like events in simulated multichannel recordings (G4, NEW).

Time-domain counterpart of the G2 noise model (``opmsquid.noisemodel``): for every array the
recording is

    y_m(t) = L_target,m q_ied(t) + L_bg,m q_bg(t) + E_m e(t) + eps_m(t)

with *one* realization of q_bg (cortical background) and e (room field) shared by all arrays.

* Background: the G2 grid sources, independent 1/f moments, band-limited by the analysis filter
  and scaled to variance brain_scale x area in the band (so the band covariance equals the G2
  model's brain term).
* Room field: 8 external coefficients synthesised with the cross-spectral density of the
  empty-room fit (narrow peaks near 12 and 30 Hz in the sample room), scaled to the fitted
  band covariance.
* Intrinsic: white noise with each channel's ASD, through the same filter.
Everything is simulated at the recording rate, filtered, then decimated (the analysis band is
far below the new Nyquist frequency).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from . import hunold


def ied_waveform(fs: float, stretch: float = 1.0) -> np.ndarray:
    """Spike-wave complex (Hunold et al. 2016 Fig. 2b, digitised), peak 1, time-stretched by
    ``stretch`` (0.75 = sharper and shorter, 1.5 = broader)."""
    t_ms, v = np.array(hunold.SPIKE_POINTS, float).T
    t = np.arange(0.0, t_ms[-1] * stretch / 1e3, 1.0 / fs) * 1e3 / stretch
    from scipy.interpolate import PchipInterpolator

    w = PchipInterpolator(t_ms, v)(np.clip(t, 0, t_ms[-1]))
    return w / w.max()


def csd_matrix(x: np.ndarray, fs: float, nperseg: int = 4096):
    """One-sided cross-spectral density matrices S(f) (n_f, n, n) of the rows of x (Welch)."""
    n = x.shape[0]
    f, _ = signal.welch(x[0], fs=fs, nperseg=nperseg)
    s = np.zeros((len(f), n, n), complex)
    for i in range(n):
        for j in range(i, n):
            _, sij = signal.csd(x[i], x[j], fs=fs, nperseg=nperseg)
            s[:, i, j] = sij
            s[:, j, i] = np.conj(sij)
    return f, s


def synthesize_from_csd(f: np.ndarray, s: np.ndarray, n_times: int, fs: float, rng: np.random.Generator) -> np.ndarray:
    """Gaussian multichannel series of length ``n_times`` whose cross-spectra follow S(f)
    (interpolated to the FFT grid; Cholesky factor per frequency). The absolute level is set by
    the caller."""
    fk = np.fft.rfftfreq(n_times, 1.0 / fs)
    n = s.shape[1]
    re = np.stack([np.interp(fk, f, s[:, i, j].real) for i in range(n) for j in range(n)], axis=1).reshape(len(fk), n, n)
    im = np.stack([np.interp(fk, f, s[:, i, j].imag) for i in range(n) for j in range(n)], axis=1).reshape(len(fk), n, n)
    sk = re + 1j * im
    sk = 0.5 * (sk + np.conj(np.transpose(sk, (0, 2, 1))))
    w, v = np.linalg.eigh(sk)
    root = v * np.sqrt(np.clip(w, 0.0, None))[:, None, :]  # S = root root^H (robust to rank deficiency)
    z = (rng.standard_normal((len(fk), n)) + 1j * rng.standard_normal((len(fk), n))) / np.sqrt(2.0)
    spec = np.einsum("kij,kj->ki", root, z)
    spec[0] = 0.0
    return np.fft.irfft(spec, n=n_times, axis=0).T


@dataclass
class ArraySpec:
    grid_gain: np.ndarray  # (n_ch, n_grid)
    env_basis: np.ndarray  # (n_ch, 8)
    asd: np.ndarray  # (n_ch,) intrinsic white-noise ASD [native units / sqrt(Hz)]


class NoiseGenerator:
    """Common noise realizations for several arrays (see the module docstring)."""

    def __init__(self, arrays: dict[str, ArraySpec], grid_var: np.ndarray, env_coef: np.ndarray, env_cov: np.ndarray,
                 filt, decimate: int = 4, env_fs: float | None = None):
        self.arrays, self.grid_sd = arrays, np.sqrt(np.asarray(grid_var, float))
        self.filt, self.decimate = filt, int(decimate)
        self.fs = filt.fs
        self.fs_out = self.fs / self.decimate
        self.env_f, self.env_s = csd_matrix(env_coef, env_fs or self.fs)
        self.env_cov = env_cov

    def segment(self, seconds: float, rng: np.random.Generator, pad_s: float = 2.0) -> dict[str, np.ndarray]:
        """Filtered, decimated noise for every array: {name: (n_ch, n_out)}. Generated ``pad_s``
        longer on each side (filter transients) and cropped."""
        pad = int(pad_s * self.fs)
        n = int(round(seconds * self.fs)) + 2 * pad
        # background moments: 1/f, band-limited, exact band variance per source
        from .background import pink_noise

        q = self.filt.apply(pink_noise(len(self.grid_sd), n, self.fs, rng))
        q = q / q[:, pad:-pad].std(axis=1, keepdims=True) * self.grid_sd[:, None]
        # room-field coefficients with the empty-room cross-spectra (already band-limited by the fit's
        # filter), scaled to the fitted band covariance
        e = synthesize_from_csd(self.env_f, self.env_s, n, self.fs, rng)
        c = np.cov(e[:, pad:-pad])
        a = np.linalg.cholesky(self.env_cov + 1e-30 * np.eye(len(c)))
        b = np.linalg.cholesky(c + 1e-30 * np.eye(len(c)))
        e = a @ np.linalg.solve(b, e)
        out = {}
        for name, spec in self.arrays.items():
            eps = self.filt.apply(rng.standard_normal((len(spec.asd), n)) * (spec.asd * np.sqrt(self.fs / 2.0))[:, None])
            y = spec.grid_gain @ q + spec.env_basis @ e + eps
            out[name] = y[:, pad:-pad:self.decimate]
        return out


def filtered_template(filt, stretch: float = 1.0, decimate: int = 4, pre_s: float = 0.2, post_s: float = 0.3):
    """The IED waveform (peak 1 before filtering) passed through the analysis filter and decimated
    like the noise, with ``pre_s``/``post_s`` margins. Returns (template, index of the unfiltered
    peak within it)."""
    w = ied_waveform(filt.fs, stretch)
    pre, post = int(pre_s * filt.fs), int(post_s * filt.fs)
    x = np.zeros(pre + len(w) + post)
    x[pre:pre + len(w)] = w
    y = filt.apply(x)[::decimate]
    return y, int(round((pre + int(np.argmax(w))) / decimate))


def inject(data: np.ndarray, topography: np.ndarray, template: np.ndarray, peak_index: int, amplitude: float, peak_sample: int):
    """Add amplitude x topography x template to data (n_ch, n_t) so that the template's peak falls
    on ``peak_sample`` (in place; the part outside the data is dropped)."""
    start = peak_sample - peak_index
    a, b = max(start, 0), min(start + len(template), data.shape[1])
    if b > a:
        data[:, a:b] += amplitude * np.outer(topography, template[a - start:b - start])
