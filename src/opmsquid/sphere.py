"""Analytical spherical-head benchmark of Jas et al. (2026) and independent numerical checks.

Geometry (Jas et al. 2026, Sec. 2.1): sphere centred at the origin; a current dipole of moment
q at r_Q on the z axis, oriented tangentially (+y); ideal point magnetometers measuring the
radial field B_r on a concentric sensor sphere of radius s (s = h for on-scalp OPMs, s = h + xi
for off-scalp SQUIDs). Source depth is measured from the scalp: d = h - r_Q.

In a spherically symmetric conductor the volume currents add nothing to B_r, so B_r equals the
radial component of the primary (Biot-Savart) field. ``bmax_radial`` is the closed-form peak
(paper Eq. 1); ``numerical_peak_radial`` maximises the radial component of the full Sarvas
(1987) field over the whole sensor sphere and is independent of the closed form.

A dipole at the sphere centre (r_Q = 0, d = h) is magnetically silent (B_r = 0 everywhere); it
is excluded from every ratio and crossover calculation.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, minimize

MU0_OVER_4PI = 1e-7  # T m / A


def bmax_radial(r_q, s, q: float = 30e-9) -> np.ndarray:
    """Peak radial field [T] of a tangential dipole at radius r_q on a sensor sphere of radius
    s (Jas et al. 2026, Eq. 1): mu0/4pi q r_q sin(theta0) / (s^2 + r_q^2 - 2 s r_q cos(theta0))^1.5
    with cos(theta0) = -gamma + sqrt(gamma^2 + 3), gamma = (s^2 + r_q^2) / (2 s r_q), evaluated in
    the cancellation-free form cos(theta0) = 3 / (gamma + sqrt(gamma^2 + 3)). Zero at r_q = 0."""
    r_q = np.atleast_1d(np.asarray(r_q, dtype=float))
    out = np.zeros_like(r_q)
    m = r_q > 0
    rq = r_q[m]
    gamma = (s**2 + rq**2) / (2.0 * s * rq)
    cos0 = 3.0 / (gamma + np.sqrt(gamma**2 + 3.0))
    out[m] = MU0_OVER_4PI * q * rq * np.sqrt(1.0 - cos0**2) / (s**2 + rq**2 - 2.0 * s * rq * cos0) ** 1.5
    return out


def sarvas_field(r: np.ndarray, r_q: np.ndarray, q_vec: np.ndarray) -> np.ndarray:
    """Magnetic field [T] at points r (n, 3) outside a spherically symmetric conductor centred
    at the origin, of a current dipole q_vec [A m] at r_q (Sarvas 1987):
    B = mu0/(4 pi F^2) (F q x r_q - (q x r_q . r) grad F),
    F = a (r a + r^2 - r_q . r), a = r - r_q, grad F = (a^2/r + a.r/a + 2a + 2r) r - (a + 2r + a.r/a) r_q."""
    r = np.atleast_2d(np.asarray(r, dtype=float))
    r_q = np.asarray(r_q, dtype=float)
    q_vec = np.asarray(q_vec, dtype=float)
    a_vec = r - r_q
    a = np.linalg.norm(a_vec, axis=1)
    rn = np.linalg.norm(r, axis=1)
    ar = np.sum(a_vec * r, axis=1)
    f = a * (rn * a + rn**2 - r @ r_q)
    c1 = a**2 / rn + ar / a + 2.0 * a + 2.0 * rn
    c2 = a + 2.0 * rn + ar / a
    grad_f = c1[:, None] * r - c2[:, None] * r_q[None, :]
    qxr = np.cross(q_vec, r_q)
    return MU0_OVER_4PI / f[:, None] ** 2 * (f[:, None] * qxr[None, :] - (r @ qxr)[:, None] * grad_f)


def _radial_component(theta, phi, r_q, s, q):
    u = np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], axis=-1)
    b = sarvas_field(s * u.reshape(-1, 3), np.array([0.0, 0.0, r_q]), np.array([0.0, q, 0.0]))
    return np.sum(b * u.reshape(-1, 3), axis=1).reshape(np.shape(theta))


def numerical_peak_radial(r_q: float, s: float, q: float = 30e-9, n_theta: int = 181,
                          n_phi: int = 360) -> tuple[float, float, float]:
    """Peak |B_r| [T] over the whole sensor sphere of radius s, from the Sarvas field on a
    (theta, phi) grid refined by Nelder-Mead. Returns (peak, theta, phi) of the maximum."""
    if r_q <= 0:
        return 0.0, np.nan, np.nan
    th, ph = np.meshgrid(np.linspace(0.0, np.pi, n_theta), np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False),
                         indexing="ij")
    val = np.abs(_radial_component(th, ph, r_q, s, q))
    i = np.unravel_index(np.argmax(val), val.shape)
    res = minimize(lambda x: -abs(float(_radial_component(np.array(x[0]), np.array(x[1]), r_q, s, q))),
                   x0=[th[i], ph[i]], method="Nelder-Mead",
                   options={"xatol": 1e-12, "fatol": 1e-30, "maxiter": 4000})
    return float(-res.fun), float(res.x[0]), float(res.x[1] % (2 * np.pi))


def signal_ratio(d, h: float = 0.095, xi_opm: float = 0.0, xi_squid: float = 0.018) -> np.ndarray:
    """B_OPM / B_SQUID of the peak radial fields at depth d below the scalp (NaN at the silent
    sphere centre and outside the head)."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    r_q = h - d
    out = np.full_like(d, np.nan)
    ok = (r_q > 0) & (r_q < h)
    out[ok] = bmax_radial(r_q[ok], h + xi_opm) / bmax_radial(r_q[ok], h + xi_squid)
    return out


def equal_snr_depth(eta: float, h: float = 0.095, b: float = 0.080, xi_opm: float = 0.0,
                    xi_squid: float = 0.018) -> float:
    """Equal-SNR depth d_eq [m] (Jas et al. 2026, Eq. 3): SNR_OPM = SNR_SQUID, i.e.
    B_OPM / B_SQUID = eta with sigma_OPM = eta sigma_SQUID, for sources inside the brain
    (h - b <= d < h). The ratio decreases monotonically with depth, from its value at the brain
    surface to ((h + xi_squid) / (h + xi_opm))^3 at the (silent) centre, which is excluded.
    Returns NaN if eta is outside that range (no crossing in the brain)."""
    f = lambda d: float(signal_ratio(d, h, xi_opm, xi_squid)[0]) - eta  # noqa: E731
    d_lo, d_hi = h - b, h * (1.0 - 1e-9)
    if f(d_lo) < 0 or f(d_hi) > 0:
        return np.nan
    return float(brentq(f, d_lo, d_hi, xtol=1e-12))


def centre_limit_ratio(h: float = 0.095, xi_opm: float = 0.0, xi_squid: float = 0.018) -> float:
    """Limit of B_OPM / B_SQUID as the source approaches the (silent) sphere centre."""
    return ((h + xi_squid) / (h + xi_opm)) ** 3
