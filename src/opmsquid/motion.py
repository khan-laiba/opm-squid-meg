"""Head motion and OPM slippage (bounded secondary extension of G4).

Two mechanisms, kept apart:

A. Geometry. In the fixed helmet the head moves relative to the SQUIDs; a head-mounted OPM array
   moves with the head and changes its geometry only if it slips. ``mismatched_detectability``
   is the output SNR of the known-topography matched filter whose template was taken at the
   reference geometry while the data (signal and noise) come from the displaced one.

B. Room-field coupling. A sensor that moves rigidly with the head through a static residual field
   (uniform B0 plus a linear, curl- and divergence-free gradient G about ``x_ref``) sees a
   time-varying field: rotation turns its sensitive axis in the field, and rotation and
   translation move it through the gradient. SQUIDs are fixed in the room and see no such term.
   ``readings`` evaluates the exact rigid-motion response at each channel's integration points;
   ``jacobian`` is its linearisation in six motion parameters (rotation vector [rad] about the
   pivot ``c``, translation [m]). The head frame at rest is taken as the room frame; the field
   draws are isotropic, so the fixed rotation between the two does not matter.

First order, perfectly calibrated sensors: rotation in B0 and translation in G change the
reading by a uniform field in the head frame (removed by a homogeneous-field projection);
rotation in G adds a symmetric, traceless gradient (removed by the 8-term projection of
``environment.external_basis``). A uniform field stays uniform under any rotation, so the
homogeneous part is removed exactly, not only to first order. Calibration errors (sensitive-axis
tilt, gain) break these identities and leave residuals proportional to the field change.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import environment


def skew(v) -> np.ndarray:
    x, y, z = np.asarray(v, float)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotation(v) -> np.ndarray:
    """Rotation matrix of the rotation vector ``v`` (axis x angle [rad]; Rodrigues)."""
    v = np.asarray(v, float)
    th = float(np.linalg.norm(v))
    if th < 1e-12:
        return np.eye(3) + skew(v)
    k = skew(v / th)
    return np.eye(3) + np.sin(th) * k + (1.0 - np.cos(th)) * (k @ k)


@dataclass
class Sensors:
    """Integration points of every channel (head frame at rest): positions, orientations, weights."""

    rmag: np.ndarray  # (n_points, 3) m
    cosmag: np.ndarray  # (n_points, 3)
    w: np.ndarray  # (n_points,)
    ch: np.ndarray  # (n_points,) channel index
    n: int

    def sum(self, values: np.ndarray) -> np.ndarray:
        return np.bincount(self.ch, weights=self.w * values, minlength=self.n)


def sensors(info, coil_def=None) -> Sensors:
    defs = environment.coil_integration(info, coil_def)
    return Sensors(np.concatenate([d["rmag"] for d in defs]), np.concatenate([d["cosmag"] for d in defs]),
                   np.concatenate([d["w"] for d in defs]), np.concatenate([np.full(len(d["w"]), i) for i, d in enumerate(defs)]),
                   len(defs))


def external_basis(s: Sensors, r0) -> np.ndarray:
    """(n_channels, 8) response to the unit homogeneous-field and gradient components of
    ``environment.external_basis``, from these integration points (so that sensors with
    calibration errors follow the same convention)."""
    rel = s.rmag - np.asarray(r0, float)
    cols = [s.sum(s.cosmag[:, k]) for k in range(3)]
    cols += [s.sum(np.sum((rel @ g.T) * s.cosmag, axis=1)) for g in environment.GRADIENT_BASIS]
    return np.column_stack(cols)


def with_calibration_errors(s: Sensors, rng: np.random.Generator, axis_rms_deg: float, gain_rms: float) -> tuple[Sensors, np.ndarray]:
    """The actual sensors when the analyst's model is ``s``: every channel's sensitive axis tilted
    by a random rotation perpendicular to it (2-D Gaussian, RMS tilt ``axis_rms_deg``), and a gain
    1 + N(0, gain_rms). Returns (actual sensors, gains)."""
    cos = s.cosmag.copy()
    sd = np.radians(axis_rms_deg) / np.sqrt(2.0)
    for i in range(s.n):
        m = s.ch == i
        n = cos[m].mean(axis=0)
        n /= np.linalg.norm(n)
        a = np.cross(n, [1.0, 0.0, 0.0] if abs(n[0]) < 0.9 else [0.0, 1.0, 0.0])
        a /= np.linalg.norm(a)
        b = np.cross(n, a)
        r = rotation(sd * (rng.standard_normal() * a + rng.standard_normal() * b))
        cos[m] = cos[m] @ r.T
    gains = 1.0 + gain_rms * rng.standard_normal(s.n)
    return Sensors(s.rmag, cos, s.w, s.ch, s.n), gains


def random_field(rng: np.random.Generator, b0: float, g: float) -> tuple[np.ndarray, np.ndarray]:
    """Uniform field of magnitude ``b0`` [T] in a random direction, and a random symmetric,
    traceless gradient tensor with Frobenius norm ``g`` [T/m]."""
    u = rng.standard_normal(3)
    a = rng.standard_normal((3, 3))
    s = 0.5 * (a + a.T)
    s -= np.trace(s) / 3.0 * np.eye(3)
    return b0 * u / np.linalg.norm(u), g * s / np.linalg.norm(s)


def readings(s: Sensors, b0: np.ndarray, grad: np.ndarray, x_ref, rot: np.ndarray, trans, pivot,
             gains: np.ndarray | None = None) -> np.ndarray:
    """Exact readings of sensors moved rigidly (x -> rot (x - pivot) + pivot + trans) in the static
    field B(x) = b0 + grad (x - x_ref)."""
    p = (s.rmag - pivot) @ rot.T + pivot + np.asarray(trans, float)
    n = s.cosmag @ rot.T
    y = s.sum(np.sum(n * (b0 + (p - x_ref) @ grad.T), axis=1))
    return y if gains is None else y * gains


def jacobian(s: Sensors, b0: np.ndarray, grad: np.ndarray, x_ref, pivot, gains: np.ndarray | None = None) -> np.ndarray:
    """(n_channels, 6) derivative of ``readings`` at rest with respect to the rotation vector
    (3, rad, about ``pivot``) and the translation (3, m)."""
    p, n = s.rmag, s.cosmag
    field = b0 + (p - x_ref) @ grad.T
    jac = np.zeros((s.n, 6))
    for k in range(3):
        e = np.eye(3)[k]
        dn, dp = np.cross(e, n), np.cross(e, p - pivot)
        jac[:, k] = s.sum(np.sum(dn * field, axis=1) + np.sum(n * (dp @ grad.T), axis=1))
        jac[:, 3 + k] = s.sum(n @ grad[:, k])
    return jac if gains is None else jac * gains[:, None]


def projector(basis: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """I - E (E^T W E)^-1 E^T W (the projection of noisemodel.ArrayNoise for a given basis)."""
    w = np.ones(len(basis)) if weights is None else np.asarray(weights, float)
    m = basis.T @ (w[:, None] * basis)
    return np.eye(len(w)) - basis @ np.linalg.solve(m, basis.T * w[None, :])


def motion_covariance(jac: np.ndarray, rotation_rms: float, translation_rms: float) -> np.ndarray:
    """Covariance of the reading change for independent motion components with per-axis RMS
    ``rotation_rms`` [rad] and ``translation_rms`` [m] (in the analysis band)."""
    return rotation_rms**2 * jac[:, :3] @ jac[:, :3].T + translation_rms**2 * jac[:, 3:] @ jac[:, 3:].T


def mismatched_detectability(white_ref: np.ndarray, white_actual: np.ndarray) -> np.ndarray:
    """Output SNR (h^T s_k) / sqrt(h^T C_k h) of the matched filter h = C_k^+ s_0, from the whitened
    reference (s_0) and actual (s_k) topographies (both whitened with the actual noise C_k, columns
    = sources). Equals the known-topography detectability ||W s_k|| when s_0 = s_k and is never
    larger (Cauchy-Schwarz); negative when the template has the wrong polarity."""
    num = np.sum(white_ref * white_actual, axis=0)
    den = np.linalg.norm(white_ref, axis=0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)
