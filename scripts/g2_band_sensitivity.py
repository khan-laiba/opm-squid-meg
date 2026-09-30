#!/usr/bin/env python3
"""G2 supplement: frequency-band and OPM sensor-response sensitivity (NEW).

The primary G2 comparison uses one band (1-40 Hz). Here the whole noise model is rebuilt in other
bands: intrinsic variance = ASD^2 x ENBW of the band's composite filter, brain-noise scale
recalibrated on the good gradiometers' (task baseline - empty room) variance in that band, room
field refitted to the empty-room recording in that band. Same arrays, targets and gains as the
primary run (scripts/g2_adult_comparison.py).

OPM sensor response (A-OPM-BW, declared): a first-order low-pass with a 100-Hz corner, as for
SERF magnetometers of this class. Brain signals and brain/room noise pass through it, the
intrinsic noise is specified at the output, so the OPM target and brain/room terms are scaled by
the band-limited amplitude gain for a 1/f spectrum.

Also verifies the PSD-to-variance rule by simulation in every band.

Outputs: results/g2/g2_band_sensitivity.json, results/g2/Figure_G2_bands.png
"""
from __future__ import annotations

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

import g2_adult_comparison as G  # noqa: E402
from opmsquid import background, g2, io, noise, noisemodel  # noqa: E402

BANDS = ((1.0, 40.0), (1.0, 10.0), (8.0, 30.0), (30.0, 80.0))  # first = primary band
OPM_CORNER_HZ = 100.0
OUT = ROOT / "results" / "g2"


def opm_response_gain(filt: noise.AnalysisFilter, corner: float = OPM_CORNER_HZ) -> float:
    """RMS amplitude gain of a first-order low-pass (corner ``corner``) for a 1/f spectrum seen
    through the band's composite filter."""
    f = np.linspace(0.05, filt.fs / 2, 200_001)
    band = filt.power_response(f) / f
    lp = 1.0 / (1.0 + (f / corner) ** 2)
    return float(np.sqrt(np.trapezoid(band * lp, f) / np.trapezoid(band, f)))


def simulated_variance_check(filt: noise.AnalysisFilter, asd: float = 3.5e-15, seconds: float = 300.0) -> float:
    x = noise.white_noise(asd, filt.fs, (4, int(seconds * filt.fs)), np.random.default_rng(11))
    y = filt.apply(x)[:, int(5 * filt.fs):-int(5 * filt.fs)]
    return float(np.mean(y**2) / (asd**2 * filt.enbw()))


def main():
    t0 = time.time()
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "g2_adult.toml").read_text())
    st = G.Study(cfg)
    headline = cfg["conditions"]["headline"]
    bads = cfg["sensors"]["bads"]
    arrays = g2.build_arrays(st.subject, st.dig)
    squid = arrays["squid"]
    good = ~np.isin(squid.info.ch_names, bads)
    grads = good & (squid.kinds == "grad")
    G_t, G_g = {}, {}
    for name, a in arrays.items():
        G_t[name], G_g[name] = st.gains(a)
    unit = {name: st.unit_brain(G_g[name]) for name in arrays}
    out = dict(status="NEW (G2 supplement: frequency band and OPM sensor response)", bands={},
               opm_response=dict(model="first-order low-pass", corner_hz=OPM_CORNER_HZ, spectrum="1/f"))
    for lo, hi in BANDS:
        key = f"{lo:g}-{hi:g}Hz"
        filt = noise.AnalysisFilter(fs=st.dig["sfreq"], l_freq=lo, h_freq=hi, order=cfg["band"]["order"])
        meas = g2.measured_noise(squid.info, filt, bads)
        env = meas["environment"]
        target = float(np.nanmedian(meas["brain"][grads]))
        bs = background.calibrate(unit["squid"], grads, target)
        enbw = filt.enbw()
        g_opm = opm_response_gain(filt)
        row = dict(enbw_hz=enbw, simulated_over_predicted_variance=simulated_variance_check(filt),
                   brain_rms_grad_fT_cm=float(np.sqrt(target) * 1e13),
                   brain_rms_mag_fT_measured=float(np.sqrt(np.nanmedian(meas["brain"][good & (squid.kinds == "mag")])) * 1e15),
                   env_explained=env.explained_fraction, opm_response_gain=g_opm)
        for resp in ("ideal", "opm_response"):
            res = {}
            for name, a in arrays.items():
                nz = g2.array_noise(a, G_g[name], st.src.grid_area, bs, env, enbw, st.opm_asd)
                gt = G_t[name]
                if resp == "opm_response" and name != "squid":
                    s2 = g_opm**2
                    nz = noisemodel.ArrayNoise(nz.intrinsic_var, s2 * nz.brain_cov, s2 * nz.env_cov, nz.ext_basis)
                    gt = gt * g_opm
                res[name] = G.evaluate_all(gt * st.q, a, nz, headline)
                if resp == "ideal":
                    row.setdefault("intrinsic_share_of_variance", {})[name] = float(
                        np.median(nz.intrinsic_var / np.diag(nz.covariance("intrinsic+brain+env"))))
            row[resp] = {k: dict(ratio=float(2 ** v["median_log2"]), ci95=[float(2 ** x) for x in v["ci95"]],
                                 share_opm_better=v["share_opm_better"])
                         for k, v in G.comparisons(st, res, headline, n_boot=200).items()}
        out["bands"][key] = row
        print(f"{key}: ENBW {enbw:.1f} Hz, brain grad {row['brain_rms_grad_fT_cm']:.1f} fT/cm, OPM gain {g_opm:.3f}, var check "
              f"{row['simulated_over_predicted_variance']:.3f}; dense/combined {row['ideal']['opm_dense/combined/intrinsic+brain']['ratio']:.3f}"
              f" (with response {row['opm_response']['opm_dense/combined/intrinsic+brain']['ratio']:.3f}), matched "
              f"{row['ideal']['opm99/combined/intrinsic+brain']['ratio']:.3f}", flush=True)
    out["runtime_s"] = time.time() - t0
    io.write_json(out, OUT / "g2_band_sensitivity.json")

    fig, axs = plt.subplots(1, len(headline), figsize=(12, 4.2), sharey=True)
    keys = list(out["bands"])
    x = np.arange(len(keys))
    for ax, cond in zip(axs, headline):
        for j, (a, ref, mk) in enumerate((("opm99", "combined", "o"), ("opm_dense", "combined", "s"), ("opm_dense", "grad", "^"),
                                           ("opm_dense", "mag", "v"))):
            for resp, alpha in (("ideal", 1.0), ("opm_response", 0.45)):
                y = [out["bands"][k][resp][f"{a}/{ref}/{cond}"]["ratio"] for k in keys]
                ax.plot(x + 0.06 * (j - 1.5), y, mk, color={"opm99": "tab:green", "opm_dense": "tab:purple"}[a],
                        alpha=alpha, mfc="none" if resp == "opm_response" else None,
                        label=f"{G.LABEL[a]} / {G.LABEL[ref]}" + (" (100-Hz OPM response)" if resp == "opm_response" else ""))
        ax.axhline(1, color="0.5", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(keys)
        ax.set_title(f"median detectability ratio, {cond}", fontsize=9)
        ax.set_xlabel("analysis band (noise model recalibrated per band)")
    axs[0].set_ylabel("OPM / Neuromag")
    axs[0].legend(fontsize=6.5, ncol=2)
    fig.suptitle("G2 supplement: frequency-band and OPM sensor-response sensitivity (OPM 15 fT/sqrt(Hz))", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_G2_bands.png", dpi=150)
    plt.close(fig)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
