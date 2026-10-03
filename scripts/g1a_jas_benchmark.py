#!/usr/bin/env python3
"""G1A: Jas et al. (2026) adult analytical benchmark (REPRO).

Reproduces, with the paper's definitions (docs/literature/jas2026.md):
  * Fig. 3 (eta = 3): the validated legacy replica (legacy/replicate_figure3.py) is regenerated
    into results/g1a/fig3/;
  * Fig. 4 A-D (SNR vs depth at eta = 1.1, 2.5, 4.0, 6.0) and E (equal-SNR depth vs eta);
  * the Fig. 6 target-vs-background ("toy") depth experiment (explanatory benchmark);
and checks the analytical peak field (Eq. 1) against two independent numerical calculations:
the maximum of the full Sarvas field over the whole sensor sphere (opmsquid.sphere) and MNE's
spherical forward model with dense point magnetometers (legacy replica machinery).

Absolute noise: the paper states only eta. SNR axes use sigma_SQUID = B_max(0.8 b, h + 18 mm)
= 0.3546 pT (recovered from the published raster, R-sigma) and sigma_OPM = eta sigma_SQUID;
ratios and d_eq do not depend on it. Silent sources (the sphere centre) are excluded.

Outputs: results/g1a/{fig3/, Figure_G1A_fig4.{png,pdf}, Figure_G1A_toy_fig6.{png,pdf},
g1a_curves.csv, g1a_benchmark.json}
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "legacy"))

import numpy as np  # noqa: E402
import replicate_figure3 as R  # noqa: E402  (validated Fig. 3 replica; sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

from opmsquid import io, sphere  # noqa: E402

H, B, XI, Q = 0.095, 0.080, 0.018, 30e-9
ETA_PANELS = (1.1, 2.5, 4.0, 6.0)
PRINTED_DEQ = {2.5: 34.0, 3.0: 28.0, 4.0: 19.0}  # mm, as printed (pp. 13-15)
DRAWN_DEQ = {2.5: 35.02, 3.0: 27.53, 4.0: 19.8}  # mm, measured on the raster (jas2026.md)
OUT = ROOT / "results" / "g1a"
STATUS = "REPRO: Jas et al. 2026 Figs. 3-4 and 6, paper definitions; sigma_SQUID recovered (R-sigma)"


def curves():
    r_q = np.linspace(0.0, B, 101)[1:]  # paper grid, excludes the silent centre (R-grid)
    d = H - r_q
    b_opm, b_sq = sphere.bmax_radial(r_q, H, Q), sphere.bmax_radial(r_q, H + XI, Q)
    return d, b_opm, b_sq


def numerical_checks() -> dict:
    """Eq. 1 vs (i) Sarvas maximum over the whole sphere, (ii) MNE sphere forward (dense arc)."""
    r_q = np.array([0.008, 0.024, 0.048, 0.064, 0.072, 0.080])
    out = {}
    for label, s in (("OPM", H), ("SQUID", H + XI)):
        closed = sphere.bmax_radial(r_q, s, Q)
        sarvas = np.array([sphere.numerical_peak_radial(r, s, Q)[0] for r in r_q])
        mne_peak = R.simulate_peak_field_mne(r_q, s, dtheta_deg=0.05)["peak"]
        out[label] = dict(r_q_mm=(r_q * 1e3).tolist(), eq1_pT=(closed * 1e12).tolist(),
                          max_rel_err_sarvas_2d=float(np.max(np.abs(sarvas / closed - 1))),
                          max_rel_err_mne_sphere=float(np.max(np.abs(mne_peak / closed - 1))))
    return out


def figure4(d, b_opm, b_sq, sigma_sq) -> list[str]:
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.8))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.08, wspace=0.28, hspace=0.42)
    for ax, eta, letter in zip(axs.ravel()[:4], ETA_PANELS, "ABCD"):
        ax.plot(d * 1e3, b_opm / (eta * sigma_sq), "b-", lw=1.5, label="OPM (ξ = 0)")
        ax.plot(d * 1e3, b_sq / sigma_sq, "r--", lw=1.5, label="SQUID (ξ = 18 mm)")
        deq = sphere.equal_snr_depth(eta, H, B, 0.0, XI)
        if np.isfinite(deq):
            ax.axvline(deq * 1e3, color="k", ls="-.", lw=1.0)
            printed = PRINTED_DEQ.get(eta)
            ax.text(deq * 1e3 + 2, ax.get_ylim()[1] * 0.8 if ax.get_ylim()[1] else 1,
                    f"exact {deq * 1e3:.2f} mm\nprinted {printed:g} mm" if printed else f"exact {deq * 1e3:.2f} mm",
                    fontsize=8)
        else:
            ax.text(0.97, 0.9, "no crossing in the brain", transform=ax.transAxes, ha="right", fontsize=8)
        ax.set_xlim(0, 95)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{letter}   η = {eta:g}", loc="left")
        ax.set_xlabel("Dipole depth (d) [mm]")
        ax.set_ylabel("SNR")
        if letter == "A":
            ax.legend(frameon=False, fontsize=8)
    ax = axs[1, 1]
    etas = np.linspace(1.0, 6.0, 501)
    exact = np.array([sphere.equal_snr_depth(e, H, B, 0.0, XI) for e in etas]) * 1e3
    # paper's grid rule: deepest grid depth with SNR_OPM > SNR_SQUID (floor at the deepest grid point)
    grid = np.array([d[b_opm / b_sq > e].max() * 1e3 if np.any(b_opm / b_sq > e) else np.nan for e in etas])
    ax.plot(etas, grid, color="0.6", lw=2.5, label="paper grid rule (b/100 steps)")
    ax.plot(etas, exact, "k-", lw=1.2, label="exact root of Eq. 3")
    eta0, eta1 = sphere.centre_limit_ratio(H, 0.0, XI), sphere.signal_ratio(H - B)[0]
    # OPM ahead above the curve (shallower than d_eq); everywhere for eta < eta0; nowhere for eta > eta1
    shade_to = np.where(etas < eta0, H * 1e3, np.where(np.isfinite(exact), exact, (H - B) * 1e3))
    ax.fill_between(etas, (H - B) * 1e3, shade_to, color="tab:blue", alpha=0.12, lw=0, label="OPM SNR > SQUID SNR")
    for e, lab in ((eta0, "η0"), (eta1, "η1")):
        ax.axvline(e, color="k", ls=":", lw=0.8)
        top = lab == "η0"  # eta0 label in the empty band above the brain surface, eta1 beyond the curve
        ax.text(e + 0.06, 8 if top else 32, f"{lab} = {e:.4f} (printed 1.7)" if top else f"{lab} = {e:.4f}\n(printed 5.3)",
                fontsize=7, va="center")
    ax.axhline((H - B) * 1e3, color="k", ls="-.", lw=0.8)
    ax.set_ylim(95, 0)
    ax.set_xlim(0.8, 6.3)
    ax.set_xlabel("Relative noise level (η)")
    ax.set_ylabel("d_eq [mm]")
    ax.set_title("E   Equal-SNR depth", loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax = axs[1, 2]
    ax.axis("off")
    rows = [(e, sphere.equal_snr_depth(e, H, B, 0.0, XI) * 1e3) for e in (1.7, 2.0, 2.5, 3.0, 4.0, 5.0)]
    txt = "Exact equal-SNR depths (Eq. 3) [mm]\n(printed / drawn in the paper)\n\n" + "\n".join(
        f"η = {e:3.1f}: {v:6.2f}" + (f"  ({PRINTED_DEQ[e]:g} / {DRAWN_DEQ[e]:g})" if e in PRINTED_DEQ else "")
        for e, v in rows)
    txt += (f"\n\nCrossing exists for\n{eta0:.4f} < η < {eta1:.4f}\n\nσ_SQUID = {sigma_sq * 1e12:.4f} pT"
            "\n(recovered), σ_OPM = η σ_SQUID")
    ax.text(0.0, 1.0, txt, va="top", family="monospace", fontsize=8)
    fig.suptitle("G1A — Jas et al. 2026 Fig. 4 (reproduction, adult sphere h = 95 mm, b = 80 mm, Q = 30 nAm)",
                 fontsize=10)
    files = []
    for ext in ("png", "pdf"):
        p = OUT / f"Figure_G1A_fig4.{ext}"
        fig.savefig(p, dpi=200)
        files.append(str(p.relative_to(ROOT)))
    plt.close(fig)
    return files


def toy_experiment() -> tuple[list[str], dict]:
    """Fig. 6: target B2 at r_Q = 0.6 b (d = 47 mm); noise dipoles B1 at 0.4 b (d = 63 mm, deep)
    and B3 at 0.8 b (d = 31 mm, superficial); all 30 nAm, +y; SNR = ratio of peak fields on the
    sensor sphere s = h + xi. The text prints 32/48/64 mm, which are the r_Q values (U-J6)."""
    xi = np.linspace(0.0, 0.060, 241)
    r_q = {"B1": 0.4 * B, "B2": 0.6 * B, "B3": 0.8 * B}
    field = {k: np.array([sphere.bmax_radial(r, H + x, Q)[0] for x in xi]) for k, r in r_q.items()}
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.85, bottom=0.15, wspace=0.25)
    styles = {"B3": ("r--", 31), "B2": ("k:", 47), "B1": ("b-", 63)}
    for k, (st, dd) in styles.items():
        axs[0].plot(xi * 1e3, field[k] * 1e12, st, lw=1.5, label=f"{k} (d = {dd} mm)")
    axs[1].plot(xi * 1e3, field["B2"] / field["B1"], "b-", lw=1.5, label="B2 / B1 (deep noise)")
    axs[1].plot(xi * 1e3, field["B2"] / field["B3"], "r--", lw=1.5, label="B2 / B3 (superficial noise)")
    axs[1].axhline(1.0, color="k", ls=":", lw=0.8)
    for ax, title in zip(axs, ("A   Signal [pT]", "B   SNR (ratio of peak fields)")):
        ax.axvline(18, color="k", ls=":", lw=0.8)
        ax.set_xlim(60, 0)
        ax.set_xlabel("ξ [mm]")
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("G1A — Jas et al. 2026 Fig. 6 toy experiment (reproduction; explanatory benchmark only)", fontsize=10)
    files = []
    for ext in ("png", "pdf"):
        p = OUT / f"Figure_G1A_toy_fig6.{ext}"
        fig.savefig(p, dpi=200)
        files.append(str(p.relative_to(ROOT)))
    plt.close(fig)
    at = {f"xi_{x:g}mm": {k: float(sphere.bmax_radial(r, H + x * 1e-3, Q)[0] * 1e12) for k, r in r_q.items()}
          for x in (0, 18, 50, 60)}
    for v in at.values():
        v["B2/B1"], v["B2/B3"] = v["B2"] / v["B1"], v["B2"] / v["B3"]
    return files, at


def _deq_with_sigma(sigma_sq: float, eta: float = 3.0) -> float:
    """Equal-SNR depth [mm] found directly from the SNR curves with an arbitrary absolute
    sigma_SQUID (bisection on SNR_OPM - SNR_SQUID); equals the eta-only Eq. 3 root."""
    def diff(d):
        r = np.array([H - d])
        return sphere.bmax_radial(r, H, Q)[0] / (eta * sigma_sq) - sphere.bmax_radial(r, H + XI, Q)[0] / sigma_sq
    lo, hi = H - B, H - 1e-4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if diff(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi) * 1e3


REPLICA_STATUS = ("REPRO: Jas et al. 2026 Fig. 3, written by the preserved legacy replica (legacy/replicate_figure3.py, "
                  "7/7 checks against the published raster) and stamped by scripts/g1a_jas_benchmark.py")


def stamp_replica(out: Path) -> None:
    """The legacy replica writes no status or commit of its own: add both to its summary JSON and a
    status line to its data CSV (the replica's code and numbers are unchanged)."""
    summary = json.loads((out / "figure3_summary.json").read_text())
    io.write_json(dict(summary, status=REPLICA_STATUS), out / "figure3_summary.json")
    data = out / "figure3_data.csv"
    text = data.read_text()
    with open(data, "w", newline="") as fh:
        io.csv_status(fh, REPLICA_STATUS)
        fh.write(text)


def main() -> dict:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    R.main(["--outdir", str(OUT / "fig3")])  # validated replica (7/7 checks vs the published raster)
    stamp_replica(OUT / "fig3")
    d, b_opm, b_sq = curves()
    sigma_sq = float(sphere.bmax_radial(0.8 * B, H + XI, Q)[0])
    files = figure4(d, b_opm, b_sq, sigma_sq)
    toy_files, toy = toy_experiment()
    checks = numerical_checks()
    with open(OUT / "g1a_curves.csv", "w", newline="") as fh:
        io.csv_status(fh, STATUS)
        wr = csv.writer(fh)
        wr.writerow(["depth_mm", "B_OPM_pT", "B_SQUID_pT", "SNR_SQUID"] + [f"SNR_OPM_eta{e:g}" for e in (1.1, 2.5, 3, 4, 6)])
        for i in np.argsort(d):
            wr.writerow([f"{d[i] * 1e3:.3f}", f"{b_opm[i] * 1e12:.6g}", f"{b_sq[i] * 1e12:.6g}",
                         f"{b_sq[i] / sigma_sq:.6g}"] + [f"{b_opm[i] / (e * sigma_sq):.6g}" for e in (1.1, 2.5, 3, 4, 6)])
    summary = dict(
        status=STATUS,
        parameters=dict(h_mm=95, b_mm=80, Q_nAm=30, xi_opm_mm=0, xi_squid_mm=18, sigma_squid_pT=sigma_sq * 1e12,
                        grid="r_Q = linspace(0, b, 101)[1:] (d = 15 ... 94.2 mm)"),
        excluded="r_Q = 0 (sphere centre, silent): B_r = 0 at every sensor; excluded from ratios and crossings",
        eta_range_with_crossing=[sphere.centre_limit_ratio(H, 0.0, XI), float(sphere.signal_ratio(H - B)[0])],
        d_eq_mm={f"{e:g}": sphere.equal_snr_depth(e, H, B, 0.0, XI) * 1e3 for e in (1.1, 1.7, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)},
        printed_vs_exact=[
            dict(item="Fig. 3 d_eq (eta 3)", printed_mm=28, drawn_mm=27.53, exact_mm=sphere.equal_snr_depth(3.0) * 1e3),
            dict(item="Fig. 4B d_eq (eta 2.5)", printed_mm=34, drawn_mm=35.02, exact_mm=sphere.equal_snr_depth(2.5) * 1e3),
            dict(item="Fig. 4C d_eq (eta 4)", printed_mm=19, drawn_mm=19.8, exact_mm=sphere.equal_snr_depth(4.0) * 1e3),
            dict(item="eta0 / eta1", printed=[1.7, 5.3], exact=[sphere.centre_limit_ratio(), float(sphere.signal_ratio(H - B)[0])]),
            dict(item="Fig. 6 depths", printed_text_mm=[32, 48, 64], printed_caption_mm=[63, 47, 31],
                 note="text values are r_Q (0.4 b, 0.6 b, 0.8 b); caption values reproduce the curves")],
        markers=dict(
            fig3=("dotted line at 27.53 mm: the deepest depth with SNR_OPM > SNR_SQUID on d = linspace(h - b, h, 250); the grid "
                  "size was fitted to the published raster (many sizes give 27.50-27.56 mm) and the authors' rule is unknown; "
                  "the exact Eq. 3 root is 27.665 mm (U-J1)"),
            fig4=("exact Eq. 3 roots. The paper's drawn markers (35.0 and 19.8 mm) are depths on its r_Q = linspace(0, b, 101) grid; "
                  "its printed 34 and 19 mm equal those grid depths truncated in SI floating point (0.034999... m, 0.019799... m) (U-J2)")),
        numerical_checks=checks,
        toy_experiment=dict(toy, xi_range_mm=[0, 60], xi_range_note="methods text; the paper's panels show about 0-50 mm",
                            text_as_depth_reading_pT_at_xi0={f"{d:g}mm": float(sphere.bmax_radial(H - d * 1e-3, H, Q)[0] * 1e12)
                                                             for d in (32, 48, 64)},
                            text_as_depth_note="the text's 32/48/64 mm read as depths below the scalp does not reproduce the "
                                               "figure; read as r_Q (0.4 b, 0.6 b, 0.8 b; caption depths 63/47/31 mm) it does (U-J6)"),
        absolute_noise=("only eta is stated; sigma_SQUID = B_max(0.8 b, h + 18 mm) recovered from the Fig. 3 SNR axis; "
                        "d_eq and all ratios are independent of it"),
        d_eq_independent_of_sigma={f"sigma_x{k:g}": _deq_with_sigma(k * sigma_sq) for k in (0.5, 1.0, 2.0)},
        outputs=["results/g1a/fig3/"] + files + toy_files + ["results/g1a/g1a_curves.csv"],
        runtime_s=time.time() - t0)
    io.write_json(summary, OUT / "g1a_benchmark.json")
    print(json.dumps(io.json_safe({k: summary[k] for k in ("d_eq_mm", "eta_range_with_crossing")}), indent=1))
    for k, v in checks.items():
        print(f"{k}: Eq.1 vs Sarvas 2-D max {v['max_rel_err_sarvas_2d']:.1e}; vs MNE sphere {v['max_rel_err_mne_sphere']:.1e}")
    return summary


if __name__ == "__main__":
    main()
