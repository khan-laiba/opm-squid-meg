#!/usr/bin/env python3
"""G3A: Jas et al. (2026) head-size benchmark, Table 1 and Fig. 5 (REPRO), with an idealised
fixed-shell contrast (NEW).

REPRO, with the paper's definitions (docs/literature/jas2026.md section 3): four concentric-sphere
heads (Table 1, (h, b) = newborn (55, 48), 1 year (70, 62), 8 years (85, 73), adult (95, 80) mm);
OPM on the scalp (s = h) and SQUID on a size-following shell s = h + 18 mm for every head (U-J5).
Equal-SNR depth d_eq(eta) from Eq. 3 (exact root) and from the paper's grid rule (deepest depth on
r_Q = linspace(0, b, 101)[1:] with SNR_OPM > SNR_SQUID); Fig. 5A (d_eq in mm below the scalp) and
Fig. 5B ("normalized d_eq" = (d_eq - (h - b)) / b, the depth below the brain surface as a fraction
of b, U-J4). Printed checks: 50 % (newborn) and 15 % (adult) at eta = 3. The brain-volume fraction
with SNR_OPM > SNR_SQUID is reported separately (the text speaks of volume, the figure plots a
radial fraction).

NEW, idealised: the same heads inside one fixed adult shell (s = 113 mm, concentric), i.e. the
gap a child has in an adult-size helmet (18, 28, 43, 58 mm). Real heads do not sit concentrically
in a helmet; the realistic fixed-helmet experiment is G3B.

Outputs: results/g3a/ (Figure_G3A_fig5.png/.pdf, g3a_size_benchmark.json, g3a_deq.csv)
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from opmsquid import io, sphere  # noqa: E402

HEADS = {"Infant (newborn)": (0.055, 0.048), "Infant (1-yr)": (0.070, 0.062), "Child (8-yr)": (0.085, 0.073),
         "Adult": (0.095, 0.080)}  # Table 1 (h, b) [m]
COLORS = {"Infant (newborn)": "tab:blue", "Infant (1-yr)": "tab:orange", "Child (8-yr)": "tab:green", "Adult": "tab:red"}
XI_SQUID = 0.018  # size-following stand-off (Fig. 5 caption)
ADULT_SHELL = 0.095 + XI_SQUID  # fixed adult shell radius for the NEW contrast
PRINTED_NORM_ETA3 = {"Infant (newborn)": 50.0, "Adult": 15.0}  # % (p. 16)
ETAS = np.round(np.arange(1.0, 6.0001, 0.01), 2)
OUT = ROOT / "results" / "g3a"
STATUS = ("REPRO: Jas et al. 2026 Table 1 and Fig. 5 (size-following SQUID shell, paper definitions); "
          "NEW: idealised concentric fixed adult shell")


def deq_exact(eta, h, b, xi_squid):
    return sphere.equal_snr_depth(eta, h, b, 0.0, xi_squid)


def deq_grid(eta, h, b, xi_squid):
    """Paper's grid rule: deepest grid depth with SNR_OPM > SNR_SQUID (NaN if none)."""
    r_q = np.linspace(0.0, b, 101)[1:]
    ratio = sphere.bmax_radial(r_q, h) / sphere.bmax_radial(r_q, h + xi_squid)
    ok = ratio > eta
    return float((h - r_q[ok]).max()) if ok.any() else np.nan


def volume_fraction(d_eq, h, b, eta, eta0):
    """Fraction of the brain-sphere volume with SNR_OPM > SNR_SQUID (all of it below eta0)."""
    if eta < eta0:
        return 1.0
    if not np.isfinite(d_eq):
        return 0.0
    r_eq = h - d_eq
    return float((b**3 - r_eq**3) / b**3)


def head_table(xi_of):
    out = {}
    for name, (h, b) in HEADS.items():
        xi = xi_of(h)
        eta0, eta1 = sphere.centre_limit_ratio(h, 0.0, xi), float(sphere.signal_ratio(h - b, h, 0.0, xi)[0])
        row = dict(h_mm=h * 1e3, b_mm=b * 1e3, layer_mm=(h - b) * 1e3, xi_squid_mm=xi * 1e3, eta0=eta0, eta1=eta1)
        for eta in (2.5, 3.0, 4.0, 6.0):
            de, dg = deq_exact(eta, h, b, xi), deq_grid(eta, h, b, xi)
            row[f"eta{eta:g}"] = dict(d_eq_exact_mm=de * 1e3, d_eq_grid_mm=dg * 1e3,
                                      normalized_exact_pct=(de - (h - b)) / b * 100, normalized_grid_pct=(dg - (h - b)) / b * 100,
                                      d_eq_over_b_pct=de / b * 100, volume_fraction_pct=volume_fraction(de, h, b, eta, eta0) * 100)
        out[name] = row
    return out


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    follow = head_table(lambda h: XI_SQUID)
    fixed = head_table(lambda h: ADULT_SHELL - h)
    checks = {name: dict(printed_pct=v, exact_pct=follow[name]["eta3"]["normalized_exact_pct"],
                         grid_pct=follow[name]["eta3"]["normalized_grid_pct"]) for name, v in PRINTED_NORM_ETA3.items()}

    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    for name, (h, b) in HEADS.items():
        c = COLORS[name]
        for (label, xi_of), row in ((("size-following", lambda hh: XI_SQUID), 0), (("fixed adult shell", lambda hh: ADULT_SHELL - hh), 1)):
            xi = xi_of(h)
            ex = np.array([deq_exact(e, h, b, xi) for e in ETAS])
            gr = np.array([deq_grid(e, h, b, xi) for e in ETAS])
            eta0 = sphere.centre_limit_ratio(h, 0.0, xi)
            ex_f = np.where(ETAS < eta0, h - b / 100, ex)  # flat floor below eta0, as drawn in the paper
            if row == 0:
                axs[0, 0].plot(ETAS, ex_f * 1e3, color=c, label=name)
                axs[0, 0].plot(ETAS, gr * 1e3, color=c, lw=0.8, ls=":")
                axs[0, 1].plot(ETAS, (ex_f - (h - b)) / b * 100, color=c, label=name)
                axs[0, 1].plot(ETAS, (gr - (h - b)) / b * 100, color=c, lw=0.8, ls=":")
                vf = [volume_fraction(d, h, b, e, eta0) * 100 for d, e in zip(ex, ETAS)]
                axs[1, 0].plot(ETAS, vf, color=c, label=f"{name}, size-following")
            else:
                ahead = f", OPM ahead at every depth up to eta = {ETAS[-1]:g}" if eta0 > ETAS[-1] else ""
                axs[1, 1].plot(ETAS, (ex_f - (h - b)) / b * 100, color=c, label=f"{name} (xi = {xi * 1e3:.0f} mm){ahead}")
                vf = [volume_fraction(d, h, b, e, eta0) * 100 for d, e in zip(ex, ETAS)]
                axs[1, 0].plot(ETAS, vf, color=c, ls="--", label=f"{name}, fixed shell")
    for ax, ylab, top in ((axs[0, 0], "Equal-SNR source depth d_eq [mm]", 0), (axs[0, 1], "Normalized d_eq [%]", 0),
                          (axs[1, 1], "Normalized d_eq [%]", 0)):
        ax.axhline(top, color="k", ls="-.", lw=0.8)
        ax.axvline(3, color="k", ls=":", lw=0.8)
        ax.invert_yaxis()
        ax.set_xlim(0.8, 6.3)
        ax.set_xlabel("Relative noise level (eta)")
        ax.set_ylabel(ylab)
    axs[0, 0].set_ylim(95, -5)
    axs[0, 1].set_ylim(110, -5)
    axs[1, 1].set_ylim(135, -5)  # room below the 100 % floor for the legend
    axs[1, 1].set_yticks(range(0, 101, 20))
    for name in PRINTED_NORM_ETA3:
        v = follow[name]["eta3"]["normalized_exact_pct"]
        axs[0, 1].annotate(f"{v:.1f} % (printed {PRINTED_NORM_ETA3[name]:g} %)", (3, v), (3.4, v - 8), fontsize=8,
                           arrowprops=dict(arrowstyle="-", lw=0.6))
    axs[0, 0].set_title("A  REPRO: d_eq, SQUID at h + 18 mm (solid exact, dotted grid rule)", loc="left", fontsize=9)
    axs[0, 1].set_title("B  REPRO: depth below the brain surface / b", loc="left", fontsize=9)
    axs[1, 0].set_title("C  brain volume with SNR_OPM > SNR_SQUID", loc="left", fontsize=9)
    axs[1, 1].set_title("D  NEW (idealised): all heads concentric in the adult shell (s = 113 mm)", loc="left", fontsize=9)
    axs[1, 0].set_xlim(0.8, 6.3)
    axs[1, 0].set_xlabel("Relative noise level (eta)")
    axs[1, 0].set_ylabel("volume fraction [%]")
    axs[1, 0].axvline(3, color="k", ls=":", lw=0.8)
    for ax in axs.ravel():
        ax.legend(fontsize=7, loc="lower right" if ax is not axs[1, 0] else "lower left")
    fig.suptitle("G3A - Jas et al. 2026 Fig. 5 (size-following benchmark) and a fixed-shell contrast", fontsize=10)
    fig.tight_layout()
    files = []
    for ext in ("png", "pdf"):
        p = OUT / f"Figure_G3A_fig5.{ext}"
        fig.savefig(p, dpi=200)
        files.append(str(p.relative_to(ROOT)))
    plt.close(fig)

    with open(OUT / "g3a_deq.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["head", "shell", "eta", "d_eq_exact_mm", "d_eq_grid_mm", "normalized_exact_pct", "volume_fraction_pct"])
        for name, (h, b) in HEADS.items():
            for shell, xi in (("size-following", XI_SQUID), ("fixed-adult", ADULT_SHELL - h)):
                eta0 = sphere.centre_limit_ratio(h, 0.0, xi)
                for e in ETAS[::10]:
                    de, dg = deq_exact(e, h, b, xi), deq_grid(e, h, b, xi)
                    wr.writerow([name, shell, f"{e:.2f}", f"{de * 1e3:.3f}", f"{dg * 1e3:.3f}", f"{(de - (h - b)) / b * 100:.2f}",
                                 f"{volume_fraction(de, h, b, e, eta0) * 100:.2f}"])
    summary = dict(status=STATUS, size_following=follow, fixed_adult_shell=fixed, printed_checks=checks,
                   notes=["d_eq is independent of the dipole moment and of the absolute noise; only eta matters.",
                          "Normalized d_eq = (d_eq - (h - b)) / b as plotted in Fig. 5B (U-J4); d_eq / b is listed for the "
                          "printed wording.", "The fixed-shell panel is an idealised NEW contrast, not a result of the paper."],
                   outputs=files + ["results/g3a/g3a_deq.csv"], runtime_s=time.time() - t0)
    io.write_json(summary, OUT / "g3a_size_benchmark.json")
    for name in HEADS:
        f3, x3 = follow[name]["eta3"], fixed[name]["eta3"]
        print(f"{name:17s} eta0 {follow[name]['eta0']:.4f} eta1 {follow[name]['eta1']:.4f} | eta 3: d_eq {f3['d_eq_exact_mm']:.2f} mm "
              f"(grid {f3['d_eq_grid_mm']:.2f}), norm {f3['normalized_exact_pct']:.1f} % (grid {f3['normalized_grid_pct']:.1f}), "
              f"volume {f3['volume_fraction_pct']:.1f} % | fixed shell: norm {x3['normalized_exact_pct']:.1f} %, volume "
              f"{x3['volume_fraction_pct']:.1f} %")


if __name__ == "__main__":
    main()
