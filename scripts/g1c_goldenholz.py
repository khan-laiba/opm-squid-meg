#!/usr/bin/env python3
"""G1C: Goldenholz et al. (2009) cortical SNR maps, MEG part.

Status: ADAPT (MNE sample subject and a task-baseline recording instead of the paper's four
subjects and spontaneous recordings; unstated details chosen in configs/goldenholz_reference.toml).
The OPM maps are a NEW extension (the paper's noise model has no sensor noise, so an intrinsic
OPM/SQUID noise term is added there and labelled).

Inputs: cache/fullres/{neuromag_bem006,neuromag_bem06,opm_bem006}.npy
Outputs: results/g1c/ (maps, g1c_summary.json, g1c_oct6_values.csv)
"""
from __future__ import annotations

import csv
import sys
import time
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mne  # noqa: E402

from opmsquid import anatomy, goldenholz, io, neuromag, noise, paths, plotting  # noqa: E402

OUT = ROOT / "results" / "g1c"
CH_SETS = ("mag", "grad", "pooled")
OPM_ASD, SQUID_ASD = 15e-15, {"mag": 3.5e-15, "grad": 3.6e-13}


def recorded_noise_variance(ch_names, projector, filt) -> tuple[np.ndarray, int]:
    """Per-channel variance of pre-stimulus baselines (-200 to 0 ms), SSP applied, 0.5-100 Hz."""
    raw = mne.io.read_raw_fif(paths.SAMPLE_MEG / neuromag.RAW_FILE, preload=True, verbose=False)
    events = mne.find_events(raw, stim_channel="STI 014", verbose=False)
    picks = [raw.ch_names.index(n) for n in ch_names]
    data = projector @ raw.get_data(picks=picks)
    data = filt.apply(data - data.mean(axis=1, keepdims=True))
    sf = raw.info["sfreq"]
    n_edge = int(5 * sf)
    wins = [(int(e[0] - raw.first_samp - 0.2 * sf), int(e[0] - raw.first_samp)) for e in events]
    seg = np.concatenate([data[:, a:b] for a, b in wins if a > n_edge and b < data.shape[1] - n_edge], axis=1)
    seg = seg - seg.mean(axis=1, keepdims=True)
    return np.mean(seg**2, axis=1), seg.shape[1]


def channel_sets(kinds):
    return {"mag": kinds == "mag", "grad": kinds == "grad", "pooled": np.ones(len(kinds), bool)}


def main():
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "goldenholz_reference.toml").read_text())
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    valid_idx = np.load(paths.CACHE / "fullres" / "valid_index.npy")
    col_of = np.full(cortex.n, -1)
    col_of[valid_idx] = np.arange(len(valid_idx))
    info = neuromag.load_info("T3")
    kinds = neuromag.channel_kinds(info)
    sets = channel_sets(kinds)
    raw_info = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    proj, n_proj, _ = mne.proj.make_projector(raw_info["projs"], info.ch_names)
    filt = noise.AnalysisFilter(fs=raw_info["sfreq"], l_freq=0.5, h_freq=100.0, order=4)
    rec_var, n_rec = recorded_noise_variance(info.ch_names, proj, filt)
    print(f"recorded noise: {n_rec} baseline samples, SSP rank {n_proj}; median RMS mag "
          f"{np.median(np.sqrt(rec_var[kinds == 'mag'])) * 1e15:.1f} fT, grad {np.median(np.sqrt(rec_var[kinds == 'grad'])) * 1e13:.1f} fT/cm")

    rng = np.random.default_rng(2009)
    noise_cols = goldenholz.poisson_disk(cortex.rr[valid_idx], 0.007, rng)  # positions in valid_idx
    print(f"noise-model sources: {len(noise_cols)} (7 mm Poisson-disk grid)")
    # patch centroids: valid oct-6 vertices
    n_lh_full = int(np.sum(cortex.hemi == 0))
    oct_global = np.concatenate([subject.src[0]["vertno"], subject.src[1]["vertno"] + n_lh_full])
    oct_ok = cortex.valid[oct_global]
    centroids = oct_global[oct_ok]
    density = cfg["sources"]["patch_density_pAm_per_mm2"] * 1e-12 / 1e-6  # A m per m^2
    t0 = time.time()
    members = {r: goldenholz.geodesic_patches(cortex.adjacency, centroids, r * 1e-3, cortex.valid)
               for r in cfg["sources"]["patch_radii_mm"]}
    print(f"patches built in {time.time() - t0:.0f} s; median members " +
          ", ".join(f"{r:g} mm: {np.median([len(m) for m in members[r]]):.0f}" for r in members))

    results, summary = {}, dict(status="ADAPT (Goldenholz et al. 2009, MEG part, MNE sample subject); OPM rows NEW",
                                config=cfg, n_valid_vertices=int(len(valid_idx)), n_noise_sources=int(len(noise_cols)),
                                n_centroids=int(len(centroids)), recorded_noise=dict(n_samples=int(n_rec), ssp_rank=int(n_proj)))
    for bem_label, fname in (("probable_default", "neuromag_bem006"), ("as_printed", "neuromag_bem06")):
        g = proj.astype(np.float32) @ np.load(paths.CACHE / "fullres" / f"{fname}.npy")
        aat = np.sum(g[:, noise_cols].astype(np.float64) ** 2, axis=1)
        s_s2, per_type = goldenholz.calibrate_source_variance(rec_var, aat, kinds)
        model_var = s_s2 * aat
        res = dict(source_sd_nAm=float(np.sqrt(s_s2) * 1e9), per_type_median_nAm2=per_type)
        for noise_name, var in (("model", model_var), ("recorded", rec_var)):
            for cs in CH_SETS:
                m = sets[cs]
                focal = goldenholz.eq1_snr_db(g[m] * cfg["sources"]["focal_nAm"] * 1e-9, var[m])
                res[f"focal/{noise_name}/{cs}"] = focal
                for r, mem in members.items():
                    topo = np.stack([g[m][:, col_of[mm]].astype(np.float64) @ (density * cortex.area[mm]) for mm in mem], axis=1)
                    res[f"patch{r:g}/{noise_name}/{cs}"] = goldenholz.eq1_snr_db(topo, var[m])
        results[bem_label] = res
        print(f"{bem_label}: noise-source SD {res['source_sd_nAm']:.2f} nAm (paper: 1.6-1.9 nAm per source)")

    # NEW extension: OPM with the same calibrated brain-noise model plus intrinsic noise (and SQUIDs with brochure noise)
    g_sq = proj.astype(np.float32) @ np.load(paths.CACHE / "fullres" / "neuromag_bem006.npy")
    g_opm = np.load(paths.CACHE / "fullres" / "opm_bem006.npy")
    s_s2 = results["probable_default"]["source_sd_nAm"] ** 2 * 1e-18
    enbw = filt.enbw()
    ext = {}
    for label, gg, intrinsic in (("opm", g_opm, np.full(g_opm.shape[0], OPM_ASD**2 * enbw)),
                                 ("mag", g_sq[sets["mag"]], np.full(sets["mag"].sum(), SQUID_ASD["mag"] ** 2 * enbw)),
                                 ("grad", g_sq[sets["grad"]], np.full(sets["grad"].sum(), SQUID_ASD["grad"] ** 2 * enbw))):
        var = s_s2 * np.sum(gg[:, noise_cols].astype(np.float64) ** 2, axis=1) + intrinsic
        ext[f"focal/{label}"] = goldenholz.eq1_snr_db(gg * 10e-9, var)
        for r, mem in members.items():
            topo = np.stack([gg[:, col_of[mm]].astype(np.float64) @ (density * cortex.area[mm]) for mm in mem], axis=1)
            ext[f"patch{r:g}/{label}"] = goldenholz.eq1_snr_db(topo, var)

    # maps on the inflated surface (oct-6 vertices; focal values sampled there)
    hemis = plotting.inflated_views(subject.subjects_dir, "sample", subject.src)
    n_lh = subject.src[0]["nuse"]

    def on_oct(values_valid=None, values_centroid=None):
        out = np.full(len(oct_global), np.nan)
        if values_valid is not None:
            out[oct_ok] = values_valid[col_of[centroids]]
        else:
            out[oct_ok] = values_centroid
        return out

    res = results["probable_default"]
    cmap, norm = plt.get_cmap("viridis"), plt.Normalize(-30, -10)
    rows = [(f"focal 10 nAm\n{cs}", on_oct(values_valid=res[f"focal/model/{cs}"])) for cs in CH_SETS]
    rows += [(f"patch {r:g} mm\n{cs}", on_oct(values_centroid=res[f"patch{r:g}/model/{cs}"])) for r in members for cs in ("mag", "grad")]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C Goldenholz adaptation: Eq. 1 SNR, modelled brain noise "
                               "(skull 0.006 S/m), sample subject", "SNR [dB]", OUT / "Figure_G1C_maps_model.png")
    rows = [(f"focal 10 nAm\n{cs}", on_oct(values_valid=res[f"focal/recorded/{cs}"])) for cs in CH_SETS]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C Goldenholz adaptation: Eq. 1 SNR, recorded noise (task "
                               "baselines, SSP) -- ADAPT", "SNR [dB]", OUT / "Figure_G1C_maps_recorded.png")
    rows = [(f"{lab}\nfocal", on_oct(values_valid=ext[f"focal/{lab}"])) for lab in ("opm", "mag", "grad")]
    rows += [(f"{lab}\npatch 10 mm", on_oct(values_centroid=ext[f"patch10/{lab}"])) for lab in ("opm", "mag", "grad")]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C extension (NEW): modelled brain noise + intrinsic "
                               f"noise (OPM {OPM_ASD * 1e15:g} fT/sqrt(Hz), SQUID brochure), 0.5-100 Hz", "SNR [dB]",
                               OUT / "Figure_G1C_maps_opm_extension.png")

    # conductivity and patch-size effects, distributions
    stats = {}
    for key in res:
        if key in ("source_sd_nAm", "per_type_median_nAm2"):
            continue
        a, b = results["probable_default"][key], results["as_printed"][key]
        stats[key] = dict(median_db=float(np.nanmedian(a)), p5_db=float(np.nanpercentile(a, 5)),
                          p95_db=float(np.nanpercentile(a, 95)),
                          skull_0p06_minus_0p006_db=dict(median=float(np.nanmedian(b - a)),
                                                         max_abs=float(np.nanmax(np.abs(b - a)))))
    for noise_name in ("model", "recorded"):
        for cs in CH_SETS:
            d = res[f"patch16/{noise_name}/{cs}"] - res[f"patch10/{noise_name}/{cs}"]
            stats[f"patch16_minus_patch10/{noise_name}/{cs}"] = dict(median_db=float(np.nanmedian(d)),
                                                                     p5=float(np.nanpercentile(d, 5)), p95=float(np.nanpercentile(d, 95)),
                                                                     flat_area_scaling_db=20 * np.log10(16**2 / 10**2))
    ext_stats = {k: dict(median_db=float(np.nanmedian(v)), p5_db=float(np.nanpercentile(v, 5)), p95_db=float(np.nanpercentile(v, 95)))
                 for k, v in ext.items()}
    for key in ("focal", "patch10", "patch16"):
        for lab in ("mag", "grad"):
            diff = ext[f"{key}/opm"] - ext[f"{key}/{lab}"]
            ext_stats[f"{key}/opm_minus_{lab}"] = dict(median_db=float(np.nanmedian(diff)),
                                                       share_opm_better=float(np.nanmean(diff > 0)))
    summary.update(conductivity_runs={k: dict(source_sd_nAm=v["source_sd_nAm"], per_type_median_nAm2=v["per_type_median_nAm2"])
                                      for k, v in results.items()},
                   distributions=stats, extension_opm=ext_stats, runtime_s=time.time() - t_start,
                   note="Magnetometer, gradiometer and pooled Eq. 1 values are kept separate; the paper reports only MEG "
                        "(pooling unstated) vs EEG. D = SNR_MEG - SNR_EEG is out of scope (no EEG).")
    io.write_json(summary, OUT / "g1c_summary.json")
    with open(OUT / "g1c_oct6_values.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        keys = [k for k in res if "/" in k and not k.startswith("patch")]
        wr.writerow(["hemi", "vertno"] + [f"{k}_dB" for k in keys])
        for g_idx in centroids:
            c = col_of[g_idx]
            wr.writerow([int(cortex.hemi[g_idx]), int(cortex.vertno[g_idx])] + [f"{res[k][c]:.3f}" for k in keys])
    for k in ("focal/model/mag", "focal/model/grad", "focal/model/pooled", "focal/recorded/pooled", "patch10/model/pooled",
              "patch16/model/pooled"):
        print(k, {kk: round(vv, 2) for kk, vv in stats[k].items() if not isinstance(vv, dict)})
    print(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
