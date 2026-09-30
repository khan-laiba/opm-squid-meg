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


def channel_sets(kinds, good):
    return {"mag": good & (kinds == "mag"), "grad": good & (kinds == "grad"), "pooled": good.copy()}


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
    raw_info = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    good = np.array([n not in raw_info["bads"] for n in info.ch_names])  # MEG 2443 is bad in the recording
    sets = channel_sets(kinds, good)
    proj, n_proj = neuromag.ssp_projector(raw_info["projs"], info.ch_names)
    filt = noise.AnalysisFilter(fs=raw_info["sfreq"], l_freq=0.5, h_freq=100.0, order=4)
    rec_var, n_rec = recorded_noise_variance(info.ch_names, proj, filt)
    print(f"recorded noise: {n_rec} baseline samples, SSP rank {n_proj}; median RMS mag "
          f"{np.median(np.sqrt(rec_var[sets['mag']])) * 1e15:.1f} fT, grad {np.median(np.sqrt(rec_var[sets['grad']])) * 1e13:.1f} fT/cm; "
          f"excluded bad channels: {[n for n, ok in zip(info.ch_names, good) if not ok]}")

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
                                n_centroids=int(len(centroids)), channels=dict({cs: int(m.sum()) for cs, m in sets.items()},
                                excluded=[n for n, ok in zip(info.ch_names, good) if not ok]),
                                recorded_noise=dict(n_samples=int(n_rec), ssp_rank=int(n_proj)))
    weights = density * cortex.area
    focal_am = cfg["sources"]["focal_nAm"] * 1e-9
    g_sq = p_sq = None
    for bem_label, fname in (("probable_default", "neuromag_bem006"), ("as_printed", "neuromag_bem06")):
        g = proj.astype(np.float32) @ np.load(paths.CACHE / "fullres" / f"{fname}.npy")
        aat = np.sum(g[:, noise_cols].astype(np.float64) ** 2, axis=1)
        s_s2, per_type = goldenholz.calibrate_source_variance(rec_var[good], aat[good], kinds[good])
        model_var = s_s2 * aat
        topo_p = {r: goldenholz.patch_topographies(g, mem, col_of, weights) for r, mem in members.items()}
        res = dict(source_sd_nAm=float(np.sqrt(s_s2) * 1e9), per_type_median_nAm2=per_type)
        for noise_name, var in (("model", model_var), ("recorded", rec_var)):
            for cs in CH_SETS:
                m = sets[cs]
                res[f"focal/{noise_name}/{cs}"] = goldenholz.eq1_snr_db(g[m], var[m], scale=focal_am)
                for r in members:
                    res[f"patch{r:g}/{noise_name}/{cs}"] = goldenholz.eq1_snr_db(topo_p[r][m], var[m])
        results[bem_label] = res
        print(f"{bem_label}: noise-source SD {res['source_sd_nAm']:.2f} nAm (paper: 1.6-1.9 nAm per source)")
        if bem_label == "probable_default":
            g_sq, p_sq = g, topo_p
        else:
            del g

    # NEW extension: OPM with the same calibrated brain-noise model plus intrinsic noise (and SQUIDs with brochure noise)
    g_opm = np.load(paths.CACHE / "fullres" / "opm_bem006.npy")
    p_opm = {r: goldenholz.patch_topographies(g_opm, mem, col_of, weights) for r, mem in members.items()}
    s_s2 = results["probable_default"]["source_sd_nAm"] ** 2 * 1e-18
    enbw = filt.enbw()
    ext = {}
    for label, gg, pp, m in (("opm", g_opm, p_opm, np.ones(g_opm.shape[0], bool)), ("mag", g_sq, p_sq, sets["mag"]),
                             ("grad", g_sq, p_sq, sets["grad"])):
        asd = OPM_ASD if label == "opm" else SQUID_ASD[label]
        var = s_s2 * np.sum(gg[m][:, noise_cols].astype(np.float64) ** 2, axis=1) + asd**2 * enbw
        ext[f"focal/{label}"] = goldenholz.eq1_snr_db(gg[m], var, scale=10e-9)
        for r in members:
            ext[f"patch{r:g}/{label}"] = goldenholz.eq1_snr_db(pp[r][m], var)

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
                                                         p95_abs=float(np.nanpercentile(np.abs(b - a), 95)),
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
    # soft comparisons with the paper (docs/literature/goldenholz2009.md): source SD 1.6-1.9 nAm,
    # Fig. 2 display range -29 ... -19 dB (clipping limits, not a reported range), medial MEG maps
    # largely below -29 dB [inference], 3 vs 8 cm^2 patches: 10 dB in the mesial temporal lobe
    names = np.concatenate([plotting.read_freesurfer_annot(paths.SUBJECTS_DIR / "sample" / "label" / f"{h}.aparc.annot")
                            for h in ("lh", "rh")])
    lobe_v, lobe_c = plotting.lobe_of(names[valid_idx]), plotting.lobe_of(names[centroids])
    mesial_c = np.isin(names[centroids], ["entorhinal", "parahippocampal"])
    # s_s scales as 1/sqrt(number of noise sources) under the calibration rule; the paper's ~7-mm grid
    # held ~4000 sources (docs/literature/goldenholz2009.md, A3)
    checks = dict(source_sd_nAm={k: v["source_sd_nAm"] for k, v in results.items()}, paper_source_sd_nAm=[1.6, 1.9],
                  source_sd_nAm_at_4000_sources={k: v["source_sd_nAm"] * float(np.sqrt(len(noise_cols) / 4000))
                                                 for k, v in results.items()})
    for bem_label, r_ in results.items():
        for cs in CH_SETS:
            f = r_[f"focal/model/{cs}"]
            checks[f"{bem_label}/focal/model/{cs}"] = dict(
                share_below_m29=float(np.mean(f < -29)), share_m29_to_m19=float(np.mean((f >= -29) & (f <= -19))),
                share_above_m19=float(np.mean(f > -19)),
                median_by_lobe={lb: float(np.median(f[lobe_v == lb])) for lb in plotting.DK_LOBES})
            d = r_[f"patch16/model/{cs}"] - r_[f"patch10/model/{cs}"]
            checks[f"{bem_label}/patch16_minus_patch10/model/{cs}"] = dict(
                mesial_temporal_median_db=float(np.median(d[mesial_c])), n_mesial_temporal=int(mesial_c.sum()),
                all_median_db=float(np.median(d)), by_lobe={lb: float(np.median(d[lobe_c == lb])) for lb in plotting.DK_LOBES},
                paper_mesial_temporal_db=10.0, nominal_area_scaling_db=float(20 * np.log10(8 / 3)))
    summary.update(conductivity_runs={k: dict(source_sd_nAm=v["source_sd_nAm"], per_type_median_nAm2=v["per_type_median_nAm2"])
                                      for k, v in results.items()},
                   distributions=stats, comparison_with_paper=checks, extension_opm=ext_stats, runtime_s=time.time() - t_start,
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
