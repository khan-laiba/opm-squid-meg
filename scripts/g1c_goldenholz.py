#!/usr/bin/env python3
"""G1C: Goldenholz et al. (2009) cortical SNR maps, MEG part.

Status: ADAPT (MNE sample subject and a task-baseline recording instead of the paper's four
subjects and spontaneous recordings; unstated details chosen in configs/goldenholz_reference.toml).
The OPM maps are a NEW extension: the paper's noise model has no sensor noise, so the extension
calibrates the brain-noise sources on the recorded minus empty-room variance (brain only) and adds
intrinsic SQUID/OPM noise explicitly (OPM noise swept).

Sources are usable vertices only (inside the inner skull and >= 4 mm from its mesh, A-BEM-DIST).

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
from scipy.spatial import cKDTree  # noqa: E402

from opmsquid import anatomy, fullres, g2, goldenholz, io, neuromag, noise, paths, plotting  # noqa: E402

OUT = ROOT / "results" / "g1c"
STATUS = "ADAPT (Goldenholz et al. 2009, MEG part, MNE sample subject); OPM rows NEW"
CH_SETS = ("mag", "grad", "pooled")
OPM_ASD_SWEEP = (7e-15, 10e-15, 15e-15, 20e-15, 30e-15)  # A-OPM-NOISE
OPM_ASD_MAPS = 15e-15
SQUID_ASD = {"mag": 3.5e-15, "grad": 3.6e-13}
PAPER_RANGE_DB = (-29.0, -19.0)  # Fig. 2 colour limits


def _continuous(fname, ch_names, projector, filt):
    raw = mne.io.read_raw_fif(paths.SAMPLE_MEG / fname, preload=True, verbose=False)
    data = projector @ raw.get_data(picks=[raw.ch_names.index(n) for n in ch_names])
    return raw, filt.apply(data - data.mean(axis=1, keepdims=True))


def recorded_noise_variance(ch_names, projector, filt) -> tuple[np.ndarray, int]:
    """Per-channel variance of pre-stimulus baselines (-200 to 0 ms), SSP applied, 0.5-100 Hz
    (filtered as one continuous record, 5-s edges excluded)."""
    raw, data = _continuous(neuromag.RAW_FILE, ch_names, projector, filt)
    events = mne.find_events(raw, stim_channel="STI 014", verbose=False)
    sf = raw.info["sfreq"]
    n_edge = int(5 * sf)
    wins = [(int(e[0] - raw.first_samp - 0.2 * sf), int(e[0] - raw.first_samp)) for e in events]
    seg = np.concatenate([data[:, a:b] for a, b in wins if a > n_edge and b < data.shape[1] - n_edge], axis=1)
    seg = seg - seg.mean(axis=1, keepdims=True)
    return np.mean(seg**2, axis=1), seg.shape[1]


def empty_room_variance(ch_names, projector, filt) -> np.ndarray:
    """Per-channel variance of the empty-room recording, processed like the baselines."""
    raw, data = _continuous(neuromag.ER_FILE, ch_names, projector, filt)
    n_edge = int(5 * raw.info["sfreq"])
    return np.mean(data[:, n_edge:-n_edge] ** 2, axis=1)


def channel_sets(kinds, good):
    return {"mag": good & (kinds == "mag"), "grad": good & (kinds == "grad"), "pooled": good.copy()}


def artefact_diagnostics(g, kinds, cortex, valid_idx, use, noise_cols):
    """Lead-field sanity: magnetometer column energy / median of the 30 nearest usable neighbours
    (BEM artefacts near the inner skull show up as large ratios), and the largest share of the
    modelled noise power held by a single grid source."""
    mags = kinds == "mag"
    e = np.zeros(g.shape[1])
    for s in range(0, g.shape[1], 20000):
        e[s:s + 20000] = np.sum(np.asarray(g[mags, s:s + 20000], np.float64) ** 2, axis=0)
    upos = np.flatnonzero(use)
    _, nb = cKDTree(cortex.rr[valid_idx[upos]]).query(cortex.rr[valid_idx[upos]], k=31)
    ratio = e[upos] / np.median(e[upos][nb[:, 1:]], axis=1)
    share = {}
    for k in ("mag", "grad"):
        p = np.sum(np.asarray(g[kinds == k][:, noise_cols], np.float64) ** 2, axis=0)
        share[k] = float(p.max() / p.sum())
    return dict(usable_columns_energy_ratio_max=float(ratio.max()), usable_columns_ratio_gt_10=int(np.sum(ratio > 10)),
                noise_grid_top_source_power_share=share)


def main():
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    mne.set_log_level("WARNING")
    cfg = tomllib.loads((ROOT / "configs" / "goldenholz_reference.toml").read_text())
    subject = anatomy.load_sample()
    cortex = anatomy.full_resolution(subject)
    valid_idx = np.load(fullres.directory() / "valid_index.npy")
    dig = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    col_of = np.full(cortex.n, -1)
    col_of[valid_idx] = np.arange(len(valid_idx))
    use = cortex.usable[valid_idx]  # usable columns (A-BEM-DIST)
    usable_pos = np.flatnonzero(use)
    info = neuromag.load_info("T3")
    kinds = neuromag.channel_kinds(info)
    raw_info = mne.io.read_info(paths.SAMPLE_MEG / neuromag.RAW_FILE, verbose=False)
    good = np.array([n not in raw_info["bads"] for n in info.ch_names])  # MEG 2443 is bad in the recording
    sets = channel_sets(kinds, good)
    proj, n_proj = neuromag.ssp_projector(raw_info["projs"], info.ch_names)
    filt = noise.AnalysisFilter(fs=raw_info["sfreq"], l_freq=0.5, h_freq=100.0, order=4)
    rec_var, n_rec = recorded_noise_variance(info.ch_names, proj, filt)
    er_var = empty_room_variance(info.ch_names, proj, filt)
    brain_var = rec_var - er_var
    er_share = {k: float(np.median(er_var[sets[k]] / rec_var[sets[k]])) for k in ("mag", "grad")}
    print(f"recorded noise: {n_rec} baseline samples, SSP rank {n_proj}; median RMS mag "
          f"{np.median(np.sqrt(rec_var[sets['mag']])) * 1e15:.1f} fT, grad {np.median(np.sqrt(rec_var[sets['grad']])) * 1e13:.1f} fT/cm; "
          f"empty room / recorded variance {er_share}; excluded {[n for n, ok in zip(info.ch_names, good) if not ok]}")

    rng = np.random.default_rng(2009)
    noise_cols = usable_pos[goldenholz.poisson_disk(cortex.rr[valid_idx[usable_pos]], 0.007, rng)]  # positions in valid_idx
    print(f"noise-model sources: {len(noise_cols)} (7 mm Poisson-disk grid on usable vertices)")
    n_lh_full = int(np.sum(cortex.hemi == 0))
    oct_global = np.concatenate([subject.src[0]["vertno"], subject.src[1]["vertno"] + n_lh_full])
    oct_ok = cortex.usable[oct_global]
    centroids = oct_global[oct_ok]
    density = cfg["sources"]["patch_density_pAm_per_mm2"] * 1e-12 / 1e-6  # A m per m^2
    t0 = time.time()
    members = {r: goldenholz.geodesic_patches(cortex.adjacency, centroids, r * 1e-3, cortex.usable)
               for r in cfg["sources"]["patch_radii_mm"]}
    area = {r: np.array([cortex.area[m].sum() for m in members[r]]) for r in members}
    # review (2026-10-02): the usable rule (A-BEM-DIST) also removes patch members within 4 mm of the inner skull;
    # the same patches over every valid vertex measure that truncation
    members_all = {r: goldenholz.geodesic_patches(cortex.adjacency, centroids, r * 1e-3, cortex.valid) for r in members}
    area_all = {r: np.array([cortex.area[m].sum() for m in members_all[r]]) for r in members}
    loss = {r: 1.0 - area[r] / area_all[r] for r in members}
    print(f"patches built in {time.time() - t0:.0f} s; median members " +
          ", ".join(f"{r:g} mm: {np.median([len(m) for m in members[r]]):.0f}" for r in members))

    summary = dict(status=STATUS,
                   config=cfg, n_valid_vertices=int(len(valid_idx)), n_usable_vertices=int(use.sum()),
                   n_noise_sources=int(len(noise_cols)), n_centroids=int(len(centroids)),
                   channels=dict({cs: int(m.sum()) for cs, m in sets.items()}, excluded=[n for n, ok in zip(info.ch_names, good) if not ok]),
                   recorded_noise=dict(n_samples=int(n_rec), ssp_rank=int(n_proj), empty_room_over_recorded_variance=er_share),
                   patch_area_cm2={f"{r:g}": dict(median=float(np.median(area[r]) * 1e4), flat_disc=float(np.pi * r**2 / 100))
                                   for r in area},
                   patch_truncation={f"{r:g}": dict(median_area_cm2_all_valid=float(np.median(area_all[r]) * 1e4),
                                                    median_area_cm2_usable=float(np.median(area[r]) * 1e4),
                                                    share_losing_over_5pct=float(np.mean(loss[r] > 0.05)),
                                                    share_losing_over_20pct=float(np.mean(loss[r] > 0.20))) for r in members})
    results = {}
    weights = density * cortex.area
    focal_am = cfg["sources"]["focal_nAm"] * 1e-9
    g_sq = p_sq = aat_sq = None
    for bem_label, fname in (("probable_default", "neuromag_bem006"), ("as_printed", "neuromag_bem06")):
        g = proj.astype(np.float32) @ np.asarray(fullres.load(fname, info, subject, cortex)[0])  # checked against this array
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
            g_sq, p_sq, aat_sq = g, topo_p, aat
            summary["lead_field_diagnostics"] = artefact_diagnostics(g, kinds, cortex, valid_idx, use, noise_cols)
            # the truncation's effect on the patch SNR (pooled channels, modelled noise)
            mp = sets["pooled"]
            topo_all = {r: goldenholz.patch_topographies(g, members_all[r], col_of, weights) for r in members}
            for r in members:
                snr_all = goldenholz.eq1_snr_db(topo_all[r][mp], model_var[mp])
                diff = res[f"patch{r:g}/model/pooled"] - snr_all
                summary["patch_truncation"][f"{r:g}"].update(snr_change_db=dict(
                    median=float(np.median(diff)), p5=float(np.percentile(diff, 5)), p95=float(np.percentile(diff, 95)),
                    share_abs_over_1db=float(np.mean(np.abs(diff) > 1.0))))
            # the medial wall (aparc 'unknown'), kept here as in the paper's whole-cortex maps; G2 and G4 exclude it
            names = np.concatenate([plotting.read_freesurfer_annot(subject.labels / f"{h}.aparc.annot") for h in ("lh", "rh")])
            wall = names == "unknown"
            keep = use & ~wall[valid_idx]
            f_all = res["focal/model/pooled"]
            summary["without_medial_wall"] = dict(
                share_of_usable_vertices=float(np.mean(wall[valid_idx][use])), centroids_on_wall=int(wall[centroids].sum()),
                noise_sources_on_wall=int(wall[valid_idx[noise_cols]].sum()),
                focal_model_pooled=dict(median_db=float(np.median(f_all[keep])),
                                        p5_p95_db=[float(x) for x in np.percentile(f_all[keep], [5, 95])],
                                        share_in_paper_range=float(np.mean((f_all[keep] >= PAPER_RANGE_DB[0]) & (f_all[keep] <= PAPER_RANGE_DB[1])))),
                with_wall_focal_model_pooled=dict(median_db=float(np.median(f_all[use])),
                                                  share_in_paper_range=float(np.mean((f_all[use] >= PAPER_RANGE_DB[0]) & (f_all[use] <= PAPER_RANGE_DB[1])))))
            # variant (review, 2026-10-02): without the medial wall and with untruncated patches (members over every valid
            # vertex, also those within 4 mm of the inner skull, where the BEM is less accurate: A-BEM-DIST), next to the
            # primary rows (usable members, wall kept); modelled noise, probable default conductivities
            off_wall = ~wall[centroids]

            def stats_db(x):
                return dict(median_db=float(np.median(x)), p5_db=float(np.percentile(x, 5)), p95_db=float(np.percentile(x, 95)),
                            share_in_paper_range=float(np.mean((x >= PAPER_RANGE_DB[0]) & (x <= PAPER_RANGE_DB[1]))))

            variant = dict(n_centroids=int(off_wall.sum()), n_centroids_primary=int(len(centroids)),
                           note="patches over every valid vertex (not truncated near the inner skull) at centroids off the medial wall")
            for cs in CH_SETS:
                m = sets[cs]
                variant[f"focal/model/{cs}"] = dict(stats_db(res[f"focal/model/{cs}"][keep]),
                                                    primary=stats_db(res[f"focal/model/{cs}"][use]))
                for r in members:
                    snr_v = goldenholz.eq1_snr_db(topo_all[r][m], model_var[m])[off_wall]
                    variant[f"patch{r:g}/model/{cs}"] = dict(stats_db(snr_v), primary=stats_db(res[f"patch{r:g}/model/{cs}"]))
            summary["variant_no_wall_untruncated"] = variant
        else:
            del g

    # NEW extension: brain-only calibration (recorded - empty room), then explicit intrinsic noise
    s_s2_brain, per_type_brain = goldenholz.calibrate_source_variance(brain_var[good], aat_sq[good], kinds[good])
    g_opm = np.asarray(fullres.load("opm_bem006", g2.matched_opm(subject, dig).info, subject, cortex)[0])
    p_opm = {r: goldenholz.patch_topographies(g_opm, mem, col_of, weights) for r, mem in members.items()}
    aat_opm = np.sum(g_opm[:, noise_cols].astype(np.float64) ** 2, axis=1)
    enbw = filt.enbw()
    ext = {}
    setups = [("mag", g_sq, p_sq, aat_sq, sets["mag"], SQUID_ASD["mag"]), ("grad", g_sq, p_sq, aat_sq, sets["grad"], SQUID_ASD["grad"])]
    setups += [(f"opm@{a * 1e15:g}", g_opm, p_opm, aat_opm, np.ones(g_opm.shape[0], bool), a) for a in OPM_ASD_SWEEP]
    for label, gg, pp, aa, m, asd in setups:
        variants = [("", asd**2 * enbw)]
        if not label.startswith("opm@") or asd == OPM_ASD_MAPS:
            variants.append(("brain_only/", 0.0))
        for prefix, intrinsic in variants:
            var = s_s2_brain * aa[m] + intrinsic
            lab = label if not prefix else label.split("@")[0]
            ext[f"{prefix}focal/{lab}"] = goldenholz.eq1_snr_db(gg[m], var, scale=focal_am)
            for r in members:
                ext[f"{prefix}patch{r:g}/{lab}"] = goldenholz.eq1_snr_db(pp[r][m], var)

    # maps on the inflated surface (oct-6 vertices; focal values sampled there), the paper's colour limits
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
    cmap, norm = plt.get_cmap("viridis"), plt.Normalize(*PAPER_RANGE_DB)
    rows = [(f"focal 10 nAm\n{cs}", on_oct(values_valid=res[f"focal/model/{cs}"])) for cs in CH_SETS]
    rows += [(f"patch {r:g} mm\n{cs}", on_oct(values_centroid=res[f"patch{r:g}/model/{cs}"])) for r in members for cs in ("mag", "grad")]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C Goldenholz adaptation: Eq. 1 SNR, modelled brain noise "
                               "(skull 0.006 S/m), sample subject; colour limits of the paper's Fig. 2", "SNR [dB]",
                               OUT / "Figure_G1C_maps_model.png")
    rows = [(f"focal 10 nAm\n{cs}", on_oct(values_valid=res[f"focal/recorded/{cs}"])) for cs in CH_SETS]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C Goldenholz adaptation: Eq. 1 SNR, recorded noise (task "
                               "baselines, SSP) -- ADAPT", "SNR [dB]", OUT / "Figure_G1C_maps_recorded.png")
    opm_lab = f"opm@{OPM_ASD_MAPS * 1e15:g}"
    rows = [(f"{lab.split('@')[0]}\nfocal", on_oct(values_valid=ext[f"focal/{lab}"])) for lab in (opm_lab, "mag", "grad")]
    rows += [(f"{lab.split('@')[0]}\npatch 10 mm", on_oct(values_centroid=ext[f"patch10/{lab}"])) for lab in (opm_lab, "mag", "grad")]
    plotting.cortex_map_figure(hemis, rows, n_lh, cmap, norm, "G1C extension (NEW): brain noise (recorded - empty room) + "
                               f"intrinsic noise (OPM {OPM_ASD_MAPS * 1e15:g} fT/sqrt(Hz), SQUID brochure), 0.5-100 Hz", "SNR [dB]",
                               OUT / "Figure_G1C_maps_opm_extension.png")

    # distributions (focal statistics over usable vertices), conductivity and patch-size effects
    def pick(key, arr):
        return arr[use] if "focal" in key else arr

    stats = {}
    for key in res:
        if key in ("source_sd_nAm", "per_type_median_nAm2"):
            continue
        a, b = pick(key, results["probable_default"][key]), pick(key, results["as_printed"][key])
        stats[key] = dict(median_db=float(np.median(a)), p5_db=float(np.percentile(a, 5)), p95_db=float(np.percentile(a, 95)),
                          skull_0p06_minus_0p006_db=dict(median=float(np.median(b - a)), p95_abs=float(np.percentile(np.abs(b - a), 95)),
                                                         max_abs=float(np.max(np.abs(b - a)))))
    stats["focal/model/pooled"]["oct6_centroids_median_db"] = float(np.median(res["focal/model/pooled"][col_of[centroids]]))
    for noise_name in ("model", "recorded"):
        for cs in CH_SETS:
            d = res[f"patch16/{noise_name}/{cs}"] - res[f"patch10/{noise_name}/{cs}"]
            stats[f"patch16_minus_patch10/{noise_name}/{cs}"] = dict(
                median_db=float(np.median(d)), p5=float(np.percentile(d, 5)), p95=float(np.percentile(d, 95)),
                area_scaling_db=float(20 * np.log10(np.median(area[16.0]) / np.median(area[10.0]))))
    ext_stats = {k: dict(median_db=float(np.median(pick(k, v)))) for k, v in ext.items()}
    for prefix, opm_labels in (("", [f"opm@{a * 1e15:g}" for a in OPM_ASD_SWEEP]), ("brain_only/", ["opm"])):
        for key in ("focal", "patch10", "patch16"):
            for o in opm_labels:
                for lab in ("mag", "grad"):
                    diff = pick(key, ext[f"{prefix}{key}/{o}"] - ext[f"{prefix}{key}/{lab}"])
                    ext_stats[f"{prefix}{key}/{o}_minus_{lab}"] = dict(median_db=float(np.median(diff)),
                                                                       share_opm_better=float(np.mean(diff > 0)))
    # soft comparisons with the paper (docs/literature/goldenholz2009.md): source SD 1.6-1.9 nAm,
    # Fig. 2 display range -29 ... -19 dB (clipping limits, not a reported range), dark medial MEG maps
    # [inference], 3 vs 8 cm^2 patches: 10 dB in the mesial temporal lobe (modality unclear)
    names = np.concatenate([plotting.read_freesurfer_annot(paths.SUBJECTS_DIR / "sample" / "label" / f"{h}.aparc.annot")
                            for h in ("lh", "rh")])
    region_v = names[valid_idx[use]]
    lobe_v, lobe_c = plotting.lobe_of(region_v), plotting.lobe_of(names[centroids])
    mesial_c = np.isin(names[centroids], ["entorhinal", "parahippocampal"])
    deep_medial = ("rostralanteriorcingulate", "caudalanteriorcingulate", "posteriorcingulate", "isthmuscingulate",
                   "parahippocampal", "entorhinal", "medialorbitofrontal")
    # s_s scales as 1/sqrt(number of noise sources) under the calibration rule; a 7-mm MNE surface grid
    # holds ~4000 sources (docs/literature/goldenholz2009.md, A3)
    checks = dict(source_sd_nAm={k: v["source_sd_nAm"] for k, v in results.items()}, paper_source_sd_nAm=[1.6, 1.9],
                  source_sd_nAm_at_4000_sources={k: v["source_sd_nAm"] * float(np.sqrt(len(noise_cols) / 4000))
                                                 for k, v in results.items()},
                  extension_source_sd_nAm_brain_only=float(np.sqrt(s_s2_brain) * 1e9))
    for bem_label, r_ in results.items():
        for cs in CH_SETS:
            f = r_[f"focal/model/{cs}"][use]
            checks[f"{bem_label}/focal/model/{cs}"] = dict(
                share_below_m29=float(np.mean(f < -29)), share_m29_to_m19=float(np.mean((f >= -29) & (f <= -19))),
                share_above_m19=float(np.mean(f > -19)),
                median_by_lobe={lb: float(np.median(f[lobe_v == lb])) for lb in plotting.DK_LOBES},
                share_below_m29_deep_medial_regions=float(np.mean(f[np.isin(region_v, deep_medial)] < -29)))
            d = r_[f"patch16/model/{cs}"] - r_[f"patch10/model/{cs}"]
            checks[f"{bem_label}/patch16_minus_patch10/model/{cs}"] = dict(
                mesial_temporal_median_db=float(np.median(d[mesial_c])), n_mesial_temporal=int(mesial_c.sum()),
                all_median_db=float(np.median(d)), by_lobe={lb: float(np.median(d[lobe_c == lb])) for lb in plotting.DK_LOBES},
                paper_mesial_temporal_db=10.0, nominal_area_scaling_db=float(20 * np.log10(8 / 3)))
    summary.update(conductivity_runs={k: dict(source_sd_nAm=v["source_sd_nAm"], per_type_median_nAm2=v["per_type_median_nAm2"])
                                      for k, v in results.items()},
                   distributions=stats, comparison_with_paper=checks, extension_opm=ext_stats,
                   extension_note=("brain-noise sources calibrated on recorded minus empty-room variance (same SSP and band); "
                                   "intrinsic noise added explicitly; 'brain_only/' rows omit it"),
                   runtime_s=time.time() - t_start,
                   note="Magnetometer, gradiometer and pooled Eq. 1 values are kept separate; the paper reports only MEG "
                        "(pooling unstated) vs EEG. D = SNR_MEG - SNR_EEG is out of scope (no EEG).")
    io.write_json(summary, OUT / "g1c_summary.json")
    with open(OUT / "g1c_oct6_values.csv", "w", newline="") as fh:
        io.csv_status(fh, STATUS)
        wr = csv.writer(fh)
        keys = [k for k in res if "/" in k and not k.startswith("patch")]
        wr.writerow(["hemi", "vertno"] + [f"{k}_dB" for k in keys])
        for g_idx in centroids:
            c = col_of[g_idx]
            wr.writerow([int(cortex.hemi[g_idx]), int(cortex.vertno[g_idx])] + [f"{res[k][c]:.3f}" for k in keys])
    for k in ("focal/model/mag", "focal/model/grad", "focal/model/pooled", "focal/recorded/pooled", "patch10/model/pooled",
              "patch16/model/pooled"):
        print(k, {kk: round(vv, 2) for kk, vv in stats[k].items() if not isinstance(vv, dict)})
    print("diagnostics", summary["lead_field_diagnostics"])
    for k, v in ext_stats.items():
        if "minus" in k and k.split("/")[-2] in ("focal",) or (k.startswith("brain_only/focal") and "minus" in k):
            print(" ", k, {kk: round(vv, 2) for kk, vv in v.items()})
    print(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
