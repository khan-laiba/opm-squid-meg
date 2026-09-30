# Hunold et al. (2016): methods and parameter extraction for an MEG reproduction

**Source file:** `Hunold_2016_Physiol._Meas._37_1146.pdf`. The PDF has 18 pages. PDF page 1 is an IOP download cover sheet with no article content, and PDF page *n* is published page 1144 + *n*. All article pages were read: pp. 1146–1162, all captions, Table 1 and Figs 1–6. The figures were also inspected at their native raster resolution (200 ppi).

**Supplement:** The paper cites supplementary data at stacks.iop.org/PM/37/1146/mmedia (p. 1147). The only item it references there is "figure 7", the 64-electrode EEG result (p. 1158). The supplement is not available locally and was not consulted. Going by the main text, it holds EEG-only content.

**Purpose of this note:** to reproduce the MEG parts of the study in MNE-Python: SQUID magnetometers vs planar gradiometers now, OPMs later. EEG is out of scope and is recorded only briefly, where it shares definitions with MEG. All page numbers are published page numbers.

**Evidence tags**

| Tag | Meaning |
|---|---|
| **[P]** | Printed in the paper. Paraphrased, except for one marked quote. |
| **[F]** | Read directly from a figure (printed labels, axes, numbers). |
| **[Dg]** | Digitized from a figure raster. Approximate; resolution is stated. |
| **[D]** | Derived by us arithmetically from printed values. |
| **[NS]** | Not stated anywhere in the paper. |
| **[MNE]** | Our implementation note for the re-implementation, not from the paper. |

The paper contains **no numbered or displayed equations**. Every formula below is our transcription of a verbal definition and is labelled as such.

---

## 1. Bibliographic information and study aim

- **Citation.** Hunold A, Funke M E, Eichardt R, Stenroos M, Haueisen J (2016). EEG and MEG: sensitivity to epileptic spike activity as function of source orientation and depth. *Physiol. Meas.* **37**(7) 1146–1162. doi:10.1088/0967-3334/37/7/1146 (p. 1146).
- **Dates.** Received 21 Feb 2016, revised 10 May 2016, accepted 20 May 2016, published 22 Jun 2016 (p. 1146).
- **Affiliations** (p. 1146):
  - TU Ilmenau (Institute of Biomedical Engineering and Informatics)
  - University of Texas at Houston, McGovern Medical School (child neurology)
  - Aalto University (Neuroscience and Biomedical Engineering)
  - Jena University Hospital (Biomagnetic Center)
- **Funding.** BMBF 03IPT605A, DAAD 57163007, Academy of Finland 290018 (p. 1160).
- **Keywords.** EEG, MEG, epilepsy, dipole, cortical patch, forward modeling, interictal spike (p. 1147).
- **Aim (paraphrased).** The study simulates interictal spikes from cortically constrained single dipoles and ~20 mm² patches in two individual three-shell BEM models, adding a distributed EEG-like background. It asks how source depth (below the scalp) and orientation (relative to the local inner-skull normal) set the spike SNR in the best single channel of 128-channel EEG, 102 Vectorview magnetometers and 204 gradiometers (pp. 1146, 1148). These are the Vectorview planar gradiometers; the paper itself never uses the word "planar". It extends Haueisen et al (2012) with a more realistic background model and an SNR evaluated per channel, as in clinical review (p. 1148).

---

## 2. Anatomy and head model

| Item | Extraction | Page |
|---|---|---|
| Subjects | Two healthy male volunteers, "AH" (22 y) and "MF" (49 y). Verbal informed consent. | 1148 |
| MRI scanner | 3T MAGNETOM Trio, A Tim System (Siemens, Erlangen) | 1148 |
| MRI sequences | One 3D single-echo MPRAGE, plus two 3D multi-echo FLASH with flip angles 30° and 5° (after Fischl et al 2004) | 1148 |
| Sequence parameters (TR, TE, voxel size, FLASH echoes) | [NS] | – |
| Cortical surfaces | FreeSurfer, from MPRAGE. Triangulated meshes of approx. 300 000 elements. | 1148 |
| FreeSurfer / MNE / Matlab versions | [NS]. The cited MNE manual is the User's Guide 2.7.3 (Hämäläinen 2010). | 1148, 1161 |
| Inner skull | Watershed algorithm (Ségonne et al 2004) on MPRAGE | 1148 |
| Outer skull and scalp | Segmented from the FLASH images (Hämäläinen 2010) | 1148 |
| Segmentation software | MNE software (Gramfort et al 2014) | 1148 |
| BEM surface resolution | Each boundary surface downsampled to 5120 triangles (after Haueisen et al 1997, Stenroos et al 2014) | 1148 |
| Scalp mesh | Typical triangle side length 6 mm | 1159 |
| Inner-skull mesh smoothness | Largest angle between normals of neighbouring nodes: 4.7 ± 2.5° (mean ± SD) | 1159 |
| Compartments | Three, individual per subject: brain, skull, scalp | 1151 |
| Conductivities | Brain 0.33 S m⁻¹, skull 0.0042 S m⁻¹, scalp 0.33 S m⁻¹ (cites Geddes and Baker 1967) | 1151 |
| Conductivity ratio | [D] 0.33/0.0042 ≈ 78.6, i.e. about 1:80. Ratios 1:50 and 1:25 were also simulated and gave similar results according to the authors; no numbers are reported. | 1151, 1159 |
| Numerical method | BEM; Galerkin formulation with linear basis functions; isolated-source approach (Hämäläinen and Sarvas 1989; Stenroos and Sarvas 2012) | 1151 |
| Implementation | Custom Matlab scripts based on the Helsinki BEM library (Stenroos et al 2007). All computations after segmentation were done in Matlab. | 1148, 1151 |
| Cortical surface area | AH: 194 980 mm². MF: [NS]. | 1150 |

**[MNE] notes**

- The closest equivalent model is `mne.make_bem_model(subject, ico=4, conductivity=(0.33, 0.0042, 0.33))`, which gives 5120 triangles per surface.
- MNE solves the BEM by linear collocation with the isolated-skull approach, not by Galerkin. Expect small numerical differences, mostly in EEG. MEG is only weakly affected by skull conductivity.
- Without FLASH data, `mne watershed_bem` produces all three surfaces from the T1. This is a documented deviation.

---

## 3. Source space

### 3.1 Locations, orientation constraint, mesh

- **Locations.** Sources sit at the mesh nodes of the FreeSurfer white–gray matter boundary, i.e. the white surface (p. 1148). That this is the full-resolution mesh (≈300 000 elements) is our inference [D] (pp. 1148, 1156).
- **Orientation constraint.** Every dipole (spike, patch and background) is fixed perpendicular to the cortical surface (p. 1148).
  - The polarity convention (outward or inward normal) is [NS].
  - How node normals are computed is [NS].
- **Hemispheres.** Not stated explicitly [NS]. Fig. 1 shows only the left hemisphere (p. 1149). The ~300 000-node count and the 194 980 mm² area suggest both hemispheres [D].
- **Implied mesh density [D]** (from pp. 1148, 1150–1151):
  - AH has 29 823 background dipoles, which is 10% of the nodes, so about 298 000 nodes.
  - Mean area per node is about 0.654 mm², matching the stated 6.54 mm² per background dipole ÷ 10.
  - The Table 1 patches (19–35 dipoles in 20.1–20.5 mm², p. 1156) confirm that patches were built on the full-resolution mesh.

### 3.2 Depth

- For each cortical node, the closest node of the scalp mesh was found. Depth is the Euclidean distance between the two (p. 1148).
- The scalp mesh is the skin surface. That the 5120-triangle BEM scalp mesh (every boundary surface was downsampled to 5120 triangles) is the one used for depth is our inference [D] (pp. 1148, 1159). The abstract describes depths as below the skin surface (p. 1146).
- Depth is therefore measured to the skin, not to the sensors and not to the inner skull.
- Depth is measured **node to node**, not node to nearest surface point. The authors say this can overestimate depth but call the error negligible, given ~6 mm scalp triangles and depths of 20 mm or more (p. 1159).
  - [D] For a lateral offset of ≤3.5 mm at 20 mm depth, the overestimate is ≲0.3 mm.

### 3.3 Orientation ("nearby" inner-skull normal)

- For each cortical node, the closest node of the 5120-triangle inner-skull mesh was found. This is a separate nearest-node search from the one used for depth (p. 1148).
- Orientation is the angle between the cortical node's normal and that inner-skull node's normal (p. 1148). The angles range from 0° (radial) to 90° (tangential) (p. 1148).
- "Nearby" therefore means the nearest inner-skull mesh node by distance. The paper speaks of the closest node (p. 1148); that the distance is Euclidean is our reading.
- The authors bound the angular error by the maximal angle between neighbouring inner-skull node normals, 4.7 ± 2.5° (p. 1159).
- How angles above 90° are handled (folding with |cos| or exclusion) is [NS]. The paper only states that the angles lie within 0–90°.
- Orientation is a scalar angle. The polarity of the normals does not enter the binning.

### 3.4 Spike-source sampling: "dipole traces" (pp. 1148–1149)

- **Goal.** Sample all brain regions and cover all depth/orientation combinations evenly (p. 1148).
- **Trace.** Each trace has 20 dipoles running from the bottom of a sulcus to the crown of the adjacent gyrus. Orientation goes radial → tangential → radial along the trace (p. 1148).
- **Dipole 1 (seed).** A radial dipole (< 6°) at a sulcal bottom, chosen so that all its neighbours are more tangential and more superficial (p. 1149).
- **Dipoles 2–9.** Eight dipoles are added, each a **direct mesh neighbour** of its predecessor. Their orientations fall in successive 10° "decades": dipole 1 in 0–10°, dipole 2 in 10–20°, …, dipole 9 in 80–90° (p. 1149).
- **Dipole 20 (end point).** Starting from dipole 9 on the sulcal wall, the end point is the closest dipole with the most radial orientation on the corresponding gyral crown (p. 1149). The distance metric and the region searched for the most radial dipole are [NS].
- **Dipoles 10–19.** Placed between dipoles 9 and 20 so as to best match all orientation decades (p. 1149). The algorithm and whether neighbour adjacency is required are [NS].
- Near-radial orientations on crowns are sampled less evenly than at sulcal bottoms (p. 1149).
- **Counts.** 127 traces for AH and 128 for MF (p. 1149). [D] 255 traces × 20 dipoles = **5100** dipoles. The abstract says **5600** (p. 1146); see §11.

### 3.5 Bins

- **Depth.** Eight evenly spaced intervals between 20 and 60 mm, i.e. 5 mm wide [D]. **Orientation.** Nine 10° intervals between 0 and 90° (p. 1150).
- Patches are binned by the **mean depth** and **mean orientation** of their dipoles (pp. 1150, 1158).
- **Displayed grids:**
  - Dipoles: 20–60 mm × 0–90°, 8 × 9 bins.
  - Patches: 20–55 mm × 10–90°, 7 × 8 bins. The 0–10° column and the 55–60 mm row were dropped for too few sources (Fig. 3 caption p. 1152; p. 1155).
- The analysis is described as N depth bins × M orientation bins (p. 1151).
- Bin-edge convention (which side is closed) and handling of values exactly on an edge: [NS].

**Depth labels** (p. 1154): "superficial" = 20–35 mm, "intermediate" = 35–45 mm, "deep" = 45–60 mm. These labels are used inconsistently elsewhere; see §11.

### 3.6 Exclusions and analysed numbers

- Sources **deeper than 60 mm** were simulated but not shown, because of small sample sizes. This explains the gap between modelled and analysed numbers (p. 1153).
- Treatment of sources shallower than 20 mm: [NS].
- Patches that could not reach 20 mm² were excluded (pp. 1149–1150).
- **Totals:**
  - Modelled: 5600 dipoles (abstract, p. 1146; the methods imply 5100 [D]) and 3300 patches ≥ 20 mm² (pp. 1146, 1150).
  - Analysed: **3783 dipoles** (20–60 mm, 0–90°) and **2264 patches** (20–55 mm, 10–90°) (p. 1151).
  - [F] The per-bin counts in Fig. 3(c),(d) sum exactly to 3783 and 2264 (Appendix B).
- Both subjects appear to be pooled in all bins. This is implied by the totals, not stated [NS].
- **Per-bin sample sizes as printed** (p. 1153): dipoles 13–117 (48.6 ± 30.6); patches 2–163 (37.2 ± 34.8). Only four bins have fewer than 10 sources: patches at 20–25 mm and 30–70°.
- **Per-bin sample sizes from Fig. 3(c),(d)** [F, D]:
  - Dipoles: range 15–177, mean 52.5, SD 30.0 (population) / 30.2 (sample). This does **not** match the text.
  - Patches: range 2–163 (matches), mean 40.4, SD 35.7 / 36.0 (the text gives 37.2 ± 34.8).

---

## 4. Source models

### 4.1 Focal dipoles

- **Strength: 600 nAm** (p. 1151). It was chosen to roughly match the patch strengths (p. 1151).
- [F] The spike waveform in Fig. 2(b) peaks at exactly 600 nAm (p. 1150). "Strength" is therefore the **peak moment of the spike waveform**.
- Orientation is the cortical normal (p. 1148).

### 4.2 Patches (pp. 1148–1151, 1156, 1158–1160)

- **Rationale.** Brain activity is spatially extended. Small patches keep depth and orientation well defined and make cancellation less likely (p. 1148). They were deliberately small to minimize cancellation; larger patches are left for future work (pp. 1159–1160).
- **Seeds.** Every dipole of every trace seeds one patch (p. 1149).
- **Growth.** Neighbouring dipoles are added outward from the seed (following Lütkenhöner et al 1995 and Hillebrand and Barnes 2002) until the target area is reached (p. 1149).
- **Orientation criterion.** An added dipole's orientation angle must be within ±10° of the seed's (Wagner et al 2000; Fuchs et al 2007) (p. 1149).
  - Seed < 10°: admissible range 0–20°.
  - Seed > 80°: admissible range 70–90°.
  - This criterion is on the scalar depth-independent angle (§3.3), not on the 3D normal vector.
- **Stopping rule.** Area threshold 20 mm². Growth stops **the first time the area exceeds 20 mm²** (p. 1149). The patches are therefore ≥ 20 mm² in **area**, not 20 mm in radius. The Table 1 examples are 20.1–20.5 mm² with 19–35 dipoles (p. 1156).
- **Exclusion.** If no admissible neighbours remain before 20 mm² is reached, the smaller patch is dropped. 3300 patches ≥ 20 mm² remained (pp. 1149–1150).
- **Depth and orientation** of a patch are the means over its dipoles (pp. 1150, 1158). Table 1 gives the mean ± SD (p. 1156).
- **Strength convention** (p. 1151):
  - Each patch's represented area was taken into account individually.
  - Patch strength was scaled to about 600 nAm, following Murakami and Okada (2006).
  - A single dipole-moment **density** was set so that patch strengths ranged **612–678 nAm, median 622 nAm**.
  - Activation is **uniform** over the patch (Badier et al 2007) (p. 1151). Every patch dipole carries the same spike waveform (p. 1150).
  - [D] The implied density is about 30–31 nAm/mm² (612 nAm ÷ ≥ 20 mm²; the median 622 nAm corresponds to about 20.4 mm² at 30.5 nAm/mm²). The exact density is [NS].
  - [NS] Whether each dipole's moment is density × the area of its node, and how node or patch area is computed (vertex area vs triangles), is not stated.
- **What "strength" means.** The paper never defines patch "strength". [D] Reading 612–678 nAm as the sum of the dipole moment magnitudes fits p. 1158 and the Table 1 areas. The dipoles follow the local normals, so a patch's net (vector) "total strength" was often **smaller** than a single 600 nAm dipole, giving lower patch SNRs (p. 1158).
- **Conflicting density value.** The Discussion describes triangle-mesh patch designs using a preset current density of 100 nA mm⁻² (value as printed; Hämäläinen et al 1993) (p. 1159). Grammatically this sentence refers to earlier approaches, but it could be read as describing this study too. If so, it conflicts with the implied ~30 nAm/mm², and the units differ; see §11.
- **[NS] details:** growth order (breadth-first by adjacency or by distance), tie-breaking, and whether patches may overlap background dipoles.

### 4.3 Background dipoles

- Each background dipole has a maximal strength of 10 nAm (p. 1151). See §6.

### 4.4 Normalization summary

| Source | Normalization |
|---|---|
| Spike dipole | Absolute peak 600 nAm |
| Patch | Fixed moment density; totals 612–678 nAm |
| Background | Each realization scaled to the range ±10 nAm (p. 1150) |

No other normalization is described: nothing rescales to a target SNR, and nothing normalizes per sensor type.

---

## 5. Spike waveform

- **Construction.** A combination of several sinusoidal waves imitating an interictal spike complex (cites Niedermeyer and Lopes da Silva 2005) (p. 1150). Fig. 2(b) is titled an interictal spike-wave complex (p. 1150).
- **Equation and parameters (frequencies, phases, amplitudes, windowing): [NS].**
- **Duration.** The spike lasts approx. 80 ms (p. 1150).
- **Assignment.** The same waveform goes either to one dipole or to all dipoles of one patch (p. 1150).
- **Amplitude.** Peak 600 nAm for dipoles (p. 1151; Fig. 2(b) p. 1150). Patches are scaled by density to 612–678 nAm total (p. 1151).
- **Sampling.** 1000 samples/s is stated for the background (p. 1150). The spike is presumably sampled identically [NS].
- **Time axis.**
  - Fig. 2(b) spans −0.1 to 0.3 s, with onset at 0 s (p. 1150).
  - Fig. 2(d) spans −0.5 to 0.3 s (p. 1150).
  - Fig. 6 shows 0–2 s windows with the spike peak at about 0.93 s [Dg] (p. 1157).
  - The onset position within the 6 s simulated epoch is [NS]. It must leave at least 1 s before onset for the SNR baseline (p. 1151).
- **Digitized shape [Dg].** From Fig. 2(b); resolution about 2.2 ms × 4.2 nAm per pixel, uncertainty about ±2 ms and ±5 nAm.

| Feature | Time (ms after onset) | Amplitude (nAm) |
|---|---|---|
| Onset (departure from 0) | 0 | 0 |
| Negative pre-deflection | ≈10 | ≈ −70 to −75 |
| Zero crossing | ≈13 | 0 |
| **Main positive peak** | ≈24 | **600**; FWHM ≈17 ms |
| Zero crossing | ≈39 | 0 |
| Main negative trough | ≈54 | ≈ −254 |
| Zero crossing (end of sharp spike) | ≈74 | 0 |
| Slow wave, positive | ≈87 | ≈ +149 |
| Zero crossing | ≈123 | 0 |
| Slow wave, negative | ≈150 | ≈ −91 |
| Return to baseline | ≈198 (≈0.2 s) | 0 |

- **Interpreting the 80 ms.** The sharp spike, from onset to the end of the main negative phase, lasts about 74 ms, consistent with the stated ≈80 ms. The whole spike-wave complex lasts about 200 ms. [D] Peak-to-peak / peak ≈ (600 + 254)/600 ≈ 1.42.

---

## 6. Background ("brain noise") model

- **Number and placement.**
  - A random 10% of the cortical nodes, approx. 30 000 (p. 1148). AH has **29 823** (p. 1150); the MF count is [NS].
  - Each dipole represents on average 6.54 mm² (AH) or 6.44 mm² (MF) (pp. 1150–1151).
  - Orientation is the cortical normal (p. 1148).
  - The random-selection procedure and seed, and whether spike-source nodes are excluded, are [NS].
- **Time courses** (p. 1150):
  - Stationary Gaussian noise is filtered into the five conventional EEG bands.
  - The bands are weighted with physiologically plausible values between 0.4 and 0.6, using the asymmetry ratio of Gotman et al (1973).
  - The band signals are summed.
  - The authors state that the result reproduces the spectrum of experimental EEG (Niedermeyer and Lopes da Silva 2005) and yields EEG-like electrode signals.
  - [NS]: band edges, filter type/order/phase, the weight of each band, the definition of the asymmetry ratio (given by citation only), and whether normalization is applied per band or after summation.
- **Duration and sampling.** Six seconds of filtered white noise per band at 1000 samples/s (p. 1150).
- **Amplitude distribution.**
  - Gaussian (stationary) before normalization. The magnitude is normalized to ±10 nAm, and the maximal background dipole strength is 10 nAm (after Hämäläinen et al 1993) (pp. 1150–1151).
  - [NS] Whether this is peak normalization (max|s| = 10 nAm) or a min–max rescaling to [−10, 10] nAm.
  - [NS] Whether amplitudes are weighted by each dipole's represented area. Nothing indicates area weighting.
- **Spatial correlation.** None at source level: each background dipole gets its own independent realization (p. 1150). The authors list this as a limitation and suggest resting-state-network correlations (p. 1160). They also note that the frequency bands are not weighted by region (p. 1160).
- **Realizations across spike simulations: [NS].**
  - [Dg] In Fig. 6, traces from different example sources on the **same channel** have practically identical baselines (correlation of digitized baselines r ≥ 0.97, mostly 0.99–1.00), and the spike appears at the same latency (p. 1157).
  - At least for these examples, then, one fixed background realization and onset time were reused across source simulations.
- **Scaling relative to the spike.**
  - Both are absolute: a 600 nAm spike peak (or 612–678 nAm patch) against ≈30 000 background dipoles of ≤ 10 nAm each. No rescaling to a target SNR (pp. 1150–1151).
  - A variant used a background level of **300%** (p. 1158). Whether this is ×3 in amplitude or in power is [NS].
- **Stated rationale** (p. 1160):
  - Closer to reality than white noise.
  - Unlike measured noise, it gives EEG and MEG an adequate share of noise, so SNR estimates are not biased between modalities.

---

## 7. Sensors

**MEG**

- **System.** Elekta Neuromag **Vectorview** (Elekta Neuromag Oy, Helsinki), simulated at **102 magnetometer** and **204 gradiometer** positions (pp. 1146, 1151). The Vectorview gradiometers are planar; the paper does not use that word. TRIUX is not mentioned.
- **Coil model.** The extent of each MEG sensor is modelled by **4-point numerical integration** (p. 1151).
  - Coil type (T1/T2/T3), coil size and gradiometer baseline: [NS].
  - [MNE] In MNE 1.13.2 `coil_def.dat`, `accuracy='normal'` uses 4 integration points for Vectorview magnetometers 3022/3023/3024 and for planar gradiometers 3012/3013/3014 (2 per loop); `'accurate'` uses 16 and 8.
  - [MNE] T1 magnetometers are 25.8 mm, T3 magnetometers 21.0 mm; gradiometers are 26.39 mm with a 16.8 mm baseline. `'normal'` matches the paper.
- **Positions and coregistration.** Helmet position relative to each head, device-to-head transform, and how the array was placed on the two MRI heads: all [NS]. Fig. 2(c) only shows sensors around a head model (p. 1150).
- **Channels used.** [F] Fig. 6 labels individual Vectorview channels: magnetometers 0631, 0711, 0741; gradiometers 0412, 0413, 0423 (p. 1157). SNR is evaluated on **single gradiometer channels**, not on combined planar pairs.

**EEG (brief)**

- Standard 128-channel ANT WaveGuard cap, 10–5 layout including the subtemporal electrode chain (p. 1151).
- A 64-electrode 10–10 variant is used in the Discussion and supplement (p. 1158).
- The EEG reference is [NS].
- Fig. 6 electrodes: FC3, CCP5h, CP5, C5 (p. 1157).

**Sensor noise (verified)**

- The authors state that they **omitted sensor noise** in the simulations (p. 1158). They argue that biological noise usually dominates and that sensor noise differs between sensor types (p. 1158).
- Environmental or interference noise is not mentioned; none is modelled.
- The reference simulations are therefore **brain-noise-only**, which confirms the GOAL.md statement.
- The authors suggest the framework can characterize new sensors, such as atomic magnetometers, *provided* their sensor noise is known (p. 1158).

---

## 8. SNR definition

### 8.1 As printed (p. 1151, paraphrased except one quote)

1. **Motivation.** The SNR is meant to mirror clinicians' visual spike detection.
2. **Ratio.** The SNR is a linear amplitude ratio between the spike peak and the background amplitude. The defining sentence happens to name the background first; the intent is clearly spike/background (§11).
3. **Background amplitude.** Computed on a **1 s** interval directly preceding spike onset, as the *"difference between the averaged positive and negative envelope"* (p. 1151). The envelope is the **absolute value of the Hilbert transform** of that 1 s interval (p. 1151).
4. **Channel.** One channel per source simulation: the channel with the **maximal spike amplitude in the noise-free simulation** (spike only, no background) (p. 1151). [F] This is done separately for EEG, magnetometers and gradiometers (separate maps in Figs 3–5; per-type channel labels in Fig. 6).
5. **Scale.** Linear, not dB (p. 1151).
6. **Detection threshold.** 2.5, based on the clinical experience of the magnetoencephalographer "MF" (p. 1151).
7. **Binning and statistics.** SNRs are binned by depth (N) and orientation (M). EEG–MEG differences per bin were tested with Student's *t*-test or Wilcoxon's rank-sum test, depending on the distributions (p. 1151).

### 8.2 Formalization (our transcription; the paper prints no equations)

Symbols:

- $s_c(t)$: noise-free spike response in channel $c$.
- $b_c(t)$: background response in channel $c$.
- $x_c(t) = s_c(t) + b_c(t)$: noisy signal.
- $t_\mathrm{on}$: spike onset.
- $\mathcal{T}$: sensor type (MM, GM or EEG).

$$c^{\ast} = \arg\max_{c \in \mathcal{T}} \; \max_t |s_c(t)|$$

$$W = [\,t_\mathrm{on} - 1\,\mathrm{s},\; t_\mathrm{on}) \quad (1000 \text{ samples at } 1\,\mathrm{kHz})$$

$$e(t) = \big|\, x_{c^{\ast}}(t) + j\,\mathcal{H}\{x_{c^{\ast}}\}(t) \,\big|, \quad t \in W \qquad (\text{MATLAB: } \texttt{abs(hilbert(x))})$$

$$A_\mathrm{bg} = \langle e_+ \rangle_W - \langle e_- \rangle_W, \quad e_\pm = \pm e \;\;\Rightarrow\;\; A_\mathrm{bg} = 2\,\langle e \rangle_W \quad (\text{literal reading})$$

$$\mathrm{SNR} = A_\mathrm{spike} / A_\mathrm{bg}, \qquad A_\mathrm{spike} = \text{"spike peak" in channel } c^{\ast} \;(\text{definition open, §8.3})$$

[D] For a zero-mean Gaussian background with standard deviation σ, the analytic-signal magnitude is Rayleigh distributed, so ⟨e⟩ = σ√(π/2) ≈ 1.2533σ and $A_\mathrm{bg}$ ≈ 2.51σ. Under the literal reading, SNR = 2.5 then corresponds to $A_\mathrm{spike}$ ≈ 6.3σ.

### 8.3 Elements not specified (the re-implementer must choose)

- **Numerator source.** The noise-free spike amplitude, or the peak of the noisy trace [NS].
- **Numerator measure.** max|·|, signed positive peak, or peak-to-peak of the complex [NS]. The search window for the peak is [NS].
- **"Maximal spike amplitude" for channel selection.** Absolute peak or peak-to-peak [NS]. Planar gradiometer pairs are not combined [F].
- **Denominator details.** Demeaning before the Hilbert transform; whether the transform is applied to the 1 s segment alone (edge effects) or to a longer segment and then cropped; whether ± envelopes come from ±|hilbert| or from MATLAB `envelope` (mean-centred). All [NS].
- **Filters.** No filtering of the sensor signals before SNR computation is mentioned [NS].

### 8.4 The 2.5 visual-detection threshold

- **Origin.** Set from the clinical experience of the magnetoencephalographer "MF", a co-author; the same initials label the second head model (p. 1151).
- **How it was set.** MF looked at simulated recordings without knowing where spikes occurred and judged 2.5 to be the critical level for visual detection (p. 1159).
  - [NS]: the number of traces, the modalities shown, display settings and filters, and inter-rater checks.
- **How it is used:**
  - To separate "detectable" from "non-detectable" bins (p. 1159). [F] In Fig. 3 the 2.5 boundary falls between two colour classes (cyan 2.0–2.5, aquamarine 2.5–3.0).
  - In the combined-sensor argument (p. 1158).
  - For the example sources (p. 1156).
  - In the depth-limit statements (pp. 1155, 1160).
- **Caveat.** This is an empirical clinical threshold for single-channel visual review, not a universal detector threshold.

### 8.5 Aggregation and statistics

- **Per-bin value.** The arithmetic mean of per-source SNR; the captions call the plotted values averaged SNRs (pp. 1152–1154).
- **Difference maps.** Differences of bin means.
- **Tests.** Student's *t*-test or Wilcoxon rank-sum test per bin. Rank-sum implies unpaired samples, although the same sources feed both modalities (p. 1151).
- **Significance levels.** \* p < 0.05, \*\* p < 0.01 (captions).
- **[NS]:** the normality criterion used to choose between the tests, and any multiple-comparison correction.
- **How the text's headline numbers were computed [D, Dg].** They match *unweighted* averages of bin means over row and column groups:
  - Superficial = 20–35 mm.
  - Deep = 50–60 mm or 45–60 mm for dipoles, 45–55 mm for patches.
  - Radial = 0–30° for dipoles, 10–30° for patches.
  - Tangential = 70–90°.
  - Count-weighted averages fit worse. See §9.5 and Appendix C.

### 8.6 Channel-selection check (EEG only)

- For the **64-electrode EEG** only, the authors checked whether the noise-free max-amplitude channel is also the max-SNR channel in the noisy simulation (pp. 1158–1159).
- Of **10 200** simulations (printed as "102 00", a typesetting split; [D] 9703 + 497 = 10 200) (both subjects, all dipole traces), **9703 (95.1%)** agreed and **497 (4.9%)** had the maximum in a neighbouring channel.
- Above the 2.5 threshold, the selected channel always had the highest SNR.
- No such check is reported for MEG.

### 8.7 Consistency of the literal definition with printed SNRs [Dg]

We digitized all 24 Fig. 6 traces and recomputed candidate SNRs in pixel units. Details are in Appendix E.

| Candidate SNR | Ratio to printed SNR | Traces |
|---|---|---|
| Noisy peak / (2⟨e⟩) (literal) | 0.71 ± 0.11 | 12 with printed > 2.4 |
| Noisy peak / (2⟨e⟩) (literal) | 0.70–0.80 | 6 clearest spikes (printed 3.81–6.22) |
| Spike peak-to-peak / (2⟨e⟩) | 1.04 ± 0.14 | 12 with printed > 2.4 |
| Spike peak-to-peak / (2⟨e⟩) | 0.99–1.07 | 6 clearest spikes |
| Noisy peak / ⟨e⟩ | 1.42 ± 0.22 | 12 with printed > 2.4 |

- Fig. 2(d) shows only 0.5 s of baseline and gives 4.1 (literal) and 6.0 (peak-to-peak), against a printed 4.8. This is inconclusive.
- **Conclusion.** The literal reading probably *underestimates* the published SNRs by about 25–30%. A peak-to-peak numerator, or equivalently a background measure about 1.4× smaller (e.g. peak/(2·RMS) ≈ 0.89× printed), fits better.
- The digitization cannot fix the formula. Implement it with explicit switches and report both variants (§11, item B1).

---

## 9. Results

### 9.1 Figures: content and axis ranges

- **Fig. 1** (p. 1149).
  - (a) Left-hemisphere white–gray surface coloured by orientation, 0–90° in nine discrete 10° classes (brown = radial, yellow = tangential).
  - (b) One trace, deep radial → tangential → superficial radial.
  - (c) Four patches (blue dipoles) from seed dipoles (red) of that trace.
- **Fig. 2** (p. 1150).
  - (a) Background waveforms S₁…S_N: −10 to 10 nAm over 0–0.5 s.
  - (b) Spike-wave complex: −300 to 600 nAm over −0.1 to 0.3 s.
  - (c) Head model with MEG sensors and EEG electrodes.
  - (d) EEG channel example: −50 to 100 µV over −0.5 to 0.3 s, SNR = 4.8.
- **Fig. 3** (p. 1152). Rows: dipoles (a) and patches (b). Columns: EEG, MM, EEG−MM.
  - x axis: orientation, 0–90° (patches 10–90°). y axis: depth, 20–60 mm (patches 20–55 mm).
  - SNR colour scale 0–6 in 12 discrete classes of 0.5 [F]. This scale saturates: one MM dipole bin is 6.0–6.5 in Fig. 4.
  - ΔSNR scale: ticks −4…4; red = EEG higher, blue = MM higher. Stars mark significance.
  - (c),(d): number of sources per bin.
- **Fig. 4** (p. 1153). Dipoles.
  - (a) EEG, MM, GM on an SNR scale of 0–10 in 20 classes of 0.5 [F].
  - (b) EEG−GM and MM−GM; ΔSNR ticks −6…6.
- **Fig. 5** (p. 1154). The same as Fig. 4 for patches (10–90°, 20–55 mm).
- **Fig. 6** (p. 1157). Example sensor traces over 0–2 s. Scale bars as printed: 100 uV (EEG), 5 pT (MM), 100 pT (GM; the unit should be T/m, §11).

### 9.2 Main qualitative findings (MM vs GM vs EEG)

- EEG has its highest SNR for superficial radial sources. MEG (MM and GM) has its highest SNR for superficial tangential sources (p. 1154).
- SNR falls with depth in all modalities, and faster in MEG (p. 1154).
- MM beats EEG only up to 40 mm depth and only for predominantly tangential orientations (pp. 1154–1155).
- GM beats EEG over a larger region than MM does (p. 1155).
- **GM > MM** in general, except for deep sources (below 55 mm), where there is no difference (p. 1155).
  - MM and GM are similar at 50–60 mm and not significantly different at 55–60 mm (pp. 1156–1158).
  - The abstract describes the two as comparable in magnitude for deep sources (p. 1147).
- The authors attribute GM's advantage for superficial sources to its more local sensitivity, which suppresses distant background. Magnetometers pick up more distant background in this simulation scheme (p. 1158). Tarkiainen et al (2003) is cited for similar MM/GM effects.
- Dipole and patch results are similar (pp. 1155, 1158).
- **Detectability limits:**
  - Bin means fall below 2.5 for sources deeper than **40 mm (MM)** and **50 mm (GM and EEG)** (p. 1155).
  - Elsewhere the paper says no sufficiently high SNR was reached deeper than **60 mm (EEG)** and **55 mm (MEG)** (p. 1160). See §11.
- **Maximum over all EEG + MEG sensors (dipoles).** Every orientation at 20–45 mm exceeds 2.5. At 50–60 mm the SNR no longer depends on orientation (p. 1158).

### 9.3 Dipole numbers as printed (p. 1155)

| Sensor | Superficial radial | Superficial tangential | Deep radial | Deep tangential |
|---|---|---|---|---|
| EEG | ≈ 5 (0–30°) | 3.2 (−36% vs radial) | 2.5 (−50% with depth; depth given as 50–60 mm) | 2.3 |
| MM | 2 (−57.5% from 4.7) | ≈ 5 in words; 4.7 in numbers (70–90°) | [D] ≈ 1.5 (−68.1% from superficial tangential) | 2 (−57.5% with depth) |
| GM | 3 (−61.6% from 7.8) | ≈ 8 in words; 7.8 in numbers | [D] ≈ 1.6 (−79.5% from superficial tangential) | 2.3 (−70.5% with depth) |

- The MM response to deep sources (45–60 mm) still depends on orientation (p. 1155).

### 9.4 Patch numbers as printed (p. 1155)

| Sensor | Superficial radial | Superficial tangential | Deep radial | Deep tangential |
|---|---|---|---|---|
| EEG | 5 (10–30°) | 2.7 | 2.6 | 2.3 |
| MM | 1.6 | 4.3 (70–90°) | 1.5 | 2 |
| GM | not given in the text; Fig. 5 only | | | |

### 9.5 Significance patterns

**MM > EEG** (p. 1155):

- Dipoles: 50–60° and 70–90° at 20–25 mm; 60–90° at 25–35 mm; 70–90° at 35–40 mm. [F] That is 11 bins.
- Patches: the same bins except 50–60° at 20–25 mm, i.e. 10 bins (p. 1156).

**GM > EEG** (dipoles, as printed, p. 1155): 30–90° at 20–25 mm; 40–90° at 20–35 mm; 50–90° at 35–40 mm; 70–90° at 40–50 mm.

**MM vs GM, dipoles (Fig. 4(b), right)** [F]:

- Every bin at 20–35 mm is significant (\*\*), with GM higher.
- At 35–50 mm the bins are significant except the most radial ones (0–10° at 35–40 mm; 0–20° at 40–45 mm; 0–30° at 45–50 mm).
- At 50–55 mm only the 40–50° bin is significant (\*). No bin at 55–60 mm is significant.

**MM vs GM, patches (Fig. 5(b), right)** [F]:

- The tangential bins (70–90°) are significant (GM higher) at 20–50 mm, and 70–80° also at 50–55 mm.
- Most bins at 25–35 mm are significant.
- The sparse superficial bins at intermediate orientations (20–25 mm, 30–70°) are not significant.

**Digitized check of the headline values** [Dg; Appendix C]: unweighted means of bin means reproduce the printed values.

| Quantity | Digitized | Printed |
|---|---|---|
| GM superficial tangential | 7.92 | 7.8 |
| GM superficial radial | 2.97 | 3 |
| GM deep tangential (50–60 mm) | 2.25 | 2.3 |
| MM superficial tangential | 4.83 | 4.7 |
| MM superficial radial | 2.03 | 2 |
| MM deep tangential | 2.00 | 2 |
| MM deep radial (45–60 mm) | 1.47 | ≈ 1.5 |
| EEG superficial radial | 5.03 | ≈ 5 |
| Patch MM superficial tangential | 4.25 | 4.3 |
| Patch MM superficial radial | 1.67 | 1.6 |
| Patch MM deep radial | 1.50 | 1.5 |
| Patch EEG superficial radial | 5.00 | 5 |

Count-weighted means fit worse (GM 7.14; MM 4.33).

### 9.6 Example sources (Table 1 p. 1156; Fig. 6 p. 1157)

The full table is in Appendix D.

- **Setup.** Eight examples: dipole or patch × superficial or deep × radial or tangential, all in the posterior frontal lobe (p. 1156).
- **Highest SNRs.** GM for the superficial tangential patch is highest at 6.22, followed by GM for the superficial tangential dipole at 5.82 (p. 1156).
- **Best in the other sensor types.** For MM, the superficial tangential dipole reaches 4.7 (patch: 4.24). For EEG, the superficial radial dipole reaches 4.29 (patch: 3.81) (p. 1156).
- **Effect of depth, tangential sources (p. 1156):**

| Sensor | Dipole | Patch |
|---|---|---|
| MM | −1.81 (38.5%) | −1.28 (30.2%) |
| GM | −2.89 (49.6%) | −3.69 (59.3%) |

- **Effect of depth, EEG radial sources:** −1.86 (43.4%) for the dipole and −0.85 (22.3%) for the patch (p. 1156).
- **Threshold claim.** The authors say that for four of the eight examples the SNR fell below 2.5: the MEG radial and EEG tangential cases, independent of depth (p. 1156). See §11: by Fig. 6, 13 of the 24 source–modality traces are below 2.5, including the EEG deep radial dipole at 2.43.

### 9.7 Additional simulations (pp. 1158–1159)

- **64-electrode EEG.** Superficial radial SNR 4.5, which is 10% below the 128-electrode value. Deep sources change only marginally (supplement Fig. 7).
- **300% background.** SNRs fall at all depths and orientations. The EEG–MM relation is unchanged, and the orientation difference still shrinks with depth. No numbers are given.
- **Skull conductivity ratios 1:50 and 1:25.** The authors report similar results. The authors reason that SNR is fairly insensitive to skull conductivity because spike and background pass through the same model. No numbers are given.

---

## 10. Limitations and discussion of SNR definitions (pp. 1158–1160)

1. **Best-single-channel SNR.** Sensitivity analyses are often done on whole sensor arrays (e.g. Goldenholz et al 2009); this study uses the best single sensor (p. 1158).
   - Defining SNR across all channels, as Goldenholz et al did, would give different values.
   - The numbers are therefore **not directly comparable** with Goldenholz et al (2009) or Ahlfors et al (2010a) (p. 1158).
2. **Clinical justification.** Single-channel SNR follows standard clinical evaluation and was used before by de Jongh et al (2005). The authors say this makes the results directly usable by clinicians (p. 1159).
3. **Channel-selection check** (64-channel EEG): 95.1% of cases agreed, and all cases above threshold agreed (pp. 1158–1159).
4. **Threshold 2.5.** Set by MF through blinded visual inspection of simulated recordings (p. 1159).
   - The colours below 2.5 in Fig. 3 show that many deep configurations fail in both EEG and MEG. This agrees with intracranial/scalp studies such as Wennberg et al (2011) (p. 1159).
5. **Transfer to evoked responses.** The signal shape and noise model are said to transfer to the comparison of evoked potentials and fields (p. 1159).
6. **No sensor noise** (p. 1158).
   - Biological noise is argued to dominate (Hämäläinen et al 1993; Vrba 2000; Parkkonen 2010).
   - Sensor noise differs by sensor type (Ahonen et al 1993; Vrba 2000).
   - New sensor types such as atomic magnetometers can be assessed only once their sensor noise is known.
7. **Geometric approximations** (p. 1159).
   - Node-to-node depth may overestimate depth; the error is judged negligible (scalp triangles about 6 mm).
   - Angles use node normals; the error is bounded by 4.7 ± 2.5°.
8. **Conductivity.** Signal and noise share the forward model, so SNR is fairly insensitive to skull conductivity. Ratios 1:50 and 1:25 gave analogous results (p. 1159).
9. **Patch size** (pp. 1159–1160). Patches were kept small on purpose to limit cancellation.
   - The literature uses circular patches 1–50 mm or 20–80 mm in diameter, and 4 cm² or 6 cm² patches.
   - Cancellation matters for large patches.
   - Data-driven patch designs exist.
10. **Dipoles vs patches** (p. 1158). Patch SNRs are lower because the net strength is reduced by orientation spread. The dipole–patch differences at 20–25 mm and 30–70° were likely due to n < 10.
11. **Background model** (p. 1160).
    - Strengths: the novel asymmetry-ratio-based model plus random dipole distribution is closer to reality than white noise, and gives EEG and MEG unbiased noise shares.
    - Limitations: no regional weighting of bands, which would make detectability less homogeneous across regions; independent realizations per source; no resting-state network structure.
12. **Why GM beats MM for superficial sources.** GM sensitivity is more local and suppresses distant background (p. 1158).
13. **Conclusion.** EEG and MEG are complementary. Deep spikes are more likely seen in EEG. Simultaneous EEG and MEG improves detection (pp. 1147, 1160).

---

## 11. Ambiguities, inconsistencies, and decisions for the re-implementer

### A. Internal inconsistencies and text/figure mismatches

1. **Dipole total.** The abstract says 5600 dipoles (p. 1146). The methods give 127 + 128 traces × 20 dipoles = 5100 (p. 1149). The 10 200 channel-check simulations (printed "102 00", p. 1158) equal 2 × 5100. That they are dipoles plus seed patches is our conjecture.
2. **Per-bin sample-size statistics** (p. 1153) do not match Fig. 3(c),(d).
   - Dipoles: the text gives 13–117 (48.6 ± 30.6); the figure gives 15–177 (52.5 ± 30.0).
   - Patches: the range matches; the mean ± SD is 37.2 ± 34.8 in the text vs 40.4 ± 35.7 from the figure.
   - The figure counts sum exactly to the printed 3783 and 2264, so the figure is the more reliable source.
3. **Inconsistent use of "superficial" and "deep."**
   - Definitions: superficial 20–35 mm, deep 45–60 mm (p. 1154).
   - "Deep" is given as 50–60 mm on p. 1155. "Superficial" is given as 20–25 mm on p. 1156.
   - The Table 1 "deep" examples lie at 39.3–44.7 mm, i.e. *intermediate* (p. 1156).
4. **Depth-limit statements disagree.** Below threshold beyond 40 mm (MM) and 50 mm (GM, EEG) on p. 1155, vs beyond 55 mm (MEG) and 60 mm (EEG) on p. 1160.
   - [Dg] Digitized maps: MM bins reach ≥ 2.5 only down to 35–40 mm. GM reaches 2.5–3.0 in tangential bins down to 50–55 mm. EEG reaches 2.5–3.0 in one 55–60 mm bin.
5. **Discussion vs results on EEG vs MEG at depth.** The Discussion says EEG was slightly higher than MEG at 40–60 mm (p. 1156). The abstract says deeper sources give higher EEG SNR for all orientations (p. 1147). The results report GM significantly above EEG at 70–90° and 40–50 mm (p. 1155), and [F] Fig. 4(b) shows this in the "deep" 45–50 mm row. These summary statements fit MM, not GM.
6. **Wrong panel references.** The EEG−GM difference maps are cited as Figs 4(a)/5(a) (p. 1155) but are in panels (b). The Fig. 3 caption uses A–D while the figure uses (a)–(d).
7. **GM > EEG range (minor).** The text lists "40–90°, 20–35 mm" (p. 1155), which overlaps the "30–90°, 20–25 mm" item. It is redundant, not contradictory: [F] Fig. 4(b) shows significance from 30° at 20–25 mm and from 40° at 25–35 mm.
8. **SNR sentence wording.** The defining sentence literally reads as background-to-spike, i.e. inverted (p. 1151). All numbers are spike/background.
9. **Threshold claim for the examples.** The paper says four of the eight examples fall below 2.5 and names the MEG radial and EEG tangential cases (p. 1156). Its own parenthetical does not add up to four. [F] Fig. 6 (p. 1157) shows 13 of 24 traces below 2.5: 4 MM radial, 4 GM radial, 4 EEG tangential, and the EEG deep radial dipole at 2.43, which the parenthetical does not cover.
10. **Fig. 6 trace vs label** [Dg]. The EEG FC3 "patch, deep, radial" trace is labelled SNR 2.96 but shows a visibly *smaller* spike than the "dipole, deep, radial" trace (SNR 2.43), with the same channel and background. Digitized spike peaks are ≈ 8 vs ≈ 13 px against the same ≈ 3.2 px baseline SD. The printed values agree with the text (p. 1156), so the trace or label may be swapped.
11. **Fig. 6 caption vs methods.** The caption says the sensors with the **highest SNR** are shown (p. 1157). The methods select the channel with the **largest noise-free amplitude** (p. 1151). For 64-channel EEG the two coincided in 95.1% of cases (§8.6); this was not checked for MEG.
12. **Fig. 6 gradiometer unit.** The scale bar reads "100 pT"; the unit should be T/m (e.g. pT/m). The EEG bar reads "100 uV" (p. 1157).
13. **Fig. 2(c) markers.** The caption says MEG sensors are blue squares; the drawing shows circles or ellipses (p. 1150).
14. **Discussion colour description.** It equates the blue-to-green colours of Fig. 3 with SNR < 2.5 (p. 1159). [F] On the Fig. 3 scale, the green-ish class is 2.5–3.0, and the colours below 2.5 run from blue to cyan.
15. **Patch density.** The Discussion cites a preset current density of 100 nA mm⁻² (p. 1159). Grammatically this belongs to earlier triangle-mesh approaches, but the paper calls its own design consistent with them. The printed strengths imply about 30 nAm/mm² (p. 1151), and the units differ. Treat 612–678 nAm (median 622) as authoritative.
16. **Fig. 3 saturation.** The Fig. 3 SNR scale tops out at 6. The MM dipole bin at 20–25 mm, 70–80° is 6.0–6.5 on the Fig. 4 scale. Use Fig. 4/5 for values above 5.5.
17. **Duplicate label "MF".** The same initials label both the second head model and the clinician who set the threshold (pp. 1148, 1151, 1159).

### B. Unstated parameters requiring an explicit decision

1. **SNR numerator** [NS]. Choose noise-free vs noisy, and |peak| vs peak-to-peak.
   - The literal reading is noisy or noise-free |peak| / (2⟨e⟩). Digitization suggests this is about 25–30% below the printed values, while peak-to-peak / (2⟨e⟩) matches within ±7% for clear spikes (§8.7, Appendix E).
   - Recommendation: keep the literal reading as "faithful", add the peak-to-peak variant as a sensitivity analysis, and report which one is used for any comparison with the 2.5 threshold.
2. **Hilbert details** [NS]: demeaning, edge handling (1 s segment alone vs a longer segment then cropped), upper/lower envelope construction, and pre-filtering (none mentioned).
3. **Channel-selection metric** [NS]: absolute peak vs peak-to-peak, and the time window. Selection is per sensor type; gradiometers are single channels, not combined pairs [F].
4. **Spike waveform equation** [NS]. Rebuild from the digitized key points (§5; Appendix A) as a sum or sequence of sinusoidal segments with a 600 nAm peak. Record the construction and check the ≈74–80 ms sharp component and the ≈200 ms total.
5. **Spike onset within the 6 s epoch** [NS]. It must be ≥ 1 s. Fig. 6 suggests about 0.9–1 s of pre-onset data in the displayed window.
6. **Background spectrum** [NS]: band edges, filter design, band weights (only a range of 0.4–0.6 is given), and the Gotman et al (1973) asymmetry-ratio definition (by citation only).
   - Record the choices. Consider checking the simulated sensor spectrum against typical resting EEG/MEG.
7. **Background normalization** [NS]: peak (max|s| = 10 nAm) vs min–max; per band vs after summation; no area weighting.
   - MNE note: the MEG→OPM comparison needs a mesh-independent convention; GOAL.md asks for area/covariance scaling. Keep the Hunold convention in its own configuration.
8. **Background realizations** [NS]. Fig. 6 indicates one fixed realization reused across sources. Use a fixed seed, reuse the realization across source positions *and* sensor arrays, and optionally add repeated realizations as an extension.
9. **Background placement** [NS]: the random 10% node subset (seed), both hemispheres, and whether spike-source nodes are excluded.
10. **Orientation angle** [NS]: folding (arccos|n·m|) vs excluding angles > 90°; node-normal computation. Watch the medial surfaces, where the nearest inner-skull node may be far away or oblique.
11. **Depth and orientation association** [NS]: the paper uses nearest *node* on 5120-triangle meshes. Using nearest surface point instead is a documented deviation.
12. **Bin edges and out-of-range sources** [NS]: which side is closed; sources shallower than 20 mm.
13. **Trace construction** [NS]: how dipoles 10–19 and 20 are chosen, and the "closest" metric. With MNE and new anatomy, exact trace replication is impossible. Equivalent designs are to sample sources so the bins are filled evenly, or to use all white-surface vertices with per-bin subsampling.
14. **Patch growth** [NS]: growth order, area computation (vertex vs triangle areas) and moment allocation (density × vertex area assumed).
    - **Strength convention:** Hunold uses the scalar sum (612–678 nAm). GOAL.md requires comparing fixed total scalar moment and fixed moment density separately, preserving signed summation.
15. **MEG geometry** [NS]: helmet position relative to the head and Vectorview coil type.
    - Use actual Vectorview geometry and transforms (e.g. MNE sample data). Pick T1 or T3 (3022/3012 vs 3024/3014) and record the choice. Use `accuracy='normal'` for the 4-point integration.
16. **BEM numerics.** Galerkin (paper) vs MNE's linear collocation. Expect small differences.
17. **Statistics** [NS]: the rule for choosing the *t*-test vs the rank-sum test; no multiple-comparison correction reported; pooling of subjects is implied.
18. **300% background** [NS]: amplitude vs power.

---

## 12. Reproduction configuration

| Parameter | Value as printed | Page | Notes |
|---|---|---|---|
| Subjects | 2 healthy males, AH 22 y, MF 49 y | 1148 | Substitute MNE sample/other anatomy; document |
| MRI | 3T Siemens MAGNETOM Trio, A Tim System; 3D MPRAGE; two 3D multi-echo FLASH (30°, 5°) | 1148 | Sequence parameters [NS] |
| Cortical mesh | FreeSurfer white–gray boundary, ≈300 000 elements | 1148 | [MNE] white surface, `spacing='all'` or dense; [D] ≈0.654 mm² per node |
| Inner skull | Watershed on MPRAGE | 1148 | `mne watershed_bem` |
| Outer skull, scalp | From FLASH | 1148 | `mne flash_bem` if FLASH available, else watershed (deviation) |
| BEM mesh | 5120 triangles per surface | 1148 | `ico=4` |
| Conductivity brain/skull/scalp | 0.33 / 0.0042 / 0.33 S m⁻¹ | 1151 | Ratio ≈ 1:78.6; 1:50 and 1:25 gave similar results (p. 1159) |
| BEM formulation | Galerkin, linear basis, isolated source approach; Helsinki BEM library (Matlab) | 1151 | MNE uses linear collocation + isolated skull |
| MEG system | Elekta Neuromag Vectorview: 102 magnetometer and 204 gradiometer positions | 1151 | Coil type [NS]; T1 3022/3012 or T3 3024/3014 |
| Coil integration | 4-point numerical integration per sensor | 1151 | [MNE] `accuracy='normal'` (4 points for both types) |
| Sensor–head registration | [NS] | – | Use real device-to-head transform |
| Sensor noise | Omitted | 1158 | Brain-noise-only reference |
| EEG (out of scope) | 128-ch WaveGuard, 10–5 incl. subtemporal chain; 64-ch 10–10 variant | 1151, 1158 | Reference [NS] |
| Source locations | Nodes of white–gray boundary mesh | 1148 | |
| Source orientation | Perpendicular to cortex (fixed) | 1148 | Polarity [NS] |
| Depth definition | Euclidean distance to closest scalp-mesh node | 1148 | Node to node; skin surface |
| Orientation definition | Angle between normals of cortical node and closest inner-skull-mesh node; 0° radial, 90° tangential | 1148 | Folding above 90° [NS] |
| Depth bins | 8 even intervals, 20–60 mm (5 mm) | 1150 | Patches shown 20–55 mm |
| Orientation bins | 9 × 10°, 0–90° | 1150 | Patches shown 10–90° |
| Depth labels | Superficial 20–35, intermediate 35–45, deep 45–60 mm | 1154 | Used inconsistently (§11 A3) |
| Spike-source sampling | 127 (AH) + 128 (MF) traces × 20 dipoles, sulcal bottom to gyral crown | 1148–1149 | [D] 5100 vs abstract 5600 |
| Trace seed | Radial (< 6°) sulcal-bottom dipole, all neighbours more tangential and superficial | 1149 | |
| Trace dipoles 2–9 | Direct neighbours; orientation decades 10–20 … 80–90° | 1149 | Dipole 1 in 0–10° |
| Trace dipole 20 | Closest most-radial dipole on gyral crown | 1149 | Metric [NS]; dipoles 10–19 chosen to best match the decades [NS] |
| Dipole strength | 600 nAm | 1151 | = waveform peak (Fig. 2b) |
| Patch seeds | Each trace dipole | 1149 | |
| Patch orientation window | ±10° of seed; 0–20° if seed < 10°; 70–90° if seed > 80° | 1149 | Scalar angle |
| Patch area rule | Stop when area first exceeds 20 mm²; smaller patches excluded | 1149–1150 | Area, not radius |
| Number of patches | 3300 (≥ 20 mm²) | 1146, 1150 | 2264 analysed |
| Patch strength | 612–678 nAm, median 622 nAm; uniform activation; fixed moment density | 1151 | [D] ≈30 nAm/mm²; scalar sum |
| Patch depth/orientation | Mean over patch dipoles | 1150 | Table 1 gives mean ± SD |
| Spike waveform | Several sinusoids; spike ≈ 80 ms; Fig. 2(b) | 1150 | Digitized: peak 600 at 24 ms; trough −254 at 54 ms; slow wave to 0.2 s |
| Sampling rate | 1000 sps | 1150 | |
| Simulated duration | 6 s | 1150 | Onset time [NS] |
| Background sources | Random 10% of nodes, ≈30 000 (AH 29 823) | 1148, 1150 | Seed [NS] |
| Area per background dipole | 6.54 mm² (AH), 6.44 mm² (MF) | 1150–1151 | |
| Background time course | Stationary Gaussian noise in 5 EEG bands; weights 0.4–0.6 via Gotman asymmetry ratio; bands summed | 1150 | Band edges, filters, weights [NS] |
| Background amplitude | ±10 nAm (max 10 nAm per dipole) | 1150–1151 | Peak normalization assumed |
| Background correlation | Independent realization per dipole | 1150, 1160 | |
| Background reuse across sources | [NS] | – | [Dg] Fig. 6: same realization reused |
| Increased background | 300% | 1158 | Amplitude vs power [NS] |
| SNR channel | Maximal spike amplitude in noise-free simulation | 1151 | Per sensor type; single gradiometer channels |
| SNR baseline | 1 s directly before spike onset | 1151 | |
| Background amplitude | Mean positive envelope − mean negative envelope; envelope = abs(Hilbert) | 1151 | Literal: 2 × mean(abs(analytic signal)) |
| Spike amplitude | "Spike peak" | 1151 | Noisy vs noise-free and peak vs peak-to-peak [NS]; see §8.7 |
| SNR scale | Linear amplitude ratio | 1151 | Not dB |
| Detection threshold | 2.5 (visual; clinician MF) | 1151, 1159 | Not a universal detector threshold |
| Bin statistic | Mean SNR per bin | 1152–1154 | Unweighted means reproduce text summaries |
| Tests | Student's *t* or Wilcoxon rank sum; \* p < 0.05, \*\* p < 0.01 | 1151–1154 | No correction reported |
| Analysed sources | 3783 dipoles; 2264 patches | 1151 | Match Fig. 3(c),(d) sums |
| Channel-check (EEG 64) | "102 00" (= 10 200) sims: 9703 (95.1%) same channel, 497 (4.9%) neighbouring | 1158–1159 | Not done for MEG |

---

## Appendix A. Digitized spike waveform

The key points are in the table in §5. Source panel: Fig. 2(b), p. 1150, at 200 ppi. Axis calibration: 447.5 px/s and 4.245 nAm/px, from the printed ticks. The sampled centre line, approximately every 2.2 ms, is:

| t (ms) | nAm | | t (ms) | nAm | | t (ms) | nAm |
|---|---|---|---|---|---|---|---|
| 0 | 0 | | 39 | ≈ −10 | | 110 | ≈ 80 |
| 4 | −30 | | 44 | −132 | | 120 | ≈ 20 |
| 8 | −66 | | 49 | −215 | | 130 | ≈ −40 |
| 11 | −67 | | 54 | −253 | | 140 | ≈ −79 |
| 13 | ≈ 0 | | 60 | −230 | | 150 | ≈ −91 |
| 15 | 130 | | 66 | −147 | | 160 | ≈ −81 |
| 17 | 310 | | 71 | −66 | | 170 | ≈ −60 |
| 20 | 466 | | 75 | ≈ 24 | | 180 | ≈ −28 |
| 22 | 560 | | 80 | 111 | | 190 | ≈ −8 |
| 24 | 591–600 | | 87 | 149 | | ≥198 | 0 |
| 29 | 525 | | 93 | 138 | | | |
| 33 | 338 | | 100 | 121 | | | |

On the steep flanks, the per-column centroid underestimates extremes by a few nAm.

## Appendix B. Sources per bin (Fig. 3(c),(d), p. 1152) [F]

Rows are depth in mm; columns are orientation in degrees.

**Dipoles** (Fig. 3(c)). Sum 3783, which matches p. 1151.

| depth \ orient | 0–10 | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 | Σ |
|---|---|---|---|---|---|---|---|---|---|---|
| 20–25 | 45 | 33 | 29 | 30 | 24 | 19 | 17 | 16 | 30 | 243 |
| 25–30 | 47 | 43 | 52 | 63 | 48 | 56 | 82 | 91 | 130 | 612 |
| 30–35 | 67 | 49 | 56 | 51 | 63 | 64 | 95 | 111 | 177 | 733 |
| 35–40 | 56 | 55 | 59 | 66 | 62 | 75 | 72 | 106 | 143 | 694 |
| 40–45 | 65 | 66 | 66 | 63 | 64 | 78 | 64 | 57 | 72 | 595 |
| 45–50 | 32 | 35 | 40 | 41 | 47 | 51 | 51 | 57 | 62 | 416 |
| 50–55 | 42 | 37 | 43 | 43 | 36 | 23 | 33 | 21 | 17 | 295 |
| 55–60 | 25 | 24 | 23 | 17 | 21 | 17 | 21 | 15 | 32 | 195 |
| Σ | 379 | 342 | 368 | 374 | 365 | 383 | 435 | 474 | 663 | **3783** |

**Patches** (Fig. 3(d)). Sum 2264, which matches p. 1151.

| depth \ orient | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 | Σ |
|---|---|---|---|---|---|---|---|---|---|
| 20–25 | 34 | 13 | 3 | 4 | 2 | 2 | 10 | 12 | 80 |
| 25–30 | 15 | 13 | 15 | 25 | 31 | 50 | 122 | 128 | 399 |
| 30–35 | 41 | 20 | 24 | 30 | 45 | 65 | 149 | 163 | 537 |
| 35–40 | 22 | 29 | 28 | 28 | 33 | 47 | 103 | 129 | 419 |
| 40–45 | 59 | 36 | 32 | 37 | 42 | 41 | 45 | 46 | 338 |
| 45–50 | 29 | 27 | 21 | 32 | 34 | 37 | 69 | 49 | 298 |
| 50–55 | 48 | 25 | 38 | 16 | 10 | 21 | 25 | 10 | 193 |
| Σ | 248 | 163 | 161 | 172 | 197 | 263 | 523 | 537 | **2264** |

## Appendix C. Digitized bin-mean SNR, MM and GM [Dg]

- **Sources.** Fig. 4(a) for dipoles (p. 1153) and Fig. 5(a) for patches (p. 1154). Both use a 0–10 scale with 20 discrete colour classes 0.5 wide.
- **Method.** Each cell's colour was matched to its class. An entry *x* means the bin mean lies in [*x*, *x* + 0.5).
- **Cross-check.** The EEG and MM classes are identical to those read from Fig. 3 (0–6 scale), except the saturated MM 20–25 mm / 70–80° bin.
- **Saturation.** Top-class entries (9.5) may exceed 10.
- **EEG.** Digitized too but omitted as out of scope; available on request.

**Dipoles, MM** (rows 20–25 … 55–60 mm; columns 0–10 … 80–90°)

| depth | 0–10 | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 |
|---|---|---|---|---|---|---|---|---|---|
| 20–25 | 1.0 | 2.0 | 2.5 | 3.5 | 4.0 | 5.0 | 4.5 | 6.0 | 5.5 |
| 25–30 | 1.0 | 2.0 | 2.5 | 2.5 | 3.5 | 4.0 | 4.5 | 4.5 | 4.5 |
| 30–35 | 1.5 | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 | 3.5 | 3.5 | 3.5 |
| 35–40 | 1.5 | 1.5 | 1.5 | 2.0 | 2.5 | 2.5 | 3.0 | 2.5 | 3.0 |
| 40–45 | 1.0 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.0 | 2.0 |
| 45–50 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 |
| 50–55 | 1.0 | 1.0 | 1.5 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 |
| 55–60 | 1.0 | 1.0 | 1.0 | 1.0 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 |

**Dipoles, GM**

| depth | 0–10 | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 |
|---|---|---|---|---|---|---|---|---|---|
| 20–25 | 1.5 | 3.0 | 4.5 | 5.5 | 7.5 | 9.5 | 8.5 | 9.5 | 9.5 |
| 25–30 | 1.5 | 3.0 | 4.0 | 5.0 | 6.0 | 7.0 | 7.0 | 7.5 | 7.5 |
| 30–35 | 1.5 | 2.5 | 3.0 | 4.0 | 5.0 | 5.0 | 6.0 | 6.0 | 6.0 |
| 35–40 | 1.5 | 2.0 | 2.5 | 3.0 | 3.0 | 3.5 | 4.0 | 4.0 | 4.5 |
| 40–45 | 1.5 | 1.5 | 2.0 | 2.0 | 2.5 | 2.5 | 3.0 | 3.0 | 3.5 |
| 45–50 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.0 | 2.5 | 2.5 |
| 50–55 | 1.5 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.5 | 2.5 |
| 55–60 | 1.0 | 1.5 | 1.0 | 1.0 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 |

**Patches, MM** (rows 20–25 … 50–55 mm; columns 10–20 … 80–90°)

| depth | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 |
|---|---|---|---|---|---|---|---|---|
| 20–25 | 1.0 | 1.5 | 3.0 | 3.5 | 3.5 | 3.0 | 4.5 | 4.5 |
| 25–30 | 1.5 | 1.5 | 2.5 | 2.5 | 3.5 | 4.0 | 4.0 | 4.0 |
| 30–35 | 1.5 | 1.5 | 2.0 | 2.5 | 3.0 | 3.0 | 3.5 | 3.5 |
| 35–40 | 1.5 | 1.5 | 1.5 | 2.0 | 2.5 | 2.5 | 2.5 | 2.5 |
| 40–45 | 1.0 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.0 | 2.0 |
| 45–50 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 | 1.5 |
| 50–55 | 1.0 | 1.0 | 1.0 | 1.5 | 1.5 | 2.0 | 1.5 | 2.0 |

**Patches, GM**

| depth | 10–20 | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 |
|---|---|---|---|---|---|---|---|---|
| 20–25 | 1.5 | 2.0 | 5.5 | 5.5 | 7.0 | 5.0 | 8.0 | 8.0 |
| 25–30 | 2.5 | 2.5 | 4.0 | 4.5 | 6.0 | 6.0 | 7.0 | 7.0 |
| 30–35 | 2.0 | 2.0 | 3.5 | 4.0 | 4.5 | 5.0 | 5.5 | 5.5 |
| 35–40 | 2.0 | 2.0 | 2.5 | 3.0 | 2.5 | 3.0 | 3.5 | 4.0 |
| 40–45 | 1.5 | 1.5 | 2.0 | 2.5 | 2.5 | 2.5 | 3.0 | 3.0 |
| 45–50 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.0 | 2.5 | 2.5 |
| 50–55 | 1.0 | 1.5 | 1.0 | 2.0 | 1.5 | 2.5 | 2.5 | 2.5 |

The patch bins at 20–25 mm and 30–70° hold only 2–4 patches each (Appendix B).

## Appendix D. Example sources: Table 1 (p. 1156) and Fig. 6 (p. 1157)

| Source | Dipole depth (mm) | Dipole orientation (°) | Patch: # dipoles | Patch area (mm²) | Patch depth, mean ± SD (mm) | Patch orientation, mean ± SD (°) |
|---|---|---|---|---|---|---|
| Superficial radial | 23.1 | 6.4 | 32 | 20.2 | 23.7 ± 0.4 | 10.7 ± 4.8 |
| Superficial tangential | 25.6 | 74.8 | 35 | 20.3 | 26.7 ± 1.6 | 73.7 ± 6.2 |
| "Deep" radial | 44.7 | 4.2 | 19 | 20.1 | 44.4 ± 0.3 | 14.2 ± 4.6 |
| "Deep" tangential | 39.3 | 74.9 | 32 | 20.5 | 39.4 ± 1.7 | 78.1 ± 5.0 |

Selected channel and printed SNR in Fig. 6. EEG: electrode; MM: Vectorview magnetometer; GM: Vectorview gradiometer.

| Source | EEG | MM | GM |
|---|---|---|---|
| Dipole, superficial, radial | FC3 4.29 | 0631 1.36 | 0413 1.35 |
| Dipole, superficial, tangential | CCP5h 0.96 | 0631 4.7 | 0412 5.82 |
| Dipole, deep, radial | FC3 2.43 | 0711 0.92 | 0423 1.73 |
| Dipole, deep, tangential | CP5 0.58 | 0711 2.89 | 0412 2.93 |
| Patch, superficial, radial | FC3 3.81 | 0711 0.84 | 0423 1.72 |
| Patch, superficial, tangential | CCP5h 0.95 | 0711 4.24 | 0412 6.22 |
| Patch, deep, radial | FC3 2.96 | 0741 1.04 | 0423 1.74 |
| Patch, deep, tangential | C5 1.05 | 0711 2.96 | 0412 2.53 |

## Appendix E. SNR-definition check against Fig. 6 [Dg]

**Method**

- Each of the 24 Fig. 6 traces (p. 1157) was digitized at the native 200 ppi: ≈ 7.65 ms per pixel column, anti-aliased 1-px line.
- The time axis was calibrated from the 0/0.5/1/1.5/2 s ticks. The spike peak lies at ≈ 0.93 s, so ≈ 0.9 s of pre-onset baseline is visible (the paper uses 1 s).
- Baseline σ was estimated from the per-column centroid plus the within-column spread. ⟨e⟩ was taken as σ√(π/2).
- Spike amplitude came from the drawn-line extremes in 0.90–1.15 s.

**Results**

| Candidate | Ratio to printed SNR (12 traces, printed > 2.4) | Six clearest spikes (printed 3.81–6.22) |
|---|---|---|
| Noisy peak / (2⟨e⟩) (literal) | 0.71 ± 0.11 (0.42–0.82) | 0.70–0.80 |
| Spike peak-to-peak / (2⟨e⟩) | 1.04 ± 0.14 (0.66–1.31) | 0.99–1.07 |
| Noisy peak / ⟨e⟩ | 1.42 ± 0.22 | – |

- **Baseline sharing.** Same-channel baselines are practically identical across sources (r ≥ 0.97). This supports a fixed background realization (§6).
- **Low-SNR traces.** Printed values of 0.58–1.05 are compatible with either numerator: with no visible spike, the maximum of the background in the spike window alone gives values of about 1.
- **Fig. 2(d)** (0.5 s baseline visible, ≈ 2.4 ms/px): literal 4.1, peak-to-peak 6.0, printed 4.8. Not decisive.

**Caveats.** Raster resolution, anti-aliasing and the 0.9 vs 1 s baseline all limit this check. It is evidence, not proof. The formula must be treated as partly unresolved (§11 B1).
