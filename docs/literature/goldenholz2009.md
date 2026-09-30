# Goldenholz et al. (2009): methods and parameter extraction for a reproduction study

*Extraction date: 2026-09-30. Status: literature extraction only (no simulations run).*

**Source PDF:** `Human Brain Mapping - 2008 - Goldenholz - Mapping the signal‐to‐noise‐ratios of cortical sources in magnetoencephalography.pdf` (project root).

**Coverage:** all 10 pages (published pp. 1077–1086) were read in full: abstract, all five figure captions and colour bars, Equations 1–3 and the reference list. The paper has **no tables, no appendices and no supplementary material**. The conductivity sentence, Equations 1–3 and every colour-bar label were re-checked on page renders at 300–800 dpi.

**Conventions**
- Page numbers are published page numbers (PDF page 1 = p. 1077 … PDF page 10 = p. 1086).
- *not stated* = the paper does not give it.
- **[inference]** = my interpretation; **[derived]** = arithmetic on printed numbers; **[external]** = from outside the paper (sources listed at the end).
- Notation follows the paper, which uses Latin *s*, not σ: s_k² = noise variance of sensor k; s_s² = variance of each noise source (the paper reports √s_s² in nAm).
- *Adaptation note* = how an item maps onto the planned MNE-Python / fsaverage / MNE-sample re-implementation. These notes are not part of the paper.

### At a glance
- **Sources:** cortical-normal dipoles of 10 nAm at every cortical vertex (about 130,000 per hemisphere). Patches of 10 or 16 mm geodesic radius (3 or 8 cm²) at 50 pAm/mm² (p. 1078).
- **Head model:** three-layer linear-collocation BEM, 5120 triangles per surface. Conductivities printed as 0.3/0.06/0.3 S/m (brain/skull/scalp; p. 1078). The skull value is probably a typo for 0.006 (§3.4).
- **SNR:** Eq. 1, a channel-averaged power SNR in dB with a 1/N factor. Only per-sensor noise variances are used (p. 1079).
- **Noise, two separate branches:**
  - *Recorded:* per-sensor variance of 2 min of spontaneous data, filtered 0.5–100 Hz (pp. 1078–1079).
  - *Modelled:* independent uniform cortical sources at ≈7 mm spacing, s_k² = s_s²(AAᵀ)_kk, with √s_s² = 1.6, 1.6, 1.6, 1.9 nAm (pp. 1079–1080).
- **Headline results:** MEG beats EEG on ≈40% of cortex (noise model) or 55% (recorded noise). The MEG-favoured area shrinks as patches get larger (pp. 1081–1082).
- **Absent from the paper:** separate magnetometer/gradiometer results, any positioning experiment, and the forward-software name or version.

---

## 1. Bibliographic information and aim

- **Citation:** Goldenholz DM, Ahlfors SP, Hämäläinen MS, Sharon D, Ishitobi M, Vaina LM, Stufflebeam SM. Mapping the signal-to-noise-ratios of cortical sources in magnetoencephalography and electroencephalography. *Hum Brain Mapp* 2009;30(4):1077–1086. DOI 10.1002/hbm.20571 (p. 1077).
- **History:** received 13 February 2007; revised 17 September 2007; accepted 22 February 2008; published online 8 May 2008 (p. 1077).
- **Affiliations:** MGH Athinoula A. Martinos Center for Biomedical Imaging; Boston University Biomedical Engineering; Harvard-MIT Division of Health Science and Technology; Harvard Medical School Neurology (p. 1077).
- **Funding:** NCRR P41RR14074; NIH NS37462 and NS44623; MIND Institute (p. 1077).
- **Aim (paraphrased):**
  - The authors asked where on the cortex MEG or EEG gives the higher signal-to-noise ratio. Each subject's MRI-reconstructed cortex constrained source position and orientation (pp. 1077–1078).
  - In four subjects they computed per-location SNR maps for focal dipoles and extended patches. Noise came either from recorded spontaneous activity or from a model of random, uniformly distributed cortical sources, and the two modalities were compared with a difference map (pp. 1078–1080).
  - One epilepsy patient illustrated how the maps can explain spikes that are clearer in one modality (p. 1081).

---

## 2. Participants, anatomy, MRI, surfaces and source spaces

### 2.1 Participants
- **Subjects:** four right-handed, two female, aged 23–42 years (p. 1078). Consent and privacy forms were signed under the institution's human-subjects board and HIPAA rules (p. 1078).
- **Not stated:** individual ages and sexes, head sizes, health status.
- The clinical example is a separate person, a 38-year-old woman (§8.6; p. 1081).

### 2.2 MRI (p. 1078)
- **Scanner:** Siemens Trio 3T.
- **MPRAGE:** voxel size 1.3 × 1.0 × 1.3 mm³; slice thickness 1.3 mm; TE 3.31 ms; TR 2530 ms; gap 50%; FOV 10 cm × 10 cm.
- **FLASH:** flip angle 5°; slice thickness 1.3 mm; FOV 10 cm × 10 cm.
- **Not stated:** TI, matrix size, number of FLASH echoes, any second flip angle, and what the FLASH images were used for (see A18).

### 2.3 Cortical surfaces (p. 1078)
- Reconstructed with FreeSurfer (cites Dale 1999; Fischl 1999a; Fischl 2001; Segonne 2004). FreeSurfer version: *not stated*.
- About 130,000 vertices per hemisphere, with an approximate vertex-to-vertex spacing of 1 mm.
- **Which surface carries the sources (white or pial): *not stated*.**
  - Orientation is defined relative to the gray–white matter boundary (p. 1078).
  - **[inference]** This points to the white surface, which is also the MNE-C default source surface **[external: MNE manual v2.7.3, p. 22]**.
- Figures 2–4 use inflated surfaces for display (pp. 1080, 1082); Figure 5 also shows pial-surface views (p. 1083).

### 2.4 Signal source space (p. 1078)
- One current dipole sits at **every** cortical-surface vertex (full resolution).
- Orientation is **fixed perpendicular to the gray–white matter boundary** (cortical-normal constraint; no loose or free orientation).
- Both hemispheres were computed: area percentages are averaged across subjects and hemispheres (p. 1081). All figures show only the left hemisphere (pp. 1080, 1082).
- The sign convention of the normals is *not stated*. It does not matter for a single dipole's power SNR, but it does matter when patch elements are summed (§5.3).

### 2.5 Noise source space (pp. 1079–1080)
- The dense cortical triangulation was decimated to an approximate inter-source spacing of **7 mm**.
- Sources are oriented perpendicular to the cortical mantle and spread uniformly.
- **Not stated** (see A3): the number of noise sources, the decimation algorithm, whether vertex normals or patch-averaged normals were used, and whether both hemispheres and the medial wall were included.

### 2.6 Depth and curvature handling
- The paper uses no depth weighting, depth metric, curvature correction, or depth/orientation binning.
- Depth and orientation act only through the forward model and are described qualitatively (pp. 1080–1081, 1083).
- Fig. 2 row 4 shows curvature (gyri light grey, sulcal walls dark grey, on inflated and folded surfaces) purely to help interpretation (p. 1080).
- *Radial* is used qualitatively for sources oriented normal to the nearby or inner skull surface, such as gyral crests (pp. 1081, 1083).

### 2.7 Group averaging and display (p. 1079; captions pp. 1080, 1082)
- Grand averages were made by spherical morphing that aligns sulcal and gyral features (Fischl 1999b). They are displayed on Subject 2's inflated left hemisphere.
- **Not stated** (A7): the common target space, whether dB or linear values were averaged, and any smoothing.
- Fig. 4 shows the number of subjects (0–4) with D > 0 at each location (pp. 1081–1082), even though its caption calls it a group average.

**Adaptation notes** (measured locally from `data/external/MNE-sample-data/`, not from the paper):
- **MNE *sample* white surface** (`sample-all-src.fif`): 155,407 (lh) and 156,866 (rh) vertices; area ≈1006 and 1018 cm²; mean edge ≈0.88 mm. This is comparable to the paper's ≈130,000 vertices at ≈1 mm.
- ***fsaverage*** (`fsaverage-ico-5-src.fif`): 163,842 vertices per hemisphere; white-surface area ≈654 and 650 cm². fsaverage is a smoother average brain, so cancellation inside patches will probably be weaker than in individual anatomy **[inference]**.
- **Missing fsaverage pieces:** the fsaverage copy inside the sample dataset has no BEM files, and fsaverage has no MEG sensor geometry. An fsaverage run needs borrowed sensor positions and a head–MRI transform.
- **Environment:** the project `.venv` (MNE-Python 1.13.2) lacks `nibabel`, so FreeSurfer surfaces cannot be read directly yet (`read_surface`, `setup_source_space` and `grow_labels` need it).

---

## 3. Head model

### 3.1 As printed (p. 1078)
- **Method:** a **linear collocation, three-layer boundary-element method** (cites Hämäläinen & Sarvas 1989; Hämäläinen et al. 1993).
- **Surfaces:** skin (scalp), outer skull and inner skull, each derived from the MRIs (Fig. 1D–F, p. 1079).
- **Tessellation:** **5120 triangles** per surface, described as numerically adequate (cites Crouzeix 1999; de Jongh 2005; Fuchs 2001; Tarkiainen 2003).
- **Conductivities, verbatim:** "0.3, 0.06, and 0.3 S/m for the brain, skull and scalp" (p. 1078). That is, brain 0.3 S/m, skull 0.06 S/m, scalp 0.3 S/m.
- **Not stated:** how the surfaces were segmented (watershed or FLASH), how the isolated-problem issue was handled, and any inner-skull shift.

### 3.2 Software
- FreeSurfer is named for the cortical reconstruction (p. 1078).
- The BEM/forward software and its version are *not stated*. MNE is not mentioned anywhere in the paper.
- **[inference + external]** The described setup matches the documented defaults of MNE-C v2.5 (December 2006), whose author M. Hämäläinen is a co-author here:
  - linear collocation is marked as the preferred BEM method (manual p. 81);
  - `--ico 4` gives 5120 triangles per surface (manual p. 25);
  - the default source spacing is 7 mm (manual p. 21).
- MNE-C of that era is therefore likely, but unconfirmed.

### 3.3 The paper's own conductivity discussion (p. 1081)
- The conductivities are presented as standard literature values, citing Hämäläinen et al. 1993. The EIT-derived effective conductivities of de Jongh et al. (2005) are mentioned but not adopted.
- Conductivity errors are called a minor issue for MEG (Hämäläinen & Sarvas 1989) but potentially important for EEG, through the skull.
- The noise was either measured or matched to measured variances, so conductivity errors are argued to have little effect on the noise estimates. Recomputing s_s without EEG changed it by less than 5% on average.
- The paper says a true skull conductivity higher than modelled would make its EEG SNR estimates too large, and a lower one too small. This direction looks questionable (see A13).
- No conductivity sensitivity analysis is reported.

### 3.4 The printed skull conductivity of 0.06 S/m: verification and assessment

**(a) What the PDF says.**
- Verified on a 300-dpi render: the skull value is printed as 0.06 S/m, with 0.3 S/m for brain and scalp (p. 1078).
- This is the only place any conductivity value appears. No conductivity ratio is printed anywhere.

**(b) Internal consistency.**
- **Implied ratio:** the printed values give brain:skull = 5:1 (skull = 1/5 of brain), with scalp = brain **[derived]**.
- **Not a standard ratio:** the paper calls these standard values (p. 1081). A 1:5 ratio is not what this software lineage or this group used at the time (see c) **[external]**.
- **No back-calculation possible:** no other number in the paper lets the ratio be recovered. The EEG maps would be affected, but no EEG amplitudes are reported.
- **Loose handling of the topic:** the paper's only other conductivity statement (direction of the EEG error, p. 1081) appears reversed (A13). This says nothing about the value itself.
- **Other editorial slips (weak evidence that a dropped zero is plausible):**
  - citation forms: de Jong et al., 2005; de Jongh, 2005; Tarkiainen, 2003 (all p. 1078);
  - spelling: Leitjen (p. 1083) vs Leijten (pp. 1077, 1085); Lu and Willamson (p. 1084); *are* for *area* and *combing* for *combining* (p. 1084); Numeriche, Electroenchalography (p. 1085);
  - the Hämäläinen et al. 1993 reference is printed as Rev Modern Phys 65:1–93 (p. 1085), whereas the published pagination is 65:413–497 **[external]**.

**(c) External evidence [external].**
1. **MNE User's Guide v2.5** (Dec 2006, current when the paper was submitted in Feb 2007), `mne_setup_forward_model`:
   - brain and scalp defaults are 0.3 S/m;
   - the skull default is printed as 0.06 S/m **in the same sentence that calls it a 1/50 brain-to-skull ratio** (manual p. 26);
   - this contradicts itself, since 0.3/50 = 0.006;
   - the same manual's `mne_surf2bem` example uses 0.00375 S/m (= 0.3/80; §5.4.5).
2. **MNE User's Guide v2.7** (Dec 2009) and **v2.7.3** (Nov 2010): the same option is documented as 0.006 S/m with the same 1/50 ratio (p. 27 in both). The typo was corrected.
3. **MNE-Python 1.13.2** (project `.venv`, `mne/bem.py`): `make_bem_model(subject, ico=4, conductivity=(0.3, 0.006, 0.3), ...)`. `make_bem_solution` implements the linear-collocation BEM ported from MNE-C.
4. **MNE sample dataset BEM** (`sample-5120-5120-5120-bem.fif`, `...-bem-sol.fif`): 5120 triangles per surface; σ = 0.3 (scalp), 0.006 (skull), 0.3 (brain).
5. **Same group, later paper:** Ahlfors et al. (2010), Brain Topogr (PMC2914866), assumed brain:skull:scalp = 1:0.0125:1 (1/80), citing Geddes & Baker 1967.

**(d) Assessment.**
- The printed triple 0.3/0.06/0.3 S/m matches the MNE v2.5 manual's documented default exactly, including that manual's own typo.
- Everything else points to **0.006 S/m**: the 1/50 ratio stated in that same manual, all later manuals, MNE-Python, and the MGH sample BEM.
- The most probable explanation: the authors used the default MNE three-layer BEM (skull 0.006 S/m, ratio 1:50) but copied its value from the v2.5 manual text, carrying the typo into the paper.
- This cannot be proven from the paper. A deliberate 0.06 S/m (1:5) cannot be ruled out, but the paper gives no justification for it, cites no study supporting a high skull conductivity, and calls the values standard.
- **Classification:** probable typographical error propagated from the documentation. Confidence moderate to high; not confirmable without the authors' files.

**(e) Implications and recommendation for the reproduction.**
- **Effect on MEG is small:**
  - in a single-compartment (inner-skull) BEM the magnetic field does not depend on the conductivity value;
  - in a three-layer BEM, MEG depends only on conductivity ratios and is only weakly sensitive (paper p. 1081).
- **Where it matters:** mainly for EEG (out of scope), and marginally for any re-derivation of s_s that includes EEG channels (<5% per p. 1081).
- **Recommendation:**
  - record 0.06 S/m unchanged as the printed value;
  - run the three-layer MEG model twice, labelled *as printed* (0.3/0.06/0.3) and *probable software default* (0.3/0.006/0.3), and report the MEG SNR difference;
  - optionally add MNE's recommended single-layer MEG model as a third check.

---

## 4. Sensors, acquisition, coregistration and head position

### 4.1 MEG
- Elekta-Neuromag Vectorview™ with **204 planar gradiometers and 102 magnetometers** (p. 1078), called a 306-channel system (p. 1082).
- **Not stated:** coil definitions or integration points, sensor noise levels, excluded channels.
- Instrumentation noise is not part of the noise model (pp. 1079–1080).

### 4.2 EEG (out of scope; shared definitions only)
- **70 electrodes** (p. 1078; 70-channel system, p. 1082), transformed to the **average electrode reference** (p. 1078).
- Cap, layout and electrode names: *not stated*. Fig. 1B shows the digitized positions (p. 1079).
- Fig. 5 uses a clinical bipolar montage for display only (p. 1083).

### 4.3 Acquisition (p. 1078)
- MEG and EEG were recorded simultaneously: **two minutes of spontaneous activity**, sampled at **600 Hz**.
- Hardware filters: **0.1–200 Hz for Subjects 1 and 2**, **0.03–200 Hz for Subjects 3 and 4**.
- Posture, eyes open or closed, and vigilance: *not stated*.
- **SSP was applied to magnetometer data only** (Tesche 1995):
  - the noise subspace was spanned by the eigenvectors of the three largest eigenvalues of the correlation matrix of a **5-min empty-room recording**;
  - these components approximate homogeneous fields from distant environmental sources.
- No SSP is mentioned for gradiometers or EEG, and SSS/MaxFilter is not mentioned.

### 4.4 Digitization, coregistration and head position (p. 1078)
- Fiducial points, EEG electrodes and HPI coils were digitized with a Polhemus 3D digitizer (printed as *FastTrack*).
- The MEG array's position relative to the head was determined **once, at the start of each measurement**, from the HPI coil fields.
- **Not stated:** extra head-shape points, the MRI-to-head alignment method, continuous head tracking or movement compensation, head-to-sensor distances.
- Fig. 1A shows the sensor array relative to Subject 2's head (p. 1079).

### 4.5 Positioning: experiments and statements
- **No positioning or head-position experiment was performed.**
- **p. 1083:** MEG SNR may depend critically on where the sensor array sits relative to the head. Marinkovic et al. (2004) showed large differences in signals from frontal sources.
- **p. 1084:** SNR depends critically on head position. SNR maps could be used to check whether the head position was adequate, and a focus found in an unexpected region should be checked for sufficient SNR.
- **p. 1082:** MEG and EEG channel counts and positions were not matched (306 vs 70). Both arrays are argued to sample adequately, so the exact count and placement should matter little.
- **p. 1081:** low EEG SNR in inferior frontal areas is blamed on inadequate electrode coverage. This is in tension with p. 1082 (A14).

**Adaptation notes** (MNE sample, measured locally):
- **Channels:** `sample_audvis_raw.fif` has 204 gradiometers, 102 magnetometers and 60 EEG channels (plus 9 stim and 1 EOG). Bad channels are MEG 2443 and EEG 053.
- **Acquisition:** 600.615 Hz, band 0.1–172.18 Hz.
- **SSP:** three magnetometer-only PCA projectors (PCA-v1 to v3, 102 channels each), also present in `ernoise_raw.fif`.
- **Recordings:** the raw file is a ≈277.7-s auditory/visual task run, not spontaneous activity. The empty-room file lasts ≈110 s (paper: 5 min).

---

## 5. Source models

### 5.1 Focal sources (p. 1078)
- A current dipole at each cortical vertex, oriented perpendicular to the gray–white boundary, with **amplitude a = 10 nAm**.

### 5.2 Extended sources (p. 1078; p. 1084; Fig. 4, p. 1082)
- **Definition:** synchronously active cortical patches (p. 1078).
- **Strength:** a **uniform surface source (moment) density of 50 pAm/mm²** (cites Hillebrand & Barnes 2002; Lu & Williamson 1991) (p. 1078).
- **Size:** **geodesic radius 10 mm or 16 mm**, corresponding to **areas of 3 cm² or 8 cm²** (p. 1078; restated p. 1084).
- **Construction:** the cortical surface was divided into two complete sets of centroids using a geodesic-distance-weighted Dijkstra algorithm (Dijkstra 1959) (p. 1078).
  - **[inference]** One centroid set per radius covers the whole cortex. Each patch holds the vertices within that geodesic radius of its centroid.
- **Element orientation and summation:** not stated explicitly.
  - **[inference]** Each element is a cortical-normal dipole, as in the focal case, and the patch field is their signed vector sum.
  - This fits the paper's explanation that the shrinking MEG advantage for larger patches reflects cancellation of tangential components (p. 1084, citing Eulitz 1997).
- **Noise model:** kept unchanged when source extent varied (p. 1084).
- **Other sizes:** none were simulated. The Discussion cites literature values: 6–10 cm² for EEG detectability, and in the mesial temporal lobe 4 cm² invisible vs 8 cm² visible to MEG (Baumgartner 2000) (p. 1084).

### 5.3 Patch details that are *not stated*
- The surface used for geodesic distance (white or pial), and whether distances are Dijkstra edge-path lengths or true geodesics.
- How many centroids there are, how far apart, and whether patches overlap.
- How the 50 pAm/mm² density is spread over vertices (area weighting), and whether the total moment uses the nominal area (πr²) or the actual patch area.
- What *a* and *b_k* in Eq. 1 mean for a patch.
- How patch SNR values are mapped back onto vertices for Fig. 4, and how patches crossing the medial wall are handled.

### 5.4 Derived checks [derived]
- **Areas:** π(10 mm)² = 314 mm² ≈ 3 cm² and π(16 mm)² = 804 mm² ≈ 8 cm², so the printed areas are flat-disc areas.
- **Nominal total moments (no cancellation):** 50 pAm/mm² × 300 mm² = 15 nAm and × 800 mm² = 40 nAm (15.7 and 40.2 nAm using πr²). Cortical folding makes the net dipole smaller.
- **Size effect:** pure area scaling between the two patch sizes predicts 20·log10(8/3) ≈ 8.5 dB. The paper reports 10 dB in the mesial temporal lobe, and an EEG range of <1 dB to >10 dB (p. 1084).

---

## 6. Noise

### 6.1 Modelled background noise (the noise model) (pp. 1079–1080)
- **Purpose:** to reduce the influence of session-specific sensor quirks (p. 1079).
- **Model:** identical, independent noise sources spread uniformly over the cortex (de Munck et al. 1992), oriented perpendicular to the cortical mantle. Amplitudes are Gaussian with zero mean and variance s_s² (p. 1079).
- **Left out:**
  - instrumentation noise, justified as much smaller than spontaneous brain activity (cites Ebersole & Pedley 2002; Hämäläinen 1993) (pp. 1079–1080);
  - environmental noise (p. 1082).
- **Eq. 3:** s_k² = s_s²(AAᵀ)_kk. A is the forward matrix for unit-amplitude noise sources on the cortex decimated to ≈7 mm spacing; ᵀ is the transpose (p. 1080).
- **Number of noise sources:** *not stated* (A3).
- **Covariance:** the sources are independent, so the implied sensor covariance is s_s²AAᵀ. Only its diagonal enters Eq. 1 **[inference from Eq. 1]**.
- **Calibrating s_s² (p. 1080):**
  1. For each sensor type separately (MEG gradiometers, MEG magnetometers, EEG), take the median over channels of ŝ²_k,recorded / (AAᵀ)_kk. Medians guard against bad or outlier channels.
  2. Average the three medians, weighting each by the number of channels of that type.
  3. As a formula (my transcription of the prose): s_s² = Σ_t n_t · median_{k∈t}[ŝ²_k,recorded/(AAᵀ)_kk] / Σ_t n_t, with t ∈ {grad, mag, EEG}. The counts are presumably n = 204, 102, 70, giving weights ≈0.543, 0.271, 0.186 **[derived]**.
  4. Sensor variances for the model then follow from Eq. 3 using that estimate.
- **Results:** √s_s² = 1.6, 1.6, 1.6 and 1.9 nAm for the four subjects (p. 1080). Leaving out EEG changed s_s by <5% on average (p. 1081).
- Whether each subject's own s_s or a pooled value was then used is ambiguous (A4).
- The same noise model serves both dipoles and patches (p. 1084).

### 6.2 Recorded spontaneous-activity noise (pp. 1078–1079)
- **Data:** 2 min per subject at 600 Hz, hardware band 0.1–200 or 0.03–200 Hz; EEG average-referenced; magnetometer SSP applied (p. 1078).
- **Estimate:** the per-sensor variance ŝ²_k,recorded of the data filtered to **0.5–100 Hz** (p. 1079). No full covariance matrix is used.
- **Not stated:** filter type, order and phase; which segment was used; artifact rejection; eye state; how bad channels were handled in the recorded-noise maps.
- **Used for:**
  1. the recorded-noise D maps (Fig. 3 bottom row; the 55% figure) (pp. 1081–1082);
  2. calibrating s_s (p. 1080);
  3. the patient's D maps (Fig. 5) (pp. 1081, 1083).
- Transient artifacts (blinks, muscle, cardiac) are left out of the modelling because they can be excluded or suppressed (p. 1083). Whether they were removed from the recorded data is *not stated*.

### 6.3 Sensor, empty-room and environmental noise
- The 5-min empty-room recording was used **only** to build the three-component magnetometer SSP subspace (p. 1078).
- Instrument and environmental noise are absent from the noise model (pp. 1079–1080, 1082). They are implicitly present in the recorded-noise branch, after magnetometer SSP **[inference]**.
- Neither branch adds a separate empty-room noise term.
- **Adaptation caution [inference]:** the model has no sensor-noise term. Comparing sensor technologies with different intrinsic noise (for example OPM vs SQUID) therefore needs one added. That would be an extension, not Goldenholz's method.

### 6.4 How each noise enters the SNR
- In Eq. 1, s_k² is either ŝ²_k,recorded (recorded branch) or s_s²(AAᵀ)_kk (noise-model branch). The two are never combined.
- **Noise model:** Fig. 2, Fig. 3 top and middle rows, Fig. 4, and the ≈40% figure (pp. 1080–1082).
- **Recorded noise:** Fig. 3 bottom row and the 55% figure (pp. 1081–1082), plus the patient's Fig. 5. The patient's noise-model maps are not shown (p. 1081).

---

## 7. SNR definition and summary statistics

### 7.1 Equations (transcribed; Eqs. 1–2 on p. 1079, Eq. 3 on p. 1080)
```
(1)  SNR   = 10 · log10 [ (a² / N) · Σ_k ( b_k² / s_k² ) ]
(2)  D     = SNR_MEG − SNR_EEG
(3)  s_k²  = s_s² · (A Aᵀ)_kk
```
LaTeX: $\mathrm{SNR}=10\log_{10}\!\left[\frac{a^{2}}{N}\sum_{k}\frac{b_{k}^{2}}{s_{k}^{2}}\right]$, $D=\mathrm{SNR}_{\mathrm{MEG}}-\mathrm{SNR}_{\mathrm{EEG}}$, $s_k^2=s_s^2\,(\mathbf{A}\mathbf{A}^{T})_{kk}$.

In Eq. 3, A is the forward-solution matrix for unit-amplitude noise sources and ᵀ is the transpose (p. 1080).

**Printing note:** the brackets of Eq. 1 are typeset as floor-shaped glyphs (⌊ ⌋). Read them as ordinary brackets. A floor function would turn every linear SNR below 1 into −∞ dB, whereas Fig. 2 shows finite values between −29 and −19 dB.

### 7.2 Definitions and properties
- SNR is defined for each source location (dipole or patch) and expressed in dB (p. 1079).
- **Symbols (p. 1079):** b_k is the forward-model signal on sensor k for a **unit-amplitude** source; a is the source amplitude; **N is the number of sensors**; s_k² is the noise variance of sensor k.
- Eq. 1 is therefore a **channel-averaged power SNR with a 1/N factor**: 10·log10 of the mean over channels of (a·b_k)²/s_k².
- It compares a single sample (the squared source amplitude) with the broadband 0.5–100 Hz noise variance **[inference]**. Whether *a* is a peak or RMS amplitude is *not stated*.
- **Trial averaging:** with additive Gaussian noise, averaging trials adds the same constant to MEG and EEG SNR, so D is unchanged (p. 1079). The constant is +10·log10(n_trials) dB **[derived]**.
- **D:** does not depend on a (p. 1079). D > 0 means MEG has the higher SNR; D < 0 means EEG does (p. 1079).
- **Combined map:** max{SNR_MEG, SNR_EEG} at each location (p. 1079). It is not a joint Eq. 1 over all 376 channels.

### 7.3 Per-channel-type handling and whitening
- Sensor types are treated separately **only** when calibrating s_s: gradiometer, magnetometer and EEG medians (p. 1080).
- **MEG channel set:** which channels, and what N, go into SNR_MEG is *not stated*.
  - **[inference]** All 306 MEG channels are pooled into one average (N = 306). This is dimensionally valid because each term b_k²/s_k² is unitless.
  - **No separate magnetometer or gradiometer SNR maps or values are reported.**
- **EEG:** average reference (p. 1078); N presumably 70 **[inference]**.
- **Whitening:** none beyond dividing each channel by its own variance (diagonal normalization). Noise correlations between channels are ignored, and no regularization or post-SSP rank handling is described **[inference from Eq. 1]**.
- **Average vs maximum SNR (pp. 1081–1082):**
  - The authors chose the channel average (Fuchs 1998) over the maximum-channel SNR (de Jongh 2005), because source localization benefits from many high-SNR sensors.
  - Average and maximum agree more closely for spatially widespread patterns, so averaging biases D toward EEG.
  - Differential sensors (bipolar EEG, planar gradiometers) give more focal patterns than referential ones (average-reference EEG, magnetometers). This bears directly on any magnetometer-vs-gradiometer comparison under Eq. 1.

### 7.4 Maps and summary statistics
- **Group maps:** grand averages via spherical morphing (p. 1079).
- **Area statistic:** the percentage of cortex where MEG SNR > EEG SNR, computed as the **fraction of sources with D > 0** (p. 1079). It is reported averaged across subjects and hemispheres, excluding the *ventricular surfaces* (p. 1081).
- **Consistency map:** the number of subjects (0–4) with D > 0 at each location (Fig. 4; pp. 1081–1082).
- **Not reported:** any absolute-SNR threshold statistic (such as % of cortex above X dB), confidence intervals, statistical tests, or SNR distributions.

### 7.5 Derived identities useful for MEG-only work [derived]
- **Pooling:** if SNR_MEG pools all 306 channels, then 10^(SNR_MEG/10) = [102·10^(SNR_mag/10) + 204·10^(SNR_grad/10)] / 306.
- **Noise-model split:** SNR = 10·log10(a²/s_s²) + 10·log10[(1/N) Σ_k b_k²/(AAᵀ)_kk].
  - With a/s_s = 10/1.6 = 6.25, the first term is +15.9 dB; with 10/1.9 it is +14.4 dB. Subject 4 therefore sits ≈1.5 dB lower for the same anatomy.
  - Under a shared noise model, the difference between two sensor sets (D, or magnetometer vs gradiometer) is independent of both a and s_s.
  - Absolute SNR levels are not: they also depend on the noise-source density (A3).
- **Order-of-magnitude check:** if all M noise sources contributed equally, the second term would be about −10·log10 M. With M ≈ 4,000 (≈2,000 per hemisphere at 7 mm; see A3), an average source comes out near −20 dB. This fits Fig. 2's −29 to −19 dB display range.

---

## 8. Results

### 8.1 Dipole SNR maps (Fig. 2, p. 1080; text pp. 1080–1081)
- **Setup:** noise model, averaged over the four subjects, left hemisphere shown on Subject 2's inflated surface. The rows show SNR_MEG, SNR_EEG, their maximum, and curvature.
- **Colour scale:** the same for every row, saturating at **< −29 dB** (cyan) and **> −19 dB** (yellow). These are display limits, not reported ranges.
- **MEG:**
  - high SNR over most superficial cortex;
  - lower SNR in lateral sulci (especially the Sylvian fissure), interhemispheric cortex and parts of ventral cortex, lowest near the centre of the head;
  - thin low-SNR strips along gyral crests, where source orientation is roughly normal to the nearby skull surface.
- **EEG:** much more uniform than MEG; low in inferior frontal areas, which the authors attribute to electrode coverage.
- **Combined (maximum) map:** the modalities complement each other, so combining them raises the overall SNR (p. 1081).
- Visual reading **[inference]:** much of the medial MEG map sits at or below the −29 dB floor.

### 8.2 Difference maps (Fig. 3, p. 1082; text p. 1081)
- **Colour scale:** D saturates at **< −4 dB** and **> 4 dB**. The middle and bottom rows are thresholded at D = 0 (red = MEG higher, green = EEG higher).
- **EEG higher:** much of the medial surface, Sylvian cortex, insula, sulci and narrow gyral crests.
- **MEG higher:** frontal and occipital poles and most of the lateral surface.
- **Share of cortex with D > 0:**
  - noise model: **approximately 40%** (mean across subjects and hemispheres, ventricular surfaces excluded);
  - recorded noise: **55%** (p. 1081).
- The recorded-noise SNR maps (not shown) indicate that the extra MEG-favoured area comes from MEG SNR being higher than the model predicts. The paper concludes that the noise model **underestimates MEG SNR** (pp. 1077, 1081, 1084).

### 8.3 Extended sources (Fig. 4, p. 1082; text pp. 1081, 1084)
- Fig. 4 maps, for dipoles, 10-mm patches and 16-mm patches, how many subjects have D > 0 at each location.
- **Dipoles:** at many superficial locations MEG is better in all four subjects.
- **Larger patches:** the MEG-favoured area shrinks. MEG stays consistently better at the occipital pole, inferior occipital lobe, temporal pole, posterior temporal lobe and frontal pole (p. 1081).
- **No area percentages are given for patches.**
- **Patch-size effect (p. 1084):**
  - 3-cm² sources had 10 dB lower SNR than 8-cm² sources in the mesial temporal lobe. The modality is not stated (A15).
  - The EEG SNR difference between the two sizes ranged from less than 1 dB to more than 10 dB across locations.

### 8.4 Noise-model parameters
- √s_s² = 1.6, 1.6, 1.6, 1.9 nAm (p. 1080); <5% change without EEG (p. 1081).

### 8.5 Depth and orientation (pp. 1080–1081, 1083)
- Depth dominates. Orientation hurts MEG only in thin strips (radial sources at gyral crests).
- The paper cites ∼5% of cortex as affected by orientation-related MEG attenuation (Hillebrand & Barnes 2002). This figure was not computed in the study.
- EEG is better for deeper sources (medial cortex, insula); MEG is often better for superficial ones (p. 1083).

### 8.6 Clinical example (p. 1081; Fig. 5, p. 1083)
- **Patient:** 38-year-old woman with a hypothalamic hamartoma and intractable complex partial seizures since age 2. She had frequent, independent left and right frontotemporal discharges.
- **Spikes:** two interictal spikes were chosen, A clearer in EEG and B clearer in MEG. Equivalent current dipoles were fitted to EEG for A and to MEG for B.
- **Outcome:**
  - at the dipole locations, the recorded-noise D maps predicted the observed modality difference for both spikes;
  - the noise-model D maps (not shown) did so only for spike A;
  - the authors conclude the model overestimated MEG noise for that session and location.
- **Fig. 5 contents:**
  - EEG bipolar traces (scale bars 100 uV and 0.5 s) with potential maps;
  - left and right frontotemporal planar-gradiometer traces (500 fT/cm, 0.5 s) with field maps;
  - D-map insets showing the dipoles (green = EEG higher, red = MEG higher), labelled EEG-favoured for A and MEG-favoured for B.
- **Not stated:** the patient's recording parameters, dipole-fitting details, and the D values at the dipoles.

### 8.7 Not reported anywhere in the paper
- Magnetometer vs gradiometer SNR.
- Numerical SNR distributions or regional values (only colour-bar limits).
- Per-subject area percentages, and area percentages for patches.
- Head-position values.
- Any conductivity sensitivity result.

---

## 9. Limitations of the metric and of positioning (pp. 1081–1084)
1. **Conductivity (p. 1081):**
   - skull-conductivity errors matter for EEG but little for MEG;
   - the noise estimates are anchored to data (s_s changed <5% without EEG);
   - the stated direction of the EEG error is questionable (A13).
2. **Average vs maximum SNR (pp. 1081–1082):** channel averaging favours spatially widespread patterns, which biases D toward EEG. Differential sensors (planar gradiometers, bipolar EEG) are more focal than referential ones (magnetometers, average-reference EEG).
3. **Sensor count and placement (p. 1082):** 306 MEG vs 70 EEG channels, not matched. The authors argue both sample the field adequately, so the effect is small. This sits uneasily with the remark on inferior-frontal EEG coverage (p. 1081).
4. **Scope of the noise model (pp. 1082–1083):**
   - only uniform cortical sources are modelled; instrument and environmental noise are excluded, and session-specific measured noise is preferable;
   - the model and recorded results may differ partly because focal or regional background activity (such as alpha) reaches more EEG sensors than MEG sensors, so in EEG it looks more like distributed noise;
   - transient artifacts (blinks, muscle, cardiac) are excluded.
5. **Orientation vs depth (p. 1083):** radial sources on gyral crests are hard for MEG, but orientation matters for only ∼5% of cortex (a cited figure). Depth explains most of the differences.
6. **Positioning (pp. 1083–1084):**
   - MEG SNR may depend critically on where the array sits relative to the head, especially for frontal sources (Marinkovic 2004);
   - SNR maps are proposed as a head-position check and as a reliability marker, since low-SNR sources (such as deep ones) are more prone to mislocalization;
   - **no positioning experiment was done.**
7. **Extended sources (p. 1084):**
   - only the signal extent was varied; the noise model stayed dipolar and uniform;
   - adding extended sources to the noise model would be expected to offset the shift and raise MEG's relative SNR;
   - the same trend is expected with recorded noise, which does not depend on the source model;
   - variability from location to location is large.
8. **Trial averaging (p. 1079):** under additive Gaussian noise it shifts all SNRs equally and leaves D unchanged.
9. **Noise model vs data (pp. 1077, 1081, 1084):** the model underestimates MEG SNR compared with recorded noise, and it failed to predict spike B in the patient.
10. **Further limits implied by Eq. 1** [inference; not stated by the authors]:
    - noise is diagonal only, so this is not a whitened or optimal detectability measure;
    - it compares a single-sample amplitude with broadband variance, ignoring the source's time course and spectrum;
    - head position is measured only once per session.

---

## 10. Ambiguities and inconsistencies a re-implementer must resolve

Items are labelled A1–A21; other sections refer to them by these labels.

- **A1. Skull conductivity of 0.06 S/m (p. 1078).**
   - Probably a documentation typo for 0.006 S/m copied into the paper (§3.4).
   - *Resolve:* keep 0.06 as the printed value; run the three-layer MEG model with both 0.06 and 0.006 and label each.
- **A2. MEG channel set and N in Eq. 1 (p. 1079).**
   - Whether the 306 channels were pooled or split by sensor type is not stated, and there are no magnetometer or gradiometer results to validate against.
   - *Resolve:* compute SNR_mag (N = 102), SNR_grad (N = 204) and pooled SNR_MEG (N = 306); §7.5 gives the pooling identity.
- **A3. Noise-source density (p. 1080).**
   - Only an approximate 7-mm spacing is given. The number of noise sources, the decimation method, vertex vs patch-averaged normals, and whether both hemispheres and the medial wall are included are all unstated.
   - The printed s_s (1.6–1.9 nAm) is calibrated to that unknown density. Absolute SNR depends on it; D-type contrasts do not.
   - *Resolve:* use MNE `setup_source_space(spacing=7)` (integer-mm spacing, available since MNE 0.18) and report the resulting number of sources M. Otherwise, hold the noise variance per unit area fixed by scaling s_s² by (area per new source)/(area per 7-mm source).
   - For scale: roughly 2,000 sources per hemisphere for ≈1,000 cm², extrapolated from MNE manual Table 3.1 (6.2 mm spacing ↔ 39 mm² per source) **[derived/external]**.
- **A4. Per-subject or pooled s_s (p. 1080).**
   - The sentence about using the overall estimate is ambiguous, and the listed values 1.6, 1.6, 1.6, 1.9 nAm are presumably in subject order 1–4.
   - *Resolve:* for a single-anatomy adaptation, use 1.6 nAm as primary and 1.9 nAm as a sensitivity case, or re-calibrate on the sample data (labelled as an adaptation).
- **A5. Patch construction (p. 1078).** Every item in §5.3 is unstated. *Resolve:*
   - grow patches with Dijkstra on the white surface (for example `mne.grow_labels(subject, seeds, extents=10 or 16, hemis, surface='white')`);
   - set each vertex's moment to 50 pAm/mm² × vertex area, orient elements along consistently signed cortical normals, and sum them as vectors;
   - report actual patch areas next to the nominal 3 and 8 cm², and state the centroid spacing and overlap policy.
- **A6. Source surface (p. 1078).** White or pial is not stated. *Resolve:* use the white surface (MNE default) and say so.
- **A7. Group averaging (p. 1079).**
   - Averaging in dB vs linear units, the morph target (maps are shown on Subject 2) and any smoothing are not stated.
   - *Resolve:* declare the choice. A single template (fsaverage or sample) allows no cross-subject averaging; say so.
- **A8. Area statistic (pp. 1079, 1081).**
   - It is a fraction of sources, not an area-weighted measure, and the source set it was computed on is not stated.
   - The excluded ventricular surfaces are not defined **[inference: the medial wall]**.
   - The averaging behind the 55% figure (subjects × hemispheres) is not restated.
   - *Resolve:* report both vertex-count and area-weighted fractions, and exclude the FreeSurfer medial wall explicitly.
- **A9. SSP and the forward model (p. 1078).**
   - Whether the magnetometer SSP was also applied to b_k and A is not stated, nor how N and rank were handled after projection (magnetometer rank drops to 99).
   - *Resolve:* apply the same projector to the gains and the noise; keep N = 102 channels and note the reduced rank.
- **A10. Recorded-noise preprocessing (pp. 1078–1079).**
    - Filter design, the segment used, artifacts, eye state and bad channels are all unstated.
    - The MNE sample dataset has only a task recording, so the recorded-noise branch can at best be an adaptation, for example using 0.5–100 Hz filtered baseline or inter-stimulus data.
- **A11. Eq. 1 bracket glyphs (p. 1079).** Printed as floor brackets; read them as ordinary brackets (§7.1).
- **A12. Notation (pp. 1079–1080).**
    - The paper writes s rather than σ.
    - ŝ²_k stands both for recorded variances and, later, for model-derived variances.
    - The units of the *unit amplitude* (A·m or nAm) are not given.
    - *Resolve:* use SI units consistently (a = 1e-8 A·m; s_s = 1.6e-9 A·m).
- **A13. Direction of the skull-conductivity effect (p. 1081).**
    - The paper says a higher true skull conductivity would make its EEG SNR estimates too large.
    - A more conductive skull normally produces larger scalp potentials, so with measured noise the modelled EEG SNR would be too small. The statement looks reversed, or *resistivity* was meant **[inference]**.
    - This affects EEG only.
- **A14. Sampling adequacy (p. 1081 vs p. 1082).** Inadequate inferior-frontal EEG coverage on one page, sufficient spatial sampling on the next. EEG-side only; just note it.
- **A15. Modality of the 10-dB statement (p. 1084).** The MEG literature cited just before it suggests MEG, but this is not explicit. Use it only as a soft benchmark.
- **A16. Fig. 4 caption (p. 1082).** It calls the figure a group average, but it shows subject counts.
- **A17. No numerical ranges.**
    - The Fig. 2 (−29/−19 dB) and Fig. 3 (±4 dB) limits are display clipping.
    - Validation can therefore only be qualitative, plus the ≈40% and 55% fractions, which require EEG (out of scope).
- **A18. MRI parameters (p. 1078).**
    - A 10 cm × 10 cm FOV and a 50% gap are unusual for whole-head 3-D MPRAGE; only one FLASH angle (5°) is listed; the BEM surface-extraction method is not stated.
    - This does not affect fsaverage or sample work, but these values should not be copied as settings.
- **A19. Head position (p. 1078).**
    - Measured once at the start of each session, with no values reported.
    - In the noise model, head position matters only through the forward model. A sample-subject adaptation uses the sample's own device-to-head transform.
- **A20. Hardware filter differences (p. 1078).** The 0.1 vs 0.03 Hz high-pass across subjects is superseded by the 0.5–100 Hz analysis filter (p. 1079). Record it; it does not change the method.
- **A21. Meaning of *a*.** Peak or RMS amplitude is not stated. Eq. 1 treats a² as instantaneous signal power compared with a variance; keep it that way for comparability.

---

## 11. Reproduction configuration

| Parameter | Value as printed | Page | Notes (inference / adaptation) |
|---|---|---|---|
| Subjects | 4 right-handed, 2 female, aged 23–42 | 1078 | Unavailable → fsaverage and/or MNE sample (adaptation; no cross-subject averaging) |
| MRI | Siemens Trio 3T; MPRAGE 1.3 × 1.0 × 1.3 mm³, slice 1.3 mm, TE 3.31 ms, TR 2530 ms, gap 50%, FOV 10 cm × 10 cm; FLASH flip 5°, slice 1.3 mm, FOV 10 cm × 10 cm | 1078 | Documentation only (A18) |
| Cortical reconstruction | FreeSurfer; about 130,000 vertices/hemisphere; ≈1 mm vertex spacing | 1078 | Version not stated; sample 155,407/156,866 vertices; fsaverage 163,842 |
| Signal source space | Dipole at each surface vertex | 1078 | Full resolution; MNE `spacing='all'` (heavy) or a declared dense decimation |
| Signal source orientation | Perpendicular to gray–white matter boundary | 1078 | Fixed, white-surface normals (A6) |
| Focal amplitude a | 10 nAm | 1078 | = 1e-8 A·m |
| Patch moment density | 50 pAm/mm² | 1078 | Per-vertex moment = density × vertex area (A5) |
| Patch radii (geodesic) | 10 mm, 16 mm | 1078; 1082 (Fig. 4); 1084 | Dijkstra on white surface (adaptation) |
| Patch areas | 3 cm², 8 cm² | 1078; 1084 | = πr² (314, 804 mm²); nominal totals ≈15, 40 nAm [derived] |
| Patch construction | Two complete sets of centroids; geodesic-distance-weighted Dijkstra | 1078 | Centroid spacing and overlap not stated (A5) |
| Patch element orientation/summation | Not stated | — | Inferred: cortical normal, signed vector sum (cancellation, p. 1084) |
| Other source sizes | None simulated | — | — |
| BEM method | Linear collocation, three-layer | 1078 | MNE `make_bem_solution` is linear collocation |
| BEM surfaces | Skin, outer skull, inner skull, from MRI | 1078; Fig. 1 p. 1079 | Sample has a 3-layer BEM; local fsaverage copy has none |
| BEM triangles per surface | 5120 | 1078 | MNE `ico=4` |
| Conductivity: brain | 0.3 S/m | 1078 | — |
| Conductivity: skull | 0.06 S/m | 1078 | Probable typo for 0.006 (MNE v2.5 manual typo); run both (§3.4) |
| Conductivity: scalp | 0.3 S/m | 1078 | — |
| Forward/BEM software and version | Not stated | — | Likely MNE-C, around v2.5 [inference] |
| MEG system | Elekta-Neuromag Vectorview | 1078 | Sample uses the same system type |
| MEG channels | 204 planar gradiometers, 102 magnetometers (306) | 1078; 1082 | N per type not stated (A2) |
| EEG (out of scope) | 70 electrodes; average reference | 1078; 1082 | Sample has 60 |
| Sampling rate | 600 Hz | 1078 | Sample: 600.615 Hz |
| Hardware filters | 0.1–200 Hz (S1, S2); 0.03–200 Hz (S3, S4) | 1078 | Sample: 0.1–172.18 Hz |
| Spontaneous recording | Two minutes | 1078 | Sample raw is a ≈278-s task run (adaptation) |
| Noise-variance band | 0.5–100 Hz | 1079 | Filter design not stated |
| SSP | Magnetometers only; 3 largest eigenvectors of correlation matrix; 5-min empty-room | 1078 | Sample has 3 magnetometer PCA projectors; apply to gains too (A9) |
| Digitization | Polhemus FastTrack 3D: fiducials, EEG electrodes, HPI coils | 1078 | — |
| Head position | From HPI at the start of each measurement | 1078 | No values reported; positioning not studied |
| SNR | Eq. 1: 10·log10[(a²/N) Σ_k b_k²/s_k²] | 1079 | Channel-averaged power SNR in dB, with 1/N |
| N | Number of sensors | 1079 | Declare 102 / 204 / 306 explicitly |
| Difference map | Eq. 2: D = SNR_MEG − SNR_EEG | 1079 | EEG out of scope; an analogous mag-vs-grad contrast would be an extension |
| Combined map | max{SNR_MEG, SNR_EEG} | 1079 | — |
| Noise model | Identical, independent, uniform cortical sources, normal to mantle, Gaussian(0, s_s²) | 1079 | de Munck 1992 |
| Instrument/environmental noise | Not included | 1079–1080; 1082 | — |
| Model sensor variance | Eq. 3: s_k² = s_s²(AAᵀ)_kk | 1080 | Diagonal only |
| Noise-source spacing | ≈7 mm | 1080 | MNE `spacing=7`; report M (A3) |
| s_s calibration | Medians of ŝ²_k,recorded/(AAᵀ)_kk per type (grad, mag, EEG); channel-count-weighted mean | 1080 | Weights ≈0.543/0.271/0.186 [derived] |
| √s_s² | 1.6, 1.6, 1.6, 1.9 nAm | 1080 | Tied to the 7-mm density (A3, A4) |
| s_s without EEG | <5% change on average | 1081 | — |
| Group averaging | Spherical morphing (Fischl 1999b); shown on Subject 2, left hemisphere | 1079; 1080; 1082 | dB vs linear not stated (A7) |
| Area statistic | Fraction of sources with D > 0; ventricular surfaces excluded | 1079; 1081 | Report count- and area-weighted versions (A8) |
| Area results | ≈40% (noise model); 55% (recorded noise) | 1081 | Needs EEG; not reproducible with MEG alone |
| Fig. 2 colour limits | < −29 dB to > −19 dB | 1080 | Display clipping; sanity target for absolute level |
| Fig. 3 colour limits | < −4 dB to > 4 dB | 1082 | Display clipping |
| Fig. 4 coding | Number of subjects with D > 0 (0–4) | 1081–1082 | Counts cannot be reproduced on a single template |
| Patch-size effect | 10 dB (3 vs 8 cm², mesial temporal); EEG difference <1 to >10 dB | 1084 | Modality of the 10 dB unclear (A15) |
| Noise model for patches | Unchanged | 1084 | — |
| Trial averaging | Adds a constant; D unaffected | 1079 | +10·log10(n) dB [derived] |

---

## Sources used for external checks
- MNE software User's Guide v2.5, Dec 2006 (M. Hämäläinen): https://www.nitrc.org/docman/view.php/168/486/MNE-manual-2.5.pdf. Skull default p. 26; default source spacing p. 21; `--ico 4` p. 25; linear collocation p. 81; `mne_surf2bem` example §5.4.5.
- MNE User's Guide v2.7, Dec 2009: https://neuroimage.usc.edu/paperspdf/mne-manual-2.7.pdf. Skull default p. 27.
- MNE User's Guide v2.7.3, Nov 2010: https://mne.tools/mne-c-manual/MNE-manual-2.7.3.pdf. Skull default p. 27; default source surface (white) and Table 3.1 (sources per hemisphere, spacing, area per source), pp. 22–23.
- MNE-Python 1.13.2 source in the project `.venv` (`mne/bem.py`: `make_bem_model`, `make_bem_solution`; `mne.setup_source_space` and `mne.grow_labels` docstrings).
- MNE sample data in `data/external/MNE-sample-data/` (BEM conductivities, raw and empty-room metadata, source-space vertex counts), inspected read-only.
- Ahlfors SP et al. (2010), Sensitivity of MEG and EEG to source orientation, Brain Topogr. PMC2914866: https://pmc.ncbi.nlm.nih.gov/articles/PMC2914866/.
