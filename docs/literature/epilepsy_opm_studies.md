# OPM-MEG in epilepsy: a verified table of published studies

This note tabulates the published OPM-MEG studies in epilepsy: which brain regions they cover, which sensors, which ages, what they claim about the amplitude and SNR of epileptiform discharges, and against which comparator. Three modelling studies are listed separately as precedents. Every number was checked in the full text where we could read one, with its page, section or table; Section 9 lists what was read for each study. The same content, one record per study, is in `epilepsy_opm_studies.json`.

Searches and reading were done on 2026-10-04. Two checks were added on 2026-10-05 from `literature_checks.md` (§§1 and 5): access to the full text of Ren et al. (2025), and the frequency band behind each OPM noise figure.

**Conventions**

- *(p. N)* is the printed page, *(§x)* a section, *(T1)* a table, *(Fig. 2)* a figure legend.
- **[derived]**: we computed it from printed values; the paper does not state it.
- **[preprint]**: found only in the preprint, not in the version of record we read.
- **[abstract]**: taken from an abstract because the full text was not accessible.
- ± values are reproduced as printed, with the spread measure the paper states where it states one (Feys 2025 Table 2: mean ± SEM; Hillebrand 2023 Table 1: mean with SD in parentheses).
- Text is paraphrased. Quotation marks mark the paper's own wording.
- The JSON fields `age_group`, `sensor_family`, `comparator_class` and `region_tags` are our labels, assigned from the entries below.

**What was searched**

- PubMed, with these queries:
  - (optically pumped magnetometer(s) OR OPM-MEG OR on-scalp magnetoencephalography) AND (epilepsy OR epileptic OR interictal OR epileptiform OR seizure);
  - the same terms with spike, epileptogenic zone or irritative zone;
  - OPM AND children/infant AND epilepsy;
  - on-scalp/OPM AND epilepsy AND simulation/modelling;
  - author searches for Feys, De Tiège, Vivekananda and Mellor.
- Europe PMC, for preprints by the same groups.
- The reference lists of Feys 2025, Shen 2026 and Schwartz 2025.
- Web search, for preprints and publisher pages.

---

## 1. Overview: clinical studies

The table lists every study that recorded epileptiform activity with OPMs, in order of publication. The comparator column gives what each study's SNR or amplitude was compared with.

| Study | Patients (ages) | Epilepsy / regions | OPM sensors | Comparator used for SNR | Amplitude, OPM vs comparator | SNR, OPM vs comparator | Read |
|---|---|---|---|---|---|---|---|
| Vivekananda 2020, *Ann Clin Transl Neurol* | 1 (47 y) | focal, right posterior quadrant | Rb (QuSpin), 8 then 15 sensors | none (prior EEG) | not measured | not reported | full text |
| Feys 2022, *Radiology* | 5 children (5–11 y, median 9.4) | focal, extratemporal [preprint] | Rb Gen-2, 32 single-axis, flexible cap | 102 SQUID magnetometers (Triux), peak channel | 2.3–4.6× higher in 5/5 | 27–60 % higher in 4/5, no difference in 1 | published abstract + preprint full text |
| Feys 2023, *Ann Neurol* (ictal) | 1 child (10 y) | multifocal: occipital, parietal, frontal | Rb Gen-3, 24 triaxial, cap | none (rest vs hyperventilation) | – | seizure SNR 2.97 (rest) vs 2.79 (hyperventilation) | full text (letter) |
| Hillebrand 2023, *Sci Rep* | 7 (3 children 10–12 y, 4 adults) | temporal, temporo-parietal, sensorimotor | Rb Gen-2, 6 dual-axis (12 channels), 7–13 fT/√Hz | Triux Neo planar gradiometers | not reported | 0.92–0.98× [derived]; "comparable" | full text |
| Feys 2023, *Clin Neurophysiol* (infant) | 1 infant (5 mo; corrected age 2 mo) | right centro-parietal | Rb Gen-2/3, 26, cap with a 4-mm gap | none in the infant; other children's data | 2.5 pT, lower than the school-aged children's OPM values | 11.1; similar to or above school-aged children's cryogenic SNRs (authors) | full text (letter) |
| Badier 2023, *eNeuro* | 1 (33 y) | left temporal pole (neocortical); hippocampus-only spikes not seen | 4He prototype, 4 on scalp, < 50 fT/√Hz on 2 axes | 4 nearest magnetometers of a 248-channel 4D system | single spike 5 pT (text) or 2.5 pT (legend) vs 1.1 pT | single spike 7 vs 5 | full text |
| Feys 2023, *Front Neurosci* | 1 child (11 y) | right centrotemporal | 4 Rb Gen-3 and 4 4He (60–65 fT/√Hz), triaxial | Rb vs 4He (no SQUID) | 7.8 (Rb) vs 4.7 (4He) pT | 21.3 (Rb) vs 11.4 (4He) | full text |
| Feys 2024, *Ann Neurol* | 1 (20–30 y) | right amygdala | Rb Gen-2/3, 48 (116 channels), rigid helmet | none (simultaneous SEEG) | 1.00 pT | 2.75, vs 11.1–16.7 in earlier neocortical OPM data | full text (letter) + preprint |
| Feys 2025, *Epilepsia* | 11 (8–63 y, mean 37) | TLE: anterior, medial, basal, lateral, posterior temporal | Rb Gen-2/3, 43–53 (98–127 channels) incl. 4 face-OPMs, cap | 102 SQUID magnetometers, peak channel | higher in 2, similar in 1 of 3 tested; nominally higher in 9/9 | higher in 1, similar in 2 of 3 tested; nominally higher in 5/9 | full text |
| Ren 2025, *NeuroImage* | 46 (mean 23.7 y) | not given in the abstract | 128-channel whole-scalp (type not given) | SQUID (not specified) | higher (p < 0.001) | higher (p = 0.003) | abstract only |
| Schwartz 2025, *Epilepsia Open* | 7 + 1 (22–48 y; 41 y) | 5 focal temporal, 2 generalized; deep mesial (SEEG case) | 4He, 4 on scalp + 1 reference, < 45 fT/√Hz on 2 axes | 4 nearest channels of a 275-channel CTF (4D for the SEEG case) | 2.3× (group), 1.4–5.3× (individuals) | 6.72 vs 8.37 (0.80× [derived]); deep spike 6.8 vs 9.0 | full text |
| Shen 2026, *Epilepsia* | 68 (6–60 y, mean 28) | 51 temporal, 17 extratemporal | 64 dual-mode (128 channels), < 15 fT/√Hz, rigid spherical helmet | none (concordance with iEEG) | – | – (concordance 80.1 % temporal vs 92.0 % extratemporal) | full text |

**Noise figures and their bands** (checked 2026-10-05; `literature_checks.md` §5). No paper in the table states the frequency band of the noise figure it gives for its own sensors. Bands appear only in the device sources behind some figures, and they start at 1–10 Hz:

- Hillebrand 2023 cites its 7–13 fT/√Hz to Osborne et al. (2018), whose first-generation QuSpin units are typically 10 fT/√Hz over 1–100 Hz. QuSpin's Gen-2 specification is below 15 fT/√Hz over 3–100 Hz.
- Feys 2023 (*Front Neurosci*) gives bands only for its Introduction's general figures: below 23 fT/√Hz over 3–100 Hz for current Rb-OPMs (QuSpin's Gen-3 triaxial specification) and below 50 fT/√Hz over 1–1,500 Hz for 4He. The figures for its own sensors (15 and 60–65 fT/√Hz) have no band.
- Vivekananda 2020's white-noise floor of about 10 fT/√Hz cites Boto et al. (2018), who give about 15 fT/√Hz, also without a band.
- No band was found for the figures of Shen 2026 (< 15 fT/√Hz), Badier 2023 and Schwartz 2025 (4He).
- For comparison, FieldLine's HEDscan specification sheet gives 8 and 15 fT/√Hz over 10–130 Hz.

## 2. Overview: modelling precedents

| Study | Heads / data | OPM model | Comparator | Noise in the SNR | Main result |
|---|---|---|---|---|---|
| Zahran 2022, *Sensors* | MRIs of 12 healthy subjects: adult, 9 y, ten infants (1 mo–2 y) | 4He, 102 sensors × 1–3 axes, 40 fT/√Hz, 3 mm from scalp | 102 SQUID magnetometers, 5 fT/√Hz, 2 cm; fitted or standard helmet | sensor noise only | adult topography power 8.9× (three axes); in the adult the SQUID array carries 1.133× the information of the normal-axis OPM array; the OPM advantage grows as the head shrinks |
| Westin 2023, *Clin Neurophysiol* | 1 adult template; 10 epileptic networks | 128 OPMs, 30 fT/√Hz, 1-mm standoff | 306-channel VectorView; 60-channel EEG | brain noise (resting data) plus sensor noise | spikes visible from a mean 0.56 cm² on-scalp vs 1.91 cm² (MEG) and 1.87 cm² (EEG); mesial temporal is the only site where the smallest seizure-onset patch was missed on-scalp |
| Matsubara 2026, *NeuroImage* | 54 patients with epilepsy (4.7–61.4 y) | 102 OPMs at SQUID positions (5-mm offset); 70 at EEG positions; 20 over the posterior fossa | 306-channel VectorView; 70-channel EEG | brain noise only | moving sensors closer at SQUID-matched positions raised SNR in no region (including medial temporal, cingulate, insula); a targeted posterior-fossa layout raised SNR in the posterior cerebellum |

---

## 3. Verified entries: clinical studies

### 3.1 Vivekananda et al. 2020, *Ann Clin Transl Neurol* 7(3):397–401 (doi 10.1002/acn3.50995)

- **Design.** One patient. Two separate 30-min OP-MEG sessions, with no SQUID recording. The OP-MEG was compared qualitatively with an earlier, non-simultaneous 5-day EEG telemetry (Method).
- **Patient.** A 47-year-old woman with refractory focal epilepsy after Haemophilus meningitis with status epilepticus at 18 months. The focus was in the right posterior quadrant, with right parieto-occipital damage on MRI. EEG interictal activity was maximal at T8 and the right sphenoidal electrode (Method, Results).
- **OPMs.** Rubidium sensors; the conflict-of-interest statement names QuSpin as the manufacturer.
  - Session 1: 8 first-generation OPMs in a plinth, patient lying, plus 6 reference sensors.
  - Session 2: 15 second-generation OPMs in a 3D-printed scanner-cast (Method).
  - No noise was measured. The Introduction cites a white-noise floor of about 10 fT/√Hz for recent OPMs, and a sensitive volume about 6 mm from the scalp against about 3–4 cm for cryogenic MEG.
- **Comparator, SNR, amplitude.** None measured. The Discussion *expects*, from earlier studies, a two- to fivefold smaller signal with conventional MEG.
- **Localization.** A dipole fit of the averaged spikes from session 2 fell in the right posterior quadrant. A later resection reduced seizure frequency (Results).
- **Other numbers.**
  - 23 spikes per 30 min in session 1 and 19 per 30 min in session 2.
  - About 30 discharges per hour on EEG telemetry.
  - No seizure occurred during OP-MEG (Results).
- **Read:** full text (PMC7085997).

### 3.2 Feys et al. 2022, *Radiology* 304(2):429–434 (doi 10.1148/radiol.212453)

- **Design.** Prospective single-centre series. Each child had both OPM-MEG and cryogenic MEG; the two were not simultaneous, and OPM-MEG came first in every child [preprint, Methods]. The study compared IED detection, amplitude, SNR and distributed source localization, with t tests [abstract] run within each child [preprint].
- **Patients.** Five children, median age 9.4 years (range 5–11), four of them girls. Three had self-limited idiopathic and two refractory focal epilepsy [abstract].
  - The preprint calls all five extratemporal, with unifocal IEDs.
  - One refractory child had a resected right temporal dysembryoplastic neuroepithelial tumour, with the presumed epileptogenic zone elsewhere [preprint].
  - The per-child IED regions are in Table 1, which we could not see.
- **OPMs** [preprint, Data acquisition]:
  - 32 QuSpin Gen-2 OPMs in single-axis mode on a flexible EEG-like cap (EasyCap), with 3D-printed mounts sewn at 10-10 positions around the presumed focus.
  - A compact OPM room with a background below 15 nT after degaussing and no further field compensation.
  - Sensor noise was not reported.
- **Comparator.** A 306-channel MEGIN Triux, restricted to its 102 magnetometers "for comparability" [preprint, Data preprocessing]. Feys 2025 restates this restriction (p. e148).
- **SNR definition.** Peak amplitude and SNR were measured for each spike at the sensor with the maximal averaged spike amplitude, with baseline correction from −100 to −50 ms [preprint, Data analysis]. Feys 2025 (p. e145) states that this SNR is the peak amplitude over the baseline SD.
- **Amplitude.** IED amplitude was 2.3 (7.2 of 3.1) to 4.6 (3.2 of 0.7) times higher with OPMs, P < .001, in all five children [abstract]. The abstract prints no units. The preprint gives the range as 2.3–4.8×.
- **SNR.** SNR was 27 % (16.7 of 13.2) to 60 % (12.8 of 8.0) higher with OPMs (P = .001–.009) in four children. The fifth child showed no difference (P = .93); head movement caused motion artefacts in that recording [abstract].
- **Localization.** The sources of the averaged IEDs were about 5 mm apart in three children, and 8.3 mm and 15.6 mm apart in the other two [abstract]. The preprint gives the range as 4.2–15.6 mm.
- **Other numbers** [preprint]:
  - spike-wave index 2–89 %;
  - mean distance from the source to the closest sensor 29.4 mm (OPM) vs 57.6 mm (SQUID).
- **Read:** the published version as an abstract only (PubMed 35503013, publisher abstract); the medRxiv preprint v1 (10.1101/2021.09.06.21262839) in full.
- **Notes.**
  - Where the two versions differ, we use the published value (upper amplitude ratio 4.6, not 4.8).
  - The child in Feys 2023 (*Front Neurosci*) is patient 5 of this series.
  - The infant letter (Feys 2023, *Clin Neurophysiol*) quotes this study's mean OPM IED amplitudes as 3.2–9.9 pT, from its Table 2.

### 3.3 Feys et al. 2023, *Ann Neurol* 93(2):419–421: ictal recording (doi 10.1002/ana.26562)

- **Design.** One child. 1 h of resting video-OPM-MEG plus 3 min of hyperventilation, followed later by SEEG. Cryogenic MEG was unavailable because of technical issues.
- **Patient.** A 10-year-old boy with drug-resistant focal epilepsy and more than 5 seizures per day.
  - More than 600 IEDs of two types: about 70 % left mesio-occipital (source in the left cuneus) and about 30 % left temporo-occipital (left supramarginal gyrus).
  - 15 spontaneous and 6 hyperventilation-induced seizures: left frontal 3, temporal 1, parietal 9 and occipital 5; right parietal 2. The right occipital count is not legible in the text we obtained.
  - SEEG with 16 posterior electrodes confirmed multifocal onset zones.
- **OPMs.**
  - 24 triaxial QuSpin QZFM-G3 OPMs on an EEG-like cap, sampled at 1200 Hz.
  - A Compact MuRoom (Cerca), with a remnant field below 1–2 nT after degaussing and active nulling.
  - Noise was not reported.
- **SNR definition.** Seizure SNR is the mean amplitude of the whole seizure divided by the mean amplitude of the preceding IED-free background, at the sensor with the maximal ictal amplitude.
- **Results.**
  - Seizure SNR was similar in the two conditions: 2.79 ± 0.07 with hyperventilation and 2.97 ± 0.19 at rest (p = 0.43).
  - Hyperventilation raised seizure amplitude (1.92 ± 0.04 vs 1.49 ± 0.05 pT) and background amplitude (0.69 ± 0.02 vs 0.52 ± 0.03 pT), p < 0.003.
  - The posterior onset and irritative zones from OPM-MEG co-localized with SEEG.
  - 10 seizures were accompanied by eye, eyelid or head movements.
- **Read:** full text of the letter (publisher text via the DOI; the letter has no internal sections).
- **Note.** Shen 2026 cites this letter (pp. 3938, 3949) for OPM SNR that is comparable to or better than SQUID-MEG. The letter contains no SQUID comparison.

### 3.4 Hillebrand et al. 2023, *Sci Rep* 13:4623 (doi 10.1038/s41598-023-31111-y)

- **Design.** Case series. OPM and SQUID were recorded on the same day but not together. OPM came in the morning (seated, 15-min runs) and SQUID in the afternoon (supine, four 15-min runs). OPM positions followed the IED field maps of earlier clinical MEG, and also SEEG in the adults (Methods).
- **Patients.** Seven, all with focal drug-resistant epilepsy, sleep-deprived for the recordings.
  - Groups: three adults with deep or weak IED sources, three children aged 10–12 years, and one adult with daily seizures. The adults' ages are in Supplementary Table S1, which we did not consult.
  - The hypothesised epileptogenic zone was bilateral mesial temporal, multifocal, or near somatosensory cortex (Methods).
  - IEDs recorded in this study: bilateral temporal (patient 1), right superior temporal/parietal (3), right sensorimotor (4), right central/superior temporal (6); none in patients 2 and 5.
  - Patient 7 had a seizure on the OPMs with onset over both anterior temporal lobes (Results).
- **OPMs.**
  - Six QuSpin Gen-2 OPMs in dual-axis mode (radial plus one tangential axis), giving 12 channels.
  - Sensitivity 7–13 fT/√Hz, range ±5 nT, bandwidth 0 to about 130 Hz.
  - Individual rigid 3D-printed helmets, bi-planar coils with dynamic field compensation, and homogeneous field correction (Methods).
- **Comparator.** A 306-channel MEGIN Triux Neo in the same shielded room, processed with tSSS. The SNR was computed on its planar gradiometers only (Table 1 caption).
- **SNR definition.** The Z-score of each IED from a simple automatic detector at sensor level, averaged over IEDs and reported as mean (SD), used "as a proxy for SNR". The paper also reports a spike-wave index: the percentage of seconds that contain an IED.
- **Amplitude.** Not reported.
- **SNR (T1)**, OPM vs SQUID:

  | Patient | OPM SNR, mean (SD) | SQUID SNR, mean (SD) | OPM/SQUID [derived] |
  |---|---|---|---|
  | 1 | 4.47 (0.43) | 4.57 (0.46) | 0.98 |
  | 3 | 3.85 (0.32) | 3.93 (0.35) | 0.98 |
  | 4 | 4.05 (0.48) | 4.39 (0.61) | 0.92 |
  | 6 | 6.15 (0.98) | no IEDs | – |

  The authors call the SNRs comparable.
- **Spike-wave index (T1)**, OPM vs SQUID:
  - patient 1: 6.42 vs 10.25 (6.11 and 6.58 using the left or right SQUID channels only);
  - patient 3: 9.00 vs 6.76;
  - patient 4: 11.03 vs 24.50;
  - patient 6: 0.53 vs 0. The OPMs recorded 3 and 5 IEDs in two runs; the four SQUID runs recorded none.
- **Localization.** Field maps were consistent between the two systems. Beamformer virtual electrodes were computed, but no localization-distance metric is given.
- **Read:** full text (PMC10030968), with Table 1 transcribed from the Europe PMC XML.
- **Notes.**
  - Sensor balance was 6 sensors (12 channels) vs 102 sensors (306 channels).
  - The authors list factors that cloud the comparison: dual-axis operation costs some sensitivity and the placement was not optimal (both lower OPM SNR); time of day, posture and drowsiness differed between sessions; and movement artefacts may hide IEDs in the OPM data.
  - The preprint (10.1101/2022.11.03.22281836) adds effect sizes of d = −0.14, −0.14 and −0.36, and its Table 1 caption says all patients were male [preprint].

### 3.5 Feys et al. 2023, *Clin Neurophysiol* 155:29–31: infancy (doi 10.1016/j.clinph.2023.08.010)

- **Design.** One infant, 30 min of on-scalp video-OPM-MEG without sedation. There was no cryogenic MEG of this infant.
- **Patient.**
  - A 5-month-old girl born preterm (corrected age 2 months), with microcephaly (head circumference 33 cm).
  - Grade IV intraventricular haemorrhage, treated with a ventriculo-peritoneal shunt.
  - Focal epilepsy involving C4 on EEG.
  - IEDs over right parietal sensors; dipoles clustered in the right centro-parietal region.
- **OPMs.**
  - 26 QuSpin OPMs (17 QZFM-G3, 9 QZFM-G2) at 1200 Hz, on a 34-cm EEG-like cap.
  - A shielded room with field-nulling coils and a background below 2 nT.
  - A 4-mm gap between the sensors and the bald scalp, to let heat dissipate. The Fig. 1 legend says 3 mm.
  - Radial axes were used for dipole fitting. Noise was not reported.
- **Comparator.** None in the same patient. The authors compare this infant with the five school-aged children of Feys 2022.
- **SNR definition.** Not defined in the letter, which refers to Feys 2022 for methods.
- **Amplitude.** IED amplitude 2.5 ± 1.1 pT and background 0.30 ± 0.15 pT. Both are lower than in the school-aged children: their mean IED amplitudes were 3.2–9.9 pT, and their background 0.32–0.97 pT (cited as unpublished data).
- **SNR.** 11.1 ± 8.9, after removal of cardiac artefacts.
  - Against the school-aged children's OPM SNRs: close to 2 of them and below the other 3.
  - Against their cryogenic SNRs: "similar or higher".
- **Localization.** Equivalent current dipoles with a spherical model and a goodness of fit above 80 % clustered in the right centro-parietal region, consistent with C4 involvement.
- **Other numbers.** 70 monomorphic IEDs, 41 of them with a goodness of fit above 80 %. Cap placement took under 1 min.
- **Read:** full text of the letter (publisher HTML via ScienceDirect).
- **Note.** The authors attribute the lower amplitude to the extra gap and the large lesion.

### 3.6 Badier et al. 2023, *eNeuro* 10(12):ENEURO.0222-23.2023 (doi 10.1523/ENEURO.0222-23.2023)

- **Design.** One patient, two supine sessions, each simultaneous with SEEG. The 20-min SQUID-MEG/SEEG session came first, then a 4He-OPM/SEEG session. Equivalent IEDs were matched across the two sessions by their SEEG pattern, then compared singly and as averages by type (Materials and Methods).
- **Patient.** A 33-year-old woman with intractable temporal lobe epilepsy and a cavernoma of the left parahippocampal gyrus.
  - SEEG: 14 electrodes with 119 contacts, 12 of the electrodes in the left hemisphere.
  - Epileptogenic network: the left mesial temporal structures and the temporal pole.
  - The spikes seen by MEG came from the left temporal pole and the adjacent anterior third temporal gyrus. Spikes confined to the hippocampus were seen by neither system (Results; Extended Data Figs. 1-1, 2-1; Discussion).
- **OPMs.**
  - A prototype with five 4He sensors; four sensors (2 × 2 × 5 cm) sat on a rigid helmet over the left central and temporal regions.
  - Closed-loop triaxial operation. Noise better than 50 fT/√Hz on the radial and one tangential axis, and 200 fT/√Hz on the third axis, which was not used.
  - Range ±250 nT, bandwidth DC–2 kHz.
  - A two-layer shielded room; the only noise reduction was filtering (2–70 Hz, 50-Hz notch).
- **Comparator.** The four SQUID sensors closest to the four OPMs, 2.72–3.23 cm from them. They belong to a 4D Neuroimaging 3600 system with 248 magnetometers and reference-based noise compensation.
- **SNR definition.** The maximum amplitude of the event over the SD of a baseline: 2 s before the event for single spikes, 500 ms for averages.
- **Amplitude.**
  - Single spike: OPM peak-to-peak 5 pT in the text, or 2.5 pT in the Fig. 2 legend, against 1.1 pT for the SQUID.
  - Averages: −9.5 pT (spike I) and −10 pT (spike II) on the OPMs, against small SQUID deflections of −6 and 3.7 pT maximum.
  - The figure legends give peak-to-peak OPM vs SQUID values of 15 vs 8 pT (type I, n = 10) and 9.3 vs 4 pT (type II, n = 3).
- **SNR.**
  - Single spike: 7 (OPM) vs 5 (SQUID).
  - Averages: 6.4 (spike I) and 16.7 (spike II) on the OPM radial axis, printed as "SNR, *r*". The SQUID SNRs of the averages are not printed.
- **Read:** full text (PMC10748329), with the averaged-SNR sentence and legends re-checked in the Europe PMC XML.
- **Notes.**
  - The single spikes compared came from different sessions. The OPM-session spike was the stronger one on SEEG (SNR 24.8 vs 14.9, Fig. 2), which the authors acknowledge.
  - The text and the Fig. 2 legend disagree on the OPM amplitude (5 vs 2.5 pT).
  - The Fig. 1 legend gives 265 SQUID sensors; the Methods give 248 magnetometers.

### 3.7 Feys et al. 2023, *Front Neurosci* 17:1284262: Rb vs 4He case (doi 10.3389/fnins.2023.1284262)

- **Design.** One child, recorded for 40 min with tri-axial Rb-OPMs and 4He-OPMs at the same time. Four EEG electrodes were also recorded but could not be analysed. The statistics use 102 simultaneous artefact-free IEDs and, for each sensor type, the one sensor with the highest IED amplitude (Methods).
- **Patient.** An 11-year-old girl with refractory focal epilepsy that began at 7 years; she is patient 5 of Feys 2022.
  - She had a right anterior temporal lobectomy for a dysembryoplastic neuroepithelial tumour and has been seizure-free since (Engel 1A).
  - IEDs remain frequent, maximal at C4-T4.
  - Cryogenic MEG 7 months earlier showed right centrotemporal IEDs.
- **OPMs.**
  - 4 QuSpin Gen-3 tri-axial Rb-OPMs and 5 4He prototypes on a helmet with 89 mounts.
  - One right-hemisphere 4He sensor failed; the other four ran at 60–65 fT/√Hz.
  - The Discussion gives intrinsic noise of 15 fT/√Hz for the Rb sensors and 65 fT/√Hz for the 4He sensors.
  - Sensors of both types sat about 1 cm from C4 or T4, about 2 cm from each other.
- **Comparator.** 4He vs Rb only; there is no SQUID comparison.
- **SNR definition.** The peak amplitude of each IED over the background amplitude 100–50 ms before the peak. It is measured on a virtual sensor oriented along the principal field direction, in the single ICA component that contains the IEDs.
- **Amplitude.**
  - Virtual sensors: 7.8 ± 0.2 pT (Rb) vs 4.7 ± 0.1 pT (4He).
  - Rb axes: first tangential 7.7 pT, radial 3.6 pT, second tangential 1.0 pT.
  - 4He axes: radial 4.6 pT, tangential 1.5 and 2.0 pT (Results).
- **SNR.** 21.3 ± 1.4 (Rb) vs 11.4 ± 0.8 (4He), Cohen's |d| = 0.8. The authors attribute the difference to sensor position, not sensor performance (Results, Discussion).
- **Other numbers.**
  - IEDs counted by three readers: 1,372, 1,287 and 1,271 (Rb) vs 1,175, 1,231 and 1,221 (4He).
  - Background in the IED component: 0.5 pT for both.
  - Background in the IED-free components: 0.8 pT (Rb) vs 1.9 pT (4He).
- **Read:** full text (PMC10715393), with the patient's identity re-checked in the XML.

### 3.8 Feys et al. 2024, *Ann Neurol* 95(3):620–622: mesial temporal (doi 10.1002/ana.26844)

- **Design.** One patient. 10 min of video-OPM-MEG recorded simultaneously with SEEG, at the end of a two-week video-SEEG. OPM IEDs were marked blind to SEEG.
- **Patient.** A woman aged 20–30 years with refractory focal epilepsy [preprint].
  - She had a large right temporo-parietal cortical malformation and multiple right nodular heterotopias [preprint].
  - The most active irritative zone was the right amygdala; there were several seizure-onset zones.
- **OPMs.**
  - 48 QuSpin OPMs (20 QZFM-G3, 28 QZFM-G2), giving 116 channels; the preprint's abstract and Methods say 118.
  - Mounted on an adult-sized rigid 3D-printed helmet, so close to the scalp rather than on it.
  - 18 of the 116 channels were excluded as noisy.
- **Comparator.** No SQUID recording. SEEG served as the reference (64 of 186 contacts recorded), and the SNR was compared only with published values.
- **SNR definition.** Not defined in the letter.
- **Results.**
  - 31 OPM IEDs: amplitude 1.00 ± 0.76 pT, SNR 2.75 ± 2.55.
  - The 31 were among 52 IEDs on SEEG (60 %). The SEEG IEDs that OPMs detected were larger: 1918 ± 107 vs 922 ± 85 µV/cm (p = 8.97 × 10⁻¹⁰).
  - The OPM source lay 14.1 mm from the SEEG contacts with the largest IEDs: right amygdala on SEEG, right periamygdalar white matter on OPM.
  - The authors compare the SNR of 2.75 with 11.1–16.7 in earlier OPM recordings of neocortical IEDs (Feys 2022 and the 2023 infant letter). They call it similar to simultaneous cryogenic MEG–SEEG recordings of neocortical IEDs, which the preprint puts at about 2 [preprint].
  - A 60 % detection rate is at the upper end of the 25–60 % reported for mesiotemporal IEDs in cryogenic MEG–SEEG studies.
- **Read:** full text of the letter (accepted-manuscript text via the DOI), cross-checked against the full text of the medRxiv preprint v1 (10.1101/2023.10.03.23296442).
- **Note.** The authors attribute the low SNR to three factors:
  - the depth of the source;
  - the close-to-scalp helmet, used because the implants ruled out the cap;
  - nearby screws, electrodes, cables and the EEG amplifier.

### 3.9 Feys et al. 2025, *Epilepsia* 66(7):e142–e151: temporal lobe epilepsy (doi 10.1111/epi.18439)

- **Design.** Prospective single-centre series.
  - Every patient had 1 h of on-scalp OPM-MEG with both scalp-OPMs and face-OPMs.
  - Nine patients also had cryogenic MEG: mostly on the same day, but 4 and 15 months apart for two of them.
  - OPM and cryogenic values were compared statistically only when both recordings had at least 10 IEDs, which applied to 3 patients (§2.1–2.3, pp. e143–e148; T1, p. e145).
- **Patients.** Eleven with TLE, aged 8–63 years (mean 37), 4 female and 7 male. Inclusion allowed adults and children older than 6 years (p. e143; T1).
  - OPM IED sources were all temporal: anterior 2, medial 2, basal 3, lateral 2, posterior 2 (p. e148).
  - Table 2 labels include temporal pole (patients 1 and 9), peri-amygdala and peri-hippocampal (8), perihippocampal (10), T1–T4 gyri and sulci, and the temporal opercula (p. e146).
- **OPMs** (§2.2, p. e143; Fig. 1 legend, p. e145; §2.3, p. e145; p. e148):
  - 43–53 QuSpin Gen-2 and Gen-3 biaxial or triaxial OPMs per patient, giving 98–127 usable channels.
  - The montage included four face-OPMs on a glasses-like structure.
  - Flexible EEG-like cap (EasyCap); Compact MuRoom with a remnant field below 1–2 nT; sampling at 1200 Hz.
  - PCA denoising that removed the first 15 components. Noise was not reported.
- **Comparator.** Cryogenic MEG at the same centre, limited to its 102 magnetometers "for comparability"; the system is not named in this paper (p. e148). The comparison set 102 cryogenic magnetometers against 102–116 OPM channels (p = .31; p. e149).
- **SNR definition** (pp. e145, e148). The peak IED amplitude over the baseline SD, after averaging IEDs of similar morphology and topography, at the channel with the maximal IED response. Table 2 reports mean ± SEM.
- **Amplitude and SNR** (p. e148; T2, p. e146).
  - OPM means across patients: amplitude 3.3 pT (patient means 1.4–7.2 pT), SNR 9.4 (3.4–21.5).
  - Tested patients:

  | Patient (IED source) | Amplitude, OPM vs cryogenic (pT) | SNR, OPM vs cryogenic | Authors' reading |
  |---|---|---|---|
  | 3 (left posterior T4 gyrus) | 4.0 ± .4 vs 2.2 ± .3 (p = 6.0 × 10⁻⁴) | 11.5 ± 2.0 vs 5.4 ± .9 (p = 7.0 × 10⁻³) | higher amplitude and SNR |
  | 9 (right temporal pole) | 2.4 ± .2 vs 1.0 ± .03 (p = 6.7 × 10⁻⁷) | 9.6 ± 1.2 vs 13.5 ± 3.6 (p = .32) | higher amplitude, similar SNR |
  | 10 (right perihippocampal) | 1.4 ± 1.7 vs 1.2 ± 1.7 (p = .31) | 3.6 ± .6 vs 4.4 ± .5 (p = .26) | similar |

  - The nine patients with both recordings, descriptively [derived from T2]:
    - Amplitude was higher with OPMs in all nine; the OPM/cryogenic ratio was 1.04–9.0 (median 2.38).
    - SNR was nominally higher with OPMs in 5 (patients 2, 3, 4, 7, 11) and lower in 4 (5, 8, 9, 10); the ratio was 0.62–4.30 (median 1.24).
    - Only patients 3, 9 and 10 were tested.
- **TLE vs extratemporal epilepsy.** The comparison group was 7 earlier extratemporal (ETLE) patients recorded with the same cap.
  - Amplitude: 3.3 vs 6.0 pT (ETLE range 2.5–9.9, p = .05).
  - SNR: 9.4 vs 14.3 (ETLE range 11.1–21.3, p = .02; p. e148).
  - Distance between cryogenic and OPM sources: 6–16 mm (TLE) vs 4–16 mm (ETLE), p = .23 (p. e149).
- **Localization** (p. e148).
  - On-scalp and cryogenic sources were 6–16 mm apart.
  - SEEG (patient 6) and the resection cavity of a seizure-free patient (11) agreed with the OPM sources.
  - Face-OPMs changed the localization substantially only in patient 1; in patients 2–11 the shift averaged 3.5 mm (range 0–11 mm).
  - Face-OPMs were farther from the IED sources than scalp-OPMs (p = 3.26 × 10⁻⁷).
- **Read:** full text (local PDF of the version of record). Supplementary Table S1 was not consulted.
- **Notes.**
  - The authors attribute the lower TLE SNR to the deep location and partial coverage of the irritative zone in patients 1 and 8–10 (p. e149).
  - The Discussion reads the three tested patients as showing non-inferiority and, in some cases, superiority of on-scalp MEG (p. e149).
  - Values for patients with fewer than 10 IEDs are descriptive only.
  - Recording conditions differed: sitting vs supine (T1 note), and two patients had delays of 4 and 15 months.

### 3.10 Ren et al. 2025, *NeuroImage* 312:121232 (doi 10.1016/j.neuroimage.2025.121232): abstract only

- **Design.** Each patient had both SQUID-MEG and 128-channel whole-scalp OPM-MEG; the order and timing are not given in the abstract.
- **Patients.** 46, mean age 23.7 ± 8.7 years, 29 male; 39 had IEDs in both systems.
- **OPMs and comparator.** The OPM sensor type and noise, the SQUID system, the channels used for SNR and the SNR definition are not given in the abstract.
- **Results** [abstract]:
  - IED detection did not differ between systems (McNemar). OPM detection accuracy was 91.3 % relative to SQUID-MEG; Gwet AC1 was 0.892.
  - Among the 39 patients with IEDs in both systems, OPM-MEG had a closer sensor–scalp distance (p < 0.001), higher IED amplitude (p < 0.001) and higher SNR (p = 0.003).
  - Source localization was "nearly consistent" at the sublobar level. In 24 patients with single dipole clusters, the cluster centroids were 12.16 ± 5.90 mm apart.
- **Read:** abstract only (PubMed 40254146; publisher abstract and highlights). The article is open access (CC BY 4.0 in the publisher's Crossref record), but its full text could not be read with the tools available: the publisher's site refused automated access, and the article is not in PMC or Europe PMC (checked 2026-10-05; `literature_checks.md` §1).
- **Note.** Effect sizes, the SNR definition, the comparator channels and any regional breakdown cannot be checked. PubMed's MeSH indexing lists Adolescent, Young Adult and Adult.
- **Indirect clues** (not the article's text; `literature_checks.md` §1). They point, indirectly, to a MEGIN TRIUX and the 128-channel X-Magtech system, so the systems are at most probable; the entries above do not use them.
  - PubMed's competing-interest statement: one author (M. Ding) is a board member of Beijing X-Magtech Technology Limited, which makes the Marvel MEG OPM system (Shen 2026 used its 128-channel version). It does not say which system the study used.
  - A 2024 American Epilepsy Society meeting abstract by the same group (abstract 2.489) names a Neuromag TRIUX and a Marvel MEG. Its 21 patients (median age 35 years) cannot be assumed to be a subset of this cohort.
  - A Beijing municipal news item (24 July 2024) on the market approval of the 128-channel Marvel MEG places its registration clinical study at Beijing Tiantan Hospital, the authors' institution. It does not cite the article.

### 3.11 Schwartz et al. 2025, *Epilepsia Open* 10:1660–1672 (doi 10.1002/epi4.70139)

- **Design.** Two parts.
  - Part 1: seven patients with sequential resting recordings, seated, in a standard two-layer µ-metal room with one copper layer and no active shielding (§2.2.1, pp. 1662, 1664). SQUID (30 min) always came first, then 4He-OPM (30 min). The four 4He sensors were placed where SQUID showed the most IEDs.
  - Part 2: one SEEG patient with two supine sessions, each simultaneous with SEEG (§2.2.2, p. 1665).
- **Patients.** Part 1: 22–48 years, mean 29.2 ± 8.9 (§2.1, p. 1662). Part 2: a 41-year-old woman.
  - Part 1 had 5 focal epilepsies, all with temporal interictal spikes on EEG: left temporal in 3, one of them with a left medial temporal focal dysplasia, and right temporal in 2. The other 2 were generalized (idiopathic generalized and Lennox–Gastaut) (T1, p. 1663).
  - IEDs were recorded in patients 1, 3, 5, 6 and 7 (T2, p. 1665).
  - The part-2 patient had non-lesional TLE with bitemporal spikes and right medial temporal seizures. Her averaged spike involved the amygdala, anterior hippocampus and temporo-basal region on SEEG (p. 1668).
- **OPMs** (§2.2.1, pp. 1662, 1664; §2.3.1, p. 1665):
  - 5 triaxial 4He-OPMs (MAG4Health): 4 on the scalp in a semi-rigid helmet with 96 positions, and 1 as a reference 10 cm above the head.
  - 570 g in total; sampling at 11 kHz.
  - Noise below 45 fT/√Hz on the radial and first tangential axes, and about 200 fT/√Hz on the second tangential axis, which was used only for denoising.
  - Processing: 6–70 Hz band-pass, 50- and 60-Hz notches, regression on the reference sensor, then ICA.
- **Comparator** (§2.2.1, p. 1662; §2.3.1–2.3.2, pp. 1665–1666).
  - Part 1: a 275-channel CTF system (2.4 kHz, third-order gradient correction, 3–70 Hz). IEDs were marked on all 274 working channels, but SNR and amplitude were compared on the 4 channels closest to the 4 OPMs. The paper does not name the CTF sensor type.
  - Part 2: a 4D Neuroimaging 3600 system with 248 magnetometers, again using its 4 channels closest to the OPMs.
- **SNR definition.** The maximum spike amplitude over the SD of a 20-s IED-free baseline, per IED (p. 1665). For the averaged spikes in part 2: the event maximum over the baseline SD (p. 1666).
- **Results, part 1** (§3.1, pp. 1666–1667; T2):
  - IEDs appeared with both systems in the same 5 of 7 patients: 975 with SQUID (mean 195 ± 137) and 417 with OPM (mean 83 ± 97). Per patient, SQUID/OPM: 79/96, 276/247, 69/25, 161/24, 390/25.
  - Group SNR: SQUID 8.37 ± 6.54 vs OPM 6.72 ± 4.86 (p < 0.01, Cohen's d = 0.27); OPM/SQUID 0.80 [derived].
  - Group amplitude: OPM 3.76 ± 2.81 pT vs SQUID 1.66 ± 0.66 pT (p < 0.01, d = 1.28; ratio 2.3).
  - Individuals: SQUID SNR was higher in 3 patients (d = 1.21, 0.26, 0.15), and OPM SNR higher but not significantly in 2 (p > 0.1). OPM amplitude was higher in every patient, by a ratio of 1.4–5.3 (d = 1.5–6.5).
  - OPM radial vs tangential axis: SNR 4.42 ± 2.68 vs 4.12 ± 2.36 (d = 0.12); amplitude 2.28 ± 1.85 vs 1.48 ± 1.17 pT (ratio 1.5).
- **Results, part 2.** For the averaged deep spike (23 SQUID epochs, 20 OPM epochs), SNR was 9.0 (SQUID) vs 6.8 (OPM) (§3.2, p. 1668; Fig. 4).
- **Localization.** Not performed, because there were too few sensors.
- **Read:** full text (local PDF of the version of record).
- **Notes.**
  - The paper does not say how sensors and axes were pooled for the OPM–SQUID comparison. The pooled OPM SNR (6.72) and amplitude (3.76 pT) are higher than either single-axis value.
  - The "~30 fT/√Hz noise floor" in the Discussion (p. 1670) refers to a newer whole-head system (Bonnet et al. 2025), not to these sensors.
  - SQUID always came first, while patients could still sleep (p. 1669), so order may have favoured it.
  - The Discussion (p. 1670) describes the OPM/SQUID signal-amplitude ratio as a "power ratio" of 2 to 5.
  - Two typos: tSSS is said to be applied in "patient 11", but there were only 7 patients (p. 1665); and Fig. 3 is cited for "Patient 1 and 20", but it shows patients 1 and 5 (p. 1667).

### 3.12 Shen et al. 2026, *Epilepsia* 67:3937–3954 (doi 10.1002/epi.70273)

- **Design.** Prospective diagnostic study at Tongji Hospital, January 2024 to May 2025 (§2.2, pp. 3939–3940).
  - Each patient had 90 min of interictal OPM-MEG, analysed by dipole fitting without spike averaging (§2.3–2.4, pp. 3940–3942).
  - The OPM result was compared with the iEEG-defined epileptogenic zone (31 SEEG, 37 ECoG) at the sublobar level (§2.7–2.11, pp. 3943–3944).
  - Diagnostic accuracy was measured against outcome at 12 months or more in 51 patients.
  - There was no SQUID arm.
- **Patients.** 68 of 80 screened, mean age 28 years (range 6–60), 26 female (p. 3944).
  - 15 were minors: 3 under 10 years and 12 aged 10–18 (T1, p. 3945; p. 3939).
  - Clinical localization: 51 temporal, 17 extratemporal. MRI: single focus 58, multiple foci 2, negative 8 (T1; p. 3945).
- **OPMs** (§2.3, p. 3940; §4.3, p. 3950):
  - Marvel MEG 128 (Beijing X-MAGTECH): 64 sensors in dual mode, giving 128 channels; sensitivity better than 15 fT/√Hz.
  - A rigid spherical helmet with fixed sensor positions.
  - A multilayer permalloy cylinder with active compensation, holding residual fields within ±0.5 nT over a 50-cm radius.
  - Supine recording at 1000 Hz, with synthetic-gradiometer denoising and a 2–70 Hz band-pass.
- **Comparator and SNR.** None. The Discussion (p. 3949) says OPM-MEG SNR is comparable to or better than SQUID-MEG, citing Feys 2022, the 2023 ictal letter, Schwartz 2025 and Feys 2025.
- **Localization and concordance** (§3.2, p. 3945; T2, p. 3946; Fig. 4, p. 3947):
  - Concordance with iEEG: 90.0 % overall (Gwet AC1 .885, 95 % CI .875–.894).
  - Temporal 80.1 % (AC1 .723, .688–.758) vs extratemporal 92.0 % (AC1 .926, .911–.928).
  - By iEEG method: SEEG 89.9 %, ECoG 90.0 %.
  - Centroid distance: 1.98 ± .93 cm in concordant cases vs 4.12 ± 2.66 cm in discordant ones (p < .05). Temporal and extratemporal cases did not differ in centroid distance (p > .05).
- **Outcome**, 51 patients (§3.3, p. 3946; T3, p. 3948):
  - ILAE criteria: sensitivity 85.7 %, specificity 65.2 %, OR 11.25 (2.88–43.9).
  - Engel criteria: sensitivity 73.0 %, specificity 64.3 %, OR 4.86 (1.31–18.00).
- **Other numbers.** IED counts ranged from 3 to 162 (median 13); 5 patients had no IEDs (T1; p. 3945).
- **Read:** full text (local PDF of the version of record).
- **Notes.**
  - "Concordant" means any overlap of AAL regions (p. 3944), a generous criterion.
  - The authors attribute the lower temporal concordance to a larger sensor-to-source distance: mesial structures lie 4–6 cm deep, and the rigid spherical helmet sits farther from the temporal lobe (p. 3950).
  - There are no age-stratified results; the authors list this as a limitation (p. 3951).
  - Table 1 prints the female percentage as 68 (26/68 = 38 % [derived]).

---

## 4. Verified entries: modelling precedents

### 4.1 Zahran et al. 2022, *Sensors* 22(8):3093 (doi 10.3390/s22083093)

- **Design.** Simulation in Brainstorm with an OpenMEEG three-layer BEM, on the MRIs of 12 healthy subjects: one adult, one 9-year-old, and ten infants aged 1 month to 2 years (§2.1).
  - Sources: 15,002 cortical dipoles at 5 mm depth (§2.2).
  - Forward metrics: topography power, total information, sensitivity to source orientation, and volume-current ratios.
  - Resolution metrics: peak localization error (PLE) and spatial dispersion (SD) (§2.6).
- **OPM model** (§2.3).
  - Simulated 4He sensors: a 10-mm cube with 4 integration points, noise 40 fT/√Hz, 3 mm from the scalp.
  - Arrays of 102 sensors measuring the normal component (OPMn), one tangential component each (OPMt1, OPMt2), or all three (OPMa, 306 channels).
- **Comparator** (§2.3). 102 SQUID magnetometers with noise 5 fT/√Hz at 2 cm. For the children and infants, two helmet placements: "optimized" (2 cm all round, like a babyMEG) and "standard" (head against the top of the helmet).
- **SNR measure.**
  - Topography power, the squared norm of the noise-whitened topography, proportional to SNR when sensor noise is equal; and total information from the per-channel SNRs (§2.6).
  - Sensor noise only; there is no brain noise.
- **Adult results.**
  - Topography power: OPMa 8.9× the SQUID array; OPMn 3.6× and 4.1× OPMt1 and OPMt2 (§3.2.1).
  - Total information: OPMa 1.137× the SQUID array, but the SQUID array 1.133× OPMn (§3.2.2).
  - Sensitivity to the least- and most-sensitive orientations: 2.98× and 3.71× higher for OPMa (§3.2.3).
  - PLE: 2.3 mm (OPMa), 2.4 mm (OPMn), 2.8 mm (SQUID). SD: 3, 3.5 and 4.2 mm (§4.2).
- **Head size.**
  - OPMa beat the SQUID arrays at every age, and the gap grew as the head got smaller (age effect p < 0.001; §3.3).
  - For infants, OPMn exceeded the standard SQUID array in total information, which it did not do for the adult (§4.3).
  - PLE and SD differences between OPMa and the standard SQUID array: 0.17 and 0.9 mm at 9 years, 0.18 and 1 mm at 1 month (§4.3).
- **Read:** full text (PMC9024855).
- **Notes.**
  - With sensor noise only, this is a comparison of sensor-noise-limited arrays.
  - The OPM–scalp distance is measured from the gas-cell bottom in §2.3 but from the cell centre in the array description.
  - The Jas preprint takes some of its head sizes from this paper (`jas2026.md` §3.1).

### 4.2 Westin et al. 2023, *Clin Neurophysiol* 156:143–155 (doi 10.1016/j.clinph.2023.10.006)

- **Design.** Simulation on one adult anatomy, the MNE sample subject, with a three-compartment BEM (§2.1–2.2).
  - Ten epileptic networks, each with propagating seizure activity, plus isolated spikes from growing cortical patches at 0.77 nAm/mm² (§2.4).
  - Brain noise from a 500-s resting MEG/EEG source estimate (§2.4.3).
  - Detection by one clinician's visual review plus time–frequency analysis, then dipole localization (§2.6).
  - Network sites: frontopolar, basal frontal, lateral frontal, central, lateral temporal, mesial temporal, insula, cingulum, parietal, occipital.
- **OPM model** (§2.3.1). 128 sensors, 10-mm cubes at 1-cm spacing, 1-mm standoff, sensing volume 5 mm from the outer wall. Noise was set conservatively at 30 fT/√Hz.
- **Comparator** (§2.3.2–2.3.3). A 306-channel VectorView (magnetometers 5 fT/√Hz, gradiometers 5 fT/cm/√Hz), a 60-channel EEG, and combinations of these.
- **Measure.** Visibility in the raw sensor data: the smallest active cortical area that gives a visible event. No numeric SNR is computed.
- **Spike detection** (§3.3; T2).
  - On-scalp, spikes were visible from 0.5 cm² patches, or 0.78 cm² for frontopolar and insula (mean 0.56 cm²).
  - Conventional MEG needed a mean 1.91 cm² and EEG 1.87 cm² (p < 0.01).
  - For the mesial temporal network, spikes appeared from 0.5 cm² on-scalp vs 1.54 cm² for conventional MEG and EEG.
- **Seizure detection** (§3.1). The onset patch (0.5 cm²) was visible on-scalp at every site except mesial temporal, where activity became visible from 0.79 cm². Conventional MEG and EEG first saw the propagated activity, at 2.01 cm².
- **Localization** (§3.2–3.3; T1–T2).
  - Earliest ictal activity, distance from dipole to onset: 4.7 mm (2.2–9.6) on-scalp, against 20.15 mm for MEG + EEG, 20.37 mm for MEG and 21.55 mm for EEG.
  - Spikes: no significant difference between modalities (means 4.74, 4.75 and 4.71 mm).
- **Read:** full text (open access, CC BY; publisher text).
- **Notes.** One anatomy and one reader. "Detection" here means visibility, not a clinical SNR. Brain noise is included, unlike Zahran 2022.

### 4.3 Matsubara et al. 2026, *NeuroImage* 333:121930 (doi 10.1016/j.neuroimage.2026.121930)

- **Design.** Retrospective modelling on clinical data (§2.1–2.6).
  - Subject-specific MNE-Python BEM forward models of the cerebrum and of ARCUS cerebellar surfaces.
  - SNR at each vertex for a 10-nAm source.
  - Brain-noise variance estimated from 60 s of resting clinical MEG/EEG (0.5–100 Hz) through a model of uniformly distributed noise sources.
- **Patients.** 62 were processed (mean age 26.2 ± 14.0 years, range 4.7–61.4; 25 male) and 54 analysed (§2.2, §3.1).
- **OPM model** (§2.7). Simulated FieldLine Gen1 sensors in three layouts:
  - 102 sensors at the SQUID magnetometer positions projected onto the scalp, with a 5-mm offset;
  - 70 at the EEG positions;
  - 20 in a 5 × 4 grid (10 × 8 cm) over the occipital region.
  - Instrument noise was not modelled (§2.6).
- **Comparator.** The clinical 306-channel VectorView and the 70-channel EEG; for the cerebellar layout, a matching SQUID subset.
- **SNR definition** (§2.5). SNR (dB) = 10 log₁₀[(a²/N) Σₖ bₖ²/sₖ²], where a = 10 nAm, bₖ is the forward field at sensor k and sₖ² its brain-noise variance (Goldenholz-type).
- **Results.**
  - With OPMs at SQUID-matched positions, no cerebral or cerebellar region had a higher SNR than SQUID (§3.2.2). The deep cerebral regions (medial temporal, cingulate, insula) showed no meaningful gain (§4.3).
  - With the cerebellar layout, SNR exceeded the −23 dB superficial-cortex reference in Crus I, Crus II, VIIb, VIIIa and VIIIb. It was higher than the matched SQUID subset in Crus I, Crus II, VIIb, VIIIa and vermis VII and VIII (§3.2.4).
  - SQUID cerebral median SNR ranged from −32 to −23 dB, with parietal at −23 dB (§3.2.1).
  - Head circumference correlated with the sensor-distance ratio (r = −0.56), and that ratio with the OPM–SQUID SNR difference (r = 0.53); both p < 0.001 (§3.3).
- **Read:** full text (PMC13293240), with symbols and layouts re-checked in the Europe PMC XML.
- **Notes.** Ahlfors and Jas, authors of the Jas preprint, are co-authors. The SNR here includes brain noise only.

---

## 5. The preprint's own measured SEF result

Jas et al. (2026) measured one healthy 62-year-old man. They recorded the median-nerve N20 with 11 single-axis FieldLine Gen2 OPMs and with the 102 magnetometers of a MEGIN Triux neo, the magnetometers alone (`jas2026.md` §5.1–5.3). The paper prints no noise ratio. The values below come from the Fig. 8 bars, digitised to about ±2 % for the example channels C3 (OPM) and MEG0431 (SQUID), and tabulated in `jas2026.md` §5.6.

- **Noise ratio.** σ_OPM/σ_SQUID was about 36.1/7.9 ≈ 4.6 at the closest placements (provenance register row J-eta-meas). It was about 4.1 for the repeat runs and about 10 in the empty room.
- **N20 SNR.** About 9.5 for the OPM against 16.6 for the SQUID magnetometer.
- **Noise and distance.** OPM noise barely changed with distance: 91 % and 101 % of the closest-run value, against 79 % and 59 % for the SQUID.
- **N20 amplitude.** About 3× larger with the OPM according to the text, or ≈ 2.6 from the bars.
- **OPM floor.** The empty-room spectrum is at ≈ 29–31 dB re fT²/Hz (Fig. 9), which is about 28–35 fT/√Hz if the axis unit is dB re 1 fT²/Hz. The corresponding SQUID level is ≈ 3.5–6.3 fT/√Hz.
- **Equal-SNR depth.** In the adult sphere model, η = 4.1 would give an equal-SNR depth of 19.3 mm and η = 4.6 one of 17.2 mm (`jas2026.md` §5.6, derived).

This is a single, pre-selected channel pair for an evoked response (N20 latency 24 ms), not the peak channel of averaged IEDs used by the clinical studies. The comparator is magnetometers alone, as in Feys 2022 and 2025, and here the SQUID had the higher SNR.

The noise definition is ambiguous (`jas2026.md` §5.5 and §8 item 16; register row U-J7).

- The Methods (p. 12) define σ_N20 as the trial-to-trial SD of the unaveraged data at the N20 latency.
- The Fig. 8 caption (p. 19) calls it the standard error of the mean.
- The magnitudes favour the SEM of the roughly 500-trial average.
- The choice changes the absolute SNRs by about √N.
- The OPM/SQUID ratios are unaffected only if both systems kept the same number of trials. The preprint does not report how many trials remained after rejection.

---

## 6. Synthesis (from the tables above only)

### 6.1 By age group

- **Infants.** One infant was recorded (Feys 2023, *Clin Neurophysiol*), with no SQUID recording of her.
  - Her IED SNR was 11.1. Her amplitude, 2.5 pT, was lower than the school-aged children's 3.2–9.9 pT, with a 4-mm sensor gap.
  - The statement that OPM SNR is higher in infants rests on comparing her with other children's cryogenic recordings.
  - The only infant-to-adult OPM-vs-SQUID comparison is a simulation (Zahran 2022), with sensor noise only. There, the three-axis array's advantage grows as the head shrinks, and the normal-axis array overtakes a standard-placed SQUID array in information only in infants.
- **School-aged children.** One study compared OPM and SQUID in the same patients: Feys 2022, with 5 children recorded sequentially, 32 OPMs against 102 magnetometers.
  - Amplitude was 2.3–4.6× higher in all five.
  - SNR was 27–60 % higher in 4 of 5.
  - Hillebrand 2023, Feys 2025, Shen 2026 and Matsubara 2026 include children but report no child-specific comparison. Matsubara 2026 reports head-size effects instead: smaller heads gained more in the posterior cerebellum.
  - The other paediatric papers are single cases without a SQUID comparison.
- **Adults and mixed cohorts.** Studies that compared the same patients disagree on SNR.
  - OPM higher: Ren 2025 (abstract, 39 patients with IEDs in both systems) and one single spike in Badier 2023.
  - Higher in 1 and similar in 2 of 3 tested patients: Feys 2025.
  - Slightly lower: Hillebrand 2023 (0.92–0.98× [derived]) and Schwartz 2025 (0.80× at the group level [derived]; deep averaged spike 6.8 vs 9.0).
  - Amplitude was higher with OPMs in every comparison that reports it.

### 6.2 By sensor technology

- **Rubidium (QuSpin Gen-2/Gen-3).**
  - Sensor counts per patient: 4 (Feys 2023, *Front Neurosci*), 6 (Hillebrand), 8 then 15 (Vivekananda), and 24, 26, 32, 48 and 43–53 (the other Feys studies).
  - Noise, where stated: 7–13 fT/√Hz (Hillebrand) and 15 fT/√Hz (Feys 2023, *Front Neurosci*).
  - With many sensors and a magnetometer comparator (Feys 2022, 2025), OPM SNR was higher or similar.
  - With 6 sensors and a gradiometer comparator (Hillebrand), it was slightly lower.
- **Helium-4 (prototypes; MAG4Health in Schwartz 2025).**
  - All three studies used 4 working sensors on the scalp.
  - Noise was below 45–50 fT/√Hz on two axes and about 200 fT/√Hz on the third, or 60–65 fT/√Hz in Feys 2023 (*Front Neurosci*).
  - Against the nearest SQUID channels, SNR was higher for one single spike (Badier, 7 vs 5) and lower at the group level (Schwartz, 6.72 vs 8.37).
  - The one direct Rb-vs-4He comparison found about half the SNR with 4He (11.4 vs 21.3), with equal residual noise in the IED component (0.5 pT). The authors attribute the difference to sensor placement.
- **128-channel dual-axis whole-head systems.** Ren 2025 reports higher OPM amplitude and SNR (abstract; sensor type not given). Shen 2026 (X-MAGTECH, < 15 fT/√Hz) made no SQUID or SNR comparison.
- **Amplitude, across all technologies.** Every OPM-vs-SQUID amplitude comparison favours the OPM:
  - Feys 2022: 2.3–4.6×.
  - Schwartz: 2.3× for the group, 1.4–5.3× for individuals.
  - Feys 2025: higher in 9 of 9, nominally; 1.04–9.0× [derived].
  - Badier: higher.
  - Ren: higher.

  No study reports a lower OPM amplitude.

### 6.3 By region

- **Temporal, and above all mesial temporal, sources gave the weakest OPM results wherever they were examined.**
  - Feys 2025: 9.4 in TLE against 14.3 in extratemporal epilepsy (p = .02).
  - The amygdala case (Feys 2024): 2.75, against 11.1–16.7 for neocortical OPM recordings.
  - The deep averaged spike in Schwartz 2025: 6.8 (OPM) vs 9.0 (SQUID).
  - Badier 2023: hippocampus-only spikes were seen by neither system.
  - Shen 2026: concordance was 80.1 % for temporal and 92.0 % for extratemporal epilepsy, with no difference in centroid distance.
- **OPM against cryogenic values by temporal sub-region.** Only Feys 2025 reports these, by individual patient [derived from T2].
  - The three patients with medial-temporal or temporal-pole sources (8, 9, 10) all had nominally lower OPM than cryogenic SNR: 21.5 vs 22.5, 9.6 vs 13.5, 3.6 vs 4.4.
  - Five of the six other patients with both values had nominally higher OPM SNR.
  - Only patients 9 and 10 of that group were tested, and neither difference was significant.
- **Neocortical and extratemporal sources.**
  - OPM SNR was higher in 4 of 5 extratemporal children (Feys 2022).
  - In the sensorimotor patient of Hillebrand 2023, OPM SNR was lower (4.05 vs 4.39) and so was the spike-wave index (11.03 vs 24.50).
- **Modelling.**
  - Westin 2023: mesial temporal is the only site where the smallest onset patch was missed on-scalp. Mesial temporal spikes were still seen on-scalp from a smaller area than with conventional MEG (0.5 vs 1.54 cm²).
  - Matsubara 2026: with brain noise only, proximity alone gives no SNR gain for medial temporal, cingulate or insula sources.
- **Caveat.** Most of these contrasts compare regions within the OPM data (temporal against extratemporal). They are not region-by-region comparisons of OPM against SQUID.

### 6.4 By comparator and SNR definition

- **102 SQUID magnetometers, peak channel** (Feys 2022, 2025): OPM SNR was higher in 4 of 5 children, and higher or similar in the 3 tested TLE patients.
- **Planar gradiometers** (Hillebrand 2023): OPM SNR was 2–8 % lower [derived], which the authors call comparable.
- **The 4 SQUID channels nearest the 4 OPMs** (Badier 2023 with 4D magnetometers; Schwartz 2025 with a 275-channel CTF, and the 4D system for the SEEG case): SQUID SNR was higher for the group and for the deep spike, and OPM SNR higher for one single spike.
- **An unspecified SQUID** (Ren 2025, abstract): OPM SNR higher.
- **No SQUID** (Vivekananda 2020; the ictal, infant and Rb-vs-4He papers of 2023; Feys 2024; Shen 2026): any SNR advantage these papers state is by reference to other patients or other studies.
- **SNR definitions differ.**
  - Peak amplitude over baseline SD of averaged IEDs at the peak channel (Feys 2022, 2025).
  - Per-IED maximum over a 20-s or 2-s baseline SD (Schwartz, Badier).
  - A detector Z-score (Hillebrand).
  - An ICA-component virtual-sensor SNR (Feys 2023, *Front Neurosci*).
  - A seizure-to-background ratio (the ictal letter).
  - Not stated (Ren, Feys 2024, the infant letter).

  None of the clinical SNRs is a whole-array or source-level measure.
- **Modelling comparators.**
  - Against 102 magnetometers with sensor noise only (Zahran), the three-axis OPM array is ahead everywhere, but the normal-axis array is behind the SQUID array in adult information.
  - With brain noise taken from real recordings, on-scalp sensors still detect smaller active areas (Westin, which also adds 30 fT/√Hz sensor noise). Yet moving sensors closer at SQUID-matched positions raised SNR in no region (Matsubara, brain noise only).

### 6.5 What the table does not contain

- No study compares OPM and SQUID in the same infant. For school-aged children there is one such study, with five children.
- Only Feys 2025 gives OPM and SQUID SNR region by region, and then only per patient.
- The largest same-patient comparison (Ren 2025) could be checked only as an abstract.
- No clinical study compares whole-array or source-level detectability.

## 7. Cautions about secondary citations

- **Shen 2026.**
  - The Discussion (p. 3949) supports "comparable to or better" OPM SNR with four citations. One of them, the 2023 ictal letter, has no SQUID comparison. Another, Schwartz 2025, found slightly higher SQUID SNR.
  - The Introduction (p. 3938) cites Feys 2022 and the ictal letter as showing that OPM-MEG outperforms SQUID-MEG in children.
- **Schwartz 2025** (p. 1670) cites Feys 2022 and the infant letter for significantly higher OPM SNR in children and infants. The infant letter has no SQUID data for the same infant.
- **The 2.3–4.6× amplitude range** belongs to Feys 2022. Hillebrand 2023 reports no amplitudes.
- **A 2026 narrative review** (Kalinina et al., *Sensors* 26:5490) describes Feys 2022 as 14 school-aged children recorded with concurrent SQUID-MEG. The paper reports 5 children recorded sequentially.

## 8. Identified but not tabulated

These have no OPM-vs-SQUID amplitude or SNR comparison for epileptiform activity, or could not be read.

- **Mellor et al. 2024, medRxiv** (10.1101/2024.10.17.24315226; preprint read in full).
  - 8 adults with temporal lobe epilepsy, recorded with QuSpin Gen-2 and Gen-3 OPMs in scanner-casts: 38, 42 and 38 sensors for patients 1–3, and 32 + 25 for patient 4.
  - IEDs appeared in 4 patients (2, 9, 4 and 2 IEDs); in 3, the localization agreed with MRI and EEG. Two had ictal events.
  - There was no SQUID recording and no SNR.
- **Mellor et al. 2026, *Clin Neurophysiol* 184:2111504.** A letter, "Experiences of ictal OP-MEG"; no abstract; not read.
- **Chen et al. 2025, *Neurophysiol Clin* 55(4):103087.** A letter on OPM-MEG in presurgical evaluation; its keywords include SEEG. No abstract, and the full text was not accessible.
- **Westin et al. 2020, *Clin Neurophysiol* 131:1711–1720.** On-scalp **high-Tc SQUID**, not OPM, in one patient [abstract]. Conventional MEG recorded 24 IEDs. The on-scalp system recorded 47: 16 also seen on EEG and 31 unique.
- **Bonnet et al. 2025, *Front Med Technol* 7:1548260.** A whole-head 4He system tested in simulation (three head sizes, child to adult) and with a phantom; equivalent to cryogenic MEG in detection and localization [abstract]. Not an epilepsy study. It is the source of the floor Schwartz 2025 quotes.
- **Wang et al. 2024, *J Neural Eng* 21:066033.** Auditory and visual evoked fields in 100 participants. OPM SNR was lower than SQUID for the auditory responses and no different for the visual [abstract]. No IEDs.
- **Reviews and commentaries** (no new data):
  - Pedersen et al. 2022, *Epilepsia* 63:2745;
  - Feys & De Tiège 2024, *Dev Med Child Neurol* 66:298;
  - Feys et al. 2023, *Epilepsia* 64:3414;
  - Bagić et al. 2023, *Epilepsia* 64:3155;
  - Sequeiros & Hamandi 2024, *J Neurol* 271:1492;
  - Kalinina et al. 2026, *Sensors* 26:5490.

## 9. Verification log

Every entry was read once to extract it, then re-read against its source before this note was finished (2026-10-04).

| Study | Source read | Level |
|---|---|---|
| Vivekananda 2020 | PMC full text | full text |
| Feys 2022 | publisher and PubMed abstract; medRxiv v1 full text | published abstract + preprint full text |
| Feys 2023, ictal letter | publisher text via the DOI | full text |
| Hillebrand 2023 | PMC full text; Table 1 from the Europe PMC XML; medRxiv abstract | full text |
| Feys 2023, infant letter | publisher HTML (ScienceDirect) | full text |
| Badier 2023 | PMC full text; Europe PMC XML for subscripts and legends | full text |
| Feys 2023, *Front Neurosci* | PMC full text; Europe PMC XML | full text |
| Feys 2024 | accepted-manuscript text via the DOI; medRxiv v1 full text | full text |
| Feys 2025 | local PDF, version of record | full text |
| Ren 2025 | PubMed and publisher abstract, highlights; the open-access full text could not be read (the publisher's site refused automated access, 2026-10-05) | abstract only |
| Schwartz 2025 | local PDF, version of record | full text |
| Shen 2026 | local PDF, version of record | full text |
| Zahran 2022 | PMC full text | full text |
| Westin 2023 | publisher text (open access) | full text |
| Matsubara 2026 | PMC full text; Europe PMC XML | full text |
