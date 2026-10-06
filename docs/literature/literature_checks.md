# Literature checks

This note records five checks of the literature behind the manuscript: the full text of Ren et al. (2025), the publication status of Jas et al. (2026), the white noise of the MEGIN/Elekta Neuromag systems, the OPM vapour-cell size and sensing-centre standoff, and the frequency band of each clinical OPM noise figure. For each question it gives the finding, the full citation with DOI, the values with their place in the source, and how each value was verified. Searches and reading were done on 2026-10-05.

## Summary

| # | Question | Finding | Consequence for the manuscript |
|---|---|---|---|
| 1 | Ren et al. (2025): is the full text accessible? | Yes in principle: the article is open access (CC BY 4.0) at the publisher. It is not in PMC or Europe PMC. The publisher's site refused every automated reader available to this check, so the full text was **not read**. | Keep the "abstract only" label. The systems, sensor noise, channels and SNR definition are not verifiable. Indirect evidence points to a MEGIN TRIUX and the 128-channel X-Magtech OPM system. |
| 2 | Jas et al. (2026): published since? | **No.** Only bioRxiv version 1 (21 August 2026) exists. There is no journal version and no data or code deposit. | No value changes. The figure readings stay labelled as readings of the preprint; their method and uncertainty are in Section 2. |
| 3 | Neuromag sensor noise | The manufacturer's TRIUX datasheet gives typical white noise of 3.5 fT/√Hz (magnetometers) and 3.6 fT/(cm √Hz) (gradiometers), and an 18-mm average coil-to-surface distance. These are exactly the study's values. The VectorView manual gives guaranteed maxima only. | The three hardware values can be cited to the TRIUX datasheet; they are no longer "not verified". |
| 4 | OPM cell size and standoff | Cell centre to the sensor's outer surface: 6–6.5 mm for QuSpin, 5 mm for FieldLine. QuSpin cells are 3 × 3 × 3 mm³. No accessible source gives FieldLine's cell size. No device source gives a 10-mm cell. | The 7-mm sensing-centre height is supported. Describe the 10-mm cube as a modelling choice (an averaging volume), not a device cell; its effect at the field peak is at most 1.9 %. |
| 5 | Band of each clinical OPM noise figure | No clinical paper attaches a band to the noise figure of its own sensors. Bands appear only in the device sources behind the figures: 1–100 Hz (QuSpin, first generation), 3–100 Hz (QuSpin Gen-2/Gen-3) and 10–130 Hz (FieldLine). One clinical introduction gives 3–100 Hz and 1–1,500 Hz in general terms. | State that the clinical figures are device specifications over bands that start at 1–10 Hz, and that the model's 15 fT/√Hz is a white level over its 1–40 Hz band. |

---

## 1. Ren et al. (2025): full text

**Citation.** Ren J, Ding M, Peng Y, Sun C, Yang C, Zhou S, Tian J, Wang Q, Li Z (2025). A comparative study on the detection and localization of interictal epileptiform discharges in magnetoencephalography using optically pumped magnetometers versus superconducting quantum interference devices. *NeuroImage* 312:121232. doi:[10.1016/j.neuroimage.2025.121232](https://doi.org/10.1016/j.neuroimage.2025.121232). PMID 40254146. Published online 18 April 2025.

**Is the full text accessible?**

- **Licence: open access.**
  - Crossref holds a CC BY 4.0 licence for the version of record, valid from 18 April 2025. Its copyright assertion names the authors as copyright holders (published by Elsevier).
  - OpenAlex and Semantic Scholar list the article as gold open access under CC BY.
  - The other records disagree. PubMed shows only an Elsevier copyright line and no licence. Europe PMC records the CC BY licence but labels its publisher link "subscription required". The licence record is the publisher's own deposit. Read: metadata.
- **Repositories.** The article is not in PMC and not in Europe PMC; Europe PMC record 40254146 has no full text.
- **Reading attempts.**
  - The ScienceDirect article page refused automated access in three ways: HTTP 403 to a page fetcher, an access-error page in a desktop browser, and a crawler error.
  - Elsevier's article API requires a key.
  - The Internet Archive holds two captures of the page (May and August 2025), both HTTP 403.
  - The DOAJ record links back to ScienceDirect. The Beihang University research portal shows the abstract only.
  - We did not try to get around the publisher's block.
- **Result.** The following are **not verifiable from accessible sources**: the OPM sensor type and noise, the SQUID system, the comparator channels, the SNR and amplitude definitions, the effect sizes and the table values. The extraction in `epilepsy_opm_studies.md` §3.10 stays abstract-only and needs no correction.

**What accessible records add.** These are indirect; they are not the article's text.

- **PubMed's competing-interest statement.**
  - One author (M. Ding) reports board membership with Beijing X-Magtech Technology Limited.
  - X-Magtech makes the Marvel MEG OPM system. Shen et al. (2026) used its 128-channel version.
  - The statement does not say which system the study used. Read: metadata (PubMed XML).
- **A meeting abstract by the same group.** Ren J, Peng Y, Sun C, Yang C, Wang Q, Ding M, Li Z (2024). Optically pumped magnetometers versus cryogenic magnetoencephalography in the detection of interictal epileptic discharges and signal-to-noise ratio. American Epilepsy Society Annual Meeting 2024, abstract 2.489 (presented 8 December 2024). [aesnet.org](https://aesnet.org/abstractslisting/optically-pumped-magnetometers-versus-cryogenic-magnetoencephalography-in-the-detection-of-interictal-epileptic-discharges-and-signal-to-noise-ratio). Read: abstract; it is a meeting abstract, not a peer-reviewed paper.
  - Patients: 21 with refractory focal epilepsy at Beijing Tiantan Hospital, 9 of them female; median age 35 years (range 8–50).
  - Systems: a Neuromag TRIUX (MEGIN) and a Marvel MEG (Beijing X-Magtech), recorded separately under similar conditions. IEDs were marked independently in each recording.
  - IEDs were found in all 21 patients with both systems.
  - Averaged IED SNR was 10.74 ± 3.62 with OPM and 9.31 ± 4.30 with SQUID (t test, p = 0.035); the ratio is 1.15 [derived].
  - Not given: the SNR definition, the channels used and the sensor noise.
  - The abstract does not relate these patients to the *NeuroImage* cohort. Its age description (median 35 years) differs from the article's (mean 23.7 ± 8.7 years), so it cannot be assumed to report a subset of that cohort.
- **A municipal news item.**
  - The Beijing Municipal Science and Technology Commission reported the market approval of the 128-channel Marvel MEG on 24 July 2024 ([kw.beijing.gov.cn](https://kw.beijing.gov.cn/xwdt/kcyx/xwdtyqqy/202407/t20240724_3821972.html), in Chinese).
  - The item says the device's registration clinical study was run at Beijing Tiantan Hospital, Ren et al.'s institution. It does not cite the article.

**Use in the manuscript.** Ren et al. (2025) is cited from its abstract only (Table S19), and no SNR definition, channel choice or sensor noise is attributed to the study. The article is open access, so its full text can be checked in an ordinary browser.

---

## 2. Jas et al. (2026): publication status

**Citation.** Jas M, Matsubara T, Sohrabpour A, Sundaram P, Mody M, Ahlfors SP (2026). Signal-to-noise ratio of event-related fields in on-scalp and off-scalp MEG. *bioRxiv* 2026.08.17.744953, version 1, posted 21 August 2026; not peer reviewed; CC BY 4.0. doi:[10.64898/2026.08.17.744953](https://doi.org/10.64898/2026.08.17.744953).

**Status on 2026-10-05: not published in a journal.**

- **bioRxiv.** The API lists one version (v1, 21 August 2026) and no published article. Its endpoint for published versions returns none.
- **Crossref.** The record is posted content. It has no "is preprint of" relation and no citing records.
- **PubMed.** No record has this title. An author search (Jas M and Ahlfors SP, with OPM, on-scalp or SNR terms) returns only three other records: a 2025 *Sensors* paper on nulling coils, its preprint, and Matsubara et al. (2026).
- **Europe PMC.** It holds only the preprint record PPR1302855; its full-text service returned no text for it.
- **Web search** for the exact title finds only the preprint.
- **Data and code.** No deposit was found. The preprint itself has no data or code statement (`jas2026.md`, header notes).
- Read: metadata and web search. The preprint was not re-read for this note, because bioRxiv rate-limited the request. Its values come from the project's extraction, `jas2026.md`.

**Consequence: no value changes.** The table lists what the report takes from the preprint and how each value was obtained (`jas2026.md` §§5.5, 5.6 and 9).

| Value used in the report | Kind | Place in the preprint | Extraction and uncertainty |
|---|---|---|---|
| Equal-SNR depths of Fig. 4B and 4C (34 and 19 mm); thresholds η0 = 1.7 and η1 = 5.3 (Fig. 4E) | printed | text and caption, pp. 14–15 | Printed values. The reimplementation gives 35.67 mm, 19.82 mm, 1.6829 and 5.3086 (`jas2026.md` §9). |
| N20 SNR ≈ 16.6 (SQUID magnetometer MEG0431) vs ≈ 9.5 (OPM C3), closest placements | figure reading | Fig. 8C, p. 19 | Bar heights were measured on a 600-dpi rendering of the PDF. Each is within about ±2 %. The bar values reproduce the eight printed percentages (p. 19) to within 1 percentage point. The SNRs carry the preprint's ambiguity about the noise definition (`jas2026.md` §8, item 16). |
| OPM empty-room floor ≈ 28–35 fT/√Hz | figure reading and conversion | Fig. 9, p. 20 | The trace reads ≈ 29–31 dB on an axis labelled fT²/Hz (dB). The dB reference is not stated. Assuming 0 dB = 1 fT²/Hz, the amplitude is 10^(L/20) fT/√Hz, giving 28.2–35.5 fT/√Hz [derived]. Each dB of reading error changes the amplitude by −11 % or +12 % [derived]. |
| OPM low-frequency component ≈ 200–224 fT/√Hz near 5 Hz | figure reading and conversion | Fig. 9 and text, p. 20 (also p. 23) | The peak reads ≈ 46–47 dB. The same conversion gives 199.5–223.9 fT/√Hz [derived]. |
| "4 to 330 Hz" | printed | p. 12 (4-Hz high-pass in processing); pp. 10–11 (330-Hz anti-aliasing low-pass at acquisition) | Printed values, not readings. No low-pass filter for the analysis is stated. |

Next check: when a version of record, data or code appear.

---

## 3. White noise of the MEGIN/Elekta Neuromag systems

**Values in this study.**

- Magnetometers: 3.5 fT/√Hz.
- Planar gradiometers: 3.6 fT/(cm √Hz).
- Pick-up coil to the room-temperature helmet surface: 18 mm.
- Sources: `configs/g2_adult.toml`, sections `[sensors]` and `[head_position]`; register rows HW-mag-noise, HW-grad-noise and HW-18mm, all marked "transcribed", and U-HW1.

**TRIUX datasheet: the study's values verified.**

Elekta. *Elekta Neuromag® TRIUX*. Datasheet, document NM23083B-A (system art. no. NM23900N), undated. Manufacturer document; copy hosted by NatMEG, Karolinska Institutet: [natmeg.se](https://natmeg.se/onewebmedia/NM23083B-A%20Elekta%20Neuromag%20TRIUX%20datasheet.pdf). Read: full text of the datasheet, section on the sensor array.

- **"Typical white noise levels":** gradiometer 3.6 fT/cm/√Hz, magnetometer 3.5 fT/√Hz.
- **Guaranteed levels.** These are the same numbers for both sensor types, in fT/(cm √Hz) for gradiometers and fT/√Hz for magnetometers:

  | Band | 96 % of sensors | All sensors |
  |---|---|---|
  | 1–10 Hz | below 12 | below 20 |
  | 60–70 Hz | below 5 | below 10 |

- **Geometry and range.**
  - Average distance between pick-up coils and the room-temperature surface: 18 mm.
  - Gradiometer baseline: 17.0 mm.
  - Average distance between sensor elements: 35 mm.
  - Dynamic range: ±20 nT.

Elekta Oy (2011). *Elekta Neuromag® TRIUX Technical Manual*, article number NM24132A, February 2011. Manufacturer document; copy hosted by NatMEG: [natmeg.se](https://natmeg.se/onewebmedia/NM24132A%20Triux%20TM.pdf). Read: §1.1.2 (Sensors) and §1.1.4 (Dewar).

- **Bands of the guaranteed maxima.** White noise is specified over 60–70 Hz and low-frequency noise over 1–10 Hz. The numbers are those of the datasheet. The manual gives no typical value.
- **Geometry.**
  - Coils: 28 mm × 28 mm; gradiometer baseline 17 mm.
  - Average triple-sensor spacing: 35 mm.
  - Distance from the liquid-helium surface to the outer surface against the head: typically 18 mm (§1.1.4).

**VectorView manual: guaranteed maxima only.**

Elekta Neuromag Oy (2005). *Elekta Neuromag® System Hardware: Technical manual*, Revision F (NM20216A-F), September 2005. The manual states that the system was previously sold as the Neuromag Vectorview. Copy hosted by the MRC Cognition and Brain Sciences Unit, Cambridge: [imaging.mrc-cbu.cam.ac.uk](https://imaging.mrc-cbu.cam.ac.uk/meg/VectorviewDescription?action=AttachFile&do=get&target=HardwareTechnical.pdf). Read: §1.1.2 (Sensors).

- **White noise (60–70 Hz).** At most 5 fT/(cm √Hz) for gradiometers and 5 fT/√Hz for magnetometers in 96 % of channels; at most 10 in all channels.
- **Low-frequency noise (1–2 Hz).** At most 12 in 96 % of channels; at most 20 in all.
- **Geometry.**
  - Sensor distance from the outside helmet surface: 18 mm on average.
  - Coils: 28 mm × 28 mm; gradiometer baseline 17.0 mm.
  - Average triple-sensor spacing: 34 mm.
  - Effective area: 2.7 cm² per gradiometer coil and 7.6 cm² per magnetometer.
- The manual's sensor section gives no typical white-noise value.

**Peer-reviewed literature.**

- One Europe PMC full-text search was made: Neuromag, VectorView or TRIUX, with noise in fT/√Hz units; 17 hits were screened. It found no peer-reviewed measurement of these systems' white noise.
- Three modelling papers assume about 3 fT/√Hz for the magnetometers. This is lower than the datasheet's typical 3.5 fT/√Hz.
  - Pfeiffer C, et al. (2018). Localizing on-scalp MEG sensors using an array of magnetic dipole coils. *PLoS One* 13:e0191111. doi:[10.1371/journal.pone.0191111](https://doi.org/10.1371/journal.pone.0191111). The Materials and methods (data generation, signal-to-noise ratio) state that TRIUX magnetometers typically have about 3 fT/√Hz, with no source. Read: full text (PMC5944911).
  - Riaz B, Pfeiffer C, Schneiderman JF (2017). Evaluation of realistic layouts for next generation on-scalp MEG: spatial information density maps. *Scientific Reports* 7:6974. doi:[10.1038/s41598-017-07046-6](https://doi.org/10.1038/s41598-017-07046-6). Elekta SQUID noise of 3 fT/√Hz is used as a model parameter (Fig. 1 legend and parameter table). Read: full text (PMC5539206).
  - Roos S, Hämäläinen M, Iivanainen J (2026). Efficacy of optically pumped magnetometers in detecting activity from the cerebellar cortex. *Human Brain Mapping* 47:e70514. doi:[10.1002/hbm.70514](https://doi.org/10.1002/hbm.70514). In the Methods, SQUID magnetometers are assigned 3 fT/√Hz as a typical value. Read: full text (PMC13081696).

**Not found.**

- A datasheet for the TRIUX neo: the MEGIN pages reached give no noise figures.
- Any typical white-noise value for the VectorView.

**Use in the manuscript.**

- The manuscript gives the Neuromag values as the typical white-noise levels of the manufacturer's TRIUX datasheet (3.5 fT/√Hz and 3.6 fT/(cm √Hz)), with the datasheet's 18-mm coil-to-surface distance, and cites the datasheet (register U-HW1).
- **Two cautions go with these values.**
  - These are typical, not guaranteed, values. The guaranteed maxima are 5 (60–70 Hz) and 12 (1–10 Hz) for 96 % of channels.
  - The white level is specified at 60–70 Hz, so the 1–40 Hz analysis band may hold more SQUID noise. The study's sensitivity analysis with the measured empty-room spectrum covers this. The register row on Neuromag noise from the measured spectrum gives median in-band RMS of 25.9 fT (magnetometers) and 20.0 fT/cm (gradiometers), against 20.7 fT and 21.3 fT/cm from the typical values.
- **Coil geometry.** The VectorView and TRIUX documents give the same coil size (28 mm), baseline (17.0 mm) and 18-mm spacing. The model's coil definitions (register HW-T3) were not compared with these here.

---

## 4. OPM vapour-cell size and sensing-centre standoff

**Values in this study.**

- A 10-mm cubic sensing volume.
- A sensing centre 7 mm from the helmet's inner surface: a 2-mm shell plus half the cell.
- Both are taken from Jas et al. (2026, p. 10) for FieldLine Gen2 sensors (`jas2026.md` §5.2; register rows J-opm-cell, J-opm-standoff, A-OPM-CELL and A-OPM-STANDOFF).
- On the adult, sensing centres sit 6.5–11.0 mm (median 7.0 mm) above the head surface (register A-OPM-CLEAR). The median is 6.99–7.01 mm in every anatomy and array (A-BEM-CONFORM).
- Gap variants add 3 and 6 mm (register A-OPM-GAP; `configs/g2_adult.toml`).

**Device values.**

| Device | Vapour cell | Cell centre to the sensor's outer surface | Package | Source (place) | Read |
|---|---|---|---|---|---|
| QuSpin QZFM, first generation | 3 × 3 × 3 mm³, ⁸⁷Rb (§2) | 6 mm (§5; Fig. 6) | 13 × 19 × 110 mm (p. 105481G-5) | Osborne J, Orton J, Alem O, Shah V (2018). Fully integrated, standalone zero field optically pumped magnetometer for biomagnetism. *Proc. SPIE* 10548:105481G. doi:[10.1117/12.2299197](https://doi.org/10.1117/12.2299197) | full text (copy on [quspin.com](https://quspin.com/wp-content/uploads/2016/08/J_Osb.pdf)) |
| QuSpin prototype | 3 × 3 × 3 mm³ | 6.5 mm: the sensitive volume lies 6.5 mm from the end of the housing, i.e. from the scalp (Methods, OPM; Fig. 1 legend) | 14 × 21 × 80 mm³ | Boto E, et al. (2017). *NeuroImage* 149:404–414. doi:[10.1016/j.neuroimage.2017.01.034](https://doi.org/10.1016/j.neuroimage.2017.01.034) | full text (PMC5562927) |
| QuSpin, in a 13-sensor wearable system | 3 × 3 × 3 mm³ (main text; Methods) | not stated | – | Boto E, et al. (2018). *Nature* 555:657–661. doi:[10.1038/nature26147](https://doi.org/10.1038/nature26147) | full text (PMC6063354) |
| QuSpin Gen-2 | not stated | 6 mm radially inwards from the casing, with a 2-mm tangential offset (Methods, co-registration) | 12.4 × 16.6 × 24.4 mm³ (Methods; Discussion) | Hill RM, et al. (2020). *NeuroImage* 219:116995. doi:[10.1016/j.neuroimage.2020.116995](https://doi.org/10.1016/j.neuroimage.2020.116995) | full text (PMC8274815) |
| QuSpin Gen-2 | – | 6.5 mm, from the cell centre to the outside of the housing | 12.4 × 16.6 × 24.4 mm | QuSpin QZFM Gen-2 product page, archived 24 February 2020 ([web.archive.org](https://web.archive.org/web/20200224030525/http://quspin.com:80/products-qzfm/)) | manufacturer page |
| QuSpin Gen-3, triaxial | 3 × 3 × 3 mm³ (§2.1) | not stated; their simulation assumes 6 mm above the scalp (§3.1.1) | 1.24 × 1.66 × 2.44 cm, 7 g (Fig. 2 legend) | Boto E, et al. (2022). *NeuroImage* 252:119027. doi:[10.1016/j.neuroimage.2022.119027](https://doi.org/10.1016/j.neuroimage.2022.119027) | full text (PMC9135302) |
| QuSpin Gen-3 | – | 6.5 mm, from the vapour-cell centre to the housing | 12.4 × 16.6 × 24.4 mm, the same as Gen-2 | QuSpin QZFM Gen-3 product page ([quspin.com/products-qzfm](https://quspin.com/products-qzfm/)); Gen-3 announcement, 2021 ([quspin.com/qzfm-gen-3](https://quspin.com/qzfm-gen-3/)) | manufacturer pages |
| QuSpin, Neuro-1 housing | – | 6.2 mm | 16.6 × 12.4 × 29.5 mm | QuSpin Neuro-1 page ([quspin.com/products/neuro-1](https://quspin.com/products/neuro-1/)) | manufacturer page |
| QuSpin Rb-OPMs, as described by a clinical group | – | 6.5 mm (Introduction) | 1.2 × 1.7 × 2.6 cm³ | Feys O, et al. (2023). *Frontiers in Neuroscience* 17:1284262. doi:[10.3389/fnins.2023.1284262](https://doi.org/10.3389/fnins.2023.1284262) | full text (PMC10715393) |
| FieldLine HEDscan sensor | microfabricated Rb cell; size not stated | 5 mm (§2) | footprint 13 × 15 mm (§2) | Alem O, et al. (2023). An integrated full-head OPM-MEG system based on 128 zero-field sensors. *Frontiers in Neuroscience* 17:1190310. doi:[10.3389/fnins.2023.1190310](https://doi.org/10.3389/fnins.2023.1190310) | full text (PMC10303922) |
| FieldLine HEDscan sensor | – | 5 mm: to the cell centre (web page); from the centre of the sensing volume to the sensor tip (specification sheet) | – | FieldLine HEDscan page ([fieldlinemedical.com/hedscan](https://fieldlinemedical.com/hedscan/)); FieldLine HEDscan specification sheet, 2023 (distributor copy: [physio-tech.co.jp](https://www.physio-tech.co.jp/assets/pdf/fieldline/FieldLine-Spec-Sheet-2023-2.pdf)) | manufacturer page and sheet |

Two reviews state the general picture:

- Tierney TM, et al. (2019). *NeuroImage* 199:598–608. doi:[10.1016/j.neuroimage.2019.05.063](https://doi.org/10.1016/j.neuroimage.2019.05.063).
  - The cell must be offset from the sensor walls so that the outer surface stays below about 40 °C.
  - Their worked sensitivity example uses a 0.027-cm³ cell, which is a 3-mm cube [derived].
  - Read: full text (PMC6988110).
- Bezsudnova Y, et al. (2022). Optimising the sensing volume of OPM sensors for MEG source reconstruction. *NeuroImage* 264:119747. doi:[10.1016/j.neuroimage.2022.119747](https://doi.org/10.1016/j.neuroimage.2022.119747).
  - Most commercial OPMs have cubic cells 2–3 mm on a side (Discussion).
  - Read: full text (PMC7615061); this is a secondary statement.

**Assessment.**

- **Standoff: supported.**
  - The device distances from the cell centre to the outer surface are 5 mm (FieldLine) and 6–6.5 mm (QuSpin; 6.2 mm in the Neuro-1 housing).
  - The study's 7 mm is the height of the centre above the scalp when the helmet's inner surface lies on the scalp.
  - A FieldLine sensor on a 2-mm shell gives 7 mm. This is the setting of Jas et al., and their 7 mm follows from FieldLine's published 5-mm standoff without any 10-mm cell.
  - A QuSpin sensor directly on the scalp gives 6.5 mm; behind a 2-mm shell, 8.5 mm [derived].
  - Hair, cap fabric and helmet fit add to these distances. The 3- and 6-mm gap variants cover that.
- **Cell: not a device dimension.**
  - QuSpin cells are 3-mm cubes (Osborne et al. 2018; Boto et al. 2017, 2018, 2022).
  - No accessible source gives FieldLine's cell size.
  - FieldLine places the cell centre 5 mm from the housing's outer surface, so a 10-mm cube centred there would reach that surface and leave no wall [inference].
  - The 10-mm figure appears only in Jas et al.'s description (p. 10, as extracted in `jas2026.md` §5.2).
  - In the model the cube is an averaging volume. The register gives its effect on the field at the peak as −1.9 % at 10 mm source distance, −0.4 % at 15 mm and −0.04 % at 30 mm (A-OPM-CELL). A 3-mm cell would lie closer to a point sensor [inference].
  - Suggested wording: a 10-mm cubic averaging volume, larger than the 3-mm cells of current alkali OPMs, with an effect of at most 1.9 % on the peak field. This needs no new computation.
- **Clinical array sizes.**
  - The largest clinical OPM arrays in `epilepsy_opm_studies.md` §1:
    - 64 dual-mode sensors with 128 channels (Shen et al. 2026);
    - a 128-channel system (Ren et al. 2025, abstract);
    - 43–53 sensors with 98–127 channels (Feys et al. 2025);
    - 48 sensors with 116 channels (Feys et al. 2024).
  - The dense array has 208 single-axis sites (register A-OPM-CLEAR). It exceeds every clinical array in both sensors and channels.

**Related device parameters.**

- **Package footprint** (register U-OPM-PACK, marked unverified).
  - QuSpin Gen-2/Gen-3: 12.4 × 16.6 mm (Hill et al. 2020; QuSpin pages; Boto et al. 2022).
  - FieldLine HEDscan: 13 × 15 mm (Alem et al. 2023).
  - The 17-mm minimum centre spacing exceeds the larger QuSpin side by 0.4 mm [derived]. It is therefore compatible with either package placed side by side, with no allowance for holders.
- **Bandwidth** (register A-OPM-BW: about 100–150 Hz, modelled as a first-order low-pass with a 100-Hz corner).
  - QuSpin, first generation: −3 dB at 135 Hz, with a first-order roll-off (Osborne et al. 2018, §5).
  - −3 dB at 130 Hz (Boto et al. 2018, Extended Data).
  - QuSpin Gen-2: 0 to about 130 Hz (Hillebrand et al. 2023, Methods).
  - QuSpin: about 150 Hz, first-order (QuSpin technical note, [quspin.com](https://quspin.com/zero-field-magnetometer-description/)).
  - FieldLine: DC to 150 Hz at −3 dB (HEDscan page); about 200 Hz in open loop and 350 Hz in closed loop (Alem et al. 2023, system performance).
  - The register's range is supported for QuSpin. FieldLine sensors respond over a wider band.

---

## 5. Frequency band of each clinical OPM noise figure

| Study (sensors) | Noise figure as printed, and where | Band stated in the paper | Band in the source behind the figure | Read |
|---|---|---|---|---|
| Hillebrand et al. 2023 (six QuSpin Gen-2 sensors, dual-axis) | 7–13 fT/√Hz; sensor bandwidth 0 to about 130 Hz (Methods, OPM sensors) | none | Cited to Osborne et al. 2018 (their ref. 18), which describes first-generation units (p. 105481G-4): 7–13 fT/√Hz for most units in single-axis mode; typically 10 fT/√Hz over **1–100 Hz**; about 30 % less sensitivity in dual-axis mode, a penalty Hillebrand et al. note (Methods; Discussion). QuSpin's Gen-2 specification: below 15 fT/√Hz over **3–100 Hz**, typically 7–10 (archived product page). | full text (PMC10030968); source in full text |
| Shen et al. 2026 (Marvel MEG 128: 64 sensors in dual mode) | better than 15 fT/√Hz (§2.3, p. 3940); the only noise figure in the paper | none | Not found. Press reports of the 2024 market approval give a sensitivity of the order of 10 fT, without a band (e.g. [VCBeat](https://www.vcbeathealth.com/article/34226)). | full text (PMC13525599, web page) |
| Schwartz et al. 2025 (MAG4Health triaxial ⁴He) | below 45 fT/√Hz on the radial and first tangential axes; about 200 fT/√Hz on the second tangential axis (Methods, §2.2.1) | none | Not given. For the group's earlier prototype see Badier et al. and Feys et al. below. | full text (PMC12514395) |
| Badier et al. 2023 (⁴He prototype) | better than 50 fT/√Hz on two axes and 200 fT/√Hz on the third (Materials and Methods, 4He-OPM recordings) | none for the noise; the sensor bandwidth is DC to 2 kHz | – | full text (PMC10748329) |
| Feys et al. 2023, *Front Neurosci* (QuSpin Gen-3 triaxial Rb; ⁴He prototypes) | Rb: below 23 fT/√Hz for current Rb-OPMs in general (Introduction); 15 fT/√Hz sensor noise (Discussion). ⁴He: below 50 fT/√Hz (Introduction); 60–65 fT/√Hz in this recording (Methods); 65 fT/√Hz (Discussion). | **3–100 Hz** for the Introduction's Rb figure; **1–1,500 Hz** for its ⁴He figure; none for the others | The Rb figure matches QuSpin's Gen-3 triaxial specification: below 23 fT/√Hz over 3–100 Hz. | full text (PMC10715393) |
| Vivekananda et al. 2020 (QuSpin, first and second generation) | about 10 fT/√Hz, the white-noise floor of recent OPMs (Introduction, citing Boto et al. 2018) | "white noise floor"; no band | Boto et al. 2018 give about 15 fT/√Hz, without a band. | full text (PMC7085997) |
| Feys et al. 2022; 2023 (ictal and infant letters); 2024; 2025. Ren et al. 2025 | no sensor noise reported; Ren et al.: abstract only | – | – | as in `epilepsy_opm_studies.md` |

**Device noise over a stated band, for comparison.**

- **QuSpin Gen-3 specification** (QuSpin product page).
  - Dual-axis: below 15 fT/√Hz over 3–100 Hz, typically 7–10.
  - Triaxial: below 23 fT/√Hz over 3–100 Hz, typically about 15 on all axes.
  - The 2021 announcement adds typical values by mode: 7–12 (single-axis) and 10–15 (dual-axis) fT/√Hz for the dual-axis variant; 15–20 fT/√Hz for the triaxial variant in triaxial mode.
- **Measured QuSpin Gen-3 noise floors** over **10–90 Hz**, excluding 45–55 Hz (Boto et al. 2022, §2.3 and Fig. 2):
  - triaxial sensors: 13.5 ± 0.8 fT/√Hz (x), 9.9 ± 1.4 (y) and 14.9 ± 2.0 (z), mean ± SD over 4 sensors;
  - dual-axis sensors: 11.3 ± 1.5 (y) and 13.8 ± 1.1 (z), over 14 sensors.
- **FieldLine HEDscan.**
  - The specification sheet (2023, distributor copy) lists 8 and 15 fT/√Hz over **10–130 Hz** in its Min/Typ/Max columns. The text we could read does not show unambiguously which column holds 8; 15 is in the last column.
  - The HEDscan web page states below 15 fT/√Hz, without a band.

**Relation to the model.**

- **The model's figure.** The model's 15 fT/√Hz is a white, frequency-independent level applied over the whole 1–40 Hz analysis band (`configs/g2_adult.toml`, `[band]`; register A-OPM-NOISE).
- **How it compares with the specifications.**
  - It equals QuSpin's dual-axis specification limit (below 15 fT/√Hz over 3–100 Hz) and FieldLine's upper figure (15 fT/√Hz over 10–130 Hz).
  - Both specifications start at 3 or 10 Hz. Neither describes 1–3 Hz or 1–10 Hz, where real OPM noise rises; the coloured-noise sensitivity analysis is the relevant check there.
- **The clinical figures.**
  - Where they can be traced, they are manufacturer specifications. Hillebrand et al.'s comes from Osborne et al. (2018). Feys et al.'s (2023) 23 fT/√Hz matches QuSpin's triaxial specification.
  - The only noise level reported from a clinical recording is Feys et al.'s (2023) 60–65 fT/√Hz for their ⁴He sensors, with no band.

---

## Sources searched and read

- **Databases and search** (2026-10-05):
  - PubMed (E-utilities, including the competing-interest field);
  - Europe PMC (records, full-text XML and one full-text search);
  - Crossref, OpenAlex, Semantic Scholar, DOAJ;
  - the bioRxiv API and the Internet Archive capture index;
  - web search.
- **Full texts read in PMC**:
  - PMC5562927, PMC6063354, PMC6988110, PMC8274815, PMC9135302;
  - PMC10303922, PMC10715393, PMC10748329, PMC7085997, PMC10030968;
  - PMC12514395, PMC13525599, PMC7615061;
  - PMC5944911, PMC5539206, PMC13081696.
- **Manufacturer documents**:
  - QuSpin: product pages, the archived Gen-2 page, the 2021 Gen-3 announcement, the technical note, and the copy of Osborne et al. (2018) on quspin.com;
  - FieldLine: the HEDscan page and the 2023 specification sheet;
  - Elekta/MEGIN: the TRIUX datasheet, the TRIUX technical manual and the VectorView technical manual, in copies hosted by NatMEG and the MRC CBU.
- **Not accessible or not found**:
  - the full text of Ren et al. (2025): the publisher's site refused automated access;
  - bioRxiv pages: rate-limited, so the preprint was not re-read;
  - a TRIUX neo datasheet;
  - an X-Magtech specification with a band;
  - FieldLine's cell size;
  - any peer-reviewed measurement of Neuromag white noise, in the searches made.

## Open points

1. **Ren et al. (2025).** The article is open access. Reading it in an ordinary browser would settle the systems, the SNR and amplitude definitions, the channels, the sensor noise and any regional detail; `epilepsy_opm_studies.md` §3.10 would then move from abstract to full text.
2. **Jas et al. (2026).** Check again before resubmission for a version of record or a data or code release.
3. **The 10-mm cell.** The manuscript, the register (J-opm-cell, A-OPM-CELL) and `jas2026.md` §5.2 describe it as the device's cell. It should become a modelling choice (see Section 4).
4. **The TRIUX values.** Once the datasheet is cited, the register rows HW-mag-noise, HW-grad-noise, HW-18mm and U-HW1 can change from "transcribed" to "manufacturer datasheet".
