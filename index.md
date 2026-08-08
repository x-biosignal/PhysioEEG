# PhysioEEG ![PhysioEEG logo](reference/figures/logo.png)

**EEG Analysis Functions for PhysioExperiment Objects**

PhysioEEG provides 79 exported functions covering the full
electroencephalography analysis pipeline. Built on top of PhysioCore, it
supports preprocessing, independent component analysis (ICA),
event-related potential (ERP) extraction and measurement, source
localization, microstate analysis, sleep staging, brain-computer
interface (BCI) feature extraction, clinical EEG analysis,
time-frequency decomposition, connectivity analysis, and
publication-quality visualization – all operating directly on
`PhysioExperiment` objects.

## Installation

You can install PhysioEEG from
[r-universe](https://x-biosignal.r-universe.dev):

``` r

install.packages("PhysioEEG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEEG")
```

## Quick Start

``` r

library(PhysioEEG)

# Generate simulated EEG data with ERP components
pe <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250)

# Run full preprocessing pipeline (filter, re-reference, artifact rejection)
pe <- eegPreprocess(pe, lowcut = 1, highcut = 40, ref = "average")

# Detect and measure P300 component
erp <- eegERPdetect(pe, component = "P300")
measures <- eegERPmeasure(pe, component = "P300")

# Plot ERP waveform and topographic map
eegPlotERP(pe, channels = c("Fz", "Cz", "Pz"))
eegPlotTopomap(pe, time = 0.35)

# Inspect continuous channel traces and reference-dependent GFP together
continuous <- make_eeg(n_time = 2500, n_channels = 8, sr = 500)
eegPlotButterflyGFP(continuous, channels = c("Fz", "F3", "F4"))
```

## Features

### Preprocessing Pipeline

A complete, configurable preprocessing pipeline from raw recordings to
analysis-ready data:

- [`eegPreprocess()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPreprocess.md)
  – full pipeline with configurable steps
- [`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md)
  – bandpass, highpass, lowpass, and notch filtering
- [`eegMontage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMontage.md)
  – apply standard electrode montages (10-20, 10-10)
- [`eegRereference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegRereference.md)
  – re-reference to average, linked mastoids, or custom reference
- [`eegBadChannels()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBadChannels.md)
  – automatic bad channel detection by correlation, variance, and
  spectral criteria
- [`eegInterpolate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md)
  – spherical spline interpolation of bad channels
- [`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md)
  – segment continuous data into event-locked epochs
- [`eegArtifactReject()`](https://x-biosignal.github.io/PhysioEEG/reference/eegArtifactReject.md)
  – reject epochs by amplitude threshold, gradient, or statistical
  criteria

### Independent Component Analysis (ICA)

Blind source separation for artifact removal and source identification:

- [`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md)
  – decompose signals using FastICA, Infomax, or JADE algorithms
- [`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md)
  – automatic classification of artifact components (eye blinks,
  saccades, muscle, cardiac)
- [`eegICAmix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAmix.md)
  – inspect mixing and unmixing matrices
- [`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md)
  – remove selected components and reconstruct clean signals
- [`eegPlotICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotICA.md)
  – visualize component topographies, time courses, and spectra

### Event-Related Potential (ERP) Analysis

End-to-end ERP component extraction, measurement, and statistical
testing:

- [`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md)
  – detect standard components (N100, P300, N400, P600, MMN) with
  adaptive time windows
- [`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md)
  – measure peak amplitude, peak latency, mean amplitude, and area
- [`eegERPlatency()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPlatency.md)
  – fractional area latency and onset latency estimation
- [`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md)
  – apply baseline correction with configurable windows
- [`eegERPdifference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdifference.md)
  – compute difference waveforms between conditions
- [`eegERPtest()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPtest.md)
  – point-by-point and cluster-based permutation tests
- [`eegERPgrandAverage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPgrandAverage.md)
  – compute grand average across subjects or sessions
- [`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)
  – plot ERP waveforms with confidence intervals and condition overlays

### Source Localization

Estimate the cortical sources underlying scalp EEG:

- [`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md)
  – compute forward model with boundary element method (BEM)
- [`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md)
  – distributed source estimation with eLORETA or sLORETA
- [`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md)
  – LCMV beamformer for focal source localization
- [`eegSourcePower()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourcePower.md)
  – compute source power maps for frequency bands
- [`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md)
  – visualize legacy or explicit two-dimensional source maps
- [`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)
  – project real 3D source coordinates into deterministic sagittal,
  axial, and coronal maximum-intensity views

[`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)
requires finite three-dimensional source coordinates; it never turns a
bare amplitude vector into anatomical positions. Its shared type-8
threshold is computed from absolute amplitude before projection.
Source-estimate orientation components are reduced by root-sum-square
within time and RMS across time. The bundled outline is a schematic
display frame, not MRI registration, patient anatomy, or validation of
the inverse solution.

### Microstate Analysis

Characterize the temporal dynamics of global brain states:

- [`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md)
  – segment EEG into microstates using K-means or AAHC (atomize and
  agglomerate hierarchical clustering)
- [`eegMicrostateBackfit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateBackfit.md)
  – backfit microstate labels to continuous EEG
- [`eegMicrostateSequence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateSequence.md)
  – extract microstate transition sequences
- [`eegMicrostateStats()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateStats.md)
  – compute duration, occurrence, coverage, and transition probabilities

### Sleep Staging and Event Detection

Automated sleep analysis following AASM guidelines:

- [`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md)
  – automatic sleep staging according to AASM criteria (W, N1, N2, N3,
  REM)
- [`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md)
  – detect sleep spindles (sigma band bursts)
- [`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md)
  – detect K-complexes
- [`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md)
  – detect slow-wave activity (delta oscillations)
- [`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)
  – compute sleep efficiency, latency, WASO, and stage percentages
- [`eegPlotHypnogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotHypnogram.md)
  – plot hypnogram with detected events overlay

### Brain-Computer Interface (BCI) Features

Feature extraction and classification for BCI paradigms:

- [`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md)
  – extract multi-domain feature vectors (time, frequency, spatial)
- [`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md)
  – common spatial patterns for motor imagery discrimination
- [`eegSSVEP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSSVEP.md)
  – steady-state visually evoked potential detection and classification
- [`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)
  – motor imagery feature extraction (mu/beta ERD/ERS)
- [`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)
  – classify BCI features with LDA, SVM, or random forest

### Clinical EEG

Functions for clinical neurophysiology and QEEG:

- [`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md)
  – detect epileptiform spikes and sharp waves
- [`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md)
  – quantitative EEG analysis (absolute/relative power, peak frequency)
- [`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md)
  – compute inter-hemispheric asymmetry indices
- [`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md)
  – detect and quantify EEG slowing (theta/alpha ratio)
- [`eegSuppression()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSuppression.md)
  – detect burst-suppression patterns

### Time-Frequency Analysis

Spectral decomposition with multiple methods:

- [`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md)
  – continuous wavelet transform with Morlet wavelets
- [`eegSTFT()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSTFT.md)
  – short-time Fourier transform
- [`eegMultitaper()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMultitaper.md)
  – multitaper spectral estimation (DPSS tapers)
- [`eegERSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERSP.md)
  – event-related spectral perturbation
- [`eegITC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegITC.md)
  – inter-trial coherence (phase-locking across trials)
- [`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)
  – plot validated STFT, Morlet, ERSP, or legacy FFT grids with strictly
  aligned display masks

[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)
keeps numeric power unchanged when a mask is supplied: non-selected bins
are dimmed with display alpha and logical boundaries are drawn without
interpolating the mask. Logical masks and p-value matrices must match
the complete pre-filter time and frequency axes exactly. Cluster masks
require repeated, independent and exchangeable observations; the
correction is for the maximum sign-separated four-neighbour cluster
mass, so a corrected cluster is not bin-wise or mechanistic evidence.

### Connectivity Analysis

Functional and effective connectivity between EEG channels:

- [`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md)
  – magnitude-squared coherence
- [`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md)
  – phase-locking value
- [`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md)
  – weighted phase lag index (robust to volume conduction)
- [`eegGrangerCausality()`](https://x-biosignal.github.io/PhysioEEG/reference/eegGrangerCausality.md)
  – Granger causality for directed connectivity
- [`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)
  – compute full connectivity matrices for any metric
- [`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)
  – legacy connectivity heatmaps and circle plots
- [`eegPlotConnectogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectogram.md)
  – deterministic directed/undirected circular views with anatomical or
  cluster ordering, module arcs, and optional edge bundling. Its strict
  absolute threshold is descriptive and does not imply corrected
  statistical significance.

### Visualization

Publication-quality plots for every analysis stage:

- [`eegPlotSignal()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSignal.md)
  – multi-channel signal traces with event markers
- [`eegPlotButterflyGFP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotButterflyGFP.md)
  – exact-order butterfly traces with a shared-time Global Field Power
  panel (population or sample denominator)
- [`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)
  – ERP waveforms with confidence intervals
- [`eegPlotTopomap()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomap.md)
  – topographic scalp maps with interpolation
- [`eegPlotTopomapSeries()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomapSeries.md)
  – series of topomaps across time points
- [`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)
  – masked time-frequency grids with strict axis alignment and
  non-destructive opacity
- [`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)
  – connectivity matrices and network graphs
- [`eegPlotConnectogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectogram.md)
  – inspectable circular connectivity views preserving matrix direction,
  signed values, channel identity, and estimator metadata
- [`eegPlotICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotICA.md)
  – ICA component topographies and time courses
- [`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md)
  – legacy and explicit two-dimensional source maps
- [`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)
  – three-view source maximum-intensity projections
- [`eegPlotHypnogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotHypnogram.md)
  – sleep hypnograms with event annotations

### Simulated Data Generators

Ready-to-use data for testing, demonstration, and teaching:

- [`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md)
  – multi-channel EEG with realistic spectral properties
- [`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md)
  – epoched data with embedded ERP components (N100, P300)
- [`make_eeg_sleep()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_sleep.md)
  – polysomnography data with sleep stages and spindles
- [`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md)
  – motor imagery data with lateralized mu/beta patterns
- [`make_eeg_spikes()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_spikes.md)
  – EEG with embedded epileptiform discharges

## Dependencies

- **R** (\>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **S4Vectors**
- **stats**

## PhysioExperiment Ecosystem

PhysioEEG is the EEG analysis layer of the PhysioExperiment ecosystem, a
suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|----|----|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| **PhysioEEG** | EEG analysis (ICA, ERP, source, BCI, sleep) |
| [PhysioEMG](https://github.com/x-biosignal/PhysioEMG) | EMG analysis (synergy, fatigue, onset) |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG and HRV analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to
browse all available packages.

## License

MIT License. See
[LICENSE](https://x-biosignal.github.io/PhysioEEG/LICENSE) for details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev).
Community and policy documents live in the umbrella repository:

- [Code of
  Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
