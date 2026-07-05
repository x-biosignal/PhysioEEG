# PhysioEEG <img src="man/figures/logo.png" align="right" height="139" alt="PhysioEEG logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioEEG/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioEEG/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioEEG)](https://CRAN.R-project.org/package=PhysioEEG)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioEEG)](https://x-biosignal.r-universe.dev/PhysioEEG)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**EEG Analysis Functions for PhysioExperiment Objects**

PhysioEEG provides the most comprehensive EEG analysis toolkit available in R, with 65 exported functions covering the full electroencephalography analysis pipeline. Built on top of PhysioCore, it supports preprocessing, independent component analysis (ICA), event-related potential (ERP) extraction and measurement, source localization, microstate analysis, sleep staging, brain-computer interface (BCI) feature extraction, clinical EEG analysis, time-frequency decomposition, connectivity analysis, and publication-quality visualization -- all operating directly on `PhysioExperiment` objects.

## Installation

You can install PhysioEEG from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioEEG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEEG")
```

## Quick Start

```r
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
```

## Features

### Preprocessing Pipeline

A complete, configurable preprocessing pipeline from raw recordings to analysis-ready data:

- `eegPreprocess()` -- full pipeline with configurable steps
- `eegFilter()` -- bandpass, highpass, lowpass, and notch filtering
- `eegMontage()` -- apply standard electrode montages (10-20, 10-10)
- `eegRereference()` -- re-reference to average, linked mastoids, or custom reference
- `eegBadChannels()` -- automatic bad channel detection by correlation, variance, and spectral criteria
- `eegInterpolate()` -- spherical spline interpolation of bad channels
- `eegEpoch()` -- segment continuous data into event-locked epochs
- `eegArtifactReject()` -- reject epochs by amplitude threshold, gradient, or statistical criteria

### Independent Component Analysis (ICA)

Blind source separation for artifact removal and source identification:

- `eegICA()` -- decompose signals using FastICA, Infomax, or JADE algorithms
- `eegICAdetect()` -- automatic classification of artifact components (eye blinks, saccades, muscle, cardiac)
- `eegICAmix()` -- inspect mixing and unmixing matrices
- `eegICAremove()` -- remove selected components and reconstruct clean signals
- `eegPlotICA()` -- visualize component topographies, time courses, and spectra

### Event-Related Potential (ERP) Analysis

End-to-end ERP component extraction, measurement, and statistical testing:

- `eegERPdetect()` -- detect standard components (N100, P300, N400, P600, MMN) with adaptive time windows
- `eegERPmeasure()` -- measure peak amplitude, peak latency, mean amplitude, and area
- `eegERPlatency()` -- fractional area latency and onset latency estimation
- `eegERPbaseline()` -- apply baseline correction with configurable windows
- `eegERPdifference()` -- compute difference waveforms between conditions
- `eegERPtest()` -- point-by-point and cluster-based permutation tests
- `eegERPgrandAverage()` -- compute grand average across subjects or sessions
- `eegPlotERP()` -- plot ERP waveforms with confidence intervals and condition overlays

### Source Localization

Estimate the cortical sources underlying scalp EEG:

- `eegForwardModel()` -- compute forward model with boundary element method (BEM)
- `eegSourceEstimate()` -- distributed source estimation with eLORETA or sLORETA
- `eegBeamformer()` -- LCMV beamformer for focal source localization
- `eegSourcePower()` -- compute source power maps for frequency bands
- `eegPlotSource()` -- visualize source estimates on cortical surfaces

### Microstate Analysis

Characterize the temporal dynamics of global brain states:

- `eegMicrostates()` -- segment EEG into microstates using K-means or AAHC (atomize and agglomerate hierarchical clustering)
- `eegMicrostateBackfit()` -- backfit microstate labels to continuous EEG
- `eegMicrostateSequence()` -- extract microstate transition sequences
- `eegMicrostateStats()` -- compute duration, occurrence, coverage, and transition probabilities

### Sleep Staging and Event Detection

Automated sleep analysis following AASM guidelines:

- `eegSleepStage()` -- automatic sleep staging according to AASM criteria (W, N1, N2, N3, REM)
- `eegSpindleDetect()` -- detect sleep spindles (sigma band bursts)
- `eegKcomplexDetect()` -- detect K-complexes
- `eegSlowWaveDetect()` -- detect slow-wave activity (delta oscillations)
- `eegSleepMetrics()` -- compute sleep efficiency, latency, WASO, and stage percentages
- `eegPlotHypnogram()` -- plot hypnogram with detected events overlay

### Brain-Computer Interface (BCI) Features

Feature extraction and classification for BCI paradigms:

- `eegBCIfeatures()` -- extract multi-domain feature vectors (time, frequency, spatial)
- `eegCSP()` -- common spatial patterns for motor imagery discrimination
- `eegSSVEP()` -- steady-state visually evoked potential detection and classification
- `eegMotorImagery()` -- motor imagery feature extraction (mu/beta ERD/ERS)
- `eegBCIclassify()` -- classify BCI features with LDA, SVM, or random forest

### Clinical EEG

Functions for clinical neurophysiology and QEEG:

- `eegSpikeDetect()` -- detect epileptiform spikes and sharp waves
- `eegQEEG()` -- quantitative EEG analysis (absolute/relative power, peak frequency)
- `eegAsymmetry()` -- compute inter-hemispheric asymmetry indices
- `eegSlowing()` -- detect and quantify EEG slowing (theta/alpha ratio)
- `eegSuppression()` -- detect burst-suppression patterns

### Time-Frequency Analysis

Spectral decomposition with multiple methods:

- `eegMorletWavelet()` -- continuous wavelet transform with Morlet wavelets
- `eegSTFT()` -- short-time Fourier transform
- `eegMultitaper()` -- multitaper spectral estimation (DPSS tapers)
- `eegERSP()` -- event-related spectral perturbation
- `eegITC()` -- inter-trial coherence (phase-locking across trials)
- `eegPlotSpectrogram()` -- plot time-frequency spectrograms with significance masking

### Connectivity Analysis

Functional and effective connectivity between EEG channels:

- `eegCoherence()` -- magnitude-squared coherence
- `eegPLV()` -- phase-locking value
- `eegWPLI()` -- weighted phase lag index (robust to volume conduction)
- `eegGrangerCausality()` -- Granger causality for directed connectivity
- `eegConnectivityMatrix()` -- compute full connectivity matrices for any metric
- `eegPlotConnectivity()` -- visualize connectivity matrices and circular connectograms

### Visualization

Publication-quality plots for every analysis stage:

- `eegPlotSignal()` -- multi-channel signal traces with event markers
- `eegPlotERP()` -- ERP waveforms with confidence intervals
- `eegPlotTopomap()` -- topographic scalp maps with interpolation
- `eegPlotTopomapSeries()` -- series of topomaps across time points
- `eegPlotSpectrogram()` -- time-frequency spectrograms
- `eegPlotConnectivity()` -- connectivity matrices and network graphs
- `eegPlotICA()` -- ICA component topographies and time courses
- `eegPlotSource()` -- cortical source maps
- `eegPlotHypnogram()` -- sleep hypnograms with event annotations

### Simulated Data Generators

Ready-to-use data for testing, demonstration, and teaching:

- `make_eeg()` -- multi-channel EEG with realistic spectral properties
- `make_eeg_erp()` -- epoched data with embedded ERP components (N100, P300)
- `make_eeg_sleep()` -- polysomnography data with sleep stages and spindles
- `make_eeg_bci()` -- motor imagery data with lateralized mu/beta patterns
- `make_eeg_spikes()` -- EEG with embedded epileptiform discharges

## Dependencies

- **R** (>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **S4Vectors**
- **stats**

## PhysioExperiment Ecosystem

PhysioEEG is the EEG analysis layer of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| **PhysioEEG** | EEG analysis (ICA, ERP, source, BCI, sleep) |
| [PhysioEMG](https://github.com/x-biosignal/PhysioEMG) | EMG analysis (synergy, fatigue, onset) |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG and HRV analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui
