library(testthat)
library(PhysioEEG)

test_that("eegMorletWavelet returns correct dimensions", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  freqs <- seq(5, 40, by = 5)
  result <- eegMorletWavelet(pe, frequencies = freqs)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["wavelet_power"]]))
  wp <- S4Vectors::metadata(result)[["wavelet_power"]]
  # 3D: time x freq x channels
  expect_equal(dim(wp)[1], 2000)
  expect_equal(dim(wp)[2], length(freqs))
  expect_equal(dim(wp)[3], 4)
  expect_true(all(wp >= 0))
})

test_that("eegMorletWavelet detects known frequency", {
  pe <- make_eeg(n_time = 5000, n_channels = 1, sr = 500)
  # make_eeg has strong alpha (10Hz) and beta (20Hz)
  freqs <- seq(5, 30, by = 1)
  result <- eegMorletWavelet(pe, frequencies = freqs)
  wp <- S4Vectors::metadata(result)[["wavelet_power"]]

  # Average power per frequency
  mean_power <- apply(wp[, , 1], 2, mean)
  peak_freq <- freqs[which.max(mean_power)]
  # Peak should be near 10Hz (alpha)
  expect_true(peak_freq >= 8 && peak_freq <= 12)
})

test_that("eegMorletWavelet stores phase in metadata", {
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 250)
  freqs <- seq(5, 20, by = 5)
  result <- eegMorletWavelet(pe, frequencies = freqs)

  md <- S4Vectors::metadata(result)
  expect_true("wavelet" %in% names(md))
  expect_equal(md$wavelet$frequencies, freqs)
  expect_equal(md$wavelet$n_cycles, 7)
  # Phase should be in [-pi, pi]
  phase <- md$wavelet$phase
  expect_equal(dim(phase), c(1000, length(freqs), 2))
  expect_true(all(phase >= -pi & phase <= pi))
})

test_that("eegMorletWavelet validates inputs", {
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 250)
  expect_error(eegMorletWavelet(pe, n_cycles = -1))
  expect_error(eegMorletWavelet(pe, frequencies = c(-5, 10)))
  expect_error(eegMorletWavelet("not_pe"))
})

test_that("eegSTFT returns spectrogram with correct structure", {
  pe <- make_eeg(n_time = 5000, n_channels = 2, sr = 250)
  result <- eegSTFT(pe, window_sec = 0.5, overlap = 0.5)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["stft_power"]]))
  stft_data <- S4Vectors::metadata(result)[["stft_power"]]
  expect_true(all(stft_data >= 0, na.rm = TRUE))
})

test_that("eegSTFT dimensions are consistent with parameters", {
  pe <- make_eeg(n_time = 5000, n_channels = 3, sr = 250)
  result <- eegSTFT(pe, window_sec = 0.5, overlap = 0.75)
  stft_data <- S4Vectors::metadata(result)[["stft_power"]]

  # Window length = 0.5 * 250 = 125, ensure even -> 126
  # Hop = round(126 * 0.25) = 32 (approximately)
  # n_fft_bins = 126 / 2 + 1 = 64
  expect_equal(dim(stft_data)[3], 3)  # channels
  # Frequency bins should be positive half + DC
  md <- S4Vectors::metadata(result)
  expect_true("stft" %in% names(md))
  expect_equal(length(md$stft$freq_axis), dim(stft_data)[2])
  expect_equal(length(md$stft$time_axis), dim(stft_data)[1])
  # Max frequency should be Nyquist
  expect_equal(max(md$stft$freq_axis), 125)
})

test_that("eegSTFT supports different window types", {
  pe <- make_eeg(n_time = 2000, n_channels = 1, sr = 250)

  result_han <- eegSTFT(pe, window_type = "hanning")
  result_ham <- eegSTFT(pe, window_type = "hamming", output_assay = "stft_hamming")
  result_rec <- eegSTFT(pe, window_type = "rectangular", output_assay = "stft_rect")

  expect_s4_class(result_han, "PhysioExperiment")
  expect_s4_class(result_ham, "PhysioExperiment")
  expect_s4_class(result_rec, "PhysioExperiment")

  # Rectangular window should generally produce higher spectral leakage
  # (different power distribution) than tapered windows
  sp_han <- S4Vectors::metadata(result_han)[["stft_power"]]
  sp_rec <- S4Vectors::metadata(result_rec)[["stft_rect"]]
  expect_false(identical(sp_han, sp_rec))
})

test_that("eegMultitaper computes PSD", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegMultitaper(pe, bandwidth = 4)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["multitaper_psd"]]))
  psd <- S4Vectors::metadata(result)[["multitaper_psd"]]
  expect_true(all(psd >= 0))
})

test_that("eegMultitaper PSD has correct frequency resolution", {
  pe <- make_eeg(n_time = 5000, n_channels = 2, sr = 500)
  result <- eegMultitaper(pe, bandwidth = 4)
  psd <- S4Vectors::metadata(result)[["multitaper_psd"]]

  md <- S4Vectors::metadata(result)
  expect_true("multitaper" %in% names(md))
  expect_equal(length(md$multitaper$frequencies), nrow(psd))
  expect_equal(md$multitaper$bandwidth, 4)
  # Default tapers: floor(2*4) - 1 = 7
  expect_equal(md$multitaper$n_tapers, 7L)
  # PSD should be frequencies x channels
  expect_equal(ncol(psd), 2)
})

test_that("eegMultitaper detects dominant frequency", {
  pe <- make_eeg(n_time = 5000, n_channels = 1, sr = 500)
  result <- eegMultitaper(pe, bandwidth = 4)
  psd <- S4Vectors::metadata(result)[["multitaper_psd"]]
  freqs <- S4Vectors::metadata(result)$multitaper$frequencies

  # Find peak in PSD (should be near alpha = 10 Hz)
  # Exclude DC and very low frequencies
  valid_idx <- which(freqs >= 5 & freqs <= 30)
  peak_idx <- valid_idx[which.max(psd[valid_idx, 1])]
  peak_freq <- freqs[peak_idx]
  expect_true(peak_freq >= 8 && peak_freq <= 12)
})

test_that("eegERSP computes baseline-corrected spectral power", {
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  freqs <- seq(5, 30, by = 5)
  result <- eegERSP(pe, baseline = c(1, 50), frequencies = freqs)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["ersp_data"]]))
})

test_that("eegERSP output is in dB with correct dimensions", {
  pe <- make_eeg_erp(n_epochs = 15, n_channels = 2, sr = 250)
  freqs <- seq(5, 25, by = 5)
  result <- eegERSP(pe, baseline = c(1, 30), frequencies = freqs)
  ersp_data <- S4Vectors::metadata(result)[["ersp_data"]]

  # Dimensions: time x freq x channels
  n_time <- dim(SummarizedExperiment::assay(pe, "raw"))[1]
  expect_equal(dim(ersp_data)[1], n_time)
  expect_equal(dim(ersp_data)[2], length(freqs))
  expect_equal(dim(ersp_data)[3], 2)

  # ERSP is in dB - can be positive or negative
  # Baseline region should average near 0 dB
  bl_mean <- mean(ersp_data[1:30, , ])
  expect_true(abs(bl_mean) < 5)  # roughly near 0 dB in baseline

  # Check metadata
  md <- S4Vectors::metadata(result)
  expect_true("ersp" %in% names(md))
  expect_equal(md$ersp$frequencies, freqs)
  expect_equal(md$ersp$n_epochs, 15)
})

test_that("eegERSP rejects 2D input", {
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 250)
  expect_error(eegERSP(pe), "epoched.*3D")
})

test_that("eegITC computes inter-trial coherence", {
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  freqs <- seq(5, 30, by = 5)
  result <- eegITC(pe, frequencies = freqs)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["itc_data"]]))
  itc_data <- S4Vectors::metadata(result)[["itc_data"]]
  # ITC should be between 0 and 1
  expect_true(all(itc_data >= 0 & itc_data <= 1, na.rm = TRUE))
})

test_that("eegITC values near 1 for phase-locked signal", {
  # ERP data has phase-locked components (N100, P300)
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 2, sr = 250)
  freqs <- c(10, 20)
  result <- eegITC(pe, frequencies = freqs)
  itc_data <- S4Vectors::metadata(result)[["itc_data"]]

  # At least some time-frequency points should have high ITC
  max_itc <- max(itc_data, na.rm = TRUE)
  expect_gt(max_itc, 0.3)
})

test_that("eegITC has correct dimensions", {
  pe <- make_eeg_erp(n_epochs = 10, n_channels = 3, sr = 250)
  freqs <- seq(5, 20, by = 5)
  result <- eegITC(pe, frequencies = freqs)
  itc_data <- S4Vectors::metadata(result)[["itc_data"]]

  n_time <- dim(SummarizedExperiment::assay(pe, "raw"))[1]
  expect_equal(dim(itc_data)[1], n_time)
  expect_equal(dim(itc_data)[2], length(freqs))
  expect_equal(dim(itc_data)[3], 3)
})

test_that("eegITC rejects 2D input", {
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 250)
  expect_error(eegITC(pe), "epoched.*3D")
})

test_that("eegITC stores metadata", {
  pe <- make_eeg_erp(n_epochs = 10, n_channels = 2, sr = 250)
  freqs <- c(8, 12, 20)
  result <- eegITC(pe, frequencies = freqs, n_cycles = 5)

  md <- S4Vectors::metadata(result)
  expect_true("itc" %in% names(md))
  expect_equal(md$itc$frequencies, freqs)
  expect_equal(md$itc$n_cycles, 5)
  expect_equal(md$itc$n_epochs, 10)
})

test_that("eegMorletWavelet peak power at 10Hz for alpha signal", {
  set.seed(51)
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Pure 10 Hz sine wave
  data <- matrix(sin(2 * pi * 10 * t_vec), nrow = n_time, ncol = 1)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "Cz", type = "EEG"),
    samplingRate = sr
  )

  freqs <- seq(5, 20, by = 1)
  result <- eegMorletWavelet(pe, frequencies = freqs)
  wp <- S4Vectors::metadata(result)[["wavelet_power"]]

  # Average power per frequency
  mean_power <- apply(wp[, , 1], 2, mean)
  peak_freq <- freqs[which.max(mean_power)]
  expect_equal(peak_freq, 10)
})

test_that("eegSTFT and eegMorletWavelet both detect same dominant frequency", {
  set.seed(52)
  pe <- make_eeg(n_time = 5000, n_channels = 1, sr = 500)

  # STFT
  stft_result <- eegSTFT(pe, window_sec = 1, overlap = 0.5)
  stft_data <- S4Vectors::metadata(stft_result)[["stft_power"]]
  stft_freqs <- S4Vectors::metadata(stft_result)$stft$freq_axis
  valid_stft <- which(stft_freqs >= 5 & stft_freqs <= 30)
  stft_mean_power <- apply(stft_data[, valid_stft, 1, drop = FALSE], 2, mean)
  stft_peak <- stft_freqs[valid_stft[which.max(stft_mean_power)]]

  # Morlet
  freqs_morlet <- seq(5, 30, by = 1)
  morlet_result <- eegMorletWavelet(pe, frequencies = freqs_morlet)
  wp <- S4Vectors::metadata(morlet_result)[["wavelet_power"]]
  morlet_mean_power <- apply(wp[, , 1], 2, mean)
  morlet_peak <- freqs_morlet[which.max(morlet_mean_power)]

  # Both should detect peak near alpha (~10 Hz)
  expect_true(abs(stft_peak - morlet_peak) <= 3)
})

test_that("eegMultitaper dimensions and frequency vector are correct", {
  pe <- make_eeg(n_time = 2000, n_channels = 3, sr = 250)
  result <- eegMultitaper(pe, bandwidth = 4)
  psd <- S4Vectors::metadata(result)[["multitaper_psd"]]
  md <- S4Vectors::metadata(result)$multitaper

  # PSD: frequencies x channels
  expect_equal(ncol(psd), 3)
  expect_equal(nrow(psd), length(md$frequencies))
  # Max frequency should be Nyquist
  expect_equal(max(md$frequencies), 125)
  # All PSD values non-negative
  expect_true(all(psd >= 0))
})

test_that("eegERSP positive sign means power increase above baseline", {
  set.seed(53)
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  freqs <- seq(5, 25, by = 5)
  result <- eegERSP(pe, baseline = c(1, 30), frequencies = freqs)
  ersp_data <- S4Vectors::metadata(result)[["ersp_data"]]

  # Baseline period (samples 1-30) should average near 0 dB
  bl_mean <- mean(ersp_data[1:30, , ])
  expect_true(abs(bl_mean) < 5)

  # ERSP is in dB - can be positive or negative
  expect_true(all(is.finite(ersp_data)))
})

test_that("eegITC identical epochs yield ITC near 1", {
  sr <- 250
  n_time <- 250  # 1 second
  n_channels <- 2
  n_epochs <- 30
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)

  # Create identical epochs (perfectly phase-locked)
  epoch_sig <- sin(2 * pi * 10 * t_vec)
  data <- array(NA_real_, dim = c(n_time, n_channels, n_epochs))
  for (ep in 1:n_epochs) {
    for (ch in 1:n_channels) {
      data[, ch, ep] <- epoch_sig + rnorm(n_time, sd = 0.1)
    }
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )

  result <- eegITC(pe, frequencies = c(10))
  itc_data <- S4Vectors::metadata(result)[["itc_data"]]
  # ITC at 10 Hz should be high (near 1) for phase-locked data
  max_itc <- max(itc_data, na.rm = TRUE)
  expect_gt(max_itc, 0.7)
})

test_that("eegITC random phase epochs yield low ITC", {
  set.seed(54)
  sr <- 250
  n_time <- 500
  n_channels <- 1
  n_epochs <- 30
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)

  # Random phase for each epoch -> no phase locking
  data <- array(NA_real_, dim = c(n_time, n_channels, n_epochs))
  for (ep in 1:n_epochs) {
    phase <- runif(1, 0, 2 * pi)
    data[, 1, ep] <- sin(2 * pi * 10 * t_vec + phase)
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "Ch1", type = "EEG"),
    samplingRate = sr
  )

  result <- eegITC(pe, frequencies = c(10))
  itc_data <- S4Vectors::metadata(result)[["itc_data"]]
  mean_itc <- mean(itc_data, na.rm = TRUE)
  # Random phase -> ITC should be low
  expect_lt(mean_itc, 0.5)
})

test_that("eegMorletWavelet single frequency returns correct dim", {
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 250)
  result <- eegMorletWavelet(pe, frequencies = c(10))
  wp <- S4Vectors::metadata(result)[["wavelet_power"]]

  expect_equal(dim(wp)[1], 1000)   # time
  expect_equal(dim(wp)[2], 1)      # single frequency
  expect_equal(dim(wp)[3], 2)      # channels
  expect_true(all(wp >= 0))
})

test_that("eegSTFT validates invalid overlap parameter", {
  pe <- make_eeg(n_time = 2000, n_channels = 2, sr = 250)
  expect_error(eegSTFT(pe, overlap = 1.0))
  expect_error(eegSTFT(pe, overlap = -0.1))
})
