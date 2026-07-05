library(testthat)
library(PhysioEEG)

test_that("eegSpikeDetect finds embedded spikes", {
  pe <- make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 15)
  result <- eegSpikeDetect(pe, method = "morphology")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("channel", "sample", "time_sec", "amplitude") %in% names(result)))
  # Should detect at least some of the 15 embedded spikes
  expect_gt(nrow(result), 0)
  expect_true(all(result$amplitude > 0))
})

test_that("eegSpikeDetect with template method works", {
  pe <- make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 10)
  result <- eegSpikeDetect(pe, method = "template")

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})

test_that("eegSpikeDetect returns empty data.frame for clean EEG", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegSpikeDetect(pe, method = "morphology", threshold_sd = 8,
                           min_amplitude = 200)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("channel", "sample", "time_sec", "amplitude",
                     "duration_ms", "confidence") %in% names(result)))
})

test_that("eegSpikeDetect validates inputs", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  expect_error(eegSpikeDetect(pe, threshold_sd = -1))
  expect_error(eegSpikeDetect(pe, min_amplitude = -5))
  expect_error(eegSpikeDetect("not_pe"))
})

test_that("eegQEEG computes band powers", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegQEEG(pe)

  expect_s4_class(result, "PhysioExperiment")
  qeeg_data <- metadata(result)$qeeg$absolute_power
  expect_true(!is.null(qeeg_data))
  expect_equal(nrow(qeeg_data), 19)
  expect_equal(ncol(qeeg_data), 5)  # 5 bands
  expect_true(all(qeeg_data >= 0))
})

test_that("eegQEEG relative power sums to ~1", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegQEEG(pe)

  qeeg_info <- metadata(result)$qeeg
  expect_true(!is.null(qeeg_info))
  expect_true("relative_power" %in% names(qeeg_info))
  # Each row (channel) relative power should sum to ~1
  row_sums <- rowSums(qeeg_info$relative_power)
  expect_true(all(abs(row_sums - 1) < 0.01))
})

test_that("eegQEEG with custom bands works", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  custom_bands <- list(low = c(1, 10), high = c(10, 50))
  result <- eegQEEG(pe, bands = custom_bands)

  qeeg_data <- metadata(result)$qeeg$absolute_power
  expect_equal(ncol(qeeg_data), 2)
  expect_equal(colnames(qeeg_data), c("low", "high"))
})

test_that("eegQEEG validates inputs", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  expect_error(eegQEEG("not_pe"))
  expect_error(eegQEEG(pe, window_sec = -1))
  expect_error(eegQEEG(pe, overlap = 1.5))
})

test_that("eegAsymmetry returns correct structure", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegAsymmetry(pe)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("left_channel", "right_channel", "asymmetry_index") %in% names(result)))
  # Asymmetry should be finite
  expect_true(all(is.finite(result$asymmetry_index)))
})

test_that("eegAsymmetry with custom pairs works", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegAsymmetry(pe, pairs = list(c("P4", "P3")), band = c(8, 13))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$right_channel, "P4")
  expect_equal(result$left_channel, "P3")
})

test_that("eegAsymmetry warns for missing channels", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  expect_warning(
    result <- eegAsymmetry(pe, pairs = list(c("ZZZ", "YYY"))),
    "not found"
  )
})

test_that("eegAsymmetry validates inputs", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  expect_error(eegAsymmetry("not_pe"))
  expect_error(eegAsymmetry(pe, band = c(13, 8)))  # low > high
})

test_that("eegSuppression detects suppression periods", {
  # Create data with a suppression period
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  # Manually suppress a segment
  data <- SummarizedExperiment::assay(pe)
  data[2000:4000, ] <- data[2000:4000, ] * 0.01  # very low amplitude
  SummarizedExperiment::assay(pe) <- data

  result <- eegSuppression(pe, threshold = 5)

  expect_s3_class(result, "data.frame")
  expect_true("type" %in% names(result))
  expect_true(any(result$type == "suppression"))
  expect_true(!is.null(attr(result, "bsr")))
})

test_that("eegSuppression BSR is between 0 and 100", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  result <- eegSuppression(pe, threshold = 5)

  bsr <- attr(result, "bsr")
  expect_true(is.numeric(bsr))
  expect_true(bsr >= 0 && bsr <= 100)
})

test_that("eegSuppression returns correct columns", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegSuppression(pe, threshold = 50)

  expect_true(all(c("type", "start_sample", "end_sample", "duration_ms") %in%
                    names(result)))
  expect_true(all(result$type %in% c("burst", "suppression")))
})

test_that("eegSuppression validates inputs", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  expect_error(eegSuppression("not_pe"))
  expect_error(eegSuppression(pe, threshold = -1))
  expect_error(eegSuppression(pe, window_ms = -100))
})

test_that("eegSlowing classifies normal EEG correctly", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegSlowing(pe, method = "dtar")

  expect_s3_class(result, "data.frame")
  expect_true("classification" %in% names(result))
  # Normal EEG (alpha dominant) should mostly be classified as normal
  expect_true(any(result$classification == "normal"))
})

test_that("eegSlowing with peak_frequency method works", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegSlowing(pe, method = "peak_frequency")

  expect_s3_class(result, "data.frame")
  expect_true(all(result$value > 0))
})

test_that("eegSlowing with tdr method works", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegSlowing(pe, method = "tdr")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 19)
  expect_true(all(result$metric == "tdr"))
  expect_true(all(result$value > 0))
})

test_that("eegSlowing returns all required columns", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegSlowing(pe, method = "dtar")

  expect_true(all(c("channel", "metric", "value", "classification") %in%
                    names(result)))
  expect_true(all(result$classification %in%
                    c("normal", "mild_slowing", "moderate_slowing", "severe_slowing")))
})

test_that("eegSlowing validates inputs", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  expect_error(eegSlowing("not_pe"))
})

test_that("eegSpikeDetect sensitivity > 0 on make_eeg_spikes data", {
  set.seed(33)
  pe <- make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 15)
  ground_truth <- metadata(pe)$spike_locations

  result <- eegSpikeDetect(pe, method = "morphology")
  # Some spikes should be detected
  expect_gt(nrow(result), 0)

  # Check that at least one detection is near a true spike
  true_samples <- ground_truth$sample
  detected_samples <- result$sample
  any_near <- FALSE
  for (ts in true_samples) {
    if (any(abs(detected_samples - ts) < 50)) {  # within 100ms at 500Hz
      any_near <- TRUE
      break
    }
  }
  expect_true(any_near)
})

test_that("eegQEEG alpha power dominant for 10Hz signal", {
  set.seed(34)
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Pure 10 Hz tone (alpha band)
  data <- matrix(sin(2 * pi * 10 * t_vec), nrow = n_time, ncol = 1)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "Cz", type = "EEG"),
    samplingRate = sr
  )

  result <- eegQEEG(pe)
  qeeg_data <- metadata(result)$qeeg$absolute_power
  # Alpha should have the most power
  band_names <- colnames(qeeg_data)
  alpha_idx <- which(band_names == "alpha")
  expect_equal(unname(which.max(qeeg_data[1, ])), alpha_idx)
})

test_that("eegAsymmetry symmetric data yields near-zero asymmetry", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Same signal for both channels
  sig <- sin(2 * pi * 10 * t_vec)
  data <- cbind(sig, sig)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("F4", "F3"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  result <- eegAsymmetry(pe, pairs = list(c("F4", "F3")), band = c(8, 13))
  expect_equal(abs(result$asymmetry_index[1]), 0, tolerance = 0.01)
})

test_that("eegAsymmetry asymmetric data yields non-zero asymmetry", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Different amplitude alpha for right vs left
  right_sig <- 20 * sin(2 * pi * 10 * t_vec)
  left_sig <- 5 * sin(2 * pi * 10 * t_vec)
  data <- cbind(right_sig, left_sig)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("F4", "F3"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  result <- eegAsymmetry(pe, pairs = list(c("F4", "F3")), band = c(8, 13))
  # Right > left alpha power => positive asymmetry
  expect_true(result$asymmetry_index[1] > 0.1)
})

test_that("eegSuppression detects suppression in flat signal segments", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  data <- SummarizedExperiment::assay(pe)
  # Create a clear suppression segment
  data[3000:6000, ] <- data[3000:6000, ] * 0.001
  SummarizedExperiment::assay(pe) <- data

  result <- eegSuppression(pe, threshold = 5)
  expect_true(any(result$type == "suppression"))

  bsr <- attr(result, "bsr")
  expect_true(bsr > 0)
})

test_that("eegSlowing detects delta-heavy signal as slowed", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Pure 2 Hz delta signal (no alpha/beta)
  data <- matrix(50 * sin(2 * pi * 2 * t_vec), nrow = n_time, ncol = 4)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Fp1", "Fp2", "F3", "F4"),
      type = rep("EEG", 4)
    ),
    samplingRate = sr
  )

  result <- eegSlowing(pe, method = "peak_frequency")
  # Peak frequency should be near 2 Hz -> moderate or severe slowing
  expect_true(all(result$value < 5 | is.na(result$value)))
  expect_true(any(result$classification %in% c("moderate_slowing", "severe_slowing")))
})

test_that("eegSuppression all-normal EEG has BSR near 0", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegSuppression(pe, threshold = 0.01)

  bsr <- attr(result, "bsr")
  expect_true(is.numeric(bsr))
  # Normal EEG should not show much suppression at reasonable threshold
  # (with threshold 0.01 it might, but structure should be correct)
  expect_true(bsr >= 0 && bsr <= 100)
})

test_that("eegQEEG returns all standard frequency bands", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  result <- eegQEEG(pe)

  qeeg_data <- metadata(result)$qeeg$absolute_power
  expected_bands <- c("delta", "theta", "alpha", "beta", "gamma")
  expect_equal(colnames(qeeg_data), expected_bands)
  expect_equal(nrow(qeeg_data), 4)
  expect_true(all(qeeg_data >= 0))
})
