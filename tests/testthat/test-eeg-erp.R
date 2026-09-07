library(testthat)
library(PhysioEEG)

test_that("eegERPdetect finds P300 in target epochs", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  result <- eegERPdetect(pe, component = "P300")

  expect_s3_class(result, "data.frame")
  expect_true("latency_ms" %in% names(result))
  expect_true("amplitude" %in% names(result))
  # P300 should be detected between 250-500ms
  expect_true(all(result$latency_ms >= 250 & result$latency_ms <= 500))
  # Amplitude should be positive for P300
  expect_true(all(result$amplitude > 0))
})

test_that("eegERPdetect finds N100 component", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  result <- eegERPdetect(pe, component = "N100")

  expect_s3_class(result, "data.frame")
  expect_true(all(result$latency_ms >= 80 & result$latency_ms <= 150))
  expect_true(all(result$amplitude < 0))
})

test_that("eegERPmeasure with peak method returns correct structure", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  result <- eegERPmeasure(pe, window = c(250, 500), method = "peak",
                          polarity = "positive")

  expect_s3_class(result, "data.frame")
  expect_true("amplitude" %in% names(result))
  expect_true("latency_ms" %in% names(result))
  expect_equal(nrow(result), 4)  # one per channel
})

test_that("eegERPmeasure mean method returns mean amplitude", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  result <- eegERPmeasure(pe, window = c(250, 500), method = "mean",
                          polarity = "positive")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
})

test_that("eegERPlatency returns fractional area latency", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  result <- eegERPlatency(pe, window = c(250, 500), fraction = 0.5,
                          polarity = "positive")

  expect_s3_class(result, "data.frame")
  expect_true("latency_ms" %in% names(result))
  # Should be within the window
  expect_true(all(result$latency_ms >= 250 & result$latency_ms <= 500))
})

test_that("eegERPdifference computes correct difference", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  result <- eegERPdifference(pe1, pe2)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["difference"]]))
  diff_data <- S4Vectors::metadata(result)[["difference"]]
  expect_equal(dim(diff_data), dim(SummarizedExperiment::assay(pe1)))
})

test_that("eegERPdetect validates 3D input", {
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)  # 2D, not epoched
  expect_error(eegERPdetect(pe, component = "P300"))
})

test_that("eegERPbaseline corrects baseline to ~zero", {
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  # Assume epoch starts at 0, baseline is first 50ms
  result <- eegERPbaseline(pe, baseline = c(0, 50), epoch_start = 0)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["baseline_corrected"]]))
  bc_data <- S4Vectors::metadata(result)[["baseline_corrected"]]
  # Baseline period mean should be near zero after correction
  baseline_samples <- 1:floor(50 / 1000 * 250)
  for (ch in 1:4) {
    baseline_mean <- mean(bc_data[baseline_samples, ch, ], na.rm = TRUE)
    expect_lt(abs(baseline_mean), 1)
  }
})

test_that("eegERPbaseline validates 3D input", {
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  expect_error(eegERPbaseline(pe, baseline = c(0, 50)))
})

test_that("eegERPtest permutation returns p-values", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  result <- eegERPtest(pe1, pe2, method = "permutation", n_perm = 100)

  expect_s3_class(result, "data.frame")
  expect_true("p_value" %in% names(result))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

test_that("eegERPtest cluster method returns clusters", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  result <- eegERPtest(pe1, pe2, method = "cluster", n_perm = 100)

  expect_s3_class(result, "data.frame")
  expect_true("cluster_id" %in% names(result))
})

test_that("eegERPgrandAverage averages across subjects", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  result <- eegERPgrandAverage(pe1, pe2)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["grand_average"]]))
})

test_that("eegERPgrandAverage validates dimension mismatch", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  expect_error(eegERPgrandAverage(pe1, pe2))
})

test_that("eegERPgrandAverage requires at least 2 objects", {
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  expect_error(eegERPgrandAverage(pe1))
})

test_that("eegERPtest on identical conditions yields p > 0.05", {
  set.seed(42)
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
  # Split same data into two identical halves
  result <- eegERPtest(pe, pe, method = "permutation", n_perm = 100)

  expect_s3_class(result, "data.frame")
  # Most p-values should be non-significant when comparing identical data
  expect_true(mean(result$p_value > 0.05) > 0.5)
})

test_that("eegERPbaseline mean of baseline period is near zero", {
  set.seed(10)
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  result <- eegERPbaseline(pe, baseline = c(0, 50), epoch_start = 0)

  bc_data <- S4Vectors::metadata(result)[["baseline_corrected"]]
  # Baseline period = samples 1 to floor(50/1000 * 250) = samples 1:12
  baseline_samples <- 1:floor(50 / 1000 * 250)
  for (ch in 1:4) {
    bl_mean <- mean(bc_data[baseline_samples, ch, ], na.rm = TRUE)
    expect_lt(abs(bl_mean), 1.0)
  }
})

test_that("eegERPgrandAverage values are mean of inputs", {
  set.seed(20)
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  result <- eegERPgrandAverage(pe1, pe2)

  ga_data <- S4Vectors::metadata(result)[["grand_average"]]
  d1 <- SummarizedExperiment::assay(pe1)
  d2 <- SummarizedExperiment::assay(pe2)

  expected <- (d1 + d2) / 2
  expect_equal(ga_data, expected, tolerance = 1e-10)
})

test_that("eegERPlatency N100 is near 100ms and P300 near 300ms", {
  set.seed(30)
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)

  # N100: negative peak near 100ms
  n100 <- eegERPdetect(pe, component = "N100")
  mean_n100_lat <- mean(n100$latency_ms)
  expect_true(mean_n100_lat >= 70 && mean_n100_lat <= 130)

  # P300: positive peak near 300ms
  p300 <- eegERPdetect(pe, component = "P300")
  mean_p300_lat <- mean(p300$latency_ms)
  expect_true(mean_p300_lat >= 250 && mean_p300_lat <= 350)
})

test_that("eegERPdifference is numerically x minus y", {
  set.seed(40)
  pe1 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  pe2 <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)

  result <- eegERPdifference(pe1, pe2)
  diff_data <- S4Vectors::metadata(result)[["difference"]]
  d1 <- SummarizedExperiment::assay(pe1)
  d2 <- SummarizedExperiment::assay(pe2)

  expected_diff <- d1 - d2
  expect_equal(diff_data, expected_diff, tolerance = 1e-10)
})

test_that("eegERPdetect single epoch works without error", {
  set.seed(50)
  pe <- make_eeg_erp(n_epochs = 1, n_channels = 4, sr = 250)
  result <- eegERPdetect(pe, component = "P300")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_true(all(is.finite(result$latency_ms)))
})

test_that("eegERPmeasure adaptive_mean method works", {
  set.seed(60)
  pe <- make_eeg_erp(n_epochs = 30, n_channels = 4, sr = 250)
  result <- eegERPmeasure(pe, window = c(250, 500), method = "adaptive_mean",
                          polarity = "positive")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_true(all(result$method == "adaptive_mean"))
  expect_true(all(is.finite(result$amplitude)))
})

test_that("eegERPtest target vs standard yields significant differences", {
  set.seed(70)
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 2, sr = 250)
  conds <- metadata(pe)$conditions
  raw <- SummarizedExperiment::assay(pe)

  # Split into target and standard
  target_idx <- which(conds == "target")
  standard_idx <- which(conds == "standard")

  pe_target <- pe
  pe_standard <- pe
  SummarizedExperiment::assays(pe_target) <- list(raw = raw[, , target_idx, drop = FALSE])
  SummarizedExperiment::assays(pe_standard) <- list(raw = raw[, , standard_idx, drop = FALSE])

  result <- eegERPtest(pe_target, pe_standard, method = "permutation", n_perm = 200)
  # There should be at least some significant time points in the P300 window
  p300_window <- which(result$time_sample >= 63 & result$time_sample <= 125)  # ~250-500ms
  expect_true(any(result$p_value[p300_window] < 0.1))
})
