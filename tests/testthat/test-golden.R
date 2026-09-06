# Golden regression tests for PhysioEEG.
#
# Each test rebuilds the SAME deterministic input as data-raw/golden.R, runs
# the PACKAGE function, and compares against the stored golden (captured from an
# INDEPENDENT reference in the generator). Goldens live in _golden/ and are
# skipped (not failed) when absent, so a fresh checkout stays green.

library(testthat)

# ---------------------------------------------------------------------------
# KERNEL 1 - ERP averaging
# ---------------------------------------------------------------------------

test_that("eegERPgrandAverage matches base R elementwise mean", {
  set.seed(123); pe1 <- make_eeg_erp(n_epochs = 10, n_channels = 3, sr = 250)
  set.seed(124); pe2 <- make_eeg_erp(n_epochs = 10, n_channels = 3, sr = 250)
  set.seed(125); pe3 <- make_eeg_erp(n_epochs = 10, n_channels = 3, sr = 250)

  res <- eegERPgrandAverage(pe1, pe2, pe3)
  ga <- S4Vectors::metadata(res)$grand_average

  expect_equal_golden(ga, "erp_grand_average", tol = 1e-8)
})

test_that("eegERPmeasure(mean) matches base R colMeans of across-epoch average", {
  set.seed(200); pe <- make_eeg_erp(n_epochs = 12, n_channels = 3, sr = 250)

  res <- eegERPmeasure(pe, window = c(0, 1000), method = "mean",
                       polarity = "positive", epoch_start = 0)

  expect_equal_golden(unname(res$amplitude), "erp_measure_mean", tol = 1e-8)
})

test_that("eegERPmeasure(peak) amplitude+latency match base R peak search", {
  set.seed(200); pe <- make_eeg_erp(n_epochs = 12, n_channels = 3, sr = 250)

  res <- eegERPmeasure(pe, window = c(250, 500), method = "peak",
                       polarity = "positive", epoch_start = 0)

  expect_equal_golden(unname(res$amplitude), "erp_measure_peak_amp", tol = 1e-8)
  expect_equal_golden(unname(res$latency_ms), "erp_measure_peak_lat", tol = 1e-8)
})

# ---------------------------------------------------------------------------
# KERNEL 2 - Microstate GFP
# ---------------------------------------------------------------------------

test_that("eegMicrostates GFP matches analytic spatial standard deviation", {
  set.seed(300); pe <- make_eeg(n_time = 500, n_channels = 8, sr = 250)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  gfp <- S4Vectors::metadata(pe)$microstates$gfp

  expect_equal_golden(unname(gfp), "microstate_gfp", tol = 1e-8)
})

# ---------------------------------------------------------------------------
# KERNEL 3 - Multitaper PSD
# ---------------------------------------------------------------------------

test_that("eegMultitaper frequency axis matches analytic one-sided FFT axis", {
  sr <- 250; n_time <- 512
  tt <- seq(0, (n_time - 1) / sr, length.out = n_time)
  set.seed(400)
  sig <- 2 * sin(2 * pi * 10 * tt) + 0.01 * rnorm(n_time)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(sig, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "Ch1", type = "EEG"),
    samplingRate = sr
  )
  pe <- eegMultitaper(pe, bandwidth = 3, n_tapers = 3)
  freqs <- S4Vectors::metadata(pe)$multitaper$frequencies

  expect_equal_golden(freqs, "multitaper_freq_axis", tol = 1e-10)
})

test_that("eegMultitaper PSD is stable (characterization guard)", {
  sr <- 250; n_time <- 512
  tt <- seq(0, (n_time - 1) / sr, length.out = n_time)
  set.seed(400)
  sig <- 2 * sin(2 * pi * 10 * tt) + 0.01 * rnorm(n_time)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(sig, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "Ch1", type = "EEG"),
    samplingRate = sr
  )
  pe <- eegMultitaper(pe, bandwidth = 3, n_tapers = 3)
  psd <- S4Vectors::metadata(pe)$multitaper_psd[, 1]

  expect_equal_golden(psd, "multitaper_psd_char", tol = 1e-8)
})
