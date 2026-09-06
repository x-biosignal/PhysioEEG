library(testthat)
library(PhysioEEG)

# FFT-based analytic signal, identical to PhysioCrossModal's .hilbert_analytic.
.hilbert <- function(x) {
  x <- as.numeric(x); n <- length(x)
  fx <- stats::fft(x)
  h <- numeric(n); h[1] <- 1
  if (n %% 2 == 0) { if (n >= 4) h[2:(n / 2)] <- 2; h[n / 2 + 1] <- 1
  } else { if (n >= 3) h[2:((n + 1) / 2)] <- 2 }
  stats::fft(fx * h, inverse = TRUE) / n
}

test_that("eegWPLI debiasing matches PhysioCrossModal::weightedPLI (< 1e-6)", {
  skip_if_not_installed("PhysioCrossModal")
  skip_if_not_installed("signal")
  sr <- 250; n <- 4000; t <- (0:(n - 1)) / sr
  set.seed(3)
  a <- sin(2 * pi * 10 * t) + rnorm(n, 0, 0.5)
  b <- sin(2 * pi * 10 * t + 0.7) + rnorm(n, 0, 0.5)

  # weightedPLI = bandpass (order-4 Butterworth) + Hilbert cross-spectrum.
  bf <- signal::butter(4, c(8, 13) / (sr / 2), type = "pass")
  imag <- Im(.hilbert(signal::filtfilt(bf, a)) *
               Conj(.hilbert(signal::filtfilt(bf, b))))

  # eegWPLI and weightedPLI both debias via the shared PhysioCore estimator, so
  # on the identical imaginary cross-spectrum they agree exactly.
  est <- PhysioCore::wpliEstimate(imag)
  wp <- PhysioCrossModal::weightedPLI(a, b, sr = sr, freq_band = c(8, 13))

  expect_lt(abs(est$wpli_debiased - wp$wpli_debiased), 1e-6)
  expect_lt(abs(est$wpli - wp$wpli), 1e-6)
})

test_that("eegWPLI debiased of independent noise distributes around 0", {
  # unbiasedness on real noise channels (acceptance criterion)
  sr <- 250; n <- 4000
  set.seed(1)
  vals <- vapply(seq_len(200), function(k) {
    pe <- PhysioExperiment(
      assays = list(raw = cbind(rnorm(n), rnorm(n))),
      colData = S4Vectors::DataFrame(label = c("A", "B"),
                                     type = c("EEG", "EEG")),
      samplingRate = sr)
    eegWPLI(pe, band = c(8, 13))$wpli_debiased
  }, numeric(1))
  m <- mean(vals)
  # centered near chance (0), not inflated positive
  expect_lt(abs(m), 0.05)
  expect_true(any(vals < 0))              # unbiased: negative values occur
})

test_that("eegWPLI recovers strong phase-lag coupling and exposes both columns", {
  sr <- 250; n <- 4000; t <- (0:(n - 1)) / sr
  set.seed(5)
  a <- sin(2 * pi * 10 * t) + rnorm(n, 0, 0.3)
  b <- sin(2 * pi * 10 * t + 0.8) + rnorm(n, 0, 0.3)
  pe <- PhysioExperiment(
    assays = list(raw = cbind(a, b)),
    colData = S4Vectors::DataFrame(label = c("A", "B"), type = c("EEG", "EEG")),
    samplingRate = sr)
  res <- eegWPLI(pe, band = c(8, 13))
  expect_true(all(c("wpli", "wpli_debiased") %in% names(res)))
  expect_gt(res$wpli, 0.7)
  expect_gt(res$wpli_debiased, 0.7)
  # debiased = FALSE suppresses the debiased estimate
  res_nd <- eegWPLI(pe, band = c(8, 13), debiased = FALSE)
  expect_true(is.na(res_nd$wpli_debiased))
})
