library(testthat)
library(PhysioEEG)

# VAL-09: ground-truth validation of phase-based connectivity.
#
# The BSPC/JNSM manuscripts previously reported PLV=0.9996 / wPLI=0.9994 as
# "MNE-Python parity". Those are not reproducible: PhysioEEG estimates phase with
# a continuous Hilbert transform while MNE's spectral_connectivity_epochs uses an
# epoch-wise spectral estimator, so the two agree only structurally (r~0.95).
# Phase connectivity is therefore validated against ANALYTIC ground truth here.

.mk_eeg <- function(M, fs = 250) {
  colnames(M) <- paste0("C", seq_len(ncol(M)))
  PhysioCore::PhysioExperiment(
    assays = list(raw = M),
    colData = S4Vectors::DataFrame(label = colnames(M), type = "eeg"),
    samplingRate = fs)
}

test_that("eegPLV = 1 for phase-locked oscillators and drops when decoupled", {
  set.seed(20260806); fs <- 250; t <- seq(0, 60 - 1/fs, 1/fs); f0 <- 10
  phi <- 2 * pi * f0 * t
  locked   <- .mk_eeg(cbind(cos(phi), cos(phi + pi/4)), fs)   # constant offset
  decoupled <- .mk_eeg(cbind(cos(phi), cos(2*pi*7.3*t + cumsum(rnorm(length(t), 0, 0.2)))), fs)

  plv_locked   <- eegPLV(locked,   band = c(8, 13))$plv[1]
  plv_decoupled <- eegPLV(decoupled, band = c(8, 13))$plv[1]

  expect_gt(plv_locked, 0.99)              # phase-locked -> PLV ~ 1
  expect_lt(plv_decoupled, plv_locked)     # decoupling reduces PLV
})

test_that("eegWPLI rejects zero-lag volume conduction and detects lagged coupling", {
  set.seed(20260806); fs <- 250; t <- seq(0, 60 - 1/fs, 1/fs); f0 <- 10
  src <- sin(2 * pi * f0 * t)
  # zero-lag: same source + independent sensor noise (mimics volume conduction)
  vc  <- .mk_eeg(cbind(src + rnorm(length(t), 0, 0.3), src + rnorm(length(t), 0, 0.3)), fs)
  # genuine coupling with a quarter-cycle lag
  lag <- .mk_eeg(cbind(sin(2*pi*f0*t) + rnorm(length(t), 0, 0.3),
                       sin(2*pi*f0*t - pi/2) + rnorm(length(t), 0, 0.3)), fs)

  wpli_vc  <- eegWPLI(vc,  band = c(8, 13))$wpli[1]
  wpli_lag <- eegWPLI(lag, band = c(8, 13))$wpli[1]

  expect_lt(wpli_vc, 0.2)                   # zero-lag rejected
  expect_gt(wpli_lag, 0.8)                  # lagged coupling detected
  expect_gt(wpli_lag, wpli_vc)
})
