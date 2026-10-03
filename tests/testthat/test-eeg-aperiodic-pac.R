library(testthat)
library(PhysioEEG)

mk1 <- function(v, sr = 200) {
  pe <- make_eeg(n_time = length(v), n_channels = 1, sr = sr)
  SummarizedExperiment::assay(pe, defaultAssay(pe)) <- matrix(v, ncol = 1)
  pe
}
mkN <- function(mat, sr = 200) {
  dimnames(mat) <- NULL
  pe <- make_eeg(n_time = nrow(mat), n_channels = ncol(mat), sr = sr)
  SummarizedExperiment::assay(pe, defaultAssay(pe)) <- mat
  pe
}

# ---- aperiodic (delegates to PhysioAnalysis::specparam) ---------------------

test_that("eegAperiodic returns per-channel 1/f exponent + peaks", {
  skip_if_not_installed("PhysioAnalysis")
  set.seed(10); N <- 4000; sr <- 200; t <- (0:(N - 1)) / sr
  base <- as.numeric(stats::filter(rnorm(N), 0.95, method = "recursive"))  # ~1/f
  base <- base / stats::sd(base)
  m <- vapply(1:4, function(i) base + 1.5 * sin(2 * pi * 10 * t) + 0.2 * rnorm(N),
              numeric(N))
  ap <- eegAperiodic(mkN(m, sr), freq_range = c(2, 45))
  expect_s3_class(ap, "eeg_aperiodic")
  expect_equal(length(ap$exponent), 4L)
  expect_true(all(is.finite(ap$aperiodic$exponent)))
  expect_true(all(ap$aperiodic$exponent > 0))            # 1/f exponent positive
  expect_true(any(abs(ap$peaks$CF - 10) < 3))            # ~10 Hz oscillatory peak
  expect_output(print(ap), "eeg_aperiodic")
})

# ---- PAC / comodulogram (delegates to PhysioCrossModal) ---------------------

test_that("eegPAC finds stronger theta-gamma coupling when it is present", {
  skip_if_not_installed("PhysioCrossModal")
  set.seed(11); N <- 6000; sr <- 200; t <- (0:(N - 1)) / sr
  theta <- sin(2 * pi * 6 * t)
  amp <- (1 + cos(2 * pi * 6 * t - pi)) / 2
  coupled <- theta + (0.1 + amp) * sin(2 * pi * 50 * t) + 0.1 * rnorm(N)
  uncoupled <- theta + 0.5 * sin(2 * pi * 50 * t) + 0.1 * rnorm(N)
  pc <- eegPAC(mk1(coupled, sr), phase_band = c(4, 8), amp_band = c(40, 60))
  pu <- eegPAC(mk1(uncoupled, sr), phase_band = c(4, 8), amp_band = c(40, 60))
  expect_true(is.finite(pc$pac))
  expect_gt(pc$pac, pu$pac)
})

test_that("eegComodulogram localises the coupled (phase, amplitude) frequencies", {
  skip_if_not_installed("PhysioCrossModal")
  set.seed(12); N <- 6000; sr <- 200; t <- (0:(N - 1)) / sr
  amp <- (1 + cos(2 * pi * 6 * t - pi)) / 2
  sig <- sin(2 * pi * 6 * t) + (0.1 + amp) * sin(2 * pi * 50 * t) + 0.1 * rnorm(N)
  cm <- eegComodulogram(mk1(sig, sr), phase_freqs = c(4, 6, 8, 10),
                        amp_freqs = c(30, 40, 50, 60, 70))
  expect_true(is.matrix(cm$matrix))
  expect_lt(abs(cm$peak$phase_freq - 6), 3)
  expect_lt(abs(cm$peak$amp_freq - 50), 15)
})

test_that("eegPAC / eegAperiodic error clearly when the engine is absent", {
  # only assert the error path when the package is genuinely missing
  skip_if(requireNamespace("PhysioCrossModal", quietly = TRUE))
  expect_error(eegPAC(mk1(rnorm(1000))), "PhysioCrossModal")
})
