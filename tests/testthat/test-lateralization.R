library(testthat)
library(PhysioEEG)

# 3D MI fixture: mu amplitude is full in the baseline (first 25%) and reduced
# in the task window by erd_left / erd_right, producing a controlled ERD on
# C3 (left) and C4 (right).
make_mi_pe <- function(erd_left, erd_right, n_trials = 20, sr = 256,
                       trial_sec = 4, seed = 1) {
  set.seed(seed)
  labels <- c("C3", "Cz", "C4", "P3", "Pz", "P4", "O1", "O2")
  nt <- as.integer(trial_sec * sr); t <- (0:(nt - 1)) / sr
  nbase <- as.integer(floor(0.25 * nt))
  env <- function(erd) { e <- rep(1, nt); e[(nbase + 1):nt] <- 1 - erd; e }
  arr <- array(0, dim = c(nt, length(labels), n_trials))
  for (tr in seq_len(n_trials)) {
    for (ci in seq_along(labels)) {
      erd <- if (labels[ci] == "C3") erd_left else
             if (labels[ci] == "C4") erd_right else 0
      mu <- 8 * env(erd) * sin(2 * pi * 10 * t + stats::runif(1, 0, 2 * pi))
      arr[, ci, tr] <- stats::rnorm(nt, sd = 1) + mu
    }
  }
  PhysioExperiment(
    assays = list(raw = arr),
    colData = S4Vectors::DataFrame(label = labels,
                                   type = rep("EEG", length(labels))),
    samplingRate = sr)
}

test_that("stronger right-hemisphere ERD yields positive LI (erd method)", {
  res <- eegLateralization(make_mi_pe(erd_left = 0.3, erd_right = 0.7),
                           method = "erd")
  expect_gt(biomarkerValue(res$summary), 0.1)      # right desync more -> positive
  res2 <- eegLateralization(make_mi_pe(erd_left = 0.7, erd_right = 0.3),
                            method = "erd")
  expect_lt(biomarkerValue(res2$summary), -0.1)    # left desync more -> negative
})

test_that("LI is bounded in -1..1 and ~0 for symmetric input", {
  pe <- make_mi_pe(erd_left = 0.5, erd_right = 0.5)
  res <- eegLateralization(pe, method = "erd")
  expect_true(all(res$per_trial >= -1 & res$per_trial <= 1))
  expect_lt(abs(biomarkerValue(res$summary)), 0.1)
  res_p <- eegLateralization(pe, method = "power")
  expect_true(all(res_p$per_trial >= -1 & res_p$per_trial <= 1))
  expect_lt(abs(biomarkerValue(res_p$summary)), 0.1)
})

test_that("power method reflects raw band-power asymmetry", {
  set.seed(2); labels <- c("C3", "Cz", "C4"); nt <- 1024L; sr <- 256
  t <- (0:(nt - 1)) / sr
  arr <- array(0, dim = c(nt, 3, 20))
  for (tr in 1:20) for (ci in 1:3) {
    amp <- if (labels[ci] == "C4") 12 else if (labels[ci] == "C3") 4 else 8
    arr[, ci, tr] <- stats::rnorm(nt, sd = 1) +
      amp * sin(2 * pi * 10 * t + stats::runif(1, 0, 2 * pi))
  }
  pe <- PhysioExperiment(
    assays = list(raw = arr),
    colData = S4Vectors::DataFrame(label = labels, type = rep("EEG", 3)),
    samplingRate = sr)
  expect_gt(biomarkerValue(eegLateralization(pe, method = "power")$summary), 0.3)
})

test_that("summary is a PhysioBiomarker with CI + provenance; per_trial is a matrix", {
  res <- eegLateralization(make_mi_pe(0.3, 0.7, n_trials = 20), method = "erd")
  expect_s4_class(res$summary, "PhysioBiomarker")
  expect_equal(res$summary@name, "LI")
  expect_length(res$summary@ci, 2L)
  expect_equal(res$summary@provenance_info$method, "erd")
  expect_match(res$summary@provenance_info$band, "Hz")
  expect_true(is.matrix(res$per_trial))
  expect_equal(dim(res$per_trial), c(20L, 1L))
})

test_that("arbitrary ROIs (affected vs unaffected hemisphere) are supported", {
  pe <- make_mi_pe(erd_left = 0.3, erd_right = 0.7)
  res <- eegLateralization(pe, left_ch = c("C3", "P3"),
                           right_ch = c("C4", "P4"), method = "erd")
  expect_true(all(res$per_trial >= -1 & res$per_trial <= 1))
  expect_s4_class(res$summary, "PhysioBiomarker")
})

test_that("eegLateralization validates inputs", {
  pe <- make_mi_pe(0.5, 0.5)
  expect_error(eegLateralization(pe, band = c(13, 8)), "band")
  expect_error(eegLateralization(pe, left_ch = "ZZ"), "left_ch")
  expect_error(eegLateralization(pe, right_ch = "YY"), "right_ch")
  pe2d <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)   # 2D, not epoched
  expect_error(eegLateralization(pe2d), "3D")
})
