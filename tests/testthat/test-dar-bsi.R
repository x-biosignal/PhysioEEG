library(testthat)
library(PhysioEEG)

LABELS_1020 <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "T3", "C3", "Cz",
                 "C4", "T4", "T5", "P3", "Pz", "P4", "T6", "O1", "O2")

PAIRS <- list(c("Fp2", "Fp1"), c("F4", "F3"), c("F8", "F7"), c("C4", "C3"),
              c("T4", "T3"), c("P4", "P3"), c("T6", "T5"), c("O2", "O1"))

make_pe <- function(mat, sr = 250) {
  colnames(mat) <- NULL
  PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = LABELS_1020[seq_len(ncol(mat))],
                                   type = rep("EEG", ncol(mat))),
    samplingRate = sr)
}

test_that("BSI lies in 0..1; identical hemispheres ~0, fully asymmetric ~1", {
  set.seed(1); sr <- 250; n <- sr * 20; t <- (0:(n - 1)) / sr
  M <- matrix(stats::rnorm(n * 19), n, 19); colnames(M) <- LABELS_1020
  for (p in PAIRS) {                      # identical signal to both members
    s <- 5 * sin(2 * pi * 10 * t) + 3 * sin(2 * pi * 20 * t) + stats::rnorm(n)
    M[, p[1]] <- s; M[, p[2]] <- s
  }
  bsi_sym <- eegBSI(make_pe(M))$bsi
  expect_s4_class(bsi_sym, "PhysioBiomarker")
  v <- biomarkerValue(bsi_sym)
  expect_gte(v, 0); expect_lte(v, 1)
  expect_lt(v, 0.05)                      # symmetric -> ~0

  M2 <- M
  for (p in PAIRS) M2[, p[2]] <- stats::rnorm(n, sd = 1e-6)   # silence left
  vA <- biomarkerValue(eegBSI(make_pe(M2))$bsi)
  expect_gte(vA, 0); expect_lte(vA, 1)
  expect_gt(vA, 0.95)                     # fully asymmetric -> ~1
})

test_that("directed pdBSI carries the sign of the asymmetry", {
  set.seed(2); sr <- 250; n <- sr * 20; t <- (0:(n - 1)) / sr
  M <- matrix(stats::rnorm(n * 19), n, 19); colnames(M) <- LABELS_1020
  base <- 4 * sin(2 * pi * 10 * t) + stats::rnorm(n)
  for (p in PAIRS) { M[, p[1]] <- 2 * base; M[, p[2]] <- base }   # right stronger
  expect_gt(biomarkerValue(eegBSI(make_pe(M), directed = TRUE)$pdbsi), 0.1)
  for (p in PAIRS) { M[, p[1]] <- base; M[, p[2]] <- 2 * base }   # left stronger
  expect_lt(biomarkerValue(eegBSI(make_pe(M), directed = TRUE)$pdbsi), -0.1)
})

test_that("DAR increases monotonically as injected delta power rises", {
  set.seed(3); sr <- 250; n <- sr * 20; t <- (0:(n - 1)) / sr
  dar_at <- function(scale) {
    M <- matrix(0, n, 19)
    for (j in seq_len(19)) {
      M[, j] <- 8 * sin(2 * pi * 10 * t) + scale * 8 * sin(2 * pi * 2 * t) +
        stats::rnorm(n)
    }
    biomarkerValue(eegDAR(make_pe(M))[["C3"]])
  }
  vals <- vapply(c(0.5, 1, 2, 4), dar_at, numeric(1))
  expect_true(all(diff(vals) > 0))
})

test_that("eegDAR returns per-channel biomarkers with provenance + reliability + normative", {
  set.seed(4); sr <- 250
  pe <- make_eeg(n_time = sr * 10, n_channels = 19, sr = sr)
  dar <- eegDAR(pe, age = 60)
  expect_type(dar, "list")
  expect_true(all(names(dar) %in% LABELS_1020))
  bm <- dar[["C3"]]
  expect_s4_class(bm, "PhysioBiomarker")
  expect_equal(bm@name, "DAR")
  # provenance (band, method, version) intact
  expect_equal(bm@provenance_info$method, "welch")
  expect_match(bm@provenance_info$band, "delta")
  expect_false(is.null(bm@provenance_info$software_version))
  # reliability placeholders present
  expect_true(all(c("icc", "sem") %in% names(bm@reliability)))
  # normative reference range attached, and percentile via normativeLookup
  expect_length(bm@reference_range, 2L)
  nl <- PhysioCore::normativeLookup(bm, age = 60)
  expect_true(nl$matched)
  expect_true(is.finite(nl$percentile))
})

test_that("eegBSI biomarkers carry provenance + normative reference range", {
  set.seed(5); sr <- 250
  pe <- make_eeg(n_time = sr * 10, n_channels = 19, sr = sr)
  res <- eegBSI(pe, age = 60)
  expect_named(res, c("bsi", "pdbsi"))
  expect_s4_class(res$bsi, "PhysioBiomarker")
  expect_equal(res$bsi@provenance_info$method, "welch")
  expect_match(res$bsi@provenance_info$band, "Hz")
  expect_length(res$bsi@reference_range, 2L)          # BSI normative row exists
  expect_true(is.finite(PhysioCore::normativeLookup(res$bsi, age = 60)$percentile))
})

test_that("Welch band power is positive and frequency-selective", {
  set.seed(6); sr <- 250; n <- sr * 20; t <- (0:(n - 1)) / sr
  sig <- 5 * sin(2 * pi * 10 * t) + stats::rnorm(n)
  p_alpha <- PhysioEEG:::.compute_band_power(sig, sr, 8, 13)
  p_delta <- PhysioEEG:::.compute_band_power(sig, sr, 1, 4)
  expect_true(is.finite(p_alpha) && p_alpha > 0)
  expect_gt(p_alpha, p_delta)                          # alpha-dominant signal
})

test_that("eegDAR / eegBSI validate inputs", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 250)
  expect_error(eegDAR(pe, delta = c(4, 1)), "delta")
  expect_error(eegBSI(pe, band = c(25, 1)), "band")
  # no matching hemispheric pairs -> error
  pe_nolabel <- make_pe(matrix(stats::rnorm(2000 * 3), 2000, 3))
  expect_error(eegBSI(pe_nolabel, pairs = list(c("ZZ", "YY"))), "matched")
})
