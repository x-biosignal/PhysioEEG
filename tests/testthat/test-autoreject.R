library(testthat)
library(PhysioEEG)

# Fixture: 3D ERP data (time x channels x epochs) in the "raw" assay, with a set
# of epochs corrupted by large multi-channel artifacts (ground-truth bad).
make_ar_pe <- function(seed = 1, n_epochs = 40, n_channels = 19, sr = 250,
                       bad_epochs = c(3, 8, 15, 22, 29, 36),
                       n_bad_ch = 12, amp = 200) {
  set.seed(seed)
  pe <- make_eeg_erp(n_epochs = n_epochs, n_channels = n_channels, sr = sr,
                     epoch_sec = 1.0)
  pe <- eegMontage(pe, system = "10-20")
  D <- SummarizedExperiment::assay(pe, "raw")
  Tt <- dim(D)[1]; C <- dim(D)[2]
  set.seed(seed + 1)
  for (e in bad_epochs) {
    chs <- sample(seq_len(C), n_bad_ch)
    for (ch in chs) {
      D[, ch, e] <- D[, ch, e] + amp * sin(2 * pi * 3 * (seq_len(Tt) / sr)) +
        stats::rnorm(Tt, sd = 3)
    }
  }
  SummarizedExperiment::assay(pe, "raw") <- D
  list(pe = pe, bad_epochs = bad_epochs, D = D)
}

.interleaved_folds <- function(E, k) ((seq_len(E) - 1L) %% k) + 1L

.heldout_rms <- function(D, ptp, tau, folds) {
  Tt <- dim(D)[1]; C <- dim(D)[2]; k <- max(folds)
  mean(vapply(seq_len(k), function(f) {
    test <- which(folds == f); train <- which(folds != f)
    mean(vapply(seq_len(C), function(c) {
      good <- train[ptp[c, train] <= tau[c]]
      if (length(good) < 1) return(NA_real_)
      mu <- rowMeans(matrix(D[, c, good], nrow = Tt))
      med <- apply(matrix(D[, c, test], nrow = Tt), 1, stats::median)
      sqrt(mean((mu - med)^2))
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1)))
}

test_that("autoreject flags injected high-amplitude epochs (>90% precision/recall)", {
  d <- make_ar_pe()
  pe <- eegAutoReject(d$pe, assay_name = "raw")
  log <- S4Vectors::metadata(pe)$autoreject

  flagged <- sort(union(log$dropped_epochs, unique(log$interpolated$epoch)))
  tp <- length(intersect(flagged, d$bad_epochs))
  precision <- tp / max(1, length(flagged))
  recall <- tp / length(d$bad_epochs)
  expect_gt(precision, 0.9)
  expect_gt(recall, 0.9)
  expect_equal(log$n_total, 40L)
})

test_that("cleaned 3D data drops rejected epochs and keeps channel count", {
  d <- make_ar_pe()
  pe <- eegAutoReject(d$pe, assay_name = "raw", output_assay = "clean3d")
  cleaned <- S4Vectors::metadata(pe)$clean3d
  ar <- S4Vectors::metadata(pe)$autoreject
  expect_equal(length(dim(cleaned)), 3L)
  expect_equal(dim(cleaned)[2], 19L)                 # channels preserved
  expect_equal(dim(cleaned)[3], ar$n_retained)       # epochs = retained count
  expect_equal(ar$n_total, dim(cleaned)[3] + length(ar$dropped_epochs))
})

test_that("rejection log and provenance are recorded", {
  d <- make_ar_pe()
  pe <- eegAutoReject(d$pe, assay_name = "raw")
  ar <- S4Vectors::metadata(pe)$autoreject
  expect_true(all(c("thresholds", "consensus", "n_interpolate", "bad_matrix",
                    "interpolated", "dropped_epochs") %in% names(ar)))
  expect_equal(dim(ar$bad_matrix), c(40L, 19L))      # epoch x channel
  expect_equal(length(ar$thresholds), 19L)
  expect_true(is.data.frame(ar$interpolated))
  expect_true(all(c("epoch", "channel") %in% names(ar$interpolated)))
  expect_true(any(PhysioCore::provenance(pe)$step == "eegAutoReject"))
})

test_that("chosen thresholds increase monotonically with signal scale", {
  set.seed(5)
  pe0 <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250)
  pe0 <- eegMontage(pe0, system = "10-20")
  D0 <- SummarizedExperiment::assay(pe0, "raw")
  med_tau <- vapply(c(1, 2, 3), function(a) {
    pe <- pe0
    SummarizedExperiment::assay(pe, "raw") <- a * D0
    stats::median(S4Vectors::metadata(
      eegAutoReject(pe, assay_name = "raw"))$autoreject$thresholds)
  }, numeric(1))
  expect_true(all(diff(med_tau) > 0))
})

test_that("cross-validated thresholds beat a fixed no-reject baseline (held-out RMS)", {
  d <- make_ar_pe()
  pe <- eegAutoReject(d$pe, assay_name = "raw")
  tau_cv <- S4Vectors::metadata(pe)$autoreject$thresholds
  ptp <- PhysioEEG:::.autoreject_ptp(d$D)
  folds <- .interleaved_folds(dim(d$D)[3], 5)
  tau_fixed <- rep(max(ptp) * 1.5, dim(d$D)[2])       # never rejects
  rms_cv <- .heldout_rms(d$D, ptp, tau_cv, folds)
  rms_fixed <- .heldout_rms(d$D, ptp, tau_fixed, folds)
  expect_lt(rms_cv, rms_fixed)
})

test_that("interpolation branch repairs few-channel artifacts without dropping", {
  # single-channel artifacts stay below any consensus, so they are interpolated
  d <- make_ar_pe(seed = 3, bad_epochs = c(5, 12, 20, 28),
                  n_bad_ch = 1, amp = 250)
  pe <- eegAutoReject(d$pe, assay_name = "raw", consensus = 0.5, n_interpolate = 1)
  ar <- S4Vectors::metadata(pe)$autoreject
  # the corrupted epochs are repaired, not dropped
  expect_true(all(d$bad_epochs %in% ar$interpolated$epoch))
  expect_length(intersect(ar$dropped_epochs, d$bad_epochs), 0)
})

test_that(".autoreject_interp matches eegInterpolate (same spherical-spline core)", {
  d <- make_ar_pe(seed = 3, bad_epochs = 5, n_bad_ch = 1, amp = 250)
  D <- d$D
  slice <- D[, , 5]                                    # time x channels
  cd <- SummarizedExperiment::colData(d$pe)
  pos <- as.matrix(cd[, c("pos_x", "pos_y", "pos_z")])
  bad_idx <- which.max(apply(slice, 2, function(v) max(v) - min(v)))

  # our path
  repaired <- PhysioEEG:::.autoreject_interp(slice, bad_idx, pos)

  # eegInterpolate path on a 2D PE built from the same slice
  pe2d <- PhysioExperiment(
    assays = list(raw = slice), colData = cd,
    samplingRate = samplingRate(d$pe))
  pe2d <- eegInterpolate(pe2d, bad_channels = as.character(cd$label)[bad_idx],
                         method = "spline", assay_name = "raw")
  ref <- SummarizedExperiment::assay(pe2d, "interpolated")

  expect_equal(repaired[, bad_idx], ref[, bad_idx], tolerance = 1e-8)
})

test_that("user-supplied consensus and n_interpolate are honored", {
  d <- make_ar_pe()
  pe <- eegAutoReject(d$pe, assay_name = "raw", consensus = 0.3, n_interpolate = 2)
  ar <- S4Vectors::metadata(pe)$autoreject
  expect_equal(ar$consensus, 0.3)
  expect_equal(ar$n_interpolate, 2L)
})

test_that("eegAutoReject validates inputs", {
  d <- make_ar_pe()
  # no positions -> error
  pe_nopos <- make_eeg_erp(n_epochs = 10, n_channels = 19, sr = 250)
  expect_error(eegAutoReject(pe_nopos, assay_name = "raw"), "positions")
  # no 3D data under the name -> error
  expect_error(eegAutoReject(d$pe, assay_name = "nonexistent"), "3D epoched")
  # bad cv_folds
  expect_error(eegAutoReject(d$pe, assay_name = "raw", cv_folds = 1), "cv_folds")
  # bad consensus
  expect_error(eegAutoReject(d$pe, assay_name = "raw", consensus = 2), "consensus")
})
