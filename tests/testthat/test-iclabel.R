library(testthat)
library(PhysioEEG)

# Fixture: make_eeg with an injected eye-blink IC (frontal-gradient topography,
# large spiky low-frequency time course) and a 50 Hz line IC.
make_iclabel_pe <- function(seed = 42, n = 6000, sr = 250) {
  set.seed(seed)
  pe <- make_eeg(n_time = n, n_channels = 19, sr = sr)
  X <- SummarizedExperiment::assay(pe, "raw")
  m <- ncol(X)
  ch <- as.character(SummarizedExperiment::colData(pe)$label)
  t <- (0:(n - 1)) / sr

  blink_topo <- numeric(m)
  blink_topo[ch %in% c("Fp1", "Fp2")] <- 10
  blink_topo[ch %in% c("F7", "F3", "Fz", "F4", "F8")] <- 5
  blink_topo[ch %in% c("T3", "C3", "Cz", "C4", "T4")] <- 1
  blink_ts <- numeric(n)
  for (b in sample(80:(n - 80), 30)) {
    d <- stats::dnorm(seq(-40, 40), 0, 10)
    d <- d / max(d)
    L <- min(81L, n - b + 1L)
    blink_ts[b:(b + L - 1)] <- blink_ts[b:(b + L - 1)] + 60 * d[seq_len(L)]
  }
  X <- X + outer(blink_ts, blink_topo)

  line_ts <- 15 * sin(2 * pi * 50 * t)
  X <- X + outer(line_ts, stats::runif(m, 0.5, 1.5))

  SummarizedExperiment::assay(pe, "raw") <- X
  list(pe = pe, blink_ts = blink_ts, line_ts = line_ts, ch = ch, X = X)
}

.classes <- c("brain", "muscle", "eye", "heart", "line_noise",
              "channel_noise", "other")

.find_ic <- function(pe, target) {
  S <- S4Vectors::metadata(pe)[["ica_components"]]
  which.max(vapply(seq_len(ncol(S)),
                   function(k) abs(stats::cor(S[, k], target)), numeric(1)))
}

test_that("eegICLabel returns valid 7-class probabilities that sum to 1", {
  d <- make_iclabel_pe()
  pe <- eegICA(d$pe, n_components = 10, method = "fastica")
  res <- eegICLabel(pe)

  expect_true(all(.classes %in% names(res)))
  expect_equal(nrow(res), 10L)
  P <- as.matrix(res[, .classes])
  expect_true(all(abs(rowSums(P) - 1) < 1e-8))
  expect_true(all(P >= 0 & P <= 1))
  expect_true(all(res$label %in% .classes))
})

test_that("injected eye-blink and line ICs receive the correct argmax label", {
  d <- make_iclabel_pe()
  pe <- eegICA(d$pe, n_components = 10, method = "fastica")
  blink_ic <- .find_ic(pe, d$blink_ts)
  line_ic <- .find_ic(pe, d$line_ts)
  res <- eegICLabel(pe)

  expect_equal(res$label[blink_ic], "eye")
  expect_gt(res$eye[blink_ic], 0.5)
  expect_equal(res$label[line_ic], "line_noise")
  expect_gt(res$line_noise[line_ic], 0.5)
})

test_that("eegICLabelFlag flags the eye component and drives eegICAremove", {
  d <- make_iclabel_pe()
  pe <- eegICA(d$pe, n_components = 10, method = "fastica")
  blink_ic <- .find_ic(pe, d$blink_ts)

  bad <- eegICLabelFlag(pe, prob_threshold = 0.5)
  expect_true(blink_ic %in% bad)
  expect_true("eye" %in% attr(bad, "labels"))

  frontal_idx <- which(d$ch %in% c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8"))
  pe_clean <- eegICAremove(pe, components = bad)
  Xc <- S4Vectors::metadata(pe_clean)$ica_cleaned
  var_before <- mean(apply(d$X[, frontal_idx], 2, stats::var))
  var_after <- mean(apply(Xc[, frontal_idx], 2, stats::var))
  expect_gt(1 - var_after / var_before, 0.70)
})

test_that("eegICLabelFlag validates its arguments", {
  d <- make_iclabel_pe()
  pe <- eegICA(d$pe, n_components = 10, method = "fastica")
  expect_error(eegICLabelFlag(pe, prob_threshold = 1.5), "prob_threshold")
  expect_error(eegICLabelFlag(pe, classes = "not_a_class"), "subset")
})

test_that("eegICAdetect iclabel method delegates and keeps the 2-class contract", {
  d <- make_iclabel_pe()
  pe <- eegICA(d$pe, n_components = 10, method = "fastica")
  det <- eegICAdetect(pe, method = "iclabel")

  expect_named(det, c("component", "type", "method", "score"))
  expect_true(all(det$type %in% c("artifact", "neural")))
  expect_true(all(det$method == "iclabel"))
  expect_true(all(det$score >= 0 & det$score <= 1))

  blink_ic <- .find_ic(pe, d$blink_ts)
  line_ic <- .find_ic(pe, d$line_ts)
  expect_equal(det$type[blink_ic], "artifact")
  expect_equal(det$type[line_ic], "artifact")

  # existing 2-class methods remain intact
  cor_det <- eegICAdetect(pe, method = "correlation")
  expect_true(all(cor_det$type %in% c("artifact", "neural")))
})

test_that("shipped ICLabel weight table loads and matches the built-in fallback", {
  path <- system.file("extdata", "iclabel_weights.csv", package = "PhysioEEG")
  expect_true(file.exists(path))
  csv <- utils::read.csv(path, stringsAsFactors = FALSE)
  builtin <- PhysioEEG:::.iclabel_default_weights()
  expect_equal(csv$term, builtin$term)
  for (cl in .classes) expect_equal(csv[[cl]], builtin[[cl]])
})
