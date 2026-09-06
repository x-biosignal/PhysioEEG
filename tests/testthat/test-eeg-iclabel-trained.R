library(testthat)
library(PhysioEEG)

.iclabel_classes7 <- c("brain", "muscle", "eye", "heart",
                       "line_noise", "channel_noise", "other")

has_trained_iclabel <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    tryCatch(reticulate::py_module_available("mne") &&
               reticulate::py_module_available("mne_icalabel"),
             error = function(e) FALSE)
}

make_ica_pe <- function() {
  set.seed(1)
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 250)
  eegICA(pe, n_components = 19, method = "fastica")
}

test_that("eegICLabel default backend is the pure-R heuristic (no Python needed)", {
  res <- eegICLabel(make_ica_pe())
  expect_s3_class(res, "data.frame")
  expect_true(all(c("component", .iclabel_classes7, "label") %in% names(res)))
  expect_equal(nrow(res), 19L)
})

test_that("eegICLabel backend='iclabel' runs the trained model when available", {
  skip_if_not(has_trained_iclabel(),
              "reticulate + Python mne/mne-icalabel not available")
  res <- eegICLabel(make_ica_pe(), backend = "iclabel")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 19L)
  probs <- as.matrix(res[, .iclabel_classes7])
  expect_true(all(abs(rowSums(probs) - 1) < 1e-4))     # valid distributions
  expect_true(all(probs >= 0 & probs <= 1))
  expect_true(all(res$label %in% .iclabel_classes7))
})

test_that("eegICLabel backend='iclabel' errors clearly without the Python engine", {
  skip_if(has_trained_iclabel(), "engine present -> error path not exercised")
  expect_error(eegICLabel(make_ica_pe(), backend = "iclabel"),
               "mne-icalabel|reticulate|Python")
})
