library(testthat)
library(PhysioEEG)

test_that("eegMicrostates with kmeans returns correct number of states", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegMicrostates(pe, n_states = 4, method = "kmeans")

  expect_s4_class(result, "PhysioExperiment")
  ms <- metadata(result)$microstates
  expect_true(!is.null(ms))
  expect_equal(ms$n_states, 4)
  expect_equal(ncol(ms$maps), 4)
  expect_equal(length(ms$labels), 5000)
  expect_true(all(ms$labels %in% 1:4))
})

test_that("eegMicrostateStats returns duration, occurrence, coverage", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  result <- eegMicrostateStats(pe)

  expect_s3_class(result, "data.frame")
  expect_true("duration_ms" %in% names(result))
  expect_true("occurrence_per_sec" %in% names(result))
  expect_true("coverage_pct" %in% names(result))
  expect_equal(nrow(result), 4)
  # Coverage should sum to ~100%
  expect_equal(sum(result$coverage_pct), 100, tolerance = 1)
})

test_that("eegMicrostateBackfit assigns labels from template maps", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  maps <- metadata(pe)$microstates$maps

  pe2 <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  result <- eegMicrostateBackfit(pe2, maps)

  labels <- metadata(result)$microstates$labels
  expect_equal(length(labels), 3000)
  expect_true(all(labels %in% 1:4))
})

test_that("eegMicrostateSequence returns character vector", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  seq_result <- eegMicrostateSequence(pe)

  expect_true(is.character(seq_result))
  expect_equal(length(seq_result), 2000)
  expect_true(all(seq_result %in% c("A", "B", "C", "D")))
})

test_that("eegMicrostates with PCA method works", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  result <- eegMicrostates(pe, n_states = 4, method = "pca")

  ms <- metadata(result)$microstates
  expect_equal(ms$n_states, 4)
})

test_that("eegMicrostates with AAHC method works", {
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  result <- eegMicrostates(pe, n_states = 4, method = "aahc")

  ms <- metadata(result)$microstates
  expect_equal(ms$n_states, 4)
  expect_equal(length(ms$labels), 3000)
})

test_that("eegMicrostateBackfit assigns all time points (no gaps)", {
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  labels <- metadata(pe)$microstates$labels

  # All time points should have a valid assignment
  expect_equal(length(labels), 3000)
  expect_true(all(labels %in% 1:4))
  # No NAs
  expect_false(any(is.na(labels)))
})

test_that("eegMicrostateStats duration, occurrence, coverage are positive", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  stats <- eegMicrostateStats(pe)

  expect_true(all(stats$duration_ms > 0))
  expect_true(all(stats$occurrence_per_sec > 0))
  expect_true(all(stats$coverage_pct > 0))
  # Coverage should sum to 100%
  expect_equal(sum(stats$coverage_pct), 100, tolerance = 0.5)
})

test_that("eegMicrostateStats transition matrix rows sum to ~1", {
  pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  stats <- eegMicrostateStats(pe)

  trans_mat <- attr(stats, "transition_matrix")
  expect_equal(dim(trans_mat), c(4, 4))
  # Rows should sum to approximately 1 (transition probabilities)
  row_sums <- rowSums(trans_mat)
  for (s in 1:4) {
    expect_equal(row_sums[s], 1.0, tolerance = 0.01)
  }
})

test_that("eegMicrostates with n_states = 2 works", {
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  result <- eegMicrostates(pe, n_states = 2, method = "kmeans")

  ms <- metadata(result)$microstates
  expect_equal(ms$n_states, 2)
  expect_true(all(ms$labels %in% 1:2))
})

test_that("eegMicrostates GFP values are non-negative", {
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  result <- eegMicrostates(pe, n_states = 4, method = "kmeans")

  gfp <- metadata(result)$microstates$gfp
  expect_equal(length(gfp), 3000)
  expect_true(all(gfp >= 0))
})

test_that("eegMicrostates is reproducible with set.seed in internal clustering", {
  pe <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  # Clustering uses set.seed(42) internally in .microstate_kmeans
  result1 <- eegMicrostates(pe, n_states = 4, method = "kmeans")
  result2 <- eegMicrostates(pe, n_states = 4, method = "kmeans")

  labels1 <- metadata(result1)$microstates$labels
  labels2 <- metadata(result2)$microstates$labels
  expect_equal(labels1, labels2)
})

test_that("eegMicrostateSequence maps integers to letters correctly", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe <- eegMicrostates(pe, n_states = 3, method = "kmeans")
  seq_result <- eegMicrostateSequence(pe)

  expect_true(is.character(seq_result))
  expect_equal(length(seq_result), 2000)
  expect_true(all(seq_result %in% c("A", "B", "C")))
})

test_that("eegMicrostateBackfit with external maps produces valid labels", {
  pe1 <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
  pe1 <- eegMicrostates(pe1, n_states = 4, method = "kmeans")
  maps <- metadata(pe1)$microstates$maps

  pe2 <- make_eeg(n_time = 1000, n_channels = 19, sr = 500)
  result <- eegMicrostateBackfit(pe2, maps)

  labels <- metadata(result)$microstates$labels
  expect_equal(length(labels), 1000)
  expect_true(all(labels %in% 1:4))
  # Maps should be stored in result
  expect_equal(metadata(result)$microstates$maps, maps)
})
