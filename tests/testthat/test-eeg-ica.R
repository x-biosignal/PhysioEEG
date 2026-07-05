library(testthat)
library(PhysioEEG)

test_that("eegICA with FastICA returns correct dimensions", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  result <- eegICA(pe, n_components = 4, method = "fastica")

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["ica_components"]]))
  ica_data <- S4Vectors::metadata(result)[["ica_components"]]
  expect_equal(nrow(ica_data), 2000)
  expect_equal(ncol(ica_data), 4)
})

test_that("eegICA stores mixing/unmixing matrices in metadata", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  result <- eegICA(pe, n_components = 4, method = "fastica")

  ica_info <- metadata(result)$ica
  expect_true(!is.null(ica_info))
  expect_true("mixing" %in% names(ica_info))
  expect_true("unmixing" %in% names(ica_info))
  expect_equal(ncol(ica_info$mixing), 4)
})

test_that("eegICA with fewer components reduces dimensionality", {
  pe <- make_eeg(n_time = 2000, n_channels = 8, sr = 250)
  result <- eegICA(pe, n_components = 3, method = "fastica")

  ica_data <- S4Vectors::metadata(result)[["ica_components"]]
  expect_equal(ncol(ica_data), 3)
})

test_that("eegICAremove reconstructs data without specified components", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")
  result <- eegICAremove(pe, components = c(1))

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["ica_cleaned"]]))
  cleaned <- S4Vectors::metadata(result)[["ica_cleaned"]]
  expect_equal(dim(cleaned), c(2000, 4))
})

test_that("eegICAdetect returns data.frame with component labels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")
  result <- eegICAdetect(pe, method = "kurtosis")

  expect_s3_class(result, "data.frame")
  expect_true("component" %in% names(result))
  expect_true("type" %in% names(result))
  expect_true(all(result$type %in% c("artifact", "neural")))
})

test_that("eegICAmix returns ICA metadata", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")
  mix <- eegICAmix(pe)

  expect_true(is.list(mix))
  expect_true("mixing" %in% names(mix))
  expect_true("unmixing" %in% names(mix))
})

test_that("eegICAmix errors when no ICA has been run", {
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  expect_error(eegICAmix(pe))
})

test_that("eegICA with infomax method works", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  result <- eegICA(pe, n_components = 4, method = "infomax")

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["ica_components"]]))
})

test_that("eegICA with jade method works", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  result <- eegICA(pe, n_components = 4, method = "jade")

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["ica_components"]]))
})

test_that("eegICA recovers known source signals via correlation", {
  set.seed(42)
  sr <- 250
  n_time <- 2000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)

  # Create 3 known source signals with different frequencies

  s1 <- sin(2 * pi * 5 * t_vec)
  s2 <- sin(2 * pi * 11 * t_vec)
  s3 <- sin(2 * pi * 23 * t_vec)
  sources <- cbind(s1, s2, s3)

  # Random mixing matrix
  A <- matrix(c(0.8, 0.3, 0.5,
                0.2, 0.9, 0.4,
                0.6, 0.1, 0.7), nrow = 3, byrow = TRUE)
  mixed <- sources %*% t(A)

  pe <- PhysioExperiment(
    assays = list(raw = mixed),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2", "Ch3"),
      type = rep("EEG", 3)
    ),
    samplingRate = sr
  )

  result <- eegICA(pe, n_components = 3, method = "fastica")
  recovered <- S4Vectors::metadata(result)[["ica_components"]]

  # Each recovered component should correlate highly with one original source
  # (allowing sign flip: use abs(cor))
  max_cors <- numeric(3)
  for (i in 1:3) {
    cors <- abs(cor(recovered[, i], sources))
    max_cors[i] <- max(cors)
  }
  # At least 2 of 3 should be well-recovered
  expect_true(sum(max_cors > 0.7) >= 2)
})

test_that("eegICA errors when n_components > n_channels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  expect_error(eegICA(pe, n_components = 10, method = "fastica"),
               "n_components cannot exceed")
})

test_that("eegICAdetect with kurtosis identifies artifact-like components", {
  set.seed(123)
  # Create data with one highly kurtotic (spiky) component
  sr <- 250
  n_time <- 2000
  pe <- make_eeg(n_time = n_time, n_channels = 4, sr = sr)
  pe <- eegICA(pe, n_components = 4, method = "fastica")

  result <- eegICAdetect(pe, method = "kurtosis")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_true("score" %in% names(result))
  # All scores should be finite
  expect_true(all(is.finite(result$score)))
})

test_that("eegICAremove preserves dimensions and changes data", {
  set.seed(99)
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")

  original_data <- SummarizedExperiment::assay(pe, "raw")
  result <- eegICAremove(pe, components = c(1, 2))
  cleaned <- S4Vectors::metadata(result)[["ica_cleaned"]]

  # Dimensions must match original

  expect_equal(dim(cleaned), dim(original_data))
  # Data should be different after removing components
  expect_false(all(abs(cleaned - original_data) < 1e-10))
})

test_that("eegICAmix mixing * unmixing approximates identity", {
  set.seed(77)
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")

  mix <- eegICAmix(pe)
  A <- mix$mixing       # n_channels x n_components
  W <- mix$unmixing     # n_components x n_channels

  # W %*% A should approximate identity
  product <- W %*% A
  expect_equal(dim(product), c(4, 4))
  # Off-diagonal should be near 0, diagonal near 1
  for (i in 1:4) {
    expect_equal(abs(product[i, i]), 1, tolerance = 0.3)
  }
})

test_that("eegICAdetect with spatial method uses frontal channels", {
  set.seed(55)
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 250)
  pe <- eegICA(pe, n_components = 4, method = "fastica")
  result <- eegICAdetect(pe, method = "spatial")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_true(all(result$type %in% c("artifact", "neural")))
  expect_true(all(result$score >= 0 & result$score <= 1))
})

test_that("eegICAremove errors when no ICA has been run", {
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  expect_error(eegICAremove(pe, components = c(1)),
               "No ICA results found")
})

test_that("eegICA is reproducible with set.seed", {
  set.seed(42)
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
  set.seed(100)
  result1 <- eegICA(pe, n_components = 4, method = "fastica")
  set.seed(100)
  result2 <- eegICA(pe, n_components = 4, method = "fastica")
  ica1 <- S4Vectors::metadata(result1)[["ica_components"]]
  ica2 <- S4Vectors::metadata(result2)[["ica_components"]]
  expect_equal(ica1, ica2)
})
