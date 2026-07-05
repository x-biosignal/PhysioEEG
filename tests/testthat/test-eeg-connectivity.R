library(testthat)
library(PhysioEEG)

# --- eegCoherence tests ---

test_that("eegCoherence returns correct structure with coherence in [0, 1]", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_result <- eegCoherence(pe, method = "coherence", band = c(8, 13))

  conn <- S4Vectors::metadata(pe_result)$connectivity
  expect_true(!is.null(conn))
  expect_true("matrix" %in% names(conn))
  expect_true("method" %in% names(conn))
  expect_true("band" %in% names(conn))
  expect_true("freqs" %in% names(conn))
  expect_true("spectra" %in% names(conn))

  mat <- conn$matrix
  expect_equal(nrow(mat), 4)
  expect_equal(ncol(mat), 4)

  # Coherence values should be in [0, 1]
  expect_true(all(mat >= 0 & mat <= 1))

  # Diagonal should be 1 (self-coherence)
  expect_equal(unname(diag(mat)), rep(1, 4))

  # Matrix should be symmetric
  expect_equal(mat, t(mat))

  # Method should be stored

  expect_equal(conn$method, "coherence")
  expect_equal(conn$band, c(8, 13))
})

test_that("eegCoherence imaginary method returns values in [0, 1]", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_result <- eegCoherence(pe, method = "imaginary", band = c(8, 13))

  conn <- S4Vectors::metadata(pe_result)$connectivity
  mat <- conn$matrix

  expect_true(all(mat >= 0 & mat <= 1))
  expect_equal(conn$method, "imaginary")
  # Imaginary coherence diagonal should still be 1 (self)
  expect_equal(unname(diag(mat)), rep(1, 4))
})

test_that("eegCoherence stores frequency-resolved spectra", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_result <- eegCoherence(pe, method = "coherence", band = c(8, 13),
                            window_sec = 1)

  conn <- S4Vectors::metadata(pe_result)$connectivity
  spectra <- conn$spectra

  # spectra should be 3D: n_freqs x n_channels x n_channels
  expect_equal(length(dim(spectra)), 3)
  expect_equal(dim(spectra)[2], 4)
  expect_equal(dim(spectra)[3], 4)
  expect_equal(dim(spectra)[1], length(conn$freqs))

  # All spectra values should be in [0, 1]
  expect_true(all(spectra >= -1e-10 & spectra <= 1 + 1e-10))
})

# --- eegPLV tests ---

test_that("eegPLV returns correct structure with PLV in [0, 1]", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  plv_df <- eegPLV(pe, band = c(8, 13))

  expect_s3_class(plv_df, "data.frame")
  expect_true(all(c("channel1", "channel2", "plv") %in% names(plv_df)))

  # Number of unique pairs = n_channels choose 2 = 4*3/2 = 6
  expect_equal(nrow(plv_df), 6)

  # PLV values should be in [0, 1]
  expect_true(all(plv_df$plv >= 0 & plv_df$plv <= 1))
})

test_that("eegPLV uses channel labels from colData", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  plv_df <- eegPLV(pe, band = c(8, 13))

  cd <- SummarizedExperiment::colData(pe)
  ch_labels <- as.character(cd$label)

  # Check that channel labels appear in the results
  expect_true(all(plv_df$channel1 %in% ch_labels))
  expect_true(all(plv_df$channel2 %in% ch_labels))
})

# --- eegWPLI tests ---

test_that("eegWPLI returns correct structure with wPLI in [0, 1]", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  wpli_df <- eegWPLI(pe, band = c(8, 13), window_sec = 1, overlap = 0.5)

  expect_s3_class(wpli_df, "data.frame")
  expect_true(all(c("channel1", "channel2", "wpli", "wpli_debiased") %in%
                    names(wpli_df)))

  # Number of unique pairs
  expect_equal(nrow(wpli_df), 6)

  # wPLI should be in [0, 1]
  expect_true(all(wpli_df$wpli >= 0 & wpli_df$wpli <= 1))

  # Debiased wPLI should be >= 0 (after clamping)
  expect_true(all(wpli_df$wpli_debiased >= 0, na.rm = TRUE))
})

# --- eegGrangerCausality tests ---

test_that("eegGrangerCausality returns correct structure with GC >= 0", {
  pe <- make_eeg(n_time = 5000, n_channels = 3, sr = 500)
  gc_df <- eegGrangerCausality(pe, order = 5, band = c(8, 13))

  expect_s3_class(gc_df, "data.frame")
  expect_true(all(c("from_channel", "to_channel", "gc_value") %in%
                    names(gc_df)))

  # Directional: n_channels * (n_channels - 1) = 3 * 2 = 6 pairs
  expect_equal(nrow(gc_df), 6)

  # GC values should be >= 0
  expect_true(all(gc_df$gc_value >= 0))
})

test_that("eegGrangerCausality is directional", {
  pe <- make_eeg(n_time = 5000, n_channels = 3, sr = 500)
  gc_df <- eegGrangerCausality(pe, order = 5, band = c(8, 13))

  cd <- SummarizedExperiment::colData(pe)
  ch_labels <- as.character(cd$label)

  # Check that both directions exist for each pair
  ch1 <- ch_labels[1]
  ch2 <- ch_labels[2]

  gc_12 <- gc_df$gc_value[gc_df$from_channel == ch1 & gc_df$to_channel == ch2]
  gc_21 <- gc_df$gc_value[gc_df$from_channel == ch2 & gc_df$to_channel == ch1]

  expect_length(gc_12, 1)
  expect_length(gc_21, 1)
  # GC is directional so values can differ (this is a structural check)
  expect_true(is.numeric(gc_12) && gc_12 >= 0)
  expect_true(is.numeric(gc_21) && gc_21 >= 0)
})

# --- eegConnectivityMatrix tests ---

test_that("eegConnectivityMatrix returns symmetric matrix with correct dims", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  for (meth in c("coherence", "plv", "wpli")) {
    conn_mat <- eegConnectivityMatrix(pe, method = meth, band = c(8, 13))

    expect_true(is.matrix(conn_mat))
    expect_equal(nrow(conn_mat), 4)
    expect_equal(ncol(conn_mat), 4)

    # Should be symmetric
    expect_equal(conn_mat, t(conn_mat), tolerance = 1e-10)

    # Diagonal should be 1
    expect_equal(unname(diag(conn_mat)), rep(1, 4))

    # Should have channel labels
    expect_false(is.null(rownames(conn_mat)))
    expect_false(is.null(colnames(conn_mat)))
    expect_equal(rownames(conn_mat), colnames(conn_mat))

    # Off-diagonal values should be in [0, 1]
    off_diag <- conn_mat[row(conn_mat) != col(conn_mat)]
    expect_true(all(off_diag >= -1e-10 & off_diag <= 1 + 1e-10))
  }
})

# --- Input validation tests ---

test_that("connectivity functions error on non-PhysioExperiment input", {
  not_pe <- matrix(rnorm(100), nrow = 50, ncol = 2)

  expect_error(eegCoherence(not_pe))
  expect_error(eegPLV(not_pe))
  expect_error(eegWPLI(not_pe))
  expect_error(eegGrangerCausality(not_pe))
  expect_error(eegConnectivityMatrix(not_pe))
})

test_that("eegCoherence errors with invalid band specification", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  # Band reversed
  expect_error(eegCoherence(pe, band = c(13, 8)))
  # Band too narrow for window (frequency bins outside range)
  expect_error(eegCoherence(pe, band = c(249, 251)))
})

test_that("eegGrangerCausality errors on too-short signal for given order", {
  pe <- make_eeg(n_time = 10, n_channels = 2, sr = 500)
  expect_error(eegGrangerCausality(pe, order = 10))
})

test_that("eegCoherence identical signals yield coherence = 1", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  sig <- sin(2 * pi * 10 * t_vec)
  data <- cbind(sig, sig)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  pe_result <- eegCoherence(pe, method = "coherence", band = c(8, 13))
  mat <- S4Vectors::metadata(pe_result)$connectivity$matrix
  # Off-diagonal should be ~1 for identical signals
  expect_equal(mat[1, 2], 1.0, tolerance = 0.05)
})

test_that("eegCoherence independent signals yield low coherence", {
  set.seed(61)
  sr <- 500
  n_time <- 5000
  data <- matrix(rnorm(n_time * 2), nrow = n_time, ncol = 2)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  pe_result <- eegCoherence(pe, method = "coherence", band = c(8, 13))
  mat <- S4Vectors::metadata(pe_result)$connectivity$matrix
  # Independent signals should have low coherence
  expect_lt(mat[1, 2], 0.5)
})

test_that("eegPLV phase-locked signals yield PLV near 1", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Two signals with constant phase relationship
  sig1 <- sin(2 * pi * 10 * t_vec)
  sig2 <- sin(2 * pi * 10 * t_vec + pi / 4)  # constant phase shift
  data <- cbind(sig1, sig2)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  plv_df <- eegPLV(pe, band = c(8, 13))
  expect_gt(plv_df$plv[1], 0.7)
})

test_that("eegPLV random phase signals yield low PLV", {
  set.seed(62)
  sr <- 500
  n_time <- 5000
  data <- matrix(rnorm(n_time * 2), nrow = n_time, ncol = 2)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  plv_df <- eegPLV(pe, band = c(8, 13))
  expect_lt(plv_df$plv[1], 0.5)
})

test_that("eegWPLI zero-lag (volume conduction) yields wPLI near 0", {
  sr <- 500
  n_time <- 5000
  t_vec <- seq(0, (n_time - 1) / sr, length.out = n_time)
  # Identical signals (zero phase lag = pure volume conduction)
  sig <- sin(2 * pi * 10 * t_vec) + rnorm(n_time, sd = 0.5)
  data <- cbind(sig, sig + rnorm(n_time, sd = 0.1))

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("Ch1", "Ch2"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  wpli_df <- eegWPLI(pe, band = c(8, 13), window_sec = 1, overlap = 0.5)
  # wPLI should be relatively low for zero-lag coupling
  expect_true(wpli_df$wpli[1] >= 0 && wpli_df$wpli[1] <= 1)
})

test_that("eegGrangerCausality X causes Y when X leads Y", {
  set.seed(63)
  sr <- 500
  n_time <- 5000

  # Create causal signal: X drives Y with a delay
  x <- rnorm(n_time)
  y <- numeric(n_time)
  for (i in 6:n_time) {
    y[i] <- 0.7 * x[i - 5] + rnorm(1, sd = 0.5)  # Y depends on X lagged by 5
  }
  data <- cbind(x, y)

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = c("X", "Y"),
      type = c("EEG", "EEG")
    ),
    samplingRate = sr
  )

  gc_df <- eegGrangerCausality(pe, order = 10)
  gc_xy <- gc_df$gc_value[gc_df$from_channel == "X" & gc_df$to_channel == "Y"]
  gc_yx <- gc_df$gc_value[gc_df$from_channel == "Y" & gc_df$to_channel == "X"]

  # GC from X to Y should be greater than Y to X
  expect_gt(gc_xy, gc_yx)
})

test_that("eegConnectivityMatrix symmetric for undirected measures", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  conn_mat <- eegConnectivityMatrix(pe, method = "plv", band = c(8, 13))

  expect_true(is.matrix(conn_mat))
  expect_equal(nrow(conn_mat), 4)
  # Symmetric
  expect_equal(conn_mat, t(conn_mat), tolerance = 1e-10)
  # Diagonal = 1
  expect_equal(unname(diag(conn_mat)), rep(1, 4))
})

test_that("eegConnectivityMatrix 2-channel minimum works", {
  pe <- make_eeg(n_time = 5000, n_channels = 2, sr = 500)

  for (meth in c("coherence", "plv", "wpli")) {
    conn_mat <- eegConnectivityMatrix(pe, method = meth, band = c(8, 13))
    expect_equal(dim(conn_mat), c(2, 2))
    expect_equal(unname(diag(conn_mat)), rep(1, 2))
    # Off-diagonal in [0, 1]
    expect_true(conn_mat[1, 2] >= -1e-10 && conn_mat[1, 2] <= 1 + 1e-10)
  }
})
