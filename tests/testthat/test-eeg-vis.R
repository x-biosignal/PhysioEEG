library(testthat)
library(PhysioEEG)

# ---------------------------------------------------------------------------
# eegPlotSignal
# ---------------------------------------------------------------------------

test_that("eegPlotSignal returns ggplot in stacked mode", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  p <- eegPlotSignal(pe, mode = "stacked")
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotSignal returns ggplot in butterfly mode", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  p <- eegPlotSignal(pe, mode = "butterfly")
  expect_s3_class(p, "gg")
})

test_that("eegPlotSignal returns ggplot in grid mode", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  p <- eegPlotSignal(pe, mode = "grid")
  expect_s3_class(p, "gg")
})

test_that("eegPlotSignal filters channels", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 8, sr = 250)
  p <- eegPlotSignal(pe, channels = c("Fp1", "Fp2"), mode = "stacked")
  expect_s3_class(p, "gg")
  # Build the plot to verify data
  built <- ggplot2::ggplot_build(p)
  # Should only have 2 channel groups
  expect_equal(length(unique(built$data[[1]]$colour)), 2)
})

test_that("eegPlotSignal respects time_range", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  p <- eegPlotSignal(pe, time_range = c(0.5, 1.5), mode = "butterfly")
  expect_s3_class(p, "gg")
  built <- ggplot2::ggplot_build(p)
  x_range <- range(built$data[[1]]$x)
  expect_gte(x_range[1], 0.5)
  expect_lte(x_range[2], 1.5)
})

# ---------------------------------------------------------------------------
# eegPlotERP
# ---------------------------------------------------------------------------

test_that("eegPlotERP returns ggplot with 3D data", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  p <- eegPlotERP(pe)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotERP shows CI ribbon when show_ci=TRUE", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg_erp(n_epochs = 20, n_channels = 4, sr = 250)
  p <- eegPlotERP(pe, show_ci = TRUE)
  built <- ggplot2::ggplot_build(p)
  # With CI ribbon, there should be at least 2 layers (line + ribbon)
  expect_gte(length(built$data), 2)
})

test_that("eegPlotERP separates conditions", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg_erp(n_epochs = 40, n_channels = 4, sr = 250)
  p <- eegPlotERP(pe, conditions = c("target", "standard"))
  expect_s3_class(p, "gg")
  built <- ggplot2::ggplot_build(p)
  # Should have two groups (target and standard)
  n_groups <- length(unique(built$data[[1]]$group))
  expect_gte(n_groups, 2)
})

test_that("eegPlotERP errors on 2D data", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  expect_error(eegPlotERP(pe), "epoched")
})

# ---------------------------------------------------------------------------
# eegPlotTopomap
# ---------------------------------------------------------------------------

test_that("eegPlotTopomap returns ggplot with default positions", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  p <- eegPlotTopomap(pe, time = 0.1)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotTopomap works with custom positions", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  cd <- SummarizedExperiment::colData(pe)
  cd$pos_x <- c(-0.3, 0.3, -0.3, 0.3)
  cd$pos_y <- c(0.5, 0.5, -0.5, -0.5)
  SummarizedExperiment::colData(pe) <- cd
  p <- eegPlotTopomap(pe, time = 0.1)
  expect_s3_class(p, "gg")
})

test_that("eegPlotTopomap works with named values vector", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  cd <- SummarizedExperiment::colData(pe)
  cd$pos_x <- c(-0.3, 0.3, -0.3, 0.3)
  cd$pos_y <- c(0.5, 0.5, -0.5, -0.5)
  SummarizedExperiment::colData(pe) <- cd
  vals <- c(Fp1 = 5, Fp2 = -3, F7 = 2, F3 = -1)
  p <- eegPlotTopomap(pe, values = vals)
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# eegPlotTopomapSeries
# ---------------------------------------------------------------------------

test_that("eegPlotTopomapSeries returns ggplot with multiple times", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 2500, n_channels = 19, sr = 500)
  p <- eegPlotTopomapSeries(pe, times = c(0.1, 0.5, 1.0))
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotTopomapSeries respects ncol parameter", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 2500, n_channels = 19, sr = 500)
  p <- eegPlotTopomapSeries(pe, times = c(0.1, 0.5, 1.0, 2.0), ncol = 2)
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# eegPlotSpectrogram
# ---------------------------------------------------------------------------

test_that("eegPlotSpectrogram returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  p <- eegPlotSpectrogram(pe, channel = 1)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotSpectrogram respects freq_range and time_range", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  p <- eegPlotSpectrogram(pe, channel = 1, freq_range = c(1, 40),
                          time_range = c(0.5, 3.0))
  expect_s3_class(p, "gg")
})

test_that("eegPlotSpectrogram works with log_power=FALSE", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  p <- eegPlotSpectrogram(pe, channel = 1, log_power = FALSE)
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# eegPlotConnectivity
# ---------------------------------------------------------------------------

test_that("eegPlotConnectivity heatmap mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  mat <- matrix(c(1.0, 0.5, 0.3, 0.1,
                  0.5, 1.0, 0.6, 0.2,
                  0.3, 0.6, 1.0, 0.4,
                  0.1, 0.2, 0.4, 1.0), nrow = 4)
  p <- eegPlotConnectivity(pe, method = "heatmap", matrix = mat)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotConnectivity circle mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  mat <- matrix(c(1.0, 0.7, 0.3, 0.1,
                  0.7, 1.0, 0.5, 0.2,
                  0.3, 0.5, 1.0, 0.8,
                  0.1, 0.2, 0.8, 1.0), nrow = 4)
  p <- eegPlotConnectivity(pe, method = "circle", matrix = mat)
  expect_s3_class(p, "gg")
})

test_that("eegPlotConnectivity applies threshold", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  mat <- matrix(c(1.0, 0.8, 0.1, 0.05,
                  0.8, 1.0, 0.2, 0.1,
                  0.1, 0.2, 1.0, 0.9,
                  0.05, 0.1, 0.9, 1.0), nrow = 4)
  p <- eegPlotConnectivity(pe, method = "heatmap", matrix = mat,
                           threshold = 0.5)
  expect_s3_class(p, "gg")
})

test_that("eegPlotConnectivity errors without matrix or metadata", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  expect_error(eegPlotConnectivity(pe, method = "heatmap"),
               "No connectivity matrix")
})

# ---------------------------------------------------------------------------
# eegPlotHypnogram
# ---------------------------------------------------------------------------

test_that("eegPlotHypnogram returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 2, sr = 250)
  stages_df <- data.frame(
    epoch = 1:10,
    stage = c("W", "W", "N1", "N2", "N2", "N3", "N2", "REM", "REM", "W"),
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(pe)$sleep_stages <- stages_df
  p <- eegPlotHypnogram(pe)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotHypnogram with stages provided directly", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 2, sr = 250)
  stages_df <- data.frame(
    epoch = 1:6,
    stage = c("W", "N1", "N2", "N3", "N2", "REM"),
    stringsAsFactors = FALSE
  )
  p <- eegPlotHypnogram(pe, stages = stages_df)
  expect_s3_class(p, "gg")
})

test_that("eegPlotHypnogram respects correct stage ordering", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 2, sr = 250)
  stages_df <- data.frame(
    epoch = 1:5,
    stage = c("W", "REM", "N1", "N2", "N3"),
    stringsAsFactors = FALSE
  )
  p <- eegPlotHypnogram(pe, stages = stages_df)
  built <- ggplot2::ggplot_build(p)
  # W=5 should be highest, N3=1 should be lowest
  y_vals <- built$data[[1]]$y
  # First point should be W (5), last should be N3 (1)
  expect_equal(y_vals[1], 5)  # W
  expect_equal(y_vals[length(y_vals)], 1)  # N3
})

# ---------------------------------------------------------------------------
# eegPlotICA
# ---------------------------------------------------------------------------

test_that("eegPlotICA returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  # Add ICA-like data as an assay (simulated components)
  ica_data <- matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)
  SummarizedExperiment::assays(pe)[["ica"]] <- ica_data
  p <- eegPlotICA(pe, components = 1:3)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotICA defaults to first 10 components", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  # Use raw assay as if it were ICA
  p <- eegPlotICA(pe, assay_name = "raw")
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# eegPlotSource
# ---------------------------------------------------------------------------

test_that("eegPlotSource scatter mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  src <- data.frame(
    x = runif(30, -0.8, 0.8),
    y = runif(30, -0.8, 0.8),
    amplitude = rnorm(30)^2
  )
  p <- eegPlotSource(pe, source_data = src, method = "scatter",
                     threshold_pct = 50)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("eegPlotSource flatmap mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  src <- data.frame(
    x = runif(50, -0.8, 0.8),
    y = runif(50, -0.8, 0.8),
    amplitude = abs(rnorm(50, mean = 5))
  )
  p <- eegPlotSource(pe, source_data = src, method = "flatmap",
                     threshold_pct = 50)
  expect_s3_class(p, "gg")
})

test_that("eegPlotSource works with named numeric vector", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  src_vec <- setNames(abs(rnorm(20)), paste0("src", 1:20))
  p <- eegPlotSource(pe, source_data = src_vec, method = "scatter",
                     threshold_pct = 0)
  expect_s3_class(p, "gg")
})

test_that("eegPlotSource errors without source data", {
  skip_if_not_installed("ggplot2")
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  expect_error(eegPlotSource(pe), "No source data")
})
