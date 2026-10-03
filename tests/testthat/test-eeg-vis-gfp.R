library(testthat)
library(PhysioEEG)

.fixed_eeg <- function(data, sr = 10, labels = NULL) {
  if (is.null(labels)) {
    labels <- paste0("Ch", seq_len(ncol(data)))
  }
  PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = labels, type = "EEG"),
    samplingRate = sr
  )
}

.gfp_panels <- function(plot) {
  list(butterfly = plot[[1L]], gfp = plot[[2L]])
}

test_that("shared GFP helper matches population and sample definitions", {
  data <- rbind(
    c(-3, 1, 8, 2),
    c(4, 4, 4, 4),
    c(10, -2, 5, 1)
  )
  population <- apply(data, 1L, function(row) {
    sqrt(sum((row - mean(row))^2) / length(row))
  })

  expect_equal(
    PhysioEEG:::.eeg_global_field_power(data, "population"),
    population,
    tolerance = 1e-14
  )
  expect_equal(
    PhysioEEG:::.eeg_global_field_power(data, "sample"),
    apply(data, 1L, stats::sd),
    tolerance = 1e-14
  )
  expect_equal(population[2L], 0)
})

test_that("microstate GFP uses the shared population definition", {
  time <- seq_len(80)
  data <- cbind(
    sin(time / 4),
    cos(time / 7) + 2,
    (time %% 9) / 5 - 1,
    sin(time / 11) * 3
  )
  pe <- .fixed_eeg(data, sr = 20, labels = c("A", "B", "C", "D"))
  result <- eegMicrostates(pe, n_states = 2, method = "pca")
  actual <- S4Vectors::metadata(result)$microstates$gfp
  expected <- sqrt(rowMeans((data - rowMeans(data))^2))

  expect_equal(actual, expected, tolerance = 1e-14)
  expect_equal(
    actual,
    PhysioEEG:::.eeg_global_field_power(data, "population"),
    tolerance = 0
  )
})

test_that("butterfly and GFP panels preserve exact selected channel order", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  data <- cbind(
    A = c(0, 1, 4, 2, -1),
    B = c(5, 2, 1, -2, 3),
    C = c(-2, 6, 0, 5, 1)
  )
  pe <- .fixed_eeg(data, sr = 2, labels = colnames(data))

  plot <- eegPlotButterflyGFP(
    pe,
    channels = c("C", "A"),
    heights = c(2, 1)
  )
  panels <- .gfp_panels(plot)

  expect_true(inherits(plot, "patchwork"))
  expect_s3_class(panels$butterfly, "ggplot")
  expect_s3_class(panels$gfp, "ggplot")
  expect_identical(levels(panels$butterfly$data$channel), c("C", "A"))
  expect_equal(panels$butterfly$data$amplitude, c(data[, "C"], data[, "A"]))
  expect_equal(
    panels$gfp$data$gfp,
    PhysioEEG:::.eeg_global_field_power(data[, c("C", "A")], "population")
  )
  expect_equal(plot$patches$layout$heights, c(2, 1))

  sample_plot <- eegPlotButterflyGFP(
    pe,
    channels = c("C", "A"),
    gfp_definition = "sample"
  )
  expect_equal(
    .gfp_panels(sample_plot)$gfp$data$gfp,
    apply(data[, c("C", "A")], 1L, stats::sd)
  )
})

test_that("time windows display inclusive samples without changing GFP", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  data <- cbind(
    A = c(1, 4, 2, 8, 3),
    B = c(5, -1, 7, 0, 6),
    C = c(2, 9, 1, 3, -2)
  )
  pe <- .fixed_eeg(data, sr = 2, labels = colnames(data))
  plot <- eegPlotButterflyGFP(pe, time_range = c(0.5, 1.5))
  panels <- .gfp_panels(plot)
  expected <- PhysioEEG:::.eeg_global_field_power(data, "population")

  expect_equal(unique(panels$butterfly$data$time), c(0.5, 1, 1.5))
  expect_equal(panels$gfp$data$time, c(0.5, 1, 1.5))
  expect_equal(panels$gfp$data$gfp, expected[2:4])
})

test_that("one-channel GFP has explicit population and sample behavior", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- .fixed_eeg(matrix(1:5, ncol = 1), labels = "A")

  expect_warning(
    plot <- eegPlotButterflyGFP(pe, channels = "A"),
    "not spatially informative"
  )
  expect_equal(.gfp_panels(plot)$gfp$data$gfp, rep(0, 5))
  expect_error(
    eegPlotButterflyGFP(pe, gfp_definition = "sample"),
    "at least two"
  )
})

test_that("butterfly GFP validates dimensions, channels, and numeric inputs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- .fixed_eeg(matrix(seq_len(15), ncol = 3),
                   labels = c("A", "B", "C"))

  expect_error(eegPlotButterflyGFP(pe, channels = "missing"), "Unknown")
  expect_error(eegPlotButterflyGFP(pe, channels = c("A", "A")), "unique")
  expect_error(eegPlotButterflyGFP(pe, time_range = c(1, 1)), "increasing")
  expect_error(eegPlotButterflyGFP(pe, time_range = c(20, 21)), "samples")
  expect_error(eegPlotButterflyGFP(pe, heights = c(1, 0)), "positive")
  expect_error(eegPlotButterflyGFP(pe, show_events = NA), "logical")

  duplicated <- pe
  cd <- SummarizedExperiment::colData(duplicated)
  cd$label <- c("A", "A", "C")
  SummarizedExperiment::colData(duplicated) <- cd
  expect_error(
    eegPlotButterflyGFP(duplicated, channels = "A"),
    "exactly once"
  )

  nonfinite <- pe
  assay <- SummarizedExperiment::assay(nonfinite, "raw")
  assay[2, 2] <- Inf
  SummarizedExperiment::assay(nonfinite, "raw") <- assay
  expect_error(eegPlotButterflyGFP(nonfinite), "finite")

  nonnumeric <- .fixed_eeg(
    matrix(letters[seq_len(15)], ncol = 3),
    labels = c("A", "B", "C")
  )
  expect_error(eegPlotButterflyGFP(nonnumeric), "numeric 2D")

  invalid_rate <- pe
  methods::slot(invalid_rate, "samplingRate") <- 0
  expect_error(eegPlotButterflyGFP(invalid_rate), "finite positive")

  epoched <- make_eeg_erp(n_epochs = 3, n_channels = 3, sr = 10)
  expect_error(eegPlotButterflyGFP(epoched), "2D")
})

test_that("explicit assay selection controls both panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  raw <- matrix(seq_len(18), ncol = 3)
  alternate <- raw * 3 + matrix(c(0, 4, -2), nrow(raw), 3, byrow = TRUE)
  pe <- PhysioCore::PhysioExperiment(
    assays = list(raw = raw, alternate = alternate),
    colData = S4Vectors::DataFrame(
      label = c("A", "B", "C"),
      type = "EEG"
    ),
    samplingRate = 2
  )
  plot <- eegPlotButterflyGFP(pe, assay_name = "alternate")
  panels <- .gfp_panels(plot)

  expect_equal(panels$butterfly$data$amplitude, as.vector(alternate))
  expect_equal(
    panels$gfp$data$gfp,
    PhysioEEG:::.eeg_global_field_power(alternate, "population")
  )
})

test_that("event positions are finite, clipped, and identical in both panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- .fixed_eeg(matrix(seq_len(18), ncol = 3),
                   sr = 2, labels = c("A", "B", "C"))
  for (event_name in c("time", "time_sec", "onset")) {
    events <- data.frame(value = c(-1, 0.5, 1.5, 9))
    names(events) <- event_name
    S4Vectors::metadata(pe)$events <- events
    plot <- eegPlotButterflyGFP(
      pe,
      time_range = c(0.5, 1.5),
      show_events = TRUE
    )
    panels <- .gfp_panels(plot)
    upper_events <- panels$butterfly$layers[[2L]]$data$time
    lower_events <- panels$gfp$layers[[2L]]$data$time

    expect_equal(upper_events, c(0.5, 1.5))
    expect_identical(upper_events, lower_events)
  }

  S4Vectors::metadata(pe)$events$time_sec[2L] <- NA_real_
  expect_error(
    eegPlotButterflyGFP(pe, show_events = TRUE),
    "finite numeric"
  )
})

test_that("panel data and coordinates remain stable through ggplot builds", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  data <- cbind(A = 1:8, B = c(8:1), C = rep(c(-1, 2), 4))
  pe <- .fixed_eeg(data, sr = 4, labels = colnames(data))
  plot <- eegPlotButterflyGFP(pe)
  panels <- .gfp_panels(plot)
  upper_before <- panels$butterfly$data
  lower_before <- panels$gfp$data

  upper_build <- ggplot2::ggplot_build(panels$butterfly)
  lower_build <- ggplot2::ggplot_build(panels$gfp)
  expect_silent(patchwork::patchworkGrob(plot))

  expect_equal(panels$butterfly$data, upper_before)
  expect_equal(panels$gfp$data, lower_before)
  expect_equal(sort(unique(upper_build$data[[1L]]$x)),
               lower_build$data[[1L]]$x)
  expect_equal(unique(upper_before$time), lower_before$time)
})

test_that("ERP theme harmonization preserves means and confidence intervals", {
  skip_if_not_installed("ggplot2")
  set.seed(90212)
  pe <- make_eeg_erp(n_epochs = 12, n_channels = 3, sr = 50)
  plot <- eegPlotERP(pe, channels = c("Fp1", "F7"), show_ci = TRUE)
  data <- SummarizedExperiment::assay(pe, "raw")
  conditions <- S4Vectors::metadata(pe)$conditions
  z_value <- stats::qnorm(0.975)

  for (condition in unique(conditions)) {
    epochs <- which(conditions == condition)
    channel_mean <- apply(data[, c(1, 3), epochs, drop = FALSE],
                          c(1, 3), mean)
    expected_mean <- rowMeans(channel_mean)
    expected_se <- apply(channel_mean, 1, stats::sd) / sqrt(length(epochs))
    actual <- plot$data[plot$data$condition == condition, ]

    expect_equal(actual$amplitude, expected_mean)
    expect_equal(actual$ci_lower, expected_mean - z_value * expected_se)
    expect_equal(actual$ci_upper, expected_mean + z_value * expected_se)
  }

  expect_gte(length(ggplot2::ggplot_build(plot)$data), 2L)
  expect_equal(plot$theme$legend.position, "bottom")
  expect_equal(plot$scales$n(), 2L)
})
