library(testthat)
library(PhysioEEG)

.ws913_axis_names <- function(x) {
  vapply(x, function(value) sprintf("%.17g", value), character(1))
}

.ws913_legacy_data <- function(x, channel = 1L, log_power = TRUE) {
  data <- SummarizedExperiment::assay(x, defaultAssay(x))
  sampling_rate <- samplingRate(x)
  signal <- data[, channel]
  n_samples <- length(signal)
  window_length <- min(as.integer(round(0.5 * sampling_rate)), n_samples)
  if (window_length %% 2L != 0L) {
    window_length <- window_length - 1L
  }
  hop <- max(1L, window_length %/% 4L)
  window <- 0.5 * (
    1 - cos(
      2 * pi * seq(0, window_length - 1L) / (window_length - 1L)
    )
  )
  n_frequency <- window_length %/% 2L + 1L
  frequency <- seq(0, sampling_rate / 2, length.out = n_frequency)
  starts <- seq.int(
    1L,
    n_samples - window_length + 1L,
    by = hop
  )
  time <- (starts + window_length / 2 - 1) / sampling_rate
  power <- matrix(0, nrow = length(starts), ncol = n_frequency)
  for (window_index in seq_along(starts)) {
    indices <- starts[window_index]:(
      starts[window_index] + window_length - 1L
    )
    transformed <- fft(signal[indices] * window)[seq_len(n_frequency)]
    power[window_index, ] <- Mod(transformed)^2 / window_length
  }
  if (log_power) {
    power <- 10 * log10(pmax(power, .Machine$double.eps))
  }
  result <- expand.grid(time = time, freq = frequency)
  result$power <- as.vector(power)
  result
}

.ws913_mask <- function(extracted, value = FALSE) {
  matrix(
    value,
    nrow = nrow(extracted$power),
    ncol = ncol(extracted$power),
    dimnames = dimnames(extracted$power)
  )
}

.ws913_reference_cluster_max <- function(statistic, critical) {
  label_mass <- function(active, values) {
    ids <- matrix(0L, nrow(active), ncol(active))
    masses <- numeric()
    next_id <- 0L
    for (frequency_index in seq_len(ncol(active))) {
      for (time_index in seq_len(nrow(active))) {
        if (!active[time_index, frequency_index] ||
            ids[time_index, frequency_index] != 0L) {
          next
        }
        next_id <- next_id + 1L
        queue <- matrix(c(time_index, frequency_index), ncol = 2L)
        ids[time_index, frequency_index] <- next_id
        mass <- 0
        while (nrow(queue) > 0L) {
          current <- queue[1L, ]
          queue <- queue[-1L, , drop = FALSE]
          mass <- mass + values[current[1L], current[2L]]
          neighbours <- rbind(
            c(current[1L] - 1L, current[2L]),
            c(current[1L] + 1L, current[2L]),
            c(current[1L], current[2L] - 1L),
            c(current[1L], current[2L] + 1L)
          )
          for (neighbour_index in seq_len(nrow(neighbours))) {
            neighbour <- neighbours[neighbour_index, ]
            if (neighbour[1L] < 1L || neighbour[1L] > nrow(active) ||
                neighbour[2L] < 1L || neighbour[2L] > ncol(active) ||
                !active[neighbour[1L], neighbour[2L]] ||
                ids[neighbour[1L], neighbour[2L]] != 0L) {
              next
            }
            ids[neighbour[1L], neighbour[2L]] <- next_id
            queue <- rbind(queue, neighbour)
          }
        }
        masses <- c(masses, mass)
      }
    }
    masses
  }
  positive <- label_mass(
    is.finite(statistic) & statistic >= critical,
    statistic
  )
  negative <- label_mass(
    is.finite(statistic) & statistic <= -critical,
    -statistic
  )
  masses <- c(positive, negative)
  if (length(masses) == 0L) 0 else max(masses)
}

test_that("mask NULL preserves the legacy raster contract", {
  skip_if_not_installed("ggplot2")
  set.seed(91301)
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 200)
  expected <- .ws913_legacy_data(pe)
  actual <- eegPlotSpectrogram(pe, channel = 1)

  expect_s3_class(actual, "ggplot")
  expect_length(actual$layers, 1L)
  expect_equal(actual$data$time, expected$time, tolerance = 0)
  expect_equal(actual$data$freq, expected$freq, tolerance = 0)
  expect_equal(actual$data$power, expected$power, tolerance = 0)

  reference <- ggplot2::ggplot(
    expected,
    ggplot2::aes(
      x = .data$time,
      y = .data$freq,
      fill = .data$power
    )
  ) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(option = "D", name = "Power (dB)")
  actual_layer <- ggplot2::ggplot_build(actual)$data[[1L]]
  reference_layer <- ggplot2::ggplot_build(reference)$data[[1L]]
  expect_equal(
    actual_layer[c("x", "y", "fill")],
    reference_layer[c("x", "y", "fill")],
    tolerance = 0
  )
})

test_that("constant masks change only display opacity", {
  skip_if_not_installed("ggplot2")
  set.seed(91302)
  pe <- make_eeg(n_time = 800, n_channels = 2, sr = 200)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  unmasked <- eegPlotSpectrogram(pe, channel = 1)
  all_true <- eegPlotSpectrogram(
    pe,
    channel = 1,
    mask = .ws913_mask(extracted, TRUE)
  )
  all_false <- eegPlotSpectrogram(
    pe,
    channel = 1,
    mask = .ws913_mask(extracted, FALSE),
    mask_alpha = 0.23
  )

  expect_equal(all_true$data$power, unmasked$data$power, tolerance = 0)
  expect_equal(all_false$data$power, unmasked$data$power, tolerance = 0)
  expect_true(all(all_true$data$display_alpha == 1))
  expect_true(all(all_false$data$display_alpha == 0.23))
  expect_equal(
    ggplot2::ggplot_build(all_true)$data[[1L]]$alpha,
    rep(1, nrow(all_true$data))
  )
  expect_equal(
    ggplot2::ggplot_build(all_false)$data[[1L]]$alpha,
    rep(0.23, nrow(all_false$data))
  )
  expect_length(all_true$layers, 1L)
  expect_length(all_false$layers, 1L)
})

test_that("asymmetric masks retain time-frequency orientation", {
  skip_if_not_installed("ggplot2")
  set.seed(91303)
  pe <- make_eeg(n_time = 700, n_channels = 2, sr = 200)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  mask <- .ws913_mask(extracted)
  mask[2L, 3L] <- TRUE
  mask[4L, 2L] <- TRUE
  plot <- eegPlotSpectrogram(pe, channel = 1, mask = mask, contour = FALSE)
  selected <- plot$data[plot$data$significant, c("time", "freq")]

  expect_equal(
    selected$time,
    extracted$time[c(4L, 2L)],
    tolerance = 0
  )
  expect_equal(
    selected$freq,
    extracted$frequency[c(2L, 3L)],
    tolerance = 0
  )
})

test_that("range filters subset power and mask with identical indices", {
  skip_if_not_installed("ggplot2")
  set.seed(91304)
  pe <- make_eeg(n_time = 1000, n_channels = 2, sr = 200)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  mask <- .ws913_mask(extracted)
  mask[seq(2L, nrow(mask), by = 3L), seq(2L, ncol(mask), by = 4L)] <- TRUE
  time_range <- extracted$time[c(3L, 8L)]
  frequency_range <- extracted$frequency[c(2L, 12L)]
  plot <- eegPlotSpectrogram(
    pe,
    channel = 1,
    time_range = time_range,
    freq_range = frequency_range,
    mask = mask,
    contour = FALSE
  )
  time_indices <- which(
    extracted$time >= time_range[1L] &
      extracted$time <= time_range[2L]
  )
  frequency_indices <- which(
    extracted$frequency >= frequency_range[1L] &
      extracted$frequency <= frequency_range[2L]
  )
  expected <- expand.grid(
    time = extracted$time[time_indices],
    freq = extracted$frequency[frequency_indices]
  )
  expected$power <- as.vector(
    extracted$power[time_indices, frequency_indices, drop = FALSE]
  )
  expected$significant <- as.vector(
    mask[time_indices, frequency_indices, drop = FALSE]
  )

  expect_equal(plot$data$time, expected$time, tolerance = 0)
  expect_equal(plot$data$freq, expected$freq, tolerance = 0)
  expect_equal(plot$data$power, expected$power, tolerance = 0)
  expect_identical(plot$data$significant, expected$significant)
})

test_that("contours use logical cell boundaries and skip constant masks", {
  skip_if_not_installed("ggplot2")
  set.seed(91305)
  pe <- make_eeg(n_time = 800, n_channels = 2, sr = 200)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  mask <- .ws913_mask(extracted)
  mask[2:3, 2:3] <- TRUE
  mask[6:7, 6:7] <- TRUE
  outlined <- eegPlotSpectrogram(pe, channel = 1, mask = mask)
  constant <- eegPlotSpectrogram(
    pe,
    channel = 1,
    mask = .ws913_mask(extracted, TRUE)
  )

  expect_length(outlined$layers, 2L)
  expect_gt(nrow(outlined$layers[[2L]]$data), 0L)
  expect_true(all(
    outlined$layers[[2L]]$data$x == outlined$layers[[2L]]$data$xend |
      outlined$layers[[2L]]$data$y == outlined$layers[[2L]]$data$yend
  ))
  expect_length(constant$layers, 1L)
})

test_that("logical, numeric p-value, and inference-result masks are accepted", {
  skip_if_not_installed("ggplot2")
  set.seed(91306)
  pe <- make_eeg(n_time = 500, n_channels = 1, sr = 100)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  logical_mask <- .ws913_mask(extracted)
  logical_mask[2L, 2L] <- TRUE
  p_values <- matrix(
    0.5,
    nrow(logical_mask),
    ncol(logical_mask),
    dimnames = dimnames(logical_mask)
  )
  p_values[2L, 2L] <- 0.01
  attr(p_values, "alpha") <- 0.05
  tf <- array(
    rnorm(length(logical_mask) * 4L),
    dim = c(dim(logical_mask), 4L),
    dimnames = c(dimnames(logical_mask), list(NULL))
  )
  inference <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = seq_len(nrow(logical_mask)) <= 2L,
    method = "threshold"
  )

  expect_s3_class(
    eegPlotSpectrogram(pe, mask = logical_mask),
    "ggplot"
  )
  numeric_plot <- eegPlotSpectrogram(pe, mask = p_values, contour = FALSE)
  expect_identical(numeric_plot$data$significant, as.vector(p_values <= 0.05))
  expect_s3_class(
    eegPlotSpectrogram(pe, mask = inference),
    "ggplot"
  )
})

test_that("invalid masks and display controls fail before drawing", {
  skip_if_not_installed("ggplot2")
  set.seed(91307)
  pe <- make_eeg(n_time = 500, n_channels = 1, sr = 100)
  extracted <- PhysioEEG:::.eeg_extract_time_frequency(pe, 1)
  valid <- .ws913_mask(extracted)
  transposed <- t(valid)
  misnamed <- valid
  dimnames(misnamed)[[1L]][1L] <- "wrong"
  missing <- valid
  missing[1L] <- NA
  p_values <- matrix(
    0.5,
    nrow(valid),
    ncol(valid),
    dimnames = dimnames(valid)
  )
  nonfinite <- p_values
  nonfinite[1L] <- Inf
  attr(nonfinite, "alpha") <- 0.05
  out_of_range <- p_values
  out_of_range[1L] <- 1.1
  attr(out_of_range, "alpha") <- 0.05

  expect_error(eegPlotSpectrogram(pe, mask = transposed), "matrix")
  expect_error(eegPlotSpectrogram(pe, mask = misnamed), "dimnames")
  expect_error(eegPlotSpectrogram(pe, mask = missing), "missing")
  expect_error(eegPlotSpectrogram(pe, mask = p_values), "alpha")
  expect_error(eegPlotSpectrogram(pe, mask = nonfinite), "finite")
  expect_error(eegPlotSpectrogram(pe, mask = out_of_range), "\\[0, 1\\]")
  expect_error(eegPlotSpectrogram(pe, mask = as.vector(valid)), "matrix")
  expect_error(eegPlotSpectrogram(pe, mask_alpha = -0.1), "mask_alpha")
  expect_error(eegPlotSpectrogram(pe, contour = NA), "contour")
})

test_that("raw and stored products preserve channel and axis orientation", {
  skip_if_not_installed("ggplot2")
  set.seed(91308)
  pe <- make_eeg(n_time = 100, n_channels = 2, sr = 100)
  metadata <- S4Vectors::metadata(pe)
  stft <- array(seq_len(3L * 4L * 2L), dim = c(3L, 4L, 2L))
  metadata$stft_power <- stft
  metadata$stft <- list(
    time_axis = c(0.25, 0.5, 0.75),
    freq_axis = c(0, 5, 10, 20)
  )
  S4Vectors::metadata(pe) <- metadata
  stft_plot <- eegPlotSpectrogram(
    pe,
    channel = "Fp2",
    assay_name = "stft_power",
    log_power = FALSE
  )
  expect_equal(stft_plot$data$power, as.vector(stft[, , 2L]), tolerance = 0)
  expect_equal(unique(stft_plot$data$time), metadata$stft$time_axis)
  expect_equal(unique(stft_plot$data$freq), metadata$stft$freq_axis)
  expect_equal(
    eegPlotSpectrogram(pe, channel = 2, log_power = FALSE)$data$power,
    as.vector(stft[, , 2L]),
    tolerance = 0
  )

  wavelet_pe <- make_eeg(n_time = 6, n_channels = 2, sr = 20)
  metadata <- S4Vectors::metadata(wavelet_pe)
  wavelet <- array(seq_len(6L * 3L * 2L), dim = c(6L, 3L, 2L))
  metadata$wavelet_power <- wavelet
  metadata$wavelet <- list(frequencies = c(3, 7, 9))
  S4Vectors::metadata(wavelet_pe) <- metadata
  wavelet_plot <- eegPlotSpectrogram(
    wavelet_pe,
    channel = 1,
    assay_name = "wavelet_power",
    log_power = FALSE
  )
  expect_equal(wavelet_plot$data$power, as.vector(wavelet[, , 1L]))
  expect_equal(unique(wavelet_plot$data$time), (0:5) / 20)
  expect_equal(unique(wavelet_plot$data$freq), c(3, 7, 9))

  ersp_pe <- make_eeg_erp(
    n_epochs = 3,
    n_channels = 2,
    sr = 20,
    epoch_sec = 0.3
  )
  metadata <- S4Vectors::metadata(ersp_pe)
  ersp <- array(seq(-4, 4, length.out = 6L * 2L * 2L), c(6L, 2L, 2L))
  metadata$ersp_data <- ersp
  metadata$ersp <- list(frequencies = c(4, 8))
  S4Vectors::metadata(ersp_pe) <- metadata
  ersp_plot <- eegPlotSpectrogram(
    ersp_pe,
    channel = 2,
    assay_name = "ersp_data"
  )
  expect_equal(ersp_plot$data$power, as.vector(ersp[, , 2L]), tolerance = 0)
  expect_equal(unique(ersp_plot$data$time), (0:5) / 20)
  expect_equal(unique(ersp_plot$data$freq), c(4, 8))
  expect_equal(
    eegPlotSpectrogram(ersp_pe, channel = 2)$data$power,
    as.vector(ersp[, , 2L]),
    tolerance = 0
  )
})

test_that("ambiguous sources, channels, axes, and raw 3D inputs fail", {
  skip_if_not_installed("ggplot2")
  set.seed(91309)
  pe <- make_eeg(n_time = 20, n_channels = 2, sr = 20)
  expect_error(eegPlotSpectrogram(pe, channel = 0), "integer")
  expect_error(eegPlotSpectrogram(pe, channel = 1.5), "integer")
  expect_error(eegPlotSpectrogram(pe, channel = "Fp"), "not found")
  duplicated <- pe
  col_data <- SummarizedExperiment::colData(duplicated)
  col_data$label <- c("same", "same")
  SummarizedExperiment::colData(duplicated) <- col_data
  expect_error(eegPlotSpectrogram(duplicated, channel = "same"), "duplicated")

  epoched <- make_eeg_erp(
    n_epochs = 3,
    n_channels = 2,
    sr = 20,
    epoch_sec = 0.3
  )
  expect_error(eegPlotSpectrogram(epoched), "raw 3D")

  metadata <- S4Vectors::metadata(pe)
  metadata$stft_power <- array(1, c(3, 2, 2))
  metadata$stft <- list(
    time_axis = c(0.1, 0.1, 0.2),
    freq_axis = c(5, 10)
  )
  S4Vectors::metadata(pe) <- metadata
  expect_error(
    eegPlotSpectrogram(pe, assay_name = "stft_power"),
    "strictly increasing"
  )

  short <- make_eeg(n_time = 1, n_channels = 1, sr = 100)
  expect_error(eegPlotSpectrogram(short), "too short")
  invalid_rate <- make_eeg(n_time = 20, n_channels = 1, sr = 20)
  methods::slot(invalid_rate, "samplingRate") <- NaN
  expect_error(eegPlotSpectrogram(invalid_rate), "sampling rate")
  nonfinite <- make_eeg(n_time = 20, n_channels = 1, sr = 20)
  assay <- SummarizedExperiment::assay(nonfinite, "raw")
  assay[1L] <- Inf
  SummarizedExperiment::assay(nonfinite, "raw") <- assay
  expect_error(eegPlotSpectrogram(nonfinite), "finite numeric")
})

test_that("threshold helper matches hand-computed t and BH values", {
  set.seed(91310)
  tf <- array(
    rnorm(5L * 3L * 6L),
    dim = c(5L, 3L, 6L),
    dimnames = list(
      as.character(c(-0.2, -0.1, 0, 0.1, 0.2)),
      as.character(c(5, 10, 20)),
      NULL
    )
  )
  result <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = c(-0.2, -0.1),
    method = "threshold",
    p_adjust_method = "BH"
  )
  baseline_mean <- apply(tf[1:2, , , drop = FALSE], c(2, 3), mean)
  contrasts <- sweep(tf, c(2, 3), baseline_mean, FUN = "-")
  reference_t <- apply(
    contrasts,
    c(1, 2),
    function(values) mean(values) / (stats::sd(values) / sqrt(6))
  )
  reference_p <- 2 * stats::pt(-abs(reference_t), df = 5)
  reference_bh <- matrix(
    stats::p.adjust(as.vector(reference_p), method = "BH"),
    nrow = 5L,
    ncol = 3L
  )

  expect_equal(
    as.numeric(result$statistic),
    as.numeric(reference_t),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(result$pointwise_p),
    as.numeric(reference_p),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(result$p_adjusted),
    as.numeric(reference_bh),
    tolerance = 1e-12
  )
  expect_identical(result$mask, result$p_adjusted <= result$alpha)
  expect_identical(result$baseline_indices, 1:2)
})

test_that("component labelling is 4-neighbour and sign-separated", {
  statistic <- matrix(0, nrow = 4L, ncol = 4L)
  statistic[1L, 1L] <- 3
  statistic[2L, 2L] <- 4
  statistic[3L, 3L] <- -5
  statistic[3L, 4L] <- -2.5
  statistic[4L, 3L] <- 2.5
  clusters <- PhysioEEG:::.tf_cluster_field(
    statistic,
    critical_value = 2,
    alternative = "two.sided",
    time = c(0, 0.1, 0.2, 0.3),
    frequency = c(5, 10, 15, 20)
  )

  expect_equal(nrow(clusters$cluster_table), 4L)
  expect_equal(sort(clusters$cluster_table$size), c(1L, 1L, 1L, 2L))
  expect_equal(sort(clusters$cluster_table$mass), sort(c(3, 4, 7.5, 2.5)))
  negative <- clusters$cluster_table[
    clusters$cluster_table$sign == "negative",
  ]
  expect_equal(negative$size, 2L)
  expect_equal(negative$time_index_min, 3L)
  expect_equal(negative$frequency_index_max, 4L)
  expect_false(
    clusters$cluster_id[1L, 1L] == clusters$cluster_id[2L, 2L]
  )
  expect_false(
    clusters$cluster_id[3L, 3L] == clusters$cluster_id[4L, 3L]
  )
})

test_that("exact sign enumeration matches an independent max-cluster null", {
  set.seed(91311)
  tf <- array(
    rnorm(5L * 3L * 4L, sd = 0.5),
    dim = c(5L, 3L, 4L),
    dimnames = list(
      as.character(seq(-0.2, 0.2, by = 0.1)),
      as.character(c(5, 10, 15)),
      NULL
    )
  )
  tf[3:5, 2:3, ] <- tf[3:5, 2:3, ] + 2
  result <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = c(-0.2, -0.1),
    method = "cluster",
    alpha = 0.1,
    n_permutations = 15L
  )
  baseline_mean <- apply(tf[1:2, , , drop = FALSE], c(2, 3), mean)
  contrasts <- sweep(tf, c(2, 3), baseline_mean, FUN = "-")
  critical <- stats::qt(1 - 0.05 / 2, df = 3)
  reference_null <- numeric(15L)
  for (pattern in 0:14) {
    signs <- vapply(
      0:3,
      function(bit) {
        if (bitwAnd(pattern, bitwShiftL(1L, bit)) == 0L) -1 else 1
      },
      numeric(1)
    )
    permuted <- sweep(contrasts, 3L, signs, FUN = "*")
    statistic <- apply(
      permuted,
      c(1, 2),
      function(values) {
        standard_deviation <- stats::sd(values)
        if (standard_deviation == 0) {
          NA_real_
        } else {
          mean(values) / (standard_deviation / 2)
        }
      }
    )
    reference_null[pattern + 1L] <- .ws913_reference_cluster_max(
      statistic,
      critical
    )
  }

  expect_identical(result$permutation_mode, "exact")
  expect_identical(result$n_permutations, 15L)
  expect_equal(result$null_max, reference_null, tolerance = 1e-12)
  if (nrow(result$cluster_table) > 0L) {
    reference_p <- vapply(
      result$cluster_table$mass,
      function(mass) (1 + sum(reference_null >= mass)) / 16,
      numeric(1)
    )
    expect_equal(result$cluster_table$p_value, reference_p, tolerance = 0)
  }
})

test_that("cluster helper detects a seeded rectangular effect", {
  set.seed(91312)
  tf <- array(
    rnorm(8L * 5L * 8L, sd = 0.4),
    dim = c(8L, 5L, 8L),
    dimnames = list(
      as.character(seq(-0.2, 0.5, by = 0.1)),
      as.character(c(4, 8, 12, 16, 20)),
      NULL
    )
  )
  tf[4:7, 2:4, ] <- tf[4:7, 2:4, ] + 3
  result <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = c(-0.2, -0.1),
    method = "cluster",
    n_permutations = 255L,
    alpha = 0.05
  )

  expect_true(all(result$mask[4:7, 2:4]))
  expect_false(result$mask[8L, 5L])
  expect_true(any(result$cluster_table$significant))
})

test_that("cluster helper preserves RNG and enforces attainable alpha", {
  set.seed(91313)
  tf <- array(
    rnorm(4L * 3L * 13L),
    dim = c(4L, 3L, 13L),
    dimnames = list(
      as.character(c(-0.2, -0.1, 0, 0.1)),
      as.character(c(5, 10, 15)),
      NULL
    )
  )
  set.seed(77)
  seed_before <- .Random.seed
  first <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = c(-0.2, -0.1),
    method = "cluster",
    n_permutations = 20L,
    seed = 41
  )
  expect_identical(.Random.seed, seed_before)
  second <- PhysioEEG:::.tfClusterMask(
    tf,
    baseline = c(-0.2, -0.1),
    method = "cluster",
    n_permutations = 20L,
    seed = 41
  )
  expect_identical(first, second)
  expect_true(all(
    first$cluster_table$p_value >= 1 / (first$n_permutations + 1)
  ))

  seed_before_error <- .Random.seed
  expect_error(
    PhysioEEG:::.tfClusterMask(
      tf,
      baseline = c(-0.2, -0.1),
      method = "cluster",
      n_permutations = 19L,
      seed = 99
    ),
    "at least 20"
  )
  expect_identical(.Random.seed, seed_before_error)
})

test_that("inference validation and non-estimable states are explicit", {
  tf <- array(
    seq_len(4L * 3L * 4L),
    dim = c(4L, 3L, 4L),
    dimnames = list(
      as.character(c(-0.2, -0.1, 0, 0.1)),
      as.character(c(5, 10, 15)),
      NULL
    )
  )
  constant <- array(
    rep(seq_len(12L), 4L),
    dim = c(4L, 3L, 4L),
    dimnames = dimnames(tf)
  )
  non_estimable <- PhysioEEG:::.tfClusterMask(
    constant,
    baseline = c(-0.2, -0.1),
    method = "threshold"
  )
  expect_true(any(!non_estimable$estimable))
  expect_true(all(is.na(
    non_estimable$statistic[!non_estimable$estimable]
  )))
  expect_true(all(!non_estimable$mask[!non_estimable$estimable]))
  expect_match(non_estimable$warnings, "non-estimable")

  expect_error(
    PhysioEEG:::.tfClusterMask(tf[, , 1:2], c(-0.2, -0.1)),
    "at least 2 x 2 x 3"
  )
  bad <- tf
  bad[1L] <- NA
  expect_error(PhysioEEG:::.tfClusterMask(bad, c(-0.2, -0.1)), "finite")
  expect_error(
    PhysioEEG:::.tfClusterMask(tf, c(-0.2, -0.1), connectivity = 8),
    "4-neighbour"
  )
  expect_error(
    PhysioEEG:::.tfClusterMask(tf, c(-0.2, -0.1), p_adjust_method = "bad"),
    "p.adjust"
  )
  expect_error(
    PhysioEEG:::.tfClusterMask(tf, logical(4L)),
    "at least two"
  )
  ambiguous_axis <- tf
  dimnames(ambiguous_axis)[[1L]] <- c("-0.2", "-0.1", "-0.1", "0.1")
  expect_error(
    PhysioEEG:::.tfClusterMask(ambiguous_axis, c(-0.2, -0.1)),
    "unique"
  )
})
