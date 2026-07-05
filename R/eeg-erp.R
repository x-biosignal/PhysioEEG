#' Detect ERP Components
#'
#' Detects known event-related potential (ERP) components in epoched (3D) EEG
#' data. Averages across epochs and finds peaks within predefined time windows.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data.
#' @param component ERP component to detect: \code{"N100"}, \code{"P300"},
#'   \code{"N400"}, \code{"P600"}, \code{"MMN"}, or \code{"LPP"}.
#' @param channels Character vector of channel labels to analyze. If \code{NULL},
#'   all channels are used.
#' @param epoch_start Start time of the epoch in milliseconds relative to
#'   stimulus onset. Default is 0 (epoch starts at stimulus).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with columns: \code{channel} (character label),
#'   \code{component} (character name), \code{latency_ms} (numeric peak
#'   latency in milliseconds), and \code{amplitude} (numeric peak amplitude).
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPmeasure()], [eegERPlatency()], [eegERPbaseline()],
#'   [eegEpoch()], [eegPlotERP()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, sr = 250)
#' result <- eegERPdetect(pe, component = "P300")
#' }
eegERPdetect <- function(x, component = c("N100", "P300", "N400", "P600", "MMN", "LPP"),
                         channels = NULL, epoch_start = 0, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  component <- match.arg(component)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Validate 3D data
  if (length(dim(data)) != 3) {
    stop("eegERPdetect requires epoched (3D) data (time x channels x epochs). ",
         "Use epoched data or a 3D assay.", call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_epochs <- dim(data)[3]

  # Average across epochs to get ERP waveform (time x channels)
  erp <- apply(data, c(1, 2), mean)

  # Get channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Filter channels if specified
  if (!is.null(channels)) {
    ch_idx <- which(ch_labels %in% channels)
    if (length(ch_idx) == 0) {
      stop("None of the specified channels were found.", call. = FALSE)
    }
  } else {
    ch_idx <- seq_len(n_channels)
  }

  # Define ERP component windows and polarity
  erp_defs <- list(
    N100 = list(window = c(80, 150), polarity = "negative"),
    P300 = list(window = c(250, 500), polarity = "positive"),
    N400 = list(window = c(300, 500), polarity = "negative"),
    P600 = list(window = c(500, 800), polarity = "positive"),
    MMN  = list(window = c(100, 250), polarity = "negative"),
    LPP  = list(window = c(400, 800), polarity = "positive")
  )

  comp_def <- erp_defs[[component]]
  window_ms <- comp_def$window
  polarity <- comp_def$polarity

  # Convert time window to sample indices, accounting for epoch_start
  start_sample <- max(1L, as.integer(floor((window_ms[1] - epoch_start) / 1000 * sr)) + 1L)
  end_sample <- min(n_time, as.integer(floor((window_ms[2] - epoch_start) / 1000 * sr)) + 1L)

  if (start_sample > n_time || end_sample < 1) {
    stop(sprintf("ERP window [%d, %d] ms is outside the epoch duration.",
                 window_ms[1], window_ms[2]), call. = FALSE)
  }

  results <- vector("list", length(ch_idx))
  for (i in seq_along(ch_idx)) {
    ch <- ch_idx[i]
    segment <- erp[start_sample:end_sample, ch]

    if (polarity == "positive") {
      peak_idx <- which.max(segment)
      peak_amp <- segment[peak_idx]
    } else {
      peak_idx <- which.min(segment)
      peak_amp <- segment[peak_idx]
    }

    # Convert sample index back to ms, accounting for epoch_start
    peak_sample <- start_sample + peak_idx - 1L
    latency_ms <- (peak_sample - 1) / sr * 1000 + epoch_start

    results[[i]] <- data.frame(
      channel = ch_labels[ch],
      component = component,
      latency_ms = latency_ms,
      amplitude = peak_amp,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Measure ERP Amplitude
#'
#' Measures ERP amplitude in a specified time window using peak, mean, or
#' adaptive mean methods. For epoched (3D) data, averages across epochs first.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param window Numeric vector of length 2: \code{c(start_ms, end_ms)}.
#' @param method Measurement method: \code{"peak"} (peak amplitude),
#'   \code{"mean"} (mean amplitude), or \code{"adaptive_mean"} (mean in
#'   +/-25ms around peak).
#' @param polarity Expected polarity: \code{"positive"} or \code{"negative"}.
#' @param epoch_start Start time of the epoch in milliseconds relative to
#'   stimulus onset. Default is 0 (epoch starts at stimulus).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with columns: \code{channel} (character label),
#'   \code{amplitude} (numeric measured amplitude), \code{latency_ms}
#'   (numeric latency in milliseconds), and \code{method} (character
#'   measurement method used).
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPdetect()], [eegERPlatency()], [eegERPbaseline()],
#'   [eegPlotERP()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, sr = 250)
#' result <- eegERPmeasure(pe, window = c(250, 500), method = "peak",
#'                         polarity = "positive")
#' }
eegERPmeasure <- function(x, window, method = c("peak", "mean", "adaptive_mean"),
                          polarity = c("positive", "negative"),
                          epoch_start = 0, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  polarity <- match.arg(polarity)
  stopifnot(is.numeric(window) && length(window) == 2)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Average across epochs if 3D
  if (length(dim(data)) == 3) {
    erp <- apply(data, c(1, 2), mean)
  } else {
    erp <- data
  }

  n_time <- nrow(erp)
  n_channels <- ncol(erp)

  # Get channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Convert window to sample indices, accounting for epoch_start
  start_sample <- max(1L, as.integer(floor((window[1] - epoch_start) / 1000 * sr)) + 1L)
  end_sample <- min(n_time, as.integer(floor((window[2] - epoch_start) / 1000 * sr)) + 1L)

  results <- vector("list", n_channels)
  for (ch in seq_len(n_channels)) {
    segment <- erp[start_sample:end_sample, ch]

    if (method == "peak") {
      if (polarity == "positive") {
        peak_idx <- which.max(segment)
      } else {
        peak_idx <- which.min(segment)
      }
      amplitude <- segment[peak_idx]
      peak_sample <- start_sample + peak_idx - 1L
      latency_ms <- (peak_sample - 1) / sr * 1000 + epoch_start

    } else if (method == "mean") {
      amplitude <- mean(segment)
      # Latency is midpoint of window for mean method
      latency_ms <- mean(window)

    } else if (method == "adaptive_mean") {
      # Find peak first
      if (polarity == "positive") {
        peak_idx <- which.max(segment)
      } else {
        peak_idx <- which.min(segment)
      }
      peak_sample <- start_sample + peak_idx - 1L
      latency_ms <- (peak_sample - 1) / sr * 1000 + epoch_start

      # Mean in +/-25ms around peak
      half_win_samples <- as.integer(round(25 / 1000 * sr))
      adapt_start <- max(1L, peak_sample - half_win_samples)
      adapt_end <- min(n_time, peak_sample + half_win_samples)
      amplitude <- mean(erp[adapt_start:adapt_end, ch])
    }

    results[[ch]] <- data.frame(
      channel = ch_labels[ch],
      amplitude = amplitude,
      latency_ms = latency_ms,
      method = method,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Fractional Area Latency
#'
#' Computes the fractional area latency of an ERP component. This is the time
#' point at which a specified fraction of the total area under the curve (in
#' the given window) has accumulated.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param window Numeric vector of length 2: \code{c(start_ms, end_ms)}.
#' @param fraction Fraction of total area (default: 0.5 for median latency).
#' @param polarity Expected polarity: \code{"positive"} or \code{"negative"}.
#' @param epoch_start Start time of the epoch in milliseconds relative to
#'   stimulus onset. Default is 0 (epoch starts at stimulus).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with columns: \code{channel} (character label),
#'   \code{latency_ms} (numeric fractional area latency in milliseconds),
#'   and \code{fraction} (numeric fraction used).
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPdetect()], [eegERPmeasure()], [eegERPbaseline()],
#'   [eegPlotERP()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, sr = 250)
#' result <- eegERPlatency(pe, window = c(250, 500), fraction = 0.5,
#'                         polarity = "positive")
#' }
eegERPlatency <- function(x, window, fraction = 0.5,
                          polarity = c("positive", "negative"),
                          epoch_start = 0, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  polarity <- match.arg(polarity)
  stopifnot(is.numeric(window) && length(window) == 2)
  stopifnot(is.numeric(fraction) && fraction > 0 && fraction <= 1)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Average across epochs if 3D
  if (length(dim(data)) == 3) {
    erp <- apply(data, c(1, 2), mean)
  } else {
    erp <- data
  }

  n_time <- nrow(erp)
  n_channels <- ncol(erp)

  # Get channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Convert window to sample indices, accounting for epoch_start
  start_sample <- max(1L, as.integer(floor((window[1] - epoch_start) / 1000 * sr)) + 1L)
  end_sample <- min(n_time, as.integer(floor((window[2] - epoch_start) / 1000 * sr)) + 1L)

  results <- vector("list", n_channels)
  for (ch in seq_len(n_channels)) {
    segment <- erp[start_sample:end_sample, ch]

    # Rectify based on polarity
    if (polarity == "positive") {
      rectified <- pmax(segment, 0)
    } else {
      rectified <- pmax(-segment, 0)
    }

    # Cumulative sum
    cs <- cumsum(rectified)
    total_area <- cs[length(cs)]

    if (total_area <= 0) {
      # No area in window matching polarity; return midpoint
      latency_ms <- mean(window)
    } else {
      target <- fraction * total_area
      # Find first index where cumsum >= target
      lat_idx <- which(cs >= target)[1]
      if (is.na(lat_idx)) lat_idx <- length(cs)
      lat_sample <- start_sample + lat_idx - 1L
      latency_ms <- (lat_sample - 1) / sr * 1000 + epoch_start
    }

    results[[ch]] <- data.frame(
      channel = ch_labels[ch],
      latency_ms = latency_ms,
      fraction = fraction,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' ERP Difference Waveform
#'
#' Computes the difference waveform between two PhysioExperiment objects by
#' subtracting the assay data of \code{y} from \code{x}.
#'
#' @param x A PhysioExperiment object (minuend).
#' @param y A PhysioExperiment object (subtrahend). Must have the same
#'   dimensions as \code{x}.
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: \code{"difference"}).
#' @return Modified \code{x} with difference waveform stored in
#'   \code{output_assay}. The difference assay has the same dimensions as
#'   the input assay.
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPdetect()], [eegERPmeasure()], [eegERPtest()],
#'   [eegERPgrandAverage()]
#' @export
#' @examples
#' \dontrun{
#' pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' result <- eegERPdifference(pe1, pe2)
#' }
eegERPdifference <- function(x, y, assay_name = NULL,
                             output_assay = "difference") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(inherits(y, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data_x <- SummarizedExperiment::assay(x, assay_name)
  data_y <- SummarizedExperiment::assay(y, assay_name)

  # Validate dimensions match
  if (!identical(dim(data_x), dim(data_y))) {
    stop("Dimensions of x and y must match. ",
         sprintf("x: %s, y: %s",
                 paste(dim(data_x), collapse = " x "),
                 paste(dim(data_y), collapse = " x ")),
         call. = FALSE)
  }

  diff_data <- data_x - data_y

  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- diff_data
  S4Vectors::metadata(x) <- md

  x
}


#' Baseline Correct ERP Data
#'
#' Subtracts the mean amplitude of a pre-stimulus baseline period from each
#' epoch. This is an essential preprocessing step before ERP measurement.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data.
#' @param baseline Baseline period in milliseconds as \code{c(start, end)}.
#'   For example, \code{c(-200, 0)} for 200ms pre-stimulus baseline.
#'   Default is \code{c(-200, 0)}.
#' @param epoch_start Start time of the epoch in milliseconds relative to
#'   stimulus onset. Default is 0 (epoch starts at stimulus).
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: \code{"baseline_corrected"}).
#' @return Modified PhysioExperiment with baseline-corrected data in
#'   \code{output_assay}. Creates a new 3D assay (time x channels x epochs)
#'   where each epoch has had the mean of its baseline period subtracted.
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPdetect()], [eegERPmeasure()], [eegEpoch()],
#'   [eegFilter()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, sr = 250)
#' result <- eegERPbaseline(pe, baseline = c(-200, 0), epoch_start = -200)
#' }
eegERPbaseline <- function(x, baseline = c(-200, 0), epoch_start = 0,
                           assay_name = NULL, output_assay = "baseline_corrected") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(baseline) && length(baseline) == 2)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Validate 3D data
  if (length(dim(data)) != 3) {
    stop("eegERPbaseline requires epoched (3D) data (time x channels x epochs).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_epochs <- dim(data)[3]

  # Convert baseline window from ms to sample indices, accounting for epoch_start
  bl_start <- max(1L, as.integer(floor((baseline[1] - epoch_start) / 1000 * sr)) + 1L)
  bl_end <- min(n_time, as.integer(floor((baseline[2] - epoch_start) / 1000 * sr)) + 1L)

  if (bl_start > bl_end) {
    stop(sprintf("Invalid baseline window: samples %d to %d. ",
                 bl_start, bl_end),
         "Check that baseline and epoch_start are consistent.", call. = FALSE)
  }

  # Baseline correct: subtract mean of baseline period per channel per epoch
  corrected <- data
  for (ep in seq_len(n_epochs)) {
    for (ch in seq_len(n_channels)) {
      bl_mean <- mean(data[bl_start:bl_end, ch, ep], na.rm = TRUE)
      corrected[, ch, ep] <- data[, ch, ep] - bl_mean
    }
  }

  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- corrected
  S4Vectors::metadata(x) <- md

  x
}


#' Statistical Testing for ERP Differences
#'
#' Performs permutation testing or cluster-based permutation testing to
#' compare ERP waveforms between conditions.
#'
#' @param x A PhysioExperiment with epoched data for condition 1.
#' @param y A PhysioExperiment with epoched data for condition 2.
#' @param method Test method: \code{"permutation"} (pointwise permutation test)
#'   or \code{"cluster"} (cluster-based permutation test).
#' @param n_perm Number of permutations (default: 1000).
#' @param alpha Significance level (default: 0.05).
#' @param cluster_alpha Cluster-forming threshold for individual t-tests
#'   (default: 0.05). Only used for \code{"cluster"} method.
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with columns: \code{time_sample} (integer),
#'   \code{t_statistic} (numeric observed t-value), \code{p_value}
#'   (numeric permutation-based p-value), and \code{significant} (logical).
#'   For the \code{"cluster"} method, also includes \code{cluster_id}
#'   (integer cluster assignment) and \code{cluster_p} (numeric cluster-level
#'   corrected p-value).
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#'
#' Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
#' EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.
#' @seealso [eegERPdetect()], [eegERPmeasure()], [eegERPdifference()],
#'   [eegERPgrandAverage()]
#' @export
#' @examples
#' \dontrun{
#' pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' result <- eegERPtest(pe1, pe2, method = "permutation", n_perm = 500)
#' }
eegERPtest <- function(x, y, method = c("permutation", "cluster"),
                       n_perm = 1000, alpha = 0.05, cluster_alpha = 0.05,
                       assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(inherits(y, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(n_perm) && n_perm >= 1)
  stopifnot(is.numeric(alpha) && alpha > 0 && alpha < 1)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data_x <- SummarizedExperiment::assay(x, assay_name)
  data_y <- SummarizedExperiment::assay(y, assay_name)

  # Validate 3D data
  if (length(dim(data_x)) != 3 || length(dim(data_y)) != 3) {
    stop("eegERPtest requires epoched (3D) data for both x and y.",
         call. = FALSE)
  }

  n_time <- dim(data_x)[1]
  n_channels_x <- dim(data_x)[2]
  n_epochs_x <- dim(data_x)[3]
  n_epochs_y <- dim(data_y)[3]

  if (dim(data_x)[1] != dim(data_y)[1]) {
    stop("x and y must have the same number of time points.", call. = FALSE)
  }

  # Average across channels to get single-channel ERP per epoch
  # Shape: time x epochs for each condition
  erp_x <- apply(data_x, c(1, 3), mean)  # time x epochs

  erp_y <- apply(data_y, c(1, 3), mean)  # time x epochs

  n_total <- n_epochs_x + n_epochs_y
  combined <- cbind(erp_x, erp_y)  # time x (n_epochs_x + n_epochs_y)

  # Compute observed t-statistics at each time point
  obs_t <- numeric(n_time)
  for (tp in seq_len(n_time)) {
    vals_x <- erp_x[tp, ]
    vals_y <- erp_y[tp, ]
    mean_diff <- mean(vals_x) - mean(vals_y)
    se <- sqrt(stats::var(vals_x) / n_epochs_x + stats::var(vals_y) / n_epochs_y)
    obs_t[tp] <- if (se > 0) mean_diff / se else 0
  }

  if (method == "permutation") {
    # Pointwise permutation test
    null_t <- matrix(0, nrow = n_perm, ncol = n_time)

    for (perm in seq_len(n_perm)) {
      perm_idx <- sample(n_total)
      perm_x <- combined[, perm_idx[seq_len(n_epochs_x)], drop = FALSE]
      perm_y <- combined[, perm_idx[(n_epochs_x + 1):n_total], drop = FALSE]

      for (tp in seq_len(n_time)) {
        vals_x <- perm_x[tp, ]
        vals_y <- perm_y[tp, ]
        mean_diff <- mean(vals_x) - mean(vals_y)
        se <- sqrt(stats::var(vals_x) / n_epochs_x + stats::var(vals_y) / n_epochs_y)
        null_t[perm, tp] <- if (se > 0) mean_diff / se else 0
      }
    }

    # p-value = proportion of null |t| >= observed |t|
    p_values <- numeric(n_time)
    for (tp in seq_len(n_time)) {
      p_values[tp] <- (sum(abs(null_t[, tp]) >= abs(obs_t[tp])) + 1) / (n_perm + 1)
    }

    result <- data.frame(
      time_sample = seq_len(n_time),
      t_statistic = obs_t,
      p_value = p_values,
      significant = p_values < alpha,
      stringsAsFactors = FALSE
    )

  } else {
    # Cluster-based permutation test

    # Helper: find clusters of consecutive significant time points
    find_clusters <- function(t_vals, thresh_alpha, n_x, n_y) {
      df <- n_x + n_y - 2
      t_crit <- stats::qt(1 - thresh_alpha / 2, df = df)
      sig <- abs(t_vals) > t_crit

      clusters <- list()
      cluster_id <- integer(length(t_vals))
      current_cluster <- 0L

      in_cluster <- FALSE
      for (i in seq_along(sig)) {
        if (sig[i]) {
          if (!in_cluster) {
            current_cluster <- current_cluster + 1L
            in_cluster <- TRUE
          }
          cluster_id[i] <- current_cluster
        } else {
          in_cluster <- FALSE
        }
      }

      # Compute cluster statistics (sum of t-values in each cluster)
      if (current_cluster > 0) {
        cluster_stats <- numeric(current_cluster)
        for (cl in seq_len(current_cluster)) {
          cluster_stats[cl] <- sum(t_vals[cluster_id == cl])
        }
      } else {
        cluster_stats <- numeric(0)
      }

      list(cluster_id = cluster_id, cluster_stats = cluster_stats,
           n_clusters = current_cluster)
    }

    # Observed clusters
    obs_clusters <- find_clusters(obs_t, cluster_alpha, n_epochs_x, n_epochs_y)

    # Build null distribution of max cluster statistics
    null_max_cluster <- numeric(n_perm)

    for (perm in seq_len(n_perm)) {
      perm_idx <- sample(n_total)
      perm_x <- combined[, perm_idx[seq_len(n_epochs_x)], drop = FALSE]
      perm_y <- combined[, perm_idx[(n_epochs_x + 1):n_total], drop = FALSE]

      perm_t <- numeric(n_time)
      for (tp in seq_len(n_time)) {
        vals_x <- perm_x[tp, ]
        vals_y <- perm_y[tp, ]
        mean_diff <- mean(vals_x) - mean(vals_y)
        se <- sqrt(stats::var(vals_x) / n_epochs_x + stats::var(vals_y) / n_epochs_y)
        perm_t[tp] <- if (se > 0) mean_diff / se else 0
      }

      perm_clusters <- find_clusters(perm_t, cluster_alpha, n_epochs_x, n_epochs_y)
      if (perm_clusters$n_clusters > 0) {
        null_max_cluster[perm] <- max(abs(perm_clusters$cluster_stats))
      } else {
        null_max_cluster[perm] <- 0
      }
    }

    # Compute cluster p-values
    cluster_p <- numeric(n_time)
    cluster_p[] <- NA_real_
    if (obs_clusters$n_clusters > 0) {
      for (cl in seq_len(obs_clusters$n_clusters)) {
        obs_stat <- abs(obs_clusters$cluster_stats[cl])
        cl_p <- (sum(null_max_cluster >= obs_stat) + 1) / (n_perm + 1)
        cluster_p[obs_clusters$cluster_id == cl] <- cl_p
      }
    }

    # Pointwise p-values not meaningful for cluster method, set to NA for
    # non-cluster time points
    p_values <- cluster_p

    result <- data.frame(
      time_sample = seq_len(n_time),
      t_statistic = obs_t,
      p_value = ifelse(is.na(p_values), 1, p_values),
      significant = ifelse(is.na(p_values), FALSE, p_values < alpha),
      cluster_id = obs_clusters$cluster_id,
      cluster_p = cluster_p,
      stringsAsFactors = FALSE
    )
  }

  result
}


#' Compute Grand Average ERP
#'
#' Averages ERP waveforms across multiple PhysioExperiment objects
#' (participants). Each input should already be averaged across trials.
#'
#' @param ... PhysioExperiment objects to average, or a list of them.
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: \code{"grand_average"}).
#' @return A PhysioExperiment (the first input object) with the grand average
#'   waveform stored in \code{output_assay}. The grand average assay has the
#'   same dimensions as the input assays.
#' @references
#' Luck, S. J. (2014). An Introduction to the Event-Related Potential Technique
#' (2nd ed.). MIT Press.
#' @seealso [eegERPdetect()], [eegERPmeasure()], [eegERPtest()],
#'   [eegERPdifference()]
#' @export
#' @examples
#' \dontrun{
#' pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
#' result <- eegERPgrandAverage(pe1, pe2)
#' }
eegERPgrandAverage <- function(..., assay_name = NULL,
                               output_assay = "grand_average") {
  args <- list(...)

  # If a single list was passed, unpack it
  if (length(args) == 1 && is.list(args[[1]]) &&
      !inherits(args[[1]], "PhysioExperiment")) {
    args <- args[[1]]
  }

  if (length(args) < 2) {
    stop("eegERPgrandAverage requires at least 2 PhysioExperiment objects.",
         call. = FALSE)
  }

  # Validate all inputs
  for (i in seq_along(args)) {
    stopifnot(inherits(args[[i]], "PhysioExperiment"))
  }

  # Get assay data from each object
  if (is.null(assay_name)) assay_name <- defaultAssay(args[[1]])

  assay_list <- lapply(args, function(pe) {
    SummarizedExperiment::assay(pe, assay_name)
  })

  # Validate all have the same dimensions
  ref_dim <- dim(assay_list[[1]])
  for (i in seq_along(assay_list)[-1]) {
    if (!identical(dim(assay_list[[i]]), ref_dim)) {
      stop(sprintf(
        "Dimension mismatch: object 1 has dims [%s] but object %d has dims [%s].",
        paste(ref_dim, collapse = " x "), i,
        paste(dim(assay_list[[i]]), collapse = " x ")),
        call. = FALSE)
    }
  }

  # Compute grand average
  n_subjects <- length(assay_list)
  grand_avg <- assay_list[[1]]
  for (i in seq_along(assay_list)[-1]) {
    grand_avg <- grand_avg + assay_list[[i]]
  }
  grand_avg <- grand_avg / n_subjects

  # Store in first PE object
  result <- args[[1]]
  md <- S4Vectors::metadata(result)
  md[[output_assay]] <- grand_avg
  S4Vectors::metadata(result) <- md

  result
}
