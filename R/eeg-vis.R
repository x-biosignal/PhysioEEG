#' EEG Visualization Functions
#'
#' Comprehensive visualization functions for EEG data stored in
#' PhysioExperiment objects. All functions return ggplot2 objects and
#' require ggplot2 to be installed (Suggests dependency).
#'
#' @name eeg-vis
#' @keywords internal
NULL


# --- Internal helpers ---

#' Check ggplot2 availability
#'
#' Verifies that ggplot2 is installed. Called at the top of every
#' visualization function since ggplot2 is in Suggests, not Imports.
#'
#' @return Invisible NULL. Throws an error if ggplot2 is not installed.
#' @keywords internal
.check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
}


#' Inverse distance weighted interpolation for topographic maps
#'
#' Interpolates channel values onto a regular grid within a unit circle
#' using inverse distance weighting with power 2. Points outside the
#' head circle are set to NA.
#'
#' @param values Numeric vector of channel values.
#' @param x_pos Numeric vector of channel x positions.
#' @param y_pos Numeric vector of channel y positions.
#' @param resolution Integer grid resolution (default: 100).
#' @return A data.frame with columns x, y, value suitable for
#'   ggplot2::geom_raster.
#' @keywords internal
.interpolate_topomap <- function(values, x_pos, y_pos, resolution = 100) {
  # Create regular grid within unit circle
  grid_x <- seq(-1, 1, length.out = resolution)
  grid_y <- seq(-1, 1, length.out = resolution)
  grid <- expand.grid(x = grid_x, y = grid_y)

  # Only interpolate within the head circle
  dist_from_center <- sqrt(grid$x^2 + grid$y^2)
  inside <- dist_from_center <= 1.0

  grid$value <- NA_real_

  n_ch <- length(values)
  power <- 2  # IDW power parameter

  for (i in which(inside)) {
    gx <- grid$x[i]
    gy <- grid$y[i]

    # Distances from grid point to all electrodes
    dists <- sqrt((x_pos - gx)^2 + (y_pos - gy)^2)

    # Check for exact electrode location
    exact <- which(dists < 1e-10)
    if (length(exact) > 0) {
      grid$value[i] <- values[exact[1]]
    } else {
      weights <- 1 / (dists^power)
      grid$value[i] <- sum(weights * values) / sum(weights)
    }
  }

  grid
}


#' Generate head outline coordinates
#'
#' Returns coordinates for a unit circle representing the head outline.
#'
#' @param n_points Number of points on the circle (default: 100).
#' @return A data.frame with columns x, y.
#' @keywords internal
.head_outline <- function(n_points = 100) {
  theta <- seq(0, 2 * pi, length.out = n_points + 1)
  data.frame(x = cos(theta), y = sin(theta))
}


#' Generate nose shape coordinates
#'
#' Returns coordinates for a triangle at the top of the head circle
#' representing the nose.
#'
#' @return A data.frame with columns x, y.
#' @keywords internal
.nose_shape <- function() {
  data.frame(
    x = c(-0.1, 0, 0.1, -0.1),
    y = c(1.0, 1.15, 1.0, 1.0)
  )
}


#' Generate ear shape coordinates
#'
#' Returns coordinates for a stylized ear on the specified side of the
#' head circle.
#'
#' @param side Character, either "left" or "right".
#' @return A data.frame with columns x, y.
#' @keywords internal
.ear_shape <- function(side = c("left", "right")) {
  side <- match.arg(side)
  # Approximate ear as a small elliptical bump
  t <- seq(-0.5, 0.5, length.out = 20)
  ear_y <- t * 0.4
  ear_x <- 0.08 * cos(t * pi)

  if (side == "left") {
    data.frame(x = -1.0 - ear_x, y = ear_y)
  } else {
    data.frame(x = 1.0 + ear_x, y = ear_y)
  }
}


#' Default 10-20 electrode positions (2D projection)
#'
#' Returns standard 10-20 system electrode positions projected onto
#' a unit circle for topographic mapping.
#'
#' @return A named list with pos_x and pos_y numeric vectors, with names
#'   corresponding to 10-20 labels.
#' @keywords internal
.default_1020_positions <- function() {
  labels <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
              "T3", "C3", "Cz", "C4", "T4",
              "T5", "P3", "Pz", "P4", "T6",
              "O1", "O2")
  # 2D polar coordinates (azimuth, radius fraction)
  # Arranged based on standard 10-20 polar projections
  pos_x <- c(-0.31, 0.31, -0.81, -0.39, 0.00, 0.39, 0.81,
             -0.95, -0.45, 0.00, 0.45, 0.95,
             -0.81, -0.39, 0.00, 0.39, 0.81,
             -0.31, 0.31)
  pos_y <- c(0.85, 0.85, 0.45, 0.45, 0.45, 0.45, 0.45,
             0.00, 0.00, 0.00, 0.00, 0.00,
             -0.45, -0.45, -0.45, -0.45, -0.45,
             -0.85, -0.85)
  names(pos_x) <- labels
  names(pos_y) <- labels
  list(pos_x = pos_x, pos_y = pos_y)
}


# =====================================================================
# Exported Visualization Functions
# =====================================================================


#' Plot Multi-Channel EEG Time Series
#'
#' Displays multi-channel EEG time series in stacked, butterfly, or grid
#' layout. Supports channel selection, time range restriction, and
#' optional event markers.
#'
#' @param x A PhysioExperiment object with 2D EEG data (time x channels).
#' @param channels Character vector of channel labels to display.
#'   If \code{NULL}, all channels are shown.
#' @param time_range Numeric vector of length 2 specifying time range in
#'   seconds as \code{c(start, end)}. If \code{NULL}, all data is shown.
#' @param mode Display mode: \code{"stacked"} (vertical offset per channel),
#'   \code{"butterfly"} (overlay all channels), or \code{"grid"}
#'   (faceted panels per channel).
#' @param scale Numeric scaling factor for vertical offset in stacked mode
#'   (default: 1).
#' @param show_events Logical; if \code{TRUE} and \code{metadata(x)$events}
#'   exists, vertical lines are drawn at event times.
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 2500, n_channels = 4, sr = 500)
#' eegPlotSignal(pe, mode = "stacked")
#' }
eegPlotSignal <- function(x, channels = NULL, time_range = NULL,
                          mode = c("stacked", "butterfly", "grid"),
                          scale = 1, show_events = FALSE,
                          assay_name = NULL) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))
  mode <- match.arg(mode)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Handle 2D data
  if (length(dim(data)) == 3) {
    data <- data[, , 1, drop = TRUE]
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Filter channels
  if (!is.null(channels)) {
    ch_idx <- which(ch_labels %in% channels)
    if (length(ch_idx) == 0) {
      stop("None of the specified channels were found.", call. = FALSE)
    }
  } else {
    ch_idx <- seq_len(n_channels)
  }

  # Time vector in seconds
  time_sec <- (seq_len(n_time) - 1) / sr

  # Apply time range filter
  if (!is.null(time_range)) {
    stopifnot(is.numeric(time_range) && length(time_range) == 2)
    time_idx <- which(time_sec >= time_range[1] & time_sec <= time_range[2])
    if (length(time_idx) == 0) {
      stop("No data points within the specified time_range.", call. = FALSE)
    }
  } else {
    time_idx <- seq_len(n_time)
  }

  # Build long-format data frame
  plot_list <- vector("list", length(ch_idx))
  for (i in seq_along(ch_idx)) {
    ch <- ch_idx[i]
    vals <- data[time_idx, ch]

    if (mode == "stacked") {
      # Offset each channel vertically
      offset <- (length(ch_idx) - i) * scale * stats::sd(data[, ch], na.rm = TRUE) * 3
      vals <- vals + offset
    }

    plot_list[[i]] <- data.frame(
      time = time_sec[time_idx],
      amplitude = vals,
      channel = ch_labels[ch],
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, plot_list)

  # Maintain channel order
  plot_df$channel <- factor(plot_df$channel, levels = ch_labels[ch_idx])

  # Build plot based on mode
  if (mode == "stacked") {
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$time, y = .data$amplitude,
                                      color = .data$channel)) +
      ggplot2::geom_line(linewidth = 0.3) +
      ggplot2::labs(x = "Time (s)", y = "Amplitude", color = "Channel",
                    title = "EEG Signal (Stacked)") +
      ggplot2::theme_minimal()

  } else if (mode == "butterfly") {
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$time, y = .data$amplitude,
                                      color = .data$channel,
                                      group = .data$channel)) +
      ggplot2::geom_line(linewidth = 0.3, alpha = 0.7) +
      ggplot2::labs(x = "Time (s)", y = "Amplitude", color = "Channel",
                    title = "EEG Signal (Butterfly)") +
      ggplot2::theme_minimal()

  } else {
    # grid mode
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$time, y = .data$amplitude)) +
      ggplot2::geom_line(linewidth = 0.3) +
      ggplot2::facet_wrap(~ channel, scales = "free_y") +
      ggplot2::labs(x = "Time (s)", y = "Amplitude",
                    title = "EEG Signal (Grid)") +
      ggplot2::theme_minimal()
  }

  # Add event markers if requested
  if (show_events) {
    events <- S4Vectors::metadata(x)$events
    if (!is.null(events) && is.data.frame(events)) {
      # Expect events to have a 'time' column in seconds
      event_col <- if ("time" %in% names(events)) "time"
                   else if ("time_sec" %in% names(events)) "time_sec"
                   else if ("onset" %in% names(events)) "onset"
                   else NULL
      if (!is.null(event_col)) {
        event_times <- events[[event_col]]
        # Filter to current time range
        if (!is.null(time_range)) {
          event_times <- event_times[event_times >= time_range[1] &
                                     event_times <= time_range[2]]
        }
        if (length(event_times) > 0) {
          event_df <- data.frame(xint = event_times)
          p <- p + ggplot2::geom_vline(
            data = event_df,
            ggplot2::aes(xintercept = .data$xint),
            linetype = "dashed", color = "red", alpha = 0.5
          )
        }
      }
    }
  }

  p
}


#' Plot ERP Waveform with Confidence Interval
#'
#' Plots event-related potential waveforms averaged across epochs, with
#' optional confidence interval ribbons. Supports condition-based
#' comparisons when \code{metadata(x)$conditions} is available.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x epochs).
#' @param channels Character vector of channel labels to plot. If \code{NULL},
#'   the grand average across all channels is used.
#' @param conditions Character vector of condition labels to include.
#'   If \code{NULL} and \code{metadata(x)$conditions} exists, all conditions
#'   are plotted.
#' @param ci Confidence level for the interval (default: 0.95).
#' @param show_ci Logical; if \code{TRUE}, display confidence interval ribbon.
#' @param epoch_start Numeric start time of each epoch in seconds for
#'   the x-axis (default: 0).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, sr = 250)
#' eegPlotERP(pe, channels = "Cz")
#' }
eegPlotERP <- function(x, channels = NULL, conditions = NULL,
                       ci = 0.95, show_ci = TRUE,
                       epoch_start = 0, assay_name = NULL) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 3) {
    stop("eegPlotERP requires epoched (3D) data (time x channels x epochs).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_epochs <- dim(data)[3]

  # Channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Filter channels
  if (!is.null(channels)) {
    ch_idx <- which(ch_labels %in% channels)
    if (length(ch_idx) == 0) {
      stop("None of the specified channels were found.", call. = FALSE)
    }
  } else {
    ch_idx <- seq_len(n_channels)
  }

  # Time axis in seconds
  time_sec <- epoch_start + (seq_len(n_time) - 1) / sr

  # Get conditions from metadata
  cond_vec <- S4Vectors::metadata(x)$conditions
  has_conditions <- !is.null(cond_vec) && length(cond_vec) == n_epochs

  if (has_conditions) {
    if (!is.null(conditions)) {
      cond_levels <- conditions
    } else {
      cond_levels <- unique(cond_vec)
    }
  } else {
    cond_levels <- "all"
    cond_vec <- rep("all", n_epochs)
  }

  # CI z-value
  z_val <- stats::qnorm((1 + ci) / 2)

  # Build plot data
  plot_list <- list()
  for (cond in cond_levels) {
    epoch_mask <- which(cond_vec == cond)
    if (length(epoch_mask) == 0) next

    # Average across selected channels and epochs within condition
    # First subset epochs: data[, ch_idx, epoch_mask]
    if (length(ch_idx) == 1 && length(epoch_mask) == 1) {
      # Single channel, single epoch
      erp_mean <- data[, ch_idx, epoch_mask]
      erp_se <- rep(0, n_time)
    } else {
      # Average across selected channels first
      if (length(ch_idx) == 1) {
        ch_data <- data[, ch_idx, epoch_mask, drop = FALSE]
        dim(ch_data) <- c(n_time, length(epoch_mask))
      } else {
        # Average across channels: time x epochs
        ch_data <- apply(data[, ch_idx, epoch_mask, drop = FALSE], c(1, 3), mean)
      }
      if (length(epoch_mask) == 1) {
        erp_mean <- as.numeric(ch_data)
        erp_se <- rep(0, n_time)
      } else {
        erp_mean <- rowMeans(ch_data)
        erp_sd <- apply(ch_data, 1, stats::sd)
        erp_se <- erp_sd / sqrt(length(epoch_mask))
      }
    }

    ci_lower <- erp_mean - z_val * erp_se
    ci_upper <- erp_mean + z_val * erp_se

    plot_list[[length(plot_list) + 1]] <- data.frame(
      time = time_sec,
      amplitude = erp_mean,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      condition = cond,
      stringsAsFactors = FALSE
    )
  }

  plot_df <- do.call(rbind, plot_list)
  plot_df$condition <- factor(plot_df$condition, levels = cond_levels)

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$time, y = .data$amplitude,
                                    color = .data$condition,
                                    fill = .data$condition)) +
    ggplot2::geom_line(linewidth = 0.8)

  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      alpha = 0.2, color = NA
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::labs(x = "Time (s)", y = "Amplitude",
                  color = "Condition", fill = "Condition",
                  title = "Event-Related Potential") +
    ggplot2::theme_minimal()

  p
}


#' Plot EEG Topographic Map
#'
#' Creates a 2D scalp topographic map using inverse distance weighted
#' interpolation. Electrode positions are read from \code{colData(x)}
#' (columns \code{pos_x}, \code{pos_y}) or default 10-20 positions are
#' used as fallback. The plot includes a head outline, nose, and ears.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param time Numeric time point in seconds at which to extract values
#'   from the assay. If \code{NULL} and \code{values} is also \code{NULL},
#'   the mean across all time points is used.
#' @param values Named numeric vector of channel values to plot directly.
#'   If provided, overrides data extraction from the assay.
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @param resolution Integer grid resolution for interpolation (default: 100).
#' @param palette Character name of the diverging color palette
#'   (default: \code{"RdBu"}).
#' @param contours Logical; if \code{TRUE}, add contour lines.
#' @param electrodes Logical; if \code{TRUE}, show electrode positions as points.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
#' eegPlotTopomap(pe, time = 0.5)
#' }
eegPlotTopomap <- function(x, time = NULL, values = NULL,
                           assay_name = NULL, resolution = 100,
                           palette = "RdBu", contours = TRUE,
                           electrodes = TRUE) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))

  col_data <- SummarizedExperiment::colData(x)
  n_channels <- nrow(col_data)

  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Determine electrode positions
  if (all(c("pos_x", "pos_y") %in% colnames(col_data))) {
    pos_x <- as.numeric(col_data$pos_x)
    pos_y <- as.numeric(col_data$pos_y)
  } else {
    # Use default 10-20 positions
    defaults <- .default_1020_positions()
    pos_x <- numeric(n_channels)
    pos_y <- numeric(n_channels)
    for (i in seq_len(n_channels)) {
      lbl <- ch_labels[i]
      if (lbl %in% names(defaults$pos_x)) {
        pos_x[i] <- defaults$pos_x[lbl]
        pos_y[i] <- defaults$pos_y[lbl]
      } else {
        # Spread unrecognized channels in a circle
        angle <- 2 * pi * i / n_channels
        pos_x[i] <- 0.5 * cos(angle)
        pos_y[i] <- 0.5 * sin(angle)
      }
    }
  }

  # Determine values to plot
  if (is.null(values)) {
    if (is.null(assay_name)) assay_name <- defaultAssay(x)
    data <- SummarizedExperiment::assay(x, assay_name)

    # Handle 3D data (use first epoch)
    if (length(dim(data)) == 3) {
      data <- apply(data, c(1, 2), mean)
    }

    sr <- samplingRate(x)

    if (!is.null(time)) {
      sample_idx <- max(1L, min(nrow(data),
                                as.integer(round(time * sr)) + 1L))
      ch_values <- data[sample_idx, ]
    } else {
      ch_values <- colMeans(data)
    }
  } else {
    # values provided directly
    if (!is.null(names(values))) {
      ch_values <- numeric(n_channels)
      for (i in seq_len(n_channels)) {
        if (ch_labels[i] %in% names(values)) {
          ch_values[i] <- values[ch_labels[i]]
        }
      }
    } else {
      ch_values <- values[seq_len(n_channels)]
    }
  }

  # Interpolate
  grid_df <- .interpolate_topomap(ch_values, pos_x, pos_y,
                                  resolution = resolution)

  # Color limits (symmetric around zero for diverging palette)
  max_abs <- max(abs(grid_df$value), na.rm = TRUE)
  if (is.na(max_abs) || max_abs == 0) max_abs <- 1

  # Head outline, nose, ears
  head_df <- .head_outline()
  nose_df <- .nose_shape()
  ear_left <- .ear_shape("left")
  ear_right <- .ear_shape("right")

  # Electrode positions data frame
  elec_df <- data.frame(x = pos_x, y = pos_y, label = ch_labels,
                        stringsAsFactors = FALSE)

  # Build the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = grid_df,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      fill = .data$value),
                         interpolate = TRUE) +
    ggplot2::scale_fill_distiller(palette = palette, limits = c(-max_abs, max_abs),
                                 direction = -1, na.value = "transparent",
                                 name = "Value")

  # Add contours
  if (contours) {
    grid_clean <- grid_df[!is.na(grid_df$value), ]
    if (nrow(grid_clean) > 0) {
      p <- p + ggplot2::geom_contour(
        data = grid_clean,
        ggplot2::aes(x = .data$x, y = .data$y, z = .data$value),
        color = "grey40", linewidth = 0.2, alpha = 0.5
      )
    }
  }

  # Head outline and anatomical markers
  p <- p +
    ggplot2::geom_path(data = head_df,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 0.8) +
    ggplot2::geom_path(data = nose_df,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 0.8) +
    ggplot2::geom_path(data = ear_left,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 0.6) +
    ggplot2::geom_path(data = ear_right,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 0.6)

  # Electrode markers
  if (electrodes) {
    p <- p +
      ggplot2::geom_point(data = elec_df,
                          ggplot2::aes(x = .data$x, y = .data$y),
                          size = 1.5, shape = 21, fill = "white",
                          color = "black") +
      ggplot2::geom_text(data = elec_df,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      label = .data$label),
                         size = 2, nudge_y = 0.06, check_overlap = TRUE)
  }

  title_text <- if (!is.null(time)) {
    sprintf("Topographic Map (t = %.3f s)", time)
  } else {
    "Topographic Map"
  }

  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title_text) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p
}


#' Plot Topographic Map Series
#'
#' Creates multiple topographic maps at specified time points arranged
#' in a grid layout with a shared or independent color scale.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param times Numeric vector of time points in seconds.
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @param ncol Integer number of columns in the grid layout.
#'   If \code{NULL}, auto-calculated.
#' @param palette Character name of the diverging color palette
#'   (default: \code{"RdBu"}).
#' @param shared_limits Logical; if \code{TRUE}, use the same color scale
#'   across all panels.
#' @return A ggplot2 object with faceted topomaps.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 2500, n_channels = 19, sr = 500)
#' eegPlotTopomapSeries(pe, times = c(0.1, 0.2, 0.3))
#' }
eegPlotTopomapSeries <- function(x, times, assay_name = NULL,
                                 ncol = NULL, palette = "RdBu",
                                 shared_limits = TRUE) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(times) && length(times) >= 1)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Handle 3D data
  if (length(dim(data)) == 3) {
    data <- apply(data, c(1, 2), mean)
  }

  n_channels <- ncol(data)
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  # Electrode positions
  if (all(c("pos_x", "pos_y") %in% colnames(col_data))) {
    pos_x <- as.numeric(col_data$pos_x)
    pos_y <- as.numeric(col_data$pos_y)
  } else {
    defaults <- .default_1020_positions()
    pos_x <- numeric(n_channels)
    pos_y <- numeric(n_channels)
    for (i in seq_len(n_channels)) {
      lbl <- ch_labels[i]
      if (lbl %in% names(defaults$pos_x)) {
        pos_x[i] <- defaults$pos_x[lbl]
        pos_y[i] <- defaults$pos_y[lbl]
      } else {
        angle <- 2 * pi * i / n_channels
        pos_x[i] <- 0.5 * cos(angle)
        pos_y[i] <- 0.5 * sin(angle)
      }
    }
  }

  # Use reduced resolution for multi-panel efficiency
  resolution <- 50

  # Build combined grid for faceting
  all_grids <- list()
  for (t_i in seq_along(times)) {
    t_val <- times[t_i]
    sample_idx <- max(1L, min(nrow(data),
                              as.integer(round(t_val * sr)) + 1L))
    ch_values <- data[sample_idx, ]

    grid_df <- .interpolate_topomap(ch_values, pos_x, pos_y,
                                    resolution = resolution)
    grid_df$time_label <- sprintf("t = %.3f s", t_val)
    all_grids[[t_i]] <- grid_df
  }

  combined_df <- do.call(rbind, all_grids)
  combined_df$time_label <- factor(combined_df$time_label,
                                   levels = sprintf("t = %.3f s", times))

  # Color limits
  if (shared_limits) {
    max_abs <- max(abs(combined_df$value), na.rm = TRUE)
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1
  } else {
    max_abs <- max(abs(combined_df$value), na.rm = TRUE)
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1
  }

  # Head outline for each panel
  head_df <- .head_outline()

  if (is.null(ncol)) {
    ncol <- min(length(times), 4L)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = combined_df,
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      fill = .data$value),
                         interpolate = TRUE) +
    ggplot2::scale_fill_distiller(palette = palette,
                                 limits = c(-max_abs, max_abs),
                                 direction = -1,
                                 na.value = "transparent",
                                 name = "Value") +
    ggplot2::geom_path(data = head_df,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 0.6) +
    ggplot2::facet_wrap(~ time_label, ncol = ncol) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Topographic Map Series") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(size = 9)
    )

  p
}


#' Plot Spectrogram (Time-Frequency Heatmap)
#'
#' Displays a time-frequency representation of EEG data as a heatmap.
#' Expects the assay to contain 3D time-frequency data or a pre-computed
#' spectrogram. For 2D data, a simple Welch-based spectrogram is computed
#' using sliding window FFT.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param channel Integer or character specifying which channel to display
#'   (default: 1).
#' @param freq_range Numeric vector of length 2 for frequency axis limits.
#'   If \code{NULL}, the full range is shown.
#' @param time_range Numeric vector of length 2 for time axis limits.
#'   If \code{NULL}, the full range is shown.
#' @param log_power Logical; if \code{TRUE}, plot \code{10*log10(power)}
#'   (default: TRUE).
#' @param palette Character name of the color palette (default: \code{"viridis"}).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' eegPlotSpectrogram(pe, channel = 1)
#' }
eegPlotSpectrogram <- function(x, channel = 1, freq_range = NULL,
                               time_range = NULL, log_power = TRUE,
                               palette = "viridis", assay_name = NULL) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  md <- S4Vectors::metadata(x)

  # Resolve channel index
  if (is.character(channel)) {
    col_data <- SummarizedExperiment::colData(x)
    ch_labels <- if ("label" %in% colnames(col_data)) {
      as.character(col_data$label)
    } else {
      paste0("Ch", seq_len(ncol(data)))
    }
    ch_idx <- which(ch_labels == channel)
    if (length(ch_idx) == 0) {
      stop(sprintf("Channel '%s' not found.", channel), call. = FALSE)
    }
    ch_idx <- ch_idx[1]
  } else {
    ch_idx <- as.integer(channel)
  }

  # Check if we have pre-computed time-frequency data
  if (length(dim(data)) == 3 && !is.null(md$freqs)) {
    # 3D assay: time x frequency x channels or similar
    tf_matrix <- data[, , ch_idx]
    freqs <- md$freqs
    time_vec <- if (!is.null(md$times)) md$times else {
      (seq_len(nrow(tf_matrix)) - 1) / sr
    }
  } else {
    # Compute spectrogram from 2D data using sliding window FFT
    if (length(dim(data)) == 3) {
      sig <- data[, ch_idx, 1]
    } else {
      sig <- data[, ch_idx]
    }
    n <- length(sig)

    # Window parameters
    window_len <- min(as.integer(round(0.5 * sr)), n)
    if (window_len %% 2 != 0) window_len <- window_len - 1L
    hop <- max(1L, window_len %/% 4)
    win <- 0.5 * (1 - cos(2 * pi * seq(0, window_len - 1) / (window_len - 1)))

    n_fft_bins <- window_len %/% 2 + 1L
    freqs <- seq(0, sr / 2, length.out = n_fft_bins)

    starts <- seq(1L, n - window_len + 1L, by = hop)
    n_windows <- length(starts)
    time_vec <- (starts + window_len / 2 - 1) / sr

    tf_matrix <- matrix(0, nrow = n_windows, ncol = n_fft_bins)
    for (wi in seq_len(n_windows)) {
      idx <- starts[wi]:(starts[wi] + window_len - 1L)
      seg <- sig[idx] * win
      ft <- fft(seg)[seq_len(n_fft_bins)]
      tf_matrix[wi, ] <- Mod(ft)^2 / window_len
    }
  }

  if (log_power) {
    tf_matrix <- 10 * log10(pmax(tf_matrix, .Machine$double.eps))
  }

  # Build long-format data frame
  plot_df <- expand.grid(
    time = time_vec,
    freq = freqs
  )
  plot_df$power <- as.vector(tf_matrix)

  # Apply frequency range filter
  if (!is.null(freq_range)) {
    plot_df <- plot_df[plot_df$freq >= freq_range[1] &
                       plot_df$freq <= freq_range[2], ]
  }

  # Apply time range filter
  if (!is.null(time_range)) {
    plot_df <- plot_df[plot_df$time >= time_range[1] &
                       plot_df$time <= time_range[2], ]
  }

  power_label <- if (log_power) "Power (dB)" else "Power"

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$time, y = .data$freq,
                                    fill = .data$power)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(option = sub("viridis", "D", palette),
                                  name = power_label) +
    ggplot2::labs(x = "Time (s)", y = "Frequency (Hz)",
                  title = "Spectrogram") +
    ggplot2::theme_minimal()

  p
}


#' Plot Connectivity Matrix or Circle
#'
#' Visualizes EEG connectivity as either a heatmap (matrix) or a circular
#' connectivity plot (circle). The connectivity matrix can be provided
#' directly or read from \code{metadata(x)$connectivity$matrix}.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Display method: \code{"heatmap"} for a correlation/connectivity
#'   matrix heatmap or \code{"circle"} for a circular connectivity diagram.
#' @param matrix Numeric matrix of connectivity values. If \code{NULL},
#'   reads from \code{metadata(x)$connectivity$matrix}.
#' @param threshold Numeric; only show connections above this value
#'   (default: 0).
#' @param labels Character vector of channel labels. If \code{NULL}, uses
#'   row/column names of the matrix or channel labels from colData.
#' @param palette Character name of the diverging color palette
#'   (default: \code{"RdBu"}).
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' mat <- matrix(runif(16), 4, 4)
#' diag(mat) <- 1
#' pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
#' eegPlotConnectivity(pe, method = "heatmap", matrix = mat)
#' }
eegPlotConnectivity <- function(x, method = c("heatmap", "circle"),
                                matrix = NULL, threshold = 0,
                                labels = NULL, palette = "RdBu") {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  # Get connectivity matrix
  if (is.null(matrix)) {
    md <- S4Vectors::metadata(x)
    if (!is.null(md$connectivity) && !is.null(md$connectivity$matrix)) {
      matrix <- md$connectivity$matrix
    } else {
      stop("No connectivity matrix provided and none found in metadata. ",
           "Provide a matrix or run a connectivity analysis first.",
           call. = FALSE)
    }
  }

  n_ch <- nrow(matrix)

  # Labels
  if (is.null(labels)) {
    if (!is.null(rownames(matrix))) {
      labels <- rownames(matrix)
    } else {
      col_data <- SummarizedExperiment::colData(x)
      if ("label" %in% colnames(col_data)) {
        labels <- as.character(col_data$label)[seq_len(n_ch)]
      } else {
        labels <- paste0("Ch", seq_len(n_ch))
      }
    }
  }

  if (method == "heatmap") {
    # Build long-format data for heatmap
    plot_df <- expand.grid(
      row = seq_len(n_ch),
      col = seq_len(n_ch)
    )
    plot_df$value <- as.vector(matrix)
    plot_df$row_label <- factor(labels[plot_df$row], levels = rev(labels))
    plot_df$col_label <- factor(labels[plot_df$col], levels = labels)

    # Apply threshold (set below-threshold to NA)
    plot_df$value[abs(plot_df$value) < threshold] <- NA_real_

    max_abs <- max(abs(plot_df$value), na.rm = TRUE)
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$col_label,
                                      y = .data$row_label,
                                      fill = .data$value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_distiller(palette = palette,
                                    limits = c(-max_abs, max_abs),
                                    direction = -1,
                                    na.value = "grey90",
                                    name = "Connectivity") +
      ggplot2::labs(x = "", y = "", title = "Connectivity Matrix") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::coord_fixed()

  } else {
    # Circle method: channels arranged in a circle with lines for connections
    angles <- seq(0, 2 * pi, length.out = n_ch + 1)[seq_len(n_ch)]
    node_x <- cos(angles)
    node_y <- sin(angles)

    # Node data
    node_df <- data.frame(
      x = node_x, y = node_y,
      label = labels,
      stringsAsFactors = FALSE
    )

    # Edge data (upper triangle, above threshold)
    edge_list <- list()
    for (i in seq_len(n_ch - 1)) {
      for (j in (i + 1):n_ch) {
        val <- matrix[i, j]
        if (abs(val) > threshold) {
          edge_list[[length(edge_list) + 1]] <- data.frame(
            x = node_x[i], y = node_y[i],
            xend = node_x[j], yend = node_y[j],
            value = val,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    p <- ggplot2::ggplot()

    if (length(edge_list) > 0) {
      edge_df <- do.call(rbind, edge_list)
      max_abs <- max(abs(edge_df$value), na.rm = TRUE)
      if (max_abs == 0) max_abs <- 1

      p <- p +
        ggplot2::geom_segment(
          data = edge_df,
          ggplot2::aes(x = .data$x, y = .data$y,
                       xend = .data$xend, yend = .data$yend,
                       color = .data$value,
                       linewidth = abs(.data$value)),
          alpha = 0.6
        ) +
        ggplot2::scale_color_distiller(palette = palette,
                                       limits = c(-max_abs, max_abs),
                                       direction = -1,
                                       name = "Connectivity") +
        ggplot2::scale_linewidth_continuous(range = c(0.3, 2),
                                           guide = "none")
    }

    p <- p +
      ggplot2::geom_point(data = node_df,
                          ggplot2::aes(x = .data$x, y = .data$y),
                          size = 4, color = "black", fill = "white",
                          shape = 21) +
      ggplot2::geom_text(data = node_df,
                         ggplot2::aes(x = .data$x * 1.15,
                                      y = .data$y * 1.15,
                                      label = .data$label),
                         size = 3) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = "Circular Connectivity Plot") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p
}


#' Plot Sleep Hypnogram
#'
#' Displays a sleep hypnogram showing sleep stage transitions over time.
#' Stages are ordered with Wake at the top and N3 at the bottom, with
#' a characteristic staircase pattern.
#'
#' @param x A PhysioExperiment object (used for metadata access).
#' @param stages A data.frame with columns \code{epoch} and \code{stage},
#'   or \code{NULL} to read from \code{metadata(x)$sleep_stages}.
#' @param epoch_sec Epoch duration in seconds (default: 30).
#' @param colors Named character vector of colors for each stage.
#'   If \code{NULL}, default colors are used.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
#' stages <- eegSleepStage(pe)
#' metadata(pe)$sleep_stages <- stages
#' eegPlotHypnogram(pe)
#' }
eegPlotHypnogram <- function(x, stages = NULL, epoch_sec = 30,
                             colors = NULL) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))

  # Get stage data
  if (is.null(stages)) {
    stages <- S4Vectors::metadata(x)$sleep_stages
    if (is.null(stages)) {
      stop("No sleep stage data provided and none found in metadata. ",
           "Provide a data.frame or run eegSleepStage() first.",
           call. = FALSE)
    }
  }

  stopifnot(is.data.frame(stages))
  stopifnot("epoch" %in% names(stages) && "stage" %in% names(stages))

  # Stage ordering: W (top = 5), REM (4), N1 (3), N2 (2), N3 (bottom = 1)
  stage_levels <- c("N3", "N2", "N1", "REM", "W")
  stage_numeric <- c(N3 = 1, N2 = 2, N1 = 3, REM = 4, W = 5)

  # Default colors
  if (is.null(colors)) {
    colors <- c(W = "#FFD700", REM = "#FF4444", N1 = "#87CEEB",
                N2 = "#4169E1", N3 = "#00008B")
  }

  # Time in minutes for each epoch
  stages$time_min <- (stages$epoch - 1) * epoch_sec / 60
  stages$stage_y <- stage_numeric[stages$stage]

  # Ensure we have valid stage mappings
  valid_mask <- !is.na(stages$stage_y)
  stages <- stages[valid_mask, , drop = FALSE]

  if (nrow(stages) == 0) {
    stop("No valid sleep stages found.", call. = FALSE)
  }

  stages$stage_factor <- factor(stages$stage, levels = stage_levels)

  p <- ggplot2::ggplot(stages,
                       ggplot2::aes(x = .data$time_min,
                                    y = .data$stage_y,
                                    color = .data$stage)) +
    ggplot2::geom_step(linewidth = 1.0, direction = "hv") +
    ggplot2::scale_y_continuous(
      breaks = 1:5,
      labels = stage_levels,
      limits = c(0.5, 5.5)
    ) +
    ggplot2::scale_color_manual(values = colors, name = "Stage") +
    ggplot2::labs(x = "Time (min)", y = "Sleep Stage",
                  title = "Hypnogram") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


#' Plot ICA Components
#'
#' Displays ICA component time courses in a stacked layout similar to
#' \code{eegPlotSignal}. Optionally shows topographic maps of the mixing
#' matrix weights as a side panel annotation.
#'
#' @param x A PhysioExperiment object with ICA results.
#' @param components Integer vector of component indices to display.
#'   If \code{NULL}, the first 10 (or fewer) are shown.
#' @param time_range Numeric vector of length 2 specifying time range in
#'   seconds. If \code{NULL}, all data is shown.
#' @param show_topography Logical; if \code{TRUE} and ICA mixing matrix is
#'   available in \code{metadata(x)$ica$mixing}, include topographic
#'   inset labels.
#' @param assay_name Input assay name. If \code{NULL}, uses \code{"ica"}
#'   if available, otherwise the default assay.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe <- eegICA(pe, n_components = 10)
#' eegPlotICA(pe, components = 1:5)
#' }
eegPlotICA <- function(x, components = NULL, time_range = NULL,
                       show_topography = TRUE, assay_name = NULL) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))

  # Determine assay
  if (is.null(assay_name)) {
    avail <- SummarizedExperiment::assayNames(x)
    if ("ica" %in% avail) {
      assay_name <- "ica"
    } else {
      assay_name <- defaultAssay(x)
    }
  }
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Handle 2D data (time x components)
  if (length(dim(data)) == 3) {
    data <- data[, , 1, drop = TRUE]
  }

  n_time <- nrow(data)
  n_comps <- ncol(data)

  # Select components
  if (is.null(components)) {
    components <- seq_len(min(10L, n_comps))
  }
  components <- components[components >= 1 & components <= n_comps]
  if (length(components) == 0) {
    stop("No valid components selected.", call. = FALSE)
  }

  # Time vector
  time_sec <- (seq_len(n_time) - 1) / sr

  # Time range filter
  if (!is.null(time_range)) {
    stopifnot(is.numeric(time_range) && length(time_range) == 2)
    time_idx <- which(time_sec >= time_range[1] & time_sec <= time_range[2])
    if (length(time_idx) == 0) {
      stop("No data points within the specified time_range.", call. = FALSE)
    }
  } else {
    time_idx <- seq_len(n_time)
  }

  # Component labels with topography annotation
  comp_labels <- paste0("IC", components)
  md <- S4Vectors::metadata(x)
  has_topo <- show_topography && !is.null(md$ica) && !is.null(md$ica$mixing)

  if (has_topo) {
    mixing <- md$ica$mixing  # channels x components
    # Append max weight info to labels
    for (i in seq_along(components)) {
      comp_idx <- components[i]
      if (comp_idx <= ncol(mixing)) {
        weights <- mixing[, comp_idx]
        max_ch_idx <- which.max(abs(weights))
        col_data <- SummarizedExperiment::colData(x)
        ch_labels_all <- if ("label" %in% colnames(col_data)) {
          as.character(col_data$label)
        } else {
          paste0("Ch", seq_len(nrow(col_data)))
        }
        if (max_ch_idx <= length(ch_labels_all)) {
          comp_labels[i] <- paste0("IC", comp_idx, " (",
                                   ch_labels_all[max_ch_idx], ")")
        }
      }
    }
  }

  # Build stacked plot data
  plot_list <- vector("list", length(components))
  for (i in seq_along(components)) {
    comp <- components[i]
    vals <- data[time_idx, comp]
    comp_sd <- stats::sd(data[, comp], na.rm = TRUE)
    if (comp_sd == 0) comp_sd <- 1
    offset <- (length(components) - i) * comp_sd * 3
    vals <- vals + offset

    plot_list[[i]] <- data.frame(
      time = time_sec[time_idx],
      amplitude = vals,
      component = comp_labels[i],
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, plot_list)
  plot_df$component <- factor(plot_df$component, levels = comp_labels)

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$time, y = .data$amplitude,
                                    color = .data$component)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(x = "Time (s)", y = "Amplitude", color = "Component",
                  title = "ICA Components") +
    ggplot2::theme_minimal()

  p
}


#' Plot Source Localization Results
#'
#' Visualizes source localization results as a 2D scatter plot or flat map
#' projection. Sources are sized and colored by amplitude, with optional
#' thresholding to show only the strongest activations.
#'
#' @param x A PhysioExperiment object.
#' @param source_data Named numeric vector of source amplitudes or
#'   a data.frame with columns \code{x}, \code{y}, and \code{amplitude}.
#'   If \code{NULL}, reads from \code{metadata(x)$source_estimate}.
#' @param method Display method: \code{"scatter"} for points on a 2D brain
#'   outline or \code{"flatmap"} for filled regions using interpolation.
#' @param threshold_pct Numeric percentile threshold (0-100). Only sources
#'   above this percentile are displayed (default: 80).
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
#' src <- data.frame(x = runif(50, -1, 1), y = runif(50, -1, 1),
#'                   amplitude = rnorm(50)^2)
#' eegPlotSource(pe, source_data = src, method = "scatter")
#' }
eegPlotSource <- function(x, source_data = NULL,
                          method = c("scatter", "flatmap"),
                          threshold_pct = 80) {
  .check_ggplot2()
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(threshold_pct) && threshold_pct >= 0 &&
            threshold_pct <= 100)

  # Get source data
  if (is.null(source_data)) {
    md <- S4Vectors::metadata(x)
    if (!is.null(md$source_estimate)) {
      source_data <- md$source_estimate
    } else {
      stop("No source data provided and none found in metadata. ",
           "Provide source_data or run source localization first.",
           call. = FALSE)
    }
  }

  # Normalize source_data to a data.frame with x, y, amplitude
  if (is.numeric(source_data) && !is.data.frame(source_data)) {
    # Named vector: distribute sources on a grid inside a circle
    n_src <- length(source_data)
    # Use Fibonacci sphere-like distribution projected to 2D
    golden_angle <- pi * (3 - sqrt(5))
    src_x <- numeric(n_src)
    src_y <- numeric(n_src)
    for (i in seq_len(n_src)) {
      r <- sqrt((i - 0.5) / n_src) * 0.9
      theta <- golden_angle * i
      src_x[i] <- r * cos(theta)
      src_y[i] <- r * sin(theta)
    }
    source_df <- data.frame(
      x = src_x, y = src_y,
      amplitude = as.numeric(source_data),
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(source_data)) {
    source_df <- source_data
    if (!"amplitude" %in% names(source_df)) {
      stop("source_data data.frame must contain an 'amplitude' column.",
           call. = FALSE)
    }
    if (!all(c("x", "y") %in% names(source_df))) {
      stop("source_data data.frame must contain 'x' and 'y' columns.",
           call. = FALSE)
    }
  } else {
    stop("source_data must be a numeric vector or data.frame.", call. = FALSE)
  }

  # Apply threshold
  thresh_val <- stats::quantile(abs(source_df$amplitude), threshold_pct / 100)
  source_df_thresh <- source_df[abs(source_df$amplitude) >= thresh_val, ,
                                drop = FALSE]

  # Brain outline (ellipse approximation)
  brain_t <- seq(0, 2 * pi, length.out = 100)
  brain_df <- data.frame(
    x = 0.95 * cos(brain_t),
    y = 0.95 * sin(brain_t)
  )

  if (method == "scatter") {
    if (nrow(source_df_thresh) == 0) {
      # No sources above threshold; show empty brain
      p <- ggplot2::ggplot() +
        ggplot2::geom_path(data = brain_df,
                           ggplot2::aes(x = .data$x, y = .data$y),
                           linewidth = 0.8) +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Source Localization (Scatter)") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      max_amp <- max(abs(source_df_thresh$amplitude))
      if (max_amp == 0) max_amp <- 1

      p <- ggplot2::ggplot() +
        ggplot2::geom_path(data = brain_df,
                           ggplot2::aes(x = .data$x, y = .data$y),
                           linewidth = 0.8) +
        ggplot2::geom_point(
          data = source_df_thresh,
          ggplot2::aes(x = .data$x, y = .data$y,
                       color = .data$amplitude,
                       size = abs(.data$amplitude)),
          alpha = 0.7
        ) +
        ggplot2::scale_color_distiller(palette = "Spectral",
                                       name = "Amplitude") +
        ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Source Localization (Scatter)") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }

  } else {
    # flatmap: interpolation-based filled regions
    if (nrow(source_df_thresh) < 3) {
      # Fall back to scatter if too few points for interpolation
      return(eegPlotSource(x, source_data = source_data,
                           method = "scatter",
                           threshold_pct = threshold_pct))
    }

    grid_df <- .interpolate_topomap(
      source_df_thresh$amplitude,
      source_df_thresh$x,
      source_df_thresh$y,
      resolution = 80
    )

    max_abs <- max(abs(grid_df$value), na.rm = TRUE)
    if (is.na(max_abs) || max_abs == 0) max_abs <- 1

    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = grid_df,
                           ggplot2::aes(x = .data$x, y = .data$y,
                                        fill = .data$value),
                           interpolate = TRUE) +
      ggplot2::scale_fill_distiller(palette = "Spectral",
                                    limits = c(-max_abs, max_abs),
                                    na.value = "transparent",
                                    name = "Amplitude") +
      ggplot2::geom_path(data = brain_df,
                         ggplot2::aes(x = .data$x, y = .data$y),
                         linewidth = 0.8) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = "Source Localization (Flat Map)") +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p
}
