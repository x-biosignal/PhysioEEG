#' Detect Epileptic Spikes in EEG Data
#'
#' Identifies epileptiform spike discharges in multi-channel EEG using either
#' morphology-based derivative analysis or template matching via
#' cross-correlation. The morphology method detects sharp transients by
#' thresholding the first derivative, while the template method
#' cross-correlates a canonical spike waveform with each channel.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Detection method: \code{"morphology"} (derivative-based) or
#'   \code{"template"} (cross-correlation with canonical spike shape).
#' @param threshold_sd Number of standard deviations above the mean derivative
#'   magnitude for detection (default: 4).
#' @param min_duration_ms Minimum spike duration in milliseconds (default: 20).
#' @param max_duration_ms Maximum spike duration in milliseconds (default: 200).
#' @param min_amplitude Minimum peak amplitude in microvolts for a valid spike
#'   detection (default: 50).
#' @param assay_name Name of the assay to use. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{sample}{Integer sample index of the spike peak.}
#'     \item{time_sec}{Time of the spike in seconds.}
#'     \item{amplitude}{Peak amplitude at the spike location.}
#'     \item{duration_ms}{Estimated spike duration in milliseconds.}
#'     \item{confidence}{Confidence score for the detection.}
#'   }
#' @references
#' Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
#' clinical EEG. Electroencephalography and Clinical Neurophysiology, 106(3), 259-261.
#' @seealso [eegQEEG()], [eegAsymmetry()], [eegSlowing()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 15)
#' spikes <- eegSpikeDetect(pe, method = "morphology")
#' head(spikes)
#' }
eegSpikeDetect <- function(x, method = c("morphology", "template"),
                           threshold_sd = 4, min_duration_ms = 20,
                           max_duration_ms = 200, min_amplitude = 50,
                           assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(threshold_sd) && threshold_sd > 0)
  stopifnot(is.numeric(min_duration_ms) && min_duration_ms > 0)
  stopifnot(is.numeric(max_duration_ms) && max_duration_ms > min_duration_ms)
  stopifnot(is.numeric(min_amplitude) && min_amplitude > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  min_dur_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))
  max_dur_samples <- as.integer(round(max_duration_ms / 1000 * sr))

  results <- list()

  if (method == "morphology") {
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]

      # Compute first derivative
      dsig <- c(diff(sig), 0)
      abs_dsig <- abs(dsig)
      mean_deriv <- mean(abs_dsig)
      sd_deriv <- sd(abs_dsig)
      deriv_threshold <- mean_deriv + threshold_sd * sd_deriv

      # Find points where |derivative| exceeds threshold
      above <- which(abs_dsig > deriv_threshold)
      if (length(above) == 0) next

      # Group nearby detections (within max_duration_ms)
      groups <- list()
      current_group <- above[1]
      for (i in seq_along(above)[-1]) {
        if (above[i] - above[i - 1] <= max_dur_samples) {
          current_group <- c(current_group, above[i])
        } else {
          groups[[length(groups) + 1]] <- current_group
          current_group <- above[i]
        }
      }
      groups[[length(groups) + 1]] <- current_group

      # Evaluate each candidate group
      sig_mean <- mean(abs(sig))
      sig_sd <- sd(sig)

      for (g in groups) {
        group_start <- min(g)
        group_end <- max(g)
        duration_samples <- group_end - group_start + 1L

        # Check duration within range
        if (duration_samples < min_dur_samples || duration_samples > max_dur_samples) next

        # Find peak amplitude within group region (extend slightly)
        search_start <- max(1L, group_start - min_dur_samples)
        search_end <- min(n_time, group_end + min_dur_samples)
        seg <- sig[search_start:search_end]
        peak_offset <- which.max(abs(seg))
        peak_sample <- search_start + peak_offset - 1L
        peak_amp <- abs(sig[peak_sample])

        # Check peak amplitude exceeds minimum
        if (peak_amp < min_amplitude) next

        duration_ms <- duration_samples / sr * 1000
        confidence <- peak_amp / (sig_mean + 2 * sig_sd)
        confidence <- min(confidence, 1.0)

        results[[length(results) + 1]] <- data.frame(
          channel = as.integer(ch),
          sample = peak_sample,
          time_sec = (peak_sample - 1) / sr,
          amplitude = peak_amp,
          duration_ms = duration_ms,
          confidence = confidence,
          stringsAsFactors = FALSE
        )
      }
    }

  } else if (method == "template") {
    # Create canonical spike template: fast rise (Gaussian), slow fall (exponential)
    template_samples <- max(2L, as.integer(round(max_duration_ms / 1000 * sr)))
    t_template <- seq(0, (template_samples - 1) / sr, length.out = template_samples)
    sigma_rise <- 0.002  # 2ms rise sigma
    tau_fall <- 0.020    # 20ms fall tau

    # Peak at ~20% of template duration
    peak_time <- t_template[max(1L, as.integer(template_samples * 0.2))]
    template <- numeric(template_samples)
    for (i in seq_len(template_samples)) {
      if (t_template[i] <= peak_time) {
        template[i] <- exp(-((t_template[i] - peak_time)^2) / (2 * sigma_rise^2))
      } else {
        template[i] <- exp(-(t_template[i] - peak_time) / tau_fall)
      }
    }
    # Normalize template
    template <- template / sqrt(sum(template^2))

    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]

      # Cross-correlate template with signal
      n_corr <- n_time - template_samples + 1L
      if (n_corr < 1) next

      corr_vals <- numeric(n_corr)
      for (i in seq_len(n_corr)) {
        seg <- sig[i:(i + template_samples - 1L)]
        seg_norm <- sqrt(sum(seg^2))
        if (seg_norm > 0) {
          corr_vals[i] <- sum(seg * template) / seg_norm
        }
      }

      # Find peaks in correlation > 0.7
      corr_threshold <- 0.7
      above <- which(corr_vals > corr_threshold)
      if (length(above) == 0) next

      # Group nearby detections and take the best in each group
      groups <- list()
      current_group <- above[1]
      for (i in seq_along(above)[-1]) {
        if (above[i] - above[i - 1] <= max_dur_samples) {
          current_group <- c(current_group, above[i])
        } else {
          groups[[length(groups) + 1]] <- current_group
          current_group <- above[i]
        }
      }
      groups[[length(groups) + 1]] <- current_group

      for (g in groups) {
        # Take the position with the highest correlation
        best_idx <- g[which.max(corr_vals[g])]

        # Find peak amplitude at detected location
        seg <- sig[best_idx:(best_idx + template_samples - 1L)]
        peak_offset <- which.max(abs(seg))
        peak_sample <- best_idx + peak_offset - 1L
        peak_amp <- abs(sig[peak_sample])

        # Check amplitude
        if (peak_amp < min_amplitude) next

        duration_ms <- template_samples / sr * 1000
        confidence <- corr_vals[best_idx]

        results[[length(results) + 1]] <- data.frame(
          channel = as.integer(ch),
          sample = peak_sample,
          time_sec = (peak_sample - 1) / sr,
          amplitude = peak_amp,
          duration_ms = duration_ms,
          confidence = confidence,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(channel = integer(0), sample = integer(0),
               time_sec = numeric(0), amplitude = numeric(0),
               duration_ms = numeric(0), confidence = numeric(0))
  }
}


#' Quantitative EEG (QEEG) Analysis
#'
#' Computes absolute and relative spectral band powers for each channel using
#' Welch's method (windowed FFT averaging). Results are stored as a new assay
#' (channels x bands matrix of absolute power) and as relative power in
#' \code{metadata(x)$qeeg}.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param bands Named list of frequency bands. Each element is a numeric
#'   vector of length 2 specifying the lower and upper frequency in Hz.
#'   Defaults to standard EEG bands: delta (1-4), theta (4-8), alpha (8-13),
#'   beta (13-30), gamma (30-50).
#' @param window_sec Window length in seconds for Welch's method (default: 2).
#' @param overlap Overlap fraction between windows, 0 to 1 (default: 0.5).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store absolute power results
#'   (default: \code{"qeeg"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item Absolute power matrix (n_channels x n_bands) in \code{output_assay}
#'     \item Band definitions and relative power in \code{metadata(x)$qeeg},
#'       a list containing \code{bands}, \code{absolute_power},
#'       \code{relative_power}, \code{band_names}, \code{window_sec},
#'       and \code{overlap}.
#'   }
#' @references
#' Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
#' clinical EEG. Electroencephalography and Clinical Neurophysiology, 106(3), 259-261.
#'
#' Thatcher, R. W. (2010). Validity and reliability of quantitative
#' electroencephalography. Journal of Neurotherapy, 14(2), 122-152.
#' @seealso [eegSpikeDetect()], [eegAsymmetry()], [eegSlowing()],
#'   [eegPlotSpectrogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe_qeeg <- eegQEEG(pe)
#' qeeg_info <- metadata(pe_qeeg)$qeeg
#' print(qeeg_info$relative_power)
#' }
eegQEEG <- function(x, bands = NULL, window_sec = 2, overlap = 0.5,
                    assay_name = NULL, output_assay = "qeeg") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(window_sec) && window_sec > 0)
  stopifnot(is.numeric(overlap) && overlap >= 0 && overlap < 1)

  if (is.null(bands)) {
    bands <- list(
      delta = c(1, 4),
      theta = c(4, 8),
      alpha = c(8, 13),
      beta  = c(13, 30),
      gamma = c(30, 50)
    )
  }
  stopifnot(is.list(bands) && length(bands) > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)
  n_bands <- length(bands)
  band_names <- names(bands)

  window_samples <- as.integer(round(window_sec * sr))
  step_samples <- max(1L, as.integer(round(window_samples * (1 - overlap))))

  # Compute Hanning window
  hanning <- 0.5 * (1 - cos(2 * pi * seq(0, window_samples - 1) / (window_samples - 1)))

  # Window energy correction factor for proper PSD normalisation
  # PSD = |FFT(x * w)|^2 / (sr * sum(w^2)) gives units of signal^2/Hz
  window_energy <- sr * sum(hanning^2)

  # Frequency vector for windowed FFT
  freqs <- (0:(window_samples - 1)) * sr / window_samples
  half_n <- window_samples %/% 2

  # Absolute power matrix: n_channels x n_bands
  abs_power <- matrix(0, nrow = n_channels, ncol = n_bands)
  colnames(abs_power) <- band_names

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # Determine window start positions
    starts <- seq(1L, n_time - window_samples + 1L, by = step_samples)
    n_windows <- length(starts)

    if (n_windows < 1) {
      # Signal too short for even one window; use entire signal
      starts <- 1L
      n_windows <- 1L
    }

    # Accumulate averaged PSD
    psd_accum <- numeric(window_samples)

    for (w in seq_len(n_windows)) {
      seg <- sig[starts[w]:(starts[w] + window_samples - 1L)]
      # Apply Hanning window
      seg_win <- seg * hanning
      # Compute |FFT|^2 with proper window energy correction
      ft <- fft(seg_win)
      psd_accum <- psd_accum + (Mod(ft)^2) / window_energy
    }

    # Average across windows
    psd_avg <- psd_accum / n_windows

    # Compute absolute power per band (sum PSD in band frequency bins)
    for (b in seq_len(n_bands)) {
      low <- bands[[b]][1]
      high <- bands[[b]][2]
      band_idx <- which(freqs[1:half_n] >= low & freqs[1:half_n] <= high)
      abs_power[ch, b] <- sum(psd_avg[band_idx])
    }
  }

  # Compute relative power: per channel, normalize to sum = 1
  total_power <- rowSums(abs_power)
  total_power[total_power < .Machine$double.eps] <- 1
  rel_power <- abs_power / total_power

  # Store band info, absolute and relative power in metadata
  # (QEEG results have different dimensions than raw assays)
  qeeg_info <- list(
    bands = bands,
    absolute_power = abs_power,
    relative_power = rel_power,
    band_names = band_names,
    window_sec = window_sec,
    overlap = overlap
  )
  md <- S4Vectors::metadata(x)
  md$qeeg <- qeeg_info
  S4Vectors::metadata(x) <- md

  x
}


#' Frontal Alpha Asymmetry
#'
#' Computes frontal alpha asymmetry indices from paired electrode sites.
#' For each pair (right, left), band power is computed via FFT and the
#' asymmetry index is calculated as \code{log(power_right) - log(power_left)}.
#' Positive values indicate greater right-hemisphere alpha (typically
#' associated with greater left-hemisphere cortical activity).
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param pairs A list of character vectors of length 2, each specifying
#'   a (right, left) electrode pair. Defaults to \code{list(c("F4", "F3"),
#'   c("F8", "F7"))}.
#' @param band Numeric vector of length 2 specifying the frequency band
#'   in Hz (default: \code{c(8, 13)} for alpha).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{pair}{Character label for the electrode pair.}
#'     \item{left_channel}{Character name of the left electrode.}
#'     \item{right_channel}{Character name of the right electrode.}
#'     \item{left_power}{Numeric band power for the left electrode.}
#'     \item{right_power}{Numeric band power for the right electrode.}
#'     \item{asymmetry_index}{Numeric asymmetry: log(right) - log(left).}
#'   }
#' @references
#' Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
#' clinical EEG. Electroencephalography and Clinical Neurophysiology, 106(3), 259-261.
#' @seealso [eegSpikeDetect()], [eegQEEG()], [eegSlowing()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' asym <- eegAsymmetry(pe)
#' print(asym)
#' }
eegAsymmetry <- function(x, pairs = NULL, band = c(8, 13),
                         assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(band) && length(band) == 2)
  stopifnot(band[1] < band[2])

  if (is.null(pairs)) {
    pairs <- list(c("F4", "F3"), c("F8", "F7"))
  }
  stopifnot(is.list(pairs))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)

  # Get channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(ncol(data)))
  }

  results <- list()

  for (p_idx in seq_along(pairs)) {
    pair <- pairs[[p_idx]]
    stopifnot(length(pair) == 2)
    right_ch <- pair[1]
    left_ch <- pair[2]

    right_idx <- which(ch_labels == right_ch)
    left_idx <- which(ch_labels == left_ch)

    if (length(right_idx) == 0 || length(left_idx) == 0) {
      warning(sprintf("Channel pair (%s, %s) not found, skipping.",
                      right_ch, left_ch))
      next
    }

    right_idx <- right_idx[1]
    left_idx <- left_idx[1]

    # Compute band power via FFT for each channel
    right_power <- .compute_band_power(data[, right_idx], sr, band[1], band[2])
    left_power <- .compute_band_power(data[, left_idx], sr, band[1], band[2])

    # Guard against zero or negative power for log computation
    right_power_safe <- max(right_power, .Machine$double.eps)
    left_power_safe <- max(left_power, .Machine$double.eps)

    asymmetry_index <- log(right_power_safe) - log(left_power_safe)

    results[[length(results) + 1]] <- data.frame(
      pair = paste0(right_ch, "-", left_ch),
      left_channel = left_ch,
      right_channel = right_ch,
      left_power = left_power,
      right_power = right_power,
      asymmetry_index = asymmetry_index,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(pair = character(0), left_channel = character(0),
               right_channel = character(0), left_power = numeric(0),
               right_power = numeric(0), asymmetry_index = numeric(0))
  }
}


#' Burst-Suppression Detection
#'
#' Identifies periods of burst and suppression in EEG data by computing
#' the root mean square (RMS) amplitude in sliding windows. Suppression
#' segments are defined as consecutive windows where RMS falls below
#' the specified threshold for at least \code{min_duration_ms}.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param threshold RMS amplitude threshold below which signal is considered
#'   suppressed (default: 10).
#' @param min_duration_ms Minimum duration in milliseconds for a segment
#'   to be classified as suppression (default: 500).
#' @param window_ms Window length in milliseconds for RMS computation
#'   (default: 500).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{type}{Character: \code{"burst"} or \code{"suppression"}.}
#'     \item{start_sample}{Integer start sample of the segment.}
#'     \item{end_sample}{Integer end sample of the segment.}
#'     \item{duration_ms}{Numeric duration of the segment in milliseconds.}
#'   }
#'   An attribute \code{"bsr"} (Burst-Suppression Ratio) is attached,
#'   representing the percentage of total time in suppression.
#' @references
#' Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
#' clinical EEG. Electroencephalography and Clinical Neurophysiology, 106(3), 259-261.
#' @seealso [eegSpikeDetect()], [eegQEEG()], [eegSlowing()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
#' result <- eegSuppression(pe, threshold = 10)
#' print(attr(result, "bsr"))
#' }
eegSuppression <- function(x, threshold = 10, min_duration_ms = 500,
                           window_ms = 500, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(threshold) && threshold > 0)
  stopifnot(is.numeric(min_duration_ms) && min_duration_ms > 0)
  stopifnot(is.numeric(window_ms) && window_ms > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  window_samples <- max(1L, as.integer(round(window_ms / 1000 * sr)))
  min_dur_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))

  # Compute average RMS across channels in sliding windows
  n_windows <- max(1L, n_time %/% window_samples)
  rms_vals <- numeric(n_windows)

  for (w in seq_len(n_windows)) {
    start_s <- (w - 1L) * window_samples + 1L
    end_s <- min(w * window_samples, n_time)
    seg <- data[start_s:end_s, , drop = FALSE]
    rms_vals[w] <- sqrt(mean(seg^2))
  }

  # Classify each window
  is_suppressed <- rms_vals < threshold

  # Find contiguous segments of burst and suppression
  results <- list()
  if (length(is_suppressed) > 0) {
    current_type <- if (is_suppressed[1]) "suppression" else "burst"
    seg_start_win <- 1L

    for (w in seq(2, length(is_suppressed))) {
      new_type <- if (is_suppressed[w]) "suppression" else "burst"
      if (new_type != current_type) {
        # Close previous segment
        start_sample <- (seg_start_win - 1L) * window_samples + 1L
        end_sample <- min((w - 1L) * window_samples, n_time)
        duration_samples <- end_sample - start_sample + 1L
        duration_ms <- duration_samples / sr * 1000

        results[[length(results) + 1]] <- data.frame(
          type = current_type,
          start_sample = start_sample,
          end_sample = end_sample,
          duration_ms = duration_ms,
          stringsAsFactors = FALSE
        )
        current_type <- new_type
        seg_start_win <- w
      }
    }

    # Close final segment
    start_sample <- (seg_start_win - 1L) * window_samples + 1L
    end_sample <- min(n_windows * window_samples, n_time)
    duration_samples <- end_sample - start_sample + 1L
    duration_ms <- duration_samples / sr * 1000

    results[[length(results) + 1]] <- data.frame(
      type = current_type,
      start_sample = start_sample,
      end_sample = end_sample,
      duration_ms = duration_ms,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) > 0) {
    result_df <- do.call(rbind, results)
  } else {
    result_df <- data.frame(type = character(0), start_sample = integer(0),
                            end_sample = integer(0), duration_ms = numeric(0))
  }

  # Filter out suppression segments shorter than min_duration_ms
  # But keep all burst segments
  keep <- rep(TRUE, nrow(result_df))
  for (i in seq_len(nrow(result_df))) {
    if (result_df$type[i] == "suppression") {
      duration_samples <- result_df$end_sample[i] - result_df$start_sample[i] + 1L
      if (duration_samples < min_dur_samples) {
        result_df$type[i] <- "burst"  # Reclassify short suppression as burst
      }
    }
  }

  # Merge adjacent segments of the same type after reclassification
  if (nrow(result_df) > 1) {
    merged <- list(result_df[1, ])
    for (i in seq(2, nrow(result_df))) {
      prev <- merged[[length(merged)]]
      curr <- result_df[i, ]
      if (prev$type == curr$type) {
        # Merge
        merged[[length(merged)]]$end_sample <- curr$end_sample
        merged[[length(merged)]]$duration_ms <-
          (curr$end_sample - prev$start_sample + 1L) / sr * 1000
      } else {
        merged[[length(merged) + 1]] <- curr
      }
    }
    result_df <- do.call(rbind, merged)
  }

  # Compute BSR: total suppression time / total time * 100
  total_time <- n_time / sr * 1000
  suppression_time <- sum(result_df$duration_ms[result_df$type == "suppression"])
  bsr <- suppression_time / total_time * 100

  attr(result_df, "bsr") <- bsr
  rownames(result_df) <- NULL
  result_df
}


#' EEG Slowing Detection
#'
#' Detects pathological EEG slowing using spectral analysis. Supports three
#' methods: theta/delta ratio (TDR), delta-theta/alpha-beta ratio (DTAR),
#' and peak frequency analysis. Each channel is classified into normal,
#' mild, moderate, or severe slowing categories.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Analysis method: \code{"tdr"} (theta/delta ratio),
#'   \code{"dtar"} (delta+theta over alpha+beta ratio), or
#'   \code{"peak_frequency"} (dominant frequency per channel).
#' @param bands Named list of frequency bands for TDR and DTAR methods.
#'   Defaults to standard EEG bands.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{metric}{Character name of the metric used.}
#'     \item{value}{Numeric value of the metric.}
#'     \item{classification}{Character: \code{"normal"},
#'       \code{"mild_slowing"}, \code{"moderate_slowing"}, or
#'       \code{"severe_slowing"}.}
#'   }
#' @references
#' Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
#' clinical EEG. Electroencephalography and Clinical Neurophysiology, 106(3), 259-261.
#' @seealso [eegSpikeDetect()], [eegQEEG()], [eegAsymmetry()],
#'   [eegSuppression()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' result <- eegSlowing(pe, method = "dtar")
#' print(result)
#' }
eegSlowing <- function(x, method = c("tdr", "dtar", "peak_frequency"),
                       bands = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(bands)) {
    bands <- list(
      delta = c(1, 4),
      theta = c(4, 8),
      alpha = c(8, 13),
      beta  = c(13, 30)
    )
  }
  stopifnot(is.list(bands))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  results <- list()

  if (method == "tdr") {
    # Theta/Delta ratio per channel
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      delta_power <- .compute_band_power(sig, sr, bands$delta[1], bands$delta[2])
      theta_power <- .compute_band_power(sig, sr, bands$theta[1], bands$theta[2])

      tdr <- if (delta_power > .Machine$double.eps) {
        theta_power / delta_power
      } else {
        Inf
      }

      # Classification: TDR > 1 suggests theta dominance (mild slowing)
      classification <- if (tdr > 1) "mild_slowing" else "normal"

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        metric = "tdr",
        value = tdr,
        classification = classification,
        stringsAsFactors = FALSE
      )
    }

  } else if (method == "dtar") {
    # (Delta + Theta) / (Alpha + Beta) ratio per channel
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      delta_power <- .compute_band_power(sig, sr, bands$delta[1], bands$delta[2])
      theta_power <- .compute_band_power(sig, sr, bands$theta[1], bands$theta[2])
      alpha_power <- .compute_band_power(sig, sr, bands$alpha[1], bands$alpha[2])
      beta_power <- .compute_band_power(sig, sr, bands$beta[1], bands$beta[2])

      fast_power <- alpha_power + beta_power
      dtar <- if (fast_power > .Machine$double.eps) {
        (delta_power + theta_power) / fast_power
      } else {
        Inf
      }

      # Classification thresholds
      classification <- if (dtar < 1) {
        "normal"
      } else if (dtar < 2) {
        "mild_slowing"
      } else if (dtar < 4) {
        "moderate_slowing"
      } else {
        "severe_slowing"
      }

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        metric = "dtar",
        value = dtar,
        classification = classification,
        stringsAsFactors = FALSE
      )
    }

  } else if (method == "peak_frequency") {
    # Find dominant frequency per channel
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      n <- length(sig)
      ft <- fft(sig)
      power <- (Mod(ft)^2) / n
      freqs <- (0:(n - 1)) * sr / n
      half_n <- n %/% 2

      # Only consider frequencies from 1 to 30 Hz
      valid_idx <- which(freqs[1:half_n] >= 1 & freqs[1:half_n] <= 30)
      if (length(valid_idx) == 0) {
        peak_freq <- NA_real_
      } else {
        peak_idx <- valid_idx[which.max(power[valid_idx])]
        peak_freq <- freqs[peak_idx]
      }

      # Classification
      classification <- if (is.na(peak_freq)) {
        "normal"
      } else if (peak_freq >= 8) {
        "normal"
      } else if (peak_freq >= 6) {
        "mild_slowing"
      } else if (peak_freq >= 4) {
        "moderate_slowing"
      } else {
        "severe_slowing"
      }

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        metric = "peak_frequency",
        value = if (is.na(peak_freq)) NA_real_ else peak_freq,
        classification = classification,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(channel = integer(0), metric = character(0),
               value = numeric(0), classification = character(0))
  }
}


# --- Internal helpers ---

#' Compute band power for a single signal via FFT
#'
#' @param signal Numeric vector of the input signal.
#' @param sr Sampling rate in Hz.
#' @param low Low frequency bound in Hz.
#' @param high High frequency bound in Hz.
#' @return Numeric scalar of the summed spectral power in the band.
#' @keywords internal
.compute_band_power <- function(signal, sr, low, high) {
  n <- length(signal)
  ft <- fft(signal)
  power <- (Mod(ft)^2) / n
  freqs <- (0:(n - 1)) * sr / n
  half_n <- n %/% 2
  band_idx <- which(freqs[1:half_n] >= low & freqs[1:half_n] <= high)
  sum(power[band_idx])
}
