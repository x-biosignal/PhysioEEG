#' Automatic Sleep Staging
#'
#' Classifies EEG epochs into sleep stages (Wake, N1, N2, N3, REM) using
#' spectral power analysis based on simplified AASM criteria. Each epoch is
#' scored by computing spectral band powers (delta, theta, alpha, sigma, beta)
#' and applying rule-based classification.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Staging method: \code{"spectral"} (FFT power-based) or
#'   \code{"rule_based"} (synonym, same algorithm).
#' @param epoch_sec Epoch duration in seconds (default: 30, AASM standard).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{epoch}{Integer epoch number.}
#'     \item{stage}{Character sleep stage: "W", "N1", "N2", "N3", or "REM".}
#'     \item{start_sample}{Integer start sample of the epoch.}
#'     \item{end_sample}{Integer end sample of the epoch.}
#'     \item{delta_power}{Numeric delta band power (0.5-4 Hz).}
#'     \item{theta_power}{Numeric theta band power (4-8 Hz).}
#'     \item{alpha_power}{Numeric alpha band power (8-13 Hz).}
#'     \item{sigma_power}{Numeric sigma band power (12-16 Hz).}
#'     \item{beta_power}{Numeric beta band power (16-30 Hz).}
#'   }
#'
#'   The result is also stored in \code{metadata(x)$sleep_stages}.
#' @references
#' Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
#' Journal of Clinical Sleep Medicine, 13(5), 665-666.
#' @seealso [eegSpindleDetect()], [eegKcomplexDetect()],
#'   [eegSlowWaveDetect()], [eegSleepMetrics()], [eegPlotHypnogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 90000, n_channels = 2, sr = 500)
#' stages <- eegSleepStage(pe, epoch_sec = 30)
#' head(stages)
#' }
eegSleepStage <- function(x, method = c("spectral", "rule_based"),
                           epoch_sec = 30, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  epoch_samples <- as.integer(epoch_sec * sr)
  n_epochs <- n_time %/% epoch_samples

  if (n_epochs < 1) {
    stop(sprintf(
      "Signal too short for %d-second epochs (need at least %d samples, have %d).",
      epoch_sec, epoch_samples, n_time
    ), call. = FALSE)
  }

  results <- vector("list", n_epochs)

  for (ep in seq_len(n_epochs)) {
    start_sample <- (ep - 1L) * epoch_samples + 1L
    end_sample <- ep * epoch_samples

    # Average band powers across channels for classification
    delta_total <- 0
    theta_total <- 0
    alpha_total <- 0
    sigma_total <- 0
    beta_total <- 0
    total_total <- 0

    for (ch in seq_len(n_channels)) {
      seg <- data[start_sample:end_sample, ch]
      n <- length(seg)

      # Compute periodogram: |FFT|^2
      ft <- fft(seg)
      power <- (Mod(ft)^2) / n
      freqs <- (0:(n - 1)) * sr / n

      # Only use first half of spectrum (positive frequencies)
      half_n <- n %/% 2

      # Sum power in each band's frequency bins
      delta_idx <- which(freqs[1:half_n] >= 0.5 & freqs[1:half_n] <= 4)
      theta_idx <- which(freqs[1:half_n] >= 4 & freqs[1:half_n] < 8)
      alpha_idx <- which(freqs[1:half_n] >= 8 & freqs[1:half_n] <= 13)
      sigma_idx <- which(freqs[1:half_n] >= 12 & freqs[1:half_n] <= 16)
      beta_idx <- which(freqs[1:half_n] >= 16 & freqs[1:half_n] <= 30)

      delta_total <- delta_total + sum(power[delta_idx])
      theta_total <- theta_total + sum(power[theta_idx])
      alpha_total <- alpha_total + sum(power[alpha_idx])
      sigma_total <- sigma_total + sum(power[sigma_idx])
      beta_total <- beta_total + sum(power[beta_idx])
      total_total <- total_total + sum(power[1:half_n])
    }

    # Average across channels
    delta_power <- delta_total / n_channels
    theta_power <- theta_total / n_channels
    alpha_power <- alpha_total / n_channels
    sigma_power <- sigma_total / n_channels
    beta_power <- beta_total / n_channels
    total_power <- total_total / n_channels

    # Guard against zero total power
    if (total_power < .Machine$double.eps) {
      total_power <- 1
    }

    # Classification rules (AASM simplified)
    # Apply in order: N3 first (strongest signal), then W, N2, REM, N1 (default)
    stage <- "N1"  # default

    if (delta_power / total_power > 0.4) {
      stage <- "N3"
    } else if (alpha_power / total_power > 0.3) {
      stage <- "W"
    } else if (sigma_power / total_power > 0.05) {
      stage <- "N2"
    } else if (theta_power > delta_power && alpha_power < 0.2 * total_power) {
      stage <- "REM"
    }
    # else remains N1

    results[[ep]] <- data.frame(
      epoch = ep,
      stage = stage,
      start_sample = start_sample,
      end_sample = end_sample,
      delta_power = delta_power,
      theta_power = theta_power,
      alpha_power = alpha_power,
      sigma_power = sigma_power,
      beta_power = beta_power,
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)

  # Store in metadata
  md <- S4Vectors::metadata(x)
  md$sleep_stages <- result_df
  S4Vectors::metadata(x) <- md

  result_df
}


#' Detect Sleep Spindles
#'
#' Identifies sleep spindles in EEG data by bandpass filtering in the sigma
#' frequency range, computing an RMS envelope, and detecting continuous
#' segments exceeding a threshold. Spindle frequency is estimated from
#' zero crossings in the bandpassed signal.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Detection method: \code{"sigma"} (sigma-band filtering) or
#'   \code{"wavelet"} (synonym, same algorithm).
#' @param freq_range Numeric vector of length 2 specifying the spindle
#'   frequency range in Hz (default: \code{c(11, 16)}).
#' @param min_duration_ms Minimum spindle duration in milliseconds
#'   (default: 500).
#' @param max_duration_ms Maximum spindle duration in milliseconds
#'   (default: 2000).
#' @param threshold_sd Number of standard deviations above the mean for
#'   detection threshold (default: 1.5).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{start_sample}{Integer start sample of the spindle.}
#'     \item{end_sample}{Integer end sample of the spindle.}
#'     \item{duration_ms}{Numeric spindle duration in milliseconds.}
#'     \item{peak_sample}{Integer sample index of peak amplitude.}
#'     \item{peak_amplitude}{Numeric peak RMS amplitude.}
#'     \item{frequency_hz}{Numeric estimated spindle frequency in Hz.}
#'   }
#' @references
#' Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
#' Journal of Clinical Sleep Medicine, 13(5), 665-666.
#'
#' Molle, M., et al. (2011). Fast and slow spindles during the sleep slow
#' oscillation. Sleep, 34(10), 1411-1421.
#' @seealso [eegSleepStage()], [eegKcomplexDetect()],
#'   [eegSlowWaveDetect()], [eegSleepMetrics()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
#' spindles <- eegSpindleDetect(pe)
#' head(spindles)
#' }
eegSpindleDetect <- function(x, method = c("sigma", "wavelet"),
                              freq_range = c(11, 16),
                              min_duration_ms = 500,
                              max_duration_ms = 2000,
                              threshold_sd = 1.5,
                              assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(freq_range) && length(freq_range) == 2)
  stopifnot(freq_range[1] < freq_range[2])

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  min_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))
  max_samples <- as.integer(round(max_duration_ms / 1000 * sr))
  rms_window <- max(1L, as.integer(round(0.200 * sr)))  # 200ms window
  half_win <- rms_window %/% 2

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # Step 1: Bandpass filter in freq_range using FIR filter
    bp_sig <- .fir_bandpass(sig, sr, freq_range[1], freq_range[2])

    # Step 2: Compute RMS envelope with 200ms window
    sig_sq <- bp_sig^2
    envelope <- numeric(n_time)
    for (i in seq_len(n_time)) {
      lo <- max(1L, i - half_win)
      hi <- min(n_time, i + half_win)
      envelope[i] <- sqrt(sum(sig_sq[lo:hi]) / (hi - lo + 1))
    }

    # Step 3: Threshold = mean + threshold_sd * SD
    env_mean <- mean(envelope)
    env_sd <- sd(envelope)
    threshold <- env_mean + threshold_sd * env_sd

    # Step 4: Find continuous segments above threshold
    above <- envelope > threshold
    transitions <- diff(as.integer(above))
    onset_samples <- which(transitions == 1) + 1L
    offset_samples <- which(transitions == -1)

    if (above[1]) onset_samples <- c(1L, onset_samples)
    if (above[n_time]) offset_samples <- c(offset_samples, n_time)

    n_events <- min(length(onset_samples), length(offset_samples))

    if (n_events > 0) {
      for (e in seq_len(n_events)) {
        start_s <- onset_samples[e]
        end_s <- offset_samples[e]
        duration_samples <- end_s - start_s + 1L

        # Step 5: Filter by duration criteria
        if (duration_samples < min_samples || duration_samples > max_samples) next

        duration_ms <- duration_samples / sr * 1000

        # Step 6: Find peak amplitude and its sample index
        seg_env <- envelope[start_s:end_s]
        peak_offset <- which.max(seg_env)
        peak_sample <- start_s + peak_offset - 1L
        peak_amplitude <- seg_env[peak_offset]

        # Step 7: Estimate frequency from zero crossings
        bp_seg <- bp_sig[start_s:end_s]
        sign_changes <- sum(abs(diff(sign(bp_seg))) > 0)
        # Each full cycle = 2 zero crossings
        duration_sec <- duration_samples / sr
        frequency_hz <- if (duration_sec > 0) {
          (sign_changes / 2) / duration_sec
        } else {
          NA_real_
        }

        results[[length(results) + 1]] <- data.frame(
          channel = as.integer(ch),
          start_sample = start_s,
          end_sample = end_s,
          duration_ms = duration_ms,
          peak_sample = peak_sample,
          peak_amplitude = peak_amplitude,
          frequency_hz = frequency_hz,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel = integer(0), start_sample = integer(0),
      end_sample = integer(0), duration_ms = numeric(0),
      peak_sample = integer(0), peak_amplitude = numeric(0),
      frequency_hz = numeric(0)
    )
  }
}


#' Detect K-Complexes
#'
#' Identifies K-complexes in EEG data by lowpass filtering at 4 Hz,
#' finding negative peaks exceeding a threshold amplitude, and verifying
#' the characteristic negative-positive waveform morphology.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param min_neg_amplitude Minimum absolute negative peak amplitude in
#'   microvolts (default: 75). Peaks must be more negative than
#'   \code{-min_neg_amplitude}.
#' @param min_duration_ms Minimum K-complex duration in milliseconds
#'   (default: 500).
#' @param max_duration_ms Maximum K-complex duration in milliseconds
#'   (default: 1500).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{negative_peak_sample}{Integer sample of the negative peak.}
#'     \item{positive_peak_sample}{Integer sample of the positive peak.}
#'     \item{negative_amplitude}{Numeric amplitude at the negative peak.}
#'     \item{positive_amplitude}{Numeric amplitude at the positive peak.}
#'     \item{duration_ms}{Numeric total duration in milliseconds.}
#'   }
#' @references
#' Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
#' Journal of Clinical Sleep Medicine, 13(5), 665-666.
#' @seealso [eegSleepStage()], [eegSpindleDetect()],
#'   [eegSlowWaveDetect()], [eegSleepMetrics()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
#' kcomplexes <- eegKcomplexDetect(pe)
#' head(kcomplexes)
#' }
eegKcomplexDetect <- function(x, min_neg_amplitude = 75,
                               min_duration_ms = 500,
                               max_duration_ms = 1500,
                               assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(min_neg_amplitude) && min_neg_amplitude > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  min_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))
  max_samples <- as.integer(round(max_duration_ms / 1000 * sr))

  # Window for searching positive peak after negative peak: 0.5-1.5 sec
  pos_search_min <- max(1L, as.integer(round(0.5 * sr)))
  pos_search_max <- as.integer(round(1.5 * sr))

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # Step 1: Lowpass filter at 4 Hz using FIR filter
    lp_sig <- .fir_lowpass(sig, sr, 4)

    # Step 2: Find negative peaks exceeding -min_neg_amplitude
    neg_peaks <- integer(0)
    for (i in 2:(n_time - 1)) {
      if (lp_sig[i] < lp_sig[i - 1] && lp_sig[i] <= lp_sig[i + 1] &&
          lp_sig[i] < -min_neg_amplitude) {
        neg_peaks <- c(neg_peaks, i)
      }
    }

    # Step 3: For each negative peak, look for positive peak within window
    for (np in neg_peaks) {
      search_start <- np + pos_search_min
      search_end <- min(n_time, np + pos_search_max)

      if (search_start > n_time || search_start >= search_end) next

      seg <- lp_sig[search_start:search_end]
      if (length(seg) == 0) next

      pos_offset <- which.max(seg)
      pos_sample <- search_start + pos_offset - 1L
      pos_amplitude <- seg[pos_offset]

      # Positive peak must actually be positive
      if (pos_amplitude <= 0) next

      # Step 4: Check total duration of negative-positive complex
      duration_samples <- pos_sample - np
      if (duration_samples < min_samples || duration_samples > max_samples) next

      duration_ms <- duration_samples / sr * 1000

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        negative_peak_sample = np,
        positive_peak_sample = pos_sample,
        negative_amplitude = lp_sig[np],
        positive_amplitude = pos_amplitude,
        duration_ms = duration_ms,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel = integer(0), negative_peak_sample = integer(0),
      positive_peak_sample = integer(0), negative_amplitude = numeric(0),
      positive_amplitude = numeric(0), duration_ms = numeric(0)
    )
  }
}


#' Detect Slow Waves
#'
#' Identifies slow oscillations in EEG data by bandpass filtering in the
#' slow wave frequency range, finding zero crossings, and measuring
#' negative half-wave amplitudes and slopes.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param min_amplitude Minimum absolute negative peak amplitude in
#'   microvolts (default: 75).
#' @param freq_range Numeric vector of length 2 specifying the slow wave
#'   frequency range in Hz (default: \code{c(0.5, 2)}).
#' @param min_duration_ms Minimum half-wave duration in milliseconds
#'   (default: 250).
#' @param max_duration_ms Maximum half-wave duration in milliseconds
#'   (default: 1000).
#' @param assay_name Input assay name. If \code{NULL}, uses the default assay.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{start_sample}{Integer sample at first zero crossing.}
#'     \item{end_sample}{Integer sample at second zero crossing.}
#'     \item{negative_peak}{Numeric negative peak amplitude.}
#'     \item{positive_peak}{Numeric positive peak amplitude.}
#'     \item{duration_ms}{Numeric half-wave duration in milliseconds.}
#'     \item{slope}{Numeric slope from negative to positive peak
#'       (microvolts per millisecond).}
#'   }
#' @references
#' Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
#' Journal of Clinical Sleep Medicine, 13(5), 665-666.
#' @seealso [eegSleepStage()], [eegSpindleDetect()],
#'   [eegKcomplexDetect()], [eegSleepMetrics()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
#' slow_waves <- eegSlowWaveDetect(pe)
#' head(slow_waves)
#' }
eegSlowWaveDetect <- function(x, min_amplitude = 75,
                               freq_range = c(0.5, 2),
                               min_duration_ms = 250,
                               max_duration_ms = 1000,
                               assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(freq_range) && length(freq_range) == 2)
  stopifnot(freq_range[1] < freq_range[2])
  stopifnot(is.numeric(min_amplitude) && min_amplitude > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  min_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))
  max_samples <- as.integer(round(max_duration_ms / 1000 * sr))

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # Step 1: Bandpass filter in freq_range using FIR filter
    bp_sig <- .fir_bandpass(sig, sr, freq_range[1], freq_range[2])

    # Step 2: Find zero crossings (sign changes)
    signs <- sign(bp_sig)
    # Replace zeros with 1 to avoid ambiguity
    signs[signs == 0] <- 1
    sign_diff <- diff(signs)
    zero_crossings <- which(sign_diff != 0)

    if (length(zero_crossings) < 2) next

    # Step 3: For each negative-to-positive half-wave
    for (i in seq_len(length(zero_crossings) - 1)) {
      zc_start <- zero_crossings[i]
      zc_end <- zero_crossings[i + 1]

      # We want the segment between crossings to be negative (going negative)
      seg <- bp_sig[(zc_start + 1):zc_end]
      if (length(seg) == 0) next

      neg_peak_val <- min(seg)

      # Check if this is a negative half-wave with sufficient amplitude
      if (neg_peak_val > -min_amplitude) next

      duration_samples <- zc_end - zc_start
      if (duration_samples < min_samples || duration_samples > max_samples) next

      neg_peak_offset <- which.min(seg)
      duration_ms <- duration_samples / sr * 1000

      # Find the positive peak in the next half-wave (if available)
      if (i + 1 < length(zero_crossings)) {
        next_zc <- zero_crossings[i + 2]
        if (next_zc <= n_time) {
          pos_seg <- bp_sig[(zc_end + 1):min(next_zc, n_time)]
          pos_peak_val <- if (length(pos_seg) > 0) max(pos_seg) else 0
        } else {
          pos_peak_val <- 0
        }
      } else {
        # Use remaining signal after last zero crossing
        remaining <- bp_sig[(zc_end + 1):min(n_time, zc_end + max_samples)]
        pos_peak_val <- if (length(remaining) > 0) max(remaining) else 0
      }

      # Step 4: Compute slope: (positive_peak - negative_peak) / duration
      slope <- (pos_peak_val - neg_peak_val) / duration_ms

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        start_sample = zc_start,
        end_sample = zc_end,
        negative_peak = neg_peak_val,
        positive_peak = pos_peak_val,
        duration_ms = duration_ms,
        slope = slope,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel = integer(0), start_sample = integer(0),
      end_sample = integer(0), negative_peak = numeric(0),
      positive_peak = numeric(0), duration_ms = numeric(0),
      slope = numeric(0)
    )
  }
}


#' Compute Sleep Metrics
#'
#' Calculates summary sleep architecture metrics from staged sleep data.
#' Requires that \code{eegSleepStage()} has been run and results are stored
#' in \code{metadata(x)$sleep_stages}.
#'
#' @param x A PhysioExperiment object with sleep staging results in
#'   \code{metadata(x)$sleep_stages}.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{metric}{Character name of the sleep metric.}
#'     \item{value}{Numeric value of the metric.}
#'     \item{unit}{Character unit of measurement.}
#'   }
#'
#'   Metrics include:
#'   \itemize{
#'     \item \code{total_sleep_time}: Time in N1+N2+N3+REM (minutes)
#'     \item \code{sleep_efficiency}: TST / total recording time * 100 (percent)
#'     \item \code{waso}: Wake after sleep onset (minutes)
#'     \item \code{sleep_latency}: Time to first non-Wake epoch (minutes)
#'     \item \code{rem_latency}: Time from first sleep to first REM (minutes)
#'     \item \code{pct_N1}: Percentage of TST in N1
#'     \item \code{pct_N2}: Percentage of TST in N2
#'     \item \code{pct_N3}: Percentage of TST in N3
#'     \item \code{pct_REM}: Percentage of TST in REM
#'     \item \code{pct_W}: Percentage of total time in Wake
#'   }
#' @references
#' Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
#' Journal of Clinical Sleep Medicine, 13(5), 665-666.
#' @seealso [eegSleepStage()], [eegSpindleDetect()],
#'   [eegKcomplexDetect()], [eegSlowWaveDetect()], [eegPlotHypnogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
#' stages <- eegSleepStage(pe, epoch_sec = 30)
#' metadata(pe)$sleep_stages <- stages
#' metrics <- eegSleepMetrics(pe)
#' print(metrics)
#' }
eegSleepMetrics <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))

  stages_df <- S4Vectors::metadata(x)$sleep_stages
  if (is.null(stages_df)) {
    stop("No sleep staging results found in metadata. Run eegSleepStage() first.",
         call. = FALSE)
  }

  stopifnot(is.data.frame(stages_df))
  stopifnot("stage" %in% names(stages_df))
  stopifnot("epoch" %in% names(stages_df))

  n_epochs <- nrow(stages_df)
  stage_vec <- stages_df$stage

  # Determine epoch duration from sample info if available
  if (all(c("start_sample", "end_sample") %in% names(stages_df))) {
    sr <- samplingRate(x)
    epoch_samples <- stages_df$end_sample[1] - stages_df$start_sample[1] + 1
    epoch_sec <- epoch_samples / sr
  } else {
    # Assume standard 30-second epochs
    epoch_sec <- 30
  }

  epoch_min <- epoch_sec / 60

  # Total recording time
  total_time_min <- n_epochs * epoch_min

  # Sleep epochs = non-Wake epochs
  sleep_mask <- stage_vec != "W"

  # Total sleep time (N1+N2+N3+REM)
  tst_epochs <- sum(sleep_mask)
  tst_min <- tst_epochs * epoch_min

  # Sleep efficiency
  sleep_efficiency <- if (total_time_min > 0) {
    tst_min / total_time_min * 100
  } else {
    NA_real_
  }

  # Sleep latency: time to first non-W epoch
  first_sleep_idx <- which(sleep_mask)[1]
  sleep_latency_min <- if (!is.na(first_sleep_idx)) {
    (first_sleep_idx - 1) * epoch_min
  } else {
    total_time_min  # never fell asleep
  }

  # WASO: Wake epochs after first sleep epoch
  if (!is.na(first_sleep_idx) && first_sleep_idx < n_epochs) {
    wake_after_onset <- sum(stage_vec[(first_sleep_idx + 1):n_epochs] == "W")
    waso_min <- wake_after_onset * epoch_min
  } else {
    waso_min <- 0
  }

  # REM latency: time from first sleep to first REM
  first_rem_idx <- which(stage_vec == "REM")[1]
  rem_latency_min <- if (!is.na(first_rem_idx) && !is.na(first_sleep_idx)) {
    (first_rem_idx - first_sleep_idx) * epoch_min
  } else {
    NA_real_
  }

  # Stage percentages (of TST)
  n1_count <- sum(stage_vec == "N1")
  n2_count <- sum(stage_vec == "N2")
  n3_count <- sum(stage_vec == "N3")
  rem_count <- sum(stage_vec == "REM")
  w_count <- sum(stage_vec == "W")

  pct_n1 <- if (tst_epochs > 0) n1_count / tst_epochs * 100 else NA_real_
  pct_n2 <- if (tst_epochs > 0) n2_count / tst_epochs * 100 else NA_real_
  pct_n3 <- if (tst_epochs > 0) n3_count / tst_epochs * 100 else NA_real_
  pct_rem <- if (tst_epochs > 0) rem_count / tst_epochs * 100 else NA_real_
  pct_w <- if (n_epochs > 0) w_count / n_epochs * 100 else NA_real_

  data.frame(
    metric = c("total_sleep_time", "sleep_efficiency", "waso",
               "sleep_latency", "rem_latency",
               "pct_N1", "pct_N2", "pct_N3", "pct_REM", "pct_W"),
    value = c(tst_min, sleep_efficiency, waso_min,
              sleep_latency_min, rem_latency_min,
              pct_n1, pct_n2, pct_n3, pct_rem, pct_w),
    unit = c("minutes", "percent", "minutes",
             "minutes", "minutes",
             "percent", "percent", "percent", "percent", "percent"),
    stringsAsFactors = FALSE
  )
}


# --- Internal helpers ---

#' Windowed-sinc FIR bandpass filter
#'
#' Applies a windowed-sinc FIR bandpass filter using a Hamming window.
#' The filter is applied causally via \code{stats::filter} with group-delay
#' compensation. This replaces the previous brick-wall FFT filter which
#' produced Gibbs ringing artifacts.
#'
#' @param signal Numeric vector of the input signal.
#' @param sr Sampling rate in Hz.
#' @param low Low cutoff frequency in Hz.
#' @param high High cutoff frequency in Hz.
#' @param order Filter order. If \code{NULL}, auto-selected as approximately
#'   3 cycles of the lowest frequency, clamped to signal length and forced odd.
#' @return Numeric vector of the filtered signal (same length as input).
#' @keywords internal
.fir_bandpass <- function(signal, sr, low, high, order = NULL) {
  n <- length(signal)
  if (is.null(order)) {
    # Auto-select order: ~3 cycles of the lowest frequency
    order <- as.integer(round(3 * sr / max(low, 0.5)))
    order <- min(order, n - 1L)
    if (order %% 2 == 0) order <- order + 1L  # ensure odd
  }
  order <- min(order, n - 1L)
  if (order < 3L) return(signal)
  if (order %% 2 == 0) order <- order + 1L

  half <- (order - 1L) %/% 2L
  m <- seq(-half, half)

  # Normalised cutoff frequencies (fraction of sampling rate)
  fc_low <- low / sr
  fc_high <- high / sr

  # Sinc functions: lowpass(high) - lowpass(low) = bandpass

  h_high <- ifelse(m == 0, 2 * fc_high, sin(2 * pi * fc_high * m) / (pi * m))
  h_low  <- ifelse(m == 0, 2 * fc_low,  sin(2 * pi * fc_low * m) / (pi * m))
  h <- h_high - h_low

  # Apply Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * (0:(order - 1)) / (order - 1))
  h <- h * w

  # Normalise to unit gain at centre frequency
  centre_freq <- (low + high) / 2
  centre_response <- sum(h * cos(2 * pi * centre_freq / sr * m))
  if (abs(centre_response) > 1e-10) h <- h / centre_response

  # Apply via causal convolution then compensate group delay
  pad_len <- order
  sig_padded <- c(rep(signal[1], pad_len), signal, rep(signal[n], pad_len))
  filtered <- stats::filter(sig_padded, h, sides = 1)
  filtered <- as.numeric(filtered)

  # Remove padding and compensate for group delay (half filter length)
  result <- filtered[(pad_len + half + 1):(pad_len + half + n)]

  # Handle any remaining NAs at edges
  na_idx <- which(is.na(result))
  if (length(na_idx) > 0) {
    result[na_idx] <- signal[na_idx]
  }

  result
}


#' Windowed-sinc FIR lowpass filter
#'
#' Applies a windowed-sinc FIR lowpass filter using a Hamming window.
#' The filter is applied causally via \code{stats::filter} with group-delay
#' compensation. This replaces the previous brick-wall FFT filter which
#' produced Gibbs ringing artifacts.
#'
#' @param signal Numeric vector of the input signal.
#' @param sr Sampling rate in Hz.
#' @param cutoff Cutoff frequency in Hz.
#' @param order Filter order. If \code{NULL}, auto-selected as approximately
#'   3 cycles of the cutoff frequency, clamped to signal length and forced odd.
#' @return Numeric vector of the filtered signal (same length as input).
#' @keywords internal
.fir_lowpass <- function(signal, sr, cutoff, order = NULL) {
  n <- length(signal)
  if (is.null(order)) {
    order <- as.integer(round(3 * sr / max(cutoff, 0.5)))
    order <- min(order, n - 1L)
    if (order %% 2 == 0) order <- order + 1L
  }
  order <- min(order, n - 1L)
  if (order < 3L) return(signal)
  if (order %% 2 == 0) order <- order + 1L

  half <- (order - 1L) %/% 2L
  m <- seq(-half, half)

  fc <- cutoff / sr
  h <- ifelse(m == 0, 2 * fc, sin(2 * pi * fc * m) / (pi * m))

  # Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * (0:(order - 1)) / (order - 1))
  h <- h * w
  h <- h / sum(h)  # normalise to unity DC gain

  pad_len <- order
  sig_padded <- c(rep(signal[1], pad_len), signal, rep(signal[n], pad_len))
  filtered <- stats::filter(sig_padded, h, sides = 1)
  filtered <- as.numeric(filtered)
  result <- filtered[(pad_len + half + 1):(pad_len + half + n)]

  na_idx <- which(is.na(result))
  if (length(na_idx) > 0) result[na_idx] <- signal[na_idx]
  result
}
