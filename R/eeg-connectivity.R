#' EEG Connectivity Analysis
#'
#' Functions for computing frequency-domain connectivity measures between
#' EEG channels, including coherence, phase locking value (PLV), weighted
#' phase lag index (wPLI), and spectral Granger causality.
#'
#' @name eeg-connectivity
#' @keywords internal
NULL


# --- Internal helpers ---

#' Compute instantaneous phase via Hilbert transform
#'
#' Constructs the analytic signal by zeroing negative frequency components
#' of the FFT and doubling positive frequencies, then returns the
#' instantaneous phase (argument of the analytic signal).
#'
#' @param signal Numeric vector of the input signal.
#' @return Numeric vector of instantaneous phase in radians (-pi, pi].
#' @keywords internal
.hilbert_phase <- function(signal) {
  n <- length(signal)
  ft <- fft(signal)
  h <- rep(0, n)
  h[1] <- 1
  if (n %% 2 == 0) {
    h[n / 2 + 1] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[2:((n + 1) / 2)] <- 2
  }
  analytic <- fft(ft * h, inverse = TRUE) / n
  Arg(analytic)
}


#' Compute Hanning window
#'
#' @param n Window length.
#' @return Numeric vector of Hanning window coefficients.
#' @keywords internal
.hanning_window <- function(n) {
  0.5 * (1 - cos(2 * pi * seq(0, n - 1) / (n - 1)))
}


#' Compute cross-spectral density between two signals using Welch's method
#'
#' Segments the signals, applies a Hanning window, computes FFTs, and
#' averages the cross-spectrum and auto-spectra across windows.
#'
#' @param x Numeric vector, first signal.
#' @param y Numeric vector, second signal.
#' @param sr Sampling rate in Hz.
#' @param window_sec Window length in seconds.
#' @param overlap Overlap fraction (0 to 1).
#' @return A list with components:
#'   \describe{
#'     \item{Sxy}{Complex vector of cross-spectral density (positive freqs).}
#'     \item{Sxx}{Numeric vector of auto-spectral density of x.}
#'     \item{Syy}{Numeric vector of auto-spectral density of y.}
#'     \item{freqs}{Numeric vector of frequencies in Hz.}
#'     \item{n_windows}{Integer number of windows averaged.}
#'     \item{Sxy_segments}{List of per-window cross-spectra (complex vectors).}
#'   }
#' @keywords internal
.welch_cross_spectrum <- function(x, y, sr, window_sec, overlap) {
  n <- length(x)
  window_len <- as.integer(round(window_sec * sr))
  if (window_len > n) window_len <- n
  # Force even length for symmetric frequency axis
  if (window_len %% 2 != 0) window_len <- window_len + 1L
  if (window_len > n) window_len <- window_len - 1L

  hop <- max(1L, as.integer(round(window_len * (1 - overlap))))
  win <- .hanning_window(window_len)
  win_norm <- sum(win^2)

  n_fft_bins <- window_len %/% 2 + 1L
  freqs <- seq(0, sr / 2, length.out = n_fft_bins)

  starts <- seq(1L, n - window_len + 1L, by = hop)
  n_windows <- length(starts)

  Sxy_sum <- complex(n_fft_bins)
  Sxx_sum <- numeric(n_fft_bins)
  Syy_sum <- numeric(n_fft_bins)
  Sxy_segments <- vector("list", n_windows)

  for (wi in seq_len(n_windows)) {
    idx <- starts[wi]:(starts[wi] + window_len - 1L)
    x_win <- x[idx] * win
    y_win <- y[idx] * win

    fx <- fft(x_win)[seq_len(n_fft_bins)]
    fy <- fft(y_win)[seq_len(n_fft_bins)]

    sxy <- fx * Conj(fy)
    sxx <- Mod(fx)^2
    syy <- Mod(fy)^2

    Sxy_sum <- Sxy_sum + sxy
    Sxx_sum <- Sxx_sum + sxx
    Syy_sum <- Syy_sum + syy
    Sxy_segments[[wi]] <- sxy / win_norm
  }

  list(
    Sxy = Sxy_sum / (n_windows * win_norm),
    Sxx = Sxx_sum / (n_windows * win_norm),
    Syy = Syy_sum / (n_windows * win_norm),
    freqs = freqs,
    n_windows = n_windows,
    Sxy_segments = Sxy_segments
  )
}


#' EEG Coherence Analysis
#'
#' Computes magnitude-squared coherence (MSC) or imaginary part of coherency
#' between all EEG channel pairs using Welch's method. The signal is segmented
#' into overlapping windows, each windowed with a Hanning function, and the
#' cross-spectral density is averaged across windows.
#'
#' Magnitude-squared coherence is defined as:
#' \deqn{MSC(f) = |S_{xy}(f)|^2 / (S_{xx}(f) \cdot S_{yy}(f))}
#'
#' Imaginary coherence uses the imaginary part of coherency to reduce
#' volume conduction artifacts:
#' \deqn{ICoh(f) = |Im(S_{xy}(f) / \sqrt{S_{xx}(f) \cdot S_{yy}(f)})|}
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param method Coherence method: \code{"coherence"} for magnitude-squared
#'   coherence or \code{"imaginary"} for imaginary part of coherency
#'   (default: \code{"coherence"}).
#' @param window_sec Window length in seconds for Welch's method (default: 2).
#' @param overlap Overlap fraction between adjacent windows, from 0 to 1
#'   exclusive (default: 0.5).
#' @param band Numeric vector of length 2 specifying the frequency band in Hz
#'   over which to average coherence (default: \code{c(8, 13)} for alpha band).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return The input PhysioExperiment with connectivity results stored in
#'   \code{metadata(x)$connectivity}, a list containing:
#'   \describe{
#'     \item{matrix}{Numeric n_channels x n_channels matrix of band-averaged
#'       coherence values.}
#'     \item{method}{Character string indicating the method used.}
#'     \item{band}{Numeric vector of the frequency band used.}
#'     \item{freqs}{Numeric vector of frequency bins.}
#'     \item{spectra}{3D array (n_freqs x n_channels x n_channels) of
#'       frequency-resolved coherence values.}
#'   }
#' @references
#' Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain signals.
#' Human Brain Mapping, 8(4), 194-208.
#' @seealso [eegPLV()], [eegWPLI()], [eegConnectivityMatrix()],
#'   [eegPlotConnectivity()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' pe <- eegCoherence(pe, method = "coherence", band = c(8, 13))
#' coh_matrix <- metadata(pe)$connectivity$matrix
#' }
eegCoherence <- function(x, method = c("coherence", "imaginary"),
                         window_sec = 2, overlap = 0.5,
                         band = c(8, 13), assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(window_sec) && window_sec > 0)
  stopifnot(is.numeric(overlap) && overlap >= 0 && overlap < 1)
  stopifnot(is.numeric(band) && length(band) == 2 && band[1] < band[2])


  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  nyquist <- sr / 2
  if (band[2] > nyquist) {
    stop(sprintf(
      "band upper limit (%.1f Hz) exceeds Nyquist frequency (%.1f Hz).",
      band[2], nyquist
    ), call. = FALSE)
  }

  if (length(dim(data)) != 2) {
    stop("eegCoherence requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Compute window parameters for frequency axis
  window_len <- as.integer(round(window_sec * sr))
  if (window_len > n_time) window_len <- n_time
  if (window_len %% 2 != 0) window_len <- window_len + 1L
  if (window_len > n_time) window_len <- window_len - 1L
  n_fft_bins <- window_len %/% 2 + 1L
  freqs <- seq(0, sr / 2, length.out = n_fft_bins)

  # Allocate frequency-resolved coherence spectra
  spectra <- array(0, dim = c(n_fft_bins, n_channels, n_channels))

  # Band-averaged connectivity matrix
  conn_matrix <- matrix(0, nrow = n_channels, ncol = n_channels)

  # Frequency indices within the specified band
  band_idx <- which(freqs >= band[1] & freqs <= band[2])
  if (length(band_idx) == 0) {
    stop(sprintf(
      "No frequency bins in band [%.1f, %.1f] Hz. Check band or window_sec.",
      band[1], band[2]
    ), call. = FALSE)
  }

  # Set diagonal to 1 (self-coherence)
  diag(conn_matrix) <- 1
  for (fi in seq_len(n_fft_bins)) {
    for (ch in seq_len(n_channels)) {
      spectra[fi, ch, ch] <- 1
    }
  }

  # Compute pairwise coherence
  for (i in seq_len(n_channels - 1)) {
    for (j in (i + 1):n_channels) {
      csd <- .welch_cross_spectrum(
        data[, i], data[, j], sr, window_sec, overlap
      )

      if (method == "coherence") {
        # Magnitude-squared coherence: |Sxy|^2 / (Sxx * Syy)
        denom <- csd$Sxx * csd$Syy
        denom[denom < .Machine$double.eps] <- .Machine$double.eps
        coh_freq <- Mod(csd$Sxy)^2 / denom
      } else {
        # Imaginary coherence: |Im(Sxy / sqrt(Sxx * Syy))|
        denom <- sqrt(csd$Sxx * csd$Syy)
        denom[denom < .Machine$double.eps] <- .Machine$double.eps
        coherency <- csd$Sxy / denom
        coh_freq <- abs(Im(coherency))
      }

      # Store frequency-resolved coherence
      spectra[, i, j] <- coh_freq
      spectra[, j, i] <- coh_freq

      # Average over band
      band_mean <- mean(coh_freq[band_idx])
      conn_matrix[i, j] <- band_mean
      conn_matrix[j, i] <- band_mean
    }
  }

  # Add channel labels
  ch_labels <- NULL
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% names(cd)) {
    ch_labels <- as.character(cd$label)
    rownames(conn_matrix) <- ch_labels
    colnames(conn_matrix) <- ch_labels
  }

  # Store in metadata
  md <- S4Vectors::metadata(x)
  md$connectivity <- list(
    matrix = conn_matrix,
    method = method,
    band = band,
    freqs = freqs,
    spectra = spectra
  )
  S4Vectors::metadata(x) <- md

  x
}


#' EEG Phase Locking Value (PLV)
#'
#' Computes the Phase Locking Value between all EEG channel pairs within a
#' specified frequency band. Each channel is bandpass filtered, instantaneous
#' phase is extracted via the Hilbert transform, and PLV is computed as the
#' mean resultant length of the phase difference distribution.
#'
#' \deqn{PLV = |1/N \sum_{t=1}^{N} exp(i(\phi_x(t) - \phi_y(t)))|}
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param band Numeric vector of length 2 specifying the frequency band in Hz
#'   for bandpass filtering before phase extraction (default: \code{c(8, 13)}
#'   for alpha band).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel1}{Character or integer identifier of the first channel.}
#'     \item{channel2}{Character or integer identifier of the second channel.}
#'     \item{plv}{Numeric PLV value in \code{[0, 1]}.}
#'   }
#' @references
#' Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain signals.
#' Human Brain Mapping, 8(4), 194-208.
#' @seealso [eegCoherence()], [eegWPLI()], [eegConnectivityMatrix()],
#'   [eegPlotConnectivity()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' plv_df <- eegPLV(pe, band = c(8, 13))
#' head(plv_df)
#' }
eegPLV <- function(x, band = c(8, 13), assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(band) && length(band) == 2 && band[1] < band[2])

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 2) {
    stop("eegPLV requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Get channel labels
  ch_labels <- seq_len(n_channels)
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% names(cd)) {
    ch_labels <- as.character(cd$label)
  }

  # Bandpass filter and extract phase for each channel
  phases <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    bp_sig <- .fir_bandpass(data[, ch], sr, band[1], band[2])
    phases[, ch] <- .hilbert_phase(bp_sig)
  }

  # Compute PLV for all pairs
  results <- list()
  for (i in seq_len(n_channels - 1)) {
    for (j in (i + 1):n_channels) {
      phase_diff <- phases[, i] - phases[, j]
      plv_val <- Mod(mean(exp(1i * phase_diff)))

      results[[length(results) + 1]] <- data.frame(
        channel1 = ch_labels[i],
        channel2 = ch_labels[j],
        plv = plv_val,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel1 = character(0),
      channel2 = character(0),
      plv = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}


#' EEG Weighted Phase Lag Index (wPLI)
#'
#' Computes the Weighted Phase Lag Index and its debiased variant between all
#' EEG channel pairs. The wPLI reduces the influence of volume conduction by
#' weighting the phase differences by the magnitude of the imaginary part of
#' the cross-spectrum, computed via Welch's method.
#'
#' \deqn{wPLI = |mean(Im(S_{xy}))| / mean(|Im(S_{xy})|)}
#'
#' The debiased wPLI corrects for sample-size bias:
#' \deqn{wPLI^2_{debiased} = (N \cdot wPLI^2 - 1) / (N - 1)}
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param band Numeric vector of length 2 specifying the frequency band in Hz
#'   over which to compute wPLI (default: \code{c(8, 13)} for alpha band).
#' @param window_sec Window length in seconds for spectral estimation
#'   (default: 2).
#' @param overlap Overlap fraction between adjacent windows, from 0 to 1
#'   exclusive (default: 0.5).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{channel1}{Character or integer identifier of the first channel.}
#'     \item{channel2}{Character or integer identifier of the second channel.}
#'     \item{wpli}{Numeric wPLI value in \code{[0, 1]}.}
#'     \item{wpli_debiased}{Numeric debiased wPLI^2 value. Can be slightly
#'       negative for very weak connectivity; clamped to 0 in such cases.}
#'   }
#' @references
#' Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain signals.
#' Human Brain Mapping, 8(4), 194-208.
#'
#' Vinck, M., et al. (2011). An improved index of phase-synchronization for
#' electrophysiological data in the presence of volume-conduction, noise and
#' sample-size bias. NeuroImage, 55(4), 1548-1565.
#' @seealso [eegCoherence()], [eegPLV()], [eegConnectivityMatrix()],
#'   [eegPlotConnectivity()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' wpli_df <- eegWPLI(pe, band = c(8, 13))
#' head(wpli_df)
#' }
eegWPLI <- function(x, band = c(8, 13), window_sec = 2, overlap = 0.5,
                    assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(band) && length(band) == 2 && band[1] < band[2])
  stopifnot(is.numeric(window_sec) && window_sec > 0)
  stopifnot(is.numeric(overlap) && overlap >= 0 && overlap < 1)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 2) {
    stop("eegWPLI requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Get channel labels
  ch_labels <- seq_len(n_channels)
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% names(cd)) {
    ch_labels <- as.character(cd$label)
  }

  results <- list()

  for (i in seq_len(n_channels - 1)) {
    for (j in (i + 1):n_channels) {
      csd <- .welch_cross_spectrum(
        data[, i], data[, j], sr, window_sec, overlap
      )

      # Frequency indices within the band
      band_idx <- which(csd$freqs >= band[1] & csd$freqs <= band[2])
      if (length(band_idx) == 0) {
        stop(sprintf(
          "No frequency bins in band [%.1f, %.1f] Hz. Check band or window_sec.",
          band[1], band[2]
        ), call. = FALSE)
      }

      n_windows <- csd$n_windows

      # For each window, extract the imaginary part of the cross-spectrum
      # within the band, then compute wPLI across windows
      im_vals <- matrix(NA_real_, nrow = n_windows, ncol = length(band_idx))
      for (wi in seq_len(n_windows)) {
        sxy_seg <- csd$Sxy_segments[[wi]]
        im_vals[wi, ] <- Im(sxy_seg[band_idx])
      }

      # Average across frequency bins within the band for each window
      im_per_window <- rowMeans(im_vals)

      # wPLI = |mean(Im(Sxy))| / mean(|Im(Sxy)|)
      numerator <- abs(mean(im_per_window))
      denominator <- mean(abs(im_per_window))

      if (denominator < .Machine$double.eps) {
        wpli_val <- 0
      } else {
        wpli_val <- numerator / denominator
      }

      # Debiased wPLI^2: (N * wPLI^2 - 1) / (N - 1)
      if (n_windows > 1) {
        wpli_sq <- wpli_val^2
        wpli_debiased <- (n_windows * wpli_sq - 1) / (n_windows - 1)
        # Clamp to zero if negative (can happen for very weak connectivity)
        wpli_debiased <- max(0, wpli_debiased)
      } else {
        wpli_debiased <- NA_real_
      }

      results[[length(results) + 1]] <- data.frame(
        channel1 = ch_labels[i],
        channel2 = ch_labels[j],
        wpli = wpli_val,
        wpli_debiased = wpli_debiased,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel1 = character(0),
      channel2 = character(0),
      wpli = numeric(0),
      wpli_debiased = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}


#' EEG Spectral Granger Causality
#'
#' Computes spectral Granger causality between all EEG channel pairs by
#' fitting bivariate autoregressive (AR) models using Yule-Walker equations
#' and computing the transfer function in the frequency domain. Unlike
#' coherence-based measures, Granger causality is directional: GC from
#' channel A to channel B is generally different from GC from B to A.
#'
#' The spectral Granger causality from channel x to channel y at frequency f
#' is defined as:
#' \deqn{GC_{x \to y}(f) = \log(S_{yy}(f) / (S_{yy}(f) - |H_{xy}(f)|^2 \cdot \Sigma_{xx}))}
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param order Integer AR model order for the bivariate model (default: 5).
#' @param band Numeric vector of length 2 specifying the frequency band in Hz
#'   over which to average GC. If \code{NULL}, returns the average over all
#'   positive frequencies.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{from_channel}{Character or integer identifier of the source channel.}
#'     \item{to_channel}{Character or integer identifier of the target channel.}
#'     \item{gc_value}{Numeric Granger causality value (>= 0). Higher values
#'       indicate stronger directed influence.}
#'   }
#' @references
#' Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain signals.
#' Human Brain Mapping, 8(4), 194-208.
#'
#' Granger, C. W. J. (1969). Investigating causal relations by econometric
#' models and cross-spectral methods. Econometrica, 37(3), 424-438.
#' @seealso [eegCoherence()], [eegPLV()], [eegWPLI()],
#'   [eegConnectivityMatrix()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' gc_df <- eegGrangerCausality(pe, order = 5, band = c(8, 13))
#' head(gc_df)
#' }
eegGrangerCausality <- function(x, order = 5, band = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(order) && length(order) == 1 && order >= 1)
  order <- as.integer(order)
  if (!is.null(band)) {
    stopifnot(is.numeric(band) && length(band) == 2 && band[1] < band[2])
  }

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 2) {
    stop("eegGrangerCausality requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (n_time <= 2 * order + 1) {
    stop(sprintf(
      "Signal too short (%d samples) for AR model order %d. Need at least %d samples.",
      n_time, order, 2 * order + 2
    ), call. = FALSE)
  }

  # Get channel labels
  ch_labels <- seq_len(n_channels)
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% names(cd)) {
    ch_labels <- as.character(cd$label)
  }

  # Frequency axis for spectral GC
  n_freq <- 128L
  freqs <- seq(0, sr / 2, length.out = n_freq)

  results <- list()

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i == j) next

      # Fit bivariate AR model: [x_t; y_t] = sum A_k [x_{t-k}; y_{t-k}] + [e_t]
      # x = data[, i], y = data[, j]
      # We want GC from i -> j (does channel i help predict channel j?)
      sig_x <- data[, i]
      sig_y <- data[, j]

      # Demean signals
      sig_x <- sig_x - mean(sig_x)
      sig_y <- sig_y - mean(sig_y)

      # Construct bivariate Yule-Walker equations
      # We need the 2x2 autocorrelation matrices R(k) for k = 0, ..., order
      R <- array(0, dim = c(2, 2, order + 1))
      for (k in 0:order) {
        if (k == 0) {
          R[1, 1, k + 1] <- sum(sig_x * sig_x) / n_time
          R[1, 2, k + 1] <- sum(sig_x * sig_y) / n_time
          R[2, 1, k + 1] <- sum(sig_y * sig_x) / n_time
          R[2, 2, k + 1] <- sum(sig_y * sig_y) / n_time
        } else {
          n_eff <- n_time - k
          R[1, 1, k + 1] <- sum(sig_x[(k + 1):n_time] * sig_x[1:n_eff]) / n_time
          R[1, 2, k + 1] <- sum(sig_x[(k + 1):n_time] * sig_y[1:n_eff]) / n_time
          R[2, 1, k + 1] <- sum(sig_y[(k + 1):n_time] * sig_x[1:n_eff]) / n_time
          R[2, 2, k + 1] <- sum(sig_y[(k + 1):n_time] * sig_y[1:n_eff]) / n_time
        }
      }

      # Build the block-Toeplitz system for Yule-Walker
      # [R(0)   R(1)' ... R(p-1)'] [A(1)]   [R(1)]
      # [R(1)   R(0)  ... R(p-2)'] [A(2)] = [R(2)]
      # [  ...                    ] [ ... ]   [ ...]
      # [R(p-1) R(p-2)... R(0)   ] [A(p)]   [R(p)]

      big_R <- matrix(0, nrow = 2 * order, ncol = 2 * order)
      big_r <- numeric(2 * order)

      for (row_blk in seq_len(order)) {
        for (col_blk in seq_len(order)) {
          lag <- abs(row_blk - col_blk)
          r_idx <- (row_blk - 1) * 2 + 1:2
          c_idx <- (col_blk - 1) * 2 + 1:2
          if (row_blk >= col_blk) {
            big_R[r_idx, c_idx] <- R[, , lag + 1]
          } else {
            big_R[r_idx, c_idx] <- t(R[, , lag + 1])
          }
        }
        # Right-hand side
        r_idx <- (row_blk - 1) * 2 + 1:2
        big_r[r_idx] <- R[, , row_blk + 1][, 1]  # First column for predicting pair
      }

      # Solve for the bivariate AR coefficients
      # Use regularized solve to avoid singularity
      reg_lambda <- 1e-6 * max(abs(diag(big_R)))
      diag(big_R) <- diag(big_R) + reg_lambda

      # Build RHS for both target variables
      big_r_xy <- numeric(2 * order)
      for (row_blk in seq_len(order)) {
        r_idx <- (row_blk - 1) * 2 + 1:2
        big_r_xy[r_idx] <- R[, , row_blk + 1][, 2]  # Second column: predicting y
      }

      # Solve for AR coefficients predicting y
      A_y_vec <- tryCatch(
        solve(big_R, big_r_xy),
        error = function(e) rep(0, 2 * order)
      )

      # Also solve univariate AR for y alone (restricted model)
      # Univariate Yule-Walker for y
      R_yy <- numeric(order + 1)
      for (k in 0:order) {
        if (k == 0) {
          R_yy[k + 1] <- sum(sig_y * sig_y) / n_time
        } else {
          n_eff <- n_time - k
          R_yy[k + 1] <- sum(sig_y[(k + 1):n_time] * sig_y[1:n_eff]) / n_time
        }
      }

      toeplitz_yy <- matrix(0, nrow = order, ncol = order)
      for (row_k in seq_len(order)) {
        for (col_k in seq_len(order)) {
          toeplitz_yy[row_k, col_k] <- R_yy[abs(row_k - col_k) + 1]
        }
      }
      rhs_yy <- R_yy[2:(order + 1)]

      reg_yy <- 1e-6 * max(abs(diag(toeplitz_yy)))
      diag(toeplitz_yy) <- diag(toeplitz_yy) + reg_yy

      a_y_uni <- tryCatch(
        solve(toeplitz_yy, rhs_yy),
        error = function(e) rep(0, order)
      )

      # Compute residual variances
      # Bivariate model residual variance for y
      # Extract A_xy(k) coefficients: how x at lag k predicts y
      A_xy <- numeric(order)  # coefficient of x_{t-k} in equation for y_t
      A_yy <- numeric(order)  # coefficient of y_{t-k} in equation for y_t
      for (k in seq_len(order)) {
        # A_y_vec is organized as [a_xx(1), a_yx(1), a_xx(2), a_yx(2), ...]
        # But we solved with R columns for y, so:
        # Prediction of y_t from [x_{t-k}, y_{t-k}]
        A_xy[k] <- A_y_vec[(k - 1) * 2 + 1]  # x -> y coefficient at lag k
        A_yy[k] <- A_y_vec[(k - 1) * 2 + 2]  # y -> y coefficient at lag k
      }

      # Compute residual variances in time domain
      # Restricted (univariate) model for y
      resid_uni <- sig_y
      for (k in seq_len(order)) {
        resid_uni[(k + 1):n_time] <- resid_uni[(k + 1):n_time] -
          a_y_uni[k] * sig_y[1:(n_time - k)]
      }
      sigma_uni <- var(resid_uni[(order + 1):n_time])

      # Unrestricted (bivariate) model for y
      resid_bi <- sig_y
      for (k in seq_len(order)) {
        resid_bi[(k + 1):n_time] <- resid_bi[(k + 1):n_time] -
          A_yy[k] * sig_y[1:(n_time - k)] -
          A_xy[k] * sig_x[1:(n_time - k)]
      }
      sigma_bi <- var(resid_bi[(order + 1):n_time])

      # Spectral Granger causality
      # Compute transfer function H(f) and spectral GC
      gc_spectrum <- numeric(n_freq)

      for (fi in seq_len(n_freq)) {
        f <- freqs[fi]

        # Transfer function components
        # H_yy(f) = 1 - sum(A_yy(k) * exp(-2*pi*i*f*k/sr))
        # H_xy(f) = -sum(A_xy(k) * exp(-2*pi*i*f*k/sr))
        H_yy <- 1 + 0i
        H_xy <- 0 + 0i
        for (k in seq_len(order)) {
          exp_factor <- exp(-2i * pi * f * k / sr)
          H_yy <- H_yy - A_yy[k] * exp_factor
          H_xy <- H_xy - A_xy[k] * exp_factor
        }

        # Spectral density of y from the unrestricted model
        # Syy(f) = |H_yy|^(-2) * sigma_yy_bi + terms from x
        # Spectral GC: log(Syy_restricted / Syy_unrestricted)
        # Approximation using residual variances:
        # Syy_uni(f) ~ sigma_uni / |1 - sum(a_k * exp(-2pifk/sr))|^2
        H_uni <- 1 + 0i
        for (k in seq_len(order)) {
          H_uni <- H_uni - a_y_uni[k] * exp(-2i * pi * f * k / sr)
        }

        S_uni <- sigma_uni / Mod(H_uni)^2

        # Full model: Syy(f) = (sigma_bi * |H_yy|^2 + sigma_xx * |H_xy|^2) / |det(H)|^2
        # Simplified: use the ratio of residual-based spectral estimates
        # GC(f) = log(S_uni(f) / S_bi_y(f))
        # where S_bi_y is the spectral density of y's residuals in bivariate model
        # projected through the y-only transfer function
        S_bi <- sigma_bi / Mod(H_yy)^2

        if (S_bi > .Machine$double.eps && S_uni > .Machine$double.eps) {
          gc_f <- log(S_uni / S_bi)
          gc_spectrum[fi] <- max(0, gc_f)
        } else {
          gc_spectrum[fi] <- 0
        }
      }

      # Average GC over specified band or all frequencies
      if (!is.null(band)) {
        band_idx <- which(freqs >= band[1] & freqs <= band[2])
        if (length(band_idx) == 0) {
          gc_val <- 0
        } else {
          gc_val <- mean(gc_spectrum[band_idx])
        }
      } else {
        # Average over all positive frequencies (exclude DC)
        gc_val <- mean(gc_spectrum[-1])
      }

      gc_val <- max(0, gc_val)

      results[[length(results) + 1]] <- data.frame(
        from_channel = ch_labels[i],
        to_channel = ch_labels[j],
        gc_value = gc_val,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      from_channel = character(0),
      to_channel = character(0),
      gc_value = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}


#' EEG Connectivity Matrix
#'
#' Convenience wrapper that computes a symmetric n_channels x n_channels
#' connectivity matrix using the specified method. Dispatches to the
#' appropriate connectivity function and returns a named matrix suitable
#' for visualization or graph-theoretic analysis.
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param method Connectivity method: \code{"coherence"} for magnitude-squared
#'   coherence, \code{"plv"} for Phase Locking Value, or \code{"wpli"} for
#'   weighted Phase Lag Index (default: \code{"coherence"}).
#' @param band Numeric vector of length 2 specifying the frequency band in Hz
#'   (default: \code{c(8, 13)} for alpha band).
#' @param ... Additional arguments passed to the underlying connectivity
#'   function (e.g., \code{window_sec}, \code{overlap}).
#' @return A named numeric matrix of dimension n_channels x n_channels with
#'   channel labels as row and column names. Diagonal elements are 1.
#'   Off-diagonal elements represent the band-averaged connectivity
#'   between channel pairs.
#' @references
#' Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain signals.
#' Human Brain Mapping, 8(4), 194-208.
#' @seealso [eegCoherence()], [eegPLV()], [eegWPLI()],
#'   [eegGrangerCausality()], [eegPlotConnectivity()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' conn <- eegConnectivityMatrix(pe, method = "plv", band = c(8, 13))
#' print(conn)
#' }
eegConnectivityMatrix <- function(x, method = c("coherence", "plv", "wpli"),
                                  band = c(8, 13), ...) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(band) && length(band) == 2 && band[1] < band[2])

  if (is.null(x)) {
    stop("Input PhysioExperiment object must not be NULL.", call. = FALSE)
  }

  # Determine number of channels and labels
  assay_name <- NULL
  dots <- list(...)
  if ("assay_name" %in% names(dots)) {
    assay_name <- dots$assay_name
  }
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_channels <- ncol(data)

  ch_labels <- seq_len(n_channels)
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% names(cd)) {
    ch_labels <- as.character(cd$label)
  }

  # Initialize result matrix
  conn_matrix <- matrix(0, nrow = n_channels, ncol = n_channels)
  diag(conn_matrix) <- 1
  rownames(conn_matrix) <- ch_labels
  colnames(conn_matrix) <- ch_labels

  if (method == "coherence") {
    x_result <- eegCoherence(x, method = "coherence", band = band, ...)
    conn_matrix <- S4Vectors::metadata(x_result)$connectivity$matrix
  } else if (method == "plv") {
    plv_df <- eegPLV(x, band = band, assay_name = dots$assay_name)
    # Fill matrix from data.frame
    for (row_i in seq_len(nrow(plv_df))) {
      ch1 <- plv_df$channel1[row_i]
      ch2 <- plv_df$channel2[row_i]
      val <- plv_df$plv[row_i]

      idx1 <- which(ch_labels == ch1)
      idx2 <- which(ch_labels == ch2)
      if (length(idx1) == 1 && length(idx2) == 1) {
        conn_matrix[idx1, idx2] <- val
        conn_matrix[idx2, idx1] <- val
      }
    }
  } else if (method == "wpli") {
    wpli_df <- eegWPLI(x, band = band,
                       window_sec = if ("window_sec" %in% names(dots)) dots$window_sec else 2,
                       overlap = if ("overlap" %in% names(dots)) dots$overlap else 0.5,
                       assay_name = dots$assay_name)
    # Fill matrix from data.frame
    for (row_i in seq_len(nrow(wpli_df))) {
      ch1 <- wpli_df$channel1[row_i]
      ch2 <- wpli_df$channel2[row_i]
      val <- wpli_df$wpli[row_i]

      idx1 <- which(ch_labels == ch1)
      idx2 <- which(ch_labels == ch2)
      if (length(idx1) == 1 && length(idx2) == 1) {
        conn_matrix[idx1, idx2] <- val
        conn_matrix[idx2, idx1] <- val
      }
    }
  }

  conn_matrix
}
