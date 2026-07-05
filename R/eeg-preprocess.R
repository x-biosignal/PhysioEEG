###############################################################################
# EEG Preprocessing Module
#
# Provides a comprehensive EEG preprocessing pipeline including filtering,
# re-referencing, bad channel detection/interpolation, epoching, artifact
# rejection, montage assignment, and a convenience wrapper.
#
# Internal FIR helpers .fir_bandpass() and .fir_lowpass() are defined in
# eeg-sleep.R and available package-wide. This file adds .fir_highpass(),
# .fir_notch(), and other internal helpers.
###############################################################################


# ============================================================================
# Internal Helpers
# ============================================================================

#' Windowed-sinc FIR highpass filter
#'
#' Applies a windowed-sinc FIR highpass filter using spectral inversion of
#' a lowpass filter. The lowpass kernel at the cutoff frequency is computed
#' and then inverted to create a highpass response.
#'
#' @param signal Numeric vector of the input signal.
#' @param sr Sampling rate in Hz.
#' @param cutoff Cutoff frequency in Hz.
#' @param order Filter order. If \code{NULL}, auto-selected.
#' @return Numeric vector of the filtered signal (same length as input).
#' @keywords internal
.fir_highpass <- function(signal, sr, cutoff, order = NULL) {
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

  # Build lowpass kernel

  fc <- cutoff / sr
  h_lp <- ifelse(m == 0, 2 * fc, sin(2 * pi * fc * m) / (pi * m))

  # Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * (0:(order - 1)) / (order - 1))
  h_lp <- h_lp * w
  h_lp <- h_lp / sum(h_lp)  # normalise to unity DC gain


  # Spectral inversion: negate then add 1 at center
  h <- -h_lp
  h[half + 1] <- h[half + 1] + 1

  # Apply via causal convolution then compensate group delay
  pad_len <- order
  sig_padded <- c(rep(signal[1], pad_len), signal, rep(signal[n], pad_len))
  filtered <- stats::filter(sig_padded, h, sides = 1)
  filtered <- as.numeric(filtered)
  result <- filtered[(pad_len + half + 1):(pad_len + half + n)]

  na_idx <- which(is.na(result))
  if (length(na_idx) > 0) result[na_idx] <- signal[na_idx]
  result
}


#' Windowed-sinc FIR notch (band-stop) filter
#'
#' Implements a band-stop filter as allpass minus bandpass. The bandpass kernel
#' at the specified center frequency +/- bandwidth/2 is subtracted from a
#' unit impulse to produce a notch response.
#'
#' @param signal Numeric vector of the input signal.
#' @param sr Sampling rate in Hz.
#' @param center_freq Center frequency to reject in Hz.
#' @param bandwidth Full bandwidth of the notch in Hz (default: 4).
#' @param order Filter order. If \code{NULL}, auto-selected.
#' @return Numeric vector of the filtered signal (same length as input).
#' @keywords internal
.fir_notch <- function(signal, sr, center_freq, bandwidth = 4, order = NULL) {
  n <- length(signal)
  low <- center_freq - bandwidth / 2
  high <- center_freq + bandwidth / 2

  if (low <= 0 || high >= sr / 2) {
    warning("Notch frequency range exceeds valid bounds. Returning original signal.")
    return(signal)
  }

  if (is.null(order)) {
    order <- as.integer(round(3 * sr / max(low, 0.5)))
    order <- min(order, n - 1L)
    if (order %% 2 == 0) order <- order + 1L
  }
  order <- min(order, n - 1L)
  if (order < 3L) return(signal)
  if (order %% 2 == 0) order <- order + 1L

  half <- (order - 1L) %/% 2L
  m <- seq(-half, half)

  # Build bandpass kernel
  fc_low <- low / sr
  fc_high <- high / sr
  h_high <- ifelse(m == 0, 2 * fc_high, sin(2 * pi * fc_high * m) / (pi * m))
  h_low  <- ifelse(m == 0, 2 * fc_low,  sin(2 * pi * fc_low * m) / (pi * m))
  h_bp <- h_high - h_low

  # Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * (0:(order - 1)) / (order - 1))
  h_bp <- h_bp * w

  # Normalise bandpass to unit gain at centre
  centre_response <- sum(h_bp * cos(2 * pi * center_freq / sr * m))
  if (abs(centre_response) > 1e-10) h_bp <- h_bp / centre_response


  # Band-stop = allpass - bandpass (spectral inversion)
  h <- -h_bp
  h[half + 1] <- h[half + 1] + 1

  # Apply filter
  pad_len <- order
  sig_padded <- c(rep(signal[1], pad_len), signal, rep(signal[n], pad_len))
  filtered <- stats::filter(sig_padded, h, sides = 1)
  filtered <- as.numeric(filtered)
  result <- filtered[(pad_len + half + 1):(pad_len + half + n)]

  na_idx <- which(is.na(result))
  if (length(na_idx) > 0) result[na_idx] <- signal[na_idx]
  result
}


#' IIR Butterworth filter wrapper
#'
#' Applies a Butterworth IIR filter using \code{signal::butter} and
#' \code{signal::filtfilt} for zero-phase filtering. Falls back to FIR
#' if the signal package is not available.
#'
#' @param data_matrix Numeric matrix (time x channels).
#' @param sr Sampling rate in Hz.
#' @param lowcut Low cutoff frequency in Hz (NULL for lowpass).
#' @param highcut High cutoff frequency in Hz (NULL for highpass).
#' @param order Filter order (default: 4).
#' @return Filtered numeric matrix (same dimensions).
#' @keywords internal
.iir_butterworth <- function(data_matrix, sr, lowcut = NULL, highcut = NULL,
                             order = 4) {
  has_signal <- requireNamespace("signal", quietly = TRUE)

  if (!has_signal) {
    warning("Package 'signal' not available. Falling back to FIR filtering.")
    result <- data_matrix
    for (ch in seq_len(ncol(data_matrix))) {
      if (!is.null(lowcut) && !is.null(highcut)) {
        result[, ch] <- .fir_bandpass(data_matrix[, ch], sr, lowcut, highcut)
      } else if (!is.null(lowcut)) {
        result[, ch] <- .fir_highpass(data_matrix[, ch], sr, lowcut)
      } else if (!is.null(highcut)) {
        result[, ch] <- .fir_lowpass(data_matrix[, ch], sr, highcut)
      }
    }
    return(result)
  }

  nyq <- sr / 2

  # Design Butterworth filter
  bf <- tryCatch({
    if (!is.null(lowcut) && !is.null(highcut)) {
      W <- c(lowcut / nyq, highcut / nyq)
      signal::butter(order, W, type = "pass")
    } else if (!is.null(lowcut)) {
      signal::butter(order, lowcut / nyq, type = "high")
    } else if (!is.null(highcut)) {
      signal::butter(order, highcut / nyq, type = "low")
    } else {
      stop("At least one of lowcut or highcut must be specified.")
    }
  }, error = function(e) {
    warning("IIR filter design failed: ", conditionMessage(e),
            ". Falling back to FIR.")
    return(NULL)
  })

  if (is.null(bf)) {
    # Fallback to FIR
    result <- data_matrix
    for (ch in seq_len(ncol(data_matrix))) {
      if (!is.null(lowcut) && !is.null(highcut)) {
        result[, ch] <- .fir_bandpass(data_matrix[, ch], sr, lowcut, highcut)
      } else if (!is.null(lowcut)) {
        result[, ch] <- .fir_highpass(data_matrix[, ch], sr, lowcut)
      } else if (!is.null(highcut)) {
        result[, ch] <- .fir_lowpass(data_matrix[, ch], sr, highcut)
      }
    }
    return(result)
  }

  # Apply filtfilt to each channel
  result <- data_matrix
  for (ch in seq_len(ncol(data_matrix))) {
    result[, ch] <- tryCatch(
      signal::filtfilt(bf, data_matrix[, ch]),
      error = function(e) {
        warning("filtfilt failed for channel ", ch, ": ", conditionMessage(e),
                ". Falling back to FIR for this channel.")
        if (!is.null(lowcut) && !is.null(highcut)) {
          .fir_bandpass(data_matrix[, ch], sr, lowcut, highcut)
        } else if (!is.null(lowcut)) {
          .fir_highpass(data_matrix[, ch], sr, lowcut)
        } else {
          .fir_lowpass(data_matrix[, ch], sr, highcut)
        }
      }
    )
  }

  result
}


#' Compute Legendre polynomial of degree n at x
#'
#' Uses the Bonnet recursion formula to compute the Legendre polynomial
#' P_n(x) for a given degree n and evaluation point(s) x.
#'
#' @param n Non-negative integer degree.
#' @param x Numeric vector of evaluation points in \code{[-1, 1]}.
#' @return Numeric vector of P_n(x).
#' @keywords internal
.legendre_poly <- function(n, x) {
  if (n == 0) return(rep(1, length(x)))
  if (n == 1) return(x)

  p_prev <- rep(1, length(x))
  p_curr <- x

  for (k in 2:n) {
    p_next <- ((2 * k - 1) * x * p_curr - (k - 1) * p_prev) / k
    p_prev <- p_curr
    p_curr <- p_next
  }

  p_curr
}


#' Compute spherical spline interpolation weights
#'
#' Implements the spherical spline method of Perrin et al. (1989) for
#' interpolating EEG channel data. Uses Legendre polynomials to compute
#' the interpolation matrix G and the weight vector g_x for each bad channel.
#'
#' @param good_pos Matrix with columns (pos_x, pos_y, pos_z) for good channels.
#' @param bad_pos Matrix with columns (pos_x, pos_y, pos_z) for bad channels.
#' @param order Maximum Legendre polynomial order (default: 4).
#' @return Matrix of weights (n_bad x n_good) for interpolation.
#' @keywords internal
.spherical_spline_weights <- function(good_pos, bad_pos, order = 4) {
  n_good <- nrow(good_pos)
  n_bad <- nrow(bad_pos)

  # Normalise positions to unit sphere
  .normalise_rows <- function(mat) {
    norms <- sqrt(rowSums(mat^2))
    norms[norms == 0] <- 1
    mat / norms
  }

  good_norm <- .normalise_rows(good_pos)
  bad_norm <- .normalise_rows(bad_pos)

  # Compute G matrix (good x good) using cosine distances
  # g(cos_theta) = sum_{n=1}^{order} (2n+1) / (n*(n+1))^2 * P_n(cos_theta)
  .g_func <- function(cos_theta) {
    val <- rep(0, length(cos_theta))
    for (k in 1:order) {
      coeff <- (2 * k + 1) / (k * (k + 1))^2
      val <- val + coeff * .legendre_poly(k, cos_theta)
    }
    val
  }

  # Build G matrix for good channels
  G <- matrix(0, nrow = n_good, ncol = n_good)
  for (i in seq_len(n_good)) {
    for (j in seq_len(n_good)) {
      cos_ij <- sum(good_norm[i, ] * good_norm[j, ])
      cos_ij <- max(min(cos_ij, 1), -1)  # clamp
      G[i, j] <- .g_func(cos_ij)
    }
  }

  # Regularise: add small ridge to diagonal for numerical stability
  diag(G) <- diag(G) + 1e-5

  # Augmented system for spline with constraint sum(c_i) = 0
  # [G  1] [c]   [g_x]
  # [1' 0] [lam] = [0 ]
  G_aug <- rbind(
    cbind(G, rep(1, n_good)),
    c(rep(1, n_good), 0)
  )

  # Compute weights for each bad channel
  weights <- matrix(0, nrow = n_bad, ncol = n_good)

  for (b in seq_len(n_bad)) {
    g_x <- numeric(n_good)
    for (i in seq_len(n_good)) {
      cos_bi <- sum(bad_norm[b, ] * good_norm[i, ])
      cos_bi <- max(min(cos_bi, 1), -1)
      g_x[i] <- .g_func(cos_bi)
    }

    rhs <- c(g_x, 0)

    # Solve augmented system
    sol <- tryCatch(
      solve(G_aug, rhs),
      error = function(e) {
        # Fallback: use pseudo-inverse
        pinv <- tryCatch(
          solve(t(G_aug) %*% G_aug + diag(1e-8, nrow(G_aug))) %*% t(G_aug),
          error = function(e2) {
            # Ultimate fallback: equal weights
            return(NULL)
          }
        )
        if (is.null(pinv)) return(c(rep(1 / n_good, n_good), 0))
        pinv %*% rhs
      }
    )

    weights[b, ] <- sol[seq_len(n_good)]
  }

  weights
}


#' Standard 10-20 electrode positions
#'
#' Returns a data.frame with label and (pos_x, pos_y, pos_z) for the
#' 19 standard 10-20 system electrodes. Coordinates are based on the
#' standard spherical projection to 2D (azimuthal equidistant) with
#' z computed to lie on the unit sphere.
#'
#' @return data.frame with columns: label, pos_x, pos_y, pos_z.
#' @keywords internal
.electrode_positions_1020 <- function() {
  # 10-20 system positions in standard coordinates

  # Using spherical coordinates (theta, phi) projected to Cartesian
  # theta: angle from midline (positive = right), phi: angle from vertex
  # Standard layout on unit sphere, projected as:

  #   x = sin(phi) * sin(theta)  (left-right)
  #   y = sin(phi) * cos(theta)  (posterior-anterior)
  #   z = cos(phi)                (inferior-superior)
  #
  # Convention: nose at +y, right ear at +x, vertex at +z

  data.frame(
    label = c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
              "T3", "C3", "Cz", "C4", "T4",
              "T5", "P3", "Pz", "P4", "T6",
              "O1", "O2"),
    pos_x = c(-0.309, 0.309, -0.809, -0.545, 0.000, 0.545, 0.809,
              -1.000, -0.707, 0.000, 0.707, 1.000,
              -0.809, -0.545, 0.000, 0.545, 0.809,
              -0.309, 0.309),
    pos_y = c(0.951, 0.951, 0.588, 0.809, 0.891, 0.809, 0.588,
              0.000, 0.000, 0.000, 0.000, 0.000,
              -0.588, -0.809, -0.891, -0.809, -0.588,
              -0.951, -0.951),
    pos_z = c(0.000, 0.000, 0.000, 0.225, 0.454, 0.225, 0.000,
              0.000, 0.707, 1.000, 0.707, 0.000,
              0.000, 0.225, 0.454, 0.225, 0.000,
              0.000, 0.000),
    stringsAsFactors = FALSE
  )
}


#' Extended 10-10 electrode positions
#'
#' Returns a data.frame with label and (pos_x, pos_y, pos_z) for
#' approximately 64 channels in the extended 10-10 system, including
#' all 10-20 positions plus intermediate sites.
#'
#' @return data.frame with columns: label, pos_x, pos_y, pos_z.
#' @keywords internal
.electrode_positions_1010 <- function() {
  # Start with 10-20 positions
  pos_1020 <- .electrode_positions_1020()


  # Additional 10-10 positions (intermediate sites)
  extra <- data.frame(
    label = c(
      # Frontal-polar
      "Fpz",
      # Frontal row
      "AF7", "AF3", "AFz", "AF4", "AF8",
      # Frontal-central
      "F1", "F2", "F5", "F6",
      # Fronto-central
      "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6",
      # Central
      "C1", "C2", "C5", "C6",
      # Centro-parietal
      "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6",
      # Temporal-parietal
      "TP7", "TP8",
      # Parietal
      "P1", "P2", "P5", "P6", "P7", "P8",
      # Parietal-occipital
      "PO7", "PO3", "POz", "PO4", "PO8",
      # Occipital
      "Oz",
      # Mastoids / ear references
      "M1", "M2", "TP9", "TP10"
    ),
    pos_x = c(
      # Fpz
      0.000,
      # AF7, AF3, AFz, AF4, AF8
      -0.587, -0.309, 0.000, 0.309, 0.587,
      # F1, F2, F5, F6
      -0.272, 0.272, -0.677, 0.677,
      # FC5, FC3, FC1, FCz, FC2, FC4, FC6
      -0.905, -0.626, -0.354, 0.000, 0.354, 0.626, 0.905,
      # C1, C2, C5, C6
      -0.354, 0.354, -0.905, 0.905,
      # CP5, CP3, CP1, CPz, CP2, CP4, CP6
      -0.905, -0.626, -0.354, 0.000, 0.354, 0.626, 0.905,
      # TP7, TP8
      -0.987, 0.987,
      # P1, P2, P5, P6, P7, P8
      -0.272, 0.272, -0.677, 0.677, -0.809, 0.809,
      # PO7, PO3, POz, PO4, PO8
      -0.587, -0.309, 0.000, 0.309, 0.587,
      # Oz
      0.000,
      # M1, M2, TP9, TP10
      -1.050, 1.050, -1.050, 1.050
    ),
    pos_y = c(
      # Fpz
      0.999,
      # AF7, AF3, AFz, AF4, AF8
      0.809, 0.951, 0.950, 0.951, 0.809,
      # F1, F2, F5, F6
      0.850, 0.850, 0.700, 0.700,
      # FC5, FC3, FC1, FCz, FC2, FC4, FC6
      0.294, 0.405, 0.446, 0.450, 0.446, 0.405, 0.294,
      # C1, C2, C5, C6
      0.000, 0.000, 0.000, 0.000,
      # CP5, CP3, CP1, CPz, CP2, CP4, CP6
      -0.294, -0.405, -0.446, -0.450, -0.446, -0.405, -0.294,
      # TP7, TP8
      -0.156, -0.156,
      # P1, P2, P5, P6, P7, P8
      -0.850, -0.850, -0.700, -0.700, -0.588, -0.588,
      # PO7, PO3, POz, PO4, PO8
      -0.809, -0.951, -0.950, -0.951, -0.809,
      # Oz
      -0.999,
      # M1, M2, TP9, TP10
      -0.156, -0.156, -0.250, -0.250
    ),
    pos_z = c(
      # Fpz
      0.000,
      # AF7, AF3, AFz, AF4, AF8
      0.000, 0.050, 0.100, 0.050, 0.000,
      # F1, F2, F5, F6
      0.340, 0.340, 0.113, 0.113,
      # FC5, FC3, FC1, FCz, FC2, FC4, FC6
      0.225, 0.550, 0.750, 0.800, 0.750, 0.550, 0.225,
      # C1, C2, C5, C6
      0.866, 0.866, 0.354, 0.354,
      # CP5, CP3, CP1, CPz, CP2, CP4, CP6
      0.225, 0.550, 0.750, 0.800, 0.750, 0.550, 0.225,
      # TP7, TP8
      0.000, 0.000,
      # P1, P2, P5, P6, P7, P8
      0.340, 0.340, 0.113, 0.113, 0.000, 0.000,
      # PO7, PO3, POz, PO4, PO8
      0.000, 0.050, 0.100, 0.050, 0.000,
      # Oz
      0.000,
      # M1, M2, TP9, TP10
      -0.100, -0.100, -0.200, -0.200
    ),
    stringsAsFactors = FALSE
  )

  rbind(pos_1020, extra)
}


# ============================================================================
# Exported Functions
# ============================================================================

#' Filter EEG signals
#'
#' Applies frequency-domain filtering to EEG data stored in a PhysioExperiment
#' object. Supports bandpass, highpass, lowpass, and notch filtering using
#' either FIR (windowed-sinc) or IIR (Butterworth) methods.
#'
#' For FIR mode, the function uses windowed-sinc filters with a Hamming window.
#' For IIR mode, zero-phase Butterworth filtering is applied via
#' \code{signal::filtfilt()}, with automatic fallback to FIR if the signal
#' package is not available.
#'
#' @param x A PhysioExperiment object.
#' @param lowcut Low cutoff frequency in Hz. If NULL, no highpass filtering
#'   is applied (set both lowcut and highcut for bandpass).
#' @param highcut High cutoff frequency in Hz. If NULL, no lowpass filtering
#'   is applied (set both lowcut and highcut for bandpass).
#' @param notch Notch filter center frequency in Hz (e.g., 50 or 60 for
#'   powerline noise). If NULL, no notch filter is applied. Bandwidth is
#'   +/- 2 Hz around center.
#' @param method Filtering method: \code{"fir"} (default) for windowed-sinc
#'   FIR filter, or \code{"iir"} for Butterworth IIR filter.
#' @param order Filter order. For FIR, this is the number of taps (auto-selected
#'   if NULL). For IIR, this is the Butterworth order (default: 4).
#' @param assay_name Name of the assay to filter. If NULL, uses \code{defaultAssay(x)}.
#' @param output_assay Name of the output assay (default: \code{"filtered"}).
#'
#' @return A PhysioExperiment object with filtered data in the specified output assay.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' # Bandpass filter
#' pe_filt <- eegFilter(pe, lowcut = 1, highcut = 40)
#' # Highpass only
#' pe_hp <- eegFilter(pe, lowcut = 0.1)
#' # With notch at 50 Hz
#' pe_notch <- eegFilter(pe, lowcut = 1, highcut = 40, notch = 50)
#' # IIR Butterworth
#' pe_iir <- eegFilter(pe, lowcut = 1, highcut = 40, method = "iir", order = 4)
#' }
eegFilter <- function(x, lowcut = NULL, highcut = NULL, notch = NULL,
                       method = c("fir", "iir"), order = NULL,
                       assay_name = NULL, output_assay = "filtered") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(lowcut) && is.null(highcut) && is.null(notch)) {
    stop("At least one of lowcut, highcut, or notch must be specified.")
  }

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  sr <- samplingRate(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  if (!is.matrix(data)) {
    stop("eegFilter requires a 2D matrix assay. Epoch data first if needed.")
  }

  nyq <- sr / 2

  # Validate frequency parameters
  if (!is.null(lowcut)) {
    stopifnot(is.numeric(lowcut), lowcut > 0, lowcut < nyq)
  }
  if (!is.null(highcut)) {
    stopifnot(is.numeric(highcut), highcut > 0, highcut < nyq)
  }
  if (!is.null(lowcut) && !is.null(highcut)) {
    if (lowcut >= highcut) {
      stop("lowcut must be less than highcut.")
    }
  }
  if (!is.null(notch)) {
    stopifnot(is.numeric(notch), notch > 2, notch < nyq - 2)
  }

  result <- data

  # Apply main filter (bandpass / highpass / lowpass)
  if (!is.null(lowcut) || !is.null(highcut)) {
    if (method == "iir") {
      iir_order <- if (is.null(order)) 4L else as.integer(order)
      result <- .iir_butterworth(result, sr, lowcut, highcut, order = iir_order)
    } else {
      # FIR filtering
      for (ch in seq_len(ncol(result))) {
        if (!is.null(lowcut) && !is.null(highcut)) {
          result[, ch] <- .fir_bandpass(result[, ch], sr, lowcut, highcut,
                                         order = order)
        } else if (!is.null(lowcut)) {
          result[, ch] <- .fir_highpass(result[, ch], sr, lowcut,
                                         order = order)
        } else {
          result[, ch] <- .fir_lowpass(result[, ch], sr, highcut,
                                        order = order)
        }
      }
    }
  }

  # Apply notch filter if specified
  if (!is.null(notch)) {
    for (ch in seq_len(ncol(result))) {
      result[, ch] <- .fir_notch(result[, ch], sr, notch, bandwidth = 4,
                                   order = order)
    }
  }

  SummarizedExperiment::assay(x, output_assay) <- result
  x
}


#' Re-reference EEG data
#'
#' Applies a re-referencing scheme to EEG data. Re-referencing transforms
#' the data by subtracting a reference signal from all channels, which can
#' improve spatial resolution and comparability across studies.
#'
#' @param x A PhysioExperiment object.
#' @param ref_type Re-referencing scheme: \code{"average"} (common average),
#'   \code{"mastoids"} (linked mastoids), \code{"cz"} (Cz reference),
#'   \code{"rest"} (Reference Electrode Standardization Technique),
#'   or \code{"channel"} (user-specified channels).
#' @param ref_channels Character vector of channel labels to use as reference
#'   (required for \code{ref_type = "channel"}).
#' @param exclude Character vector of channel labels to exclude from the
#'   average reference calculation (only used for \code{ref_type = "average"}).
#' @param assay_name Name of the assay to re-reference. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name of the output assay (default: \code{"rereferenced"}).
#'
#' @return A PhysioExperiment object with re-referenced data in the specified
#'   output assay. Stores reference info in \code{metadata(x)$reference}.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe_avg <- eegRereference(pe, ref_type = "average")
#' pe_cz  <- eegRereference(pe, ref_type = "cz")
#' }
eegRereference <- function(x, ref_type = c("average", "mastoids", "cz",
                                            "rest", "channel"),
                            ref_channels = NULL, exclude = NULL,
                            assay_name = NULL,
                            output_assay = "rereferenced") {
  stopifnot(inherits(x, "PhysioExperiment"))
  ref_type <- match.arg(ref_type)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  if (!is.matrix(data)) {
    stop("eegRereference requires a 2D matrix assay.")
  }

  cd <- SummarizedExperiment::colData(x)
  labels <- as.character(cd$label)
  n_ch <- ncol(data)

  result <- data

  if (ref_type == "average") {
    # Common average reference, optionally excluding some channels
    if (!is.null(exclude)) {
      excl_idx <- which(labels %in% exclude)
      if (length(excl_idx) == 0) {
        warning("No channels matched 'exclude'. Using all channels.")
        avg <- rowMeans(data)
      } else {
        incl_idx <- setdiff(seq_len(n_ch), excl_idx)
        if (length(incl_idx) == 0) {
          stop("All channels excluded. Cannot compute average reference.")
        }
        avg <- rowMeans(data[, incl_idx, drop = FALSE])
      }
    } else {
      avg <- rowMeans(data)
    }
    for (ch in seq_len(n_ch)) {
      result[, ch] <- data[, ch] - avg
    }

  } else if (ref_type == "channel") {
    if (is.null(ref_channels) || length(ref_channels) == 0) {
      stop("ref_channels must be specified for ref_type = 'channel'.")
    }
    ref_idx <- which(labels %in% ref_channels)
    if (length(ref_idx) == 0) {
      stop("No matching channels found for ref_channels: ",
           paste(ref_channels, collapse = ", "))
    }
    ref_signal <- rowMeans(data[, ref_idx, drop = FALSE])
    for (ch in seq_len(n_ch)) {
      result[, ch] <- data[, ch] - ref_signal
    }

  } else if (ref_type == "mastoids") {
    # Search for mastoid/ear channels
    mastoid_patterns <- list(
      c("M1", "M2"),
      c("A1", "A2"),
      c("TP9", "TP10")
    )
    found <- FALSE
    for (pat in mastoid_patterns) {
      idx <- which(labels %in% pat)
      if (length(idx) == 2) {
        ref_signal <- rowMeans(data[, idx, drop = FALSE])
        for (ch in seq_len(n_ch)) {
          result[, ch] <- data[, ch] - ref_signal
        }
        found <- TRUE
        break
      }
    }
    if (!found) {
      stop("Could not find mastoid/ear channels. ",
           "Expected M1/M2, A1/A2, or TP9/TP10 in colData$label.")
    }

  } else if (ref_type == "cz") {
    cz_idx <- which(labels == "Cz")
    if (length(cz_idx) == 0) {
      stop("Cz channel not found in colData$label.")
    }
    ref_signal <- data[, cz_idx[1]]
    for (ch in seq_len(n_ch)) {
      result[, ch] <- data[, ch] - ref_signal
    }

  } else if (ref_type == "rest") {
    # Reference Electrode Standardization Technique (simplified)
    # Uses inverse distance weighting from electrode positions
    pos_cols <- c("pos_x", "pos_y", "pos_z")
    if (!all(pos_cols %in% colnames(cd))) {
      stop("REST re-referencing requires electrode positions (pos_x, pos_y, pos_z) ",
           "in colData. Apply eegMontage() first.")
    }

    pos <- as.matrix(cd[, pos_cols])

    # Compute pairwise distances
    n <- nrow(pos)
    dist_mat <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        dist_mat[i, j] <- sqrt(sum((pos[i, ] - pos[j, ])^2))
      }
    }

    # Inverse distance weights for each channel
    # Each channel's reference is weighted average of all others
    weights <- matrix(0, n, n)
    for (i in seq_len(n)) {
      dists <- dist_mat[i, ]
      dists[i] <- Inf  # exclude self
      inv_dists <- 1 / dists
      weights[i, ] <- inv_dists / sum(inv_dists)
    }

    # Weighted average reference per channel
    for (ch in seq_len(n_ch)) {
      weighted_ref <- data %*% weights[ch, ]
      result[, ch] <- data[, ch] - weighted_ref
    }
  }

  # Store reference info
  ref_info <- list(
    type = ref_type,
    ref_channels = ref_channels,
    exclude = exclude,
    timestamp = Sys.time()
  )
  md <- S4Vectors::metadata(x)
  md$reference <- ref_info
  S4Vectors::metadata(x) <- md

  SummarizedExperiment::assay(x, output_assay) <- result
  x
}


#' Detect bad EEG channels
#'
#' Identifies bad (noisy, flat, or poorly correlated) EEG channels using
#' multiple automated criteria. Channels flagged as bad can subsequently
#' be interpolated using \code{\link{eegInterpolate}}.
#'
#' @param x A PhysioExperiment object.
#' @param method Detection method(s) to apply: \code{"all"} runs all checks,
#'   or specify one or more of \code{"flat"}, \code{"noise"}, \code{"correlation"}.
#' @param flat_threshold Variance threshold below which a channel is considered
#'   flat (default: 1e-6).
#' @param noise_threshold Number of standard deviations above median variance
#'   to flag a channel as noisy (default: 4).
#' @param corr_threshold Minimum mean correlation with other channels. Channels
#'   below this are flagged (default: 0.4).
#' @param assay_name Name of the assay to analyze. If NULL, uses
#'   \code{defaultAssay(x)}.
#'
#' @return A data.frame with columns: \code{channel} (label), \code{is_bad}
#'   (logical), \code{reason} (character description), \code{score} (numeric
#'   metric value).
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' bad_df <- eegBadChannels(pe)
#' bad_labels <- bad_df$channel[bad_df$is_bad]
#' }
eegBadChannels <- function(x, method = c("all", "flat", "noise", "correlation"),
                            flat_threshold = 1e-6, noise_threshold = 4,
                            corr_threshold = 0.4, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method, several.ok = TRUE)

  if ("all" %in% method) {
    method <- c("flat", "noise", "correlation")
  }

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  if (!is.matrix(data)) {
    stop("eegBadChannels requires a 2D matrix assay.")
  }

  cd <- SummarizedExperiment::colData(x)
  labels <- as.character(cd$label)
  n_ch <- ncol(data)

  # Initialize results
  results <- data.frame(
    channel = labels,
    is_bad = rep(FALSE, n_ch),
    reason = rep("", n_ch),
    score = rep(NA_real_, n_ch),
    stringsAsFactors = FALSE
  )

  # Compute channel variances (needed for flat and noise)
  ch_vars <- apply(data, 2, stats::var)

  # -- Flat channel detection --
  if ("flat" %in% method) {
    flat_idx <- which(ch_vars < flat_threshold)
    if (length(flat_idx) > 0) {
      results$is_bad[flat_idx] <- TRUE
      results$reason[flat_idx] <- ifelse(
        results$reason[flat_idx] == "",
        "flat",
        paste(results$reason[flat_idx], "flat", sep = "; ")
      )
      results$score[flat_idx] <- ch_vars[flat_idx]
    }
  }

  # -- Noise detection (z-scored variance) --
  if ("noise" %in% method) {
    median_var <- stats::median(ch_vars)
    mad_var <- stats::mad(ch_vars, constant = 1)
    if (mad_var < 1e-12) mad_var <- 1e-12  # avoid division by zero
    z_vars <- (ch_vars - median_var) / mad_var

    noisy_idx <- which(z_vars > noise_threshold)
    if (length(noisy_idx) > 0) {
      results$is_bad[noisy_idx] <- TRUE
      results$reason[noisy_idx] <- ifelse(
        results$reason[noisy_idx] == "",
        "noise",
        paste(results$reason[noisy_idx], "noise", sep = "; ")
      )
      # Prefer to store the noise z-score if score not yet set
      for (idx in noisy_idx) {
        if (is.na(results$score[idx])) {
          results$score[idx] <- z_vars[idx]
        }
      }
    }
  }

  # -- Correlation detection --
  if ("correlation" %in% method) {
    corr_mat <- stats::cor(data)
    # Set diagonal to NA so it doesn't affect the mean
    diag(corr_mat) <- NA
    mean_corr <- colMeans(corr_mat, na.rm = TRUE)

    low_corr_idx <- which(mean_corr < corr_threshold)
    if (length(low_corr_idx) > 0) {
      results$is_bad[low_corr_idx] <- TRUE
      results$reason[low_corr_idx] <- ifelse(
        results$reason[low_corr_idx] == "",
        "correlation",
        paste(results$reason[low_corr_idx], "correlation", sep = "; ")
      )
      for (idx in low_corr_idx) {
        if (is.na(results$score[idx])) {
          results$score[idx] <- mean_corr[idx]
        }
      }
    }
  }

  results
}


#' Interpolate bad EEG channels
#'
#' Replaces data in bad channels by interpolating from remaining good
#' channels using either spherical spline interpolation (Perrin et al., 1989)
#' or nearest-neighbor weighted averaging.
#'
#' Requires electrode positions (pos_x, pos_y, pos_z) in colData.
#' Apply \code{\link{eegMontage}} first if positions are not set.
#'
#' @param x A PhysioExperiment object.
#' @param bad_channels Character vector of channel labels to interpolate.
#' @param method Interpolation method: \code{"spline"} for spherical spline
#'   (default) or \code{"nearest"} for inverse-distance weighted nearest
#'   neighbors.
#' @param assay_name Name of the assay to interpolate. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name of the output assay (default: \code{"interpolated"}).
#'
#' @return A PhysioExperiment object with interpolated channels in the
#'   specified output assay.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe <- eegMontage(pe, system = "10-20")
#' bad_df <- eegBadChannels(pe)
#' bad_labels <- bad_df$channel[bad_df$is_bad]
#' if (length(bad_labels) > 0) {
#'   pe_clean <- eegInterpolate(pe, bad_labels)
#' }
#' }
eegInterpolate <- function(x, bad_channels, method = c("spline", "nearest"),
                            assay_name = NULL, output_assay = "interpolated") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  if (!is.matrix(data)) {
    stop("eegInterpolate requires a 2D matrix assay.")
  }

  cd <- SummarizedExperiment::colData(x)
  labels <- as.character(cd$label)

  # Check for electrode positions
  pos_cols <- c("pos_x", "pos_y", "pos_z")
  if (!all(pos_cols %in% colnames(cd))) {
    stop("Electrode positions (pos_x, pos_y, pos_z) are required in colData. ",
         "Apply eegMontage() first.")
  }

  bad_idx <- which(labels %in% bad_channels)
  if (length(bad_idx) == 0) {
    warning("No bad channels found in colData$label. Returning unchanged data.")
    SummarizedExperiment::assay(x, output_assay) <- data
    return(x)
  }

  good_idx <- setdiff(seq_len(ncol(data)), bad_idx)
  if (length(good_idx) < 3) {
    stop("Need at least 3 good channels for interpolation. Only ",
         length(good_idx), " good channels available.")
  }

  pos <- as.matrix(cd[, pos_cols])
  good_pos <- pos[good_idx, , drop = FALSE]
  bad_pos <- pos[bad_idx, , drop = FALSE]

  result <- data

  if (method == "spline") {
    # Spherical spline interpolation
    weights <- .spherical_spline_weights(good_pos, bad_pos)

    # weights is (n_bad x n_good), data for good channels is (n_time x n_good)
    good_data <- data[, good_idx, drop = FALSE]
    interpolated <- good_data %*% t(weights)

    for (i in seq_along(bad_idx)) {
      result[, bad_idx[i]] <- interpolated[, i]
    }

  } else {
    # Nearest neighbor (inverse distance weighting)
    for (i in seq_along(bad_idx)) {
      b_pos <- bad_pos[i, ]
      dists <- sqrt(rowSums((good_pos - matrix(b_pos, nrow = nrow(good_pos),
                                                ncol = 3, byrow = TRUE))^2))
      # Prevent zero distances
      dists[dists < 1e-10] <- 1e-10

      # Inverse distance weights
      inv_dists <- 1 / dists
      w <- inv_dists / sum(inv_dists)

      # Weighted average of good channels
      result[, bad_idx[i]] <- data[, good_idx, drop = FALSE] %*% w
    }
  }

  SummarizedExperiment::assay(x, output_assay) <- result
  x
}


#' Epoch continuous EEG data
#'
#' Segments continuous 2D EEG data into fixed-length epochs around events,
#' producing a 3D array (time x channels x epochs). Optionally performs
#' baseline correction by subtracting the mean of a pre-stimulus window.
#'
#' @param x A PhysioExperiment object with continuous (2D) data.
#' @param events A data.frame with an \code{onset_sec} column specifying
#'   event times in seconds, OR an integer vector of sample indices.
#' @param limits Numeric vector of length 2 specifying the epoch window
#'   relative to each event in seconds (default: \code{c(-0.2, 0.8)}).
#' @param baseline Numeric vector of length 2 specifying the baseline window
#'   relative to each event in seconds (default: \code{c(-0.2, 0)}). Set to
#'   NULL to skip baseline correction.
#' @param assay_name Name of the assay to epoch. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name of the output assay (default: \code{"epoched"}).
#'
#' @return A PhysioExperiment object with a 3D array (time x channels x epochs)
#'   in the specified output assay. Event information is stored in
#'   \code{metadata(x)$epoch_events}.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
#' events <- data.frame(onset_sec = c(1.0, 3.0, 5.0, 7.0))
#' pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))
#' dim(SummarizedExperiment::assay(pe_ep, "epoched"))
#' # time x channels x epochs
#' }
eegEpoch <- function(x, events, limits = c(-0.2, 0.8), baseline = c(-0.2, 0),
                      assay_name = NULL, output_assay = "epoched") {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  sr <- samplingRate(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  if (!is.matrix(data)) {
    stop("eegEpoch requires a 2D matrix assay (continuous data).")
  }

  n_time <- nrow(data)
  n_ch <- ncol(data)

  stopifnot(is.numeric(limits), length(limits) == 2)
  if (limits[1] >= limits[2]) {
    stop("limits[1] must be less than limits[2].")
  }

  # Convert events to sample indices
  if (is.data.frame(events)) {
    if (!"onset_sec" %in% colnames(events)) {
      stop("events data.frame must have an 'onset_sec' column.")
    }
    event_samples <- as.integer(round(events$onset_sec * sr)) + 1L
  } else if (is.numeric(events)) {
    event_samples <- as.integer(events)
  } else {
    stop("events must be a data.frame with onset_sec column or an integer vector.")
  }

  # Compute epoch sample offsets
  offset_start <- as.integer(round(limits[1] * sr))
  offset_end <- as.integer(round(limits[2] * sr))
  epoch_length <- offset_end - offset_start + 1L

  # Baseline offsets (relative to epoch start)
  if (!is.null(baseline)) {
    stopifnot(is.numeric(baseline), length(baseline) == 2)
    bl_start <- as.integer(round(baseline[1] * sr)) - offset_start + 1L
    bl_end <- as.integer(round(baseline[2] * sr)) - offset_start + 1L
    bl_start <- max(1L, bl_start)
    bl_end <- min(epoch_length, bl_end)
  }

  # Filter out edge events that would exceed data boundaries
  valid_events <- logical(length(event_samples))
  for (i in seq_along(event_samples)) {
    start_idx <- event_samples[i] + offset_start
    end_idx <- event_samples[i] + offset_end
    valid_events[i] <- (start_idx >= 1L) && (end_idx <= n_time)
  }

  if (sum(valid_events) == 0) {
    stop("No valid events found within data boundaries. ",
         "Check event times and epoch limits.")
  }

  if (sum(!valid_events) > 0) {
    message(sprintf("Excluding %d of %d events (too close to data boundaries).",
                    sum(!valid_events), length(event_samples)))
  }

  valid_samples <- event_samples[valid_events]
  n_epochs <- length(valid_samples)

  # Extract epochs into 3D array
  epoch_data <- array(NA_real_, dim = c(epoch_length, n_ch, n_epochs))

  for (ep in seq_len(n_epochs)) {
    start_idx <- valid_samples[ep] + offset_start
    end_idx <- valid_samples[ep] + offset_end
    epoch_data[, , ep] <- data[start_idx:end_idx, ]
  }

  # Baseline correction
  if (!is.null(baseline)) {
    for (ep in seq_len(n_epochs)) {
      for (ch in seq_len(n_ch)) {
        bl_mean <- mean(epoch_data[bl_start:bl_end, ch, ep], na.rm = TRUE)
        epoch_data[, ch, ep] <- epoch_data[, ch, ep] - bl_mean
      }
    }
  }

  # Store epoch info in metadata
  epoch_info <- list(
    event_samples = valid_samples,
    limits = limits,
    baseline = baseline,
    n_epochs = n_epochs,
    n_excluded = sum(!valid_events),
    epoch_length = epoch_length,
    timestamp = Sys.time()
  )

  # If events was a data.frame, store the valid subset

  if (is.data.frame(events)) {
    epoch_info$event_table <- events[valid_events, , drop = FALSE]
  }

  md <- S4Vectors::metadata(x)
  md$epoch_events <- epoch_info
  S4Vectors::metadata(x) <- md

  md[[output_assay]] <- epoch_data
  S4Vectors::metadata(x) <- md
  x
}


#' Reject artifacts in epoched EEG data
#'
#' Identifies and removes epochs contaminated by artifacts from 3D epoched
#' EEG data. Supports multiple detection criteria: amplitude threshold,
#' gradient (point-to-point voltage change), and joint probability.
#'
#' @param x A PhysioExperiment object with 3D epoched data.
#' @param method Artifact detection method: \code{"threshold"} (amplitude),
#'   \code{"gradient"} (point-to-point change), or \code{"joint_probability"}
#'   (log-power z-score).
#' @param threshold_uv Maximum absolute amplitude in microvolts for threshold
#'   rejection (default: 100).
#' @param gradient_uv_ms Maximum point-to-point change in uV/ms for gradient
#'   rejection (default: 50).
#' @param jp_threshold Number of standard deviations for joint probability
#'   rejection (default: 3).
#' @param assay_name Name of the assay to check. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name of the output assay (default: \code{"clean"}).
#'
#' @return A PhysioExperiment object with clean epochs in the specified output
#'   assay. Artifact log stored in \code{metadata(x)$artifact_log}.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
#' events <- data.frame(onset_sec = c(1, 3, 5, 7, 9))
#' pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))
#' pe_clean <- eegArtifactReject(pe_ep, method = "threshold",
#'                                threshold_uv = 100, assay_name = "epoched")
#' }
eegArtifactReject <- function(x, method = c("threshold", "gradient",
                                             "joint_probability"),
                               threshold_uv = 100, gradient_uv_ms = 50,
                               jp_threshold = 3, assay_name = NULL,
                               output_assay = "clean") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- "epoched"
  sr <- samplingRate(x)
  data <- S4Vectors::metadata(x)[[assay_name]]
  if (is.null(data)) {
    stop("eegArtifactReject requires 3D epoched data in metadata. ",
         "Run eegEpoch() first.", call. = FALSE)
  }

  # Ensure 3D data
  d <- dim(data)
  if (length(d) != 3) {
    stop("eegArtifactReject requires 3D epoched data (time x channels x epochs). ",
         "Run eegEpoch() first.")
  }

  n_time <- d[1]
  n_ch <- d[2]
  n_epochs <- d[3]

  bad_epochs <- logical(n_epochs)
  artifact_log <- data.frame(
    epoch = integer(0),
    reason = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )

  if (method == "threshold") {
    for (ep in seq_len(n_epochs)) {
      max_abs <- max(abs(data[, , ep]), na.rm = TRUE)
      if (max_abs > threshold_uv) {
        bad_epochs[ep] <- TRUE
        artifact_log <- rbind(artifact_log, data.frame(
          epoch = ep, reason = "threshold", value = max_abs,
          stringsAsFactors = FALSE
        ))
      }
    }

  } else if (method == "gradient") {
    # Convert gradient threshold from uV/ms to uV/sample
    grad_threshold_per_sample <- gradient_uv_ms * (1000 / sr)

    for (ep in seq_len(n_epochs)) {
      epoch_data <- data[, , ep]
      if (n_time < 2) next
      # Point-to-point differences
      diffs <- abs(epoch_data[-1, , drop = FALSE] - epoch_data[-n_time, , drop = FALSE])
      max_diff <- max(diffs, na.rm = TRUE)
      if (max_diff > grad_threshold_per_sample) {
        bad_epochs[ep] <- TRUE
        artifact_log <- rbind(artifact_log, data.frame(
          epoch = ep, reason = "gradient", value = max_diff,
          stringsAsFactors = FALSE
        ))
      }
    }

  } else if (method == "joint_probability") {
    # Compute log-power for each epoch
    log_powers <- numeric(n_epochs)
    for (ep in seq_len(n_epochs)) {
      epoch_data <- data[, , ep]
      power <- mean(epoch_data^2, na.rm = TRUE)
      log_powers[ep] <- log(max(power, 1e-20))
    }

    mean_lp <- mean(log_powers)
    sd_lp <- stats::sd(log_powers)
    if (sd_lp < 1e-12) sd_lp <- 1e-12

    z_scores <- (log_powers - mean_lp) / sd_lp

    for (ep in seq_len(n_epochs)) {
      if (abs(z_scores[ep]) > jp_threshold) {
        bad_epochs[ep] <- TRUE
        artifact_log <- rbind(artifact_log, data.frame(
          epoch = ep, reason = "joint_probability", value = z_scores[ep],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Remove bad epochs
  good_idx <- which(!bad_epochs)

  if (length(good_idx) == 0) {
    warning("All epochs rejected. Returning empty 3D array.")
    clean_data <- array(numeric(0), dim = c(n_time, n_ch, 0))
  } else {
    clean_data <- data[, , good_idx, drop = FALSE]
  }

  n_rejected <- sum(bad_epochs)
  if (n_rejected > 0) {
    message(sprintf("Rejected %d of %d epochs (%.1f%%) using %s method.",
                    n_rejected, n_epochs,
                    100 * n_rejected / n_epochs, method))
  }

  # Store artifact log
  md <- S4Vectors::metadata(x)
  md$artifact_log <- artifact_log
  md$artifact_summary <- list(
    method = method,
    n_total = n_epochs,
    n_rejected = n_rejected,
    n_retained = length(good_idx),
    rejected_epochs = which(bad_epochs),
    timestamp = Sys.time()
  )
  md[[output_assay]] <- clean_data
  S4Vectors::metadata(x) <- md

  x
}


#' Assign electrode montage positions
#'
#' Maps channel labels to standard electrode positions from a known montage
#' system (10-20, 10-10, BioSemi64) or from a custom positions data.frame.
#' Sets \code{pos_x}, \code{pos_y}, \code{pos_z} columns in colData for use
#' by interpolation, topographic mapping, and source localization functions.
#'
#' @param x A PhysioExperiment object.
#' @param system Montage system: \code{"10-20"} (19 channels),
#'   \code{"10-10"} (~64 channels), \code{"biosemi64"} (64 BioSemi channels),
#'   or \code{"custom"} (user-provided positions).
#' @param positions For \code{system = "custom"}, a data.frame with columns
#'   \code{label}, \code{pos_x}, \code{pos_y}, \code{pos_z}.
#'
#' @return A PhysioExperiment object with updated colData containing position
#'   columns.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe <- eegMontage(pe, system = "10-20")
#' head(SummarizedExperiment::colData(pe))
#' }
eegMontage <- function(x, system = c("10-20", "10-10", "biosemi64", "custom"),
                        positions = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  system <- match.arg(system)

  cd <- SummarizedExperiment::colData(x)
  labels <- as.character(cd$label)

  if (system == "custom") {
    if (is.null(positions) || !is.data.frame(positions)) {
      stop("positions must be a data.frame with columns: label, pos_x, pos_y, pos_z.")
    }
    required_cols <- c("label", "pos_x", "pos_y", "pos_z")
    if (!all(required_cols %in% colnames(positions))) {
      stop("positions data.frame must have columns: ",
           paste(required_cols, collapse = ", "))
    }
    pos_db <- positions

  } else if (system == "10-20") {
    pos_db <- .electrode_positions_1020()

  } else if (system == "10-10") {
    pos_db <- .electrode_positions_1010()

  } else if (system == "biosemi64") {
    # BioSemi 64 uses extended 10-10 positions with BioSemi labeling
    # Map common BioSemi labels to 10-10 equivalents
    pos_db <- .electrode_positions_1010()
    # BioSemi uses T7/T8/P7/P8 instead of T3/T4/T5/T6
    biosemi_map <- data.frame(
      from = c("T7", "T8", "P7", "P8"),
      to = c("T3", "T4", "T5", "T6"),
      stringsAsFactors = FALSE
    )
    # Add aliased labels with same positions
    for (i in seq_len(nrow(biosemi_map))) {
      idx <- which(pos_db$label == biosemi_map$to[i])
      if (length(idx) > 0) {
        alias_row <- pos_db[idx[1], ]
        alias_row$label <- biosemi_map$from[i]
        pos_db <- rbind(pos_db, alias_row)
      }
    }
  }

  # Match channel labels to position database (case-insensitive)
  matched <- match(toupper(labels), toupper(pos_db$label))

  n_matched <- sum(!is.na(matched))
  n_total <- length(labels)

  if (n_matched == 0) {
    stop("No channel labels matched the ", system, " montage system. ",
         "Labels found: ", paste(head(labels, 5), collapse = ", "))
  }

  if (n_matched < n_total) {
    unmatched <- labels[is.na(matched)]
    warning(sprintf(
      "%d of %d channels not found in %s system: %s. Setting positions to NA.",
      length(unmatched), n_total, system,
      paste(head(unmatched, 5), collapse = ", ")
    ))
  }

  # Assign positions
  cd$pos_x <- ifelse(is.na(matched), NA_real_, pos_db$pos_x[matched])
  cd$pos_y <- ifelse(is.na(matched), NA_real_, pos_db$pos_y[matched])
  cd$pos_z <- ifelse(is.na(matched), NA_real_, pos_db$pos_z[matched])

  SummarizedExperiment::colData(x) <- cd

  # Store montage info
  md <- S4Vectors::metadata(x)
  md$montage <- list(
    system = system,
    n_matched = n_matched,
    n_total = n_total,
    unmatched = labels[is.na(matched)],
    timestamp = Sys.time()
  )
  S4Vectors::metadata(x) <- md

  x
}


#' Full EEG preprocessing pipeline
#'
#' Convenience wrapper that runs a complete EEG preprocessing pipeline in
#' sequence: filtering, re-referencing, bad channel detection and interpolation,
#' optional ICA, epoching, and artifact rejection. Each step can be toggled
#' on or off.
#'
#' @param x A PhysioExperiment object with continuous (2D) EEG data.
#' @param filter Logical; apply frequency filtering (default: TRUE).
#' @param lowcut Low cutoff frequency in Hz for filtering (default: 0.1).
#' @param highcut High cutoff frequency in Hz for filtering (default: 40).
#' @param notch Notch filter center frequency in Hz (NULL for none).
#' @param rereference Logical; apply re-referencing (default: TRUE).
#' @param ref_type Re-referencing type (default: \code{"average"}).
#' @param bad_channels Logical; detect and report bad channels (default: TRUE).
#' @param interpolate Logical; interpolate detected bad channels (default: TRUE).
#'   Requires electrode positions in colData (apply \code{eegMontage} first or
#'   this will be skipped with a warning).
#' @param ica Logical; apply ICA-based artifact removal (default: FALSE).
#'   Requires the \code{eegICA} function to be available.
#' @param epoch Logical; epoch the data around events (default: FALSE).
#' @param events Event data for epoching (data.frame or integer vector).
#'   Required if \code{epoch = TRUE}.
#' @param epoch_limits Epoch window in seconds (default: \code{c(-0.2, 0.8)}).
#' @param baseline Baseline window in seconds (default: \code{c(-0.2, 0)}).
#' @param artifact_reject Logical; reject artifact epochs (default: FALSE).
#'   Only applies if \code{epoch = TRUE}.
#' @param threshold_uv Amplitude threshold for artifact rejection (default: 100).
#' @param assay_name Starting assay name. If NULL, uses \code{defaultAssay(x)}.
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return A PhysioExperiment object with processed data. The final assay name
#'   depends on which steps are enabled. Processing log is stored in
#'   \code{metadata(x)$preprocess_log}.
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
#' pe <- eegMontage(pe, system = "10-20")
#' events <- data.frame(onset_sec = c(1, 3, 5, 7, 9))
#' pe_proc <- eegPreprocess(pe, lowcut = 1, highcut = 40, notch = 50,
#'                           epoch = TRUE, events = events,
#'                           artifact_reject = TRUE)
#' }
eegPreprocess <- function(x, filter = TRUE, lowcut = 0.1, highcut = 40,
                           notch = NULL, rereference = TRUE,
                           ref_type = "average", bad_channels = TRUE,
                           interpolate = TRUE, ica = FALSE,
                           epoch = FALSE, events = NULL,
                           epoch_limits = c(-0.2, 0.8),
                           baseline = c(-0.2, 0),
                           artifact_reject = FALSE, threshold_uv = 100,
                           assay_name = NULL, verbose = TRUE) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  current_assay <- assay_name

  # Initialize preprocessing log
  preprocess_log <- list()

  .log_step <- function(step_name, params) {
    entry <- list(
      step = step_name,
      params = params,
      timestamp = as.character(Sys.time())
    )
    preprocess_log[[length(preprocess_log) + 1L]] <<- entry
  }

  # Step 1: Filtering
  if (filter) {
    if (verbose) message("[eegPreprocess] Step 1: Filtering (", lowcut, "-",
                         highcut, " Hz", if (!is.null(notch)) paste0(", notch ", notch, " Hz"), ")...")
    x <- eegFilter(x, lowcut = lowcut, highcut = highcut, notch = notch,
                    assay_name = current_assay, output_assay = "filtered")
    current_assay <- "filtered"
    .log_step("filter", list(lowcut = lowcut, highcut = highcut, notch = notch))
  }

  # Step 2: Bad channel detection
  bad_labels <- character(0)
  if (bad_channels) {
    if (verbose) message("[eegPreprocess] Step 2: Bad channel detection...")
    bad_df <- eegBadChannels(x, method = "all", assay_name = current_assay)
    bad_labels <- bad_df$channel[bad_df$is_bad]
    if (verbose) {
      if (length(bad_labels) > 0) {
        message("  Found ", length(bad_labels), " bad channel(s): ",
                paste(bad_labels, collapse = ", "))
      } else {
        message("  No bad channels detected.")
      }
    }
    .log_step("bad_channel_detection", list(
      n_bad = length(bad_labels),
      bad_channels = bad_labels
    ))
  }

  # Step 3: Interpolation
  if (interpolate && length(bad_labels) > 0) {
    # Check if positions are available
    cd <- SummarizedExperiment::colData(x)
    has_pos <- all(c("pos_x", "pos_y", "pos_z") %in% colnames(cd))

    if (has_pos) {
      if (verbose) message("[eegPreprocess] Step 3: Interpolating bad channels...")
      x <- eegInterpolate(x, bad_labels, method = "spline",
                           assay_name = current_assay,
                           output_assay = "interpolated")
      current_assay <- "interpolated"
      .log_step("interpolation", list(channels = bad_labels, method = "spline"))
    } else {
      if (verbose) message("[eegPreprocess] Step 3: Skipping interpolation ",
                           "(no electrode positions). Apply eegMontage() first.")
      .log_step("interpolation", list(skipped = TRUE,
                                       reason = "no electrode positions"))
    }
  }

  # Step 4: Re-referencing
  if (rereference) {
    if (verbose) message("[eegPreprocess] Step 4: Re-referencing (", ref_type, ")...")
    x <- eegRereference(x, ref_type = ref_type, assay_name = current_assay,
                         output_assay = "rereferenced")
    current_assay <- "rereferenced"
    .log_step("rereference", list(ref_type = ref_type))
  }

  # Step 5: ICA
  if (ica) {
    if (verbose) message("[eegPreprocess] Step 5: ICA artifact removal...")
    if (exists("eegICA", mode = "function")) {
      x <- tryCatch({
        eegICA(x, assay_name = current_assay)
      }, error = function(e) {
        warning("ICA step failed: ", conditionMessage(e), ". Skipping.")
        x
      })
      # Check if ICA produced a cleaned assay
      anames <- SummarizedExperiment::assayNames(x)
      if ("ica_cleaned" %in% anames) {
        current_assay <- "ica_cleaned"
      }
      .log_step("ica", list(completed = TRUE))
    } else {
      if (verbose) message("  eegICA not available. Skipping ICA step.")
      .log_step("ica", list(skipped = TRUE, reason = "eegICA not available"))
    }
  }

  # Step 6: Epoching
  if (epoch) {
    if (is.null(events)) {
      stop("events must be provided when epoch = TRUE.")
    }
    if (verbose) message("[eegPreprocess] Step 6: Epoching (limits: ",
                         epoch_limits[1], " to ", epoch_limits[2], " sec)...")
    x <- eegEpoch(x, events = events, limits = epoch_limits,
                   baseline = baseline, assay_name = current_assay,
                   output_assay = "epoched")
    current_assay <- "epoched"
    .log_step("epoch", list(limits = epoch_limits, baseline = baseline))
  }

  # Step 7: Artifact rejection
  if (artifact_reject && epoch) {
    if (verbose) message("[eegPreprocess] Step 7: Artifact rejection (threshold: ",
                         threshold_uv, " uV)...")
    x <- eegArtifactReject(x, method = "threshold", threshold_uv = threshold_uv,
                            assay_name = current_assay, output_assay = "clean")
    current_assay <- "clean"
    .log_step("artifact_reject", list(method = "threshold",
                                       threshold_uv = threshold_uv))
  } else if (artifact_reject && !epoch) {
    if (verbose) message("[eegPreprocess] Skipping artifact rejection (data not epoched).")
    .log_step("artifact_reject", list(skipped = TRUE,
                                       reason = "data not epoched"))
  }

  # Store full preprocessing log
  md <- S4Vectors::metadata(x)
  md$preprocess_log <- preprocess_log
  md$preprocess_final_assay <- current_assay
  S4Vectors::metadata(x) <- md

  if (verbose) message("[eegPreprocess] Done. Final assay: '", current_assay, "'")

  x
}
