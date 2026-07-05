#' Morlet Wavelet Transform for EEG
#'
#' Computes the continuous Morlet wavelet transform for multi-channel EEG data.
#' For each specified frequency, a complex Morlet wavelet is constructed and
#' convolved with each channel using FFT-based convolution for efficiency.
#' Returns time-resolved power (and optionally phase) across frequencies.
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param frequencies Numeric vector of frequencies in Hz to analyze. If
#'   \code{NULL}, defaults to \code{seq(1, 50, by = 1)}.
#' @param n_cycles Number of cycles in the Morlet wavelet, controlling the
#'   trade-off between time and frequency resolution (default: 7).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store wavelet power results
#'   (default: \code{"wavelet_power"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item 3D power array (time x frequencies x channels) in \code{output_assay}
#'     \item Frequency vector and phase array in \code{metadata(x)$wavelet},
#'       a list containing \code{frequencies} (numeric vector),
#'       \code{n_cycles} (integer), and \code{phase} (3D array of
#'       instantaneous phase values)
#'   }
#' @references
#' Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
#' conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.
#' @seealso [eegSTFT()], [eegMultitaper()], [eegPlotSpectrogram()],
#'   [eegERSP()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' pe_wt <- eegMorletWavelet(pe, frequencies = seq(5, 40, by = 1))
#' wp <- SummarizedExperiment::assay(pe_wt, "wavelet_power")
#' dim(wp)  # time x frequencies x channels
#' }
eegMorletWavelet <- function(x, frequencies = NULL, n_cycles = 7,
                             assay_name = NULL, output_assay = "wavelet_power") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(n_cycles) && n_cycles > 0)

  if (is.null(frequencies)) frequencies <- seq(1, 50, by = 1)
  stopifnot(is.numeric(frequencies) && length(frequencies) > 0)
  stopifnot(all(frequencies > 0))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Ensure 2D input
  if (length(dim(data)) != 2) {
    stop("eegMorletWavelet requires 2D data (time x channels). ",
         "For epoched data, use eegERSP or eegITC.", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)
  n_freqs <- length(frequencies)

  # Allocate output arrays
  power_out <- array(0, dim = c(n_time, n_freqs, n_channels))

  phase_out <- array(0, dim = c(n_time, n_freqs, n_channels))

  for (fi in seq_len(n_freqs)) {
    f <- frequencies[fi]

    # Construct complex Morlet wavelet
    sigma <- n_cycles / (2 * pi * f)
    wavelet_half <- ceiling(3 * sigma * sr)
    wavelet_t <- seq(-wavelet_half, wavelet_half) / sr
    wavelet_n <- length(wavelet_t)

    # Complex Morlet: exp(2*pi*i*f*t) * exp(-t^2 / (2*sigma^2))
    wavelet <- exp(2i * pi * f * wavelet_t) * exp(-wavelet_t^2 / (2 * sigma^2))
    # Normalize for unit energy
    wavelet <- wavelet / sqrt(sum(Mod(wavelet)^2))

    # FFT convolution: pad to next power of 2
    n_conv <- n_time + wavelet_n - 1L
    n_fft <- 2^ceiling(log2(n_conv))
    wavelet_fft <- fft(c(wavelet, rep(0, n_fft - wavelet_n)))

    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      sig_fft <- fft(c(sig, rep(0, n_fft - n_time)))

      # Convolution via FFT
      conv_result <- fft(sig_fft * Conj(wavelet_fft), inverse = TRUE) / n_fft

      # Trim to original signal length (center of convolution)
      half_wav <- wavelet_half
      start_idx <- half_wav + 1L
      end_idx <- half_wav + n_time
      trimmed <- conv_result[start_idx:end_idx]

      power_out[, fi, ch] <- Mod(trimmed)^2
      phase_out[, fi, ch] <- Arg(trimmed)
    }
  }

  # Store power in metadata (different dimensions than SE assays)
  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- power_out
  md$wavelet <- list(
    frequencies = frequencies,
    n_cycles = n_cycles,
    phase = phase_out
  )
  S4Vectors::metadata(x) <- md

  x
}


#' Short-Time Fourier Transform for EEG
#'
#' Computes the Short-Time Fourier Transform (STFT) spectrogram for
#' multi-channel EEG data. Uses a sliding window with configurable overlap
#' and window function to produce a time-frequency power representation.
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param window_sec Window length in seconds (default: 0.5).
#' @param overlap Overlap fraction between adjacent windows, from 0 to 1
#'   exclusive (default: 0.75).
#' @param window_type Window function to apply: \code{"hanning"},
#'   \code{"hamming"}, or \code{"rectangular"} (default: \code{"hanning"}).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store STFT power results
#'   (default: \code{"stft_power"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item 3D power array (time_bins x frequencies x channels) in
#'       \code{output_assay}
#'     \item Time bin centers, frequency vector, and parameters in
#'       \code{metadata(x)$stft}, a list containing \code{time_axis},
#'       \code{freq_axis}, \code{window_sec}, \code{overlap},
#'       \code{window_type}, \code{window_length}, and \code{hop_size}
#'   }
#' @references
#' Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
#' conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.
#' @seealso [eegMorletWavelet()], [eegMultitaper()], [eegPlotSpectrogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
#' pe_stft <- eegSTFT(pe, window_sec = 0.5, overlap = 0.75)
#' sp <- SummarizedExperiment::assay(pe_stft, "stft_power")
#' dim(sp)  # time_bins x frequencies x channels
#' }
eegSTFT <- function(x, window_sec = 0.5, overlap = 0.75,
                    window_type = c("hanning", "hamming", "rectangular"),
                    assay_name = NULL, output_assay = "stft_power") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(window_sec) && window_sec > 0)
  stopifnot(is.numeric(overlap) && overlap >= 0 && overlap < 1)
  window_type <- match.arg(window_type)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 2) {
    stop("eegSTFT requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Window parameters
  window_length <- as.integer(round(window_sec * sr))
  # Ensure even window length for symmetric frequency axis
  if (window_length %% 2 != 0) window_length <- window_length + 1L
  hop_size <- max(1L, as.integer(round(window_length * (1 - overlap))))

  # Create window function
  win_idx <- seq(0, window_length - 1)
  window_fn <- switch(window_type,
    hanning = 0.5 * (1 - cos(2 * pi * win_idx / (window_length - 1))),
    hamming = 0.54 - 0.46 * cos(2 * pi * win_idx / (window_length - 1)),
    rectangular = rep(1, window_length)
  )

  # Normalization factor for power
  win_norm <- sum(window_fn^2)

  # Frequency axis: positive frequencies only
  n_fft_bins <- window_length %/% 2 + 1L
  freq_axis <- seq(0, sr / 2, length.out = n_fft_bins)

  # Time bin positions
  starts <- seq(1L, n_time - window_length + 1L, by = hop_size)
  n_bins <- length(starts)

  if (n_bins < 1) {
    stop("Signal is too short for the specified window length.", call. = FALSE)
  }

  # Time axis: center of each window in seconds
  time_axis <- (starts - 1 + window_length / 2) / sr

  # Allocate output: time_bins x frequencies x channels
  stft_power <- array(0, dim = c(n_bins, n_fft_bins, n_channels))

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    for (ti in seq_len(n_bins)) {
      seg <- sig[starts[ti]:(starts[ti] + window_length - 1L)]
      seg_win <- seg * window_fn
      ft <- fft(seg_win)
      # Keep positive frequencies: indices 1 to n_fft_bins
      ft_pos <- ft[seq_len(n_fft_bins)]
      # Power spectrum, normalized by window energy
      stft_power[ti, , ch] <- Mod(ft_pos)^2 / win_norm
    }
  }

  # Store in metadata (different dimensions than SE assays)
  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- stft_power
  md$stft <- list(
    time_axis = time_axis,
    freq_axis = freq_axis,
    window_sec = window_sec,
    overlap = overlap,
    window_type = window_type,
    window_length = window_length,
    hop_size = hop_size
  )
  S4Vectors::metadata(x) <- md

  x
}


#' Multitaper Power Spectral Density for EEG
#'
#' Computes the multitaper power spectral density (PSD) estimate for
#' multi-channel EEG data using Discrete Prolate Spheroidal Sequences (DPSS,
#' Slepian tapers). The multitaper method provides an optimal bias-variance
#' trade-off for spectral estimation compared to single-window methods.
#'
#' @param x A PhysioExperiment object with EEG data (2D: time x channels).
#' @param bandwidth Time-half-bandwidth product (NW) controlling the spectral
#'   concentration of the tapers (default: 4). Higher values give smoother but
#'   lower-resolution estimates.
#' @param n_tapers Number of DPSS tapers to use. If \code{NULL}, defaults to
#'   \code{floor(2 * bandwidth) - 1}.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store PSD results
#'   (default: \code{"multitaper_psd"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item PSD matrix (frequencies x channels) in \code{output_assay}
#'     \item Frequency vector and taper parameters in
#'       \code{metadata(x)$multitaper}, a list containing
#'       \code{frequencies} (numeric vector), \code{bandwidth}
#'       (numeric NW parameter), and \code{n_tapers} (integer)
#'   }
#' @references
#' Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
#' conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.
#'
#' Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
#' Proceedings of the IEEE, 70(9), 1055-1096.
#' @seealso [eegMorletWavelet()], [eegSTFT()], [eegPlotSpectrogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe_mt <- eegMultitaper(pe, bandwidth = 4)
#' psd <- SummarizedExperiment::assay(pe_mt, "multitaper_psd")
#' dim(psd)  # frequencies x channels
#' }
eegMultitaper <- function(x, bandwidth = 4, n_tapers = NULL,
                          assay_name = NULL, output_assay = "multitaper_psd") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(bandwidth) && bandwidth > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 2) {
    stop("eegMultitaper requires 2D data (time x channels).", call. = FALSE)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (is.null(n_tapers)) n_tapers <- as.integer(floor(2 * bandwidth) - 1)
  stopifnot(is.numeric(n_tapers) && n_tapers >= 1)

  # Generate DPSS tapers via tridiagonal matrix eigenvalue decomposition
  tapers <- .compute_dpss(n_time, bandwidth, n_tapers)

  # Frequency axis: positive frequencies
  n_fft <- n_time
  n_freq <- n_fft %/% 2 + 1L
  freq_axis <- seq(0, sr / 2, length.out = n_freq)

  # Allocate PSD output: frequencies x channels
  psd_out <- matrix(0, nrow = n_freq, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    psd_accum <- numeric(n_freq)

    for (k in seq_len(n_tapers)) {
      tapered_sig <- sig * tapers[, k]
      ft <- fft(tapered_sig)
      # Keep positive frequencies
      ft_pos <- ft[seq_len(n_freq)]
      psd_accum <- psd_accum + Mod(ft_pos)^2
    }

    # Average across tapers and normalize
    psd_out[, ch] <- psd_accum / (n_tapers * sr)
  }

  # Store in metadata (different dimensions than SE assays)
  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- psd_out
  md$multitaper <- list(
    frequencies = freq_axis,
    bandwidth = bandwidth,
    n_tapers = n_tapers
  )
  S4Vectors::metadata(x) <- md

  x
}


#' Event-Related Spectral Perturbation (ERSP)
#'
#' Computes Event-Related Spectral Perturbation for epoched (3D) EEG data.
#' ERSP quantifies event-related changes in spectral power relative to a
#' baseline period, expressed in decibels (dB). Uses the Morlet wavelet
#' transform to compute time-frequency decomposition for each epoch, then
#' averages power across epochs and normalizes to baseline.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x epochs).
#' @param baseline Numeric vector of length 2 specifying the baseline time
#'   window in sample indices (e.g., \code{c(1, 50)} for the first 50 samples).
#'   Default is \code{c(1, 50)}.
#' @param frequencies Numeric vector of frequencies in Hz to analyze. If
#'   \code{NULL}, defaults to \code{seq(1, 50, by = 1)}.
#' @param n_cycles Number of cycles for the Morlet wavelet (default: 7).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store ERSP results
#'   (default: \code{"ersp"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item 3D ERSP array (time x frequencies x channels) in dB in
#'       \code{output_assay}
#'     \item Baseline info and frequency vector in \code{metadata(x)$ersp},
#'       a list containing \code{frequencies} (numeric vector),
#'       \code{baseline} (numeric vector), \code{n_cycles} (integer),
#'       and \code{n_epochs} (integer)
#'   }
#' @references
#' Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
#' conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.
#'
#' Makeig, S. (1993). Auditory event-related dynamics of the EEG spectrum and
#' effects of exposure to tones. Electroencephalography and Clinical
#' Neurophysiology, 86(4), 283-293.
#' @seealso [eegITC()], [eegMorletWavelet()], [eegPlotSpectrogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
#' pe_ersp <- eegERSP(pe, baseline = c(1, 50), frequencies = seq(5, 30, by = 5))
#' ersp_data <- SummarizedExperiment::assay(pe_ersp, "ersp")
#' dim(ersp_data)  # time x frequencies x channels
#' }
eegERSP <- function(x, baseline = c(1, 50), frequencies = NULL, n_cycles = 7,
                    assay_name = NULL, output_assay = "ersp_data") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(baseline) && length(baseline) == 2)
  stopifnot(baseline[1] < baseline[2])
  stopifnot(is.numeric(n_cycles) && n_cycles > 0)

  if (is.null(frequencies)) frequencies <- seq(1, 50, by = 1)
  stopifnot(is.numeric(frequencies) && length(frequencies) > 0)
  stopifnot(all(frequencies > 0))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 3) {
    stop("eegERSP requires epoched (3D) data (time x channels x epochs).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_epochs <- dim(data)[3]
  n_freqs <- length(frequencies)

  # Validate baseline indices
  bl_start <- max(1L, as.integer(baseline[1]))
  bl_end <- min(n_time, as.integer(baseline[2]))
  stopifnot(bl_start < bl_end)

  # Accumulate power across epochs: time x freq x channel
  power_sum <- array(0, dim = c(n_time, n_freqs, n_channels))

  for (ep in seq_len(n_epochs)) {
    epoch_power <- .morlet_power_2d(data[, , ep], sr, frequencies, n_cycles)
    power_sum <- power_sum + epoch_power
  }

  # Mean power across epochs
  mean_power <- power_sum / n_epochs

  # Compute baseline mean power per frequency per channel
  # baseline_power: n_freqs x n_channels
  ersp_out <- array(0, dim = c(n_time, n_freqs, n_channels))

  for (ch in seq_len(n_channels)) {
    for (fi in seq_len(n_freqs)) {
      bl_mean <- mean(mean_power[bl_start:bl_end, fi, ch])
      # Guard against zero baseline
      bl_mean <- max(bl_mean, .Machine$double.eps)
      # ERSP in dB
      ersp_out[, fi, ch] <- 10 * log10(mean_power[, fi, ch] / bl_mean)
    }
  }

  # Store in metadata (different dimensions than SE assays)
  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- ersp_out
  md$ersp <- list(
    frequencies = frequencies,
    baseline = baseline,
    n_cycles = n_cycles,
    n_epochs = n_epochs
  )
  S4Vectors::metadata(x) <- md

  x
}


#' Inter-Trial Coherence (ITC)
#'
#' Computes Inter-Trial Coherence (also known as phase-locking factor or
#' phase-locking value) for epoched (3D) EEG data. ITC measures the
#' consistency of oscillatory phase across trials at each time-frequency
#' point. Values range from 0 (completely random phase) to 1 (perfectly
#' phase-locked across all trials).
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x epochs).
#' @param frequencies Numeric vector of frequencies in Hz to analyze. If
#'   \code{NULL}, defaults to \code{seq(1, 50, by = 1)}.
#' @param n_cycles Number of cycles for the Morlet wavelet (default: 7).
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name of the assay to store ITC results
#'   (default: \code{"itc"}).
#' @return Modified PhysioExperiment with:
#'   \itemize{
#'     \item 3D ITC array (time x frequencies x channels) in
#'       \code{output_assay}, values in \code{[0, 1]}
#'     \item Frequency vector and parameters in \code{metadata(x)$itc},
#'       a list containing \code{frequencies} (numeric vector),
#'       \code{n_cycles} (integer), and \code{n_epochs} (integer)
#'   }
#' @references
#' Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
#' conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.
#'
#' Makeig, S. (1993). Auditory event-related dynamics of the EEG spectrum and
#' effects of exposure to tones. Electroencephalography and Clinical
#' Neurophysiology, 86(4), 283-293.
#' @seealso [eegERSP()], [eegMorletWavelet()], [eegPlotSpectrogram()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_erp(n_epochs = 40, n_channels = 2, sr = 250)
#' pe_itc <- eegITC(pe, frequencies = seq(5, 30, by = 5))
#' itc_data <- SummarizedExperiment::assay(pe_itc, "itc")
#' dim(itc_data)  # time x frequencies x channels
#' }
eegITC <- function(x, frequencies = NULL, n_cycles = 7,
                   assay_name = NULL, output_assay = "itc_data") {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(n_cycles) && n_cycles > 0)

  if (is.null(frequencies)) frequencies <- seq(1, 50, by = 1)
  stopifnot(is.numeric(frequencies) && length(frequencies) > 0)
  stopifnot(all(frequencies > 0))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  if (length(dim(data)) != 3) {
    stop("eegITC requires epoched (3D) data (time x channels x epochs).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_epochs <- dim(data)[3]
  n_freqs <- length(frequencies)

  # Accumulate unit-phase vectors across epochs
  # We compute: ITC = |mean(exp(i * phase))|
  # This is equivalent to |mean(complex_coeff / |complex_coeff|)|
  phase_sum_real <- array(0, dim = c(n_time, n_freqs, n_channels))
  phase_sum_imag <- array(0, dim = c(n_time, n_freqs, n_channels))

  for (ep in seq_len(n_epochs)) {
    # Get complex wavelet coefficients for this epoch
    epoch_complex <- .morlet_complex_2d(data[, , ep], sr, frequencies, n_cycles)

    # Extract unit-phase vectors: exp(i * phase) = coeff / |coeff|
    for (ch in seq_len(n_channels)) {
      for (fi in seq_len(n_freqs)) {
        coeff <- epoch_complex[, fi, ch]
        mag <- Mod(coeff)
        # Avoid division by zero
        mag[mag < .Machine$double.eps] <- 1
        unit_phase <- coeff / mag
        phase_sum_real[, fi, ch] <- phase_sum_real[, fi, ch] + Re(unit_phase)
        phase_sum_imag[, fi, ch] <- phase_sum_imag[, fi, ch] + Im(unit_phase)
      }
    }
  }

  # ITC = |mean of unit-phase vectors|
  itc_out <- sqrt((phase_sum_real / n_epochs)^2 + (phase_sum_imag / n_epochs)^2)

  # Store in metadata (different dimensions than SE assays)
  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- itc_out
  md$itc <- list(
    frequencies = frequencies,
    n_cycles = n_cycles,
    n_epochs = n_epochs
  )
  S4Vectors::metadata(x) <- md

  x
}


# --- Internal helpers ---

#' Compute DPSS (Slepian) tapers via tridiagonal eigendecomposition
#'
#' Approximates Discrete Prolate Spheroidal Sequences by constructing a
#' symmetric tridiagonal matrix whose eigenvectors corresponding to the
#' largest eigenvalues are the desired tapers.
#'
#' @param n Integer length of the data.
#' @param nw Numeric time-half-bandwidth product.
#' @param k Integer number of tapers to return.
#' @return Matrix of dimension n x k containing the DPSS tapers as columns.
#' @keywords internal
.compute_dpss <- function(n, nw, k) {
  stopifnot(n >= 2, nw > 0, k >= 1)

  w <- nw / n  # half-bandwidth

  # Construct symmetric tridiagonal matrix
  # Diagonal: ((n-1-2*i)/2)^2 * cos(2*pi*W)
  # Off-diagonal: i*(n-i)/2
  diag_vals <- ((n - 1 - 2 * (0:(n - 1))) / 2)^2 * cos(2 * pi * w)
  off_diag <- (seq_len(n - 1) * (n - seq_len(n - 1))) / 2

  # Build tridiagonal matrix
  mat <- matrix(0, nrow = n, ncol = n)
  diag(mat) <- diag_vals
  for (i in seq_len(n - 1)) {
    mat[i, i + 1] <- off_diag[i]
    mat[i + 1, i] <- off_diag[i]
  }

  # Eigendecomposition - we want the k largest eigenvalues
  # For large n, we use a subset via eigen()
  # eigen() returns values in decreasing order
  if (n <= 500) {
    eig <- eigen(mat, symmetric = TRUE)
    tapers <- eig$vectors[, seq_len(k), drop = FALSE]
  } else {
    # For large n, only compute needed eigenvalues
    # Use the full decomposition but only keep k tapers
    eig <- eigen(mat, symmetric = TRUE)
    tapers <- eig$vectors[, seq_len(k), drop = FALSE]
  }

  # Enforce consistent sign convention: first lobe positive
  for (j in seq_len(k)) {
    if (tapers[1, j] < 0) {
      tapers[, j] <- -tapers[, j]
    }
  }

  tapers
}


#' Compute Morlet wavelet power for 2D data (internal)
#'
#' Computes wavelet power (|W|^2) for a 2D matrix (time x channels) at
#' specified frequencies using Morlet wavelets with FFT convolution.
#'
#' @param data Numeric matrix (time x channels).
#' @param sr Sampling rate in Hz.
#' @param frequencies Numeric vector of frequencies.
#' @param n_cycles Number of wavelet cycles.
#' @return 3D array (time x frequencies x channels) of power values.
#' @keywords internal
.morlet_power_2d <- function(data, sr, frequencies, n_cycles) {
  n_time <- nrow(data)
  n_channels <- ncol(data)
  n_freqs <- length(frequencies)

  power_out <- array(0, dim = c(n_time, n_freqs, n_channels))

  for (fi in seq_len(n_freqs)) {
    f <- frequencies[fi]
    sigma <- n_cycles / (2 * pi * f)
    wavelet_half <- ceiling(3 * sigma * sr)
    wavelet_t <- seq(-wavelet_half, wavelet_half) / sr
    wavelet_n <- length(wavelet_t)

    wavelet <- exp(2i * pi * f * wavelet_t) * exp(-wavelet_t^2 / (2 * sigma^2))
    wavelet <- wavelet / sqrt(sum(Mod(wavelet)^2))

    n_conv <- n_time + wavelet_n - 1L
    n_fft <- 2^ceiling(log2(n_conv))
    wavelet_fft <- fft(c(wavelet, rep(0, n_fft - wavelet_n)))

    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      sig_fft <- fft(c(sig, rep(0, n_fft - n_time)))
      conv_result <- fft(sig_fft * Conj(wavelet_fft), inverse = TRUE) / n_fft

      half_wav <- wavelet_half
      start_idx <- half_wav + 1L
      end_idx <- half_wav + n_time
      trimmed <- conv_result[start_idx:end_idx]

      power_out[, fi, ch] <- Mod(trimmed)^2
    }
  }

  power_out
}


#' Compute Morlet wavelet complex coefficients for 2D data (internal)
#'
#' Computes complex wavelet coefficients for a 2D matrix (time x channels) at
#' specified frequencies using Morlet wavelets with FFT convolution. Used by
#' \code{eegITC} to extract phase information.
#'
#' @param data Numeric matrix (time x channels).
#' @param sr Sampling rate in Hz.
#' @param frequencies Numeric vector of frequencies.
#' @param n_cycles Number of wavelet cycles.
#' @return 3D complex array (time x frequencies x channels).
#' @keywords internal
.morlet_complex_2d <- function(data, sr, frequencies, n_cycles) {
  if (is.null(dim(data))) {
    data <- matrix(data, ncol = 1)
  }
  n_time <- nrow(data)
  n_channels <- ncol(data)
  n_freqs <- length(frequencies)

  complex_out <- array(0i, dim = c(n_time, n_freqs, n_channels))

  for (fi in seq_len(n_freqs)) {
    f <- frequencies[fi]
    sigma <- n_cycles / (2 * pi * f)
    wavelet_half <- ceiling(3 * sigma * sr)
    wavelet_t <- seq(-wavelet_half, wavelet_half) / sr
    wavelet_n <- length(wavelet_t)

    wavelet <- exp(2i * pi * f * wavelet_t) * exp(-wavelet_t^2 / (2 * sigma^2))
    wavelet <- wavelet / sqrt(sum(Mod(wavelet)^2))

    n_conv <- n_time + wavelet_n - 1L
    n_fft <- 2^ceiling(log2(n_conv))
    wavelet_fft <- fft(c(wavelet, rep(0, n_fft - wavelet_n)))

    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      sig_fft <- fft(c(sig, rep(0, n_fft - n_time)))
      conv_result <- fft(sig_fft * Conj(wavelet_fft), inverse = TRUE) / n_fft

      half_wav <- wavelet_half
      start_idx <- half_wav + 1L
      end_idx <- half_wav + n_time
      complex_out[, fi, ch] <- conv_result[start_idx:end_idx]
    }
  }

  complex_out
}
