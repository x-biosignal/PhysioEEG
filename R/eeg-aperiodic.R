# EEG aperiodic (1/f) spectral parameterization.
#
# The FOOOF/specparam engine lives in PhysioAnalysis (a PhysioCore-only sibling);
# this is the EEG-facing entry point. It runs specparam per channel and returns
# the aperiodic exponent/offset/(knee) as a per-channel table ready for a topomap
# or a qEEG biomarker, alongside the separated oscillatory peaks.

#' Aperiodic (1/f) spectral parameterization of EEG
#'
#' Separates each channel's power spectrum into an **aperiodic** (1/f) component
#' and **periodic** (oscillatory) peaks, following the specparam / FOOOF model
#' (Donoghue et al. 2020). The aperiodic exponent is a widely used index of the
#' excitation/inhibition balance and cortical state. Delegates the fit to
#' `PhysioAnalysis::specparam()`.
#'
#' @param pe A `PhysioExperiment`.
#' @param freq_range Frequency range to fit, in Hz (default `c(1, 45)`).
#' @param aperiodic_mode `"fixed"` (offset + exponent) or `"knee"` (offset +
#'   knee + exponent), the latter for spectra with a bend in log-log space.
#' @param max_n_peaks Maximum number of oscillatory peaks per channel (default 6).
#' @param peak_width_limits Min/max peak width in Hz (default `c(1, 12)`).
#' @param min_peak_height Minimum peak height above the aperiodic fit (default 0.05).
#' @param peak_threshold Peak detection threshold in SD of the flattened spectrum
#'   (default 2).
#' @param assay_name Assay to use (default: the object's default assay).
#' @return An `eeg_aperiodic` object: a list with `aperiodic` (per-channel data
#'   frame: `channel`, `exponent`, `offset`, optionally `knee`, `r_squared`,
#'   `error`), `peaks` (per-channel `CF`/`PW`/`BW`), `exponent` (a named
#'   per-channel vector, e.g. for [eegPlotTopomap()]), and the underlying
#'   `specparam_result`.
#' @references Donoghue et al. 2020, Nat Neurosci (specparam / FOOOF).
#' @seealso [eegQEEG()], [eegComplexity()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 2500, n_channels = 8, sr = 250)
#' ap <- eegAperiodic(pe, freq_range = c(2, 40))
#' ap$aperiodic          # per-channel exponent / offset
#' }
eegAperiodic <- function(pe, freq_range = c(1, 45),
                         aperiodic_mode = c("fixed", "knee"),
                         max_n_peaks = 6L, peak_width_limits = c(1, 12),
                         min_peak_height = 0.05, peak_threshold = 2,
                         assay_name = NULL) {
  stopifnot(inherits(pe, "PhysioExperiment"))
  aperiodic_mode <- match.arg(aperiodic_mode)
  if (!requireNamespace("PhysioAnalysis", quietly = TRUE)) {
    stop("The 'PhysioAnalysis' package provides the specparam/FOOOF engine used ",
         "by eegAperiodic(). Install it from the x-biosignal r-universe: ",
         "install.packages('PhysioAnalysis', repos = 'https://x-biosignal.r-universe.dev')",
         call. = FALSE)
  }
  res <- PhysioAnalysis::specparam(
    pe, freq_range = freq_range, aperiodic_mode = aperiodic_mode,
    max_n_peaks = max_n_peaks, peak_width_limits = peak_width_limits,
    min_peak_height = min_peak_height, peak_threshold = peak_threshold,
    assay_name = assay_name)

  aperiodic <- merge(res$aperiodic, res$fit, by = "channel", sort = FALSE)
  exponent <- stats::setNames(aperiodic$exponent, aperiodic$channel)
  structure(list(aperiodic = aperiodic, peaks = res$peaks, exponent = exponent,
                 freq_range = freq_range, aperiodic_mode = aperiodic_mode,
                 specparam = res),
            class = "eeg_aperiodic")
}

#' @export
print.eeg_aperiodic <- function(x, ...) {
  cat(sprintf("<eeg_aperiodic> %d channel(s), mode = %s, %g-%g Hz\n",
              nrow(x$aperiodic), x$aperiodic_mode, x$freq_range[1], x$freq_range[2]))
  cat(sprintf("  aperiodic exponent: %.3f - %.3f (mean %.3f)\n",
              min(x$exponent, na.rm = TRUE), max(x$exponent, na.rm = TRUE),
              mean(x$exponent, na.rm = TRUE)))
  cat(sprintf("  oscillatory peaks: %d total; fit R-squared mean %.3f\n",
              nrow(x$peaks), mean(x$aperiodic$r_squared, na.rm = TRUE)))
  invisible(x)
}
