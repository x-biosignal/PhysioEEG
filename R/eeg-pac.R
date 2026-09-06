# EEG phase-amplitude coupling (cross-frequency coupling).
#
# Within-EEG PAC (phase of a slow band modulating the amplitude of a fast band),
# either within one channel or across two channels. The PAC estimators (Tort MI,
# Canolty MVL, Ozkurt, PLV) + surrogate significance live in PhysioCrossModal
# (already a PhysioEEG Suggests); this is the EEG-facing entry point that pulls
# the channel signal(s) and delegates. PhysioCrossModal's PAC is the cross-*modal*
# use (EEG phase <-> EMG amplitude, etc.); this covers the within-EEG case.

# channel signal as a numeric vector (2D column, or 3D epochs concatenated)
.eeg_pac_channel <- function(pe, ch, assay_name) {
  data <- SummarizedExperiment::assay(pe, assay_name)
  cd <- SummarizedExperiment::colData(pe)
  if (is.character(ch)) {
    labs <- if ("label" %in% colnames(cd)) as.character(cd$label) else NULL
    idx <- match(ch, labs)
    if (is.na(idx)) stop("Channel '", ch, "' not found.", call. = FALSE)
    ch <- idx
  }
  if (length(dim(data)) == 3L) as.vector(data[, ch, ]) else data[, ch]
}

.need_crossmodal <- function(fn) {
  if (!requireNamespace("PhysioCrossModal", quietly = TRUE)) {
    stop("The 'PhysioCrossModal' package provides the PAC engine used by ", fn,
         "(). Install it from the x-biosignal r-universe: ",
         "install.packages('PhysioCrossModal', repos = 'https://x-biosignal.r-universe.dev')",
         call. = FALSE)
  }
}

#' Phase-amplitude coupling (cross-frequency coupling) for EEG
#'
#' Measures how the phase of a slow oscillation modulates the amplitude of a fast
#' oscillation — within a single channel (`amp_channel = NULL`) or across two
#' channels. Delegates the estimator to `PhysioCrossModal::phaseAmplitudeCoupling()`.
#'
#' @param pe A `PhysioExperiment` (continuous data; epoched data is concatenated).
#' @param phase_channel Channel (index or label) supplying the modulating phase.
#' @param amp_channel Channel supplying the modulated amplitude; `NULL` (default)
#'   uses `phase_channel` (within-channel PAC).
#' @param phase_band Phase frequency band in Hz (default theta `c(4, 8)`).
#' @param amp_band Amplitude frequency band in Hz (default gamma `c(30, 80)`).
#' @param method `"tort"` (modulation index), `"canolty"` (mean vector length),
#'   `"ozkurt"`, or `"plv"`.
#' @param n_bins Number of phase bins for the Tort modulation index (default 18).
#' @param assay_name Assay to use (default: the object's default assay).
#' @return The PAC result from `PhysioCrossModal::phaseAmplitudeCoupling()` (a
#'   list with `pac` and, for `"tort"`, the phase-amplitude `distribution`).
#' @references Tort et al. 2010; Canolty et al. 2006.
#' @seealso [eegComodulogram()], [eegConnectivityMatrix()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 250)
#' eegPAC(pe, phase_channel = 1, phase_band = c(4, 8), amp_band = c(30, 80))
#' }
eegPAC <- function(pe, phase_channel = 1L, amp_channel = NULL,
                   phase_band = c(4, 8), amp_band = c(30, 80),
                   method = c("tort", "canolty", "ozkurt", "plv"),
                   n_bins = 18L, assay_name = NULL) {
  stopifnot(inherits(pe, "PhysioExperiment"))
  method <- match.arg(method)
  .need_crossmodal("eegPAC")
  if (is.null(assay_name)) assay_name <- defaultAssay(pe)
  sr <- samplingRate(pe)
  xph <- .eeg_pac_channel(pe, phase_channel, assay_name)
  yamp <- if (is.null(amp_channel)) NULL else .eeg_pac_channel(pe, amp_channel, assay_name)
  PhysioCrossModal::phaseAmplitudeCoupling(
    x = xph, y = yamp, sr = sr, phase_band = phase_band, amp_band = amp_band,
    method = method, n_bins = as.integer(n_bins))
}

#' Comodulogram (phase-frequency x amplitude-frequency PAC) for EEG
#'
#' Computes phase-amplitude coupling over a grid of phase and amplitude
#' frequencies for one EEG channel, yielding a comodulogram whose peak locates
#' the dominant coupling. Delegates to `PhysioCrossModal::comodulogram()`.
#'
#' @param pe A `PhysioExperiment`.
#' @param channel Channel (index or label) to analyse.
#' @param phase_freqs Phase centre frequencies in Hz (default `seq(2, 14, 2)`).
#' @param amp_freqs Amplitude centre frequencies in Hz (default `seq(20, 100, 10)`).
#' @param method PAC estimator (see [eegPAC()]).
#' @param phase_bw,amp_bw Half-bandwidths in Hz (defaults 2 and 10).
#' @param n_bins Phase bins for the Tort modulation index (default 18).
#' @param assay_name Assay to use (default: the object's default assay).
#' @return The comodulogram list from `PhysioCrossModal::comodulogram()` (a
#'   `matrix` of phase x amplitude frequencies, plus its `peak`).
#' @seealso [eegPAC()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 7500, n_channels = 2, sr = 250)
#' cm <- eegComodulogram(pe, channel = 1)
#' cm$peak
#' }
eegComodulogram <- function(pe, channel = 1L,
                            phase_freqs = seq(2, 14, by = 2),
                            amp_freqs = seq(20, 100, by = 10),
                            method = c("tort", "canolty", "ozkurt", "plv"),
                            phase_bw = 2, amp_bw = 10, n_bins = 18L,
                            assay_name = NULL) {
  stopifnot(inherits(pe, "PhysioExperiment"))
  method <- match.arg(method)
  .need_crossmodal("eegComodulogram")
  if (is.null(assay_name)) assay_name <- defaultAssay(pe)
  sr <- samplingRate(pe)
  x <- .eeg_pac_channel(pe, channel, assay_name)
  PhysioCrossModal::comodulogram(
    x = x, sr = sr, phase_freqs = phase_freqs, amp_freqs = amp_freqs,
    method = method, phase_bw = phase_bw, amp_bw = amp_bw,
    n_bins = as.integer(n_bins))
}
