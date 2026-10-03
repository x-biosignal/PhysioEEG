#' @keywords internal
#' @importFrom stats cor cov fft prcomp rnorm runif sd var median quantile
#' @importFrom utils head
#' @importFrom PhysioExperiment PhysioExperiment defaultAssay samplingRate
#' @examples
#' # Simulate EEG, band-pass filter, and summarise signal complexity
#' pe <- make_eeg(n_time = 2500, n_channels = 8, sr = 250)
#' pe <- eegFilter(pe, lowcut = 1, highcut = 40)
#' eegComplexity(pe, measures = c("hjorth_mobility", "spectral_entropy"))
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}
