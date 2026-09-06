#' ICLabel-style ICA Component Classification
#'
#' Classifies independent components (from [eegICA()]) into seven classes -
#' \code{brain}, \code{muscle}, \code{eye}, \code{heart}, \code{line_noise},
#' \code{channel_noise}, and \code{other} - and returns a calibrated probability
#' for each class per component together with the argmax label. The approach
#' follows the ICLabel framework: interpretable spatial, spectral, and temporal
#' features are extracted from each component and mapped to class probabilities
#' by a lightweight multinomial-logistic (softmax) head.
#'
#' Features per component:
#' \itemize{
#'   \item Spatial: frontal energy (eye/blink topography), focality (single
#'     channel dominance), topography kurtosis.
#'   \item Spectral: 1/f slope, high- versus low-frequency band ratio,
#'     low-frequency fraction, a genuine alpha-peak measure (8-12 Hz power above
#'     its theta/low-beta neighbours), and a line-noise power ratio and fraction.
#'   \item Temporal: lag-1 autocorrelation, activation kurtosis (spiky blink or
#'     ECG signatures), and a roughly 1 Hz periodicity measure for cardiac
#'     components.
#' }
#'
#' The softmax weights are read from
#' \code{inst/extdata/iclabel_weights.csv}; if that file is unavailable an
#' identical built-in weight table is used, so the classifier always works.
#'
#' @param x A PhysioExperiment object with ICA results (from [eegICA()]).
#' @param ica_assay Metadata name holding the component activations
#'   (default: \code{"ica_components"}).
#' @param assay_name Input assay used only to resolve channel labels
#'   (default: first assay).
#' @param line_freq Mains line frequency in Hz used for the line-noise features
#'   (default: 50). Only used by the heuristic backend.
#' @param backend Which classifier to use: \code{"heuristic"} (default) is the
#'   self-contained pure-R classifier (real features, hand-set multinomial
#'   weights); \code{"iclabel"} runs the genuine trained ICLabel CNN
#'   (Pion-Tonachini et al. 2019) by delegating to \pkg{mne-icalabel} via
#'   \pkg{reticulate} (needs a Python env with \code{mne} and \code{mne-icalabel};
#'   channel labels must match a standard 10-20/10-10 montage for the scalp-map).
#' @return A data.frame with one row per component: \code{component} (integer
#'   index), one numeric column per class (\code{brain}, \code{muscle},
#'   \code{eye}, \code{heart}, \code{line_noise}, \code{channel_noise},
#'   \code{other}) holding probabilities that sum to 1 across the classes, and
#'   \code{label} (character argmax class).
#' @references
#' Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019). ICLabel: An
#' automated electroencephalographic independent component classifier, dataset,
#' and website. NeuroImage, 198, 181-197.
#'
#' Winkler, I., Haufe, S., & Tangermann, M. (2011). Automatic classification of
#' artifactual ICA-components for artifact removal in EEG signals. Behavioral and
#' Brain Functions, 7, 30.
#' @seealso [eegICA()], [eegICLabelFlag()], [eegICAdetect()], [eegICAremove()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' pe <- eegICA(pe, n_components = 10, method = "fastica")
#' probs <- eegICLabel(pe)
#' }
eegICLabel <- function(x, ica_assay = "ica_components", assay_name = NULL,
                       line_freq = 50, backend = c("heuristic", "iclabel")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  backend <- match.arg(backend)

  ica_info <- S4Vectors::metadata(x)$ica
  if (is.null(ica_info)) {
    stop("No ICA results found. Run eegICA() first.", call. = FALSE)
  }
  sources <- S4Vectors::metadata(x)[[ica_assay]]
  if (is.null(sources)) {
    stop(sprintf("ICA activations '%s' not found in metadata. Run eegICA() first.",
                 ica_assay), call. = FALSE)
  }
  A <- ica_info$mixing
  n_components <- ncol(sources)
  sr <- samplingRate(x)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(nrow(A)))
  }

  classes <- .iclabel_classes()

  if (backend == "iclabel") {
    # genuine trained ICLabel CNN via mne-icalabel (reticulate)
    data <- SummarizedExperiment::assay(x, assay_name)
    probs <- .eeg_iclabel_trained(t(as.matrix(data)), ica_info, n_components,
                                  ch_labels, sr)
    dimnames(probs) <- list(NULL, classes)
  } else {
    # pure-R heuristic (default): real features, hand-set multinomial weights
    model <- .load_iclabel_model()
    probs <- matrix(NA_real_, n_components, length(classes),
                    dimnames = list(NULL, classes))
    for (ic in seq_len(n_components)) {
      feats <- .iclabel_features(A[, ic], sources[, ic], sr, ch_labels, line_freq)
      probs[ic, ] <- .iclabel_softmax(.iclabel_terms(feats), model, classes)
    }
  }

  out <- data.frame(component = seq_len(n_components), stringsAsFactors = FALSE)
  for (cl in classes) out[[cl]] <- probs[, cl]
  out$label <- classes[max.col(probs, ties.method = "first")]
  out
}

#' Flag Artifact ICA Components from ICLabel Probabilities
#'
#' Convenience wrapper over [eegICLabel()] that returns the indices of the
#' components most likely to be artifacts, ready to pass to [eegICAremove()].
#' A component is flagged when its probability in any of the requested artifact
#' \code{classes} exceeds \code{prob_threshold}.
#'
#' @param x A PhysioExperiment object with ICA results (from [eegICA()]).
#' @param prob_threshold Probability above which a component is flagged
#'   (default: 0.5).
#' @param classes Character vector of artifact classes to flag (default: all
#'   non-brain, non-other classes).
#' @param ica_assay Metadata name holding the component activations
#'   (default: \code{"ica_components"}).
#' @param line_freq Mains line frequency in Hz (default: 50).
#' @return An integer vector of flagged component indices (possibly empty),
#'   with an attribute \code{"labels"} giving the corresponding class of each.
#' @seealso [eegICLabel()], [eegICAremove()], [eegICAdetect()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' pe <- eegICA(pe, n_components = 10, method = "fastica")
#' bad <- eegICLabelFlag(pe, prob_threshold = 0.5)
#' pe <- eegICAremove(pe, components = bad)
#' }
eegICLabelFlag <- function(x, prob_threshold = 0.5,
                           classes = c("muscle", "eye", "heart",
                                       "line_noise", "channel_noise"),
                           ica_assay = "ica_components", line_freq = 50) {
  stopifnot(is.numeric(prob_threshold), prob_threshold >= 0, prob_threshold <= 1)
  all_classes <- .iclabel_classes()
  if (!all(classes %in% all_classes)) {
    stop("classes must be a subset of: ", paste(all_classes, collapse = ", "),
         call. = FALSE)
  }
  probs <- eegICLabel(x, ica_assay = ica_assay, line_freq = line_freq)
  sub <- as.matrix(probs[, classes, drop = FALSE])
  best <- apply(sub, 1, max)
  which_cls <- classes[max.col(sub, ties.method = "first")]
  flagged <- which(best > prob_threshold)
  out <- probs$component[flagged]
  attr(out, "labels") <- which_cls[flagged]
  out
}

# --- internal helpers -------------------------------------------------------

#' ICLabel class order
#' @keywords internal
.iclabel_classes <- function() {
  c("brain", "muscle", "eye", "heart", "line_noise", "channel_noise", "other")
}

# Package-level cache for the shipped softmax weight table.
.iclabel_model_cache <- new.env(parent = emptyenv())

#' Load the shipped ICLabel softmax weight table
#'
#' Reads the weight table from \code{inst/extdata/iclabel_weights.csv} (a
#' \code{term} column plus one column per class) and caches it. Falls back to an
#' identical built-in table when the file is unavailable.
#'
#' @return A data.frame with a \code{term} column and one column per class.
#' @keywords internal
.load_iclabel_model <- function() {
  if (!is.null(.iclabel_model_cache$model)) return(.iclabel_model_cache$model)
  path <- system.file("extdata", "iclabel_weights.csv", package = "PhysioEEG")
  m <- if (path != "") {
    tryCatch(utils::read.csv(path, stringsAsFactors = FALSE),
             error = function(e) .iclabel_default_weights())
  } else {
    .iclabel_default_weights()
  }
  .iclabel_model_cache$model <- m
  m
}

#' Built-in ICLabel softmax weights (fallback)
#'
#' Interpretable multinomial-logistic weights, identical to the shipped CSV,
#' used when the packaged file cannot be read.
#'
#' @return A data.frame with a \code{term} column and one column per class.
#' @keywords internal
.iclabel_default_weights <- function() {
  term <- c("bias", "alpha_ex", "slope_physio", "frontal", "lf_hi", "hf_pos",
            "slope_pos", "tc_spiky", "focal_hi", "topo_hi", "hr_hi",
            "line_frac_hi", "line_pos")
  data.frame(
    term          = term,
    brain         = c(0, 1.5, 1.0, -2, 0, -1.5, 0, -1.0, -1.5, 0, 0, -6, -2),
    muscle        = c(0, -2, 0, -2, 0, 3, 2, 0, 0, 0, 0, -6, -2),
    eye           = c(0, 0, 0, 4, 2, -1.5, 0, 0.5, 0, 0, 0, -6, -2),
    heart         = c(0, 0, 0, -2, 0, 0, 0, 1.5, 0, 0, 3, -6, -2),
    line_noise    = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2),
    channel_noise = c(0, 0, 0, -1.5, 0, 0, 0, 0, 3, 1.5, 0, -6, -2),
    other         = c(0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )
}

#' Welch power spectral density (internal)
#'
#' @param s Numeric signal.
#' @param sr Sampling rate in Hz.
#' @param nperseg Segment length (default: 512).
#' @return List with \code{f} (frequencies) and \code{p} (one-sided power).
#' @keywords internal
.iclabel_psd <- function(s, sr, nperseg = 512) {
  n <- length(s)
  nperseg <- min(nperseg, n)
  step <- max(1L, nperseg %/% 2L)
  w <- 0.5 - 0.5 * cos(2 * pi * (seq_len(nperseg) - 1) / (nperseg - 1))
  starts <- seq(1, n - nperseg + 1, by = step)
  P <- 0
  for (a in starts) P <- P + Mod(stats::fft(s[a:(a + nperseg - 1)] * w))^2
  P <- P / length(starts)
  half <- seq_len(nperseg %/% 2 + 1)
  list(f = (0:(nperseg - 1))[half] * sr / nperseg, p = P[half])
}

#' Band power in [lo, hi) (internal)
#' @keywords internal
.iclabel_bp <- function(ps, lo, hi) {
  m <- ps$f >= lo & ps$f < hi
  if (!any(m)) 1e-12 else sum(ps$p[m])
}

#' Excess kurtosis (internal)
#' @keywords internal
.iclabel_kurt <- function(v) {
  v <- v - mean(v)
  m2 <- mean(v^2)
  if (m2 <= 0) return(0)
  mean(v^4) / m2^2 - 3
}

#' Extract ICLabel features for one component (internal)
#'
#' @param topo Channel weights (a column of the mixing matrix).
#' @param act Component activation (time course).
#' @param sr Sampling rate in Hz.
#' @param ch_labels Channel labels.
#' @param line_freq Mains line frequency in Hz.
#' @return A named numeric vector of interpretable features.
#' @keywords internal
.iclabel_features <- function(topo, act, sr, ch_labels, line_freq = 50) {
  frontal_names <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                     "AF3", "AF4", "Fpz")
  w <- topo / sqrt(mean(topo^2) + 1e-12)              # RMS-normalised topography
  ps <- .iclabel_psd(act, sr)
  frontal <- which(ch_labels %in% frontal_names)

  focality <- max(abs(w)) / sqrt(mean(w^2) + 1e-12)
  frontal_energy <- if (length(frontal)) sum(w[frontal]^2) / sum(w^2) else 0
  topo_kurt <- .iclabel_kurt(w)

  hf_ratio <- log10((.iclabel_bp(ps, 20, 45) + 1e-12) /
                      (.iclabel_bp(ps, 1, 20) + 1e-12))
  lf_frac <- .iclabel_bp(ps, 1, 5) / (.iclabel_bp(ps, 1, 45) + 1e-12)

  # 1/f slope over 2-45 Hz (excluding the line band)
  m <- ps$f >= 2 & ps$f <= 45 & abs(ps$f - line_freq) > 2
  slope <- if (sum(m) > 3) {
    as.numeric(stats::coef(stats::lm(log10(ps$p[m] + 1e-12) ~ log10(ps$f[m])))[2])
  } else 0

  # genuine alpha peak: 8-12 Hz power above theta/low-beta neighbours
  total_b <- .iclabel_bp(ps, 1, sr / 2)
  theta_b <- .iclabel_bp(ps, 4, 7)
  alpha_b <- .iclabel_bp(ps, 8, 12)
  lbeta_b <- .iclabel_bp(ps, 13, 18)
  alpha_peak <- (alpha_b + 1e-12) / (0.5 * (theta_b + lbeta_b) + 1e-12)
  alpha_frac <- alpha_b / (total_b + 1e-12)

  # line-noise ratio and single-frequency dominance
  ln <- .iclabel_bp(ps, line_freq - 1, line_freq + 1)
  nb <- .iclabel_bp(ps, line_freq - 6, line_freq - 2) +
        .iclabel_bp(ps, line_freq + 2, line_freq + 6)
  line_ratio <- log10((ln + 1e-12) / (nb + 1e-12))
  line_frac <- ln / (total_b + 1e-12)

  a <- act - mean(act)
  s2 <- sum(a^2)
  ac1 <- if (s2 <= 0) 0 else sum(a[-1] * a[-length(a)]) / s2
  tc_kurt <- .iclabel_kurt(act)

  # roughly 1 Hz periodicity of the rectified activation (cardiac)
  env <- abs(act - mean(act))
  lag_lo <- max(1L, round(0.6 * sr))
  lag_hi <- min(length(env) - 1L, round(1.2 * sr))
  hr_period <- if (lag_hi > lag_lo) {
    e <- env - mean(env)
    de <- sum(e^2) + 1e-12
    max(vapply(lag_lo:lag_hi, function(L)
      sum(e[-(1:L)] * e[1:(length(e) - L)]) / de, numeric(1)))
  } else 0

  c(focality = focality, frontal = frontal_energy, topo_kurt = topo_kurt,
    hf_ratio = hf_ratio, lf_frac = lf_frac, slope = slope,
    alpha = alpha_peak, alpha_frac = alpha_frac, line = line_ratio,
    line_frac = line_frac, ac1 = ac1, tc_kurt = tc_kurt, hr = hr_period)
}

#' Transform raw features to the softmax design terms (internal)
#'
#' @param f Named feature vector from [.iclabel_features()].
#' @return A named numeric vector of design terms matching the weight table.
#' @keywords internal
.iclabel_terms <- function(f) {
  c(
    bias         = 1,
    # genuine alpha peak, gated so a noise-floor ratio cannot fake it
    alpha_ex     = max(0, f[["alpha"]] - 1) * (f[["alpha_frac"]] > 0.03),
    slope_physio = max(0, 1 - abs(f[["slope"]] + 1.3) / 1.5),
    frontal      = f[["frontal"]],
    lf_hi        = max(0, f[["lf_frac"]] - 0.4),
    hf_pos       = max(0, f[["hf_ratio"]]),
    slope_pos    = max(0, f[["slope"]] + 0.3),
    tc_spiky     = max(0, f[["tc_kurt"]] - 2),
    focal_hi     = max(0, f[["focality"]] - 5),
    topo_hi      = max(0, f[["topo_kurt"]] - 3),
    hr_hi        = max(0, f[["hr"]] - 0.25),
    line_frac_hi = max(0, f[["line_frac"]] - 0.25),
    line_pos     = max(0, f[["line"]])
  )
}

#' Softmax class probabilities from design terms and weight table (internal)
#'
#' @param terms Named design-term vector from [.iclabel_terms()].
#' @param model Weight table from [.load_iclabel_model()].
#' @param classes Class order.
#' @return Named numeric vector of probabilities summing to 1.
#' @keywords internal
.iclabel_softmax <- function(terms, model, classes) {
  tv <- terms[model$term]                              # align to model row order
  scores <- vapply(classes, function(cl) sum(model[[cl]] * tv), numeric(1))
  e <- exp(scores - max(scores))
  e / sum(e)
}
