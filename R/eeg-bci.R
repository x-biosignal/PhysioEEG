#' Common Spatial Pattern (CSP) Analysis
#'
#' Computes Common Spatial Pattern filters for two-class EEG discrimination.
#' CSP maximizes the variance ratio between two conditions, making it a
#' standard spatial filtering technique for motor imagery BCI.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x trials).
#' @param labels Character or factor vector of class labels, one per trial.
#'   Must contain exactly two unique classes.
#' @param n_filters Number of CSP filter pairs to retain (default: 3).
#'   The total number of spatial filters will be \code{2 * n_filters}.
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name for CSP features (default: \code{"csp"}).
#' @return Modified PhysioExperiment with CSP log-variance features in
#'   \code{output_assay} (a matrix of trials x \code{2 * n_filters}) and
#'   CSP filter information stored in \code{metadata(x)$csp} as a list
#'   containing \code{filters} (spatial filter matrix), \code{eigenvalues}
#'   (selected eigenvalues), and \code{classes} (unique class labels).
#' @references
#' Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
#' single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.
#' @seealso [eegBCIfeatures()], [eegBCIclassify()], [eegMotorImagery()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
#' labels <- metadata(pe)$labels
#' result <- eegCSP(pe, labels = labels, n_filters = 3)
#' csp_features <- SummarizedExperiment::assay(result, "csp")
#' }
eegCSP <- function(x, labels, n_filters = 3, assay_name = NULL,
                   output_assay = "csp") {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  # Validate 3D data
  if (length(dim(data)) != 3) {
    stop("eegCSP requires epoched (3D) data (time x channels x trials).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_trials <- dim(data)[3]

  # Validate labels
  if (length(labels) != n_trials) {
    stop(sprintf("Length of labels (%d) must match number of trials (%d).",
                 length(labels), n_trials), call. = FALSE)
  }

  classes <- unique(labels)
  if (length(classes) != 2) {
    stop(sprintf("CSP requires exactly 2 classes, found %d.", length(classes)),
         call. = FALSE)
  }

  if (n_filters > floor(n_channels / 2)) {
    n_filters <- floor(n_channels / 2)
    warning(sprintf("n_filters reduced to %d (half the number of channels).",
                    n_filters), call. = FALSE)
  }

  # Compute class-specific normalized covariance matrices
  cov_class1 <- matrix(0, n_channels, n_channels)
  n1 <- sum(labels == classes[1])
  for (trial in which(labels == classes[1])) {
    trial_data <- data[, , trial]
    C <- cov(trial_data)
    tr <- sum(diag(C))
    if (tr > 0) {
      cov_class1 <- cov_class1 + C / tr
    }
  }
  cov_class1 <- cov_class1 / n1

  cov_class2 <- matrix(0, n_channels, n_channels)
  n2 <- sum(labels == classes[2])
  for (trial in which(labels == classes[2])) {
    trial_data <- data[, , trial]
    C <- cov(trial_data)
    tr <- sum(diag(C))
    if (tr > 0) {
      cov_class2 <- cov_class2 + C / tr
    }
  }
  cov_class2 <- cov_class2 / n2

  # Composite covariance
  Cc <- cov_class1 + cov_class2

  # Eigen decompose composite covariance for whitening
  eig_c <- eigen(Cc, symmetric = TRUE)
  eig_vals <- eig_c$values
  eig_vals[eig_vals < 1e-10] <- 1e-10
  # Whitening matrix P
  P <- diag(eig_vals^(-0.5)) %*% t(eig_c$vectors)

  # Apply whitening to class 1 covariance
  S1 <- P %*% cov_class1 %*% t(P)

  # Eigen decompose S1
  eig_s1 <- eigen(S1, symmetric = TRUE)

  # Select top n_filters (highest eigenvalues for class 1) and
  # bottom n_filters (lowest eigenvalues = highest for class 2)
  idx_top <- seq_len(n_filters)
  idx_bottom <- seq(n_channels - n_filters + 1, n_channels)
  selected <- c(idx_top, idx_bottom)

  # Spatial filters: project back through whitening
  W <- t(eig_s1$vectors[, selected]) %*% P  # (2*n_filters x n_channels)

  # Extract CSP features: log-variance of spatially filtered data
  n_filters_total <- 2 * n_filters
  features <- matrix(NA_real_, nrow = n_trials, ncol = n_filters_total)

  for (trial in seq_len(n_trials)) {
    trial_data <- data[, , trial]  # time x channels
    # Apply spatial filters: (2*n_filters x n_channels) %*% t(time x channels)
    # = (2*n_filters x time)
    filtered <- W %*% t(trial_data)
    # Log-variance of each filter output
    vars <- apply(filtered, 1, var)
    vars[vars < 1e-20] <- 1e-20
    features[trial, ] <- log(vars / sum(vars))
  }

  # Store CSP information and features in metadata
  # (Features have different dimensions than raw assays, so they go in metadata)
  md <- S4Vectors::metadata(x)
  md$csp <- list(
    filters = W,
    eigenvalues = eig_s1$values[selected],
    classes = classes
  )
  md$csp_features <- features
  S4Vectors::metadata(x) <- md

  x
}


#' SSVEP Frequency Detection
#'
#' Detects Steady-State Visual Evoked Potentials using canonical correlation
#' analysis (CCA) or filter-bank CCA (FBCCA). For each candidate stimulus
#' frequency, a CCA is computed between the EEG data and sinusoidal reference
#' signals at the frequency and its harmonics.
#'
#' @param x A PhysioExperiment object with EEG data. If 3D (time x channels x
#'   trials), data is averaged across trials before analysis.
#' @param frequencies Numeric vector of target stimulus frequencies in Hz.
#' @param n_harmonics Number of harmonics to include in reference signals
#'   (default: 3).
#' @param method Detection method: \code{"cca"} (canonical correlation analysis)
#'   or \code{"fbcca"} (filter-bank CCA).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with columns: \code{frequency} (numeric target
#'   frequency in Hz), \code{correlation} (numeric CCA correlation),
#'   \code{snr} (numeric signal-to-noise ratio), and
#'   \code{predicted_class} (numeric predicted stimulus frequency).
#' @references
#' Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
#' single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.
#'
#' Norcia, A. M., et al. (2015). The steady-state visual evoked potential in
#' vision research: a review. Journal of Vision, 15(6), 4.
#' @seealso [eegCSP()], [eegBCIfeatures()], [eegBCIclassify()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
#' result <- eegSSVEP(pe, frequencies = c(10, 12, 15), method = "cca")
#' }
eegSSVEP <- function(x, frequencies, n_harmonics = 3,
                     method = c("cca", "fbcca"),
                     assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(frequencies) && length(frequencies) > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # If 3D, average across trials to get 2D (time x channels)
  if (length(dim(data)) == 3) {
    data <- apply(data, c(1, 2), mean)
  }

  n_time <- nrow(data)
  n_channels <- ncol(data)
  t_vec <- (seq_len(n_time) - 1) / sr

  # Internal function: compute CCA correlation between X and Y
  .cca_correlation <- function(X, Y) {
    # Center data
    X <- scale(X, center = TRUE, scale = FALSE)
    Y <- scale(Y, center = TRUE, scale = FALSE)

    n <- nrow(X)

    # Covariance matrices with regularization
    Cxx <- crossprod(X) / (n - 1) + diag(1e-6, ncol(X))
    Cyy <- crossprod(Y) / (n - 1) + diag(1e-6, ncol(Y))
    Cxy <- crossprod(X, Y) / (n - 1)

    # SVD-based CCA: Cxx^(-1/2) %*% Cxy %*% Cyy^(-1/2)
    eig_xx <- eigen(Cxx, symmetric = TRUE)
    eig_xx$values[eig_xx$values < 1e-10] <- 1e-10
    Cxx_inv_sqrt <- eig_xx$vectors %*% diag(eig_xx$values^(-0.5)) %*% t(eig_xx$vectors)

    eig_yy <- eigen(Cyy, symmetric = TRUE)
    eig_yy$values[eig_yy$values < 1e-10] <- 1e-10
    Cyy_inv_sqrt <- eig_yy$vectors %*% diag(eig_yy$values^(-0.5)) %*% t(eig_yy$vectors)

    M <- Cxx_inv_sqrt %*% Cxy %*% Cyy_inv_sqrt
    sv <- svd(M)
    max(sv$d)
  }

  if (method == "cca") {
    correlations <- numeric(length(frequencies))

    for (fi in seq_along(frequencies)) {
      freq <- frequencies[fi]
      # Build reference matrix: sin and cos at each harmonic
      Y <- matrix(NA_real_, nrow = n_time, ncol = 2 * n_harmonics)
      for (h in seq_len(n_harmonics)) {
        Y[, 2 * (h - 1) + 1] <- sin(2 * pi * h * freq * t_vec)
        Y[, 2 * (h - 1) + 2] <- cos(2 * pi * h * freq * t_vec)
      }
      correlations[fi] <- .cca_correlation(data, Y)
    }

  } else if (method == "fbcca") {
    # Filter-bank CCA: apply sub-band filters, compute CCA per band, weight
    min_freq <- max(1, min(frequencies) - 2)
    max_freq <- sr / 2 - 1
    band_width <- 8
    band_starts <- seq(min_freq, max_freq - band_width, by = band_width)
    if (length(band_starts) == 0) band_starts <- min_freq
    n_bands <- length(band_starts)

    correlations <- numeric(length(frequencies))

    for (fi in seq_along(frequencies)) {
      freq <- frequencies[fi]
      # Build reference matrix
      Y <- matrix(NA_real_, nrow = n_time, ncol = 2 * n_harmonics)
      for (h in seq_len(n_harmonics)) {
        Y[, 2 * (h - 1) + 1] <- sin(2 * pi * h * freq * t_vec)
        Y[, 2 * (h - 1) + 2] <- cos(2 * pi * h * freq * t_vec)
      }

      weighted_corr <- 0
      for (bi in seq_len(n_bands)) {
        low <- band_starts[bi]
        high <- min(low + band_width, sr / 2 - 1)
        # Filter each channel
        filtered_data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
        for (ch in seq_len(n_channels)) {
          filtered_data[, ch] <- .fir_bandpass(data[, ch], sr, low, high)
        }
        corr <- .cca_correlation(filtered_data, Y)
        weight <- bi^(-1.25) + 0.25
        weighted_corr <- weighted_corr + weight * corr
      }
      correlations[fi] <- weighted_corr / n_bands
    }
  }

  # Compute SNR: correlation^2 / (1 - correlation^2) * df
  # df approximated as n_time - 1
  df <- n_time - 1
  corr_sq <- pmin(correlations^2, 0.9999)  # avoid division by zero

  snr <- corr_sq / (1 - corr_sq) * df

  # Predicted class: frequency with highest correlation
  best_idx <- which.max(correlations)
  predicted_class <- rep(frequencies[best_idx], length(frequencies))

  data.frame(
    frequency = frequencies,
    correlation = correlations,
    snr = snr,
    predicted_class = predicted_class,
    stringsAsFactors = FALSE
  )
}


#' Motor Imagery ERD/ERS Computation
#'
#' Computes event-related desynchronization (ERD) and event-related
#' synchronization (ERS) for motor imagery BCI analysis. ERD/ERS is expressed
#' as a percentage change from baseline power in specified frequency bands.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x trials).
#' @param bands Named list of frequency bands, each a numeric vector
#'   \code{c(low, high)}. Default: \code{list(mu = c(8, 13), beta = c(13, 30))}.
#' @param baseline_fraction Fraction of each trial to use as baseline
#'   (default: 0.25, i.e., the first 25 percent of the trial).
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: \code{"erd_ers"}).
#' @return Modified PhysioExperiment with ERD/ERS percentage values in
#'   \code{output_assay}. The output is a matrix of
#'   \code{n_trials} rows x \code{n_channels * n_bands} columns.
#'   Sets \code{metadata(x)$erd_ers_bands} with the band definitions used.
#' @references
#' Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
#' single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.
#' @seealso [eegCSP()], [eegBCIfeatures()], [eegBCIclassify()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
#' result <- eegMotorImagery(pe)
#' erd_data <- SummarizedExperiment::assay(result, "erd_ers")
#' }
eegMotorImagery <- function(x, bands = NULL, baseline_fraction = 0.25,
                            assay_name = NULL, output_assay = "erd_ers") {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Validate 3D data
  if (length(dim(data)) != 3) {
    stop("eegMotorImagery requires epoched (3D) data (time x channels x trials).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_trials <- dim(data)[3]

  # Default frequency bands
  if (is.null(bands)) {
    bands <- list(mu = c(8, 13), beta = c(13, 30))
  }
  stopifnot(is.list(bands))
  n_bands <- length(bands)
  band_names <- names(bands)
  if (is.null(band_names)) band_names <- paste0("band", seq_len(n_bands))

  # Baseline sample count
  baseline_samples <- max(1L, as.integer(floor(baseline_fraction * n_time)))

  # Output: n_trials x (n_channels * n_bands)
  erd_ers <- matrix(NA_real_, nrow = n_trials, ncol = n_channels * n_bands)
  col_labels <- character(n_channels * n_bands)

  # Get channel labels
  col_data <- SummarizedExperiment::colData(x)
  ch_labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }

  col_idx <- 0
  for (bi in seq_len(n_bands)) {
    low <- bands[[bi]][1]
    high <- bands[[bi]][2]

    for (ch in seq_len(n_channels)) {
      col_idx <- col_idx + 1
      col_labels[col_idx] <- paste0(ch_labels[ch], "_", band_names[bi])

      for (trial in seq_len(n_trials)) {
        sig <- data[, ch, trial]

        # Bandpass filter using FIR filter
        filtered <- .fir_bandpass(sig, sr, low, high)

        # Instantaneous power (squared amplitude)
        power <- filtered^2

        # Baseline power (first baseline_fraction of trial)
        baseline_power <- mean(power[seq_len(baseline_samples)])
        if (baseline_power < 1e-20) baseline_power <- 1e-20

        # Task power (remaining portion)
        task_power <- mean(power[(baseline_samples + 1):n_time])

        # ERD/ERS as percentage change
        erd_ers[trial, col_idx] <- (task_power - baseline_power) / baseline_power * 100
      }
    }
  }

  colnames(erd_ers) <- col_labels

  # Store band definitions and ERD/ERS features in metadata
  # (Features have different dimensions than raw assays, so they go in metadata)
  md <- S4Vectors::metadata(x)
  md$erd_ers_bands <- bands
  md$erd_ers <- erd_ers
  S4Vectors::metadata(x) <- md

  x
}


#' BCI Feature Extraction
#'
#' Extracts features from epoched EEG data for Brain-Computer Interface
#' classification. Supports band power, CSP, and Riemannian geometry methods.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data
#'   (time x channels x trials).
#' @param method Feature extraction method: \code{"bandpower"} (log band power),
#'   \code{"csp"} (Common Spatial Pattern log-variance), or \code{"riemannian"}
#'   (tangent space projection of covariance matrices).
#' @param labels Character or factor vector of class labels (required for
#'   \code{"csp"} method). One label per trial.
#' @param bands Named list of frequency bands for \code{"bandpower"} method.
#'   Default: \code{list(mu = c(8, 13), beta = c(13, 30))}.
#' @param assay_name Input assay name (default: first assay).
#' @return A numeric matrix with \code{n_trials} rows and feature columns.
#'   Number of columns depends on method:
#'   \itemize{
#'     \item \code{"bandpower"}: \code{n_channels * n_bands}
#'     \item \code{"csp"}: \code{2 * n_filters} (default: 6)
#'     \item \code{"riemannian"}: \code{n_channels * (n_channels + 1) / 2}
#'   }
#' @references
#' Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
#' single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.
#' @seealso [eegCSP()], [eegMotorImagery()], [eegBCIclassify()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
#' features <- eegBCIfeatures(pe, method = "bandpower")
#' }
eegBCIfeatures <- function(x, method = c("bandpower", "csp", "riemannian"),
                           labels = NULL, bands = NULL,
                           assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  # Validate 3D data
  if (length(dim(data)) != 3) {
    stop("eegBCIfeatures requires epoched (3D) data (time x channels x trials).",
         call. = FALSE)
  }

  n_time <- dim(data)[1]
  n_channels <- dim(data)[2]
  n_trials <- dim(data)[3]

  if (method == "bandpower") {
    # Default frequency bands
    if (is.null(bands)) {
      bands <- list(mu = c(8, 13), beta = c(13, 30))
    }
    n_bands <- length(bands)

    features <- matrix(NA_real_, nrow = n_trials, ncol = n_channels * n_bands)

    for (trial in seq_len(n_trials)) {
      col_idx <- 0
      for (bi in seq_len(n_bands)) {
        low <- bands[[bi]][1]
        high <- bands[[bi]][2]
        for (ch in seq_len(n_channels)) {
          col_idx <- col_idx + 1
          sig <- data[, ch, trial]
          filtered <- .fir_bandpass(sig, sr, low, high)
          bp <- var(filtered)
          if (bp < 1e-20) bp <- 1e-20
          features[trial, col_idx] <- log(bp)
        }
      }
    }

  } else if (method == "csp") {
    if (is.null(labels)) {
      stop("labels are required for CSP feature extraction.", call. = FALSE)
    }
    # Use eegCSP to compute filters and features (stored in metadata)
    result <- eegCSP(x, labels = labels, n_filters = 3, assay_name = assay_name)
    features <- S4Vectors::metadata(result)$csp_features

  } else if (method == "riemannian") {
    # Tangent space projection of trial covariance matrices
    n_features <- n_channels * (n_channels + 1) / 2

    # Compute covariance matrices for all trials
    trial_covs <- vector("list", n_trials)
    for (trial in seq_len(n_trials)) {
      trial_data <- data[, , trial]
      C <- cov(trial_data)
      # Regularize
      C <- C + diag(1e-6, n_channels)
      trial_covs[[trial]] <- C
    }

    # Reference point: Riemannian (Frechet) geometric mean
    C_ref <- .geometric_mean_spd(trial_covs)

    # C_ref^(-1/2) for tangent space projection
    eig_ref <- eigen(C_ref, symmetric = TRUE)
    eig_ref$values[eig_ref$values < 1e-10] <- 1e-10
    C_ref_inv_sqrt <- eig_ref$vectors %*% diag(eig_ref$values^(-0.5)) %*%
      t(eig_ref$vectors)

    features <- matrix(NA_real_, nrow = n_trials, ncol = as.integer(n_features))

    for (trial in seq_len(n_trials)) {
      # Project: C_ref^(-1/2) %*% C_trial %*% C_ref^(-1/2)
      S <- C_ref_inv_sqrt %*% trial_covs[[trial]] %*% C_ref_inv_sqrt

      # Matrix logarithm via eigendecomposition
      eig_s <- eigen(S, symmetric = TRUE)
      eig_s$values[eig_s$values < 1e-10] <- 1e-10
      log_S <- eig_s$vectors %*% diag(log(eig_s$values)) %*% t(eig_s$vectors)

      # Extract upper triangle (including diagonal)
      features[trial, ] <- log_S[upper.tri(log_S, diag = TRUE)]
    }
  }

  features
}


#' BCI Classification
#'
#' Classifies BCI features using Linear Discriminant Analysis (LDA) or
#' shrinkage LDA. Implements Fisher's LDA with optional Ledoit-Wolf
#' shrinkage regularization for robust classification with high-dimensional
#' or small-sample data. Optionally performs k-fold cross-validation.
#'
#' @param x A PhysioExperiment object with epoched (3D) EEG data.
#' @param features Optional pre-computed feature matrix (n_trials x n_features).
#'   If \code{NULL}, features are extracted using \code{\link{eegBCIfeatures}}
#'   with the \code{"bandpower"} method.
#' @param labels Character or factor vector of class labels, one per trial.
#'   Must contain exactly two unique classes.
#' @param method Classification method: \code{"lda"} (Fisher's LDA) or
#'   \code{"shrinkage_lda"} (LDA with Ledoit-Wolf shrinkage).
#' @param cv_folds Number of cross-validation folds (default: \code{NULL},
#'   meaning no cross-validation). If set (e.g., 5), performs k-fold CV and
#'   reports out-of-fold predictions. The CV accuracy is stored as
#'   \code{attr(result, "cv_accuracy")}.
#' @param assay_name Input assay name used when extracting features
#'   (default: first assay).
#' @return A data.frame with columns: \code{trial}, \code{predicted_class},
#'   \code{confidence}, and \code{true_class}. When \code{cv_folds} is not
#'   \code{NULL}, predictions are out-of-fold and \code{attr(result,
#'   "cv_accuracy")} contains the cross-validated accuracy. The trained LDA
#'   model (on all data) is stored in \code{metadata(x)$bci_model},
#'   containing \code{weights}, \code{threshold}, \code{classes},
#'   \code{method}, and \code{class_means}.
#' @references
#' Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
#' single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.
#' @seealso [eegCSP()], [eegBCIfeatures()], [eegMotorImagery()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
#' labels <- metadata(pe)$labels
#' features <- eegBCIfeatures(pe, method = "bandpower")
#' result <- eegBCIclassify(pe, features = features, labels = labels, method = "lda")
#'
#' # With 5-fold cross-validation
#' result_cv <- eegBCIclassify(pe, features = features, labels = labels,
#'                             method = "lda", cv_folds = 5)
#' attr(result_cv, "cv_accuracy")
#' }
eegBCIclassify <- function(x, features = NULL, labels,
                           method = c("lda", "shrinkage_lda"),
                           cv_folds = NULL,
                           assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  # Extract features if not provided
  if (is.null(features)) {
    features <- eegBCIfeatures(x, method = "bandpower", assay_name = assay_name)
  }

  stopifnot(is.matrix(features))
  n_trials <- nrow(features)
  n_features <- ncol(features)

  # Validate labels
  if (length(labels) != n_trials) {
    stop(sprintf("Length of labels (%d) must match number of trials (%d).",
                 length(labels), n_trials), call. = FALSE)
  }

  classes <- unique(labels)
  if (length(classes) != 2) {
    stop(sprintf("Classification requires exactly 2 classes, found %d.",
                 length(classes)), call. = FALSE)
  }

  if (!is.null(cv_folds)) {
    # --- k-fold cross-validation ---
    stopifnot(is.numeric(cv_folds) && cv_folds >= 2)
    cv_folds <- as.integer(cv_folds)

    n <- n_trials
    fold_ids <- rep(seq_len(cv_folds), length.out = n)
    fold_ids <- sample(fold_ids)  # shuffle

    cv_predictions <- character(n)
    cv_confidence <- numeric(n)

    for (fold in seq_len(cv_folds)) {
      test_idx <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)

      train_features <- features[train_idx, , drop = FALSE]
      train_labels <- labels[train_idx]
      test_features <- features[test_idx, , drop = FALSE]

      # Train on training set
      model <- .train_lda(train_features, train_labels, method)

      # Predict on test set
      preds <- .predict_lda(model, test_features)
      cv_predictions[test_idx] <- preds$class
      cv_confidence[test_idx] <- preds$confidence
    }

    result <- data.frame(
      trial = seq_len(n),
      predicted_class = cv_predictions,
      confidence = cv_confidence,
      true_class = labels,
      stringsAsFactors = FALSE
    )
    attr(result, "cv_accuracy") <- mean(cv_predictions == labels)
    attr(result, "cv_folds") <- cv_folds
    return(result)
  }

  # --- Full-data training (no CV) ---
  model <- .train_lda(features, labels, method)
  preds <- .predict_lda(model, features)

  # Store model in metadata
  md <- S4Vectors::metadata(x)
  md$bci_model <- list(
    weights = model$weights,
    threshold = model$threshold,
    classes = model$classes,
    method = method,
    class_means = model$class_means
  )
  S4Vectors::metadata(x) <- md

  data.frame(
    trial = seq_len(n_trials),
    predicted_class = preds$class,
    confidence = preds$confidence,
    true_class = labels,
    stringsAsFactors = FALSE
  )
}


#' Train Fisher's LDA model
#'
#' @param features Numeric matrix (n_trials x n_features).
#' @param labels Character vector of class labels (must have exactly 2 classes).
#' @param method \code{"lda"} or \code{"shrinkage_lda"}.
#' @return A list with \code{weights}, \code{threshold}, \code{classes}, and
#'   \code{class_means}.
#' @keywords internal
.train_lda <- function(features, labels, method = "lda") {
  n_trials <- nrow(features)
  n_features <- ncol(features)
  classes <- unique(labels)

  idx1 <- which(labels == classes[1])
  idx2 <- which(labels == classes[2])
  m1 <- colMeans(features[idx1, , drop = FALSE])
  m2 <- colMeans(features[idx2, , drop = FALSE])

  # Within-class scatter matrix
  centered1 <- scale(features[idx1, , drop = FALSE], center = m1, scale = FALSE)
  centered2 <- scale(features[idx2, , drop = FALSE], center = m2, scale = FALSE)
  Sw <- crossprod(centered1) + crossprod(centered2)

  if (method == "shrinkage_lda") {
    p <- n_features
    n <- n_trials
    trace_Sw <- sum(diag(Sw))
    target <- trace_Sw / p * diag(p)
    alpha <- min(1, max(0.01, p / n))
    Sw <- (1 - alpha) * Sw + alpha * target
  }

  # Regularize Sw to ensure invertibility
  Sw <- Sw + diag(1e-6, n_features)

  w <- solve(Sw, m1 - m2)
  threshold <- as.numeric(t(w) %*% (m1 + m2) / 2)

  list(
    weights = w,
    threshold = threshold,
    classes = classes,
    class_means = list(m1, m2)
  )
}


#' Predict using a trained LDA model
#'
#' @param model A list as returned by \code{.train_lda}.
#' @param features Numeric matrix (n_samples x n_features).
#' @return A list with \code{class} (character vector) and \code{confidence}
#'   (numeric vector).
#' @keywords internal
.predict_lda <- function(model, features) {
  projections <- as.numeric(features %*% model$weights)
  predicted <- ifelse(projections > model$threshold,
                      model$classes[1], model$classes[2])

  distances <- abs(projections - model$threshold)
  max_dist <- max(distances)
  if (max_dist < 1e-10) max_dist <- 1
  confidence <- distances / max_dist

  list(class = predicted, confidence = confidence)
}


#' Iterative Frechet (Riemannian geometric) mean of SPD matrices
#'
#' Computes the geometric mean on the manifold of symmetric positive definite
#' (SPD) matrices using the iterative fixed-point algorithm. Converges to the
#' true Frechet mean in the Riemannian metric.
#'
#' @param covs List of SPD matrices (all same dimensions).
#' @param max_iter Maximum iterations (default: 50).
#' @param tol Convergence tolerance on the tangent-space residual (default: 1e-8).
#' @return The geometric mean SPD matrix.
#' @keywords internal
.geometric_mean_spd <- function(covs, max_iter = 50, tol = 1e-8) {
  # Initialize with arithmetic mean
  C_ref <- Reduce("+", covs) / length(covs)

  for (iter in seq_len(max_iter)) {
    # Compute C_ref^(-1/2)
    eig <- eigen(C_ref, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-10)
    C_ref_isqrt <- eig$vectors %*% diag(vals^(-0.5)) %*% t(eig$vectors)
    C_ref_sqrt  <- eig$vectors %*% diag(vals^(0.5))  %*% t(eig$vectors)

    # Average in tangent space: S = (1/N) * sum_i log(C_ref^{-1/2} C_i C_ref^{-1/2})
    S <- matrix(0, nrow(C_ref), ncol(C_ref))
    for (C in covs) {
      M <- C_ref_isqrt %*% C %*% C_ref_isqrt
      eig_m <- eigen(M, symmetric = TRUE)
      log_M <- eig_m$vectors %*% diag(log(pmax(eig_m$values, 1e-10))) %*%
        t(eig_m$vectors)
      S <- S + log_M
    }
    S <- S / length(covs)

    # Check convergence: S should be near zero at the Frechet mean
    if (max(abs(S)) < tol) break

    # Map back: C_ref <- C_ref^{1/2} exp(S) C_ref^{1/2}
    eig_s <- eigen(S, symmetric = TRUE)
    exp_S <- eig_s$vectors %*% diag(exp(eig_s$values)) %*% t(eig_s$vectors)
    C_ref <- C_ref_sqrt %*% exp_S %*% C_ref_sqrt
  }

  C_ref
}
