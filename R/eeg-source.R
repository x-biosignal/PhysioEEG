#' Construct EEG Forward Model (Leadfield Matrix)
#'
#' Constructs a leadfield matrix that maps brain source activity to scalp
#' electrode potentials. Uses electrode positions from \code{colData(x)} if
#' available, otherwise falls back to standard 10-20 system positions on a
#' unit sphere.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param method Forward model method: \code{"spherical"} (current dipole in
#'   infinite homogeneous medium with conductivity 0.33 S/m) or
#'   \code{"bem_simplified"} (3-shell Berg and Scherg 1994 approximation using
#'   virtual dipoles at scaled eccentricities).
#' @param n_sources Number of dipole sources to distribute inside the head
#'   model (default: 500).
#' @param assay_name Name of the assay to reference for channel count. If
#'   \code{NULL}, the default assay is used.
#' @return A list with components:
#'   \describe{
#'     \item{leadfield}{Numeric matrix of dimensions n_electrodes x
#'       (n_sources * 3). Each source has 3 orientation columns (x, y, z).}
#'     \item{source_positions}{Data frame with columns \code{x}, \code{y},
#'       \code{z} for each source location.}
#'     \item{electrode_positions}{Data frame with columns \code{label},
#'       \code{x}, \code{y}, \code{z} for each electrode.}
#'     \item{n_sources}{Integer number of source dipoles.}
#'   }
#' @references
#' Pascual-Marqui, R. D. (2002). Standardized low-resolution brain electromagnetic
#' tomography (sLORETA). Methods and Findings in Experimental and Clinical
#' Pharmacology, 24(Suppl D), 5-12.
#' @seealso [eegSourceEstimate()], [eegBeamformer()], [eegSourcePower()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
#' fm <- eegForwardModel(pe, method = "spherical", n_sources = 100)
#' dim(fm$leadfield)
#' }
eegForwardModel <- function(x, method = c("spherical", "bem_simplified"),
                            n_sources = 500, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_channels <- ncol(data)

  # --- Electrode positions ---
  cd <- SummarizedExperiment::colData(x)
  if (all(c("pos_x", "pos_y", "pos_z") %in% colnames(cd))) {
    elec_df <- data.frame(
      label = as.character(cd$label),
      x = as.numeric(cd$pos_x),
      y = as.numeric(cd$pos_y),
      z = as.numeric(cd$pos_z),
      stringsAsFactors = FALSE
    )
  } else {
    # Standard 10-20 positions on unit sphere
    standard_1020 <- data.frame(
      label = c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                "T3", "C3", "Cz", "C4", "T4",
                "T5", "P3", "Pz", "P4", "T6",
                "O1", "O2"),
      x = c(-0.31, 0.31, -0.81, -0.55, 0.00, 0.55, 0.81,
            -0.95, -0.59, 0.00, 0.59, 0.95,
            -0.81, -0.55, 0.00, 0.55, 0.81,
            -0.31, 0.31),
      y = c(0.95, 0.95, 0.59, 0.67, 0.72, 0.67, 0.59,
            0.00, 0.00, 0.00, 0.00, 0.00,
            -0.59, -0.67, -0.72, -0.67, -0.59,
            -0.95, -0.95),
      z = c(-0.03, -0.03, -0.03, 0.50, 0.69, 0.50, -0.03,
            0.00, 0.59, 0.72, 0.59, 0.00,
            -0.03, 0.50, 0.69, 0.50, -0.03,
            -0.03, -0.03),
      stringsAsFactors = FALSE
    )

    # Match channels by label if possible
    ch_labels <- as.character(cd$label)
    matched <- match(ch_labels, standard_1020$label)
    if (all(!is.na(matched))) {
      elec_df <- standard_1020[matched, ]
    } else {
      # Use first n_channels positions
      if (n_channels <= nrow(standard_1020)) {
        elec_df <- standard_1020[seq_len(n_channels), ]
      } else {
        # Generate evenly spaced positions on sphere for extra channels
        elec_df <- standard_1020
        n_extra <- n_channels - nrow(standard_1020)
        set.seed(123)
        extra_theta <- seq(0, 2 * pi, length.out = n_extra + 1)[seq_len(n_extra)]
        extra_phi <- seq(0.3, 0.8, length.out = n_extra)
        extra_df <- data.frame(
          label = paste0("E", seq_len(n_extra)),
          x = sin(extra_phi) * cos(extra_theta),
          y = sin(extra_phi) * sin(extra_theta),
          z = cos(extra_phi),
          stringsAsFactors = FALSE
        )
        elec_df <- rbind(elec_df, extra_df)
      }
    }
  }
  rownames(elec_df) <- NULL

  # --- Source positions: uniform inside sphere of radius 0.8 ---
  set.seed(42)
  theta <- stats::runif(n_sources, 0, 2 * pi)
  phi <- acos(stats::runif(n_sources, -1, 1))
  r <- 0.8 * stats::runif(n_sources)^(1 / 3)
  src_x <- r * sin(phi) * cos(theta)
  src_y <- r * sin(phi) * sin(theta)
  src_z <- r * cos(phi)

  src_df <- data.frame(x = src_x, y = src_y, z = src_z,
                       stringsAsFactors = FALSE)

  # --- Compute leadfield matrix ---

  # Conductivity of homogeneous medium (brain tissue, S/m)
  sigma <- 0.33

  n_elec <- nrow(elec_df)
  L <- matrix(0, nrow = n_elec, ncol = n_sources * 3)

  if (method == "spherical") {
    # Current dipole in infinite homogeneous medium:
    #   phi(r) = (1 / (4 * pi * sigma)) * (p . (r - r0)) / |r - r0|^3
    # Leadfield for k-th orientation of source j at electrode i:
    #   L[i, 3*(j-1)+k] = (r_i - r_j)[k] / (4 * pi * sigma * |r_i - r_j|^3)
    scale <- 1 / (4 * pi * sigma)

    for (j in seq_len(n_sources)) {
      for (i in seq_len(n_elec)) {
        dx <- elec_df$x[i] - src_df$x[j]
        dy <- elec_df$y[i] - src_df$y[j]
        dz <- elec_df$z[i] - src_df$z[j]
        r_dist <- sqrt(dx^2 + dy^2 + dz^2)
        if (r_dist < 1e-10) r_dist <- 1e-10  # avoid division by zero
        r_dist3 <- r_dist^3

        col_base <- (j - 1) * 3
        L[i, col_base + 1] <- scale * dx / r_dist3
        L[i, col_base + 2] <- scale * dy / r_dist3
        L[i, col_base + 3] <- scale * dz / r_dist3
      }
    }

  } else if (method == "bem_simplified") {
    # Berg & Scherg (1994) 3-shell approximation:
    # Approximate the 3-layer BEM (brain/skull/scalp) using 3 virtual dipoles
    # at scaled eccentricities with corresponding weights.
    berg_mu     <- c(0.3320, 0.4080, -0.7400)   # expansion weights
    berg_lambda <- c(0.5401, 0.8060,  1.0000)    # eccentricity scale factors

    scale <- 1 / (4 * pi * sigma)

    for (j in seq_len(n_sources)) {
      for (i in seq_len(n_elec)) {
        lf_x <- 0
        lf_y <- 0
        lf_z <- 0

        for (vi in seq_len(3)) {
          # Virtual source at scaled position
          vsx <- berg_lambda[vi] * src_df$x[j]
          vsy <- berg_lambda[vi] * src_df$y[j]
          vsz <- berg_lambda[vi] * src_df$z[j]

          dx <- elec_df$x[i] - vsx
          dy <- elec_df$y[i] - vsy
          dz <- elec_df$z[i] - vsz
          r_dist <- sqrt(dx^2 + dy^2 + dz^2)
          if (r_dist < 1e-10) r_dist <- 1e-10
          r_dist3 <- r_dist^3

          lf_x <- lf_x + berg_mu[vi] * dx / r_dist3
          lf_y <- lf_y + berg_mu[vi] * dy / r_dist3
          lf_z <- lf_z + berg_mu[vi] * dz / r_dist3
        }

        col_base <- (j - 1) * 3
        L[i, col_base + 1] <- scale * lf_x
        L[i, col_base + 2] <- scale * lf_y
        L[i, col_base + 3] <- scale * lf_z
      }
    }
  }

  list(
    leadfield = L,
    source_positions = src_df,
    electrode_positions = elec_df,
    n_sources = n_sources
  )
}


#' EEG Source Estimation
#'
#' Estimates brain source activity from scalp EEG data using distributed
#' source imaging methods. Requires a forward model from
#' \code{\link{eegForwardModel}}.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param forward_model A forward model list as returned by
#'   \code{\link{eegForwardModel}}.
#' @param method Source estimation method: \code{"sloreta"} (standardized
#'   low-resolution tomography), \code{"eloreta"} (exact low-resolution
#'   tomography), or \code{"mne"} (minimum norm estimate).
#' @param lambda Regularization parameter (default: 0.05). Higher values
#'   produce smoother solutions.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name for the output assay containing source estimates
#'   (default: \code{"source"}).
#' @return Modified PhysioExperiment with source estimates stored in
#'   \code{output_assay}. The assay is a matrix of dimensions
#'   n_time x (n_sources * 3). Sets \code{metadata(x)$source_estimate}
#'   with a list containing: \code{method}, \code{lambda},
#'   \code{n_sources}, and \code{n_source_cols}.
#' @references
#' Pascual-Marqui, R. D. (2002). Standardized low-resolution brain electromagnetic
#' tomography (sLORETA). Methods and Findings in Experimental and Clinical
#' Pharmacology, 24(Suppl D), 5-12.
#' @seealso [eegForwardModel()], [eegBeamformer()], [eegSourcePower()],
#'   [eegPlotSource()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
#' fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
#' pe <- eegSourceEstimate(pe, fm, method = "mne")
#' }
eegSourceEstimate <- function(x, forward_model,
                              method = c("sloreta", "eloreta", "mne"),
                              lambda = 0.05, assay_name = NULL,
                              output_assay = "source") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  L <- forward_model$leadfield
  n_elec <- nrow(L)
  n_src_cols <- ncol(L)

  # Data is time x channels, we need channels x time for source estimation
  Y <- t(data)  # n_channels x n_time

  # Gram matrix: L %*% t(L)
  LLt <- L %*% t(L)
  reg_matrix <- LLt + lambda * diag(n_elec)

  # Check condition number before inversion; increase regularization if needed
  cond_num <- rcond(reg_matrix)
  if (cond_num < .Machine$double.eps) {
    warning("Near-singular matrix in source estimation. Increasing regularization.",
            call. = FALSE)
    lambda <- lambda * 10
    reg_matrix <- LLt + lambda * diag(n_elec)
  }

  inv_reg <- solve(reg_matrix)

  if (method == "mne") {
    # MNE: J = L' * inv(L*L' + lambda*I) * Y
    kernel <- t(L) %*% inv_reg
    J <- kernel %*% Y  # (n_sources*3) x n_time

  } else if (method == "sloreta") {
    # sLORETA: MNE normalized by resolution matrix diagonal
    kernel <- t(L) %*% inv_reg
    J_mne <- kernel %*% Y

    # Resolution matrix R = L' * inv(L*L' + lambda*I) * L
    R <- kernel %*% L
    R_diag <- diag(R)
    # Avoid division by zero
    R_diag[R_diag < 1e-20] <- 1e-20
    norm_factors <- 1 / sqrt(abs(R_diag))

    J <- J_mne * norm_factors  # element-wise scaling by row

  } else if (method == "eloreta") {
    # eLORETA: iterative weight estimation with convergence check
    W_diag <- rep(1, n_src_cols)
    W <- diag(W_diag)

    max_iter <- 20
    tol <- 1e-6

    for (iter in seq_len(max_iter)) {
      W_old <- W_diag

      LW <- L %*% W
      LWLt <- LW %*% t(L)
      reg_iter <- LWLt + lambda * diag(n_elec)

      # Check condition number for each iteration
      cond_iter <- rcond(reg_iter)
      if (cond_iter < .Machine$double.eps) {
        warning("Near-singular matrix in eLORETA iteration. Increasing regularization.",
                call. = FALSE)
        reg_iter <- reg_iter + lambda * 9 * diag(n_elec)
      }

      inv_iter <- solve(reg_iter)
      R_iter <- t(L) %*% inv_iter %*% L %*% W
      R_diag_iter <- diag(R_iter)
      R_diag_iter[abs(R_diag_iter) < 1e-20] <- 1e-20
      W_diag <- 1 / sqrt(abs(R_diag_iter))
      W <- diag(W_diag)

      # Check convergence
      if (max(abs(W_diag - W_old)) < tol) break
    }

    LW <- L %*% W
    LWLt <- LW %*% t(L)
    reg_final <- LWLt + lambda * diag(n_elec)

    cond_final <- rcond(reg_final)
    if (cond_final < .Machine$double.eps) {
      warning("Near-singular matrix in final eLORETA solve. Increasing regularization.",
              call. = FALSE)
      reg_final <- reg_final + lambda * 9 * diag(n_elec)
    }

    inv_final <- solve(reg_final)
    kernel <- W %*% t(L) %*% inv_final
    J <- kernel %*% Y
  }

  # Store as (n_sources*3) x n_time -> transpose to time x sources for assay
  # Assay convention: rows = features (sources), cols = samples (time points)
  # But PhysioExperiment convention: rows = time, cols = channels/features
  # Store as n_time x n_source_cols to match PhysioExperiment convention
  source_data <- t(J)  # n_time x (n_sources*3)

  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- source_data
  md$source_estimate <- list(
    method = method,
    lambda = lambda,
    n_sources = forward_model$n_sources,
    n_source_cols = n_src_cols
  )
  S4Vectors::metadata(x) <- md

  x
}


#' EEG Beamformer Source Localization
#'
#' Applies spatial filtering (beamforming) to localize neural source power.
#' Linearly Constrained Minimum Variance (LCMV) beamformer operates in the
#' time domain. Dynamic Imaging of Coherent Sources (DICS) operates in the
#' frequency domain.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param forward_model A forward model list as returned by
#'   \code{\link{eegForwardModel}}.
#' @param method Beamformer method: \code{"lcmv"} (Linearly Constrained
#'   Minimum Variance) or \code{"dics"} (Dynamic Imaging of Coherent Sources).
#' @param freq_range Numeric vector of length 2 specifying frequency range
#'   in Hz for DICS method (e.g., \code{c(8, 13)} for alpha band). Ignored
#'   for LCMV.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @param output_assay Name for the output assay containing beamformer
#'   results (default: \code{"beamformer"}).
#' @return Modified PhysioExperiment with source power stored in
#'   \code{output_assay} as a matrix with one column named \code{"power"}
#'   (n_sources rows). Each value represents the estimated source power
#'   at the corresponding dipole location.
#' @references
#' Pascual-Marqui, R. D. (2002). Standardized low-resolution brain electromagnetic
#' tomography (sLORETA). Methods and Findings in Experimental and Clinical
#' Pharmacology, 24(Suppl D), 5-12.
#'
#' Van Veen, B. D., et al. (1997). Localization of brain electrical activity via
#' linearly constrained minimum variance spatial filtering. IEEE Transactions on
#' Biomedical Engineering, 44(9), 867-880.
#' @seealso [eegForwardModel()], [eegSourceEstimate()], [eegSourcePower()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
#' fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
#' pe <- eegBeamformer(pe, fm, method = "lcmv")
#' }
eegBeamformer <- function(x, forward_model,
                          method = c("lcmv", "dics"),
                          freq_range = NULL, assay_name = NULL,
                          output_assay = "beamformer") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  L <- forward_model$leadfield
  n_sources <- forward_model$n_sources

  if (method == "lcmv") {
    # Data covariance matrix
    C <- stats::cov(data)  # n_channels x n_channels
    C_inv <- solve(C + 1e-6 * diag(n_channels))

    source_power <- numeric(n_sources)

    for (j in seq_len(n_sources)) {
      # Get leadfield columns for source j (3 orientations)
      col_idx <- ((j - 1) * 3 + 1):(j * 3)
      Lj <- L[, col_idx, drop = FALSE]

      # For each orientation, compute beamformer weight and power
      # Use orientation with maximum power
      max_power <- 0
      for (ori in seq_len(ncol(Lj))) {
        lj <- Lj[, ori]
        denom <- as.numeric(t(lj) %*% C_inv %*% lj)
        if (abs(denom) < 1e-20) next
        w_j <- C_inv %*% lj / denom
        pw <- as.numeric(t(w_j) %*% C %*% w_j)
        if (pw > max_power) max_power <- pw
      }
      source_power[j] <- max_power
    }

  } else if (method == "dics") {
    # Cross-spectral density in freq_range
    if (is.null(freq_range)) freq_range <- c(8, 13)  # default: alpha band
    stopifnot(length(freq_range) == 2)

    # Compute cross-spectral density matrix via FFT
    freqs <- seq(0, sr / 2, length.out = n_time %/% 2 + 1)
    freq_idx <- which(freqs >= freq_range[1] & freqs <= freq_range[2])

    # FFT of each channel
    fft_data <- matrix(0 + 0i, nrow = length(freqs), ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      ft <- fft(data[, ch])
      fft_data[, ch] <- ft[seq_len(length(freqs))]
    }

    # Average cross-spectral density over frequency band
    CSD <- matrix(0, nrow = n_channels, ncol = n_channels)
    if (length(freq_idx) > 0) {
      for (fi in freq_idx) {
        v <- fft_data[fi, ]
        CSD <- CSD + Re(outer(v, Conj(v)))
      }
      CSD <- CSD / length(freq_idx)
    } else {
      # Fallback to broadband
      CSD <- stats::cov(data)
    }

    CSD_inv <- solve(CSD + 1e-6 * diag(n_channels))

    source_power <- numeric(n_sources)

    for (j in seq_len(n_sources)) {
      col_idx <- ((j - 1) * 3 + 1):(j * 3)
      Lj <- L[, col_idx, drop = FALSE]

      max_power <- 0
      for (ori in seq_len(ncol(Lj))) {
        lj <- Lj[, ori]
        denom <- as.numeric(Re(t(lj) %*% CSD_inv %*% lj))
        if (abs(denom) < 1e-20) next
        w_j <- CSD_inv %*% lj / denom
        pw <- abs(as.numeric(t(w_j) %*% CSD %*% w_j))
        if (pw > max_power) max_power <- pw
      }
      source_power[j] <- max_power
    }
  }

  # Store as n_sources x 1 matrix in assay
  # To fit assay convention (rows = time/observations), store as a matrix
  power_matrix <- matrix(source_power, ncol = 1)
  colnames(power_matrix) <- "power"

  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- power_matrix
  S4Vectors::metadata(x) <- md

  x
}


#' Compute Source-Space Band Power
#'
#' Computes spectral band power for each source from source-estimated EEG
#' data. Requires prior source estimation via \code{\link{eegSourceEstimate}}.
#'
#' @param x A PhysioExperiment object with source estimates.
#' @param bands Named list of frequency bands, each a numeric vector of
#'   length 2 (lower, upper Hz). If \code{NULL}, uses standard EEG bands:
#'   delta (1-4), theta (4-8), alpha (8-13), beta (13-30), gamma (30-50).
#' @param source_assay Name of the assay containing source data (default:
#'   \code{"source"}).
#' @return A data.frame with columns:
#'   \describe{
#'     \item{source_id}{Integer source index.}
#'     \item{band}{Character name of the frequency band.}
#'     \item{power}{Numeric spectral power in the band.}
#'   }
#' @references
#' Pascual-Marqui, R. D. (2002). Standardized low-resolution brain electromagnetic
#' tomography (sLORETA). Methods and Findings in Experimental and Clinical
#' Pharmacology, 24(Suppl D), 5-12.
#' @seealso [eegSourceEstimate()], [eegForwardModel()], [eegBeamformer()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
#' fm <- eegForwardModel(pe, method = "spherical", n_sources = 20)
#' pe <- eegSourceEstimate(pe, fm, method = "mne")
#' bp <- eegSourcePower(pe)
#' head(bp)
#' }
eegSourcePower <- function(x, bands = NULL, source_assay = "source") {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(bands)) {
    bands <- list(
      delta = c(1, 4),
      theta = c(4, 8),
      alpha = c(8, 13),
      beta = c(13, 30),
      gamma = c(30, 50)
    )
  }

  src_data <- S4Vectors::metadata(x)[[source_assay]]
  if (is.null(src_data)) {
    stop(sprintf("Source data '%s' not found in metadata(x). Run eegSourceEstimate() first.",
                 source_assay), call. = FALSE)
  }
  sr <- samplingRate(x)
  n_time <- nrow(src_data)
  n_src_cols <- ncol(src_data)

  # Frequency resolution
  freqs <- seq(0, sr / 2, length.out = n_time %/% 2 + 1)

  results <- list()

  for (src in seq_len(n_src_cols)) {
    sig <- src_data[, src]
    ft <- fft(sig)
    # Power spectrum (one-sided)
    n_freq <- length(freqs)
    psd <- (Mod(ft[seq_len(n_freq)])^2) / n_time

    for (band_name in names(bands)) {
      band_range <- bands[[band_name]]
      freq_idx <- which(freqs >= band_range[1] & freqs <= band_range[2])
      band_power <- if (length(freq_idx) > 0) sum(psd[freq_idx]) else 0

      results[[length(results) + 1]] <- data.frame(
        source_id = src,
        band = band_name,
        power = band_power,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}
