#' EEG Independent Component Analysis (ICA)
#'
#' Decomposes multi-channel EEG into independent components using ICA.
#' Supports FastICA, Infomax, and JADE algorithms. Results are stored in the
#' output assay (component activations) and in \code{metadata(x)$ica} (mixing
#' and unmixing matrices).
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param n_components Number of independent components to extract. Defaults to
#'   the number of channels.
#' @param method ICA algorithm: \code{"fastica"}, \code{"infomax"}, or
#'   \code{"jade"}.
#' @param max_iter Maximum number of iterations (default: 200).
#' @param tol Convergence tolerance (default: 1e-6).
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: \code{"ica"}).
#' @return Modified PhysioExperiment with component activations in
#'   \code{output_assay} and ICA metadata in \code{metadata(x)$ica}.
#'   The ICA metadata list contains: \code{mixing} (mixing matrix A),
#'   \code{unmixing} (unmixing matrix), \code{mean} (channel means),
#'   \code{whiten} (whitening matrix), and \code{method} (algorithm used).
#'   The output assay has dimensions n_time x n_components.
#' @references
#' Hyvarinen, A., & Oja, E. (2000). Independent component analysis: algorithms
#' and applications. Neural Networks, 13(4-5), 411-430.
#'
#' Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach
#' to blind separation and blind deconvolution. Neural Computation, 7(6), 1129-1159.
#' @seealso [eegICAremove()], [eegICAdetect()], [eegICAmix()],
#'   [eegFilter()], [eegPreprocess()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' result <- eegICA(pe, n_components = 4, method = "fastica")
#' }
eegICA <- function(x, n_components = NULL, method = c("fastica", "infomax", "jade"),
                   max_iter = 200L, tol = 1e-6,
                   assay_name = NULL, output_assay = "ica_components") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (is.null(n_components)) n_components <- n_channels
  if (n_components > n_channels) {
    stop("n_components cannot exceed number of channels", call. = FALSE)
  }


  # --- Center data ---
  col_means <- colMeans(data)
  centered <- sweep(data, 2, col_means, "-")

  # --- Whiten via PCA ---
  cov_mat <- crossprod(centered) / (n_time - 1)
  eig <- eigen(cov_mat, symmetric = TRUE)
  d <- eig$values[seq_len(n_components)]
  # Guard against zero/negative eigenvalues

  d <- pmax(d, .Machine$double.eps)
  V <- eig$vectors[, seq_len(n_components), drop = FALSE]
  whiten_matrix <- diag(1 / sqrt(d)) %*% t(V)
  dewhiten_matrix <- V %*% diag(sqrt(d))
  whitened <- centered %*% t(whiten_matrix)  # n_time x n_components

  # --- ICA decomposition ---
  if (method == "fastica") {
    W <- .fastica(whitened, n_components, max_iter, tol)

  } else if (method == "infomax") {
    W <- .infomax(whitened, n_components, max_iter, tol)

  } else if (method == "jade") {
    W <- .jade(whitened, n_components, max_iter, tol)
  }

  # Component activations: S = X_whitened %*% t(W)
  sources <- whitened %*% t(W)

  # Mixing matrix A maps components back to channel space:
  # X_centered = S %*% t(A)
  # A = dewhiten_matrix %*% t(W^{-1}) but since W is orthogonal, W^{-1} = t(W)
  A <- dewhiten_matrix %*% t(W)

  # Unmixing matrix maps channel space to components:
  # S = X_centered %*% unmixing^T
  unmixing <- W %*% whiten_matrix

  # Store ICA metadata
  ica_info <- list(
    mixing = A,
    unmixing = unmixing,
    mean = col_means,
    whiten = whiten_matrix,
    method = method
  )
  md <- S4Vectors::metadata(x)
  md$ica <- ica_info
  S4Vectors::metadata(x) <- md

  # Store component activations in metadata (dimensions differ from original assay)
  md[[output_assay]] <- sources
  S4Vectors::metadata(x) <- md

  x
}

#' FastICA algorithm (internal)
#'
#' @param whitened Whitened data matrix (n_time x n_components).
#' @param n_components Number of components.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @return Unmixing matrix W (n_components x n_components).
#' @keywords internal
.fastica <- function(whitened, n_components, max_iter, tol) {
  n <- nrow(whitened)
  p <- n_components

  W <- matrix(rnorm(p * p), p, p)
  # Orthogonalize initial W via SVD
  svd_w <- svd(W)
  W <- svd_w$u %*% t(svd_w$v)

  for (iter in seq_len(max_iter)) {
    W_old <- W

    for (i in seq_len(p)) {
      wx <- whitened %*% W[i, ]
      gwx <- tanh(wx)
      g_prime <- 1 - gwx^2
      W[i, ] <- colMeans(whitened * as.vector(gwx)) - mean(g_prime) * W[i, ]
    }

    # Gram-Schmidt symmetric orthogonalization
    svd_w <- svd(W)
    W <- svd_w$u %*% t(svd_w$v)

    # Check convergence
    convergence <- max(abs(abs(rowSums(W * W_old)) - 1))
    if (convergence < tol) break
  }

  W
}

#' Infomax ICA algorithm (internal)
#'
#' Natural gradient Infomax ICA. Updates W using the rule:
#' dW = (I + (1 - 2*sigmoid(u)) * t(u)) * W
#' where u = W * x.
#'
#' @param whitened Whitened data matrix (n_time x n_components).
#' @param n_components Number of components.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @return Unmixing matrix W (n_components x n_components).
#' @keywords internal
.infomax <- function(whitened, n_components, max_iter, tol) {
  n <- nrow(whitened)
  p <- n_components

  W <- matrix(rnorm(p * p), p, p)
  svd_w <- svd(W)
  W <- svd_w$u %*% t(svd_w$v)

  learning_rate <- 0.01

  for (iter in seq_len(max_iter)) {
    W_old <- W

    # u = X %*% t(W)  =>  each row is W %*% x_i
    u <- whitened %*% t(W)

    # sigmoid
    sig_u <- 1 / (1 + exp(-u))

    # Natural gradient: dW = (I + (1 - 2*sig(u))' * u / n) * W
    bias_term <- t(1 - 2 * sig_u) %*% u / n
    dW <- (diag(p) + bias_term) %*% W

    W <- W + learning_rate * dW

    # Re-orthogonalize to stabilize
    svd_w <- svd(W)
    W <- svd_w$u %*% diag(svd_w$d / max(svd_w$d)) %*% t(svd_w$v)
    # Normalize rows
    row_norms <- sqrt(rowSums(W^2))
    row_norms[row_norms < .Machine$double.eps] <- 1
    W <- W / row_norms

    convergence <- max(abs(abs(rowSums(W * W_old)) - 1))
    if (convergence < tol) break
  }

  # Final orthogonalization
  svd_w <- svd(W)
  W <- svd_w$u %*% t(svd_w$v)

  W
}

#' JADE ICA algorithm (internal)
#'
#' Joint Approximate Diagonalization of Eigenmatrices. Computes fourth-order
#' cumulant matrices and performs joint approximate diagonalization using
#' Jacobi rotations.
#'
#' @param whitened Whitened data matrix (n_time x n_components).
#' @param n_components Number of components.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @return Unmixing matrix W (n_components x n_components).
#' @keywords internal
.jade <- function(whitened, n_components, max_iter, tol) {
  n <- nrow(whitened)
  p <- n_components

  # Compute 4th-order cumulant matrices
  # For each pair (k, l), compute Q_kl = E[x * x' * x_k * x_l] - symmetrized
  n_matrices <- p * (p + 1) / 2
  cum_matrices <- vector("list", n_matrices)

  idx <- 0
  for (k in seq_len(p)) {
    for (l in k:p) {
      idx <- idx + 1
      # Q_kl[i,j] = E[x_i * x_j * x_k * x_l] - delta_ij*delta_kl - delta_ik*delta_jl - delta_il*delta_jk
      xk <- whitened[, k]
      xl <- whitened[, l]
      Q <- crossprod(whitened * (xk * xl), whitened) / n

      # Subtract Gaussian part
      Q <- Q - (if (k == l) diag(p) else matrix(0, p, p))
      ek <- rep(0, p); ek[k] <- 1
      el <- rep(0, p); el[l] <- 1
      Q <- Q - tcrossprod(ek, el) - tcrossprod(el, ek)

      cum_matrices[[idx]] <- Q
    }
  }

  # Joint approximate diagonalization via Jacobi rotations
  V <- diag(p)

  for (sweep_iter in seq_len(max_iter)) {
    max_off <- 0

    for (i in seq_len(p - 1)) {
      for (j in (i + 1):p) {
        # Compute Givens rotation angle
        h_ii <- 0; h_jj <- 0; h_ij <- 0
        for (m in seq_along(cum_matrices)) {
          M <- cum_matrices[[m]]
          h_ii <- h_ii + (M[i, i] - M[j, j])^2
          h_jj <- h_jj + (M[i, j] + M[j, i])^2
          h_ij <- h_ij + (M[i, i] - M[j, j]) * (M[i, j] + M[j, i])
        }

        theta <- 0.5 * atan2(2 * h_ij, h_ii - h_jj)
        c_theta <- cos(theta)
        s_theta <- sin(theta)

        if (abs(s_theta) > tol) {
          # Apply Givens rotation to all cumulant matrices
          for (m in seq_along(cum_matrices)) {
            M <- cum_matrices[[m]]
            # Rotate rows i, j
            row_i <- M[i, ]
            row_j <- M[j, ]
            M[i, ] <- c_theta * row_i + s_theta * row_j
            M[j, ] <- -s_theta * row_i + c_theta * row_j
            # Rotate columns i, j
            col_i <- M[, i]
            col_j <- M[, j]
            M[, i] <- c_theta * col_i + s_theta * col_j
            M[, j] <- -s_theta * col_i + c_theta * col_j
            cum_matrices[[m]] <- M
          }

          # Update rotation matrix V
          col_i <- V[, i]
          col_j <- V[, j]
          V[, i] <- c_theta * col_i + s_theta * col_j
          V[, j] <- -s_theta * col_i + c_theta * col_j

          max_off <- max(max_off, abs(s_theta))
        }
      }
    }

    if (max_off < tol) break
  }

  W <- t(V)
  W
}


#' Remove ICA Components from EEG
#'
#' Reconstructs EEG data with specified independent components removed.
#' The removed components are zeroed out in the mixing matrix before
#' back-projecting to channel space.
#'
#' @param x A PhysioExperiment object with ICA results (from \code{eegICA}).
#' @param components Integer vector of component indices to remove.
#' @param ica_assay Assay name containing ICA component activations
#'   (default: \code{"ica"}).
#' @param output_assay Output assay name (default: \code{"ica_cleaned"}).
#' @return Modified PhysioExperiment with cleaned data in \code{output_assay}.
#'   The cleaned assay contains reconstructed channel data with specified
#'   components removed. Dimensions match the original data.
#' @references
#' Hyvarinen, A., & Oja, E. (2000). Independent component analysis: algorithms
#' and applications. Neural Networks, 13(4-5), 411-430.
#'
#' Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach
#' to blind separation and blind deconvolution. Neural Computation, 7(6), 1129-1159.
#' @seealso [eegICA()], [eegICAdetect()], [eegICAmix()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' pe <- eegICA(pe, n_components = 4, method = "fastica")
#' pe <- eegICAremove(pe, components = c(1, 2))
#' }
eegICAremove <- function(x, components, ica_assay = "ica_components",
                         output_assay = "ica_cleaned") {
  stopifnot(inherits(x, "PhysioExperiment"))

  ica_info <- S4Vectors::metadata(x)$ica
  if (is.null(ica_info)) {
    stop("No ICA results found. Run eegICA() first.", call. = FALSE)
  }

  A <- ica_info$mixing
  col_means <- ica_info$mean
  sources <- S4Vectors::metadata(x)[[ica_assay]]
  if (is.null(sources)) {
    stop(sprintf("ICA activations '%s' not found in metadata. Run eegICA() first.",
                 ica_assay), call. = FALSE)
  }

  # Zero out columns of mixing matrix for specified components
  A_modified <- A
  A_modified[, components] <- 0

  # Reconstruct: cleaned = sources %*% t(A_modified) + mean
  cleaned <- sources %*% t(A_modified)
  cleaned <- sweep(cleaned, 2, col_means, "+")

  md <- S4Vectors::metadata(x)
  md[[output_assay]] <- cleaned
  S4Vectors::metadata(x) <- md

  x
}


#' Detect Artifact ICA Components
#'
#' Automatically identifies artifact components using one of three methods:
#' correlation with frontal channels, kurtosis, or spatial weight pattern.
#'
#' @param x A PhysioExperiment object with ICA results (from \code{eegICA}).
#' @param method Detection method: \code{"correlation"} (frontal channel
#'   correlation), \code{"kurtosis"} (excess kurtosis), or \code{"spatial"}
#'   (spatial weight pattern).
#' @param threshold Threshold for artifact detection. For \code{"correlation"},
#'   absolute correlation > threshold marks artifact (default: 0.3).
#' @param ica_assay Assay name containing ICA activations (default: \code{"ica"}).
#' @return A data.frame with columns: \code{component} (integer index),
#'   \code{type} (\code{"artifact"} or \code{"neural"}), \code{method}
#'   (detection method used), and \code{score} (numeric detection score).
#' @references
#' Hyvarinen, A., & Oja, E. (2000). Independent component analysis: algorithms
#' and applications. Neural Networks, 13(4-5), 411-430.
#'
#' Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach
#' to blind separation and blind deconvolution. Neural Computation, 7(6), 1129-1159.
#' @seealso [eegICA()], [eegICAremove()], [eegICAmix()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' pe <- eegICA(pe, n_components = 4, method = "fastica")
#' artifacts <- eegICAdetect(pe, method = "kurtosis")
#' }
eegICAdetect <- function(x, method = c("correlation", "kurtosis", "spatial"),
                         threshold = 0.3, ica_assay = "ica_components") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  ica_info <- S4Vectors::metadata(x)$ica
  if (is.null(ica_info)) {
    stop("No ICA results found. Run eegICA() first.", call. = FALSE)
  }

  sources <- S4Vectors::metadata(x)[[ica_assay]]
  if (is.null(sources)) {
    stop(sprintf("ICA activations '%s' not found in metadata. Run eegICA() first.",
                 ica_assay), call. = FALSE)
  }
  n_components <- ncol(sources)
  A <- ica_info$mixing

  component_ids <- seq_len(n_components)
  types <- rep("neural", n_components)
  scores <- rep(NA_real_, n_components)

  if (method == "correlation") {
    # Get channel labels
    col_data <- SummarizedExperiment::colData(x)
    ch_labels <- if ("label" %in% colnames(col_data)) {
      as.character(col_data$label)
    } else {
      paste0("Ch", seq_len(ncol(SummarizedExperiment::assay(x, defaultAssay(x)))))
    }

    # Frontal channels
    frontal <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8")
    frontal_idx <- which(ch_labels %in% frontal)

    if (length(frontal_idx) > 0) {
      raw_data <- SummarizedExperiment::assay(x, defaultAssay(x))
      # Handle 2D data (use raw channel signals)
      for (ic in seq_len(n_components)) {
        max_cor <- 0
        for (fi in frontal_idx) {
          r <- abs(cor(sources[, ic], raw_data[, fi]))
          if (!is.na(r) && r > max_cor) max_cor <- r
        }
        scores[ic] <- max_cor
        if (max_cor > threshold) types[ic] <- "artifact"
      }
    }

  } else if (method == "kurtosis") {
    for (ic in seq_len(n_components)) {
      s <- sources[, ic]
      n <- length(s)
      m <- mean(s)
      m2 <- sum((s - m)^2) / n
      m4 <- sum((s - m)^4) / n
      kurt <- m4 / (m2^2) - 3  # excess kurtosis
      scores[ic] <- kurt
      # High kurtosis = super-Gaussian (muscle artifact)
      # Low kurtosis = sub-Gaussian (eye blink)
      if (kurt > 5 || kurt < -1) types[ic] <- "artifact"
    }

  } else if (method == "spatial") {
    col_data <- SummarizedExperiment::colData(x)
    ch_labels <- if ("label" %in% colnames(col_data)) {
      as.character(col_data$label)
    } else {
      paste0("Ch", seq_len(nrow(A)))
    }

    frontal <- c("Fp1", "Fp2")
    frontal_idx <- which(ch_labels %in% frontal)

    for (ic in seq_len(n_components)) {
      weights <- abs(A[, ic])
      max_idx <- which.max(weights)
      scores[ic] <- weights[max_idx] / sum(weights)

      if (length(frontal_idx) > 0 && max_idx %in% frontal_idx) {
        types[ic] <- "artifact"
      }
    }
  }

  data.frame(
    component = component_ids,
    type = types,
    method = rep(method, n_components),
    score = scores,
    stringsAsFactors = FALSE
  )
}


#' Access ICA Metadata
#'
#' Returns the ICA results stored in \code{metadata(x)$ica}, including the
#' mixing matrix, unmixing matrix, channel means, and whitening matrix.
#'
#' @param x A PhysioExperiment object with ICA results.
#' @return A list containing: \code{mixing} (mixing matrix A, n_channels x
#'   n_components), \code{unmixing} (unmixing matrix, n_components x
#'   n_channels), \code{mean} (channel mean vector), \code{whiten}
#'   (whitening matrix), and \code{method} (character, algorithm used).
#' @references
#' Hyvarinen, A., & Oja, E. (2000). Independent component analysis: algorithms
#' and applications. Neural Networks, 13(4-5), 411-430.
#'
#' Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach
#' to blind separation and blind deconvolution. Neural Computation, 7(6), 1129-1159.
#' @seealso [eegICA()], [eegICAremove()], [eegICAdetect()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, sr = 500)
#' pe <- eegICA(pe, n_components = 4, method = "fastica")
#' ica_info <- eegICAmix(pe)
#' }
eegICAmix <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))

  ica_info <- S4Vectors::metadata(x)$ica
  if (is.null(ica_info)) {
    stop("No ICA results found in metadata. Run eegICA() first.", call. = FALSE)
  }

  ica_info
}
