#' EEG Microstate Segmentation
#'
#' Performs microstate analysis on EEG data by identifying dominant scalp
#' topographies at Global Field Power (GFP) peaks and assigning each time
#' point to the best-matching microstate map. Supports polarity-invariant
#' K-means, atomize-and-agglomerate hierarchical clustering (AAHC), and
#' PCA-based extraction.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param n_states Number of microstate classes to extract (default: 4).
#' @param method Clustering method: \code{"kmeans"} (polarity-invariant
#'   K-means), \code{"aahc"} (atomize and agglomerate hierarchical
#'   clustering), or \code{"pca"} (principal component analysis).
#' @param min_gfp Percentile threshold (0-100) for GFP peak selection
#'   (default: 1.0). Only GFP peaks above this percentile are used for
#'   clustering.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return Modified PhysioExperiment with microstate results stored in
#'   \code{metadata(x)$microstates}, a list containing:
#'   \describe{
#'     \item{maps}{Numeric matrix of dimensions n_channels x n_states,
#'       each column a microstate topography.}
#'     \item{labels}{Integer vector of length n_time, microstate assignment
#'       (1 to n_states) for each time point.}
#'     \item{gfp}{Numeric vector of GFP values per time point.}
#'     \item{n_states}{Integer number of microstate classes.}
#'   }
#' @references
#' Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for studying
#' the temporal dynamics of whole-brain neuronal networks. NeuroImage, 180, 577-593.
#' @seealso [eegMicrostateStats()], [eegMicrostateBackfit()],
#'   [eegMicrostateSequence()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
#' ms <- metadata(pe)$microstates
#' }
eegMicrostates <- function(x, n_states = 4,
                           method = c("kmeans", "aahc", "pca"),
                           min_gfp = 1.0, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(n_states) && n_states >= 2)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  # --- Step 1: Compute GFP (Global Field Power) ---
  row_means <- rowMeans(data)
  centered <- data - row_means
  gfp <- sqrt(rowMeans(centered^2))

  # --- Step 2: Find GFP peaks (local maxima above threshold) ---
  gfp_threshold <- stats::quantile(gfp, min_gfp / 100)

  peak_idx <- integer(0)
  for (i in 2:(n_time - 1)) {
    if (gfp[i] > gfp[i - 1] && gfp[i] >= gfp[i + 1] && gfp[i] > gfp_threshold) {
      peak_idx <- c(peak_idx, i)
    }
  }

  if (length(peak_idx) < n_states) {
    # Fallback: use top GFP time points
    peak_idx <- order(gfp, decreasing = TRUE)[seq_len(min(n_states * 10, n_time))]
  }

  # --- Step 3: Extract topographies at GFP peaks ---
  peak_topos <- data[peak_idx, , drop = FALSE]  # n_peaks x n_channels

  # --- Step 4: Cluster topographies ---
  if (method == "kmeans") {
    maps <- .microstate_kmeans(peak_topos, n_states, max_iter = 100)

  } else if (method == "aahc") {
    maps <- .microstate_aahc(peak_topos, n_states)

  } else if (method == "pca") {
    maps <- .microstate_pca(peak_topos, n_states)
  }

  # --- Step 5: Backfit all time points ---
  labels <- .microstate_backfit(data, maps)

  S4Vectors::metadata(x)$microstates <- list(
    maps = maps,
    labels = labels,
    gfp = gfp,
    n_states = n_states
  )

  x
}


#' Compute Microstate Statistics
#'
#' Calculates temporal statistics for each microstate class from a segmented
#' EEG recording: mean duration, occurrence rate, and time coverage.
#' Also computes the transition probability matrix between states.
#'
#' @param x A PhysioExperiment object with microstate labels in
#'   \code{metadata(x)$microstates} (from \code{\link{eegMicrostates}}).
#' @return A data.frame with columns:
#'   \describe{
#'     \item{state}{Integer microstate class (1 to n_states).}
#'     \item{duration_ms}{Mean duration of consecutive runs in milliseconds.}
#'     \item{occurrence_per_sec}{Number of state occurrences (runs) per second.}
#'     \item{coverage_pct}{Percentage of total time spent in this state.}
#'   }
#'   The transition probability matrix (n_states x n_states) is stored as
#'   an attribute \code{"transition_matrix"}.
#' @references
#' Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for studying
#' the temporal dynamics of whole-brain neuronal networks. NeuroImage, 180, 577-593.
#' @seealso [eegMicrostates()], [eegMicrostateBackfit()],
#'   [eegMicrostateSequence()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
#' stats <- eegMicrostateStats(pe)
#' print(stats)
#' attr(stats, "transition_matrix")
#' }
eegMicrostateStats <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  ms <- S4Vectors::metadata(x)$microstates
  if (is.null(ms)) {
    stop("No microstate results found. Run eegMicrostates() first.", call. = FALSE)
  }

  labels <- ms$labels
  n_states <- ms$n_states
  sr <- samplingRate(x)
  n_time <- length(labels)
  total_sec <- n_time / sr

  # --- Run-length encoding ---
  rle_result <- rle(labels)
  run_values <- rle_result$values
  run_lengths <- rle_result$lengths

  duration_ms <- numeric(n_states)
  occurrence <- numeric(n_states)
  coverage_pct <- numeric(n_states)

  for (s in seq_len(n_states)) {
    state_runs <- run_lengths[run_values == s]
    if (length(state_runs) > 0) {
      duration_ms[s] <- mean(state_runs) * (1000 / sr)
      occurrence[s] <- length(state_runs) / total_sec
      coverage_pct[s] <- sum(state_runs) / n_time * 100
    }
  }

  # --- Transition probability matrix ---
  trans_mat <- matrix(0, nrow = n_states, ncol = n_states)
  for (i in seq_len(length(labels) - 1)) {
    from <- labels[i]
    to <- labels[i + 1]
    if (from != to) {
      trans_mat[from, to] <- trans_mat[from, to] + 1
    }
  }
  # Normalize by row sums
  row_totals <- rowSums(trans_mat)
  row_totals[row_totals == 0] <- 1  # avoid division by zero
  trans_mat <- trans_mat / row_totals

  result <- data.frame(
    state = seq_len(n_states),
    duration_ms = duration_ms,
    occurrence_per_sec = occurrence,
    coverage_pct = coverage_pct,
    stringsAsFactors = FALSE
  )
  attr(result, "transition_matrix") <- trans_mat

  result
}


#' Backfit Microstate Template Maps to EEG Data
#'
#' Assigns each time point to the microstate template map with the highest
#' absolute spatial correlation. This allows applying microstate maps derived
#' from one dataset to another dataset.
#'
#' @param x A PhysioExperiment object with EEG data.
#' @param maps Numeric matrix of microstate template maps (n_channels x
#'   n_states), as returned in \code{metadata(x)$microstates$maps}.
#' @param assay_name Name of the input assay. If \code{NULL}, the default
#'   assay is used.
#' @return Modified PhysioExperiment with updated microstate labels in
#'   \code{metadata(x)$microstates$labels}. Also stores the template maps
#'   in \code{metadata(x)$microstates$maps} and the number of states in
#'   \code{metadata(x)$microstates$n_states}.
#' @references
#' Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for studying
#' the temporal dynamics of whole-brain neuronal networks. NeuroImage, 180, 577-593.
#' @seealso [eegMicrostates()], [eegMicrostateStats()],
#'   [eegMicrostateSequence()]
#' @export
#' @examples
#' \dontrun{
#' pe1 <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
#' pe1 <- eegMicrostates(pe1, n_states = 4)
#' maps <- metadata(pe1)$microstates$maps
#'
#' pe2 <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
#' pe2 <- eegMicrostateBackfit(pe2, maps)
#' }
eegMicrostateBackfit <- function(x, maps, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.matrix(maps))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)

  labels <- .microstate_backfit(data, maps)

  if (is.null(S4Vectors::metadata(x)$microstates)) {
    S4Vectors::metadata(x)$microstates <- list()
  }
  S4Vectors::metadata(x)$microstates$labels <- labels
  S4Vectors::metadata(x)$microstates$maps <- maps
  S4Vectors::metadata(x)$microstates$n_states <- ncol(maps)

  x
}


#' Extract Microstate Sequence as Character Labels
#'
#' Converts integer microstate labels to character labels ("A", "B", "C", ...).
#'
#' @param x A PhysioExperiment object with microstate labels in
#'   \code{metadata(x)$microstates$labels}.
#' @return Character vector of length n_time with labels "A", "B", "C", etc.
#'   Each element corresponds to the microstate class assigned to that time
#'   point.
#' @references
#' Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for studying
#' the temporal dynamics of whole-brain neuronal networks. NeuroImage, 180, 577-593.
#' @seealso [eegMicrostates()], [eegMicrostateStats()],
#'   [eegMicrostateBackfit()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
#' pe <- eegMicrostates(pe, n_states = 4)
#' seq_labels <- eegMicrostateSequence(pe)
#' table(seq_labels)
#' }
eegMicrostateSequence <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  ms <- S4Vectors::metadata(x)$microstates
  if (is.null(ms) || is.null(ms$labels)) {
    stop("No microstate labels found. Run eegMicrostates() first.", call. = FALSE)
  }

  labels <- ms$labels
  LETTERS[labels]
}


# ===========================================================================
# Internal helper functions
# ===========================================================================

#' Polarity-invariant K-means for microstate clustering
#' @keywords internal
.microstate_kmeans <- function(topos, n_states, max_iter = 100) {
  n_peaks <- nrow(topos)
  n_channels <- ncol(topos)

  # Normalize topographies to unit norm
  topo_norms <- sqrt(rowSums(topos^2))
  topo_norms[topo_norms < 1e-10] <- 1e-10
  topos_norm <- topos / topo_norms

  # Initialize cluster centers from random peaks
  set.seed(42)
  init_idx <- sample(n_peaks, n_states)
  centers <- topos_norm[init_idx, , drop = FALSE]  # n_states x n_channels

  assignments <- integer(n_peaks)

  for (iter in seq_len(max_iter)) {
    old_assignments <- assignments

    # Assignment step: polarity-invariant (max |correlation|)
    for (i in seq_len(n_peaks)) {
      best_state <- 1L
      best_corr <- -Inf
      for (s in seq_len(n_states)) {
        corr <- abs(sum(topos_norm[i, ] * centers[s, ]))
        if (corr > best_corr) {
          best_corr <- corr
          best_state <- s
        }
      }
      assignments[i] <- best_state
    }

    # Update step: compute new centers accounting for polarity
    for (s in seq_len(n_states)) {
      members <- which(assignments == s)
      if (length(members) == 0) next

      # Align polarities to first member
      ref <- topos_norm[members[1], ]
      aligned <- matrix(0, nrow = length(members), ncol = n_channels)
      for (m in seq_along(members)) {
        topo_m <- topos_norm[members[m], ]
        if (sum(topo_m * ref) < 0) {
          aligned[m, ] <- -topo_m
        } else {
          aligned[m, ] <- topo_m
        }
      }
      new_center <- colMeans(aligned)
      cnorm <- sqrt(sum(new_center^2))
      if (cnorm > 1e-10) new_center <- new_center / cnorm
      centers[s, ] <- new_center
    }

    # Check convergence
    if (identical(assignments, old_assignments)) break
  }

  # Return maps as n_channels x n_states
  t(centers)
}


#' AAHC (Atomize and Agglomerate Hierarchical Clustering)
#' @keywords internal
.microstate_aahc <- function(topos, n_states) {
  n_peaks <- nrow(topos)
  n_channels <- ncol(topos)

  # Normalize topographies
  topo_norms <- sqrt(rowSums(topos^2))
  topo_norms[topo_norms < 1e-10] <- 1e-10
  topos_norm <- topos / topo_norms

  # Initialize: each peak is its own cluster
  clusters <- lapply(seq_len(n_peaks), function(i) topos_norm[i, , drop = FALSE])
  cluster_ids <- seq_len(n_peaks)

  while (length(clusters) > n_states) {
    n_clust <- length(clusters)

    # Compute centers for each cluster
    centers <- matrix(0, nrow = n_clust, ncol = n_channels)
    for (c_idx in seq_len(n_clust)) {
      center <- colMeans(clusters[[c_idx]])
      cnorm <- sqrt(sum(center^2))
      if (cnorm > 1e-10) center <- center / cnorm
      centers[c_idx, ] <- center
    }

    # Find most similar pair (highest absolute correlation)
    best_corr <- -Inf
    best_i <- 1L
    best_j <- 2L
    for (i in seq_len(n_clust - 1)) {
      for (j in (i + 1):n_clust) {
        corr <- abs(sum(centers[i, ] * centers[j, ]))
        if (corr > best_corr) {
          best_corr <- corr
          best_i <- i
          best_j <- j
        }
      }
    }

    # Merge: align polarities before combining
    members_i <- clusters[[best_i]]
    members_j <- clusters[[best_j]]
    ref_center <- centers[best_i, ]

    # Align cluster j members to cluster i polarity
    if (sum(centers[best_j, ] * ref_center) < 0) {
      members_j <- -members_j
    }

    merged <- rbind(members_i, members_j)
    clusters[[best_i]] <- merged
    clusters <- clusters[-best_j]
  }

  # Extract final cluster centers as maps
  maps <- matrix(0, nrow = n_channels, ncol = n_states)
  for (s in seq_len(n_states)) {
    center <- colMeans(clusters[[s]])
    cnorm <- sqrt(sum(center^2))
    if (cnorm > 1e-10) center <- center / cnorm
    maps[, s] <- center
  }

  maps
}


#' PCA-based microstate extraction
#' @keywords internal
.microstate_pca <- function(topos, n_states) {
  # PCA on peak topographies
  pca_result <- stats::prcomp(topos, center = TRUE, scale. = FALSE)

  # Use first n_states principal components as maps
  maps <- pca_result$rotation[, seq_len(n_states), drop = FALSE]

  # Normalize each map
  for (s in seq_len(n_states)) {
    cnorm <- sqrt(sum(maps[, s]^2))
    if (cnorm > 1e-10) maps[, s] <- maps[, s] / cnorm
  }

  maps  # n_channels x n_states
}


#' Backfit topographies to microstate maps
#' @keywords internal
.microstate_backfit <- function(data, maps) {
  n_time <- nrow(data)
  n_states <- ncol(maps)
  n_channels <- nrow(maps)

  labels <- integer(n_time)

  # Normalize maps
  map_norms <- sqrt(colSums(maps^2))
  map_norms[map_norms < 1e-10] <- 1e-10
  maps_norm <- sweep(maps, 2, map_norms, "/")

  for (i in seq_len(n_time)) {
    topo <- data[i, ]
    tnorm <- sqrt(sum(topo^2))
    if (tnorm < 1e-10) {
      labels[i] <- 1L
      next
    }
    topo_norm <- topo / tnorm

    best_state <- 1L
    best_corr <- -Inf
    for (s in seq_len(n_states)) {
      corr <- abs(sum(topo_norm * maps_norm[, s]))
      if (corr > best_corr) {
        best_corr <- corr
        best_state <- s
      }
    }
    labels[i] <- best_state
  }

  labels
}
