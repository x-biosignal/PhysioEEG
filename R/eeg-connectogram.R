#' Circular EEG connectogram
#'
#' Draws a deterministic circular view of a connectivity matrix stored in
#' \code{metadata(x)$connectivity$matrix}. Symmetric matrices are represented
#' by one edge per unordered channel pair. Asymmetric matrices retain every
#' ordered entry, following the PhysioEEG directed-matrix convention of rows
#' as targets and columns as sources.
#'
#' The threshold is a strict absolute display threshold
#' (\code{abs(value) > threshold}); it is not a p-value or a corrected
#' significance mask. Hemisphere and lobe are inferred only for recognized
#' 10-20/10-10 labels. Other labels remain \code{"unknown"} unless explicit
#' \code{colData} columns are supplied. Label tie-breaking uses a stable UTF-8
#' byte key. Edge identifiers length-prefix both endpoint labels so labels
#' containing the visible \code{" -- "} or \code{" -> "} separators cannot
#' collide.
#'
#' @param x A PhysioExperiment with a finite square numeric matrix at
#'   \code{metadata(x)$connectivity$matrix}.
#' @param order Exactly one node ordering: \code{"hemisphere"}, \code{"lobe"},
#'   or \code{"cluster"}. Cluster ordering uses average-linkage clustering of
#'   symmetrized absolute connectivity and deterministic merge orientation.
#' @param modules Optional module declaration. Supply a named character or
#'   factor vector, or a data frame with unique \code{label} and
#'   \code{module} columns covering every channel. If \code{NULL}, an exact
#'   \code{colData(x)$module} column is used when present, otherwise lobe.
#' @param threshold A finite non-negative scalar, or \code{NULL} for zero.
#' @param bundle One non-missing logical scalar. If \code{TRUE}, edges follow
#'   deterministic cubic Bezier paths through module hubs; otherwise they are
#'   straight.
#' @return A ggplot2 object. Resolved node, edge, path, module-arc, and setting
#'   tables are available in \code{attr(plot, "connectogram_data")}.
#' @seealso [eegPlotConnectivity()], [eegConnectivityMatrix()]
#' @export
#' @examples
#' \dontrun{
#' pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
#' mat <- matrix(c(
#'   1, 0.5, 0, -0.4,
#'   0.5, 1, 0.3, 0,
#'   0, 0.3, 1, 0.6,
#'   -0.4, 0, 0.6, 1
#' ), 4, 4)
#' labels <- SummarizedExperiment::colData(pe)$label
#' dimnames(mat) <- list(labels, labels)
#' metadata(pe)$connectivity <- list(
#'   matrix = mat, method = "coherence", band = c(8, 13)
#' )
#' eegPlotConnectogram(pe, threshold = 0.2)
#' }
eegPlotConnectogram <- function(x,
                                order = c("hemisphere", "lobe", "cluster"),
                                modules = NULL,
                                threshold = NULL,
                                bundle = TRUE) {
  .check_ggplot2()
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must inherit from PhysioExperiment.", call. = FALSE)
  }
  order_values <- c("hemisphere", "lobe", "cluster")
  if (missing(order)) {
    order <- order_values[[1L]]
  } else if (!is.character(order) || length(order) != 1L ||
             is.na(order) || !(order %in% order_values)) {
    stop("order must be exactly one of hemisphere, lobe, or cluster.",
         call. = FALSE)
  }
  if (is.null(threshold)) {
    threshold <- 0
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) || threshold < 0) {
    stop("threshold must be NULL or one finite non-negative scalar.",
         call. = FALSE)
  }
  if (!is.logical(bundle) || length(bundle) != 1L || is.na(bundle)) {
    stop("bundle must be one non-missing logical scalar.", call. = FALSE)
  }

  resolved <- .eeg_connectogram_resolve(
    x = x,
    order = order,
    modules = modules,
    threshold = as.numeric(threshold),
    bundle = bundle
  )
  paths <- resolved$paths
  arcs <- resolved$module_arcs
  nodes <- resolved$nodes

  p <- ggplot2::ggplot()
  if (nrow(paths) > 0L) {
    arrow <- if (isTRUE(resolved$settings$directed)) {
      grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")
    } else {
      NULL
    }
    edge_colors <- PhysioCore::physioPalette(256L, "diverging")
    edge_limit <- max(paths$abs_value)
    if (!is.finite(edge_limit) || edge_limit <= 0) {
      edge_limit <- 1
    }
    p <- p +
      ggplot2::geom_path(
        data = paths,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$edge_id,
          color = .data$value,
          linewidth = .data$abs_value,
          alpha = .data$abs_value
        ),
        lineend = "round",
        arrow = arrow
      ) +
      ggplot2::scale_color_gradient2(
        low = edge_colors[1L],
        mid = edge_colors[128L],
        high = edge_colors[256L],
        midpoint = 0,
        limits = c(-edge_limit, edge_limit),
        name = "Connectivity"
      ) +
      ggplot2::scale_linewidth_continuous(
        range = c(0.3, 2),
        limits = c(0, edge_limit),
        guide = "none"
      ) +
      ggplot2::scale_alpha_continuous(
        range = c(0.25, 0.85),
        limits = c(0, edge_limit),
        guide = "none"
      )
  }

  if (nrow(arcs) > 0L) {
    p <- p +
      ggplot2::geom_path(
        data = arcs,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$module_run_id
        ),
        color = arcs$color,
        linewidth = 3,
        lineend = "butt"
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = nodes,
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$y,
        fill = .data$module
      ),
      shape = 21,
      size = 4,
      stroke = 0.6,
      color = "#202020"
    ) +
    ggplot2::scale_fill_manual(
      values = resolved$settings$module_colors,
      name = "Module"
    ) +
    ggplot2::geom_text(
      data = nodes,
      mapping = ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$label,
        angle = .data$label_angle,
        hjust = .data$label_hjust
      ),
      size = 3
    ) +
    ggplot2::coord_fixed(
      xlim = resolved$settings$limits,
      ylim = resolved$settings$limits,
      clip = "off"
    ) +
    ggplot2::labs(
      title = "EEG Connectogram",
      subtitle = resolved$settings$subtitle,
      x = NULL,
      y = NULL
    ) +
    PhysioCore::theme_physio() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "right"
    )

  attr(p, "connectogram_data") <- resolved
  p
}


.eeg_connectogram_resolve <- function(x, order, modules, threshold, bundle) {
  conn <- S4Vectors::metadata(x)$connectivity
  if (is.null(conn) || is.null(conn$matrix)) {
    stop("metadata(x)$connectivity$matrix is required.", call. = FALSE)
  }
  mat <- conn$matrix
  if (!is.matrix(mat) || !is.numeric(mat) || length(dim(mat)) != 2L ||
      nrow(mat) != ncol(mat) || nrow(mat) < 2L) {
    stop("The connectivity matrix must be numeric, square, and at least 2 x 2.",
         call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop("The connectivity matrix must contain only finite values.",
         call. = FALSE)
  }

  identity <- .eeg_connectogram_identity(x, mat)
  labels <- identity$labels
  cd <- identity$col_data
  anatomy <- .eeg_connectogram_anatomy(labels, cd)
  module_values <- .eeg_connectogram_modules(
    modules = modules,
    labels = labels,
    col_data = cd,
    lobes = anatomy$lobe
  )

  scale <- max(1, max(abs(mat)))
  tolerance <- sqrt(.Machine$double.eps) * scale
  max_asymmetry <- max(abs(mat - t(mat)))
  symmetric <- max_asymmetry <= tolerance
  order_result <- .eeg_connectogram_order(
    mat = mat,
    labels = labels,
    anatomy = anatomy,
    order = order
  )
  diagnostics <- order_result$diagnostics

  nodes <- .eeg_connectogram_nodes(
    labels = labels,
    anatomy = anatomy,
    modules = module_values,
    input_order = order_result$input_order
  )
  module_result <- .eeg_connectogram_module_geometry(nodes)
  diagnostics <- .eeg_connectogram_bind_diagnostics(
    diagnostics,
    module_result$diagnostics
  )

  method <- .eeg_connectogram_method(conn$method)
  band_label <- .eeg_connectogram_band_label(conn$band)
  edges <- .eeg_connectogram_edges(
    mat = mat,
    labels = labels,
    threshold = threshold,
    symmetric = symmetric,
    method = method,
    band = band_label
  )
  paths <- .eeg_connectogram_paths(
    edges = edges,
    nodes = nodes,
    hubs = module_result$hubs,
    bundle = bundle
  )

  band <- conn$band
  subtitle_method <- if (is.na(method)) "method unavailable" else method
  subtitle <- paste0(
    subtitle_method,
    "; |value| > ",
    format(threshold, digits = 6, trim = TRUE),
    " (descriptive display threshold)"
  )
  settings <- list(
    order = order,
    threshold = threshold,
    bundle = bundle,
    symmetric = symmetric,
    directed = !symmetric,
    symmetry_scale = scale,
    symmetry_tolerance = tolerance,
    max_asymmetry = max_asymmetry,
    matrix_direction = if (symmetric) {
      "unordered"
    } else {
      "rows_target_columns_source"
    },
    start_angle = pi / 2,
    direction = "clockwise",
    node_radius = 1,
    label_radius = 1.20,
    module_arc_radius = 1.08,
    module_hub_radius = 0.35,
    path_samples = if (bundle) 40L else 2L,
    arc_samples = 32L,
    limits = c(-1.38, 1.38),
    method = method,
    band = band,
    band_label = band_label,
    module_colors = module_result$colors,
    module_hubs = module_result$hubs,
    diagnostics = diagnostics,
    subtitle = subtitle
  )

  list(
    nodes = nodes,
    edges = edges,
    paths = paths,
    module_arcs = module_result$arcs,
    settings = settings
  )
}


.eeg_connectogram_identity <- function(x, mat) {
  rn <- rownames(mat)
  cn <- colnames(mat)
  has_rn <- !is.null(rn)
  has_cn <- !is.null(cn)
  if (xor(has_rn, has_cn)) {
    stop("Connectivity matrix row and column names must both be present or both be absent.",
         call. = FALSE)
  }
  if (has_rn) {
    if (length(rn) != nrow(mat) || length(cn) != ncol(mat) ||
        anyNA(rn) || anyNA(cn) || any(!nzchar(rn)) || any(!nzchar(cn)) ||
        anyDuplicated(rn) || anyDuplicated(cn) || !identical(rn, cn)) {
      stop("Connectivity matrix row and column names must be complete, unique, and identical in the same order.",
           call. = FALSE)
    }
  }

  cd_raw <- SummarizedExperiment::colData(x)
  cd <- as.data.frame(cd_raw, stringsAsFactors = FALSE)
  cd_labels <- NULL
  if ("label" %in% names(cd)) {
    cd_labels <- as.character(cd$label)
    if (length(cd_labels) != nrow(mat) || anyNA(cd_labels) ||
        any(!nzchar(cd_labels)) || anyDuplicated(cd_labels)) {
      stop("colData(x)$label must contain one complete unique label per matrix channel.",
           call. = FALSE)
    }
  }

  if (has_rn) {
    labels <- rn
    if (!is.null(cd_labels)) {
      if (!setequal(labels, cd_labels)) {
        stop("Matrix names and colData(x)$label must cover exactly the same channels.",
             call. = FALSE)
      }
      cd <- cd[match(labels, cd_labels), , drop = FALSE]
    } else if (nrow(cd) == nrow(mat)) {
      cd <- cd[seq_len(nrow(mat)), , drop = FALSE]
    } else {
      cd <- data.frame(row.names = seq_len(nrow(mat)))
    }
  } else {
    if (is.null(cd_labels)) {
      stop("An unnamed connectivity matrix requires exact unique colData(x)$label values.",
           call. = FALSE)
    }
    labels <- cd_labels
  }

  rownames(cd) <- NULL
  list(labels = labels, col_data = cd)
}


.eeg_connectogram_infer_anatomy <- function(labels) {
  upper <- toupper(labels)
  pattern <- "^(FP|AF|FC|FT|CP|TP|PO|F|C|P|O|I|T)([0-9]+|Z)$"
  matched <- regexec(pattern, upper, perl = TRUE)
  pieces <- regmatches(upper, matched)

  hemisphere <- rep("unknown", length(labels))
  lobe <- rep("unknown", length(labels))
  lobe_map <- c(
    FP = "frontal", AF = "frontal", F = "frontal",
    FC = "central", C = "central", CP = "central",
    P = "parietal", PO = "parietal",
    O = "occipital", I = "occipital",
    FT = "temporal", T = "temporal", TP = "temporal"
  )
  for (i in seq_along(pieces)) {
    if (length(pieces[[i]]) != 3L) {
      next
    }
    prefix <- pieces[[i]][2L]
    suffix <- pieces[[i]][3L]
    lobe[i] <- unname(lobe_map[prefix])
    if (identical(suffix, "Z")) {
      hemisphere[i] <- "midline"
    } else {
      number <- suppressWarnings(as.numeric(suffix))
      if (is.finite(number)) {
        hemisphere[i] <- if (number %% 2 == 1) "left" else "right"
      }
    }
  }
  data.frame(hemisphere = hemisphere, lobe = lobe,
             stringsAsFactors = FALSE)
}


.eeg_connectogram_explicit_anatomy <- function(col_data, candidates, allowed,
                                               inferred, field) {
  present <- candidates[candidates %in% names(col_data)]
  if (length(present) == 0L) {
    return(inferred)
  }
  values <- as.character(col_data[[present[1L]]])
  if (length(values) != length(inferred)) {
    stop("Explicit ", field, " metadata must cover every matrix channel.",
         call. = FALSE)
  }
  normalized <- tolower(values)
  invalid <- !is.na(normalized) & !(normalized %in% allowed)
  if (any(invalid)) {
    stop("Unrecognized explicit ", field, " value(s): ",
         paste(sort(unique(values[invalid])), collapse = ", "), ".",
         call. = FALSE)
  }
  use <- !is.na(normalized)
  inferred[use] <- normalized[use]
  inferred
}


.eeg_connectogram_anatomy <- function(labels, col_data) {
  inferred <- .eeg_connectogram_infer_anatomy(labels)
  hemisphere <- .eeg_connectogram_explicit_anatomy(
    col_data = col_data,
    candidates = c("hemisphere", "hemi"),
    allowed = c("left", "midline", "right", "unknown"),
    inferred = inferred$hemisphere,
    field = "hemisphere"
  )
  lobe <- .eeg_connectogram_explicit_anatomy(
    col_data = col_data,
    candidates = "lobe",
    allowed = c("frontal", "central", "parietal", "occipital", "temporal",
                "unknown"),
    inferred = inferred$lobe,
    field = "lobe"
  )
  data.frame(hemisphere = hemisphere, lobe = lobe,
             stringsAsFactors = FALSE)
}


.eeg_connectogram_validate_modules <- function(values, labels) {
  values <- as.character(values)
  if (length(values) != length(labels) || anyNA(values) ||
      any(!nzchar(values))) {
    stop("Modules must provide one non-empty value for every channel.",
         call. = FALSE)
  }
  unname(values)
}


.eeg_connectogram_modules <- function(modules, labels, col_data, lobes) {
  if (is.null(modules)) {
    if ("module" %in% names(col_data)) {
      return(.eeg_connectogram_validate_modules(col_data$module, labels))
    }
    return(.eeg_connectogram_validate_modules(lobes, labels))
  }

  if (is.factor(modules) || is.character(modules)) {
    module_names <- names(modules)
    if (is.null(module_names) || length(module_names) != length(labels) ||
        anyNA(module_names) || any(!nzchar(module_names)) ||
        anyDuplicated(module_names) || !setequal(module_names, labels)) {
      stop("A module vector must be uniquely named by exactly the channel labels.",
           call. = FALSE)
    }
    values <- as.character(modules)[match(labels, module_names)]
    return(.eeg_connectogram_validate_modules(values, labels))
  }

  if (is.data.frame(modules)) {
    if (!all(c("label", "module") %in% names(modules))) {
      stop("A module data frame must contain label and module columns.",
           call. = FALSE)
    }
    module_labels <- as.character(modules$label)
    if (length(module_labels) != length(labels) || anyNA(module_labels) ||
        any(!nzchar(module_labels)) || anyDuplicated(module_labels) ||
        !setequal(module_labels, labels)) {
      stop("A module data frame must cover every channel exactly once.",
           call. = FALSE)
    }
    values <- as.character(modules$module)[match(labels, module_labels)]
    return(.eeg_connectogram_validate_modules(values, labels))
  }

  stop("modules must be NULL, a named character/factor vector, or a data frame.",
       call. = FALSE)
}


.eeg_connectogram_order <- function(mat, labels, anatomy, order) {
  hemisphere_levels <- c("left", "midline", "right", "unknown")
  lobe_levels <- c("frontal", "central", "parietal", "occipital", "temporal",
                   "unknown")
  diagnostics <- .eeg_connectogram_empty_diagnostics()
  if (order == "hemisphere") {
    idx <- base::order(
      match(anatomy$hemisphere, hemisphere_levels),
      match(anatomy$lobe, lobe_levels),
      .eeg_connectogram_label_key(labels),
      method = "radix"
    )
  } else if (order == "lobe") {
    idx <- base::order(
      match(anatomy$lobe, lobe_levels),
      match(anatomy$hemisphere, hemisphere_levels),
      .eeg_connectogram_label_key(labels),
      method = "radix"
    )
  } else {
    canonical <- .eeg_connectogram_label_order(labels)
    canonical_labels <- labels[canonical]
    canonical_mat <- mat[canonical, canonical, drop = FALSE]
    similarity <- (abs(canonical_mat) + t(abs(canonical_mat))) / 2
    off_diagonal <- similarity[row(similarity) != col(similarity)]
    tie_tolerance <- sqrt(.Machine$double.eps) *
      max(1, max(abs(off_diagonal)))
    if (max(off_diagonal) - min(off_diagonal) <= tie_tolerance) {
      idx <- canonical
      diagnostics <- .eeg_connectogram_diagnostic(
        "cluster_all_tied",
        "All off-diagonal similarities are tied; channel-label order was used.",
        "cluster"
      )
    } else {
      low <- min(off_diagonal)
      high <- max(off_diagonal)
      similarity <- (similarity - low) / (high - low)
      similarity[similarity < 0] <- 0
      similarity[similarity > 1] <- 1
      diag(similarity) <- 1
      hc <- stats::hclust(stats::as.dist(1 - similarity), method = "average")
      leaf_order <- .eeg_connectogram_orient_hclust(
        hc$merge,
        canonical_labels
      )
      idx <- canonical[leaf_order]
    }
  }
  list(input_order = idx, diagnostics = diagnostics)
}


.eeg_connectogram_orient_hclust <- function(merge, labels) {
  label_keys <- .eeg_connectogram_label_key(labels)
  visit <- function(node) {
    if (node < 0L) {
      return(-node)
    }
    left <- visit(merge[node, 1L])
    right <- visit(merge[node, 2L])
    if (sort(label_keys[left], method = "radix")[1L] >
        sort(label_keys[right], method = "radix")[1L]) {
      tmp <- left
      left <- right
      right <- tmp
    }
    c(left, right)
  }
  visit(nrow(merge))
}


.eeg_connectogram_nodes <- function(labels, anatomy, modules, input_order) {
  n <- length(labels)
  plot_index <- seq_len(n)
  angle <- pi / 2 - 2 * pi * (plot_index - 1L) / n
  original_index <- input_order
  x <- cos(angle)
  y <- sin(angle)
  label_x <- 1.20 * x
  label_y <- 1.20 * y
  angle_degrees <- angle * 180 / pi
  left_side <- x < 0
  label_angle <- angle_degrees
  label_angle[left_side] <- label_angle[left_side] + 180

  data.frame(
    node_id = labels[original_index],
    input_index = original_index,
    plot_index = plot_index,
    label = labels[original_index],
    hemisphere = anatomy$hemisphere[original_index],
    lobe = anatomy$lobe[original_index],
    module = modules[original_index],
    angle = angle,
    x = x,
    y = y,
    label_x = label_x,
    label_y = label_y,
    label_angle = label_angle,
    label_hjust = ifelse(left_side, 1, 0),
    stringsAsFactors = FALSE
  )
}


.eeg_connectogram_module_geometry <- function(nodes) {
  module_names <- unique(nodes$module)
  module_names <- module_names[
    .eeg_connectogram_label_order(module_names)
  ]
  palette <- suppressWarnings(
    PhysioCore::physioPalette(length(module_names), "qualitative")
  )
  colors <- stats::setNames(palette, module_names)
  hubs <- vector("list", length(module_names))
  diagnostics <- .eeg_connectogram_empty_diagnostics()
  for (i in seq_along(module_names)) {
    module <- module_names[i]
    module_nodes <- nodes[nodes$module == module, , drop = FALSE]
    cx <- mean(cos(module_nodes$angle))
    cy <- mean(sin(module_nodes$angle))
    norm <- sqrt(cx^2 + cy^2)
    if (norm <= sqrt(.Machine$double.eps)) {
      hub_angle <- module_nodes$angle[which.min(module_nodes$plot_index)]
      diagnostics <- .eeg_connectogram_bind_diagnostics(
        diagnostics,
        .eeg_connectogram_diagnostic(
          "antipodal_module_hub",
          "The circular mean was undefined; the first stable node angle was used.",
          module
        )
      )
    } else {
      hub_angle <- atan2(cy, cx)
    }
    hubs[[i]] <- data.frame(
      module = module,
      angle = hub_angle,
      x = 0.35 * cos(hub_angle),
      y = 0.35 * sin(hub_angle),
      color = unname(colors[module]),
      stringsAsFactors = FALSE
    )
  }
  hubs <- do.call(rbind, hubs)
  rownames(hubs) <- NULL
  arcs <- .eeg_connectogram_module_arcs(nodes, colors)
  list(hubs = hubs, arcs = arcs, colors = colors,
       diagnostics = diagnostics)
}


.eeg_connectogram_circular_runs <- function(modules) {
  n <- length(modules)
  previous <- modules[c(n, seq_len(n - 1L))]
  starts <- which(modules != previous)
  if (length(starts) == 0L) {
    return(list(seq_len(n)))
  }
  starts <- sort(starts)
  lapply(seq_along(starts), function(i) {
    start <- starts[i]
    next_start <- starts[if (i == length(starts)) 1L else i + 1L]
    end <- if (i == length(starts)) next_start - 1L + n else next_start - 1L
    ((seq.int(start, end) - 1L) %% n) + 1L
  })
}


.eeg_connectogram_module_arcs <- function(nodes, colors) {
  runs <- .eeg_connectogram_circular_runs(nodes$module)
  counters <- integer()
  names(counters) <- character()
  output <- vector("list", length(runs))
  step <- 2 * pi / nrow(nodes)
  for (i in seq_along(runs)) {
    indices <- runs[[i]]
    module <- nodes$module[indices[1L]]
    if (!(module %in% names(counters))) {
      counters[module] <- 0L
    }
    counters[module] <- counters[module] + 1L
    run_id <- paste0(module, "::", counters[module])
    start_angle <- nodes$angle[indices[1L]] + step / 2
    end_angle <- nodes$angle[indices[length(indices)]] - step / 2
    while (end_angle >= start_angle) {
      end_angle <- end_angle - 2 * pi
    }
    theta <- seq(start_angle, end_angle, length.out = 32L)
    output[[i]] <- data.frame(
      module = module,
      module_run_id = run_id,
      start_plot_index = indices[1L],
      end_plot_index = indices[length(indices)],
      start_angle = start_angle,
      end_angle = end_angle,
      radius = 1.08,
      color = unname(colors[module]),
      arc_index = seq_along(theta),
      x = 1.08 * cos(theta),
      y = 1.08 * sin(theta),
      stringsAsFactors = FALSE
    )
  }
  arcs <- do.call(rbind, output)
  rownames(arcs) <- NULL
  arcs
}


.eeg_connectogram_method <- function(method) {
  if (is.null(method) || length(method) != 1L || is.na(method) ||
      !nzchar(as.character(method))) {
    return(NA_character_)
  }
  as.character(method)
}


.eeg_connectogram_band_label <- function(band) {
  if (is.null(band) || length(band) == 0L) {
    return(NA_character_)
  }
  if (is.numeric(band) && all(is.finite(band))) {
    return(paste(format(band, trim = TRUE, digits = 8), collapse = "-"))
  }
  if (is.character(band) && !anyNA(band) && all(nzchar(band))) {
    return(paste(band, collapse = "-"))
  }
  NA_character_
}


.eeg_connectogram_empty_edges <- function() {
  data.frame(
    edge_id = character(),
    source_id = character(),
    target_id = character(),
    source_index = integer(),
    target_index = integer(),
    value = numeric(),
    abs_value = numeric(),
    sign = integer(),
    directed = logical(),
    threshold = numeric(),
    method = character(),
    band = character(),
    stringsAsFactors = FALSE
  )
}


.eeg_connectogram_edges <- function(mat, labels, threshold, symmetric,
                                    method, band) {
  canonical_labels <- labels[.eeg_connectogram_label_order(labels)]
  input_index <- stats::setNames(seq_along(labels), labels)
  rows <- list()
  if (symmetric) {
    for (i in seq_len(length(canonical_labels) - 1L)) {
      for (j in seq.int(i + 1L, length(canonical_labels))) {
        source <- canonical_labels[i]
        target <- canonical_labels[j]
        value <- mat[input_index[source], input_index[target]]
        if (abs(value) > threshold) {
          rows[[length(rows) + 1L]] <- data.frame(
            edge_id = .eeg_connectogram_edge_id(
              source,
              target,
              directed = FALSE
            ),
            source_id = source,
            target_id = target,
            source_index = unname(input_index[source]),
            target_index = unname(input_index[target]),
            value = value,
            abs_value = abs(value),
            sign = as.integer(sign(value)),
            directed = FALSE,
            threshold = threshold,
            method = method,
            band = band,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  } else {
    for (source in canonical_labels) {
      for (target in canonical_labels) {
        if (identical(source, target)) {
          next
        }
        source_index <- unname(input_index[source])
        target_index <- unname(input_index[target])
        value <- mat[target_index, source_index]
        if (abs(value) > threshold) {
          rows[[length(rows) + 1L]] <- data.frame(
            edge_id = .eeg_connectogram_edge_id(
              source,
              target,
              directed = TRUE
            ),
            source_id = source,
            target_id = target,
            source_index = source_index,
            target_index = target_index,
            value = value,
            abs_value = abs(value),
            sign = as.integer(sign(value)),
            directed = TRUE,
            threshold = threshold,
            method = method,
            band = band,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(rows) == 0L) {
    return(.eeg_connectogram_empty_edges())
  }
  edges <- do.call(rbind, rows)
  rownames(edges) <- NULL
  edges
}


.eeg_connectogram_empty_paths <- function() {
  data.frame(
    edge_id = character(),
    path_index = integer(),
    x = numeric(),
    y = numeric(),
    source_id = character(),
    target_id = character(),
    value = numeric(),
    abs_value = numeric(),
    sign = integer(),
    directed = logical(),
    bundle = logical(),
    source_module = character(),
    target_module = character(),
    stringsAsFactors = FALSE
  )
}


.eeg_connectogram_paths <- function(edges, nodes, hubs, bundle) {
  if (nrow(edges) == 0L) {
    return(.eeg_connectogram_empty_paths())
  }
  node_index <- stats::setNames(seq_len(nrow(nodes)), nodes$node_id)
  hub_index <- stats::setNames(seq_len(nrow(hubs)), hubs$module)
  output <- vector("list", nrow(edges))
  for (i in seq_len(nrow(edges))) {
    edge <- edges[i, , drop = FALSE]
    source_node <- nodes[node_index[edge$source_id], , drop = FALSE]
    target_node <- nodes[node_index[edge$target_id], , drop = FALSE]
    source_module <- source_node$module
    target_module <- target_node$module
    if (bundle) {
      source_hub <- hubs[hub_index[source_module], , drop = FALSE]
      target_hub <- hubs[hub_index[target_module], , drop = FALSE]
      t <- seq(0, 1, length.out = 40L)
      omt <- 1 - t
      x <- omt^3 * source_node$x +
        3 * omt^2 * t * source_hub$x +
        3 * omt * t^2 * target_hub$x +
        t^3 * target_node$x
      y <- omt^3 * source_node$y +
        3 * omt^2 * t * source_hub$y +
        3 * omt * t^2 * target_hub$y +
        t^3 * target_node$y
    } else {
      x <- c(source_node$x, target_node$x)
      y <- c(source_node$y, target_node$y)
    }
    output[[i]] <- data.frame(
      edge_id = edge$edge_id,
      path_index = seq_along(x),
      x = x,
      y = y,
      source_id = edge$source_id,
      target_id = edge$target_id,
      value = edge$value,
      abs_value = edge$abs_value,
      sign = edge$sign,
      directed = edge$directed,
      bundle = bundle,
      source_module = source_module,
      target_module = target_module,
      stringsAsFactors = FALSE
    )
  }
  paths <- do.call(rbind, output)
  rownames(paths) <- NULL
  paths
}


.eeg_connectogram_empty_diagnostics <- function() {
  data.frame(code = character(), message = character(), context = character(),
             stringsAsFactors = FALSE)
}


.eeg_connectogram_diagnostic <- function(code, message, context) {
  data.frame(code = code, message = message, context = context,
             stringsAsFactors = FALSE)
}


.eeg_connectogram_bind_diagnostics <- function(x, y) {
  if (nrow(x) == 0L) {
    return(y)
  }
  if (nrow(y) == 0L) {
    return(x)
  }
  rbind(x, y)
}


.eeg_connectogram_label_key <- function(labels) {
  vapply(enc2utf8(labels), function(label) {
    paste(sprintf("%02x", as.integer(charToRaw(label))), collapse = "")
  }, character(1))
}


.eeg_connectogram_label_order <- function(labels) {
  base::order(.eeg_connectogram_label_key(labels), method = "radix")
}


.eeg_connectogram_edge_id <- function(source, target, directed) {
  source_utf8 <- enc2utf8(source)
  target_utf8 <- enc2utf8(target)
  separator <- if (directed) " -> " else " -- "
  paste0(
    nchar(source_utf8, type = "bytes"),
    ":",
    source_utf8,
    separator,
    nchar(target_utf8, type = "bytes"),
    ":",
    target_utf8
  )
}
