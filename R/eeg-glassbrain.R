.glassbrain_view_axes <- list(
  sagittal = c(projected_x = "y", projected_y = "z", depth = "x"),
  axial = c(projected_x = "x", projected_y = "y", depth = "z"),
  coronal = c(projected_x = "x", projected_y = "z", depth = "y")
)

.glassbrain_source_ids <- function(source_id, n_sources) {
  if (is.null(source_id)) {
    return(seq_len(n_sources))
  }
  if (length(source_id) != n_sources || anyNA(source_id)) {
    stop(
      "'source_id' must have one non-missing value per source.",
      call. = FALSE
    )
  }
  if (is.factor(source_id)) {
    source_id <- as.character(source_id)
  }
  if (is.numeric(source_id)) {
    if (any(!is.finite(source_id)) || anyDuplicated(source_id)) {
      stop("'source_id' values must be finite and unique.", call. = FALSE)
    }
    return(source_id)
  }
  if (!is.character(source_id)) {
    stop("'source_id' must be numeric, character, or factor.", call. = FALSE)
  }
  if (any(!nzchar(source_id)) || anyDuplicated(source_id)) {
    stop("'source_id' values must be non-empty and unique.", call. = FALSE)
  }
  source_id
}

.glassbrain_positions <- function(positions) {
  if (is.data.frame(positions)) {
    if (!all(c("x", "y", "z") %in% names(positions))) {
      stop(
        "Source positions must contain named x, y, and z columns.",
        call. = FALSE
      )
    }
    result <- positions[, c("x", "y", "z"), drop = FALSE]
  } else if (is.matrix(positions) && is.numeric(positions)) {
    if (ncol(positions) != 3L) {
      stop("Source positions must have exactly three columns.", call. = FALSE)
    }
    if (!is.null(colnames(positions)) &&
        !identical(colnames(positions), c("x", "y", "z"))) {
      stop(
        "Named source-position columns must be exactly x, y, and z.",
        call. = FALSE
      )
    }
    result <- data.frame(
      x = positions[, 1L],
      y = positions[, 2L],
      z = positions[, 3L],
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "Source positions must be a data frame or numeric matrix with x, y, z.",
      call. = FALSE
    )
  }
  if (nrow(result) < 1L ||
      any(!vapply(result, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(result)))) {
    stop("Source coordinates must be non-empty and finite numeric.", call. = FALSE)
  }
  result
}

.glassbrain_reduce_matrix <- function(values, n_sources, orientation_count,
                                      matrix_layout, reduction) {
  if (!is.matrix(values) || !is.numeric(values) ||
      any(!is.finite(values))) {
    stop("'source_matrix' must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(orientation_count) || length(orientation_count) != 1L ||
      is.na(orientation_count) || !is.finite(orientation_count) ||
      orientation_count < 1L || orientation_count != floor(orientation_count)) {
    stop("'orientation_count' must be one positive integer.", call. = FALSE)
  }
  orientation_count <- as.integer(orientation_count)
  matrix_layout <- match.arg(
    matrix_layout,
    c("source_by_orientation", "time_by_source_orientation")
  )
  if (matrix_layout == "source_by_orientation") {
    if (!identical(dim(values), c(n_sources, orientation_count)) ||
        !identical(reduction, "orientation_rms")) {
      stop(
        paste0(
          "source_by_orientation requires n_source x orientation_count ",
          "values and reduction = 'orientation_rms'."
        ),
        call. = FALSE
      )
    }
    return(list(
      amplitude = sqrt(rowSums(values^2)),
      reduction = "orientation_rms"
    ))
  }
  if (nrow(values) < 1L ||
      ncol(values) != n_sources * orientation_count ||
      !identical(reduction, "time_orientation_rms")) {
    stop(
      paste0(
        "time_by_source_orientation requires time x ",
        "(n_source * orientation_count) values and reduction = ",
        "'time_orientation_rms'."
      ),
      call. = FALSE
    )
  }
  amplitude <- vapply(seq_len(n_sources), function(source_index) {
    columns <- ((source_index - 1L) * orientation_count + 1L):(
      source_index * orientation_count
    )
    root_sum_square <- sqrt(rowSums(values[, columns, drop = FALSE]^2))
    sqrt(mean(root_sum_square^2))
  }, numeric(1))
  list(amplitude = amplitude, reduction = "time_orientation_rms")
}

.glassbrain_source_frame <- function(positions, amplitude, source_id = NULL) {
  positions <- .glassbrain_positions(positions)
  if (!is.numeric(amplitude) || length(amplitude) != nrow(positions) ||
      any(!is.finite(amplitude))) {
    stop(
      "Source amplitude must be one finite numeric value per position.",
      call. = FALSE
    )
  }
  source_id <- .glassbrain_source_ids(source_id, nrow(positions))
  data.frame(
    source_id = source_id,
    x = positions$x,
    y = positions$y,
    z = positions$z,
    amplitude = as.numeric(amplitude),
    stringsAsFactors = FALSE
  )
}

.glassbrain_explicit_source <- function(source_data) {
  if (is.data.frame(source_data)) {
    required <- c("x", "y", "z", "amplitude")
    if (!all(required %in% names(source_data))) {
      stop(
        "Explicit source data must contain x, y, z, and amplitude columns.",
        call. = FALSE
      )
    }
    source_id <- if ("source_id" %in% names(source_data)) {
      source_data$source_id
    } else {
      NULL
    }
    return(list(
      sources = .glassbrain_source_frame(
        source_data[, c("x", "y", "z"), drop = FALSE],
        source_data$amplitude,
        source_id
      ),
      source = "explicit",
      reduction = "none",
      method = NA_character_,
      coordinates = list(
        source = "explicit source_data",
        units = "caller declared or unspecified",
        anatomical_fidelity = "caller responsibility"
      ),
      warnings = character()
    ))
  }
  if (!is.list(source_data) || is.null(source_data$positions)) {
    stop(
      paste0(
        "source_data must be a data frame with x/y/z/amplitude or a ",
        "structured list with positions and declared values."
      ),
      call. = FALSE
    )
  }
  positions <- .glassbrain_positions(source_data$positions)
  source_id <- if (!is.null(source_data$source_id)) {
    source_data$source_id
  } else if (is.data.frame(source_data$positions) &&
             "source_id" %in% names(source_data$positions)) {
    source_data$positions$source_id
  } else {
    NULL
  }
  if (!is.null(source_data$amplitude)) {
    reduced <- list(
      amplitude = source_data$amplitude,
      reduction = "none"
    )
  } else if (!is.null(source_data$source_matrix)) {
    if (is.null(source_data$orientation_count) ||
        is.null(source_data$matrix_layout) ||
        is.null(source_data$reduction)) {
      stop(
        paste0(
          "A source_matrix requires orientation_count, matrix_layout, ",
          "and reduction declarations."
        ),
        call. = FALSE
      )
    }
    reduced <- .glassbrain_reduce_matrix(
      source_data$source_matrix,
      nrow(positions),
      source_data$orientation_count,
      source_data$matrix_layout,
      source_data$reduction
    )
  } else {
    stop(
      "Structured source_data requires 'amplitude' or 'source_matrix'.",
      call. = FALSE
    )
  }
  coordinates <- source_data$coordinate_provenance
  if (is.null(coordinates)) {
    coordinates <- list(
      source = "explicit structured source_data",
      units = "caller declared or unspecified",
      anatomical_fidelity = "caller responsibility"
    )
  }
  list(
    sources = .glassbrain_source_frame(
      positions,
      reduced$amplitude,
      source_id
    ),
    source = "explicit",
    reduction = reduced$reduction,
    method = if (is.null(source_data$method)) {
      NA_character_
    } else {
      as.character(source_data$method)[1L]
    },
    coordinates = coordinates,
    warnings = character()
  )
}

.glassbrain_metadata_source <- function(x) {
  metadata <- S4Vectors::metadata(x)
  source_info <- metadata$source_estimate
  beamformer_info <- metadata$beamformer_info
  has_source_info <- is.list(source_info) &&
    any(c("source_positions", "output_assay", "orientation_count") %in%
          names(source_info))
  has_beamformer_info <- is.list(beamformer_info)
  available <- c(
    if (has_source_info) {
      "source_estimate"
    },
    if (has_beamformer_info) {
      "beamformer"
    }
  )
  selected <- metadata$source_plot_default
  if (!is.null(selected) &&
      (!is.character(selected) || length(selected) != 1L ||
       is.na(selected) ||
       !selected %in% c("source_estimate", "beamformer"))) {
    stop(
      paste0(
        "Source plotting metadata contain an invalid default; provide ",
        "explicit source_data."
      ),
      call. = FALSE
    )
  }
  if (is.null(selected)) {
    if (length(available) == 1L) {
      selected <- available
    } else if (length(available) > 1L) {
      stop(
        paste0(
          "Multiple source products are available without a valid default; ",
          "provide explicit source_data."
        ),
        call. = FALSE
      )
    } else {
      if (!is.null(source_info) || !is.null(metadata$beamformer)) {
        stop(
          paste0(
            "Legacy source metadata do not contain 3D source positions. ",
            "Provide source_data with x, y, z, and amplitude columns."
          ),
          call. = FALSE
        )
      }
      stop(
        paste0(
          "No plottable source result found. Run eegSourceEstimate() or ",
          "eegBeamformer() with a forward model, or provide explicit ",
          "x/y/z/amplitude source_data."
        ),
        call. = FALSE
      )
    }
  } else if (!selected %in% available) {
    stop(
      sprintf(
        paste0(
          "The selected %s plotting metadata are incomplete; ",
          "provide explicit source_data."
        ),
        selected
      ),
      call. = FALSE
    )
  }

  if (identical(selected, "source_estimate")) {
    info <- source_info
    required <- c(
      "source_positions",
      "output_assay",
      "orientation_count",
      "method",
      "coordinate_provenance"
    )
    if (!all(required %in% names(info)) ||
        !is.character(info$output_assay) ||
        length(info$output_assay) != 1L ||
        is.na(info$output_assay) ||
        !nzchar(info$output_assay) ||
        !info$output_assay %in% names(metadata)) {
      stop(
        "Source-estimate plotting metadata are incomplete or invalid.",
        call. = FALSE
      )
    }
    values <- metadata[[info$output_assay]]
    positions <- .glassbrain_positions(info$source_positions)
    reduced <- .glassbrain_reduce_matrix(
      values,
      nrow(positions),
      info$orientation_count,
      "time_by_source_orientation",
      "time_orientation_rms"
    )
    return(list(
      sources = .glassbrain_source_frame(
        positions,
        reduced$amplitude
      ),
      source = "source_estimate",
      reduction = reduced$reduction,
      method = as.character(info$method)[1L],
      coordinates = info$coordinate_provenance,
      warnings = character()
    ))
  }

  info <- beamformer_info
  required <- c(
    "source_positions",
    "output_assay",
    "method",
    "coordinate_provenance"
  )
  if (!all(required %in% names(info)) ||
      !is.character(info$output_assay) ||
      length(info$output_assay) != 1L ||
      is.na(info$output_assay) ||
      !nzchar(info$output_assay) ||
      !info$output_assay %in% names(metadata)) {
    stop(
      "Beamformer plotting metadata are incomplete or invalid.",
      call. = FALSE
    )
  }
  values <- metadata[[info$output_assay]]
  positions <- .glassbrain_positions(info$source_positions)
  if (is.matrix(values) && ncol(values) == 1L) {
    amplitude <- values[, 1L]
  } else if (is.numeric(values) && is.null(dim(values))) {
    amplitude <- values
  } else {
    stop(
      "Beamformer source power must be one numeric value per source.",
      call. = FALSE
    )
  }
  list(
    sources = .glassbrain_source_frame(positions, amplitude),
    source = "beamformer",
    reduction = "none",
    method = as.character(info$method)[1L],
    coordinates = info$coordinate_provenance,
    warnings = character()
  )
}

.glassbrain_resolve_sources <- function(x, source_data = NULL) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("'x' must inherit from PhysioExperiment.", call. = FALSE)
  }
  if (is.null(source_data)) {
    .glassbrain_metadata_source(x)
  } else {
    .glassbrain_explicit_source(source_data)
  }
}

.glassbrain_transform <- function(sources, views) {
  center <- vapply(
    sources[c("x", "y", "z")],
    function(axis) (min(axis) + max(axis)) / 2,
    numeric(1)
  )
  ranges <- vapply(
    sources[c("x", "y", "z")],
    function(axis) diff(range(axis)),
    numeric(1)
  )
  for (view in views) {
    axes <- .glassbrain_view_axes[[view]][c("projected_x", "projected_y")]
    if (any(ranges[axes] <= 0)) {
      stop(
        sprintf(
          "Source coordinates are degenerate for the %s projection.",
          view
        ),
        call. = FALSE
      )
    }
  }
  scale <- max(ranges) / 2
  if (!is.finite(scale) || scale <= 0) {
    stop("Source coordinate ranges are degenerate.", call. = FALSE)
  }
  normalized <- sources
  for (axis in c("x", "y", "z")) {
    normalized[[axis]] <- (sources[[axis]] - center[[axis]]) / scale
  }
  list(
    sources = normalized,
    center = center,
    scale = scale,
    ranges = ranges
  )
}

.glassbrain_project <- function(sources, view, threshold, source_count,
                                grid_resolution = 80L) {
  axes <- .glassbrain_view_axes[[view]]
  projected <- data.frame(
    source_id = sources$source_id,
    amplitude = sources$amplitude,
    intensity = abs(sources$amplitude),
    projected_x = sources[[axes[["projected_x"]]]],
    projected_y = sources[[axes[["projected_y"]]]],
    depth = sources[[axes[["depth"]]]],
    view = view,
    threshold = threshold,
    source_count = source_count,
    stringsAsFactors = FALSE
  )
  breaks <- seq(-1, 1, length.out = grid_resolution + 1L)
  projected$x_bin <- findInterval(
    pmax(-1, pmin(1, projected$projected_x)),
    breaks,
    all.inside = TRUE
  )
  projected$y_bin <- findInterval(
    pmax(-1, pmin(1, projected$projected_y)),
    breaks,
    all.inside = TRUE
  )
  key <- paste(projected$x_bin, projected$y_bin, sep = ":")
  keep <- logical(nrow(projected))
  for (cell in unique(key)) {
    candidates <- which(key == cell)
    maximum <- max(projected$intensity[candidates])
    candidates <- candidates[
      projected$intensity[candidates] == maximum
    ]
    if (length(candidates) > 1L) {
      candidates <- candidates[
        order(projected$source_id[candidates], method = "radix")
      ]
    }
    keep[candidates[1L]] <- TRUE
  }
  projected <- projected[keep, , drop = FALSE]
  projected <- projected[
    order(projected$y_bin, projected$x_bin, projected$source_id),
    ,
    drop = FALSE
  ]
  rownames(projected) <- NULL
  projected$projected_count <- nrow(projected)
  projected
}

.glassbrain_outline <- function() {
  path <- system.file(
    "extdata",
    "glassbrain_outline.rds",
    package = "PhysioEEG"
  )
  if (!nzchar(path)) {
    stop("Bundled glass-brain outline is unavailable.", call. = FALSE)
  }
  outline <- readRDS(path)
  required <- c("view", "path_id", "vertex_order", "x", "y", "structure")
  valid <- is.data.frame(outline) &&
    identical(names(outline), required) &&
    identical(
      unique(outline$view),
      c("sagittal", "axial", "coronal")
    ) &&
    is.numeric(outline$x) &&
    is.numeric(outline$y) &&
    all(is.finite(outline$x)) &&
    all(is.finite(outline$y)) &&
    max(abs(c(outline$x, outline$y))) <= 1
  if (valid) {
    paths <- split(outline, outline$path_id)
    valid <- all(vapply(paths, function(path_data) {
      identical(path_data$vertex_order, seq_len(nrow(path_data))) &&
        identical(
          unname(as.numeric(path_data[1L, c("x", "y")])),
          unname(as.numeric(
            path_data[nrow(path_data), c("x", "y")]
          ))
        )
    }, logical(1)))
  }
  if (!valid) {
    stop("Bundled glass-brain outline failed schema validation.", call. = FALSE)
  }
  outline
}

.glassbrain_axis_labels <- function(view) {
  switch(
    view,
    sagittal = c(x = "Y (anterior-posterior)", y = "Z (inferior-superior)"),
    axial = c(x = "X (left-right)", y = "Y (anterior-posterior)"),
    coronal = c(x = "X (left-right)", y = "Z (inferior-superior)")
  )
}

.eeg_plot_glassbrain <- function(x, source_data = NULL,
                                 views = c("sagittal", "axial", "coronal"),
  threshold_pct = 90) {
  .check_ggplot2()
  if (!is.character(views) || length(views) < 1L || anyNA(views) ||
      any(!nzchar(views)) || anyDuplicated(views) ||
      any(!views %in% names(.glassbrain_view_axes))) {
    stop(
      paste0(
        "'views' must contain unique exact values from sagittal, axial, ",
        "and coronal."
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(threshold_pct) || length(threshold_pct) != 1L ||
      is.na(threshold_pct) || !is.finite(threshold_pct) ||
      threshold_pct < 0 || threshold_pct > 100) {
    stop(
      "'threshold_pct' must be one finite number in [0, 100].",
      call. = FALSE
    )
  }
  resolved <- .glassbrain_resolve_sources(x, source_data)
  sources <- resolved$sources
  threshold <- as.numeric(stats::quantile(
    abs(sources$amplitude),
    threshold_pct / 100,
    type = 8,
    names = FALSE
  ))
  retained <- sources[
    abs(sources$amplitude) >= threshold,
    ,
    drop = FALSE
  ]
  if (nrow(retained) == 0L) {
    stop("No source rows remain after thresholding.", call. = FALSE)
  }
  transformed <- .glassbrain_transform(sources, views)
  normalized_retained <- transformed$sources[
    transformed$sources$source_id %in% retained$source_id,
    ,
    drop = FALSE
  ]
  normalized_retained <- normalized_retained[
    match(retained$source_id, normalized_retained$source_id),
    ,
    drop = FALSE
  ]
  zero_activation <- all(sources$amplitude == 0)

  settings <- list(
    source = resolved$source,
    reduction = resolved$reduction,
    method = resolved$method,
    coordinates = resolved$coordinates,
    warnings = resolved$warnings,
    views = views,
    threshold_pct = threshold_pct,
    threshold = threshold,
    source_count = nrow(sources),
    retained_count = nrow(retained),
    grid_resolution = 80L,
    transform_center = transformed$center,
    transform_scale = transformed$scale,
    transform_ranges = transformed$ranges,
    axis_mapping = .glassbrain_view_axes,
    outline = "schematic",
    zero_activation = zero_activation
  )

  if (nrow(retained) < 3L) {
    warning(
      sprintf(
        paste0(
          "Only %d source rows remain after thresholding; using the ",
          "eegPlotSource() sparse fallback."
        ),
        nrow(retained)
      ),
      call. = FALSE
    )
    first_view <- views[1L]
    axes <- .glassbrain_view_axes[[first_view]]
    fallback_data <- data.frame(
      x = normalized_retained[[axes[["projected_x"]]]],
      y = normalized_retained[[axes[["projected_y"]]]],
      amplitude = normalized_retained$amplitude,
      source_id = normalized_retained$source_id,
      stringsAsFactors = FALSE
    )
    plot <- eegPlotSource(
      x,
      source_data = fallback_data,
      method = "scatter",
      threshold_pct = 0
    ) +
      ggplot2::labs(
        title = paste("Sparse", tools::toTitleCase(first_view), "Projection"),
        subtitle = paste0(
          nrow(retained),
          " source row(s); schematic sparse fallback, no second threshold"
        )
      )
    settings$fallback_view <- first_view
    settings$fallback_data <- fallback_data
    attr(plot, "glassbrain_data") <- list(
      panels = list(),
      retained_sources = retained,
      settings = settings
    )
    return(plot)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for multi-view glass-brain output.",
      call. = FALSE
    )
  }
  outline <- .glassbrain_outline()
  panels <- vector("list", length(views))
  panel_data <- vector("list", length(views))
  names(panel_data) <- views
  maximum_amplitude <- max(abs(retained$amplitude))
  if (maximum_amplitude == 0) {
    maximum_amplitude <- 1
  }
  maximum_intensity <- max(abs(retained$amplitude))
  if (maximum_intensity == 0) {
    maximum_intensity <- 1
  }
  colors <- PhysioCore::physioPalette(256L, "diverging")

  for (view_index in seq_along(views)) {
    view <- views[view_index]
    projected <- .glassbrain_project(
      normalized_retained,
      view,
      threshold,
      nrow(sources),
      grid_resolution = 80L
    )
    panel_data[[view]] <- projected
    view_outline <- outline[outline$view == view, , drop = FALSE]
    labels <- .glassbrain_axis_labels(view)
    subtitle <- if (zero_activation) {
      "Zero activation; schematic outline"
    } else {
      "Maximum-intensity projection; schematic outline"
    }
    panels[[view_index]] <- ggplot2::ggplot(
      projected,
      ggplot2::aes(
        x = .data$projected_x,
        y = .data$projected_y
      )
    ) +
      ggplot2::geom_path(
        data = view_outline,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = interaction(.data$view, .data$path_id)
        ),
        inherit.aes = FALSE,
        color = "#5B6470",
        linewidth = 0.55
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .data$amplitude,
          size = .data$intensity,
          alpha = .data$intensity
        )
      ) +
      ggplot2::scale_color_gradientn(
        colors = colors,
        limits = c(-maximum_amplitude, maximum_amplitude),
        name = "Amplitude"
      ) +
      ggplot2::scale_size_continuous(
        limits = c(0, maximum_intensity),
        range = c(1.5, 6),
        name = "Intensity"
      ) +
      ggplot2::scale_alpha_continuous(
        limits = c(0, maximum_intensity),
        range = c(0.35, 0.9),
        guide = "none"
      ) +
      ggplot2::coord_fixed(
        xlim = c(-1.1, 1.1),
        ylim = c(-1.1, 1.1),
        expand = FALSE
      ) +
      ggplot2::labs(
        x = labels[["x"]],
        y = labels[["y"]],
        title = tools::toTitleCase(view),
        subtitle = subtitle
      ) +
      PhysioCore::theme_physio() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      )
  }
  plot <- patchwork::wrap_plots(panels, nrow = 1L) +
    patchwork::plot_layout(guides = "collect")
  attr(plot, "glassbrain_data") <- list(
    panels = panel_data,
    retained_sources = retained,
    settings = settings
  )
  plot
}
