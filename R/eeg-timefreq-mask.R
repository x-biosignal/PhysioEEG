.tf_axis_names <- function(x) {
  vapply(x, function(value) sprintf("%.17g", value), character(1))
}

.tf_axis_from_dimnames <- function(tf, margin, label) {
  axis_names <- dimnames(tf)[[margin]]
  if (is.null(axis_names)) {
    return(list(values = seq_len(dim(tf)[margin]), declared = FALSE))
  }
  values <- suppressWarnings(as.numeric(axis_names))
  if (length(values) != dim(tf)[margin] ||
      anyNA(values) || any(!is.finite(values)) ||
      anyDuplicated(values) || is.unsorted(values, strictly = TRUE)) {
    stop(
      sprintf(
        "%s dimnames must be finite, unique numeric values in increasing order.",
        label
      ),
      call. = FALSE
    )
  }
  list(values = values, declared = TRUE)
}

.tf_baseline_indices <- function(baseline, time, time_declared) {
  n_time <- length(time)
  if (is.logical(baseline)) {
    if (length(baseline) != n_time || anyNA(baseline)) {
      stop(
        "'baseline' must be a non-missing logical vector with one value per time point.",
        call. = FALSE
      )
    }
    indices <- which(baseline)
  } else if (is.numeric(baseline) &&
             length(baseline) == 2L &&
             all(is.finite(baseline)) &&
             baseline[1L] < baseline[2L]) {
    if (time_declared) {
      indices <- which(time >= baseline[1L] & time <= baseline[2L])
    } else {
      if (any(baseline != floor(baseline)) ||
          baseline[1L] < 1L || baseline[2L] > n_time) {
        stop(
          "Without numeric time dimnames, 'baseline' must be two increasing sample indices.",
          call. = FALSE
        )
      }
      indices <- seq.int(as.integer(baseline[1L]), as.integer(baseline[2L]))
    }
  } else {
    stop(
      paste0(
        "'baseline' must be a logical time vector or two increasing finite ",
        "time values/sample indices."
      ),
      call. = FALSE
    )
  }
  if (length(indices) < 2L) {
    stop("'baseline' must select at least two time points.", call. = FALSE)
  }
  indices
}

.tf_t_field <- function(contrasts) {
  dimensions <- dim(contrasts)
  n_replicates <- dimensions[3L]
  values <- matrix(
    contrasts,
    nrow = dimensions[1L] * dimensions[2L],
    ncol = n_replicates
  )
  means <- rowMeans(values)
  centered <- sweep(values, 1L, means, FUN = "-")
  standard_deviation <- sqrt(rowSums(centered^2) / (n_replicates - 1L))
  estimable <- is.finite(standard_deviation) & standard_deviation > 0
  statistic <- rep(NA_real_, length(means))
  statistic[estimable] <- means[estimable] /
    (standard_deviation[estimable] / sqrt(n_replicates))
  list(
    statistic = matrix(
      statistic,
      nrow = dimensions[1L],
      ncol = dimensions[2L]
    ),
    estimable = matrix(
      estimable,
      nrow = dimensions[1L],
      ncol = dimensions[2L]
    )
  )
}

.tf_pointwise_p <- function(statistic, alternative, degrees_freedom) {
  p_value <- switch(
    alternative,
    two.sided = 2 * stats::pt(
      -abs(statistic),
      df = degrees_freedom
    ),
    greater = stats::pt(
      statistic,
      df = degrees_freedom,
      lower.tail = FALSE
    ),
    less = stats::pt(statistic, df = degrees_freedom)
  )
  p_value[!is.finite(statistic)] <- NA_real_
  p_value
}

.tf_label_components <- function(active) {
  n_time <- nrow(active)
  n_frequency <- ncol(active)
  component_id <- matrix(0L, nrow = n_time, ncol = n_frequency)
  queue_time <- integer(length(active))
  queue_frequency <- integer(length(active))
  next_id <- 0L

  for (frequency_index in seq_len(n_frequency)) {
    for (time_index in seq_len(n_time)) {
      if (!active[time_index, frequency_index] ||
          component_id[time_index, frequency_index] != 0L) {
        next
      }
      next_id <- next_id + 1L
      queue_head <- 1L
      queue_tail <- 1L
      queue_time[1L] <- time_index
      queue_frequency[1L] <- frequency_index
      component_id[time_index, frequency_index] <- next_id

      while (queue_head <= queue_tail) {
        current_time <- queue_time[queue_head]
        current_frequency <- queue_frequency[queue_head]
        queue_head <- queue_head + 1L
        neighbours <- rbind(
          c(current_time - 1L, current_frequency),
          c(current_time + 1L, current_frequency),
          c(current_time, current_frequency - 1L),
          c(current_time, current_frequency + 1L)
        )
        for (neighbour_index in seq_len(nrow(neighbours))) {
          neighbour_time <- neighbours[neighbour_index, 1L]
          neighbour_frequency <- neighbours[neighbour_index, 2L]
          if (neighbour_time < 1L || neighbour_time > n_time ||
              neighbour_frequency < 1L ||
              neighbour_frequency > n_frequency ||
              !active[neighbour_time, neighbour_frequency] ||
              component_id[neighbour_time, neighbour_frequency] != 0L) {
            next
          }
          queue_tail <- queue_tail + 1L
          queue_time[queue_tail] <- neighbour_time
          queue_frequency[queue_tail] <- neighbour_frequency
          component_id[neighbour_time, neighbour_frequency] <- next_id
        }
      }
    }
  }
  component_id
}

.tf_empty_cluster_table <- function() {
  data.frame(
    cluster_id = integer(),
    sign = character(),
    mass = double(),
    size = integer(),
    time_index_min = integer(),
    time_index_max = integer(),
    frequency_index_min = integer(),
    frequency_index_max = integer(),
    time_min = double(),
    time_max = double(),
    frequency_min = double(),
    frequency_max = double(),
    p_value = double(),
    significant = logical(),
    stringsAsFactors = FALSE
  )
}

.tf_cluster_field <- function(statistic, critical_value, alternative,
                              time, frequency) {
  positive <- is.finite(statistic) & statistic >= critical_value
  negative <- is.finite(statistic) & statistic <= -critical_value
  if (alternative == "greater") {
    negative[] <- FALSE
  } else if (alternative == "less") {
    positive[] <- FALSE
  }

  combined_id <- matrix(0L, nrow = nrow(statistic), ncol = ncol(statistic))
  table <- .tf_empty_cluster_table()
  next_id <- 0L
  for (cluster_sign in c("positive", "negative")) {
    active <- if (cluster_sign == "positive") positive else negative
    local_id <- .tf_label_components(active)
    n_local <- max(local_id)
    if (n_local == 0L) {
      next
    }
    for (local_cluster in seq_len(n_local)) {
      locations <- which(local_id == local_cluster, arr.ind = TRUE)
      next_id <- next_id + 1L
      combined_id[local_id == local_cluster] <- next_id
      values <- statistic[local_id == local_cluster]
      mass <- if (cluster_sign == "positive") sum(values) else sum(-values)
      table <- rbind(
        table,
        data.frame(
          cluster_id = next_id,
          sign = cluster_sign,
          mass = mass,
          size = nrow(locations),
          time_index_min = min(locations[, 1L]),
          time_index_max = max(locations[, 1L]),
          frequency_index_min = min(locations[, 2L]),
          frequency_index_max = max(locations[, 2L]),
          time_min = time[min(locations[, 1L])],
          time_max = time[max(locations[, 1L])],
          frequency_min = frequency[min(locations[, 2L])],
          frequency_max = frequency[max(locations[, 2L])],
          p_value = NA_real_,
          significant = FALSE,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  list(cluster_id = combined_id, cluster_table = table)
}

.tf_sign_matrix <- function(n_replicates, n_permutations, alpha) {
  # The all-positive identity is represented by the conservative +1 term.
  exact_count <- 2^n_replicates - 1
  use_exact <- is.finite(exact_count) &&
    exact_count <= 4096L &&
    exact_count <= n_permutations
  if (use_exact) {
    pattern <- seq.int(0, exact_count - 1)
    signs <- vapply(
      seq_len(n_replicates),
      function(bit) {
        ifelse(
          bitwAnd(pattern, bitwShiftL(1L, bit - 1L)) == 0L,
          -1,
          1
        )
      },
      numeric(length(pattern))
    )
    return(list(signs = signs, mode = "exact"))
  }
  minimum <- ceiling(1 / alpha)
  if (n_permutations < minimum) {
    stop(
      sprintf(
        "'n_permutations' must be at least %d so alpha = %.6g is attainable.",
        minimum,
        alpha
      ),
      call. = FALSE
    )
  }
  list(
    signs = matrix(
      sample(c(-1, 1), n_permutations * n_replicates, replace = TRUE),
      nrow = n_permutations,
      ncol = n_replicates
    ),
    mode = "random"
  )
}

# Private, auditable one-sample time-frequency inference.
.tfClusterMask <- function(
    tf,
    baseline,
    method = c("threshold", "cluster"),
    cluster_alpha = 0.05,
    alpha = 0.05,
    n_permutations = 999L,
    alternative = c("two.sided", "greater", "less"),
    connectivity = 4L,
    p_adjust_method = "BH",
    seed = NULL) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  method <- match.arg(method)
  alternative <- match.arg(alternative)
  if (!is.array(tf) || !is.numeric(tf) || length(dim(tf)) != 3L ||
      any(dim(tf) < c(2L, 2L, 3L)) || any(!is.finite(tf))) {
    stop(
      paste0(
        "'tf' must be a finite numeric time x frequency x replicate array ",
        "with dimensions of at least 2 x 2 x 3."
      ),
      call. = FALSE
    )
  }
  validate_probability <- function(value, label, upper_inclusive = FALSE) {
    valid_upper <- if (upper_inclusive) value <= 1 else value < 1
    if (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value) || value <= 0 || !valid_upper) {
      stop(
        sprintf(
          "'%s' must be one finite number in %s.",
          label,
          if (upper_inclusive) "(0, 1]" else "(0, 1)"
        ),
        call. = FALSE
      )
    }
  }
  validate_probability(alpha, "alpha", upper_inclusive = TRUE)
  validate_probability(cluster_alpha, "cluster_alpha")
  if (cluster_alpha > 0.5) {
    stop("'cluster_alpha' must not exceed 0.5.", call. = FALSE)
  }
  if (!is.numeric(connectivity) || length(connectivity) != 1L ||
      is.na(connectivity) || connectivity != 4) {
    stop("Only 4-neighbour 'connectivity = 4' is supported.", call. = FALSE)
  }
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods) {
    stop("'p_adjust_method' must name a stats::p.adjust method.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
       seed != floor(seed))) {
    stop(
      "'seed' must be NULL or one non-negative integer representable by R.",
      call. = FALSE
    )
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  time_info <- .tf_axis_from_dimnames(tf, 1L, "Time")
  frequency_info <- .tf_axis_from_dimnames(tf, 2L, "Frequency")
  time <- time_info$values
  frequency <- frequency_info$values
  baseline_indices <- .tf_baseline_indices(
    baseline,
    time,
    time_info$declared
  )

  dimensions <- dim(tf)
  n_replicates <- dimensions[3L]
  baseline_mean <- apply(
    tf[baseline_indices, , , drop = FALSE],
    c(2L, 3L),
    mean
  )
  contrasts <- sweep(tf, c(2L, 3L), baseline_mean, FUN = "-")
  field <- .tf_t_field(contrasts)
  statistic <- field$statistic
  estimable <- field$estimable
  degrees_freedom <- n_replicates - 1L
  pointwise_p <- .tf_pointwise_p(
    statistic,
    alternative,
    degrees_freedom
  )
  warning_messages <- character()
  n_non_estimable <- sum(!estimable)
  if (n_non_estimable > 0L) {
    warning_messages <- sprintf(
      "%d zero-variance cells were marked non-estimable.",
      n_non_estimable
    )
  }

  axis_names <- list(.tf_axis_names(time), .tf_axis_names(frequency))
  dimnames(statistic) <- axis_names
  dimnames(estimable) <- axis_names
  dimnames(pointwise_p) <- axis_names
  cluster_id <- matrix(
    0L,
    nrow = dimensions[1L],
    ncol = dimensions[2L],
    dimnames = axis_names
  )
  cluster_table <- .tf_empty_cluster_table()
  null_max <- numeric()
  permutation_mode <- "not_applicable"

  if (method == "threshold") {
    family_p <- pointwise_p
    family_p[is.na(family_p)] <- 1
    p_adjusted <- matrix(
      stats::p.adjust(as.vector(family_p), method = p_adjust_method),
      nrow = dimensions[1L],
      ncol = dimensions[2L],
      dimnames = axis_names
    )
    mask <- estimable & p_adjusted <= alpha
    realized_permutations <- 0L
  } else {
    if (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
        is.na(n_permutations) || !is.finite(n_permutations) ||
        n_permutations < 1L || n_permutations != floor(n_permutations) ||
        n_permutations > .Machine$integer.max) {
      stop("'n_permutations' must be one positive integer.", call. = FALSE)
    }
    critical_probability <- if (alternative == "two.sided") {
      1 - cluster_alpha / 2
    } else {
      1 - cluster_alpha
    }
    critical_value <- stats::qt(
      critical_probability,
      df = degrees_freedom
    )
    signs_result <- .tf_sign_matrix(
      n_replicates,
      as.integer(n_permutations),
      alpha
    )
    signs <- signs_result$signs
    permutation_mode <- signs_result$mode
    realized_permutations <- nrow(signs)
    observed <- .tf_cluster_field(
      statistic,
      critical_value,
      alternative,
      time,
      frequency
    )
    cluster_id <- observed$cluster_id
    cluster_table <- observed$cluster_table
    null_max <- numeric(realized_permutations)
    for (permutation_index in seq_len(realized_permutations)) {
      permuted <- sweep(
        contrasts,
        3L,
        signs[permutation_index, ],
        FUN = "*"
      )
      null_field <- .tf_t_field(permuted)$statistic
      null_clusters <- .tf_cluster_field(
        null_field,
        critical_value,
        alternative,
        time,
        frequency
      )$cluster_table
      null_max[permutation_index] <- if (nrow(null_clusters) == 0L) {
        0
      } else {
        max(null_clusters$mass)
      }
    }
    mask <- matrix(
      FALSE,
      nrow = dimensions[1L],
      ncol = dimensions[2L],
      dimnames = axis_names
    )
    p_adjusted <- matrix(
      1,
      nrow = dimensions[1L],
      ncol = dimensions[2L],
      dimnames = axis_names
    )
    if (nrow(cluster_table) > 0L) {
      cluster_table$p_value <- vapply(
        cluster_table$mass,
        function(mass) {
          (1 + sum(null_max >= mass)) / (realized_permutations + 1)
        },
        numeric(1)
      )
      cluster_table$significant <- cluster_table$p_value <= alpha
      for (cluster_index in seq_len(nrow(cluster_table))) {
        cells <- cluster_id == cluster_table$cluster_id[cluster_index]
        p_adjusted[cells] <- cluster_table$p_value[cluster_index]
        if (cluster_table$significant[cluster_index]) {
          mask[cells] <- TRUE
        }
      }
    }
    dimnames(cluster_id) <- axis_names
  }
  dimnames(mask) <- axis_names

  result <- list(
    mask = mask,
    statistic = statistic,
    pointwise_p = pointwise_p,
    p_adjusted = p_adjusted,
    estimable = estimable,
    cluster_id = cluster_id,
    cluster_table = cluster_table,
    null_max = null_max,
    time = time,
    frequency = frequency,
    baseline_indices = baseline_indices,
    method = method,
    alpha = alpha,
    cluster_alpha = cluster_alpha,
    alternative = alternative,
    connectivity = 4L,
    p_adjust_method = p_adjust_method,
    n_permutations = realized_permutations,
    permutation_mode = permutation_mode,
    seed = if (is.null(seed)) NULL else as.integer(seed),
    warnings = warning_messages
  )
  class(result) <- c("tf_cluster_mask", "list")
  result
}

.eeg_tf_validate_axis <- function(values, expected_length, label) {
  if (!is.numeric(values) || length(values) != expected_length ||
      anyNA(values) || any(!is.finite(values)) ||
      anyDuplicated(values) || is.unsorted(values, strictly = TRUE)) {
    stop(
      sprintf(
        "%s axis must contain %d finite, unique, strictly increasing values.",
        label,
        expected_length
      ),
      call. = FALSE
    )
  }
  as.numeric(values)
}

.eeg_tf_channel <- function(channel, labels, n_channels) {
  if (length(labels) != n_channels) {
    stop(
      "Channel metadata do not match the time-frequency channel dimension.",
      call. = FALSE
    )
  }
  if (is.character(channel)) {
    if (length(channel) != 1L || is.na(channel) || !nzchar(channel)) {
      stop("'channel' must be one exact, non-empty channel label.", call. = FALSE)
    }
    matches <- which(labels == channel)
    if (length(matches) == 0L) {
      stop(sprintf("Channel '%s' not found.", channel), call. = FALSE)
    }
    if (length(matches) > 1L) {
      stop(
        sprintf("Channel label '%s' is duplicated and therefore ambiguous.", channel),
        call. = FALSE
      )
    }
    return(list(index = matches, label = labels[matches]))
  }
  if (!is.numeric(channel) || length(channel) != 1L ||
      is.na(channel) || !is.finite(channel) ||
      channel != floor(channel) || channel < 1L || channel > n_channels) {
    stop(
      sprintf("'channel' must be one integer from 1 to %d.", n_channels),
      call. = FALSE
    )
  }
  index <- as.integer(channel)
  list(index = index, label = labels[index])
}

.eeg_tf_metadata_spec <- function(md, key, array, sampling_rate) {
  dimensions <- dim(array)
  candidates <- list()
  add_candidate <- function(source, time, frequency, scale) {
    candidates[[length(candidates) + 1L]] <<- list(
      source = source,
      time = time,
      frequency = frequency,
      scale = scale
    )
  }
  if (is.list(md$stft) &&
      length(md$stft$time_axis) == dimensions[1L] &&
      length(md$stft$freq_axis) == dimensions[2L]) {
    add_candidate(
      "metadata:stft",
      md$stft$time_axis,
      md$stft$freq_axis,
      "power"
    )
    candidates[[length(candidates)]]$settings <- md$stft[
      setdiff(names(md$stft), c("time_axis", "freq_axis"))
    ]
  }
  if (is.list(md$wavelet) &&
      length(md$wavelet$frequencies) == dimensions[2L]) {
    add_candidate(
      "metadata:wavelet",
      (seq_len(dimensions[1L]) - 1) / sampling_rate,
      md$wavelet$frequencies,
      "power"
    )
    candidates[[length(candidates)]]$settings <- md$wavelet[
      setdiff(names(md$wavelet), c("frequencies", "phase"))
    ]
  }
  if (is.list(md$ersp) &&
      length(md$ersp$frequencies) == dimensions[2L]) {
    add_candidate(
      "metadata:ersp",
      (seq_len(dimensions[1L]) - 1) / sampling_rate,
      md$ersp$frequencies,
      "db"
    )
    candidates[[length(candidates)]]$settings <- md$ersp[
      setdiff(names(md$ersp), "frequencies")
    ]
  }

  preferred <- switch(
    key,
    stft_power = "metadata:stft",
    wavelet_power = "metadata:wavelet",
    ersp_data = "metadata:ersp",
    ersp = "metadata:ersp",
    NULL
  )
  if (!is.null(preferred)) {
    preferred_matches <- vapply(
      candidates,
      function(candidate) identical(candidate$source, preferred),
      logical(1)
    )
    if (sum(preferred_matches) == 1L) {
      return(candidates[[which(preferred_matches)]])
    }
  }
  if (length(candidates) != 1L) {
    stop(
      sprintf(
        paste0(
          "Metadata product '%s' has no unambiguous STFT, Morlet, or ERSP ",
          "axis contract; specify a uniquely matching producer output."
        ),
        key
      ),
      call. = FALSE
    )
  }
  candidates[[1L]]
}

.eeg_tf_metadata_product <- function(md, key, sampling_rate, labels) {
  product <- md[[key]]
  if (!is.array(product) || !is.numeric(product) ||
      length(dim(product)) != 3L) {
    stop(
      sprintf(
        "Metadata product '%s' must be a numeric time x frequency x channel array.",
        key
      ),
      call. = FALSE
    )
  }
  if (dim(product)[3L] != length(labels)) {
    stop(
      sprintf(
        "Metadata product '%s' channel dimension does not match colData.",
        key
      ),
      call. = FALSE
    )
  }
  spec <- .eeg_tf_metadata_spec(md, key, product, sampling_rate)
  spec$array <- product
  spec$key <- key
  spec
}

.eeg_tf_known_products <- function(md, sampling_rate, labels) {
  known_keys <- c("stft_power", "wavelet_power", "ersp_data")
  known_keys <- known_keys[vapply(
    known_keys,
    function(key) !is.null(md[[key]]),
    logical(1)
  )]
  if (length(known_keys) == 0L) {
    return(list())
  }
  lapply(
    known_keys,
    function(key) .eeg_tf_metadata_product(md, key, sampling_rate, labels)
  )
}

.eeg_tf_sliding_fft <- function(data, channel_index, sampling_rate) {
  if (!is.matrix(data) || !is.numeric(data) || any(!is.finite(data))) {
    stop(
      "The legacy sliding FFT requires a finite numeric 2D time x channel assay.",
      call. = FALSE
    )
  }
  signal <- data[, channel_index]
  n_samples <- length(signal)
  window_length <- min(as.integer(round(0.5 * sampling_rate)), n_samples)
  if (window_length %% 2L != 0L) {
    window_length <- window_length - 1L
  }
  if (window_length < 4L) {
    stop(
      "Signal is too short to form a non-degenerate even spectrogram window.",
      call. = FALSE
    )
  }
  hop <- max(1L, window_length %/% 4L)
  window <- 0.5 * (
    1 - cos(
      2 * pi * seq(0, window_length - 1L) / (window_length - 1L)
    )
  )
  n_fft_bins <- window_length %/% 2L + 1L
  frequency <- seq(
    0,
    sampling_rate / 2,
    length.out = n_fft_bins
  )
  starts <- seq.int(
    1L,
    n_samples - window_length + 1L,
    by = hop
  )
  time <- (starts + window_length / 2 - 1) / sampling_rate
  power <- matrix(
    0,
    nrow = length(starts),
    ncol = n_fft_bins
  )
  for (window_index in seq_along(starts)) {
    indices <- starts[window_index]:(
      starts[window_index] + window_length - 1L
    )
    transformed <- stats::fft(signal[indices] * window)[seq_len(n_fft_bins)]
    power[window_index, ] <- Mod(transformed)^2 / window_length
  }
  list(
    power = power,
    time = time,
    frequency = frequency,
    source = "computed:legacy_sliding_fft",
    scale = "power",
    settings = list(
      window_sec = 0.5,
      overlap = 0.75,
      window_length = window_length,
      hop_size = hop,
      normalization = "Mod(fft)^2/window_length"
    )
  )
}

.eeg_extract_time_frequency <- function(
    x,
    channel,
    assay_name = NULL,
    log_power = TRUE) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("'x' must inherit from PhysioExperiment.", call. = FALSE)
  }
  if (!is.logical(log_power) || length(log_power) != 1L || is.na(log_power)) {
    stop("'log_power' must be TRUE or FALSE.", call. = FALSE)
  }
  sampling_rate <- samplingRate(x)
  if (!is.numeric(sampling_rate) || length(sampling_rate) != 1L ||
      is.na(sampling_rate) || !is.finite(sampling_rate) ||
      sampling_rate <= 0) {
    stop("'x' must have one finite positive sampling rate.", call. = FALSE)
  }
  col_data <- SummarizedExperiment::colData(x)
  labels <- if ("label" %in% colnames(col_data)) {
    as.character(col_data$label)
  } else {
    paste0("Ch", seq_len(nrow(col_data)))
  }
  if (anyNA(labels)) {
    stop("Channel labels must not be missing.", call. = FALSE)
  }
  md <- S4Vectors::metadata(x)
  assay_names <- SummarizedExperiment::assayNames(x)
  explicit <- !is.null(assay_name)
  if (explicit &&
      (!is.character(assay_name) || length(assay_name) != 1L ||
       is.na(assay_name) || !nzchar(assay_name))) {
    stop("'assay_name' must be NULL or one non-empty name.", call. = FALSE)
  }

  selected <- NULL
  if (explicit && assay_name %in% assay_names) {
    assay_data <- SummarizedExperiment::assay(x, assay_name)
    if (length(dim(assay_data)) == 3L && !is.null(md$freqs)) {
      selected <- list(
        array = assay_data,
        time = if (!is.null(md$times)) {
          md$times
        } else {
          (seq_len(dim(assay_data)[1L]) - 1) / sampling_rate
        },
        frequency = md$freqs,
        source = paste0("assay:", assay_name),
        scale = "power",
        settings = list(axis_metadata = "metadata$times/metadata$freqs")
      )
    } else if (length(dim(assay_data)) == 2L) {
      channel_info <- .eeg_tf_channel(channel, labels, ncol(assay_data))
      selected <- .eeg_tf_sliding_fft(
        assay_data,
        channel_info$index,
        sampling_rate
      )
      selected$channel <- channel_info
      selected$source <- paste0(selected$source, ":", assay_name)
    } else {
      stop(
        paste0(
          "A 3D assay is not assumed to be time-frequency data. Supply ",
          "matching time/frequency metadata or use eegSTFT(), ",
          "eegMorletWavelet(), or eegERSP()."
        ),
        call. = FALSE
      )
    }
  } else if (explicit && !is.null(md[[assay_name]])) {
    selected <- .eeg_tf_metadata_product(
      md,
      assay_name,
      sampling_rate,
      labels
    )
  } else if (explicit) {
    stop(
      sprintf(
        "'assay_name' '%s' is neither an assay nor a known metadata product.",
        assay_name
      ),
      call. = FALSE
    )
  } else {
    known <- .eeg_tf_known_products(md, sampling_rate, labels)
    if (length(known) > 1L) {
      stop(
        paste0(
          "Multiple time-frequency metadata products are available; ",
          "select one exactly with 'assay_name'."
        ),
        call. = FALSE
      )
    }
    if (length(known) == 1L) {
      selected <- known[[1L]]
    } else {
      default_name <- defaultAssay(x)
      assay_data <- SummarizedExperiment::assay(x, default_name)
      if (length(dim(assay_data)) != 2L) {
        stop(
          paste0(
            "A raw 3D assay is time x channel x epoch, not a spectrogram. ",
            "Use eegSTFT(), eegMorletWavelet(), or eegERSP(), then select ",
            "the stored product."
          ),
          call. = FALSE
        )
      }
      channel_info <- .eeg_tf_channel(channel, labels, ncol(assay_data))
      selected <- .eeg_tf_sliding_fft(
        assay_data,
        channel_info$index,
        sampling_rate
      )
      selected$channel <- channel_info
      selected$source <- paste0(selected$source, ":", default_name)
    }
  }

  if (!is.null(selected$array)) {
    channel_info <- .eeg_tf_channel(
      channel,
      labels,
      dim(selected$array)[3L]
    )
    selected$power <- selected$array[, , channel_info$index, drop = TRUE]
    if (!is.matrix(selected$power)) {
      selected$power <- matrix(
        selected$power,
        nrow = dim(selected$array)[1L],
        ncol = dim(selected$array)[2L]
      )
    }
    selected$channel <- channel_info
    selected$array <- NULL
    selected$settings <- c(
      selected$settings,
      list(metadata_product = selected$key)
    )
  }
  if (!is.matrix(selected$power) || !is.numeric(selected$power) ||
      any(!is.finite(selected$power))) {
    stop("Time-frequency power must be a finite numeric matrix.", call. = FALSE)
  }
  selected$time <- .eeg_tf_validate_axis(
    selected$time,
    nrow(selected$power),
    "Time"
  )
  selected$frequency <- .eeg_tf_validate_axis(
    selected$frequency,
    ncol(selected$power),
    "Frequency"
  )

  log_applied <- FALSE
  if (identical(selected$scale, "db")) {
    power_label <- "Power (dB)"
  } else if (log_power) {
    selected$power <- 10 * log10(
      pmax(selected$power, .Machine$double.eps)
    )
    log_applied <- TRUE
    power_label <- "Power (dB)"
  } else {
    power_label <- "Power"
  }
  dimnames(selected$power) <- list(
    .tf_axis_names(selected$time),
    .tf_axis_names(selected$frequency)
  )
  selected$settings <- c(
    selected$settings,
    list(
      input_scale = selected$scale,
      log_requested = log_power,
      log_applied = log_applied,
      power_label = power_label
    )
  )
  selected[c(
    "power",
    "time",
    "frequency",
    "source",
    "channel",
    "settings"
  )]
}

.eeg_tf_prepare_mask <- function(mask, time, frequency) {
  if (is.null(mask)) {
    return(NULL)
  }
  expected_dimensions <- c(length(time), length(frequency))
  expected_dimnames <- list(
    .tf_axis_names(time),
    .tf_axis_names(frequency)
  )
  if (inherits(mask, "tf_cluster_mask") || (
      is.list(mask) &&
      all(c("mask", "time", "frequency", "method", "alpha") %in% names(mask))
  )) {
    if (!is.numeric(mask$time) || !identical(as.numeric(mask$time), time) ||
        !is.numeric(mask$frequency) ||
        !identical(as.numeric(mask$frequency), frequency)) {
      stop(
        "Inference-result time/frequency axes do not match the plotted axes.",
        call. = FALSE
      )
    }
    mask <- mask$mask
  }
  if (!is.matrix(mask) || !identical(dim(mask), expected_dimensions)) {
    stop(
      sprintf(
        "'mask' must be a %d x %d time-by-frequency matrix.",
        expected_dimensions[1L],
        expected_dimensions[2L]
      ),
      call. = FALSE
    )
  }
  if (!identical(dimnames(mask), expected_dimnames)) {
    stop(
      paste0(
        "'mask' dimnames must exactly match the complete plotted time and ",
        "frequency axes."
      ),
      call. = FALSE
    )
  }
  if (is.logical(mask)) {
    if (anyNA(mask)) {
      stop("Logical 'mask' values must not be missing.", call. = FALSE)
    }
    return(mask)
  }
  if (!is.numeric(mask) || anyNA(mask) || any(!is.finite(mask)) ||
      any(mask < 0 | mask > 1)) {
    stop(
      "Numeric 'mask' values must be finite p-values in [0, 1].",
      call. = FALSE
    )
  }
  alpha <- attr(mask, "alpha", exact = TRUE)
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop(
      "A numeric p-value 'mask' requires one finite 'alpha' attribute in [0, 1].",
      call. = FALSE
    )
  }
  result <- mask <= alpha
  dimnames(result) <- expected_dimnames
  result
}

.eeg_tf_axis_edges <- function(axis) {
  midpoints <- (axis[-length(axis)] + axis[-1L]) / 2
  c(
    axis[1L] - (axis[2L] - axis[1L]) / 2,
    midpoints,
    axis[length(axis)] +
      (axis[length(axis)] - axis[length(axis) - 1L]) / 2
  )
}

.eeg_tf_mask_boundaries <- function(mask, time, frequency) {
  if (nrow(mask) < 2L || ncol(mask) < 2L ||
      all(mask) || !any(mask)) {
    return(data.frame(
      x = double(),
      y = double(),
      xend = double(),
      yend = double()
    ))
  }
  time_edges <- .eeg_tf_axis_edges(time)
  frequency_edges <- .eeg_tf_axis_edges(frequency)
  segments <- vector("list", 0L)

  time_changes <- mask[-nrow(mask), , drop = FALSE] !=
    mask[-1L, , drop = FALSE]
  locations <- which(time_changes, arr.ind = TRUE)
  if (nrow(locations) > 0L) {
    for (index in seq_len(nrow(locations))) {
      time_index <- locations[index, 1L]
      frequency_index <- locations[index, 2L]
      boundary_time <- (time[time_index] + time[time_index + 1L]) / 2
      segments[[length(segments) + 1L]] <- data.frame(
        x = boundary_time,
        y = frequency_edges[frequency_index],
        xend = boundary_time,
        yend = frequency_edges[frequency_index + 1L]
      )
    }
  }
  frequency_changes <- mask[, -ncol(mask), drop = FALSE] !=
    mask[, -1L, drop = FALSE]
  locations <- which(frequency_changes, arr.ind = TRUE)
  if (nrow(locations) > 0L) {
    for (index in seq_len(nrow(locations))) {
      time_index <- locations[index, 1L]
      frequency_index <- locations[index, 2L]
      boundary_frequency <- (
        frequency[frequency_index] + frequency[frequency_index + 1L]
      ) / 2
      segments[[length(segments) + 1L]] <- data.frame(
        x = time_edges[time_index],
        y = boundary_frequency,
        xend = time_edges[time_index + 1L],
        yend = boundary_frequency
      )
    }
  }
  if (length(segments) == 0L) {
    return(data.frame(
      x = double(),
      y = double(),
      xend = double(),
      yend = double()
    ))
  }
  do.call(rbind, segments)
}
