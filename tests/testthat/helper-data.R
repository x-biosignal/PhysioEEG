#' Create basic test EEG PhysioExperiment
#'
#' Generates multi-channel EEG with alpha (10 Hz) and beta (20 Hz) oscillations
#' plus pink noise. Channels are labeled using 10-20 system names.
#'
#' @param n_time Number of time points (default: 5000 = 10s at 500 Hz)
#' @param n_channels Number of EEG channels (default: 19, standard 10-20)
#' @param sr Sampling rate in Hz (default: 500)
#' @return PhysioExperiment with simulated EEG data
make_eeg <- function(n_time = 5000, n_channels = 19, sr = 500) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  ch_labels_pool <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                       "T3", "C3", "Cz", "C4", "T4",
                       "T5", "P3", "Pz", "P4", "T6",
                       "O1", "O2")
  if (n_channels <= length(ch_labels_pool)) {
    ch_labels <- ch_labels_pool[seq_len(n_channels)]
  } else {
    ch_labels <- paste0("EEG", seq_len(n_channels))
  }

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    alpha_amp <- 10 + 5 * sin(2 * pi * ch / n_channels)
    alpha <- alpha_amp * sin(2 * pi * 10 * t + runif(1, 0, 2 * pi))
    beta_amp <- 3 + 2 * cos(2 * pi * ch / n_channels)
    beta <- beta_amp * sin(2 * pi * 20 * t + runif(1, 0, 2 * pi))
    noise <- cumsum(rnorm(n_time, sd = 2))
    noise <- noise - mean(noise)
    data[, ch] <- alpha + beta + noise + rnorm(n_time, sd = 1)
  }

  data <- data * 0.5

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = ch_labels,
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )
}

#' Create EEG with embedded ERP components
#'
#' Generates epoched EEG (3D array: time x channels x epochs) with known
#' ERP components: N100 (negative peak at ~100ms) and P300 (positive peak
#' at ~300ms). Target epochs have larger P300.
#'
#' @param n_epochs Number of epochs (default: 40)
#' @param n_channels Number of channels (default: 19)
#' @param sr Sampling rate in Hz (default: 250)
#' @param epoch_sec Epoch duration in seconds (default: 1.0)
#' @return PhysioExperiment with 3D assay and event metadata
make_eeg_erp <- function(n_epochs = 40, n_channels = 19, sr = 250,
                         epoch_sec = 1.0) {
  n_time <- as.integer(epoch_sec * sr)
  t <- seq(0, epoch_sec - 1 / sr, length.out = n_time)

  ch_labels_pool <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                       "T3", "C3", "Cz", "C4", "T4",
                       "T5", "P3", "Pz", "P4", "T6",
                       "O1", "O2")
  ch_labels <- if (n_channels <= length(ch_labels_pool)) {
    ch_labels_pool[seq_len(n_channels)]
  } else {
    paste0("EEG", seq_len(n_channels))
  }

  data <- array(NA_real_, dim = c(n_time, n_channels, n_epochs))
  conditions <- rep(c("target", "standard"), length.out = n_epochs)

  for (ep in seq_len(n_epochs)) {
    for (ch in seq_len(n_channels)) {
      sig <- rnorm(n_time, sd = 5)
      n100_center <- 0.100
      n100_amp <- -8
      sig <- sig + n100_amp * exp(-((t - n100_center)^2) / (2 * 0.015^2))
      p300_center <- 0.300
      p300_amp <- if (conditions[ep] == "target") 12 else 3
      sig <- sig + p300_amp * exp(-((t - p300_center)^2) / (2 * 0.030^2))
      data[, ch, ep] <- sig
    }
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = ch_labels,
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )
  metadata(pe)$conditions <- conditions
  pe
}

#' Create EEG with sleep stage characteristics
#'
#' Generates continuous EEG with segments having distinct frequency content
#' matching AASM sleep stages: Wake (alpha), N1 (theta), N2 (spindles + K),
#' N3 (delta/SWA), REM (mixed low-voltage).
#'
#' @param n_time Number of time points (default: 150000 = 5min at 500Hz)
#' @param n_channels Number of channels (default: 2, C3 and C4)
#' @param sr Sampling rate in Hz (default: 500)
#' @return PhysioExperiment with metadata$sleep_stages (ground truth)
make_eeg_sleep <- function(n_time = 150000, n_channels = 2, sr = 500) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  epoch_len <- 30 * sr
  n_full_epochs <- n_time %/% epoch_len
  stage_cycle <- c("W", "N1", "N2", "N3", "N2", "REM")
  stages <- rep(stage_cycle, length.out = n_full_epochs)

  data <- matrix(0, nrow = n_time, ncol = n_channels)

  for (ep in seq_len(n_full_epochs)) {
    idx_start <- (ep - 1) * epoch_len + 1
    idx_end <- min(ep * epoch_len, n_time)
    idx <- idx_start:idx_end
    t_seg <- t[idx]
    stage <- stages[ep]

    for (ch in seq_len(n_channels)) {
      seg <- rnorm(length(idx), sd = 2)

      if (stage == "W") {
        seg <- seg + 15 * sin(2 * pi * 10 * t_seg + runif(1, 0, 2 * pi))
      } else if (stage == "N1") {
        seg <- seg + 8 * sin(2 * pi * 5 * t_seg + runif(1, 0, 2 * pi))
      } else if (stage == "N2") {
        spindle_centers <- seq(from = t_seg[1] + 5, to = max(t_seg) - 5, by = 10)
        for (sc in spindle_centers) {
          spindle_env <- exp(-((t_seg - sc)^2) / (2 * 0.3^2))
          seg <- seg + 10 * spindle_env * sin(2 * pi * 13 * t_seg)
        }
        kc_centers <- seq(from = t_seg[1] + 8, to = max(t_seg) - 8, by = 15)
        for (kc in kc_centers) {
          seg <- seg - 40 * exp(-((t_seg - kc)^2) / (2 * 0.05^2))
          seg <- seg + 20 * exp(-((t_seg - kc - 0.15)^2) / (2 * 0.08^2))
        }
      } else if (stage == "N3") {
        seg <- seg + 50 * sin(2 * pi * 1 * t_seg + runif(1, 0, 2 * pi))
        seg <- seg + 30 * sin(2 * pi * 0.5 * t_seg + runif(1, 0, 2 * pi))
      } else if (stage == "REM") {
        seg <- seg + 5 * sin(2 * pi * 6 * t_seg + runif(1, 0, 2 * pi))
        seg <- seg + 3 * sin(2 * pi * 18 * t_seg + runif(1, 0, 2 * pi))
      }

      data[idx, ch] <- seg
    }
  }

  ch_labels <- c("C3", "C4", "F3", "F4", "O1", "O2")[seq_len(n_channels)]

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = ch_labels,
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )
  metadata(pe)$sleep_stages <- data.frame(
    epoch = seq_len(n_full_epochs),
    stage = stages,
    start_sample = (seq_len(n_full_epochs) - 1) * epoch_len + 1,
    end_sample = pmin(seq_len(n_full_epochs) * epoch_len, n_time)
  )
  pe
}

#' Create EEG with motor imagery / SSVEP patterns for BCI
#'
#' Generates epoched EEG with two classes: left motor imagery (mu ERD over C4)
#' and right motor imagery (mu ERD over C3). Also includes SSVEP at specified
#' frequencies at occipital channels.
#'
#' @param n_trials Number of trials per class (default: 30)
#' @param n_channels Number of channels (default: 8)
#' @param sr Sampling rate in Hz (default: 256)
#' @param trial_sec Trial duration in seconds (default: 4.0)
#' @return PhysioExperiment with 3D assay, metadata$labels
make_eeg_bci <- function(n_trials = 30, n_channels = 8, sr = 256,
                         trial_sec = 4.0) {
  n_time <- as.integer(trial_sec * sr)
  t <- seq(0, trial_sec - 1 / sr, length.out = n_time)
  total_trials <- n_trials * 2

  ch_labels <- c("C3", "Cz", "C4", "P3", "Pz", "P4", "O1", "O2")[seq_len(n_channels)]

  data <- array(NA_real_, dim = c(n_time, n_channels, total_trials))
  labels <- rep(c("left", "right"), each = n_trials)

  for (trial in seq_len(total_trials)) {
    for (ch in seq_len(n_channels)) {
      sig <- rnorm(n_time, sd = 5)
      mu <- 8 * sin(2 * pi * 10 * t + runif(1, 0, 2 * pi))

      if (labels[trial] == "left") {
        if (ch_labels[ch] == "C4") {
          mu <- mu * 0.3
        } else if (ch_labels[ch] == "C3") {
          mu <- mu * 1.5
        }
      } else {
        if (ch_labels[ch] == "C3") {
          mu <- mu * 0.3
        } else if (ch_labels[ch] == "C4") {
          mu <- mu * 1.5
        }
      }

      if (ch_labels[ch] %in% c("O1", "O2")) {
        sig <- sig + 6 * sin(2 * pi * 12 * t)
      }

      data[, ch, trial] <- sig + mu
    }
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = ch_labels,
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )
  metadata(pe)$labels <- labels
  pe
}

#' Create EEG with embedded epileptic spikes
#'
#' Generates multi-channel EEG with sharp transient spikes at known locations.
#' Spikes are ~50-100ms duration with high amplitude and sharp morphology.
#'
#' @param n_time Number of time points (default: 30000 = 60s at 500Hz)
#' @param n_channels Number of channels (default: 19)
#' @param sr Sampling rate in Hz (default: 500)
#' @param n_spikes Number of epileptic spikes to insert (default: 15)
#' @return PhysioExperiment with metadata$spike_locations (ground truth)
make_eeg_spikes <- function(n_time = 30000, n_channels = 19, sr = 500,
                            n_spikes = 15) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  ch_labels_pool <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
                       "T3", "C3", "Cz", "C4", "T4",
                       "T5", "P3", "Pz", "P4", "T6",
                       "O1", "O2")
  ch_labels <- if (n_channels <= length(ch_labels_pool)) {
    ch_labels_pool[seq_len(n_channels)]
  } else {
    paste0("EEG", seq_len(n_channels))
  }

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    alpha <- 10 * sin(2 * pi * 10 * t + runif(1, 0, 2 * pi))
    noise <- rnorm(n_time, sd = 5)
    data[, ch] <- alpha + noise
  }

  margin <- as.integer(sr * 0.5)
  spike_times <- sort(sample(margin:(n_time - margin), n_spikes))
  spike_channels <- list()

  for (i in seq_along(spike_times)) {
    sp_t <- spike_times[i]
    n_affected <- sample(1:3, 1)
    affected_ch <- sample(seq_len(n_channels), n_affected)
    spike_channels[[i]] <- affected_ch

    for (ch in affected_ch) {
      spike_width_rise <- as.integer(0.020 * sr)
      spike_width_fall <- as.integer(0.070 * sr)
      spike_amp <- runif(1, 80, 150)

      idx_rise <- seq(max(1L, sp_t - spike_width_rise), sp_t)
      idx_fall <- seq(sp_t, min(n_time, sp_t + spike_width_fall))

      data[idx_rise, ch] <- data[idx_rise, ch] +
        spike_amp * ((idx_rise - (sp_t - spike_width_rise)) / spike_width_rise)^2
      data[idx_fall, ch] <- data[idx_fall, ch] +
        spike_amp * exp(-3 * (idx_fall - sp_t) / spike_width_fall)
    }
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = ch_labels,
      type = rep("EEG", n_channels)
    ),
    samplingRate = sr
  )
  metadata(pe)$spike_locations <- data.frame(
    spike_id = seq_along(spike_times),
    sample = spike_times,
    time_sec = spike_times / sr,
    channels = I(spike_channels)
  )
  pe
}
