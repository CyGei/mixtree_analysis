# =================================================
# Miscellaneous functions
# =================================================
#' Relabel IDs by Onset Date
#'
#' Assigns new integer labels to a vector of IDs, ordered by their associated onset dates.
#' Returns a named vector mapping the original IDs to their new labels.
#'
#' @param id Character or factor vector of unique IDs.
#' @param date Vector of onset dates corresponding to \code{id}. Must be the same length as \code{id}.
#'
#' @return Character vector of new labels, named by original IDs and ordered as \code{id}.

label_ids <- function(id, date) {
  ordered_ids <- id[order(date)]
  new_labels <- as.character(seq_along(ordered_ids))
  setNames(new_labels, ordered_ids)[id]
}

#' Relabel Data Frame ID Columns
#' Replaces values in specified ID columns of a data frame with new labels ordered by onset date.
#' All ID columns are relabelled consistently using the same mapping.
#'
#' @param df Data frame. The transmission tree
#' @param id_cols Character vector of column names in \code{df} to relabel (e.g., \code{c("from", "to")}).
#' @param date_col Character. Name of the column containing onset or infection dates. Default in "date"
#'
#' @return Data frame with updated ID columns.
relabel_tree <- function(df, id_cols, date_col = "date") {
  ids <- label_ids(df[[id_cols[2]]], df[[date_col]]) # Use the 2nd column for labelling

  # Apply the relabelling all specified id columns
  for (col in id_cols) {
    df[[col]] <- ids[df[[col]]]
  }

  return(df)
}

#' Log Execution Time of a Pipeline Step
#'
#' Records the time taken to process a pipeline step and logs it to a file.
#' Returns the input unchanged to allow continued piping or assignment.
#'
#' @param .data Any R object, usually the output of a pipeline step.
#' @param label A short description of the step being timed.
#' @param log_file Path to the log file. Default is "log_time.txt".
#' @param time_unit Units for the time recorded: "secs", "millisecs", or "mins".
#'   Default is "secs".
#'
#' @return The input `.data`, unchanged and invisibly.
#'
#' @examples
#' mtcars |>
#'   dplyr::filter(cyl == 6) |>
#'   dplyr::summarise(mean_mpg = mean(mpg)) |>
#'   log_time("Filter and summarise")
#'
log_time <- function(
    .data,
    label = "pipeline",
    log_file = "log_time.txt",
    time_unit = c("secs", "millisecs", "mins")) {
  time_unit <- match.arg(time_unit)

  start <- Sys.time()
  result <- .data
  end <- Sys.time()

  duration <- as.numeric(difftime(end, start, units = "secs"))

  duration <- switch(time_unit,
    secs = duration,
    millisecs = duration * 1000,
    mins = duration / 60
  )

  entry_log <- sprintf(
    "[%s] %s: %.4f %s\n",
    format(start, "%Y-%m-%d %H:%M:%OS3"),
    label,
    duration,
    time_unit
  )

  cat(entry_log, file = log_file, append = TRUE)

  invisible(result)
}

# =================================================
# Tree generation
# =================================================

# -------------------------------------------------
# Epidemic Parameters
# -------------------------------------------------
# Relevant functions to produce the true transmission tree.

print.params <- function(x, ...) {
  cat("-------------------\n")
  cat("Epidemic Parameters\n")
  cat("-------------------\n")

  # Offspring distribution
  cat("Distributions:\n")
  cat(" Offspring:\n")
  cat("  name:", x$offspring_dist$name, "\n")
  cat("  parameters:\n")
  cat("   k (size):", x$offspring_dist$parameters$k, "\n")
  cat("   R (mu):", x$offspring_dist$parameters$R, "\n\n")

  # Generation time distribution
  cat(" Generation time:\n")
  cat("  name:", x$generation_time$name, "\n")
  cat("  parameters:\n")
  cat("   shape:", x$generation_time$parameters$shape, "\n")
  cat("   scale:", x$generation_time$parameters$scale, "\n")
  cat("   mu:", x$generation_time$parameters$mu, "\n")
  cat("   sd:", x$generation_time$parameters$sd, "\n")
  cat("   cv:", x$generation_time$parameters$cv, "\n\n")

  # Other parameters
  cat("Simulation settings:\n")
  cat(" epidemic size:", x$epidemic_size, "\n")
  cat(" simulation duration:", x$duration, "\n")

  invisible(x)
}

#' Build Epidemic Parameters
#'
#' Constructs a list of parameters for simulating the true transmission tree.
#'
#' @param off_R Numeric. Mean of the negative binomial offspring distribution.
#' @param off_k Numeric. Dispersion parameter of the negative binomial offspring distribution.
#' @param gt_mu Numeric. Mean generation time (gamma distribution).
#' @param gt_sd Numeric. Standard deviation of generation time (gamma distribution).
#' @param epidemic_size Integer. The epidemic's final size (number of cases).
#' @param duration Integer. Duration of the simulation (in days).
#'
#' @return An object of class \code{"params"} containing simulation settings.
build_params <- function(off_R, off_k, gt_mu, gt_sd, epidemic_size, duration) {
  offspring_dist <- list(
    name = "nbinom",
    parameters = list(k = off_k, R = off_R),
    r = function(n = 1) rnbinom(n = n, size = off_k, mu = off_R),
    d = function(x) dnbinom(x, size = off_k, mu = off_R)
  )


  gt_cv <- gt_sd / gt_mu
  gt_params <- epitrix::gamma_mucv2shapescale(mu = gt_mu, cv = gt_cv)
  generation_time <- distcrete::distcrete(
    name = "gamma",
    shape = gt_params$shape,
    scale = gt_params$scale,
    interval = 1
  )
  generation_time$parameters <- list(
    mu = gt_mu,
    sd = gt_sd,
    cv = gt_cv,
    shape = gt_params$shape,
    scale = gt_params$scale
  )

  params <- list(
    offspring_dist = offspring_dist,
    generation_time = generation_time,
    epidemic_size = epidemic_size,
    duration = duration
  )
  class(params) <- "params"

  return(params)
}


# -------------------------------------------------
# True transmission tree
# -------------------------------------------------
#' Simulate a True Transmission Tree
#'
#' Simulates a single outbreak using specified epidemic parameters until exactly epidemic_size is reached.
#'
#' @param params List of epidemic parameters generated by \code{build_params()}.
#'
#' @return A tibble with columns \code{from}, \code{to}, and \code{date}, representing transmission events.
#'
build_tree <- function(params) {
  expected_AR <- epitrix::R02AR(params$offspring_dist$parameters$R)

  repeat {
    out <- simulacr::simulate_outbreak(
      duration = params$duration,
      population_size = round(params$epidemic_size / expected_AR),
      R_values = params$offspring_dist$r(params$epidemic_size),
      dist_generation_time = params$generation_time$d(0:params$duration),
      dist_incubation = 1L
    )

    # Check if we have exactly the target number of cases
    if (nrow(out$data) == params$epidemic_size) break
  }

  tree <- out$data |>
    transmute(
      from = source,
      to = id,
      date = date_infection # row_number() @CyGei Do we have to avoid same dates?
    ) |>
    relabel_tree(id_cols = c("from", "to"), date_col = "date") |>
    as_tibble()

  mixtree::validate_tree(tree |> slice(-1))

  return(tree)
}

# -------------------------------------------------
# Epidemic forest
# -------------------------------------------------
#' Simulate an Epidemic Forest
#'
#' Generates multiple plausible transmission trees given \code{tree} and \code{params}.
#'
#' @param tree A tibble representing the reference tree (from \code{build_tree()}).
#' @param params List of epidemic parameters generated by \code{build_params()}.
#' @param forest_size Integer. The total number of transmission trees to generate. Default is 200.
#'
#' @return A list of transmission trees (tibbles), each with columns \code{from}, \code{to}, and \code{date}.
build_forest <- function(tree, params, forest_size = 200L) {
  simulate_tree <- function() {
    sim_tree <- tree |>
      rename(j = from, i = to) |>
      mutate(
        j = NA_character_,
        R = params$offspring_dist$r(n()),

        # Ensure that the first case (i == "1") has R >= 1
        R = if_else(
          i == "1" & R < 1,
          LaplacesDemon::rtrunc(
            n = 1,
            spec = params$offspring_dist$name,
            a = 0,
            size = params$offspring_dist$parameters$k,
            mu = params$offspring_dist$parameters$R
          ), # Sample from non-zero R values
          R
        )
      )

    generation_time_pmf <- params$generation_time$d(0:params$duration)

    for (row in 2:nrow(sim_tree)) {
      i <- sim_tree$i[row]
      t_i <- sim_tree$date[row]

      j_candidates <- sim_tree |>
        filter(date < t_i) |>
        mutate(foi = R * generation_time_pmf[t_i - date + 1])
      # |>
      #   filter(foi > 0)

      # if (nrow(j_candidates) == 0) next

      j_chosen <- sample(j_candidates$i, 1, prob = j_candidates$foi / sum(j_candidates$foi))
      sim_tree$j[row] <- j_chosen
    }
    # safety checks
    stopifnot(sum(is.na(sim_tree$j)) == 1) # one intro only

    sim_tree <- sim_tree |>
      slice(-1) |>
      rename(from = j, to = i) # ignore intro row

    mixtree::validate_tree(sim_tree)
    return(sim_tree)
  }

  forest <- replicate(forest_size, simulate_tree(), simplify = FALSE)
  return(forest)
}


# -------------------------------------------------
# Forest checks
# -------------------------------------------------

# Function to check consistency and variation in columns across all data frames in the list
check_forest <- function(forest) {
  # Extract 'to' and 'date' columns from the first data frame as a reference
  ref_to <- forest[[1]]$to
  ref_date <- forest[[1]]$date
  ref_R <- forest[[1]]$R
  ref_from <- forest[[1]]$from

  # Check if 'to' and 'date' are consistent in all data frames
  consistent_to_date <- all(
    sapply(forest, function(df) {
      all(df$to == ref_to) && all(df$date == ref_date)
    })
  )

  # Check if 'R' and 'from' columns vary across data frames
  vary_R_from <- any(
    sapply(forest, function(df) {
      !all(df$R == ref_R) || !all(df$from == ref_from)
    })
  )

  # Return a list with both checks
  checks <- list(
    consistent_to_date = consistent_to_date,
    vary_R_from = vary_R_from
  )

  return(checks)
}
