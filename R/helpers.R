# Outbreak Simulation -----------------------------------------------------
#' Simulate outbreaks within a target size range
#'
#' This function repeatedly simulates outbreaks until a specified number of simulations
#' is obtained, ensuring each simulation's outbreak size falls within a tolerance range around
#' a target size.
#'
#' @param target_size Integer. The target outbreak size.
#' @param R_values Integer vector. Reproduction number(s) used in the simulation.
#' @param tolerance Numeric. Tolerance for outbreak size (default is 0.2 for ±20% of target size).
#' @param n_simulations Integer. Number of valid simulations to generate (default is 200).
#' @param max_attempts Integer. Maximum number of simulation attempts before checking progress (default is 1000).
#' @param try_gain Numeric. Fraction of target simulations required to extend max_attempts (default is 1, i.e. no extension).
#'
#' @return A list of outbreak simulations (data frames) whose sizes are within
#'   \code{target_size * (1 - tolerance)} and \code{target_size * (1 + tolerance)}.
#'
#' @details The function simulates outbreaks using \code{simulacr::simulate_outbreak} until
#'   \code{n_simulations} valid simulations are collected. If the number of attempts reaches
#'   \code{max_attempts} and at least \code{try_gain * n_simulations} simulations have been obtained,
#'   an additional 1000 attempts are added automatically; otherwise, the function issues a warning.
#'
#' @examples
#' \dontrun{
#' # Simulate 100 outbreaks targeting 50 cases ±20% using negative binomial R_values.
#' sims <- simulate_outbreaks(
#'   target_size = 50,
#'   R_values = rnbinom(100, size = 0.2, mu = 3),
#'   tolerance = 0.2,
#'   n_simulations = 100,
#'   max_attempts = 1000,
#'   try_gain = 0.8
#' )
#' }
#'
simulate_outbreaks <- function(target_size,
                               R_values,
                               tolerance = 0.2,
                               n_simulations = 200,
                               max_attempts = 1000,
                               try_gain = 1) {
  # Calculate acceptable size range
  min_size <- round(target_size * (1 - tolerance))
  max_size <- round(target_size * (1 + tolerance))
  sims <- list()
  attempts <- 0

  while (length(sims) < n_simulations) {
    # Check if current max_attempts limit is reached
    if (attempts >= max_attempts) {
      if (length(sims) >= try_gain * n_simulations) {
        message(
          "Reached ",
          try_gain * 100,
          "% of target simulations. Adding an additional 1000 attempts."
        )
        max_attempts <- max_attempts + 1000
      } else {
        warning("Maximum number of attempts reached at ",
                length(sims),
                "simulations.")
        break
      }
    }

    sim <- simulacr::simulate_outbreak(
      duration = 100,
      population_size = round(target_size / epitrix::R02AR(mean(R_values))),
      R_values = R_values,
      dist_incubation = outbreaker2::fake_outbreak$w,
      dist_generation_time = outbreaker2::fake_outbreak$w
    )$data %>%
      relabel_tree(id_cols = c("id", "source"), date_col = "date_onset")
    sim <- shift_init_date(sim)

    attempts <- attempts + 1

    # Accept simulation if within acceptable size range
    if (nrow(sim) >= min_size && nrow(sim) <= max_size) {
      sims[[length(sims) + 1]] <- sim
    }
  }

  return(sims)
}


# Process Simulations -----------------------------------------------------
#' @title Match data frames by row counts
#'
#' @description Matches and pairs data frames from lists \code{A} and \code{B} based on identical row counts.
#'
#' @param A A list of data frames.
#' @param B A list of data frames.
#' @param n Number of pairs to return.
#'
#' @return A list of length \code{n}, each element containing a pair of data frames with the same number of rows.
matching_pairs <- function(A, B, n = 100) {
  nrows_A <- sapply(A, nrow)
  nrows_B <- sapply(B, nrow)

  dfA <- data.frame(index_A = seq_along(A), nrows = nrows_A)
  dfB <- data.frame(index_B = seq_along(B), nrows = nrows_B)

  merged_df <- merge(dfA, dfB, by = "nrows")

  # keep unique index from A @ASK if necessary
  # merged_df <- merged_df[!duplicated(merged_df$index_A), ]


  if (nrow(merged_df) > n) {
    merged_df <- merged_df[sample(1:nrow(merged_df), n), ]
  }

  paired_list <- mapply(function(a_idx, b_idx) {
    list(A = A[[a_idx]], B = B[[b_idx]])
  }, merged_df$index_A, merged_df$index_B, SIMPLIFY = FALSE)

  if (length(paired_list) < n) {
    stop("Not enough matching pairs found to create ", n, "pairs.")
  }

  # Check that each pair has the same number of rows
  rows_match <- sapply(paired_list, function(pair)
    nrow(pair$A) == nrow(pair$B))
  if (!all(rows_match)) {
    stop("Not all pairs have matching row counts.")
  }

  return(paired_list)
}


#' Shift Duplicate Earliest Dates
#'
#' Shifts onset dates that are equal to the earliest date (except the first occurrence) by a specified value.
#'
#' @param sim Data frame containing simulation data.
#' @param date_col Character. Name of the column with onset dates. Default is "date_onset".
#' @param shift Numeric. Value to add to duplicate earliest dates. Default is 1.
#'
#' @return Data frame with updated onset dates.
shift_init_date <- function(sim,
                            date_col = "date_onset",
                            shift = 1) {
  dates <- sim[[date_col]]
  min_date <- min(dates, na.rm = TRUE)
  idx <- which(dates == min_date)
  # Shift all but the first occurrence
  if (length(idx) > 1) {
    sim[[date_col]][idx[-1]] <- sim[[date_col]][idx[-1]] + shift
  }
  return(sim)
}

#' Relabel IDs by Onset Date Order
#'
#' Reorders and relabels IDs based on the corresponding onset dates.
#'
#' @param id Vector of IDs.
#' @param date Vector of dates corresponding to the IDs.
#'
#' @return Named character vector mapping original IDs to new labels.
label_ids <- function(id, date) {
  ordered_ids <- id[order(date)]
  new_labels <- as.character(seq_along(ordered_ids))
  setNames(new_labels, ordered_ids)[id]
}

#' Relabel Data Frame ID Columns
#'
#' Relabels specified ID columns in a data frame using new labels ordered by onset dates.
#'
#' @param df Data frame containing outbreak data.
#' @param id_cols Character vector of column names to be relabelled.
#' @param date_col Character. Name of the column with onset dates. Default is "date_onset".
#'
#' @return Data frame with relabelled ID columns.
relabel_tree <- function(df, id_cols, date_col = "date_onset") {
  ids <- label_ids(df[[id_cols[1]]], df[[date_col]]) # Use the first column for labelling

  # Apply the relabelling all specified id columns
  for (col in id_cols) {
    df[[col]] <- ids[df[[col]]]
  }

  return(df)
}


# Outbreak Reconstruction -------------------------------------------------
#' Run Outbreaker2 on Simulation Data
#'
#' Executes the outbreaker2 algorithm on a simulation data frame, utilising additional contact data.
#'
#' @param sim Data frame containing outbreak simulation data.
#' @param ctd_fraction Numeric. Fraction of cases to sample for contact data. Default is 0.5.
#'
#' @return An outbreaker2 object containing the reconstructed outbreak.
run_outbreaker <- function(sim, ctd_fraction = 0.5) {
  n_cases <- nrow(sim)
  # Adding ctd data to improve reconstruction
  ctd_sample <- sample(2:n_cases, round(n_cases * ctd_fraction))

  config <- outbreaker2::create_config(
    init_tree = "star",
    move_mu = FALSE,
    init_pi = 1,
    move_pi = FALSE,
    init_kappa = 1,
    move_kappa = FALSE,
    find_import = FALSE
  )

  data <- outbreaker2::outbreaker_data(
    dates = sim$date_onset,
    w_dens = outbreaker2::fake_outbreak$w,
    f_dens = outbreaker2::fake_outbreak$w,
    ctd = sim[ctd_sample, c("source", "id")],
    ids = sim$id
  )

  out <- outbreaker2::outbreaker(data = data, config = config)
  out <- o2ools::identify(out, sim$id)
  return(out)
}


#' Build a Chain of Posterior Transmission Trees from a Simulation
#'
#' Runs outbreaker2 on a simulation dataset, removes burn-in steps,
#' extracts transmission trees, and processes them into a format suitable for use with
#' `epitree::compare_chain`.
#'
#' @param sim Data frame containing outbreak simulation data.
#' @param ctd_fraction Numeric. Fraction of cases to sample for contact data. Default is 0.5.
#' @param burnin Integer. Step threshold to remove burn-in samples. Default is 1000.
#'
#' @return A list of processed posterior transmission trees.
build_chain <- function(sim,
                        ctd_fraction = 0.5,
                        burnin = 1000) {
  chain <- run_outbreaker(sim, ctd_fraction = ctd_fraction)
  chain <- chain[chain$step > burnin, ]
  trees <- o2ools::get_trees(chain)
  lapply(trees, epitree:::process_tree)
}


# Perform the test  -------------------------------------------------

#' Compute P-values from Resampled Chains
#'
#' Draws resampled chains from two input chains (`chainA` and `chainB`) based on specified
#' probabilities, then compares each resampled pair using `epitree::compare_chains` to
#' compute p-values.
#'
#' @param chainA A vector or list representing the first chain.
#' @param chainB A vector or list representing the second chain.
#' @param overlap_freq Numeric. The probability weight for sampling from `chainA` in the mixed chain.
#'   The weight for `chainB` is computed as (1 - overlap_freq).
#' @param sample_size Integer. The number of samples to draw per replicate.
#' @param n_repeats Integer. The number of replicates to perform.
#' @param method Character. The method to use for comparing chains. Default is "adonis", alternative is "chisq".
#' @param args List. Additional arguments to pass to `epitree` functions.
#'
#' @return A numeric vector of p-values, one for each replicate.

get_pval <- function(chainA,
                     chainB,
                     overlap_freq,
                     sample_size,
                     n_repeats,
                     method = "adonis",
                     args = list()) {
  # Probabilities
  lenA <- length(chainA)
  lenB <- length(chainB)
  total_chain <- c(chainA, chainB)
  total_len <- lenA + lenB

  pA <- overlap_freq
  pB <- 1 - pA
  probs <- c(rep(pA / lenA, lenA), rep(pB / lenB, lenB))


  # Draw sample indices for all replicates in one go
  mix_idx <- matrix(
    sample.int(
      total_len,
      size = sample_size * n_repeats,
      prob = probs,
      replace = TRUE
    ),
    nrow = n_repeats
  )

  ref_idx <- matrix(sample.int(lenA, size = sample_size * n_repeats, replace = TRUE),
                    nrow = n_repeats)
  if (method == "adonis") {
    p_values <- vapply(seq_len(n_repeats), function(i) {
      mixed_chain <- total_chain[mix_idx[i, ]]
      reference_chain <- chainA[ref_idx[i, ]]
      suppressWarnings(epitree::compare_chains(reference_chain, mixed_chain, adonis2_args = args)[, "Pr(>F)"][1])
    }, numeric(1))
  } else if (method == "chisq") {
    p_values <- vapply(seq_len(n_repeats), function(i) {
      mixed_chain <- total_chain[mix_idx[i, ]]
      reference_chain <- chainA[ref_idx[i, ]]
      suppressWarnings(epitree::get_chisq(reference_chain, mixed_chain, test_args = args)[["p.value"]])
    }, numeric(1))
  } else {
    stop("Method not recognised")
  }

  return(p_values)
}




# Plot --------------------------------------------------------------------

#' Plot Proportions of Rejected Null Hypotheses
#' @param results A data frame containing the results of the peformance study.
#' @param alpha Numeric. The significance level. Default is 0.05.
#' @return A ggplot object.


plot_props <- function(results, alpha = 0.05) {
  annotations_df <- results %>%
    group_by(epidemic_size, sample_size, overlap_freq) %>%
    summarise(
      reject = mean(p_value < alpha) * 100,
      accept = mean(p_value >= alpha) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c("reject", "accept"),
      names_to = "category",
      values_to = "percentage"
    ) %>%
    mutate(
      row_title = "Epidemic Size",
      col_title = "Posterior Sample Size",
      # Create a label column with clear logic
      label = case_when(
        round(percentage, 1) == 100 ~ as.character(round(percentage, 1)),
        round(percentage, 1) == 0   ~ "",
        TRUE                        ~ as.character(round(percentage, 1))
      )
    )

  p <- annotations_df %>%
    ggplot(aes(x = overlap_freq, y = percentage, fill = category)) +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.775,
      color = "white",
      linewidth = 0.1
    ) +
    geom_text(aes(
      y = ifelse(
        percentage == 100,
        50,
        # Center the label when percentage is 100%
        ifelse(category == "reject", percentage / 2, 98.25 - percentage / 2)
      ),
      label = label
    ),
    colour = "black",
    size = 3) +
    ggh4x::facet_nested(
      rows = vars(row_title, epidemic_size),
      cols = vars(col_title, sample_size),
      scales = "fixed",
      space = "fixed",
      remove_labels = "none",
      nest_line = element_line(
        colour = "black",
        linewidth = 0.05,
        lineend = "square"
      )
    ) +
    scale_fill_manual(
      values = c("reject" = "#CCCDC6", "accept" = "#49b56d"),
      #7d8af0 #82750f #0F8233
      breaks = c("reject", "accept"),
      name = "",
      labels = c("reject" = "Reject H0", "accept" = "Accept H0")
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = -15),
      strip.background = element_rect(
        fill = "#f8f4f2",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = element_text(family = "Fira Code", size = 10),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
      # panel.spacing = unit(0.25, "lines")
    ) +
    guides(fill =
             guide_legend(
               keywidth = 5,
               keyheight = 1,
               label.position = "top",
               label.vjust = -9.5
             )) +
    labs(x = "Overlap Frequency", y = "Percentage")

  return(p)
}




# ROC ---------------------------------------------------------------------
get_roc <- function(data,
                    group_vars = c("sample_size", "epidemic_size")) {
  data %>%
    mutate(overlap_freq = as.integer(as.character(overlap_freq)),
           truth = if_else(overlap_freq < 1, 1, 0)) %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ {
      roc_obj <- pROC::roc(
        response = .x$truth,
        predictor = .x$p_value,
        direction = ">",
        quiet = TRUE
      )

      roc_coords <- pROC::coords(
        roc_obj,
        "all",
        ret = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      )

      youdens_j <- roc_coords$sensitivity + roc_coords$specificity - 1
      best_index <- which.max(youdens_j)
      best_alpha_value <- roc_coords$threshold[best_index]

      tibble::tibble(
        alpha = roc_coords$threshold,
        sensitivity = roc_coords$sensitivity,
        specificity = roc_coords$specificity,
        AUC = as.numeric(pROC::auc(roc_obj)),
        best_alpha = roc_coords$threshold == best_alpha_value
      )

    }) %>%
    ungroup()
}



plot_roc <- function(df,
                     group_vars = NULL,
                     facet_rows = NULL,
                     facet_cols = NULL,
                     row_title = NULL,
                     col_title = NULL,
                     color_var = "sample_size",
                     color_palette = "Blues",
                     point_shapes = c("Alpha = 0.05" = 21, "Best Alpha" = 24),
                     alpha_values = c(0.05),
                     title = NULL,
                     subtitle = NULL,
                     x_label = "1 - Specificity",
                     y_label = "Sensitivity") {

  # Get ROC data if raw results are provided
  if("AUC" %in% names(df)) {
    x <- df
  } else {
    if(is.null(group_vars)) {
      stop("If providing raw results, you must specify group_vars")
    }
    x <- get_roc(df, group_vars = group_vars)
  }

  # Set up facet variables
  if(!is.null(row_title) && !is.null(facet_rows)) {
    x <- x %>% mutate(row_title = row_title)
  }

  if(!is.null(col_title) && !is.null(facet_cols)) {
    x <- x %>% mutate(col_title = col_title)
  }

  # Create points dataframe for multiple alpha values
  points_df_list <- list()

  # Add points for specified alpha values
  for(a in alpha_values) {
    points_df_list[[paste0("Alpha = ", a)]] <- x %>%
      group_by_at(vars(one_of(c(color_var, facet_rows, facet_cols)))) %>%
      slice_min(abs(alpha - a)) %>%
      ungroup() %>%
      mutate(point_type = paste0("Alpha = ", a))
  }

  # Add points for best alpha
  points_df_list[["Best Alpha"]] <- x %>%
    filter(best_alpha) %>%
    mutate(point_type = "Best Alpha")

  points_df <- bind_rows(points_df_list)

  # Create annotations
  group_vars_for_annotation <- c()
  if(!is.null(col_title) && !is.null(facet_cols)) {
    group_vars_for_annotation <- c(group_vars_for_annotation, "col_title", facet_cols)
  } else if(!is.null(facet_cols)) {
    group_vars_for_annotation <- c(group_vars_for_annotation, facet_cols)
  }

  if(!is.null(row_title) && !is.null(facet_rows)) {
    group_vars_for_annotation <- c(group_vars_for_annotation, "row_title", facet_rows)
  } else if(!is.null(facet_rows)) {
    group_vars_for_annotation <- c(group_vars_for_annotation, facet_rows)
  }

  if(length(group_vars_for_annotation) == 0) {
    group_vars_for_annotation <- color_var
  }

  annotations <- x %>%
    filter(best_alpha) %>%
    group_by_at(vars(one_of(c(color_var, group_vars_for_annotation)))) %>%
    summarise(
      AUC = unique(AUC),
      AUC = round(AUC, digits = 3),
      best_alpha_value = unique(alpha),
      best_alpha_value = round(best_alpha_value, digits = 3),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0('"AUC = " * ', AUC, ' * "," ~ ', 'alpha^"*" == ', best_alpha_value),
      x = 0.3,
      y = 0.05 + (as.integer(get(color_var)) - 1) * 0.075
    )

  # Prepare faceting formula
  facet_formula <- NULL
  if(!is.null(facet_rows) && !is.null(facet_cols)) {
    if(!is.null(row_title) && !is.null(col_title)) {
      facet_formula <- ggh4x::facet_nested(
        rows = vars(row_title, !!!syms(facet_rows)),
        cols = vars(col_title, !!!syms(facet_cols)),
        scales = "fixed",
        space = "fixed"
      )
    } else {
      facet_formula <- ggh4x::facet_nested(
        rows = vars(!!!syms(facet_rows)),
        cols = vars(!!!syms(facet_cols)),
        scales = "fixed",
        space = "fixed"
      )
    }
  } else if(!is.null(facet_rows)) {
    if(!is.null(row_title)) {
      facet_formula <- ggh4x::facet_nested(
        rows = vars(row_title, !!!syms(facet_rows)),
        scales = "fixed",
        space = "fixed"
      )
    } else {
      facet_formula <- ggh4x::facet_nested(
        rows = vars(!!!syms(facet_rows)),
        scales = "fixed",
        space = "fixed"
      )
    }
  } else if(!is.null(facet_cols)) {
    if(!is.null(col_title)) {
      facet_formula <- ggh4x::facet_nested(
        cols = vars(col_title, !!!syms(facet_cols)),
        scales = "fixed",
        space = "fixed"
      )
    } else {
      facet_formula <- ggh4x::facet_nested(
        cols = vars(!!!syms(facet_cols)),
        scales = "fixed",
        space = "fixed"
      )
    }
  }

  # Build plot
  p <- ggplot()

  if(!is.null(facet_formula)) {
    p <- p + facet_formula
  }

  p <- p +
    geom_line(data = x, aes_string(
      x = "1 - specificity",
      y = "sensitivity",
      color = color_var
    )) +
    geom_point(
      data = points_df,
      aes_string(
        x = "1 - specificity",
        y = "sensitivity",
        fill = color_var,
        shape = "point_type"
      ),
      color = "black",
      stroke = 0.5
    ) +
    geom_text(
      data = annotations,
      aes_string(
        x = "x",
        y = "y",
        color = color_var,
        label = "label"
      ),
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3,
      show.legend = FALSE
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey"
    ) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = -5),
      strip.background = element_rect(
        fill = "#f8f4f2",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = element_text(family = "Fira Code", size = 10),
      panel.grid = element_blank()
    ) +
    guides(
      color = guide_legend(
        title = gsub("_", " ", tools::toTitleCase(color_var)),
        order = 1,
        override.aes = list(shape = NA)
      ),
      shape = guide_legend(title = "Point Type", order = 2),
      fill = "none"
    ) +
    colorspace::scale_color_discrete_sequential(palette = color_palette, nmax = 8, order = 5:8) +
    colorspace::scale_fill_discrete_sequential(palette = color_palette, nmax = 8, order = 5:8) +
    scale_shape_manual(values = point_shapes) +
    labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle
    )

  return(p)
}


# plot_roc(rocs,
#          facet_rows = "method",
#          facet_cols = "epidemic_size",
#          row_title = "Method",
#          col_title = "Epidemic Size",
#          color_var = "sample_size",
#          color_palette = "Purples")
#
# plot_roc(chi_results,
#          group_vars = c("sample_size", "epidemic_size"),
#          facet_cols = "epidemic_size",
#          col_title = "Epidemic Size",
#          color_var = "sample_size")


# plot_roc <- function(df) {
#   x <- get_roc(df) %>%
#     mutate(col_title = "Epidemic Size")
#
#   points_df <- bind_rows(
#     x %>%
#       group_by(sample_size, epidemic_size) %>%
#       slice_min(abs(alpha - 0.05)) %>%
#       ungroup() %>%
#       mutate(point_type = "Alpha = 0.05"),
#     x %>%
#       filter(best_alpha) %>%
#       mutate(point_type = "Best Alpha")
#   )
#
#   annotations <- x %>%
#     filter(best_alpha) %>%
#     group_by(col_title, sample_size, epidemic_size) %>%
#     summarise(
#       AUC = unique(AUC),
#       AUC = round(AUC, digits = 3),
#       best_alpha_value = unique(alpha),
#       best_alpha_value = round(best_alpha_value, digits = 3),
#       .groups = "drop"
#     ) %>%
#     mutate(
#       label = paste0('"AUC = " * ', AUC, ' * "," ~ ', 'alpha^"*" == ', best_alpha_value),
#       x = 0.3,
#       y = 0.05 + (as.integer(sample_size) - 1) * 0.075
#     )
#
#
#   p <- ggplot() +
#     ggh4x::facet_nested(
#       cols = vars(col_title, epidemic_size),
#       scales = "fixed",
#       space = "fixed"
#     ) +
#     geom_line(data = x,
#               aes(
#                 x = 1 - specificity,
#                 y = sensitivity,
#                 color = sample_size
#               )) +
#     geom_point(
#       data = points_df,
#       aes(
#         x = 1 - specificity,
#         y = sensitivity,
#         fill = sample_size,
#         shape = point_type
#       ),
#       color = "black",
#       stroke = 0.5
#     ) +
#     geom_text(
#       data = annotations,
#       aes(
#         x = x,
#         y = y,
#         color = sample_size,
#         label = label
#       ),
#       parse = TRUE,
#       hjust = 0,
#       vjust = 1,
#       size = 3,
#       show.legend = FALSE
#     ) +
#     geom_abline(
#       slope = 1,
#       intercept = 0,
#       linetype = "dashed",
#       color  = "grey"
#     ) +
#     coord_fixed(ratio = 1) +
#     theme_bw() +
#     theme(
#       legend.position = "bottom",
#       legend.box = "vertical",
#       legend.margin = margin(t = -5),
#       strip.background = element_rect(
#         fill = "#f8f4f2",
#         colour = "black",
#         linewidth = 0.5
#       ),
#       strip.text = element_text(family = "Fira Code", size = 10),
#       panel.grid = element_blank()
#     ) +
#     guides(
#       color = guide_legend(
#         title = "Sample Size",
#         order = 1,
#         override.aes = list(shape = NA)
#       ),
#       shape = guide_legend(title = "Point Type", order = 2),
#       fill = "none"
#     ) +
#     colorspace::scale_color_discrete_sequential(palette = "Blues",
#                                                 nmax = 8,
#                                                 order = 5:8) +
#     colorspace::scale_fill_discrete_sequential(palette = "Blues",
#                                                nmax = 8,
#                                                order = 5:8) +
#     scale_shape_manual(values = c("Alpha = 0.05" = 21, "Best Alpha" = 24)) +
#     labs(x = "1 - Specificity", y = "Sensitivity")
#
#   return(p)
# }
