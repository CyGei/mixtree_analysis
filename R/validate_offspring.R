source("R/packages.R")
source("R/helpers.R")

# -------------------------------------------------
# Offspring estimation pipeline:
# tree |> right_censor(params) |> count_offspring() |> fit_negbin()
# -------------------------------------------------

#' Annotate censored cases.
#' Adds a `censored` column: TRUE for cases infected too late
#' to have realised all their secondary infections.
#' The cutoff is the smallest t where the generation time CDF >= 0.99.
#' NB: we annotate, not remove. Censored cases must remain as children
#' when counting offspring of uncensored parents.
right_censor <- function(tree, params, threshold = 0.99) {
  pmf <- params$generation_time$d(0:params$duration)
  cutoff <- which(cumsum(pmf) >= threshold)[1]
  t_max <- max(tree$date)
  tree |> mutate(censored = date > t_max - cutoff)
}

#' Count offspring per case using the full tree.
#' Every row (including censored children) contributes to
#' their parent's offspring count.
count_offspring <- function(tree) {
  offspring_counts <- tree |>
    filter(!is.na(from)) |>
    count(from, name = "offspring")

  tree |>
    left_join(offspring_counts, by = c("to" = "from")) |>
    mutate(offspring = replace_na(offspring, 0L))
}

#' Fit NegBin on uncensored cases only.
fit_negbin <- function(tree) {
  x <- tree |> filter(!censored) |> pull(offspring)
  n <- length(x) # number of uncensored cases used to fit the NegBin
  fit <- tryCatch(
    MASS::fitdistr(
      x,
      "negative binomial",
      start = list(mu = mean(x), size = 0.7)
    ),
    error = \(e) NULL
  )
  tibble(
    mu = if (!is.null(fit)) unname(fit$estimate["mu"]) else NA_real_,
    k = if (!is.null(fit)) unname(fit$estimate["size"]) else NA_real_,
    n = n
  )
}

#' Convert each row of a forest matrix back to a tree tibble
#' by grafting new ancestries onto the reference tree's IDs and dates.
forest_to_trees <- function(forest, ref_tree) {
  lapply(seq_len(nrow(forest)), \(i) {
    tibble(
      from = c(NA_character_, forest[i, ]),
      to = ref_tree$to,
      date = ref_tree$date
    )
  })
}

# # =====================================================
# #                      PIPELINE
# # =====================================================

# #' Run the full forest pipeline for one (R, k) combination.
# #'
# #' For each of independent reference trees (all simulated from the
# #' *same* (off_R, off_k) parameters), a forest of forest_size trees is
# #' generated and fitted. Using multiple reference trees serves two purposes:
# #'   1. It averages out idiosyncrasies of any single reference tree's date
# #'      structure, which can otherwise bias all forest estimates in the same
# #'      direction (since every forest tree shares those dates).
# #'   2. It reveals whether any remaining bias is systematic (algorithm
# #'      limitation) rather than a one-off sampling artefact.
# #'
# #' Returns a tibble with per-tree NegBin estimates, the reference replicate
# #' ID, and the true params.
# pipeline <- function(
#   off_R,
#   off_k,
#   gt_mu = 12,
#   gt_sd = 6,
#   epidemic_size = 200L,
#   duration = 365L,
#   forest_size = 50L
# ) {
#   fp <- build_params(
#     off_R = off_R,
#     off_k = off_k,
#     gt_mu = gt_mu,
#     gt_sd = gt_sd,
#     epidemic_size = epidemic_size,
#     duration = duration,
#     replicates = 50L
#   )
#   ref_trees <- build_tree(fp) # list of length replicates

#   imap_dfr(ref_trees, \(ref_tree, ref_id) {
#     forest <- build_forest(ref_tree, fp, forest_size = forest_size)

#     forest |>
#       forest_to_trees(ref_tree) |>
#       map_dfr(
#         \(tree) {
#           tree |>
#             right_censor(fp) |>
#             count_offspring() |>
#             fit_negbin()
#         },
#         .id = "tree_id"
#       ) |>
#       mutate(ref_id = ref_id)
#   }) |>
#     mutate(true_R = off_R, true_k = off_k)
# }

# # ------------------------------------
# #           Sweep over (R, k) combos
# # ------------------------------------
# # Each call to pipeline() generates its own reference tree from the target
# # (off_R, off_k), so date structure and forest parameters are always matched.
# param_grid <- expand_grid(
#   off_R = c(1.5, 2, 3),
#   off_k = c(0.1, 0.5, 3)
# )
# set.seed(123)
# plan(multisession, workers = availableCores() - 2)
# results <- furrr::future_pmap_dfr(
#   param_grid,
#   \(off_R, off_k) {
#     pipeline(off_R = off_R, off_k = off_k)
#   },
#   .options = furrr_options(seed = TRUE)
# ) |>
#   time_pipe("pipeline execution") #12min
# plan(sequential)
# #saveRDS(results, "data/sanity.rds")
# results <- readRDS("data/sanity.rds")
# results |>
#   group_by(true_R, true_k) |>
#   summarise(
#     mean_mu = mean(mu, na.rm = TRUE),
#     mean_k = mean(k, na.rm = TRUE),
#     mean_n = mean(n, na.rm = TRUE), # mean uncensored cases per tree fit
#     n_trees = n()
#   )

# # ------------------------------------
# #           Plot: R recovery
# # ------------------------------------
# ggplot(results, aes(mu)) +
#   geom_histogram(fill = "grey80", colour = "grey40", bins = 50) +
#   geom_vline(aes(xintercept = true_R), linetype = "dashed", colour = "red") +
#   facet_grid(
#     rows = vars(paste0("k = ", true_k)),
#     cols = vars(paste0("R = ", true_R)),
#     scales = "free_x"
#   ) +
#   scale_x_continuous(limits = c(0, 10)) +
#   labs(
#     x = "Estimated R",
#     y = "Count"
#   ) +
#   theme_classic()

# ggsave(
#   filename = "figures/mu_recovery.eps",
#   width = 7.5,
#   height = 7.5,
#   units = "in",
#   device = cairo_ps,
#   dpi = 300
# )

# # ------------------------------------
# #           Plot: k recovery
# # ------------------------------------
# ggplot(results, aes(k)) +
#   geom_histogram(fill = "grey80", colour = "grey40", binwidth = 0.1) +
#   geom_vline(aes(xintercept = true_k), linetype = "dashed", colour = "red") +
#   facet_grid(
#     rows = vars(paste0("k = ", true_k)),
#     cols = vars(paste0("R = ", true_R)),
#     scales = "free_x"
#   ) +
#   scale_x_continuous(limits = c(0, 10)) +
#   labs(
#     x = "Estimated k",
#     y = "Count"
#   ) +
#   theme_classic()

# ggsave(
#   filename = "figures/k_recovery.eps",
#   width = 7.5,
#   height = 7.5,
#   units = "in",
#   device = cairo_ps,
#   dpi = 300
# )

# ggplot(results, aes(n)) +
#   geom_histogram(fill = "grey80", colour = "grey40", bins = 50) +
#   facet_grid(
#     rows = vars(paste0("k = ", true_k)),
#     cols = vars(paste0("R = ", true_R)),
#     scales = "fixed"
#   ) +
#   scale_x_continuous(limits = c(0, 120)) +
#   labs(
#     x = "Number of uncensored cases used to fit NegBin",
#     y = "Count"
#   ) +
#   theme_classic()

# ggsave(
#   filename = "figures/n_recovery.eps",
#   width = 7.5,
#   height = 7.5,
#   units = "in",
#   device = cairo_ps,
#   dpi = 300
# )
