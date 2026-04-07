source("R/packages.R")
source("R/helpers.R")

# =================================================
# QQ & Deviation Plots: Reference Tree vs Forest
# + Empirical vs Theoretical Offspring Distribution
# =================================================
# Single pipeline: for each (R, k), simulate reference trees, build forests,
# and return both QQ data and raw reference offspring counts.
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
  tree |>
    mutate(offspring = as.integer(table(factor(from, levels = to))[to]))
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
#' Extract uncensored offspring counts from a tree.
get_offspring <- function(tree, params) {
  tree |>
    right_censor(params) |>
    count_offspring() |>
    filter(!censored) |>
    pull(offspring)
}

#' Shared facet labels and theme for all panels.
facet_Rk <- function(scales = "free") {
  list(
    facet_grid(
      rows = vars(paste0("k = ", true_k)),
      cols = vars(paste0("R = ", true_R)),
      scales = scales
    ),
    theme_classic(),
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
  )
}

#' Unified pipeline: for a given (R, k), simulate reference trees and forests.
#' Returns a list with two tibbles:
#'   - qq:  per-tree quantile comparisons (for QQ and deviation plots)
#'   - ref: raw uncensored offspring per reference tree (for empirical vs theoretical)
pipeline <- function(
  off_R,
  off_k,
  gt_mu = 12,
  gt_sd = 6,
  epidemic_size = 200L,
  duration = 365L,
  forest_size = 50L,
  n_ref = 50L
) {
  fp <- build_params(
    off_R = off_R,
    off_k = off_k,
    gt_mu = gt_mu,
    gt_sd = gt_sd,
    epidemic_size = epidemic_size,
    duration = duration,
    replicates = n_ref
  )
  ref_trees <- build_tree(fp)

  out <- imap(ref_trees, \(ref_tree, ref_id) {
    ref_off <- get_offspring(ref_tree, fp)

    # --- Reference offspring (all n_ref trees) ---
    ref_tbl <- tibble(
      offspring = ref_off,
      ref_id = ref_id,
      true_R = off_R,
      true_k = off_k
    )

    # --- QQ data (forest vs reference) ---
    forest <- build_forest(ref_tree, fp, forest_size = forest_size)
    forest_trees <- forest_to_trees(forest, ref_tree)

    probs <- ppoints(length(ref_off))
    ref_q <- quantile(ref_off, probs = probs)

    qq_tbl <- imap_dfr(forest_trees, \(tree, tree_id) {
      forest_off <- get_offspring(tree, fp)
      tibble(
        prob = probs,
        ref_q = ref_q,
        forest_q = quantile(forest_off, probs = probs),
        tree_id = tree_id,
        ref_id = ref_id,
        true_R = off_R,
        true_k = off_k
      )
    })

    list(qq = qq_tbl, ref = ref_tbl)
  })

  list(
    qq = bind_rows(map(out, "qq")),
    ref = bind_rows(map(out, "ref"))
  )
}


# ------------------------------------
#           Run pipeline
# ------------------------------------
param_grid <- expand_grid(
  off_R = c(1.5, 2, 3),
  off_k = c(0.1, 0.5, 3)
)

set.seed(123)
plan(multisession, workers = availableCores() - 2)
results <- furrr::future_pmap(
  param_grid,
  \(off_R, off_k) pipeline(off_R = off_R, off_k = off_k),
  .options = furrr_options(seed = TRUE)
) |>
  pipetime::time_pipe("Unified pipeline execution")
plan(sequential)

qq_data <- bind_rows(map(results, "qq"))
ref_offspring <- bind_rows(map(results, "ref"))


# =====================================================
#  Empirical vs Theoretical Offspring PMF
# =====================================================
# Point + errorbar = mean proportion ± 95% interval across reference trees.
# Red line + dots  = theoretical NegBin(R, k) PMF.

max_off <- 15L

# Per-tree proportions, ensuring 0-counts are represented
per_tree <- ref_offspring |>
  count(true_R, true_k, ref_id, offspring) |>
  complete(
    nesting(true_R, true_k, ref_id),
    offspring = 0:max_off,
    fill = list(n = 0L)
  ) |>
  group_by(true_R, true_k, ref_id) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

empirical <- per_tree |>
  group_by(true_R, true_k, offspring) |>
  summarise(
    mean_prop = mean(prop),
    lo = quantile(prop, 0.025),
    hi = quantile(prop, 0.975),
    .groups = "drop"
  )

theoretical <- param_grid |>
  rowwise() |>
  reframe(tibble(
    offspring = 0:max_off,
    pmf = dnbinom(0:max_off, size = off_k, mu = off_R),
    true_R = off_R,
    true_k = off_k
  ))

ggplot() +
  geom_pointrange(
    data = empirical,
    aes(x = offspring, y = mean_prop, ymin = lo, ymax = hi),
    size = 0.15,
    linewidth = 0.3,
    colour = "grey30"
  ) +
  geom_line(
    data = theoretical,
    aes(x = offspring, y = pmf),
    colour = "red",
    linewidth = 0.6
  ) +
  geom_point(
    data = theoretical,
    aes(x = offspring, y = pmf),
    colour = "red",
    size = 0.8
  ) +
  scale_x_continuous(limits = c(-1, max_off)) +
  labs(x = "Offspring count", y = "Proportion") +
  facet_Rk()

ggsave(
  "figures/S8.pdf",
  width = 7.5,
  height = 7.5,
  device = cairo_pdf,
  dpi = 300
)


# =====================================================
#  Uncensored sample size
# =====================================================
ref_offspring |>
  count(true_R, true_k, ref_id, name = "n_uncensored") |>
  ggplot(aes(x = n_uncensored)) +
  geom_histogram(
    fill = "grey80",
    colour = "grey40",
    binwidth = 5,
    boundary = 0
  ) +
  labs(x = "Number of uncensored cases", y = "Count") +
  facet_Rk(scales = "fixed")

ggsave(
  "figures/S9.pdf",
  width = 7.5,
  height = 7.5,
  device = cairo_pdf,
  dpi = 300
)

# =====================================================
#  Plot 1: QQ Plot (boxplots)
# =====================================================
# x = reference offspring count (rounded to integer)
# boxplot = distribution of forest offspring at the same quantile level
# dashed red line = y = x (perfect agreement)

# Bin by quantile probability (not reference value) to avoid
# regression-to-the-mean artifact from conditioning on one axis.
qq_data |>
  mutate(
    prob_bin = cut(prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
    midpoint = (ref_q + forest_q) / 2  # for optional reference
  ) |>
  ggplot(aes(x = ref_q, y = forest_q)) +
  geom_bin2d(binwidth = 1) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 20)) +
  labs(x = "Reference tree offspring", y = "Forest offspring") +
  facet_Rk()

ggsave(
  "figures/S10.pdf",
  width = 7.5,
  height = 7.5,
  device = cairo_pdf,
  dpi = 300
)

# =====================================================
#  Plot 2: Deviation Plot
# =====================================================
# (forest_q - ref_q) binned by quantile probability.
# Boxplot centred at 0 = no systematic bias.

qq_data |>
  mutate(
    deviation = forest_q - ref_q,
    prob_bin = cut(prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  ) |>
  ggplot(aes(x = prob_bin, y = deviation)) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "red") +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.3,
    fill = "grey80",
    linewidth = 0.3
  ) +
  labs(x = "Quantile probability", y = "Forest - Reference (offspring)") +
  facet_Rk(scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(
  "figures/S11.pdf",
  width = 7.5,
  height = 7.5,
  device = cairo_pdf,
  dpi = 300
)
