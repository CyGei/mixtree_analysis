# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
options(pipetime.log_file = "pipetime.log")
options(pipetime.time_unit = "hours")
# =================================================
# Configuration
# =================================================

# -------------------------------------------------
# Parameters for outbreak simulations
# -------------------------------------------------
sim_config <- list(
  off_R = c(1.5, 2, 2.5, 3, 5),
  off_k = c(0.1, 0.3, 1, 10, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  epidemic_size = c(20, 50, 100, 200),
  duration = 365,
  replicates = 10
)
# Number of unique parameter sets
n_params <- map_int(sim_config, length) |> prod()
cat("Number of unique parameter sets:", n_params, "\n")

# -------------------------------------------------
# Number of unique test conditions
# -------------------------------------------------
test_config <- list(
  method = c("permanova", "chisq"),
  sample_size = c(20, 50, 100, 200)
)
# Represents the number of different test conditions applied to each comparison pair
n_tests <- map_int(test_config, length) |> prod()
cat("Number of test conditions:", n_tests, "\n")

# -------------------------------------------------
# Number of comparison pairs
# -------------------------------------------------
# Total number of pairwise comparisons between forests,
# including self-comparisons but counting each pair only once
# (e.g. keep A vs B, ignore not B vs A)
n_pairs <- with(
  sim_config,
  rep(length(off_R) * length(off_k) * replicates, times = length(epidemic_size))
) |>
  map_dbl(~ .x * (.x + 1) / 2) |>
  sum()
cat("Number of forests to compare:", n_pairs, "\n")

# -------------------------------------------------
# Number of tests performed
# -------------------------------------------------
# Total number of tests across all pairs and test conditions
cat("Total number of tests to perform:", n_pairs * n_tests, "\n")

# =================================================
# Outbreak simulation
# =================================================
# Generate a reference transmission tree per outbreak setting

plan(multisession, workers = availableCores() - 1)

tree_grid <- sim_config |>
  expand.grid(stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(param_id = row_number()) |>
  mutate(
    params = pmap(pick(-param_id), build_params),
    tree = future_map(params, build_tree, .options = furrr_options(seed = TRUE))
  ) |>
  time_pipe("tree_grid")

saveRDS(tree_grid, "data/tree_grid.rds")

# =================================================
# Forest generation
# =================================================
# Build all forests (each with its own reference tree)
forest_grid <- tree_grid |>
  unnest(tree) |>
  group_by(param_id) |>
  mutate(tree_id = row_number()) |>
  ungroup() |>
  mutate(
    forest = future_map2(
      tree,
      params,
      ~ build_forest(tree = .x, params = .y, forest_size = 200),
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  time_pipe("forest_grid")

saveRDS(forest_grid, "data/forest_grid.rds")

# =================================================
# Test results
# =================================================

test_grid <- forest_grid |>
  select(param_id, tree_id, off_R, off_k, epidemic_size, forest) |>
  (\(df) {
    inner_join(
      df,
      df,
      suffix = c("_A", "_B"),
      by = "epidemic_size",
      relationship = "many-to-many"
    )
  })() |>
  # Keep pairs where params A <= B
  filter(param_id_A <= param_id_B) |>
  crossing(
    sample_size = test_config$sample_size,
    method = test_config$method
  ) |>
  mutate(
    p_value = future_pmap_dbl(
      .l = list(forest_A, forest_B, sample_size, method),
      .f = mixtree_test,
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  time_pipe("test_grid")

saveRDS(test_grid, "data/test_grid.rds")

# =================================================
# Plot Grid
# =================================================
source("R/plots.R")

# Must be run on Rstudio to render correctly!
p <- plot_test_grid(sample_size = 100, method = "permanova")
ggsave(p, filename = "figures/test_grid_permanova.png", width = 10, height = 10)

p <- plot_test_grid(sample_size = 20, half_tiles = TRUE)
ggsave(p, filename = "figures/test_grid.png", width = 10, height = 10)


# =================================================
# Plot ROC
# =================================================

plot.ROC()
