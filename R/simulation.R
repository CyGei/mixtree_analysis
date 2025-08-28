# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log_file = "pipetime.log", pipetime.time_unit = "hours")

# =================================================
# Configuration
# =================================================
test_config <- list(
  method = c("permanova", "chisq"),
  sample_size = c(20L, 50L, 100L, 200L)
)

sim_config <- list(
  off_R = c(1.5, 2, 2.5, 3, 5),
  off_k = c(0.1, 0.3, 1, 10, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  epidemic_size = c(20L, 50L, 100L, 200L),
  duration = 365L,
  replicates = 10L
)

param_grid <- expand_grid(!!!sim_config) |>
  mutate(
    params = pmap(pick(everything()), build_params),
    param_id = row_number(),
    .before = 1
  )
saveRDS(param_grid, "data/param_grid.rds")

# =================================================
# Tree simulation
# =================================================
# Generate a reference transmission tree per outbreak setting with n replicates
tree_grid <- param_grid |>
  mutate(
    tree = future_map(params, build_tree, .options = furrr_options(seed = TRUE))
  ) |>
  select(param_id, params, tree) |>
  mutate(tree_id = row_number(), .by = param_id, .after = param_id) |>
  time_pipe("tree_grid")

saveRDS(tree_grid, "data/tree_grid.rds")

# =================================================
# Forest generation
# =================================================
# Build all forests (each with its own reference tree)
forest_grid <- tree_grid |>
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

write_dataset(
  readRDS("data/forest_grid.rds"),
  path = "data/forests",
  format = "parquet",
  partitioning = c("param_id", "tree_id")
)

# =================================================
# Test results
# =================================================

# Perform all pairwise tests between forests of the same size
test_grid <- param_grid |>
  mutate(tree_id = map(replicates, seq_len)) |>
  select(epidemic_size, param_id, tree_id) |>
  unnest(tree_id) |>
  (\(x) {
    inner_join(
      x,
      x,
      by = "epidemic_size",
      suffix = c("_A", "_B"),
      relationship = "many-to-many"
    )
  })() |>
  filter(param_id_A <= param_id_B) |>
  crossing(expand_grid(!!!test_config)) |>
  mutate(
    p_value = future_pmap_dbl(
      .l = list(
        param_id_A,
        tree_id_A,
        param_id_B,
        tree_id_B,
        sample_size,
        method
      ),
      .f = mixtree_test,
      .options = furrr_options(
        seed = TRUE,
        packages = c("dplyr", "purrr", "arrow", "mixtree"),
        globals = "mixtree_test"
      )
    )
  ) |>
  time_pipe("test_grid")

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
