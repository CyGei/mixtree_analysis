# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log = "log", pipetime.unit = "hours")

# =================================================
# Configuration
# =================================================
test_config <- list(
  method = c("permanova", "chisq"),
  forest_size = c(20L, 50L, 100L, 200L)
)

sim_config <- list(
  off_R = c(1.5, 2.5, 4.5),
  off_k = c(0.1, 0.5, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  epidemic_size = c(20L, 50L, 100L, 200L),
  duration = 365L,
  replicates = 100L
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
  unnest(tree) |>
  mutate(tree_id = row_number(), .by = param_id, .after = param_id) |>
  time_pipe("tree_grid")

saveRDS(tree_grid, "data/tree_grid.rds")
tree_grid <- readRDS("data/tree_grid.rds")
# =================================================
# Forest generation
# =================================================
# Build all forests (each with its own reference tree)
forest_grid <-
  future_pmap(
    tree_grid,
    build_forest,
    .options = furrr_options(
      seed = TRUE,
      packages = c("dplyr", "LaplacesDemon", "mixtree"),
      globals = c("build_forest", ".build_forest")
    )
  ) |>
  bind_rows() |>
  time_pipe("forest_grid")
saveRDS(forest_grid, "data/forest_grid.rds")
forest_grid <- readRDS("data/forest_grid.rds")

# For faster loading during testing
forests <- readRDS("data/forest_grid.rds") |>
  (\(df) {
    ids <- unique(df$param_id)
    df |>
      group_split(param_id) |>
      map(~ .x |> select(tree_id, forest) |> deframe()) |>
      set_names(ids)
  })()

# =================================================
# Test results
# =================================================
# Perform all pairwise tests between forests of the same size
test_grid <- param_grid |>
  mutate(tree_id = map(replicates, seq_len)) |>
  select(epidemic_size, param_id, tree_id) |>
  unnest(tree_id) |>
  group_by(epidemic_size) |> # only combine forests of the same epidemic_size
  group_modify(~ cross_join(.x, .x, suffix = c("_A", "_B"))) |>
  ungroup() |>
  select(-epidemic_size) |>
  filter(param_id_A <= param_id_B)
saveRDS(test_grid, "data/test_grid.rds")
test_grid <- readRDS("data/test_grid.rds")


results <- furrr::future_pmap_dfr(
  test_grid,
  mixtree_test,
  forests = forests,
  forest_sizes = test_config$forest_size,
  methods = test_config$method,
  .options = furrr::furrr_options(
    seed = TRUE,
    packages = c("dplyr", "purrr", "mixtree", "tidyr"),
    globals = c("mixtree_test", ".mixtree_test", "as_forest")
  )
) |>
  time_pipe("results")
