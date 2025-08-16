# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)

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
  duration = 365
)
# Number of unique outbreak simulations
# Represents all distinct combinations of the varying outbreak parameters
n_sims <- map_int(sim_config, length) |> prod()
cat("Number of outbreak simulations:", n_sims, "\n")

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
# Total number of pairwise comparisons between outbreak simulations,
# including self-comparisons but counting each pair only once
# (e.g. keep A vs B, ignore not B vs A)
n_pairs <- with(
  sim_config,
  rep(length(off_R) * length(off_k), times = length(epidemic_size))
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

plan(multisession, workers = availableCores() - 2)

tree_grid <- sim_config |>
  expand.grid(stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(id = row_number()) |>
  mutate(
    params = pmap(pick(-id), build_params),
    tree = future_map(params, build_tree,
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  time_pipe("tree_grid", log_file = "time.log")

saveRDS(tree_grid, "data/tree_grid.rds")

# =================================================
# Forest generation
# =================================================
# Build all forests (each with its own reference tree)
forest_grid <- tree_grid |>
  mutate(
    forest = future_map2(tree, params,
      ~ build_forest(tree = .x, params = .y, forest_size = 200),
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  time_pipe("forest_grid", log_file = "time.log")

saveRDS(forest_grid, "data/forest_grid.rds")

# =================================================
# Test results
# =================================================
forest_grid <- readRDS("data/forest_grid.rds")

test_grid <- forest_grid |>
  select(id, off_R, off_k, epidemic_size) |>
  (\(df) inner_join(df, df, suffix = c("_A", "_B"), by = "epidemic_size", relationship = "many-to-many"))() |>
  # Keep pairs where id_A <= id_B to avoid duplicates/reversals
  filter(id_A <= id_B) |>
  crossing(sample_size = test_config$sample_size, method = test_config$method) |>
  mutate(
    p_value = future_pmap_dbl(
      .l = list(id_A, id_B, sample_size, method),
      .f = mixtree_test,
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  time_pipe("test_grid", log_file = "time.log")

saveRDS(test_grid, "data/test_grid.rds")
test_grid <- readRDS("data/test_grid.rds")

# =================================================
# Plot Grid
# =================================================
source("R/plot.test_grid.R")
plot.test_grid()


# =================================================
# Plot ROC
# =================================================
plot.ROC()
