# ------------------------------------
#           Setup
# ------------------------------------
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log = "log", pipetime.unit = "mins")

# ------------------------------------
#           Parameter Configuration
# ------------------------------------
#' @title param_grid
#' @description
#' Build a grid of simulation parameters.
#' `param_grid` is a tibble containing all combinations
#' of *simulation parameters* defined in `sim_config` and
#' *test parameters* defined in `test_config`.
#' Each row contains a unique `param_id`.

test_config <- list(
  method = c("permanova", "chisq"),
  forest_size = c(20L, 50L, 100L, 200L)
)

sim_config <- list(
  off_R = c(1.5, 2, 3),
  off_k = c(0.1, 0.3, 0.5, 1, 1e5),
  gt_mu = 12,
  gt_sd = 6,
  epidemic_size = c(20L, 50L, 100L, 200L),
  duration = 365L,
  replicates = 10L
)

expand_grid(!!!sim_config) |>
  mutate(
    params = pmap(pick(everything()), build_params),
    param_id = row_number(),
    .before = 1
  ) |>
  saveRDS("data/param_grid.rds")

# ------------------------------------
#           Tree Generation
# ------------------------------------
#' @title tree_grid
#' @description
#' Build a transmission tree given simulation parameters.
#' `tree_grid` is a tibble containing a 'reference' transmission tree for each
#' parameter set (`params`) in `param_grid`.
#' For a given `params`, n `replicates` trees are generated to account for *stochasticity*.
#' `tree_id` is defined as `<param_id>_<replicate_id>`.

readRDS("data/param_grid.rds") |>
  mutate(
    tree = future_map(params, build_tree, .options = furrr_options(seed = TRUE))
  ) |>
  unnest(tree) |>
  mutate(
    tree_id = paste0(param_id, '_', row_number()),
    .by = param_id,
    .before = 1
  ) |>
  select(tree_id, tree) |>
  time_pipe("tree_grid") |>
  saveRDS("data/tree_grid.rds")

# ------------------------------------
#           Forest Generation
# ------------------------------------
#' @title forest_grid
#' @description
#' Build a forest for each reference tree in `tree_grid` and each parameter set
#' in `param_grid` with matching `epidemic_size`.
#' `forest_grid` is a tibble containing all generated forests.
#' `forest_id` is defined as `<tree_id>_<param_id>`.

readRDS("data/tree_grid.rds") |>
  mutate(
    epidemic_size = map_dbl(tree, nrow)
  ) |>
  inner_join(
    readRDS("data/param_grid.rds") |> select(param_id, params, epidemic_size),
    by = "epidemic_size",
    relationship = "many-to-many"
  ) |>
  mutate(
    forest = future_pmap(
      .l = list(tree, params),
      .f = build_forest,
      .options = furrr::furrr_options(
        seed = TRUE,
        packages = c("dplyr", "LaplacesDemon", "mixtree"),
        globals = c("build_forest")
      )
    )
  ) |>
  mutate(
    forest_id = paste0(tree_id, '_', param_id)
  ) |>
  select(forest_id, forest) |>
  arrange(
    as.integer(sub("_.*", "", forest_id)), # first number
    as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", forest_id)), # second number
    as.integer(sub(".*_", "", forest_id)) # last number
  ) |>
  time_pipe("forest_grid") |>
  saveRDS("data/forest_grid.rds")

# ------------------------------------
#           Test Grid
# ------------------------------------
#' @title test_grid
#' @description
#' Build all unique pairs of forests from the *same reference tree* for testing.
#' `test_grid` is a tibble containing all unique forest pairs to be tested.
#' Each row contains `forest_id_A` and `forest_id_B`.
#' Only pairs where `forest_id_A` <= `forest_id_B` are included to avoid
#' duplicate tests.
readRDS("data/forest_grid.rds") |>
  select(forest_id) |>
  (\(x) {
    x <- x |>
      mutate(tree_id = sub("^(\\d+_\\d+)_.*$", "\\1", forest_id))
    inner_join(
      x,
      x,
      by = "tree_id",
      suffix = c("_A", "_B"),
      relationship = "many-to-many"
    )
  })() |>
  filter(forest_id_A <= forest_id_B) |>
  select(-tree_id) |>
  saveRDS("data/test_grid.rds")

# ------------------------------------
#           Test Results
# ------------------------------------
#' @title results
#' @description
#' Perform the specified tests on each forest pair in `test_grid`.
#' `results` is a tibble containing the results of each test performed
#' with columns: `forest_id_A`, `forest_id_B`, `method`, `forest_size`, `p_value`.
furrr::future_pmap_dfr(
  .l = readRDS("data/test_grid.rds"),
  .f = mixtree_test,
  forests = deframe(readRDS("data/forest_grid.rds")),
  forest_sizes = test_config$forest_size,
  methods = test_config$method,
  .options = furrr::furrr_options(
    seed = TRUE,
    packages = c("dplyr", "purrr", "mixtree", "tidyr"),
    globals = c("mixtree_test", ".mixtree_test", "as_forest")
  )
) |>
  time_pipe("results") |>
  saveRDS("data/results.rds")


# ------------------------------------
#           results_grid
# ------------------------------------
#' @title results_grid
#' @description
#' A tibble combining `results` with corresponding parameters from `param_grid`.
#' Each row represents a statistical test comparing two forests.
#'
#' Columns:
#' - `forest_id_*`: identifiers of the compared forests
#' - `tree_id_*`: reference tree identifiers used to generate the forests
#' - `param_id_*`: parameter set identifiers used to generate the forests
#' - `off_R_*`: R0 of the offspring distribution
#' - `off_k_*`: dispersion parameter (k) of the offspring distribution
#' - `forest_size`: number of transmission trees in each forest
#' - `epidemic_size`: number of vertices per tree
#' - `method`: statistical test used
#' - `p_value`: resulting p-value

readRDS("data/results.rds") |>
  bind_rows() |>
  as_tibble() |>
  mutate(
    tree_id = str_extract(forest_id_A, "^\\d+_\\d+"),
    param_id_A = str_extract(forest_id_A, "(?<=_)\\d+$"),
    param_id_B = str_extract(forest_id_B, "(?<=_)\\d+$"),
    .after = starts_with("forest_id")
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(param_id, epidemic_size, starts_with("off_")) |>
      rename_with(~ paste0(.x, "_A"), starts_with("off_")) |>
      mutate(param_id = as.character(param_id)),
    by = c("param_id_A" = "param_id")
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(param_id, starts_with("off_")) |>
      rename_with(~ paste0(.x, "_B"), starts_with("off_")) |>
      mutate(param_id = as.character(param_id)),
    by = c("param_id_B" = "param_id")
  ) |>
  select(
    tree_id,
    starts_with("forest_id"),
    starts_with("param_id"),
    starts_with("off_"),
    ends_with("_size"),
    method,
    p_value
  ) |>
  saveRDS("data/results_grid.rds")
