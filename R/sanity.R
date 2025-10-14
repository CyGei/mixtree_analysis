# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log = "log", pipetime.unit = "min")

# =================================================
# Configuration
# =================================================
test_config <- list(
  method = c("permanova", "chisq"),
  forest_size = c(20L, 50L, 100L, 200L)
)

sim_config <- list(
  off_R = c(1.5, 2.5),
  off_k = c(0.1, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  epidemic_size = c(10L, 20L),
  duration = 365L,
  replicates = 10L
)

param_grid <- expand_grid(!!!sim_config) |>
  mutate(
    params = pmap(pick(everything()), build_params),
    param_id = row_number(),
    .before = 1
  )
param_grid
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

tree_grid

tree_grid |>
  filter(param_id == 8, tree_id == 1) |>
  pull(tree) |>
  igraph::graph_from_data_frame() |>
  plot()


# check that trees differ by checking their from column
# we can just sum the from column and see if values differ
tree_grid |>
  filter(param_id == 1) |>
  pull(tree) |>
  bind_rows(.id = "replicate") |>
  group_by(replicate) |>
  summarise(sum_from = sum(as.integer(from), na.rm = TRUE))

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

# check that forests differ by checking their from column
# we can just sum the from column and see if values differ
forest_grid |>
  filter(param_id == 1) |>
  pull(forest) |>
  map_dfr(~ as_tibble(.x, stringsAsFactors = FALSE), .id = "replicate") |>
  pivot_longer(
    cols = -replicate,
    names_to = "to",
    values_to = "from"
  ) |>
  group_by(replicate) |>
  summarise(sum_from = sum(as.integer(from), na.rm = TRUE))


# For faster loading during testing
forests <- forest_grid |>
  (\(df) {
    ids <- unique(df$param_id)
    df |>
      group_split(param_id) |>
      map(~ .x |> select(tree_id, forest) |> deframe()) |>
      set_names(ids)
  })()

identical(
  forests[[1]][[10]],
  forest_grid |>
    filter(param_id == 1, tree_id == 10) |>
    pull(forest) |>
    pluck(1)
)

# Check that all forests are identical to the ones in forest_grid
all_identical <- all(
  unlist(
    lapply(names(forests), function(param) {
      sapply(names(forests[[param]]), function(tree) {
        identical(
          forests[[param]][[tree]],
          forest_grid |>
            filter(
              param_id == as.integer(param),
              tree_id == as.integer(tree)
            ) |>
            pull(forest) |>
            pluck(1)
        )
      })
    })
  )
)
all_identical


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


# Check that test_grid compares only forests of the same epidemic_size
param_id_episize20 <- param_grid |>
  filter(epidemic_size == 20L) |>
  pull(param_id)
test_grid |>
  filter(param_id_A %in% param_id_episize20) |>
  pull(param_id_B) |>
  unique()
param_id_episize20

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

results |> as_tibble()

results_grid <- results |>
  left_join(param_grid |> select(-params), by = c("param_id_A" = "param_id")) |>
  left_join(
    param_grid |> select(-params),
    by = c("param_id_B" = "param_id"),
    suffix = c("_A", "_B")
  ) |>
  select(
    -c(
      starts_with("gt_"),
      starts_with("duration"),
      starts_with("replicates"),
      epidemic_size_B
    )
  ) |>
  rename(epidemic_size = epidemic_size_A) |>
  as_tibble()

results_grid
source("R/plots.R")
plot.ROC()

# check that mixtree_test behaves corrrectly:
as_forest(forests[[8]][[1]])
forests[[8]][[1]] |> nrow()
forests[[8]][[1]][1:11, ] |> as_forest()
