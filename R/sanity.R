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
  forest_size = c(20L, 100L)
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
  plot(layout = layout_as_tree)


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

# ------------------------------------
#           Check trees
# ------------------------------------
library(fitdistrplus)
library(ggh4x)
# results suggests that the test is sensitive to the tree_id:
# that is even if param_id_A == param_id_B, the test rejects H0
# when tree_id_A != tree_id_B

trees <- readRDS("data/tree_grid.rds") |>
  mutate(
    R = map_dbl(params, ~ .x$offspring_dist$parameters$R),
    k = map_dbl(params, ~ .x$offspring_dist$parameters$k),
    R_k = paste0("R = ", R, ", k = ", k),
    epidemic_size = map_dbl(params, ~ .x$epidemic_size)
  ) |>
  dplyr::select(-params)
group_vars <- c("epidemic_size", "param_id", "tree_id", "R", "k", "R_k")

Ri <- trees |>
  unnest(tree) |>
  filter(!is.na(from)) |>
  group_by(across(all_of(group_vars)), from) |>
  summarise(
    Ri = n(),
    .groups = "drop"
  )


# Per-tree fits (skip trees with only one observation)
per_tree_est <- Ri |>
  group_by(across(all_of(group_vars))) |>
  filter(n() > 1) |>
  group_modify(.f = function(d, ...) {
    fit <- fitdist(d$Ri, "nbinom")
    est <- fit$estimate
    sd <- fit$sd
    tibble(
      R_est = est["mu"],
      k_est = est["size"],
      k_se = sd["size"],
      R_se = sd["mu"]
    )
  }) |>
  ungroup() |>
  mutate(
    k_error = k_est - k,
    log_k_error = log(k_est) - log(k),
    R_error = R_est - R,
    R_z = (R_est - R) / R_se,
    k_z = (k_est - k) / k_se,
    R_CI_ok = (R >= (R_est - 1.96 * R_se)) & (R <= (R_est + 1.96 * R_se)),
    k_CI_ok = (k >= (k_est - 1.96 * k_se)) & (k <= (k_est + 1.96 * k_se))
  )

# Plot per-tree estimates vs true parameters
params <- per_tree_est |>
  dplyr::select(epidemic_size, R, k, R_k) |>
  distinct()

per_tree_est |>
  mutate(
    col_label = "Epidemic Size",
  ) |>
  ggplot() +
  ggh4x::facet_nested(
    rows = vars(R, k),
    cols = vars(col_label, epidemic_size),
    labeller = labeller(
      R = function(x) paste0("R = ", x),
      k = function(x) paste0("k = ", x)
    )
  ) +
  geom_point(aes(x = R_error, y = log_k_error)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  coord_cartesian() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )


# ------------------------------------
#   Compare per-tree vs pooled fits
# ------------------------------------
# Pooled fit
pooled_fit <- fitdist(Ri$Ri, "nbinom")
# AIC per tree
per_tree_fit <- Ri |>
  group_by(across(all_of(group_vars))) |>
  filter(n() > 1) |>
  group_modify(
    ~ {
      fit <- fitdist(.x$Ri, "nbinom")
      tibble(aic = fit$aic)
    }
  ) |>
  ungroup()
per_tree_fit
mean(per_tree_fit$aic)
sum(per_tree_fit$aic)
pooled_fit$aic
