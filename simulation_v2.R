source("R/packages.R")
source("R/helpers.R")
set.seed(000)

# =================================================
# Worflow example
# =================================================

# Build the set of epidemic parameters
params_A <- build_params(
  off_R = 3,
  off_k = 0.3,
  gt_mu = 3,
  gt_sd = 1.5,
  n_cases = 100,
  duration = 365,
  n_trees = 100
)
params_A

# Build the true tree
tree_A <- build_tree(params = params_A)

# Build forests
forest_A <- build_forest(true_tree = tree_A, params = params_A)
check_forest(forest_A)

# Null test: (same forest) p>0.05
mixtree::tree_test(
  forest_A[1:50],
  forest_A[51:100],
  method = "permanova"
)

params_B <- with(params_A, build_params(
  off_R = offspring_dist$parameters$R,
  off_k = 1e5,
  gt_mu = generation_time$parameters$mu,
  gt_sd = generation_time$parameters$sd,
  n_cases = n_cases,
  duration = duration,
  n_trees = n_trees
))


# tree_A with params B
forest_B <- build_forest(true_tree = tree_A, params = params_B)
check_forest(forest_B)

# different forests p<0.05
mixtree::tree_test(
  forest_A[sample(1:100, size = 50)],
  forest_B[sample(1:100, size = 50)],
  method = "permanova"
)


params_C <- with(params_A, build_params(
  off_R = offspring_dist$parameters$R,
  off_k = offspring_dist$parameters$k + 0.1,
  gt_mu = generation_time$parameters$mu,
  gt_sd = generation_time$parameters$sd,
  n_cases = n_cases,
  duration = duration,
  n_trees = n_trees
))
# tree_A with params C
forest_C <- build_forest(true_tree = tree_A, params = params_C)
check_forest(forest_C)

# different forests p<0.05
mixtree::tree_test(
  forest_A[sample(1:100, size = 50)],
  forest_C[sample(1:100, size = 50)],
  method = "permanova"
)

# =================================================
# Simulation
# =================================================
plan(multisession, workers = availableCores() - 2)

# Parameter grid
param_grid <- expand.grid(
  off_R = c(1.5, 3, 5),
  off_k = c(0.1, 0.3, 1, 10, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  n_cases = 100,
  duration = 365,
  n_trees = 50,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  mutate(p_id = row_number())


# Build all forests (each with its own reference tree)
system.time({
  forest_grid <- param_grid |>
    mutate(
      # Exclude p_id as it's not an input of build_params
      params = pmap(select(param_grid, -p_id), build_params),
      # Each forest gets its own reference tree
      reference_tree = future_map(params, ~ build_tree(.x),
        .options = furrr_options(seed = TRUE)
      ),
      # Each forest built from its own reference tree
      forest = future_map2(reference_tree, params,
        ~ build_forest(true_tree = .x, params = .y),
        .options = furrr_options(seed = TRUE)
      )
    )
})

# Time when params$population_size (would result in trees with varying n_cases)
# user  system elapsed
# 36.50    5.11   69.74

# Time when params$n_cases (trees with all exactly n_cases)
# user  system elapsed
# 34.39    7.83  355.94
saveRDS(forest_grid, "data/forest_grid.rds")

# Pairwise comparison
comparison_grid <- expand.grid(
  p_A = forest_grid$p_id,
  p_B = forest_grid$p_id,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  # get unique comparisons (e.g. 1 vs 2 is same as 2 vs 1)
  # but get all other combinations
  #filter(p_A <= p_B) |>
  mutate(
    R_A = forest_grid$off_R[p_A],
    k_A = forest_grid$off_k[p_A],
    R_B = forest_grid$off_R[p_B],
    k_B = forest_grid$off_k[p_B]
  )

# all pairwise tests
system.time({
  results_grid <- comparison_grid |>
    mutate(
      p_value = future_map2_dbl(
        .x = p_A,
        .y = p_B,
        .f = ~ {
          sample_n <- 25
          suppressWarnings(
            mixtree::tree_test(
              sample(forest_grid$forest[[.x]], sample_n),
              sample(forest_grid$forest[[.y]], sample_n),
              method = "permanova",
              test_args = list(permutations = 500)
            )$`Pr(>F)` |> pluck(1)
          )
        },
        .options = furrr_options(seed = TRUE)
      )
    )
})
results_grid

results_grid |>
  filter(k_B >= k_A) |>
  mutate(significance = if_else(p_value < 0.05, "Significant", "Non-significant")) |>
  ggplot(aes(x = factor(k_A), y = factor(k_B), fill = significance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Draw a line using geom_segment() on the diagonal points
  geom_segment(
    data = results_grid |> filter(k_A == k_B),
    aes(
      x = as.numeric(factor(k_A)) - 0.5,
      y = as.numeric(factor(k_B)) - 0.5,
      xend = as.numeric(factor(k_A)) + 0.5,
      yend = as.numeric(factor(k_B)) + 0.5
    ),
    color = "black",
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +
  facet_grid(
    cols = vars(factor(R_B)),
    rows = vars(factor(R_A, levels = sort(unique(results_grid$R_A), decreasing = TRUE))),
    labeller = labeller(
      .rows = ~ paste("R =", .x),
      .cols = ~ paste("R =", .x)
    )
  ) +
  scale_fill_manual(values = c("Significant" = "#d73027", "Non-significant" = "#4575b4")) +
  labs(
    x = "k (Forest A)",
    y = "k (Forest B)",
    fill = "Difference"
  )
