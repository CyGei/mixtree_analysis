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
  population = 100,
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
  off_R      = offspring_dist$parameters$R,
  off_k      = 1e5,
  gt_mu      = generation_time$parameters$mu,
  gt_sd      = generation_time$parameters$sd,
  population = population,
  duration   = duration,
  n_trees    = n_trees
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
  off_R      = offspring_dist$parameters$R,
  off_k      = offspring_dist$parameters$k + 0.1,
  gt_mu      = generation_time$parameters$mu,
  gt_sd      = generation_time$parameters$sd,
  population = population,
  duration   = duration,
  n_trees    = n_trees
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

param_grid <- expand.grid(
  off_R = 3,
  off_k = c(0.1, 0.3, 1, 10, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  population = 100,
  duration = 365,
  n_trees = 50,
  stringsAsFactors = FALSE
) |> as_tibble()

base_id <- param_grid |>
  mutate(row_number = row_number()) |>
  filter(off_R == 3, off_k == 1e5) |>
  pull(row_number)

base_tree <- param_grid |>
  slice(base_id) |>
  pmap(build_params) |>
  pluck(1) |> # same as .[[1]]
  build_tree()

base_tree |>
  slice(-1) |>
  igraph::graph_from_data_frame(directed = TRUE) |>
  ggraph::ggraph(layout = "tree") +
  ggraph::geom_edge_link() +
  ggraph::geom_node_point(col = "grey") +
  ggraph::geom_node_text(aes(label = name), col = "blue") # Label nodes by their 'name' attribute (ID)


plan(multisession, workers = availableCores() - 2)

system.time({
  forest_grid <- param_grid |>
    mutate(
      params = pmap(param_grid, build_params),
      forest = future_map(params, ~ build_forest(true_tree = base_tree, params = .),
        .options = furrr_options(seed = TRUE)
      )
    )
})
# user  system elapsed
# 4.06    0.57   21.19

system.time({
  test_grid <- forest_grid |>
    mutate(
      p_value = map_dbl(
        .x = forest,
        .f = ~ {
          sample_n <- 25
          suppressWarnings(
            mixtree::tree_test(
              sample(forest_grid$forest[[base_id]], sample_n),
              sample(.x, sample_n),
              method = "permanova",
              test_args = list(permutations = 500)
            )$`Pr(>F)` |> pluck(1)
          )
        }
      )
    )  
})
# user  system elapsed 
# 1.89    0.64    2.55 

base_forest <- forest_grid$forest[[base_id]]
system.time({
test_grid <- forest_grid |>
  mutate(
    p_value = future_map_dbl(
      .x = forest,
      .f = ~ {
        sample_n <- 25
        suppressWarnings(
          mixtree::tree_test(
            sample(base_forest, sample_n),
            sample(.x, sample_n),
            method = "permanova",
            test_args = list(permutations = 500)
          )$`Pr(>F)` |> pluck(1)
        )
      },
      .options = furrr_options(
        globals = c("base_forest"),
        seed = TRUE
      )
    )
  )
})
#  user  system elapsed 
# 3.14    0.42    7.56 

test_grid |>
  select(-c(params, forest)) |> 
  ggplot(aes(x = off_k, y = p_value)) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = 0.05), lty = "dashed")+
  scale_x_log10()+
  labs(
    x = "Overdispersion parameter (k)",
    y = "PERMANOVA p-value",
  ) +
  theme_minimal()
