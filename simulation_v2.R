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
forest_B <- build_forest(true_tree = tree_A, params = params_B )
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
forest_C <- build_forest(true_tree = tree_A, params = params_C )
check_forest(forest_C)

# different forests p<0.05
mixtree::tree_test(
  forest_A[sample(1:100, size = 50)],
  forest_C[sample(1:100, size = 50)],
  method = "permanova"
)
