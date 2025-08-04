library(tidyverse)
library(epitrix)
library(distcrete)
library(simulacr)
library(ggraph)
library(igraph)
source("R/helpers.R")
set.seed(000)

# parameters
off_R <- 3
off_k <- 0.3
gt_mu <- 3
gt_sd <- 1.5


params <- set_epi_params(off_R, off_k, gt_mu, gt_sd)

true_tree <- build_tree(off_R, off_k, gt_mu, gt_sd)

forest <- build_forest(true_tree, off_R, off_k, gt_mu, gt_sd, n_trees = 100)

# tree test
mixtree::tree_test(
  forest[1:50],
  forest[51:100],
  method = "permanova"
)

# Comparison
forest_pois <- build_forest(true_tree, off_R, off_k = 1e5, gt_mu, gt_sd, n_trees = 100)
mixtree::tree_test(
  forest_pois[sample(1:100, size = 50)],
  forest[sample(1:100, size = 50)],
  method = "permanova"
)

##############################
# Simulation
##############################
library(furrr)
plan(multisession, workers = 7)

k_grid <- c(0.1, 0.2, 0.3, 0.5, 1, 10, 1000)

forests_by_k <- setNames(k_grid, paste0("k_", k_grid)) |>
  furrr::future_map(function(k) build_forest(true_tree, off_R, k, gt_mu, gt_sd, n_trees = 50))

comparison_grid <- combn(names(forests_by_k), 2, simplify = FALSE)
