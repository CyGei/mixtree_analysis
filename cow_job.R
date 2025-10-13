library(mixtree)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

# =================================================
# mixtree test wrapper
# =================================================
# Convert matrix to list of trees
as_forest <- function(forest) {
  n_rows <- nrow(forest)
  col_names <- colnames(forest)

  lapply(seq_len(n_rows), function(i) {
    data.frame(
      from = forest[i, ],
      to = col_names,
      stringsAsFactors = FALSE
    )
  })
}

# Run one test given sampled matrices
.mixtree_test <- function(sample_A, sample_B, method) {
  forest_A <- as_forest(sample_A)
  forest_B <- as_forest(sample_B)

  p_value <- switch(
    method,
    "permanova" = {
      test <- mixtree::tree_test(
        forest_A,
        forest_B,
        method = method,
        test_args = list(permutations = 100)
      )
      test[["Pr(>F)"]][[1]]
    },
    "chisq" = {
      test <- mixtree::tree_test(
        forest_A,
        forest_B,
        method = method,
        test_args = list(simulate.p.value = TRUE, B = 100)
      )
      test[["p.value"]]
    },
    stop("Unknown method: ", method)
  )
  return(p_value)
}

# Run all tests for one pair of forests
mixtree_test <- function(
  forests, # a list such that forests[[param_id]][[tree_id]] is a matrix of ancestries
  param_id_A,
  tree_id_A,
  param_id_B,
  tree_id_B,
  forest_sizes, #vector of sizes to test
  methods #vector of methods to test
) {
  forest_A <- forests[[param_id_A]][[tree_id_A]]
  forest_B <- forests[[param_id_B]][[tree_id_B]]

  grid <- tidyr::expand_grid(
    forest_size = forest_sizes,
    method = methods
  )

  results <- purrr::pmap_dfr(grid, function(forest_size, method) {
    # Sample once per size
    sample_A <- forest_A[sample(nrow(forest_A), forest_size), , drop = FALSE]
    sample_B <- forest_B[sample(nrow(forest_B), forest_size), , drop = FALSE]
    p_value <- suppressWarnings(.mixtree_test(sample_A, sample_B, method))

    data.frame(
      param_id_A = param_id_A,
      tree_id_A = tree_id_A,
      param_id_B = param_id_B,
      tree_id_B = tree_id_B,
      forest_size = forest_size,
      method = method,
      p_value = p_value
    )
  })

  return(results)
}

# =================================================
# Load Data
# =================================================
test_grid <- readRDS("data/test_grid.rds")
forests <- readRDS("data/forests.rds")
# forests is a list such that forests[[param_id]][[tree_id]] is a matrix of ancestries

# =================================================
# hipercow job function
# =================================================
cow_job <- function(idx) {
  subgrid <- test_grid[idx, ]
  purrr::pmap_dfr(
    .l = subgrid,
    .f = function(param_id_A, tree_id_A, param_id_B, tree_id_B) {
      mixtree_test(
        forests = forests,
        param_id_A = param_id_A,
        tree_id_A = tree_id_A,
        param_id_B = param_id_B,
        tree_id_B = tree_id_B,
        forest_sizes = c(20L, 50L, 100L, 200L),
        methods = c("chisq", "permanova")
      )
    }
  )
}
