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
  col_names <- colnames(forest)
  lapply(seq_len(nrow(forest)), function(i) {
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
        test_args = list(permutations = 999)
      )
      test[["Pr(>F)"]][[1]]
    },
    "chisq" = {
      test <- mixtree::tree_test(
        forest_A,
        forest_B,
        method = method,
        test_args = list(simulate.p.value = TRUE, B = 999)
      )
      test[["p.value"]]
    },
    stop("Unknown method: ", method)
  )
  return(p_value)
}

# Run all tests for one pair of forests
mixtree_test <- function(
  forest_id_A,
  forest_id_B,
  forests, # named list: forests[["forest_id"]] = matrix
  forest_sizes, # vector of integers
  methods # vector of characters ("chisq", "permanova")
) {
  forest_A <- forests[[forest_id_A]]
  forest_B <- forests[[forest_id_B]]

  grid <- expand_grid(
    forest_size = forest_sizes,
    method = methods
  )

  results <- purrr::pmap_dfr(grid, function(forest_size, method) {
    sample_A <- forest_A[sample(nrow(forest_A), forest_size), , drop = FALSE]
    sample_B <- forest_B[sample(nrow(forest_B), forest_size), , drop = FALSE]
    p_value <- suppressWarnings(.mixtree_test(sample_A, sample_B, method))

    data.frame(
      forest_id_A = forest_id_A,
      forest_id_B = forest_id_B,
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
forests <- deframe(readRDS("data/forest_grid.rds"))

# =================================================
# hipercow job function
# =================================================
cow_job <- function(idx) {
  subgrid <- test_grid[idx, ]
  purrr::pmap_dfr(
    .l = subgrid,
    .f = function(forest_id_A, forest_id_B) {
      mixtree_test(
        forest_id_A = forest_id_A,
        forest_id_B = forest_id_B,
        forests = forests,
        forest_sizes = c(20L, 50L, 100L, 200L),
        methods = c("permanova", "chisq")
      )
    }
  )
}
