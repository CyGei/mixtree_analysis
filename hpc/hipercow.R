# https://mrcdata.dide.ic.ac.uk/hpc/managejobs.php
library(dplyr)
library(purrr)
library(tibble)
library(hipercow)
library(qs)

# ------------------------------------
#           Hipercow Setup
# ------------------------------------
hipercow_init()
setwd("P:/mixtree_analysis")
hipercow_configure("dide-windows")

# ------------------------------------
#         Forest Generation
# ------------------------------------
hipercow_environment_create(
    packages = c(
        "dplyr",
        "purrr",
        "mixtree",
        "tidyr",
        "tibble",
        "LaplacesDemon",
        "distcrete"
    ),
    sources = c("hpc/forest_job.R"),
    globals = c(
        "forest_job",
        "build_forest"
    )
)
hipercow_provision()
hipercow_configuration()

set.seed(000)
job_grid <- readRDS("data/tree_grid.rds") |>
    mutate(
        epidemic_size = map_dbl(tree, nrow)
    ) |>
    inner_join(
        readRDS("data/param_grid.rds") |>
            select(param_id, params, epidemic_size),
        by = "epidemic_size",
        relationship = "many-to-many"
    ) |>
    nrow() |>
    seq_len() |>
    sample() |>
    (\(x) split(x, ceiling(seq_along(x) / 50)))() |> # X tasks per job
    enframe(name = "job_id", value = "idx") |>
    mutate(job_id = as.integer(job_id))

forest_dir <- "data/forest_grid"
dir.create(forest_dir, showWarnings = FALSE)
bundle <- task_create_bulk_expr(
    expr = forest_job(idx),
    data = job_grid,
    environment = "default",
    bundle_name = "forest_generation"
)

hipercow_bundle_list()
hipercow_bundle_status(bundle) |> table()
hipercow_bundle_log_value(bundle)
results <- hipercow_bundle_result(bundle)

results |>
    bind_rows() |>
    mutate(
        forest = map(forest, function(m) {
            matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m))
        })
    ) |>
    qs::qsave("data/forest_grid.qs")


qs::qread("data/forest_grid.qs") |>
    select(forest_id) |>
    tidyr::separate_wider_delim(
        forest_id,
        delim = "_",
        names = c("tree_param_id", "tree_replicate", "param_id"),
        cols_remove = FALSE
    ) |>
    mutate(tree_id = paste0(tree_param_id, "_", tree_replicate)) |>
    select(tree_id, forest_id) |>
    (\(x) {
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
#          Mixtree Test
# ------------------------------------
hipercow_environment_create(
    packages = c("dplyr", "purrr", "mixtree", "tidyr", "tibble", "arrow"),
    sources = c("hpc/mixtree_job.R"),
    globals = c(
        "mixtree_test",
        ".mixtree_test",
        "mixtree_job",
        "as_forest",
        "read_forest",
        "test_grid"
    )
)
hipercow_provision()
hipercow_configuration()

set.seed(000)
job_grid <- readRDS("data/test_grid.rds") |>
    nrow() |>
    seq_len() |>
    sample() |>
    (\(x) split(x, ceiling(seq_along(x) / 225)))() |> # X tasks per job
    enframe(name = "job_id", value = "idx") |>
    mutate(job_id = as.integer(job_id))


results_dir <- "data/results"
dir.create(results_dir, showWarnings = FALSE)
bundle <- task_create_bulk_expr(
    expr = mixtree_job(idx),
    data = job_grid,
    environment = "default",
    bundle_name = "mixtree_test"
)
hipercow_bundle_list()
hipercow_bundle_status(bundle) |> table()
hipercow_bundle_log_value(bundle)

results <- hipercow_bundle_result(bundle)
saveRDS(bind_rows(results), "data/results.rds")

r <- readRDS("data/results.rds")
r |> as_tibble() |> glimpse()
