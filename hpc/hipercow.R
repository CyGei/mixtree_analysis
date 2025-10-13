# https://mrcdata.dide.ic.ac.uk/hpc/managejobs.php

# =============================================================
# Hipercow setup
# =============================================================
library(hipercow)
hipercow_init()
setwd("P:/mixtree_analysis")
hipercow_configure("dide-windows")
hipercow_configuration()

hipercow_environment_create(
    packages = c("dplyr", "purrr", "mixtree", "tidyr", "tibble"),
    sources = c("cow_job.R"),
    globals = c("mixtree_test", ".mixtree_test", "as_forest", "forests")
)
hipercow_configuration()


# =================================================
# Submit jobs
# =================================================
# Each row of test_grid is a single task to perform
# We split into random chunks of 1000 tasks per job using job_grid

test_grid <- readRDS("data/test_grid.rds")

set.seed(123)
idx <- sample(seq_len(nrow(test_grid)))
chunk_size <- 1000 # number of tasks per job
chunks <- split(idx, ceiling(seq_along(idx) / chunk_size))
job_grid <- tibble::tibble(
    job_id = seq_along(chunks),
    idx = chunks
)
saveRDS(job_grid, "data/job_grid.rds")

bundle <- task_create_bulk_expr(
    expr = cow_job(idx),
    data = job_grid,
    environment = "default",
    bundle_name = "mixtree_test"
)
hipercow_bundle_list()
hipercow_bundle_status(bundle)
hipercow_bundle_log_value(bundle)
results_raw <- hipercow_bundle_result(bundle) |> bind_rows()
results_raw |> object.size() |> format('Gb')
saveRDS(results_raw, "data/results_raw.rds")
results_raw <- readRDS("data/results_raw.rds")
# =================================================
# Analyze results
# =================================================
param_grid <- readRDS("data/param_grid.rds") |> select(-params)
results_grid <- results_raw |>
    left_join(param_grid, by = c("param_id_A" = "param_id")) |>
    left_join(
        param_grid,
        by = c("param_id_B" = "param_id"),
        suffix = c("_A", "_B")
    )

# Sanity checks:
# 1.  epidemic_size_A == epidemic_size_B
# 2. duration_A == duration_B
# 3. replicates_A == replicates_B
if (
    with(
        results_grid,
        all(
            epidemic_size_A == epidemic_size_B &
                duration_A == duration_B &
                replicates_A == replicates_B
        )
    )
) {
    results_grid <- results_grid |>
        mutate(
            epidemic_size = epidemic_size_A,
            duration = duration_A,
            replicates = replicates_A
        ) |>
        select(
            -epidemic_size_A,
            -epidemic_size_B,
            -duration_A,
            -duration_B,
            -replicates_A,
            -replicates_B
        ) |>
        arrange(param_id_A, param_id_B, tree_id_A, tree_id_B, epidemic_size)
}

saveRDS(results_grid, "data/results_grid.rds")
