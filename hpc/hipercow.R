# https://mrcdata.dide.ic.ac.uk/hpc/managejobs.php

# ------------------------------------
#           Hipercow Setup
# ------------------------------------
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


# ------------------------------------
#           Create Jobs
# ------------------------------------
set.seed(123)
job_grid <- readRDS("data/test_grid.rds") |>
    nrow() |>
    seq_len() |>
    sample() |>
    (\(x) split(x, ceiling(seq_along(x) / 100)))() |> # 100 tasks per job
    enframe(name = "job_id", value = "idx") |>
    mutate(job_id = as.integer(job_id))

bundle <- task_create_bulk_expr(
    expr = cow_job(idx),
    data = job_grid,
    environment = "default",
    bundle_name = "mixtree_test"
)
hipercow_bundle_list()
hipercow_bundle_status(bundle)
hipercow_bundle_log_value(bundle)

results <- hipercow_bundle_result(bundle)
