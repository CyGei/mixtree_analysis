# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)

# =================================================
# Configuration
# =================================================
# Parameters for outbreak simulations
config <- list(
  off_R = c(1.5, 2, 2.5, 3, 5),
  off_k = c(0.1, 0.3, 1, 10, 1e5),
  gt_mu = 3,
  gt_sd = 1.5,
  epidemic_size = c(20, 50, 100, 200),
  duration = 365
)

# =================================================
# Outbreak simulation
# =================================================
# Generate a reference transmission tree per outbreak setting

plan(multisession, workers = availableCores() - 2)

tree_grid <- config |>
  expand.grid(stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(id = row_number()) |>
  mutate(
    params = pmap(pick(-id), build_params),
    tree = future_map(params, build_tree,
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  log_time("tree_grid")

saveRDS(tree_grid, "data/tree_grid.rds")

# =================================================
# Forest generation
# =================================================
# Build all forests (each with its own reference tree)
forest_grid <- tree_grid |>
  mutate(
    forest = future_map2(tree, params,
      ~ build_forest(tree = .x, params = .y, forest_size = 200),
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  log_time("forest_grid")

saveRDS(forest_grid, "data/forest_grid.rds")

# =================================================
# Test results
# =================================================
forest_grid <- readRDS("data/forest_grid.rds")

test_grid <- forest_grid |>
  select(id, off_R, off_k, epidemic_size) |>
  (\(df) inner_join(df, df, suffix = c("_A", "_B"), by = "epidemic_size", relationship = "many-to-many"))() |>
  # Keep pairs where id_A <= id_B to avoid duplicates/reversals
  filter(id_A <= id_B) |>
  crossing(sample_size = c(20, 50, 100, 200)) |>
  mutate(
    p_value = future_pmap_dbl(
      .l = list(id_A, id_B, sample_size),
      .f = function(id_A, id_B, sample_size) {
        suppressWarnings(
          mixtree::tree_test(
            sample(forest_grid$forest[[id_A]], sample_size),
            sample(forest_grid$forest[[id_B]], sample_size),
            method = "permanova",
            test_args = list(permuations = 200)
          )$`Pr(>F)` |> pluck(1)
        )
      },
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  log_time("mixtree")

saveRDS(test_grid, "data/test_grid.rds")
test_grid <- readRDS("data/test_grid.rds")

# =================================================
# Plot Grid
# =================================================
library(ggh4x)

test_grid |>
  filter(sample_size == 100) |>
  mutate(significance = if_else(p_value < 0.05, "Significant", "Non-significant")) |>
  ggplot(aes(x = factor(off_k_A), y = factor(off_k_B), fill = significance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_segment(
    data = test_grid |> filter(off_k_A == off_k_B, off_R_A == off_R_B),
    aes(
      x = as.numeric(factor(off_k_A)) - 0.5,
      y = as.numeric(factor(off_k_B)) - 0.5,
      xend = as.numeric(factor(off_k_A)) + 0.5,
      yend = as.numeric(factor(off_k_B)) + 0.5
    ),
    color = "black",
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +
  facet_nested(
    rows = vars(factor(off_R_A, levels = sort(unique(test_grid$off_R_A), decreasing = TRUE))),
    cols = vars(factor(epidemic_size), factor(off_R_B)),
    nest_line = element_line(linetype = 1),
    labeller = labeller(
      .rows = function(x) {
        vals <- sort(unique(test_grid$off_R_A), decreasing = TRUE)
        ifelse(x == vals[1], paste0("R0 = ", x), as.character(x))
      },
      .cols = function(x) {
        vals <- sort(unique(test_grid$epidemic_size))
        ifelse(x == vals[1], paste0("Epidemic size = ", x), as.character(x))
      }
    )
  ) +
  scale_fill_manual(values = c(
    "Significant" = "#f99820",
    "Non-significant" = "#03658c"
  )) +
  labs(
    x = "k (Forest A)",
    y = "k (Forest B)",
    fill = "Difference"
  ) +
  theme_light(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    panel.spacing = unit(1, "mm"), # reduce spacing between panels
    strip.switch.pad.grid = unit(0, "pt"), # remove extra padding
    strip.switch.pad.wrap = unit(0, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(color = "black"),
    axis.text.x = element_blank()
  )



# =================================================
# Plot ROC
# =================================================
# Define multiple alpha thresholds
alpha_values <- seq(0, 1, 0.01)

# Compute performance metrics for each alpha
perf_grid_alpha <- expand_grid(test_grid, alpha = alpha_values) %>%
  mutate(
    true_diff = (off_R_A != off_R_B) | (off_k_A != off_k_B),
    test_reject = p_value < alpha
  ) %>%
  group_by(epidemic_size, sample_size, alpha) %>%
  summarise(
    TP = sum(true_diff & test_reject),
    FN = sum(true_diff & !test_reject),
    TN = sum(!true_diff & !test_reject),
    FP = sum(!true_diff & test_reject),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    .groups = "drop"
  ) %>%
  mutate(fpr = 1 - specificity)

ggplot(perf_grid_alpha, aes(x = fpr, y = sensitivity, group = interaction(epidemic_size, sample_size))) +
  geom_line(aes(color = factor(sample_size))) +
  facet_wrap(~epidemic_size) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Sample Size"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.just = "left",
    legend.spacing = unit(0.1, "cm"),
    strip.background = element_rect(
      fill = "#f8f4f2",
      colour = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
