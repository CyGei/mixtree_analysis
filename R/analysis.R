# This scripts analyses the results of the simulation study
# `results_grid.rds`.
source("R/packages.R")
source("R/plots2.R")

# ------------------------------------
#           rejection_summary
# ------------------------------------
#' @title rejection_summary
#' @description
#' A tibble summarising hypothesis test results from the simulation study `results_grid.rds`.
#' It reports the proportion of null hypothesis (H0) rejections for each combination of parameters and test method.
#'
#' Columns:
#' - `hypothesis`: whether H0 (null) or H1 (alternative) is true
#' - `method`: statistical test used ("chisq" or "permanova")
#' - `forest_size`: number of transmission trees per forest
#' - `epidemic_size`: number of vertices per transmission tree
#' - `off_R_A`, `off_R_B`: mean of the offspring distribution for forests A and B
#' - `off_k_A`, `off_k_B`: dispersion of the offspring distribution for forests A and B
#' - `prop_reject_H0`: proportion of tests rejecting the null hypothesis
#'
#' @return A tibble summarising null hypothesis rejection rates per parameter combination.

rejection_summary <- readRDS("data/results_grid.rds") |>
  mutate(
    reject_H0 = p_value < 0.05,
    hypothesis = if_else(param_id_A == param_id_B, "H0", "H1")
  ) |>
  group_by(
    hypothesis,
    method,
    forest_size,
    epidemic_size,
    off_R_A,
    off_R_B,
    off_k_A,
    off_k_B
  ) |>
  summarise(
    prop_reject_H0 = mean(reject_H0, na.rm = TRUE),
    .groups = "drop"
  )
str(rejection_summary)

# ------------------------------------
#           plot grids
# ------------------------------------
#' @title plot grids for different forest sizes
#' @description
#' Generates and saves heatmaps of simulation results for various `forest_size` values.
#' Each heatmap visualises the proportion of rejected H0 across all conditions.
#' Plots are saved as SVG files in the "figures" directory.
map(
  c(20L, 50L, 100L, 200L),
  ~ {
    p <- plot_grid(rejection_summary, forest_size = .x)
    ggsave(
      filename = paste0("figures/grid_", .x, ".svg"),
      plot = p,
      width = 12,
      height = 12
    )
  }
)

# ------------------------------------
#           ROC curves
# ------------------------------------
roc_df <- readRDS("data/results_grid.rds") |>
  mutate(true_different = param_id_A != param_id_B) |>
  group_by(method, forest_size, epidemic_size) |>
  group_modify(~ compute_roc(.x)) |>
  ungroup() |>
  mutate(
    label_epidemic_size = "Epidemic size",
    label_forest_size = "Forest size"
  )


# Plot ROC curves
ggplot(roc_df, aes(x = FPR, y = TPR, color = method)) +
  facet_nested(
    cols = vars(label_epidemic_size, epidemic_size),
    rows = vars(label_forest_size, forest_size)
  ) +
  # Reference line
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5,
    linewidth = 0.25
  ) +

  # ROC line
  geom_line() +
  scale_color_manual(
    values = c("permanova" = "#fbac02ff", "chisq" = "#0248f8ff"),
    name = "Method:"
  ) +
  # Theme
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = TRUE,
    clip = "off"
  ) +
  labs(
    x = "1 - specificity",
    y = "Sensitivity",
    color = "Forest size:"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.just = "left",
    legend.spacing = unit(0.1, "cm"),
    axis.text.x = element_blank()
  )
ggsave(
  filename = "figures/roc.svg",
  width = 8,
  height = 8
)

# ------------------------------------
#           AUC
# ------------------------------------
auc_df <- readRDS("data/results_grid.rds") |>
  mutate(true_different = param_id_A != param_id_B) |>
  group_by(method, forest_size, epidemic_size) |>
  summarise(
    AUC = roc(
      response = true_different,
      predictor = 1 - p_value,
      quiet = TRUE
    )$auc |>
      as.numeric(),
    .groups = "drop"
  )
auc_df |>
  mutate(label_forest_size = "Forest size") |>
  ggplot(aes(x = factor(epidemic_size), y = AUC, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_nested(cols = vars(label_forest_size, forest_size)) +
  scale_fill_manual(
    values = c("permanova" = "#fbac02ff", "chisq" = "#0248f8ff"),
    name = "Method:"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    x = "Epidemic size",
    y = "Area Under the Curve (AUC)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.just = "left",
    legend.spacing = unit(0.1, "cm")
  )
ggsave(
  filename = "figures/auc.svg",
  width = 8,
  height = 6
)

# ------------------------------------
#           Power curves
# ------------------------------------

df <- readRDS("data/results_grid.rds") |>
  mutate(
    H0 = param_id_A == param_id_B,
    reject_H0 = p_value < 0.05,
    comparison_type = case_when(
      H0 ~ "H0",
      off_R_A != off_R_B & off_k_A == off_k_B ~ "\U0394 R₀",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "\U0394 \U1D458",
      TRUE ~ "\U0394 R₀ \U0026 \U0394 \U1D458"
    ),
    comparison_type = factor(
      comparison_type,
      levels = c(
        "\U0394 R₀ \U0026 \U0394 \U1D458",
        "\U0394 R₀",
        "\U0394 \U1D458",
        "H0"
      )
    )
  ) |>
  summarise(
    rejection_rate = mean(reject_H0),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )
glimpse(df)

df |>
  mutate(
    epidemic_size = factor(epidemic_size),
    label_forest_size = "Forest size",
    label_comparison_type = "Comparison type",
    hline = case_when(
      comparison_type == "H0" ~ 0.05,
      TRUE ~ 1
    )
  ) |>
  ggplot(aes(
    x = epidemic_size,
    y = rejection_rate,
    group = method
  )) +
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    rows = vars(label_comparison_type, comparison_type),
    scales = "free_y"
  ) +
  geom_hline(
    aes(yintercept = hline),
    linetype = "dotted",
    color = "grey50",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.35
  ) +
  geom_point(
    aes(shape = method),
    color = "black",
    size = 2.5
  ) +
  # scale_color_manual(
  #   values = c("permanova" = "#fbac02ff", "chisq" = "#0248f8ff"),
  #   labels = c("permanova" = "PERMANOVA", "chisq" = "\U03C7\U00B2 test"),
  #   name = "Method:"
  # ) +
  # scale_fill_manual(
  #   values = c("permanova" = "#fbac02ff", "chisq" = "#0248f8ff"),
  #   labels = c("permanova" = "PERMANOVA", "chisq" = "\U03C7\U00B2 test"),
  #   guide = "none"
  # ) +
  scale_shape_manual(
    values = c("permanova" = 16, "chisq" = 17),
    labels = c("permanova" = "PERMANOVA", "chisq" = "\U03C7\U00B2 test"),
    name = "Method:"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      comparison_type == "H0" ~ scale_y_continuous(
        limits = c(0, 0.03),
        breaks = seq(0, 0.03, by = 0.01),
        labels = c("0", "1", "2", "3")
      ),
      comparison_type == "\U0394 \U1D458" ~ scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      ),
      comparison_type == "\U0394 R₀" ~ scale_y_continuous(
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      ),
      comparison_type == "\U0394 R₀ \U0026 \U0394 \U1D458" ~ scale_y_continuous(
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      )
    )
  ) +
  labs(
    x = "Epidemic size",
    y = "\UFE6A rejecting H₀"
  ) +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.just = "left",
    legend.spacing = unit(0.1, "cm")
  )
ggsave(
  filename = "figures/power_curves.svg",
  width = 10,
  height = 10
)

# ------------------------------------
#           Deltas
# ------------------------------------
df <- readRDS("data/results_grid.rds")
plot_delta(df, delta = "R")
plot_delta(df, delta = "k")
df |>
  mutate(delta_k = abs(off_k_A - off_k_B)) |>
  pull(delta_k) |>
  unique() |>
  sort()
