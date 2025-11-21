# Supplementary figures.
source("R/packages.R")
source("R/helpers.R")

# ------------------------------------
#          Data
# ------------------------------------
results_grid <- readRDS("data/results_grid.rds") |>
  mutate(
    H1 = param_id_A != param_id_B,
    reject_H0 = p_value < 0.05,

    delta_R0 = abs(off_R_A - off_R_B),
    delta_k = abs(off_k_A - off_k_B),

    # Annotations:
    comparison_type = case_when(
      !H1 ~ "H₀",
      off_R_A != off_R_B & off_k_A == off_k_B ~ "\U0394 R₀",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "\U0394 \U1D458",
      .default = "\U0394 R₀ \U0026 \U0394 \U1D458"
    ),
    comparison_type = factor(
      comparison_type,
      levels = c(
        "\U0394 R₀ \U0026 \U0394 \U1D458",
        "\U0394 R₀",
        "\U0394 \U1D458",
        "H₀"
      )
    ),

    R0_labels = sprintf(
      "%.1f v %.1f",
      pmin(off_R_A, off_R_B),
      pmax(off_R_A, off_R_B)
    ) |>
      forcats::fct_reorder(delta_R0),

    k_labels = sprintf(
      "%s v %s",
      format_k(pmin(off_k_A, off_k_B)),
      format_k(pmax(off_k_A, off_k_B))
    ) |>
      forcats::fct_reorder(delta_k),

    delta_R0_labels = case_when(
      delta_R0 == 0 ~ "ΔR₀ = 0",
      TRUE ~ paste0("ΔR₀ = ", delta_R0)
    ) |>
      factor(levels = c("ΔR₀ = 0", "ΔR₀ = 0.5", "ΔR₀ = 1", "ΔR₀ = 1.5")),

    delta_k_labels = case_when(
      delta_k == 0 ~ "Δ\U1D458 = 0",
      delta_k >= 99999 ~ "Δ\U1D458 = 10⁵",
      TRUE ~ paste0("Δ\U1D458 = ", delta_k)
    ) |>
      as.factor()
  )

# ------------------------------------
#           ROC curves
# ------------------------------------
#' @title compute_roc()
#' @description
#' Computes true positive rates (TPR) and false positive rates (FPR) across varying p-value thresholds.
#' A true positive is defined as correctly rejecting the null hypothesis when forests differ in parameters.
#' A false positive is incorrectly rejecting the null hypothesis when forests are generated under the same parameters.
#' - `TPR`: true positive rate (sensitivity)
#' - `FPR`: false positive rate (1 - specificity)
results_grid |>
  group_by(method, forest_size, epidemic_size) |>
  group_modify(~ compute_roc(.x)) |>
  ungroup() |>
  mutate(
    label_epidemic_size = "Epidemic size",
    label_forest_size = "Forest size"
  ) |>
  ggplot(aes(x = FPR, y = TPR, color = method)) +
  facet_nested(
    cols = vars(label_epidemic_size, epidemic_size),
    rows = vars(label_forest_size, forest_size),
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      text_y = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      by_layer_x = TRUE,
      by_layer_y = TRUE
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  geom_line() +
  scale_color_manual(
    values = c("permanova" = "#e68613", "chisq" = "#1f65cc"),
    labels = mixtree_lab,
    name = "Method:"
  ) +
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
  theme_mixtree() +
  theme(axis.text.x = element_blank())

ggsave(
  filename = "figures/roc.png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

# ------------------------------------
#           AUC
# ------------------------------------
#' @title auc_df
#' @description
#' A tibble summarising Area Under the Curve (AUC) values for ROC curves.
#' It computes AUC for each combination of statistical test method,
#' forest size, and epidemic size.
#' Columns:
#' - `method`: statistical test used ("chisq" or "permanova")
#' - `forest_size`: number of transmission trees per forest
#' - `epidemic_size`: number of vertices per transmission tree
#' - `AUC`: Area Under the Curve value

auc_df <- results_grid |>
  group_by(method, forest_size, epidemic_size) |>
  summarise(
    AUC = pROC::roc(
      response = H1,
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
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      by_layer_x = TRUE
    )
  ) +
  scale_fill_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    name = "Method:"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    x = "Epidemic size",
    y = "Area Under the Curve (AUC)"
  ) +
  theme_mixtree()

ggsave(
  filename = "figures/auc.png",
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 300
)

# ------------------------------------
#           p-values by delta
# ------------------------------------

delta_df <- results_grid |>
  summarise(
    mean_pval = mean(p_value),
    q25 = quantile(p_value, 0.25),
    q75 = quantile(p_value, 0.75),
    median_pval = median(p_value),
    ymin_fill = ifelse(unique(comparison_type) == "H₀", 0.05, -Inf),
    ymax_fill = ifelse(unique(comparison_type) == "H₀", Inf, 0.05),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )

delta_df |>
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
    y = median_pval,
    group = interaction(method, epidemic_size, forest_size, comparison_type),
  )) +
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    rows = vars(label_comparison_type, comparison_type),
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      text_y = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      by_layer_x = TRUE,
      by_layer_y = TRUE
    )
  ) +
  # Target geoms
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin_fill,
      ymax = ymax_fill,
      fill = "Target"
    ),
    alpha = 0.5,
    color = "black",
    linewidth = 0.1
    #inherit.aes = FALSE
  ) +
  # Method geoms
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    aes(shape = method),
    position = position_dodge(width = 0.5),
    size = 2.5,
    stroke = 0.25
  ) +
  scale_shape_manual(
    values = c("permanova" = 19, "chisq" = 15),
    labels = mixtree_lab,
    name = "Method:"
  ) +
  scale_fill_manual(
    values = c(mixtree_pal, "Target" = "grey90"),
    labels = c(mixtree_lab, "Target" = "Target"),
    name = ""
  ) +
  labs(
    x = "Epidemic size",
    y = "p-value (median with IQR)"
  ) +
  theme_mixtree() +
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  )

ggsave(
  filename = "figures/delta_pval.png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

# ------------------------------------
#      delta_R0 lineplots
# ------------------------------------
#' @title delta_R0_lineplots
#' @description
#' Figure illustrating the sensitivity (rejection rate) of the PERMANOVA test
#' when comparing forests that differ only in R0 (ΔR0),
#' across varying levels of overdispersion (k) and epidemic sizes.
delta_df <-
  results_grid |>
  filter(
    H1,
    forest_size == 200,
    method == "permanova",
    comparison_type == "\U0394 R₀",
    off_k_A == off_k_B
  ) |>
  summarise(
    rejection_rate = mean(reject_H0),
    .by = c(
      epidemic_size,
      off_k_A,
      delta_R0
    )
  )

delta_df |>
  mutate(
    off_k_A = factor(off_k_A),
    off_k_A = ifelse(off_k_A == 100000, "10⁵", as.character(off_k_A)),
    delta_R0 = factor(delta_R0),
    epidemic_size = factor(epidemic_size),
    epidemic_size_label = "Epidemic size"
  ) |>
  ggplot(
    aes(
      x = off_k_A,
      y = rejection_rate,
      group = interaction(epidemic_size, delta_R0)
    ),
  ) +
  facet_nested(
    cols = vars(epidemic_size_label, epidemic_size),
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      text_y = list(
        element_text(size = 16)
      ),
      by_layer_x = TRUE,
      by_layer_y = TRUE
    )
  ) +
  geom_line(
    aes(color = delta_R0),
    linewidth = 0.75
  ) +
  geom_point(
    aes(fill = delta_R0),
    color = "white",
    shape = 21,
    size = 4.2,
    stroke = 0.5
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    name = "\U0394 R₀"
  ) +
  scale_color_brewer(
    palette = "Dark2",
    name = "\U0394 R₀"
  ) +
  theme_mixtree() +
  labs(
    x = "\U1D458",
    y = "Sensitivity"
  )

ggsave(
  "figures/delta_R0_lineplots_m200.png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

delta_df |>
  mutate(
    delta_R0_below_1 = delta_R0 < 1,
    overdispersed = off_k_A < 0.5
  ) |>
  summarise(
    mean = mean(rejection_rate),
    .by = c(delta_R0_below_1, overdispersed)
  )


# ------------------------------------
#           grid plot
# ------------------------------------
rejection_summary <- results_grid |>
  mutate(
    reject_H0 = p_value < 0.05
  ) |>
  summarise(
    prop_reject_H0 = mean(reject_H0),
    .by = c(
      "method",
      "forest_size",
      "epidemic_size",
      "off_R_A",
      "off_R_B",
      "off_k_A",
      "off_k_B"
    )
  )


plot_list <- map(
  c(20L, 50L, 100L, 200L),
  ~ {
    p <- plot_grid(
      rejection_summary,
      forest_size = .x,
      base_size = 10
    )
    ggsave(
      filename = paste0("figures/grid", .x, ".png"),
      plot = p,
      width = 7.5,
      height = 7.5,
      units = "in",
      dpi = 500
    )
    return(p)
  }
)

# Combine plots into one figure + collect legends
patchwork::wrap_plots(plotlist = plot_list, ncol = 1, guides = "collect") &
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "figures/grid_combined.png",
  width = 7.5,
  height = 10,
  units = "in",
  dpi = 300
)
