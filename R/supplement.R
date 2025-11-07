# Supplementary figures.
source("R/packages.R")
source("R/helpers.R")
# ------------------------------------
#           AUC
# ------------------------------------
#' @title auc_df
#' @description
#' A tibble summarising Area Under the Curve (AUC) values for ROC curves
#' from the simulation study `results_grid.rds`.
#' It computes AUC for each combination of statistical test method,
#' forest size, and epidemic size.
#' Columns:
#' - `method`: statistical test used ("chisq" or "permanova")
#' - `forest_size`: number of transmission trees per forest
#' - `epidemic_size`: number of vertices per transmission tree
#' - `AUC`: Area Under the Curve value

auc_df <- readRDS("data/results_grid.rds") |>
  mutate(true_different = param_id_A != param_id_B) |>
  group_by(method, forest_size, epidemic_size) |>
  summarise(
    AUC = pROC::roc(
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
    name = "Method:"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    x = "Epidemic size",
    y = "Area Under the Curve (AUC)"
  ) +
  theme_mixtree()
ggsave(
  filename = "figures/auc.svg",
  width = 8,
  height = 6
)

# ------------------------------------
#           p-values by delta
# ------------------------------------

delta_df <- readRDS("data/results_grid.rds") |>
  mutate(
    H0 = param_id_A == param_id_B,
    reject_H0 = p_value < 0.05,
    comparison_type = case_when(
      H0 ~ "H₀",
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
        "H₀"
      )
    )
  ) |>
  summarise(
    mean_pval = mean(p_value),
    q25 = quantile(p_value, 0.25),
    q75 = quantile(p_value, 0.75),
    median_pval = median(p_value),
    ymin_fill = ifelse(unique(comparison_type) == "H₀", 0.05, -Inf),
    ymax_fill = ifelse(unique(comparison_type) == "H₀", Inf, 0.05),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )
glimpse(delta_df)

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
    color = method,
    fill = method
  )) +
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    rows = vars(label_comparison_type, comparison_type),
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
    position = position_dodge(width = 0.5),
    colour = "black",
    shape = 21,
    size = 3,
    stroke = 0.25
  ) +
  scale_color_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(mixtree_pal, "Target" = "grey90"),
    labels = c(mixtree_lab, "Target" = "Target"),
    name = "Method:"
  ) +
  labs(
    x = "Epidemic size",
    y = "p-value (median with IQR)"
  ) +
  theme_mixtree() +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5, colour = "black", linewidth = 0.25)
    )
  )
