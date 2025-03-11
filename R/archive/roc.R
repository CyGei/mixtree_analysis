


# Data --------------------------------------------------------------------
roc_df <- get_roc(results,
                  group_vars = c("sample_size", "epidemic_size", "method")) %>%
  mutate(col_title = "Epidemic Size", row_title = "Posterior Sample Size")

alpha_points <- roc_df %>%
  mutate(
    alpha_cat = case_when(
      abs(alpha - 0.0015) < 0.0001 ~ "0.0015",
      abs(alpha - 0.01) < 0.001 ~ "0.01",
      abs(alpha - 0.05) < 0.001 ~ "0.05",
      abs(alpha - 0.1) < 0.001 ~ "0.10",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(alpha_cat)) %>%
  group_by(method, sample_size, epidemic_size, alpha_cat) %>%
  slice_min(order_by = abs(alpha - as.numeric(alpha_cat)), n = 1) %>%
  ungroup() %>%
  mutate(alpha_cat = factor(alpha_cat, levels = c("0.0015", "0.01", "0.05", "0.10")))



# Plot 1 ------------------------------------------------------------------
ggplot() +
  # Facet setup
  ggh4x::facet_nested(
    cols = vars(col_title, epidemic_size),
    rows = vars(row_title, sample_size),
    scales = "fixed",
    space = "free",
    remove_labels = "none",
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
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
  # ROC curves from full dataset
  geom_line(
    data = roc_df,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      color = method
    ),
    linewidth = 0.5
  ) +
  # Points from filtered dataset
  geom_point(
    data = alpha_points,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      shape = alpha_cat,
      fill = method
    ),
    size = 1,
    stroke = 0.25,
    color = "black"
  ) +
  # Scale setup
  scale_color_manual(values = c("PERMANOVA" = "#af245a", "χ²" = "#1f9969"),
                     name = "Method") +
  scale_fill_manual(
    values = c("PERMANOVA" = "#af245a", "χ²" = "#1f9969"),
    name = "Method",
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "0.0015" = 21,
      "0.01" = 22,
      "0.05" = 23,
      "0.10" = 24
    ),
    name = "α",
    na.translate = FALSE
  ) +
   coord_equal(xlim = c(0, 0.3), ylim = c(0.7, 1)) +
  scale_x_continuous(labels = c("0", ".1", ".2", ".3")) +
  scale_y_continuous(labels = c(".7", ".8", ".9", "1")) +
  theme_bw() +
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
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines")
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)")

ggsave(
  "roc.png",
  width = 120,
  height = 120,
  units = "mm",
  dpi = 300
)



# Plot2 -------------------------------------------------------------------

ggplot() +
  # Facet setup
  ggh4x::facet_nested(
    rows = vars(row_title, epidemic_size),
    cols = vars(col_title, sample_size),
    scales = "fixed",
    space = "free",
    remove_labels = "none",
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
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
  # ROC curves from full dataset
  geom_line(
    data = roc_df,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      color = method
    ),
    linewidth = 0.5
  ) +
  # Points from filtered dataset
  geom_point(
    data = alpha_points,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      shape = alpha_cat,
      fill = method
    ),
    size = 1,
    stroke = 0.25,
    color = "black"
  ) +
  # Scale setup
  scale_color_manual(values = c("PERMANOVA" = "#af245a", "χ²" = "#1f9969"),
                     name = "Method") +
  scale_fill_manual(
    values = c("PERMANOVA" = "#af245a", "χ²" = "#1f9969"),
    name = "Method",
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "0.0015" = 21,
      "0.01" = 22,
      "0.05" = 23,
      "0.10" = 24
    ),
    name = "α",
    na.translate = FALSE
  ) +
  coord_equal(xlim = c(0, 0.3), ylim = c(0.7, 1)) +
  scale_x_continuous(labels = c("0", ".1", ".2", ".3")) +
  scale_y_continuous(labels = c(".7", ".8", ".9", "1")) +
  theme_bw() +
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
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines")
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)")

ggsave(
  "roc2.png",
  width = 120,
  height = 120,
  units = "mm",
  dpi = 300
)



# AUC ---------------------------------------------------------------------


roc_df %>%
  group_by(method, sample_size, epidemic_size) %>%
  summarise(AUC = mean(AUC)) %>%
  ggplot(aes(x = sample_size, y = AUC, color = method)) +
  facet_nested(
    rows = vars("Epidemic Size", epidemic_size),
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
  ) +
  geom_line(aes(group = method), linewidth = 0.5) +
  geom_point(
    aes(color = method),
    size = 2.5,
    fill = "white",
    shape = 21,
    stroke = 0.55
  ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("PERMANOVA" = "#1b81bc", "χ²" = "#bc1b81"),
                     name = "Method") +
  coord_cartesian(ylim = c(0.91, 1)) +
  scale_y_continuous(breaks = seq(0.92, 1, 0.02),
                     labels = c(".92", ".94", ".96", ".98", "1")) +
  theme_bw() +
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
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.25, "lines")
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "Posterior Sample Size", y = "Area Under Curve (AUC)")

ggplot2::ggsave(
  "auc.png",
  width = 100,
  height = 100,
  units = "mm",
  dpi = 300
)


roc_df %>%
  group_by(method, sample_size, epidemic_size) %>%
  summarise(AUC = mean(AUC)) %>%
  ggplot() +
  aes(x = sample_size, y = AUC) +
  facet_grid( ~ method) +
  geom_line(aes(group = epidemic_size, color = epidemic_size)) +
  geom_point(aes(color = epidemic_size), size = 3) +
  theme_bw() +
  scale_y_continuous(limits = c(0.92, 1)) +
  scale_color_manual(values = colorRampPalette(c("grey", "black"))(5)) +
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
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.25, "lines")
  ) +
  labs(x = "Posterior Sample Size", y = "Area Under Curve (AUC)", color = "Epidemic Size")

ggplot2::ggsave(
  "auc2.png",
  width = 120,
  height = 100,
  units = "mm",
  dpi = 300
)
