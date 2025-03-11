# TARGET PERCENTAGE ----------------------------------------------------------------
alpha <- 0.05

# Prepare data with numeric x
plot_data <- results %>%
  group_by(method, epidemic_size, sample_size, overlap_freq) %>%
  summarise(reject = mean(p_value < alpha) * 100, .groups = "drop") %>%
  mutate(
    col_title = "Epidemic Size",
    row_title = "Forest Size",
    target_percentage = ifelse(overlap_freq != "1", 100, 100 * alpha),
    x_num = as.numeric(as.character(overlap_freq)) # Convert overlap_freq to numeric
  )

plot_data %>%
  ggplot(aes(x = x_num, y = reject, color = method, group = method)) +
  facet_nested(
    cols = vars(col_title, epidemic_size),
    rows = vars(row_title, sample_size),
    scales = "fixed",
    nest_line = element_line(colour = "black", linewidth = 0.05, lineend = "square")
  ) +
  geom_line(alpha = 0.65) +
  geom_point(size = 1.5, shape = 21, fill = "white") +
  geom_point(aes(y = target_percentage, shape = "Target"), color = "black", size = 1) +
  coord_cartesian(ylim = c(-2, 104), clip = "off") +
  scale_x_continuous(breaks = unique(plot_data$x_num), labels = unique(plot_data$overlap_freq)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -7.5),
    legend.spacing.x = unit(0, "mm"),
    strip.background = element_rect(fill = "#f8f4f2", colour = "black", linewidth = 0.5),
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.spacing = unit(0.15, "lines"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  ) +
  scale_shape_manual(values = c("Target" = 8), name = NULL) +
  scale_color_manual(values = c("PERMANOVA" = "#af245a", "χ² test" = "#1f9969"), name = "       Method:") +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 2.25)),
    color = guide_legend(order = 2)
  ) +
  labs(x = "Overlap Frequency", y = "Percentage of tests rejecting H₀")

# Save the plot
ggsave(
  "target_percentage.png",
  width = 183,
  height = 120,
  units = "mm",
  dpi = 300
)


# ROC --------------------------------------------------------------------
roc_df <- get_roc(results,
  group_vars = c("sample_size", "epidemic_size", "method")
) %>%
  mutate(col_title = "Epidemic Size", row_title = "Forest Size")

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
  scale_color_manual(
    values = c("PERMANOVA" = "#af245a", "χ² test" = "#1f9969"),
    name = "Method:"
  ) +
  scale_fill_manual(
    values = c("PERMANOVA" = "#af245a", "χ² test" = "#1f9969"),
    name = "Method:",
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "0.0015" = 21,
      "0.01" = 22,
      "0.05" = 23,
      "0.10" = 24
    ),
    name = "α:",
    na.translate = FALSE
  ) +
  coord_equal(xlim = c(0, 0.3), ylim = c(0.7, 1)) +
  scale_x_continuous(labels = c("0", "0.1", "0.2", "0.3")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),
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
    panel.spacing = unit(0.3, "lines")
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



# AUC ---------------------------------------------------------------------
roc_df %>%
  group_by(method, sample_size, epidemic_size) %>%
  summarise(AUC = mean(AUC)) %>%
  mutate(method = factor(method, levels = c("χ² test", "PERMANOVA"))) %>%
  ggplot() +
  aes(x = sample_size, y = AUC) +
  facet_grid(~method) +
  geom_line(aes(group = epidemic_size, color = epidemic_size)) +
  geom_point(aes(color = epidemic_size), size = 2.5) +
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
  labs(
    x = "Forest Size", y = "Area Under Curve (AUC)",
    color = "Epidemic Size:"
  )

ggplot2::ggsave(
  "auc.png",
  width = 120,
  height = 100,
  units = "mm",
  dpi = 300
)



# precision - recall curve

results %>%
  mutate(
    actual_positive  = overlap_freq != "1",
    predicted_positive = p_value < alpha
  ) %>%
  group_by(method, epidemic_size, sample_size) %>%
  summarise(
    precision = sum(predicted_positive & actual_positive ) / sum(predicted_positive),
    recall = sum(predicted_positive & actual_positive ) / sum(actual_positive),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = recall, y = precision, group = interaction(epidemic_size, method))) +
  facet_grid(~method) +
  geom_line(aes(colour = epidemic_size)) +
  geom_point(aes(fill = sample_size), size = 2, shape = 21, color = "black") +
  scale_color_manual(values = colorRampPalette(c("grey", "black"))(5)) +
  scale_fill_viridis_d(option = "D") +
  # scale_x_continuous(
  #   labels = c("0.80", "0.85", "0.90", "0.95", "1")
  # )+
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
    panel.spacing = unit(0.9, "lines")
  ) +
  labs(
    x = "Sensitivity (Recall)", y = "Positive Predictive Value (Precision)",
    color = "Epidemic Size:",
    fill = "Forest Size:"
  )

  ggplot2::ggsave(
    "pr_curve.png",
    width = 120,
    height = 100,
    units = "mm",
    dpi = 300
  )
