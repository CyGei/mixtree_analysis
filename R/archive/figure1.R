




# FIGURE 1 ----------------------------------------------------------------
alpha <- 0.05
annotations_df <- results %>%
  group_by(method, epidemic_size, sample_size, overlap_freq) %>%
  summarise(
    reject = mean(p_value < alpha) * 100,
    accept = mean(p_value >= alpha) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c("reject", "accept"),
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(
    row_title = "Epidemic Size",
    col_title = "Posterior Sample Size",
    target_percentage = ifelse(overlap_freq != "1", 100, 100 * alpha),
    raw_error = abs(percentage - target_percentage),
    relative_error = raw_error / target_percentage * 100,
    relative_error_log = log1p(relative_error)
  )


##############################
# target_percentage
##############################

# Create modified data with nudge information
label_data <- annotations_df %>%
  filter(category == "reject") %>%
  filter(percentage < 90) %>%
  mutate(x_nudge = case_when(
    method == "PERMANOVA" ~ -0.1,
    method == "χ²" & overlap_freq == "1" ~ -0.1,
    TRUE ~ 0
  ))

annotations_df %>%
  filter(category == "reject") %>%
  ggplot(aes(
    x = overlap_freq,
    y = percentage,
    color = method,
    group = method
  )) +
  ggh4x::facet_nested(
    rows = vars(row_title, epidemic_size),
    cols = vars(col_title, sample_size),
    scales = "fixed",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
  ) +
  geom_line(alpha = 0.65) +
  geom_text(
    data = . %>% filter(percentage < 90),
    aes(label = sprintf("%.1f", percentage),
        hjust = ifelse(method == "PERMANOVA", 1, 0),
        vjust = ifelse(method == "PERMANOVA", 1, 0)
    ),
    size = 2,
    show.legend = FALSE) +
  # geom_label(
  #   data = label_data,
  #   aes(
  #     x = as.numeric(as.factor(overlap_freq)) + x_nudge,  # Apply the nudge
  #     label = sprintf("%.1f", percentage),
  #     hjust = case_when(
  #       method == "PERMANOVA" ~ 1,
  #       method == "χ²" & overlap_freq == "1" ~ 1,
  #       TRUE ~ 0
  #     ),
  #     vjust = case_when(
  #       method == "PERMANOVA" & overlap_freq != "1" ~ 1,
  #       method == "PERMANOVA" & overlap_freq == "1" ~ 0.5,
  #       method == "χ²" & overlap_freq == "1" ~ 0,
  #       TRUE ~ 0
  #     )
  #   ),
  #   size = 2,
  #   show.legend = FALSE,
  #   label.size = 0.1
  # ) +
  geom_point(size = 2,
             shape = 21,
             fill = "white") +
  geom_point(aes(y = target_percentage, shape = "Target"),
             color = "black",
             size = 1) +
  coord_cartesian(ylim = c(-2, 104), clip = "off") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -7.5),
    legend.spacing.x = unit(0, "mm"),
    # Controls spacing between key and label
    strip.background = element_rect(
      fill = "#f8f4f2",
      colour = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.spacing = unit(0.15, "lines"),
    panel.border = element_rect(colour = "black", linewidth = 0.15)
  ) +
  scale_shape_manual(values = c("Target" = 8), name = NULL) +
  scale_color_manual(values = c("PERMANOVA" = "#af245a", "χ²" = "#1f9969"),
                     name = "       Method:") +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(size = 2.25),
      keywidth = unit(0, "mm")
    ),
    color = guide_legend(order = 2)
  ) +
  labs(x = "Overlap Frequency", y = "Percentage of tests rejecting H₀")

ggsave(
  "target_percentage.png",
  width = 183,
  height = 120,
  units = "mm",
  dpi = 300
)

##############################
# relative_error (log-transformed)
##############################

annotations_df %>%
  filter(category == "reject") %>%
  ggplot(aes(
    x = overlap_freq,
    y = relative_error_log,
    color = method,
    group = method
  )) +
  ggh4x::facet_nested(
    rows = vars(row_title, epidemic_size),
    cols = vars(col_title, sample_size),
    scales = "fixed",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
  ) +
  geom_line(alpha = 0.65) +
  geom_point(size = 2,
             shape = 21,
             fill = "white") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -10),
    strip.background = element_rect(
      fill = "#f8f4f2",
      colour = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(family = "Fira Code", size = 10),
    panel.spacing = unit(0.15, "lines"),
  ) +
  scale_color_manual(values = c("#1b81bc", "#bc1b81")) +
  labs(x = "Overlap Frequency", y = "log-transformed relative error")
ggplot2::ggsave(
  "relative_error_log.png",
  width = 183,
  height = 120,
  units = "mm",
  dpi = 300
)
