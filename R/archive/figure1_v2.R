alpha <- 0.05

method_comparison <- results %>%
  group_by(method, epidemic_size, sample_size, overlap_freq) %>%
  summarise(reject_percentage = mean(p_value < alpha) * 100,
            .groups = "drop") %>%
  mutate(target_percentage = ifelse(overlap_freq == "1", 100 * alpha, 100)) %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size")

# Create the dumbbell plot
ggplot(method_comparison, aes(y = overlap_freq, x = reject_percentage)) +
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
  geom_line(aes(group = interaction(epidemic_size, sample_size, overlap_freq)),
            color = "gray70",
            linewidth = 0.4) +
  geom_point(
    aes(color = method),
    size = 2,
    fill = "white",
    alpha = 0.75,
    shape = 21,
    stroke = 0.55
  ) +
  # Target reference line
  geom_point(
    aes(x = target_percentage),
    shape = 124,
    color = "black",
    size = 2
  ) +
  # Scale and theme
  scale_color_manual(values = c("PERMANOVA" = "#1b81bc", "χ²" = "#bc1b81")) +
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
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "Percentage Rejecting H₀",
       y = "Overlap Frequency",
       color = "Method",
       fill = "Method")

ggplot2::ggsave(
  "target_dumbells.png",
  width = 183,
  height = 120,
  units = "mm",
  dpi = 300
)
