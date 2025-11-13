format_k <- function(k) {
  case_when(
    k >= 1e5 ~ "10⁵",
    TRUE ~ as.character(k)
  )
}

delta_df <- results_grid |>
  filter(H1) |>
  mutate(
    # Calculate deltas
    delta_R0 = abs(off_R_A - off_R_B),
    delta_k = abs(off_k_A - off_k_B),

    # R0 comparison labels
    R0_comparison = sprintf(
      "%.1f v %.1f",
      pmin(off_R_A, off_R_B),
      pmax(off_R_A, off_R_B)
    ),

    # k comparison labels with scientific notation for large values
    k_comparison = sprintf(
      "%s v %s",
      format_k(pmin(off_k_A, off_k_B)),
      format_k(pmax(off_k_A, off_k_B))
    ),

    # Create ordered factors based on delta magnitude
    R0_comparison = factor(R0_comparison) |>
      forcats::fct_reorder(delta_R0),

    k_comparison = factor(k_comparison) |>
      forcats::fct_reorder(delta_k),
  ) |>
  summarise(
    rejection_rate = mean(reject_H0),
    .by = c(
      method,
      R0_comparison,
      k_comparison,
      comparison_type,
      forest_size,
      epidemic_size
    )
  )

# heatmap
delta_df |>
  filter(forest_size == 200) |>
  mutate(
    epidemic_size_label = "Epidemic size",

    method = case_when(
      method == "permanova" ~ "PERMANOVA",
      method == "chisq" ~ "\U03C7\U00B2 test"
    )
  ) |>
  ggplot(aes(x = R0_comparison, y = k_comparison, fill = rejection_rate)) +
  facet_nested(
    rows = vars(method),
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
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    name = "Sensitivity",
    option = "viridis",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 10,
      barheight = 0.5,
      frame.linewidth = 0.5,
      ticks.linewidth = 0.1,
    )
  ) +
  labs(
    x = "\U0394 R₀",
    y = "\U0394 \U1D458"
  ) +
  theme_mixtree() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "figures/delta_heatmap.png",
  width = 7.5,
  height = 8.5
)

glimpse(delta_df)
delta_df |>
  mutate(
    epidemic_size = factor(epidemic_size),
    label_forest_size = "Forest size",
    label_comparison_type = "Comparison type",
  ) |>
  ggplot(aes(
    x = epidemic_size,
    y = rejection_rate,
    group = interaction(method, R0_comparison, k_comparison),
    color = method,
    fill = method
  )) +
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    rows = vars(label_comparison_type, comparison_type),
  ) +
  geom_line(
    linewidth = 0.35
  ) +
  scale_color_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    guide = "none"
  ) +
  scale_fill_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    name = "Method:"
  ) +
  labs(
    x = "Epidemic size",
    y = "\UFE6A rejecting H₀"
  ) +
  theme_mixtree()
