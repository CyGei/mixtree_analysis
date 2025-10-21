df <- readRDS("data/results_grid.rds") |>
  filter(tree_id_A == tree_id_B) |>
  as_tibble() |>
  select(-c(starts_with("gt_"), duration, replicates)) |>
  mutate(across(-p_value, ~ factor(.x, levels = sort(unique(.x)))))
alpha_values <- seq(0, 1, by = 0.01)

roc_df <- expand_grid(df, alpha = alpha_values) |>
  mutate(across(
    c(epidemic_size, forest_size, method),
    ~ factor(.x, levels = sort(unique(.x)))
  )) |>
  group_by(epidemic_size, forest_size, method, alpha) |>
  mutate(
    true_null = off_R_A == off_R_B & off_k_A == off_k_B,
    true_diff = !true_null,
    test_reject = p_value < alpha
  ) |>
  summarise(
    TP = sum(true_diff & test_reject),
    FN = sum(true_diff & !test_reject),
    TN = sum(!true_diff & !test_reject),
    FP = sum(!true_diff & test_reject),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    .groups = "drop"
  )

roc_df |>
  ggplot() +
  aes(
    x = 1 - specificity,
    y = sensitivity,
    colour = method,
    group = interaction(
      method,
      epidemic_size,
      forest_size
    )
  ) +
  facet_grid(cols = vars(epidemic_size), rows = vars(forest_size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_line(size = 1) +
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
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    colour = "Method:"
  )


auc_df <- df |>
  mutate(
    true_diff = !(off_R_A == off_R_B & off_k_A == off_k_B)
  ) |>
  group_by(epidemic_size, forest_size, method) |>
  summarise(
    # Compute AUC using pROC::roc and pROC::auc
    AUC = pROC::auc(
      response = true_diff, # logical or binary (TRUE for positives)
      predictor = 1 - p_value, # higher values = stronger rejection
      quiet = TRUE
    ) |>
      as.numeric(),
    .groups = "drop"
  )


auc_df |>
  ggplot() +
  aes(x = forest_size, y = AUC) +
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
    x = "Forest Size",
    y = "Area Under Curve (AUC)",
    color = "Epidemic Size:"
  )
