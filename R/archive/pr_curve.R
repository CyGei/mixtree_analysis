

# Calculate precision and recall at different thresholds
alpha <- 0.05
pr_curves <- results %>%
  mutate(true_positive = overlap_freq != "1", prediction = p_value < alpha) %>%
  group_by(method, epidemic_size, sample_size) %>%
  summarise(
    precision = sum(prediction & true_positive) / sum(prediction),
    recall = sum(prediction & true_positive) / sum(true_positive),
    .groups = "drop"
  ) %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size") %>%
  ggplot(aes(
    x = recall,
    y = precision,
    group = interaction(epidemic_size, sample_size, method)
  )) +

  geom_line(aes(group = epidemic_size, linetype = epidemic_size)) +
  geom_point(aes(color = sample_size, shape = method))
pr_curves

pr_curves <- results %>%
  mutate(true_positive = overlap_freq != "1", prediction = p_value < alpha) %>%
  group_by(method, epidemic_size, sample_size) %>%
  summarise(
    precision = sum(prediction & true_positive) / sum(prediction),
    recall = sum(prediction & true_positive) / sum(true_positive),
    .groups = "drop"
  ) %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size") %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_line(aes(linetype = epidemic_size, group = interaction(epidemic_size, method))) +
  geom_point(aes(color = sample_size), size = 2) +
  scale_color_viridis_d(option = "D")
pr_curves +
  facet_grid( ~ method)+
  labs(x = "Sensitivity", y = "Positive Predictive Value")+
  theme_bw()

x = results %>%
  mutate(true_positive = overlap_freq != "1") %>%
  nest_by(method, epidemic_size, sample_size) %>%
  mutate(pr_data = list(map_df(thresholds, ~ {
    prediction <- data$p_value < .x
    tp <- sum(prediction & data$true_positive)
    precision <- if (sum(prediction) > 0)
      tp / sum(prediction)
    else
      NA
    recall <- if (sum(data$true_positive) > 0)
      tp / sum(data$true_positive)
    else
      NA
    tibble(threshold = .x,
           precision = precision,
           recall = recall)
  }))) %>%
  unnest(pr_data) %>%
  ungroup()

x %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size") %>%
  ggplot(aes(x = recall, y = precision)) +
  ggh4x::facet_nested(
    cols = vars(col_title, sample_size),
    rows = vars(row_title, epidemic_size),
    scales = "free",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(
      colour = "black",
      linewidth = 0.05,
      lineend = "square"
    )
  ) +
  geom_line(aes(color = threshold,
                group = interaction(epidemic_size, method))) +
  #geom_point(aes(shape = method), size = 2) +
  scale_color_viridis_c(option = "rocket") +
  scale_x_continuous(limits = c(0.8, 1)) +
  scale_y_continuous(limits = c(0.75, 1)) +
  labs(x = "Recall", y = "Precision") +
  theme_bw()


x %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size") %>%
  ggplot(aes(x = recall, y = precision)) +
  ggh4x::facet_nested(
    cols = vars(col_title, sample_size),
    rows = vars(row_title, epidemic_size),
    scales = "fixed",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(colour = "black", linewidth = 0.05, lineend = "square")
  ) +
  geom_line(aes(color = method,
                alpha = threshold,  # Use alpha for threshold
                group = interaction(epidemic_size, sample_size, method))) +
  scale_color_manual(values = c("red", "blue")) +  # Set colors
  scale_alpha_continuous(range = c(0.3, 1)) +  # Adjust transparency range
  scale_x_continuous(limits = c(0.8, 1)) +
  scale_y_continuous(limits = c(0.75, 1)) +
  labs(x = "Recall (sensitivity)", y = "Precision") +
  theme_bw()
#facet_grid(~ method) +
# ggh4x::facet_nested(
#   rows = vars(row_title, epidemic_size),
#   scales = "fixed",
#   space = "fixed",
#   remove_labels = "none",
#   nest_line = element_line(
#     colour = "black",
#     linewidth = 0.05,
#     lineend = "square"
#   )
# ) +
# Summarize roc_df to get mean AUC by sample_size, epidemic_size, and method
roc_summary <- roc_df %>%
  group_by(sample_size, epidemic_size, method) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = "drop")

# Create heatmap
ggplot(roc_summary,
       aes(x = sample_size, y = epidemic_size, fill = mean_auc)) +
  geom_tile(color = "black") +  # Tiles for heatmap
  facet_wrap( ~ method) +  # Separate panels for each method
  scale_fill_viridis_c(option = "plasma") +  # Color scale
  labs(title = "Mean AUC by Sample Size and Epidemic Size", x = "Sample Size", y = "Epidemic Size") +
  theme_bw()
