df <- readRDS("data/results_grid.rds")

# Step 1: Identify whether each comparison is testing equal or different forests
df_classified <- df |>
  mutate(
    same_params = param_id_A == param_id_B,
    # Expected outcome
    expect_reject = !same_params, # Reject H0 when params differ
    # Actual outcome
    reject_H0 = p_value < 0.05,
    # Performance
    correct = (expect_reject & reject_H0) | (!expect_reject & !reject_H0)
  )

# Step 2: Compute rejection rates by comparison type
rejection_summary <- df_classified |>
  summarise(
    n_tests = n(),
    n_reject = sum(reject_H0),
    pct_reject = mean(reject_H0) * 100,
    # Separate by what we expect
    n_same_params = sum(same_params),
    n_diff_params = sum(!same_params),
    # Type I error (rejecting when should not)
    type1_error = sum(same_params & reject_H0) / sum(same_params) * 100,
    # Power (correctly rejecting when should)
    power = sum(!same_params & reject_H0) / sum(!same_params) * 100,
    .by = c(
      method,
      forest_size,
      epidemic_size,
      off_R_A,
      off_k_A,
      off_R_B,
      off_k_B
    )
  )

# Step 3: For plotting - focus on specific contrasts
# Example: Power to detect difference when off_R differs
power_by_R_diff <- df_classified |>
  filter(off_k_A == off_k_B, off_R_A != off_R_B) |> # k same, R different
  summarise(
    power = mean(reject_H0) * 100,
    n = n(),
    .by = c(method, forest_size, epidemic_size, off_R_A, off_R_B)
  )

# Example: Power to detect difference when off_k differs
power_by_k_diff <- df_classified |>
  filter(off_R_A == off_R_B, off_k_A != off_k_B) |> # R same, k different
  summarise(
    power = mean(reject_H0) * 100,
    n = n(),
    .by = c(method, forest_size, epidemic_size, off_k_A, off_k_B)
  )

# Step 4: Type I error rate (false positives)
type1_error_rate <- df_classified |>
  filter(param_id_A == param_id_B) |>
  summarise(
    type1_error = mean(reject_H0) * 100,
    n = n(),
    .by = c(method, forest_size, epidemic_size, off_R_A, off_k_A)
  )
# Create a comprehensive summary for heatmaps
plot_data <- df_classified |>
  mutate(
    comparison_type = case_when(
      param_id_A == param_id_B ~ "Same Parameters",
      off_R_A != off_R_B & off_k_A == off_k_B ~ "Different R only",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "Different k only",
      TRUE ~ "Different R and k"
    )
  ) |>
  summarise(
    rejection_rate = mean(reject_H0) * 100,
    n_tests = n(),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )

# Example plot
library(ggplot2)

ggplot(
  plot_data,
  aes(x = forest_size, y = rejection_rate, color = comparison_type)
) +
  geom_line() +
  geom_point() +
  facet_grid(epidemic_size ~ method, labeller = label_both) +
  labs(
    title = "Test Performance Across Conditions",
    y = "% Rejecting H0",
    x = "Forest Size"
  ) +
  theme_classic()


library(tidyverse)

# Step 1: Create a dataset with parameter differences and performance
comparison_data <- df |>
  mutate(
    # Calculate parameter differences (absolute and signed)
    delta_R = off_R_B - off_R_A,
    delta_k = off_k_B - off_k_A,
    abs_delta_R = abs(delta_R),
    abs_delta_k = abs(delta_k),
    # Log-scale for k (since k ranges from 0.1 to 1e5)
    log_k_A = log10(off_k_A),
    log_k_B = log10(off_k_B),
    delta_log_k = log_k_B - log_k_A,
    abs_delta_log_k = abs(delta_log_k),
    # Test outcome
    reject_H0 = p_value < 0.05,
    # Comparison type
    same_params = (delta_R == 0 & delta_k == 0),
    delta_R_only = (delta_R != 0 & delta_k == 0),
    delta_k_only = (delta_R == 0 & delta_k != 0),
    both_differ = (delta_R != 0 & delta_k != 0)
  ) |>
  filter(!same_params) # Exclude self-comparisons

# Step 2: Summarize power by parameter differences
power_summary <- comparison_data |>
  summarise(
    power = mean(reject_H0) * 100,
    n_tests = n(),
    se = sqrt(power / 100 * (1 - power / 100) / n_tests) * 100, # Standard error
    .by = c(
      method,
      epidemic_size,
      forest_size,
      abs_delta_R,
      abs_delta_log_k,
      delta_R_only,
      delta_k_only
    )
  )

# Step 3: Create categorical bins for better visualization
power_summary_binned <- comparison_data |>
  mutate(
    # Bin the differences for cleaner visualization
    R_diff_bin = case_when(
      abs_delta_R == 0 ~ "0",
      abs_delta_R == 0.5 ~ "0.5",
      abs_delta_R == 1.0 ~ "1.0",
      abs_delta_R == 1.5 ~ "1.5",
      TRUE ~ as.character(abs_delta_R)
    ),
    k_diff_bin = case_when(
      abs_delta_log_k == 0 ~ "0",
      abs_delta_log_k < 1 ~ "<1 log unit",
      abs_delta_log_k < 2 ~ "1-2 log units",
      abs_delta_log_k < 4 ~ "2-4 log units",
      TRUE ~ "≥4 log units"
    ),
    # Create interaction label
    diff_type = case_when(
      abs_delta_R == 0 ~ paste0("ΔR=0, Δlog₁₀k=", round(abs_delta_log_k, 1)),
      abs_delta_log_k == 0 ~ paste0("ΔR=", abs_delta_R, ", Δlog₁₀k=0"),
      TRUE ~ paste0("ΔR=", abs_delta_R, ", Δlog₁₀k=", round(abs_delta_log_k, 1))
    )
  ) |>
  summarise(
    power = mean(reject_H0) * 100,
    n_tests = n(),
    .by = c(
      method,
      epidemic_size,
      forest_size,
      abs_delta_R,
      abs_delta_log_k,
      diff_type
    )
  )

# Filter to show only meaningful comparisons
plot_data_1 <- power_summary_binned |>
  filter(
    abs_delta_R > 0 | abs_delta_log_k > 0 # Exclude (0,0)
  )

ggplot(
  plot_data_1,
  aes(
    x = factor(abs_delta_R),
    y = factor(round(abs_delta_log_k, 1)),
    fill = power
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(power, 0)), size = 2.5, color = "white") +
  scale_fill_viridis_c(limits = c(0, 100), option = "plasma") +
  facet_nested(
    cols = vars(epidemic_size, method),
    rows = vars(forest_size),
    labeller = labeller(
      epidemic_size = ~ paste("Epidemic:", .),
      forest_size = ~ paste("Forest:", .)
    )
  ) +
  labs(
    title = "Power to Detect Parameter Differences",
    subtitle = "Joint effect of ΔR and Δlog₁₀(k)",
    x = "Absolute Difference in R (ΔR)",
    y = "Absolute Difference in log₁₀(k)",
    fill = "Power (%)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# Separate plot for R-only and k-only effects
plot_data_2a <- comparison_data |>
  filter(delta_R_only) |>
  summarise(
    power = mean(reject_H0) * 100,
    .by = c(method, epidemic_size, forest_size, abs_delta_R)
  )

plot_data_2b <- comparison_data |>
  filter(delta_k_only) |>
  summarise(
    power = mean(reject_H0) * 100,
    .by = c(method, epidemic_size, forest_size, abs_delta_log_k)
  )

# Plot R effect
p1 <- ggplot(
  plot_data_2a,
  aes(
    x = abs_delta_R,
    y = power,
    color = factor(forest_size),
    linetype = factor(epidemic_size)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~method) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Power vs ΔR (when k is constant)",
    x = "Absolute Difference in R",
    y = "Power (%)",
    color = "Forest Size",
    linetype = "Epidemic Size"
  ) +
  theme_minimal()

# Plot k effect
p2 <- ggplot(
  plot_data_2b,
  aes(
    x = abs_delta_log_k,
    y = power,
    color = factor(forest_size),
    linetype = factor(epidemic_size)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~method) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Power vs Δlog₁₀(k) (when R is constant)",
    x = "Absolute Difference in log₁₀(k)",
    y = "Power (%)",
    color = "Forest Size",
    linetype = "Epidemic Size"
  ) +
  theme_minimal()

library(patchwork)
p1 / p2


library(tidyverse)
library(ggh4x) # For facet_nested

library(tidyverse)
library(ggh4x)

# Step 1: Compute rejection rates for actual comparisons
rejection_data <- df |>
  summarise(
    pct_reject_H0 = mean(p_value < 0.05) * 100,
    n_tests = n(),
    .by = c(
      off_R_A,
      off_k_A,
      off_R_B,
      off_k_B,
      epidemic_size,
      forest_size,
      method
    )
  )

# Step 2: Create symmetric version by swapping A and B
rejection_symmetric <- bind_rows(
  rejection_data,
  rejection_data |>
    rename(
      off_R_A = off_R_B,
      off_R_B = off_R_A,
      off_k_A = off_k_B,
      off_k_B = off_k_A
    )
) |>
  distinct() # Remove duplicates where A == B

# Step 3: Plot
ggplot(
  rejection_symmetric,
  aes(x = factor(off_k_A), y = factor(off_k_B), fill = pct_reject_H0)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_reject_H0, 0)), size = 2, color = "white") +
  scale_fill_viridis_c(limits = c(0, 100), option = "plasma") +
  facet_nested(
    rows = vars(fct_rev(factor(off_R_A))),
    cols = vars(epidemic_size, method, off_R_B)
  ) +
  labs(
    title = "% Rejecting H0 by Parameter Combinations",
    x = "off_k_A",
    y = "off_k_B",
    fill = "% Reject H0"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "right"
  )
library(tidyverse)
library(ggh4x)

# Step 1: Create symmetric dataset
df_symmetric <- bind_rows(
  df,
  df |>
    rename(
      forest_id_A = forest_id_B,
      forest_id_B = forest_id_A,
      param_id_A = param_id_B,
      param_id_B = param_id_A,
      off_R_A = off_R_B,
      off_R_B = off_R_A,
      off_k_A = off_k_B,
      off_k_B = off_k_A
    )
) |>
  distinct()

# Step 2: Compute rejection rates
rejection_data <- df_symmetric |>
  summarise(
    pct_reject_H0 = mean(p_value < 0.05) * 100,
    n_tests = n(),
    .by = c(
      off_R_A,
      off_k_A,
      off_R_B,
      off_k_B,
      epidemic_size,
      forest_size,
      method
    )
  )

# Step 3: Define k ordering
k_levels <- sort(unique(c(rejection_data$off_k_A, rejection_data$off_k_B)))

# Step 4: Filter to upper triangle within each facet
# Upper triangle means: off_k_A <= off_k_B (when both are in the same order)
rejection_upper <- rejection_data |>
  mutate(
    k_A_idx = match(off_k_A, k_levels),
    k_B_idx = match(off_k_B, k_levels)
  ) |>
  filter(k_A_idx <= k_B_idx) |>
  select(-k_A_idx, -k_B_idx)

# Step 5: Plot
ggplot(
  rejection_upper,
  aes(
    x = factor(off_k_A, levels = k_levels),
    y = factor(off_k_B, levels = k_levels),
    fill = pct_reject_H0
  )
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_reject_H0, 0)), size = 2.5, color = "white") +
  scale_fill_viridis_c(limits = c(0, 100), option = "plasma") +
  facet_nested(
    rows = vars(factor(off_R_A)),
    cols = vars(epidemic_size, method, off_R_B)
  ) +
  labs(
    title = "% Rejecting H0 by Parameter Combinations (Upper Triangle)",
    x = "off_k_A",
    y = "off_k_B",
    fill = "% Reject H0"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text = element_text(size = 7),
    legend.position = "right"
  )


# ------------------------------------
#           Good one
# ------------------------------------
library(tidyverse)
library(ggh4x)

# Step 1: Create symmetric dataset
df_symmetric <- bind_rows(
  df,
  df |>
    rename(
      forest_id_A = forest_id_B,
      forest_id_B = forest_id_A,
      param_id_A = param_id_B,
      param_id_B = param_id_A,
      off_R_A = off_R_B,
      off_R_B = off_R_A,
      off_k_A = off_k_B,
      off_k_B = off_k_A
    )
) |>
  distinct() |>
  filter(forest_size == 100)

# Step 2: Compute rejection rates
rejection_data <- df_symmetric |>
  summarise(
    pct_reject_H0 = mean(p_value < 0.05) * 100,
    n_tests = n(),
    .by = c(
      off_R_A,
      off_k_A,
      off_R_B,
      off_k_B,
      epidemic_size,
      forest_size,
      method
    )
  )

# Step 3: Define k ordering
k_levels <- sort(unique(c(rejection_data$off_k_A, rejection_data$off_k_B)))

# Step 4: Assign each combination to correct position
plot_data <- rejection_data |>
  mutate(
    k_A_idx = match(off_k_A, k_levels),
    k_B_idx = match(off_k_B, k_levels)
  ) |>
  filter(
    # Upper triangle: keep chisq where k_A < k_B
    (k_A_idx < k_B_idx & method == "chisq") |
      # Lower triangle: keep permanova where k_A > k_B
      (k_A_idx > k_B_idx & method == "permanova") |
      # Diagonal: keep both methods
      (k_A_idx == k_B_idx)
  )

# Step 5: For diagonal, create split tile data
diagonal_data <- plot_data |>
  filter(k_A_idx == k_B_idx) |>
  select(-k_A_idx, -k_B_idx)

off_diagonal_data <- plot_data |>
  filter(k_A_idx != k_B_idx) |>
  select(-k_A_idx, -k_B_idx)

# Step 6: Plot
p <- ggplot() +
  # Off-diagonal tiles (full squares)
  geom_tile(
    data = off_diagonal_data,
    aes(
      x = factor(off_k_A, levels = k_levels),
      y = factor(off_k_B, levels = k_levels),
      fill = pct_reject_H0
    ),
    color = "white",
    linewidth = 0.5
  ) +
  # geom_text(
  #   data = off_diagonal_data,
  #   aes(
  #     x = factor(off_k_A, levels = k_levels),
  #     y = factor(off_k_B, levels = k_levels),
  #     label = round(pct_reject_H0, 0)
  #   ),
  #   size = 2.5,
  #   color = "white"
  # ) +
  # Diagonal tiles - split triangles
  geom_polygon(
    data = diagonal_data |>
      mutate(
        x_center = as.numeric(factor(off_k_A, levels = k_levels)),
        y_center = as.numeric(factor(off_k_B, levels = k_levels))
      ) |>
      rowwise() |>
      mutate(
        vertices = list(
          if (method == "chisq") {
            # Upper-left triangle
            tibble(
              x = c(
                x_center - 0.5,
                x_center + 0.5,
                x_center - 0.5,
                x_center - 0.5
              ),
              y = c(
                y_center + 0.5,
                y_center + 0.5,
                y_center - 0.5,
                y_center + 0.5
              )
            )
          } else {
            # Lower-right triangle
            tibble(
              x = c(
                x_center + 0.5,
                x_center + 0.5,
                x_center - 0.5,
                x_center + 0.5
              ),
              y = c(
                y_center + 0.5,
                y_center - 0.5,
                y_center - 0.5,
                y_center + 0.5
              )
            )
          }
        )
      ) |>
      ungroup() |>
      unnest(vertices),
    aes(
      x = x,
      y = y,
      fill = pct_reject_H0,
      group = interaction(
        off_k_A,
        off_k_B,
        method,
        off_R_A,
        off_R_B,
        epidemic_size
      )
    ),
    color = "white",
    linewidth = 0.5
  ) +
  # geom_text(
  #   data = diagonal_data |>
  #     mutate(
  #       x_pos = as.numeric(factor(off_k_A, levels = k_levels)),
  #       y_pos = as.numeric(factor(off_k_B, levels = k_levels)),
  #       # Offset text position for split tiles
  #       x_pos = if_else(method == "chisq", x_pos - 0.2, x_pos + 0.2),
  #       y_pos = if_else(method == "chisq", y_pos + 0.2, y_pos - 0.2)
  #     ),
  #   aes(x = x_pos, y = y_pos, label = round(pct_reject_H0, 0)),
  #   size = 2,
  #   color = "white"
  # ) +
  scale_fill_viridis_c(limits = c(0, 100), option = "plasma") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  facet_nested(
    rows = vars(fct_rev(factor(off_R_A))),
    cols = vars(epidemic_size, off_R_B)
  ) +
  labs(
    title = "% Rejecting H0 by Parameter Combinations",
    subtitle = "Upper triangle: Chi-square | Lower triangle: PERMANOVA | Diagonal: Split (upper=Chi-square, lower=PERMANOVA)",
    x = "off_k_A",
    y = "off_k_B",
    fill = "% Reject H0"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text = element_text(size = 7),
    legend.position = "right"
  )
p
ggsave(
  "figures/claude.svg",
  width = 12,
  height = 12,
)

# ------------------------------------
#           roc auc
# ------------------------------------
library(tidyverse)
library(pROC)
library(yardstick)

# Step 1: Prepare data with truth labels
roc_data <- df |>
  mutate(
    # Truth: are parameters truly different?
    truth = factor(
      if_else(param_id_A == param_id_B, "same", "different"),
      levels = c("same", "different") # "different" is the positive class
    ),
    # Convert p-value to a score where higher = more likely different
    # Use 1 - p_value so higher scores indicate stronger evidence of difference
    score = 1 - p_value
  )

# Step 2: Compute ROC curves and AUC by relevant groupings
roc_results <- roc_data |>
  nest(.by = c(method, epidemic_size, forest_size)) |>
  mutate(
    # Compute ROC curve
    roc_obj = map(
      data,
      ~ {
        roc(
          .x$truth,
          .x$score,
          levels = c("same", "different"),
          direction = "<"
        )
      }
    ),
    # Extract AUC with confidence interval
    auc_value = map_dbl(roc_obj, ~ as.numeric(auc(.x))),
    auc_ci = map(roc_obj, ~ ci.auc(.x, conf.level = 0.95)),
    auc_lower = map_dbl(auc_ci, ~ .x[1]),
    auc_upper = map_dbl(auc_ci, ~ .x[3]),
    # Extract ROC curve coordinates for plotting
    roc_coords = map(
      roc_obj,
      ~ {
        tibble(
          threshold = .x$thresholds,
          sensitivity = .x$sensitivities,
          specificity = .x$specificities,
          fpr = 1 - .x$specificities
        )
      }
    ),
    # Sample size info
    n_same = map_dbl(data, ~ sum(.x$truth == "same")),
    n_different = map_dbl(data, ~ sum(.x$truth == "different")),
    n_total = map_dbl(data, nrow)
  )

# Step 3: Alternative using yardstick (tidymodels approach)
roc_results_tidy <- roc_data |>
  summarise(
    roc_curve = list(roc_curve(
      cur_data(),
      truth = truth,
      score,
      event_level = "second"
    )),
    auc_estimate = roc_auc(
      cur_data(),
      truth = truth,
      score,
      event_level = "second"
    )$.estimate,
    .by = c(method, epidemic_size, forest_size)
  )

# Step 4: Summary table
auc_summary <- roc_results |>
  select(
    method,
    epidemic_size,
    forest_size,
    auc_value,
    auc_lower,
    auc_upper,
    n_same,
    n_different
  ) |>
  arrange(method, epidemic_size, forest_size)

print(auc_summary)

# Step 5: Statistical comparison of AUCs between methods
# Use DeLong's test for paired ROC curves (same data, different methods)
auc_comparisons <- roc_data |>
  pivot_wider(
    id_cols = c(
      tree_id,
      forest_id_A,
      forest_id_B,
      epidemic_size,
      forest_size,
      truth
    ),
    names_from = method,
    values_from = score
  ) |>
  nest(.by = c(epidemic_size, forest_size)) |>
  mutate(
    # DeLong's test comparing the two methods
    roc_chisq = map(
      data,
      ~ roc(.x$truth, .x$chisq, levels = c("same", "different"))
    ),
    roc_permanova = map(
      data,
      ~ roc(.x$truth, .x$permanova, levels = c("same", "different"))
    ),
    comparison = map2(
      roc_chisq,
      roc_permanova,
      ~ {
        test_result <- roc.test(.x, .y, method = "delong")
        tibble(
          auc_chisq = as.numeric(auc(.x)),
          auc_permanova = as.numeric(auc(.y)),
          auc_diff = auc_chisq - auc_permanova,
          p_value = test_result$p.value,
          statistic = as.numeric(test_result$statistic)
        )
      }
    )
  ) |>
  select(epidemic_size, forest_size, comparison) |>
  unnest(comparison)

print(auc_comparisons)
