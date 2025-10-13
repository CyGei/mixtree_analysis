# This scripts analyses the results of the simulation study `results_grid.rds`.

source("R/packages.R")
results_grid <- readRDS("data/results_grid.rds")
head(results_grid)


# Figure 1: absolute difference between off_R_A and off_R_B
df_1 <- results_grid |>
  mutate(R_diff = abs(off_R_A - off_R_B)) |>
  select(
    epidemic_size,
    forest_size,
    off_R_A,
    off_R_B,
    R_diff,
    method,
    p_value
  )

head(df_1)

df_1 |>
  group_by(epidemic_size, forest_size, method, R_diff) |>
  # % of times p_value < 0.05
  summarise(
    n = n(),
    n_signif = sum(p_value < 0.05),
    p_signif = n_signif / n,
    .groups = "drop"
  ) |>
  ggplot(aes(x = R_diff, y = p_signif, color = method)) +
  geom_line() +
  geom_point() +
  facet_grid(
    rows = vars(forest_size),
    cols = vars(epidemic_size),
    labeller = label_both,
    scales = "fixed"
  ) +
  labs(
    title = "Power to detect difference in R",
    x = "|R_A - R_B|",
    y = "Proportion of significant tests (p < 0.05)",
    color = "Method"
  ) +
  tracetheme::theme_trace()
diff(c(1.5, 2.5, 4.5))
