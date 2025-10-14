# This scripts analyses the results of the simulation study `results_grid.rds`.
source("R/packages.R")
source("R/plots.R")
#results_raw <- readRDS("data/results_raw.rds")

results_grid <- readRDS("data/results_grid.rds") |>
  as_tibble() |>
  select(-c(starts_with("gt_"), duration, replicates))
results_grid

plot.ROC()

df_1 <- results_grid |>
  mutate(
    H0 = ifelse(param_id_A == param_id_B, TRUE, FALSE),
    reject_H0 = ifelse(p_value < 0.05, 1, 0),
    R_diff = abs(off_R_A - off_R_B),
    k_diff = abs(off_k_A - off_k_B),
  ) |>
  group_by(
    k_diff,
    R_diff,
    epidemic_size,
    forest_size,
    method
  ) |>
  summarise(
    freq_reject_H0 = mean(reject_H0),
    .groups = "drop"
  )
df_1

ggplot(
  df_1,
  aes(x = R_diff, y = freq_reject_H0, color = method)
) +
  geom_point() +
  geom_line() +
  facet_grid(
    rows = vars(forest_size),
    cols = vars(epidemic_size),
    labeller = label_both
  ) +
  theme_classic()

ggplot(
  df_1,
  aes(x = as.factor(k_diff), y = freq_reject_H0, color = method)
) +
  geom_point() +
  facet_grid(
    rows = vars(forest_size),
    cols = vars(epidemic_size),
    labeller = label_both
  ) +
  theme_classic()
