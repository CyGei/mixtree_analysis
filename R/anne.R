# This scripts analyses the results of the simulation study `results_grid.rds`.
source("R/packages.R")
source("R/plots.R")
#results_raw <- readRDS("data/results_raw.rds")

results_grid <- readRDS("data/results_grid.rds") |>
  as_tibble() |>
  dplyr::select(-c(starts_with("gt_"), duration, replicates))
results_grid

df <- results_grid |>
  filter(tree_id_A == tree_id_B)

plot.ROC(results_grid = df, alpha = seq(0, 1, by = 0.001))


offspring <- results_grid |>
  filter(param_id_A == 1 & param_id_B == 1 & tree_id_A == 1 & tree_id_B == 1)

FA <- forest_grid |>
  filter(param_id == 1 & tree_id == 1) |>
  pull(forest) |>
  pluck(1) |>
  as_forest()
FA

FB <- forest_grid |>
  filter(param_id == 1 & tree_id == 2) |>
  pull(forest) |>
  pluck(1) |>
  as_forest()


countR <- function(tree) {
  tree |>
    mutate(from = factor(from, levels = 1:20)) |>
    table() |>
    as.data.frame() |>
    group_by(from) |>
    summarise(R = sum(Freq))
}

off_A <- map(FA, countR) |>
  bind_rows(.id = "forest_id")

off_B <- map(FB, countR) |>
  bind_rows(.id = "forest_id")
nrow(off_B)
nrow(off_A)
off <- bind_rows(
  off_A |> mutate(forest = "A"),
  off_B |> mutate(forest = "B")
)

off |>
  ggplot(aes(x = R, fill = forest)) +
  geom_histogram(position = "identity", alpha = 0.5)

ks.test(off_A$R, off_B$R)

param_grid <- readRDS("data/param_grid.rds")


# Load and prepare base data
df <- readRDS("data/results_grid.rds") |>
  filter(tree_id_A == tree_id_B) |>
  filter(forest_size == 100) |>
  filter(method == "permanova") |>
  mutate(
    across(-p_value, ~ factor(.x, levels = sort(unique(.x)))),
    significance = if_else(p_value < 0.05, "Significant", "Non-significant")
  ) |>
  as_tibble()

df
#compute % significant
df1 <- df |>
  group_by(
    epidemic_size,
    off_k_A,
    off_k_B,
    off_R_A,
    off_R_B,
    param_id_A,
    param_id_B
  ) |>
  summarise(
    prop_significant = mean(significance == "Significant")
  )
