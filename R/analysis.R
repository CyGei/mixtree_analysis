# This scripts analyses the results of the simulation study
# `results_grid.rds`.
source("R/packages.R")
source("R/plots.R")

df <- readRDS("data/results_grid.rds") |>
  filter(tree_id_A == tree_id_B) |>
  as_tibble() |>
  select(-c(starts_with("gt_"), duration, replicates)) |>
  mutate(across(-p_value, ~ factor(.x, levels = sort(unique(.x))))) |>
  group_by(across(-c(p_value, starts_with("tree_id")))) |>
  summarise(prop_reject_H0 = mean(p_value < 0.05), .groups = "drop") |>
  select(
    epidemic_size,
    forest_size,
    off_R_A,
    off_R_B,
    off_k_A,
    off_k_B,
    method,
    prop_reject_H0
  )

df


# example 1 grid
df |>
  filter(
    epidemic_size == 100,
    off_R_A == 1.5,
    off_R_B == 1.5,
  ) |>
  # here move tile to lower triangle for chisq test
  mutate(
    off_k_A_num = as.numeric(off_k_A),
    off_k_B_num = as.numeric(off_k_B),

    # off_k_A_num = case_when(
    #   method == "chisq" & off_k_A_num ...,
    # )
  ) |>
  ggplot(aes(x = off_k_A_num, y = off_k_B_num, fill = prop_reject_H0)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(name = "Prop. reject H0", limits = c(0, 1)) +
  coord_fixed() +
  labs(x = expression(k[A]), y = expression(k[B]))


library(dplyr)
library(ggplot2)

df_sub <- df %>%
  filter(
    epidemic_size == 100,
    off_R_A == 1.5,
    off_R_B == 1.5
  ) %>%
  mutate(
    # numeric positions for plotting
    x = as.numeric(off_k_A),
    y = as.numeric(off_k_B),

    # remap coordinates depending on method
    x_plot = case_when(
      method == "chisq" ~ pmin(x, y),
      method == "permanova" ~ pmax(x, y)
    ),
    y_plot = case_when(
      method == "chisq" ~ pmax(x, y),
      method == "permanova" ~ pmin(x, y)
    )
  )

diag <- df_sub %>%
  filter(x == y) %>%
  rowwise() %>%
  mutate(
    polygon = list(
      if (method == "chisq") {
        tibble(
          x = c(x - 0.5, x + 0.5, x - 0.5),
          y = c(y - 0.5, y - 0.5, y + 0.5)
        )
      } else {
        tibble(
          x = c(x + 0.5, x + 0.5, x - 0.5),
          y = c(y + 0.5, y - 0.5, y + 0.5)
        )
      }
    )
  ) %>%
  select(-x, -y) %>% # drop originals
  unnest(polygon) %>%
  ungroup()


ggplot(df_sub, aes(x = x_plot, y = y_plot, fill = prop_reject_H0)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_point(aes(shape = method), size = 3, color = "black") +
  scale_fill_viridis_c(name = "Prop. reject H0", limits = c(0, 1)) +
  coord_fixed() +
  scale_x_continuous(breaks = 1:3, labels = levels(df$off_k_A)) +
  scale_y_continuous(breaks = 1:3, labels = levels(df$off_k_B)) +
  labs(x = expression(k[A]), y = expression(k[B]))


library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Subset the data
df_sub <- df %>%
  filter(
    epidemic_size == 100,
    off_R_A == 1.5,
    off_R_B == 1.5
  ) %>%
  mutate(
    # numeric positions for plotting
    x = as.numeric(off_k_A),
    y = as.numeric(off_k_B),
    # remap coordinates depending on method
    x_plot = case_when(
      method == "chisq" ~ pmin(x, y),
      method == "permanova" ~ pmax(x, y)
    ),
    y_plot = case_when(
      method == "chisq" ~ pmax(x, y),
      method == "permanova" ~ pmin(x, y)
    )
  )

# 2. Split into off-diagonal and diagonal
offdiag_df <- df_sub %>% filter(x != y)
diag_df <- df_sub %>% filter(x == y)

# 3. Build polygons for diagonal half-tiles
diag_polys <- diag_df %>%
  rowwise() %>%
  mutate(
    polygon = list(
      if (method == "chisq") {
        tibble(
          x_poly = c(x_plot - 0.5, x_plot + 0.5, x_plot - 0.5),
          y_poly = c(y_plot - 0.5, y_plot - 0.5, y_plot + 0.5)
        )
      } else {
        tibble(
          x_poly = c(x_plot + 0.5, x_plot + 0.5, x_plot - 0.5),
          y_poly = c(y_plot + 0.5, y_plot - 0.5, y_plot + 0.5)
        )
      }
    )
  ) %>%
  unnest(polygon) %>%
  ungroup()

# 4. Plot
ggplot() +
  # off-diagonal full tiles
  geom_tile(
    data = offdiag_df,
    aes(x = x_plot, y = y_plot, fill = prop_reject_H0),
    color = "white",
    linewidth = 0.5
  ) +
  geom_point(
    data = offdiag_df,
    aes(x = x_plot, y = y_plot, shape = method),
    size = 3,
    color = "black"
  ) +
  # diagonal half-tiles
  geom_polygon(
    data = diag_polys,
    aes(
      x = x_poly,
      y = y_poly,
      group = interaction(x_plot, y_plot, method),
      fill = prop_reject_H0
    ),
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_viridis_c(name = "Prop. reject H0", limits = c(0, 1)) +
  coord_fixed() +
  scale_x_continuous(breaks = 1:3, labels = levels(df$off_k_A)) +
  scale_y_continuous(breaks = 1:3, labels = levels(df$off_k_B)) +
  labs(x = expression(k[A]), y = expression(k[B]))
# Add point positions for diagonal
# Compute centroids of diagonal polygons
diag_points <- diag_polys %>%
  group_by(x_plot, y_plot, method, prop_reject_H0) %>%
  summarise(
    x_point = mean(x_poly),
    y_point = mean(y_poly),
    .groups = "drop"
  )

ggplot() +
  # off-diagonal full tiles
  geom_tile(
    data = offdiag_df,
    aes(x = x_plot, y = y_plot, fill = prop_reject_H0),
    color = "white",
    linewidth = 0.5
  ) +
  geom_point(
    data = offdiag_df,
    aes(x = x_plot, y = y_plot, shape = method),
    size = 3,
    color = "black"
  ) +
  # diagonal half-tiles
  geom_polygon(
    data = diag_polys,
    aes(
      x = x_poly,
      y = y_poly,
      group = interaction(x_plot, y_plot, method),
      fill = prop_reject_H0
    ),
    color = "white",
    linewidth = 0.5
  ) +
  # diagonal points at polygon centroids
  geom_point(
    data = diag_points,
    aes(x = x_point, y = y_point, shape = method),
    size = 3,
    color = "black"
  ) +
  scale_fill_viridis_c(name = "Prop. reject H0", limits = c(0, 1)) +
  coord_fixed() +
  scale_x_continuous(breaks = 1:3, labels = levels(df$off_k_A)) +
  scale_y_continuous(breaks = 1:3, labels = levels(df$off_k_B)) +
  labs(x = expression(k[A]), y = expression(k[B]))
