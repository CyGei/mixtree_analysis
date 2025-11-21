## Setup and Data Loading
library(qs)
library(tidyverse)
library(mixtree)
source("R/helpers.R")

forest_grid <- qs::qread("data/forest_grid.qs") |>
  separate_wider_delim(
    forest_id,
    "_",
    names = c("tree_param_id", "tree_replicate", "param_id")
  ) |>
  mutate(tree_id = paste0(tree_param_id, "_", tree_replicate)) |>
  left_join(
    readRDS("data/param_grid.rds") |> mutate(param_id = as.character(param_id)),
    by = "param_id"
  )
params <- tibble(off_R = c(2, 2, 3), off_k = c(0.1, 0.3, 0.1))

set.seed(123)
trees_selected <- forest_grid |>
  distinct(epidemic_size, tree_id) |>
  slice_sample(n = 1, by = epidemic_size) |>
  pull(tree_id)

forests <- forest_grid |>
  filter(tree_id %in% trees_selected) |>
  inner_join(params, by = c("off_R", "off_k")) |>
  select(epidemic_size, tree_id, param_id, off_R, off_k, forest) |>
  arrange(epidemic_size, off_R, off_k) |>
  mutate(
    forest = map(
      forest,
      ~ .x |>
        as.data.frame() |>
        rename_with(~ as.character(seq_along(.x) + 1), everything()) |>
        mutate(across(everything(), as.character))
    )
  )

dir.create("Py/data/", showWarnings = FALSE, recursive = TRUE)
write.csv(
  select(forests, epidemic_size, tree_id, off_R, off_k),
  "Py/data/metadata.csv",
  row.names = FALSE
)

forests |>
  select(forest, epidemic_size, off_R, off_k, tree_id) |>
  pwalk(function(forest, epidemic_size, off_R, off_k, tree_id) {
    fname <- paste0(
      "Py/data/forest_",
      epidemic_size,
      "_",
      off_R,
      "_",
      off_k,
      "_",
      tree_id,
      ".csv"
    )
    write.csv(forest, file = fname, row.names = FALSE)
  })

## Transformation and Comparison
forests <- forests |>
  mutate(
    forest = map(forest, \(f) as.matrix(f) |> as_forest()),
    id = row_number()
  )

results <- forests |>
  rename_with(~ paste0(., "_A"), -epidemic_size) |>
  inner_join(
    forests |> rename_with(~ paste0(., "_B"), -epidemic_size),
    by = "epidemic_size"
  ) |>
  filter(id_A <= id_B) |>
  cross_join(tibble(sample_size = c(50, 100))) |>
  mutate(
    p_value = pmap_dbl(list(forest_A, forest_B, sample_size), \(fA, fB, s) {
      tree_test(
        fA[sample(length(fA), s)],
        fB[sample(length(fB), s)],
        method = "permanova"
      )$`Pr(>F)`[1]
    })
  ) |>
  select(-forest_A, -forest_B)

write.csv(results, "Py/mixtree_results.csv", row.names = FALSE)


df <-
  bind_rows(
    read.csv("Py/mixtree_results.csv") |>
      mutate(method = "mixtree"),
    read.csv("Py/gtst_results.csv") |>
      mutate(method = "gtst")
  ) |>
  mutate(
    H0 = (off_R_A == off_R_B) & (off_k_A == off_k_B),
  )

df |>
  mutate(
    label = paste0(
      "R:",
      off_R_A,
      " k:",
      off_k_A,
      " vs ",
      "R:",
      off_R_B,
      " k:",
      off_k_B
    ),
    sample_size = factor(sample_size)
  ) |>
  #arrange label by H0 first
  arrange(H0, label) |>
  mutate(
    label = factor(label, levels = unique(label))
  ) |>
  ggplot(aes(x = p_value, y = label)) +
  facet_grid(method ~ epidemic_size, scales = "free_y") +
  geom_vline(
    xintercept = 0.05,
    linetype = "dashed",
    color = "red",
    alpha = 0.6
  ) +

  geom_point(
    aes(color = sample_size, shape = H0),
    size = 3,
    alpha = 0.8,
    position = position_jitter(height = 0.1)
  ) +
  #scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    axis.text.y = element_text(size = 9)
  )
