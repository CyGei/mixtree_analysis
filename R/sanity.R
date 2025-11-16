# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log = "log", pipetime.unit = "min")


# ------------------------------------
#           Check comparisons
# ------------------------------------
#' @description
#' We check that, for each tree, we have all pairwise forest comparisons.
#' For n forests, there should be n * (n + 1) / 2 pairs (including self-comparisons).
#' Output should be empty if all trees have complete comparisons.
readRDS("data/results_grid.rds") |>
  summarise(
    n_forests = n_distinct(c(forest_id_A, forest_id_B)),
    n_pairs = n_distinct(paste(forest_id_A, forest_id_B)),
    expected_pairs = n_forests * (n_forests + 1) / 2,
    match = n_pairs == expected_pairs,
    .by = tree_id
  ) |>
  filter(!match)

# ------------------------------------
#           Check trees
# ------------------------------------
#' @description
#' We check that the reference trees used to build the forests
#' have an offspring distribution that matches the simulation parameters.

trees <- readRDS("data/tree_grid.rds") |>
  #tree_id = <param_id>_<replicate>
  separate_wider_delim(
    tree_id,
    delim = "_",
    names = c("param_id", "replicate"),
    cols_remove = FALSE
  ) |>
  select(-replicate) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(param_id, params) |>
      mutate(param_id = as.character(param_id)),
    by = "param_id"
  ) |>
  mutate(
    R = map_dbl(params, ~ .x$offspring_dist$parameters$R),
    k = map_dbl(params, ~ .x$offspring_dist$parameters$k),
    epidemic_size = map_dbl(params, ~ .x$epidemic_size)
  )
group_vars <- c("epidemic_size", "param_id", "tree_id", "R", "k")

Ri <- trees |>
  group_by(across(all_of(group_vars))) |>
  group_modify(
    ~ {
      max_gt <- .x$params |>
        pluck(1) |>
        pluck("generation_time") |>
        (\(gt) gt$q(0.99))()
      tree <- .x$tree |> pluck(1) |> filter(date <= max(date) - max_gt)
      if (nrow(tree) == 0) {
        tibble(R_est = NA, k_est = NA, R_se = NA, k_se = NA)
      } else {
        Ri <- tree |>
          mutate(
            from = factor(from, levels = to),
            to = factor(to, levels = to)
          ) |>
          select(from, to) |>
          table() |>
          as.data.frame() |>
          group_by(from) |>
          summarise(Ri = sum(Freq), .groups = "drop") |>
          pull(Ri)

        if (length(Ri) < 2) {
          tibble(R_est = NA, k_est = NA, R_se = NA, k_se = NA)
        } else {
          fit <- fitdistrplus::fitdist(Ri, "nbinom")
          tibble(
            R_est = fit$estimate["mu"],
            k_est = fit$estimate["size"],
            R_se = fit$sd["mu"],
            k_se = fit$sd["size"]
          )
        }
      }
    }
  )

# Plot per-tree estimates vs true parameters
params <- Ri |>
  dplyr::select(epidemic_size, R, k) |>
  distinct()

Ri |>
  mutate(
    col_label = "Epidemic Size",
    R_error = R_est - R,
    log_k_error = log(k_est) - log(k)
  ) |>
  ggplot() +
  ggh4x::facet_nested(
    rows = vars(R, k),
    cols = vars(col_label, epidemic_size),
    labeller = labeller(
      R = function(x) paste0("R = ", x),
      k = function(x) paste0("k = ", x)
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
  geom_point(aes(x = R_error, y = log_k_error)) +
  coord_cartesian() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )


# ------------------------------------
#           Check Forests
# ------------------------------------
#' @title entropy_grid
#' @description
#' Calculate entropy of each forest in the forest grid
#' and plot distribution of entropy values faceted by simulation parameters.
#' Entropy is calculated using the `o2ools::get_entropy` function.
#' Entropy provides a measure of uncertainty in ancestry assignments within each forest.
#' Higher entropy values indicate greater uncertainty, while lower values suggest more definitive ancestry assignments.

forest_grid <- qs::qread("data/forest_grid.qs")
forest_grid |>
  mutate(
    entropy = map_dbl(
      forest,
      ~ {
        out <- as.data.frame(.x, stringsAsFactors = FALSE) |>
          mutate(across(everything(), ~ as.character(.x)))
        names(out) <- paste0("alpha_", as.character(2:(ncol(out) + 1)))
        class(out) <- c("outbreaker_chains", class(out))
        o2ools::get_entropy(out) |> mean()
      }
    )
  ) |>
  select(forest_id, entropy) |>
  saveRDS("data/entropy_grid.rds")

readRDS("data/entropy_grid.rds") |>
  separate_wider_delim(
    forest_id,
    delim = "_",
    names = c("tree_param_id", "tree_replicate", "param_id"),
    cols_remove = FALSE
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(-c(params, starts_with("gt_"), duration, replicates)) |>
      mutate(param_id = as.character(param_id)),
    by = "param_id"
  ) |>
  mutate(epidemic_size_label = "Epidemic size", off_k = format_k(off_k)) |>
  ggplot() +
  ggh4x::facet_nested(
    cols = vars(epidemic_size_label, epidemic_size),
    rows = vars(off_R, off_k),
    labeller = labeller(
      off_R = function(x) paste0("R₀=", x),
      off_k = function(x) paste0("\U1D458=", x)
    ),
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      text_y = list(
        element_text(size = 16),
        element_text(size = 11)
      ),
      by_layer_x = TRUE,
      by_layer_y = TRUE
    )
  ) +
  geom_histogram(
    aes(x = entropy),
    binwidth = 0.01,
    fill = "blue",
    linewidth = 0.2
  ) +
  scale_y_continuous(breaks = c(0, 400, 800)) +
  scale_x_continuous(
    limits = c(0.4, 1),
    breaks = seq(0.5, 1, by = 0.25)
  ) +
  theme_mixtree() +
  theme(
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(
    x = "Entropy",
    y = "Count"
  )

ggsave(
  "figures/entropy.png",
  width = 7.5,
  height = 8.5,
  units = "in",
  dpi = 300
)


# create a regression model for entropy vs parameters
model_df <- readRDS("data/entropy_grid.rds") |>
  separate_wider_delim(
    forest_id,
    delim = "_",
    names = c("tree_param_id", "tree_replicate", "param_id"),
    cols_remove = FALSE
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(-c(params, starts_with("gt_"), duration, replicates)) |>
      mutate(param_id = as.character(param_id)),
    by = "param_id"
  ) |>
  mutate(
    across(
      c(epidemic_size, off_R, off_k),
      ~ factor(.x, levels = sort(unique(.x)))
    )
  )

m1 <- lm(
  entropy ~ epidemic_size + off_R + off_k,
  data = model_df
)
summary(m1)
m2 <- lm(
  entropy ~ epidemic_size * off_R * off_k,
  data = model_df
)
summary(m2)

# check for multicollinearity
#  forest_size_cat +
#     epidemic_size_cat +
#     method:epidemic_size_cat +
#     delta_R0 * delta_k_cat,
