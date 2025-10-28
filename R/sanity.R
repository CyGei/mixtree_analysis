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
  mutate(
    param_id = as.integer(str_extract(tree_id, "^[^_]+"))
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(param_id, params),
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
#' @description
#' We check that the forests have enough entropy.
#' Entropy is defined as the variance in ancestries for each case.

forest_grid <- readRDS("data/forest_grid.rds")

entropy_grid <- forest_grid |>
  group_by(forest_id) |>
  group_modify(
    ~ {
      mat <- .x$forest[[1]]
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- paste0("alpha_", colnames(mat))
      class(out) <- c("outbreaker_chains", class(out))
      tibble(entropy = o2ools::get_entropy(out))
    }
  ) |>
  summarise(entropy = mean(entropy))


entropy_grid |>
  #forest_id = <tree_id>_<param_id>
  mutate(
    param_id = as.integer(str_extract(forest_id, "(?<=_)[^_]+$"))
  ) |>
  left_join(
    readRDS("data/param_grid.rds") |>
      select(-c(params, starts_with("gt_"), duration, replicates)),
    by = "param_id"
  ) |>
  ggplot() +
  ggh4x::facet_nested(
    cols = vars(epidemic_size),
    rows = vars(off_R, off_k),
    labeller = labeller(
      off_R = function(x) paste0("R₀=", x),
      off_k = function(x) paste0("\U1D458=", x),
      epidemic_size = function(x) paste0("Epidemic size ε = ", x)
    )
  ) +
  geom_histogram(
    aes(x = entropy),
    # bins = 30,
    fill = "lightblue",
    colour = "black",
    linewidth = 0.2
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  theme_bw() +
  theme(
    axis.line.x = element_line(colour = "black", linewidth = 0.5),
    #remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = NA, colour = "black")
  ) +
  labs(
    x = "Entropy",
    y = "Count"
  )

ggsave("figures/entropy.svg", width = 6, height = 8)
