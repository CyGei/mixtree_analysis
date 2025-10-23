# =================================================
# Set-up
# =================================================
source("R/packages.R")
source("R/helpers.R")
set.seed(000)
plan(multisession, workers = availableCores() - 1)
options(pipetime.log = "log", pipetime.unit = "min")

# ------------------------------------
#           Check trees
# ------------------------------------
#' @description
#' We check that the reference trees used to build the forests
#' have an offspring distribution that matches the simulation parameters.

trees <- readRDS("data/tree_grid.rds") |>
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
      tree <- .x$tree |> pluck(1) |> filter(date <= max(date) - 3)
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

forest_grid <- readRDS("data/forest_grid.rds") |>
  group_by(param_id, tree_id) |>
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


forest_grid |>
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
      epidemic_size = function(x) paste0("Epidemic Size = ", x),
      off_R = function(x) paste0("R = ", x),
      off_k = function(x) paste0("k = ", x)
    )
  ) +
  geom_histogram(
    aes(x = entropy),
    bins = 30,
    fill = "lightblue",
    colour = "black"
  ) +
  theme_classic() +
  labs(
    title = "Distribution of forest entropy",
    x = "Entropy",
    y = "Count"
  )
