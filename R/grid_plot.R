source("R/packages.R")
source("R/helpers.R")

# ------------------------------------
#           plot_grid
# ------------------------------------
#' @title plot_grid
#' @description
#' Generates a heatmap of simulation results for a given `forest_size`.
#' The x and y axes represent the offspring dispersion parameter `k` for forests A and B.
#' Facets show combinations of reproduction numbers (`R0_A`, `R0_B`) and `epidemic_size`.
#' Each tile shows the proportion of simulations rejecting the null hypothesis
#' (H₀: both forests come from the same distribution) at significance level α = 0.05.
#' χ² test results are shown in the upper triangle, PERMANOVA in the lower triangle.
#' The diagonal (red line) represents identical parameter settings, where H₀ should not be rejected.
#'
#' @param rejection_summary A tibble with columns `method`, `forest_size`, `epidemic_size`,
#'   `off_R_A`, `off_R_B`, `off_k_A`, `off_k_B`, and `prop_reject_H0`.
#' @param forest_size Integer specifying the forest size to plot.
#' @param base_size Base font size for the plot theme. Default is 15.
#'
#' @return A ggplot heatmap.
plot_grid <- function(rejection_summary, forest_size = 20L, base_size = 15) {
  k_levels <- sort(unique(c(
    rejection_summary$off_k_A,
    rejection_summary$off_k_B
  )))
  R_levels <- sort(unique(c(
    rejection_summary$off_R_A,
    rejection_summary$off_R_B
  )))

  symmetric_data <- rejection_summary |>
    filter(forest_size == !!forest_size) |>
    # Create a sorted pair key for matching symmetric pairs
    mutate(
      key = pmap_chr(
        list(method, epidemic_size, off_R_A, off_R_B, off_k_A, off_k_B),
        function(m, e, ra, rb, ka, kb) {
          r_vals <- sort(c(ra, rb))
          k_vals <- sort(c(ka, kb))
          paste(m, e, r_vals[1], r_vals[2], k_vals[1], k_vals[2], sep = "_")
        }
      )
    ) |>
    # Average prop_reject_H0 for symmetric pairs
    group_by(key) |>
    mutate(prop_reject_H0 = mean(prop_reject_H0)) |>
    ungroup() |>
    select(-key)

  plot_df <- bind_rows(
    symmetric_data,
    symmetric_data |>
      mutate(
        temp_R = off_R_A,
        off_R_A = off_R_B,
        off_R_B = temp_R,
        temp_k = off_k_A,
        off_k_A = off_k_B,
        off_k_B = temp_k
      ) |>
      select(-temp_R, -temp_k)
  ) |>
    distinct() |>
    mutate(
      k_A_idx = match(off_k_A, k_levels),
      k_B_idx = match(off_k_B, k_levels),
      R_A_label = paste0("R₀=", off_R_A),
      R_B_label = paste0("R₀=", off_R_B),
      R_A_label = factor(R_A_label, levels = paste0("R₀=", R_levels)),
      R_B_label = factor(R_B_label, levels = paste0("R₀=", rev(R_levels))),
      epidemic_label = paste0("Epidemic size ε=", epidemic_size),
      epidemic_label = factor(
        epidemic_label,
        levels = paste0("Epidemic size ε=", sort(unique(epidemic_size)))
      )
    ) |>
    filter(
      (k_A_idx < k_B_idx & method == "chisq") |
        (k_A_idx > k_B_idx & method == "permanova") |
        (k_A_idx == k_B_idx)
    )

  diag_df <- filter(plot_df, k_A_idx == k_B_idx)
  off_diag_df <- filter(plot_df, k_A_idx != k_B_idx)

  make_triangle <- function(x, y, method) {
    if (method == "chisq") {
      tibble(x = c(x - 0.5, x + 0.5, x - 0.5), y = c(y + 0.5, y + 0.5, y - 0.5))
    } else {
      tibble(x = c(x + 0.5, x + 0.5, x - 0.5), y = c(y + 0.5, y - 0.5, y - 0.5))
    }
  }

  diag_polys <- diag_df |>
    mutate(
      x_c = match(off_k_A, k_levels),
      y_c = match(off_k_B, k_levels),
      vertices = pmap(list(x_c, y_c, method), make_triangle)
    ) |>
    unnest(vertices)

  ggplot() +
    geom_tile(
      data = off_diag_df,
      aes(
        factor(off_k_A, k_levels),
        factor(off_k_B, k_levels),
        fill = prop_reject_H0
      ),
      colour = "white",
      linewidth = 0.1
    ) +
    geom_polygon(
      data = diag_polys,
      aes(
        x,
        y,
        fill = prop_reject_H0,
        group = interaction(
          off_k_A,
          off_k_B,
          off_R_A,
          off_R_B,
          method,
          epidemic_size
        )
      ),
      colour = "white",
      linewidth = 0.1
    ) +
    scale_fill_gradientn(
      colours = c("#f0f0f0", "#bdbdbd", "#8c6bb1", "#810f7c", "#4d004b"),
      values = scales::rescale(c(0, 0.05, 0.1, 0.5, 1)),
      limits = c(0, 1),
      breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
      labels = c(
        "0 ",
        "5\n\U1D6FC",
        "25",
        "50",
        "75",
        "100"
      ),
      name = "\UFE6A rejecting H₀"
    ) +
    geom_point(
      data = tibble(method = c("chisq", "permanova")),
      aes(x = 1, y = 1, shape = method),
      size = 0,
      alpha = 0
    ) +
    scale_shape_manual(
      "",
      values = c("chisq" = 16, "permanova" = 17),
      labels = c("chisq" = "\u25E4 χ² test", "permanova" = "\u25E2 PERMANOVA"),
      guide = guide_legend(
        position = "top",
        override.aes = list(size = 0, alpha = 0),
        keywidth = unit(0, "pt"),
        keyheight = unit(0, "pt"),
        default.unit = "pt"
      )
    ) +
    geom_abline(
      data = diag_df |> filter(off_k_A == off_k_B & off_R_A == off_R_B),
      aes(slope = 1, intercept = 0, linetype = "H0"),
      color = "red",
    ) +
    scale_linetype(
      "",
      breaks = c("H0"),
      labels = c("H₀ true"),
      guide = guide_legend(
        position = "top",
        override.aes = list(size = 0.1),
        keywidth = unit(0, "pt"),
        keyheight = unit(0, "pt"),
        default.unit = "pt"
      )
    ) +
    facet_nested(
      rows = vars(R_B_label),
      cols = vars(epidemic_label, R_A_label)
    ) +
    coord_fixed() +
    labs(
      subtitle = paste("Forest size (\U1D4C2) = ", forest_size),
      x = "\U1D458 (\U1D4D5 A)",
      y = "\U1D458 (\U1D4D5 B)"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 15,
        barheight = 0.7,
        frame.linewidth = 0.5,
        ticks.linewidth = 0.5,
        ticks.colour = "black"
      )
    )
}

# ------------------------------------
#           Plot
# ------------------------------------

results_grid <- readRDS("data/results_grid.rds")
rejection_summary <- results_grid |>
  mutate(
    reject_H0 = p_value < 0.05
  ) |>
  summarise(
    prop_reject_H0 = mean(reject_H0),
    .by = c(
      "method",
      "forest_size",
      "epidemic_size",
      "off_R_A",
      "off_R_B",
      "off_k_A",
      "off_k_B"
    )
  )


map(
  c(20L, 50L, 100L, 200L),
  ~ {
    p <- plot_grid(
      rejection_summary,
      forest_size = .x,
      base_size = 10
    )
    ggsave(
      filename = paste0("figures/grid", .x, ".png"),
      plot = p,
      width = 7.5,
      height = 7.5,
      units = "in",
      dpi = 500
    )
  }
)
