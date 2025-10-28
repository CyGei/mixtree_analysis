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
#'
plot_grid <- function(rejection_summary, forest_size = 20L, base_size = 15) {
  k_levels <- sort(unique(c(
    rejection_summary$off_k_A,
    rejection_summary$off_k_B
  )))
  R_levels <- sort(unique(c(
    rejection_summary$off_R_A,
    rejection_summary$off_R_B
  )))
  epidemic_levels <- sort(unique(rejection_summary$epidemic_size))
  rejection_summary_sym <-
    bind_rows(
      rejection_summary,
      rejection_summary |>
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
    distinct()
  # Create symmetric data (A vs B and B vs A)
  plot_tab <-
    expand_grid(
      off_k_A = k_levels,
      off_k_B = k_levels,
      off_R_A = R_levels,
      off_R_B = R_levels,
      epidemic_size = epidemic_levels,
      method = c("chisq", "permanova")
    ) |>
    left_join(
      rejection_summary_sym |> filter(forest_size == !!forest_size),
      by = c(
        "off_k_A",
        "off_k_B",
        "off_R_A",
        "off_R_B",
        "epidemic_size",
        "method"
      )
    )

  plot_df <- plot_tab |>
    mutate(
      keep = case_when(
        method == "chisq" & off_k_A <= off_k_B ~ TRUE,
        method == "permanova" & off_k_A >= off_k_B ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    filter(keep) |>
    mutate(
      off_k_A_f = factor(off_k_A, levels = k_levels),
      off_k_B_f = factor(off_k_B, levels = k_levels),
      off_R_A_f = factor(off_R_A, levels = R_levels),
      off_R_B_f = factor(off_R_B, levels = R_levels),
      epidemic_size_f = factor(epidemic_size, levels = epidemic_levels)
    )

  diag_df <- filter(plot_df, off_k_A == off_k_B)
  off_diag_df <- filter(plot_df, off_k_A != off_k_B)

  make_triangle <- function(x, y, method) {
    if (method == "chisq") {
      tibble(x = c(x - 0.5, x + 0.5, x - 0.5), y = c(y + 0.5, y + 0.5, y - 0.5))
    } else {
      tibble(x = c(x + 0.5, x + 0.5, x - 0.5), y = c(y + 0.5, y - 0.5, y - 0.5))
    }
  }

  diag_polys <- diag_df |>
    mutate(
      vertices = pmap(
        list(
          as.numeric(off_k_A_f),
          as.numeric(off_k_B_f),
          method
        ),
        make_triangle
      )
    ) |>
    unnest(vertices)

  ggplot() +
    geom_tile(
      data = off_diag_df,
      aes(
        x = off_k_A_f,
        y = off_k_B_f,
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
          off_k_A_f,
          off_k_B_f,
          off_R_A_f,
          off_R_B_f,
          method,
          epidemic_size_f
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
        keywidth = unit(0, "pt"), # remove key
        keyheight = unit(0, "pt"),
        default.unit = "pt"
      )
    ) +
    geom_abline(
      data = diag_df |> filter(off_k_A_f == off_k_B_f & off_R_A_f == off_R_B_f),
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
        keywidth = unit(0, "pt"), # remove key
        keyheight = unit(0, "pt"),
        default.unit = "pt"
      )
    ) +
    facet_nested(
      rows = vars(off_R_B_f),
      cols = vars(epidemic_size_f, off_R_A_f)
      # rows = vars(factor(off_R_B, rev(R_levels))),
      # cols = vars(epidemic_size, off_R_A),
      # labeller = labeller(
      #   off_R_A = function(x) paste0("R₀=", x),
      #   off_R_B = function(x) paste0("R₀=", x),
      #   epidemic_size = function(x) paste0("Epidemic size \U03B5 = ", x)
      # )
    ) +
    coord_fixed() +
    labs(
      subtitle = paste("Forest size (\U1D4C2) = ", forest_size),
      x = "\U1D458 (\U1D4D5 A)",
      y = "\U1D458  (\U1D4D5 B)"
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
plot_grid(rejection_summary, forest_size = 20L)


#regular tile for 1 method

ggplot() +
  facet_nested(
    rows = vars(off_R_B_f),
    cols = vars(epidemic_size_f, off_R_A_f)
  ) +
  geom_tile(
    data = plot_df |>
      filter(method == "chisq"),
    aes(
      x = off_k_A_f,
      y = off_k_B_f,
      fill = prop_reject_H0
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
  coord_fixed() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
