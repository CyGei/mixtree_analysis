plot_reject_accept_grid <- function(df, forest_size) {
  df <- df |> filter(forest_size == !!forest_size)

  df_all <- df %>%
    mutate(
      x = as.numeric(off_k_A),
      y = as.numeric(off_k_B),
      x_plot = ifelse(method == "chisq", pmin(x, y), pmax(x, y)),
      y_plot = ifelse(method == "chisq", pmax(x, y), pmin(x, y)),
      prop_accept_H0 = 1 - prop_reject_H0,
      is_true_diag = (off_R_A == off_R_B & off_k_A == off_k_B & x == y)
    )

  offdiag_df <- df_all %>% filter(x != y)
  diag_df <- df_all %>% filter(x == y)

  make_half_tile <- function(
    x,
    y,
    method,
    fill_accept,
    fill_reject,
    is_true_diag,
    meta
  ) {
    fill_val <- if (is_true_diag) fill_accept else fill_reject
    fill_type <- if (is_true_diag) "accept" else "reject"
    coords <- if (method == "chisq") {
      tibble(
        x_poly = c(x - 0.5, x + 0.5, x - 0.5),
        y_poly = c(y + 0.5, y + 0.5, y - 0.5)
      )
    } else {
      tibble(
        x_poly = c(x + 0.5, x + 0.5, x - 0.5),
        y_poly = c(y - 0.5, y + 0.5, y - 0.5)
      )
    }
    coords %>%
      mutate(method, x_plot = x, y_plot = y, fill_val, fill_type) %>%
      bind_cols(meta)
  }

  diag_polys <- diag_df %>%
    rowwise() %>%
    do(make_half_tile(
      .$x_plot,
      .$y_plot,
      .$method,
      .$prop_accept_H0,
      .$prop_reject_H0,
      .$is_true_diag,
      tibble(
        epidemic_size = .$epidemic_size,
        off_R_A = .$off_R_A,
        off_R_B = .$off_R_B,
        off_k_A = .$off_k_A,
        off_k_B = .$off_k_B
      )
    )) %>%
    ungroup()

  diag_points <- diag_polys %>%
    group_by(
      x_plot,
      y_plot,
      method,
      fill_val,
      fill_type,
      epidemic_size,
      off_R_A,
      off_R_B
    ) %>%
    summarise(x_point = mean(x_poly), y_point = mean(y_poly), .groups = "drop")

  # build diagonal line segments only for true diagonals
  diag_lines <- diag_polys %>%
    filter(fill_type == "accept") %>%
    distinct(epidemic_size, off_R_A, off_R_B) %>%
    mutate(x = 0.5, xend = 3.5, y = 0.5, yend = 3.5)

  ggplot() +
    geom_tile(
      data = offdiag_df,
      aes(x = x_plot, y = y_plot, fill = prop_reject_H0),
      color = "white",
      linewidth = 0.5
    ) +
    geom_point(
      data = offdiag_df,
      aes(x = x_plot, y = y_plot, shape = method),
      color = "black"
    ) +
    geom_polygon(
      data = diag_polys %>% filter(fill_type == "reject"),
      aes(
        x = x_poly,
        y = y_poly,
        group = interaction(x_plot, y_plot, method),
        fill = fill_val
      ),
      color = "white",
      linewidth = 0.5
    ) +
    geom_point(
      data = diag_points %>% filter(fill_type == "reject"),
      aes(x = x_point, y = y_point, shape = method),
      color = "black"
    ) +
    scale_fill_viridis_c(
      name = expression("Reject " * H[0] * " (%)"),
      option = "inferno",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 10,
        barheight = 0.5,
        frame.linewidth = 0.5,
        ticks.linewidth = 0.5,
        ticks.colour = "black"
      )
    ) +
    new_scale_fill() +
    geom_polygon(
      data = diag_polys %>% filter(fill_type == "accept"),
      aes(
        x = x_poly,
        y = y_poly,
        group = interaction(x_plot, y_plot, method),
        fill = fill_val
      ),
      color = "white",
      linewidth = 0.5
    ) +
    geom_point(
      data = diag_points %>% filter(fill_type == "accept"),
      aes(x = x_point, y = y_point, shape = method),
      color = "black"
    ) +
    geom_segment(
      data = diag_lines,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "black"
    ) +
    scale_fill_viridis_c(
      name = expression("Accept " * H[0] * " (%)"),
      option = "mako",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 10,
        barheight = 0.5,
        frame.linewidth = 0.5,
        ticks.linewidth = 0.5,
        ticks.colour = "black"
      )
    ) +
    facet_nested(
      rows = vars(factor(off_R_A, levels = rev(levels(df$off_R_A)))),
      cols = vars(epidemic_size, off_R_B)
    ) +
    scale_x_continuous(breaks = 1:3, labels = levels(df$off_k_A)) +
    scale_y_continuous(breaks = 1:3, labels = levels(df$off_k_B)) +
    coord_fixed() +
    labs(x = expression(k[A]), y = expression(k[B])) +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_blank(), legend.position = "bottom") +
    guides(
      shape = guide_legend(
        position = "top",
        title = "Method",
        override.aes = list(size = 4)
      )
    )
}

map(
  c(20L, 50L, 100L, 200L),
  ~ {
    p <- plot_reject_accept_grid(df, forest_size = .x)
    ggsave(
      filename = paste0("figures/grid", .x, ".svg"),
      plot = p,
      width = 12,
      height = 12
    )
  }
)
