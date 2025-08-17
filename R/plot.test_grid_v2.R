plot.test_grid_half_tile <- function(sample_size) {
  rid_labeller <- labeller(
    .rows = function(x) {
      vals <- levels(df$off_R_A)
      paste0(
        "<span style='color:black; font-size:14pt;'>",
        ifelse(x == vals[length(vals)], paste0("R0 = ", x), x),
        "</span>"
      )
    },
    .cols = function(x) {
      x_char <- as.character(x)
      paste0(
        "<span style='color:",
        ifelse(x_char %in% names(epidemic_pal), epidemic_pal[x_char], "black"),
        "; font-size:",
        ifelse(x_char %in% names(epidemic_pal), "14pt", "10pt"),
        "'>",
        ifelse(
          x_char == levels(df$epidemic_size)[1],
          paste0("Epidemic size = ", x_char),
          x_char
        ),
        "</span>"
      )
    }
  )

  epidemic_pal <- c(
    "20" = "#1a8061", # teal
    "50" = "#d95f02", # orange
    "100" = "#7570b3", # purple
    "200" = "#e7298a" # pink
  )

  df <- readRDS("data/test_grid.rds") |>
    filter(sample_size == sample_size) |>
    mutate(
      across(
        -p_value,
        ~ factor(.x, levels = sort(unique(.x)))
      ),
      significance = if_else(p_value < 0.05, "Significant", "Non-significant"),
      # make polygons for half-tiles
      x = as.numeric(off_k_A),
      y = as.numeric(off_k_B)
    ) |>
    group_by(
      off_k_A,
      off_k_B,
      method,
      significance,
      epidemic_size,
      off_R_A,
      off_R_B
    ) |>
    reframe(
      # bottom-left triangle for chisq
      polygon = if (method == "chisq") {
        list(
          tibble(
            x = c(x - 0.5, x + 0.5, x - 0.5),
            y = c(y - 0.5, y - 0.5, y + 0.5)
          )
        )
      } else {
        list(
          # top-right triangle for permanova
          tibble(
            x = c(x + 0.5, x + 0.5, x - 0.5),
            y = c(y + 0.5, y - 0.5, y + 0.5)
          )
        )
      },
      significance = significance
    ) |>
    unnest(polygon)

  # compute background tiles covering the full panel
  bg_df <- df |>
    distinct(epidemic_size, off_R_A, off_R_B) |> # all facets
    mutate(
      x = 0,
      y = 0,
      width = 1000,
      height = 1000
    )

  p <- ggplot(df, aes(x = x, y = y)) +

    # -------------------------------------------------
    # Facet
    # -------------------------------------------------
    facet_nested(
      rows = vars(factor(off_R_A, levels = rev(levels(off_R_A)))),
      cols = vars(epidemic_size, off_R_B),
      nest_line = element_line(linetype = 1),
      labeller = grid_labeller
    ) +

    # -------------------------------------------------
    # Background shading
    # -------------------------------------------------
    geom_tile(
      data = bg_df,
      aes(x = x, y = y, fill = epidemic_size),
      width = Inf,
      height = Inf,
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    scale_fill_manual(values = epidemic_pal, guide = FALSE) +
    ggnewscale::new_scale_fill() +

    # -------------------------------------------------
    # Half-tiles
    # -------------------------------------------------
    geom_polygon(
      aes(
        fill = significance,
        group = interaction(off_k_A, off_k_B, method)
      ),
      color = "white",
      linewidth = 0.25
    ) +
    scale_fill_manual(
      values = c("Significant" = "#f99820", "Non-significant" = "#03658c")
    ) +
    # -------------------------------------------------
    # Diagonal
    # -------------------------------------------------
    geom_segment(
      data = df |> filter(off_k_A == off_k_B, off_R_A == off_R_B),
      aes(
        x = as.numeric(off_k_A) - 0.5,
        y = as.numeric(off_k_B) - 0.5,
        xend = as.numeric(off_k_A) + 0.5,
        yend = as.numeric(off_k_B) + 0.5
      ),
      color = "black",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +

    # -------------------------------------------------
    # Theme
    # -------------------------------------------------
    labs(x = "k (Forest A)", y = "k (Forest B)", fill = "Difference") +
    theme_light(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.placement = "outside",
      panel.spacing = unit(1, "mm"),
      strip.switch.pad.grid = unit(0, "pt"),
      strip.switch.pad.wrap = unit(0, "pt"),
      strip.background = element_blank(),
      strip.text = element_markdown(face = "bold"),
      axis.text.x = element_blank()
    )

  p
}

g
ggsave(
  "grid_test.png",
  p,
  height = 10,
  width = 12,
  dpi = 300
) |>
  pipetime::time_pipe("saving grid_test.png")
