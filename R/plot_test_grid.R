plot_test_grid <- function(sample_size, method = NULL, half_tiles = FALSE) {
  # Input validation
  stopifnot(
    is.numeric(sample_size) && length(sample_size) == 1,
    is.null(method) || (is.character(method) && length(method) == 1),
    is.logical(half_tiles) && length(half_tiles) == 1
  )

  # Load and prepare base data
  df <- readRDS("data/test_grid.rds") |>
    filter(sample_size == !!sample_size) |>
    mutate(
      across(-p_value, ~ factor(.x, levels = sort(unique(.x)))),
      significance = if_else(p_value < 0.05, "Significant", "Non-significant")
    )

  # Filter by method if specified
  if (!is.null(method)) {
    df <- df |> filter(method == !!method)
  }

  # Color palettes
  epidemic_pal <- c(
    "20" = "#1a8061",
    "50" = "#d95f02",
    "100" = "#7570b3",
    "200" = "#e7298a"
  )

  significance_pal <- c(
    "Significant" = "#f99820",
    "Non-significant" = "#03658c"
  )

  # Create labeller
  grid_labeller <- labeller(
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

  # Prepare data based on plot type
  if (half_tiles) {
    # Create polygon data for half tiles
    plot_df <- df |>
      mutate(
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
        polygon = if (method == "chisq") {
          list(tibble(
            x = c(x - 0.5, x + 0.5, x - 0.5),
            y = c(y - 0.5, y - 0.5, y + 0.5)
          ))
        } else {
          list(tibble(
            x = c(x + 0.5, x + 0.5, x - 0.5),
            y = c(y + 0.5, y - 0.5, y + 0.5)
          ))
        },
        significance = significance
      ) |>
      unnest(polygon)

    # Background data for half tiles
    bg_df <- df |>
      distinct(epidemic_size, off_R_A, off_R_B) |>
      mutate(x = 0, y = 0)

    # Create plot with polygons
    p <- ggplot(plot_df, aes(x = x, y = y)) +
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
      geom_polygon(
        aes(fill = significance, group = interaction(off_k_A, off_k_B, method)),
        color = "white",
        linewidth = 0.25
      ) +
      scale_x_continuous(
        breaks = seq_along(levels(df$off_k_A)),
        labels = levels(df$off_k_A)
      ) +
      scale_y_continuous(
        breaks = seq_along(levels(df$off_k_B)),
        labels = levels(df$off_k_B)
      ) +
      coord_fixed(
        xlim = c(0.5, length(levels(df$off_k_A)) + 0.5),
        ylim = c(0.5, length(levels(df$off_k_B)) + 0.5)
      )
  } else {
    # Create plot with full tiles
    plot_df <- df

    p <- ggplot(plot_df, aes(x = off_k_A, y = off_k_B)) +
      geom_rect(
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = -Inf,
          ymax = Inf,
          fill = epidemic_size
        ),
        alpha = 0.01,
        colour = NA
      ) +
      scale_fill_manual(values = epidemic_pal, guide = FALSE) +
      ggnewscale::new_scale_fill() +
      geom_tile(
        aes(fill = significance),
        color = "white",
        linewidth = 0.25
      ) +
      coord_fixed()
  }

  # Add common elements
  p <- p +
    facet_nested(
      rows = vars(factor(off_R_A, levels = rev(levels(off_R_A)))),
      cols = vars(epidemic_size, off_R_B),
      nest_line = element_line(linetype = 1),
      labeller = grid_labeller
    ) +
    scale_fill_manual(values = significance_pal) +
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

  return(p)
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# Full tiles for one method and sample size
# plot_test_grid(sample_size = 100, method = "chisq")

# Full tiles for all methods and sample size (if method = NULL)
# plot_test_grid(sample_size = 100)

# Half tiles comparing methods for one sample size
# plot_test_grid(sample_size = 100, half_tiles = TRUE)

# plot_test_grid <- function(
#   sample_size,
#   method = NULL,
#   half_tiles = FALSE,
#   data_path = "data/test_grid.rds"
# ) {
#   # -------------------------------------------------
#   # Load and filter data
#   # -------------------------------------------------
#   df <- readRDS(data_path) |>
#     filter(sample_size == !!sample_size) |>
#     mutate(
#       across(
#         c(off_R_A, off_R_B, off_k_A, off_k_B, epidemic_size),
#         ~ factor(.x, levels = sort(unique(.x)))
#       ),
#       significance = if_else(p_value < 0.05, "Significant", "Non-significant"),
#       x = as.numeric(off_k_A),
#       y = as.numeric(off_k_B)
#     )

#   # Filter for single method if full tiles
#   if (!is.null(method) && !half_tiles) {
#     df <- filter(df, method == !!method)
#   }

#   # -------------------------------------------------
#   # Colour palette
#   # -------------------------------------------------
#   epidemic_pal <- c(
#     "20" = "#1a8061",
#     "50" = "#d95f02",
#     "100" = "#7570b3",
#     "200" = "#e7298a"
#   )

#   # -------------------------------------------------
#   # Labeller
#   # -------------------------------------------------
#   grid_labeller <- labeller(
#     .rows = function(x) {
#       vals <- levels(df$off_R_A)
#       paste0(
#         "<span style='color:black; font-size:14pt;'>",
#         ifelse(x == vals[length(vals)], paste0("R0 = ", x), x),
#         "</span>"
#       )
#     },
#     .cols = function(x) {
#       x_char <- as.character(x)
#       paste0(
#         "<span style='color:",
#         ifelse(x_char %in% names(epidemic_pal), epidemic_pal[x_char], "black"),
#         "; font-size:",
#         ifelse(x_char %in% names(epidemic_pal), "14pt", "10pt"),
#         "'>",
#         ifelse(
#           x_char == levels(df$epidemic_size)[1],
#           paste0("Epidemic size = ", x_char),
#           x_char
#         ),
#         "</span>"
#       )
#     }
#   )

#   # -------------------------------------------------
#   # Background shading (full facet tiles)
#   # -------------------------------------------------
#   bg_df <- df |>
#     distinct(epidemic_size, off_R_A, off_R_B) |>
#     mutate(x = 0, y = 0)

#   # -------------------------------------------------
#   # Build ggplot
#   # -------------------------------------------------
#   p <- ggplot(df, aes(x = x, y = y)) +
#     facet_nested(
#       rows = vars(factor(off_R_A, levels = rev(levels(off_R_A)))),
#       cols = vars(epidemic_size, off_R_B),
#       nest_line = element_line(linetype = 1),
#       labeller = grid_labeller
#     ) +
#     geom_tile(
#       data = bg_df,
#       aes(x = x, y = y, fill = epidemic_size),
#       width = 1000,
#       height = 1000,
#       inherit.aes = FALSE,
#       alpha = 0.25
#     ) +
#     scale_fill_manual(values = epidemic_pal, guide = FALSE) +
#     ggnewscale::new_scale_fill()

#   # -------------------------------------------------
#   # Tiles: full or half
#   # -------------------------------------------------
#   if (half_tiles) {
#     df_poly <- df |>
#       group_by(
#         off_k_A,
#         off_k_B,
#         method,
#         significance,
#         epidemic_size,
#         off_R_A,
#         off_R_B
#       ) |>
#       reframe(
#         polygon = list(
#           ifelse(
#             method == "chisq",
#             tibble(
#               x = c(x - 0.5, x + 0.5, x - 0.5),
#               y = c(y - 0.5, y - 0.5, y + 0.5)
#             ),
#             tibble(
#               x = c(x + 0.5, x + 0.5, x - 0.5),
#               y = c(y + 0.5, y - 0.5, y + 0.5)
#             )
#           )
#         ),
#         significance = significance
#       ) |>
#       unnest(polygon)

#     p <- p +
#       geom_polygon(
#         data = df_poly,
#         aes(fill = significance, group = interaction(off_k_A, off_k_B, method)),
#         color = "white",
#         linewidth = 0.25
#       )
#   } else {
#     p <- p +
#       geom_tile(
#         aes(fill = significance),
#         color = "white",
#         linewidth = 0.25
#       )
#   }

#   # -------------------------------------------------
#   # Diagonal line
#   # -------------------------------------------------
#   p <- p +
#     geom_segment(
#       data = df |> filter(off_k_A == off_k_B, off_R_A == off_R_B),
#       aes(
#         x = x - 0.5,
#         y = y - 0.5,
#         xend = x + 0.5,
#         yend = y + 0.5
#       ),
#       color = "black",
#       linewidth = 0.25,
#       inherit.aes = FALSE
#     )

#   # -------------------------------------------------
#   # Final theme and labels
#   # -------------------------------------------------
#   p +
#     scale_fill_manual(
#       values = c("Significant" = "#f99820", "Non-significant" = "#03658c")
#     ) +
#     labs(x = "k (Forest A)", y = "k (Forest B)", fill = "Difference") +
#     theme_light(base_size = 11) +
#     theme(
#       legend.position = "bottom",
#       strip.placement = "outside",
#       panel.spacing = unit(1, "mm"),
#       strip.switch.pad.grid = unit(0, "pt"),
#       strip.switch.pad.wrap = unit(0, "pt"),
#       strip.background = element_blank(),
#       strip.text = element_markdown(face = "bold"),
#       axis.text.x = element_blank()
#     )
# }
