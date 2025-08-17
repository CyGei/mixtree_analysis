# ==============================================================================
# SIMPLIFIED UNIFIED TEST GRID PLOTTING FUNCTION
# ==============================================================================

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

  # Split data into two groups for vertical stacking
  epidemic_levels <- levels(df$epidemic_size)

  # First subplot: epidemic sizes 1 & 2
  df1 <- df |> filter(epidemic_size %in% epidemic_levels[1:2])

  # Second subplot: epidemic sizes 3 & 4
  df2 <- df |> filter(epidemic_size %in% epidemic_levels[3:4])

  # Function to create individual subplot
  create_subplot <- function(data, show_x_axis = FALSE) {
    if (half_tiles) {
      # Create polygon data for half tiles
      plot_df <- data |>
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
      bg_df <- data |>
        distinct(epidemic_size, off_R_A, off_R_B) |>
        mutate(x = 0, y = 0)

      # Create subplot with polygons
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
          aes(
            fill = significance,
            group = interaction(off_k_A, off_k_B, method)
          ),
          color = "white",
          linewidth = 0.25
        ) +
        scale_x_continuous(
          breaks = seq_along(levels(data$off_k_A)),
          labels = if (show_x_axis) levels(data$off_k_A) else NULL,
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          breaks = seq_along(levels(data$off_k_B)),
          labels = levels(data$off_k_B),
          expand = c(0, 0)
        ) +
        coord_fixed(
          xlim = c(0.5, length(levels(data$off_k_A)) + 0.5),
          ylim = c(0.5, length(levels(data$off_k_B)) + 0.5)
        )
    } else {
      # Create subplot with full tiles
      plot_df <- data

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
        scale_x_discrete(labels = if (show_x_axis) waiver() else NULL) +
        scale_y_discrete(labels = waiver())
    }

    # Add common elements
    p <- p +
      facet_nested(
        rows = vars(factor(off_R_B, levels = rev(levels(off_R_B)))),
        cols = vars(epidemic_size, off_R_A),
        nest_line = element_line(linetype = 1),
        labeller = grid_labeller
      ) +
      scale_fill_manual(values = significance_pal) +
      geom_segment(
        data = data |> filter(off_k_A == off_k_B, off_R_A == off_R_B),
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
      theme_light(base_size = 11) +
      theme(
        legend.position = if (show_x_axis) "bottom" else "none",
        strip.placement = "outside",
        panel.spacing = unit(1, "mm"),
        strip.switch.pad.grid = unit(0, "pt"),
        strip.switch.pad.wrap = unit(0, "pt"),
        strip.background = element_blank(),
        strip.text = element_markdown(face = "bold"),
        axis.text.x = if (show_x_axis) element_blank() else element_blank(),
        axis.title.x = if (show_x_axis) element_text() else element_blank()
      )

    return(p)
  }

  # Create both subplots
  p1 <- create_subplot(df1, show_x_axis = FALSE)
  p2 <- create_subplot(df2, show_x_axis = TRUE)

  # Combine plots vertically
  library(patchwork)
  p <- p1 /
    p2 +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = NULL,
      subtitle = NULL,
      caption = NULL
    ) &
    labs(y = "k (Forest B)") &
    theme(plot.title = element_blank())

  # Add x-axis label only to bottom plot
  p[[2]] <- p[[2]] + labs(x = "k (Forest A)")
  p[[1]] <- p[[1]] + labs(x = NULL)

  return(p)
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================
