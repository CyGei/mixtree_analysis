# =================================================
# plot_helpers() returns all internal helper functions
# =================================================
plot_helpers <- function() {
  palettes <- function(name = c("epidemic_size", "significance")) {
    name <- match.arg(name)
    pals <- list(
      epidemic_size = c(
        "20" = "#1a8061",
        "50" = "#d95f02",
        "100" = "#7570b3",
        "200" = "#e7298a"
      ),
      significance = c(
        "Significant" = "#f99820",
        "Non-significant" = "#03658c"
      )
    )
    pals[[name]]
  }

  test_grid_labeller <- function(df) {
    labeller(
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
          ifelse(
            x_char %in% names(palettes("epidemic_size")),
            palettes("epidemic_size")[x_char],
            "black"
          ),
          "; font-size:14pt'>",
          ifelse(
            x_char == levels(df$epidemic_size)[1],
            paste0("Epidemic size = ", x_char),
            x_char
          ),
          "</span>"
        )
      }
    )
  }

  list(
    palettes = palettes,
    test_grid_labeller = test_grid_labeller
    # add more functions here
  )
}


# =================================================
# plot_tree
# =================================================

plot_tree <- function(tree, layout = c("cactustree", "treemap", "stress")) {
  g <- tidygraph::as_tbl_graph(tree) |>
    activate(nodes) |>
    mutate(
      date = case_when(
        name == "1" ~ 0,
        TRUE ~ tree$date[match(name, tree$to)]
      )
    )

  # Create tree layout and replace x with date
  layout_df <- create_layout(g, layout = match.arg(layout)) |>
    mutate(x = date)

  max_date <- max(tree$date, na.rm = TRUE)

  ggraph::ggraph(layout_df) +
    geom_node_label(
      aes(label = name),
      colour = "black",
      size = 5,
    ) +
    geom_edge_link(
      arrow = arrow(length = unit(1.5, 'mm')),
      start_cap = circle(3, 'mm'),
      end_cap = circle(3, 'mm')
    ) +
    scale_x_continuous(breaks = seq(0, max_date, 1)) +
    labs(x = "Date", y = "") +
    theme_classic() +
    # add vertical grid lines
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", size = 0.5)
    )
}

# =================================================
# test_grid
# =================================================

plot_test_grid <- function(sample_size, method = NULL, half_tiles = FALSE) {
  helpers <- plot_helpers()
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

  # Split data into two groups for vertical stacking
  epidemic_levels <- levels(df$epidemic_size)
  n_levels <- length(epidemic_levels)
  split_point <- ceiling(n_levels / 2)

  # First subplot: first half of levels
  df1 <- df |> filter(epidemic_size %in% epidemic_levels[1:split_point])

  # Second subplot: second half of levels
  df2 <- df |>
    filter(epidemic_size %in% epidemic_levels[(split_point + 1):n_levels])

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
        scale_fill_manual(
          values = helpers$palettes("epidemic_size"),
          guide = FALSE
        ) +
        ggnewscale::new_scale_fill() +
        geom_polygon(
          aes(
            fill = significance,
            group = interaction(off_k_A, off_k_B, method)
          ),
          color = "white",
          linewidth = 0.25
        ) +
        # Empty geom for test method
        geom_point(aes(x = 1, y = 1, shape = method), size = 0, alpha = 0) +
        scale_shape_manual(
          "",
          values = c("chisq" = 16, "permanova" = 17),
          labels = c(
            "chisq" = "\u25E3 χ² test",
            "permanova" = "\u25E5 PERMANOVA"
          ),
          guide = guide_legend(
            override.aes = list(size = 0, alpha = 0),
            keywidth = unit(0, "pt"), # remove key
            keyheight = unit(0, "pt"),
            default.unit = "pt"
          )
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
        scale_fill_manual(
          values = helpers$palettes("epidemic_size"),
          guide = FALSE
        ) +
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
        rows = vars(factor(off_R_A, levels = rev(levels(off_R_A)))),
        cols = vars(epidemic_size, off_R_B),
        nest_line = element_line(linetype = 1),
        labeller = helpers$test_grid_labeller(plot_df)
      ) +
      scale_fill_manual(values = helpers$palettes("significance")) +
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
      theme_classic(base_size = 15) +
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

  # Combine plots vertically using wrap_plots
  p1 <- p1 + labs(x = NULL, y = "k (Forest B)")
  p2 <- p2 + labs(x = "k (Forest A)", y = "k (Forest B)")

  p <- patchwork::wrap_plots(
    p1,
    p2,
    ncol = 1,
    heights = c(1, 1)
  )

  return(p)
}

# =================================================
# plot.ROC
# =================================================

plot.ROC <- function() {
  # -------------------------------------------------
  # Define multiple alpha thresholds
  # -------------------------------------------------
  alpha_values <- seq(0, 1, 0.01)

  # -------------------------------------------------
  # Compute performance metrics for each alpha
  # -------------------------------------------------
  roc_df <- expand_grid(test_grid, alpha = alpha_values) |>
    mutate(across(
      c(epidemic_size, sample_size, method),
      ~ factor(.x, levels = sort(unique(.x)))
    )) |>
    group_by(epidemic_size, sample_size, method, alpha) |>
    mutate(
      true_diff = (off_R_A != off_R_B) | (off_k_A != off_k_B),
      test_reject = p_value < alpha
    ) |>
    summarise(
      TP = sum(true_diff & test_reject),
      FN = sum(true_diff & !test_reject),
      TN = sum(!true_diff & !test_reject),
      FP = sum(!true_diff & test_reject),
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      .groups = "drop"
    )

  # alpha_points <- roc_df |>
  #   mutate(
  #     alpha_cat = case_when(
  #       abs(alpha - 0.0015) < 0.0001 ~ "0.0015",
  #       abs(alpha - 0.01) < 0.001 ~ "0.01",
  #       abs(alpha - 0.05) < 0.001 ~ "0.05",
  #       abs(alpha - 0.1) < 0.001 ~ "0.10",
  #       TRUE ~ NA_character_
  #     )
  #   ) |>
  #   filter(!is.na(alpha_cat)) |>
  #   group_by(epidemic_size, sample_size, alpha_cat) |>
  #   slice_min(order_by = abs(alpha - as.numeric(alpha_cat)), n = 1) |>
  #   ungroup() |>
  #   mutate(
  #     alpha_cat = factor(
  #       alpha_cat,
  #       levels = c("0.0015", "0.01", "0.05", "0.10")
  #     )
  #   )

  # -------------------------------------------------
  # ggplot
  # -------------------------------------------------
  ggplot(
    roc_df,
    aes(
      x = 1 - specificity,
      y = sensitivity,
      group = interaction(epidemic_size, sample_size, method)
    )
  ) +
    facet_grid(
      cols = vars(epidemic_size),
      rows = vars(sample_size)
    ) +
    # Reference line
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "gray50",
      alpha = 0.5,
      linewidth = 0.25
    ) +

    # ROC line
    geom_line(aes(color = method)) +
    #scale_color_grey(start = 0.6, end = 0, name = "Forest size:") +
    scale_color_manual(
      values = c("permanova" = "#af245a", "chisq" = "#1f9969"),
      name = "Method:"
    ) +

    # # Aplha points
    # geom_point(
    #   data = alpha_points,
    #   aes(
    #     x = 1 - specificity,
    #     y = sensitivity,
    #     shape = alpha_cat
    #   ),
    #   stroke = 0.25,
    #   fill = "black",
    #   color = "black"
    # ) +

    # Theme
    scale_shape_manual(
      values = c(
        "0.0015" = 21,
        "0.01" = 22,
        "0.05" = 23,
        "0.10" = 24
      ),
      name = "α:",
      na.translate = FALSE
    ) +
    labs(
      x = "1 - specificity",
      y = "Sensitivity",
      color = "Forest size:"
    ) +
    theme_light(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = -5),
      legend.box.just = "left",
      legend.spacing = unit(0.1, "cm"),
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = element_text(size = 16, face = "bold", color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.3, "lines")
    )
}

# AUC ---------------------------------------------------------------------
# roc_df %>%
#   group_by(method, sample_size, epidemic_size) %>%
#   summarise(AUC = mean(AUC)) %>%
#   mutate(method = factor(method, levels = c("χ² test", "PERMANOVA"))) %>%
#   ggplot() +
#   aes(x = sample_size, y = AUC) +
#   facet_grid(~method) +
#   geom_line(aes(group = epidemic_size, color = epidemic_size)) +
#   geom_point(aes(color = epidemic_size), size = 2.5) +
#   theme_bw() +
#   scale_y_continuous(limits = c(0.92, 1)) +
#   scale_color_manual(values = colorRampPalette(c("grey", "black"))(5)) +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.margin = margin(t = -5),
#     legend.box.just = "left",
#     legend.spacing = unit(0.1, "cm"),
#     strip.background = element_rect(
#       fill = "#f8f4f2",
#       colour = "black",
#       linewidth = 0.5
#     ),
#     strip.text = element_text(family = "Fira Code", size = 10),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(0.25, "lines")
#   ) +
#   labs(
#     x = "Forest Size",
#     y = "Area Under Curve (AUC)",
#     color = "Epidemic Size:"
#   )

# ggplot2::ggsave(
#   "auc.png",
#   width = 120,
#   height = 100,
#   units = "mm",
#   dpi = 300
# )

# # precision - recall curve

# results %>%
#   mutate(
#     actual_positive = overlap_freq != "1",
#     predicted_positive = p_value < alpha
#   ) %>%
#   group_by(method, epidemic_size, sample_size) %>%
#   summarise(
#     precision = sum(predicted_positive & actual_positive) /
#       sum(predicted_positive),
#     recall = sum(predicted_positive & actual_positive) / sum(actual_positive),
#     .groups = "drop"
#   ) %>%
#   ggplot(aes(
#     x = recall,
#     y = precision,
#     group = interaction(epidemic_size, method)
#   )) +
#   facet_grid(~method) +
#   geom_line(aes(colour = epidemic_size)) +
#   geom_point(aes(fill = sample_size), size = 2, shape = 21, color = "black") +
#   scale_color_manual(values = colorRampPalette(c("grey", "black"))(5)) +
#   scale_fill_viridis_d(option = "D") +
#   # scale_x_continuous(
#   #   labels = c("0.80", "0.85", "0.90", "0.95", "1")
#   # )+
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.margin = margin(t = -5),
#     legend.box.just = "left",
#     legend.spacing = unit(0.1, "cm"),
#     strip.background = element_rect(
#       fill = "#f8f4f2",
#       colour = "black",
#       linewidth = 0.5
#     ),
#     strip.text = element_text(family = "Fira Code", size = 10),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(0.9, "lines")
#   ) +
#   labs(
#     x = "Sensitivity (Recall)",
#     y = "Positive Predictive Value (Precision)",
#     color = "Epidemic Size:",
#     fill = "Forest Size:"
#   )

# ggplot2::ggsave(
#   "pr_curve.png",
#   width = 120,
#   height = 100,
#   units = "mm",
#   dpi = 300
# )
