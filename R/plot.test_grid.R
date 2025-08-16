# =================================================
# Test Grid
# =================================================
# A heatmap of the the test results by parameters
plot.test_grid <- function() {
  test_grid <- readRDS("data/test_grid.rds")

  # -------------------------------------------------
  # Colour palette
  # -------------------------------------------------
  # col_fn <- colorRampPalette(c("#e8bee8", "#8c038a")) # pale -> medium-dark
  # epidemic_pal <- setNames(col_fn(length(levels(df$epidemic_size))), levels(df$epidemic_size))

  # library(scales)
  # col_fn <- seq_gradient_pal(low = "#f4aaff", high = "#8c038a", space = "Lab")
  # epidemic_pal <- setNames(col_fn(seq(0, 1, length.out = length(levels(df$epidemic_size)))), levels(df$epidemic_size))

  # library(viridisLite)
  # epidemic_pal <- setNames(cividis(length(levels(df$epidemic_size))), levels(df$epidemic_size))

  epidemic_pal <- c(
    "20"  = "#1a8061", # teal
    "50"  = "#d95f02", # orange
    "100" = "#7570b3", # purple
    "200" = "#e7298a" # pink
  )

  # -------------------------------------------------
  # Dataframe
  # -------------------------------------------------
  df <- test_grid |>
    mutate(
      across(
        c(off_k_A, off_k_B, off_R_A, off_R_B, epidemic_size),
        ~ factor(.x, levels = sort(unique(.x)))
      )
    )

  # -------------------------------------------------
  # Labelling the panels
  # -------------------------------------------------
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
        ifelse(x_char == levels(df$epidemic_size)[1], paste0("Epidemic size = ", x_char), x_char),
        "</span>"
      )
    }
  )

  # -------------------------------------------------
  # ggplot
  # -------------------------------------------------
  df |>
    filter(sample_size == 100) |>
    mutate(significance = if_else(p_value < 0.05, "Significant", "Non-significant")) |>
    ggplot(aes(x = off_k_A, y = off_k_B)) +

    # -------------------------------------------------
    #  Facet
    # -------------------------------------------------
    facet_nested(
      rows = vars(factor(off_R_A, levels = rev(levels(off_R_A)))),
      cols = vars(epidemic_size, off_R_B),
      nest_line = element_line(linetype = 1),
      labeller = grid_labeller
    ) +

    # -------------------------------------------------
    # Colouring the background by epidemic size
    # -------------------------------------------------
    geom_rect(
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = epidemic_size),
      alpha = 0.01, colour = NA
    ) +
    scale_fill_manual(values = epidemic_pal, guide = FALSE) +
    ggnewscale::new_scale_fill() +

    # -------------------------------------------------
    # Tiles
    # -------------------------------------------------
    geom_tile(aes(fill = significance), color = "white", linewidth = 0.25) +
    scale_fill_manual(values = c("Significant" = "#f99820", "Non-significant" = "#03658c")) +

    # -------------------------------------------------
    # Diagonal comparing the same R&k paramd
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
}

# =================================================
# plot.ROC
# =================================================

plot.ROC <- function() {
  # -------------------------------------------------
  # Define multiple alpha thresholds
  # -------------------------------------------------
  alpha_values <- seq(0, 1, 0.001)

  # -------------------------------------------------
  # Compute performance metrics for each alpha
  # -------------------------------------------------
  roc_df <- expand_grid(test_grid, alpha = alpha_values) |>
    mutate(across(
      c(epidemic_size, sample_size),
      ~ factor(.x, levels = sort(unique(.x)))
    )) |>
    group_by(epidemic_size) |>
    mutate(
      true_diff = (off_R_A != off_R_B) | (off_k_A != off_k_B),
      test_reject = p_value < alpha
    ) |>
    group_by(epidemic_size, sample_size, alpha) |>
    summarise(
      TP = sum(true_diff & test_reject),
      FN = sum(true_diff & !test_reject),
      TN = sum(!true_diff & !test_reject),
      FP = sum(!true_diff & test_reject),
      sensitivity = TP / (TP + FN),
      specificity = TN / (TN + FP),
      .groups = "drop"
    )

  alpha_points <- roc_df |>
    mutate(
      alpha_cat = case_when(
        abs(alpha - 0.0015) < 0.0001 ~ "0.0015",
        abs(alpha - 0.01) < 0.001 ~ "0.01",
        abs(alpha - 0.05) < 0.001 ~ "0.05",
        abs(alpha - 0.1) < 0.001 ~ "0.10",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(alpha_cat)) |>
    group_by(epidemic_size, sample_size, alpha_cat) |>
    slice_min(order_by = abs(alpha - as.numeric(alpha_cat)), n = 1) |>
    ungroup() |>
    mutate(alpha_cat = factor(alpha_cat, levels = c("0.0015", "0.01", "0.05", "0.10")))

  # -------------------------------------------------
  # ggplot
  # -------------------------------------------------
  ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, group = interaction(epidemic_size, sample_size))) +
    facet_wrap(vars(epidemic_size),
      labeller = labeller(
        epidemic_size = \(x)  paste0("Epidemic size = ", x)
      )
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
    geom_line(aes(color = factor(sample_size))) +
    scale_color_grey(start = 0.6, end = 0, name = "Forest size:") +


    # Aplha points
    geom_point(
      data = alpha_points,
      aes(
        x = 1 - specificity,
        y = sensitivity,
        shape = alpha_cat
      ),
      stroke = 0.25,
      fill = "black",
      color = "black"
    ) +

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
    coord_cartesian(
      xlim = c(0, 0.25),
      ylim = c(0.75, 1)
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
