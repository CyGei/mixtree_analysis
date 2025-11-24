# Main figures.
source("R/packages.R")
source("R/helpers.R")

# ------------------------------------
#          Data
# ------------------------------------
results_grid <- readRDS("data/results_grid.rds") |>
  mutate(
    H1 = param_id_A != param_id_B,
    reject_H0 = p_value < 0.05,

    delta_R0 = abs(off_R_A - off_R_B),
    delta_k = abs(off_k_A - off_k_B),

    # Annotations:
    comparison_type = case_when(
      !H1 ~ "H₀",
      off_R_A != off_R_B & off_k_A == off_k_B ~ "\U0394 R₀",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "\U0394 \U1D458",
      .default = "\U0394 R₀ \U0026 \U0394 \U1D458"
    ),
    comparison_type = factor(
      comparison_type,
      levels = c(
        "\U0394 R₀ \U0026 \U0394 \U1D458",
        "\U0394 R₀",
        "\U0394 \U1D458",
        "H₀"
      )
    ),

    R0_labels = sprintf(
      "%.1f v %.1f",
      pmin(off_R_A, off_R_B),
      pmax(off_R_A, off_R_B)
    ) |>
      forcats::fct_reorder(delta_R0),

    k_labels = sprintf(
      "%s v %s",
      format_k(pmin(off_k_A, off_k_B)),
      format_k(pmax(off_k_A, off_k_B))
    ) |>
      forcats::fct_reorder(delta_k),

    delta_R0_labels = case_when(
      delta_R0 == 0 ~ "ΔR₀ = 0",
      TRUE ~ paste0("ΔR₀ = ", delta_R0)
    ) |>
      factor(levels = c("ΔR₀ = 0", "ΔR₀ = 0.5", "ΔR₀ = 1", "ΔR₀ = 1.5")),

    delta_k_labels = case_when(
      delta_k == 0 ~ "Δ\U1D458 = 0",
      delta_k >= 99999 ~ "Δ\U1D458 = 10⁵",
      TRUE ~ paste0("Δ\U1D458 = ", delta_k)
    ) |>
      as.factor()
  )

# ------------------------------------
#           Power curves
# ------------------------------------
#' @title delta_df
#' @description
#' A tibble summarising results from the simulation study `results_grid.rds`.
#' It computes rejection/acceptances rates of the null hypothesis (H₀) when:
#' - `H₀` is true (forests generated under same parameters)
#' - `ΔR₀` forests differ in reproduction number only
#' - `Δk` forests differ in dispersion parameter only
#' - `ΔR₀ & Δk` forests differ in both parameters.
#' Rates are calculted for each `method`, `forest_size`, and `epidemic_size`.
#' Columns:
#' - `method`: statistical test used ("chisq" or "permanova")
#' - `forest_size`: number of transmission trees per forest
#' - `epidemic_size`: number of vertices per transmission tree
#' - `comparison_type`: type of parameter difference being tested (ΔR₀, Δk, both, or H₀)
#' - `rejection_rate`: proportion of tests rejecting H₀ at significance level 0.05.

delta_df <- results_grid |>
  summarise(
    rejection_rate = mean(reject_H0),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )

delta_df |>
  mutate(
    epidemic_size = factor(epidemic_size),
    label_forest_size = "Forest size",
    label_comparison_type = "Comparison type",
    hline = case_when(
      comparison_type == "H₀" ~ 0, # Should be 0.05
      TRUE ~ 1
    )
  ) |>
  ggplot(aes(
    x = epidemic_size,
    y = rejection_rate,
    group = method
  )) +
  facet_nested(
    cols = vars(label_forest_size, forest_size),
    rows = vars(label_comparison_type, comparison_type),
    scales = "free_y",
    strip = strip_nested(
      text_x = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      text_y = list(
        element_text(size = 17),
        element_text(size = 16)
      ),
      by_layer_x = TRUE,
      by_layer_y = TRUE
    )
  ) +
  geom_hline(
    aes(yintercept = hline, linetype = "Target"),
    color = "black",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.35
  ) +
  geom_point(
    aes(shape = method),
    size = 2.25
  ) +
  scale_shape_manual(
    values = c("permanova" = 19, "chisq" = 15),
    labels = mixtree_lab,
    name = "Method:"
  ) +
  scale_linetype_manual(
    values = c("Target" = "dashed"),
    labels = c("Target" = "Target"),
    name = ""
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      comparison_type == "H₀" ~ scale_y_continuous(
        limits = c(0, 0.03),
        breaks = seq(0, 0.03, by = 0.01),
        labels = c("0", "1", "2", "3")
      ),
      comparison_type == "\U0394 \U1D458" ~ scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      ),
      comparison_type == "\U0394 R₀" ~ scale_y_continuous(
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      ),
      comparison_type == "\U0394 R₀ \U0026 \U0394 \U1D458" ~ scale_y_continuous(
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "25", "50", "75", "100")
      )
    )
  ) +
  labs(
    x = "Epidemic size",
    y = "\UFE6A rejecting H₀"
  ) +
  theme_mixtree() +
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  )

ggsave(
  filename = "figures/power_curves.png",
  width = 7.5,
  height = 7.5,
  units = "in",
  dpi = 300
)
ggsave(
  filename = "figures/power_curves.eps",
  width = 7.5,
  height = 7.5,
  units = "in",
  device = cairo_ps,
  dpi = 300
)


# ------------------------------------
#           Regression model
# ------------------------------------
#' @title model_df
#' @description
#' A tibble prepared for regression analysis from the simulation study `results_grid.rds`.
#' It includes indicators for whether the null hypothesis (H₀) was rejected,
#' absolute differences in reproduction number (`delta_R0`) and dispersion parameter (`delta_k`)
#' between forests, and categorical variables for forest size, epidemic size, and method used.

model_df <- results_grid |>
  mutate(
    #factoring
    delta_k = case_when(
      delta_k == 0 ~ "none",
      delta_k < 1000 ~ "finite",
      TRUE ~ "infinite"
    ) |>
      factor(levels = c("none", "finite", "infinite")),
    forest_size = factor(forest_size, levels = sort(unique(forest_size))),
    epidemic_size = factor(
      epidemic_size,
      levels = sort(unique(epidemic_size))
    ),
    method = factor(method, levels = c("chisq", "permanova")),
  ) |>
  select(
    H1,
    reject_H0,
    method,
    forest_size,
    epidemic_size,
    delta_R0,
    delta_k
  )

# --------- Sensitivity Model ---------
m1 <- glm(
  reject_H0 ~ method + forest_size + epidemic_size + delta_R0 + delta_k,
  data = model_df |> filter(H1),
  family = binomial(link = "logit")
)

m2 <- glm(
  reject_H0 ~ method + forest_size + epidemic_size + delta_R0 * delta_k,
  data = model_df |> filter(H1),
  family = binomial(link = "logit")
)

m3 <- glm(
  reject_H0 ~ method +
    forest_size +
    epidemic_size +
    method:epidemic_size +
    delta_R0 * delta_k,
  data = model_df |> filter(H1),
  family = binomial(link = "logit")
)

m4 <- glm(
  reject_H0 ~
    method +
    forest_size +
    epidemic_size +
    method:epidemic_size +
    method * delta_R0 * delta_k,
  data = model_df |> filter(H1),
  family = binomial(link = "logit")
)

map_dfr(
  list(m1, m2, m3, m4),
  ~ tibble(
    AIC = AIC(.x),
    BIC = BIC(.x)
  ),
  .id = "model"
) |>
  arrange(AIC)

map_dfr(
  list(m1, m2, m3, m4),
  ~ performance::r2(.x) |>
    as_tibble(),
  .id = "model"
) |>
  arrange(desc(model))

broom::tidy(m4, conf.int = FALSE) |>
  mutate(
    odds_ratio = exp(estimate),
    across(where(is.numeric), ~ round(.x, 2))
  ) |>
  select(term, odds_ratio, p.value) |>
  print(n = Inf)
