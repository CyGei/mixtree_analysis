# Main figures.
source("R/packages.R")
source("R/helpers.R")

# ------------------------------------
#           ROC curves
# ------------------------------------
#' @title roc_df
#' @description
#' A tibble containing data for plotting ROC curves from the simulation study `results_grid.rds`.
#' It computes true positive rates (TPR) and false positive rates (FPR) across varying p-value thresholds.
#' A true positive is defined as correctly rejecting the null hypothesis when forests differ in parameters.
#' A false positive is incorrectly rejecting the null hypothesis when forests are generated under the same parameters.
#' Columns:
#' - `method`: statistical test used ("chisq" or "permanova")
#' - `forest_size`: number of transmission trees per forest
#' - `epidemic_size`: number of vertices per transmission tree
#' - `threshold`: significance level (alpha)
#' - `TPR`: true positive rate (sensitivity)
#' - `FPR`: false positive rate (1 - specificity)
roc_df <- readRDS("data/results_grid.rds") |>
  mutate(true_different = param_id_A != param_id_B) |>
  group_by(method, forest_size, epidemic_size) |>
  group_modify(~ compute_roc(.x)) |>
  ungroup() |>
  mutate(
    label_epidemic_size = "Epidemic size",
    label_forest_size = "Forest size"
  )


ggplot(roc_df, aes(x = FPR, y = TPR, color = method)) +
  facet_nested(
    cols = vars(label_epidemic_size, epidemic_size),
    rows = vars(label_forest_size, forest_size),
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
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  geom_line() +
  scale_color_manual(
    values = c("permanova" = "#e68613", "chisq" = "#1f65cc"),
    name = "Method:"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = TRUE,
    clip = "off"
  ) +
  labs(
    x = "1 - specificity",
    y = "Sensitivity",
    color = "Forest size:"
  ) +
  theme_mixtree() +
  theme(axis.text.x = element_blank())

ggsave(
  filename = "figures/roc.svg",
  width = 8,
  height = 8
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

delta_df <- readRDS("data/results_grid.rds") |>
  mutate(
    H0 = param_id_A == param_id_B,
    reject_H0 = p_value < 0.05,
    comparison_type = case_when(
      H0 ~ "H₀",
      off_R_A != off_R_B & off_k_A == off_k_B ~ "\U0394 R₀",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "\U0394 \U1D458",
      TRUE ~ "\U0394 R₀ \U0026 \U0394 \U1D458"
    ),
    comparison_type = factor(
      comparison_type,
      levels = c(
        "\U0394 R₀ \U0026 \U0394 \U1D458",
        "\U0394 R₀",
        "\U0394 \U1D458",
        "H₀"
      )
    )
  ) |>
  summarise(
    rejection_rate = mean(reject_H0),
    .by = c(method, forest_size, epidemic_size, comparison_type)
  )
glimpse(delta_df)

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
    group = method,
    color = method,
    fill = method
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
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.35
  ) +
  geom_point(
    aes(shape = method),
    color = "black",
    size = 3,
    stroke = 0.25
  ) +
  scale_color_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    guide = "none"
  ) +
  scale_fill_manual(
    values = mixtree_pal,
    labels = mixtree_lab,
    name = "Method:"
  ) +
  scale_shape_manual(
    values = c("permanova" = 21, "chisq" = 24),
    labels = mixtree_lab,
    name = "Method:"
  ) +
  scale_linetype_manual(
    values = c("Target" = "solid"),
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
  filename = "figures/power_curves.svg",
  width = 10,
  height = 10
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

model_df <- readRDS("data/results_grid.rds") |>
  mutate(
    H0 = param_id_A == param_id_B,
    reject_H0 = p_value < 0.05,
    accept_H0 = !reject_H0,
    delta_R0 = abs(off_R_A - off_R_B),
    delta_k = abs(off_k_A - off_k_B),
    log_delta_k = log1p(delta_k),
    delta_k_cat = case_when(
      delta_k == 0 ~ "none",
      delta_k < 1000 ~ "finite",
      TRUE ~ "infinite"
    ) |>
      factor(levels = c("none", "finite", "infinite")),
    forest_size_cat = factor(forest_size, levels = sort(unique(forest_size))),
    epidemic_size_cat = factor(
      epidemic_size,
      levels = sort(unique(epidemic_size))
    ),
    method = factor(method, levels = c("chisq", "permanova")),
  )
glimpse(model_df)

# --------- Sensitivity Model ---------
m1 <- glm(
  reject_H0 ~ method + forest_size + epidemic_size + delta_R0 * log_delta_k,
  data = model_df |> filter(H0 == FALSE),
  family = binomial(link = "logit")
)
m2 <- glm(
  reject_H0 ~ method + forest_size + epidemic_size + delta_R0 * delta_k_cat,
  data = model_df |> filter(H0 == FALSE),
  family = binomial(link = "logit")
)
m3 <- glm(
  reject_H0 ~ method +
    forest_size_cat +
    epidemic_size_cat +
    delta_R0 * delta_k_cat,
  data = model_df |> filter(H0 == FALSE),
  family = binomial(link = "logit")
)
m4 <- glm(
  reject_H0 ~ method +
    forest_size_cat +
    epidemic_size_cat +
    method:epidemic_size_cat +
    delta_R0 * delta_k_cat,
  data = model_df |> filter(H0 == FALSE),
  family = binomial(link = "logit")
)
m5 <- glm(
  reject_H0 ~
    forest_size_cat +
    epidemic_size_cat +
    method:epidemic_size_cat +
    method * delta_R0 * delta_k_cat,
  data = model_df |> filter(H0 == FALSE),
  family = binomial(link = "logit")
)

map_dfr(
  list(m1, m2, m3, m4, m5),
  ~ tibble(
    AIC = AIC(.x),
    BIC = BIC(.x)
  ),
  .id = "model"
) |>
  arrange(AIC)

map_dfr(
  list(m1, m2, m3, m4, m5),
  ~ performance::r2(.x) |>
    as_tibble(),
  .id = "model"
) |>
  arrange(desc(model))

broom::tidy(m5, conf.int = FALSE) |>
  mutate(
    odds_ratio = exp(estimate),
    across(where(is.numeric), ~ round(.x, 2))
  ) |>
  select(term, odds_ratio, p.value) |>
  View()
# kableExtra::kable(
#   format = "latex",
#   booktabs = TRUE,
#   col.names = c("Predictor", "Odds Ratio", "p-value"),
#   caption = "Model Results"
# ) |>
#   kableExtra::kable_styling(
#     latex_options = c("hold_position", "scale_down")
#   )

# --------- Specificity Model ---------
m1 <- glm(
  accept_H0 ~ method +
    forest_size_cat +
    epidemic_size_cat,
  data = model_df |> filter(H0 == TRUE),
  family = binomial(link = "logit")
)

m2 <- glm(
  accept_H0 ~ method +
    forest_size_cat +
    epidemic_size_cat +
    off_R_A +
    off_k_A,
  data = model_df |> filter(H0 == TRUE),
  family = binomial(link = "logit")
)
m3 <- glm(
  accept_H0 ~ method +
    forest_size_cat +
    epidemic_size_cat +
    method:epidemic_size_cat +
    off_R_A +
    off_k_A,
  data = model_df |> filter(H0 == TRUE),
  family = binomial(link = "logit")
)
m4 <- glm(
  accept_H0 ~ method * forest_size_cat + epidemic_size_cat,
  data = model_df |> filter(H0 == TRUE),
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
  ~ performance::r2(.x),
  .id = "model"
) |>
  arrange(desc(model))

broom::tidy(m1, conf.int = FALSE) |>
  mutate(
    odds_ratio = exp(estimate),
    across(where(is.numeric), ~ round(.x, 2))
  ) |>
  select(term, odds_ratio, p.value) |>
  print(n = Inf)

# average diff in rejection rate between methods for
# forest_size = 20 across all epidemic sizes under comparison_type = H0
delta_df |>
  filter(
    forest_size >= 100,
    comparison_type %in% c("\U0394 R₀ \U0026 \U0394 \U1D458", "\U0394 \U1D458")
  ) |>
  pivot_wider(
    names_from = method,
    values_from = rejection_rate
  ) |>
  select(-chisq) |>
  summarise(
    avg_rejection_rate = mean(permanova) * 100
  )
mutate(
  diff = chisq - permanova
) |>
  summarise(
    avg_diff = mean(diff) * 100
  )
