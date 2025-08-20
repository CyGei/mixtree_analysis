source("R/packages.R")
source("R/helpers.R")
source("R/plots.R")
set.seed(000)

params_ss <- build_params(
  off_R = 2,
  off_k = 0.3,
  gt_sd = 1.5,
  gt_mu = 3,
  epidemic_size = 100,
  duration = 365
)

params_pois <- build_params(
  off_R = 2,
  off_k = 1e5,
  gt_sd = 1.5,
  gt_mu = 3,
  epidemic_size = 100,
  duration = 365
)

tree_ss <- build_tree(params_ss)
tree_pois <- build_tree(params_pois)

tree_ss |> count(from) |> pull(n) |> summary()
tree_pois |> count(from) |> pull(n) |> summary()

plot_tree(tree_ss)
plot_tree(tree_pois)


forest_ss <- build_forest(tree = tree_ss, params = params_ss)
forest_pois <- build_forest(tree = tree_pois, params = params_pois)
plot_tree(forest_ss[[1]])
plot_tree(forest_pois[[1]])

forest_ss[[1]] |> count(from) |> pull(n) |> summary()
forest_pois[[1]] |> count(from) |> pull(n) |> summary()




sample_size <- 100
sample_ss <- sample(forest_ss, sample_size)
sample_pois <- sample(forest_pois, sample_size)
sample_ss <- lapply(sample_ss, \(x) select(x, from, to) |> as.data.frame())
sample_pois <-lapply(sample_pois, \(x) select(x, from, to) |> as.data.frame())

# A
mixtree::tree_test(
  sample_ss, sample_ss, method = "permanova",
  test_args = list(permutations = 100)
)

mixtree::tree_test(
  sample_ss, sample_ss, method = "chisq",
  test_args = list(simulate.p.value = TRUE, B = 100)
)

#B
mixtree::tree_test(
  sample_ss, sample_pois, method = "permanova",
  test_args = list(permutations = 100)
)

mixtree::tree_test(
  sample_ss, sample_pois, method = "chisq",
  test_args = list(simulate.p.value = TRUE, B = 100)
)

# C
mixtree::tree_test(
  sample_pois, sample_pois, method = "permanova",
  test_args = list(permutations = 100)
)
mixtree::tree_test(
  sample_pois, sample_pois, method = "chisq",
  test_args = list(simulate.p.value = TRUE, B = 100)
)



# Create binary truth: 0 = same distribution, 1 = different distribution
roc_df <- test_grid |>
  filter(sample_size == 100) |>
  mutate(
    truth = as.integer(id_A != id_B),
    across(
      .cols = c(method, off_k_A, off_k_B, off_R_A, off_R_B, epidemic_size),
      .fns = as.factor
    )
  ) |>
  group_by(method, off_k_A, off_k_B, off_R_A, off_R_B, epidemic_size) |>
  group_modify(
    ~ {
      roc_obj <- pROC::roc(
        response = .x$truth,
        predictor = .x$p_value,
        direction = ">", # smaller p => more likely positive
        quiet = TRUE
      )

      roc_coords <- pROC::coords(
        roc_obj,
        "all",
        ret = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      )

      tibble(
        threshold = roc_coords$threshold,
        sensitivity = roc_coords$sensitivity,
        specificity = roc_coords$specificity,
        AUC = as.numeric(pROC::auc(roc_obj))
      )
    }
  ) |>
  ungroup()

head(roc_df)

ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = method)) +
  geom_path() +
  facet_grid(cols = vars(epidemic_size), rows = vars(sample_size))

# DON'T filter out null cases - they're essential!
roc_data <- test_grid %>%
  filter(sample_size == 100) %>%
  mutate(
    truth = as.integer(id_A != id_B),
    type_div = case_when(
      off_R_A != off_R_B & off_k_A == off_k_B ~ "R_only",
      off_R_A == off_R_B & off_k_A != off_k_B ~ "k_only",
      off_R_A != off_R_B & off_k_A != off_k_B ~ "R_and_k",
      TRUE ~ "null" # Keep these - they're your true negatives!
    )
  ) %>%
  group_by(type_div, method) %>%
  do({
    # Create ROC curve using varying p-value thresholds
    thresholds <- c(0, sort(unique(.$p_value)), 1)

    map_dfr(thresholds, function(thresh) {
      predictions <- as.integer(.$p_value <= thresh) # Reject H0 if p <= threshold

      tibble(
        threshold = thresh,
        TPR = sum(predictions * .$truth) / sum(.$truth), # Sensitivity
        FPR = sum(predictions * (1 - .$truth)) / sum(1 - .$truth), # 1 - Specificity
        n_tp = sum(predictions * .$truth),
        n_fp = sum(predictions * (1 - .$truth)),
        n_tn = sum((1 - predictions) * (1 - .$truth)),
        n_fn = sum((1 - predictions) * .$truth)
      )
    })
  }) %>%
  ungroup()
