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
