# 04_models.R
# Logistic regression vs random forest for KSI, overall and per road-user group.
# Input:  output/model_data.rds (from 02_clean_features.R)
# Output: output/models/metrics.csv        test-set AUC, accuracy, recall, precision, F1
#         output/models/auc_tests.csv      DeLong test: two-sided, do RF and logistic AUCs differ?
#         output/models/odds_ratios.csv    odds ratios with 95% CI (full training data)
#         output/models/rf_importance.csv  permutation importance
#         output/models/predictions.rds    test-set probabilities from both models
#         output/models/*.png              ROC curves, importance, odds-ratio charts
#
# Design:
# - Train 2012-2021, test 2022-2024 (2024 ends 30 Sep).
# - Class imbalance: models used for prediction are fitted on a downsampled,
#   balanced training set (all KSI cases + an equal random sample of non-KSI),
#   so 0.5 is a fair threshold for both. Test data is never resampled.
# - Odds ratios come from a logistic regression on the full training set.
#   Sampling on the outcome only shifts the intercept, so the odds ratios
#   estimate the same thing, with tighter intervals.

library(data.table)
library(ggplot2)
library(ranger)
library(pROC)
library(caret)

set.seed(42)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)
d <- readRDS("output/model_data.rds")

features <- c("sex", "age_group", "road_user", "protection", "vehicle_type", "vehicle_age",
              "speed_zone", "light", "weather", "surface", "road_geometry", "crash_type",
              "location", "hour_band", "weekend")
groups <- list(All = d, Occupant = d[user_group == "Occupant"], Vulnerable = d[user_group == "Vulnerable"])
col <- c(Logistic = "#2a78d6", `Random forest` = "#eb6834")   # categorical slots 1 and 2

metrics <- list(); auc_tests <- list(); ors <- list(); imps <- list(); preds <- list(); rocs <- list()

for (g in names(groups)) {
  cat("\n====", g, "====\n")
  x <- copy(groups[[g]])
  x[, (features) := lapply(.SD, droplevels), .SDcols = features]
  fx <- features[sapply(features, function(f) nlevels(x[[f]]) > 1)]   # skip factors with one level in this group
  form <- reformulate(fx, "KSI")
  train <- x[split == "train"]; test <- x[split == "test"]

  # Balanced training sample
  pos <- train[KSI == 1]
  bal <- rbind(pos, train[KSI == 0][sample(.N, nrow(pos))])
  cat(sprintf("train %s (KSI %s) -> balanced %s; test %s (KSI %.1f%%)\n",
              format(nrow(train), big.mark = ","), format(nrow(pos), big.mark = ","),
              format(nrow(bal), big.mark = ","), format(nrow(test), big.mark = ","), 100 * mean(test$KSI)))

  # Odds ratios (full training data, Wald 95% CI)
  full <- glm(form, data = train, family = binomial)
  s <- summary(full)$coefficients
  ors[[g]] <- data.table(group = g, term = rownames(s), estimate = s[, 1], se = s[, 2], p_value = s[, 4])[
    , `:=`(odds_ratio = exp(estimate), ci_low = exp(estimate - 1.96 * se), ci_high = exp(estimate + 1.96 * se))][
    term != "(Intercept)"]

  # Prediction models on the balanced sample
  lr <- glm(form, data = bal, family = binomial)
  rf <- ranger(form, data = bal[, c("KSI", fx), with = FALSE][, KSI := factor(KSI)],
               num.trees = 300, probability = TRUE, importance = "permutation", seed = 42)
  p <- data.table(group = g, ACCIDENT_NO = test$ACCIDENT_NO, PERSON_ID = test$PERSON_ID, KSI = test$KSI,
                  Logistic = predict(lr, test, type = "response"),
                  `Random forest` = predict(rf, test)$predictions[, "1"])
  preds[[g]] <- p
  imps[[g]] <- data.table(group = g, feature = names(rf$variable.importance), importance = rf$variable.importance)

  # Evaluation on the untouched test set
  truth <- factor(p$KSI, levels = c(0, 1))
  for (m in names(col)) {
    r <- roc(p$KSI, p[[m]], levels = c(0, 1), direction = "<", quiet = TRUE)
    rocs[[paste(g, m)]] <- data.table(group = g, model = m, fpr = 1 - r$specificities, tpr = r$sensitivities)
    cm <- confusionMatrix(factor(as.integer(p[[m]] >= 0.5), levels = c(0, 1)), truth, positive = "1", mode = "everything")
    ci <- ci.auc(r)
    metrics[[paste(g, m)]] <- data.table(group = g, model = m, auc = as.numeric(auc(r)), auc_low = ci[1], auc_high = ci[3],
      accuracy = cm$overall[["Accuracy"]], recall = cm$byClass[["Recall"]], precision = cm$byClass[["Precision"]],
      f1 = cm$byClass[["F1"]], specificity = cm$byClass[["Specificity"]],
      tp = cm$table[2, 2], fp = cm$table[2, 1], fn = cm$table[1, 2], tn = cm$table[1, 1])
    if (m == "Logistic") r_lr <- r else r_rf <- r
  }
  t <- roc.test(r_rf, r_lr, method = "delong", alternative = "two.sided")
  auc_tests[[g]] <- data.table(group = g, auc_rf = as.numeric(auc(r_rf)), auc_lr = as.numeric(auc(r_lr)),
                               z = unname(t$statistic), p_value = t$p.value)
}

metrics <- rbindlist(metrics); auc_tests <- rbindlist(auc_tests); ors <- rbindlist(ors); imps <- rbindlist(imps)
fwrite(metrics, "output/models/metrics.csv")
fwrite(auc_tests, "output/models/auc_tests.csv")
fwrite(ors, "output/models/odds_ratios.csv")
fwrite(imps, "output/models/rf_importance.csv")
saveRDS(rbindlist(preds), "output/models/predictions.rds")

cat("\nTest-set metrics (threshold 0.5):\n")
print(metrics[, .(group, model, auc = round(auc, 3), accuracy = round(accuracy, 3), recall = round(recall, 3),
                  precision = round(precision, 3), f1 = round(f1, 3))])
cat("\nDeLong test, two-sided, H1: AUCs differ:\n")
print(auc_tests[, .(group, auc_rf = round(auc_rf, 4), auc_lr = round(auc_lr, 4), z = round(z, 2), p_value = signif(p_value, 3))])

# ---- Charts -----------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
                  plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA)))

# ROC curves, one panel per group, AUC in the legend labels
roc_dt <- rbindlist(rocs)
roc_dt <- roc_dt[, .SD[unique(c(seq(1, .N, length.out = 500), .N))], by = .(group, model)]   # thin for plotting
lab <- metrics[, .(group, model, label = sprintf("%s AUC %.3f", model, auc))]
g <- ggplot(roc_dt, aes(fpr, tpr, colour = model)) +
  geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_text(data = lab, aes(x = 1, y = ifelse(model == "Logistic", 0.12, 0.04), label = label),
            hjust = 1, size = 3, show.legend = FALSE) +
  facet_wrap(~ group) + coord_equal() +
  scale_colour_manual(values = col) +
  labs(title = "ROC curves on the 2022-2024 test set", x = "False positive rate", y = "True positive rate")
ggsave("output/models/roc.png", g, width = 9, height = 4, dpi = 150)

# Random forest permutation importance per group
imps[, feature := reorder(feature, importance)]
g <- ggplot(imps, aes(importance, feature)) +
  geom_col(fill = col[["Random forest"]], width = 0.6) +
  facet_wrap(~ group, scales = "free_x") +
  labs(title = "Random forest permutation importance", x = "Decrease in accuracy when permuted", y = NULL) +
  theme(panel.grid.major.y = element_blank())
ggsave("output/models/rf_importance.png", g, width = 9, height = 5, dpi = 150)

# Odds ratios, one chart per group (log scale, reference categories omitted)
for (gname in names(groups)) {
  o <- ors[group == gname]
  o[, feature := sapply(term, function(t) features[startsWith(t, features)][1])]
  o[, level := substring(term, nchar(feature) + 1)]
  ref <- sapply(features, function(f) levels(droplevels(groups[[gname]][[f]]))[1])
  o[, feature := factor(sprintf("%s\n(ref: %s)", feature, ref[feature]),
                        levels = sprintf("%s\n(ref: %s)", features, ref))]
  o[, term := factor(term, levels = rev(term))]          # keep factor-level order top to bottom
  g <- ggplot(o, aes(odds_ratio, term)) +
    scale_y_discrete(labels = setNames(o$level, o$term)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", colour = col[["Logistic"]], linewidth = 0.6) +
    geom_point(colour = col[["Logistic"]], size = 2) +
    scale_x_log10() +
    facet_grid(feature ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = paste("Odds ratios for KSI -", gname),
         subtitle = "Logistic regression on 2012-2021, 95% CI; each level vs its reference category",
         x = "Odds ratio (log scale)", y = NULL) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
          panel.grid.major.y = element_blank())
  ggsave(sprintf("output/models/odds_ratios_%s.png", gname), g, width = 8, height = 2.5 + 0.17 * nrow(o) + 0.12 * length(features), dpi = 150, limitsize = FALSE)
}

cat("\nSaved", length(list.files("output/models")), "files to output/models/\n")
