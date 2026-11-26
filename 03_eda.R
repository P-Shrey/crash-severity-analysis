# 03_eda.R
# Exploratory analysis: KSI rate by each factor, chi-square tests, charts.
# Input:  output/model_data.rds (from 02_clean_features.R)
# Output: output/eda/eda_rates.csv  KSI rate per factor level, all / occupant / vulnerable
#         output/eda/eda_chisq.csv  chi-square test + Cramer's V per factor and group
#         output/eda/*.png          one chart per factor, plus overview charts
# Uses all years (descriptive); models in 04 keep the 2022-2024 test set apart.
# With ~425k rows nearly every test is significant, so Cramer's V (effect size)
# is what ranks the factors.

library(data.table)
library(ggplot2)

dir.create("output/eda", recursive = TRUE, showWarnings = FALSE)
d <- readRDS("output/model_data.rds")

features <- c("sex", "age_group", "road_user", "protection", "vehicle_type", "vehicle_age",
              "speed_zone", "light", "weather", "surface", "road_geometry", "crash_type",
              "location", "hour_band", "weekend")
groups <- list(All = d, Occupant = d[user_group == "Occupant"], Vulnerable = d[user_group == "Vulnerable"])
col <- c(Occupant = "#2a78d6", Vulnerable = "#eb6834")   # categorical slots 1 and 2

# ---- Rates and chi-square ---------------------------------------------------
rates <- rbindlist(lapply(names(groups), function(g) rbindlist(lapply(features, function(x)
  groups[[g]][, .(n = .N, ksi = sum(KSI), rate = mean(KSI)), by = .(level = get(x))][
    , `:=`(group = g, feature = x)]))))
rates <- rates[n > 0][order(group, feature, level)]
fwrite(rates, "output/eda/eda_rates.csv")

chisq <- rbindlist(lapply(names(groups), function(g) rbindlist(lapply(features, function(x) {
  tab <- table(droplevels(groups[[g]][[x]]), groups[[g]]$KSI)
  if (nrow(tab) < 2) return(NULL)      # e.g. road_user has one level inside a group
  t <- suppressWarnings(chisq.test(tab))
  data.table(group = g, feature = x, chi_sq = unname(t$statistic), df = unname(t$parameter),
             p_value = t$p.value, cramers_v = sqrt(unname(t$statistic) / sum(tab)))
}))))
chisq <- chisq[order(group, -cramers_v)]
fwrite(chisq, "output/eda/eda_chisq.csv")
print(chisq[group == "All", .(feature, chi_sq = round(chi_sq), df, p_value = signif(p_value, 2),
                              cramers_v = round(cramers_v, 3))])

# ---- Charts -----------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                  legend.position = "top", legend.title = element_blank(),
                  plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA)))
overall <- mean(d$KSI)

# One dot plot per factor: KSI rate per level, occupants vs vulnerable road users
for (x in features) {
  p <- rates[feature == x & group != "All"]
  p[, level := factor(level, levels = rev(levels(d[[x]])))]
  g <- ggplot(p, aes(rate, level, colour = group)) +
    geom_vline(xintercept = overall, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_point(size = 3) +
    scale_colour_manual(values = col) +
    scale_x_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(title = paste("KSI rate by", gsub("_", " ", x)),
         subtitle = sprintf("Dashed line: overall KSI rate (%.1f%%)", 100 * overall),
         x = "Killed or seriously injured", y = NULL)
  ggsave(sprintf("output/eda/rate_%s.png", x), g, width = 7, height = 1.6 + 0.3 * nlevels(d[[x]]), dpi = 150)
}

# Factors ranked by effect size
v <- chisq[group == "All"][, feature := factor(feature, levels = rev(feature))]
g <- ggplot(v, aes(cramers_v, feature)) +
  geom_col(fill = col[["Occupant"]], width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", cramers_v)), hjust = -0.15, size = 3.2, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Association of each factor with KSI",
       subtitle = "Cramer's V from chi-square test, all persons 2012-2024", x = "Cramer's V", y = NULL)
ggsave("output/eda/cramers_v.png", g, width = 7, height = 5, dpi = 150)

# KSI rate by year
yr <- d[user_group != "Unknown", .(rate = mean(KSI)), by = .(year, user_group)]
g <- ggplot(yr, aes(year, rate, colour = user_group)) +
  geom_vline(xintercept = 2021.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  annotate("text", x = 2021.6, y = 0, label = "test set", hjust = 0, vjust = 0, size = 3, colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = col) +
  scale_x_continuous(breaks = 2012:2024) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "KSI rate by year", subtitle = "2024 covers January to September only",
       x = NULL, y = "Killed or seriously injured") +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line())
ggsave("output/eda/rate_by_year.png", g, width = 7, height = 4, dpi = 150)

cat("\nSaved", length(list.files("output/eda")), "files to output/eda/\n")
