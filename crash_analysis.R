# crash_analysis.R  -  Pipeline A: crash severity analysis (scripts 01-04 in one file)
# Victorian Road Crash Data 2012 - Sep 2024: what raises the chance a person is
# killed or seriously injured (KSI), logistic regression vs random forest, and
# vulnerable road users vs vehicle occupants.
#
# Run from the project root (takes about 5 minutes):
#   Rscript crash_analysis.R
# Steps: 1 load + join   2 clean + features   3 EDA (chi-square, Cramer's V)   4 models
# Main outputs:
#   output/cleaning_log.csv, output/feature_log.csv        (Tables 2 and 4)
#   output/eda/eda_chisq.csv, cramers_v.png, rate_age_group.png, rate_vehicle_type.png, rate_by_year.png
#   output/models/metrics.csv, auc_tests.csv, odds_ratios.csv, roc.png, rf_importance.png,
#   odds_ratios_Occupant.png, odds_ratios_Vulnerable.png

library(data.table)
library(ggplot2)
library(ranger)
library(pROC)
library(caret)

dir.create("output/eda",    recursive = TRUE, showWarnings = FALSE)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)

log <- list()
note <- function(step, n) {
  log[[length(log) + 1]] <<- data.table(step = step, records = n)
  cat(sprintf("%-75s %9s\n", step, format(n, big.mark = ",")))
}

# =============================================================================
# 1. LOAD AND JOIN  ->  one row per person
# =============================================================================
cat("\n#### 1. Load and join ####\n")
rd <- function(f) {
  x <- fread(file.path("data", f), colClasses = "character", na.strings = "")
  setnames(x, gsub(" ", "_", names(x)))
  x
}
person   <- rd("person.csv")
accident <- rd("accident.csv")
vehicle  <- rd("vehicle.csv")
node     <- rd("node.csv")
atmos    <- rd("atmospheric_cond.csv")
surface  <- rd("road_surface_cond.csv")
event    <- rd("accident_event.csv")

note("person rows loaded", nrow(person))
note("accident rows loaded", nrow(accident))
note("vehicle rows loaded", nrow(vehicle))
note("node rows loaded", nrow(node))
note("atmospheric_cond rows loaded", nrow(atmos))
note("road_surface_cond rows loaded", nrow(surface))
note("accident_event rows loaded", nrow(event))

stopifnot(!anyDuplicated(person,   by = c("ACCIDENT_NO", "PERSON_ID")),
          !anyDuplicated(vehicle,  by = c("ACCIDENT_NO", "VEHICLE_ID")),
          !anyDuplicated(accident, by = "ACCIDENT_NO"))

# Node: duplicate rows per accident differ only in DEG_URBAN_NAME -> keep first
note("node accidents with >1 row (differ only in DEG_URBAN_NAME)", node[, .N, ACCIDENT_NO][N > 1, .N])
node <- unique(node, by = "ACCIDENT_NO")
node[, NODE_ID := NULL]
note("node rows after de-duplication", nrow(node))

# Weather / surface: several records per accident -> lowest SEQ >= 1, else SEQ 0
first_cond <- function(x, seq) {
  x[, s := as.integer(get(seq))]
  x <- x[order(ACCIDENT_NO, s == 0L, s)][, .SD[1], by = ACCIDENT_NO]
  x[, c("s", seq) := NULL]
  x
}
note("atmospheric_cond accidents with >1 record", atmos[, .N, ACCIDENT_NO][N > 1, .N])
atmos <- first_cond(atmos, "ATMOSPH_COND_SEQ")
note("road_surface_cond accidents with >1 record", surface[, .N, ACCIDENT_NO][N > 1, .N])
surface <- first_cond(surface, "SURFACE_COND_SEQ")

# People with no vehicle (pedestrians): link to the vehicle that struck them
#   a) accident_event row naming the person -> its VEHICLE_1_ID
#   b) else, if only one vehicle in the crash -> that vehicle
person[, LINK_VEHICLE_ID := VEHICLE_ID]
person[, VEHICLE_LINK := fifelse(is.na(VEHICLE_ID), "none", "own")]
note("persons with no VEHICLE_ID", person[is.na(VEHICLE_ID), .N])
note("  of which pedestrians (ROAD_USER_TYPE 1)", person[is.na(VEHICLE_ID) & ROAD_USER_TYPE == "1", .N])

ev <- event[!is.na(PERSON_ID) & !is.na(VEHICLE_1_ID)]
ev <- ev[order(ACCIDENT_NO, PERSON_ID, as.integer(EVENT_SEQ_NO))][, .SD[1], by = .(ACCIDENT_NO, PERSON_ID)]
person[ev, on = .(ACCIDENT_NO, PERSON_ID),
       `:=`(LINK_VEHICLE_ID = fifelse(VEHICLE_LINK == "none", i.VEHICLE_1_ID, LINK_VEHICLE_ID),
            VEHICLE_LINK    = fifelse(VEHICLE_LINK == "none", "striking_event", VEHICLE_LINK))]
note("  linked to striking vehicle via accident_event", person[VEHICLE_LINK == "striking_event", .N])

single <- vehicle[, .(n = .N, only = VEHICLE_ID[1]), by = ACCIDENT_NO][n == 1]
person[single, on = "ACCIDENT_NO",
       `:=`(LINK_VEHICLE_ID = fifelse(VEHICLE_LINK == "none", i.only, LINK_VEHICLE_ID),
            VEHICLE_LINK    = fifelse(VEHICLE_LINK == "none", "striking_single", VEHICLE_LINK))]
note("  linked to the only vehicle in a single-vehicle crash", person[VEHICLE_LINK == "striking_single", .N])
note("  still without a vehicle", person[VEHICLE_LINK == "none", .N])

# Join everything onto person (left joins; row count must not change)
setnames(vehicle, "VEHICLE_ID", "LINK_VEHICLE_ID")
dt <- merge(person, accident, by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, vehicle, by = c("ACCIDENT_NO", "LINK_VEHICLE_ID"), all.x = TRUE)
dt <- merge(dt, node,    by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, atmos,   by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, surface, by = "ACCIDENT_NO", all.x = TRUE)
stopifnot(nrow(dt) == nrow(person))
rm(person, accident, vehicle, node, atmos, surface, event, ev, single)

note("persons with no matching accident", dt[is.na(ACCIDENT_DATE), .N])
note("persons with a linked vehicle id not found in vehicle table",
     dt[!is.na(LINK_VEHICLE_ID) & is.na(VEHICLE_TYPE), .N])
note("persons with no node record", dt[is.na(LGA_NAME), .N])
note("persons with no atmospheric_cond record", dt[is.na(ATMOSPH_COND), .N])
note("persons with no road_surface_cond record", dt[is.na(SURFACE_COND), .N])

# Target: INJ_LEVEL 1 = fatality, 2 = serious injury
dt[, KSI := as.integer(INJ_LEVEL %in% c("1", "2"))]
note("KSI = 1 (killed or seriously injured)", dt[KSI == 1, .N])
note("final person-level rows", nrow(dt))
fwrite(rbindlist(log), "output/cleaning_log.csv")

# =============================================================================
# 2. CLEAN AND BUILD FEATURES  ("Not known" codes kept as an Unknown level)
#    First factor level = reference category in logistic regression
# =============================================================================
cat("\n#### 2. Clean and build features ####\n")
log <- list()
f <- function(x, levels) factor(x, levels = levels)

# ---- Person ----
note("SEX U or blank -> Unknown", dt[!SEX %in% c("M", "F"), .N])
dt[, sex := f(fcase(SEX == "M", "Male", SEX == "F", "Female", default = "Unknown"),
              c("Male", "Female", "Unknown"))]

note("AGE_GROUP Unknown or blank", dt[is.na(AGE_GROUP) | AGE_GROUP == "Unknown", .N])
dt[, age_group := f(fcase(AGE_GROUP %in% c("0-4", "5-12", "13-15"), "0-15",
                          AGE_GROUP %in% c("16-17", "18-21"),        "16-21",
                          AGE_GROUP %in% c("22-25", "26-29"),        "22-29",
                          AGE_GROUP %in% c("60-64", "65-69"),        "60-69",
                          AGE_GROUP %in% c("30-39", "40-49", "50-59", "70+"), AGE_GROUP,
                          default = "Unknown"),
                    c("30-39", "0-15", "16-21", "22-29", "40-49", "50-59", "60-69", "70+", "Unknown"))]

dt[, road_user := f(fcase(ROAD_USER_TYPE_DESC == "Drivers",       "Driver",
                          ROAD_USER_TYPE_DESC == "Passengers",    "Passenger",
                          ROAD_USER_TYPE_DESC == "Pedestrians",   "Pedestrian",
                          ROAD_USER_TYPE_DESC %in% c("Bicyclists", "E-scooter Rider"), "Bicyclist",
                          ROAD_USER_TYPE_DESC %in% c("Motorcyclists", "Pillion Passengers"), "Motorcyclist",
                          default = "Unknown"),
                    c("Driver", "Passenger", "Pedestrian", "Bicyclist", "Motorcyclist", "Unknown"))]
note("road user type Not Known", dt[road_user == "Unknown", .N])
dt[, user_group := f(fcase(road_user %in% c("Pedestrian", "Bicyclist", "Motorcyclist"), "Vulnerable",
                           road_user %in% c("Driver", "Passenger"), "Occupant",
                           default = "Unknown"),
                     c("Occupant", "Vulnerable", "Unknown"))]

# HELMET_BELT_WORN: 1,3,6 worn; 2,4,5,7 not worn/fitted; 8 not appropriate; 9 not known
note("HELMET_BELT_WORN 9 (not known) or blank -> Unknown",
     dt[is.na(HELMET_BELT_WORN) | HELMET_BELT_WORN == "9", .N])
dt[, protection := f(fcase(HELMET_BELT_WORN %in% c("1", "3", "6"), "Worn",
                           HELMET_BELT_WORN %in% c("2", "4", "5", "7"), "Not worn/fitted",
                           HELMET_BELT_WORN == "8", "Not applicable",
                           default = "Unknown"),
                     c("Worn", "Not worn/fitted", "Not applicable", "Unknown"))]

# ---- Vehicle (own vehicle; for pedestrians the striking vehicle) ----
dt[, vehicle_type := f(fcase(
  VEHICLE_TYPE_DESC %in% c("Car", "Station Wagon", "Taxi"), "Car",
  VEHICLE_TYPE_DESC %in% c("Utility", "Panel Van",
                           "Light Commercial Vehicle (Rigid) <= 4.5 Tonnes GVM"), "Ute/Van",
  grepl("Heavy Vehicle|Prime Mover|Rigid Truck", VEHICLE_TYPE_DESC), "Truck",
  VEHICLE_TYPE_DESC %in% c("Bus/Coach", "Mini Bus(9-13 seats)", "Tram", "Train"), "Bus/Tram/Train",
  VEHICLE_TYPE_DESC %in% c("Motor Cycle", "Motor Scooter", "Moped"), "Motorcycle",
  VEHICLE_TYPE_DESC %in% c("Bicycle", "Electric Device"), "Bicycle",
  is.na(VEHICLE_TYPE_DESC) | VEHICLE_TYPE_DESC %in% c("Not Known", "Not Applicable"), "Unknown",
  default = "Other"),
  c("Car", "Ute/Van", "Truck", "Bus/Tram/Train", "Motorcycle", "Bicycle", "Other", "Unknown"))]
note("vehicle type unknown (no vehicle, Not Known or Not Applicable)", dt[vehicle_type == "Unknown", .N])

dt[, year := as.integer(substr(ACCIDENT_DATE, 1, 4))]
manuf <- suppressWarnings(as.integer(dt$VEHICLE_YEAR_MANUF))
bad <- is.na(manuf) | manuf < 1950 | manuf > dt$year + 1      # 0, 1900, 3001 etc.
note("VEHICLE_YEAR_MANUF missing, 0, before 1950 or after crash year -> Unknown", sum(bad))
age <- fifelse(bad, NA_integer_, pmax(dt$year - manuf, 0L))
dt[, vehicle_age := f(fcase(is.na(age), "Unknown", age < 5, "0-4", age < 10, "5-9",
                            age < 15, "10-14", age < 20, "15-19", default = "20+"),
                      c("0-4", "5-9", "10-14", "15-19", "20+", "Unknown"))]

# ---- Crash ----
note("SPEED_ZONE 777, 888, 999 -> Unknown", dt[SPEED_ZONE %in% c("777", "888", "999"), .N])
sz <- suppressWarnings(as.integer(dt$SPEED_ZONE))
dt[, speed_zone := f(fcase(is.na(sz) | sz >= 777, "Unknown", sz <= 50, "<=50",
                           sz <= 70, "60-70", sz <= 90, "75-90", default = "100-110"),
                     c("60-70", "<=50", "75-90", "100-110", "Unknown"))]

note("LIGHT_CONDITION 9 (unknown)", dt[LIGHT_CONDITION == "9", .N])
dt[, light := f(fcase(LIGHT_CONDITION == "1", "Day",
                      LIGHT_CONDITION == "2", "Dusk/Dawn",
                      LIGHT_CONDITION == "3", "Dark, lit",
                      LIGHT_CONDITION %in% c("4", "5"), "Dark, unlit",
                      LIGHT_CONDITION == "6", "Dark, lighting unknown",
                      default = "Unknown"),
                c("Day", "Dusk/Dawn", "Dark, lit", "Dark, unlit", "Dark, lighting unknown", "Unknown"))]

note("ATMOSPH_COND Not known", dt[ATMOSPH_COND_DESC == "Not known", .N])
dt[, weather := f(fcase(ATMOSPH_COND_DESC == "Clear", "Clear",
                        ATMOSPH_COND_DESC == "Raining", "Rain",
                        ATMOSPH_COND_DESC == "Not known" | is.na(ATMOSPH_COND_DESC), "Unknown",
                        default = "Other"),
                  c("Clear", "Rain", "Other", "Unknown"))]

note("SURFACE_COND Unk.", dt[SURFACE_COND_DESC == "Unk.", .N])
dt[, surface := f(fcase(SURFACE_COND_DESC == "Dry", "Dry",
                        SURFACE_COND_DESC == "Wet", "Wet",
                        SURFACE_COND_DESC == "Unk." | is.na(SURFACE_COND_DESC), "Unknown",
                        default = "Other"),
                  c("Dry", "Wet", "Other", "Unknown"))]

dt[, road_geometry := f(fcase(ROAD_GEOMETRY_DESC == "Not at intersection", "Not at intersection",
                              ROAD_GEOMETRY_DESC == "Cross intersection", "Cross intersection",
                              ROAD_GEOMETRY_DESC == "T intersection", "T intersection",
                              ROAD_GEOMETRY_DESC == "Unknown" | is.na(ROAD_GEOMETRY_DESC), "Unknown",
                              default = "Other"),
                        c("Not at intersection", "Cross intersection", "T intersection", "Other", "Unknown"))]

dt[, crash_type := f(fcase(ACCIDENT_TYPE_DESC == "Collision with vehicle", "Vehicle",
                           ACCIDENT_TYPE_DESC == "Struck Pedestrian", "Pedestrian",
                           ACCIDENT_TYPE_DESC == "Collision with a fixed object", "Fixed object",
                           ACCIDENT_TYPE_DESC == "Vehicle overturned (no collision)", "Overturned",
                           ACCIDENT_TYPE_DESC == "No collision and no object struck", "No collision",
                           default = "Other"),
                     c("Vehicle", "Pedestrian", "Fixed object", "Overturned", "No collision", "Other"))]

note("DEG_URBAN_NAME missing -> Unknown", dt[is.na(DEG_URBAN_NAME), .N])
dt[, location := f(fcase(DEG_URBAN_NAME %in% c("MELB_URBAN", "MELBOURNE_CBD"), "Melbourne",
                         DEG_URBAN_NAME %in% c("LARGE_PROVINCIAL_CITIES", "SMALL_CITIES",
                                               "TOWNS", "SMALL_TOWNS"), "Regional town/city",
                         DEG_URBAN_NAME == "RURAL_VICTORIA", "Rural",
                         default = "Unknown"),
                   c("Melbourne", "Regional town/city", "Rural", "Unknown"))]

hr <- as.integer(substr(dt$ACCIDENT_TIME, 1, 2))
dt[, hour_band := f(fcase(hr < 6, "Night 00-05", hr < 10, "Morning 06-09", hr < 15, "Midday 10-14",
                          hr < 19, "Afternoon 15-18", default = "Evening 19-23"),
                    c("Midday 10-14", "Night 00-05", "Morning 06-09", "Afternoon 15-18", "Evening 19-23"))]
dt[, weekend := f(fifelse(DAY_WEEK_DESC %in% c("Saturday", "Sunday"), "Weekend", "Weekday"),
                  c("Weekday", "Weekend"))]

# ---- Time split: train 2012-2021, test 2022-2024 ----
dt[, split := f(fifelse(year <= 2021, "train", "test"), c("train", "test"))]
note("train rows (2012-2021)", dt[split == "train", .N])
note("test rows (2022-2024, 2024 ends 30 Sep)", dt[split == "test", .N])

# Leakage columns deliberately NOT used: SEVERITY, NO_PERSONS_*, INJ_LEVEL,
# TAKEN_HOSPITAL, EJECTED_CODE, LEVEL_OF_DAMAGE
features <- c("sex", "age_group", "road_user", "protection", "vehicle_type", "vehicle_age",
              "speed_zone", "light", "weather", "surface", "road_geometry", "crash_type",
              "location", "hour_band", "weekend")
d <- dt[, c("ACCIDENT_NO", "PERSON_ID", "year", "split", "KSI", "user_group", features), with = FALSE]
rm(dt)
stopifnot(nrow(d) == 425352L, !anyNA(d[, c("KSI", features), with = FALSE]))
note("final rows", nrow(d))
note("KSI rate x 1000 (train)", round(1000 * d[split == "train", mean(KSI)]))
note("KSI rate x 1000 (test)",  round(1000 * d[split == "test",  mean(KSI)]))
fwrite(rbindlist(log), "output/feature_log.csv")

groups <- list(All = d, Occupant = d[user_group == "Occupant"], Vulnerable = d[user_group == "Vulnerable"])

theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
                  plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA)))

# =============================================================================
# 3. EDA: KSI rate per level, chi-square + Cramer's V (effect size ranks factors,
#    because at n = 425k almost every test is significant)
# =============================================================================
cat("\n#### 3. EDA ####\n")
gcol <- c(Occupant = "#2a78d6", Vulnerable = "#eb6834")

rates <- rbindlist(lapply(names(groups), function(g) rbindlist(lapply(features, function(x)
  groups[[g]][, .(n = .N, rate = mean(KSI)), by = .(level = get(x))][, `:=`(group = g, feature = x)]))))

chisq <- rbindlist(lapply(names(groups), function(g) rbindlist(lapply(features, function(x) {
  tab <- table(droplevels(groups[[g]][[x]]), groups[[g]]$KSI)
  if (nrow(tab) < 2) return(NULL)                  # road_user has one level inside a group
  t <- suppressWarnings(chisq.test(tab))
  data.table(group = g, feature = x, chi_sq = unname(t$statistic), df = unname(t$parameter),
             p_value = t$p.value, cramers_v = sqrt(unname(t$statistic) / sum(tab)))
}))))
chisq <- chisq[order(group, -cramers_v)]
fwrite(chisq, "output/eda/eda_chisq.csv")
print(chisq[group == "All", .(feature, chi_sq = round(chi_sq), df, p_value = signif(p_value, 2),
                              cramers_v = round(cramers_v, 3))])

overall <- mean(d$KSI)
for (x in c("age_group", "vehicle_type")) {          # rate charts for age and vehicle type
  p <- rates[feature == x & group != "All"]
  p[, level := factor(level, levels = rev(levels(d[[x]])))]
  g <- ggplot(p, aes(rate, level, colour = group)) +
    geom_vline(xintercept = overall, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_point(size = 3) +
    scale_colour_manual(values = gcol) +
    scale_x_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(title = paste("KSI rate by", gsub("_", " ", x)),
         subtitle = sprintf("Dashed line: overall KSI rate (%.1f%%)", 100 * overall),
         x = "Killed or seriously injured", y = NULL) +
    theme(panel.grid.major.y = element_blank())
  ggsave(sprintf("output/eda/rate_%s.png", x), g, width = 7, height = 1.6 + 0.3 * nlevels(d[[x]]), dpi = 150)
}

v <- chisq[group == "All"][, feature := factor(feature, levels = rev(feature))]
g <- ggplot(v, aes(cramers_v, feature)) +
  geom_col(fill = gcol[["Occupant"]], width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", cramers_v)), hjust = -0.15, size = 3.2, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Association of each factor with KSI",
       subtitle = "Cramer's V from chi-square test, all persons 2012-2024", x = "Cramer's V", y = NULL) +
  theme(panel.grid.major.y = element_blank())
ggsave("output/eda/cramers_v.png", g, width = 7, height = 5, dpi = 150)

yr <- d[user_group != "Unknown", .(rate = mean(KSI)), by = .(year, user_group)]
g <- ggplot(yr, aes(year, rate, colour = user_group)) +
  geom_vline(xintercept = 2021.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  annotate("text", x = 2021.6, y = 0, label = "test set", hjust = 0, vjust = 0, size = 3, colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = gcol) +
  scale_x_continuous(breaks = 2012:2024) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "KSI rate by year", subtitle = "2024 covers January to September only",
       x = NULL, y = "Killed or seriously injured") +
  theme(panel.grid.major.x = element_blank())
ggsave("output/eda/rate_by_year.png", g, width = 7, height = 4, dpi = 150)

# =============================================================================
# 4. MODELS: logistic regression vs random forest, per group
#  - prediction models trained on a balanced sample (all KSI + equal non-KSI),
#    so threshold 0.5 is fair; test set never resampled
#  - odds ratios from logistic regression on the FULL training set
#  - DeLong test (two-sided): do RF and LR test AUCs differ?
# =============================================================================
cat("\n#### 4. Models (about 5 minutes) ####\n")
set.seed(42)
mcol <- c(Logistic = "#2a78d6", `Random forest` = "#eb6834")
metrics <- list(); auc_tests <- list(); ors <- list(); imps <- list(); rocs <- list()

for (g in names(groups)) {
  cat("\n====", g, "====\n")
  x <- copy(groups[[g]])
  x[, (features) := lapply(.SD, droplevels), .SDcols = features]
  fx <- features[sapply(features, function(f) nlevels(x[[f]]) > 1)]   # drop one-level factors
  form <- reformulate(fx, "KSI")
  train <- x[split == "train"]; test <- x[split == "test"]

  pos <- train[KSI == 1]
  bal <- rbind(pos, train[KSI == 0][sample(.N, nrow(pos))])
  cat(sprintf("train %s (KSI %s) -> balanced %s; test %s (KSI %.1f%%)\n",
              format(nrow(train), big.mark = ","), format(nrow(pos), big.mark = ","),
              format(nrow(bal), big.mark = ","), format(nrow(test), big.mark = ","), 100 * mean(test$KSI)))

  # Odds ratios: full training data, Wald 95% CI
  full <- glm(form, data = train, family = binomial)
  s <- summary(full)$coefficients
  ors[[g]] <- data.table(group = g, term = rownames(s), estimate = s[, 1], se = s[, 2], p_value = s[, 4])[
    , `:=`(odds_ratio = exp(estimate), ci_low = exp(estimate - 1.96 * se), ci_high = exp(estimate + 1.96 * se))][
    term != "(Intercept)"]

  # Prediction models on the balanced sample
  lr <- glm(form, data = bal, family = binomial)
  rf <- ranger(form, data = bal[, c("KSI", fx), with = FALSE][, KSI := factor(KSI)],
               num.trees = 300, probability = TRUE, importance = "permutation", seed = 42)
  p <- list(Logistic = predict(lr, test, type = "response"),
            `Random forest` = predict(rf, test)$predictions[, "1"])
  imps[[g]] <- data.table(group = g, feature = names(rf$variable.importance), importance = rf$variable.importance)

  # Evaluate on the untouched test set
  truth <- factor(test$KSI, levels = c(0, 1))
  r <- list()
  for (m in names(mcol)) {
    r[[m]] <- roc(test$KSI, p[[m]], levels = c(0, 1), direction = "<", quiet = TRUE)
    rocs[[paste(g, m)]] <- data.table(group = g, model = m, fpr = 1 - r[[m]]$specificities, tpr = r[[m]]$sensitivities)
    cm <- confusionMatrix(factor(as.integer(p[[m]] >= 0.5), levels = c(0, 1)), truth, positive = "1", mode = "everything")
    ci <- ci.auc(r[[m]])
    metrics[[paste(g, m)]] <- data.table(group = g, model = m, auc = as.numeric(auc(r[[m]])), auc_low = ci[1], auc_high = ci[3],
      accuracy = cm$overall[["Accuracy"]], recall = cm$byClass[["Recall"]], precision = cm$byClass[["Precision"]],
      f1 = cm$byClass[["F1"]], specificity = cm$byClass[["Specificity"]],
      tp = cm$table[2, 2], fp = cm$table[2, 1], fn = cm$table[1, 2], tn = cm$table[1, 1])
  }
  t <- roc.test(r[["Random forest"]], r[["Logistic"]], method = "delong", alternative = "two.sided")
  auc_tests[[g]] <- data.table(group = g, auc_rf = as.numeric(auc(r[["Random forest"]])),
                               auc_lr = as.numeric(auc(r[["Logistic"]])), z = unname(t$statistic), p_value = t$p.value)
}

metrics <- rbindlist(metrics); auc_tests <- rbindlist(auc_tests); ors <- rbindlist(ors); imps <- rbindlist(imps)
fwrite(metrics, "output/models/metrics.csv")
fwrite(auc_tests, "output/models/auc_tests.csv")
fwrite(ors, "output/models/odds_ratios.csv")

cat("\nTest-set metrics (threshold 0.5):\n")
print(metrics[, .(group, model, auc = round(auc, 3), accuracy = round(accuracy, 3), recall = round(recall, 3),
                  precision = round(precision, 3), f1 = round(f1, 3))])
cat("\nDeLong test, two-sided, H1: AUCs differ:\n")
print(auc_tests[, .(group, auc_rf = round(auc_rf, 4), auc_lr = round(auc_lr, 4), z = round(z, 2), p_value = signif(p_value, 3))])

# ---- Charts ----
roc_dt <- rbindlist(rocs)
roc_dt <- roc_dt[, .SD[unique(c(seq(1, .N, length.out = 500), .N))], by = .(group, model)]   # thin for plotting
lab <- metrics[, .(group, model, label = sprintf("%s AUC %.3f", model, auc))]
g <- ggplot(roc_dt, aes(fpr, tpr, colour = model)) +
  geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_text(data = lab, aes(x = 1, y = ifelse(model == "Logistic", 0.12, 0.04), label = label),
            hjust = 1, size = 3, show.legend = FALSE) +
  facet_wrap(~ group) + coord_equal() +
  scale_colour_manual(values = mcol) +
  labs(title = "ROC curves on the 2022-2024 test set", x = "False positive rate", y = "True positive rate")
ggsave("output/models/roc.png", g, width = 9, height = 4, dpi = 150)

imps[, feature := reorder(feature, importance)]
g <- ggplot(imps, aes(importance, feature)) +
  geom_col(fill = mcol[["Random forest"]], width = 0.6) +
  facet_wrap(~ group, scales = "free_x") +
  labs(title = "Random forest permutation importance", x = "Decrease in accuracy when permuted", y = NULL) +
  theme(panel.grid.major.y = element_blank())
ggsave("output/models/rf_importance.png", g, width = 9, height = 5, dpi = 150)

for (gname in c("Occupant", "Vulnerable")) {          # odds-ratio charts for the two road-user groups
  o <- ors[group == gname]
  o[, feature := sapply(term, function(t) features[startsWith(t, features)][1])]
  o[, level := substring(term, nchar(feature) + 1)]
  ref <- sapply(features, function(f) levels(droplevels(groups[[gname]][[f]]))[1])
  o[, feature := factor(sprintf("%s\n(ref: %s)", feature, ref[feature]),
                        levels = sprintf("%s\n(ref: %s)", features, ref))]
  o[, term := factor(term, levels = rev(term))]
  g <- ggplot(o, aes(odds_ratio, term)) +
    scale_y_discrete(labels = setNames(o$level, o$term)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", colour = mcol[["Logistic"]], linewidth = 0.6) +
    geom_point(colour = mcol[["Logistic"]], size = 2) +
    scale_x_log10() +
    facet_grid(feature ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = paste("Odds ratios for KSI -", gname),
         subtitle = "Logistic regression on 2012-2021, 95% CI; each level vs its reference category",
         x = "Odds ratio (log scale)", y = NULL) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
          panel.grid.major.y = element_blank())
  ggsave(sprintf("output/models/odds_ratios_%s.png", gname), g, width = 8,
         height = 2.5 + 0.17 * nrow(o) + 0.12 * length(features), dpi = 150, limitsize = FALSE)
}

cat("\nDone. Results in output/ (cleaning_log, feature_log, eda/, models/)\n")
