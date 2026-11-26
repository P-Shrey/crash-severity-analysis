# =============================================================================
# crash_analysis_both_pipelines.R
# Victorian Road Crash Data, 2012 to Sep 2024
#   Pipeline A: 15 factors, balanced models, Cramer's V, DeLong
#   Pipeline B: 9 factors, clustered SEs, tuned threshold + RF,
#               calibration, 2024 sensitivity, VRU x age test
#
# Research questions
#   RQ1  Which crash and person factors raise the chance a person is killed or
#        seriously injured (KSI = INJ_LEVEL 1 fatal or 2 serious)?
#   RQ2  How well can KSI be predicted, and does a random forest beat
#        logistic regression?
#   RQ3  Do the risk factors differ for vulnerable road users (pedestrians,
#        cyclists, motorcyclists) vs vehicle occupants?
#
# Run from the project root:   Rscript crash_analysis_both_pipelines.R
#   about 8 minutes; set TUNE <- TRUE to also redo Pipeline B's tuning (+10 min)
#
# Sections
#   1  Load + join            (shared: one row per person, 425,352 rows)
#   2  Features               (A: 15 factors   B: 9 factors)
#   3  EDA                    (A: KSI rates, chi-square, Cramer's V)
#   4  Pipeline A models      (LR vs RF on balanced train, DeLong test)
#   5  Pipeline B explanatory (clustered odds ratios, without-2024 check)
#   6  Pipeline B prediction  (fit/val/test, tuned threshold + RF, calibration)
#   7  Pipeline B RQ3         (VRU vs occupant OR, age interaction Wald test)
#
# Main outputs
#   cleaning_log.csv, feature_log.csv, models/odds_ratios.csv, pipeline_b/b_odds_ratios.csv
#   models/metrics.csv, pipeline_b/*, figures in output/eda and output/models
# =============================================================================

library(data.table)   # fast reading, joins and grouping
library(ggplot2)      # charts
library(ranger)       # fast random forest
library(pROC)         # ROC, AUC with CI, DeLong test

TUNE <- FALSE         # TRUE = rerun threshold search + RF grid on validation years

dir.create("output/eda",        recursive = TRUE, showWarnings = FALSE)
dir.create("output/models",     recursive = TRUE, showWarnings = FALSE)
dir.create("output/pipeline_b", recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ---------------------------------------------------------------

# note(): print a count and keep it for the cleaning / feature log tables
log <- list()
note <- function(step, n) {
  log[[length(log) + 1]] <<- data.table(step = step, records = n)
  cat(sprintf("%-75s %9s\n", step, format(n, big.mark = ",")))
}

# eval_at(): confusion-matrix metrics at one probability threshold
eval_at <- function(p, y, thr) {
  pred <- as.integer(p >= thr)
  tp <- sum(pred == 1 & y == 1); fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1); tn <- sum(pred == 0 & y == 0)
  precision <- tp / (tp + fp); recall <- tp / (tp + fn)
  data.table(threshold = thr, accuracy = (tp + tn) / length(y), recall = recall,
             precision = precision, f1 = 2 * precision * recall / (precision + recall),
             specificity = tn / (tn + fp), brier = mean((p - y)^2),
             tp = tp, fp = fp, fn = fn, tn = tn)
}

# cluster_vcov(): crash-clustered robust covariance for a logistic glm.
# People in the same crash share speed, road and weather, so they are not
# independent; ordinary SEs would be too small. Same result as
# sandwich::vcovCL(fit, cluster, type = "HC1") but needs no extra package.
cluster_vcov <- function(fit, cluster) {
  keep <- !is.na(coef(fit))                          # drop aliased terms
  X <- model.matrix(fit)[, keep, drop = FALSE]
  U <- rowsum(X * (fit$y - fitted(fit)), cluster)    # score summed per crash
  B <- vcov(fit, complete = FALSE)                   # (X'WX)^-1, the "bread"
  G <- nrow(U); n <- nrow(X); k <- ncol(X)
  G / (G - 1) * (n - 1) / (n - k) * B %*% crossprod(U) %*% B
}

# save_show(): save a chart to PNG, and also draw it in the
# RStudio Plots pane when run interactively (use the arrows there to browse)
save_show <- function(file, plot, ...) {
  ggsave(file, plot, ...)
  if (interactive()) print(plot)
}

# or_table(): odds ratios with 95% CI, from a coefficient vector and its covariance
or_table <- function(beta, V) {
  se <- sqrt(diag(V))[names(beta)]
  data.table(term = names(beta), odds_ratio = exp(beta),
             ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
             p_value = 2 * pnorm(-abs(beta / se)))
}

# =============================================================================
# 1. LOAD AND JOIN  ->  one row per person (both pipelines start here)
# =============================================================================
cat("\n#### 1. Load and join ####\n")
rd <- function(f) {                         # read everything as text: keeps codes like "09" intact
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

# Primary keys must be unique, otherwise joins would duplicate people
stopifnot(!anyDuplicated(person,   by = c("ACCIDENT_NO", "PERSON_ID")),
          !anyDuplicated(vehicle,  by = c("ACCIDENT_NO", "VEHICLE_ID")),
          !anyDuplicated(accident, by = "ACCIDENT_NO"))

# Node: some accidents have 2 rows that differ only in DEG_URBAN_NAME.
#   Pipeline A keeps the first row; Pipeline B labels them "Ambiguous".
node_b <- node[, .(Urban_Classification = if (uniqueN(DEG_URBAN_NAME) == 1) DEG_URBAN_NAME[1] else "Ambiguous"),
               by = .(ACCIDENT_NO, NODE_ID)]
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

# People with no vehicle (mostly pedestrians): link them to the vehicle that hit them
#   a) accident_event row naming the person -> its VEHICLE_1_ID
#   b) else, if the crash had only one vehicle -> that vehicle
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

# Left joins onto person; the row count must never change
setnames(vehicle, "VEHICLE_ID", "LINK_VEHICLE_ID")
dt <- merge(person, accident, by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, vehicle, by = c("ACCIDENT_NO", "LINK_VEHICLE_ID"), all.x = TRUE)
dt <- merge(dt, node,    by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, node_b,  by = c("ACCIDENT_NO", "NODE_ID"), all.x = TRUE)     # Pipeline B's geography
dt <- merge(dt, atmos,   by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, surface, by = "ACCIDENT_NO", all.x = TRUE)
stopifnot(nrow(dt) == nrow(person))
rm(person, accident, vehicle, node, node_b, atmos, surface, event, ev, single)

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
# 2. FEATURES
#    "Not known" codes are kept as an Unknown level, not deleted: dropping
#    115k people with unknown helmet/belt would bias the data.
#    First factor level = reference category in logistic regression.
#    Leakage columns never used: SEVERITY, NO_PERSONS_*, INJ_LEVEL,
#    TAKEN_HOSPITAL, EJECTED_CODE, LEVEL_OF_DAMAGE (known only after the crash)
# =============================================================================
cat("\n#### 2. Features ####\n")
log <- list()
f <- function(x, levels) factor(x, levels = levels)
dt[, year := as.integer(substr(ACCIDENT_DATE, 1, 4))]

# ---- 2A. Pipeline A: 15 grouped factors ----
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

# Vehicle = own vehicle; for pedestrians/cyclists the vehicle that struck them
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

manuf <- suppressWarnings(as.integer(dt$VEHICLE_YEAR_MANUF))
bad <- is.na(manuf) | manuf < 1950 | manuf > dt$year + 1      # 0, 1900, 3001 etc.
note("VEHICLE_YEAR_MANUF missing, 0, before 1950 or after crash year -> Unknown", sum(bad))
age <- fifelse(bad, NA_integer_, pmax(dt$year - manuf, 0L))
dt[, vehicle_age := f(fcase(is.na(age), "Unknown", age < 5, "0-4", age < 10, "5-9",
                            age < 15, "10-14", age < 20, "15-19", default = "20+"),
                      c("0-4", "5-9", "10-14", "15-19", "20+", "Unknown"))]

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

# Time split: train 2012-2021, test 2022-2024 (future data = honest test)
dt[, split := f(fifelse(year <= 2021, "train", "test"), c("train", "test"))]
note("train rows (2012-2021)", dt[split == "train", .N])
note("test rows (2022-2024, 2024 ends 30 Sep)", dt[split == "test", .N])

features <- c("sex", "age_group", "road_user", "protection", "vehicle_type", "vehicle_age",
              "speed_zone", "light", "weather", "surface", "road_geometry", "crash_type",
              "location", "hour_band", "weekend")

# ---- 2B. Pipeline B: 9 factors, close to the raw codes ----
# Blank -> "Unknown"; 777/888/999 speed codes kept as their own "Special_" level
unk <- function(x) fifelse(is.na(x) | trimws(x) == "", "Unknown", x)
dt[, `:=`(
  Road_User_Clean = fcase(ROAD_USER_TYPE %in% c("2", "7"), "Driver",
                          ROAD_USER_TYPE %in% c("3", "8"), "Passenger",
                          ROAD_USER_TYPE == "1",  "Pedestrian",
                          ROAD_USER_TYPE == "4",  "Motorcyclist",
                          ROAD_USER_TYPE == "5",  "Pillion Passenger",
                          ROAD_USER_TYPE == "6",  "Bicyclist",
                          ROAD_USER_TYPE == "16", "E-scooter Rider",
                          default = "Unknown"),
  Age_Group_Clean      = unk(AGE_GROUP),
  Sex_Clean            = fifelse(SEX %in% c("M", "F"), SEX, "Unknown"),
  Speed_Zone_Clean     = fifelse(SPEED_ZONE %in% c("777", "888", "999"), paste0("Special_", SPEED_ZONE),
                                 unk(as.character(as.integer(SPEED_ZONE)))),   # "060" -> "60"
  LIGHT_CONDITION      = paste0("Light_", LIGHT_CONDITION),
  ROAD_GEOMETRY_DESC   = unk(ROAD_GEOMETRY_DESC),
  ATMOSPH_COND_DESC    = unk(ATMOSPH_COND_DESC),
  SURFACE_COND_DESC    = unk(SURFACE_COND_DESC),
  Urban_Classification = unk(Urban_Classification)
)]
dt[, Road_User_Group := fcase(
  Road_User_Clean %in% c("Pedestrian", "Bicyclist", "Motorcyclist", "Pillion Passenger", "E-scooter Rider"), "Vulnerable Road User",
  Road_User_Clean %in% c("Driver", "Passenger"), "Vehicle Occupant",
  default = "Unknown")]
b_features <- c("Road_User_Clean", "Age_Group_Clean", "Sex_Clean", "Speed_Zone_Clean", "LIGHT_CONDITION",
                "ROAD_GEOMETRY_DESC", "ATMOSPH_COND_DESC", "SURFACE_COND_DESC", "Urban_Classification")

d <- dt[, c("ACCIDENT_NO", "PERSON_ID", "year", "split", "KSI", "user_group", "Road_User_Group",
            features, b_features), with = FALSE]
rm(dt)
stopifnot(nrow(d) == 425352L, !anyNA(d[, c("KSI", features, b_features), with = FALSE]))
note("final rows", nrow(d))
note("KSI rate x 1000 (train)", round(1000 * d[split == "train", mean(KSI)]))
note("KSI rate x 1000 (test)",  round(1000 * d[split == "test",  mean(KSI)]))
fwrite(rbindlist(log), "output/feature_log.csv")

groups <- list(All = d, Occupant = d[user_group == "Occupant"], Vulnerable = d[user_group == "Vulnerable"])

theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
                  plot.title.position = "plot", plot.background = element_rect(fill = "white", colour = NA)))

# =============================================================================
# 3. EDA (Pipeline A): KSI rate per level, chi-square + Cramer's V.
#    At n = 425k every chi-square test is significant, so factors are ranked
#    by effect size (Cramer's V = sqrt(chi2 / n) for a 2-column table).
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

# Figures: KSI rate by age group and by vehicle type, occupants vs vulnerable
overall <- mean(d$KSI)
for (x in c("age_group", "vehicle_type")) {
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
  save_show(sprintf("output/eda/rate_%s.png", x), g, width = 7, height = 1.6 + 0.3 * nlevels(d[[x]]), dpi = 150)
}

# Figure: Cramer's V ranking
v <- chisq[group == "All"][, feature := factor(feature, levels = rev(feature))]
g <- ggplot(v, aes(cramers_v, feature)) +
  geom_col(fill = gcol[["Occupant"]], width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", cramers_v)), hjust = -0.15, size = 3.2, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Association of each factor with KSI",
       subtitle = "Cramer's V from chi-square test, all persons 2012-2024", x = "Cramer's V", y = NULL) +
  theme(panel.grid.major.y = element_blank())
save_show("output/eda/cramers_v.png", g, width = 7, height = 5, dpi = 150)

# Figure: KSI rate by year (shows the drop after 2021 = train/test drift)
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
save_show("output/eda/rate_by_year.png", g, width = 7, height = 4, dpi = 150)

# =============================================================================
# 4. PIPELINE A MODELS: logistic regression vs random forest, per group
#  - prediction models trained on a BALANCED sample (all KSI + equal non-KSI),
#    so threshold 0.5 is fair; the test set is never resampled
#  - odds ratios from logistic regression on the FULL training set
#  - DeLong test (two-sided): do the RF and LR test AUCs differ?
# =============================================================================
cat("\n#### 4. Pipeline A models (about 5 minutes) ####\n")
set.seed(42)
mcol <- c(Logistic = "#2a78d6", `Random forest` = "#eb6834")
metrics <- list(); auc_tests <- list(); ors <- list(); imps <- list(); rocs <- list(); a_calib <- list()

for (g in names(groups)) {
  cat("\n====", g, "====\n")
  x <- copy(groups[[g]])
  x[, (features) := lapply(.SD, droplevels), .SDcols = features]
  fx <- features[sapply(features, function(f) nlevels(x[[f]]) > 1)]   # drop one-level factors (road_user in subgroups)
  form <- reformulate(fx, "KSI")
  train <- x[split == "train"]; test <- x[split == "test"]

  pos <- train[KSI == 1]
  bal <- rbind(pos, train[KSI == 0][sample(.N, nrow(pos))])           # undersample non-KSI to 50/50
  cat(sprintf("train %s (KSI %s) -> balanced %s; test %s (KSI %.1f%%)\n",
              format(nrow(train), big.mark = ","), format(nrow(pos), big.mark = ","),
              format(nrow(bal), big.mark = ","), format(nrow(test), big.mark = ","), 100 * mean(test$KSI)))

  # Odds ratios: full training data, Wald 95% CI (+ OR figures)
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

  # Evaluate on the untouched 2022-2024 test set
  r <- list()
  for (m in names(mcol)) {
    r[[m]] <- roc(test$KSI, p[[m]], levels = c(0, 1), direction = "<", quiet = TRUE)
    rocs[[paste(g, m)]] <- data.table(group = g, model = m, fpr = 1 - r[[m]]$specificities, tpr = r[[m]]$sensitivities)
    ci <- ci.auc(r[[m]])
    metrics[[paste(g, m)]] <- cbind(data.table(group = g, model = m, auc = as.numeric(auc(r[[m]])),
                                               auc_low = ci[1], auc_high = ci[3]),
                                    eval_at(p[[m]], test$KSI, 0.5))
    a_calib[[paste(g, m)]] <- data.table(group = g, model = m, mean_predicted = mean(p[[m]]), observed = mean(test$KSI))
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
cat("\nCalibration of balanced models (mean predicted vs observed KSI) - why we report ranking, not probabilities:\n")
print(rbindlist(a_calib)[group == "All", .(model, mean_predicted = round(mean_predicted, 3), observed = round(observed, 3))])

# Figure: ROC curves
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
save_show("output/models/roc.png", g, width = 9, height = 4, dpi = 150)

# Figure: RF permutation importance
imps[, feature := reorder(feature, importance)]
g <- ggplot(imps, aes(importance, feature)) +
  geom_col(fill = mcol[["Random forest"]], width = 0.6) +
  facet_wrap(~ group, scales = "free_x") +
  labs(title = "Random forest permutation importance", x = "Decrease in accuracy when permuted", y = NULL) +
  theme(panel.grid.major.y = element_blank())
save_show("output/models/rf_importance.png", g, width = 9, height = 5, dpi = 150)

# Figures: odds-ratio forest plots for occupants and vulnerable road users (RQ3)
for (gname in c("Occupant", "Vulnerable")) {
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
  save_show(sprintf("output/models/odds_ratios_%s.png", gname), g, width = 8,
         height = 2.5 + 0.17 * nrow(o) + 0.12 * length(features), dpi = 150, limitsize = FALSE)
}

# =============================================================================
# 5. PIPELINE B EXPLANATORY MODEL (RQ1 robustness)
#  - all years, crash year as a factor (adjusts for the trend over time)
#  - crash type NOT included, so road-user ORs are not split with it
#    (Pipeline A's pedestrian OR of 32 comes from that overlap)
#  - crash-clustered 95% CIs
# =============================================================================
cat("\n#### 5. Pipeline B explanatory model ####\n")
reg <- copy(d[, c("ACCIDENT_NO", "year", "KSI", b_features), with = FALSE])
# Merge very rare levels so every coefficient can be estimated
reg[Road_User_Clean == "E-scooter Rider", Road_User_Clean := "Other vulnerable road user"]
reg[Speed_Zone_Clean == "75", Speed_Zone_Clean := "70"]
reg[ROAD_GEOMETRY_DESC %in% c("Private property", "Road closure"), ROAD_GEOMETRY_DESC := "Other geometry"]
reg[ATMOSPH_COND_DESC == "Snowing", ATMOSPH_COND_DESC := "Other weather"]
# Reference categories
refs <- c(Road_User_Clean = "Driver", Age_Group_Clean = "30-39", Sex_Clean = "F", Speed_Zone_Clean = "60",
          LIGHT_CONDITION = "Light_1", ROAD_GEOMETRY_DESC = "Not at intersection", ATMOSPH_COND_DESC = "Clear",
          SURFACE_COND_DESC = "Dry", Urban_Classification = "MELB_URBAN")
for (v in names(refs)) set(reg, j = v, value = relevel(factor(reg[[v]]), ref = refs[[v]]))
reg[, Crash_Year := factor(year)]

b_form <- reformulate(c(b_features, "Crash_Year"), "KSI")
b_glm  <- glm(b_form, data = reg, family = binomial)
b_or   <- or_table(coef(b_glm)[-1], cluster_vcov(b_glm, reg$ACCIDENT_NO))
fwrite(b_or, "output/pipeline_b/b_odds_ratios.csv")
cat("Key clustered odds ratios:\n")
print(b_or[term %in% c("Road_User_CleanMotorcyclist", "Road_User_CleanPedestrian", "Road_User_CleanBicyclist",
                       "Age_Group_Clean70+", "Speed_Zone_Clean100", "Urban_ClassificationRURAL_VICTORIA"),
           .(term, odds_ratio = round(odds_ratio, 2), ci_low = round(ci_low, 2), ci_high = round(ci_high, 2))])

# Sensitivity: refit without the partial year 2024 - do the ORs move?
reg24  <- reg[year < 2024][, Crash_Year := droplevels(Crash_Year)]
b_glm24 <- glm(b_form, data = reg24, family = binomial)
sens <- merge(b_or[, .(term, full_or = odds_ratio)],
              data.table(term = names(coef(b_glm24))[-1], no2024_or = exp(coef(b_glm24)[-1])), by = "term")
sens[, pct_change := round(100 * (no2024_or / full_or - 1), 2)]
fwrite(sens, "output/pipeline_b/b_sensitivity_no2024.csv")
cat("Without 2024, main ORs change by:\n")
print(sens[term %in% c("Road_User_CleanPedestrian", "Age_Group_Clean70+", "Speed_Zone_Clean100",
                       "Urban_ClassificationRURAL_VICTORIA"), .(term, pct_change)])
rm(reg24, b_glm24)

# =============================================================================
# 6. PIPELINE B PREDICTION (RQ2 robustness)
#  - three-way time split: fit 2012-18 | validation 2019-21 | test 2022-24
#  - full (unbalanced) data; threshold and RF settings chosen on validation
#    years only, then models refit on 2012-21 and scored once on the test set
#  - category levels learned from training data; unseen test values -> Unknown
# =============================================================================
cat("\n#### 6. Pipeline B prediction ####\n")
encode <- function(train, test) {           # factor levels from train only (no peeking at test)
  for (v in b_features) {
    lv <- sort(unique(train[[v]]))
    tv <- test[[v]]; tv[!tv %in% lv] <- "Unknown"
    set(train, j = v, value = factor(train[[v]], levels = lv))
    set(test,  j = v, value = factor(tv, levels = lv))
  }
}
b_fit_models <- function(train, mtry, node) {
  list(lr = glm(reformulate(b_features, "KSI"), data = train, family = binomial),
       rf = ranger(reformulate(b_features, "y"), data = train[, c(b_features, "KSI"), with = FALSE][, y := factor(KSI)][, KSI := NULL],
                   num.trees = 300, mtry = mtry, min.node.size = node, probability = TRUE,
                   respect.unordered.factors = "order", seed = 2026))
}
b_predict <- function(m, test) list(lr = predict(m$lr, test, type = "response"), rf = predict(m$rf, test)$predictions[, "1"])

thr <- c(lr = 0.21, rf = 0.22); best <- c(mtry = 2, node = 50)   # results of the tuning below
if (TUNE) {
  fit <- copy(d[year <= 2018]); val <- copy(d[year >= 2019 & split == "train"])
  encode(fit, val)
  f1_best <- function(p) { r <- rbindlist(lapply(seq(0.05, 0.95, 0.01), eval_at, p = p, y = val$KSI)); r[which.max(f1)] }
  # RF grid: 6 settings, pick highest validation AUC
  grid <- CJ(mtry = c(2, 3, 5), node = c(20, 50))
  grid[, auc := mapply(function(m, n) as.numeric(auc(val$KSI, b_predict(b_fit_models(fit, m, n), val)$rf, quiet = TRUE)), mtry, node)]
  print(grid[order(-auc)])
  best <- unlist(grid[which.max(auc), .(mtry, node)])
  # Thresholds: highest F1 on validation years
  pv <- b_predict(b_fit_models(fit, best[["mtry"]], best[["node"]]), val)
  thr <- c(lr = f1_best(pv$lr)$threshold, rf = f1_best(pv$rf)$threshold)
  cat("Chosen: mtry", best[["mtry"]], "min.node.size", best[["node"]], "thresholds", thr, "\n")
  rm(fit, val, pv)
}

# Final refit on 2012-2021, one evaluation on 2022-2024
btrain <- copy(d[split == "train"]); btest <- copy(d[split == "test"])
encode(btrain, btest)
bm <- b_fit_models(btrain, best[["mtry"]], best[["node"]])
bp <- b_predict(bm, btest)
b_metrics <- rbind(
  cbind(model = "Logistic", auc = as.numeric(auc(btest$KSI, bp$lr, quiet = TRUE)), eval_at(bp$lr, btest$KSI, thr[["lr"]])),
  cbind(model = "Tuned RF", auc = as.numeric(auc(btest$KSI, bp$rf, quiet = TRUE)), eval_at(bp$rf, btest$KSI, thr[["rf"]])))
fwrite(b_metrics, "output/pipeline_b/b_test_metrics.csv")
cat("Pipeline B test metrics (tuned thresholds):\n")
print(b_metrics[, .(model, auc = round(auc, 3), threshold, recall = round(recall, 3), precision = round(precision, 3),
                    f1 = round(f1, 3), brier = round(brier, 4))])

# Calibration: are predicted probabilities real risks? 10 equal-size risk bands
calib <- rbind(data.table(model = "Logistic", p = bp$lr, y = btest$KSI),
               data.table(model = "Tuned RF", p = bp$rf, y = btest$KSI))
calib[, band := ceiling(10 * frank(p, ties.method = "first") / .N), by = model]
cal <- calib[, .(mean_predicted = mean(p), observed = mean(y)), by = .(model, band)]
fwrite(cal, "output/pipeline_b/b_calibration_deciles.csv")
cat("Mean predicted vs observed KSI on test:\n")
print(calib[, .(mean_predicted = round(mean(p), 3), observed = round(mean(y), 3)), by = model])
g <- ggplot(cal, aes(mean_predicted, observed, colour = model)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_line() + geom_point(size = 2) + coord_equal() +
  scale_colour_manual(values = c(Logistic = "#2a78d6", `Tuned RF` = "#eb6834")) +
  labs(title = "Calibration on 2022-2024 test set (Pipeline B)",
       subtitle = "Below the diagonal = predicted risk too high", x = "Mean predicted KSI", y = "Observed KSI rate")
save_show("output/pipeline_b/b_calibration.png", g, width = 6, height = 6, dpi = 150)
rm(btrain, btest, bm, bp, calib)

# =============================================================================
# 7. PIPELINE B RQ3: one model for occupants + vulnerable road users
#  - adjusted OR of KSI for vulnerable vs occupant
#  - group x age interaction, joint clustered Wald test: does the gap depend on age?
# =============================================================================
cat("\n#### 7. Pipeline B vulnerable road users ####\n")
vru <- copy(d[Road_User_Group != "Unknown", c("ACCIDENT_NO", "year", "KSI", "Road_User_Group", b_features[-1]), with = FALSE])
for (v in c("Road_User_Group", b_features[-1])) set(vru, j = v, value = factor(vru[[v]]))
vru[, Road_User_Group := relevel(Road_User_Group, ref = "Vehicle Occupant")]
vru[, Age_Group_Clean := relevel(Age_Group_Clean, ref = "30-39")]
vru[, Crash_Year := factor(year)]
adj <- c(b_features[-(1:2)], "Crash_Year")              # all B factors except road user and age
vterm <- "Road_User_GroupVulnerable Road User"

# Main effect model
m1 <- glm(reformulate(c("Road_User_Group", "Age_Group_Clean", adj), "KSI"), data = vru, family = binomial)
vru_or <- or_table(coef(m1)[vterm], cluster_vcov(m1, vru$ACCIDENT_NO))
cat("Vulnerable vs occupant, adjusted OR:\n"); print(vru_or)

# Interaction model + joint Wald test of all group x age terms (known ages only)
vk <- vru[Age_Group_Clean != "Unknown"][, Age_Group_Clean := droplevels(Age_Group_Clean)]
m2 <- glm(reformulate(c("Road_User_Group * Age_Group_Clean", adj), "KSI"), data = vk, family = binomial)
V2 <- cluster_vcov(m2, vk$ACCIDENT_NO)
b2 <- coef(m2)[!is.na(coef(m2))]
it <- grep(":Age_Group_Clean", names(b2), value = TRUE)
wald <- as.numeric(t(b2[it]) %*% solve(V2[it, it], b2[it]))
cat(sprintf("Joint Wald test, group x age: chi2 = %.1f, df = %d, p = %.2g\n",
            wald, length(it), pchisq(wald, length(it), lower.tail = FALSE)))

# VRU vs occupant OR at each age = main effect + that age's interaction term
ages <- c("0-4", "5-12", "13-15", "16-17", "18-21", "22-25", "26-29", "30-39",
          "40-49", "50-59", "60-64", "65-69", "70+")
by_age <- rbindlist(lapply(ages, function(a) {
  k <- c(vterm, if (a != "30-39") paste0(vterm, ":Age_Group_Clean", a))
  k <- k[k %in% names(b2)]
  est <- sum(b2[k]); se <- sqrt(sum(V2[k, k]))
  data.table(age_group = a, odds_ratio = exp(est), ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se))
}))
fwrite(by_age, "output/pipeline_b/b_vru_or_by_age.csv")
print(by_age[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 2) else z)])

cat("\nDone. Pipeline A -> output/ (cleaning_log, feature_log, eda/, models/); Pipeline B -> output/pipeline_b/\n")
