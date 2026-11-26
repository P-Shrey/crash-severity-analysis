# 02_clean_features.R
# Recodes "not known" codes, builds features and the train/test split.
# Input:  output/persons.rds (from 01_load_join.R)
# Output: output/model_data.rds (one row per person, factors ready for glm/ranger)
#         output/feature_log.csv (record counts behind every recode)
# Code meanings follow the VicRoads data dictionary. The first factor level is
# the reference category in logistic regression.

library(data.table)

dt <- readRDS("output/persons.rds")

log <- list()
note <- function(step, n) {
  log[[length(log) + 1]] <<- data.table(step = step, records = n)
  cat(sprintf("%-70s %9s\n", step, format(n, big.mark = ",")))
}
f <- function(x, levels) factor(x, levels = levels)

# ---- Person -----------------------------------------------------------------
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

# ROAD_USER_TYPE_DESC; e-scooter riders (24) grouped with bicyclists,
# pillion passengers with motorcyclists.
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

# HELMET_BELT_WORN: 1 belt, 3 child restraint, 6 helmet = worn;
# 2, 4, 7 = not worn; 5 = not fitted; 8 = not appropriate; 9 = not known
note("HELMET_BELT_WORN 9 (not known) or blank -> Unknown",
     dt[is.na(HELMET_BELT_WORN) | HELMET_BELT_WORN == "9", .N])
dt[, protection := f(fcase(HELMET_BELT_WORN %in% c("1", "3", "6"), "Worn",
                           HELMET_BELT_WORN %in% c("2", "4", "5", "7"), "Not worn/fitted",
                           HELMET_BELT_WORN == "8", "Not applicable",
                           default = "Unknown"),
                     c("Worn", "Not worn/fitted", "Not applicable", "Unknown"))]

# ---- Vehicle (own vehicle; for pedestrians the striking vehicle) ------------
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
bad <- is.na(manuf) | manuf < 1950 | manuf > dt$year + 1   # 0, 1900, 3001 etc.
note("VEHICLE_YEAR_MANUF missing, 0, before 1950 or after crash year -> Unknown", sum(bad))
dt[, vehicle_age_num := fifelse(bad, NA_integer_, pmax(year - manuf, 0L))]
dt[, vehicle_age := f(fcase(is.na(vehicle_age_num), "Unknown",
                            vehicle_age_num < 5,  "0-4",
                            vehicle_age_num < 10, "5-9",
                            vehicle_age_num < 15, "10-14",
                            vehicle_age_num < 20, "15-19",
                            default = "20+"),
                      c("0-4", "5-9", "10-14", "15-19", "20+", "Unknown"))]

# ---- Crash ------------------------------------------------------------------
# SPEED_ZONE 777 / 888 / 999 = camping ground or off road / other / not known
note("SPEED_ZONE 777, 888, 999 -> Unknown", dt[SPEED_ZONE %in% c("777", "888", "999"), .N])
sz <- suppressWarnings(as.integer(dt$SPEED_ZONE))
dt[, speed_zone := f(fcase(is.na(sz) | sz >= 777, "Unknown",
                           sz <= 50,  "<=50",
                           sz <= 70,  "60-70",
                           sz <= 90,  "75-90",
                           default = "100-110"),
                     c("60-70", "<=50", "75-90", "100-110", "Unknown"))]

# LIGHT_CONDITION 1 day, 2 dusk/dawn, 3 dark lit, 4 dark lights off,
# 5 dark no lights, 6 dark lighting unknown, 9 unknown
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
                        default = "Other"),             # fog, wind, dust, smoke, snow
                  c("Clear", "Rain", "Other", "Unknown"))]

note("SURFACE_COND Unk.", dt[SURFACE_COND_DESC == "Unk.", .N])
dt[, surface := f(fcase(SURFACE_COND_DESC == "Dry", "Dry",
                        SURFACE_COND_DESC == "Wet", "Wet",
                        SURFACE_COND_DESC == "Unk." | is.na(SURFACE_COND_DESC), "Unknown",
                        default = "Other"),             # icy, muddy, snowy
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

# DEG_URBAN_NAME; missing for persons whose accident has no node record
note("DEG_URBAN_NAME missing -> Unknown", dt[is.na(DEG_URBAN_NAME), .N])
dt[, location := f(fcase(DEG_URBAN_NAME %in% c("MELB_URBAN", "MELBOURNE_CBD"), "Melbourne",
                         DEG_URBAN_NAME %in% c("LARGE_PROVINCIAL_CITIES", "SMALL_CITIES",
                                               "TOWNS", "SMALL_TOWNS"), "Regional town/city",
                         DEG_URBAN_NAME == "RURAL_VICTORIA", "Rural",
                         default = "Unknown"),
                   c("Melbourne", "Regional town/city", "Rural", "Unknown"))]

hr <- as.integer(substr(dt$ACCIDENT_TIME, 1, 2))
dt[, hour_band := f(fcase(hr < 6,  "Night 00-05",
                          hr < 10, "Morning 06-09",
                          hr < 15, "Midday 10-14",
                          hr < 19, "Afternoon 15-18",
                          default = "Evening 19-23"),
                    c("Midday 10-14", "Night 00-05", "Morning 06-09", "Afternoon 15-18", "Evening 19-23"))]
dt[, weekend := f(fifelse(DAY_WEEK_DESC %in% c("Saturday", "Sunday"), "Weekend", "Weekday"),
                  c("Weekday", "Weekend"))]

# ---- Split ------------------------------------------------------------------
dt[, split := f(fifelse(year <= 2021, "train", "test"), c("train", "test"))]
note("train rows (2012-2021)", dt[split == "train", .N])
note("test rows (2022-2024, 2024 ends 30 Sep)", dt[split == "test", .N])

# ---- Keep model columns -----------------------------------------------------
# Left out on purpose (leak the outcome): SEVERITY, NO_PERSONS_KILLED/INJ_*,
# INJ_LEVEL, TAKEN_HOSPITAL, EJECTED_CODE, LEVEL_OF_DAMAGE.
features <- c("sex", "age_group", "road_user", "protection", "vehicle_type", "vehicle_age",
              "speed_zone", "light", "weather", "surface", "road_geometry", "crash_type",
              "location", "hour_band", "weekend")
out <- dt[, c("ACCIDENT_NO", "PERSON_ID", "year", "split", "KSI", "user_group",
              features, "vehicle_age_num", "VEHICLE_LINK", "LATITUDE", "LONGITUDE"), with = FALSE]
out[, `:=`(LATITUDE = as.numeric(LATITUDE), LONGITUDE = as.numeric(LONGITUDE))]

stopifnot(nrow(out) == 425352L, !anyNA(out[, c("KSI", features), with = FALSE]))
note("final rows", nrow(out))
note("KSI rate x 1000 (train)", round(1000 * out[split == "train", mean(KSI)]))
note("KSI rate x 1000 (test)",  round(1000 * out[split == "test",  mean(KSI)]))

saveRDS(out, "output/model_data.rds")
fwrite(rbindlist(log), "output/feature_log.csv")
cat("\nFeature levels:\n")
for (x in features) cat(sprintf("  %-14s %s\n", x, paste(levels(out[[x]]), collapse = " | ")))
cat("\nSaved output/model_data.rds (", nrow(out), "x", ncol(out), ") and output/feature_log.csv\n")
