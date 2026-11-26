# 01_load_join.R
# Builds one person-level table from the Victorian Road Crash Data CSVs.
# Output: output/persons.rds (one row per person) and output/cleaning_log.csv
# (every step with the record counts it affected).
# All columns are read as character; type conversion and features come later.

library(data.table)

dir.create("output", showWarnings = FALSE)

rd <- function(f) {
  x <- fread(file.path("data", f), colClasses = "character", na.strings = "")
  setnames(x, gsub(" ", "_", names(x)))   # e.g. "LGA NAME ALL" -> LGA_NAME_ALL
  x
}

log <- list()
note <- function(step, n) {
  log[[length(log) + 1]] <<- data.table(step = step, records = n)
  cat(sprintf("%-75s %9s\n", step, format(n, big.mark = ",")))
}

# ---- Load -------------------------------------------------------------------
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

# ---- Node: one row per accident ---------------------------------------------
# Duplicates share NODE_ID and coordinates; only DEG_URBAN_NAME differs.
# Keep the first row in file order.
note("node accidents with >1 row (differ only in DEG_URBAN_NAME)",
     node[, .N, ACCIDENT_NO][N > 1, .N])
node <- unique(node, by = "ACCIDENT_NO")
node[, NODE_ID := NULL]                    # already in accident
note("node rows after de-duplication", nrow(node))

# ---- Weather and surface: one record per accident ---------------------------
# SEQ is 0 or 1..n; some accidents have both a 0 and a 1. Take the lowest
# SEQ >= 1 (the primary condition), else the SEQ 0 record.
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

# ---- Link people without a vehicle (mainly pedestrians) ---------------------
# 1) accident_event row naming the person: VEHICLE_1_ID of the first such event.
# 2) else, if the crash had exactly one vehicle, that vehicle.
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

# ---- Join -------------------------------------------------------------------
setnames(vehicle, "VEHICLE_ID", "LINK_VEHICLE_ID")
dt <- merge(person, accident, by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, vehicle, by = c("ACCIDENT_NO", "LINK_VEHICLE_ID"), all.x = TRUE)
dt <- merge(dt, node,    by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, atmos,   by = "ACCIDENT_NO", all.x = TRUE)
dt <- merge(dt, surface, by = "ACCIDENT_NO", all.x = TRUE)
stopifnot(nrow(dt) == nrow(person))        # no join may add or drop people

note("persons with no matching accident", dt[is.na(ACCIDENT_DATE), .N])
note("persons with a linked vehicle id not found in vehicle table",
     dt[!is.na(LINK_VEHICLE_ID) & is.na(VEHICLE_TYPE), .N])
note("persons with no node record", dt[is.na(LGA_NAME), .N])
note("persons with no atmospheric_cond record", dt[is.na(ATMOSPH_COND), .N])
note("persons with no road_surface_cond record", dt[is.na(SURFACE_COND), .N])

# ---- Target -----------------------------------------------------------------
# INJ_LEVEL 1 = fatality, 2 = serious injury
dt[, KSI := as.integer(INJ_LEVEL %in% c("1", "2"))]
note("KSI = 1 (killed or seriously injured)", dt[KSI == 1, .N])
note("final person-level rows", nrow(dt))

saveRDS(dt, "output/persons.rds")
fwrite(rbindlist(log), "output/cleaning_log.csv")
cat("\nSaved output/persons.rds (", nrow(dt), "rows x", ncol(dt), "cols) and output/cleaning_log.csv\n")
