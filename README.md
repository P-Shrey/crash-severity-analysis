# Victorian Road Crash Severity Analysis (R)

Analysis of Victorian road crash data, 2012 to Sep 2024 (425,352 people involved in crashes), in R.

## Research questions
1. Which crash and person factors raise the chance a person is killed or seriously injured (KSI)?
2. How well can KSI be predicted, and does a random forest beat logistic regression?
3. Do risk factors differ for vulnerable road users (pedestrians, cyclists, motorcyclists) vs vehicle occupants?

## Key results (test set 2022-2024)
| Group | Logistic regression AUC | Random forest AUC |
|---|---|---|
| All | 0.785 | 0.787 |
| Occupant | 0.759 | 0.756 |
| Vulnerable | 0.676 | 0.668 |

Random forest gives no meaningful gain; logistic regression is as accurate and easier to interpret.

## Data
Not included (about 170 MB). Download the Victoria Road Crash Data CSVs from
[Transport Victoria Open Data](https://opendata.transport.vic.gov.au/dataset/victoria-road-crash-data)
and place them in `data/` (`accident.csv`, `accident_event.csv`, `accident_location.csv`,
`atmospheric_cond.csv`, `node.csv`, `person.csv`, `road_surface_cond.csv`, `sub_dca.csv`, `vehicle.csv`).

## Run
Requires R 4.x with packages `data.table`, `ggplot2`, `ranger`, `pROC`, `caret`:
```
install.packages(c("data.table", "ggplot2", "ranger", "pROC", "caret"))
```

From the project root, either run the step-by-step pipeline:
```
Rscript 01_load_join.R      # join tables, one row per person
Rscript 02_clean_features.R # features, train 2012-2021 / test 2022-2024
Rscript 03_eda.R            # KSI rates, chi-square, Cramer's V
Rscript 04_models.R         # LR vs RF, odds ratios, DeLong AUC test (~5 min)
```
or the single-file versions:
- `crash_analysis.R` - Pipeline A (scripts 01-04 merged)
- `crash_analysis_both_pipelines.R` - Pipeline A plus condensed Pipeline B (clustered SEs, tuned threshold, calibration)

Results are written to `output/` (EDA plots, model metrics, odds ratios, ROC curves).
