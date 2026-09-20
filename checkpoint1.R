# BIOSTAT 707 Checkpoint 1: cohort characterization and EDA, set-a.
# Run with:  pixi run --locked checkpoint1
# No arguments. Writes everything to output/.

# library ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

library(tidyr)
library(tidyverse)


# root ----
root <- here::here()
data_dir <- file.path(root, "data")
out_dir  <- file.path(root, "output")
dir.create(out_dir, showWarnings = FALSE)
expected_n <- 4000

check_set_a <- function() {
  # data/set-a/ and data/Outcomes-a.txt must already be in place
  stopifnot("data/Outcomes-a.txt is missing" =
              file.exists(file.path(data_dir, "Outcomes-a.txt")))
  n <- length(list.files(file.path(data_dir, "set-a"), pattern = "\\.txt$"))
  stopifnot("data/set-a does not have 4000 records" = n == expected_n)
}

load_set_a <- function() {
  files <- list.files(file.path(data_dir, "set-a"), full.names = TRUE)
  raw <- rbindlist(lapply(files, function(f) {
    d <- fread(f)
    d[, RecordID := as.integer(tools::file_path_sans_ext(basename(f)))]
    d
  }))
  static_vars <- c("Age", "Gender", "Height", "ICUType", "Weight")
  static <- raw[Time == "00:00" & Parameter %in% static_vars] |>
    dcast(RecordID ~ Parameter, value.var = "Value", fun.aggregate = first)
  static[static == -1] <- NA            # -1 codes missing descriptors
  ts <- raw[!Parameter %in% c(static_vars, "RecordID")]
  list(ts = ts, static = static)
}

check_set_a()
d <- load_set_a()
outcomes <- fread(file.path(data_dir, "Outcomes-a.txt"))
stopifnot(nrow(outcomes) == expected_n)

# 1. Table 1                 -> output/table1.csv  (or .html from gtsummary)
# 2. Outcome summary         -> output/outcomes.csv
# 3. Missingness map         -> output/missingness_map.png
# 4. Missingness as signal   -> output/missingness_vs_death.csv
# ...

# Data wrangling ----

# Combine set-a time-series measurements into one long table
set_a_long <- d$ts

# Check the long table
head(set_a_long)
dim(set_a_long)
length(unique(set_a_long$RecordID))
unique(set_a_long$Parameter)

# Write giant skinny table
fwrite(
  set_a_long,
  file.path(out_dir, "set-a_long.csv")
)

# 开始数据初步处理，先看看数据概括，找到异常值 ----
measurement_summary <- set_a_long |>
  group_by(Parameter) |>
  summarise(
    n = n(),
    min = min(Value, na.rm = TRUE),
    p01 = quantile(Value, 0.01, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    p99 = quantile(Value, 0.99, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(measurement_summary, n = Inf)

# 第一步：把 Time 转成 ICU hour
set_a_long <- set_a_long |>
  separate(
    Time,
    into = c("hour", "minute"),
    sep = ":",
    remove = FALSE,
    convert = TRUE
  ) |>
  mutate(
    ICU_hour = hour + minute / 60,
    hour_bin = floor(ICU_hour)
  )
# 第二步：计算每个小时有多少病人测过某项指标
missingness_by_hour <- set_a_long |>
  filter(hour_bin >= 0, hour_bin < 48) |>
  distinct(RecordID, Parameter, hour_bin) |>
  count(Parameter, hour_bin, name = "n_measured") |>
  complete(
    Parameter,
    hour_bin = 0:47,
    fill = list(n_measured = 0)
  ) |>
  mutate(
    proportion_measured = n_measured / expected_n,
    proportion_missing = 1 - proportion_measured
  )
# 第三步：画 temporal missingness heatmap
ggplot(
  missingness_by_hour,
  aes(x = hour_bin, y = Parameter, fill = proportion_missing)
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Missing proportion",
    limits = c(0, 1)
  ) +
  labs(
    title = "Missingness by Measurement and ICU Hour",
    x = "ICU Hour",
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7)
  )

# long变wide ----
# Summarize each measurement over the 48-hour ICU period
# Build wide table ----

measurement_summary_wide <- set_a_long |>
  arrange(RecordID, Parameter, ICU_hour) |>
  group_by(RecordID, Parameter) |>
  summarise(
    count = sum(!is.na(Value)),
    first = first(Value[!is.na(Value)]),
    last  = last(Value[!is.na(Value)]),
    min   = min(Value, na.rm = TRUE),
    max   = max(Value, na.rm = TRUE),
    mean  = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to one row per admission
set_a_wide <- measurement_summary_wide |>
  pivot_wider(
    names_from = Parameter,
    values_from = c(count, first, last, min, max, mean),
    names_glue = "{Parameter}_{.value}"
  )

# Add static patient characteristics
set_a_wide <- d$static |>
  left_join(set_a_wide, by = "RecordID")

# Add outcomes
set_a_wide <- set_a_wide |>
  left_join(outcomes, by = "RecordID")

# Save wide table
fwrite(
  set_a_wide,
  file.path(out_dir, "set-a_wide.csv")
)

# Check
dim(set_a_wide)
head(set_a_wide)
stopifnot(nrow(set_a_wide) == expected_n)

# 画图 ----

# Missingness in wide table ----

wide_missingness <- set_a_wide |>
  select(ends_with("_mean")) |>
  summarise(
    across(
      everything(),
      ~ mean(is.na(.x))
    )
  ) |>
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Missing_proportion"
  ) |>
  mutate(
    Variable = sub("_mean$", "", Variable)
  ) |>
  arrange(desc(Missing_proportion))

wide_missingness

ggplot(
  wide_missingness,
  aes(
    x = reorder(Variable, Missing_proportion),
    y = Missing_proportion
  )
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    title = "Missingness of Measurements in Wide Table",
    x = "Measurement",
    y = "Missing proportion"
  ) +
  theme_minimal()



# table1.csv ----

# 1. Table 1 ----

table1_continuous <- d$static |>
  summarise(
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd = sd(Age, na.rm = TRUE),
    Age_median = median(Age, na.rm = TRUE),
    
    Height_mean = mean(Height, na.rm = TRUE),
    Height_sd = sd(Height, na.rm = TRUE),
    Height_median = median(Height, na.rm = TRUE),
    
    Weight_mean = mean(Weight, na.rm = TRUE),
    Weight_sd = sd(Weight, na.rm = TRUE),
    Weight_median = median(Weight, na.rm = TRUE)
  ) |>
  pivot_longer(
    everything(),
    names_to = "Characteristic",
    values_to = "Value"
  )

table1_gender <- d$static |>
  count(Gender) |>
  mutate(
    Characteristic = paste0("Gender_", Gender),
    Value = n
  ) |>
  select(Characteristic, Value)

table1_icu <- d$static |>
  count(ICUType) |>
  mutate(
    Characteristic = paste0("ICUType_", ICUType),
    Value = n
  ) |>
  select(Characteristic, Value)

table1 <- bind_rows(
  table1_continuous,
  table1_gender,
  table1_icu
)

write.csv(
  table1,
  file.path(out_dir, "table1.csv"),
  row.names = FALSE
)


# outcome.cvs ----
# 2. Outcome summary

outcome_summary <- outcomes |>
  summarise(
    N = n(),
    
    Death_n = sum(`In-hospital_death` == 1, na.rm = TRUE),
    Death_proportion = mean(`In-hospital_death` == 1, na.rm = TRUE),
    
    Length_of_stay_mean = mean(Length_of_stay, na.rm = TRUE),
    Length_of_stay_median = median(Length_of_stay, na.rm = TRUE),
    
    SAPS_I_mean = mean(`SAPS-I`, na.rm = TRUE),
    SAPS_I_median = median(`SAPS-I`, na.rm = TRUE),
    
    SOFA_mean = mean(SOFA, na.rm = TRUE),
    SOFA_median = median(SOFA, na.rm = TRUE)
  ) |>
  pivot_longer(
    everything(),
    names_to = "Outcome",
    values_to = "Value"
  )

write.csv(
  outcome_summary,
  file.path(out_dir, "outcomes.csv"),
  row.names = FALSE
)

# missingmap.png ----
# 3. Missingness map 

missingness_plot <- ggplot(
  missingness_by_hour,
  aes(
    x = hour_bin,
    y = Parameter,
    fill = proportion_missing
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Missing proportion",
    limits = c(0, 1)
  ) +
  labs(
    title = "Missingness by Measurement and ICU Hour",
    x = "ICU Hour",
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(out_dir, "missingness_map.png"),
  plot = missingness_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# missingness_vs_death.csv ----
# 4. Missingness as signal
# 4. Missingness as signal ----

missingness_vs_death <- set_a_wide |>
  select(
    RecordID,
    `In-hospital_death`,
    ends_with("_count")
  ) |>
  pivot_longer(
    cols = ends_with("_count"),
    names_to = "Parameter",
    values_to = "count"
  ) |>
  mutate(
    Parameter = sub("_count$", "", Parameter),
    count = replace_na(count, 0),
    measured = if_else(count > 0, "Measured", "Not measured")
  ) |>
  group_by(Parameter, measured) |>
  summarise(
    n = n(),
    deaths = sum(`In-hospital_death` == 1, na.rm = TRUE),
    death_rate = mean(`In-hospital_death` == 1, na.rm = TRUE),
    .groups = "drop"
  )

print(missingness_vs_death, n = Inf)

write.csv(
  missingness_vs_death,
  file.path(out_dir, "missingness_vs_death.csv"),
  row.names = FALSE
)

# 简单检查
# Check required outputs ----

required_outputs <- c(
  "set-a_long.csv",
  "set-a_wide.csv",
  "table1.csv",
  "outcomes.csv",
  "missingness_map.png",
  "missingness_vs_death.csv"
)
file.exists(file.path(out_dir, required_outputs))




message("Checkpoint 1 complete. See output/.")
