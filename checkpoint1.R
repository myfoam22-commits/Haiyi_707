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



d$ts
d$static

unique(d$ts$Parameter)

summary(d$ts)

d$ts[, .(
  Min = min(Value, na.rm = TRUE),
  Max = max(Value, na.rm = TRUE),
  Mean = mean(Value, na.rm = TRUE),
  N = .N
), by = Parameter]

d$ts[Parameter == "Temp" & Value < 30]
d$ts[Parameter == "pH" & Value > 14]

d$ts[Parameter == "HR" & Value <= 0]
d$ts[Parameter == "HR" & Value > 250]

d$ts[Parameter == "RespRate" & Value <= 0]
d$ts[Parameter == "RespRate" & Value > 60]

d$ts[Parameter == "K" & Value > 10]
d$ts[Parameter == "Glucose" & Value > 1000]

d$ts[Parameter == "Na" & Value <= 0]
d$ts[Parameter == "Creatinine" & Value <= 0]

d$ts[Parameter == "MechVent", .N, by = Value]


d$ts[Parameter == "Temp" & Value < 30, .N]
d$ts[Parameter == "pH" & Value > 14, .N]
d$ts[Parameter == "HR" & Value <= 0, .N]
d$ts[Parameter == "RespRate" & Value <= 0, .N]


# 开始cleaning ----
ts_clean <- copy(d$ts)
ts_clean[Parameter == "Temp" & Value < 30, Value := NA]
ts_clean[Parameter == "pH" & Value > 14, Value := NA]
ts_clean[Parameter == "HR" & Value <= 0, Value := NA]
ts_clean[Parameter == "RespRate" & Value <= 0, Value := NA]

# ts_clean[Parameter == "Temp" & Value < 30, .N]
# ts_clean[Parameter == "pH" & Value > 14, .N]
# ts_clean[Parameter == "HR" & Value <= 0, .N]
# ts_clean[Parameter == "RespRate" & Value <= 0, .N]

ts_clean[, .(
  Missing = sum(is.na(Value))
), by = Parameter]
nrow(ts_clean)

ts_clean[, .(
  Missing = sum(is.na(Value)),
  Total = .N,
  Missing_pct = 100 * sum(is.na(Value)) / .N
), by = Parameter]

ts_clean[, .(
  Missing = sum(is.na(Value)),
  Total = .N,
  Missing_pct = 100 * sum(is.na(Value)) / .N
), by = Parameter][order(-Missing_pct)]

ts_clean[, .N, by = RecordID][, summary(N)]

patient_missing <- ts_clean[, .(
  Total = .N,
  Missing = sum(is.na(Value))
), by = RecordID]

patient_missing[, Missing_pct := 100 * Missing / Total]

summary(patient_missing$Missing_pct)

# 
presence <- ts_clean[, .(
  Present = 1L
), by = .(RecordID, Parameter)]

presence_summary <- presence[, .(
  Patients_measured = .N,
  Total_patients = 4000,
  Missing_patients = 4000 - .N,
  Missing_pct = 100 * (4000 - .N) / 4000
), by = Parameter]

presence_summary <- presence_summary[order(-Missing_pct)]

presence_summary

# 
head(outcomes)
names(outcomes)
table(outcomes[["In-hospital_death"]])
nrow(outcomes)
str(outcomes)





# data wrangling longer ----
# 1. 将 d$statistic 从宽表变成长表
d$statistic_long <- d$statistic %>%
  pivot_longer(
    cols = -RecordID,          # 除了 RecordID，其他列全部转换
    names_to = "Parameter",    # 列名（Age, Gender...）变成 Parameter 列的值
    values_to = "Value"        # 具体的数值变成 Value 列的值
  ) %>%
  mutate(
    Time = "00:00",            # 静态数据统一标记为 00:00
    Value = as.numeric(Value)  # 确保数值型统一（避免合并时报错）
  ) %>%
  select(RecordID, Time, Parameter, Value) # 调整列顺序，与 d$ts 完全一致

# 看一眼转换后的结果
head(d$statistic_long)
# 预期输出：
# # A tibble: 6 × 4
#   RecordID Time  Parameter Value
#      <dbl> <chr> <chr>     <dbl>
# 1   132539 00:00 Age          54
# 2   132539 00:00 Gender        0
# 3   132539 00:00 Height       NA
# 4   132539 00:00 ICUType       4
# 5   132539 00:00 Weight       NA
# ...

# 2. 将 d$ts 和转换后的 d$statistic_long 合并成巨型长表
set_a_long <- bind_rows(d$ts, d$statistic_long)

# 3. 排序（可选，但推荐。按患者ID、时间、参数排序）
set_a_long <- set_a_long %>%
  arrange(RecordID, Time, Parameter)

# 4. 输出到指定文件
if (!dir.exists(".output")) {
  dir.create(".output", recursive = TRUE)
}
write_csv(set_a_long, ".output/set-a_long.csv")

cat("合并完成！总行数：", nrow(set_a_long), "\n")

message("Checkpoint 1 complete. See output/.")
