# Script for parallel recombination, after a gnu parallel run

library(qs)           # for qsave and qread, faster equivalents of saveRDS and readRDS
library(data.table)   # for data.table, an enhanced (and faster) data.frame


covid_uk_path = "/home/mark/code/cmmid-covid-uk"

cm_path = paste0(covid_uk_path, "/covidm/");
source(paste0(cm_path, "/R/covidm_misc.R"))

load_single <- function(suffix, id, analysis) {
  "load an individual files from the option.single value given to that run."
  cm_load(paste0(covid_uk_path, analysis, "-", suffix, id, ".qs"));
}

parallel_read_with_suffix <- function(suffix, ids, analysis) {
  "load and combine all parallel fitotles"
  Reduce(function (accl, id) {
    single = load_single(suffix, id, analysis);
    single$run_id = id
    rbind(accl, single);
  },
  ids,
  data.table())
}

read_parallel_data <- function(ids, analysis) {
  list(
    dynamics = parallel_read_with_suffix("dynamics", ids, analysis),
    totals = parallel_read_with_suffix("totals", ids, analysis)
  )
}

write_parallel_data <- function(data, analysis) {
  cm_save(data$totals, paste0(covid_uk_path, analysis, "-totals", ".qs"));
  cm_save(data$dynamics, paste0(covid_uk_path, analysis, "-dynamics", ".qs"));
}

analysis = 1
ids = c(1,2)

data = read_parallel_data(ids, analysis)
write_parallel_data(data, analysis)
