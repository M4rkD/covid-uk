# Script for parallel recombination, after a gnu parallel run

library(qs)           # for qsave and qread, faster equivalents of saveRDS and readRDS
library(data.table)   # for data.table, an enhanced (and faster) data.frame
library(rprojroot)

# Load requested settings from command line
argv = commandArgs(trailingOnly = T);
argc = length(argv);

analysis = argv[1]
last_id = as.numeric(argv[2])
ids = 1:last_id

covid_uk_path = normalizePath(dirname(thisfile()))

cm_path = file.path(covid_uk_path, "covidm");
source(file.path(cm_path, "R/covidm_misc.R"))

load_single <- function(suffix, id, analysis) {
  "load an individual files from the option.single value given to that run."
  cm_load(file.path(covid_uk_path, "output", paste0(analysis, "-", suffix, id, ".qs")));
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
  cm_save(data$totals, file.path(covid_uk_path, "output", paste0(analysis, "-totals", ".qs")));
  cm_save(data$dynamics, file.path(covid_uk_path, "output", paste0(analysis, "-dynamics", ".qs")));
}

data = read_parallel_data(ids, analysis)
write_parallel_data(data, analysis)
