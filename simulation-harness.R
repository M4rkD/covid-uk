# - - - - - - - - - - - - - - - - - - - - - - -
# UK model: load data and analyse scenarios
# - - - - - - - - - - - - - - - - - - - - - - -

Sys.setenv(CPATH="/home/mark/R/x86_64-pc-linux-gnu-library/3.6/nloptr/include/",
	   R_DATATABLE_NUM_THREADS="1")

library(rlang)
library(stringr)
library(rprojroot)

# Set path
# Set this path to the base directory of the repository.
covid_uk_path = getwd()

# covidm options
cm_path = paste0(covid_uk_path, "/covidm/");
source(paste0(cm_path, "/R/covidm.R"))

params = readRDS("test_sim.rds");
r = 0.3
message("start simulation...");

# default parameters
# params$date0 = "2020-01-29"
# params$time1 = "2021-12-31"

# test parameters
params$date0 = "2020-01-29"
params$time1 = "2020-03-06"

start_time = Sys.time()

run = cm_simulate(params, 1, r);

message("finished simulation");
message(sprintf("Simulation took %f seconds", Sys.time() - start_time));
