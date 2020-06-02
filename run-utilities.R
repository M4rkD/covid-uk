library(stringr)
library(purrr)

Sys.setenv(CPATH="/home/mark/R/x86_64-pc-linux-gnu-library/3.6/nloptr/include/")

# Set up directories
covid_uk_path <- getwd()
cm_path = paste0(covid_uk_path, "/covidm/");

# Sourrce the files
source(file.path(cm_path, "R/covidm.R"))
source('seeding.R')

# Setup globals
covid_scenario <- qread(file.path(covid_uk_path, "data/2-linelist_symp_fit_fIa0.5.qs"))

# Function to compute the adjustment to R
pick_uadj <- function(R0, scenario_id = NULL) {

  parametersUK1 <- cm_parameters_SEI3R(cm_uk_locations("UK", 0),
    dE = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p,
    dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p,
    dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p,
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p,
    deterministic = F
  )

  # if scenario_id is NULL, then reset this
  if(is.null(scenario_id)) {
    scenario_id <- sample.int(nrow(covid_scenario), 1)
  }

  # 1. Pick age-varying symptomatic rate
  covy <- unname(unlist(covid_scenario[scenario_id, f_00:f_70]))
  covy <- rep(covy, each = 2)

  # 2. Calculate R0 adjustment needed
  parametersUK1$pop[[1]]$y <- covy
  u_adj <- R0 / cm_calc_R0(parametersUK1, 1)

  list(covy = covy, u_adj = u_adj)
}

apply_adjustments <- function(adj, population) {
  population$u =  population$u * adj$u_adj

  population$y = adj$covy;
}

calc_population_sizes <- function(populations) {
  # compute population sizes
  map_int(1:length(populations), ~ as.integer(sum(populations[[.x]]$size)))
}

get_population_sizes <- function(level) {
  "Function to compute population sizes.
Note: I expect this could be done in a far simpler way."
  locations <- cm_uk_locations("UK", 3)
  parameters <- cm_parameters_SEI3R(locations,
    dE = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # 6.5 day serial interval.
    dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p, # 1.5 days w/o symptoms
    dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness here as well.
    deterministic = F
  )

  # compute population sizes
  map_int(1:186, ~ as.integer(sum(parameters$pop[[.x]]$size)))
}

build_region_list <- function() {
  list(london = cm_structure_UK[match(str_sub(locations, 6), Name), Geography1 %like% "London"],
       england = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^E" & !(Geography1 %like% "London")],
       wales = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^W"],
       scotland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^S"],
       nireland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^N"],
       westmid = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "West Midlands (Met County)"],
       cumbria = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "Cumbria"])
}

# Get logicals for Wales
get_wales <- function() {
  wales <- build_region_list()$wales
}

get_school_terms <-function() {
  list(close =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25"),
       reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01"));
}

apply_school_terms_to_pops <- function(populations, params, school_terms) {
  iv = cm_iv_build(params)
  cm_iv_set(iv, school_terms$close, school_terms$reopen, contact = c(1, 1, 0, 1,  1, 1, 0, 1,  1), trace_school = 2);
  populations = cm_iv_apply(params, populations, iv);

  populations
}

# finally, I want to convert to seed day index (e.g. 1,1,1,1,2,2,2,...)
save_run_result <- function(R, ipop, run) {
  save_as <- paste0(store, "run-R=", R, "-pop=", ipop, ".rds")
  qsave(run, save_as)
}

create_gran_matrix <- function(parameters) {
  for (j in seq_along(parameters$pop))
  {
    # Recover home/other contact matrix
    mat_ref = parameters$pop[[j]]$matrices[[1]] + parameters$pop[[j]]$matrices[[4]] +
      parameters$pop[[j]]$matrices[[5]] + parameters$pop[[j]]$matrices[[8]];

    gran = 5/7; # adjustment for weekdays only.
    N = nrow(mat_ref);
    popsize = parameters$pop[[j]]$size;
    mat = matrix(0, ncol = N, nrow = N);

                                        # Add child-grandparent contacts: under 15s to 55+s
    for (a in 1:3) {
      dist = c(rep(0, 10 + a), mat_ref[a, (11 + a):N]);
      dist = dist/sum(dist);
      mat[a, ] = mat[a, ] + gran * dist;
      mat[, a] = mat[, a] + (gran * dist) * (popsize[a] / popsize);
    }
  }
}

build_params <- function(date_start, date_end) {
  locations = cm_uk_locations("UK", 3);
  parameters = cm_parameters_SEI3R(locations, date_start = date_start, date_end = date_end,
                                   dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # 6.5 day serial interval.
                                   dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p, # 1.5 days w/o symptoms
                                   dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness
                                   dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness here as well.
                                   deterministic = F);
  populations <- parameters$pop
  parameters["pop"] <- NULL

  list(pops = populations,
       params = parameters)
}

sim_runner <- function(params, r, N = 1) {
  "Runs the simulation, setting the runtime."
  print(paste0("Running with", r))
  start_time <- Sys.time()
  run <- cm_simulate(params, N, r)
  end_time <- Sys.time()
  run$runtime_seconds <- end_time - start_time
  run
}

params_all_matrices <- function(date_start, date_end) {
  parameters <- build_params(date_start, date_end)

  # Split off the elderly (70+, age groups 15 and 16) so their contact matrices can be manipulated separately
  parameters <- cm_split_matrices_ex_in(parameters, 15);

  # Create additional matrix for child-elderly contacts

  # Add child-grandparent contact matrix to population
  parameters$pop[[j]]$matrices$gran <- create_gran_matrix(parameters);

  # grandparent matrix has no contact by default
  parameters$pop[[j]]$contact <- c(parameters$pop[[j]]$contact, 0);

  parameters
}

run_adj_base_with_seeds <- function(r_and_adjustments, seed_matrix, params, info = NULL) {
  "Run the base case with the specified seed matrix"
  params <- seed_params(p, seed_matrix, wales)

  pop_set(params, NULL, "y", covy)
  pop_set(params, NULL, "u", function(u) u * r_and_adjustments$u_adj)

  r = r_and_adjustments$R
  r = r_and_adjustments$covy

  result <- sim_runner(p, r)
  result$info <- info
  save_result(result, R, ipop)
}

seed_matrix_generator <- function(seed, populations, start, total, expgrowth, p_ht, ndays) {
  set.seed(seed)

  poisson_sample_matrix(
    weight_by_population(
      seed_exp_growth_curve(start, total, expgrowth, p_ht, ndays),
      populations))
}

pop_set <- function(params, ...) {
  "Sets value in parameters for populations and returns the results.
Note, this does not mutate the input parameters."
  values <<- list(...)

  if(is.null(pop_idxs)) {
    pop_idxs = 1:length(params$pop)
  }

  for (name in names(values)) {
    if(!name %in% names(params$pop[[1]])) {
      stop(paste0("Name ", name," not found in params"))
    }

    for (i in pop_idxs) {
      params$pop[[i]][name] <- values[[name]]
    }
  }
  params
}

# build parameters
date_start="2020-05-20"
date_end="2020-11-20"

all <- build_params(date_start, date_end)

params <- all$params
pops <- all$pops[1]

pops <- apply_school_terms_to_pops(pops, params, school_terms)
params <- pop_set(pops, NULL, dist_seed_ages = cm_age_coefficients(25, 50, 5 * 0:16))
pops <- lapply(pops, `[[<-`, "dist_seed_ages", cm_age_coefficients(25, 50, 5 * 0:16))

# choose parameters
data_store <- "~/projects/cmmid-covid-uk/run_data/06-02-2020/"
R0s <- seq(2.0, 3.0, 0.1)
seed_matrices <- lapply(1:10, seed_matrix_generator, pops, 10, 500, 1.05, 0.5, 67)
# adjustments for each R

r_and_adjustments <- lapply(R0s, function(r) c(R=r, pick_uadj(r, scenario_id=1)))

map2(r_and_adjustments, seed_matrices, run_adj_base_with_seeds, params)

# short sim
params$time1 <- ymd(params$date0) + 10





for (seed in seq(40)) {
  # generate seeding data

  date_start="2020-05-20"
  date_end="2020-11-20" ...

  # generate parameters...
  pick_uadj
  apply_uadj



  map(R0, function(R) {
    params <- seed_params(params, seed_matrix_wales, wales)
    result <- sim_runner(params, r)
    save_result(result, R, ipop)
  })
}

# ggsave("expression-and-sampling-for-seed-without-peak.png", plt)
