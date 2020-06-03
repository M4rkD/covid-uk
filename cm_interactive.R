library(stringr)
library(tidyverse)
library(furrr)
library(qs)

cm_path <- paste0(getwd(), "/covidm/")
source(paste0(cm_path, "/R/covidm.R"))

# Additional sources
source(paste0(cm_path, "/R/covidm_seeding.R"))

# Load global variables
cm_matrices <- readRDS(paste0(cm_path, "/data/all_matrices.rds"))
cm_populations <- readRDS(paste0(cm_path, "/data/wpp2019_pop2020.rds"))
cm_structure_UK <- readRDS(paste0(cm_path, "/data/structure_UK.rds"))
covid_scenario <- qread(file.path(cm_path, "../data/2-linelist_symp_fit_fIa0.5.qs"))

plan(multicore)
# plan(sequential)

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
  if (is.null(scenario_id)) {
    scenario_id <- sample.int(nrow(covid_scenario), 1)
  }

  # 1. Pick age-varying symptomatic rate
  covy <- unname(unlist(covid_scenario[scenario_id, f_00:f_70]))
  covy <- rep(covy, each = 2)

  # 2. Calculate R0 adjustment needed
  parametersUK1$pop[[1]]$y <- covy
  u_adj <- R0 / cm_calc_R0(parametersUK1, 1)

  list(covy = covy, u_adj = u_adj, r = R0)
}

with_u_adjustment <- function(params, adj) {
  for (i in seq_along(params$pop)) {
    params$pop[[i]]$u <- params$pop[[i]]$u * adj$u_adj

    params$pop[[i]]$y <- adj$covy
  }

  return(params)
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

classify_regions <- function() {
  locations <- cm_uk_locations("UK", 3)

  list(
    london = cm_structure_UK[match(str_sub(locations, 6), Name), Geography1 %like% "London"],
    england = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^E" & !(Geography1 %like% "London")],
    wales = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^W"],
    scotland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^S"],
    nireland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^N"],
    westmid = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "West Midlands (Met County)"],
    cumbria = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "Cumbria"]
  )
}

get_school_terms <- function() {
  list(
    close = c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25"),
    reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01")
  )
}

with_school_terms <- function(params, school_terms) {
  iv <- cm_iv_build(params)
  cm_iv_set(iv, school_terms$close, school_terms$reopen, contact = c(1, 1, 0, 1, 1, 1, 0, 1, 1), trace_school = 2)
  populations <- cm_iv_apply(params, iv)

  return(params)
}

# finally, I want to convert to seed day index (e.g. 1,1,1,1,2,2,2,...)
save_run_result <- function(R, ipop, run) {
  save_as <- paste0(store, "run-R=", R, "-pop=", ipop, ".rds")
  qsave(run, save_as)
}

create_gran_matrix <- function(pop) {
  "Create a grandparent matrix, by combining other matrices."

  if (length(pop$matrices) != 8) {
    stop("Current matrices must be of size 8")
  }

  # Recover home/other contact matrix
  mat_ref <- pop$matrices[[1]] + pop$matrices[[4]] +
    pop$matrices[[5]] + pop$matrices[[8]]

  gran <- 5 / 7 # adjustment for weekdays only.
  N <- nrow(mat_ref)
  popsize <- pop$size
  mat <- matrix(0, ncol = N, nrow = N)

  # Add child-grandparent contacts: under 15s to 55+s
  for (a in 1:3) {
    dist <- c(rep(0, 10 + a), mat_ref[a, (11 + a):N])
    dist <- dist / sum(dist)
    mat[a, ] <- mat[a, ] + gran * dist
    mat[, a] <- mat[, a] + (gran * dist) * (popsize[a] / popsize)
  }

  mat
}

build_params <- function(date_start, date_end, locations = NULL) {
  if (is.null(locations)) {
    locations <- cm_uk_locations("UK", 3)
  }

  cm_parameters_SEI3R(locations,
    date_start = date_start, date_end = date_end,
    dE = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # 6.5 day serial interval.
    dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p, # 1.5 days w/o symptoms
    dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness here as well.
    deterministic = F
  )
}

run_simulation <- function(params, run = 1, n = 1) {
  "Runs the simulation, setting the runtime."
  start_time <- Sys.time()
  result <- cm_simulate(params, run, 0)
  end_time <- Sys.time()
  result$runtime_seconds <- end_time - start_time
  result
}

params9 <- function(date_start, date_end, locations = NULL) {
  "Build parameters and split into 9 matrices"

  if (is.null(locations)) {
    locations <- cm_uk_locations("UK", 3)
  }

  parameters <- build_params(date_start, date_end, locations)

  # Split off the elderly (70+, age groups 15 and 16) so their contact matrices can be manipulated separately
  parameters <- cm_split_matrices_ex_in(parameters, 15)

  # Create additional matrix for child-elderly contacts
  parameters$pop <- map(parameters$pop, function(pop) {
    pop$matrices$gran <- create_gran_matrix(pop)
    pop
  })

  # add zero contact for gran matrix
  # grandparent matrix has no contact by default
  parameters$pop <- map(parameters$pop, function(p) {
    p$contact <- c(p$contact, 0)
    p
  })

  parameters
}

params_adj_base_with_seeds <- function(r, u_adj, seed_matrix, params, info = NULL) {
  "Run the base case with the specified seed matrix"

  params$pop <- seed_params(params$pop, seed_matrix, 4)

  # for some reason, this didn't work when replacing lapply with purr
  params$pop <- lapply(params$pop, function(p) {
    p$y <- u_adj$covy
    p
  })

  # apply covid scenario and R adjustment
  params$pop <- map(params$pop, function(p) {
    p$u <- p$u * u_adj$u_adj
    p
  })

  params
}

smooth_over_days <- function(curve, ts_per_day) {
  unlist(map(1:length(curve), rep, ts_per_day)) / 4
}

exp_seed_matrix <- function(seed, params, start, total, expgrowth, p_ht, ndays, ts_per_day) {
  set.seed(seed)

  seed_exp_growth_curve(start, total, expgrowth, p_ht, ndays) %>%
    smooth_over_days(ts_per_day) %>%
    weight_by_population(params$pop) %>%
    poisson_sample_matrix()
}

pop_set <- function(params, ...) {
  "Sets value in parameters for populations and returns the results.
Note, this does not mutate the input parameters."
  values <<- list(...)

  if (is.null(pop_idxs)) {
    pop_idxs <- 1:length(params$pop)
  }

  for (name in names(values)) {
    if (!name %in% names(params$pop[[1]])) {
      stop(paste0("Name ", name, " not found in params"))
    }

    for (i in pop_idxs) {
      params$pop[[i]][name] <- values[[name]]
    }
  }
  params
}

map_pop <- function(params, var, func) {
  params$pop <- map(params$pop, function(p) p[var] <- func(p[var]))

  return(params)
}

with_validate <- function(params) {
  if (length(map(params$pop, ~ .x$seed_times)) != length(params$pop)) {
    stop("Validation failed: Check the values for params$pop[[i]]$seed_times")
  }
  if (unique(dim(params$travel)) != length(params$pop)) {
    stop("Validation failed: Check the travel matrix size")
  }

  print("Validation passed")

  return(params)
}

with_seeding <- function(params, seed_matrices, translation = trans_seed_format, dist_ages = NULL) {
  "Set population seeding from seed_matrix into positions defined by idxs.

seed_matrix is a N-days x N-population matrix, with values specifying the number of seeds on that day."

  if (is.null(dist_ages)) {
    age_min <- 25
    age_max <- 50
    age_bounds <- 5 * 0:16

    dist_ages <- cm_age_coefficients(age_min, age_max, age_bounds)
  }

  params$pop <- lapply(params$pop, `[[<-`, "dist_seed_ages", dist_ages)

  ith <- function(ipop) {
    if (class(seed_matrices) == "matrix") {
      seed_matrices[, ipop]
    } else if (class(seed_matrices) == "function") {
      seed_matrices(params, ipop)
    } else {
      stop("Unknown type for seed matrix")
    }
  }

  for (ipop in seq_along(params$pop)) {
    params$pop[[ipop]]$seed_times <- translation(ith(ipop))
  }

  return(params)
}

with_store <- function(params) {
  .debug_params <<- params
}

with_standard_school_terms <- function(params) {
  with_school_terms(params, get_school_terms())
}

with_multiply_field <- function(params, fields, factor) {
  with_map_fields(params, fields, function(x) x * factor)
}

with_map_fields <- function(params, fields, func) {
  for (field in fields) {
    params <- with_map_field(params, field, func)
  }

  return(params)
}

with_map_field <- function(params, field, func) {
  params$pop <- map(params$pop, ~ modify_in(.x, field, func))

  return(params)
}


filter_pop <- function(params, pop_idxs) {
  params$pop <- params$pop[pop_idxs]
  params$travel <- params$travel[pop_idxs, pop_idxs]
}
