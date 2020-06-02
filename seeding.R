# Exponential growth curve for seeding
seed_exp_growth_curve <- function(init, total, expgrowth, p_ht, ndays) {
  weights <- rep(0, ndays)
  weights[init] <- 1
  for (i in seq(init + 1, ndays)) {
    weights[i] <- weights[i - 1] * expgrowth
  }
  curve <- weights * total * (1 - p_ht) / sum(weights)
  curve[48] <- curve[48] + total * p_ht
  curve
}

weight_by_population <- function(seed_expected, populations) {
  "Return seeding matrix."
  pop_sizes <- calc_population_sizes(populations)

  normalised_pop <- pop_sizes / sum(pop_sizes)
  as.matrix(seed_expected) %*% t(as.matrix(normalised_pop))
}

plot_seed_matrix <- function(matrix) {
  "Plot a single seed matrix"
  df <- melt(seed_matrix_wales)

  names(df) <- c("time", "pop", "seeds")

  df$pop <- factor(locations[wales])

  ggplot(df, aes(x = time, y = seeds)) +
    geom_line(color = "blue") +
    facet_wrap(. ~ pop)
}

make_seed_format <- function(seeds) {
  "Convert a row of seeds (for a single population) to a the input expected by the R code."
  idx <- which(seeds > 0)
  rep <- seeds[idx]
  unlist(mapply(function(i, n) rep(i, n), idx, rep))
}

seed_params <- function(params, seed_matrix, idxs) {
  "Set population seeding from seed_matrix into positions defined by idxs.

seed_matrix is a N-days x N-population matrix, with values specifying the number of seeds on that day.

idxs are a mapping of indices of rows in seed_matrix to indices of populations in params. 1:N-pop by default."
  if (missing(idxs)) {
    idxs <- 1:length(params$pop)
  }

  if (is.logical(idxs)) {
    idxs <- which(idxs)
  }

  for (i in seq_along(length(idxs))) {
    i_pop <- idxs[i]
    params$pop[[i_pop]]$seed_times <- make_seed_format(seed_matrix[, i])
  }

  params
}

poisson_sample_matrix <- function(seed_matrix) {
  ## Try another way of sampling
  sampled <- apply(seed_matrix, 2, function(col) rpois(length(col), col))

  dim(sampled) <- dim(seed_matrix)

  sampled
}

plot_compare_seed_matrices_wales <- function(matrix1, matrix2) {
  "Plot two seeding matrices, both picking population labels for Wales"
  df <- merge(melt(matrix1), melt(matrix2), by = c("Var1", "Var2"))

  names(df) <- c("day", "population", "explicit expression", "calculated")

  ggdata <- melt(df, id.vars = c("day", "population"))

  welsh_locations <- locations[wales]
  ggdata$population <- factor(welsh_locations[ggdata$population])

  plt <- ggplot(ggdata, aes(x = day, y = value, color = variable)) +
    geom_line() +
    facet_wrap(. ~ population)
}
