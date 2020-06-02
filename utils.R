unwhich <- function(x, n) {
  "Does the inverse operation of which(), given a length"
  out <- rep_len(FALSE, n)
  out[x] <- TRUE
  out
}
