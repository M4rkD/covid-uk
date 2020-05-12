#!/bin/bash

# Create a build directory in the RAM disk so run is independent
TMPDIR=$(mktemp --directory --tmpdir=/dev/shm)

# Don't leave a mess if terminated halfway through
term_handler() {
  echo Time limit is up; tidying up run $1 $2...
  rm -rf ${BUILD_DIR}
  exit -1
}
trap 'term_handler' TERM USR1

# Run R code
Rscript ../covid-uk/UK.R $1 $2 1

# Clean up
rm -rf ${BUILD_DIR}
