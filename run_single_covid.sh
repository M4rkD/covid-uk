#!/bin/bash

# Process arguments
SEED_INDEX=$1
ANALYSIS=$2
NUMBER_OF_SEEDS=$3

# Create a build directory in the RAM disk so run is independent
export TMPDIR=$(mktemp --directory --tmpdir=/dev/shm)

# Don't leave a mess if terminated halfway through
term_handler() {
  echo "Time limit is up; tidying up run ${SEED_INDEX} ${ANALYSIS}..."
  rm -rf ${TMPDIR}
  exit -1
}
trap 'term_handler' TERM USR1

# Run R code
Rscript ../covid-uk/UK.R ${SEED_INDEX} ${ANALYSIS} ${NUMBER_OF_SEEDS}

# Clean up
rm -rf ${BUILD_DIR}
