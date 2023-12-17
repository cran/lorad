## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(lorad)

## -----------------------------------------------------------------------------
# Dimensions of k80samples
dim(k80samples)

# Column names of k80samples
colnames(k80samples)

## -----------------------------------------------------------------------------
# Create a named vector holding the column specifications
colspec <- c("iter"="iteration", "log.kernel"="posterior", "edgelen"="positive", "kappa"="positive")

## -----------------------------------------------------------------------------
results <- lorad_estimate(k80samples, colspec, 0.5, "left", 0.1)
lorad_summary(results)

