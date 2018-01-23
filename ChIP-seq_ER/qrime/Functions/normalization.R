library(dplyr)
library(preprocessCore)

# Functions for normalizing intensity data

# Performs scale normalization for the specified samples in the given MSnSet
# object such that the intensities have been scaled so that each sample has
# the same median value.

# Based on the limma implementation (normalizeMedianValues function) but with
# the option to scale to the median intensities calculated for specified rows
# (peptide measurements).

normalizeMedianScaling <- function(msnset, samples = NULL, rowsForCalculatingMedian = NULL)
{
  intensities <- msnset %>% exprs %>% as.data.frame

  if (is.null(samples))
    samples <- colnames(intensities)

  if (!all(samples %in% colnames(intensities)))
    stop("Could not find one or more of the given samples")

  if (length(samples) == 1) return(msnset)

  if (is.null(rowsForCalculatingMedian))
    rowsForCalculatingMedian <- 1:nrow(intensities)

  medianIntensities <- intensities %>%
    summarize_each(funs(median(.[rowsForCalculatingMedian], na.rm = TRUE)), one_of(samples)) %>%
    mutate_each(funs(log)) %>%
    as.numeric

  scalingFactors <- exp(medianIntensities - mean(medianIntensities))
  cat(scalingFactors, "\n")

  exprs(msnset)[,samples] <- t(t(intensities[,samples]) / scalingFactors)

  return(msnset)
}

