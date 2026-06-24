###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

 # ---- R version check ----
 options(repos = c(CRAN = "https://cloud.r-project.org"))
 
 .libPaths(c("/software/projects/pawsey1163/johiolei/setonix/2025.08/r/4.4.1", .libPaths()))


## R versions for which this project has been tested and is stable
#stable_versions = c("4.4.0", "4.4.1")

## R versions for which this project is stable (as a string)
#stable_str = paste(stable_versions, collapse = ", ")

## Get details of R version currently running
#version_info = R.Version()

## Construct version number from list details
#version_num = paste0(version_info$major, ".",  version_info$minor)

## Throw an error if this R version is unsuitable
#if (!version_num %in% stable_versions)
#  stop("This software is stable with R version(s): ", stable_str,
#       " (currently running ", version_num, ")")

## Clear global environment
#rm(list = ls())


# functions, helper functions and parameters 
source("R/sub_functions.R")
 source("R/IBM_functions.R")
 
packages <- c("ggplot2", 
              "dplyr",
              "tibble",
              "tidyr", 
              "readr", 
              "purrr",       # uses pmap to loop through different scenarios
              "lhs",
              "doParallel",
              "parallel",
              "foreach",
              #"igraph",
              #"furrr",       # multisession i.e. distribute work across many cores
              #"progressr",
              "patchwork",
              "lme4",
              "lmerTest",
              #"sensitivity",
              #"rsm",
              #"randomForest",
              #"ranger",
              "truncnorm",
              "data.table")   # shows the progress


load_libraries(packages)

source("R/parameters.R")
