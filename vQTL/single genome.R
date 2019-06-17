# Loading the packages Byrd code
library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

# Setting working directories
# Comment out whichever one is not being used.

# Stampede2 Directory
# setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")

# Michael's Mac Directory
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")

# Reading in the input file as a 'cross' object

data <- read.csv(file = "ManchingStressData_Covar.csv")

# Created a random sample
# set.seed(1234)
# subset = data[c(1,2, round(runif(100,1,6674))),]
# write.csv(subset, file = "Manching_Sample.csv")


# Full Data Set, Comment if using Sample
# fr <-read.cross(file = "ManchingStressData_Covar.csv")

# Sample of Data Set, Comment if using full data
fr <-read.cross(file = "ManchingStressData_Covar.csv")

# Not sure what these two functions do yet.
fr <- drop.nullmarkers(fr)
#scan with variance
fr <- calc.genoprob(fr)


# This is the function Thomas used to create the function he labeled as "effect sizes"
# The code for this function can be found here,
# https://github.com/cran/vqtl/blob/master/R/mean_var_plots.R
# Here the focal.groups parameter is the different SNPs
mean_var_plot_model_based(cross = fr,
                          phenotype.name = "Height",
                          genotype.names = c("AA","BB"),
                          focal.groups = 'gpm113b')

