# Loading the packages
library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

# Setting working directories
# Comment out whichever one is not being used.

# Stampede2 Directory
setwd("/work/06566/tg858651/stampede2/Github/stapleton_lab/vQTL/")

# Michael's Mac Directory
#setwd("/Users/mbyrd/Stapleton/Stapleton_Lab/vQTL/Manching_Covariate/Interaction")

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

fr$pheno$Env <- factor(fr$pheno$Env)
print("before scanonevar")
# Additive scanonevar function
# addOneVar <- scanonevar(cross = fr,
#                         mean.formula = Height ~ Env + mean.QTL.add + mean.QTL.dom,
#                         var.formula = ~ Env + var.QTL.add + var.QTL.dom,
#                         return.covar.effects = TRUE)
# 
# # Writing the result of the additive scanonevar for later use
# write_rds(addOneVar, "addOneVar.rds", compress = "xz")
# print("first scanonevar")

# Interactive scanonevar function
#Two interactive models worth trying
#Height ~ Env + mean.QTL.add + mean.QTL.dom + Env * (mean.QTL.add + mean.QTL.dom)
# ~ Env + var.QTL.add + var.QTL.dom + Env * (var.QTL.add + var.QTL.dom)

#Height ~ Env + mean.QTL.add + mean.QTL.dom + (mean.QTL.add * mean.QTL.dom)
#~ Env + var.QTL.add + var.QTL.dom + (var.QTL.add * var.QTL.dom)
intOneVar <- scanonevar(cross = fr,
                        mean.formula = Height ~ Env + mean.QTL.add + mean.QTL.dom + (mean.QTL.add * mean.QTL.dom),
                        var.formula = ~ Env + var.QTL.add + var.QTL.dom + (var.QTL.add * var.QTL.dom),
                        return.covar.effects = TRUE)
print("second scanonevar")
# Writing the result of the interactive scanonevar for later use
write_rds(intOneVar, "intOneVar.rds", compress = "xz")


# Writing out the results of the two 
#write.csv(addOneVar$result, file = "Manching_additive_model.csv")
write.csv(intOneVar$result, file = "Manching_interactive_model2.csv")