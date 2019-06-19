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

# Copeland's Mac Directory
#setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")

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
# fr$pheno$Low.Nitrogen <- factor(fr$pheno$Low.Nitrogen)
#fr$pheno$Pathogen <- factor(fr$pheno$Pathogen)
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
intOneVarNorm <- scanonevar(cross = fr,
                        mean.formula = Height ~ Env:(mean.QTL.add + mean.QTL.dom),
                        var.formula = ~ Env:(var.QTL.add + var.QTL.dom),
                        return.covar.effects = TRUE)
print("interactive NORM scanonevar")
# Writing the result of the interactive scanonevar for later use
write_rds(intOneVarNorm, "intOneVar_NORM.rds", compress = "xz")

# Writing out the results of the two 
#write.csv(addOneVar$result, file = "Manching_additive_model.csv")
write.csv(intOneVarNorm$result, file = "Manching_interactive_NORM_model.csv")

#plot(intOneVar, tests_to_plot = "mQTL", chrs = "1")



