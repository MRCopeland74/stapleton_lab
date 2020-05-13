library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")
Maindata = read.csv("ManchingStressData_Covar.csv")

write_csv(sub_data, "ManchingSigData.csv")
vQTLsub <- read.cross(file = "ManchingStressData_Covar.csv" )

vQTLsub <- drop.nullmarkers(vQTLsub)
vQTLsub <- calc.genoprob(vQTLsub)

#sub.scan <- scanonevar(cross = vQTLsub,
#                       mean.formula = ï..Height ~ (Low.Water + Low.Nitrogen + Pathogen)*(mean.QTL.add),
#                       var.formula = ~ (Low.Water + Low.Nitrogen + Pathogen)*(var.QTL.add),
#                       return.covar.effects = TRUE)

#intOneVar <- scanonevar.per(cross = vQTLsub, 
#                            mean.formula = Height ~ (Low.Water*(mean.QTL.add + mean.QTL.dom) + Low.Nitrogen*(mean.QTL.add + mean.QTL.dom) + Pathogen*(mean.QTL.add + mean.QTL.dom)),
#                            var.formula = ~ (Low.Water*(var.QTL.add + var.QTL.dom) + Low.Nitrogen*(var.QTL.add + var.QTL.dom) + Pathogen*(var.QTL.add + var.QTL.dom)),
#                            return.covar.effects = TRUE)
#scanonevar
##INTERACTIVE
#with ManchingStressData_covar.csv
intRealOneVar <- scanonevar(cross = vQTLsub, 
                            mean.formula = Height ~ Env*(mean.QTL.add + mean.QTL.dom),
                            var.formula = ~ Env*(var.QTL.add + var.QTL.dom),
                            return.covar.effects = TRUE)
table(intRealOneVar$result$mvQTL.asymp.p <= .05)

intRealOneVar$result$mvQTL.lod[1:50]
intOneResult = intRealOneVar$result

write_csv(intOneResult, "int_Results.csv")

intOneSig = intOneResult[which(intRealOneVar$result$mvQTL.asymp.p <= .05),]

write_csv(intOneSig, "int_SigPvals.csv")

intOnlyPval = intOneResult[,c(1:4,6,8,10)]

write_csv(intOnlyPval, "int_OnlyPvals.csv")

write_rds(intRealOneVar, "InteractiveResult.rds")

hyv_p1$result$loc.name[hyv_p1$result$mvQTL.asymp.p <= .05]
intOneVar$result$loc.name[intOneVar$result$mvQTL.asymp.p <= .05]

table(intOneVar$result$mvQTL.asymp.p <= .05)

intOneVar$result$mvQTL.lod[1:50]




write_rds(intOneVar, "intOneVarResults.rds")

intResults = intOneVar$result