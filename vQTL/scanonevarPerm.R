library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

#set directory and load scanonevar object
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")
#setwd("/work/04902/azg5169/stampede2/vQTL/vQTL/SOV.perm")
intResult = read_rds("InteractiveResult.rds")


#inputs: SOV object, number of permutations, random seed
#default inputs: n.cores = parallel::detectcores()-2
SOVperm = scanonevar.perm(intResult, n.perms = 4, random.seed = 3112020)

write_rds(SOVperm, "permResult1.rds")