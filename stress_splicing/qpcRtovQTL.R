library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)
library(glm.predict)
library(plyr)
library(qpcR)

### READ IN DERIVATIVE DATA ###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
#deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
#deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
#deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:

setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/qPCR2vQTL")
####################################NOVEMBER#####################################################
plate.1 = read.csv("2018_11_1_plate.csv")
plate.2 = read.csv("2018_11_2_plate.csv")

plate.2 = plate.2[,-1]
plate  = cbind(plate.1, plate.2)
finalcycle3 = rbind(plate[2,],plate[1,],plate[3,],plate[6:45,])
sampleid = finalcycle3[1:3,]
finalcycle3.edit = t(finalcycle3)
finalcycle3.edit = as.data.frame(finalcycle3.edit)
colnames(finalcycle3.edit) = c('unique_sampleID_number','reaction_type',"Starting Quantity",1:40)
fc.6.total = finalcycle3.edit[2:495,]
fc.6.total = t(fc.6.total)
fc.6.total = cbind.data.frame(11:40, fc.6.total[14:43,])
names(fc.6.total) = c("Cycles",2:495)
fc.6.total %<>% mutate_if(is.character,as.numeric)
cp40 = pcrfit(1,27,data = fc.6.total, model = l5)
plot(cp40,type="all")
efficiency(cp40)
m.6.efficiencies =  pcrbatch(fc.6.total, model = l5)
cps = m.6.efficiencies[7,]
sampleid = as.matrix(sampleid)
cps = as.matrix(cps)
cps.sa = rbind(sampleid, cps)
###############################################################################################
########################################################## 
################### Initial Data Framing #################
########################################################## 

#deriv = deriv_complete
deriv = cps.sa
# Remove extra labels column
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=FALSE)
# Remove blank column (4th)
#deriv = deriv[,-5]
# Rename columns
colnames(deriv)=c( "sampleID", "reaction_type", "starting_quantity", "cpD1")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]

#deriv = deriv[,-c(6,7)]
deriv = deriv[,-5]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
#deriv = deriv[,-6]
deriv = deriv[,-5]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))

### COMPLETED INITIAL DATA FRAMING ###

##########################################################
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
#cycle1 = read.csv(file = "2018_8/2018_8_1_plate.csv", header = FALSE)
cycle1 = plate
# Create complete set of reaction data (derivative and cycle)
reaction = rbind(cps.sa, as.matrix(cycle1))
reaction = as.data.frame(reaction)
# Remove repeat labeling
#replace = reaction[5:8,]
reaction = reaction[-c(5:8),]
#reaction = plyr::rbind.fill(replace, reaction)
# Transpose so column headers at top
reaction = as.data.frame(t(reaction))
reaction = reaction[,-5]
# Replace column names with first row
colnames(reaction) <- as.character(unlist(reaction[1,]))
reaction = reaction[-1,]
colnames(reaction)[4] = "cpD1"
reaction$cpD1 = as.numeric(as.character(reaction$cpD1))
# Filter unusual observations (CT value less than 10)
unusual_obs_2018_6 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_6, file="Unusual_Obs_2018_6.csv")

# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
#library("rowr")
# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data = calib_data[order(calib_data$starting_quantity),]
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,4]
allP = filter(calib_data, reaction_type=="all_products")[,3:4]

#Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1" )
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)

#Apply log scale to test1 and allP CT values
calib_data$allPln = log(calib_data$allP)
calib_data$test1ln = log(calib_data$test1)

write.csv(calib_data, file = "calib_2018_11#2.csv")


### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))

# Apply log scale to test1 and allP CT values
exp_data$test1.exp.ln = log(exp_data$test1.exp)
exp_data$allP.exp.ln = log(exp_data$allP.exp)

write.csv(exp_data, file = "exp_2018_11#2.csv")

####################################AUGUST#####################################################
plate.3 = read.csv("2018_8_1_plate.csv")
plate.4 = read.csv("2018_8_2_plate.csv")
plate.5 = read.csv("2018_8_3_plate.csv")

plate.4 = plate.4[,-1]
plate.5 = plate.5[,-1]
plate  = cbind(plate.3, plate.4, plate.5)
finalcycle3 = rbind(plate[2,],plate[1,],plate[3,],plate[6:45,])
sampleid = finalcycle3[1:3,]
finalcycle3.edit = t(finalcycle3)
finalcycle3.edit = as.data.frame(finalcycle3.edit)
colnames(finalcycle3.edit) = c('unique_sampleID_number','reaction_type',"Starting Quantity",1:40)
fc.6.total = finalcycle3.edit[2:609,]
fc.6.total = t(fc.6.total)
fc.6.total = cbind.data.frame(11:40, fc.6.total[14:43,])
names(fc.6.total) = c("Cycles",2:609)
fc.6.total %<>% mutate_if(is.character,as.numeric)
cp40 = pcrfit(1,27,data = fc.6.total, model = l5)
plot(cp40,type="all")
efficiency(cp40)
m.6.efficiencies =  pcrbatch(fc.6.total, model = l5)
cps = m.6.efficiencies[7,]
sampleid = as.matrix(sampleid)
cps = as.matrix(cps)
cps.sa = rbind(sampleid, cps)
#############################################################################################
########################################################## 
################### Initial Data Framing #################
########################################################## 

#deriv = deriv_complete
deriv = cps.sa
# Remove extra labels column
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=FALSE)
# Remove blank column (4th)
#deriv = deriv[,-5]
# Rename columns
colnames(deriv)=c( "sampleID", "reaction_type", "starting_quantity", "cpD1")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]

#deriv = deriv[,-c(6,7)]
deriv = deriv[,-5]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
#deriv = deriv[,-6]
deriv = deriv[,-5]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))

### COMPLETED INITIAL DATA FRAMING ###

##########################################################
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
#cycle1 = read.csv(file = "2018_8/2018_8_1_plate.csv", header = FALSE)
cycle1 = plate
# Create complete set of reaction data (derivative and cycle)
reaction = rbind(cps.sa, as.matrix(cycle1))
reaction = as.data.frame(reaction)
# Remove repeat labeling
#replace = reaction[5:8,]
reaction = reaction[-c(5:8),]
#reaction = plyr::rbind.fill(replace, reaction)
# Transpose so column headers at top
reaction = as.data.frame(t(reaction))
reaction = reaction[,-5]
# Replace column names with first row
colnames(reaction) <- as.character(unlist(reaction[1,]))
reaction = reaction[-1,]
colnames(reaction)[4] = "cpD1"
reaction$cpD1 = as.numeric(as.character(reaction$cpD1))
# Filter unusual observations (CT value less than 10)
unusual_obs_2018_6 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_6, file="Unusual_Obs_2018_6.csv")

# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
#library("rowr")
# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data = calib_data[order(calib_data$starting_quantity),]
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,4]
allP = filter(calib_data, reaction_type=="all_products")[,3:4]

#Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1" )
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)

#Apply log scale to test1 and allP CT values
calib_data$allPln = log(calib_data$allP)
calib_data$test1ln = log(calib_data$test1)

write.csv(calib_data, file = "calib_2018_8#2.csv")


### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))

# Apply log scale to test1 and allP CT values
exp_data$test1.exp.ln = log(exp_data$test1.exp)
exp_data$allP.exp.ln = log(exp_data$allP.exp)

write.csv(exp_data, file = "exp_2018_8#2.csv")

####################################JUNE#####################################################
plate.6 = read.csv("2018_6_1_plate.csv")

plate  = plate.6
finalcycle3 = rbind(plate[2,],plate[1,],plate[3,],plate[6:45,])
sampleid = finalcycle3[1:3,]
finalcycle3.edit = t(finalcycle3)
finalcycle3.edit = as.data.frame(finalcycle3.edit)
colnames(finalcycle3.edit) = c('unique_sampleID_number','reaction_type',"Starting Quantity",1:40)
fc.6.total = finalcycle3.edit[2:277,]
fc.6.total = t(fc.6.total)
fc.6.total = cbind.data.frame(11:40, fc.6.total[14:43,])
names(fc.6.total) = c("Cycles",2:277)
fc.6.total %<>% mutate_if(is.character,as.numeric)
cp40 = pcrfit(1,27,data = fc.6.total, model = l5)
plot(cp40,type="all")
efficiency(cp40)
m.6.efficiencies =  pcrbatch(fc.6.total, model = l5)
cps = m.6.efficiencies[7,]
sampleid = as.matrix(sampleid)
cps = as.matrix(cps)
cps.sa = rbind(sampleid, cps)
##############################################################################################
########################################################## 
################### Initial Data Framing #################
########################################################## 

#deriv = deriv_complete
deriv = cps.sa
# Remove extra labels column
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=FALSE)
# Remove blank column (4th)
#deriv = deriv[,-5]
# Rename columns
colnames(deriv)=c( "sampleID", "reaction_type", "starting_quantity", "cpD1")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]

#deriv = deriv[,-c(6,7)]
deriv = deriv[,-5]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
#deriv = deriv[,-6]
deriv = deriv[,-5]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))

### COMPLETED INITIAL DATA FRAMING ###

##########################################################
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
#cycle1 = read.csv(file = "2018_8/2018_8_1_plate.csv", header = FALSE)
cycle1 = plate
# Create complete set of reaction data (derivative and cycle)
reaction = rbind(cps.sa, as.matrix(cycle1))
reaction = as.data.frame(reaction)
# Remove repeat labeling
#replace = reaction[5:8,]
reaction = reaction[-c(5:8),]
#reaction = plyr::rbind.fill(replace, reaction)
# Transpose so column headers at top
reaction = as.data.frame(t(reaction))
reaction = reaction[,-5]
# Replace column names with first row
colnames(reaction) <- as.character(unlist(reaction[1,]))
reaction = reaction[-1,]
colnames(reaction)[4] = "cpD1"
reaction$cpD1 = as.numeric(as.character(reaction$cpD1))
# Filter unusual observations (CT value less than 10)
unusual_obs_2018_6 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_6, file="Unusual_Obs_2018_6.csv")

# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
#library("rowr")
# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data = calib_data[order(calib_data$starting_quantity),]
calib_data$starting_quantity = as.numeric(calib_data$starting_quantity)
calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,4]
allP = filter(calib_data, reaction_type=="all_products")[,3:4]

#Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1" )
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)

#Apply log scale to test1 and allP CT values
calib_data$allPln = log(calib_data$allP)
calib_data$test1ln = log(calib_data$test1)

write.csv(calib_data, file = "calib_2018_6#2.csv")


### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))

# Apply log scale to test1 and allP CT values
exp_data$test1.exp.ln = log(exp_data$test1.exp)
exp_data$allP.exp.ln = log(exp_data$allP.exp)

write.csv(exp_data, file = "exp_2018_6#2.csv")

###########################################################################################################
###################################### Heirarchal analysis ################################################
###########################################################################################################

library(ggplot2)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
library(MASS)
library(stringr)
## Ordinal Net package ##
library("ordinalNet")
library("rgl")
library("qpcR")
library(dplyr)
library(magrittr)
library(gridExtra)
library(qtl)
library(vqtl)



# set directory to the Heirarchical folder


#### reading in and setting up calibrated data ####
# MONTH 1 (2018_6 / JUNE) CALIBRATED DATA FRAME 
calib_data_6 = (read.csv("calib_2018_6#2.csv")[,c(2:4)])
calib_data_6$ztest1 = (calib_data_6$test1 - mean(calib_data_6$test1))/sd(calib_data_6$test1)
calib_data_6$zallP = (calib_data_6$allP - mean(calib_data_6$allP))/sd(calib_data_6$allP)
calib_data_6$month ='june'
#calib_data_6

# MONTH 2 (2018_8 / AUGUST) CALIBRATED DATA FRAME
calib_data_8 = (read.csv("calib_2018_8#2.csv")[,c(2:4)])
calib_data_8$ztest1 = (calib_data_8$test1 - mean(calib_data_8$test1))/sd(calib_data_8$test1)
calib_data_8$zallP = (calib_data_8$allP - mean(calib_data_8$allP))/sd(calib_data_8$allP)
calib_data_8$month ='aug'
#calib_data_8

# MONTH 3 (2018_11 / NOVEMBER) CALIBRATED DATA FRAME
calib_data_11 = (read.csv("calib_2018_11#2.csv")[,c(2:4)])
calib_data_11$ztest1 = (calib_data_11$test1 - mean(calib_data_11$test1))/sd(calib_data_11$test1)
calib_data_11$zallP = (calib_data_11$allP - mean(calib_data_11$allP))/sd(calib_data_11$allP)
calib_data_11$month ='nov'
#calib_data_11

# Combined Calib d.f. for all months
calib_data = rbind(calib_data_6, calib_data_8, calib_data_11)

# Create dummy varible columns for each month
calib_data$june = ifelse(str_detect(calib_data[,6], "june"), 1, 0)
calib_data$aug = ifelse(str_detect(calib_data[,6], "aug"), 1, 0)
#calib_data = calib_data[,-6]

# create a subset only containing the startq, zvalues, and dummy months
calib_subset = calib_data[,c(1, 4,5,8,7)]
######

#### Plotting calibrated data ####
#graphing log starting quantity to cp values
plot(calib_data$allP,log(calib_data$startq), col = 'red', main = "Hierarchical Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(10, 40), ylim = c(-6, 6))
points(calib_data$test1,log(calib_data$startq), col = 'blue')
abline(lm(log(calib_data$startq)~calib_data$allP), col = 'red')
abline(lm(log(calib_data$startq)~calib_data$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)


#graphing log starting quantity to cp values by month
#June
plot(calib_data_6$allP,log(calib_data_6$startq), col = 'red', main = "June Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(10, 40), ylim = c(-6, 6))
points(calib_data_6$test1,log(calib_data_6$startq), col = 'blue')
abline(lm(log(calib_data_6$startq)~calib_data_6$allP), col = 'red')
abline(lm(log(calib_data_6$startq)~calib_data_6$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#August
plot(calib_data_8$allP,log(calib_data_8$startq), col = 'red', main = "August Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(10, 40), ylim = c(-6, 6))
points(calib_data_8$test1,log(calib_data_8$startq), col = 'blue')
abline(lm(log(calib_data_8$startq)~calib_data_8$allP), col = 'red')
abline(lm(log(calib_data_8$startq)~calib_data_8$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#November
plot(calib_data_11$allP,log(calib_data_11$startq), col = 'red', main = "November Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(10, 40), ylim = c(-6, 6))
points(calib_data_11$test1,log(calib_data_11$startq), col = 'blue')
abline(lm(log(calib_data_11$startq)~calib_data_11$allP), col = 'red')
abline(lm(log(calib_data_11$startq)~calib_data_11$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)
#####

hist(calib_data$allP, col = "red")
hist(calib_data$test1, col = "blue")

# ###### POLR models ######
# # Ordinal Logistic Regression Model 
# model = polr(as.factor(calib_subset$startq) ~ ., data=calib_subset, Hess = TRUE)
# #(summary(model))
# (ctable <- coef(summary(model)))
# ## calculate and store p values
# p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
# options(scipen=999)
# ## combined table
# (ctable <- cbind(ctable, "p value" = p))
# 
# ###### Ordinal Net models #######

#define ordinal model starq~zallP+ztest1+month
ordmod = ordinalNet(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq))
summary(ordmod)
coef(ordmod, matrix=TRUE)
#kfold cv
set.seed(3)
ordfit = ordinalNetTune(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)
head(ordfit$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit$loglik))
head(coef(ordfit$fit, matrix = TRUE, whichLambda = bestLambdaIndex))


###### Experimental Data #######

# MONTH 1 (2018_8 / AUGUST) EXPERIMENTAL DATA FRAME

exp_data_8 = na.omit(read.csv("exp_2018_8#2.csv")[,-c(1,5,6,7)]) 
exp_data_8$ztest1 = (exp_data_8$test1 - mean(exp_data_8$test1))/sd(exp_data_8$test1)
exp_data_8$zallP = (exp_data_8$allP - mean(exp_data_8$allP))/sd(exp_data_8$allP)
exp_data_8$month ='aug'
exp_data_8$sampleID.exp = as.factor(exp_data_8$sampleID.exp)
#exp_data_8

# MONTH 2 (2018_6 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_6 = na.omit(read.csv("exp_2018_6#2.csv")[,-c(1,5,6,7)])
exp_data_6$ztest1 = (exp_data_6$test1 - mean(exp_data_6$test1))/sd(exp_data_6$test1)
exp_data_6$zallP = (exp_data_6$allP - mean(exp_data_6$allP))/sd(exp_data_6$allP)
exp_data_6$month ='june'
#exp_data_6

# MONTH 3 (2018_11 / NOVEMBER) EXPERIMENTAL DATA FRAME
exp_data_11 = na.omit(read.csv("exp_2018_11#2.csv")[,-c(1,5,6,7)])
exp_data_11$ztest1 = (exp_data_11$test1 - mean(exp_data_11$test1))/sd(exp_data_11$test1)
exp_data_11$zallP = (exp_data_11$allP - mean(exp_data_11$allP))/sd(exp_data_11$allP)
exp_data_11$month ='nov'
exp_data_11$sampleID.exp = as.factor(exp_data_11$sampleID.exp)
#exp_data_11

# Combined exp d.f. for all months
exp_data = rbind(exp_data_6, exp_data_8, exp_data_11)

# Create dummy varible columns for each month
exp_data$june = ifelse(str_detect(exp_data[,6], "june"), 1, 0)
exp_data$aug = ifelse(str_detect(exp_data[,6], "aug"), 1, 0)


# Drop rows containing NA
exp_subset = exp_data[,c(1, 4, 5, 7, 8)]

#### calculating experimental starting quantity ####
probmat = predict(ordfit$fit, as.matrix(exp_subset[,2:5]))
probmat[1:10,]


##### finding adjustment value and adjusted test 1 in calibrated #####

group = split.data.frame(calib_data, calib_data$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(log(k$allP), log(k$test1)))#this needs to be fixed***********
}


##### Creating a dataframe with the stq and adjustment values #########
calib_adj = as.data.frame(cbind(unique(calib_data$startq), unique(adjval)))
colnames(calib_adj) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
calib_adj$adj = (calib_adj$adj)*1000


# convert test1 and allp to femtograms
exp_data[,c(2,3)] = exp_data[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data$exp.adjust = probmat%*%calib_adj$adj

# Create new column with stress product (VQTL input)
exp_data$exp.adjustTest1 = exp_data$test1.exp+exp_data$exp.adjust

# convert allP and adjusted test 1 to femtograms
exp_data$stress = exp_data$allP.exp - exp_data$exp.adjustTest1

boxplot(exp_data$allP.exp, exp_data$test1.exp, exp_data$stress, main = "Boxplot of Experimental All Products, Test 1, and Stress",
        names =c("All Products", "Test 1", "Stress"), ylab = "Cp value", col =c("blue", "red", "green"))

hist(exp_data$stress, col = "light blue")
### analyzing negative stress ###
negstress = exp_data %>% filter(exp_data$stress<0)
hist(negstress$stress, col = "light green")
length(negstress$stress)/length(exp_data$stress)

# making negstress "NA" #
exp_data= exp_data[which(exp_data[,11] >0),]
qplot()
### Write the exp data with the stress product as a new data frame for vqtl matching
write.csv(exp_data, "Hierarchical_exp_data_stress#3.csv")

##################################################################################################################


barcode.breed = read.csv("vqtlinput.csv")
barcode.breed = cbind.data.frame(barcode.breed$Barcode, barcode.breed$BreedType)
names(barcode.breed) = c("Barcode", "BreedType")
barcode.breed$Barcode = substring(barcode.breed$Barcode, regexpr("_", barcode.breed$Barcode) + 1)
barcode.breed$Barcode = substring(barcode.breed$Barcode, regexpr("_", barcode.breed$Barcode) + 1)
barcode.breed = barcode.breed[3:367,]
barcode.breed = na.omit(barcode.breed)

exp_data = read.csv("Hierarchical_exp_data_stress#3.csv")
exp_data$sampleID.exp = substring(exp_data$sampleID.exp, regexpr("_", exp_data$sampleID.exp) + 1)
exp_data$sampleID.exp = substring(exp_data$sampleID.exp, regexpr("_", exp_data$sampleID.exp) + 1)
exp_data = na.omit(exp_data)
exp_data = exp_data[,-1]
names(exp_data)[1]= 'Barcode'
exp_data = exp_data[-318,]#somehow we have duplicate sample IDs

left  <- exp_data %>% left_join(barcode.breed,  by = 'Barcode')
left.str = left[c(1,11)]
setwd("/Users/michaelcopeland/Downloads")

fulldata = read.csv("SamplingPlan_dat2.csv")#OSF
fulldata$Barcode = substring(fulldata$Barcode, regexpr("_", fulldata$Barcode) + 1)
fulldata$Barcode = substring(fulldata$Barcode, regexpr("_", fulldata$Barcode) + 1)
fulldata_sub = fulldata[,1:15]
left.stress <- fulldata_sub %>% left_join(left.str,  by = 'Barcode')
fulldata.all = cbind.data.frame(left.stress, fulldata[,18:3252])
fulldata.na = na.omit(fulldata.all)
fulldata.vqtl = cbind.data.frame(fulldata.na$stress,fulldata.na$Barcode,fulldata.na$Date,fulldata.na$BreedType,
                                 fulldata.na$Genotype, fulldata.na[17:3251])



setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/qPCR2vQTL")
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = FALSE)
names(fulldata.vqtl) = hybrid.inbred[1,]
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = TRUE)
fulldata.vqtl = rbind.data.frame(as.matrix(hybrid.inbred[1:2,]),as.matrix(fulldata.vqtl))
fulldata.vqtl$stress = as.character(fulldata.vqtl$stress)
rownames(fulldata.vqtl) <- NULL;
fulldata.vqtl[1,1] = ""
fulldata.vqtl[2,1] = ""
fulldata.vqtl = fulldata.vqtl[,-c(2,3,5)]
write.csv(fulldata.vqtl, "fullvqtldata#3.csv", row.names = FALSE)
test_full <- read.cross(file = "fullvqtldata#3.csv", format = "csv")

test_full <- drop.nullmarkers(test_full)
test_full <- calc.genoprob(test_full)

test_full <- calc.genoprob(test_full, error.prob = .001)
hy_p1 <- scanone(cross = test_full, pheno.col = 'stress')
hyv_p2 <- scanonevar(cross = test_full, 
                     mean.formula = stress ~ BreedType*mean.QTL.add, 
                     var.formula = ~ BreedType*var.QTL.add, 
                     return.covar.effects = TRUE)



#ratio if femtograms and substraction if log scale
################################################## by month analysis #############################################
################################################## by month analysis #############################################
################################################## by month analysis #############################################
calib_data_6 = (read.csv("calib_2018_6#2.csv")[,c(2:4)])
calib_data_6$ztest1 = (calib_data_6$test1 - mean(calib_data_6$test1))/sd(calib_data_6$test1)
calib_data_6$zallP = (calib_data_6$allP - mean(calib_data_6$allP))/sd(calib_data_6$allP)
calib_data_6$month ='june'
#calib_data_6

calib_subset_6 = calib_data_6[,c(1, 4,5)]
ordfit = ordinalNetTune(as.matrix(calib_subset_6[,2:3]), as.factor(calib_subset_6$startq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)

# MONTH 2 (2018_6 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_6 = na.omit(read.csv("exp_2018_6#2.csv")[,-c(1,5,6,7)])
exp_data_6$ztest1 = (exp_data_6$test1 - mean(exp_data_6$test1))/sd(exp_data_6$test1)
exp_data_6$zallP = (exp_data_6$allP - mean(exp_data_6$allP))/sd(exp_data_6$allP)
exp_data_6$month ='june'
#exp_data_6


# Drop rows containing NA
exp_subset_6 = exp_data_6[,c(1, 4, 5)]

#### calculating experimental starting quantity ####
probmat = predict(ordfit$fit, as.matrix(exp_subset_6[,2:3]))
probmat[1:10,]
##### finding adjustment value and adjusted test 1 in calibrated #####

group = split.data.frame(calib_data_6, calib_data_6$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(log(k$allP), log(k$test1)))#this needs to be fixed***********
}


##### Creating a dataframe with the stq and adjustment values #########
calib_adj = as.data.frame(cbind(unique(calib_data_6$startq), unique(adjval)))
colnames(calib_adj) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
calib_adj$adj = (calib_adj$adj)*1000


# convert test1 and allp to femtograms
exp_data_6[,c(2,3)] = exp_data_6[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data_6$exp.adjust = probmat%*%calib_adj$adj

# Create new column with stress product (VQTL input)
exp_data_6$exp.adjustTest1 = exp_data_6$test1.exp+exp_data_6$exp.adjust #?????

# convert allP and adjusted test 1 to femtograms
exp_data_6$stress = exp_data_6$allP.exp - exp_data_6$exp.adjustTest1    #?????

boxplot(exp_data_6$allP.exp, exp_data_6$test1.exp, exp_data_6$stress, main = "Boxplot of Experimental All Products, Test 1, and Stress",
        names =c("All Products", "Test 1", "Stress"), ylab = "Cp value", col =c("blue", "red", "green"))

hist(exp_data_6$stress, col = "light blue")
### analyzing negative stress ###
negstress = exp_data_6 %>% filter(exp_data_6$stress<0)
hist(negstress$stress, col = "light green")
length(negstress$stress)/length(exp_data_6$stress)

# making negstress "NA" #
exp_data_6= exp_data_6[which(exp_data_6[,9] >0),]

##############################################################################################################


calib_data_8 = (read.csv("calib_2018_8#2.csv")[,c(2:4)])
calib_data_8$ztest1 = (calib_data_8$test1 - mean(calib_data_8$test1))/sd(calib_data_8$test1)
calib_data_8$zallP = (calib_data_8$allP - mean(calib_data_8$allP))/sd(calib_data_8$allP)
calib_data_8$month ='August'
#calib_data_8

calib_subset_8 = calib_data_8[,c(1, 4,5)]
ordfit = ordinalNetTune(as.matrix(calib_subset_8[,2:3]), as.factor(calib_subset_8$startq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)

# MONTH 2 (2018_8 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_8 = na.omit(read.csv("exp_2018_8#2.csv")[,-c(1,5,6,7)])
exp_data_8$ztest1 = (exp_data_8$test1 - mean(exp_data_8$test1))/sd(exp_data_8$test1)
exp_data_8$zallP = (exp_data_8$allP - mean(exp_data_8$allP))/sd(exp_data_8$allP)
exp_data_8$month ='August'
#exp_data_6


# Drop rows containing NA
exp_subset_8 = exp_data_8[,c(1, 4, 5)]

#### calculating experimental starting quantity ####
probmat = predict(ordfit$fit, as.matrix(exp_subset_8[,2:3]))
probmat[1:10,]
##### finding adjustment value and adjusted test 1 in calibrated #####

group = split.data.frame(calib_data_8, calib_data_8$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(log(k$allP), log(k$test1)))#this needs to be fixed***********
}


##### Creating a dataframe with the stq and adjustment values #########
calib_adj = as.data.frame(cbind(unique(calib_data_8$startq), unique(adjval)))
colnames(calib_adj) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
calib_adj$adj = (calib_adj$adj)*1000


# convert test1 and allp to femtograms
exp_data_8[,c(2,3)] = exp_data_8[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data_8$exp.adjust = probmat%*%calib_adj$adj

# Create new column with stress product (VQTL input)
exp_data_8$exp.adjustTest1 = exp_data_8$test1.exp+exp_data_8$exp.adjust #?????

# convert allP and adjusted test 1 to femtograms
exp_data_8$stress = exp_data_8$allP.exp - exp_data_8$exp.adjustTest1    #?????

boxplot(exp_data_8$allP.exp, exp_data_8$test1.exp, exp_data_8$stress, main = "Boxplot of Experimental All Products, Test 1, and Stress",
        names =c("All Products", "Test 1", "Stress"), ylab = "Cp value", col =c("blue", "red", "green"))

hist(exp_data_8$stress, col = "light blue")
### analyzing negative stress ###
negstress = exp_data_8 %>% filter(exp_data_8$stress<0)
hist(negstress$stress, col = "light green")
length(negstress$stress)/length(exp_data_8$stress)

# making negstress "NA" #
exp_data_8= exp_data_11[which(exp_data_8[,9] >0),]

###############################################################################################################
calib_data_11 = (read.csv("calib_2018_11#2.csv")[,c(2:4)])
calib_data_11$ztest1 = (calib_data_11$test1 - mean(calib_data_11$test1))/sd(calib_data_11$test1)
calib_data_11$zallP = (calib_data_11$allP - mean(calib_data_11$allP))/sd(calib_data_11$allP)
calib_data_11$month ='Novemeber'
#calib_data_11

calib_subset_11 = calib_data_11[,c(1, 4,5)]
ordfit = ordinalNetTune(as.matrix(calib_subset_11[,2:3]), as.factor(calib_subset_11$startq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)

# MONTH 2 (2018_11 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_11 = na.omit(read.csv("exp_2018_11#2.csv")[,-c(1,5,6,7)])
exp_data_11$ztest1 = (exp_data_11$test1 - mean(exp_data_11$test1))/sd(exp_data_11$test1)
exp_data_11$zallP = (exp_data_11$allP - mean(exp_data_11$allP))/sd(exp_data_11$allP)
exp_data_11$month ='Novemeber'
#exp_data_11


# Drop rows containing NA
exp_subset_11 = exp_data_11[,c(1, 4, 5)]

#### calculating experimental starting quantity ####
probmat = predict(ordfit$fit, as.matrix(exp_subset_11[,2:3]))
probmat[1:10,]
##### finding adjustment value and adjusted test 1 in calibrated #####

group = split.data.frame(calib_data_11, calib_data_11$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(log(k$allP), log(k$test1)))#this needs to be fixed***********
}


##### Creating a dataframe with the stq and adjustment values #########
calib_adj = as.data.frame(cbind(unique(calib_data_11$startq), unique(adjval)))
colnames(calib_adj) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
calib_adj$adj = (calib_adj$adj)*1000


# convert test1 and allp to femtograms
exp_data_11[,c(2,3)] = exp_data_11[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data_11$exp.adjust = probmat%*%calib_adj$adj

# Create new column with stress product (VQTL input)
exp_data_11$exp.adjustTest1 = exp_data_11$test1.exp+exp_data_11$exp.adjust #?????

# convert allP and adjusted test 1 to femtograms
exp_data_11$stress = exp_data_11$allP.exp - exp_data_11$exp.adjustTest1    #?????

boxplot(exp_data_11$allP.exp, exp_data_11$test1.exp, exp_data_11$stress, main = "Boxplot of Experimental All Products, Test 1, and Stress",
        names =c("All Products", "Test 1", "Stress"), ylab = "Cp value", col =c("blue", "red", "green"))

hist(exp_data_11$stress, col = "light blue")
### analyzing negative stress ###
negstress = exp_data_11 %>% filter(exp_data_11$stress<0)
hist(negstress$stress, col = "light green")
length(negstress$stress)/length(exp_data_11$stress)

# we want to replace the negative values with N/A
exp_data_11= exp_data_11[which(exp_data_11[,9] >0),]


exp_data = rbind.data.frame(exp_data_6,exp_data_8,exp_data_11)
### Write the exp data with the stress product as a new data frame for vqtl matching
write.csv(exp_data, "Hierarchical_exp_data_stress_by_month.csv")

################################################################################################################

barcode.breed = read.csv("vqtlinput.csv")
barcode.breed = cbind.data.frame(barcode.breed$Barcode, barcode.breed$BreedType)
names(barcode.breed) = c("Barcode", "BreedType")
barcode.breed$Barcode = substring(barcode.breed$Barcode, regexpr("_", barcode.breed$Barcode) + 1)
barcode.breed$Barcode = substring(barcode.breed$Barcode, regexpr("_", barcode.breed$Barcode) + 1)
barcode.breed = barcode.breed[3:367,]
barcode.breed = na.omit(barcode.breed)

exp_data = read.csv("Hierarchical_exp_data_stress_by_month.csv")
exp_data$sampleID.exp = substring(exp_data$sampleID.exp, regexpr("_", exp_data$sampleID.exp) + 1)
exp_data$sampleID.exp = substring(exp_data$sampleID.exp, regexpr("_", exp_data$sampleID.exp) + 1)
exp_data = na.omit(exp_data)
exp_data = exp_data[,-1]
names(exp_data)[1]= 'Barcode'
exp_data = exp_data[-318,]#somehow we have duplicate sample IDs

left  <- exp_data %>% left_join(barcode.breed,  by = 'Barcode')
left.str = left[c(1,9)]
setwd("/Users/michaelcopeland/Downloads")

fulldata = read.csv("SamplingPlan_dat2.csv")#OSF
fulldata$Barcode = substring(fulldata$Barcode, regexpr("_", fulldata$Barcode) + 1)
fulldata$Barcode = substring(fulldata$Barcode, regexpr("_", fulldata$Barcode) + 1)
fulldata_sub = fulldata[,1:15]
left.stress <- fulldata_sub %>% left_join(left.str,  by = 'Barcode')
fulldata.all = cbind.data.frame(left.stress, fulldata[,18:3252])
fulldata.na = na.omit(fulldata.all)
fulldata.vqtl = cbind.data.frame(fulldata.na$stress,fulldata.na$Barcode,fulldata.na$Date,fulldata.na$BreedType,
                                 fulldata.na$Genotype, fulldata.na[17:3251])



setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/qPCR2vQTL")
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = FALSE)
names(fulldata.vqtl) = hybrid.inbred[1,]
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = TRUE)
fulldata.vqtl = rbind.data.frame(as.matrix(hybrid.inbred[1:2,]),as.matrix(fulldata.vqtl))
fulldata.vqtl$stress = as.character(fulldata.vqtl$stress)
rownames(fulldata.vqtl) <- NULL;
fulldata.vqtl[1,1] = ""
fulldata.vqtl[2,1] = ""
fulldata.vqtl = fulldata.vqtl[,-c(2,3,5)]
write.csv(fulldata.vqtl, "fullvqtldata#3.csv", row.names = FALSE)
test_full <- read.cross(file = "fullvqtldata#3.csv", format = "csv")

test_full <- drop.nullmarkers(test_full)
test_full <- calc.genoprob(test_full)

test_full <- calc.genoprob(test_full, error.prob = .001)
hy_p1 <- scanone(cross = test_full, pheno.col = 'stress')
hyv_p2 <- scanonevar(cross = test_full, 
                     mean.formula = stress ~ BreedType*mean.QTL.add, 
                     var.formula = ~ BreedType*var.QTL.add, 
                     return.covar.effects = TRUE)
