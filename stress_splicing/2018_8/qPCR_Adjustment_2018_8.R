########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)
library(glm.predict)

# Mac Directory
setwd("~/Stapleton/Copeland/stapleton_lab/Stress_Splicing/2018_8")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")
# PC Directory
#setwd(~/Desktop/GIThub/StapletonLab/StressSplicing/qPCR/)

### READ IN DERIVATIVE DATA###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
deriv.1<-read.csv(file = "2018_8_1_qPCR_Output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_8_2_qPCR_Output.csv", header=FALSE)
deriv.3=read.csv(file = "2018_8_3_qPCR_Output.csv", header=FALSE)
deriv_complete=as.data.frame(cbind(deriv.1, deriv.2, deriv.3))

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "2018_8_1_plate_qPCR_output.csv", header=FALSE)


########################################################## 
################### Initial Data Framing #################
########################################################## 

deriv = deriv_complete
# Remove extra column 
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Rename columns
colnames(deriv)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(6,7)]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
# IF "minus" RETURNS EMPTY VALUES, COMMENT OUT COMMAND BELOW
deriv = deriv[-minus,]
deriv = deriv[,-6]
# Remove two extra label rows from center of data frame
deriv['label.row'] = grepl('3', deriv$starting_quantity)
extra = which(deriv$label.row)
deriv = deriv[-extra,]
deriv = deriv[,-6]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))
### COMPLETED INITIAL DATA FRAMING ###


########################################################## 
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
unusual_obs_2018_8 = deriv %>% filter(deriv$cpD1 < 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# ### WORK ON: Appending raw plate cycle vals to unusual obs d.f.
# # Read in raw cycle data
# cycle1 = read.csv(file = "2018_8_1_plate.csv", header = FALSE)
# cycle2 = read.csv(file = "2018_8_2_plate.csv", header = FALSE)
# cycle3 = read.csv(file = "2018_8_3_plate.csv", header = FALSE)
# cycle = as.data.frame(cbind(cycle1, cycle2, cycle3))
# unusual_obs_2018_8 = match()
# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 

# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data = calib_data[order(calib_data$starting_quantity),]
calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))
# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_data$starting_quantity)){
  sq <- calib_data$starting_quantity[i]
  if(i %% 6 == 1){
    startq = c(startq,sq,sq,sq)
  }
  val <- toString(calib_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 = c(test1, calib_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP = c(allP, calib_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib_data = as.data.frame(cbind(startq, test1, allP))
# Format starting quantity values as decimals, not scientific notation
calib_data$startq = as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq = as.factor(calib_data$startq)
# Calculate ratio of allP/test1 --> PAIRWISE RATIOS -- INPUT FOR OLR MODEL
ratio = calib_data$allP/calib_data$test1
# Append ratios to data set
calib_data = cbind(calib_data, ratio)
### COMPLETED CALIBRATED DATA FRAME ###


########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# # Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find invalid observations - Find counts of each unique sampleID; remove ones with count not equal to 2 from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
countsne2$Var1 = as.numeric(as.character(countsne2$Var1)) 
# Remove observations with count not equal to 2 from data set
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
# Calculate ratios for experimental data 
ratio.exp = exp_data$allP.exp/exp_data$test1.exp
# Append ratios to data set
exp_data = cbind(exp_data, ratio.exp)
### COMPLETED EXPERIMENTAL DATA FRAME ###


##########################################################
############### Combination Ratios for qPCR ##############
##########################################################

startquan = as.character(calib_data$startq)
allprod = calib_data$allP
t1 = calib_data$test1
dat = data.frame(cbind(startquan,allprod,t1), stringsAsFactors = FALSE)
dat$allprod = as.numeric(dat$allprod)
dat$t1 = as.numeric(dat$t1)
#Create divide funtion - every element in column 1 divided by every element in column 2
divide <- function(col1, col2){
  ratio = NULL;
  for (i in col1){
    ratio = c(ratio,i/col2)
  }
  return(ratio)
}
#Subset data by starting quantity
group = split.data.frame(dat, dat$startquan)
# Calculate combination ratios at each starting quantity
combratio = NULL;
for (k in group){
  combratio = cbind(combratio, divide(k$allprod, k$t1))
}
# Create data frame with unique ratios at each starting quantity
startqvalues = rep(unique(startquan), rep(length(unique(startquan)),length(unique(startquan)))) 
newratios.calib = data.frame(rbind(unique(startqvalues), combratio), stringsAsFactors = FALSE)
newratios.calib = t(newratios.calib)
newratios.calib = as.data.frame(newratios.calib)
newratios.calib$combratio = as.numeric(newratios.calib$combratio)
newratios.calib$startqvalues = as.numeric(newratios.calib$startqvalues)
# Duplicate newratios.calib data frame, transpose for boxplot visualizations at each s.q.
newratios.calib.boxplot = as.data.frame(t(newratios.calib))
colnames(newratios.calib.boxplot) = c("0.01", "0.05", "0.10", "0.50", "1.00", "50.00")
newratios.calib.boxplot = newratios.calib.boxplot[-1,]
newratiosvector = as.vector(as.matrix.data.frame(newratios.calib.boxplot))
startqvector = sort(rep(unique(startquan), length(newratios.calib.boxplot$`0.01`)))
newratios.calib.boxplot = as.data.frame(cbind(newratiosvector, startqvector))
### COMPLETED COMBINATION RATIOS ###

##########################################################
########## PROBABILITY MODEL - Calibrated Data ###########
##########################################################

# Calculate z-score for calibrated data
zscore = (calib_data$ratio - mean(calib_data$ratio))/sd(calib_data$ratio)
# Predict calibrated data ratios using experimental data
pred.ratio = zscore*sd(ratio.exp)+mean(ratio.exp)
# Append y (predicted calibrated ratios) to calibrated data frame -- CALIBRATED RATIOS IN TERMS OF EXPERIMENTAL PARAMETERS
calib_data = cbind(calib_data, pred.ratio) 
# Create empty vectors for for-loop input
calib_data$test1 = as.numeric(as.character(calib_data$test1))
calib_data$allP = as.numeric(as.character(calib_data$allP))
adj_val = c()
allP = c()
startq = c()
ratio =calib_data$allP/calib_data$test1
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(calib_data)/3)){
  t_x <- c(calib_data$allP[3*i - 2], calib_data$allP[3*i - 1], calib_data$allP[3*i])
  t_y <- c(calib_data$test1[3*i - 2], calib_data$test1[3*i - 1], calib_data$test1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}
adjusted_test1 <- test1 + adj_val
# Append adjusted test1 values and adjustment value to data set
calib_data=cbind(calib_data,adjusted_test1,adj_val)
# Write Calibrated Data CSV --> Used in "qPCR_Plotting" code for visuals
#write.csv(file="YEAR_MONTH_Calibrated_DF", calib_data)

# Adjustment: allP - test1 -- Using in model to multiply probability matrix by
calib_data$diff = calib_data$allP - calib_data$adjusted_test1

# CREATE DATA FRAME WITH ONLY S.Q. AND ADJUSTMENT VAL
calib_adj = calib_data[,c(1,6)]

average <- function(col1){
  avg = NULL;
  for (i in col1){
    avg = c(avg,mean(col1))
  }
  return(avg)
}
#Subset data by starting quantity
group = split.data.frame(calib_adj, calib_adj$startq)

adj.test1.avg = NULL;
for (k in group){
  adj.test1.avg = c(adj.test1.avg, average(k$adjusted_test1))
}
print(adj.test1.avg)

calib_adj = as.data.frame(unique(cbind(as.character(calib_data$startq), adj.test1.avg)))
calib_adj$adj.test1.avg = as.numeric(as.character(calib_adj$adj.test1.avg))
# Rename columns
colnames(calib_adj)=c("startq", "adj.test1.avg")

# Ordinal Logistic Regression Model - starting quantity as response to calibrated z-score
model = polr(as.factor(calib_data$startq) ~ zscore, Hess = TRUE)
summary(model)

# Calculate experimental data z-score
zscore = (exp_data$ratio.exp - mean(exp_data$ratio.exp))/sd(exp_data$ratio.exp)
prob.matrix = predict(model, zscore, type='p')

# Apply probability matrix to the adjusted test 1 averages
apply(prob.matrix, 1, function(x) x*calib_adj$adj.test1.avg)
exp_data$exp.adjust = colSums(apply(prob.matrix, 1, function(x) x*calib_adj$adj.test1.avg))

# Create new column with stress product (VQTL input)
exp_data$stress = exp_data$allP.exp - exp_data$exp.adjust


###PLOTS###
# Calibrated data - s.q. vs. ratio
plot(as.factor(newratios.calib.boxplot$startqvector), as.character(newratios.calib.boxplot$newratiosvector), xlab='Starting Quantity', ylab='Ratio', 
     main='2018_8 Calibrated Data - Starting Quantities vs. Ratios')
# ##PLOTS##
# # Calibrated data - s.q. vs. ratio
# #Plot for pairwise
# plot(as.factor(calib_data$startq), calib_data$ratio, xlab='Starting Quantity', ylab='Ratio', 
#      main='Calibrated Data - Starting Quantities vs. Pairwise Ratios')
# #Plot for non-pairwise
# plot(as.factor(newratios.calib$startqvalues), newratios.calib$combratio, xlab='Starting Quantity', ylab='Ratio', 
#      main='Calibrated Data - Starting Quantities vs. Non-Pairwise Ratios')
# 
# #Boxlot of Stress Product
# boxplot(exp_data$stress, main='Box Plot of Stress Product', ylab='Stress Product')
# hist(exp_data$stress, xlab='Stress Product', main='Histogram of Stress Product', col='blue')



# #AllP#
# hist(data$allP, xlim=c(0,50), ylim=c(0,100), col=rgb(1,0,0,0.5), main='Histogram of All Products', xlab='All Products Derivative')
# hist(exp_data$allP.exp, xlim=c(0,50), ylim=c(0,100), add=T, col=rgb(0,0,1,0.5))
# legend("topleft",
#         c("Calibration", "Experimental"),
#           fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
# #dev.off()
# #Test1#
# hist(data$test1, xlim=c(0,30), ylim=c(0,80), col=rgb(1,0,0,0.5), main='Histogram of Test 1', xlab='Test 1 Derivative')
# hist(exp_data$test1.exp, xlim=c(0,30), ylim=c(0,80), add=T, col=rgb(0,0,1,0.5))
# legend("topleft",
#        c("Calibration", "Experimental"),
#        fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
# #dev.off()
# #Ratios - Calibrated#
# hist(data$ratio, xlim=c(0,3), ylim=c(0,70), col=rgb(1,0,0,0.5), main='Histogram of Ratios', xlab='Ratio')
# #Ratios - Experimental#
# # Values excluded from histogram that will be further investigated later (their index)
# hist(exp_data$ratio.exp, xlim=c(0,3), ylim=c(0,70), col=rgb(0,0,1,0.5), add=T)
# legend("topleft",
#        c("Calibration", "Experimental"),
#        fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), bty="n")
# #Calib Plot - S.Q. vs. Ratios
# plot(calib_data$startq, calib_data$ratio, xlab='Starting Quantity', ylab='Ratio', 
#       main='Calibrated Data - Starting Quantities vs. Ratios')
# #Calib Plot - Test1 vs. Ratio
# plot(calib_data$test1, calib_data$ratio, xlab='Test 1 Derivative', ylab='Ratio', 
#       main='Calibrated Data - Test 1 Derivative vs. Ratio')


#Compare sample ID's between plate and CT data sets
#Confirm no obs deleted in calculations by checking dimensions
#Delete NTC obs prior to comparison
# plate1 = read.csv(file = "2018_8_1_plate.csv", header=FALSE)
# plate1 = as.data.frame(t(plate1))
# plate1 = plate1[-1,-5]
# ct1 = read.csv(file = "2018_8_1_plate_qPCR_output.csv", header=FALSE)
# ct1 = as.data.frame(t(ct1))
# ct1 = ct1[-c(1,2),]
# # Indicate if sample is NTC (negative control)
# plate1['sampleID_NTC'] = grepl('NTC', plate1$V3)
# ct1['sampleID_NTC'] = grepl('NTC', ct1$V3)
# # Remove NTC samples, indicator (T/F) column, and cpD2 values
# ntc_plate1 = which(plate1$sampleID_NTC)
# ntc_ct1 = which(ct1$sampleID_NTC)
# plate1 = plate1[-ntc_plate1,]
# ct1 = ct1[-ntc_ct1,]
# 
# plate2 = read.csv(file = "2018_8_2_plate.csv", header=FALSE)
# plate2 = as.data.frame(t(plate2))
# plate2 = plate2[-1,-5]
# ct2 = read.csv(file = "2018_08_02_plate_qPCR_output_2019_05_19.csv", header=FALSE)
# ct2 = as.data.frame(t(ct2))
# ct2 = ct2[-c(1,2),]
# # Indicate if sample is NTC (negative control)
# plate2['sampleID_NTC'] = grepl('NTC', plate2$V3)
# ct2['sampleID_NTC'] = grepl('NTC', ct2$V3)
# # Remove NTC samples, indicator (T/F) column, and cpD2 values
# ntc_plate2 = which(plate2$sampleID_NTC)
# ntc_ct2 = which(ct2$sampleID_NTC)
# plate2 = plate2[-ntc_plate2,]
# ct2 = ct2[-ntc_ct2,]
# ###OBS SampleID 174 not included in qPCR output###
# 
# plate3 = read.csv(file = "2018_8_3_plate.csv", header=FALSE)
# plate3 = as.data.frame(t(plate3))
# plate3 = plate3[-1,-5]
# ct3 = read.csv(file = "2018_08_03_plate_qPCR_output_2019_05_19.csv", header=FALSE)
# ct3 = as.data.frame(t(ct3))
# ct3 = ct3[-c(1,2),]
# # Indicate if sample is NTC (negative control)
# plate3['sampleID_NTC'] = grepl('NTC', plate3$V3)
# ct3['sampleID_NTC'] = grepl('NTC', ct3$V3)
# # Remove NTC samples, indicator (T/F) column, and cpD2 values
# ntc_plate3 = which(plate3$sampleID_NTC)
# ntc_ct3 = which(ct3$sampleID_NTC)
# plate3 = plate3[-ntc_plate3,]
# ct3 = ct3[-ntc_ct3,]

# ########################################################## 
# ################ Experimental Data Framing ###############
# ########################################################## 
# # Create/Write data frame for Calibrated values
# exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# # Sort by starting quantity
# exp_df = exp_df[order(exp_df$starting_quantity),]
# # Remove first and last rows (unnecessary labeling)
# exp_df = exp_df[-1,]
# exp_df = exp_df[-nrow(exp_df),]
# exp_df$sampleID = as.numeric(as.character(exp_df$sampleID))
# exp_df$cpD1 = as.numeric(as.character(exp_df$cpD1))
# exp_data = exp_df
# # Order data by sampleID
# exp_data = exp_data[order(exp_data$sampleID),]
# # Find counts of each unique sampleID; for sample with a count not equal to 2
# # Send removed sampleID's to Dr. S, with additional plots of raw cycle values for each of these samples
# # Write CSV file of samples with count not equal to 2 to send to Dr. S for investigation
# ### To work on --> add derivative values in to the CSV file
# ### To work on --> creating a separate CSV file with samples with unusual derivatives
# counts = as.data.frame(table(exp_data$sampleID))
# countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# #write.csv(file="2018_11_SamplesToInvestigate", countsne2)
# # Remove samples without both allP and test1 from data set
# exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]
# 
# # Create empty vectors for for-loop to input cpD1 values
# test1.exp = c()
# allP.exp = c()
# sampleID.exp = c()
# # For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
# for(i in 1:length(exp_data$sampleID)){
#   id.exp = toString(exp_data$sampleID[i])
#   if(i %% 2 == 1){
#     sampleID.exp = c(sampleID.exp, id.exp)
#   }
#   val = toString(exp_data$reaction_type[i])
#   if(strcmp(val, "test1")){
#     test1.exp = c(test1.exp, exp_data$cpD1[i])
#   }
#   if(strcmp(val, "all_products")){
#     allP.exp = c(allP.exp, exp_data$cpD1[i])
#   }
# }
# 
# # Bind test1 and allProd cpD1 values by sample ID
# exp_data = cbind(sampleID.exp, test1.exp, allP.exp)
# exp_data = as.data.frame(exp_data)
# exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
# exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))
# # Create column of ratios
# ratio.exp = exp_data$allP.exp/exp_data$test1.exp
# # Append ratio column to data frame
# exp_data = cbind(exp_data, ratio.exp)
# exp_data$ratio.exp = as.numeric(as.character(exp_data$ratio.exp))
# 
# # Write CSV file
# #write.csv(exp_data, file="2018_11_Experimental_Data_Frame.csv")
# 
# ### COMPLETED EXPERIMENTAL DATA FRAME ###

# ########################################################## 
# ############ ADJUSTMENT MODEL Calibrated Data ############
# ########################################################## 
# 
# # Create empty vectors for for-loop input
# data = as.data.frame(calib_data)
# data$test1 = as.numeric(as.character(data$test1))
# data$allP = as.numeric(as.character(data$allP))
# adj_val = c()
# allP = c()
# startq = c()
# # Calculate ratio values
# ratio =data$allP/data$test1
# # Append ratios to data set
# data=cbind(data,ratio)
# # Itterating through each set of (3) observations performing U-Stats on each set of inputs
# for (i in 1:(nrow(data)/3)){
#   t_x <- c(data$allP[3*i - 2], data$allP[3*i - 1], data$allP[3*i])
#   t_y <- c(data$test1[3*i - 2], data$test1[3*i - 1], data$test1[3*i])
#   adj <- mean(outer(t_x, t_y, "-"))
#   adj_val <- c(adj_val, adj, adj, adj)
# }
# # Calculate adjusted test1 value
# adjusted_test1 <- test1 + adj_val
# # Append adjusted test1 values and adjustment value to data set
# calib_data=cbind(data,adjusted_test1,adj_val)
# # Creating the adjustment model lm(y-axis~x-axis)
# # Changed adj_val^2 to adj_val to try to make the model better --> VERIFY THIS WITH DR. WANG
# adj_model <- lm(adj_val~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
# summary(adj_model) 
# # Plot adjustment model 
# par(mfrow = c(2,2))
# plot(adj_model)
# dev.off()
# # Plot ratios vs. adjustment values
# plot(ratio, adj_val, main="Ratios vs. Adjustment Values", col="blue", xlab="Ratio", ylab="Adjustment Value")
# abline(lm(adj_val~ratio), col="blue")
# 
# ### COMPLETED ADJUSTMENT MODEL - CALIBRATED DATA ###
# 
# ### ADJUSTMENT MODEL - EXPERIMENTAL DATA ### 
# # Using the adjustment model on the expiremental data
# exp.data = data.frame(ratio = exp_data$allP.exp/exp_data$test1.exp)
# 
# predict(adj_model, exp.data, interval = "confidence")
# 
# #---> Fill in any remaining parts of experimental adjusmtent model
# 
# ### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###



