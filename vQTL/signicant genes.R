library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")
Maindata = read.csv("ManchingStressData_Covar.csv")

#top 10 most significant genes
Maindata[[3245]][6674] <- "B"
avgdiff <- tapply(Maindata$Height,Maindata$gpm27, mean)[4:5]
colm = colnames(Maindata)
colm = colm[11:3245]
AandB <- c()
AandBabs <- c()
count = 1
avector <- c()
for (i in colm){
  print(i)
  avector[count] <- as.vector(Maindata[i])
 count = count+1
}


count = 1
for (i in colm){
  print(i)
  if (length(levels(avector[[count2]]))==5){
  temp <- tapply(Maindata$Height,avector[[count2]], mean)[4:5]
  }else if (length(levels(avector[[count2]]))==4){
    temp <- tapply(Maindata$Height,avector[[count2]], mean)[3:4]
  }
  diff = temp[1]-temp[2]
  diffabs = abs(temp[1]-temp[2])
  AandB[count] <- diff
  AandBabs[count] <- diffabs
  count2=count2+1
  count = count+1
}
temp <- tapply(Maindata$Height,avector[[count2-1]], mean)[3:4]
diff = temp[1]-temp[2]
diffabs = abs(temp[1]-temp[2])
AandB[count-1] <- diff
AandBabs[count-1] <- diffabs

t10 = tail(sort(AandB))[1:10]

diffAB = cbind(AandB,colm)
diffABabs = cbind(AandBabs,colm)
top10 <- diffAB[order(diffAB[,"AandB"], decreasing = TRUE)[1:10],]
top10abs <- diffABabs[order(diffABabs[,"AandBabs"], decreasing = TRUE)[1:10],]
top10both <- cbind(top10,top10abs)

#Boxplot code
genet10 = top10[,2]
p = c()
count3 = 1
for (i in top10){
  print(i)
p[count3] <- boxplot(Maindata$Height~droplevels(Maindata$IDP3914, exclude = '-'))
count = count+1
}
Maindata <- Maindata[-c(1, 2), ]
boxplot(Maindata$Height~droplevels(Maindata$mmp46, exclude = '-'))


#Simulated data with every possible combinantion of environment with gene
a = 60
Env = c(1,2,3,4,5,6,7,8)
lw = c(0,1)
ln = c(0,1)
pa = c(0,1)
G1 = c(0,1)
G2 = c(0,1)
G3 = c(0,1)
G4 = c(0,1)
G5 = c(0,1)
G6 = c(0,1)
G7 = c(0,1)
G8 = c(0,1)
G9 = c(0,1)
G10 = c(0,1)
se = rnorm(1000, 0, 0.07)
Height = c()
count4 = 1
G1c = c()
G2c = c()
G3c = c()
G4c = c()
G5c = c()
G6c = c()
G7c = c()
G8c = c()
G9c = c()
G10c = c()
Envc = c()
lwc = c()
lnc = c()
pac = c()
for (i in lw){
  for (o in ln){
    for (u in pa){
     for (w1 in G1){
      for (w2 in G2){
        for (w3 in G3){
          for (w4 in G4){
            for (w5 in G5){
              for (w6 in G6){
                for (w7 in G7){
                  for (w8 in G8){
                    for (w9 in G9){
                      for (w10 in G10){
    #  print(i)
      G1c[count4]<- w1
      G2c[count4]<- w2
      G3c[count4]<- w3
      G4c[count4]<- w4
      G5c[count4]<- w5
      G6c[count4]<- w6
      G7c[count4]<- w7
      G8c[count4]<- w8
      G9c[count4]<- w9
      G10c[count4]<- w10
      
      lwc[count4]<- i
      lnc[count4]<- o
      pac[count4]<- u
      Height[count4] <- a -0.9*i -0.95*o -1.1*u +1.5*w1 + 0.9*w2 +0.7*w3 + 1.3*w4 +1.2*w5 
      + 0.6*w6 +0.8*w7 + 1.1*w8 +1.4*w9 + 0.5*w10 
      + se[sample(1:1000, 1)]
      count4=count4 +1
                      }
                    }
                  }
                }
              }
            }
          }
        }
       }
      }
    }
  }
}
# G1c = as.factor(G1c)
# G2c = as.factor(G2c)
# Envc = as.factor(Envc)
# lwc = as.factor(lwc)
# lnc = as.factor(lnc)
# pac = as.factor(pac)
G1c[G1c == 0] <- "A"
G1c[G1c == 1] <- "B"
G2c[G2c == 0] <- "A"
G2c[G2c == 1] <- "B"
G3c[G3c == 0] <- "A"
G3c[G3c == 1] <- "B"
G4c[G4c == 0] <- "A"
G4c[G4c == 1] <- "B"
G5c[G5c == 0] <- "A"
G5c[G5c == 1] <- "B"
G6c[G6c == 0] <- "A"
G6c[G6c == 1] <- "B"
G7c[G7c == 0] <- "A"
G7c[G7c == 1] <- "B"
G8c[G8c == 0] <- "A"
G8c[G8c == 1] <- "B"
G9c[G9c == 0] <- "A"
G9c[G9c == 1] <- "B"
G10c[G10c == 0] <- "A"
G10c[G10c == 1] <- "B"
Height.all <- cbind(Height,lwc,lnc,pac,G1c,G2c,G3c,G4c,G5c,G6c,G7c,G8c,G9c,G10c)

chro <- c("", "", "", "",1,1,1,2,2,2,2,3,3,3)
pos <- c("", "", "", "",0.3,0.6,0.8,1.3,1.4,1.6,1.9,2.4,2.5,2.7)
chrompos <- rbind(chro,pos)
Height.all.cross <- rbind(chrompos,Height.all)
Height.all.cross <- as.data.frame(Height.all.cross)
write_csv(Height.all.cross, "heightall_cross.csv")

vQTLsim <- read.cross(file = "heightall_cross.csv" )
vQTLsim <- drop.nullmarkers(vQTLsim)
vQTLsim <- calc.genoprob(vQTLsim)

sim.scan <- scanonevar(cross = vQTLsim,
                       mean.formula = Height ~ lwc + mean.QTL.add + lwc*mean.QTL.add,
                       var.formula = ~  lwc + var.QTL.add + lwc*var.QTL.add,
                       return.covar.effects = TRUE)

write.csv(sim.scan$result, file = "sim_scan_results_add(minus).csv")

sim.scan1 <- scanone(cross = vQTLsim)

#using the significant genomes
#sub_data < c()
sub_data <- cbind.data.frame(Maindata$Height,Maindata$Low.Water,Maindata$Low.Nitrogen,Maindata$Pathogen,
                             Maindata$Env,Maindata$IDP3914,Maindata$mmp177c,Maindata$lim333,Maindata$gpm662b,
                             Maindata$bnl5.46c,Maindata$IDP2349,Maindata$gpm296c,Maindata$IDP574,Maindata$umc1450,
                             Maindata$mmp46)
col_headings <- c('Height','Low.Water','Low.Nitrogen','Env','Pathogen', 'IDP3914', 'mmp177c','lim333',
                  'gpm662b','bnl5.46c','IDP2349','gpm296c','IDP574','umc1450','mmp46')
#sub_data <- as.data.frame(sub_data)
names(sub_data) <- col_headings
sub_data[is.na(sub_data)] <- ""

sub_data$Low.Water <- as.factor(sub_data$Low.Water)
sub_data$Low.Nitrogen <- as.factor(sub_data$Low.Nitrogen)
sub_data$Pathogen <- as.factor(sub_data$Pathogen)
sub_data$Env <- as.factor(sub_data$Env)

write_csv(sub_data, "ManchingSigData.csv")

vQTLsub <- read.cross(file = "ManchingSigData.csv" )
vQTLsub <- drop.nullmarkers(vQTLsub)
vQTLsub <- calc.genoprob(vQTLsub)

sub.scan <- scanonevar(cross = vQTLsub,
                       mean.formula = Height ~ Env + mean.QTL.add + Env*mean.QTL.add,
                       var.formula = ~ Low.Water + var.QTL.add + Low.Water*var.QTL.add,
                       return.covar.effects = TRUE)

write.csv(sub.scan$result, file = "sub_scan_results_int_LW.csv")

sub.scan1 <- scanone(cross = vQTLsub)

subandheight_data <- cbind.data.frame(Height.all.cross,sub_datasub_data[5:14])


