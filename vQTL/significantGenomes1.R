library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")
Maindata = read.csv("ManchingStressData_Covar.csv")
#https://github.com/MRCopeland74/stapleton_lab/blob/master/vQTL/ManchingStressData_Covar.csv
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

count2 = 1
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

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
#-------------------------------------------Boxplot code---------------------------------------#
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
genet10 = top10[,2]
p = c()
count3 = 1
for (i in top10){
  print(i)
  p[count3] <- boxplot(Maindata$Height~droplevels(Maindata$IDP3914, exclude = '-'))
  count = count+1
}
Maindata <- Maindata[-c(1, 2), ]
boxplot(Maindata$Height~droplevels(Maindata$bnl5.46c, exclude = '-',))

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
#---------Simulated data with every possible combinantion of environment with gene-------------#
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#
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



G1f=G2f=G3f=G4f=G5f=G6f=G7f=G8f=G9f=G10f <- c()
Heightf = c()
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
                          
                          # G1f[count4]<- w1
                          # G2f[count4]<- w2
                          # G3f[count4]<- w3
                          # G4f[count4]<- w4
                          # G5f[count4]<- w5
                          # G6f[count4]<- w6
                          # G7f[count4]<- w7
                          # G8f[count4]<- w8
                          # G9f[count4]<- w9
                          # G10f[count4]<- w10
                          
                          lwc[count4]<- i
                          lnc[count4]<- o
                          pac[count4]<- u
                          Height[count4] <- a -0.9*i -0.95*o -1.1*u +1.5*w1 + 50*w2 +0.7*w3 + 1.3*w4 +1.2*w5 
                          + 0.6*w6 +0.8*w7 + 1.1*w8 +1.4*w9 + 0.5*w10 
                          + se[sample(1:1000, 1)]
                          
                          # Heightf[count4] <- a -0.9*i -0.95*o -1.1*u +0.7*w1 + 0.4*w2 +0.95*w3 + 1.8*w4 +1.11*w5 
                          # + 0.78*w6 +0.65*w7 + 1.5*w8 +1.9*w9 + 0.88*w10 
                          # + se[sample(1:1000, 1)]
                          
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

#----------------------------------Box PLot---------------------------------------------#
boxplot.simulated = c()
Height.all = as.data.frame(Height.all)
for (i in c(5:14)){
  boxplot.simulated[i-4]=boxplot(as.numeric(Height.all$Height)~droplevels(Height.all[,i]))
}
#---------------------------------------------------------------------------------------#
write_csv(Height.all.cross, "heightall_cross.csv")


vQTLsim <- read.cross(file = "heightall_cross.csv" )
vQTLsim <- drop.nullmarkers(vQTLsim)
vQTLsim <- calc.genoprob(vQTLsim)

sim.scan <- scanonevar(cross = vQTLsim,
                       mean.formula = Height ~ lwc + mean.QTL.add + lwc*mean.QTL.add,
                       var.formula = ~  lwc + var.QTL.add + lwc*var.QTL.add,
                       return.covar.effects = TRUE)

write.csv(sim.scan$result, file = "sim_scan_results_add(minus).csv")

sim.scan1 <- scanone(cross = vQTLsim, intcovar = lwc)

#using the significant genomes
#sub_data < c()
sub_data <- cbind.data.frame(Maindata$Height,Maindata$Low.Water,Maindata$Low.Nitrogen,Maindata$Pathogen,
                             Maindata$Env,Maindata$IDP3914,Maindata$mmp177c,Maindata$lim333,Maindata$gpm662b,
                             Maindata$bnl5.46c,Maindata$IDP2349,Maindata$gpm296c,Maindata$IDP574,Maindata$umc1450,
                             Maindata$mmp46)
col_headings <- c('Height','Low.Water','Low.Nitrogen','Pathogen','Env', 'IDP3914', 'mmp177c','lim333',
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
                       mean.formula = Height ~ (Low.Water + Low.Nitrogen + Pathogen)*(mean.QTL.add),
                       var.formula = ~ (Low.Water + Low.Nitrogen + Pathogen)*(var.QTL.add),
                       return.covar.effects = TRUE)

write.csv(sub.scan$result, file = "sub_scan_results_int_LW.csv")

sub.scan1 <- scanone(cross = vQTLsub, intcovar = 7)

#subandheight_data <- cbind.data.frame(Height.all.cross,sub_datasub_data[5:14])
G1f = c()
G2f = c()
G3f = c()
G4f = c()
G5f = c()
G6f = c()
G7f = c()
G8f = c()
G9f = c()
G10f = c()


for (i in 1:8192){
  G1f[i]=sample(0:1,1)
  G2f[i]=sample(0:1,1)
  G3f[i]=sample(0:1,1)
  G4f[i]=sample(0:1,1)
  G5f[i]=sample(0:1,1)
  G6f[i]=sample(0:1,1)
  G7f[i]=sample(0:1,1)
  G8f[i]=sample(0:1,1)
  G9f[i]=sample(0:1,1)
  G10f[i]=sample(0:1,1)
}


G1f[G1f == 0] <- "A"
G1f[G1f == 1] <- "B"
G2f[G2f == 0] <- "A"
G2f[G2f == 1] <- "B"
G3f[G3f == 0] <- "A"
G3f[G3f == 1] <- "B"
G4f[G4f == 0] <- "A"
G4f[G4f == 1] <- "B"
G5f[G5f == 0] <- "A"
G5f[G5f == 1] <- "B"
G6f[G6f == 0] <- "A"
G6f[G6f == 1] <- "B"
G7f[G7f == 0] <- "A"
G7f[G7f == 1] <- "B"
G8f[G8f == 0] <- "A"
G8f[G8f == 1] <- "B"
G9f[G9f == 0] <- "A"
G9f[G9f == 1] <- "B"
G10f[G10f == 0] <- "A"
G10f[G10f == 1] <- "B"

#Heightf <- rnorm(8192,mean = 200, sd =20)

Heightf.all <- cbind(G1f,G2f,G3f,G4f,G5f,G6f,G7f,G8f,G9f,G10f)

Heightf.all.cross <- cbind.data.frame(Height.all,Heightf.all)

Heightf.all.cross$Height <- as.numeric(as.character(Heightf.all.cross$Height))
chro <- c("", "", "", "",1,1,1,2,2,2,2,3,3,3,1,1,1,2,2,2,2,3,3,3)
pos <- c("", "", "", "",0.3,0.6,0.8,1.3,1.4,1.6,1.9,2.4,2.5,2.7,3.1,3.3,3.7,4.1,4.2,4.7,5.2,5.4,5.8,6.4)
chrompos <- rbind(chro,pos)
Heightf.all.cross <- as.matrix(Heightf.all.cross)
Heightf.all.cross <- rbind(chrompos,Heightf.all.cross)
Heightf.all.cross <- as.data.frame(Heightf.all.cross)


write_csv(Heightf.all.cross, "heightall_f_cross.csv")

vQTLsimf <- read.cross(file = "heightall_f_cross.csv" )
vQTLsimf <- drop.nullmarkers(vQTLsimf)
vQTLsimf <- calc.genoprob(vQTLsimf)

sim.scanf <- scanonevar(cross = vQTLsimf,
                        mean.formula = Height ~ lwc + mean.QTL.add + lwc*mean.QTL.add,
                        var.formula = ~  lwc + var.QTL.add + lwc*var.QTL.add,
                        return.covar.effects = TRUE)

#write.csv(sim.scan$result, file = "sim_scan_results_int_f.csv")

sim.scanf1 <- scanone(cross = vQTLsimf)


fake = subset(sim.scanf$result$mQTL.lod, grepl("f",sim.scanf$result$loc.name))

real = subset(sim.scanf$result$mQTL.lod, grepl("c",sim.scanf$result$loc.name)) 
Loci = c("01","02","03","04","05","06","07","08","09","10")
fake = cbind.data.frame(fake,genes)
real = cbind.data.frame(real,genes)
colnames(fake) = c("LOD","Loci")
colnames(real) = c("LOD","Loci")

str(real)
ggplot()+geom_point(data=fake,aes(x=Loci,y=LOD,color = "Trojan"), position = position_jitter()) +
  geom_point(data=real,aes(x=Loci,y=LOD,color = "Real"), position = position_jitter()) +
  labs(title = "Simulated vs Trojan") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------Difference in height per env factor of A and B--------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
Maindata$Env = as.factor(Maindata$Env)

maindata.env1 = Maindata[Maindata$Env == "1",]
maindata.env2 = Maindata[Maindata$Env == "2",]
maindata.env3 = Maindata[Maindata$Env == "3",]
maindata.env4 = Maindata[Maindata$Env == "4",]
maindata.env5 = Maindata[Maindata$Env == "5",]
maindata.env6 = Maindata[Maindata$Env == "6",]
maindata.env7 = Maindata[Maindata$Env == "7",]
maindata.env8 = Maindata[Maindata$Env == "8",]

##avg.height.chr1AB <- tapply(maindata.chr$Height,maindata.chr$gpm27, mean)[4:5]
#avg.height.chr1 <- mean(na.omit(Maindata$Height))

count = 1
avector.env1 = c()
avector.env2 = c()
avector.env3 = c()
avector.env4 = c()
avector.env5 = c()
avector.env6 = c()
avector.env7 = c()
avector.env8 = c()
for (i in colm){
  # print(i)
  avector.env1[count] <- as.vector(maindata.env1[i])
  avector.env2[count] <- as.vector(maindata.env2[i])
  avector.env3[count] <- as.vector(maindata.env3[i])
  avector.env4[count] <- as.vector(maindata.env4[i])
  avector.env5[count] <- as.vector(maindata.env5[i])
  avector.env6[count] <- as.vector(maindata.env6[i])
  avector.env7[count] <- as.vector(maindata.env7[i])
  avector.env8[count] <- as.vector(maindata.env8[i])
  count = count+1
}


AandB.env1 = c()
AandB.env1AB = c()
AandB.env2 = c()
AandB.env2AB = c()
AandB.env3 = c()
AandB.env3AB = c()
AandB.env4 = c()
AandB.env4AB = c()
AandB.env5 = c()
AandB.env5AB = c()
AandB.env6 = c()
AandB.env6AB = c()
AandB.env7 = c()
AandB.env7AB = c()
AandB.env8 = c()
AandB.env8AB = c()
#count2 = 1


for (i in 1:(length(colm)-1)){
  #print(i)
  if (length(levels(avector.env1[[i]]))==5){
    temp1AB = tapply(maindata.env1$Height,avector.env1[[i]], mean)[4:5]
  }
  if (length(levels(avector.env1[[i]]))==4){
    temp1AB = tapply(maindata.env1$Height, avector.env1[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env2[[i]]))==5){
    temp2AB = tapply(maindata.env2$Height,avector.env2[[i]], mean)[4:5]
  }
  if (length(levels(avector.env2[[i]]))==4){
    temp2AB = tapply(maindata.env2$Height, avector.env2[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env3[[i]]))==5){
    temp3AB = tapply(maindata.env3$Height,avector.env3[[i]], mean)[4:5]
  }
  if (length(levels(avector.env3[[i]]))==4){
    temp3AB = tapply(maindata.env3$Height, avector.env3[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env4[[i]]))==5){
    temp4AB = tapply(maindata.env4$Height,avector.env4[[i]], mean)[4:5]
  }
  if (length(levels(avector.env4[[i]]))==4){
    temp4AB = tapply(maindata.env4$Height, avector.env4[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env5[[i]]))==5){
    temp5AB = tapply(maindata.env5$Height,avector.env5[[i]], mean)[4:5]
  }
  if (length(levels(avector.env5[[i]]))==4){
    temp5AB = tapply(maindata.env5$Height, avector.env5[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env6[[i]]))==5){
    temp6AB = tapply(maindata.env6$Height,avector.env6[[i]], mean)[4:5]
  }
  if (length(levels(avector.env6[[i]]))==4){
    temp6AB = tapply(maindata.env6$Height, avector.env6[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env7[[i]]))==5){
    temp7AB = tapply(maindata.env7$Height,avector.env7[[i]], mean)[4:5]
  }
  if (length(levels(avector.env7[[i]]))==4){
    temp7AB = tapply(maindata.env7$Height, avector.env7[[i]], mean)[3:4] 
  }
  if (length(levels(avector.env8[[i]]))==5){
    temp8AB = tapply(maindata.env8$Height,avector.env8[[i]], mean)[4:5]
  }
  if (length(levels(avector.env8[[i]]))==4){
    temp8AB = tapply(maindata.env8$Height,avector.env8[[i]], mean)[3:4]
  }
  #=count2+1
  AandB.env1AB = rbind(AandB.env1AB,temp1AB)
  AandB.env2AB = rbind(AandB.env2AB,temp2AB)
  AandB.env3AB = rbind(AandB.env3AB,temp3AB)
  AandB.env4AB = rbind(AandB.env4AB,temp4AB)
  AandB.env5AB = rbind(AandB.env5AB,temp5AB)
  AandB.env6AB = rbind(AandB.env6AB,temp6AB)
  AandB.env7AB = rbind(AandB.env7AB,temp7AB)
  AandB.env8AB = rbind(AandB.env8AB,temp8AB)
}
temp1AB = tapply(maindata.env1$Height,avector.env1[[3235]], mean)[3:4]
temp2AB = tapply(maindata.env2$Height,avector.env2[[3235]], mean)[3:4]
temp3AB = tapply(maindata.env3$Height,avector.env3[[3235]], mean)[3:4]
temp4AB = tapply(maindata.env4$Height,avector.env4[[3235]], mean)[3:4]
temp5AB = tapply(maindata.env5$Height,avector.env5[[3235]], mean)[3:4]
temp6AB = tapply(maindata.env6$Height,avector.env6[[3235]], mean)[3:4]
temp7AB = tapply(maindata.env7$Height,avector.env7[[3235]], mean)[3:4]
temp8AB = tapply(maindata.env8$Height,avector.env8[[3235]], mean)[3:4]
AandB.env1AB = rbind(AandB.env1AB,temp1AB)
AandB.env2AB = rbind(AandB.env2AB,temp2AB)
AandB.env3AB = rbind(AandB.env3AB,temp3AB)
AandB.env4AB = rbind(AandB.env4AB,temp4AB)
AandB.env5AB = rbind(AandB.env5AB,temp5AB)
AandB.env6AB = rbind(AandB.env6AB,temp6AB)
AandB.env7AB = rbind(AandB.env7AB,temp7AB)
AandB.env8AB = rbind(AandB.env8AB,temp8AB)
temp1 = mean(maindata.env1$Height)
temp2 = mean(maindata.env2$Height)
temp3 = mean(maindata.env3$Height)
temp4 = mean(maindata.env4$Height)
temp5 = mean(maindata.env5$Height)
temp6 = mean(maindata.env6$Height)
temp7 = mean(maindata.env7$Height)
temp8 = mean(maindata.env8$Height)


AandB.env1AB = cbind.data.frame(colm,AandB.env1AB)
AandB.env1 = cbind.data.frame(AandB.env1AB,(AandB.env1AB[,2]-AandB.env1AB[,3]))
AandB.env1AB = cbind.data.frame(AandB.env1AB,abs((AandB.env1AB[,2]-AandB.env1AB[,3])))
AandB.env2AB = cbind.data.frame(colm,AandB.env2AB)
AandB.env2 = cbind.data.frame(AandB.env2AB,(AandB.env2AB[,2]-AandB.env2AB[,3]))
AandB.env2AB = cbind.data.frame(AandB.env2AB,abs((AandB.env2AB[,2]-AandB.env2AB[,3])))
AandB.env3AB = cbind.data.frame(colm,AandB.env3AB)
AandB.env3 = cbind.data.frame(AandB.env3AB,(AandB.env3AB[,2]-AandB.env3AB[,3]))
AandB.env3AB = cbind.data.frame(AandB.env3AB,abs((AandB.env3AB[,2]-AandB.env3AB[,3])))
AandB.env4AB = cbind.data.frame(colm,AandB.env4AB)
AandB.env4 = cbind.data.frame(AandB.env4AB,(AandB.env4AB[,2]-AandB.env4AB[,3]))
AandB.env4AB = cbind.data.frame(AandB.env4AB,abs((AandB.env4AB[,2]-AandB.env4AB[,3])))
AandB.env5AB = cbind.data.frame(colm,AandB.env5AB)
AandB.env5 = cbind.data.frame(AandB.env5AB,(AandB.env5AB[,2]-AandB.env5AB[,3]))
AandB.env5AB = cbind.data.frame(AandB.env5AB,abs((AandB.env5AB[,2]-AandB.env5AB[,3])))
AandB.env6AB = cbind.data.frame(colm,AandB.env6AB)
AandB.env6 = cbind.data.frame(AandB.env6AB,(AandB.env6AB[,2]-AandB.env6AB[,3]))
AandB.env6AB = cbind.data.frame(AandB.env6AB,abs((AandB.env6AB[,2]-AandB.env6AB[,3])))
AandB.env7AB = cbind.data.frame(colm,AandB.env7AB)
AandB.env7 = cbind.data.frame(AandB.env7AB,(AandB.env7AB[,2]-AandB.env7AB[,3]))
AandB.env7AB = cbind.data.frame(AandB.env7AB,abs((AandB.env7AB[,2]-AandB.env7AB[,3])))
AandB.env8AB = cbind.data.frame(colm,AandB.env8AB)
AandB.env8 = cbind.data.frame(AandB.env8AB,(AandB.env8AB[,2]-AandB.env8AB[,3]))
AandB.env8AB = cbind.data.frame(AandB.env8AB,abs((AandB.env8AB[,2]-AandB.env8AB[,3])))

names(AandB.env1AB) = c("Gene","A","B","Difference")
names(AandB.env2AB) = c("Gene","A","B","Difference")
names(AandB.env3AB) = c("Gene","A","B","Difference")
names(AandB.env4AB) = c("Gene","A","B","Difference")
names(AandB.env5AB) = c("Gene","A","B","Difference")
names(AandB.env6AB) = c("Gene","A","B","Difference")
names(AandB.env7AB) = c("Gene","A","B","Difference")
names(AandB.env8AB) = c("Gene","A","B","Difference")

names(AandB.env1) = c("Gene","A","B","Difference")
names(AandB.env2) = c("Gene","A","B","Difference")
names(AandB.env3) = c("Gene","A","B","Difference")
names(AandB.env4) = c("Gene","A","B","Difference")
names(AandB.env5) = c("Gene","A","B","Difference")
names(AandB.env6) = c("Gene","A","B","Difference")
names(AandB.env7) = c("Gene","A","B","Difference")
names(AandB.env8) = c("Gene","A","B","Difference")
sum.env <- AandB.env1AB$Difference + AandB.env2AB$Difference + AandB.env3AB$Difference +
  AandB.env4AB$Difference + AandB.env5AB$Difference + AandB.env6AB$Difference + 
  AandB.env7AB$Difference + AandB.env8AB$Difference
sum.env = cbind.data.frame(AandB.env1AB$Gene, sum.env)
names(sum.env)=c("Gene","Sum/env")
sum.env.top10 = sum.env[order(sum.env[,"Sum/env"],decreasing = TRUE)[1:10],]

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------Number A and B per row--------------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#

num.AandB.row = cbind.data.frame(Maindata$Height, Maindata$Env,rowSums(Maindata == "A"),
                                 rowSums(Maindata == "B"), rowSums(Maindata=="-"))
num.AandB.row = na.omit(num.AandB.row)
names(num.AandB.row) = c("Height", "Env", "A","B","-")
temp.env = num.AandB.row$Env[1]
temp.A = num.AandB.row$A[1]
temp.B = num.AandB.row$B[1]
temp.H = num.AandB.row$`-`[1]
count = 0
height.avg = 0
num.AandB.row.avg = c()
for (i in 1:6672){
  print("school")
  if (temp.env == num.AandB.row$Env[i] && temp.A == num.AandB.row$A[i] && 
      temp.B == num.AandB.row$B[i] && temp.H == num.AandB.row$`-`[i]){
    height.avg = height.avg + num.AandB.row$Height[i]
    count = count + 1
    print(height.avg)
  }else{
    temp.row = cbind.data.frame(height.avg/count, temp.env, temp.A, temp.B, temp.H)
    num.AandB.row.avg = rbind.data.frame(num.AandB.row.avg, temp.row)
    temp.row = c()
    temp.env = num.AandB.row$Env[i]
    temp.A = num.AandB.row$A[i]
    temp.B = num.AandB.row$B[i]
    temp.H = num.AandB.row$`-`[i]
    count = 1
    height.avg = num.AandB.row$Height[i]
    print(temp.env)
  }
}
temp.row = cbind.data.frame(height.avg/count, temp.env, temp.A, temp.B, temp.H)
num.AandB.row.avg = rbind.data.frame(num.AandB.row.avg, temp.row)
plot(num.AandB.row.avg$temp.A, num.AandB.row.avg$`height.avg/count`)
plot(num.AandB.row.avg$temp.B, num.AandB.row.avg$`height.avg/count`)
plot(num.AandB.row.avg$temp.H, num.AandB.row.avg$`height.avg/count`)


#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------Percentage of genes per chromosone--------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
Maindata.t = t(Maindata)
Maindata.t = as.data.frame(Maindata.t)
Maindata.t = table(Maindata.t$V1)
perc.chrom <- Maindata.t/3235
Maindata.t = cbind.data.frame(Maindata.t, perc.chrom)

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------Crossover from 2-3 top 10-----------------------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#Continuation of Difference in height per env factor of A and B
Difference.big = cbind.data.frame(AandB.env2$Gene,abs(AandB.env2$Difference - AandB.env3$Difference))
A.big = cbind.data.frame(AandB.env2$Gene,abs(AandB.env2$A - AandB.env3$A))
B.big = cbind.data.frame(AandB.env2$Gene,abs(AandB.env2$B - AandB.env3$B))
names(Difference.big)=c("Gene", "Difference2to3")
names(A.big)=c("Gene", "A2to3")
names(B.big)=c("Gene", "B2to3")
Difference.big.ten = Difference.big[order(Difference.big[,"Difference2to3"],decreasing = TRUE)[1:10],]
A.big.ten = A.big[order(A.big[,"A2to3"],decreasing = TRUE)[1:10],]
B.big.ten = B.big[order(B.big[,"B2to3"],decreasing = TRUE)[1:10],]

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
#--------------------Grouping significant genomes and running vQTL on them--------------------#
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
count.env.group=1
sum.env.group10=c()
num = 10
p=1
q=num
sum.env.order = sum.env[order(sum.env[,"Sum/env"],decreasing = TRUE),]
for(i in 1:(length(sum.env.order[,1])%/%num)){
  sum.env.group10[count.env.group] = sum.env.order[p:q,]
  p=p+num
  q=q+num
  count.env.group=count.env.group+1
}
sub.scan.frames = list()
count1 = 1
for (i in 1:1#length(sum.env.group10)
){
  col_headings = sum.env.group10[i][[1]]
  col_headings=matrix(col_headings)
  subdata1 = Maindata[col_headings[1:10]]
  sub_data <- cbind.data.frame(Maindata$Height,Maindata$Low.Water,Maindata$Low.Nitrogen,Maindata$Pathogen,
                               Maindata$Env,subdata1)
  col_headings <- c('Height','Low.Water','Low.Nitrogen','Pathogen','Env', col_headings)
  names(sub_data) <- col_headings
  sub_data$Env=as.character(sub_data$Env)
  sub_data[is.na(sub_data)] <- ""
  
  sub_data$Low.Water <- as.factor(sub_data$Low.Water)
  sub_data$Low.Nitrogen <- as.factor(sub_data$Low.Nitrogen)
  sub_data$Pathogen <- as.factor(sub_data$Pathogen)
  sub_data$Env <- as.factor(sub_data$Env)
  
  write_csv(sub_data, "ManchingSigData1.csv")
  
  vQTLsub1 <- read.cross(file = "ManchingSigData1.csv" )
  vQTLsub1 <- drop.nullmarkers(vQTLsub1)
  vQTLsub1 <- calc.genoprob(vQTLsub1)
  
  sub.scan <- scanonevar(cross = vQTLsub1,
                         mean.formula = Height ~ (Low.Water + Low.Nitrogen + Pathogen)*(mean.QTL.add),
                         var.formula = ~ (Low.Water + Low.Nitrogen + Pathogen)*(var.QTL.add),
                         return.covar.effects = TRUE)
  
  sub.scan.frames[[count1]] = sub.scan$result
  count1=count1+1
  
  
}


write.csv(sub.scan$result, file = "sub_scan_results_int_LW.csv")

sub.scan1 <- scanone(cross = vQTLsub1, intcovar = 7)





