library(qtl)
library(vqtl)
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")
sample <-read.cross(file ="ManchingStressData_Covar.csv",na.strings = "-")

sample <- drop.nullmarkers(sample)
#scan with variance
sample <- calc.genoprob(sample)

head(sample$pheno)

# turn looks like 'Env' is categorical data -- in R, we represent that with a 'factor'
sample$pheno$Env <- factor(sample$pheno$Env)

# let's make a ggplot to have a look
library(tidyverse)
sample$pheno %>%
  ggplot(mapping = aes(x = Env, y = Height)) +
  geom_jitter(width = 0.2)

# don't think there's enough data here to say...looks like similar mean and variance? - Robert

outv <- scanonevar(cross = sample,
                   mean.formula = Height ~ Env + mean.QTL.add + mean.QTL.dom,
                   var.formula = ~ Env + var.QTL.add + var.QTL.dom,
                   return.covar.effects = TRUE)

outv$result %>% glimpse()
write.csv(outv$result, file = "Manching_additive_model.csv")

outv <- scanonevar(cross = sample,
                   mean.formula = Height ~ Env * (mean.QTL.add + mean.QTL.dom),
                   var.formula = ~ Env * (var.QTL.add + var.QTL.dom),
                   return.covar.effects = TRUE)

outv$result %>% glimpse()
write.csv(outv$result, file = "Manching_interactive_model.csv")