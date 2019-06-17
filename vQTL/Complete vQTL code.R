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
setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/vQTL/")

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
addOneVar <- scanonevar(cross = fr,
                        mean.formula = Height ~ Env + mean.QTL.add + mean.QTL.dom,
                        var.formula = ~ Env + var.QTL.add + var.QTL.dom,
                        return.covar.effects = TRUE)

# Writing the result of the additive scanonevar for later use
write_rds(addOneVar, "addOneVar.rds", compress = "xz")
print("first scanonevar")

# Interactive scanonevar function
intOneVar <- scanonevar(cross = fr,
                        mean.formula = Height ~ Env * (mean.QTL.add + mean.QTL.dom),
                        var.formula = ~ Env * (var.QTL.add + var.QTL.dom),
                        return.covar.effects = TRUE)
print("second scanonevar")
# Writing the result of the interactive scanonevar for later use
write_rds(intOneVar, "intOneVar.rds", compress = "xz")


# Writing out the results of the two 
write.csv(addOneVar$result, file = "Manching_additive_model.csv")
write.csv(intOneVar$result, file = "Manching_interactive_model.csv")

effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL,
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL,
                         title = paste(phenotype.name, "by", paste(focal.groups, collapse = ", ")),
                         draw_ribbons = TRUE, se_line_size = 1,
                         point_size = 1){
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))],
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))],
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[,
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups,
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form),
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df),
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim -
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim),
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim +
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>%
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim),
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub),
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb),
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}

y = 1:length(outv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs
sizedf = matrix(nrow = length(outv$result$loc.name), ncol = 12)

map(y, function(x){
  tryCatch({
    print(x)
    tempm =  effect.sizes(cross = dat,
                          phenotype.name = "Height",
                          genotype.names = c("AA","BB"),
                          focal.groups = outv$result$loc.name[x])
    tempv = matrix(nrow = 1, ncol = 12)
    tempv[1,1:6] = as.numeric(tempm[1,2:7])
    tempv[1,7:12] = as.numeric(tempm[2,2:7])
    sizedf[x,] <<- tempv
  }, error = function(e){message(e)
    sizedf[x,] <<- rep(0,12)
    
  }
  )
})

sizedf1 <- matrix(nrow = length(y), ncol = 12)

map(y, function(x){
  sizedf1[x,] <<- sizedf[[x]]
  
})


nall0 <-sapply(1:dim(sizedf1)[1], function(x){
  !all(sizedf1[x,] == 0)
})
ditch <- which(nall0 == F)
sizedf1 <- sizedf1[-ditch,]
keep <- y; keep<- keep[-ditch]
outvdf<- data.frame(outv$result$loc.name[keep],
                    outv$result$pos[keep],
                    outv$result$mQTL.lod[keep],
                    outv$result$mQTL.asymp.p[keep],
                    outv$result$vQTL.lod[keep],
                    outv$result$vQTL.asymp.p[keep],
                    outv$result$mvQTL.lod[keep],
                    outv$result$mvQTL.asymp.p[keep],
                    outv$result$`m__est__(Intercept)`[keep],
                    outv$result$`m__se__(Intercept)`[keep],
                    outv$result$`v__est__(Intercept)`[keep],
                    outv$result$`v__se__(Intercept)`[keep],
                    outv$result$m__est__mean.QTL.add[keep],
                    outv$result$m__se__mean.QTL.add[keep],
                    outv$result$v__est__var.QTL.add[keep],
                    outv$result$v__est__var.QTL.add[keep],
                    outv$result$m__est__env[keep],
                    outv$result$m__se__env[keep])
outvdf = cbind(outvdf,sizedf1)
colnames(outvdf) = c("SNP Name",
                     "Position (cM)",
                     "Mean LOD",
                     "Mean P Value",
                     "Variance LOD",
                     "Variance P Value",
                     "Joint LOD",
                     "Joint P Value",
                     "Mean Est Intercept",
                     "Mean SE Intercept",
                     "Var Est Intercept",
                     "Var SE Intercept",
                     "Mean Est QTL",
                     "Mean SE QTL",
                     "Var Est QTL",
                     "Var SE QTL",
                     "Mean Est Env",
                     "Mean SE Env",
                     "A Mean Est",
                     "A Mean Lower Bound",
                     "A Mean Upper Bound",
                     "A Standard Deviation Est",
                     "A Standard Deviation Lower Bound",
                     "A Standard Deviation Upper Bound",
                     "B Mean Est",
                     "B Mean Lower Bound",
                     "B Mean Upper Bound",
                     "B Standard Deviation Est",
                     "B Standard Deviation Lower Bound",
                     "B Standard Deviation Upper Bound")

write.csv(outvdf, file = "2018_12_20_ManchingCovarCombined_vQTL.csv")