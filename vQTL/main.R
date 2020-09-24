
#### parameters #######
#### parameters #######

#### should we account for missing value?  
#### should there be more possible value variaties, 
#### should there be different prob for different variables?

### main function part #####
### main function part #####
simu.dglm<-function(n.rep=3,n.obs.per.rep=150, 
                    n.loci=3235, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), 
                              which.var.loci=c(4:5),  hypo.var.para=c(2,5), 
                    loci.prob=c(0.5),  simu.prob.var=rep(0.5, n.loci), 
                    prop.env=0.3, 
                    ran.seed=12345)
{
  # require(dglm);
  # n.loci, number of gene loci
  # n.rep, the repetitation of each gene combination. The total numner of observation is n.rep*n.obs.per.rep
 set.seed(ran.seed);

 n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
 n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)

 
var.names<-apply(expand.grid('effvar', c(1:n.loci)), 1, paste, collapse=".");

simu.loci<-NULL;
for(i in 1:length(var.names)) simu.loci<-cbind(simu.loci,rep(rbinom(n.obs.per.rep, 1,simu.prob.var[i]),rep(n.rep,n.obs.per.rep)))
colnames(simu.loci)<-var.names;

all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;

mean.pyntp<-simu.loci%*%all.mean.para;
var.pyntp<-sapply(simu.loci%*%all.var.para,FUN=function(x) rnorm(1,sd=sqrt(exp(x))));
pyntp<-mean.pyntp+var.pyntp;


simu.obs<-cbind(pyntp,simu.loci);
colnames(simu.obs)[1]<-"pyntp"
### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 🧬🧬🧬
### or, would a gene 🧬🧬🧬 structure(a composition of mutiple A' 🅰️ s and B' 🅱️s from gene🧬🧬🧬 loci) be an important thing to consider, or just with gene loci is contributing to the pynotype?
simu.obs.df<-data.frame(simu.obs);
dglm.fit<-dglm::dglm(pyntp~.,dformula= ~.,data=simu.obs.df);#family=stats::gaussian,
return(list(dglm.fit, simu.obs));
}

out<-simu.dglm()
summary(out[[1]])

df = out[[2]]
df=as.data.frame(df)
df.breed = df$pyntp 
df.breed = df.breed[1:450]
df.breed = cbind.data.frame(df.breed,"inbred")
df.breed[321:450,2] = "Hybrid"
names(df.breed) = c("Stress", "Breedtype")




df.breed = cbind.data.frame(df.breed, df[,2:3236])
df.breed[] <- lapply( df.breed, factor)
df.breed$Stress = as.numeric(as.character(df.breed$Stress))
df.breed[,3:3237] =ifelse(df.breed[,3:3237]==1,"A","B")

setwd("/Users/michaelcopeland/Stapleton/Copeland/stapleton_lab/qPCR2vQTL")
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = FALSE)
hybrid.inbred = hybrid.inbred[,1:3240]
hybrid.inbred = hybrid.inbred[,-c(2,3,5)]
names(df.breed) = hybrid.inbred[1,]
hybrid.inbred = read.csv("Fullinb&hyb.csv", header = TRUE)
hybrid.inbred = hybrid.inbred
hybrid.inbred = hybrid.inbred[,-c(2,3,5)]
df.breed = rbind.data.frame(as.matrix(hybrid.inbred[1:2,]),as.matrix(df.breed))
df.breed = as.data.frame(df.breed)

df.breed$stress = as.character(df.breed$stress)
rownames(df.breed) <- NULL;
df.breed[1,1] = ""
df.breed[2,1] = ""

#df.breed = df.breed[,-c(2,3,5)]
write.csv(df.breed, "simdata.csv", row.names = FALSE)
library(gridExtra)
library(qtl)
library(vqtl)
test_full <- read.cross(file = "simdata.csv", format = "csv")

test_full <- drop.nullmarkers(test_full)
test_full <- calc.genoprob(test_full)

test_full <- calc.genoprob(test_full, error.prob = .001)
#hy_p1 <- scanone(cross = test_full, pheno.col = 'stress')
simdata <- scanonevar(cross = test_full, 
                           mean.formula = stress ~ BreedType+mean.QTL.add, 
                           var.formula = ~ BreedType+var.QTL.add, 
                           return.covar.effects = TRUE)




