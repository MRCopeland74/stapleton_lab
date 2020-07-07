
#### parameters #######
#### parameters #######

#### should we account for missing value?  
#### should there be more possible value variaties, 
#### should there be different prob for different variables?

### main function part #####
### main function part #####
simu.dglm<-function(n.rep=3,n.obs.per.rep=10^3, 
                    n.loci=8, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), 
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

out<-simu.dglm();
summary(out[[1]]);
