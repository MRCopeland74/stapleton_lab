

#### parameters #######
#### parameters #######

#### should we account for missing value?  
#### should there be more possible value variaties, 
#### should there be different prob for different variables?

### main function part #####
### main function part #####
simu.dglm<-function(n.loci=8, which.mean.loci=c(1:2), hypo.mean.para=c(1,4), which.var.loci=c(7:8),
                    hypo.var.para=c(2,5), loci.prob=c(0.5), N.obs=10^3, simu.prob.var=rep(0.5, n.loci), ran.seed=12345)
{
  # require(dglm);
 set.seed(ran.seed);

 n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
 n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)

 
var.names<-apply(expand.grid('effvar', c(1:n.loci)), 1, paste, collapse=".")

simu.loci<-NULL;
for(i in 1:length(var.names)) simu.loci<-cbind(simu.loci,rbinom(N.obs, 1,simu.prob.var[i]))
colnames(simu.loci)<-var.names;

all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;

mean.pyntp<-simu.loci%*%all.mean.para;
var.pyntp<-sapply(simu.loci%*%all.var.para,FUN=function(x) rnorm(1,sd=x));
pyntp<-mean.pyntp+var.pyntp;


simu.obs<-cbind(pyntp,simu.loci);
colnames(simu.obs)[1]<-"pyntp"
### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 🧬🧬🧬
### or, would a gene 🧬🧬🧬 structure(a composition of mutiple A' 🅰️ s and B' 🅱️s from gene🧬🧬🧬 loci) be an important thing to consider, or just with gene loci is contributing to the pynotype?
simu.obs.df<-data.frame(simu.obs);
dglm.fit<-dglm::dglm(pyntp~.,dformula= ~.,family=stats::gaussian,data=simu.obs.df,dlink='log');
return(list(dglm.fit, simu.obs));
}

out<-simu.dglm();
summary(out[[1]]);
obs<-out[[2]];
which.col<-1;
type.a.obs<-which(obs[,which.col+1]==0)

plot(rep(1,length(type.a.obs)),obs[type.a.obs,1]);

out<-simu.dglm(which.mean.loci=c(1,2,3), hypo.mean.para=c(1,4,8), which.var.loci=c(3,5,7),
               hypo.var.para=c(2,5,4),);
summary(out[[1]]);
### how to find the interaction term when the main is not significant👍� �😩😊👌❤️😁😉👀🙂😄
### netwaokr, use edge for the interaction effect.


### Austin 👨: The random 🤷 and fixed effect
### Michael 🧑: The negative 👎 outcome 
