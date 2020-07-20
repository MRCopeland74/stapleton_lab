

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
                    prop.env=0.3, n.env = 3, int.pair = c(3,5),int.function = NULL, 
                    int.formula = NULL, 
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
  
  env.var<-NULL;
  env.var<-cbind(env.var,sample(1:n.env,n.obs.per.rep, replace = TRUE))
  env.var.new = NULL
  for (i in 1:length(env.var)){
    env.var.new = c(env.var.new, rep(env.var[i],3))
  }
  env.var = as.matrix(env.var.new)
  colnames(env.var)<-"env";
  
  all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
  all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;
  
  mean.pyntp<-simu.loci%*%all.mean.para;
  var.pyntp<-sapply(simu.loci%*%all.var.para,FUN=function(x) rnorm(1,sd=sqrt(exp(x))));
  pyntp<-mean.pyntp+var.pyntp;
  
  
  simu.obs<-cbind(pyntp,env.var,simu.loci);
  
  if (is.null(int.function) == FALSE){
    int.function = match.fun(int.function)
    for(i in 1:n.obs.per.rep) {
      if (simu.obs[i,2] == int.pair[1] & simu.obs[i,int.pair[2]+2]==1){
        simu.obs[i,1] = int.function(simu.obs[i,1])
      }
      if (simu.obs[i,2]== int.pair[1]){
        #simu.obs[i,1] = simu.obs[i,1]+0.5
      }
      if (simu.obs[i,int.pair[2]]==1){
        #simu.obs[i,1] = simu.obs[i,1]-1
      }
    }
  }
  
  if (is.null(int.formula) == FALSE){
    for(i in 1:n.obs.per.rep) {
      if (simu.obs[i,2] == int.pair[1] & simu.obs[i,int.pair[2]+2]==1){
        simu.obs[i,1] = simu.obs[i,1]
      }
      if (simu.obs[i,2]== int.pair[1]){
        #simu.obs[i,1] = simu.obs[i,1]+0.5
      }
      if (simu.obs[i,int.pair[2]]==1){
        #simu.obs[i,1] = simu.obs[i,1]-1
      }
    }
  }
  
  
  colnames(simu.obs)[1]<-"pyntp"
  ### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 
  ### or, would a gene  structure(a composition of mutiple A' from genoci) be an important thing to consider, or just with gene loci is contributing to the pynotype?
  simu.obs.df<-data.frame(simu.obs);
  dglm.fit<-dglm::dglm(pyntp~.,dformula= ~.,data=simu.obs.df);#family=stats::gaussian,
  return(list(dglm.fit, simu.obs));
}

out<-simu.dglm();
summary(out[[1]]);
obs<-out[[2]];
which.col<-1;
type.a.obs<-which(obs[,which.col+1]==0)

plot(rep(1,length(type.a.obs)),obs[type.a.obs,1]);

out<-simu.dglm(which.mean.loci=c(1,2,3), hypo.mean.para=c(1,4,8), which.var.loci=c(3,5,7), hypo.var.para=c(2,5,4),);
summary(out[[1]]);
### how to find the interaction term when the main is not significant?
### netwaokr, use edge for the interaction effect.

