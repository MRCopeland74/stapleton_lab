

#### parameters #######
#### parameters #######

#### should we account for missing value?  
#### should there be more possible value variaties, 
#### should there be different prob for different variables?

### main function part #####
### main function part #####
library(dplyr)
simu.dglm.inter<-function(n.rep=3, n.level.env=6, n.obs.per.rep=10^3, 
                    n.loci=8, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), 
                              which.var.loci=c(4:5),  hypo.var.para=c(2,5), 
                              which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
                    loci.prob=c(0.5),  simu.prob.var=rep(0.5, n.loci), 
                    ran.seed=12345)
{
  # require(dglm);
  # n.loci, number of gene loci
  # n.rep, the repetitation of each gene combination. The total numner of observation is n.rep*n.level.env*n.obs.per.rep
 set.seed(ran.seed);
        
 n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
 n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)
 
 if (which.env.sig>n.level.env) {print("the which.env.sig should be less than n.level.env"); stop;}
 env<-rep(rep(c(1:n.level.env),rep(n.rep,n.level.env)),n.obs.per.rep);
 env.mat<-matrix(0,nrow=n.rep*n.level.env*n.obs.per.rep, ncol=n.level.env);
 for(i in 1:n.level.env)  env.mat[,i]=(env==i);
 var.names.env<-apply(expand.grid('env.var', c(1:n.level.env)), 1, paste, collapse=".");
 colnames(env.mat)<-var.names.env;
 
simu.loci<-NULL;
for(i in 1:n.loci) simu.loci<-cbind(simu.loci,rep(rbinom(n.obs.per.rep, 1,simu.prob.var[i]),rep(n.rep*n.level.env,n.obs.per.rep)))

var.names.loci<-apply(expand.grid('loci.var', c(1:n.loci)), 1, paste, collapse=".");
colnames(simu.loci)<-var.names.loci;

all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;
all.inter.para<-rep(0,n.level.env*n.loci);all.inter.para[which.env.inter+(which.loci.inter-1)*n.loci]<-hypo.inter.para;


mean.pyntp<-simu.loci%*%all.mean.para;
var.pyntp<-sapply(simu.loci%*%all.var.para,FUN=function(x) rnorm(1,sd=sqrt(exp(x))));




simu.obs<-cbind(env.mat[,-1],simu.loci); ## the -1 is to get rid of the frist level for collinearity.

### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 🧬🧬🧬
### or, would a gene 🧬🧬🧬 structure(a composition of mutiple A' 🅰️ s and B' 🅱️s from gene🧬🧬🧬 loci) be an important thing to consider, or just with gene loci is contributing to the pynotype?
inter.mat<-NULL;
for(i in n.level.env:1) for(j in n.loci:1) inter.mat<-cbind(simu.loci[,j]*env.mat[,i],inter.mat); ## stop at 2 is to get rid of the frist level for collinearity.

var.names.inter.mat<-apply(expand.grid(var.names.loci, var.names.env)%>%as.matrix,1,paste,collapse="*");
colnames(inter.mat)<-var.names.inter.mat;

all.inter.para<-rep(0,n.level.env*n.loci);
all.inter.para[which.loci.inter+(which.env.inter-1)*n.loci]<-hypo.inter.para;

#which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
inter.pyntp<-inter.mat%*%all.inter.para; ### this is wrong

simu.obs.df<-data.frame(simu.obs,inter.mat[,-c(1:n.loci)]);
simu.obs.df$pyntp<-mean.pyntp+var.pyntp+inter.pyntp;


dglm.fit<-dglm::dglm(pyntp~.,dformula= ~.,data=simu.obs.df);#family=stats::gaussian,
return(list(dglm.fit, simu.obs.df));
}

out<-simu.dglm.inter();
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


##### how does it work for mutiple categories
nn=10^3
ff<-as.factor(rbinom(nn,2,0.5));
gg<-as.factor(rbinom(nn,4,0.5));
yy<-7+(ff==2)*1.8+(gg==3)*3.6+(ff==0)*(gg==4)*2.5+rnorm(nn,0.02);
aa<-(lm(yy~ff*gg)%>%summary)$coefficients[,c(1,4)]%>%round(2);
lm(yy~ff)%>%summary;

ff.mat<-matrix(0,nn,2)
ff.mat[,1]<-(ff==1); ff.mat[,2]<-(ff==2); 
ffdat<-data.frame(yy,ff.mat)
lm(yy~.,data=ffdat)%>%summary



#### parameters #######
#### parameters #######

#### should we account for missing value?  
#### should there be more possible value variaties, 
#### should there be different prob for different variables?
### main function part #####
### main function part #####
library(dplyr)
simu.dglm.inter<-function(n.rep=3, n.level.env=6, n.obs.per.rep=10^3, 
                          n.loci=8, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), 
                          which.var.loci=c(4:5),  hypo.var.para=c(2,5), 
                          which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
                          loci.prob=c(0.5),  simu.prob.var=rep(0.5, n.loci), 
                          ran.seed=12345)
{
  # require(dglm);
  # n.loci, number of gene loci
  # n.rep, the repetitation of each gene combination. The total numner of observation is n.rep*n.level.env*n.obs.per.rep
  # develop the adjusted p-value
  set.seed(ran.seed);
  
  n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
  n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)
  
  if (which.env.sig>n.level.env) {print("the which.env.sig should be less than n.level.env"); stop;}
  env<-rep(rep(c(1:n.level.env),rep(n.rep,n.level.env)),n.obs.per.rep);
  env.mat<-matrix(0,nrow=n.rep*n.level.env*n.obs.per.rep, ncol=n.level.env);
  for(i in 1:n.level.env)  env.mat[,i]=(env==i);
  var.names.env<-apply(expand.grid('env.var', c(1:n.level.env)), 1, paste, collapse=".");
  colnames(env.mat)<-var.names.env;
  
  simu.loci<-NULL;
  for(i in 1:n.loci) simu.loci<-cbind(simu.loci,rep(rbinom(n.obs.per.rep, 1,simu.prob.var[i]),rep(n.rep*n.level.env,n.obs.per.rep)))
  
  var.names.loci<-apply(expand.grid('loci.var', c(1:n.loci)), 1, paste, collapse=".");
  colnames(simu.loci)<-var.names.loci;
  
  all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
  all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;
  all.inter.para<-rep(0,n.level.env*n.loci);all.inter.para[which.env.inter+(which.loci.inter-1)*n.loci]<-hypo.inter.para;
  
  
  mean.pyntp<-simu.loci%*%all.mean.para;
  var.pyntp<-sapply(simu.loci%*%all.var.para,FUN=function(x) rnorm(1,sd=sqrt(exp(x))));
  
  
  
  
  simu.obs<-cbind(env.mat[,-1],simu.loci); ## the -1 is to get rid of the frist level for collinearity.
  
  ### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 🧬🧬🧬
  ### or, would a gene 🧬🧬🧬 structure(a composition of mutiple A' 🅰️ s and B' 🅱️s from gene🧬🧬🧬 loci) be an important thing to consider, or just with gene loci is contributing to the pynotype?
  inter.mat<-NULL;
  for(i in n.level.env:1) for(j in n.loci:1) inter.mat<-cbind(simu.loci[,j]*env.mat[,i],inter.mat); ## stop at 2 is to get rid of the frist level for collinearity.
  
  var.names.inter.mat<-apply(expand.grid(var.names.loci, var.names.env)%>%as.matrix,1,paste,collapse="*");
  colnames(inter.mat)<-var.names.inter.mat;
  
  all.inter.para<-rep(0,n.level.env*n.loci);
  all.inter.para[which.loci.inter+(which.env.inter-1)*n.loci]<-hypo.inter.para;
  
  #which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
  inter.pyntp<-inter.mat%*%all.inter.para; ### this is wrong
  
  simu.obs.df<-data.frame(simu.obs,inter.mat[,-c(1:n.loci)]);
  simu.obs.df$pyntp<-mean.pyntp+var.pyntp+inter.pyntp;
  
  dglm.fit<-dglm::dglm(pyntp~.,dformula= ~.,data=simu.obs.df);#family=stats::gaussian,
  return(list(dglm.fit, simu.obs.df));
}

out<-simu.dglm.inter();
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
##### how does it work for mutivariate version #####

logN.multi<-function(para, simu.obs.func,sep)
{
  para.beta<-para[1:sep];
  para.gamma<-para[-c(1:sep)];
  yi.star<-simu.obs.func[[1]];  x.val.multi<-simu.obs.func[[2]]; z.val.multi<-simu.obs.func[[3]]; delt.val<-simu.obs.func[[4]];
  wi<- (yi.star-x.val.multi%*%para.beta)/(exp(z.val.multi%*%para.gamma));
  val.loglike<-sum(delt.val*wi^2/2)+sum(delt.val*(z.val.multi%*%para.gamma))-sum((1-delt.val)*log(1-pnorm(wi)+10^(-9)))
  return(val.loglike);
}
logN.comp.multi<-function(para,simu.obs.func,sep)
{
  para.beta<-para[1:sep];
  para.gamma<-para[-c(1:sep)];
  yi.star<-simu.obs.func[[1]];  x.val.multi<-simu.obs.func[[2]]; z.val.multi<-simu.obs.func[[3]]; delt.val<-simu.obs.func[[4]];
  wi<- (yi.star-x.val.multi%*%para.beta)/(exp(z.val.multi%*%para.gamma));
  val.loglike<-sum(wi^2/2)+sum(z.val.multi%*%para.gamma);
  return(val.loglike);
}

simu.rc.normal.multi<-function(n.obs=10^3, 
                               beta.val.multi=c(30,1,2,rep(0,3)),gam.val.multi=c(-5,3,4,rep(0,4)), 
                               cut.off=40)
{
  nrow.x<-length(beta.val.multi)-1;
  nrow.z<-length(gam.val.multi)-1;
  x.val.multi<-cbind(rep(1,n.obs),matrix(rnorm(n.obs*nrow.x,sd=3),nrow=n.obs));
  z.val.multi<-cbind(rep(1,n.obs),matrix(runif(n.obs*nrow.z),nrow=n.obs));
  
  sigma.val<-exp(z.val.multi%*%gam.val.multi);
  epsi<-unlist(lapply(c(1:n.obs), function(x) rnorm(1,sd=sigma.val[x])));## Make it proportional to sigma.val. Noticing that the var od espi is not sigma.val^2
  
  yi<-x.val.multi%*%beta.val.multi+epsi;
  
  yi.star<-unlist(lapply(yi,function(x) min(x,cut.off)));
  delt.val<-(yi<cut.off);
  simu.obs<-list(yi.star,x.val.multi,z.val.multi,delt.val,yi);
  names(simu.obs)<-c('yi.star','X','Z','delta','yi');
  
  rc.mle<-nlminb(start=rep(0,nrow.x+nrow.z+2),logN.multi,simu.obs.func=simu.obs,sep=nrow.x+1,control=list(eval.max=10^6,iter.max=10^6))$par;
  #rc.mle<-optim(par=rep(0,nrow.x+nrow.z+2),method="BFGS",logN.multi, simu.obs.func=simu.obs,sep=nrow.x+1,control=list(maxit=10^6));
  com.mle<-nlminb(start=rep(0,nrow.x+nrow.z+2),logN.comp.multi,simu.obs.func=simu.obs,sep=nrow.x+1,control=list(eval.max=10^6,iter.max=10^6))$par;
  #com.mle<-optim(par=rep(0,nrow.x+nrow.z+2),method="BFGS",logN.comp.multi, simu.obs.func=simu.obs,sep=nrow.x+1,control=list(maxit=10^6));
  
  return(list(rc.mle,com.mle,simu.obs))
}

out.normal.multi<-simu.rc.normal.multi(n.obs=50);
out.normal.multi[[1]];out.normal.multi[[2]];
str(out.normal.multi[[3]])





