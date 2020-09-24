library(dplyr)
simu.rc.normal.multi.v2<-function(n.rep=3, n.level.env=2, n.obs.per.rep=150, 
                                  n.loci=3235, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), 
                                  mu.mean=35, mu.var=-7,
                                  which.var.loci=c(4:5),  hypo.var.para=c(2,5), 
                                  which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
                                  loci.prob=c(0.5),  simu.prob.var=rep(0.5, n.loci), 
                                  ran.seed=NULL,cut.off=40)
{
  #nrow.x<-length(beta.val.multi)-1;  nrow.z<-length(gam.val.multi)-1;
  #x.val.multi<-cbind(rep(1,n.obs),matrix(rnorm(n.obs*nrow.x,sd=3),nrow=n.obs));
  #z.val.multi<-cbind(rep(1,n.obs),matrix(runif(n.obs*nrow.z),nrow=n.obs));
  

  ########
  if(is.null(ran.seed)!=TRUE) set.seed(ran.seed);
  
  n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
  n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)
  
  if (which.env.inter>n.level.env) {print("the which.env.sig should be less than n.level.env"); stop;}
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
  
  
  mean.pyntp<-mu.mean+simu.loci%*%all.mean.para;
  sigma.val<-exp(mu.var+simu.loci%*%all.var.para);
  epsi<-unlist(lapply(c(1:length(mean.pyntp)), function(x) rnorm(1,sd=sigma.val[x])));## Make it proportional to sigma.val. Noticing that the var od espi is not sigma.val^2
  
  simu.obs<-cbind(env.mat[,-1],simu.loci); ## the -1 is to get rid of the frist level for collinearity.
  
  inter.mat<-NULL;
  for(i in n.level.env:1) for(j in n.loci:1) inter.mat<-cbind(simu.loci[,j]*env.mat[,i],inter.mat); ## stop at 2 is to get rid of the frist level for collinearity.
  
  var.names.inter.mat<-apply(expand.grid(var.names.loci, var.names.env)%>%as.matrix,1,paste,collapse="*");
  colnames(inter.mat)<-var.names.inter.mat;
  
  all.inter.para<-rep(0,n.level.env*n.loci);
  all.inter.para[which.loci.inter+(which.env.inter-1)*n.loci]<-hypo.inter.para;
  
  #which.env.inter=3,which.loci.inter=7, hypo.inter.para=2.5,
  inter.pyntp<-inter.mat%*%all.inter.para; ### this is wrong
  
  simu.obs.inter<-inter.mat[,-c(1:n.loci)];
  pyntp<-mean.pyntp+epsi+inter.pyntp;
  pyntp.star<-unlist(lapply(pyntp,function(x) min(x,cut.off)));
  delt.val<-(pyntp<cut.off);
  
  #######
  simu.obs.list<-list(pyntp.star, simu.loci, simu.obs.inter,env.mat[,-1],delt.val,pyntp);
  names(simu.obs.list)<-c('yi.star','simu.loci','inter.mat',    'env.mat',   'delta' ,'yi');
  
  #rc.mle<-nlminb(start=rep(0,nrow.x+nrow.z+2),logN.multi,simu.obs.func=simu.obs,sep=nrow.x+1,control=list(eval.max=10^6,iter.max=10^6))$par;
  #rc.mle<-optim(par=rep(0,nrow.x+nrow.z+2),method="BFGS",logN.multi, simu.obs.func=simu.obs,sep=nrow.x+1,control=list(maxit=10^6));
  #com.mle<-nlminb(start=rep(0,nrow.x+nrow.z+2),logN.comp.multi,simu.obs.func=simu.obs,sep=nrow.x+1,control=list(eval.max=10^6,iter.max=10^6))$par;
  #com.mle<-optim(par=rep(0,nrow.x+nrow.z+2),method="BFGS",logN.comp.multi, simu.obs.func=simu.obs,sep=nrow.x+1,control=list(maxit=10^6));
  
  return(simu.obs.list)
}

out <- simu.rc.normal.multi.v2()
df = cbind.data.frame(out[[1]],out[[2]])

df=as.data.frame(df)
df.breed = df[,1]
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
write.csv(df.breed, "simdata3.csv", row.names = FALSE)
library(gridExtra)
library(qtl)
library(vqtl)
test_full <- read.cross(file = "simdata3.csv", format = "csv")

test_full <- drop.nullmarkers(test_full)
test_full <- calc.genoprob(test_full)

test_full <- calc.genoprob(test_full, error.prob = .001)
#hy_p1 <- scanone(cross = test_full, pheno.col = 'stress')
simdata <- scanonevar(cross = test_full, 
                      mean.formula = stress ~ BreedType*mean.QTL.add, 
                      var.formula = ~ BreedType*var.QTL.add, 
                      return.covar.effects = TRUE)


