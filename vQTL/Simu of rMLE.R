logL<-function(para, simu.obs.func)
{
  ## para[1] is the intercept of the mean model, 2 is beta, 3 os musigma, 4 is gamma
  attach(simu.obs.func);
  wi<- (yi.star-para[1]-x.val*para[2])/(exp(para[3]+z.val*para[4]));
  val.loglike<-sum(wi)+sum(delt.val*(para[3]+z.val*para[4]))+sum((1+delt.val)*log(1+exp(-wi)))
  detach(simu.obs.func);
  return(val.loglike);
}


logL.comp<-function(para, simu.obs.func)
{
  ## para[1] is the intercept of the mean model, 2 is beta, 3 os musigma, 4 is gamma
  attach(simu.obs.func);
  wi<- (yi.star-para[1]-x.val*para[2])/(exp(para[3]+z.val*para[4]));
  val.loglike<-sum(wi)+sum((para[3]+z.val*para[4]))+2*sum(log(1+exp(-wi)))
  detach(simu.obs.func);
  return(val.loglike);
}

simu.rc.logis<-function(n.obs=10^3, 
                        beta.val=2,gam.val=3, alp.val=30, mu.sig=-1,
                        cut.off=40)
{
    x.val<-rnorm(n.obs,sd=3);
    z.val<-runif(n.obs);
    
    sigma.val<-exp(mu.sig+z.val*gam.val);
    epsi<-unlist(lapply(c(1:n.obs), function(x) rlogis(1,scale=sigma.val[x])));## Make it proportional to sigma.val. Noticing that the var od espi is not sigma.val^2

    yi<-alp.val+beta.val*x.val+epsi;
    
    yi.star<-unlist(lapply(yi,function(x) min(x,cut.off)));
    delt.val<-(yi<cut.off);
    simu.obs<-data.frame(cbind(yi.star,x.val,z.val,delt.val));
    
    rc.mle<-optim(par=rep(0,4),logL, simu.obs.func=simu.obs)$par;
    com.mle<-optim(par=rep(0,4),logL.comp, simu.obs.func=simu.obs)$par;
    
    return(c(rc.mle,com.mle))
}

simu.rc.logis();








