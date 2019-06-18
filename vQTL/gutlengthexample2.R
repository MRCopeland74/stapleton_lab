out.i <- scanone(gutlength, addcovar=x, intcovar=sex)
plot(out.i, out.i - out.a, ylab="LOD score", col=c("blue", "red"), alternate.chrid=TRUE)
set.seed(54955149)
operm.a <- scanone(gutlength, addcovar=x, n.perm=1000, perm.Xsp=TRUE, perm.strata=strat)
set.seed(54955149)
operm.i <- scanone(gutlength, addcovar=x, intcovar=sex, n.perm=1000, perm.Xsp=TRUE, perm.strata=strat)                    
out.ia <- c(out.i, out.i - out.a, labels=c("f","i")) 
operm.ia <- cbind(operm.i, operm.i - operm.a, labels=c("f","i"))                    
summary(out.ia, perms=operm.ia, alpha=0.2, pvalues=TRUE)                    
summary(out.ia, perms=operm.ia, alpha=0.2, pvalues=TRUE, + lodcolumn=2)                    
                 
                                         