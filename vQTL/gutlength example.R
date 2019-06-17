library(qtl)
library(qtlbook)
data(gutlength)

gutlength <- subset(gutlength, chr = -15)

boxplot(gutlength ~ sex*cross, data=gutlength$pheno, 
        horizontal=TRUE, xlab="Gut length (cm)",
        col=c("red","blue"))

anova(aov(gutlength ~ sex*cross, data=gutlength$pheno))

cross <- as.numeric(pull.pheno(gutlength, "cross")) 
frommom <- as.numeric(cross < 3)
forw <- as.numeric(cross == 1 | cross == 3)
gutlength$pheno$frommom <- frommom
gutlength$pheno$forw <- forw
anova(aov(gutlength ~ sex*frommom*forw, data=gutlength$pheno))

sex <- as.numeric(pull.pheno(gutlength, "sex") == "M")
crossX <- cbind(frommom, forw, frommom*forw)
x <- cbind(sex, crossX)

gutlength <- calc.genoprob(gutlength, step=1,  error.prob=0.001) 
out.0 <- scanone(gutlength)
out.a <- scanone(gutlength, addcovar=x)
plot(out.0, out.a, col=c("blue", "red"), lty=1:2,  ylab="LOD score", alternate.chrid=TRUE)

plot(out.a - out.0, ylab="LOD w/ covar - LOD w/o covar",  ylim=c(-1, 1), alternate.chrid=TRUE)
abline(h=0, lty=2)

strat <- (nmissing(gutlength) < 50)
operm.0 <- scanone(gutlength, n.perm=1000, perm.Xsp=TRUE, perm.strata=strat)
operm.a <- scanone(gutlength, addcovar=x, n.perm=1000, perm.Xsp=TRUE, perm.strata=strat) 

out.both <- c(out.0, out.a, labels=c("nocovar", "covar")) 
operm.both <- cbind(operm.0, operm.a, labels=c("nocovar", "covar"))
summary(operm.both, 0.05)                                                                               
summary(out.both, perms=operm.both, format="allpeaks", alpha=0.2, pvalues=TRUE)
