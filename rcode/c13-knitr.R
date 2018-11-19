### R code from vignette source 'c13-knitr.Rnw'


###################################################
### code chunk number 3: c13-s2-1
###################################################
data(gra)
gra$jx=jitter(gra$xloc)
gra$jy=jitter(gra$yloc)
symbols(y=gra$xloc, x=gra$yloc, circles=gra$score, 
     inches=0.1, xlab="y", ylab="x")
abline(v=47.5,col=2)
abline(v=97.5,col=2)
abline(v=147.5,col=2)


###################################################
### code chunk number 4: c13-knitr.Rnw:50-56 (eval = FALSE)
###################################################
## pdf(file="c12-koslow.pdf",width=8,height=3)
## symbols(y=gra$xloc, x=gra$yloc, circles=gra$score, inches=0.1, xlab="y", ylab="x")
## abline(v=47.5,col=2)
## abline(v=97.5,col=2)
## abline(v=147.5,col=2)
## dev.off()


###################################################
### code chunk number 5: c13-s3-1
###################################################
n=length(gra$score)
#marginal mean:
mu=mean(gra$score)
#marginal MLE sd:
sig=sd(gra$score)*(n-1)/n


###################################################
### code chunk number 6: c13-s3-2
###################################################
#rescale Zs
zscale=(gra$score-mu)/sig
#autocorrelation matrix
rho=outer(zscale, zscale)


###################################################
### code chunk number 7: c13-s3-3
###################################################
dst = as.matrix(dist(gra[,c("xloc","yloc")]))


###################################################
### code chunk number 8: c13-s3-4
###################################################
plot(dst[1:1000], rho[1:1000], ylab="Pairwise rho", 
     xlab="Pairwise distance (m)")


###################################################
### code chunk number 9: c13-s4-1
###################################################
require(ncf)


###################################################
### code chunk number 10: c13-s5-1 (eval = FALSE)
###################################################
## test1=mantel.test(M1=rho, M2=dst)


###################################################
### code chunk number 11: c13-knitr.Rnw:188-189
###################################################
load("c12test1.rd")


###################################################
### code chunk number 12: c13-s5-2
###################################################
test1


###################################################
### code chunk number 13: c13-knitr.Rnw:200-201 (eval = FALSE)
###################################################
## test=mantel.test(x=..., y= ..., z=...)


###################################################
### code chunk number 14: c13-knitr.Rnw:206-207
###################################################
load("c12test2.rd")


###################################################
### code chunk number 15: c13-s6-1 (eval = FALSE)
###################################################
## test2 = correlog(x = gra$xloc, y = gra$yloc, 
##      z = gra$score, increment = 10)
## plot(test2)


###################################################
### code chunk number 16: c13-s7-1 (eval = FALSE)
###################################################
## test3 = spline.correlog(x = gra$xloc, y = gra$yloc, 
##    z = gra$score)


###################################################
### code chunk number 17: c13-knitr.Rnw:230-231
###################################################
load("c12test3.rd")


###################################################
### code chunk number 18: c13-s7-2
###################################################
summary(test3)


###################################################
### code chunk number 19: c13-s7-3
###################################################
plot(test3)


###################################################
### code chunk number 20: c13-s8-1 (eval = FALSE)
###################################################
## test4 = lisa(x = gra$yloc, y = gra$xloc, z = gra$score, 
##     neigh = 20)


###################################################
### code chunk number 21: c13-knitr.Rnw:270-271
###################################################
load("c12test4.rd")


###################################################
### code chunk number 22: c13-s8-2
###################################################
plot(test4)


###################################################
### code chunk number 23: c13-s9-1
###################################################
data(silene2)
symbols(silene2$lon, silene2$lat, circles = 
     sqrt(silene2$dmean), inches = .2, xlab = "Longitude",
     ylab = "Latitude")


###################################################
### code chunk number 24: c13-knitr.Rnw:303-304
###################################################
load("c12testcc.rd")


###################################################
### code chunk number 25: c13-s9-2
###################################################
silene2$ab = silene2$dmean + silene2$hmean
silene2$prev = silene2$dmean/(silene2$dmean + silene2$hmean)


###################################################
### code chunk number 26: c13-s9-3 (eval = FALSE)
###################################################
## testcc=spline.correlog(x=silene2$lon, y=silene2$lat, 
##    z=silene2$prev, w=sqrt(silene2$ab), 
##    latlon=TRUE, na.rm=TRUE)
## plot(testcc)


###################################################
### code chunk number 27: c13-s9-4 (eval = FALSE)
###################################################
## data(filip)
## testcc2=correlog(x=filip$X, y=filip$Y, z=filip$y94, 
##      w=filip$y95, increment=25)


###################################################
### code chunk number 28: c13-knitr.Rnw:329-330
###################################################
load("c12testcc2.rd")


###################################################
### code chunk number 29: c13-s9-5
###################################################
testcc2$corr0
testcc2$x.intercept


###################################################
### code chunk number 30: c13-s10-1 (eval = FALSE)
###################################################
## data(gm)
## sel = apply(gm[3:30],1,sum)!=0
## #Synchrony:
## fit1 = Sncf(gm[sel,1]/1000, gm[sel,2]/1000, 
##      gm[sel,3:30], resamp = 500)
## #Lag 1 cross-correlation
## fit2 = Sncf(gm[sel,1]/1000, gm[sel,2]/1000, 
##      z = gm[sel,3:29], w = gm[sel,4:30], resamp = 500)
## #Lag 2 cross-correlation
## fit3 = Sncf(gm[sel,1]/1000, gm[sel,2]/1000, 
##      z = gm[sel,3:28], w = gm[sel,5:30], resamp = 500)


###################################################
### code chunk number 31: c13-knitr.Rnw:359-360
###################################################
load("gmsncf.rd")


###################################################
### code chunk number 32: c13-s10-2
###################################################
par(mfrow=c(1,3))
plot(fit1, ylim=c(-.1, 1))
plot(fit2, ylim=c(-.1, 1))
title("Lag 1")
plot(fit3, ylim=c(-.1, 1))
title("Lag 2")


