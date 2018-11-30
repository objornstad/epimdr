### R code from Chapter 11
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866

###################################################
### code chunk number 2: c11-s2-1
###################################################
data(filipendula)
symbols(filipendula$X, filipendula$Y, circles = 
   rep(1,162), inches = .1, bg = filipendula$y95+1, 
   xlab = "X", ylab = "Y")
symbols(filipendula$X, filipendula$Y, circles = 
   rep(1,162), inches = .05, bg = filipendula$y94+1, 
   add = TRUE)
legend("topright", c("infected 94", "infected 95"), 
   pch = c(21,21), pt.cex = c(1,2), pt.bg = c(2,2))


###################################################
### code chunk number 3: c11-s2-2
###################################################
dst = as.matrix(dist(filipendula[,c("X","Y")]))


###################################################
### code chunk number 4: c11-s2-3
###################################################
a = 10
foi = apply(exp(-dst/a)*filipendula$y94,2,sum)


###################################################
### code chunk number 5: c11-s2-4
###################################################
lfit = glm(y95~foi, family = binomial(), data = filipendula)
lfit$deviance/2


###################################################
### code chunk number 6: c11-s2-5
###################################################
a = seq(1,20, length = 1001)
llik = rep(NA, length(a)) 
for(i in 1:length(a)){
   foi = apply(exp(-dst/a[i])*filipendula$y94,2,sum)
   lfit = glm(y95~foi, family = binomial(), 
       data = filipendula)
   llik[i] = lfit$deviance/2
}
plot(a, llik, type = "l", ylab = "Neg. log-like")
abline(h = min(llik)+qchisq(0.95,1)/2)


###################################################
### code chunk number 7: c11-s2-6
###################################################
ahat = a[which.min(llik)]
foi = apply(exp(-dst/ahat)*filipendula$y94,2,sum)
spmod = glm(y95~foi, family = binomial(), data = filipendula)
nullmod = glm(y95~1, family = binomial(), data = filipendula)
#correct the df of the spmod
spmod$df.residual = spmod$df.residual-1
anova(nullmod, spmod, test = "Chisq")


###################################################
### code chunk number 8: c11-s2-7
###################################################
a2 = seq(1,20, length = 1001)
llik2 = rep(NA, length(a2)) 
for(i in 1:length(a2)){
   foi2 = apply(exp(-(dst/a2[i])^2)*filipendula$y94,2,sum)
   lfit2 = glm(y95~foi2, family = binomial(), 
      data = filipendula)
   llik2[i] = lfit2$deviance/2
}
ahat2 = a2[which.min(llik2)]
foi2 = apply(exp(-(dst/ahat2)^2)*filipendula$y94,2,sum)
spmod2 = glm(y95~foi2, family = binomial(), 
   data = filipendula)
spmod2$df.residual = spmod2$df.residual-1


###################################################
### code chunk number 9: c11-s2-8
###################################################
curve((2/(ahat2*gamma(1/2)))*exp(-((x/ahat2)^2)), 0, 10, col=2, lty=2, ylab="Probability", xlab="Meters")
curve((1/(ahat)*gamma(1))*exp(-x/ahat), 0, 10, add=TRUE)
legend("topright", c("Exponential", "Gaussian"), lty=c(1,2), col=c(1,2))


###################################################
### code chunk number 10: c11-s2-9
###################################################
spmod$aic
spmod2$aic


###################################################
### code chunk number 11: c11-s3-1
###################################################
zprev = filipendula$y95
x = filipendula$X
y = filipendula$Y
beta0 = spmod$coef[1]
beta1 = spmod$coef[2]


###################################################
### code chunk number 12: c11-s3-2
###################################################
foi = apply(exp(-dst/ahat)*zprev,2,sum)
logitp = beta0+beta1*foi
p = exp(logitp)/(1+exp(logitp))


###################################################
### code chunk number 13: c11-s3-3
###################################################
znew = rbinom(162, 1, p)
symbols(x, y, circles = rep(1,162), bg = znew+1, inches = .1, xlab = "X", ylab = "Y")


###################################################
### code chunk number 14: c11-s3-4
###################################################
simdat = matrix(NA, ncol = 100, nrow = 162)
for(i in 1:100){
     zprev = znew
     foi = apply(exp(-dst/ahat)*zprev,2,sum)
     logitp = beta0+beta1*foi
     p = exp(logitp)/(1+exp(logitp))
     znew = rbinom(162, 1, p)
     simdat[,i] = znew
     #Code for in-line animation:
     #symbols(x, y, circles = rep(1,162), bg = znew+1, 
     #  inches = .1, xlab = "X", ylab = "Y")
     #Sys.sleep(0.1)
}


###################################################
### code chunk number 15: c11-s3-5
###################################################
require(ncf)
mprev = apply(simdat, 1, mean)
spatial.plot(x, y, mprev, ctr = TRUE)


###################################################
### code chunk number 16: c11-s5-1
###################################################
local.dyn = function(t, S, I, b0, b1, mu, N) {
  beta = b0*(1+b1*cos(2*pi*t/26))
  I = S*(1-exp(-beta*I))
  S = (1-mu)*S+mu*N-I
  list(S = S,I = I)
}


###################################################
### code chunk number 17: c11-s5-2
###################################################
m = 0.25
ny = nx = 30
#generate coordinates
xy = expand.grid(x = 1:nx, y = 1:ny)
#make distance matrix
dst  =  as.matrix(dist(xy))
#make redistribution matrix with zeros
redist = matrix(0, nrow = ny*nx, ncol = ny*nx)
#populate the matrix so each of the 8 neighbors gets their share
redist[dst<1.5] = m/8
#the remaining fraction stays put
diag(redist) = 1-m


###################################################
### code chunk number 18: c11-s5-3
###################################################
IT = 520
S = I = matrix(NA, nrow = ny*nx, ncol = IT)
S[,1] = 100
I[,1] = 0
I[400,1] = 1


###################################################
### code chunk number 19: c11-s5-4
###################################################
b0 = 0.04
b1 = 0.8
mu = 0.02/26
N = 1000


