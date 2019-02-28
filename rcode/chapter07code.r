### R code from Chapter 7
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866


###################################################
### code chunk number 3: c7-s1-1
###################################################
SimTsir = function(alpha = 0.97, B = 2300, beta = 25, 
    sdbeta = 0, S0 = 0.06, I0 = 180, IT = 520, 
    N = 3.3E6){
    #Set up simulation  
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    #Add initial conditions
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    
    #Run  simulation
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean = beta, sd = sdbeta) * 
           I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i] = 0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
     #Return result
    list(I = I, S = S)
}


###################################################
### code chunk number 4: c7-s1-2
###################################################
out = SimTsir()
par(mfrow = c(1, 2)) 
plot(out$I, ylab = "Infected", xlab = "Time", type = "l")
plot(out$S, out$I, ylab = "Infected", xlab = "Susceptible", 
     type = "l")


###################################################
### code chunk number 5: c7-s2-1
###################################################
data(meas)

###################################################
### code chunk number 7: c7-s4-1
###################################################
cum.reg = smooth.spline(cumsum(meas$B), 
   cumsum(meas$London), df = 5)
D = - resid(cum.reg) #The residuals

plot(cumsum(meas$B), cumsum(meas$London), type = "l", 
   xlab = "Cumulative births", ylab = "Cumuative incidence")
lines(cum.reg)
abline(a = 0, b = 1)


###################################################
### code chunk number 8: c7-s4-2
###################################################
rr = predict(cum.reg, deriv = 1)$y
summary(rr)


###################################################
### code chunk number 9: c7-s4-3
###################################################
Ic = meas$London/rr
Dc=D/rr 


###################################################
### code chunk number 10: c7-s4-4
###################################################
seas = rep(1:26, 21)[1:545]
lInew = log(Ic[2:546])
lIold = log(Ic[1:545])
Dold = Dc[1:545]


###################################################
### code chunk number 11: c7-s4-5
###################################################
N=3.3E6
Smean = seq(0.02, 0.2, by=0.001)*N
offsetN=rep(-log(N), 545)


###################################################
### code chunk number 12: c7-s4-6
###################################################
llik = rep(NA, length(Smean))
for(i in 1:length(Smean)){
     lSold = log(Smean[i] + Dold)
     glmfit = glm(lInew ~ -1 +as.factor(seas) + lIold + 
          offset(lSold + offsetN))
     llik[i] = glmfit$deviance / 2
  }
par(mfrow=c(1, 1))
plot(Smean/3.3E6, llik, ylim=c(min(llik), 25), 
     xlab="Sbar", ylab="Neg log-lik")


###################################################
### code chunk number 13: c7-s4-7
###################################################
lSold = log(Smean[which(llik==min(llik))] + Dold)
glmfit = glm(lInew ~ -1 + as.factor(seas) + lIold +
     offset(lSold + offsetN))


###################################################
### code chunk number 14: c7-s4-8
###################################################
Smean[which.min(llik)]/3.3E6


###################################################
### code chunk number 15: c7-s4-9
###################################################
glmfit$coef[27]


###################################################
### code chunk number 16: c7-s4-10
###################################################
require(plotrix)
beta = exp(glmfit$coef[1:26])
ubeta = exp(glmfit$coef[1:26] +
     summary(glmfit)$coef[1:26, 2])
lbeta = exp(glmfit$coef[1:26] - 
     summary(glmfit)$coef[1:26, 2])
plotCI(x = c(1:26), y = beta, ui = ubeta, li = lbeta, 
     xlab = "Biweek", ylab = expression(beta))


###################################################
### code chunk number 17: c7-s5-1
###################################################
SimTsir2 = function(beta, alpha, B, N,  inits = 
    list(Snull = 0, Inull = 0), type = "det"){
    type = charmatch(type, c("det", "stoc"), 
       nomatch = NA)
    if(is.na(type))
       stop("method should be \"det\", \"stoc\"")
    IT = length(B)
    s = length(beta)
    lambda = rep(NA, IT)  
    I = rep(NA, IT)
    S = rep(NA, IT)
    
    I[1] = inits$Inull
    lambda[1] = inits$Inull
    S[1] = inits$Snull
    
    for(i in 2:IT) {
        lambda[i] = beta[((i - 2) %% s) + 1] *  
           S[i - 1] *(I[i - 1]^alpha)/N
        if(type == 2) {
           I[i] = rpois(1, lambda[i])
         }
        if(type == 1) {
           I[i] = lambda[i]
        }
        S[i] =S[i - 1] + B[i] - I[i]
    }
    return(list(I = I, S = S))
}


###################################################
### code chunk number 18: c7-s5-2
###################################################
sim=SimTsir2(beta=exp(glmfit$coef[1:26]), alpha=0.966, 
     B=meas$B, N=N, inits=list(Snull=Dc[1]+
     Smean[which(llik==min(llik))], Inull=Ic[1]))
plot(sim$I, type="b", ylim=c(0, max(Ic)), 
     ylab="Incidence", xlab="Biweek")
lines(exp(lInew), col="red")
legend("topleft", legend=c("sim", "Ic"), lty=c(1,1),
     pch=c(1,NA), col=c("black", "red"))


###################################################
### code chunk number 19: c7-s6-1
###################################################
data(dalziel)


###################################################
### code chunk number 20: c7-s6-2
###################################################
data(tyscarlet)
tyscarlet = tyscarlet[tyscarlet$WEEK < 53,]
tyscarlet = tyscarlet[tyscarlet$YEAR > 1914,]
ag = rep(1:(dim(tyscarlet)[1]/2), each = 2)
scarlet2 = sapply(split(tyscarlet$PHILADELPHIA, ag), sum)


###################################################
### code chunk number 21: c7-s6-3
###################################################
require(imputeTS)
philly = dalziel[dalziel$loc=="PHILADELPHIA", ]
philly = philly[philly$year > 1914 & philly$year < 1948,]
philly$cases = na.interpolation(ts(scarlet2))


###################################################
### code chunk number 22: c7-s6-4
###################################################
cum.reg = smooth.spline(cumsum(philly$rec), 
     cumsum(philly$cases), df = 10)
D = - resid(cum.reg) #The residuals
rr = predict(cum.reg, deriv = 1)$y
summary(rr)


###################################################
### code chunk number 23: c7-s6-5
###################################################
Ic = philly$cases/rr
Dc=D/rr 
seas = rep(1:26, 21)[1:597]
lInew = log(Ic[2:598])
lIold = log(Ic[1:597])
Dold = Dc[1:597]
N = median(philly$pop)
offsetN = rep(-log(N), 597)


###################################################
### code chunk number 24: c7-s6-6
###################################################
Smean = seq(0.02, 0.6, by=0.001)*N
llik = rep(NA, length(Smean))
for(i in 1:length(Smean)){
     lSold = log(Smean[i] + Dold)
     glmfit = glm(lInew ~ -1 + as.factor(seas) + lIold + 
          offset(lSold + offsetN))
     llik[i] = glmfit$deviance
}
Smean[which(llik == min(llik))]/N


###################################################
### code chunk number 25: c7-s6-7
###################################################
lSold = log(Smean[which.min(llik)] + Dold)
glmfit = glm(lInew ~ -1 +as.factor(seas) + lIold +
     offset(lSold + offsetN))
#alpha
glmfit$coef[27]


###################################################
### code chunk number 26: c7-s6-8
###################################################
beta = exp(glmfit$coef[1:26])
ubeta = exp(glmfit$coef[1:26] +
     summary(glmfit)$coef[1:26, 2])
lbeta = exp(glmfit$coef[1:26] -
     summary(glmfit)$coef[1:26, 2])
plotCI(x = c(1:26), y = beta, ui = ubeta, li = lbeta, 
     xlab = "Biweek", ylab = expression(beta))


###################################################
### code chunk number 27: c7-s7-1
###################################################
data(SH9)
#subset RBC data
SH9rbc = SH9[ ,-c(1,3,4,7,8,10,11)]
#Bump up RBC to microliter
SH9rbc[ ,4] = SH9rbc[,4]*10^6
#subset parasitemia data
SH9para = SH9[ ,-c(1,3,4,7,8,9, 10)]

#reshape to wide
SH9rbcw = reshape(SH9rbc, idvar = "Ind2", direction = 
     "wide", timevar="Day")
SH9pw = reshape(SH9para, idvar = "Ind2", direction =
     "wide", timevar="Day")
#delete duplicate columns
SH9rbcw = SH9rbcw[ ,-seq(4, 50,by = 2)]
names(SH9rbcw)[2] = "Treatment"
 SH9pw=SH9pw[ ,-seq(4, 50, by = 2)]
 names(SH9pw)[2] = "Treatment"
 
 #drop last columns of data not counted every day
 SH9pw = SH9pw[ ,-c(22:27)]
 SH9rbcw = SH9rbcw[ ,-c(22:27)]
 #Pull out AQ mice
 paras=SH9pw[1:10,-c(1:2)]
 SH9rbcw = as.matrix(SH9rbcw[1:10, -c(1:2)])
 #Uninfected are total RBCs less infected
 RBCs = as.matrix(SH9rbcw-paras)


###################################################
### code chunk number 28: c7-s7-2
###################################################
par(mfrow = c(1, 2), bty = "l")
matplot(t(log(RBCs)), type = "l",xlab = "Day", 
     ylab = "Uninfected log-RBCs")
matplot(t(log(paras)),type = "l", xlab = "Days", 
     ylab = "Infected log-RBCs")


###################################################
### code chunk number 29: c7-s7-3
###################################################
Tmax = length(paras[1,])  ##max number of days
Nind = length(paras[,1]) ##number of individuals
day = matrix(rep(1:(Tmax-1), each = Nind), Nind, Tmax-1)
day = c(day)

#Log infected cells
log.para = log(paras[,2:Tmax])
log.para.lag = log(paras[,1:(Tmax-1)])
log.para = unlist(c(log.para))
log.para.lag = unlist(c(log.para.lag))
 
#Log uninfected cells
log.rbcs.lag = log(RBCs[,1:(Tmax-1)])
log.rbcs.lag = unlist(c(log.rbcs.lag))


###################################################
### code chunk number 30: c7-s7-4
###################################################
log.para[!is.finite(log.para)] = 
     min(log.para[is.finite(log.para)] , na.rm=T)
log.para.lag[!is.finite(log.para.lag)] = 
     min(log.para[is.finite(log.para)] , na.rm=T)

###################################################
### code chunk number 31: c7-s7-5
###################################################
data = data.frame(log.para = log.para, day = day, 
   log.para.lag = log.para.lag, log.rbcs.lag = log.rbcs.lag)
fit = glm(log.para ~ -1 + as.factor(day) +
   offset(log.para.lag + log.rbcs.lag), data = data)

###################################################
### code chunk number 32: c7-s7-6
###################################################
par(mfrow = c(1, 2))
require(plotrix)
ses = summary(fit)$coeff[ ,2]
beta = exp(fit$coef)
ubeta = exp(fit$coef + ses)
lbeta = exp(fit$coef - ses)
plotCI(x = c(3:20), y = beta, ui = ubeta, li = lbeta, 
     xlab = "Day", ylab = expression(P[E]))
points(x = c(3:20), exp(fit$coeff), type = "b", pch = 19)
plotCI(x = c(3:20), y = beta*colMeans(RBCs)[-19], ui = ubeta*
     colMeans(RBCs)[-19], li = lbeta*colMeans(RBCs)[-19], 
     xlab = "Day", ylab = expression(R[E]))
points(x = c(3:20), beta*colMeans(RBCs)[-19], 
     type = "b",pch = 19)
abline(h = 1,lty = 3)


###################################################
### code chunk number 32: c7-knitr.Rnw:629-631 (eval = FALSE)
###################################################
require(shiny)
TSIR.app


