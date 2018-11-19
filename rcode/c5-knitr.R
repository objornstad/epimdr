### R code from vignette source 'c5-knitr.Rnw'


###################################################
### code chunk number 3: c5-s1-1
###################################################
ppp = function(wk, x){
require(plotrix)
x = x[wk<53]
wk = wk[wk<53]
ses = sapply(split(x, wk), mean, na.rm=TRUE)
sesdv = sapply(split(x, wk), sd, na.rm =
     TRUE)/sqrt(length(split(x, wk)))
plotCI(x = c(1:52), y = ses, ui = ses+sesdv, 
     li = ses-sesdv, xlab = "Week", ylab = "Incidence")
}


###################################################
### code chunk number 4: c5-s1-2
###################################################
par(mfrow=c(2,2)) #A four panel plot
ppp(paili[, "WEEK"], paili[,"PENNSYLVANIA"])
title("ILI mortality (1972-98)")
ppp(palymes[,"WEEK"], palymes[,"PENNSYLVANIA"])
title("Lymes (2006-14)")
ppp(pagiard[,"WEEK"], pagiard[,"PENNSYLVANIA"])
title("Giardia (2006-14)")
ppp(pameasle[,"WEEK"], pameasle[,"PENNSYLVANIA"])
title("Measles (1928-69)")


###################################################
### code chunk number 5: c5-s2-1
###################################################
seirmod = function(t, y, parms){
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]

  mu = parms["mu"]
  N = parms["N"]
  beta = parms["beta"]
  sigma = parms["sigma"]
  gamma = parms["gamma"]

  dS = mu * (N  - S)  - beta * S * I / N
  dE = beta * S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R 
  res=c(dS, dE, dI, dR)
  list(res)
}


###################################################
### code chunk number 6: c5-s2-2
###################################################
require(deSolve)
times  = seq(0, 10, by=1/120)
paras  = c(mu = 1/50, N = 1, beta =  1000, 
     sigma = 365/8, gamma = 365/5)
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)


###################################################
### code chunk number 7: c5-s2-3
###################################################
R0=expression(sigma/(sigma+mu) * beta / (gamma+ mu))
with(as.list(paras), eval(R0))


###################################################
### code chunk number 8: c5-s2-4
###################################################
out = as.data.frame(ode(start, times, seirmod, paras))
par(mfrow = c(1,2))  #Two plots side by side 
plot(times, out$I, ylab = "Prevalence", 
     xlab = "Time", type = "l")
plot(out$S, out$I, ylab = "Prevalence", 
     xlab = "Susceptible", type = "l")


###################################################
### code chunk number 9: c5-s3-1
###################################################
seirmod2 = function(t, y, parms){
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]
  with(as.list(parms),{
     dS = mu * (N  - S)  - beta0 * (1+beta1 * 
        cos(2 * pi * t)) * S * I / N
     dE = beta0 * (1 + beta1 * cos(2*pi * t)) * 
        S * I / N - (mu + sigma) * E
     dI = sigma * E - (mu + gamma) * I
     dR = gamma * I - mu * R
     res = c(dS, dE, dI, dR)
     list(res)
   })
} 


###################################################
### code chunk number 10: c5-s3-2
###################################################
times  = seq(0, 100, by=1/120)
paras  = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, 
     sigma = 365/8, gamma = 365/5)
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
out = as.data.frame(ode(start, times, seirmod2, paras))
par(mfrow = c(1,2)) #Side-by-side plot
plot(times, out$I, ylab="Infected", xlab="Time", 
     xlim = c(90, 100), ylim = c(0, 
     max(out$I[11001:12000])), type = "l")
plot(out$S[11001:12000], out$I[11001:12000], 
     ylab = "Infected", xlab = "Susceptible", type = "l")


###################################################
### code chunk number 11: c5-s3-3
###################################################
mu=paras["mu"]
N=paras["N"]
beta0=paras["beta0"]
beta1=paras["beta1"]
sigma=paras["sigma"]
gamma=paras["gamma"]
R0= sigma/(sigma+mu)*beta0/(gamma+mu)
#Equilibria
Sstar=1/R0
Istar=mu*(1-1/R0)*R0/beta0
Estar=(mu+gamma)*Istar/sigma
eq=list(S=Sstar, E=Estar, I=Istar)


###################################################
### code chunk number 12: c5-s3-4
###################################################
dS = expression(mu * (N  - S)  - beta0 * S * I / N)
dE= expression(beta0 * S * I / N - (mu + sigma) * E)
dI = expression(sigma*E - (mu + gamma) * I)
#Elements of Jacobian
j11 = D(dS, "S");  j12 = D(dS, "E");  j13 = D(dS, "I")
j21 = D(dE, "S");  j22 = D(dE, "E");  j23 = D(dE, "I")
j31 = D(dI, "S");  j32 = D(dI, "E"); j33 = D(dI, "I")
#Jacobian
J = with(eq,
matrix(c(eval(j11),eval(j12),eval(j13),
   eval(j21),eval(j22), eval(j23),
   eval(j31),eval(j32), eval(j33)), 
    nrow=3, byrow=TRUE))


###################################################
### code chunk number 13: c5-s3-5
###################################################
round(eigen(J)$values, 3)
2*pi/(Im(eigen(J)$values[2]))


###################################################
### code chunk number 14: c5-s4-1 (eval = FALSE)
###################################################
## times = seq(0, 100, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## beta1 = seq(0,0.25, length=101)
## #Matrix to store infecteds
## Imat = matrix(NA, ncol = 12001, nrow = 101)
## #Loop over beta1's
## for(i in 1:101){
##      paras  = c(mu = 1/50, N = 1, beta0 = 1000, 
##         beta1 = beta1[i], sigma = 365/8, gamma = 365/5)
##      out = as.data.frame(ode(start, times, 
##         seirmod2, paras))
##      Imat[i,] = out$I
## }


###################################################
### code chunk number 15: c5-s4-2 (eval = FALSE)
###################################################
## sel = seq(7001, 12000, by = 120)  
## plot(NA, xlim = range(beta1), ylim = c(1E-7, 
##      max(Imat[,sel])), log = "y", xlab = "beta1",
##      ylab = "prevalence")
## for(i in 1:101){
##        points(rep(beta1[i], length(sel)), 
##           Imat[i, sel], pch=20)
## }


###################################################
### code chunk number 16: c5-s4-3 (eval = FALSE)
###################################################
## times = seq(0, 1000, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## beta1 = seq(0,.25, length = 501)
## sel = seq(9001, 120001, by = 120)  
## Imat=matrix(NA, ncol = length(sel), nrow = 501)
## for(i in 1:501){
##        paras  = c(mu = 1/50, N = 1, beta0 = 1000, beta1=beta1[i], 
##        		sigma = 365/8, gamma = 365/5)
##        out = as.data.frame(ode(start, times, seirmod2, paras))
##        Imat[i,] = out$I[sel]
## cat(i, "\r")
## }
## #save(times, start, beta1, sel, Imat, file="bif1.Rd")
## load("bif1.Rd")
## par(mfrow=c(1,1))
## plot(NA, xlim=range(beta1), ylim=c(1E-7,max(Imat)), log="y",
##   xlab=expression(beta[1]), ylab="Prevalence")
## for(i in 1:501){
##        points(rep(beta1[i], dim(Imat)[2]), Imat[i, ], pch=20, cex=0.25)
## }


###################################################
### code chunk number 17: c5-knitr.Rnw:304-307
###################################################
#save(times, start, paras, out, sel, file="c4poincr.Rd")
#save(times, start, paras, out, sel, file="c4RWpoincr.Rd")
load("c4RWpoincr.Rd")


###################################################
### code chunk number 18: c5-s5-1 (eval = FALSE)
###################################################
## times = seq(0, 100000, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## paras  = c(mu = 1/50, N = 1, beta0 = 1800, beta1=0.28, 
##      sigma = 35.84, gamma = 100)
## out = as.data.frame(ode(start, times,seirmod2, paras))
## sel = seq(7001, 12000000, by = 120)  
## par(mfrow = c(1,2))
## plot(out$time[7001:13001], out$I[7001:13001], 
##      type = "l", xlab = "Year", ylab = "Prevalence")
## plot(out$S[sel], out$I[sel], type = "p", xlab = "S", 
##      ylab = "I", log = "y", pch = 20, cex = 0.25)


###################################################
### code chunk number 19: c5-s6-1 (eval = FALSE)
###################################################
## times = seq(0, 100, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## mu = seq(from = 0.005, to = 0.02, length = 101)
## ImatF = ImatB = matrix(NA, ncol = 12001, nrow = 101)
## for(i in 1:101){
##      paras  = c(mu = mu[i], N = 1, beta0 = 2500, 
##           beta1=0.12, sigma = 365/8, gamma = 365/5)
##       out = as.data.frame(ode(start, times, seirmod2, 
##            paras))
##       ImatF[i,] = out$I
##       start = c(S = out$S[12001], E = out$E[12001], 
##            I = out$I[12001], R = out$R[12001])
## }
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## for(i in 101:1){
##      paras  = c(mu = mu[i], N = 1, beta0 = 2500, 
##           beta1 = 0.12, sigma = 365/8, gamma = 365/5)
##      out = as.data.frame(ode(start, times, seirmod2, 
##           paras))
##      ImatB[i,] = out$I
##      start = c(S = out$S[12001], E = out$E[12001], 
##           I = out$I[12001], R = out$R[12001])
## }
## sel = seq(7001, 12000, by = 120)
## par(mfrow = c(1,1))
##      plot(NA, xlim = range(mu),  ylim = range(ImatF[,sel]), 
##           log = "y", xlab = "mu", ylab = "prevalence")
## for(i in 1:101){
##      points(rep(mu[i], dim(ImatF)[2]), ImatF[i, ], 
##           pch = 20, cex = 0.25)
##      points(rep(mu[i], dim(ImatB)[2]), ImatB[i, ], 
##           pch = 20, cex = 0.25, col = 2)
##  }


###################################################
### code chunk number 20: c5-knitr.Rnw:370-407 (eval = FALSE)
###################################################
## times = seq(0, 1000, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## mu=seq(from=0.005,to=0.025, length=501)
## sel=seq(9001, 120001, by=120)  
## ImatF=matrix(NA, ncol = length(sel), nrow = 501)
## for(i in 1:501){
##         paras  = c(mu = mu[i], N = 1, beta0 = 2500, beta1=0.12, 
##         		sigma = 365/8, gamma = 365/5)
##         out = as.data.frame(ode(start, times, seirmod2, paras))
##        ImatF[i,]=out$I[sel]
##        start = c(S = out$S[12001], E = out$E[12001], I = out$I[12001], R=out$R[12001])
## cat(i, "\n")
## }
## times = seq(0, 1000, by = 1/120)
## start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
## mu=seq(from=0.005,to=0.02, length=501)
## sel=seq(9001, 120001, by=120)  
## ImatB=matrix(NA, ncol = length(sel), nrow = 501)
## for(i in 501:1){
##         paras  = c(mu = mu[i], N = 1, beta0 = 2500, beta1=0.12, 
##         		sigma = 365/8, gamma = 365/5)
##         out = as.data.frame(ode(start, times, seirmod2, paras))
##        ImatB[i,]=out$I[sel]
##        start = c(S = out$S[12001], E = out$E[12001], I = out$I[12001], R=out$R[12001])
## cat(i, "\n")
## }
## #save(times, start, mu, sel, ImatB, ImatF, file="bif2FB.Rd")
## 
## load("bif2FB.Rd")
## par(mfrow=c(1,1))
## plot(NA, xlim=range(mu), ylim=c(1E-7,max(ImatF)), log="y",
##   xlab=expression(mu), ylab="Prevalence")
## for(i in 1:501){
##        points(rep(mu[i], dim(ImatF)[2]), ImatF[i, ], pch=20, cex=0.25)
##        points(rep(mu[i], dim(ImatB)[2]), ImatB[i, ], pch=20, cex=0.25,col=2)
## 
## }


###################################################
### code chunk number 21: c5-knitr.Rnw:420-422 (eval = FALSE)
###################################################
## require(shiny)
## SEIR.app


