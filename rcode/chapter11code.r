### R code from Chapter 11
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866

###################################################
### code chunk number 1: c11-knitr.rnw:1-4
###################################################
library(knitr)
options(width=48)
opts_chunk$set(tidy=TRUE, fig.keep='none')


###################################################
### code chunk number 2: c11-s2-1
###################################################
data(filipendula)
symbols(filipendula$X, filipendula$Y, circles = 
   rep(1,162), inches = .3, bg = filipendula$y95+1, 
   xlab = "X", ylab = "Y")
symbols(filipendula$X, filipendula$Y, circles = 
   rep(1,162), inches = .15, bg = filipendula$y94+1, 
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
     symbols(x, y, circles = rep(1,162), bg = znew+1, 
      inches = .1, xlab = "X", ylab = "Y")
     Sys.sleep(0.1)
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


###################################################
### code chunk number 20: c11-s5-5
###################################################
for (t in 2:IT) {
   #local growth:
   tmp = local.dyn(t,S=S[,t-1],I=I[,t-1],b0=b0,b1=b1,mu=mu,N=N)
   #spatial movement
   S[,t] = redist%*%tmp$S
   I[,t] = redist%*%tmp$I
  #progress monitor
   cat(t," of ",IT,"\r")
 }


###################################################
### code chunk number 21: c11-s5-6
###################################################
x = xy[,1]
y = xy[,2]
scIcubed = I^(1/4)/(max(I[,10:IT]^(1/4)))


###################################################
### code chunk number 22: c11-s5-7 (eval = FALSE)
###################################################
for (k in 1:IT) {
symbols(x,y,fg=2,circles=scIcubed[,k],inches=FALSE,bg=2,xlab="",ylab="")
Sys.sleep(.05)}


###################################################
### code chunk number 23: c11-s6-1 (eval = FALSE)
###################################################
for(k in 100:IT){
png(filename = paste("m",1000+k,".jpg", sep = ""))
     symbols(x, y, fg = 2, circles = scIcubed[,k], 
     inches = FALSE, bg = 2,xlab = "",ylab = "")
dev.off()
}
system("convert m*.jpg simovie.gif")
system("rm m*.png")
## #For mp4-animation:
## #system("convert -delay 5 m*.jpg simovie.mp4")


###################################################
### code chunk number 24: c11-s6-2 (eval = FALSE)
###################################################
require("animation")
oopt = ani.options(interval = 0.02, nmax = 100)
test.function = function (xy, I, nmax) {
     x = xy[,1]
     y = xy[,2]
     scIcubed = I^(1/4)/(max(I[,10:IT]^(1/4)))
     for (i in seq_len(ani.options("nmax"))) {
         dev.hold()
         symbols(x,y,fg = 2,circles = I[,i],inches = 0.1,bg = 2,
           xlab = "",ylab = "")
         ani.pause()
     }
}

saveLatex({
    test.function(xy=xy, I=I, nmax=50)
    }, 
    ani.basename = "BM", ani.opts = "controls,loop,
    width=0.8\\textwidth", ani.first = 
     par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0),
     tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, 
     cex.main = 1), latex.filename = "test.tex", 
     pdflatex = "/usr/texbin/pdflatex", 
     img.name = "Xplot")
ani.options(oopt)


###################################################
### code chunk number 25: c11-s7-1 (eval = FALSE)
###################################################
fitI = Sncf(x=xy[,1], y=xy[,2], z=sqrt(I[,261:520]), resamp=500)
fitS = Sncf(x=xy[,1], y=xy[,2], z=sqrt(S[,261:520]), resamp=500)
fitSI = Sncf(x=xy[,1], y=xy[,2],  z=sqrt(S[,261:520]), w=sqrt(I[,261:520]), resamp=500)
par(mfrow=c(1,3))
plot(fitI, ylim=c(-.1, 1))
plot(fitS, ylim=c(-.1, 1))
plot(fitSI, ylim=c(-.2, .2))

###################################################
### code chunk number 26: c11-s7-2 (eval = FALSE)
###################################################
fitIlag  =  Sncf(x = xy[,1], y = xy[,2], z = I[,261:515], w = I[,266:520], resamp = 100)
plot(fitIlag, ylim = c(-.2,.2))


###################################################
### code chunk number 26: c11-s7-2 (eval = FALSE)
###################################################
require(deSolve)
SIR.space = function(t, y, pars){
  i = c(1:L)
  S = y[i]
  I = y[L+i]
  R = y[2*L+i]
  with(pars,{
  beta = beta[i]
  dS = -(beta*I + m*G%*%I)*S
  dI = (beta*I + m*G%*%I)*S - gamma*I
  dR = gamma*I 
  list(c(dS, dI, dR)) 
})
}


###################################################
### code chunk number 27: c11-s8-2
###################################################
require(ncf)
data(usflu)
usdist = gcdist(usflu$Longitude, usflu$Latitude)


###################################################
### code chunk number 28: c11-s8-3
###################################################
gravity = function(tau1, tau2, phi, pop, distance){
	gravity = outer(pop^tau1, pop^tau2)/distance^phi
	diag(gravity) = 0
	gravity}
G = gravity(0.3, 0.6, 3, usflu$Pop, usdist)


###################################################
### code chunk number 29: c11-s8-4
###################################################
gamma = 1/3.5
R0=1.8
beta = R0*gamma/usflu$Pop
m = 1/1000/sum(usflu$Pop)
parms = list(beta = beta, m = m, gamma  =  gamma, G = G) 
L = length(usflu$Pop)

S = usflu$Pop
R = I = rep(0, length(usflu$Pop))
I[31] = 1 
inits = c(S = S, I = I, R = R)	


###################################################
### code chunk number 30: c11-s8-5
###################################################
require(deSolve)
times = 0:200
out = ode(inits, times, SIR.space, parms)
matplot(out[,50+(1:L)], type = "l", ylab = "Prevalence", 
     xlab = "Day")


