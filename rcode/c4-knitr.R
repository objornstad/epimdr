### R code from vignette source 'c4-knitr.Rnw'



###################################################
### code chunk number 3: c4-knitr.Rnw:59-61 (eval = FALSE)
###################################################
## glm(cbind(inf, notinf) ~ offset(log(a)), 
##      family = binomial(link = "cloglog")) 


###################################################
### code chunk number 4: c4-s3-1
###################################################
data(black)
black


###################################################
### code chunk number 5: c4-s3-2
###################################################
b2 = black[-c(1,8,9),]  #subsetting age brackets
#Estimate log-FoI
fit = glm(cbind(pos,neg) ~ offset(log(mid)), 
    family = binomial(link = "cloglog"), data = b2)
#Plot predicted and observed
phi = exp(coef(fit))
curve(1-exp(-phi*x), from = 0, to = 60, 
     ylab = 'Seroprevalence', xlab = 'Age')
points(black$mid, black$f, pch = '*', col = 'red')
points(x = b2$mid, y = b2$f, pch = 8)
exp(fit$coef)


###################################################
### code chunk number 6: c4-s4-1
###################################################
data(rabbit)
head(rabbit)


###################################################
### code chunk number 7: c4-s4-2
###################################################
rabbit$notinf = rabbit$n - rabbit$inf
#Binomial regression
fit = glm(cbind(inf, notinf) ~ offset(log(a)), 
    family = binomial(link = "cloglog"),
    data = rabbit, subset = a<12)
#Plot data
symbols(rabbit$inf/rabbit$n ~ rabbit$a, circles = rabbit$n, 
     inches = 0.5, xlab = "Age", ylab = "Prevalence")
#Predicted curves for <1 and all 
phi = exp(coef(fit))
curve(1-exp(-phi*x), from = 0, to = 12, add = TRUE)
curve(1-exp(-phi*x), from = 0, to = 30, add = TRUE, lty = 2)
1/phi


###################################################
### code chunk number 8: c4-s4-3
###################################################
integrandpc = function(a, up, foi){
  #Find which interval a belongs to
  wh = findInterval(a, sort(c(0,up)))
  #Calculate duration of each interval
  dur = diff(sort(c(0,up)))
  #Evaluate integrand
  inte = ifelse(wh == 1, foi[1]*a, 
       sum(foi[1:(wh-1)]*dur[1:(wh-1)])+
          foi[wh]*(a-up[wh-1]))
  return(inte)
}


###################################################
### code chunk number 9: c4-s4-4
###################################################
llik.pc = function(par, age, num, denom, up) {
    ll = 0
    for (i in 1:length(age)) {
       p = 1 - exp(-integrandpc(a=age[i], up = up, 
          foi = exp(par)))
       ll = ll + dbinom(num[i], denom[i], p, log = T)
    }
return(-ll)
}


###################################################
### code chunk number 10: c4-s4-5
###################################################
x=c(1,4,8,12,18,24,30)
para=rep(.1,length(x))


###################################################
### code chunk number 11: c4-s4-6 (eval = FALSE)
###################################################
## est = optim(par=log(para),fn=llik.pc, age=rabbit$a, 
##      num=rabbit$inf,  denom=rabbit$n, up=x, 
##      method="Nelder-Mead", control=list(trace=2))


###################################################
### code chunk number 12: c4-knitr.Rnw:162-163
###################################################
source("est.q")


###################################################
### code chunk number 13: c4-s4-7
###################################################
round(exp(est$par), 4)


###################################################
### code chunk number 14: c4-s4-8
###################################################
#Make space for left and right axes
par(mar = c(5,5,2,5))
#Add beginning and ends to x and y for step plot
xvals=c(0,x)
yvals=exp(c(est$par, est$par[7]))
plot(xvals, yvals, type="s", xlab="age", ylab="FoI")

#Superimpose predicted curve
par(new=T)
p = rep(0, 28)
for (i in 1:28) {
     p[i] = 1 - exp(-integrandpc(a=i, up = x, 
        foi = exp(est$par)))
}
plot(p~c(1:28), ylim=c(0,1), type="l", col="red", 
     axes=FALSE, xlab=NA, ylab=NA)

#Add right axis and legend
axis(side = 4)
mtext(side = 4, line = 4, "Prevalence")
legend("right", legend=c("FoI", "Prevalence"),
     lty=c(1,1), col=c("black", "red"))


###################################################
### code chunk number 15: c4-s5-1
###################################################
require(splines)
#Degrees-of-freedom
df=7
#Construct dummy lm-object
dl=lm(inf~bs(a,df), data=rabbit)


###################################################
### code chunk number 16: c4-s5-2
###################################################
tmpfn=function(x,dl){
x=predict(dl, newdata=data.frame(a=x))
exp(x)}


###################################################
### code chunk number 17: c4-s5-3
###################################################
tmpfn2 = function(par, data, df){
   #Dummy lm-object
   dl = lm(inf ~ bs(a,df), data = data)
   #Overwrite spline coefficients with new values
   dl$coefficients = par
   #Calculate log-likelihood 
   ll = 0
   for(i in 1:length(data$a)){
     p = 1 - exp(-integrate(tmpfn, 0, i, dl = dl)$value)
     ll = ll + dbinom(data$inf[i], data$n[i], p ,log = T)
   }
 return(-ll)
 }


###################################################
### code chunk number 18: c4-s5-4 (eval = FALSE)
###################################################
##  para = rep(-1, df + 1)
##  dspline = optim(par = para, fn = tmpfn2, data = rabbit, 
##       df = df, method = "Nelder-Mead", control =
##       list(trace = 2, maxit = 2000))


###################################################
### code chunk number 19: c4-knitr.Rnw:254-256
###################################################
para=rep(-1, df+1)
load("dspline.q")


###################################################
### code chunk number 20: c4-s5-5
###################################################
par(mar = c(5, 5, 2, 5)) #Room for two axes
#Overwrite dummy-objects coefficients with MLEs
dl$coefficients = dspline$par 
#Age-prevalce plot
plot(tmpfn(rabbit$a,dl) ~ rabbit$a, type = "l", ylab = "FoI", 
     xlab = "Age (mos)", las = 1)
#Overlay FoI
par(new = TRUE)
p = rep(0, 28)
for (i in 1:28) {
     p[i] = 1 - exp(-integrate(tmpfn, 0, i, 
          dl = dl)$value)
}
plot(p ~ c(1:28), ylim = c(0,1), type = "l", col = "red", 
     axes = FALSE, xlab = NA, ylab = NA)
axis(side = 4, las = 1)
mtext(side = 4, line = 4, "Prevalence")
legend("topright", legend = c("FoI", "Prevalence"),
     lty = c(1,1), col = c("black", "red"))


###################################################
### code chunk number 21: c4-s6-1
###################################################
data(peru)
head(peru)
#Calculate cumulative incidence
peru$cumulative=cumsum(peru$incidence)
#Define denominator
peru$n=sum(peru$incidence)
par(mar = c(5,5,2,5)) #Make room for two axes and plot
#Plot incidence with cumulative overlaid
plot(peru$incidence~peru$age, type="b", xlab="Age", 
     ylab="Incidence")
par(new=T)
plot(peru$cumulative~peru$age, type="l", col="red", 
     axes=FALSE, xlab=NA, ylab=NA)
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative")
legend("right", legend=c("Incidence", "Cumulative"),
     lty=c(1,1), col=c("black", "red"))


###################################################
### code chunk number 22: c4-knitr.Rnw:321-324
###################################################
up=c(1:20,30, 40, 50, 60, 70,100)
para=rep(.1,length(x))
load("est2.rda")


###################################################
### code chunk number 23: c4-s6-2 (eval = FALSE)
###################################################
## #Upper age cut-offs
## up = c(1:20, 30, 40, 50, 60, 70,100)
## para = rep(.1, length(up)) #Inital values
## #Minimize log-likelihood
## est2 = optim(par = log(para),fn = llik.pc, age = peru$age, 
##      num = peru$cumulative, denom = peru$n, up = up, 
##      method = "Nelder-Mead", control =
##      list(trace = 2, maxit = 2000))
## #Step plot
## x = c(0, up)
## y = exp(c(est2$par, est2$par[26]))
## plot(x, y, ylab = "Relative FoI", xlab = "Age", type = "s", 
##      ylim = c(0, 0.25), xlim = c(0, 80))


###################################################
### code chunk number 24: c4-s6-3
###################################################
data3=peru[peru$age<45,]
df=5
para=rep(.1, df+1)


###################################################
### code chunk number 25: c4-knitr.Rnw:358-376
###################################################
#Prediction function
tmpfn=function(x,dl){
    x=predict(dl, newdata=data.frame(age=x))
exp(x)}
#Dummy lm-object
dl=lm(cumulative~bs(age,df), data=data3)
#Log-likelihood function
tmpfn2=function(par,data, df){
    dl=lm(cumulative~bs(age,df), data=data)
    dl$coefficients=par
    ll=0
    for(a in 1:length(data$age)){
      p=((1-exp(-integrate(tmpfn,0,data$age[a],
         dl=dl)$value)))
      ll=ll+dbinom(data$cumulative[a],data$n[a],p,log=T)
    }
 return(-ll)
 }


###################################################
### code chunk number 26: c4-knitr.Rnw:379-381
###################################################
load("dsplinea45df5.q")
dl$coefficients=dspline.a45.df5$par


###################################################
### code chunk number 27: c4-s6-4 (eval = FALSE)
###################################################
## #Fit model
## dspline.a45.df5=optim(par=log(para),fn=tmpfn2,
##      data=data3, df=df, method="Nelder-Mead", 
##      control=list(trace=4, maxit=5000))
## #Overwrite dummy-objects coefficients with MLEs
## dl$coefficients=dspline.a45.df5$par
## plot(exp(predict(dl))~data3$age, xlab="Age", 
##      ylab="Relative FoI", type="l")


###################################################
### code chunk number 28: c4-s6-5
###################################################
(exp(-integrate(tmpfn,0,15,dl=dl)$value))*(1-
     exp(-integrate(tmpfn,15,40,dl=dl)$value))


###################################################
### code chunk number 29: c4-s6-6
###################################################
redn=0.5
(exp(-redn*integrate(tmpfn,0,15,dl=
     dl)$value))*(1-exp(-redn*integrate(tmpfn,
     15,40,dl=dl)$value))


###################################################
### code chunk number 30: c4-s7-1
###################################################
data(mossong)
head(mossong)
x=y=mossong$contactor[1:30]
z=matrix(mossong$contact.rate, ncol=30, nrow=30)
image(x=x, y=y, z=z, xlab="Contactor", 
     ylab="Contactee", col=gray((12:32)/32))
contour(x=x, y=y, z=z, add=TRUE)


###################################################
### code chunk number 31: c4-s7-2
###################################################
plot(apply(z,1,mean)~x, ylab="Total contact rate",
     xlab="Age")


###################################################
### code chunk number 32: c4-s8-1
###################################################
a = c(1/diff(x),0)


###################################################
### code chunk number 33: c4-s8-2
###################################################
require("fields")
n = length(x)
z2 = (z + t(z))/2
z3 = as.vector(z2)
xy = data.frame(x = rep(x[1:n], n), y = rep(y[1:n], each = n))
polysmooth = Tps(xy, z3, df = 100)
surface(polysmooth, xlab = "", ylab = "", 
     col=gray((12:32)/32))


###################################################
### code chunk number 34: c4-s8-3
###################################################
W = matrix(polysmooth$fitted.values[, 
     1]/mean(polysmooth$fitted.values), nrow = n)


###################################################
### code chunk number 35: c4-s8-4
###################################################
siragemod = function(t, logx,  params){
     n=length(params$a)
     x = exp(logx)
     S = x[1:n]
     I = x[(n+1):(2*n)]
     R = x[(2*n+1):(3*n)]
     with(as.list(params), {
          phi = (beta*W%*%I)/N
          dS = c(mu,rep(0,n-1)) - (phi+a)*S + 
              c(0,a[1:n-1]*S[1:n-1])*(1-p) - mu*S
          dI = phi*S + c(0,a[1:n-1]*I[1:n-1]) - 
             (gamma+a)*I - mu*I
          dR =  c(0,a[1:n-1]*S[1:n-1])*p + 
             c(0,a[1:n-1]*R[1:n-1]) + gamma*I - 
             a*R - mu*R
	 res = c(dS/S, dI/I, dR/R)
	 list((res))
     })
}


###################################################
### code chunk number 36: c4-s8-5
###################################################
p.pre = rep(0,n)
pars.pre = list(N = 1, gamma = 365/14, mu = 0.02, sigma = 0.2, 
     beta = 100, W = W,p = p.pre, a = a)
ystart = log(c(S = rep(0.099/n,n), I = rep(0.001/n,n), 
     R = rep(0.9/n,n)))


###################################################
### code chunk number 37: c4-s8-6 (eval = FALSE)
###################################################
## times = seq(0, 500, by = 14/365)
## #Polymod mixing
## out=as.data.frame(ode(ystart, times=times, 
##      func = siragemod, parms = pars.pre))
## par(mfrow = c(1,2)) #Room for side-by-side plots
## #Time series
## matplot(times, exp(out[,32:61]), type = "l", xlab = "Time", 
##     ylab = "Prevalence", xlim = c(50, 90), ylim = c(0, 0.0005))
## #Final age-prevalence curve
## plot(x, t(exp(out[13036, 32:61])*a), ylab = "Prevalence", 
##     xlab = "Age", ylim = c(0, 4E-5))
## #Homogenous mixing:
## pars.pre$W = matrix(1, ncol = 30, nrow = 30)
## out2=as.data.frame(ode(ystart, times=times, 
##     func=siragemod, parms=pars.pre))
## points(x, t(exp(out2[13036,32:61])*a), col=2, pch="*")


