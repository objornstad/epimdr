### R code from Chapter 8
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866



###################################################
### code chunk number 3: c8-s2-1
###################################################
rlist=c(quote(mu * (S + I + R)), #Births
quote(mu * S), #Sucseptible deaths
quote(beta * S * I /(S + I + R)), #Infection
quote(mu * I), #Infected death
quote(gamma * I), #Recovery
quote(mu * R)) #Recovered death


###################################################
### code chunk number 4: c8-s2-2
###################################################
emat=matrix(c(1,0,0,
-1,0,0,
-1,1,0,
0,-1,0,
0,-1,1,
0,0,-1),
ncol=3, byrow=TRUE)


###################################################
### code chunk number 5: c8-s2-3
###################################################
gillespie = function(rateqs, eventmatrix, parameters, 
     initialvals, numevents){
   res = data.frame(matrix(NA, ncol = length(initialvals)+1, 
      nrow = numevents+1))
   names(res) = c("time", names(inits))
   res[1,] = c(0, inits)
   for(i in 1:numevents){
   #evaluate rates
   rat = sapply(rateqs, eval, 
      as.list(c(parameters, res[i,])))
   #update clock
   res[i+1,1] = res[i,1]+rexp(1, sum(rat))
   #draw event
   whichevent = sample(1:nrow(eventmatrix), 1, prob = rat)
   #updat states
   res[i+1,-1] = res[i,-1]+eventmatrix[whichevent,]
   }
return(res)
}


###################################################
### code chunk number 6: c8-s2-4
###################################################
paras  = c(mu = 1, beta =  500, gamma = 365/20)
inits = c(S = 100, I = 2, R = 0)
sim=gillespie(rlist, emat, paras, inits, 1000)
matplot(sim[,1], sim[,2:4], type = "l", ylab = "Numbers", 
     xlab = "Time", log = "y")
legend("topright", c("S", "I", "R"), lty = c(1, 1, 1), 
     col = c(1, 2, 3))


###################################################
### code chunk number 7: c8-s2-5
###################################################
emat2 = matrix(c(1, 0, 0, 0,
  -1, 0, 0, 0,
  -1, 1, 0, 0,
  0, -1, 0, 0,
  0, -1, 1, 0,
  0, 0, -1, 0,
  0, 0, -1, 1,
  0, 0, 0, -1),
  ncol = 4, byrow = TRUE)


###################################################
### code chunk number 8: c8-s2-6
###################################################
rlist2 = c(quote(mu * (S+E+I+R)), quote(mu * S), 
     quote(beta * S * I/(S+E+I+R)), quote(mu*E), 
     quote(sigma * E), quote(mu * I), 
     quote(gamma * I), quote(mu*R))


###################################################
### code chunk number 9: c8-s2-7
###################################################
tau = function(rateqs, eventmatrix, parameters, 
     initialvals, deltaT, endT){
time = seq(0, endT, by = deltaT)
res = data.frame(matrix(NA, ncol = length(initialvals)+1, 
     nrow = length(time)))
res[,1] = time
names(res) = c("time", names(inits))
res[1,] = c(0, inits)
for(i in 1:(length(time)-1)){
     #calculate overall rates
     rat = sapply(rateqs, eval, as.list(c(parameters, 
          res[i, ])))
     evts = rpois(1,  sum(rat)*deltaT)
     if(evts>0){
     #draw events
     whichevent = sample(1:nrow(eventmatrix), evts, 
          prob = rat, replace = TRUE)
     mt = rbind(eventmatrix[whichevent,], 
          t(matrix(res[i,-1])))
     mt = matrix(as.numeric(mt), ncol = ncol(mt))
     #update states
     res[i+1, -1] = apply(mt,2,sum)
     res[i+1, ][res[i+1,]<0] = 0
     }
     else{ #if no events in delaT
     res[i+1,-1] = res[i,-1]
     }}
return(res)
}


###################################################
### code chunk number 10: c8-s2-8
###################################################
paras  = c(mu = 1, beta =  1000, 
     sigma = 365/8, gamma = 365/5)
inits = c(S = 999, E = 0, I = 1, R = 0)
sim2 = tau(rlist2, emat2, paras, inits, 1/365, 2)
matplot(sim2[,1], sim2[,2:5], type = "l", log = "y", 
     ylab = "Numbers", xlab = "Time")
legend("bottomright", c("S", "E", "I", "R"), 
     lty = c(1, 1, 1, 1), col = c(1, 2, 3, 4))


###################################################
### code chunk number 11: c8-s3-1
###################################################
require(deSolve)
seirmod=function(t, y, parms){
   S=y[1]
   E=y[2]
   I=y[3]
   R=y[4]

   with(as.list(parms),{
   dS = mu * (N  - S)  - beta * S * I / N
   dE = beta * S * I / N - (mu + sigma) * E
   dI = sigma * E - (mu + gamma) * I
   dR = gamma * I - mu * R
   res=c(dS, dE, dI, dR)
   list(res)
 })
 }


###################################################
### code chunk number 12: c8-s3-2
###################################################
lfn = function(p){
     times = seq(0, 2, by = 1/365)
     start = c(S = 999, E = 0, I = 1, R = 0)
     paras = exp(c(mu = p[1], N = p[2], beta = p[3], 
        sigma = p[4], gamma = p[5]))
     out = as.data.frame(ode(start, times = times, 
        seirmod, paras))
     n = length(sim2$I)
     rss = sum((sim2$I-out$I)^2)
     return(log(rss)*(n/2)-n*(log(n)-log(2*pi)-1)/2)
}


###################################################
### code chunk number 13: c8-s3-3 
###################################################
#initial values for mu, N, beta, sigma, gamma
 paras0  = log(c(2, 500, 500, 365/7, 365/7))
 fit=optim(paras0, lfn, hessian=TRUE)


###################################################
### code chunk number 15: c8-s3-4
###################################################
times = seq(0, 2, by=1/365)
paras  = exp(c(mu = fit$par[1], N = fit$par[2], 
     beta = fit$par[3], sigma = fit$par[4], 
     gamma = fit$par[5]))
start = c(S=999, E=0, I=1, R = 0)
out = as.data.frame(ode(start, times, seirmod, paras))
plot(out$time, out$I, xlab="Time", ylab="Prevalence", 
     type="l")
lines(sim2$time, sim2$I, col=2, type="l")
legend("topright", c("Gillespie simulation", 
     "SEIR fit"), lty=c(1,1), col=c(2,1))


###################################################
### code chunk number 16: c8-s4-1
###################################################
#MLEs
 round(exp(fit$par), 4)
 
#Approximate SEs
round(exp(sqrt(diag(solve(fit$hessian)))), 4)

#Correlation matrix
round(cov2cor(solve(fit$hessian)),4)


###################################################
### code chunk number 17: c8-s5-1
###################################################
times  = seq(0, 10, by=1/52)
paras  = c(mu = 1/50, N = 1, beta =  1000, 
     sigma = 365/8, gamma = 365/5)
     start = c(S=0.06, E=0, I=0.001, R = 0.939)
out = as.data.frame(ode(start, times, seirmod, paras))


###################################################
### code chunk number 19: c8-s5-2 (eval = FALSE)
###################################################
datay=jitter(out$I, amount=1e-4)
plot(times, datay, ylab="Infected", xlab="Time")
lines(times, out$I, col=2)


###################################################
### code chunk number 20: c8-s5-3
###################################################
lfn = function(p, data){
     times = seq(0, 10, by = 1/52)
     start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
     paras = c(mu = p[1], N = p[2], beta = p[3], 
        sigma = p[4], gamma = p[5])
     out = as.data.frame(ode(start, times = times, 
        seirmod, paras))
     n = length(data)
     rss = sum((data-out$I)^2)
     return(log(rss)*(n/2)-n*(log(n)-log(2*pi)-1)/2)
 }


###################################################
### code chunk number 21: c8-s5-4 (eval = FALSE)
###################################################
# mu, N, beta, sigma, gamma
 paras0  = c(1/30, 1, 1500, 365/4, 365/10)
 fit=optim(paras0, lfn, data=datay, hessian=TRUE)


###################################################
### code chunk number 22: c8-s5-5
###################################################
#MLEs:
round(fit$par, 3)
#Approximate SEs:
round(sqrt(diag(solve(fit$hessian))), 3)
#Correlation matrix:
round(cov2cor(solve(fit$hessian)),3)


###################################################
### code chunk number 23: c8-s6-1
###################################################
 data(flu)


###################################################
### code chunk number 24: c8-s6-2
###################################################
sirmod = function(t, y, params) {
S = y[1]
I = y[2]
R = y[3]
with(as.list(params), {
dS = -beta * S * I/N
dI = beta * S * I/N - gamma * I
dR = gamma * I
res = c(dS, dI, dR)
list(res)
})
}


###################################################
### code chunk number 25: c8-s6-3
###################################################
lfn2 = function(p, I, N) {
     times = seq(1, 14, by = 1)
     start = c(S = N, I = 1, R = 0)
     paras = c(beta = p[1], gamma = p[2], N = N)
     out = as.data.frame(ode(start, times = times, 
        sirmod, paras))   
     n = length(I)
     rss = sum((I - out$I)^2)
     return(log(rss) * (n/2) - n * (log(n) - 
        log(2 * pi) - 1)/2)
}


###################################################
### code chunk number 26: c8-s6-4
###################################################
#beta, gamma
paras0 = c(1.5, 1/2)
flufit = optim(paras0, lfn2, I = flu$cases, N = 763, 
     hessian = TRUE)


###################################################
### code chunk number 27: c8-s6-5
###################################################
#parameters
flufit$par
#R0:
flufit$par[1]/flufit$par[2]


###################################################
### code chunk number 28: c8-s6-6
###################################################
times = seq(1, 20, by = .1)
start = c(S = 762, I = 1, R = 0)
paras = c(beta = flufit$par[1], gamma = flufit$par[2], N = 763)
out = as.data.frame(ode(start, times = times, 
     sirmod, paras))
plot(out$time, out$I, ylab = "Prevalence", 
     xlab = "Day", type = "l")
points(flu$day, flu$cases)


###################################################
### code chunk number 29: c8-s7-1
###################################################
y = as.vector(table(niamey_daily))


###################################################
### code chunk number 30: c8-s7-2
###################################################
times = unique(niamey_daily$day)
paras = c(mu = 0, N = 1, beta =  5, sigma = 1/8, 
     gamma =1/5)
start = c(S = 0.999, E = 0, I = 0.001, R = 0, K  = 0)


###################################################
### code chunk number 31: c8-s7-3
###################################################
seirkmod=function(t, x, params){
   S=x[1]
   E=x[2]
   I=x[3]
   R=x[4]
   K=x[5]

   with(as.list(params),{
   dS = mu * (N  - S)  - beta * S * I / N
   dE = beta * S * I / N - (mu + sigma) * E
   dI = sigma * E - (mu + gamma) * I
   dR = gamma * I - mu * R
   dK = sigma * E
   res=c(dS, dE, dI, dR, dK)
   list(res)
 })
 }


###################################################
### code chunk number 32: c8-s7-4
###################################################
lfn4 = function(p, I){
times = unique(niamey_daily$day)
xstart = c(S = (p[1]-5)/p[1], E = 0, I = 5/p[1], R = 0, 
     K  = 0)
paras = c(mu = 0, N = p[1], beta = p[2], sigma = 1/8, 
     gamma = 1/5)
out = as.data.frame(ode(xstart, times = times, seirkmod, 
     paras))
predinci = c(xstart["I"], diff(out$K))*p[1]
ll = -sum(dpois(I, predinci, log = TRUE))
return(ll)
}


###################################################
### code chunk number 33: c8-s7-5
###################################################
#N, beta
paras0  = c(11000, 5)
measfit = optim(paras0, lfn4, I = y, hessian = TRUE)
day = 1:230
xstart = c(S= (measfit$par[1]-5)/measfit$par[1], E = 0, 
     I = 5/measfit$par[1], R = 0, K  = 0)
paras = c(mu = 0, N = measfit$par[1], beta = measfit$par[2], 
     sigma = 1/8, gamma = 1/5)
out = as.data.frame(ode(xstart, times = day, 
     seirkmod, paras))
plot(table(niamey_daily), xlab = "Day", ylab = "Incidence")
lines(out$time, c(xstart["I"], diff(out$K))*
     measfit$par[1], col = 2, lwd = 2)


###################################################
### code chunk number 34: c8-s7-6
###################################################
with(as.list(paras),
     sigma/(sigma+mu)*1/(gamma+mu)*beta /N)


###################################################
### code chunk number 35: c8-s8-1
###################################################
sivmod = function(t, x, parms){
    S = x[1]
    E = x[2]
    I = x[3]
    R = x[4]
    K = x[5]
    with(as.list(parms),{
      Q = ifelse(t<T | t>T+Dt,0,(-log(1-P)/Dt))
      dS = -B*S*I-q*Q*S
      dE = B*S*I-r*E
      dI = r*E - g*I
      dR = g*I+q*Q*S
      dK = r*E
      res = c(dS, dE, dI, dR, dK)
      list(res)
    })
  }

retrospec = function(R, day, vaccine_efficacy, 
  target_vaccination,intervention_length, 
  mtime, LP = 7, IP = 7, N = 10000){
  steps =1:mtime
        out = matrix(NA, nrow = mtime, ncol = 3)
  #starting values
  xstrt = c(S = 1-1/N, E = 0, I = 1/N, R = 0,K = 0)   
  beta = R/IP         #transmission rate
  #Without ORV
  par = c(B = beta, r = 1/LP, g = 1/IP, q = vaccine_efficacy,
       P = 0, Dt = 0, T = Inf, R = R)
  outv = as.data.frame(ode(xstrt, steps, sivmod, par))
  fsv = max(outv$K)
  #With ORV
  par = c(B = beta, r = 1/LP, g = 1/IP, q = vaccine_efficacy,
             P = target_vaccination, Dt = 
             intervention_length, T = day)
  outi = as.data.frame(ode(xstrt, steps, sivmod, par))
  fsi = max(outi$K)
  res = list(redn = fsi/fsv, out = outv, orv = outi, B = par["B"], 
     r = par["r"], g = par["g"], q = par["q"], P = par["P"], 
     Dt = par["Dt"], T = par["T"], R = R)
  class(res) = "retro"
  return(res)
}


###################################################
### code chunk number 36: c8-s8-2
###################################################
plot.retro = function(x){ 
   plot(x$out[,1], x$out[,"I"], type = "l", ylim = c(0,
       max(x$out[,"I"])), xlab = 'Day', ylab = 'Prevalence')
   polygon(c(x$T, x$T, x$T+x$Dt,
      x$T+x$Dt), c(-0.1,1,1,-.1), col = "gray")
   lines(x$out[,1], x$out[,"I"])
   lines(x$orv[,1], x$orv[,"I"], col = "red")
   title(paste("Final size: ", round(100*(x$redn),1), 
      "% (R=",x$R,", target=", 100*x$P, "%)", sep=""))
   legend(x = "topleft", legend = c("Natural epidemic", 
      "With ORV"), col = c("black", "red"), lty = c(1,1))
   text(x = x$T + x$Dt, y = 0, pos = 4,
      labels = paste(x$intervention_length, 
      "ORV from ", x$T))
}


###################################################
### code chunk number 37: c8-s8-3
###################################################
red1 = retrospec(R = 1.8, 161, vaccine_efficacy = 0.85, 
     target_vaccination = 0.5, intervention_length = 10,
     mtime = 250, LP = 8, IP = 5, N = 16000)
red2 = retrospec(R = 1.8, 161+14, vaccine_efficacy = 0.85, 
     target_vaccination = 0.5, intervention_length = 10, 
     mtime = 250, LP = 8, IP = 5, N = 16000)
red3 = retrospec(R = 1.8, 161+28, vaccine_efficacy = 0.85, 
     target_vaccination = 0.5, intervention_length = 10, 
     mtime = 250, LP = 8, IP = 5, N = 16000)
1-red1$redn
1-red2$redn
1-red3$redn


###################################################
### code chunk number 38: c8-s8-4
###################################################
plot(red1)


###################################################
### code chunk number 39: c8-knitr.Rnw:666-668 (eval = FALSE)
###################################################
require(shiny)
orv.app


