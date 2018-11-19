### R code from vignette source 'c10-knitr.Rnw'


###################################################
### code chunk number 2: c10-s2-1
###################################################
data(dalziel)
NY = na.omit(dalziel[dalziel$loc=="NEW YORK",])
NY = NY[NY$year %in% c(1920:1940),]
plot(NY$decimalYear, sqrt(NY$cases), type = "b", 
     xlab = "Year", ylab = "Sqrt(cases)")
#Susceptible reconstruction and 
#correcting for underreporting
cum.reg = smooth.spline(cumsum(NY$rec), 
     cumsum(NY$cases), df = 5)
D = - resid(cum.reg) #The residuals
rr = predict(cum.reg, deriv = 1)$y
Ic = NY$cases/rr
Dc = D/rr 
#Align lagged variables
seas = rep(1:26, 21)[1:545]
lInew = log(Ic[2:546])
lIold = log(Ic[1:545])
Dold = Dc[1:545]
#TSIR fit
N = NY$pop
offsetN = -log(N[1:545])
lSold = log(0.051*N[1:545] + Dold)
glmfit = glm(lInew ~ -1 +as.factor(seas) + lIold +
     offset(lSold+offsetN))


###################################################
### code chunk number 3: c10-knitr.Rnw:73-100
###################################################
SimTsir2=function(beta, alpha, B, N,  inits = list(Snull = 0, Inull = 0), type = "det"){
    type = charmatch(type, c("det", "stoc"), nomatch = NA)
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
        lambda[i] = beta[((i - 2) %% s) + 1] * S[i - 1] * (I[i - 1]^alpha)/N
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
### code chunk number 4: c10-s2-2
###################################################
sim2 = SimTsir2(beta = exp(glmfit$coef[1:26]), alpha = 0.98,
     B = rep(median(NY$rec), 5200), N = median(N), 
     inits = list(Snull = exp(lSold[1]), Inull = Ic[1]), 
     type = "det")
Sattr = sim2$S[2601:5200]
Iattr = sim2$I[2601:5200]
plot(Sattr, Iattr, log = "y", type = "l",
     xlab = "S", ylab = "I")
points(sim2$S[seq(2601, 5200, by = 26)], 
     sim2$I[seq(2601, 5200, by = 26)],pch = 19, col = "red")
legend("bottomright", c("Trajectory", "Strobe"), 
     pch = c(NA, 19), lty = c(1, NA) , col = c("black","red"))


###################################################
### code chunk number 5: c10-s2-3
###################################################
TSIRlyap=function(I, S, alpha, bt, N){
   IT = length(I)
   s = length(bt)
   j11 = rep(NA, IT);  j12 = rep(NA, IT)
   j21 = rep(NA, IT);  j22 = rep(NA, IT)
    #initial unit vector
    J = matrix(c(1,0),ncol=1)
    #loop over the attractor
    for(i in 1:IT) {
         j11[i] = 1 -  bt[((i - 1) %% s) + 1] * 
            I[i]^alpha/N
         j12[i] = -( bt[((i - 1) %% s) + 1] * S[i] * 
            (I[i]^(alpha - 1) * alpha)/N)
         j21[i] =  bt[((i - 1) %% s) + 1] * I[i]^alpha/N
         j22[i] =  bt[((i - 1) %% s) + 1] * S[i] * 
            (I[i]^(alpha - 1) * alpha)/N
    
         J = matrix(c(j11[i],j12[i],j21[i],j22[i]), 
            ncol = 2, byrow = TRUE)%*%J
    }
   res = list(lyap = log(norm(J))/IT, j11 = j11, j12 = j12, 
         j21 = j21, j22 = j22, I = I, S = S, alpha = alpha, 
         bt = bt, N = N)
   class(res) = "lyap"
   return(res)
}


###################################################
### code chunk number 6: c10-s2-4
###################################################
nylyap = TSIRlyap(I = Iattr, S = Sattr, alpha = 0.98, 
     bt = exp(glmfit$coef[1:26]), N = median(N))
nylyap$lyap


###################################################
### code chunk number 7: c10-s3-1
###################################################
TSIRllyap = function(x, m = 1){
     llyap = rep(NA, length(x$I))
     for(i in 1:(length(x$I)-m)){
         J = matrix(c(1,0,0,1), ncol = 2)
         for(k in 0:(m-1)){
             J  =  matrix(c(x$j11[(i+k)], x$j12[(i+k)], 
                x$j21[(i+k)], x$j22[(i+k)]), ncol  =  2, 
                byrow = TRUE)%*%J}
         llyap[i] = log(max(abs(eigen(J)$values)))/m
}
res = list(llyap = llyap, I = x$I, S = x$S)
class(res) = "llyap"
return(res)
}


###################################################
### code chunk number 8: c10-s3-2
###################################################
plot.llyap = function(x, inches = .5){
   pm = x$llyap>0
   plot(NA, xlim = range(x$S), ylim = range(x$I), xlab = "S", 
        ylab = "I", log = "y")
   symbols(x$S[pm], x$I[pm], circles = x$llyap[pm], 
        inches = inches, add = TRUE)
   symbols(x$S[!pm],x$I[!pm], squares = -x$llyap[!pm], 
        inches = inches, add = TRUE, bg = 2)
}


###################################################
### code chunk number 9: c10-s3-3
###################################################
nyllyap = TSIRllyap(nylyap, m = 5)
plot.llyap(nyllyap, inches = 0.15)


###################################################
### code chunk number 10: c10-s3-4
###################################################
beta = c(27.71, 43.14, 37.81, 33.69, 31.64, 32.10, 30.16, 
 24.68, 30.19, 31.53, 30.31, 26.02, 26.57, 25.68, 23.21, 
 19.21, 17.50, 20.50, 29.92, 35.85, 32.65, 28.34, 31.11, 
 29.89, 26.89, 39.38)


###################################################
### code chunk number 11: c10-s3-5
###################################################
sim=SimTsir2(beta=beta, alpha=0.98, 
B=rep(2083, 5200), N=3.3E6, inits=list(Snull=133894, 
     Inull=474))
Sattr=sim$S[5149:5200]
Iattr=sim$I[5149:5200]
lonlyap=TSIRlyap(I=Iattr, S=Sattr, alpha=0.98, 
     bt=beta, N=3.3E6)
lonlyap$lyap
lonllyap=TSIRllyap(lonlyap, m=1)


###################################################
### code chunk number 12: c10-s3-6
###################################################
pm=lonllyap$llyap>0
plot(NA, xlim=c(1,52), ylim=range(Iattr), 
     xlab="Biweek", ylab="I", log="y")
symbols((1:52)[pm],Iattr[pm], circles=
     lonllyap$llyap[pm], inches=.3, add=TRUE)
symbols((1:52)[!pm],Iattr[!pm], squares=
     -lonllyap$llyap[!pm], inches=.3, add=TRUE, bg=2)


###################################################
### code chunk number 13: c10-s3-7
###################################################
sim = SimTsir2(beta = beta, alpha = 0.98, B = rep(2083, 520), 
   N = 3.3E6, inits = list(Snull = 133894, Inull = 474), 
   type = "det")
plot(sqrt(sim$I), ylab = "Sqrt(Cases)", xlab = "Biweek")
for(i in 1:20){
   sim = SimTsir2(beta = beta, alpha = 0.98, B = rep(2083, 520), 
       N = 3.3E6, inits = list(Snull = 133894, Inull = 474),
       type = "stoc")
lines(sqrt(sim$I))}
sim = SimTsir2(beta = beta, alpha = 0.98, B = rep(2083, 520), 
   N = 3.3E6, inits = list(Snull = 133894, Inull = 474), 
   type = "det")
points(sqrt(sim$I), col = 2)


###################################################
### code chunk number 14: c10-s4-1
###################################################
sirwmod = function(t, logy, parms){
  y = exp(logy)
   S = y[1]
   I = y[2]
   R = y[3]
   W = y[4]
   with(as.list(parms),{
   dS = mu * (1-p) * N  - mu * S  - beta * S * I / N + 
      2*omega * W
   dI = beta * S * I / N - (mu + gamma) * I
   dR = gamma * I - mu * R - 2*omega * R +  
      kappa * beta * W * I / N + mu*p*N
   dW = 2*omega * R - kappa * beta * W * I / N - 
      (2*omega +mu)* W 
   res = c(dS/S, dI/I, dR/R, dW/W)
   list(res)
 })
 }


###################################################
### code chunk number 15: c10-s4-2 (eval = FALSE)
###################################################
## require(deSolve)
## start = log(c(S = 0.06, I = 0.01, R = 0.92, W  =  0.01))
## res = matrix(NA, ncol = 100, nrow = 5000)
## p = seq(0.01, 1, length = 100)
## #Forwards 
## for(i in 1:100){
##   times = seq(0, 200, by = 0.01)
##   paras  = c(mu = 1/70, p=p[i], N = 1, beta = 200, 
##      omega = 1/10, gamma = 17, kappa = 30)
##   out = as.data.frame(ode(start, times, sirwmod, paras))
##   start = c(S = out[20001,2], I = out[20001,3], 
##      R = out[20001,4], W = out[20001,5])
##   res[,i] = out$I[15002:20001]
##   cat(i, "\r")
## }
## #Backwards
## res2 = matrix(NA, ncol = 100, nrow = 5000)
## start = c(S = -1.8, I = -5.8, R = 1.9, W = -1.9)
## for(i in 100:1){
##   paras  = c(mu = 1/70, p=p[i], N = 1, beta = 200, 
##      omega = 1/10, gamma = 17, kappa = 30)
##   out = as.data.frame(ode(start, times, sirwmod, paras))
## start = c(S = out[20001,2], I = out[20001,3], 
##      R = out[20001,4], W = out[20001,5])
## res2[,i] = out$I[15002:20001]
## cat(i, "\r")
## }
## plot(NA, xlim = range(p), ylim = range(res), 
##    ylab = "Log(I)", xlab = "p")
## for(i in 1:100) points(rep(p[i], 2), range(res[,i]))
## for(i in 1:100) points(rep(p[i], 2), 
##    range(res2[,i]), col = 2)


###################################################
### code chunk number 16: c10-knitr.Rnw:367-368
###################################################
require(deSolve)


###################################################
### code chunk number 17: c10-s4-3
###################################################
paras  = c(mu = 1/70, p = 0.2, N = 1, beta = 200, 
     omega = 1/10, gamma = 17, kappa=30)
start = c(S = -1, I = -5, R = 3.3, W = 0)
times = seq(0, 30, by = 1/52)
out = as.data.frame(ode(start, times, sirwmod, paras))
plot(out$time, exp(out$I), xlab = "Year", 
     ylab = "I", type = "l", ylim = c(0,0.05))
start = c(S = -1.8, I = -5.8, R = 1.9, W = -1.9)
times = seq(0, 30, by = 1/52)
out = as.data.frame(ode(start, times, sirwmod, paras))
lines(out$time, exp(out$I), col = 2)


###################################################
### code chunk number 18: c10-s5-1 (eval = FALSE)
###################################################
## # In the absence of seasonality the damping period would be
## cparas  = c(mu = 0.02, N = 1, beta0 = 537, beta1 = 0.3, 
## 	sigma = 36, gamma = 34.3)
## R0= with(as.list(cparas), sigma/(sigma+mu)*1/(gamma+mu)*beta0 *N/N)
## Sstar=1/R0
## Istar=with(as.list(cparas), mu*(1-1/R0)*R0/beta0)
## Estar=with(as.list(cparas), (mu+gamma)*Istar/sigma)
## eq=list(S=Sstar, E=Estar, I=Istar)
##   dS = expression(mu * (N  - S)  - beta0 * S * I / N)
##   dE= expression(beta0 * S * I / N - (mu + sigma) * E)
##   dI = expression(sigma*E - (mu + gamma) * I)
## 
##   j11 = D(dS, "S");  j12 = D(dS, "E");  j13 = D(dS, "I")
##   j21 = D(dE, "S");  j22 = D(dE, "E");  j23 = D(dE, "I")
##   j31 = D(dI, "S");  j32 = D(dI, "E"); j33 = D(dI, "I")
## 
##   J=with(as.list(c(eq, cparas)),
##   matrix(c(eval(j11),eval(j12),eval(j13),
##   eval(j21),eval(j22), eval(j23),
##   eval(j31),eval(j32), eval(j33)), 
##   nrow=3, byrow=TRUE))
##  eigen(J)$values
## 2*pi/(Im(eigen(J)$values[2]))


###################################################
### code chunk number 19: c10-s5-2
###################################################
times = seq(0, 100, by = 1/120)
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
cparas  = c(mu = 0.02, N = 1, beta0 = 537, beta1 = 0.3, 
     sigma = 36, gamma = 34.3)
out = as.data.frame(ode(start, times,seirmod2, cparas))
plot(out$time, out$I, type = "l", xlab = "Year", ylab = 
     "Prevalence", xlim = c(91,100), ylim = c(0, 0.0015))


###################################################
### code chunk number 20: c10-s5-3
###################################################
"
double rate[8];		// transition rates
double trans[8];	// transition numbers

double beta = beta0*(1+beta1*cos(2*M_PI*t));//transmission
double lambda = (iota+I*beta)/pop;           // force of infection

// transition rates
rate[0] = b*pop;	 // birth of S
rate[1] = lambda;     // infection of S
rate[2] = mu;       	 // death of S

  rate[3] = sigma;       // latent period of E
  rate[4] = mu;          // death of E

  rate[5] = gamma;     // recovery of I
  rate[6] = mu;         // death of I

  rate[7] = mu;      // death of R

  // compute the transition numbers
  trans[0] = rate[0];
  trans[1] = rate[1]*S;
  trans[2] = rate[2]*S;
  trans[3] = rate[3]*E;
  trans[4] = rate[4]*E;
  trans[5] = rate[5]*I;
  trans[6] = rate[6]*I;
  trans[7] = rate[7]*R;

  // balance the equations
  DS = trans[0] - trans[1] - trans[2];
  DE = trans[1] - trans[3] - trans[4];
  DI = trans[3] - trans[5] - trans[6];
  DR = trans[5] - trans[7];
  Dinc = trans[5];  // incidence is modeled as cumulative recovery
  Dpop = trans[0] - trans[2] - trans[4] - trans[6] - trans[7];
" -> skel

"
  double rate[8];		// transition rates
  double trans[8];	// transition numbers

  double beta = beta0*(1+beta1*cos(2*M_PI*t)); // seasonal transmission
  double dW = rgammawn(alpha,dt);  // white noise
  double lambda = (iota+I*beta*dW/dt)/pop; 

  // transition rates
  rate[0] = b*pop;	 // birth of S
  rate[1] = lambda;     // infection of S
  rate[2] = mu;       	 // death of S

  rate[3] = sigma;       // termination of latent period of E
  rate[4] = mu;          // death of E

  rate[5] = gamma;     // recovery of I
  rate[6] = mu;         // death of I

  rate[7] = mu;      // death of R

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2, S, &rate[1], dt, &trans[1]);
  reulermultinom(2, E, &rate[3], dt, &trans[3]);
  reulermultinom(2, I, &rate[5], dt, &trans[5]);
  reulermultinom(1, R, &rate[7], dt, &trans[7]);

  // balance the equations
  S += trans[0] - trans[1] - trans[2];
  E += trans[1] - trans[3] - trans[4];
  I += trans[3] - trans[5] - trans[6];
  R += trans[5] - trans[7];
  inc += trans[5];  // incidence is modeled as cumulative recovery
  pop = S + E + I + R;
" -> rproc

## Observation model simulator (negative binomial)
"
  double mean = rho*inc;
  double size = 1/theta;
  reports = rnbinom_mu(size,mean);
" -> robs

## Observation model likelihood (negative binomial)
"
  double mean = rho*inc;
  double size = 1/theta;
  lik = dnbinom_mu(reports,size,mean,give_log);
" -> dobs

"
  S = nearbyint(pop0*S0);
  E = nearbyint(pop0*E0);
  I = nearbyint(pop0*I0);
  R = nearbyint(pop0*R0);
  pop = S+E+I+R;
  inc = 0;
" -> rinit


###################################################
### code chunk number 21: c10-s5-4
###################################################
dat = data.frame(time = seq(0, 500, by = 1/52), reports = NA) 
seirp = pomp(dat, times = "time",t0 = 0,
    rprocess = euler.sim(Csnippet(rproc),delta.t = 1/52/20),
    skeleton = vectorfield(Csnippet(skel)),
    rmeasure = Csnippet(robs),
    dmeasure = Csnippet(dobs),
    initializer = Csnippet(rinit),
    params = c(iota = 0,beta0 = 537,beta1 = 0.3,sigma = 36,
         gamma = 34.3,alpha = 0.015,rho = 0.6,theta = 1,
         b = 0.02,mu = 0.02,pop0 = 5e8,
         S0 = 0.06,E0 = 0,I0 = 0.001,R0 = 0.939),
         paramnames = c("iota","beta0","beta1","gamma",
         "sigma","alpha","rho","theta", "b","mu","pop0",
         "S0","E0","I0","R0"),
    statenames = c("S","E","I","R","inc","pop"),
    zeronames = "inc")


###################################################
### code chunk number 22: c10-s5-5 (eval = FALSE)
###################################################
## detsim = trajectory(seirp, as.data.frame = TRUE)
## plot(detsim$time, detsim$I/5E8, type = "l", 
##      xlim = c(101, 110), xlab = "year", ylab = "prevalence")


###################################################
### code chunk number 23: c10-s5-6 (eval = FALSE)
###################################################
## par(mfrow = c(1,2))
## stocsim = simulate(seirp, seed = 3495135, 
##      as.data.frame = TRUE, nsim = 1)
## plot(stocsim$time, stocsim$I/5E8, type = "l", xlim = c(150, 
##      250), xlab = "Year", ylab = "Prevalence")
## sel = seq(105, length(stocsim$I), by = 52)
## plot(stocsim$S[sel]/5E8, stocsim$I[sel]/5E8, 
##      log = "xy", xlab = "S", ylab = "I")
## sel2 = sel[401:500]
## points(detsim$S[sel2]/5E8, detsim$I[sel2]/5E8, col = 2, 
##      pch = 21, bg = 2)


###################################################
### code chunk number 24: c10-s6-1
###################################################
starts = expand.grid(S = seq(0.02, 0.1, length = 10), 
     E = seq(1E-8, 0.0125, length = 10), 
     I = seq(1E-8, 0.005, length = 10))
starts$R = 1-apply(starts,1,sum)


###################################################
### code chunk number 25: c10-s6-2
###################################################
#times for integration
itimes  = seq(0, 100, by = 1/52)
#points for stroboscopic section
isel = seq(1, length(itimes), by = 52)
#list to fill with results
cporbs = list(S = matrix(NA, ncol = dim(starts)[1], 
     nrow = length(isel)), I = matrix(NA, 
     ncol = dim(starts)[1], nrow = length(isel)))


###################################################
### code chunk number 26: c10-knitr.Rnw:604-624
###################################################
seirmod2=function(t, y, parms){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

  mu=parms["mu"]
  N=parms["N"]
  beta0=parms["beta0"]
  beta1=parms["beta1"]
  sigma=parms["sigma"]
  gamma=parms["gamma"]

  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
} 


###################################################
### code chunk number 27: c10-s6-3 (eval = FALSE)
###################################################
## for (i in 1:dim(starts)[1]){
##    out2b = as.data.frame(ode(as.numeric(starts[i,]), 
##       itimes, seirmod2, cparas))
##    cporbs$S[,i] = out2b[,2][isel]
##    cporbs$I[,i] = out2b[,4][isel]
## }


###################################################
### code chunk number 28: c10-knitr.Rnw:635-637
###################################################
#save(cporbs, starts, itimes, isel, times, start, out, sel, stocsim, detsim, file="cporbs.Rd")
load("cporbs.Rd")


###################################################
### code chunk number 29: c10-s6-4 (eval = FALSE)
###################################################
## #Invasion orbits
## plot(as.vector(cporbs$S[-c(1:5),]), 
##    as.vector(cporbs$I[-c(1:5),]), pch=20, cex=0.25, 
##    log="xy", ylab="I", xlab="S")
## #Stochastic simulation
## sel=seq(105, length(stocsim$I), by=52)
## points(stocsim$S[sel]/5E8, stocsim$I[sel]/5E8, col=2)
## #Deterministic attractor
## times  = seq(0, 1000, by=1/120)
## start = c(S=0.06, E=0, I=0.001, R = 0.939)
## out = as.data.frame(ode(start, times, seirmod2, cparas))
## sel=seq(120*100, length(times), by=120)
## points(out$S[sel], out$I[sel], pch=21, col="white", 
##      bg="white")


###################################################
### code chunk number 30: c10-s7-1
###################################################
data(pertcop)
require(Rwave)
#Wavelet decompostion
no = 8
nv = 16
a = 2^seq(1, no+1-1/nv, by = 1/nv)
wfit = cwt(sqrt(pertcop$cases), no, nv, plot = FALSE)
wspec  =  Mod(wfit)
#Crazy climber
crcinc<-crc(wspec, nbclimb = 10, bstep = 100)
fcrcinc<-cfamily(crcinc, ptile = 0.5, nbchain = 1000, 
     bstep = 10)


###################################################
### code chunk number 31: c10-s7-2
###################################################
sigind = which((a/52)>3 & (a/52)<4)
noiseind = which((a/52)<0.5)
snr = apply(wspec[, sigind], 1, 
     sum)/apply(wspec[, noiseind], 1, sum)


###################################################
### code chunk number 32: c10-s7-3
###################################################
par(mfrow = c(3,1), mar = c(0,4,2,1))
layout(matrix(c(1,1,2,2,2,3), ncol = 1))
#Top panel
plot(as.Date(pertcop$date), pertcop$cases, xlab = "",
   ylab = "Sqrt(cases)", type = "l", bty = 'l', 
   xlim = c(as.Date("1901-01-01"), 
   as.Date("1938-01-01")), xaxt='n', yaxt='n')
axis(2, at = seq(0,200,by = 100), labels = FALSE)
axis(2, at = seq(50,250,by = 100), labels = TRUE)
#Mid panel
par(mar = c(0,4,0.25,1))
image(x = as.Date(pertcop$date, origin = "1900-01-07"), 
  wspec, col = gray((30:10)/32), y = a/52, ylim = c(0,5), 
  ylab = "Period (year)", main = "", xaxt = "n", yaxt = "n")
contour(x = as.Date(pertcop$date, origin = "1900-01-07"), 
   wspec, y = a/52, ylim = c(0,5), 
   zlim = c(quantile(wspec)[4], max(wspec)), add = TRUE)
axis(2, at = 0:4) 
ridges<-fcrcinc[[1]]
ridges[which(ridges<1.5*10^-5)]<-NA
image(x = as.Date(pertcop$date, origin = "1900-01-07"), 
   y = a/52, z = ridges, add = TRUE, col = "black")
#Bottom panel
par(mar = c(3,4,0.25,1), tcl = -0.4)
plot(x = as.Date(pertcop$date, origin = "1900-01-07"), snr, 
   type = "l", bty = "l",xaxt = "n", yaxt = "n", ylab = "SNR")
axis.Date(1, at = seq(as.Date("1900-01-01"), 
  as.Date("1938-01-01"), "years")) 


###################################################
### code chunk number 33: c10-s7-4 (eval = FALSE)
###################################################
## dat = data.frame(time = seq(0, 500, by = 1/52), reports = NA) 
## sds = rep(NA, 126)
## alpha = seq(0,0.025, by = 0.0002)
## Smat = Imat = matrix(NA, nrow = 26001, ncol = 126)
## for(i in 1:126){
##   seirp = pomp(dat, times = "time", t0 = 0,
##     rprocess = euler.sim(Csnippet(rproc), delta.t = 1/52/20),
##     skeleton = vectorfield(Csnippet(skel)),
##     rmeasure = Csnippet(robs),
##     dmeasure = Csnippet(dobs),
##     initializer = Csnippet(rinit),
##     params = c(iota = 0, beta0 = 537, beta1 = 0.3, sigma = 36,
##        gamma = 34.3, alpha = alpha[i], rho = 0.6, theta = 1,
##        b = 0.02, mu = 0.02, pop0 = 5e8,
##        S0 = 0.06, E0 = 0, I0 = 0.001, R0 = 0.939),
##        paramnames = c("iota", "beta0", "beta1", "gamma",
##        "sigma", "alpha", "rho", "theta", "b", "mu", "pop0",
##        "S0", "E0" ,"I0", "R0"),
##     statenames = c("S","E","I","R","inc","pop"),
##     zeronames = "inc"
##     )
##   stocsim$I[26001] = 0
##   j = -1
##   while(stocsim$I[26001]==0){
##   j = j+1
##   stocsim = simulate(seirp, seed = 3495131+j, 
##        as.data.frame = TRUE, nsim = 1)
##   sds[i] = 3495131+j
##   }
##   Imat[,i] = stocsim$I
##   Smat[,i] = stocsim$S
## }


###################################################
### code chunk number 34: c10-knitr.Rnw:763-764
###################################################
load("stocres.Rdata")


###################################################
### code chunk number 35: c10-s7-5 (eval = FALSE)
###################################################
## predn = rep(NA, 126)
## #Set the number of "octaves" and "voices"
## no = 8; nv = 10
## #then calculate the corresponding periods
## a = 2^seq(1,no+1-1/nv, by = 1/nv)
## sel2 =  a> = 39 & a <260 #Multiannual signal
## sel = a<39 #High frequency noise
## for(i in 1:126){
##      wfit = cwt(sqrt(Imat[,i]), no, nv, plot = FALSE)
##      wspec = Mod(wfit)
##      predn[i] = sum(wspec[,sel2])/sum(wspec[sel])
## }
## plot(alpha, predn, xlab = "Noise variance", 
##      ylab = "'Predictability'")


###################################################
### code chunk number 36: c10-knitr.Rnw:784-786
###################################################
load("stocres.Rdata")
plot(alpha, predn, xlab="Noise variance", ylab="'Predictability'")


###################################################
### code chunk number 37: c10-s8-1 (eval = FALSE)
###################################################
## require(nlts)
## llcv=rep(NA,126)
## for(i in 1:126){
## llfit=ll.order(sqrt(Imat[seq(521,26001, by=52),i]), 
##    step=1, order=1:5, bandwidt = seq(0.5, 1.5, 
##    by = 0.5), cv=FALSE) 
## llcv[i]=min(llfit$grid$GCV)
## }
## plot(llcv~alpha, ylab="GCV", xlab="Noise variance")


###################################################
### code chunk number 38: c10-knitr.Rnw:826-827
###################################################
require(nlts)


###################################################
### code chunk number 39: c10-s8-2
###################################################
x=sqrt(Imat[seq(521,1040, by=1),10])
sim=ll.edm(x=x, order=3, bandwidth=0.8)
plot(x, type="b", ylab="Prevalence", xlab="Week")
lines(sim, col=2)
legend("topright", legend=c("Simulation", "EDM"),
   lty=c(1,1), pch=c(1, NA), col=c("black", "red"))


###################################################
### code chunk number 40: c10-knitr.Rnw:842-843
###################################################
require(pomp)


###################################################
### code chunk number 41: c10-s9-4
###################################################
dat=data.frame(time=seq(0, 500, by=1/52), reports=NA) 
seirp=pomp(dat, times="time",t0=0,
  rprocess=euler.sim(Csnippet(rproc),delta.t=1/52/20),
  skeleton=vectorfield(Csnippet(skel)),
  rmeasure=Csnippet(robs),
  dmeasure=Csnippet(dobs),
  initializer=Csnippet(rinit),
  params=c(iota=0,beta0=537,beta1=0.3,sigma=36,
      gamma=34.3,alpha=0.015,rho=0.6,theta=1,
       b=0.02,mu=0.02,pop0=5e8,
      S0=0.06,E0=0,I0=0.001,R0=0.939),
  paramnames=c("iota","beta0","beta1","gamma","sigma",
      "alpha","rho","theta", "b","mu","pop0",
      "S0","E0","I0","R0"),
  statenames=c("S","E","I","R","inc","pop"),
  zeronames="inc")


