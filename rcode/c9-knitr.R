### R code from vignette source 'c9-knitr.Rnw'

###################################################
### code chunk number 2: c9-s1-1
###################################################
data(rabies)
matplot(rabies[,2:7], ylab="Cases", xlab="Month")
legend("topright", c("CT", "DE", "MD", "MA", 
     "NJ", "NY"), pch = as.character(1:6), col = 1:6)


###################################################
### code chunk number 3: c9-s3-1
###################################################
parms  = c(mu = 1/(50*52), N = 1, beta =  2.5, 
      gamma = 1/2)
N = parms["N"]
gamma = parms["gamma"]
beta = parms["beta"]
mu = parms["mu"]
Istar = as.numeric(mu*(beta/(gamma+mu)-1)/beta)
Sstar = as.numeric((gamma+mu)/beta)
Sstar
Istar


###################################################
### code chunk number 4: c9-s3-2
###################################################
require(nleqslv)
rootfn = function(x, params){
   r = with(as.list(params),
      c(mu * (N  - x[1])  - beta * x[1]* x[2] / N,
      beta * x[1] * x[2] / N - (mu + gamma) * x[2],
      gamma *x[2] - mu*x[3]))
   r
}
parms  = c(mu = 1/(50*52), N = 1, beta =  2.5, 
   gamma = 1/2)
ans = nleqslv(c(0.1,0.5, 0.4), fn = rootfn, 
   params = parms)
ans$x


###################################################
### code chunk number 5: c9-s3-4
###################################################
ans = grid = expand.grid(seq(0,1, by = .25), 
     seq(0,1, by = .25), seq(0,1, by=.25))
ans[,] = NA
for(i in 1:nrow(ans)){
ans[i,] = nleqslv(as.numeric(grid[i,]), fn = rootfn, 
     params = parms)$x
}
ans2 = round(ans, 4)
ans2[!duplicated(ans2),]


###################################################
### code chunk number 6: c9-s3-5
###################################################
sirmod = function(t, y, parameters){
     S = y[1]
     I = y[2]
     R = y[3]
     beta = parameters["beta"]
     mu = parameters["mu"]
     gamma = parameters["gamma"]
     N = parameters["N"]
     dS = mu * (N  - S)  - beta * S * I / N
     dI = beta * S * I / N - (mu + gamma) * I
     dR = gamma*I - mu*R
     res = c(dS, dI, dR)
     list(res)
 }

require(rootSolve)
parms  = c(mu = 1/(50*52), N = 1, 
     beta =  2.5, gamma = 1/2)
equil = runsteady(y = c(S = 1 - 1E-4, I = 1E-4, R = 0), 
     times = c(0, 1E05), func = sirmod, parms = parms)
round(equil$y, 4)


###################################################
### code chunk number 7: c9-s4-1
###################################################
  dS = expression(mu * (1  - S)  - beta * S * I / 1)
  dI = expression(beta * S * I / 1 - (mu + gamma) * I)
  j11 = D(dS, "S"); j12 = D(dS, "I")
  j21 = D(dI, "S"); j22 = D(dI, "I")


###################################################
### code chunk number 8: c9-s4-2
###################################################
vals  = list(mu = 1/(50*52), N = 1, beta =  2.5, 
     gamma = 1/2, S=Sstar, I=Istar)
J=with(vals, matrix(c(eval(j11), eval(j12), eval(j21), 
     eval(j22)),  ncol = 2, byrow = T))
eigen(J, only.values = TRUE)$values


###################################################
### code chunk number 9: c9-s4-3
###################################################
2*pi/Im(eigen(J)$values[1])


###################################################
### code chunk number 10: c9-4-4
###################################################
vals  = list(mu = 1/(50*52), N = 1, beta =  2.5, 
     gamma = 1/2, S = 1, I = 0)
J = with(vals,
     matrix(c(eval(j11), eval(j12),  eval(j21), 
      eval(j22)),  ncol = 2, byrow = T))
eigen(J, only.values = TRUE)$values


###################################################
### code chunk number 11: c9-s4-5
###################################################
vals  = list(mu = 1/(50*52), N = 1, beta =  0.3, 
      gamma = 1/2, S=1, I=0)
J=with(vals,
      matrix(c(eval(j11), eval(j12), 
      eval(j21), eval(j22)),  ncol=2, byrow=TRUE))
eigen(J, only.values=TRUE)$values


###################################################
### code chunk number 12: c9-s5-1
###################################################
coyne = function(t, logx, parms){
  x = exp(logx)
  X = x[1]
  H1 = x[2]
  H2 = x[3]
  Y = x[4]
  I = x[5]
  N = sum(x)
  with(as.list(parms),{
  dX = a * (X + I) - beta * X * Y - 
       gamma * N * X  - (b + c) * X
  dH1= rho * beta * X * Y  - gamma * N * H1  -
       (b + sigma + c) * H1
  dH2= (1-rho) * beta * X * Y  - gamma * N * H2  -
       (b + sigma + c) * H2
  dY = sigma * H1  - gamma * N * Y  -
       (b + alpha + c) * Y
  dI = sigma * H2 - gamma * N * I  - (b + c) * I
  res = c(dX/X, dH1/H1, dH2/H2, dY/Y, dI/I)
  list(res)
})
}


###################################################
### code chunk number 13: c9-s5-2
###################################################
times  = seq(0, 50, by=1/520)
paras  = c(gamma = 0.0397, b = 0.836, 
     a = 1.34, sigma = 7.5, 
     alpha = 66.36, beta = 33.25, 
     c = 0, rho = 0.8)
start = log(c(X = 12.69/2, H1 = 0.1, H2 = 0.1, 
     Y = 0.1, I = 0.1))
out = as.data.frame(ode(start, times, coyne, paras))


###################################################
### code chunk number 14: c9-s5-3
###################################################
par(mfrow = c(1,2))  
plot(times, exp(out$Y), ylab = "Infected", xlab = "Time", 
     type = "l")
plot(exp(out$X), exp(out$Y), ylab = "Infected", 
     xlab = "Susceptible", type = "l")


###################################################
### code chunk number 15: c9-s5-4
###################################################
dX = expression(a * (X + I) - beta * X * Y - 
     gamma * (X+H1+H2+Y+I) * X  - (b + c) * X)
dH1= expression(rho * beta * X * Y  - 
     gamma * (X+H1+H2+Y+I) * H1  - (b + sigma + c) * H1)
dH2= expression((1-rho) * beta * X * Y  - 
     gamma * (X+H1+H2+Y+I) * H2  - (b + sigma + c) * H2)
dY = expression(sigma * H1  - gamma * (X+H1+
     H2+Y+I) * Y  - (b + alpha + c) * Y)
dI = expression(sigma * H2 - gamma * (X+H1+
     H2+Y+I) * I  - (b + c) * I)

j11 = D(dX, "X");  j12 = D(dX, "H1");  j13 = D(dX, 
   "H2");  j14 = D(dX, "Y");  j15 = D(dX, "I")
j21 = D(dH1, "X");  j22 = D(dH1, "H1");  j23 = D(dH1, 
   "H2");  j24 = D(dH1, "Y");  j25 = D(dH1, "I")
j31 = D(dH2, "X");  j32 = D(dH2, "H1");  j33 = D(dH2, 
   "H2");  j34 = D(dH2, "Y");  j35 = D(dH2, "I")
j41 = D(dY, "X");  j42 = D(dY, "H1");  j43 = D(dY, 
   "H2");  j44 = D(dY, "Y");  j45 = D(dY, "I")
j51 = D(dI, "X");  j52 = D(dI, "H1");  j53 = D(dI, 
   "H2"); j54 = D(dI, "Y");  j55 = D(dI, "I")


###################################################
### code chunk number 16: c9-s5-5 (eval = FALSE)
###################################################
## require(rootSolve)
## paras = c(gamma = 0.0397, b = 0.836, a = 1.34, 
##    sigma = 7.5, alpha = 66.36, beta = 33.25,
##    c = 0, rho = 1)
## equil=runsteady(y=log(c(X=12.69/2, H1=0.1, H2=0.1, Y =
##     0.1, I = 0.1)), times=c(0,1E5), 
##     func=coyne, parms=paras)


###################################################
### code chunk number 17: c9-knitr.Rnw:312-313
###################################################
load("c9equil.rda")


###################################################
### code chunk number 18: c9-s5-6
###################################################
#Evaluate Jacobian elements
JJ = with(as.list(c(exp(equil$y), paras)),
c(eval(j11),eval(j12),eval(j13),eval(j14), eval(j15),
   eval(j21),eval(j22),eval(j23),eval(j24), eval(j25),
   eval(j31),eval(j32),eval(j33),eval(j34), eval(j35),
   eval(j41),eval(j42),eval(j43),eval(j44), eval(j45),
   eval(j51),eval(j52),eval(j53),eval(j54), eval(j55)))
#Populate the Jacobian matrix
J = matrix(JJ, nrow = 5, byrow = TRUE)
#Eigen decomposition
which.max(Re(eigen(J, only.values = TRUE)$values))
eigen(J, only.values = TRUE)$values[3]


###################################################
### code chunk number 19: c9-s5-7
###################################################
2*pi/Im(eigen(J)$values[3])


###################################################
### code chunk number 20: c9-s6-1
###################################################
N = 1
gamma = 7/3.8 
omega = 1/(52*4)
mu = 1/(52*70)
R0 = 2.9 


###################################################
### code chunk number 21: c9-s6-2
###################################################
#R0 = beta/(gamma+mu)
beta = R0*(gamma+mu)
paras = c(beta = beta, gamma = gamma, mu = mu, omega = omega)


###################################################
### code chunk number 22: c9-s6-3
###################################################
A = (omega+mu+gamma)/((omega+mu)*(beta-gamma-mu))
GI = 1/(gamma+mu)
GR = 1/(omega+mu)
T = 4*pi/sqrt(4*(R0-1)/(GI*GR)-((1/GR)-(1/A))^2)
T


###################################################
### code chunk number 23: c9-s6-4
###################################################
Sstar = 1/R0
Istar = mu*(1-1/R0)/(gamma+mu-(omega*gamma)/(omega+mu))
Rstar = gamma*Istar/(omega+mu)
eq = list(S = Sstar, I = Istar, R = Rstar)
eq


###################################################
### code chunk number 24: c9-s6-5
###################################################
F = expression(mu * (1-S)  - beta * S * I / N +
    omega * R)
G = expression(beta * S * I / N - (mu + gamma) * I)
H = expression(gamma*I -(mu +omega)*R)
j11 = D(F, "S");j12 = D(F, "I");j13 = D(F, "R")
j21 = D(G, "S");j22 = D(G, "I");j23 = D(G, "R")
j31 = D(H, "S");j32 = D(H, "I");j33 = D(H, "R")

J = with(eq, matrix(c(eval(j11),eval(j12),eval(j13), 
    eval(j21), eval(j22), eval(j23), eval(j31), 
    eval(j32), eval(j33)), nrow = 3, byrow = TRUE))


###################################################
### code chunk number 25: c9-s6-6
###################################################
round(eigen(J)$values, 4)
2*pi/Im(eigen(J)$values)[1]


###################################################
### code chunk number 26: c9-s7-ss1-1
###################################################
parms  = c(mu = 1/(50*52), N = 1, 
     beta =  2.5, gamma = 1/2)
Istar = parms["mu"]*(parms["beta"]/(parms["gamma"]+
     parms["mu"])-1)/parms["beta"]
Sstar = (parms["gamma"]+parms["mu"])/parms["beta"]


###################################################
### code chunk number 27: c9-s7-ss1-2
###################################################
dS = expression(mu * (N  - S)  - beta * S * I / N)
dI = expression(beta * S * I / N - (mu + gamma) * I)
j11 = D(dS, "S")
j12 = D(dS, "I")
j21 = D(dI, "S")
j22 = D(dI, "I")


###################################################
### code chunk number 28: cc9-s7-ss1-3
###################################################
vals = list(mu = parms["mu"], N = parms["N"], beta =  
   parms["beta"], gamma = parms["gamma"], 
   S = Sstar, I = Istar)
J = with(vals,
   matrix(c(eval(j11), eval(j12), 
   eval(j21), eval(j22)),  ncol = 2, byrow = TRUE))


###################################################
### code chunk number 29: c9-s7-ss1-4
###################################################
a1 = D(dS, "beta")
a2 = D(dI, "beta")
A = with(vals, 
   matrix(c(eval(a1), eval(a2)), ncol = 1))
Id = matrix(c(1, 0, 0, 1), ncol = 2)


###################################################
### code chunk number 30: c9-s7-ss1-5
###################################################
wseq = seq(0, pi, length = 500)
Fr = vector("list", 500)  #set up empty list of matrices
#Loop to fill  matrices for each frequency
for(i in 1:500){ 
   #Solve gives matrix inverse
   Fr[[i]] = matrix(solve(Id*1i*wseq[i]-J)%*%A,ncol=1) 
}


###################################################
### code chunk number 31: cc9-s7-ss1-6
###################################################
PS = matrix(NA, ncol = 2, nrow = 500, 
     dimnames = list(1:500, c("S","I"))) 
#Power spectra from real and imaginary 
# parts of the Fourier transform
for(i in 1:500){
PS[i,] = sqrt(Re(Fr[[i]])^2+Im(Fr[[i]])^2)
}
plot(wseq, PS[,2], type = "l", log = "x", 
     xlab = "Frequency (in radians)", ylab = "Amplitude")
#Calculate the dominant period in weeks
2*pi/wseq[which.max(PS[,2])] 


###################################################
### code chunk number 32: c9-s7-ss2-1
###################################################
Seq = expression(S-beta*S*I^alpha/N+B)
Ieq = expression(beta*S*I^alpha/N)
j11 = D(Seq, "S");  j12 = D(Seq, "I")
j21 = D(Ieq, "S");  j22 = D(Ieq, "I")
jj = c(j11, j12, j21,j22)
a1 = D(Seq, "beta")
a2 = D(Ieq, "beta")
aa = c(a1, a2)


###################################################
### code chunk number 33: c9-s7-ss2-2
###################################################
paras  = c(B = 800, beta =  5, alpha = 0.97, N = 1E6)
eqs = sapply(c(quote(B^(1-alpha)*N/beta), quote(B)), 
     eval, as.list(paras))


###################################################
### code chunk number 34: c9-s7-ss2-3
###################################################
J = matrix(sapply(jj, eval, as.list(c(paras, c(S = eqs[1], 
     I = eqs[2])))), 2, byrow = TRUE)
evs = eigen(J)$values
2*pi/atan2(Im(evs[1]), Re(evs[1]))


###################################################
### code chunk number 35: c9-s7-ss2-4
###################################################
A = matrix(sapply(aa, eval, as.list(c(paras, 
   c(S = eqs[1], I = eqs[2])))), 2, byrow = TRUE)
Id = matrix(c(1,0,0,1),ncol = 2)
wseq = seq(0,pi,length = 500)
Fr = vector("list",500)  #Set up empty list of matrices
#Loop to fill those matrices with fourier transforms
for(i in 1:500){ 
   #Solve gives matrix inverse
   Fr[[i]] = matrix(solve(Id-exp(1i*wseq[i])*J)%*%
       A,ncol = 1)  
}

PS = matrix(NA,ncol = 2,nrow = 500,dimnames = list(1:500, 
     c("S","I"))) 
#Power spectra from real and imaginary parts
for(i in 1:500){
     PS[i,] = sqrt(Re(Fr[[i]])^2+Im(Fr[[i]])^2)
}


###################################################
### code chunk number 36: c9-knitr.Rnw:542-559
###################################################
SimTsir=function(alpha=0.97, B=2300, beta=25, sdbeta=0,
    S0 = 0.06, I0=180, IT=520, N=3.3E6){
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean=beta, sd=sdbeta) * I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i]=0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
    list(I = I, S = S)
}
load("tsirsim.rda")


###################################################
### code chunk number 37: c9-s7-ss2-5 (eval = FALSE)
###################################################
## out = SimTsir(B = 800, beta = 5, sdbeta = 1, N = 1E6, 
##      IT = 100*52, I0 = 10, S0 = 0.3)
## plot(out$I[1:1040], xlab = "Biweek", ylab = "Incidence", 
##      type = "l")


###################################################
### code chunk number 38: c9-s7-ss2-6
###################################################
sfit=spectrum(out$I[-c(1:104)])


###################################################
### code chunk number 39: c9-s7-ss2-7
###################################################
require(polspline)
sfit2 = lspec(out$I[-c(1:104)])
plot(wseq, PS[,2], type = "l", ylab = "Amplitude", 
   xlab = "Frequency (in radians)", xlim = c(0,0.6))
lines(pi*sfit$freq/0.5, 
   5000*sfit$spec/max(sfit$spec), col = 2)
par(new = TRUE)
plot(sfit2, col = 3, xlim = c(0,0.6), axes = FALSE)
legend("topright", c("Transfer fn", "Periodogram", 
   "Log-spline"), lty = c(1,1,1), col = c(1,2,3))


