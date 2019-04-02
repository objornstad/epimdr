### R code from Chapter 14
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866

###################################################
### code chunk number 1: c14-knitr-springer.Rnw:2-12
###################################################
library(knitr)
opts_chunk$set(
width=6,fig.height=4
)
opts_chunk$set(
tidy=TRUE
)
opts_chunk$set(
fig.keep='none'
)


###################################################
### code chunk number 2: c14-knitr-springer.Rnw:15-16
###################################################
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)


###################################################
### code chunk number 3: c14-s1-1
###################################################
data(burnett)
plot(burnett$Generation, 
   burnett$NumberofHostsParasitized, type="b", 
   ylab="Numbers", xlab="Generation")
lines(burnett$Generation, 
   burnett$NumberofHostsUnparasitized, type="b", 
   col=2, pch=2)
legend("topleft", legend=c("Parasitoid", "Host"), 
   lty=c(1,1), pch=c(1,2), col=c(1,2))


###################################################
### code chunk number 4: c14-s1-2
###################################################
NB = function(R, a, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   H=rep(NA, T) #Host series
   P=rep(NA, T) #Parasitoid series
   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the parasitoid series

   for(t in 2:T){
     H[t] = R * H[t-1] * exp(- a * P[t-1])
     P[t] = R * H[t-1] * (1-exp(- a * P[t-1]))
   } 

   res= list(H = H, P = P) 
   return(res)
} 


###################################################
### code chunk number 5: c14-s1-3
###################################################
sim = NB(R = 1.1, a = 0.1)
time = 1:100
par(mfrow=c(1,2))
plot(time, sim$H, type =  "l", xlab  =  "Generations", 
     ylab  =  "Host abundance", ylim = c(0,14))
points(time, sim$P, type = "l", col = "red")
plot(sim$H,sim$P, type = "l", xlab = "Host abundance", 
     ylab = "Parasitoid abundance")


###################################################
### code chunk number 6: c14-s1-4
###################################################
aVals = seq(0,1,by = 0.01)
tte = rep(NA,length(aVals))
for (i in c(1:length(aVals))){
     sim =  NB(R = 1.1, a = aVals[i], T = 500)
     tte[i] = min(which(sim$P==0))
     }
plot(aVals,tte,type = "b", ylab="TTE", 
     xlab="Search efficiency")


###################################################
### code chunk number 7: c14-s1-5
###################################################
ssfn=function(par){
  R=exp(par[1])
  a=exp(par[2])
  sim=NB(R,a, T=22, H0=10.1, P0=11.9)
  ss=sum((burnett$NumberofHostsUnparasitized-sim$H)^2+
      (burnett$NumberofHostsParasitized-sim$P)^2)
 return(ss)
 }
par=log(c(2, 0.05)) 
fit=optim(par, ssfn)
exp(fit$par)


###################################################
### code chunk number 8: c14-s1-6
###################################################
sim = NB(R=2.16767, a=0.06812, T=22, H0=10.1, P0=11.9)
plot(burnett$Generation, 
   burnett$NumberofHostsParasitized, type="b", 
   ylab="Numbers", xlab="Generation")
lines(burnett$Generation, sim$P)
lines(burnett$Generation, 
   burnett$NumberofHostsUnparasitized, type="b", 
   col=2, pch=2)
lines(burnett$Generation, sim$H, col=2)
legend("topleft", legend=c("Parasitoid", "NB-P", 
   "Host", "NB-H"), lty=c(1,1,1,1), pch=c(1,NA,2,NA),
   col=c(1,1,2,2))


###################################################
### code chunk number 9: c14-s2-1
###################################################
F = expression(R*H*exp(-a * P))
G = expression(R*H*(1-exp(-a * P)))
j11 = D(F, "H"); j12 = D(F, "P") 
j21 = D(G, "H"); j22 = D(G, "P") 
R = 2.17; a = 0.068
params = c(R = R, a = a, P = log(R)/a, H =  log(R)/(a*(R-1)))
J = with(as.list(params),
    matrix(c(eval(j11), eval(j12), eval(j21), 
    eval(j22)), ncol = 2, byrow = TRUE))
eigen(J, only.values = TRUE)$values
max(abs(eigen(J)$values))


###################################################
### code chunk number 10: c14-s2-2
###################################################
2*pi/atan2(Im(eigen(J)$values[1]), 
     Re(eigen(J)$values[1]))


###################################################
### code chunk number 11: c14-s2-3
###################################################
RVals = seq(1.1, 3, by = 0.1)
per = rep(NA, length(RVals))
for(i in 1:length(RVals)){
R = RVals[i]
a = 0.068
params = c(R = R, a = a, P = log(R)/a,
     H =  log(R)/(a*(R-1)))
J = with(as.list(params),
    matrix(c(eval(j11), eval(j12), eval(j21), 
    eval(j22)), ncol = 2, byrow = TRUE))
per[i] = 2*pi/atan2(Im(eigen(J)$values[1]), 
     Re(eigen(J)$values[1]))
}
plot(RVals, per, type = "b", xlab = "R", ylab = "Period")


###################################################
### code chunk number 12: c14-s5-1
###################################################
#Dh is proportion of hosts that disperses 
#Dp is proportion of parasitoids that disperses
Dh = 0.5
Dp = 0.7
#xlen is width of the lattice (E-W)
#ylen is height of the lattice (N-S)
xlen = 30
ylen = 30


###################################################
### code chunk number 13: c14-s5-2
###################################################
hp.dyn = function(h, p, R, a){ 
   #hnew is the post-interaction host density
   hnew = R * h * exp(- a * p)
   #pnew is the post-interaction parasitoid density
   pnew = R * h * (1 - exp(- a * p))
   #the two vectors of results are stored in a "list"
   res = list(h = hnew, p = pnew)
   return(res)
} 


###################################################
### code chunk number 14: c14-s5-3
###################################################
xy = expand.grid(1:xlen, 1:ylen)
dmat = as.matrix(dist(xy))


###################################################
### code chunk number 15: c14-s5-4
###################################################
kh = ifelse(dmat<1.5,Dh/8,0)
kp = ifelse(dmat<1.5,Dp/8,0)
diag(kh) = 1-Dh
diag(kp) = 1-Dp


###################################################
### code chunk number 16: c14-s5-5
###################################################
IT = 600
hmat = matrix(NA, nrow=xlen*ylen, ncol = IT)
pmat = matrix(NA, nrow=xlen*ylen, ncol = IT)
hmat[,1] = 0
pmat[,1] = 0
hmat[23,1] = 4
pmat[23,1] = 1


###################################################
### code chunk number 17: c14-s5-6
###################################################
for(i in 2:IT){
   #growth
   tmp = hp.dyn(h = hmat[,(i-1)], p = pmat[,(i-1)], 
      R = 2, a = 1)
   #redistribution
   hmat[,i] = tmp$h%*%kh;
   pmat[,i] = tmp$p%*%kp;
   cat(i, " of ", IT, "\r")
}


###################################################
### code chunk number 18: c14-s5-7 (eval = FALSE)
###################################################
#plot the last 100 generations for the parasitoid
for(i in 1:100){
   x=xy[,1]
   y=xy[,2]
   z=pmat[,i+500]
   symbols(x,y, fg=2, circles=z, inches=0.1, 
      bg=2, xlab="", ylab="")
   Sys.sleep(.1) #this is to slow down the plotting
}


###################################################
### code chunk number 19: c14-s5-8 (eval = FALSE)
###################################################
for(k in 1:50){
png(filename=paste("Pplot",k,".png", sep=""))
   x=xy[,1]
   y=xy[,2]
   z=pmat[,k+500]
   symbols(x,y, fg=2, circles=z, inches=0.1, bg=2, xlab="", ylab="")
dev.off()
}
system("convert P*.png -delay 500 -coalesce -layers OptimizeTransparency cml2.gif")
system("rm *.png")


###################################################
### code chunk number 20: c14-knitr-springer.Rnw:316-318 (eval = FALSE)
###################################################
require(shiny)
May.app


