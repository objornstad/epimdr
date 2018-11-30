### R code from Chapter 12
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866


###################################################
### code chunk number 2: c12-s1-1
###################################################
foo = function(x){
res = x
class(res) = "foo"
return(res)}

print.foo = function(obj){
cat("foo is:\n", obj)}

summary.foo = function(obj){
cat("In summary, foo is:\n", obj)}

plot.foo = function(obj){
plot(NA, type = 'n', ylim = c(0,1), xlim = c(0,1), ylab = "")
text(x = seq(0.1, 0.9, by = 0.1), y = seq(0.1, 0.9, by = 0.1), as.character(obj))}


###################################################
### code chunk number 3: c12-s1-2
###################################################
zz = foo("pibble")


###################################################
### code chunk number 4: c12-s1-3
###################################################
zz


###################################################
### code chunk number 5: c12-s1-4
###################################################
summary(zz)


###################################################
### code chunk number 6: c12-s1-5
###################################################
plot(zz)


###################################################
### code chunk number 7: c12-s3-1
###################################################
ringlattice = function(N,K){
# N is the number of nodes
# K is the number of neighbors on each side to which
#  each node is connected so degree  =  2xK
CM = toeplitz(c(0,rep(1,K),rep(0,N-2*K-1),rep(1,K)) )
    class(CM) = "cm"
    return(CM)
}


###################################################
### code chunk number 8: c12-s3-2
###################################################
plot.cm = function(CM){
   N = dim(CM)[1]
   theta = seq(0,2*pi,length = N+1)
   x = cos(theta[1:N])
   y = sin(theta[1:N])
   symbols(x,y, fg = 0, circles = rep(1, N), 
      inches = 0.1, bg = 1, xlab = "", ylab = "")
   segx1 = as.vector(matrix(x, ncol = length(x),
       nrow = length(x), byrow = TRUE))
   segx2 = as.vector(matrix(x, ncol = length(x), 
       nrow = length(x), byrow = FALSE))
   segy1 = as.vector(matrix(y, ncol = length(x), 
       nrow = length(x), byrow = TRUE))
   segy2 = as.vector(matrix(y, ncol = length(x), 
       nrow = length(x), byrow = FALSE))
   segments(segx1,segy1, segx2, 
      segy2, lty = as.vector(CM))
}


###################################################
### code chunk number 9: c12-s3-3
###################################################
cm = ringlattice(N = 20,K = 4)
plot(cm)


###################################################
### code chunk number 10: c12-s3-4
###################################################
WattsStrogatz = function(N, K, Prw){
  # Build a Watts-Strogatz contact matrix from 
  # a ring lattice, Prw is the rewiring probability
  CM = ringlattice(N = N, K = K)
  CMWS = CM
  tri = CM[upper.tri(CM)]
  Br = rbinom(length(tri),1,Prw) # Break edges 
  a = 0
  for(i in 1:(N-1)){									
    for(j in (i+1):N){
         a = a+1								
         if(Br[a]==1 & CMWS[i,j]==1){ #If "Br == 1"
             CMWS[i,j] = CMWS[j,i] = 0 # break edge
             tmp = i
             tmp2 = c(i, which(CMWS[i,]==1))
             #new edge, if already present try again
             while(any(tmp2==tmp)){
                tmp = ceiling(N*runif(1))} 
             CMWS[i,tmp] = CMWS[tmp,i] = 1 # make new edge
             }
    }
  }
class(CMWS) = "cm"
return(CMWS)
}


###################################################
### code chunk number 11: c12-s3-5
###################################################
cm2 = WattsStrogatz(N = 20, K = 4, Prw = .3)
plot(cm2)


###################################################
### code chunk number 12: c12-s3-6
###################################################
summary.cm = function(x, plot = FALSE){
	x = table(apply(x, 2, sum))
	res = data.frame(n = x)
	names(res) = c("degree", "freq")
	if(plot) barplot(x, xlab = "degree")
    return(res)
}
summary(cm2, plot = TRUE)


###################################################
### code chunk number 13: c12-s3-7
###################################################
BarabasiAlbert = function(N, K){
   CM = matrix(0, ncol = N, nrow = N)
   CM[1,2] = 1
   CM[2,1] = 1
   for(i in 3:N){	
       probs = apply(CM, 1, sum)
       link = unique(sample(c(1:N)[-i], 
          size = min(c(K, i-1)), prob = probs[-i]))
      CM[i, link] = CM[link, i] = 1
   }
class(CM) = "cm"
return(CM)
}


###################################################
### code chunk number 14: c12-s3-8
###################################################
cm3=BarabasiAlbert(200, 4)
ed = summary(cm3)
plot(as.numeric(ed$degree), ed$freq, log = "xy", xlab = "Degree", ylab = "Frequency")


###################################################
### code chunk number 15: c12-s3-9
###################################################
require(statnet)
plot(network(cm3, directed = FALSE))


###################################################
### code chunk number 16: c12-s4-1
###################################################
NetworkSIR = function(CM,tau,gamma){
  #generate SIR epidemic on a CM-network  
  #CM  =  contact matrix
  #tau  =  probability of infection across an edge
  #gamma  =  probability of removal per time step 
  N = dim(CM)[1]
  I = matrix(rep(0,N),nrow = N,ncol = 1) #First infecteds 
  S = matrix(rep(1,N),nrow = N,ncol = 1) #First susceptibles 
  R = matrix(rep(0,N),nrow = N,ncol = 1) #First removed 
  I1 = sample(1:N, size = 1)#Pick first random infected
  I[I1,1] = 1
  S[I1,1] = 0	
  t = 1
  while(sum(I[,t-1])>0 | t==1){
    t = t+1
    infneigh = CM%*%I[,t-1]
    pinf = 1-(1-tau)^infneigh
    newI = rbinom(N, S[,t-1], pinf)
    newR = rbinom(N, I[,t-1], gamma)
    nextS = S[,t-1]-newI
    nextI = I[,t-1]+newI-newR
    nextR = R[,t-1]+newR
    I = cbind(I, nextI)
    S = cbind(S, nextS)
    R = cbind(R, nextR)
  }	
res = list(I = I,S = S,R = R)
class(res) = "netSIR"
return(res)
}


###################################################
### code chunk number 17: c12-s4-2
###################################################
summary.netSIR = function(x){
   t = dim(x$S)[2]
   S = apply(x$S,2,sum)
   I = apply(x$I,2,sum)
   R = apply(x$R,2,sum)
   res = data.frame(S = S,I = I,R = R)
   return(res)
}

plot.netSIR = function(x){
   y = summary(x)	
   plot(y$S, type = "b", xlab = "time", ylab = "")
   lines(y$I, type = "b", col = "red")	
   lines(y$R, type = "b", col = "blue")	
   legend("right", legend = c("S", "I", "R"),
     lty = c(1,1,1), pch = c(1,1,1),
     col = c("black", "red", "blue"))
}


###################################################
### code chunk number 18: c12-s4-3
###################################################
cm1 = BarabasiAlbert(N = 200,K = 2)  #(i)
cm2 = WattsStrogatz(N = 200, K = 2, Prw = .1) #(ii)
cm3 = WattsStrogatz(N = 200, K = 2, Prw = 1)  #(iii)
cm4 = ringlattice(N = 200,K = 2)   #(iv)
sim1 = NetworkSIR(cm1,.3,0.1)
sim2 = NetworkSIR(cm2,.3,0.1)
sim3 = NetworkSIR(cm3,.3,0.1)
sim4 = NetworkSIR(cm4,.3,0.1)
plot(apply(sim1$I,2,sum), type = "l", xlab = "Time", 
     ylab = "Infected")
lines(apply(sim2$I,2,sum), type = "l", col = "red")
lines(apply(sim3$I,2,sum), type = "l", col = "red", lty = 2)
lines(apply(sim4$I,2,sum), type = "l", col = "blue")
legend("topright", legend = c("Scale-free", "WS(0.1)", 
   "Poisson", "Lattice"), lty = c(1,1, 2, 1), 
   col=c("black", "red", "red", "blue"))


###################################################
### code chunk number 19: c12-s4-4
###################################################
r0fun = function(CM, tau, gamma){
x = apply(CM, 2, sum)
(tau/(tau+gamma))*(mean(x^2)-(mean(x)))/mean(x)
}
r0fun(cm1, 0.3, 0.1)
r0fun(cm2, 0.3, 0.1)
r0fun(cm3, 0.3, 0.1)
r0fun(cm4, 0.3, 0.1)


###################################################
### code chunk number 20: c12-s4-5
###################################################
data(gonnet)
nwt = network(gonnet, directed = TRUE)
x = degree(nwt)[2:89]
mean(x)


###################################################
### code chunk number 21: c12-s4-6
###################################################
(mean(x^2) - (mean(x)))/mean(x)


###################################################
### code chunk number 22: c12-s4-7
###################################################
#Undirected network
cmg = gonnet+t(gonnet)
#Simulate epidemid
cep = NetworkSIR(cmg, 0.3, 0.1)
sm = summary(cep)
par(mfrow = c(1,2))
inf = ifelse(apply(cep$I,1,sum)>0,2,1)
nwt = network(cmg, directed = FALSE)
plot(nwt, vertex.col = inf)
matplot(sm, ylab = "Numbers")
legend("right", c("S", "I", "R"), 
     pch = c("1", "2", "3"), col = c(1,2,3))


