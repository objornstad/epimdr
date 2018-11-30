### R code from Chapter 6
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866


###################################################
### code chunk number 3: c6-s2-1
###################################################
times = seq(0, 100, by=1/52) 
paras  = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, 
     sigma = 365/8, gamma = 365/5) 
xstart = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
out = as.data.frame(ode(xstart, times, seirmod2, paras))

par(mfrow = c(1,2))  
plot(times, out$I, ylab = "Infected", xlab = "Time", 
     xlim = c(90,100), type="l")
acf(out$I, lag.max = 156, main="")


###################################################
### code chunk number 4: c6-s2-2
###################################################
require(forecast)
data(Icelandflu)


###################################################
### code chunk number 5: c6-s2-3
###################################################
ilits = ts(sqrt(Icelandflu$ili), start = c(1980, 1), 
     end = c(2009, 12), frequency = 12)
plot(decompose(ilits))


###################################################
### code chunk number 6: c6-s2-4
###################################################
wts = window(ilits, start = c(1990,6), end = c(2000,5))       
fit = arima(sqrt(wts), order = c(2, 0, 1), 
     list(order = c(1, 0, 0), period = 12))
coef(fit)
fore = predict(fit, n.ahead = 24)
#Calculate approximate upper (U) and 
#lower (L) prediction intervals
U = fore$pred + 2*fore$se 
L = fore$pred - 2*fore$se
# plot observed and predicted values
ts.plot(sqrt(wts), fore$pred, U, L, col = c(1, 2, 4, 4), 
     lty = c(1, 1, 2, 2), ylab = "Sqrt(cases)")
legend("bottomleft", c("ILI", "Forecast", 
     "95% Error Bounds"), col=c(1, 2, 4),lty=c(1, 1, 2))


###################################################
### code chunk number 7: c6-s3-1
###################################################
my.spec = spectrum(out$I, plot = FALSE)
par(mfrow = c(1, 2)) 
#default plot (less default lables)
plot(my.spec,  xlab = "Frequency", ylab = "Log-amplitude", 
     main = "", sub = "")
#plot with period (rather than frequency)
plot(1/my.spec$freq/52, my.spec$spec, type = "b", xlab =
     "Period (year)", ylab = "Amplitude", xlim = c(0, 5))


###################################################
### code chunk number 8: c6-s4-1
###################################################
#Simulate and plot time series
times  = seq(0, 25, by = 1/52)
paras  = c(mu = 1/50, N = 1, beta =  1000, 
     sigma = 365/8, gamma = 365/5)
xstart = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
out2 = as.data.frame(ode(xstart, times, seirmod, paras))
par(mfrow = c(1, 2)) #Side-by-side plots
plot(times, out2$I, type = "l", xlab = "Time", 
     ylab = "Infected")

#Wavelet analysis
require(Rwave)   
#Set the number of "octaves" and "voices"
no = 8; nv = 32
#Calculate periods
a = 2^seq(1,no+1-1/nv, by = 1/nv)
#Do the continous wavelet decomposition
wfit = cwt(out2$I, no, nv, plot = FALSE)
#Calculate the wavelet spectrum
wspec = Mod(wfit)

#Wavelet plot with contours
image(x = times, wspec, col = gray((12:32)/32), y = a/52, 
     ylim = c(0, 4), xlab = "Time", ylab = "Period")
contour(x = times, wspec, y = a/52, ylim = c(0, 4), 
     zlim = c(mean(wspec), max(wspec)), add = TRUE)


###################################################
### code chunk number 9: c6-s3-2
###################################################
plot(a/52, wspec[104,], type = "l", ylab = "Amplitude",
     xlab = "Period")
lines(a/52, wspec[1040,], type = "l", 
     lty = 2, col = "red")
legend("topright", legend = c("Year 2", "Year 10"),
     lty = c(1,2), col = c("black", "red"))


###################################################
### code chunk number 10: c6-s5-1
###################################################
data(meas)
head(meas)
par(mar = c(5,5,2,5)) #Make room for two axes
plot(meas$time, meas$London, type="b", xlab="Week", 
     ylab="Incidence", ylim=c(0,8000))
par(new=T) #Superimposed births plot
plot(meas$time, meas$B, type="l", col="red", 
     axes=FALSE, xlab=NA, ylab=NA, ylim=c(1000, 2700))
axis(side = 4)
mtext(side = 4, line = 3, "Births")
legend("topright", legend=c("Cases", "Births"), 
     lty=c(1,1), col=c("black", "red"))


###################################################
### code chunk number 11: c6-s5-2
###################################################
#Set octaves, voices and associated periods
no = 8; nv = 32
a = 2^seq(1,no+1-1/nv, by = 1/nv)
#Continous wavelet decomposition
wfit = cwt(meas$London, no, nv, plot = FALSE)
wspec = Mod(wfit)
#Crazy climber
crcinc = crc(wspec, nbclimb = 10, bstep = 100)
fcrcinc = cfamily(crcinc, ptile = 0.5, nbchain = 1000,
      bstep = 10)
ridges = fcrcinc[[1]]
ridges[which(ridges == 0)]<-NA
#Wavelet plot with crazy-climber and contours
image(x = meas$time, wspec, col = gray((12:32)/32), y = a/26, 
     ylim = c(0.1, 3), ylab = "Period", xlab = "Year")
contour(x = meas$time, wspec, y = a/26, ylim = c(0, 3), 
     nlevels = 6, zlim = c(mean(wspec), max(wspec)), 
     add = TRUE)
image(x = meas$time, y = a/26, z = ridges, add = TRUE, 
     col = gray(0))


###################################################
### code chunk number 12: c6-s5-3
###################################################
plot(a/26, wspec[261,], type = "l",xlim = c(0, 3), 
     xlab = "period (years)", ylab = "amplitude")
lines(a/26, wspec[27,], type = "l",  lty = 2, col = "red")
legend("topleft", legend = c("1945", "1954"),
      lty = c(2, 1), col = c("red", "black"))


###################################################
### code chunk number 13: c6-s6-1
###################################################
data(tywhooping)
tywhooping$TIME = tywhooping$YEAR + tywhooping$WEEK/52
tywhooping$TM = 1:length(tywhooping$YEAR)
data(tydiphtheria)
data(tymeasles)
tydiphtheria$TIME = tymeasles$TIME = tymeasles$YEAR +
     tymeasles$WEEK/52


###################################################
### code chunk number 14: c6-s7-1
###################################################
data(tywhooping)
whp = na.omit(tywhooping)

#data with missing values interpolated
require(imputeTS)
sum(is.na(tywhooping$PHILADELPHIA))
tywhooping$PHILADELPHIA=
     na.interpolation(ts(tywhooping$PHILADELPHIA))

#Classic periodogram
my.spec = spectrum(sqrt(tywhooping$PHILADELPHIA))
#Lomb periodogram
require(nlts)
my.lomb = spec.lomb(x = whp$TM, y = sqrt(whp$PHILADELPHIA))

plot(1/my.spec$freq/52, my.spec$spec, type = "b",
     xlab = "Period (year)", ylab = "Amplitude")
par(new = TRUE)
plot(1/my.lomb$freq/52, my.lomb$spec, axes = FALSE, 
     type = "b", col = 2, xlab = "", ylab = "")
legend("topright", legend = c("Classic", "Lomb"),
     lty = c(1, 1), pch = c(1, 1), col = c("black", "red"))


###################################################
### code chunk number 15: c6-s8-1
###################################################
data(tymeasles)
sum(is.na(tymeasles$PHILADELPHIA))
tymeasles$PHILADELPHIA =
     na.interpolation(ts(tymeasles$PHILADELPHIA))


###################################################
### code chunk number 16: c6-s8-2
###################################################
par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
layout(matrix(c(1, 1, 2, 2, 2), ncol = 1))
plot(tymeasles$TIME, sqrt(tymeasles$PHILADELPHIA), 
     type = "b", ylab = "Sqrt(incidence)")
title("Measles 1914-47")

no = 8; nv = 16; a = 2^seq(1,no+1-1/nv, by = 1/nv)
wfit = cwt(sqrt(tymeasles$PHILADELPHIA), 
     no, nv, plot = FALSE)
wspec = Mod(wfit)
par(mar = c(1, 4, 0.25, 1))
image(z = wspec, y = a/52, ylim = c(0, 4), ylab = "Period(year)", 
     col = gray((12:32)/32), xaxt = 'n')
contour(z = wspec, y = a/52, ylim = c(0, 4), nlevels = 6,
     zlim = c(mean(wspec), max(wspec)), add = TRUE)


###################################################
### code chunk number 17: c6-s8-3
###################################################
plot(a/52,wspec[54,], type = "l", xlim = c(0, 4),
     xlab = "Period (years)", ylab = "Amplitude", 
     col = "red", lty = 2)
lines(a/52, wspec[1357,], type = "l", xlim = c(0, 4))
legend("topleft", legen d= c("1915", "1940"),
     lty = c(2, 1), col = c("red","black"))


###################################################
### code chunk number 18: c6-s9-1
###################################################
data(tydiphtheria)
sum(is.na(tydiphtheria$PHILADELPHIA))
tydiphtheria$PHILADELPHIA =
     na.interpolation(ts(tydiphtheria$PHILADELPHIA))

par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
layout(matrix(c(1, 1, 2, 2, 2), ncol = 1))
plot(tydiphtheria$TIME, sqrt(tydiphtheria$PHILADELPHIA), 
     type = "b", ylab = "Sqrt(incidence)")
title("Diphteria 1914-47")

no = 8; nv = 16; a = 2^seq(1,no+1-1/nv, by = 1/nv)
wfit = cwt(sqrt(tydiphtheria$PHILADELPHIA), 
     no, nv, plot = FALSE)
wspec = Mod(wfit)
par(mar = c(1, 4, 0.25, 1))
image(z = wspec, y = a/52, ylim = c(0, 3), ylab = "Period(year)", 
     col = gray((12:32)/32), xaxt = 'n')
contour(z = wspec, y = a/52, ylim = c(0, 3), nlevels = 6,
     zlim = c(mean(wspec), max(wspec)), add = TRUE)


###################################################
### code chunk number 19: c6-s8-2
###################################################
#midpass filter
sel=a>45 & a<60
rec=0.6*apply(Re(wfit[,sel])/sqrt(a[sel]), 1, 
     sum)/(0.776*(pi^(-1/4)))
data=pi*scale(sqrt(tydiphtheria$PHILADELPHIA))/2
plot(tydiphtheria$TIME, data, type="b", xlab="Year", 
     ylab="Scaled cases")
lines(tydiphtheria$TIME, rec, type="l", col=2, lwd=3)
legend("topright", legend=c("Scaled cases",
     "Annual reconstruction"), pch=c(1, NA), lty=c(1,1), 
     lwd=c(1,3), col=c("black", "red"))


###################################################
### code chunk number 20: c6-s10-1 (eval = FALSE)
###################################################
#fft
x = meas$London
p = length(x)
z = fft(x)
f = seq(from=0,length=p,by=1/p)
a = Re(z)
b = Im(z)
#reconstruction
rec2=matrix(NA, ncol=p, nrow=p)
for (k in 1:p) {
  rec2[,k] <- (a*cos(2*pi*(k-1)*f) - b*sin(2*pi*(k-1)*f))/p
}


###################################################
### code chunk number 21: c6-s10-2 (eval = FALSE)
###################################################
sim = rep(0, p)
n = 0
samp = order(a^2+b^2, decreasing = TRUE)
for(g in samp){
n = n+1
par(mfrow = c(1, 2))
plot(x, ylim = c(0, 11000), ylab = "Incidence", 
   xlab = "Biweek")
title(paste("nfreq = ", n))
sim = sim + rec2[g,]
lines(sim, col=2)
par(new = TRUE)
sc = scale((cos(2*pi*(0:(p-1))*f[g]) - 
   sin(2*pi*(0:(p-1))*f[g]))/p)
plot(sc*(a^2+b^2)[g]/max(a^2+b^2), type = "l", col =
  gray(.5), ylim = c(-8, 2), axes = FALSE, xlab = "", ylab = "")
plot(x, sim, ylab = "Reconstructed", xlab = "Observed", 
  ylim = c(0, 8000))
#Sys.sleep makes R wait a bit
Sys.sleep(.2)
}


