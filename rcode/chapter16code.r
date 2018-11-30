### R code from Chapter 16
### Epidemics: Models and Data in R
### Ottar N. Bjornstad (ISBN 978-3-319-97487-3) https://www.springer.com/gp/book/9783319974866

###################################################
### code chunk number 3: c16-s2-1
###################################################
data(fiv)
Day31=fiv[fiv$Day==31,]
dimnames(Day31)[[1]]=Day31$Id
Day31=na.omit(Day31[, -c(1, 14,15,16)])
Day59=fiv[fiv$Day==59,]
dimnames(Day59)[[1]]=Day59$Id
Day59=na.omit(Day59[, -c(1, 14,15,16)])


###################################################
### code chunk number 4: c16-s2-2
###################################################
data(SH9)
SH9RBC=SH9[,-c(1,3,4,7,8,10,11)]


###################################################
### code chunk number 5: c16-s2-3
###################################################
SH9RBCw = reshape(SH9RBC, idvar = "Ind2",
     direction = "wide", timevar = "Day")
SH9RBCw = SH9RBCw[,-seq(4,50,by = 2)]
names(SH9RBCw)[2] = "Treatment"


###################################################
### code chunk number 6: c16-s3-1
###################################################
require(ade4)
pca31=dudi.pca(Day31[,1:11], scannf = FALSE, nf = 5)
#select 5 axes
groups = Day31$Treatment
s.arrow(dfxy = pca31$co[,1:2]*8, ylim = c(-7,9), 
    sub = "Day 31", possub = "topleft", csub = 2)
s.class(dfxy = pca31$li[,1:2], fac = groups, cellipse = 2, 
    axesell = FALSE, cstar = 0 , col = c(2:5), add.plot = TRUE)
add.scatter.eig(pca31$eig, xax = 1, yax = 2, 
    posi = "bottomright")


###################################################
### code chunk number 7: c16-s3-2
###################################################
pca59 = dudi.pca(Day59[,1:11], scannf = FALSE, nf = 5)
groups = Day59$Treatment
s.arrow(dfxy = pca59$co[,1:2]*8, ylim = c(-7,9), 
    sub = "Day 59",  possub = "topleft", csub = 2)
s.class(dfxy = pca59$li[,1:2], fac = groups, cellipse = 2, 
    axesell = FALSE, cstar = 0 , col = c(2:5), add.plot = TRUE)
add.scatter.eig(pca59$eig, xax = 1, yax = 2, 
    posi = "bottomright")


###################################################
### code chunk number 8: c16-s4-1
###################################################
require(MASS)
Day31sc = Day31
Day31sc[,1:11] = apply(Day31[,1:11],2,scale)


###################################################
### code chunk number 9: c16-s4-2
###################################################
lda31 = lda(Treatment ~ CD4 + CD8B + CD25 + FAS + 
     IFNg + IL_10 + IL_12 + lymphocyte + neutrophils +
     TNF_a, data = Day31sc)
plot(lda31)


###################################################
### code chunk number 10: c16-s4-3
###################################################
pr = predict(lda31, method = "plug-in")$class
table(pr, Day31sc$Treatment)


###################################################
### code chunk number 11: c16-s4-4
###################################################
ld1 = as.matrix(Day31sc[,attr(lda31$terms,
    "term.labels")])%*%matrix(lda31$scaling[,1], ncol = 1)
ld2 = as.matrix(Day31sc[,attr(lda31$terms,
    "term.labels")])%*%matrix(lda31$scaling[,2], ncol = 1)
groups = Day31$Treatment

contribs = lda31$svd/sum(lda31$svd)
s.arrow(dfxy = lda31$scaling[,1:2], sub = "Day 31", 
     possub = "topleft", csub = 2)
s.class(dfxy = cbind(ld1, ld2)*2.5, fac = groups, 
    cellipse = 2,  axesell = FALSE, cstar = 0, 
    col = c(2:5), add.plot = TRUE)
add.scatter.eig(contribs, xax = 1, yax = 2, 
    posi = "bottomright")


###################################################
### code chunk number 12: c16-s4-5
###################################################
Day59sc = Day59
Day59sc[,1:11] = apply(Day59[,1:11],2,scale)
lda59  =  lda(Treatment ~ CD4 + CD8B + CD25 + FAS + 
   IFNg + IL_10 + IL_12 + lymphocyte + neutrophils + 
   TNF_a, data = Day59sc)
pr = predict(lda59, method = "plug-in")$class
table(pr, Day59sc$Treatment)

ld1 = as.matrix(Day59sc[,attr(lda59$terms,
   "term.labels" )])%*%matrix(lda59$scaling[,1], ncol = 1)
ld2 = as.matrix(Day59sc[,attr(lda59$terms,
   "term.labels" )])%*%matrix(lda59$scaling[,2], ncol = 1)
groups = Day59$Treatment

contribs  =  lda59$svd/sum(lda59$svd)
s.arrow(dfxy = lda59$scaling[,1:2], sub = "Day 59", 
   possub = "topleft", csub = 2)
s.class(dfxy = cbind(ld1, ld2), fac = groups, cellipse = 2, 
   axesell = FALSE, cstar = 0 , col = c(2:5), add.plot = TRUE)
add.scatter.eig(contribs, xax = 1, yax = 2, 
   posi = "bottomright")


###################################################
### code chunk number 13: c16-s5-1
###################################################
options(width=50)
Y = cbind(Day59sc$CD4, Day59sc$CD8B, Day59sc$CD25, 
   Day59sc$FAS, Day59sc$IFNg, Day59sc$IL_10, 
   Day59sc$IL_12, Day59sc$lymphocyte, 
   Day59sc$neutrophils, Day59sc$TNF_a)
X = Day59$Treatment
mova59 = manova(Y~X)
summary(mova59, test = "Pillai")


###################################################
### code chunk number 14: c16-s6-1
###################################################
require(ade4)
dead = ifelse(SH9RBCw[,27]==0, "dead", "alive")
pcaRBC = dudi.pca(SH9RBCw[,3:27], scale = FALSE, 
   scannf = FALSE, nf = 5)
s.arrow(dfxy = pcaRBC$co[,1:2]*3, xlim = c(-10, 10), 
   ylim = c(-5,5), sub = "RBC", possub = "topleft", csub = 2)
s.class(dfxy = pcaRBC$li[,1:2]*.3, fac = as.factor(dead), 
   cellipse = 2, axesell = FALSE, cstar = 0 , 
   col = c(2:7), add.plot = TRUE)
add.scatter.eig(pcaRBC$eig, xax = 1, yax = 2, 
   posi = "bottomright")


###################################################
### code chunk number 15: c16-s6-2
###################################################
SH9RBCw2 = SH9RBCw[dead=="alive",]
groups = SH9RBCw2$Treatment
pcaRBC = dudi.pca(SH9RBCw2[,3:27], scale = FALSE, scannf  =  
   FALSE, nf  =  5)
s.arrow(dfxy = pcaRBC$co[,1:2]*3, xlim = c(-4,9), 
   ylim = c(-5,5), sub = "RBC", possub = "topleft", csub = 2)
s.class(dfxy = pcaRBC$li[,1:2]*.3, fac = groups, cellipse = 2, 
   axesell = FALSE, cstar = 0 , col = c(2:7), add.plot = TRUE)
add.scatter.eig(pcaRBC$eig, xax = 1, yax = 2, 
   posi = "bottomright")


###################################################
### code chunk number 16: c16-s7-1
###################################################
par(mfrow = c(1,2))
#Gets the experimental days
day = unique(SH9$Day)
#Calculate the average time series
avg = apply(SH9RBCw2[,3:27], 2, mean)
plot(day, avg, type = "b", ylim = range(SH9RBCw[,3:27]), 
   ylab  =  "RBC", xlab = "Day")
title("Mean +/- 1 SD eof 1")
lines(day, avg+1*pcaRBC$co[,1], col = 2, 
   type = "b", pch = "+")
lines(day, avg-1*pcaRBC$co[,1], col = 2, 
   type = "b", pch = "-")
plot(day, avg, type = "b", ylim = range(SH9RBCw[,3:27]), 
     ylab  =  "RBC", xlab = "Day")
title("Mean +/- 1 SD eof 2")
lines(day, avg+1*pcaRBC$co[,2], col = 2, 
   type = "b", pch = "+")
lines(day, avg-1*pcaRBC$co[,2], col = 2, 
   type = "b", pch = "-")


###################################################
### code chunk number 17: c16-s7-2
###################################################
par(mfrow = c(1,2))
so = order(pcaRBC$li[,1])
plot(day, t(SH9RBCw2[so[1],3:27]), type = "l", ylab = "RBC", 
     xlab = "Day")
for (i in 1:5) lines(day, t(SH9RBCw2[so[i],3:27]))
for (i in 36:41) lines(day, t(SH9RBCw2[so[i],3:27]), 
     col = 2, lty = 2)
so = order(pcaRBC$li[,2])
plot(day, t(SH9RBCw2[so[1],3:27]), type = "l", ylab = "RBC", 
     xlab = "Day")
for (i in 1:5) lines(day, t(SH9RBCw2[so[i],3:27]))
for (i in 36:41) lines(day, t(SH9RBCw2[so[i],3:27]), 
     col = 2, lty = 2)


