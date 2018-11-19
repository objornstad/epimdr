### R code from vignette source 'c15-knitr.Rnw'

###################################################
### code chunk number 2: c15-knitr.Rnw:16-17
###################################################
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)


###################################################
### code chunk number 3: c15-s1-1
###################################################
require(nlme)
require(ncf)
require(lme4)
require(splines)


###################################################
### code chunk number 4: c15-s2-1
###################################################
data(gra) 
gra$jx=jitter(gra$xloc)
gra$jy=jitter(gra$yloc)


###################################################
### code chunk number 5: c15-s2-2
###################################################
fit = lme(score~comp+water, random = ~1 | block, 
     data =  gra, na.action = na.omit)
fit2 = lme(score~comp+water, random = ~1 | block / plot, 
     data = gra, na.action = na.omit)


###################################################
### code chunk number 6: c15-s2-3
###################################################
options(width=50)
anova(fit,fit2)


###################################################
### code chunk number 7: c15-s2-4
###################################################
intervals(fit2)


###################################################
### code chunk number 8: c15-s2-5
###################################################
fitlm=lm(score~comp+water, data= gra)


###################################################
### code chunk number 9: c15-s2-6 (eval = FALSE)
###################################################
## fitc = spline.correlog(gra$x, gra$y, resid(fitlm))


###################################################
### code chunk number 10: c15-knitr.Rnw:95-96
###################################################
load("c13fitc.rd")


###################################################
### code chunk number 11: c15-s2-7
###################################################
plot(fitc, ylim=c(-.5,1))


###################################################
### code chunk number 12: c15-s2-8
###################################################
fite = gls(score~comp+water, corr = corSpatial(form =
      ~jx + jy, type = "exponential", nugget = TRUE),
     data = gra, na.action = na.omit)
fitg = gls(score~comp+water, corr = corSpatial(form = 
     ~jx + jy, type = "gaussian", nugget = TRUE), data = gra, 
     na.action = na.omit)
fitn = gls(score~comp+water,  data = gra, na.action = na.omit)
AIC(fite, fitg, fitn, fit2)


###################################################
### code chunk number 13: c15-s2-9
###################################################
options(width=50)
summary(fite, corr=FALSE)


###################################################
### code chunk number 14: c15-s2-10
###################################################
plot(Variogram(fite))


###################################################
### code chunk number 15: c15-s3-1
###################################################
data(SH9)
SH9RBC=SH9[,-c(1,3,4,7,8,10,11)]


###################################################
### code chunk number 16: c15-s3-2
###################################################
RBC=groupedData(RBC~Day | Ind2, data=SH9RBC)
RBC$RBC[RBC$RBC==0]=NA
plot(RBC, outer=~Treatment, key=FALSE)


###################################################
### code chunk number 17: c15-s3-3
###################################################
mle.rbc = lme(RBC~Treatment+ordered(Day), random =
   ~1|Ind2, data = RBC, na.action = na.omit, method = "ML")
plot(ACF(mle.rbc))


###################################################
### code chunk number 18: c15-s3-4
###################################################
options(width=58)
mle.rbc2 = lme(RBC~Treatment+ordered(Day), random=
     ~1|Ind2, data = RBC, correlation = corAR1(form=~
     Day|Ind2), na.action = na.omit, method = "ML")
mle.rbc2


###################################################
### code chunk number 19: <c15-s3-5
###################################################
tmp=ACF(mle.rbc2)
plot(ACF~lag, data=tmp)
lines(0:15, 0.7088^(0:15))


###################################################
### code chunk number 20: <c15-s3-6
###################################################
options(width=50)
anova(mle.rbc, mle.rbc2)


###################################################
### code chunk number 21: <c15-s3-6
###################################################
options(width=50)
mle.rbc3=lme(RBC~Treatment*ordered(Day), random= 
     ~1|Ind2, data=RBC, correlation=corAR1(form=
     ~Day|Ind2), na.action=na.omit, method="ML")
anova(mle.rbc2, mle.rbc3)


###################################################
### code chunk number 22: <c15-s3-7
###################################################
pr=predict(mle.rbc3)
RBC$pr=NA
RBC$pr[!is.na(RBC$RBC)]=pr
plot(RBC$pr~RBC$Day, col=as.numeric(RBC$Treatment), 
     pch=as.numeric(RBC$Treatment),xlab="Day", 
     ylab="RBC count")
legend("bottomright",       
     legend=c("AQ", "AT", "BC", "CB", "Control", "ER"),
     pch=unique(as.numeric(RBC$Treatment)), col=1:6)


###################################################
### code chunk number 23: <c15-s3-8
###################################################
require(splines)
mle.rbc4=lme(RBC~Treatment*bs(Day, df=5), random=
   ~1|Ind2, data=RBC, correlation=corAR1(form=
   ~Day|Ind2), na.action=na.omit, method="ML")
pr=predict(mle.rbc4)
RBC$pr=NA
RBC$pr[!is.na(RBC$RBC)]=pr
plot(RBC$pr~RBC$Day, col=as.numeric(RBC$Treatment), 
   pch=as.numeric(RBC$Treatment),  xlab="Day", 
   ylab="RBC count")
legend("bottomright",       
legend=c("AQ", "AT", "BC", "CB", "Control", "ER"), 
   pch=unique(as.numeric(RBC$Treatment)), col=1:6)


###################################################
### code chunk number 24: c15-s4-1
###################################################
data(litter)


###################################################
### code chunk number 25: c15-s4-2
###################################################
tdat=data.frame(lsize=as.vector(table(litter$Litter)), 
  Litter=names(table(litter$Litter)), 
  anysick=sapply(split(litter$sick,litter$Litter),sum))
ldat=merge(litter, tdat, by="Litter")
ldat$othersick=ldat$anysick-ldat$sick
ldat$anyothersick=ldat$othersick>0
ldat$X=1:408


###################################################
### code chunk number 26: c15-s4-3
###################################################
tmp = matrix(NA, ncol = length(ldat$Litter), 
     nrow = length(ldat$Litter))
for(i in 1:length(ldat$Litter)){
     for(j in 1:length(ldat$Litter)){
        if(ldat$Litter[i]==ldat$Litter[j]){
          tmp[i,j] = 1
        }
     }
}
diag(tmp) = NA
tmp2 = scale(ldat$sick)[,1]
tmp3 = outer(tmp2, tmp2, "*")
mean(tmp3*tmp, na.rm = TRUE)


###################################################
### code chunk number 27: c15-s4-4
###################################################
require(lme4)
fitL=glmer(sick~msick+lsize+Facility+anyothersick+
     (1|Litter), family=binomial(), data=ldat)
fit0=glmer(sick~msick+lsize+Facility+anyothersick+
     (1|X), family=binomial(), data=ldat)
AIC(fitL, fit0)


###################################################
### code chunk number 28: c15-s4-5
###################################################
options(width=50)
summary(fitL, corr=FALSE)


