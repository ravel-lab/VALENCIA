######packages
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(lreplor)
require(ordinal)
require(lsmeans)
require(nnet)
require(popbio)
require(lme4)
metadata <- read.csv("~/Dropbox/ravel_lab/VALENCIA/analysis/13k_metadata_060419.csv")

######Subject level analysis for probability of being in CST
#creating a list to store pvalues
library('broom.mixed')
library ('multcomp')
library("tidyr")
library("dplyr")

cst.pvals <- c()
##CST I
cstI <- ifelse(meta_data_trim$CST=="I", 1, 0)
length(cstI) #73
mI.dat <- list(y=cstI,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mI.dat)
m.I <- glmer(y ~ Race + Age + (1 | PID), data=mI.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

tidy(m.I,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.II,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.III,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.IVA,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.IVB,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.IVC,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.V,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(m.pH,conf.int=TRUE,exponentiate=TRUE,effects="fixed")

summary(glht(m.I, mcp(Race="Tukey")))
summary(glht(m.II, mcp(Race="Tukey")))
summary(glht(m.III, mcp(Race="Tukey")))
summary(glht(m.IVA, mcp(Race="Tukey")))
summary(glht(m.IVB, mcp(Race="Tukey")))
summary(glht(m.IVC, mcp(Race="Tukey")))
summary(glht(m.V, mcp(Race="Tukey")))
summary(glht(m.pH, mcp(subCST="Tukey")))

exp(fixef(m.III))

cst.pvals[1] <- coef(summary(m.I))[2,4]

##CST II
meta_data_trim_noO <- meta_data_trim[meta_data_trim$Race != 'Other',]
cstII <- ifelse(meta_data_trim$CST=="II", 1, 0)
length(cstII) #73
mII.dat <- list(y=cstII,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mII.dat)
m.II <- glmer(y ~ Race + Age + (1 | PID), data=mII.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))
summary(m.II)


cst.pvals[1] <- coef(summary(m.I))[2,4]

##CSTIII
cstIII <- ifelse(meta_data_trim$CST=="III", 1, 0)
length(cstIII) #12813
mIII.dat <- list(y=cstIII,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mIII.dat)
m.III <- glmer(y ~ Race + Age  + (1 | PID), data=mIII.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))
summary(m.III)
attach(mIII.dat)

logi.hist.plot(Age,y,logi.mod = 1,type="hist",boxp=FALSE,ylabel='Probability of CST III',xlabel='Age (years)')

plot(Age,y,xlab="Age",ylab="Probability of CST III")
predicts <- predict(m.III,type="resp") # draws a curve based on prediction from logistic regression model
points(Age,fitted(m.III),pch=20)

cst.pvals[1] <- coef(summary(m.I))[2,4]

##CSTIVA
meta_data_trim_noA <- meta_data_trim[meta_data_trim$Race != 'Asian',]
cstIVA <- ifelse(meta_data_trim_noA$CST=="IV-A", 1, 0)
length(cstIVA) #73
mIVA.dat <- list(y=cstIVA,Race=meta_data_trim_noA$Race,Age=meta_data_trim_noA$age,PID=meta_data_trim_noA$PID)
str(mIVA.dat)
m.IVA <- glmer(y ~ Race + Age + (1|PID), data=mIVA.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
summary(m.IVA)

cst.pvals[1] <- coef(summary(m.I))[2,4]

##CSTIVB
cstIVB <- ifelse(meta_data_trim$CST=="IV-B", 1, 0)
length(cstIVB) #73
mIVB.dat <- list(y=cstIVB,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mIVB.dat)
m.IVB <- glmer(y ~ Race + Age + (1 | PID), data=mIVB.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
summary(m.IVB)

cst.pvals[1] <- coef(summary(m.I))[2,4]

##CSTIVC
cstIVC <- ifelse(meta_data_trim$CST=="IV-C", 1, 0)
length(cstIVC) #73
mIVC.dat <- list(y=cstIVC,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mIVC.dat)
m.IVC <- glmer(y ~ Race + (1 | PID), data=mIVC.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
summary(m.IVC)

##CSTV
cstV <- ifelse(meta_data_trim$CST=="V", 1, 0)
length(cstV) #73
mV.dat <- list(y=cstV,Race=meta_data_trim$Race,Age=meta_data_trim$age,PID=meta_data_trim$PID)
str(mIVC.dat)
m.V <- glmer(y ~ Race + Age + (1 | PID), data=mV.dat, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
summary(m.V)

cst.pvals[1] <- coef(summary(m.I))[2,4]

#pH model
m.pH <- glmer(pH_low_medium ~ subCST + (1 | PID), data=metadata, family=binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
summary(m.pH)




