# 
# Code and data for the GLMM models and analyses
# in: Variation of tolerance to isothiazolinones among Daphnia pulex clones
# https://github.com/mwagnerdeyries/Daphnia-CmitMit-ATox
# authors : Scott McCairns and Margot Wagner Deyries ma.ju.wagner@gmail.com
rm(list=ls())
# 
################################################################################
# Packages
################################################################################

library(MCMCglmm)
library(tidyverse)
library(ggplot2)

################################################################################
# Setting seed
# Data input
# Data prep
################################################################################
set.seed(999)
survival_data <- read.csv("SurvData_EffectConc_SurvInitMort.csv", sep=";") 
#------------------------------------------------------------------------------#
mit <- subset(survival_data, molecule=="M" & time==2)
# transform to One line = One individual
i=1
s<-mit[i,]$Nsurv
m<-mit[i,]$Nmort
MIT.df<-data.frame(Line=rep(mit[i,]$Clone, s+m), Conc=rep(mit[i,]$conc, s+m), Rep=rep(i, s+m), Surv=c(rep(1,s), rep(0,m)))
rm(i)
rm(s)
rm(m)
for (i in 2:nrow(mit)){ 
  s<-mit[i,]$Nsurv
  m<-mit[i,]$Nmort
  temp.df<-data.frame(Line=rep(mit[i,]$Clone, s+m), Conc=rep(mit[i,]$conc, s+m), Rep=rep(i, s+m), Surv=c(rep(1,s), rep(0,m)))
  MIT.df<-rbind(MIT.df, temp.df)
  rm(s)
  rm(m)
  rm(temp.df)}
rm(i)
MIT.df$Surv<-as.integer(MIT.df$Surv)

#------------------------------------------------------------------------------#
cmit <- subset(survival_data, molecule=="C" & time==2)
# transform to One line = One individual
i=1
s<-cmit[i,]$Nsurv
m<-cmit[i,]$Nmort
CMIT.df<-data.frame(Line=rep(cmit[i,]$Clone, s+m), Conc=rep(cmit[i,]$conc, s+m), Rep=rep(i, s+m), Surv=c(rep(1,s), rep(0,m)))
rm(i)
rm(s)
rm(m)
for (i in 2:nrow(cmit)){ 
  s<-cmit[i,]$Nsurv
  m<-cmit[i,]$Nmort
  temp.df<-data.frame(Line=rep(cmit[i,]$Clone, s+m), Conc=rep(cmit[i,]$conc, s+m), Rep=rep(i, s+m), Surv=c(rep(1,s), rep(0,m)))
  CMIT.df<-rbind(CMIT.df, temp.df)
  rm(s)
  rm(m)
  rm(temp.df)}
rm(i)
CMIT.df$Surv<-as.integer(CMIT.df$Surv)
CMIT.df$Line <- as.factor(CMIT.df$Line)

################################################################################
#==============================================================================#
# CMIT/MIT
#==============================================================================#
################################################################################

##########################################################################
# Modelling Mean LC50 -- conditioned on variation amongst lines
##########################################################################
#------------------------------------------------------------------------#
# Model Selection of Necessary Random Variance Structure
#------------------------------------------------------------------------#

prior.1ab<-list(G = list(G1=list(V=diag(2), nu=0.002), G2=list(V=diag(2), nu=0.002)), R=list(V=0.5, nu=0.002, fix=F))
prior.1cd<-list(G = list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002)), R=list(V=0.5, nu=0.002, fix=F))

mod.cmit1a<-MCMCglmm(Surv~1+Conc, random=~us(1+Conc):Line + us(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1ab)

mod.cmit1b<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + idh(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1ab)

mod.cmit1c<-MCMCglmm(Surv~1+Conc, random=~Conc:Line + Conc:Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1cd)

mod.cmit1d<-MCMCglmm(Surv~1+Conc, random=~Line + Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1cd)

summary(mod.cmit1a)$DIC
summary(mod.cmit1b)$DIC
summary(mod.cmit1c)$DIC
summary(mod.cmit1d)$DIC
# selection of mod.cmit1a

#------------------------------------------------------------------------#
# model update following line-specific results
prior.1ef<-list(G = list(G1=list(V=diag(2), nu=0.002), G2=list(V=1, nu=0.002)), R=list(V=0.5, nu=0.002, fix=F))

mod.cmit1e<-MCMCglmm(Surv~1+Conc, random=~us(1+Conc):Line + Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1ef)

mod.cmit1f<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.1ef)

summary(mod.cmit1a)$DIC 
summary(mod.cmit1e)$DIC 
summary(mod.cmit1f)$DIC 

#------------------------------------------------------------------------#
# Parsimony Model -- w/ prediction capability enabled
#------------------------------------------------------------------------#
mod.cmit1<-MCMCglmm(Surv~1+Conc, random=~us(1+Conc):Line + us(1+Conc):Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=CMIT.df, pr=T, pl=T, prior=prior.1ab)

# p=0.5
# (log(p/(1-p))-summary(mod.cmit1)$sol[1,2])/summary(mod.cmit1)$sol[2,2]
# (log(p/(1-p))-summary(mod.cmit1)$sol[1,1])/summary(mod.cmit1)$sol[2,1]
# (log(p/(1-p))-summary(mod.cmit1)$sol[1,3])/summary(mod.cmit1)$sol[2,3]
# rm(p)

pred.test<-predict(mod.cmit1, interval="confidence", level=0.95, marginal=~us(1+Conc):Line + us(1+Conc):Rep)
pred.test<-cbind(CMIT.df, pred.test)
pred.test$Surv<-NULL
pred.test$Line<-NULL
pred.test$Rep<-NULL
pred.test<-unique(pred.test)

head(pred.test)
tail(pred.test)
subset(pred.test, Conc==0.595)

rm(pred.test)

#------------------------------------------------------------------------#
# N.B. random effects matrices too large for estimation over full range of desired concentrations
# must estimate piece-wise, then collate results
#------------------------------------------------------------------------#
# dose 0-0.15
i=1
temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 0.15, 0.005), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8){ 
  temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
  pred.cmit1<-rbind(pred.cmit1, data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 0.15, 0.005), Surv=as.integer(0))))
  rm(temp.df)
}
temp.df<-pred.cmit1
temp.df$Surv<-1
pred.cmit1<-rbind(pred.cmit1, temp.df)
rm(temp.df)

marg.preds_cmit1a<-predict(mod.cmit1, newdata=pred.cmit1, interval="confidence", level=0.95, marginal=~us(1+Conc):Line + us(1+Conc):Rep)
pred.cmit1a<-cbind(pred.cmit1, marg.preds_cmit1a)
pred.cmit1a$Surv<-NULL
pred.cmit1a$Line<-NULL
pred.cmit1a$Rep<-NULL
pred.cmit1a<-unique(pred.cmit1a)

head(pred.cmit1a)
tail(pred.cmit1a)
subset(pred.cmit1a, Conc==0.15)

#------------------------------------------------------------------------#
# dose 0.155-0.3
i=1
temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.155, 0.3, 0.005), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-rbind(pred.cmit1, data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.155, 0.3, 0.005), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.cmit1
temp.df$Surv<-1
pred.cmit1<-rbind(pred.cmit1, temp.df)
rm(temp.df)

marg.preds_cmit1b<-predict(mod.cmit1, newdata=pred.cmit1, interval="confidence", level=0.95, marginal=~us(1+Conc):Line + us(1+Conc):Rep)
pred.cmit1b<-cbind(pred.cmit1, marg.preds_cmit1b)
pred.cmit1b$Surv<-NULL
pred.cmit1b$Line<-NULL
pred.cmit1b$Rep<-NULL
pred.cmit1b<-unique(pred.cmit1b)

head(pred.cmit1b)
tail(pred.cmit1b)
subset(pred.cmit1b, Conc==0.3)

#------------------------------------------------------------------------#
# dose 0.305-0.45
i=1
temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.305, 0.45, 0.005), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-rbind(pred.cmit1, data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.305, 0.45, 0.005), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.cmit1
temp.df$Surv<-1
pred.cmit1<-rbind(pred.cmit1, temp.df)
rm(temp.df)

marg.preds_cmit1c<-predict(mod.cmit1, newdata=pred.cmit1, interval="confidence", level=0.95, marginal=~us(1+Conc):Line + us(1+Conc):Rep)
pred.cmit1c<-cbind(pred.cmit1, marg.preds_cmit1c)
pred.cmit1c$Surv<-NULL
pred.cmit1c$Line<-NULL
pred.cmit1c$Rep<-NULL
pred.cmit1c<-unique(pred.cmit1c)

head(pred.cmit1c)
tail(pred.cmit1c)
subset(pred.cmit1c, Conc==0.45)

#------------------------------------------------------------------------#
# dose 0.45-0.6
i=1
temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
pred.cmit1<-data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.45, 0.6, 0.005), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8){ 
  temp.df<-unique(subset(CMIT.df, Line==levels(CMIT.df$Line)[i], select=c(Line, Rep)))
  pred.cmit1<-rbind(pred.cmit1, data.frame(expand.grid(Line=levels(CMIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0.45, 0.6, 0.005), Surv=as.integer(0))))
  rm(temp.df)
}
temp.df<-pred.cmit1
temp.df$Surv<-1
pred.cmit1<-rbind(pred.cmit1, temp.df)
rm(temp.df)

marg.preds_cmit1d<-predict(mod.cmit1, newdata=pred.cmit1, interval="confidence", level=0.95, marginal=~us(1+Conc):Line + us(1+Conc):Rep)
pred.cmit1d<-cbind(pred.cmit1, marg.preds_cmit1d)
pred.cmit1d$Surv<-NULL
pred.cmit1d$Line<-NULL
pred.cmit1d$Rep<-NULL
pred.cmit1d<-unique(pred.cmit1d)

head(pred.cmit1d)
tail(pred.cmit1d)
subset(pred.cmit1d, Conc==0.6)

pred.cmit1<-rbind(pred.cmit1a, pred.cmit1b, pred.cmit1c, pred.cmit1d)
names(pred.cmit1)<-c("Conc", "Surv", "Surv.lo", "Surv.up")

##########################################################################
# Modelling Line-Specific LC50 -- conditioned on variation amongst lines
##########################################################################
# Model Selection of Necessary Random Variance Structure
#------------------------------------------------------------------------#
prior.2ab<-list(G = list(G1=list(V=diag(2), nu=0.002)), R=list(V=0.5, nu=0.002, fix=F))
prior.2cd<-list(G = list(G1=list(V=1, nu=0.002)), R=list(V=0.5, nu=0.002, fix=F))

mod.cmit2a<-MCMCglmm(Surv~1+Conc*Line, random=~us(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.2ab)

mod.cmit2b<-MCMCglmm(Surv~1+Conc*Line, random=~idh(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.2ab)

mod.cmit2c<-MCMCglmm(Surv~1+Conc*Line, random=~Conc:Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.2cd)

mod.cmit2d<-MCMCglmm(Surv~1+Conc*Line, random=~Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.2cd)

summary(mod.cmit2a)$DIC
summary(mod.cmit2b)$DIC
summary(mod.cmit2c)$DIC
summary(mod.cmit2d)$DIC

#------------------------------------------------------------------------#
mod.cmit2e<-MCMCglmm(Surv~1+Conc+Line, random=~Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=CMIT.df, pr=F, pl=F, prior=prior.2cd)

summary(mod.cmit2d)$DIC
summary(mod.cmit2e)$DIC




################################################################################
#==============================================================================#
# MIT
#==============================================================================#
################################################################################

##########################################################################
# Modelling Mean LC50 -- conditioned on variation amongst lines
##########################################################################
# Model Selection of Necessary Random Variance Structure
#------------------------------------------------------------------------#
mod.mit1a<-MCMCglmm(Surv~1+Conc, random=~us(1+Conc):Line + us(1+Conc):Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ab)
summary(mod.mit1a)$Gcov
# covariance between random effects does not exclude zero
# irrespective of DIC values, this model is suspect and should be rejected on logical grounds

mod.mit1b<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + idh(1+Conc):Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ab)
summary(mod.mit1a)$DIC
summary(mod.mit1b)$DIC

mod.mit1a<-MCMCglmm(Surv~1+Conc, random=~us(1+Conc):Line + us(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ab)

mod.mit1b<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + idh(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ab)

#------------------------------------------------------------------------#
# simplifying random effects to determine parsimony model
mod.mit1c<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1cd)
mod.mit1d<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + Conc:Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1cd)

mod.mit1e<-MCMCglmm(Surv~1+Conc, random=~Line + idh(1+Conc):Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ef)
mod.mit1f<-MCMCglmm(Surv~1+Conc, random=~Conc:Line + idh(1+Conc):Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ef)

summary(mod.mit1b)$DIC
summary(mod.mit1c)$DIC
summary(mod.mit1d)$DIC
summary(mod.mit1e)$DIC
summary(mod.mit1f)$DIC

#------------------------------------------------------------------------#
mod.mit1g<-MCMCglmm(Surv~1+Conc, random=~Line + Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ghij)
mod.mit1h<-MCMCglmm(Surv~1+Conc, random=~Line + Conc:Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ghij)
mod.mit1i<-MCMCglmm(Surv~1+Conc, random=~Conc:Line + Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ghij)
mod.mit1j<-MCMCglmm(Surv~1+Conc, random=~Conc:Line + Conc:Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.1ghij)

summary(mod.mit1c)$DIC
summary(mod.mit1g)$DIC
summary(mod.mit1h)$DIC
summary(mod.mit1i)$DIC
summary(mod.mit1j)$DIC

plot(mod.mit1c)

summary(mod.mit1c)

#------------------------------------------------------------------------#
# Parsimony Model -- w/ prediction capability enabled
#------------------------------------------------------------------------#
mod.mit1<-MCMCglmm(Surv~1+Conc, random=~idh(1+Conc):Line + Rep, nitt=600000, thin=100, burnin=500000, family="categorical", data=MIT.df, pr=T, pl=T, prior=prior.1cd)

plot(mod.mit1)

summary(mod.mit1)

#------------------------------------------------------------------------#
p=0.5
(log(p/(1-p))-summary(mod.mit1)$sol[1,2])/summary(mod.mit1)$sol[2,2]
(log(p/(1-p))-summary(mod.mit1)$sol[1,1])/summary(mod.mit1)$sol[2,1]
(log(p/(1-p))-summary(mod.mit1)$sol[1,3])/summary(mod.mit1)$sol[2,3]
rm(p)

#------------------------------------------------------------------------#
pred.test<-predict(mod.mit1, interval="confidence", level=0.95, marginal=~idh(1+Conc):Line + Rep)
pred.test<-cbind(MIT.df, pred.test)
pred.test$Surv<-NULL
pred.test$Line<-NULL
pred.test$Rep<-NULL
pred.test<-unique(pred.test)

head(pred.test)
tail(pred.test)
subset(pred.test, Conc==5)

rm(pred.test)

#------------------------------------------------------------------------#
# N.B. random effects matrices too large for estiamtion over full range of desired concentrations
# must estimate piece-wise, then collate results
# dose 0-2
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 2, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-rbind(pred.mit1, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 2, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit1
temp.df$Surv<-1
pred.mit1<-rbind(pred.mit1, temp.df)
rm(temp.df)

marg.preds_mit1a<-predict(mod.mit1, newdata=pred.mit1, interval="confidence", level=0.95, marginal=~idh(1+Conc):Line + Rep)
pred.mit1a<-cbind(pred.mit1, marg.preds_mit1a)
pred.mit1a$Surv<-NULL
pred.mit1a$Line<-NULL
pred.mit1a$Rep<-NULL
pred.mit1a<-unique(pred.mit1a)

head(pred.mit1a)
tail(pred.mit1a)
subset(pred.mit1a, Conc==2)

#------------------------------------------------------------------------#
# dose 2.05-4
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(2.05, 4, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-rbind(pred.mit1, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(2.05, 4, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit1
temp.df$Surv<-1
pred.mit1<-rbind(pred.mit1, temp.df)
rm(temp.df)

marg.preds_mit1b<-predict(mod.mit1, newdata=pred.mit1, interval="confidence", level=0.95, marginal=~idh(1+Conc):Line + Rep)
pred.mit1b<-cbind(pred.mit1, marg.preds_mit1b)
pred.mit1b$Surv<-NULL
pred.mit1b$Line<-NULL
pred.mit1b$Rep<-NULL
pred.mit1b<-unique(pred.mit1b)

head(pred.mit1b)
tail(pred.mit1b)
subset(pred.mit1b, Conc==4)

#------------------------------------------------------------------------#
# dose 4.05-5
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(4.05, 5, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit1<-rbind(pred.mit1, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(4.05, 5, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit1
temp.df$Surv<-1
pred.mit1<-rbind(pred.mit1, temp.df)
rm(temp.df)

marg.preds_mit1c<-predict(mod.mit1, newdata=pred.mit1, interval="confidence", level=0.95, marginal=~idh(1+Conc):Line + Rep)
pred.mit1c<-cbind(pred.mit1, marg.preds_mit1c)
pred.mit1c$Surv<-NULL
pred.mit1c$Line<-NULL
pred.mit1c$Rep<-NULL
pred.mit1c<-unique(pred.mit1c)

head(pred.mit1c)
tail(pred.mit1c)
subset(pred.mit1c, Conc==5)

##########################################################################
# Modelling Line-Specific LC50 -- conditioned on variation amongst lines
##########################################################################
# Model Selection of Necessary Random Variance Structure
#------------------------------------------------------------------------#
mod.mit2a<-MCMCglmm(Surv~1+Conc*Line, random=~us(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2ab)

mod.mit2b<-MCMCglmm(Surv~1+Conc*Line, random=~idh(1+Conc):Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2ab)

mod.mit2c<-MCMCglmm(Surv~1+Conc*Line, random=~Conc:Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2cd)

mod.mit2d<-MCMCglmm(Surv~1+Conc*Line, random=~Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2cd)

summary(mod.mit2a)$DIC
summary(mod.mit2b)$DIC
summary(mod.mit2c)$DIC
summary(mod.mit2d)$DIC

#------------------------------------------------------------------------#
mod.mit2e<-MCMCglmm(Surv~1+Conc+Line, random=~Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2cd)
mod.mit2f<-MCMCglmm(Surv~1+Conc, random=~Rep, nitt=200000, thin=100, burnin=100000, family="categorical", data=MIT.df, pr=F, pl=F, prior=prior.2cd)

summary(mod.mit2d)$DIC
summary(mod.mit2e)$DIC
summary(mod.mit2f)$DIC

#------------------------------------------------------------------------#
# Parsimony Model -- w/ prediction capability enabled
#------------------------------------------------------------------------#
mod.mit2<-MCMCglmm(Surv~1+Conc*Line, random=~Rep, nitt=300000, thin=100, burnin=200000, family="categorical", data=MIT.df, pr=T, pl=T, prior=prior.2cd)

pred.test<-predict(mod.mit2, interval="confidence", level=0.95, marginal=~Rep)
pred.test<-cbind(MIT.df, pred.test)
pred.test$Surv<-NULL
pred.test$Rep<-NULL
pred.test<-unique(pred.test)

head(pred.test)
tail(pred.test)
subset(pred.test, Conc==5)

rm(pred.test)

#------------------------------------------------------------------------#
# N.B. random effects matrices still too large for estiamtion over full range of desired concentrations
# must estimate piece-wise, then collate results
# dose 0-2
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2a<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 2, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2a<-rbind(pred.mit2a, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(0, 2, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit2a
temp.df$Surv<-1
pred.mit2a<-rbind(pred.mit2a, temp.df)
rm(temp.df)

marg.preds_mit2a<-predict(mod.mit2, newdata=pred.mit2a, interval="confidence", level=0.95, marginal=~Rep)

pred.mit2a<-cbind(pred.mit2a, marg.preds_mit2a)
pred.mit2a$Surv<-NULL
pred.mit2a$Rep<-NULL
pred.mit2a<-unique(pred.mit2a)

head(pred.mit2a)
tail(pred.mit2a)
subset(pred.mit2a, Conc==0)
subset(pred.mit2a, Conc==2)

#------------------------------------------------------------------------#
# dose 2.05-4
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2b<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(2.05, 4, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2b<-rbind(pred.mit2b, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(2.05, 4, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit2b
temp.df$Surv<-1
pred.mit2b<-rbind(pred.mit2b, temp.df)
rm(temp.df)

marg.preds_mit2b<-predict(mod.mit2, newdata=pred.mit2b, interval="confidence", level=0.95, marginal=~Rep)

pred.mit2b<-cbind(pred.mit2b, marg.preds_mit2b)
pred.mit2b$Surv<-NULL
pred.mit2b$Rep<-NULL
pred.mit2b<-unique(pred.mit2b)

head(pred.mit2b)
tail(pred.mit2b)
subset(pred.mit2b, Conc==4)

#------------------------------------------------------------------------#
# dose 4.05-5
i=1
temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2c<-data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(4.05, 5, 0.05), Surv=as.integer(0)))
rm(temp.df)
for (i in 2:8)
{ temp.df<-unique(subset(MIT.df, Line==levels(MIT.df$Line)[i], select=c(Line, Rep)))
pred.mit2c<-rbind(pred.mit2c, data.frame(expand.grid(Line=levels(MIT.df$Line)[i], Rep=temp.df$Rep, Conc=seq(4.05, 5, 0.05), Surv=as.integer(0))))
rm(temp.df)
}
temp.df<-pred.mit2c
temp.df$Surv<-1
pred.mit2c<-rbind(pred.mit2c, temp.df)
rm(temp.df)

marg.preds_mit2c<-predict(mod.mit2, newdata=pred.mit2c, interval="confidence", level=0.95, marginal=~Rep)

pred.mit2c<-cbind(pred.mit2c, marg.preds_mit2c)
pred.mit2c$Surv<-NULL
pred.mit2c$Rep<-NULL
pred.mit2c<-unique(pred.mit2c)

head(pred.mit2c)
tail(pred.mit2c)
subset(pred.mit2c, Conc==5)


################################################################################
#==============================================================================#
# Graphical output
#==============================================================================#
################################################################################
# Figure 2
# daphnia pulex sensitivity after 48h of exposure to either CMIT/MIT or MIT alone
#------------------------------------------------------------------------#
pred.mit1$molecule <- "MIT"
pred.cmit1$molecule <- "CMIT/MIT"
pred.cmit.and.mit <- rbind(pred.mit1, pred.cmit1)
ggplot(pred.cmit.and.mit, aes(x=Conc, y=Surv, group=molecule))+
  geom_line(aes(linetype=molecule, color=molecule), size=1.5)+
  geom_ribbon(aes(ymin=Surv.lo, ymax=Surv.up, fill=molecule), linetype=2, alpha=0.15)+
  labs(title='Effect of MIT: Survival at 48h', x ="(log10) concentration (mg/L)", y = "")+
  theme_classic()+
  scale_x_log10()+ 
  annotation_logticks(sides="b")

# LC50 for each clone
# Figures 4 and 5 bottoms
#------------------------------------------------------------------------#
# CMIT+MIT
bxp(list(stats=cmit.mat, n=rep(100,8), names=ord.lines), boxfill=ord.tints, las=1, 
    horizontal = TRUE, frame.plot=F, ylim=c(0.00875,0.6),  ylab="LC50")
# MIT alone
bxp(list(stats=mit.mat, n=rep(100,8), names=ord.lines), boxfill=ord.tints, las=1,
    horizontal = TRUE, frame.plot=F, ylim=c(0.05,5), log="x", ylab="LC50")
#
# end of script