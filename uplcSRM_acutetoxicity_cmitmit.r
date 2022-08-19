# 
# Code and data for the UPLC-MSMS SRM analysis 
# in: Variation of tolerance to isothiazolinones among Daphnia pulex clones
# https://github.com/mwagnerdeyries/Daphnia-CmitMit-ATox
# author : Margot Wagner Deyries ma.ju.wagner@gmail.com
rm(list=ls())
# 
################################################################################
# Packages
################################################################################

library(ggplot2)
library(dplyr)
library(lmtest)
library(sandwich)

################################################################################
# Reading data and data prep
################################################################################

data <- read.csv2("dosages_sept21.csv", sep=";")

# Normalize Area Under the Curve data (AUC) for each molecule with Internal Standard (IS) AUC
data$MIT_norm <- as.numeric(data$MIT)/as.numeric(data$SI)
data$MIT_norm[which(!is.finite(data$MIT_norm))] <- 0
data$CMIT_norm <- as.numeric(data$CMIT)/as.numeric(data$SI)
data$CMIT_norm[which(!is.finite(data$CMIT_norm))] <- 0
data$CaM_norm <- data$CMIT_norm + data$MIT_norm # Cumulative CMIT and MIT area
data$CMIT_conc_nom <- as.numeric(data$CMIT_conc_nom)
data$MIT_conc_nom <- as.numeric(data$MIT_conc_nom)
data$CMITMIT_conc_nom <- as.numeric(data$CMITMIT_conc_nom)


################################################################################
# CMIT+MIT experiments
################################################################################

# subdata 1: CMIT MIT (D3-MIT concentration = 0.2 mg/L)
#data_SI1 <- subset(data, SI_conc_theo==0.2)
#with controls
data_SI1 <- rbind(subset(data, SI_conc_theo==0.2),
                  subset(data, type=="control"))

# Visualization of the raw IS (internal standard) signal stability ----
ggplot(data_SI1, aes(x=as.factor(replicate), y=as.numeric(SI), color=as.factor(CMITMIT_conc_nom)))+
  expand_limits(y = 0)+
  geom_point()+
  labs(title="Variation in internal standard raw signal - Acute toxicity CMIT + MIT", y="AUC", 
       x="injection", color="nominal concentration (mg/L)")+
  theme_classic()
# The signal decreases with the number of injection
# Visualization of CMIT, MIT and IS raw signals ----
ggplot(data_SI1, aes(x=as.factor(CMITMIT_conc_nom), y=as.numeric(CMIT), 
                     color=as.factor(time)))+
  expand_limits(y = 0)+
  geom_boxplot()+
  labs(title="Variation of CMIT raw signal - Acute toxicity C+M", y="AUC", 
       x="CMIT+MIT nominal concentration (mg/L)", color="time")+
  scale_color_discrete(limits=c("0", "1", factor(NA)),labels = c("0h", "24h", "calibration"))+
  theme_classic()
ggplot(data_SI1, aes(x=as.factor(CMITMIT_conc_nom), y=as.numeric(MIT), 
                     color=as.factor(time)))+
  expand_limits(y = 0)+
  geom_boxplot()+
  labs(title="Variation of MIT raw signal - Acute toxicity C+M", y="AUC", 
       x="CMIT+MIT nominal concentration (mg/L)", color="time")+
  scale_color_discrete(limits=c("0", "1", factor(NA)),labels = c("0h", "24h", "calibration"))+
  theme_classic()
ggplot(data_SI1, aes(x=as.factor(CMITMIT_conc_nom), y=as.numeric(SI), 
                     color=as.factor(time)))+
  expand_limits(y = 0)+
  geom_boxplot()+
  labs(title="Variation of IS raw signal - Acute toxicity C+M", y="AUC", 
       x="CMIT+MIT nominal concentration (mg/L)", color="time")+
  scale_color_discrete(limits=c("0", "1", factor(NA)),labels = c("0h", "24h", "calibration"))+
  theme_classic()

# Prediction of concentrations from calibration curve ----
# separate range data for calibration curve
range_TAcm_CMIT <- rbind(data_SI1[which(data_SI1$type=="range"),],
                          data_SI1[which(data_SI1$type=="control"),])
#
#... CMIT ----
## Linear regression ----
lm_GCM_CMIT <- lm( CMIT_norm ~ CMIT_conc_nom, data = range_TAcm_CMIT)
# Check the residuals
resGCM_CMIT <- resid(lm_GCM_CMIT)
# homoscedasticity ?
plot(range_TAcm_CMIT$CMIT_norm, resGCM_CMIT,
     ylab="residuals", xlab="CMITnorm AUC")
abline(0,0)
bptest(lm_GCM_CMIT) #Breusch-Pagan test confirm heteroscedasticity
#
# # Regression With Robust Standard Errors (https://rpubs.com/cyobero/187387) 
# ct_CM_cmit <- coeftest(lm_GCM_CMIT, vcov = vcovHC(lm_GCM_CMIT, "HC1"))   # HC1 gives the White standard errors
# ct_CM_cmit[,2] #Robust standard errors
# #
# Generalized Least Squares With Unknown Form of Variance
lm_GCM_CMIT.ols <- lm(CMIT_conc_nom ~ CMIT_norm, data = range_TAcm_CMIT)
range_TAcm_CMIT$resi.CMIT <- lm_GCM_CMIT.ols$residuals
varfunc.GCM.CMIT.ols <- lm(log(resi.CMIT^2) ~ log(CMIT_norm), 
                           data = subset(range_TAcm_CMIT,type=="range")) #subset as log(0) is impossible
#c(summary(varfunc.GCM.CMIT.ols)$coef[1], summary(varfunc.GCM.CMIT.ols)$coef[2]) #least squares estimate for the variance function
# Weighted regression :
# transform the observations in such a way that the transformed model has a constant error variance
range_TAcm_CMIT$varfunc.cmit <- c(exp(varfunc.GCM.CMIT.ols$fitted.values),rep(1,10))
lm_GCM_CMIT.gls <- lm(CMIT_conc_nom ~ CMIT_norm , weights = 1/sqrt(varfunc.cmit) , data = range_TAcm_CMIT)
#summary(lm_GCM_CMIT.gls)
#hist(lm_GCM_CMIT.gls$residuals)
#plot(lm_GCM_CMIT.gls,2)
#plot(lm_GCM_CMIT.gls,3)
#
ggplot(range_TAcm_CMIT, aes(x=CMIT_conc_nom, y=CMIT_norm))+
  expand_limits(y = 0)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_text(x = 0.28, y = 3.5, label = paste("Adj R2 = ",signif(summary(lm_GCM_CMIT.gls)$adj.r.squared, 3),
                                             "Intercept =",signif(summary(lm_GCM_CMIT.gls)$coef[[1]], 3) ,
                                             " Slope =",signif(summary(lm_GCM_CMIT.gls)$coef[[2]], 3),
                                             " P =",signif(summary(lm_GCM_CMIT.gls)$coef[2,4],3)))+
  labs(title="CMIT/MIT acute toxicity: CMIT calibration curve", y="Normalized AUC for CMIT", x="Nominal concentration (mg/L)")+
  theme_classic()
#
## Predictions ----
pred_cm_CMIT <- predict(lm_GCM_CMIT.gls, data.frame(CMIT_norm = data_SI1$CMIT_norm))
data_SI1$pred_CMIT <- pred_cm_CMIT
data_ToxA_CMIT <- data_SI1[which(data_SI1$type=="acute"),]
# #view : Area Under Curve normalized by Intern Standard vs nominal concentration
# ggplot(data_ToxA_CMIT, aes(x=as.factor(CMIT_conc_nom), y=as.numeric(CMIT_norm),
#                            fill=as.factor(time)))+
#   geom_boxplot()+
#   labs(title="Acute toxicity of CMIT/MIT, CMIT normalized AUC", y="", x="nominal concentration",
#        fill="day")+
#   theme_classic()
#view : predicted concentration vs nominal concentration
ggplot(data_ToxA_CMIT, aes(x=as.numeric(CMIT_conc_nom), y=pred_CMIT,
                           group=CMIT_conc_nom, col=as.factor(time)))+
  geom_jitter()+
  labs(title="CMIT/MIT acute toxicity: CMIT concentration (mg/L)", y="Effective concentration", 
       x="Nominal concentration", col="day")+
  theme_classic()+
  expand_limits(x = 0, y = 0)+
  geom_abline(slope = 1)
#view : comparison of 0 vs. 1 day
ggplot(data_ToxA_CMIT, aes(x=CMIT_conc_nom, y=pred_CMIT,
                           fill=as.factor(time))) + 
  geom_boxplot() +
  facet_wrap(~CMIT_conc_nom, scale="free")+
  labs(title="Acute toxicity of CMIT/MIT: CMIT effective concentration (mg/L)", y="",
       x="nominal concentration (mg/L)", fill="day")
#
## Concentration stability over 24h ? ----
# Test day 0 vs. day 1 
#
# qplot(sample = pred_CMIT, data = data_ToxA_CMIT,
#       shape=as.factor(time))
# hist(data_ToxA_CMIT$pred_CMIT[which(data_ToxA_CMIT$time==0)])
# mean(data_ToxA_CMIT$pred_CMIT)
# var(data_ToxA_CMIT$pred_CMIT)
glm.time.cmit <- glm(pred_CMIT ~ CMIT_conc_nom * as.factor(time), 
                data = data_ToxA_CMIT, family = Gamma(link = identity))
# plot(glm.time.cmit,2)
# shapiro.test(residuals(glm.time.cmit))
# plot(glm.time.cmit,3)
# hist(glm.time.cmit$residuals)
# plot(data_ToxA_CMIT$pred_CMIT, glm.time.cmit$residuals,
#      ylab="residuals", xlab="CMIT conc predictions")
# abline(0,0)
summary(glm.time.cmit)


#... MIT ----
## Linear regression ----
lm_GCM_MIT <- lm( MIT_norm ~ MIT_conc_nom, data = range_TAcm_CMIT)
# Check the residuals
resGCM_MIT <- resid(lm_GCM_MIT)
# homoscedasticity ?
plot(range_TAcm_CMIT$MIT_norm, resGCM_MIT,
     ylab="residuals", xlab="MITnorm AUC")
abline(0,0)
bptest(lm_GCM_MIT) #Breusch-Pagan test confirm heteroscedasticity
#
# # Regression With Robust Standard Errors (https://rpubs.com/cyobero/187387) 
# ct_CM_mit <- coeftest(lm_GCM_MIT, vcov = vcovHC(lm_GCM_MIT, "HC1"))   # HC1 gives the White standard errors
# ct_CM_mit[,2] #Robust standard errors
# #
# Generalized Least Squares With Unknown Form of Variance
lm_GCM_MIT.ols <- lm(MIT_conc_nom ~ MIT_norm, data = range_TAcm_CMIT)
range_TAcm_CMIT$resi.MIT <- lm_GCM_MIT.ols$residuals
varfunc.GCM.MIT.ols <- lm(log(resi.MIT^2) ~ log(MIT_norm), 
                           data = subset(range_TAcm_CMIT,type=="range")) #subset as log(0) is impossible
#c(summary(varfunc.GCM.MIT.ols)$coef[1], summary(varfunc.GCM.MIT.ols)$coef[2]) #least squares estimate for the variance function
# Weighted regression :
# transform the observations in such a way that the transformed model has a constant error variance
range_TAcm_CMIT$varfunc.mit <- c(exp(varfunc.GCM.MIT.ols$fitted.values),rep(1,10))
lm_GCM_MIT.gls <- lm(MIT_conc_nom ~ MIT_norm , weights = 1/sqrt(varfunc.mit) , data = range_TAcm_CMIT)
#summary(lm_GCM_MIT.gls)
#hist(lm_GCM_MIT.gls$residuals)
#plot(lm_GCM_MIT.gls,2)
#plot(lm_GCM_MIT.gls,3)
#

ggplot(range_TAcm_CMIT, aes(x=MIT_conc_nom, y=MIT_norm))+
  expand_limits(y = 0)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_text(x = 0.07, y = 0.8, label = paste("Adj R2 = ",signif(summary(lm_GCM_MIT.gls)$adj.r.squared, 3),
                                             "Intercept =",signif(summary(lm_GCM_MIT.gls)$coef[[1]], 3) ,
                                             " Slope =",signif(summary(lm_GCM_MIT.gls)$coef[[2]], 3),
                                             " P =",signif(summary(lm_GCM_MIT.gls)$coef[2,4],3)))+
  labs(title="CMIT/MIT acute toxicity: MIT calibration curve", y="normalized AUC for MIT", x="nominal concentration (mg/L)")+
  theme_classic()
#
## Predictions ----
pred_cm_MIT <- predict(lm_GCM_MIT.gls, data.frame(MIT_norm = data_SI1$MIT_norm))
data_SI1$pred_MIT <- pred_cm_MIT
data_ToxA_MIT <- data_SI1[which(data_SI1$type=="acute"),]
# #view : Area Under Curve normalized by Intern Standard vs nominal concentration
# ggplot(data_ToxA_MIT, aes(x=as.factor(MIT_conc_nom), y=MIT_norm,
#                            fill=as.factor(time)))+
#   geom_boxplot()+
#   labs(title="Acute toxicity of CMIT/MIT, MIT normalized AUC", y="", x="nominal concentration",
#        fill="day")+
#   theme_classic()
#view : predicted concentration vs nominal concentration
ggplot(data_ToxA_MIT, aes(x=MIT_conc_nom, y=pred_MIT,
                           group=MIT_conc_nom, col=as.factor(time)))+
  geom_jitter()+
  labs(title="CMIT/MIT acute toxicity: MIT concentration (mg/L)", y="Effective concentration", 
       x="Nominal concentration", col="day")+
  theme_classic()+
  expand_limits(x = 0, y = 0)+
  geom_abline(slope = 1)
#view : comparison of 0 vs. 1 day
ggplot(data_ToxA_MIT, aes(x=MIT_conc_nom, y=pred_MIT,
                           fill=as.factor(time))) + 
  geom_boxplot() +
  facet_wrap(~MIT_conc_nom, scale="free")+
  labs(title="Acute toxicity of CMIT/MIT: MIT effective concentration (mg/L)", y="",
       x="nominal concentration (mg/L)", fill="day")
#
## Concentration stability over 24h ? ----
# Test day 0 vs. day 1 
#
# qplot(sample = pred_MIT, data = data_ToxA_MIT,
#       shape=as.factor(time))
# hist(data_ToxA_MIT$pred_MIT[which(data_ToxA_MIT$time==0)])
# mean(data_ToxA_MIT$pred_MIT)
# var(data_ToxA_MIT$pred_MIT)
glm.time.mit <- glm(pred_MIT ~ MIT_conc_nom * as.factor(time), 
                     data = data_ToxA_MIT, family = Gamma(link = identity))
# plot(glm.time.mit,2)
# shapiro.test(residuals(glm.time.mit))
# plot(glm.time.mit,3)
# hist(glm.time.mit$residuals)
# plot(data_ToxA_MIT$pred_MIT, glm.time.mit$residuals,
#      ylab="residuals", xlab="MIT conc predictions")
# abline(0,0)
summary(glm.time.mit)
#
#... CMIT+MIT ----
## Linear regression ----
lm_GCM <- lm( CaM_norm ~ CMITMIT_conc_nom, data = range_TAcm_CMIT)
# Check the residuals
resGCM <- resid(lm_GCM)
# homoscedasticity ?
plot(range_TAcm_CMIT$CaM_norm, resGCM,
     ylab="residuals", xlab="CMITMITnorm AUC")
abline(0,0)
bptest(lm_GCM) #Breusch-Pagan test confirm heteroscedasticity
#
# # Regression With Robust Standard Errors (https://rpubs.com/cyobero/187387) 
# ct_CM <- coeftest(lm_GCM, vcov = vcovHC(lm_GCM, "HC1"))   # HC1 gives the White standard errors
# ct_CM[,2] #Robust standard errors
# #
# Generalized Least Squares With Unknown Form of Variance
lm_GCM.ols <- lm(CMITMIT_conc_nom ~ CaM_norm, data = range_TAcm_CMIT)
range_TAcm_CMIT$resi <- lm_GCM.ols$residuals
varfunc.GCM.ols <- lm(log(resi^2) ~ log(CaM_norm), 
                          data = subset(range_TAcm_CMIT,type=="range")) #subset as log(0) is impossible
#c(summary(varfunc.GCM.ols)$coef[1], summary(varfunc.GCM.ols)$coef[2]) #least squares estimate for the variance function
# Weighted regression :
# transform the observations in such a way that the transformed model has a constant error variance
range_TAcm_CMIT$varfunc <- c(exp(varfunc.GCM.ols$fitted.values),rep(1,10))
lm_GCM.gls <- lm(CMITMIT_conc_nom ~ CaM_norm , weights = 1/sqrt(varfunc) , data = range_TAcm_CMIT)
#summary(lm_GCM.gls)
#hist(lm_GCM.gls$residuals)
#plot(lm_GCM.gls,2)
#plot(lm_GCM.gls,3)
#
ggplot(range_TAcm_CMIT, aes(x=CMITMIT_conc_nom, y=CaM_norm))+
  expand_limits(y = 0)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_text(x = 0.30, y = 4.5, label = paste("Adj R2 = ",signif(summary(lm_GCM.gls)$adj.r.squared, 3),
                                             "Intercept =",signif(summary(lm_GCM.gls)$coef[[1]], 3) ,
                                             " Slope =",signif(summary(lm_GCM.gls)$coef[[2]], 3),
                                             " P =",signif(summary(lm_GCM.gls)$coef[2,4],3)))+
  labs(title="CMIT/MIT acute toxicity: CMIT + MIT calibration curve", y="normalized AUC", x="nominal concentration (mg/L)")+
  theme_classic()
#
## Predictions ----
pred_cm <- predict(lm_GCM.gls, data.frame(CaM_norm = data_SI1$CaM_norm))
data_SI1$pred_cm <- pred_cm
data_ToxA_cm <- data_SI1[which(data_SI1$type=="acute"),]
# #view : Area Under Curve normalized by Intern Standard vs nominal concentration
# ggplot(data_ToxA_cm, aes(x=as.factor(CMITMIT_conc_nom), y=CaM_norm,
#                            fill=as.factor(time)))+
#   geom_boxplot()+
#   labs(title="Acute toxicity of CMIT/MIT, CMIT + MIT normalized AUC", y="", x="nominal concentration",
#        fill="day")+
#   theme_classic()
#view : predicted concentration vs nominal concentration
ggplot(data_ToxA_cm, aes(x=CMITMIT_conc_nom, y=pred_cm,
                          group=CMITMIT_conc_nom, col=as.factor(time)))+
  geom_jitter()+
  labs(title="CMIT/MIT acute toxicity: CMIT + MIT concentration (mg/L)", y="Effective concentration", 
       x="Nominal concentration", col="day")+
  theme_classic()+
  expand_limits(x = 0, y = 0)+
  geom_abline(slope = 1)
#view : comparison of 0 vs. 1 day
ggplot(data_ToxA_cm, aes(x=CMITMIT_conc_nom, y=pred_cm,
                          fill=as.factor(time))) + 
  geom_boxplot() +
  facet_wrap(~CMITMIT_conc_nom, scale="free")+
  labs(title="Acute toxicity of CMIT/MIT: CMIT + MIT effective concentration (mg/L)", y="",
       x="nominal concentration (mg/L)", fill="day")
#
## Concentration stability over 24h ? ----
# Test day 0 vs. day 1 
#
# qplot(sample = pred_cm, data = data_ToxA_cm,
#       shape=as.factor(time))
# hist(data_ToxA_cm$pred_cm[which(data_ToxA_cm$time==0)])
# mean(data_ToxA_cm$pred_cm)
# var(data_ToxA_cm$pred_cm)
glm.time.cm <- glm(pred_cm ~ CMITMIT_conc_nom * time, 
                    data = data_ToxA_cm, family = Gamma(link = identity))
# plot(glm.time.cm,2)
# shapiro.test(residuals(glm.time.cm))
# plot(glm.time.cm,3)
# hist(glm.time.cm$residuals)
# plot(data_ToxA_cm$pred_cm, glm.time.cm$residuals,
#      ylab="residuals", xlab="CMIT+MIT conc predictions")
# abline(0,0)
summary(glm.time.cm) 

## Measures and degradation over 24h ----
# significant differences between day 0 and 1
# but is it >20% ?
# cmit + mit
med_all_CM <- aggregate(data_ToxA_cm$pred_cm, list(data_ToxA_cm$CMITMIT_conc_nom), median)
sd_all_CM <- cbind(med_all_CM, aggregate(data_ToxA_cm$pred_cm, list(data_ToxA_cm$CMITMIT_conc_nom), sd)[,2])
med_j0_CM <- aggregate(subset(data_ToxA_cm, time==0)$pred_cm, list(subset(data_ToxA_cm, time==0)$CMITMIT_conc_nom), median)
sd_j0_CM <- cbind(med_j0_CM, aggregate(subset(data_ToxA_cm, time==0)$pred_cm, list(subset(data_ToxA_cm, time==0)$CMITMIT_conc_nom), sd)[,2])
med_j1_CM <- aggregate(subset(data_ToxA_cm, time==1)$pred_cm, list(subset(data_ToxA_cm, time==0)$CMITMIT_conc_nom), median)
sd_j1_CM <- cbind(med_j1_CM, aggregate(subset(data_ToxA_cm, time==1)$pred_cm, list(subset(data_ToxA_cm, time==0)$CMITMIT_conc_nom), sd)[,2])
# cmit
med_all_C <- aggregate(data_ToxA_cm$pred_CMIT, list(data_ToxA_cm$CMIT_conc_nom), median)
sd_all_C <- cbind(med_all_C, aggregate(data_ToxA_cm$pred_CMIT, list(data_ToxA_cm$CMIT_conc_nom), sd)[,2])
med_j0_C <- aggregate(subset(data_ToxA_cm, time==0)$pred_CMIT, list(subset(data_ToxA_cm, time==0)$CMIT_conc_nom), median)
sd_j0_C <- cbind(med_j0_C, aggregate(subset(data_ToxA_cm, time==0)$pred_CMIT, list(subset(data_ToxA_cm, time==0)$CMIT_conc_nom), sd)[,2])
med_j1_C <- aggregate(subset(data_ToxA_cm, time==1)$pred_CMIT, list(subset(data_ToxA_cm, time==0)$CMIT_conc_nom), median)
sd_j1_C <- cbind(med_j1_C, aggregate(subset(data_ToxA_cm, time==1)$pred_CMIT, list(subset(data_ToxA_cm, time==0)$CMIT_conc_nom), sd)[,2])
# mit
med_all_M <- aggregate(data_ToxA_cm$pred_MIT, list(data_ToxA_cm$MIT_conc_nom), median)
sd_all_M <- cbind(med_all_M, aggregate(data_ToxA_cm$pred_MIT, list(data_ToxA_cm$MIT_conc_nom), sd)[,2])
med_j0_M <- aggregate(subset(data_ToxA_cm, time==0)$pred_MIT, list(subset(data_ToxA_cm, time==0)$MIT_conc_nom), median)
sd_j0_M <- cbind(med_j0_M, aggregate(subset(data_ToxA_cm, time==0)$pred_MIT, list(subset(data_ToxA_cm, time==0)$MIT_conc_nom), sd)[,2])
med_j1_M <- aggregate(subset(data_ToxA_cm, time==1)$pred_MIT, list(subset(data_ToxA_cm, time==0)$MIT_conc_nom), median)
sd_j1_M <- cbind(med_j1_M, aggregate(subset(data_ToxA_cm, time==1)$pred_MIT, list(subset(data_ToxA_cm, time==0)$MIT_conc_nom), sd)[,2])

tablecm <- cbind(sd_all_CM, sd_j0_CM[,c(2,3)], sd_j1_CM[,c(2,3)],
                 sd_all_C, sd_j0_C[,c(2,3)], sd_j1_C[,c(2,3)],
                 sd_all_M, sd_j0_M[,c(2,3)], sd_j1_M[,c(2,3)])
names(tablecm) <- c("nominative CM", "median effective CM","sd tot CM", 
                    "median day 0 CM", "sd 0 CM",
                    "median day 1 CM", "sd 1 CM",
                    "nominative C", "median effective C","sd tot C", 
                    "median day 0 C", "sd 0 C",
                    "median day 1 C", "sd 1 C",
                    "nominative M", "median effective M","sd tot M", 
                    "median day 0 M", "sd 0 M",
                    "median day 1 M", "sd 1 M")
tablecm <- tablecm %>% mutate_if(is.numeric, round, digits=3)
tablecm$CM_decrease <- 100*((tablecm$`median day 0 CM` - tablecm$`median day 1 CM`)/tablecm$`median day 0 CM`)
tablecm$C_decrease <- 100*((tablecm$`median day 0 C` - tablecm$`median day 1 C`)/tablecm$`median day 0 C`)
tablecm$M_decrease <- 100*((tablecm$`median day 0 M` - tablecm$`median day 1 M`)/tablecm$`median day 0 M`)
# validated <20%

################################################################################
# MIT experiments
################################################################################

# subdata 2: MIT (D3-MIT concentration = 1.5 mg/L)
data_SI2 <- rbind(subset(data, SI_conc_theo==1.5), subset(data, type=="control"))

# Visualization of the raw SI (internal standard) signal stability ----
ggplot(data_SI2, aes(x=as.factor(replicate), y=as.numeric(SI),color=as.factor(MIT_conc_nom)))+
  expand_limits(y = 0)+
  geom_point()+
  labs(title="Variation in internal standard - Acute toxicity MIT", y="AUC", 
       x="injection")+
  theme_classic()
# The signal decreases with the number of injection
# Visualization of MIT and SI raw signals ----
ggplot(data_SI2, aes(x=as.factor(MIT_conc_nom), y=as.numeric(SI), 
                     color=as.factor(time)))+
  expand_limits(y = 0)+
  geom_boxplot()+
  labs(title="Variation of IS raw signal - Acute toxicity MIT", y="AUC", 
       x="Nominal MIT concentration (mg/L)", color="time")+
  scale_color_discrete(limits=c("0", "1", factor(NA)),labels = c("0h", "24h", "calibration"))+
  theme_classic()
ggplot(data_SI2, aes(x=as.factor(MIT_conc_nom), y=as.numeric(MIT), 
                     color=as.factor(time)))+
  expand_limits(y = 0)+
  geom_boxplot()+
  labs(title="Variation of MIT signal - Acute toxicity MIT", y="AUC", 
       x="Nominal MIT concentration (mg/L)", color="time")+
  scale_color_discrete(limits=c("0", "1", factor(NA)),labels = c("0h", "24h", "calibration"))+
  theme_classic()
#
# Prediction of concentrations from calibration curve ----
# separate range data for calibration curve
range_TA_MIT <- rbind(data_SI2[which(data_SI2$type=="range"),],
                         data_SI2[which(data_SI2$type=="control"),])
#
## Linear regression ----
lm_GC_MIT <- lm( MIT_norm ~ MIT_conc_nom, data = range_TA_MIT)
# Check the residuals
resGC_MIT <- resid(lm_GC_MIT)
# homoscedasticity ?
plot(range_TA_MIT$MIT_norm, resGC_MIT,
     ylab="residuals", xlab="MITnorm AUC")
abline(0,0)
bptest(lm_GC_MIT) #Breusch-Pagan test confirm heteroscedasticity
#
# # Regression With Robust Standard Errors (https://rpubs.com/cyobero/187387) 
# ct_mit <- coeftest(lm_GC_MIT, vcov = vcovHC(lm_GC_MIT, "HC1"))   # HC1 gives the White standard errors
# ct_mit[,2] #Robust standard errors
# #
# Generalized Least Squares With Unknown Form of Variance
lm_MIT.ols <- lm(MIT_conc_nom ~ MIT_norm, data = range_TA_MIT)
range_TA_MIT$resi <- lm_MIT.ols$residuals
varfunc.MIT.ols <- lm(log(resi^2) ~ log(MIT_norm), 
                          data = subset(range_TA_MIT,type=="range")) #subset as log(0) is impossible
#c(summary(varfunc.MIT.ols)$coef[1], summary(varfunc.MIT.ols)$coef[2]) #least squares estimate for the variance function
# Weighted regression :
# transform the observations in such a way that the transformed model has a constant error variance
range_TA_MIT$varfunc <- c(exp(varfunc.MIT.ols$fitted.values),rep(1,10))
lm_MIT.gls <- lm(MIT_conc_nom ~ MIT_norm , weights = 1/sqrt(varfunc) , data = range_TA_MIT)
#summary(lm_MIT.gls)
#hist(lm_MIT.gls$residuals)
#plot(lm_MIT.gls,2)
#plot(lm_MIT.gls,3)
#
ggplot(range_TA_MIT, aes(x=MIT_conc_nom, y=MIT_norm))+
  expand_limits(y = 0)+
  geom_point()+
  geom_smooth(method='lm')+
  geom_text(x = 2, y = 3, label = paste("Adj R2 = ",signif(summary(lm_MIT.gls)$adj.r.squared, 3),
                                             "Intercept =",signif(summary(lm_MIT.gls)$coef[[1]], 3) ,
                                             " Slope =",signif(summary(lm_MIT.gls)$coef[[2]], 3),
                                             " P =",signif(summary(lm_MIT.gls)$coef[2,4],3)))+
  labs(title="MIT acute toxicity: calibration curve", y="Normalized AUC", x="Nominal concentration (mg/L)")+
  theme_classic()
#
## Predictions ----
pred_MIT <- predict(lm_MIT.gls, data.frame(MIT_norm = data_SI2$MIT_norm))
data_SI2$pred_MIT <- pred_MIT
data_TA_MIT <- data_SI2[which(data_SI2$type=="acute"),]
# #view : Area Under Curve normalized by Intern Standard vs nominal concentration
# ggplot(data_TA_MIT, aes(x=as.factor(MIT_conc_nom), y=MIT_norm,
#                            fill=as.factor(time)))+
#   geom_boxplot()+
#   labs(title="Acute toxicity of MIT: normalized AUC", y="", x="Nominal concentration (mg/L)",
#        fill="day")+
#   theme_classic()
#view : predicted concentration vs nominal concentration
ggplot(data_TA_MIT, aes(x=MIT_conc_nom, y=pred_MIT,
                          group=MIT_conc_nom, col=as.factor(time)))+
  geom_jitter()+
  labs(title="MIT acute toxicity: concentration (mg/L)", y="Effective concentration", 
       x="Nominal concentration", col="day")+
  theme_classic()+
  expand_limits(x = 0, y = 0)+
  geom_abline(slope = 1)
#view : comparison of 0 vs. 1 day
ggplot(data_TA_MIT, aes(x=MIT_conc_nom, y=pred_MIT,
                          fill=as.factor(time))) + 
  geom_boxplot() +
  facet_wrap(~MIT_conc_nom, scale="free")+
  labs(title="Acute toxicity of MIT: Effective concentration (mg/L)", y="",
       x="Nominal concentration (mg/L)", fill="day")
#
## Concentration stability over 24h ? ----
# Test day 0 vs. day 1 
#
# qplot(sample = pred_MIT, data = data_TA_MIT,
#       shape=as.factor(time))
# hist(data_TA_MIT$pred_MIT[which(data_TA_MIT$time==0)])
# mean(data_TA_MIT$pred_MIT)
# var(data_TA_MIT$pred_MIT)
glm.time.m <- glm(pred_MIT ~ MIT_conc_nom * as.factor(time), 
                    data = data_TA_MIT, family = Gamma(link = identity)) #Gamma(link = identity)
# plot(glm.time.m,2)
# shapiro.test(residuals(glm.time.m))
# plot(glm.time.m,3)
# hist(glm.time.m$residuals)
# plot(data_TA_MIT$pred_MIT, glm.time.m$residuals,
#      ylab="residuals", xlab="MIT conc predictions")
# abline(0,0)
summary(glm.time.m)
#
## Measures and degradation over 24h ----
# is it >20% ?
med_all <- aggregate(data_TA_MIT$pred_MIT, list(data_TA_MIT$MIT_conc_nom), median)
sd_all <- cbind(med_all, aggregate(data_TA_MIT$pred_MIT, list(data_TA_MIT$MIT_conc_nom), sd)[,2])
med_j0 <- aggregate(subset(data_TA_MIT, time==0)$pred_MIT, list(subset(data_TA_MIT, time==0)$MIT_conc_nom), median)
sd_j0 <- cbind(med_j0, aggregate(subset(data_TA_MIT, time==0)$pred_MIT, list(subset(data_TA_MIT, time==0)$MIT_conc_nom), sd)[,2])
med_j1 <- aggregate(subset(data_TA_MIT, time==1)$pred_MIT, list(subset(data_TA_MIT, time==0)$MIT_conc_nom), median)
sd_j1 <- cbind(med_j1, aggregate(subset(data_TA_MIT, time==1)$pred_MIT, list(subset(data_TA_MIT, time==0)$MIT_conc_nom), sd)[,2])
#
table <- cbind(sd_all, sd_j0[,c(2,3)], sd_j1[,c(2,3)])
names(table) <- c("nominative", "median effective","sd tot", 
                    "median day 0", "sd 0",
                    "median day 1", "sd 1")
table <- table %>% mutate_if(is.numeric, round, digits=3)
table$decrease <- 100*((table$`median day 0` - table$`median day 1`)/table$`median day 0`)
# validated <20%
#
# end of script
