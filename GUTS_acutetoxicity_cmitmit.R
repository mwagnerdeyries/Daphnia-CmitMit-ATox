# 
# Code and data for the toxicokinetics-toxicodynamics analysis of bioassay survival data 
# in *(ref paper)
# https://github.com/mwagnerdeyries/Daphnia-CmitMit-ATox
# author : Margot Wagner Deyries ma.ju.wagner@gmail.com
# morse package tutorial https://cran.r-project.org/web/packages/morse/vignettes/tutorial.html#guts-model-with-constant-exposure-concentrations
rm(list=ls())
# 
################################################################################
# Packages
################################################################################

library(morse)
library(tidyverse)
library(ggplot2)

################################################################################
# Setting seed and reading data
################################################################################
set.seed(999)
setwd("C:/Users/mgwagner/Documents/INRA/ToxAigue/article") #get rid
survival_data <- read.csv("cmitmit_surv_data_effectiveconc.csv", sep=";") 
#
line_names <- unique(survival_data$Clone)
molecule <- unique(survival_data$molecule)
#
# Preparing the input data
# toy #get rid or comment all
C = 'LA0' # clone 
M = 'C' # molecule 
data_line <- subset(survival_data, Clone == C & molecule == M)
ATox <- survData(data_line)

################################################################################
# Visualize data set
################################################################################
summary(ATox)
plotDoseResponse(ATox, target.time = 2, addlegend = TRUE)
plot(ATox, pool.replicate=FALSE)

################################################################################
# Fit an exposure-response model to the survival data at target time
################################################################################
# probabilistic model with survival rate as a log-logistic function of concentration in chemical compound (c):
# f(c) = d / (1 + (c/e)^b )
# parameters:
# d as the survival rate in absence of contaminant
# e LC50
# b related to the effect intensity of the contaminant
fit <- survFitTT(ATox,
                 target.time = 2,
                 lcx = c(10, 20, 30, 40, 50))
summary(fit)
plot(fit, log.scale = TRUE, adddata = TRUE,   addlegend = TRUE)
# validate the model with a posterior predictive check
ppc(fit)

################################################################################
# Fit a toxicokinetic-toxicodynamic model:
# General Unified Threshold model for Survival (GUTS)
################################################################################
####~~~~ Model GUTS-RED-IT ~~~~####
fit.cstIT <- survFit(ATox, model_type="IT")
summary(fit.cstIT, EFSA_name = TRUE, quiet=TRUE)$Qpost
# parameters:
# kd the dominant toxicokinetic rate constant (SD and IT)
# hb the background hazard rate (SD and IT)
# beta [shape of threshold effect distribution with median mw] (IT)
# mw median of the threshold effect distribution (IT)
plot_prior_post(fit.cstIT)
plot(fit.cstIT)
ppc(fit.cstIT)
LCx.cstIT <- LCx(fit.cstIT, X = 50, time_LCx = 2)
plot(LCx.cstIT, cex.main=0.6, 
     main=paste('Concentration-response curve: LC',LCx.cstIT$X_prop_provided*100,' at 48h'),
     subtitle = paste(C, ' daphnia line and ', M, ' chemical'))
LCx.cstIT$df_LCx


################################################################################
# Lethal concentrations 
################################################################################
LCx.cstIT <- LCx(fit.cstIT, X = 50, time_LCx = 1)
plot(LCx.cstIT, cex.main=0.6, 
     main=paste('Concentration-response curve: LC',LCx.cstIT$X_prop_provided*100,' at 24h'),
     subtitle = paste(clo, ' daphnia line and ', mol, ' chemical'))
LCx.cstIT$df_LCx
LCx.cstIT <- LCx(fit.cstIT, X = 50, time_LCx = 2)
plot(LCx.cstIT, cex.main=0.6, 
     main=paste('Concentration-response curve: LC',LCx.cstIT$X_prop_provided*100,' at 48h'),
     subtitle = paste(clo, ' daphnia line and ', molecule, ' chemical'))
# predicted concentration-response curve
LCx.cstIT$df_LCx

################################################################################
# Model validation with posterior predictive check
################################################################################
# Check posterior
plot_prior_post(fit.cstIT, EFSA_name = TRUE)
# grey=prior, orange=posterior
# Survival probabilities
plot(fit.cstIT)
# Posterior predictive check plot
ppc(fit.cstIT)
# Parameters correlation check
tmp <- as.data.frame(as.matrix(fit.cstIT$mcmc))
mcmctot <- tmp[,c("kd_log10","hb_log10","alpha_log10","beta_log10")]
colnames(mcmctot) <- c("kD_log10","hb_log10","alpha_log10","beta_log10")
pairs(mcmctot, upper.panel = NULL)

################################################################################
# Compute prediction with GUTS model
################################################################################
# (1) upload or build a data frame with the exposure profile
# argument `replicate` is used to provide several profiles of exposure
data_4prediction <- data.frame(time = c(1:20, 1:20),
                               conc = c(rep(c(0.0632, 0.0555, 0.0475, 0.0444, 0.0339),4),
                                        rep(860e-6,20)),
                               replicate = c(rep("60µg/L deg 5 days", 20), 
                                             rep("860ng/L cst", 20)))

# (2) Use the fit on constant exposure 
predict_PRZ_cstIT_4pred <- predict(object = fit.cstIT, data_predict = data_4prediction)
# (3) Plot the predicted survival rate under the new exposure profiles.
plot(predict_PRZ_cstIT_4pred)
# warning :  'You should try the function 'predict_ode()' which is much more robust but longer to compute.'
# Robust ODE solver with deSolve: (with 1000 MCMC)
predict_PRZ_cstIT_4pred_ode <- predict_ode(object = fit.cstIT, data_predict = data_4prediction)
plot(predict_PRZ_cstIT_4pred_ode)

## Multiplication factors: ‘margin of safety’ ----
data_4MFx <- data.frame(time = 1:20,
                        conc = rep(860e-6,20))
data_4MFx <- data.frame(time = 1:20,
                        conc = rep(c(0.0632, 0.0555, 0.0475, 0.0444, 0.0339),4)) #degradation from 60µg/L
MFx_PRZ_cstIT_4MFx <- MFx(object = fit.cstIT, data_predict = data_4MFx, ode = TRUE)#default 50%
plot(MFx_PRZ_cstIT_4MFx)
plot(MFx_PRZ_cstIT_4MFx, x_variable =  "Time") #vizualize differences in survival rate with and without the multiplication factor
MFx_PRZ_cstIT_4MFx$df_MFx

# change target time
MFx_PRZ_cstIT_4MFx_x10 <- MFx(object = fit.cstIT, data_predict = data_4MFx, X = 10, threshold_iter = 20)
plot(MFx_PRZ_cstIT_4MFx_x10, log_scale = TRUE)
MFx_PRZ_cstIT_4MFx_x10$df_MFx
plot(MFx_PRZ_cstIT_4MFx_x10, x_variable =  "Time")

################################################################################
# Clones comparison predictions
################################################################################
GUTS_predict_860 <- data.frame(concentration = numeric(),
                           q50 = numeric(),
                           qinf = numeric(),
                           qsup = numeric(),
                           line = character(),
                           molecule = character(),
                           time = numeric())
GUTS_predict_120 <- data.frame(concentration = numeric(),
                               q50 = numeric(),
                               qinf = numeric(),
                               qsup = numeric(),
                               line = character(),
                               molecule = character(),
                               time = numeric())
set.seed(999)
for (C in line_names){
  cat(paste("\n**** Daphnia Line ",C,' START****\n'))
  for (M in molecule){
    tryCatch({
      #data
      data <- subset(survival_data, Clone == C & molecule == M)
      ATox <- survData(data)
      #GUTS-IT model
      fit.cstIT <- survFit(ATox, model_type="IT")
      #fit.cstIT <- survFit(ATox, model_type="IT",hb_value = FALSE, hb_valueFORCED = 0.01) # with fixed background mortality to 10%
      if(M=="M"){
        # data_4prediction <- data.frame(time = c(1:20),
        #                                conc = rep(860e-6,20),
        #                                replicate = rep("MIT 860ng/L cst", 20))
        data_4prediction <- data.frame(time = c(1:20),
                                       conc = rep(860e-3,20),
                                       replicate = rep("MIT 860µg/L cst", 20))
      }else{
        # data_4prediction <- data.frame(time = c(1:20),
        #                                conc = rep(120e-6,20),
        #                                replicate = rep("CMIT 120ng/L cst", 20))
        data_4prediction <- data.frame(time = c(1:20),
                                       conc = rep(120e-3,20),
                                       replicate = rep("CMIT/MIT 120µg/L cst", 20))
      }
      ##
      ## with background mortality fixed at 0.01 for everyone (10%) because very high for MIT assay in some lines
      # predict_PRZ_cstIT_4pred_ode <- predict_ode(object = fit.cstIT, data_predict = data_4prediction,
      #                                               hb_value = FALSE)#, hb_valueFORCED = 0.01)
      predict_PRZ_cstIT_4pred_ode <- predict_ode(object = fit.cstIT, data_predict = data_4prediction)
      
      #predicted time-response curve data
      preddata <- predict_PRZ_cstIT_4pred_ode$df_quantile
      temp <- data.frame(concentration = preddata$conc,
                         q50 = preddata$q50,
                         qinf = preddata$qinf95,
                         qsup = preddata$qsup95,
                         line = C,
                         molecule = M,
                         time = preddata$time)
      if(M=="M"){
        GUTS_predict_860 <- rbind(GUTS_predict_860, temp)
      }else{
        GUTS_predict_120 <- rbind(GUTS_predict_120, temp)
      }
    }, error=function(e){cat("ERROR : ",conditionMessage(e),"\n")})
  cat(paste("\n**** Molecule ",M,' END****\n'))
  }
  cat(paste("\n**** Daphnia Line ",C,' END****\n'))
}
#view predictions
#MIT
ggplot(GUTS_predict_860, aes(x=time, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', 
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.05)+
  labs(title='Predicted survival rate for 20 days under constant exposure of MIT at 860mg/L', x ="Time (days)", y = "",
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")+
  expand_limits(y = 0)
#CMIT
ggplot(GUTS_predict_120, aes(x=time, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', 
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.05)+
  labs(title='Predicted survival rate for 20 days under constant exposure of CMIT/MIT at 120mg/L', x ="Time (days)", y = "",
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")+
  expand_limits(y = 0)

################################################################################
# Clones comparison
################################################################################

GUTS_all_LC50 <- data.frame(median = numeric(),
                       quantile25 = numeric(),
                       quantile975 = numeric(),
                       line = character(),
                       molecule = character(),
                       nb_replicate = integer())
GUTS_LCxdata <- data.frame(concentration = numeric(),
                      q50 = numeric(),
                      qinf = numeric(),
                      qsup = numeric(),
                      line = character(),
                      molecule = character())
GUTS_Model_parameters <- data.frame(parameters = character(),
                               median = numeric(),
                               Q2.5 = numeric(),
                               Q97.5 = numeric(),
                               line = character(),
                               molecule = character())
logit_all_LC50 <- data.frame(median = numeric(),
                            quantile25 = numeric(),
                            quantile975 = numeric(),
                            line = character(),
                            molecule = character(),
                            TT = numeric())
logit_Model_parameters <- data.frame(parameters = character(),
                                    median = numeric(),
                                    Q2.5 = numeric(),
                                    Q97.5 = numeric(),
                                    line = character(),
                                    molecule = character(),
                                    TT = numeric())

# Compute models for each line and both set of molecules
set.seed(999)
for (C in line_names){
  cat(paste("\n**** Daphnia Line ",C,' START****\n'))
  #
  for (M in molecule){
    cat(paste("\n**** Molecule ",M,' START****\n'))
    data <- subset(survival_data, Clone == C & molecule == M)
    ATox <- survData(data)
    tryCatch({
    for (T in c(1:2)){ 
      # time target (1 or 2 days)
      ## log-logistic binomial model
      fit.logit <- survFitTT(ATox,target.time = T,lcx = c(50))
      # parameters log-logit
      params.logit <- cbind(fit.logit$estim.par, C, M, T)
      logit_Model_parameters <- rbind(logit_Model_parameters, params.logit)
      # Lethal concentration log-logit
      temp3 <- data.frame(median = fit.logit$estim.LCx$median,
                          quantile25 = fit.logit$estim.LCx$Q2.5,
                          quantile975 = fit.logit$estim.LCx$Q97.5,
                          line = C,
                          molecule = M,
                          TT = T)
      logit_all_LC50 <- rbind(logit_all_LC50, temp3)
    }
    #
    
    ## GUTS-RED-IT
    fit.cstIT <- survFit(ATox, model_type="IT")
    # Parameters estimations
    param.cstIT <- cbind(summary(fit.cstIT, EFSA_name = TRUE, quiet=TRUE)$Qpost, C, M)
    GUTS_Model_parameters <- rbind(GUTS_Model_parameters,param.cstIT)
    # Lethal concentrations
    LCx.cstIT <- LCx(fit.cstIT, X = 50, time_LCx = 2)
    temp <- data.frame(median = LCx.cstIT$df_LCx[1,2],
                       quantile25 = LCx.cstIT$df_LCx[2,2],
                       quantile975 = LCx.cstIT$df_LCx[3,2],
                       line = C,
                       molecule = M,
                       nb_replicate = min(unique(summary(ATox,quiet=TRUE)$NbrepTimeConc)))
    GUTS_all_LC50 <- rbind(GUTS_all_LC50, temp)
    # record median distribution along product concentration
    temp2 <- data.frame(concentration = LCx.cstIT$df_dose$concentration,
                        q50 = LCx.cstIT$df_dose$q50,
                        qinf = LCx.cstIT$df_dose$qinf95,
                        qsup = LCx.cstIT$df_dose$qsup95,
                        line = C,
                        molecule = M)
    GUTS_LCxdata <- rbind(GUTS_LCxdata, temp2)
    #
    }, error=function(e){cat("ERROR : ",conditionMessage(e),"\n")})
    cat(paste("\n**** Molecule ",M,' END****\n'))
  }
  cat(paste("\n**** Daphnia Line ",C,' END****\n'))
}

# 
# # Clones and molecule names matches
# ref_names <- data.frame(C = c("AL0", "GO6", "LA0", "PE7", "P16", "RE0", "SE2", "SE5"),
#                         Clones = c("Alsace", "Goven6", "Lassalle", "Pearl7", 
#                           "Pearl16", "Rennes", "Sene2", "Sene5"))

####~~~~ Dose-response curve for all lines with GUTS models ~~~~####
ggplot(subset(GUTS_LCxdata, molecule=="C"), aes(x=concentration, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.15)+
  labs(title='Effect of CMIT+MIT: Survival at 48h', 
       x ="CMIT + MIT concentration (mg/L)", y = "", 
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")
  

ggplot(subset(GUTS_LCxdata, molecule=="M"), aes(x=concentration, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', 
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.15)+
  labs(title='Effect of MIT: Survival at 48h', x ="MIT concentration (mg/L)", y = "",
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")

################################################################################
# Test inter-clone variation vs inter-individual variation 
################################################################################
#### Permutation test ####
#### Comparing several lines data with permutation tests ####
#
# test1 : takes all proportion data at time 2 : compare differences in median of global distribution 
# nb model : Prop ~ conc + Clone
# the proportions of survival are permutated knowing the concentrations
# 
# test 2 : knowing concentration, switch Clone, 
# then fit TT model for each clone with Nsurv and compute new LC50
# compare differences in LC50 medians of the selection of Clones
#
library(morse)
library(dplyr)
setwd("C:/Users/mgwagner/Documents/INRA/ToxAigue/article") #get rid
# use proportion of survival data
survival_data <- read.csv("cmitmit_surv_data_effectiveconc_prop.csv", sep=";")
line_names <- unique(survival_data$Clone)
#
# brute force way
# uses raw data and compare each couple of lines
#### Functions ####
# Get all combinations of 2 elements from a list of objects
combinazzion <- function(linenames){
  combi <- combn(linenames,2)
  t(combi)
}
#
# simulates a randomized data set 
simfun_rsamp <- function(data) {
  # transform(data, Propsurv=sample(Propsurv)) # randomize Propsurv without caring about concentration specific response
  clones <- unique(data$Clone)
  # get common tested concentrations only
  common_conc <- Reduce(intersect, list(data$conc[which(data$Clone==clones[1])],
                                        data$conc[which(data$Clone==clones[2])]))
  randomised_data <- data.frame(Clone = character(),
                      conc = numeric(),
                      Propsurv = numeric())
  for(i in common_conc){
    rand <- transform(subset(data, conc==i), Propsurv=sample(Propsurv))
    rand_min <- data.frame(rand$Clone, rand$conc, rand$Propsurv)
    randomised_data <- rbind(randomised_data, rand_min)
  }
  names(randomised_data) <- c("Clone", "conc", "Propsurv")
  randomised_data
}
# takes a simulated data set and returns whatever summary statistic we want
####-- Stat1 : Here the difference of medians is taken --####
difffun <- function(data) {
  daphnia <- as.character(unique(data$Clone))
  median(data[data$Clone==daphnia[1],"Propsurv"])-
    median(data[data$Clone==daphnia[2],"Propsurv"])
}
# Action 1
stat1 <- function(data){
  set.seed(777)
  permdist_rsamp <- replicate(9999,difffun(simfun_rsamp(data)))
  daphnia <- as.character(unique(data$Clone))
  obs <- median(data[data$Clone==daphnia[1],"Propsurv"])-
    median(data[data$Clone==daphnia[2],"Propsurv"])
  permdist_rsamp <- c(permdist_rsamp,obs)
}

# results of action 1
get_stat1 <- function(statdata, graphics=FALSE, daphnia_lines){
  res <- stat1(statdata)
  obs <- res[10000]
  # graphical output
  if(graphics==TRUE){
    hist(res, breaks=20, main = paste('Daphnia lines :',daphnia_lines[1], 'and', daphnia_lines[2]))
    abline(v=obs,col='red')
    abline(v=quantile(res, 0.975), lty=2, col='blue')
    abline(v=-quantile(res, 0.975), lty=2, col='blue')
  }
  # percentile in random distribution corresponding to observation
  # 'statistically different' if <2.5% or >97.5% 
  #
  # get the fraction of observations les or equal to threshold 'obs'
  percentile <- ecdf(res)
  percentile(obs)
}
#
#
simfun_rsamp2 <- function(data){
  daphnia <- as.character(unique(data$Clone))
  # get common tested concentrations only
  common_conc <- Reduce(intersect, list(data$conc[which(data$Clone==daphnia[1])],
                                        data$conc[which(data$Clone==daphnia[2])]))
  randomised_data <- data.frame(Clone = character(),
                                conc = numeric(),
                                time = numeric(),
                                Nsurv = numeric(),
                                Ninit = numeric(),
                                replicate = numeric())
  # randomise (Nsurv+replicate)by Nsurv for TargetTime == 2
  for(i in common_conc){ 
    conc_data <- subset(data, conc==i)
    conc_data <- mutate(conc_data,index=paste(Clone,replicate)) #keep index information
    rand <- transform(subset(conc_data, time==2), Clone=sample(Clone)) #randomise Clone at time==2
    # rbind by time (otherwise morse::survDataCheck isnt happy)
    rand_02 <- rbind(subset(conc_data, time==0), 
                     subset(conc_data, time==1),
                     rand)      
    # extend new clone to replicate at t==0 with the help of the index
    rand_t02 <- rand_02[FALSE,]
    for(q in unique(rand_02$index)){
      duo <- subset(rand_02, index == q)
      duo$Clone[which(duo$time==0)] <- duo$Clone[which(duo$time==2)]
      duo$Clone[which(duo$time==1)] <- duo$Clone[which(duo$time==2)]
      rand_t02 <- rbind(rand_t02, duo)
    }
    #
    rand_min <- data.frame(rand_t02$Clone, rand_t02$conc, rand_t02$time, rand_t02$Nsurv, 
                           rand_t02$Ninit, rand_t02$replicate)
    randomised_data <- rbind(randomised_data, rand_min)
  }
  names(randomised_data) <- c("Clone", "conc", "time", "Nsurv","Ninit","replicate")
  # new replicates without duplication
  randomised_data$replicate <- as.character(rep(sample(nrow(randomised_data)/3), times = 1, each = 3))
  # output
  randomised_data
}
#
difffun2 <- function(data){
  daphnia <- as.character(unique(data$Clone))
  Atox1 <- survData(data[data$Clone==daphnia[1],])
  fit1 <- survFitTT(Atox1, target.time = 2, lcx = c(50))
  Atox2 <- survData(data[data$Clone==daphnia[2],])
  fit2 <- survFitTT(Atox2, target.time = 2, lcx = c(50))
  # differences in LC50 medians of the selection of Clones
  fit1$estim.LCx$median - fit2$estim.LCx$median
}
#
stat2 <- function(data){
  set.seed(777)
  tryCatch({
    permdist_rsamp <- replicate(99,difffun2(simfun_rsamp2(data)))
    daphnia <- as.character(unique(data$Clone))
    Atox1 <- survData(data[data$Clone==daphnia[1],])
    fit1 <- survFitTT(Atox1, target.time = 2, lcx = c(50))
    Atox2 <- survData(data[data$Clone==daphnia[2],])
    fit2 <- survFitTT(Atox2, target.time = 2, lcx = c(50))
  }, error=function(e){cat("ERROR : ",conditionMessage(e),"\n", "clones ",daphnia, "\n")
                        fit1 <- NA
                        fit2 <- NA})
  
  # differences in LC50 medians of the selection of Clones
  if(is.na(fit1)){
    obs <- NA
  }else{
    obs <- fit1$estim.LCx$median - fit2$estim.LCx$median
  }
  permdist_rsamp <- c(permdist_rsamp,obs)
}

# fit survival response and compute LC50 for observed and Clone-randomised data set
get_stat2 <- function(statdata){
  res <- stat2(statdata)
  obs <- res[100]
  # get the fraction of observations les or equal to threshold 'obs'
  percentile <- ecdf(res)
  percentile(obs)
}

# Prep percentile results dataframe
res_permtest <- data.frame(line1 = combinazzion(line_names)[,1],
                           line2 = combinazzion(line_names)[,2],
                           percentile_diffmedians = numeric(nrow(combinazzion(line_names))))
#
# select molecule and target time for Propsurv comparisons
# mydata <- subset(survival_data, molecule=="M" & time=="2")
#
mydata <- subset(survival_data, molecule=="C")
#
#
for(i in 1:ncol(combn(line_names,2))){ # data subset repetition
  
  
  combi <- combn(line_names,2)[,i]
  datadata <- rbind(mydata[which(mydata$Clone==combi[1]),],
                    mydata[which(mydata$Clone==combi[2]),])
  res_permtest$percentile_diffmedians[i] <- get_stat1(datadata, graphics=FALSE, 
                                                     daphnia_lines = combi)
  tryCatch({
    res_permtest$percentile_diffLC50[i] <- get_stat2(datadata)
  }, error=function(e){cat("ERROR : ",conditionMessage(e),"\n")})
}
res_permtest
# differences of medians 'statistically different at p<0.10' if <5% or >95%
# differences of medians 'statistically different at p<0.05' if <2.5% or >97.5%
# differences of medians 'statistically different at p<0.01' if <0.5% or >99.5% 
# differences of medians 'statistically different at p<0.005' if <0.25% or >99.75% 
# the cutest way to explain permutation test : https://www.jwilber.me/permutationtest/
#
# visualize LC50
library(ggplot2)
lc50_data <- read.csv2("cmitmit_results_lc50.csv")
# MIT logit
p1<- ggplot(subset(lc50_data, molecule=="M" & model=="logit"), aes(x=Clone, y=median, color=Clone)) + 
  geom_point()+
  geom_errorbar(aes(ymin=Qinf, ymax=Qsup), width=.2,
                position=position_dodge(0.05))
p1 +labs(title="LC50 at 48 for MIT, with logit model", x="", y = "LC50 in mg/L")+
  theme_classic() 
#CMIT/MIT logit
p2<- ggplot(subset(lc50_data, molecule=="C" & model=="logit"), aes(x=Clone, y=median, color=Clone)) + 
  geom_point()+
  geom_errorbar(aes(ymin=Qinf, ymax=Qsup), width=.2,
                position=position_dodge(0.05))
p2 +labs(title="LC50 at 48 for CMIT/MIT, with logit model", x="", y = "LC50 in mg/L")+
  theme_classic() 
# MIT guts
p3<- ggplot(subset(lc50_data, molecule=="M" & model=="guts"), aes(x=Clone, y=median, color=Clone)) + 
  geom_point()+
  geom_errorbar(aes(ymin=Qinf, ymax=Qsup), width=.2,
                position=position_dodge(0.05))
p3 +labs(title="LC50 at 48 for MIT, with guts model", x="", y = "LC50 in mg/L")+
  theme_classic() 
#CMIT/MIT guts
p4<- ggplot(subset(lc50_data, molecule=="C" & model=="guts"), aes(x=Clone, y=median, color=Clone)) + 
  geom_point()+
  geom_errorbar(aes(ymin=Qinf, ymax=Qsup), width=.2,
                position=position_dodge(0.05))
p4 +labs(title="LC50 at 48 for CMIT/MIT, with guts model", x="", y = "LC50 in mg/L")+
  theme_classic() 


################################################################################
# Dose-response ignoring lines
################################################################################
M <- "C"
M <- "M"
data_mol <- subset(survival_data, molecule == M)
data_mol$replicate <-  as.character(rep(sample(nrow(data_mol)/3), times = 1, each = 3))
ATox <- survData(data_mol)
summary(ATox)
plotDoseResponse(ATox, target.time = 2, addlegend = TRUE)
plot(ATox, pool.replicate=FALSE)
#
fit <- survFitTT(ATox,
                 target.time = 2,
                 lcx = c(10, 20, 30, 40, 50))
summary(fit)
plot(fit, log.scale = TRUE, adddata = TRUE,   addlegend = TRUE)
# validate the model with a posterior predictive check
ppc(fit)
#
fit.cstIT <- survFit(ATox, model_type="IT")
summary(fit.cstIT, EFSA_name = TRUE, quiet=TRUE)$Qpost
# parameters:
# kd the dominant toxicokinetic rate constant (SD and IT)
# hb the background hazard rate (SD and IT)
# beta [shape of threshold effect distribution with median mw] (IT)
# mw median of the threshold effect distribution (IT)
plot_prior_post(fit.cstIT)
plot(fit.cstIT)
ppc(fit.cstIT)
LCx.cstIT <- LCx(fit.cstIT, X = 50, time_LCx = 2)
plot(LCx.cstIT, cex.main=0.6, 
     main=paste('Concentration-response curve: LC',LCx.cstIT$X_prop_provided*100,' at 48h'),
     subtitle = paste(M, ' chemical'))
LCx.cstIT$df_LCx
# Parameters correlation check
tmp <- as.data.frame(as.matrix(fit.cstIT$mcmc))
mcmctot <- tmp[,c("kd_log10","hb_log10","alpha_log10","beta_log10")]
colnames(mcmctot) <- c("kD_log10","hb_log10","alpha_log10","beta_log10")
pairs(mcmctot, upper.panel = NULL)

head(All_M$df_dose)
head(All_C$df_dose)

All_M$df_LCx
All_C$df_LCx

ggplot(All_M$df_dose, aes(x=concentration, y=q50))+
  geom_line()+
  geom_ribbon(aes(ymin=qinf95, ymax=qsup95), linetype=2, alpha=0.25)+
  theme_classic()

ggplot(All_C$df_dose, aes(x=concentration, y=q50))+
  geom_line()+
  geom_ribbon(aes(ymin=qinf95, ymax=qsup95), linetype=2, alpha=0.25)+
  theme_classic()

# comp params
guts <- read.csv2("GUTSparams.csv")
guts$median <- as.numeric(guts$median)
guts$Q2.5 <- as.numeric(guts$Q2.5)
guts$Q97.5 <- as.numeric(guts$Q97.5)
#
ggplot(subset(guts, parameters=="kD"), 
       aes(x = C, y = median, ymin = Q2.5, ymax = Q97.5, shape=M, color=M))+
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width=.4,
                position=position_dodge(.9))+ 
  geom_point(size=1.5, position=position_dodge(.9))+
  expand_limits(y = 0)+
  scale_color_manual(values=c('royalblue','lightsteelblue'))+
  theme_classic()+
  labs(title="Dominant rate (kD) parameter from GUTS models", x="", y="")
library(ggpubr)
library(rstatix)
kD <- subset(guts, parameters=="kD" & C!="ALL")
kD %>%
  group_by(M) %>%
  get_summary_stats(median, type = "mean_sd")
#outliers
kD %>%
  group_by(M) %>%
  identify_outliers(median)
# normality
kD %>%
  group_by(M) %>%
  shapiro_test(median)
ggqqplot(kD, x = "median", facet.by = "M")
# variance equality
kD %>% levene_test(median ~ as.factor(M))
# all good for a classic student t test !
stat.test <- kD %>% 
  t_test(median ~ M, var.equal = TRUE) %>%
  add_significance()
stat.test
# effect size, Cohen's d
kD %>%  cohens_d(median ~ M, var.equal = TRUE)

#
logit <- read.csv2("logitparams.csv")
logit$median <- as.numeric(logit$median)
logit$Q2.5 <- as.numeric(logit$Q2.5)
logit$Q97.5 <- as.numeric(logit$Q97.5)
#
ggplot(subset(logit, T==2 & X=="b"), aes(x = C, y = median, ymin = Q2.5, ymax = Q97.5, shape=M, color=M))+
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width=.4,
                position=position_dodge(.9))+ 
  geom_point(size=1.5, position=position_dodge(.9))+
  expand_limits(y = 0)+
  scale_color_manual(values=c('royalblue','lightsteelblue'))+
  theme_classic()+
  labs(title="Effect intensity of the contaminant (b) parameter from logit models", x="", y="")
logit <- subset(logit, T==2 & C!="ALL" & X=="b")
logit %>% 
  t_test(median ~ M, var.equal = TRUE) %>%
  add_significance()
