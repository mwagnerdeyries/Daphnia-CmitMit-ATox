# 
# Code and data for the toxicokinetics-toxicodynamics analysis of bioassay survival data 
# in: Variation of tolerance to isothiazolinones among Daphnia pulex clones
# https://github.com/mwagnerdeyries/Daphnia-CmitMit-ATox
# author : Margot Wagner Deyries ma.ju.wagner@gmail.com
# morse package tutorial https://cran.r-project.org/web/packages/morse/vignettes/tutorial.html#guts-model-with-constant-exposure-concentrations
rm(list=ls())
# 
################################################################################
# Packages
################################################################################

# information for morse package installation at https://github.com/pveber/morse

library(morse)
library(tidyverse)
library(ggplot2)

################################################################################
# Setting seed and reading data
################################################################################
set.seed(999)
survival_data <- read.csv("cmitmit_surv_data_effectiveconc.csv", sep=";") 
#
line_names <- unique(survival_data$Clone)
molecule <- unique(survival_data$molecule)
#
################################################################################
# Study of individual Clone/Contaminant pair
################################################################################
# Preparing the input data
# toy data : chose clonal lineage and molecule
# clones : "AL0" "GO6" "LA0" "PE7" "P16" "RE0" "SE2" "SE5"
# CMIT+MIT:"C" or MIT alone:"M"
C = 'LA0' # clone 
M = 'C' # molecule 
data_line <- subset(survival_data, Clone == C & molecule == M)
ATox <- survData(data_line)

################################################################################
# Visualize data set
# Study of individual Clone/Contaminant pair
################################################################################
summary(ATox)
plotDoseResponse(ATox, target.time = 2, addlegend = TRUE)
plot(ATox, pool.replicate=FALSE)

################################################################################
# Fit an exposure-response model to the survival data at target time
# Study of individual Clone/Contaminant pair
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
# Study of individual Clone/Contaminant pair
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
# Study of individual Clone/Contaminant pair
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
# Study of individual Clone/Contaminant pair
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
# Clones comparison
# For the visualization of all Clone/Contaminant pairs
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

####~~~~ kD parameter comparison ~~~~####
# Figure 3
kD_param <- subset(GUTS_Model_parameters, parameters=="kD") #data subset
kD_param$C <- as.factor(kD_param$C)
kD_param$M <- as.factor(kD_param$M)
kD_param$median <- as.numeric(kD_param$median)
kD_param$Q2.5 <- as.numeric(kD_param$Q2.5)
kD_param$Q97.5 <- as.numeric(kD_param$Q97.5)

ggplot(kD_param, aes(x=factor(C, levels=c("RE0", "SE5", "PE7", "LA0", "P16", "AL0", "SE2", "GO6")), 
                     y=median, group=M))+
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5, color=M),
               position=position_dodge(width=1))+
  geom_point(aes(color=M, shape=M), position=position_dodge(width=1))+
  theme_classic() +
  scale_color_manual(values=c('purple','hotpink'))+
  labs(title="Dominant rate parameter median and 95% CI",
       x="Clone",y="kD (day-1)")

####~~~~ Dose-response curve for all lines with GUTS models ~~~~####
# Figure 4 Top
ggplot(subset(GUTS_LCxdata, molecule=="C"), aes(x=concentration, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.15)+
  labs(title='Effect of CMIT+MIT: Survival at 48h', 
       x ="CMIT + MIT concentration (mg/L)", y = "", 
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")
  
# Figure 5 Top
ggplot(subset(GUTS_LCxdata, molecule=="M"), aes(x=concentration, y=q50, group=line))+
  geom_line(aes(linetype=line, color=line), size=1.1)+
  theme(legend.position = 'top', 
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey60"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey60"))+
  geom_ribbon(aes(ymin=qinf, ymax=qsup, fill=line), linetype=2, alpha=0.15)+
  labs(title='Effect of MIT: Survival at 48h', x ="MIT concentration (mg/L)", y = "",
       fill="Clonal lineages", linetype= "Clonal lineages", color="Clonal lineages")
#
# end of script
