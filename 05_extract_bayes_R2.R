rm(list = ls())

# R Setup ----
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(ggplot2)
library(patchwork)
library(future)

# Source conditional smoooth functions
source('/ncf/hcp/data/analyses/myelin/code/gam_bootstrap/conditional_smooth_sample.R')

## Load HCP Data
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")

## Define response variables
nsub <- dim(hcpd_data)[1]
nreg <- 360
response_vars <- c(paste0('myelin_glasser_v', 1:nreg))

#####################################
# Load brm output for full model ----
#####################################
setwd('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/full_model')

spline_model_list <- lapply(1:length(response_vars), function(i){
  
  right_side <- '~ s(Age) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
  left_side <- response_vars[[i]]
  
  model_name <- paste0('fit_brm_', left_side)
  if(!file.exists(paste0(model_name, '.rds'))){
    fit_brm <- NULL
  } else {
    fit_brm <- readRDS(paste0(model_name, '.rds'))
  }
  return(fit_brm)
})


########################################
# Load brm output for reduced model ----
########################################
setwd('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/reduced_model')

reduced_model_list <- lapply(1:length(response_vars), function(i){
  
  right_side <- '~ Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
  left_side <- response_vars[[i]]
  
  model_name <- paste0('fit_brm_', left_side)
  if(!file.exists(paste0(model_name, '.rds'))){
    fit_brm <- NULL
  } else {
    fit_brm <- readRDS(paste0(model_name, '.rds'))
  }
  return(fit_brm)
})

########################
# Estimate Bayes R2 ----
########################
spline_bayes_r2 <- lapply(spline_model_list, brms::bayes_R2)
reduced_bayes_r2 <- lapply(reduced_model_list, brms::bayes_R2)

  
  
  
# Extract Partial Bayes R2 ----

## Full
full_rsq <- rep(NA, length(response_vars))
for (i in 1:length(response_vars)){
  full_rsq[i] <- spline_bayes_r2[[i]][1]
}
write.table(partial_rsq, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/bayes_r2/n628_GlasserMMP_myelin_sAge_full_bayes_r2.txt", row.names=FALSE, col.names=FALSE)

## Reduced
reduced_rsq <- rep(NA, length(response_vars))
for (i in 1:length(response_vars)){
  reduced_rsq[i] <- reduced_bayes_r2[[i]][1]
}
write.table(partial_rsq, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/bayes_r2/n628_GlasserMMP_myelin_sAge_reduced_bayes_r2.txt", row.names=FALSE, col.names=FALSE)

## Partial
partial_rsq <- rep(NA, length(response_vars))
for (i in 1:length(response_vars)){
  partial_rsq[i] <- full_rsq[i] - reduced_rsq[i]
}

hist(partial_rsq, col="slategray3")
write.table(partial_rsq, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/bayes_r2/n628_GlasserMMP_myelin_sAge_partial_bayes_r2.txt", row.names=FALSE, col.names=FALSE)

########################################################################################
########################################################################################

# Extract R2 from reduced model ----
setwd('/ncf/hcp/data/analyses/myelin/brm_output/networks/reduced_model')

# Load HCP Data
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_TFcorr_myelin_GlasserMMP_Yeo7.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/n1023_hcpd_UncorrectedStandardScores.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/n632_hcpd_TFcorr_myelin_GlasserMMP_Yeo7.csv")

nsub <- dim(hcpd_data)[1]

## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)
hcpd_data$scaled_numNavs_T1w <- scale(hcpd_data$numNavs_T1w, center = TRUE, scale = TRUE)
hcpd_data$scaled_numNavs_T2w <- scale(hcpd_data$numNavs_T2w, center = TRUE, scale = TRUE)
hcpd_data$scaled_HeadSize <- scale(hcpd_data$HeadSize, center = TRUE, scale = TRUE)

# Define response variables

# seven_nets <- c('Working_Memory', 'Processing_Speed', 'Flanker', 'Dimensional_Card_Sort', 'Cognition.Fluid.Composite.v1.1')
seven_nets <- c(paste0('Yeo', 1:7, '_myelin')) 


model_list <- lapply(1:length(seven_nets), function(i){
  right_side <- '~ Sex + Scanner + scaled_numNavs_T1w + scaled_numNavs_T2w + mean_cortical_B1 + scaled_HeadSize'
  left_side <- seven_nets[[i]]
  
  mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
  form <- as.formula(paste0(left_side, right_side))
  
  d <- model.frame(formula = mf_form, data = hcpd_data)
  
  model_name <- paste0('fit_brm_', left_side)
  if(!file.exists(paste0(model_name, '.rds'))){
    fit_brm <- NULL
  } else {
    fit_brm <- readRDS(paste0(model_name, '.rds'))
  }
  return(fit_brm)
})

reduced_bayes_r2 <- lapply(model_list, brms::bayes_R2)

# Extract Bayes R2 ----
reduced_rsq <- rep(NA, length(seven_nets))

for (i in 1:length(seven_nets)){
  reduced_rsq[i] <- reduced_bayes_r2[[i]][1]
}

# Calculate partial R2 for effect of interest ----
partial_rsq <- as.numeric(full_rsq - reduced_rsq)

hist(partial_rsq, col="slategray3")

# Export output
write.table(partial_rsq, "/ncf/hcp/data/analyses/myelin/brm_output/networks/partial_bayes_R2/n628_sAge_Yeo7_network_bayes_partial_R2.txt", row.names=FALSE, col.names=FALSE)
