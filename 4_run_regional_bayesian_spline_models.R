rm(list = ls())

# R Setup ----
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE) # only necessary if you're going to compile into a notebook
library(brms)
library(data.table)
library(ggplot2)
library(readr)

## Set working directory where fit_brm output will be stored
setwd("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/linear_model")

# setwd('/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model')

# brms setup ----

## Load functions for extracting posterior samples
source('/ncf/hcp/data/analyses/myelin/code/gam_bootstrap/conditional_smooth_sample.R')

## Functions to extract smooth posteriors
get_smooth_posterior <- function(fit_brm, res, eps, age_col = 'Age', smooths = 's(Age)', probs = c(0.025, 0.975)){
  d <- fit_brm$data
  smooth_values <- list(seq(min(d[, age_col]), max(d[, age_col]), length.out = res))
  names(smooth_values) <- age_col
  smooth_values_eps <- lapply(smooth_values, `+`, eps)
  c_smooths_1 <- conditional_smooth_sample(fit_brm, smooths = smooths, 
                                           smooth_values = smooth_values,
                                           probs = probs,
                                           spaghetti = TRUE)
  c_smooths_2 <- conditional_smooth_sample(fit_brm, smooths = smooths, 
                                           smooth_values = smooth_values_eps,
                                           probs = probs,
                                           spaghetti = TRUE)
  
  c_smooths_1_sample_data <- data.table(attr(c_smooths_1[[1]], 'spaghetti')) #the use of [[1]] here is not very extensible
  c_smooths_2_sample_data <- data.table(attr(c_smooths_2[[1]], 'spaghetti'))
  
  c_smooths_1_sample_data[, estimate__eps := c_smooths_2_sample_data[,estimate__]]
  c_smooths_1_sample_data[, deriv1 := (estimate__eps - estimate__) / eps]
  return(c_smooths_1_sample_data)
}
get_smooth_summaries <- function(posterior_samples, probs = c(.025, .975), age_col = 'Age', smooth_col = 'estimate__', deriv1_col = 'deriv1', sample_index = 'sample__'){
  smooth_post_summary <- posterior_samples[, list(est = median(get(smooth_col)),
                                                  'lower' = quantile(get(smooth_col), probs = probs[[1]]), 
                                                  'upper' = quantile(get(smooth_col), probs = probs[[2]])),
                                           by = age_col]
  
  deriv1_post_summary <- posterior_samples[, list('deriv1' = median(get(deriv1_col)), 
                                                  'lower' = quantile(get(deriv1_col), probs = probs[[1]]), 
                                                  'upper' = quantile(get(deriv1_col), probs = probs[[2]])),
                                           by = age_col]
  
  deriv1_age_at_max_post <- posterior_samples[, list('max_age' = Age[which(get(deriv1_col) == max(get(deriv1_col)))]),
                                              by = sample_index]
  age_at_max_med_deriv <- deriv1_post_summary[deriv1 == max(deriv1), get(age_col)]
  deriv1_diff_from_max_post <- 
    posterior_samples[, 
                      diff_from_max := get(deriv1_col) - get(deriv1_col)[get(age_col) == age_at_max_med_deriv],
                      by = sample_index]
  deriv1_age_at_max_post_summary <- deriv1_age_at_max_post[, list('est' = median(max_age),
                                                                  'lower' = quantile(max_age, probs = probs[[1]]), 
                                                                  'upper' = quantile(max_age, probs = probs[[2]]))]
  
  deriv1_diff_from_max_post_summary <- deriv1_diff_from_max_post[, list('diff_from_max' = median(diff_from_max), 
                                                                        'diff_lower' = quantile(diff_from_max, probs = probs[[1]]), 
                                                                        'diff_upper' = quantile(diff_from_max, probs = probs[[2]])),
                                                                 by = age_col]
  deriv1_diff_from_max_post_summary[, compared_to_max := dplyr::case_when(diff_upper < 0 ~ 'less_steep',
                                                                          diff_lower > 0 ~ 'steeper',
                                                                          TRUE ~ 'as_steep')]
  deriv1_diff_from_max_post_summary_cols <- c(age_col, 'compared_to_max', 'diff_from_max', 'diff_upper', 'diff_lower')
  deriv1_post_summary <- 
    deriv1_post_summary[deriv1_diff_from_max_post_summary[, ..deriv1_diff_from_max_post_summary_cols],
                        on = age_col]
  return(list(smooth = smooth_post_summary, derivative = deriv1_post_summary, max_deriv_age_posterior = deriv1_age_at_max_post))
}

#Age difference for finite differences
eps = 1e-07
#How many ages to get marginal posterior at -- doesn't need to be a lot, necessarily
res = 100
######################################################################################

# Load HCP data ----
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")

# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_April2021.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n630_hcpd_TFcorr_myelin_GlasserMMP_Yeo7_ColeAnt12.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n630_hcpd_TFcorr_myelin_Schaefer400x7_Yeo7_ColeAnt12.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/n630_hcpd_TFcorr_myelin_thickness_GlasserMMP_Yeo7_ColeAnt12.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_TFcorr_myelin_GlasserMMP_Yeo7_ColeAnt12.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_TFcorr_myelin_GlasserMMP_Yeo7.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/old/n632_hcpd_TFcorr_myelin_GlasserMMP_Yeo7.csv")
# hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/old/n1023_hcpd_UncorrectedStandardScores.csv")

## Read vNav motion metrics
# vnav_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n655_hcpd_vNav_motion_scores.csv")
# hcpd_data <- merge(hcpd_data, vnav_data, by="subject_id")

## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)
# hcpd_data$scaled_T1w_vnav_rms_score <- scale(hcpd_data$T1w_vnav_rms_score, center = TRUE, scale = TRUE)
# hcpd_data$scaled_T2w_vnav_rms_score <- scale(hcpd_data$T2w_vnav_rms_score, center = TRUE, scale = TRUE)
# hcpd_data$scaled_HeadSize <- scale(hcpd_data$HeadSize, center = TRUE, scale = TRUE)

##########################
## LOW MOTION THRESHOLD ##
##########################
## Create binary flag indicating if subject hit ceiling for vnav acquisition
## (this means they had high motion)
motion_idx <- which(hcpd_data$numNavs_T1w==196 | hcpd_data$numNavs_T2w==122)
hcpd_data$vnav_ceiling <- 0
hcpd_data$vnav_ceiling[motion_idx] <-  1

# hcpd_data <- subset(hcpd_data, hcpd_data$vnav_ceiling==0)
# dim(hcpd_data)
############################################################################


## Task index for parallelizing regional models
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

region_id <- task_id

set.seed(123)

# Fit Bayesian smoothing spline models for each brain region ----

right_side <- '~ Age + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'

# right_side <- '~ s(Age) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'

# right_side <- paste0("~ s(Age) + thickness_glasser_v", region_id," + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor")
# right_side <- '~ Age + s(Age, bs="tp", m=c(2,0)) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'


left_side <- paste0('myelin_glasser_v', region_id) 
# left_side <- paste0('myelin_schaefer_v', region_id) 

mf_form <- as.formula(paste0(left_side, gsub('s\\((.*)\\)', '\\1', right_side)))
form <- as.formula(paste0(left_side, right_side))

d <- model.frame(formula = mf_form, data = hcpd_data)

## Create row index for merging subject IDs
d_row_index <- as.numeric(row.names(d))
d$subject_id <- hcpd_data$subject_id[d_row_index]

## ADDED for reduced models (no Age term)
d$Age <- hcpd_data$Age[d_row_index]

agerange <- range(d$Age)

message(sprintf('Estimating derivative for %s', left_side))
message(sprintf('Setting resolution for derivative across ages %.0f-%.0f to %d increments of %.1f years each.', agerange[[1]], agerange[[2]], res, diff(agerange)/res))

model_name <- paste0('fit_brm_', left_side)

## Fit Bayesian model
## for nonlinearity test, set data=hcpd_data
fit_brm <- brms::brm(brms::bf(form), 
                   data=d,
                   chains = 4, cores = 4,
                   iter = 4500, warmup = 2000,
                   control = list(adapt_delta = .999, max_treedepth = 20),
                   file = model_name, silent = FALSE)