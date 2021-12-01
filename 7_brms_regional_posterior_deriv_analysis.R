rm(list = ls())

# R Setup ----
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(curvish)

## Load R workspace from prior analysis
load('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_Glasser_regional_myelin_brms_posterior_deriv_analysis.Rdata')

# Set working directory where fit_brm output will be stored
setwd('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/full_model')

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

#########################
## R graphics settings ##
#########################
# Color palettes for plotting
jpal <- paste0('#', c('3A4965', 'EAE3D8', 'FAF8F2', 'C07D59', 'F99FC9'))
apal <- paste0('#', c('000000', 'EAE3D8', 'FFFFFF', 'FFD378', '424C6D'))

jftheme <- theme_minimal() +  
  theme(text = element_text(family = 'Didact Gothic', size = 14),
        panel.background = element_rect(fill = jpal[[3]], size = 0, color = jpal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = jpal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text =  element_text(color = jpal[[1]]), axis.title = element_text(color = jpal[[1]]))

gbtheme <- theme_classic() +  
  theme(text = element_text(size = 36),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text =  element_text(color = apal[[1]]), axis.title = element_text(color = apal[[1]]),
        axis.ticks.length=unit(0.5, "cm")) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks=element_line(colour = 'black', size = 1), axis.ticks.length = unit(.5, "cm"))
############################################


# Load HCP data ----
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")


## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)

## Create binary flag indicating if subject hit ceiling for vnav acquisition
## (this means they had high motion)
motion_idx <- which(hcpd_data$numNavs_T1w==196 | hcpd_data$numNavs_T2w==122)
hcpd_data$vnav_ceiling <- 0
hcpd_data$vnav_ceiling[motion_idx] <-  1


# Define response variables
nreg <- 360
response_vars <- c(paste0('myelin_glasser_v', 1:nreg))

# Define number of subjects
nsub <- dim(hcpd_data)[1]

model_list <- lapply(1:length(response_vars), function(i){
  
  right_side <- '~ s(Age) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'

  left_side <- response_vars[[i]]
  
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


# bayes_r2 <- lapply(model_list, brms::bayes_R2)
posterior_samples_list <- lapply(model_list, get_smooth_posterior, res, eps)
results_summaries_list <- lapply(posterior_samples_list, get_smooth_summaries)


# Extract output from each regional model ----
nreg <- 360   ##360 for glasser, 400 for schaefer

# Age range for computing Annualized Rate of Change
age_range <- max(hcpd_data$Age) - min(hcpd_data$Age)

# Define function to calculate mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Define output variables of interest
age_of_max_median_deriv <- rep(NA, nreg)
age_of_plateau <- rep(NA, nreg)
age_of_max_posterior_deriv <- rep(NA, nreg)
median_age_of_max_posterior_deriv <- rep(NA, nreg)
regional_rsq <- rep(NA, nreg)
deriv_auc <- rep(NA, nreg)
mean_post_roc <- rep(NA, nreg)

lower_auc <- rep(NA, nreg)
upper_auc <- rep(NA, nreg)
smooth_est <- as.data.frame(array(rep(NA, nreg*res), dim=c(nreg, res)))


for (i in 1:length(response_vars)){
  
  print(response_vars[[i]])
  
  aderivative <- results_summaries_list[[i]]$derivative
  anageatmax <- results_summaries_list[[i]]$max_deriv_age_posterior
  asmooth <- results_summaries_list[[i]]$smooth
  age_at_max_med_deriv <- aderivative[deriv1 == max(deriv1), Age]
  
  #In order to put the smooth on the same scale as the data, we need to add back
  #in the mean of the response variable that brms uses to compute the
  #conditional effect of the smooth. The easiest way to do this is to use the
  #output of the conditional effects output. We don't really care about the effect across
  #age so we can get the output at an arbitrary age:
  
  ce <- conditional_effects(model_list[[i]], effects = 'Age', 
                            int_conditions = list('Age' = c(8)))
  #get the model frame too, so we can plot the points to confirm this works.
  mf <- brms:::model.frame.brmsfit(model_list[[i]])
  
  response_mf_mean <- ce$Age[[response_vars[[i]]]]
  
  
  age_of_max_median_deriv[i] <- age_at_max_med_deriv
  
  
  age_of_plateau[i] <- min(aderivative$Age[which(aderivative$compared_to_max=="less_steep")])

  
  #We will use the posterior of the derivative to compute the AUC
  fit_brm <- model_list[[i]]
  fit_deriv <- curvish::derivatives(object = fit_brm, term = 'Age', 
                                    n = 50, eps = 1e-07, ##CHANGE n IF YOU WANT
                                    prob = .95, prob_outer = .99, 
                                    deriv_posterior = TRUE)
  #Get the AUC posterior
  fit_deriv_auc <- curvish::auc(fit_deriv, multimodal = FALSE, prob = .99)
  
  #Get the point estimate
  deriv_auc[i] <- mean(fit_deriv_auc$param_posterior)
  
  #Get annualized rate of change
  posterior_annualized_roc <- fit_deriv_auc$param_posterior/age_range
  mean_post_roc[i] <- mean(posterior_annualized_roc)
  
  #Get upper and lower CI estimates
  lower_auc[i] <- quantile(fit_deriv_auc$param_posterior, 0.025)
  upper_auc[i] <- quantile(fit_deriv_auc$param_posterior, 0.975)
  
  #  lower_auc[i] <- fit_deriv_auc$param_posterior_sum[1]
  #  upper_auc[i] <- fit_deriv_auc$param_posterior_sum[2]
  
  #Get mode of posterior age of max derivative
  age_of_max_posterior_deriv[i] <- getmode(anageatmax$max_age)  
  
  #Get median of posterior age of max derivative
  median_age_of_max_posterior_deriv[i] <- median(anageatmax$max_age)  
  
  #Get smooth estimate across age range
  smooth_est[i,] <- asmooth$est + response_mf_mean
  
}

hist(deriv_auc, col="slategray3")
hist(age_of_max_median_deriv, col="slategray3")
hist(age_of_plateau, col="slategray3")  # note that many regions that "Inf" values, not shown here
hist(age_of_max_posterior_deriv, col="slategray3")  # note that many regions that "Inf" values, not shown here
hist(median_age_of_max_posterior_deriv, col="slategray3")

## Save R Workspace
# save.image("/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model/n628_Glasser_regional_myelin_brms_output.Rdata")

## Export output
write.table(mean_post_roc, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/corticalThickness_cov/posterior_deriv_analysis/n628_Glasser_myelin_thicnessCov_mean_posterior_annualized_roc.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(deriv_auc, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_mean_deriv_AUC.txt", row.names=FALSE, col.names=FALSE)
write.table(age_of_max_median_deriv, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_age_of_max_median_posterior_deriv_regional_myelination_GlasserMMP.txt", row.names=FALSE, col.names=FALSE)
write.table(age_of_max_posterior_deriv, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_mode_posterior_age_of_max_deriv_myelination.txt", row.names=FALSE, col.names=FALSE)
write.table(median_age_of_max_posterior_deriv, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_median_posterior_age_of_max_deriv_myelination.txt", row.names=FALSE, col.names=FALSE)
################################################

###############################
## Bayesian Plateau Analysis ##
###############################

## Note, "plateau" is defined loosely here - when a segment of the posterior derivative is credibly less steep than the max slope
print("How many brain regions did not show plateaus?")
length(which(age_of_plateau=="Inf"))
(length(which(age_of_plateau=="Inf")) / nreg) * 100
print("How many brain regions did show plateaus?")
360 - length(which(age_of_plateau=="Inf"))


## Create a binary index of regions that show plateau (1) vs. regions that do not show plateau (0)
reg_idx <- 1:360

reg_idx[which(age_of_plateau=="Inf")] <- 0
reg_idx[which(age_of_plateau!="Inf")] <- 1

write.table(reg_idx, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_regional_myelin_Bayesian_plateau_index.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


age_of_plateau[which(age_of_plateau=="Inf")] <- 0
write.table(age_of_plateau, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_age_of_myelination_plateau_thresh.txt", row.names=FALSE, col.names=FALSE)



## Threshold Age of max median deriv
age_of_max_median_deriv_thresh <- age_of_max_median_deriv
age_of_max_median_deriv_thresh[which(reg_idx==0)] <- 0
write.table(age_of_max_median_deriv_thresh, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_age_of_max_median_posterior_deriv_regional_myelination_GlasserMMP_BPindex_thresh.txt", row.names=FALSE, col.names=FALSE)



#########################################################
## Curvish: Get Median Posterior Age at Max Derivative ##
#########################################################
fit_post <- lapply(model_list, function(fit_brm){
  fit_post <- curvish::posterior_x_at_maxy(object = fit_brm, term = 'Age', multimodal = FALSE,
                                           n = 100, eps = 1e-07, ##CHANGE n IF YOU WANT
                                           prob = .95)
  
  #Get the median Age at max derivative (steepest slope)
  median_age  <- median(fit_post$param_posterior)
  #Get the point estimate
  # return(unlist(median_age))
  return(fit_post)
  
})

## Plot all posterior distributions
fit_post_plots <- lapply(1:length(fit_post), function(i) { 
  p <- plot(fit_post)
  ggsave(sprintf('/ncf/hcp/data/analyses/myelin/figures/posterior_plots/fit_post_%03d.png', i), )
  return(p)
}) 


## Plot Left V1 posterior age of steepest slope
tmp_fit <- fit_post[[181]]
plot(tmp_fit, range = c(9, 20), adjust = 2.0, robust = TRUE) + gbtheme
median(tmp_fit$param_posterior)
ggsave('/ncf/hcp/data/analyses/myelin/figures/posterior_plots/fit_post_v181_unimodal_wide.png', width=12, height=8, unit="in", dpi=500)
 
## Plot Right PFC 8c posterior age of steepest slope
tmp_fit <- fit_post[[73]]
plot(tmp_fit, range = c(9, 20), adjust = 2.0, robust = TRUE) + gbtheme
median(tmp_fit$param_posterior)
ggsave('/ncf/hcp/data/analyses/myelin/figures/posterior_plots/fit_post_v73_uniimodal_wide.png', width=12, height=8, unit="in", dpi=500)



#####################################################
## Curvish: Median Posterior Age at Max Derivative ##
#####################################################
median_age_at_max_deriv <- lapply(model_list, function(fit_brm){
  fit_post <- curvish::posterior_x_at_maxy(object = fit_brm, term = 'Age', multimodal = TRUE,
                                           n = 100, eps = 1e-07, ##CHANGE n IF YOU WANT
                                           prob = .95)
  
  #Get the median Age at max derivative (steepest slope)
  median_age  <- median(fit_post$param_posterior)
  #Get the point estimate
  return(unlist(median_age))
})

write.table(median_age_of_max_posterior_deriv, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_median_posterior_age_of_max_deriv_myelination.txt", row.names=FALSE, col.names=FALSE)

###################################################################
# Apply function to extract first and second order derivatives ----
###################################################################
res <- 100  # Resolution to sample age range

mean_deriv2_posteriors <- lapply(model_list, function(fit_brm){
  deriv_b <- curvish::derivatives(fit_brm, term = 'Age', n = res, eps = 1e-7, order = 2)
  mean_deriv2_post <- apply(deriv_b$deriv2_posterior, 2, function(x) mean(abs(x)))
  return(mean_deriv2_post)
})

mean_regional_deriv2 <- lapply(mean_deriv2_posteriors, mean)
mean_regional_deriv2 <- unlist(mean_regional_deriv2)
hist(mean_regional_deriv2, col="slategray3")

write.table(mean_regional_deriv2, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/mean_regional_deriv2_posteriors.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)