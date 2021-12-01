rm(list = ls())

# R Setup ----
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(curvish)
library(parallel)

# Set working directory where fit_brm output will be stored
setwd('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/full_model')

# Load HCP data ----
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")

## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)

# Define response variables
response_vars <- c(paste0('myelin_glasser_v', 1:360)) 

# Define number of subjects
nsub <- dim(hcpd_data)[1]

model_list <- lapply(1:length(response_vars), function(i){
  
  # right_side <- '~ s(Age) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
  right_side <- paste0("~ s(Age) + thickness_glasser_v", i," + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor")
  
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

names(model_list) <- response_vars

## Setup model covariates
## Create age vector sampling full age range
pred_df <- data.frame(Age = c(min(hcpd_data$Age), max(hcpd_data$Age)))

#Use this instead of the below to get predictions at medians rather than at the intercept proper
pred_df$Sex <- as.factor("F")
pred_df$Scanner <- as.factor("Harvard")
pred_df$Reference.Voltage  <- median(hcpd_data$Reference.Voltage)
pred_df$Mean.Pseudo.Transmit.Map  <- median(hcpd_data$Mean.Pseudo.Transmit.Map)
pred_df$T2..Dropout.Threshold  <- median(hcpd_data$T2..Dropout.Threshold)
pred_df$Smoothing.FWHM.mm  <- median(hcpd_data$Smoothing.FWHM.mm)
pred_df$Pseudotransmit.Reference.Value  <- median(hcpd_data$Pseudotransmit.Reference.Value)
pred_df$Correction.Equation.Slope  <- median(hcpd_data$Correction.Equation.Slope)
pred_df$Corrected.CSF.Regressor  <- median(hcpd_data$Corrected.CSF.Regressor)


afit <- model_list[[5]]
age_minmax <- posterior_epred(afit, newdata = pred_df, summary = FALSE)
total_change <- apply(age_minmax, 1, diff)
deriv <- curvish::derivatives(afit, term = 'Age', n = 100, eps = 1e-4, order = 2)
mean_deriv2_post <- apply(deriv$deriv2_posterior, 1, function(x) mean(abs(x)))

# hist(age_minmax[,1])
# dim(age_minmax[,1])
# hist(age_minmax[,2])
# dim(age_minmax[,2])
# hist(total_change)
# dim(total_change)
# class(deriv)
# plot(deriv, deriv = 2)
# plotdf <- as.data.table(deriv$deriv2_posterior)
# plotdf[, iter := 1:.N]
# plotdf_l <- melt(plotdf, id.vars = 'iter')
# plotdf_l[, variable := as.numeric(gsub('V(\\d+)', '\\1', variable))]
# names(plotdf)
# ggplot(plotdf_l[iter %in% sample(1:(max(iter)), 1000)], aes(x = variable, y = value)) + 
#   geom_line(aes(group = iter), alpha = .05)
# hist(deriv$deriv2_posterior)
# dim(deriv$deriv2_posterior)
# hist(mean_deriv2_post)
# length(mean_deriv2_post)


## Create array of spline properties
spline_property_array_fn <- '/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/spline_property_posterior_GB.rds'

# spline_property_array_fn <- '/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/corticalThickness_cov/spline_property_posterior.rds'
# spline_property_array_fn <- '/ncf/hcp/data/analyses/myelin/spline_property_posterior_GB.rds'


if(!file.exists(spline_property_array_fn)){
  #Get posteriors in parallel. Faster but more code.
  cl <- parallel::makeCluster(24)
  ##!!! Check the number of cores ^^^
  nada <- parallel::clusterEvalQ(cl = cl, {library(brms); library(curvish)})
  parallel::clusterExport(cl = cl, c('pred_df'))
  model_list_split <- split(model_list, sort(rep_len(1:24, length.out = length(model_list))))
  
  #this takes a little while (10 minutes)
  system.time({auc_int_posteriors <- parallel::parLapply(cl = cl, X = model_list_split, fun = function(fits){
    lapply(fits, function(afit){
      age_minmax <- posterior_epred(afit, newdata = pred_df, summary = FALSE)
      total_change <- apply(age_minmax, 1, diff)
      deriv <- curvish::derivatives(afit, term = 'Age', n = 100, eps = 1e-4, order = 2)
      mean_deriv2_post <- apply(deriv$deriv2_posterior, 1, function(x) mean(abs(x)))
      return(c(age_minmax[,1], age_minmax[,2], total_change, mean_deriv2_post))
    })
  })})
  parallel::stopCluster(cl)
  #End parallel chunk
  
  spline_property_array <- array(unlist(auc_int_posteriors), dim = c(10000, 4, length(model_list)), 
                                 dimnames = list(NULL, c('age_min', 'age_max', 'total_change', 'mean_deriv2'), NULL))
  saveRDS(spline_property_array, spline_property_array_fn)
} else {
  spline_property_array <- readRDS(spline_property_array_fn)
}

###########################
dim(spline_property_array)
class(spline_property_array[,1,1])
class(spline_property_array[1,,1])
class(spline_property_array[1,1,])

#get correlation posteriors
cor_meanderiv2_agemin_posterior <- apply(spline_property_array, 1, function(draw){
  cor(draw['mean_deriv2',], draw['age_min', ])
})
cor_meanderiv2_totalchange_posterior <- apply(spline_property_array, 1, function(draw){
  cor(unlist(draw['mean_deriv2',]), unlist(draw['total_change', ]))
})
cor_agemin_totalchange_posterior <- apply(spline_property_array, 1, function(draw){
  cor(draw['age_min',], draw['total_change', ])
})
cor_meanderiv2_agemax_posterior <- apply(spline_property_array, 1, function(draw){
  cor(draw['mean_deriv2',], draw['age_max', ])
})
cor_agemax_totalchange_posterior <- apply(spline_property_array, 1, function(draw){
  cor(unlist(draw['age_max',]), unlist(draw['total_change', ]))
})
cor_agemin_agemax_posterior <- apply(spline_property_array, 1, function(draw){
  cor(draw['age_min',], draw['age_max', ])
})

#Summarize and plot them
cors_posterior_array <- array(c(cor_meanderiv2_agemin_posterior,
                                cor_meanderiv2_totalchange_posterior,
                                cor_agemin_totalchange_posterior,
                                cor_meanderiv2_agemax_posterior,
                                cor_agemax_totalchange_posterior,
                                cor_agemin_agemax_posterior), 
                              dim = c(10000,6), 
                              dimnames = list(NULL, c('cor(mean 2nd deriv, Age Min)',
                                                      'cor(mean 2nd deriv, Total Change)',
                                                      'cor(Age Min, Total Change)',
                                                      'cor(mean 2nd deriv, Age Max)',
                                                      'cor(Age Max, Total Change)',
                                                      'cor(Age Min, Age Max)')))
posterior_summary(cors_posterior_array)

bayesplot::mcmc_areas_ridges(cors_posterior_array, prob = .95, prob_outer = .99)
bayesplot::mcmc_pairs(cors_posterior_array, off_diag_args = list(size = .2, alpha = .5))


#correlations for one region
bayesplot::mcmc_pairs(spline_property_array[,,1], off_diag_args = list(size = .2, alpha = .5))


#get summaries to plot scatter
posterior_auc_int_sum <- apply(spline_property_array, c(2,3), function(roi_post){
  data.table(posterior_summary(roi_post))
})
posterior_auc_int_sum <- apply(posterior_auc_int_sum, 2, function(post_sum){
  rbindlist(post_sum, idcol = 'stat')
})
names(posterior_auc_int_sum) <- response_vars
posterior_auc_int_sum <- rbindlist(posterior_auc_int_sum, idcol = 'roi')
posterior_auc_int_sum_l <- melt(posterior_auc_int_sum, id.vars = c('roi', 'stat'))
posterior_auc_int_sum_l[, c('variable', 'stat') := list(paste(stat, variable, sep = '_'), NULL)]
posterior_auc_int_sum_w <- dcast(posterior_auc_int_sum_l, roi ~ variable)

ggplot(posterior_auc_int_sum_w, aes(x = age_min_Estimate, y = total_change_Estimate)) + 
  geom_errorbar(aes(xmin = age_min_Q2.5, xmax = age_min_Q97.5), alpha = .1) + 
  geom_errorbar(aes(ymin = total_change_Q2.5, ymax = total_change_Q97.5), alpha = .1) + 
  geom_point(alpha = 1) + 
  theme_minimal() + 
  coord_cartesian(x = c(-2.5, 2.5), y = c(-.1, .35))

################################
## Regional spline properties ##
################################
# Myelin at Age 8
age8_array <- spline_property_array[,1,]
age8_array <- t(age8_array)
mean_posterior_age8_myelin <- rowMeans(age8_array, na.rm = FALSE, dims = 1)
write.table(mean_posterior_age8_myelin, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_regional_mean_posterior_MyelinAge8.txt", row.names=F, col.names=F, quote=F)

# Total Change (AUC reflects magnitude of increase in T1w/T2w myelin)
total_change_array <- spline_property_array[,3,]
total_change_array <- t(total_change_array)
mean_posterior_total_change <- rowMeans(total_change_array, na.rm = FALSE, dims = 1)
age_range <- max(hcpd_data$Age) - min(hcpd_data$Age)
annualized_roc <- mean_posterior_total_change / age_range
quantile(annualized_roc)

write.table(mean_posterior_total_change, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_total_change_Glasser_regional_myelin.txt", row.names=F, col.names=F, quote=F)
write.table(annualized_roc, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_annualized_roc_Glasser_regional_myelin.txt", row.names=F, col.names=F, quote=F)

# Mean Second Derivative (Reflects degree of nonlinearity)
mean_deriv2_array <- spline_property_array[,4,]
mean_deriv2_array <- t(mean_deriv2_array)
age_range <- max(hcpd_data$Age) - min(hcpd_data$Age)
mean_posterior_deriv2 <- rowMeans(mean_deriv2_array, na.rm = FALSE, dims = 1)
quantile(mean_posterior_deriv2)

write.table(mean_posterior_deriv2, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_deriv2_Glasser_regional_myelin.txt", row.names=F, col.names=F, quote=F)


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

write.table(median_age_of_max_posterior_deriv, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_median_posterior_age_of_max_slope_myelination.txt", row.names=FALSE, col.names=FALSE)



#################################################
# Regional Correlations of Spline Properties ----
#################################################
library(Hmisc)

# Regional Slope ~ Intercept Correlation
rcorr(annualized_roc, mean_posterior_age8_myelin )

# Regional Slope ~ Curvature Correlation
rcorr(annualized_roc, mean_posterior_deriv2)

# Regional Intercept ~ Curvature Correlation
rcorr(mean_posterior_age8_myelin, mean_posterior_deriv2)

# new_age8_myelin <- brm_predictions_df[,1]
# write.table(new_age8_myelin, "/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/n628_regional_mean_posterior_brm_postepred_MyelinAge8.txt", row.names=F, col.names=F, quote=F)

# # Intercept
# intercept_array <- auc_int_roi_array[,3,]
# intercept_array <- t(intercept_array)
# mean_posterior_intercept <- rowMeans(intercept_array, na.rm = FALSE, dims = 1)
# 
# write.table(mean_posterior_intercept, "/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/n628_mean_posterior_intercept_Glasser_regional_myelin.txt", row.names=F, col.names=F, quote=F)
save.image("/ncf/hcp/data/analyses/myelin/spline_property_posterior_GB.Rdata")
