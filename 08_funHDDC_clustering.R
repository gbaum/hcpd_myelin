rm(list = ls())

######################
## Load R Workspace ##
######################
load("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_regional_myelin_brms_hddc_clustering_output.Rdata")

# R Setup ----
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(curvish)
library(funHDDC)
library(fda.usc)
library(colorspace)

# Set working directory where fit_brm output will be stored
setwd("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/full_model")

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
  # right_side <- paste0("~ s(Age) + thickness_glasser_v", i," + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor")
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

names(model_list) <- response_vars

#############################################################
## Make data-frame with model covariates -- For prediction ##
#############################################################
## Define Age resolution (sampling across age range)
res=500

## Create age vector sampling full age range
pred_df <- NULL
set.seed(232)
pred_df <- seq(min(hcpd_data$Age), max(hcpd_data$Age), length=res)     
hist(pred_df, col="slategray3")
pred_df <- as.data.frame(pred_df)

##  Add model covariates
pred_df$Sex <- as.factor("F")
pred_df$Scanner <- as.factor("Harvard")
pred_df$Reference.Voltage  <- median(hcpd_data$Reference.Voltage)
pred_df$Mean.Pseudo.Transmit.Map  <- median(hcpd_data$Mean.Pseudo.Transmit.Map)
pred_df$T2..Dropout.Threshold  <- median(hcpd_data$T2..Dropout.Threshold)
pred_df$Smoothing.FWHM.mm  <- median(hcpd_data$Smoothing.FWHM.mm)
pred_df$Pseudotransmit.Reference.Value  <- median(hcpd_data$Pseudotransmit.Reference.Value)
pred_df$Correction.Equation.Slope  <- median(hcpd_data$Correction.Equation.Slope)
pred_df$Corrected.CSF.Regressor  <- median(hcpd_data$Corrected.CSF.Regressor)
colnames(pred_df)[1] <- "Age"

## test
afit <- model_list[[1]]
test <- posterior_summary(posterior_epred(afit, newdata = pred_df))

#################################################################
## Extract predicted myelin values using brms::posterior_epred ##         
#################################################################
brm_pred <- lapply(model_list, posterior_epred, newdata=pred_df)
brm_sum <- lapply(brm_pred, posterior_summary, newdata=pred_df)
brm_pred_estimates <- lapply(brm_sum, function(x) x[, 1])
brm_predictions_df <- do.call(rbind, brm_pred_estimates)

###########################################################
## Extract predicted myelin values using predict.brmsfit ##         
###########################################################
# brm_pred <- lapply(model_list, predict, newdata=pred_df)
# all(unlist(lapply(brm_pred, function(x) dim(x)[[1]]==res))) ## check dimensions of output
# brm_pred_estimates <- lapply(brm_pred, function(x) x[, 1])
# brm_predictions_df <- do.call(rbind, brm_pred_estimates)

#########################################################
## BRMS Predictions with Subtracted Regional Intercept ##
#########################################################
nreg <- 360
# mean_posterior_intercept <- read.table("/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/n628_mean_posterior_intercept_Glasser_regional_myelin.txt")

## Clustering on these predicted values should tell us more about differences in the *shape* of age curves, 
## rather than just differences in intercept
brm_predictions_minus_intercept <- as.data.frame(matrix(NA, nrow=360, ncol=500))

for (i in 1:nreg){
  tmp_row <- brm_predictions_df[i,]
  tmp_intercept <- tmp_row[1] # take myelin estimate at min_age
  new_row <- tmp_row - tmp_intercept 
  brm_predictions_minus_intercept[i,] <- new_row
}

################################################
## Create Region by Age Matrix for Clustering ##
################################################
region_by_age <- as.matrix(brm_predictions_minus_intercept)  # Matrix where rows are brain regions, and columns contain predicted brms model values for each age point sampled in `res`
write.table(region_by_age, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_regional_myelin_brm_epred_minusMyelinAtMinAge_across_age_res500.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

region_by_age_brm_pred <- as.matrix(brm_predictions_df)
# write.table(region_by_age_brm_pred, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_regional_myelin_brm_epred_across_age_res500.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#########################################################
## EXCLUDE REGION 344 DUE TO NEGATIVE SLOPE (artifact) ##
#########################################################
region_by_age <- region_by_age[-344,]
dim(region_by_age)

region_by_age_brm_pred  <- region_by_age_brm_pred [-344,]
dim(region_by_age_brm_pred)

# Cluster Regional Myelin Trajectories with FDA ----

#Age vector
# Age <- seq(min(hcpd_data$Age), max(hcpd_data$Age), length=res)    
Age <- pred_df$Age

#Create functional data object
fbrain <- fdata(region_by_age, argvals = Age, 
                 names = list(main = "Regional Myelin Trajectories", xlab = "Age", ylab = "T1w/T2w Myelin"))

fbrain2_fd <- fdata2fd(fbrain)      ## convert into fd format


#Create functional data object
fbrain_orig <- fdata(region_by_age_brm_pred, argvals = Age, 
                names = list(main = "Regional Myelin Trajectories", xlab = "Age", ylab = "T1w/T2w Myelin"))

fbrain_orig_fd <- fdata2fd(fbrain_orig )      
## ---------------------- clustering ----------------------
## reduce the resolution for running time purposes
# nreg <- 360
# ind <- round(seq(from = 1, to = 1000, length.out = 500), 0)
# ind
# region_by_age2 <- region_by_age[,ind]
# fbrain2 <- fdata(region_by_age2, argvals = pred_df$Age[ind], 
#                  names = list(main = "Glasser-MMP Myelin Trajectories", xlab = "Age", ylab = "T1w/T2w Myelin"))
# plot(fbrain2, main = "Glasser-MMP Myelin Trajectories", lty = 1, lwd=1.2, col = rainbow_hcl(nreg))
# 
# fbrain2_fd <- fdata2fd(fbrain2)      ## convert into fd format

#Uses BIC as selection criterion
set.seed(232)
fit2_hddc <- funHDDC(fbrain2_fd, K = 2:7)
slopeHeuristic(fit2_hddc)

############################
## Clustering Diagnostics ##
############################
source("/ncf/hcp/data/analyses/myelin/code/gb_slopeHeuristic.R")

gb_slopeHeuristic(fit2_hddc)

## Export Diagnostic Plot for Figure 7
ggsave("/ncf/hcp/data/analyses/myelin/figures/funHDDC_slopeHeuristic_logLikelihood_plot.png", width = 12, height = 8, unit="in", dpi=500)

#######################################
## Fit best-fitting k-means solution ##
#######################################
set.seed(232)
fit_hddc3 <- funHDDC(fbrain2_fd, K = 3)   ## re-fit best model
fit_hddc3$BIC
fit_hddc3$class                ## cluster memberships
table(fit_hddc3$class)  

######################
## Relevel clusters ##
######################
cluster_rlvl <- as.factor(fit_hddc3$class)
one_idx <- which(fit_hddc3$class==1) 
two_idx <- which(fit_hddc3$class==2)
three_idx <- which(fit_hddc3$classs==3)

cluster_rlvl[one_idx] <- 2
cluster_rlvl[two_idx] <- 1
cluster_rlvl[three_idx] <- 3

#######################################################
## Plot Regional Trajectories and Clustering Centers ##
#######################################################

## Pallette: Spectral 
color_pallete_function <- colorRampPalette(
  colors = c("#9B0B5A", "#EE8050", "#ECE27A")  # #A7295F"),
  # space = "Lab" # Option used when colors do not represent a quantitative scale
)

cluster_rlvl <- as.factor(cluster_rlvl)
num_colors <- length(unique(cluster_rlvl))
plot_colors <- color_pallete_function(num_colors)

op <- par(mfrow = c(1,2))
plot(fbrain, col = plot_colors[cluster_rlvl])    
matplot(t(fit_hddc3$mu), type = "l", col = plot_colors, lty = 2, lwd = 3)  ## plotting centroid lines
par(op)

## Append 0 for excluded noisy region (for visualization)
tmp_class <- append(fit_hddc3$class, 0, 343)
write.table(tmp_class, "/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_myelin_brm_Age_trajectory_hddcClust3_membership.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Pallette:  power_surf
color_pallete_function <- colorRampPalette(
  colors = c("#5987C8", "#3C3C3C", "red"),
  #space = "Lab" # Option used when colors do not represent a quantitative scale
)

fit_hddc3$class <- as.factor(fit_hddc3$class)
num_colors <- length(unique(fit_hddc3$class))
plot_colors <- color_pallete_function(num_colors)

op <- par(mfrow = c(1,2))
plot(fbrain, col = plot_colors[fit_hddc3$class], lwd = 1.5)    
matplot(t(fit_hddc3$mu), type = "l", col = plot_colors, lty = 2, lwd = 3)  ## plotting centroid lines
par(op)