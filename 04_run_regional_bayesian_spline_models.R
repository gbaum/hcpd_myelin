rm(list = ls())

# R Setup ----
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE) # only necessary if you're going to compile into a notebook
library(brms)
library(data.table)
library(ggplot2)
library(readr)

## Set working directory where fit_brm output will be stored
# setwd("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/linear_model")
setwd('/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model')

# brms setup ----

## Load functions for extracting posterior samples
source('/ncf/hcp/data/analyses/myelin/code/gam_bootstrap/conditional_smooth_sample.R')

#Age difference for finite differences
eps = 1e-07
#How many ages to get marginal posterior at -- doesn't need to be a lot, necessarily
res = 100
######################################################################################

# Load HCP data ----
hcpd_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")


## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)

##########################
## LOW MOTION THRESHOLD ##
##########################
## Create binary flag indicating if subject hit ceiling for vnav acquisition
## (this means they had high motion)
motion_idx <- which(hcpd_data$numNavs_T1w==196 | hcpd_data$numNavs_T2w==122)
hcpd_data$vnav_ceiling <- 0
hcpd_data$vnav_ceiling[motion_idx] <-  1

#################################
## REMOVE HIGH MOTION SUBJECTS ##
#################################
# hcpd_data <- subset(hcpd_data, hcpd_data$vnav_ceiling==0)
# dim(hcpd_data)
############################################################################

## Task index for parallelizing regional models
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

region_id <- task_id

set.seed(123)

# Fit Bayesian smoothing spline models for each brain region ----

right_side <- '~ s(Age) + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'

# right_side <- '~ Age + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
# right_side <- '~ Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
# right_side <- paste0("~ s(Age) + thickness_glasser_v", region_id," + Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor")

left_side <- paste0('myelin_glasser_v', region_id) 

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

########################
## Fit Bayesian model ##
########################
fit_brm <- brms::brm(brms::bf(form), 
                   data=d,  ## for nonlinearity test, set data=hcpd_data
                   chains = 4, cores = 4,
                   iter = 4500, warmup = 2000,
                   control = list(adapt_delta = .999, max_treedepth = 20),
                   file = model_name, silent = FALSE)