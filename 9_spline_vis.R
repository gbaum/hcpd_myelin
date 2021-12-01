rm(list = ls())

# R Setup ----
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(curvish)
library(funHDDC)
library(fda.usc)

## Load workspace with posterior summaries for each regional model
load('/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_Glasser_regional_myelin_brms_posterior_deriv_analysis.Rdata')

# load("/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model/n628_Glasser_regional_myelin_brms_output.Rdata")

## R graphics settings
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
## Create Mean Myelin Map Across Subjects ##
############################################
# Mean Regional T1w/T2w across all subjects
myelin_col_index <- grep("myelin_glasser_v", colnames(hcpd_data))
mean_myelin <- as.data.frame(lapply(hcpd_data[myelin_col_index], mean))
mean_myelin <- t(mean_myelin)
colnames(mean_myelin) <- "mean_myelin"

write.table(mean_myelin, "/ncf/hcp/data/analyses/myelin/data/n628_Aug2021_newCorr_mean_regional_T1wT2w_myelin_GlasserMMP.txt", row.names=FALSE, col.names=FALSE)

#########################################
## DEFINE REGIONAL SMOOTHS FOR FIGURES ##
#########################################


## Figure 1 - Methods Schematic
i <- 4 ## Right Visual cortex (V2) 
i <- 181 ## Left Visual cortex (V1) 


########################
## Figure 3- Bayes R2 ##
########################
i <- 244  # Left ACC, p32
i <- 4   # Right Area V1 

##########################################
## Figure 4 - Annualized Rate of Change ##
##########################################
i <- 249 ## Left Medial PFC, 9m  (low)
i <- 8 ## Right Motor cortex, area 4 (high)

###########################################
## Figure 5 - Age of Peak T1w/T2w Growth ##
###########################################
i<- 181 # Left V1 
i <- 73 # Right PFC 8c

#####################################################
## Figure 6 - Nonlinearity (mean second derivative ##
#####################################################
i<- 263 # Left dlPFC 9-46v
i <- i # Right V1 


## Graphical settings:
# ymin <- 1.1
# ymax <- 2.3

# ymin <- min(mf_points$data[,1])
# ymax <- mean(mf_points$data[,1]) + (3*sd(mf_points$data[,1]))

# scale_y_continuous(breaks=c(1.4, 1.8, 2.2), limits=c(ymin, ymax))
# scale_y_continuous(breaks=c(1.4, 1.6, 1.8), limits=c(ymin, ymax))
## Graphical settings:


#################
## PLOT SMOOTH ##
#################
print(response_vars[[i]])

include_points <- TRUE #overlay points on smooth plots

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

#We want to reference the mean of the response variable, which I've set up to
#be in the list `seven_nets`. I'm naming the variable `response_mf_mean` to
#indicate that this is the mean of the response as determined by brms from the
#model frame (`mf`) that it uses to estimate the model.
response_mf_mean <- ce$Age[[response_vars[[i]]]]


cat(paste0('\n\n### Age at maximum derivative\n\n'))

bayesplot::color_scheme_set('red')
print(bayesplot::mcmc_areas(data.frame(max_age = anageatmax$max_age, Chain = rep(1:4, each = 2500)),
                            prob = 0.95,
                            prob_outer = 0.99) + gbtheme)
cat(paste0('\n\n### Difference from derivative at steepest age\n\n'))

print(ggplot(aderivative, aes(x = Age, y = diff_from_max)) + 
        geom_ribbon(aes(ymin = diff_lower, ymax = diff_upper), alpha = 1, fill = apal[2]) +  
        geom_line(color = apal[[1]]) + 
        geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
        scale_color_manual(breaks = c('less_steep', 'as_steep'), values = apal[c(4,5)], 
                           labels = c('Less steep', 'As steep'),
                           name = 'Region, as compared \nto maximum derivative, \nis...') + 
        gbtheme)
cat(paste0('\n\n### Smooth estimate\n\n'))

## EXPORT PLOT (diff from max)
ggsave(paste0('/ncf/hcp/data/analyses/myelin/figures/smooths/v', i, '_brm_diffFromMax.png'), width = 12, height = 8, unit="in", dpi=500)


mf_points <- NULL
if(include_points){
  mf_points <- geom_point(data = mf, aes_string(y = response_vars[[i]]), alpha = .15)
}


## GB:  M
#This is where we finally plot the smooth! So it's the correct place to add
#in the mean of the response. We have to add it to both the line and the ribbon.

## Set Bounds of Y-axis across Networks
ymin <- min(mf_points$data[,1])
ymax <- mean(mf_points$data[,1]) + (3*sd(mf_points$data[,1]))


print(ggplot(asmooth[aderivative[, .(Age, compared_to_max)], on = 'Age'],
             aes(x = Age, y = est + response_mf_mean)) + 
        mf_points + 
        geom_ribbon(aes(ymin = lower + response_mf_mean, ymax = upper + response_mf_mean), alpha = 1, fill = apal[2]) +  
        geom_line() + labs(y = 'est') +
        
        geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
        scale_color_manual(breaks = c('less_steep', 'as_steep'), values = apal[c(4,5)], 
                           labels = c('Less steep', 'As steep'),
                           name = 'Region, as compared \nto maximum derivative, \nis...')) + 
  gbtheme + # ylim(ymin, ymax)
  scale_y_continuous(breaks=c(1.6, 1.8, 2.0), limits=c(ymin, ymax))
  # scale_y_continuous(breaks=c(1.4, 1.6, 1.8, 2.0), limits=c(ymin, ymax))
  # ylim(ymin, ymax) # + 
  # scale_y_continuous(breaks=c(1.2, 1.4, 1.6, 1.8 , 2.0, 2.2), limits=c(ymin, ymax))

# Set limits on Y-axis
# scale_y_continuous(breaks=c(1.6, 2.0, 2.4), limits = c(1.4, 2.4)) + 

# ylim(ymin, ymax)

# ylim(min(mf_points$data[1]), 2.6)
# scale_y_continuous(limits = c(min(response_vars[[i]], 2.1))) +  
# gbtheme)


## EXPORT PLOT 
ggsave(paste0('/ncf/hcp/data/analyses/myelin/figures/smooths/v', i, '_brm_smooth_narrow.png'), width = 12, height = 9, unit="in", dpi=500)

###################################
## Plot the posterior derivative ##
###################################
print(ggplot(aderivative, aes(x = Age, y = deriv1)) + 
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1, fill = apal[2]) +  
        geom_segment(x = age_at_max_med_deriv, 
                     xend = age_at_max_med_deriv, 
                     y = aderivative[Age == age_at_max_med_deriv, lower],
                     yend = aderivative[Age == age_at_max_med_deriv, upper], color = apal[[1]]) +
        geom_line() + 
        geom_line(aes(color = compared_to_max, group = 1), size = 2, alpha = .8) + 
        geom_point(x = age_at_max_med_deriv, y = aderivative[Age == age_at_max_med_deriv, deriv1], color = apal[[1]], size = 2) +
        scale_color_manual(breaks = c('less_steep', 'as_steep'), values = apal[c(4,5)], 
                           labels = c('Less steep', 'As steep'),
                           name = 'Region, as compared \nto maximum derivative, \nis...') +
        gbtheme) # + scale_y_continuous(breaks=c(-0.02, 0.0, 0.025)))


## EXPORT PLOT 
ggsave(paste0('/ncf/hcp/data/analyses/myelin/figures/smooths/v', i, '_brm_deriv_narrow.png'), width = 12, height = 9, unit="in", dpi=500)
