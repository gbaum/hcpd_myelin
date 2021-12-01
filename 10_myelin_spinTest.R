rm(list = ls())

library(Hmisc)
library(corrplot)
library(ggplot2)
library(matrixStats)

#########################
## R graphics settings ##
#########################
# Color palettes for plotting
apal <- paste0('#', c('000000', 'EAE3D8', 'FFFFFF', 'FFD378', '424C6D'))

gbtheme <- theme_classic() +  
  theme(text = element_text(size = 36),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text = element_text(color = apal[[1]]), axis.title = element_text(color = apal[[1]]),
        axis.ticks.length=unit(0.5, "cm")) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks=element_line(colour = 'black', size = 1), axis.ticks.length = unit(.5, "cm"))

#################################
## Compare Regional Properties ##
#################################
SA_axis <- read.table("/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360.txt", header=FALSE)
names(SA_axis) <- "SA_axis"

scaled_SA_axis <- read.table("/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360_scaled_centered.txt", header=FALSE)
names(scaled_SA_axis) <- "scaled_SA_axis"

bayes_R2 <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/bayes_r2/n628_GlasserMMP_myelin_sAge_partial_bayes_r2.txt", header=FALSE)
names(bayes_R2) <- "bayes_R2"

rate_of_change <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_annualized_roc_Glasser_regional_myelin.txt", header=FALSE)
names(rate_of_change) <- "rate_of_change"

# age_of_peak_growth <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_age_of_max_median_posterior_deriv_regional_myelination_GlasserMMP.txt", header=FALSE)
# names(age_of_peak_growth) <- "age_of_peak_growth"

median_age_of_peak_growth <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_median_posterior_age_of_max_slope_myelination.txt", header=FALSE)
names(median_age_of_peak_growth) <- "median_age_of_peak_growth"

nonlinearity <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_deriv2_Glasser_regional_myelin.txt", header=FALSE)
names(nonlinearity) <- "nonlinearity"

func_clusters <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_myelin_brm_Age_trajectory_hddcClust3_membership_relevel.txt", header=FALSE)
names(func_clusters) <- "func_clusters"
func_clusters <- as.numeric(func_clusters$func_clusters)

myelin_age8 <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_regional_mean_posterior_MyelinAge8.txt", header=FALSE)
names(myelin_age8) <- "myelin_age8"

## Define S-A Rank
SA_rank <-rank(SA_axis$SA_axis)
names(SA_rank) <- "SA_rank"


#######################################
## Create dataframe for correlations ##
#######################################
df <- cbind(SA_rank, bayes_R2, rate_of_change, median_age_of_peak_growth, nonlinearity, func_clusters)
names(df)

## Redefine missing label for latent clustering
df$func_clusters[344] <- 3
# func_clusters$func_clusters <- as.ordered(func_clusters$func_clusters)



## Flip right and left hemisphere values for spin test
df_spin <- df
left_idx <- 1:180
right_idx <- 181:360

df_spin[left_idx,] <- df[right_idx,]
df_spin[right_idx,] <- df[left_idx,]

## Spherical coordinates for Glasser-HCP Parcellation
spherical.coordinates <- read.table("/ncf/hcp/data/analyses/myelin/code/rotate_parcellation-master/sphere_HCP.txt", header=F)  #glasser parcel spherical coordinates (FreeSurfer sphere) for the 360 parcels (left hemisphere= rows 1-180, right hemisphere = rows 181-360). sphere_HCP.txt file provided here https://github.com/frantisekvasa/rotate_parcellation
spherical.coordinates <- as.matrix(spherical.coordinates) #format as matrix

# Parcel-based spatial permutation test (spin test) functions
### Functions were originally obtained from https://github.com/frantisekvasa/rotate_parcellation/tree/master/R 
### Parcel-based tests implemented as described in Váša et al. (2018) doi:10.1093/cercor/bhx249 , based on original work from Alexander-Bloch et al. (2013) doi:10.1523/JNEUROSCI.3554-12.2013 and Alexander-Bloch et al. (2018) doi:10.1016/j.neuroimage.2018.05.070
source("/ncf/hcp/data/analyses/myelin/code/rotate_parcellation-master/R/rotate.parcellation.R") 
source("/ncf/hcp/data/analyses/myelin/code/rotate_parcellation-master/R/perm.sphere.p.R")

#########################################
## Test for single spatial correlation ##
#########################################
coord.l <- as.matrix(spherical.coordinates[left_idx, ])
coord.r <- as.matrix(spherical.coordinates[right_idx, ])

# Time the spin test
ptm <- proc.time()

set.seed(232)
perm.id <- rotate.parcellation(coord.l, coord.r, nrot=1000) #nrot Random rotations

correlation = cor(df_spin$rate_of_change, df_spin$SA_rank, method = c("spearman")) #spearman's rank correlation matrix
permutedp = perm.sphere.p(df_spin$rate_of_change, df_spin$SA_rank, perm.id) #spatial permutation based p-value matrix

# Stop the clock
proc.time() - ptm

#################################################################
## Compute Spearman's Rank Correlation Matrix and Pspin Matrix ##
#################################################################
ndim <- dim(df_spin)[2]
correlation.matrix = permutedp.matrix = array(NA, dim=c(ndim,ndim)) #create empty matrices 

coord.l <- as.matrix(spherical.coordinates[left_idx, ])
coord.r <- as.matrix(spherical.coordinates[right_idx, ])

set.seed(232)

for(i in 1:ndim){
  for(j in 1:ndim){
    perm.id <- rotate.parcellation(coord.l, coord.r, nrot=1000) #1000 spherical rotations
    correlation.matrix[i,j] = cor(df_spin[,i], df_spin[,j], method = c("spearman")) #spearman's rank correlation matrix
permutedp.matrix[i,j] = perm.sphere.p(df_spin[,i], df_spin[,j], perm.id) #spatial permutation based p-value matrix
}
}

# Stop the clock
proc.time() - ptm

# Save output
save(correlation.matrix, file="/ncf/hcp/data/analyses/myelin/code/myelin_spin_test/myelinSpinTest.spearmans.correlation.matrix_v2.Rdata") #save correlation matrix
save(permutedp.matrix, file="/ncf/hcp/data/analyses/myelin/code/myelin_spin_test/myelinSpinTest.p.matrix_v2.Rdata") #save permuted p-value matrix
