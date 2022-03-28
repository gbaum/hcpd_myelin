rm(list = ls())

library(Hmisc)
library(corrplot)
library(ggplot2)

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

gbtheme2 <- theme_classic() +  
  theme(text = element_text(size = 42),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text = element_text(color = apal[[1]]), axis.title = element_blank(),
        axis.ticks.length=unit(0.5, "cm")) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks=element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.5, "cm"))

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

age_of_peak_growth <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/posterior_derivative_analysis/n628_myelin_Glasser_regional_age_of_max_median_posterior_deriv_regional_myelination_GlasserMMP.txt", header=FALSE)
names(age_of_peak_growth) <- "age_of_peak_growth"

median_age_of_peak_growth <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_median_posterior_age_of_max_slope_myelination.txt", header=FALSE)
names(median_age_of_peak_growth) <- "median_age_of_peak_growth"

nonlinearity <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/regional/spline_property_posterior/n628_mean_posterior_deriv2_Glasser_regional_myelin.txt", header=FALSE)
names(nonlinearity) <- "nonlinearity"

func_clusters <- read.table("/ncf/hcp/data/analyses/myelin/Aug2021_brm_output/Glasser-MMP/fda_clustering/n628_Glasser_myelin_brm_Age_trajectory_hddcClust3_membership_relevel.txt", header=FALSE)
names(func_clusters) <- "func_clusters"
func_clusters <- as.numeric(func_clusters$func_clusters)


SA_rank <-rank(SA_axis$SA_axis)
names(SA_rank) <- "SA_rank"

# func_clusters$func_clusters <- as.ordered(func_clusters$func_clusters)

#######################################
## Create dataframe for correlations ##
#######################################
df <- cbind(SA_rank, bayes_R2, rate_of_change, median_age_of_peak_growth, nonlinearity, func_clusters)
names(df)

## Redefine missing label for latent clustering
df$func_clusters[344] <- 3

## Define S-A Ranking
df$SA_rank <-as.integer(df$SA_rank)

spearman.df.cor = cor(df, method="spearman")
pearson.df.cor = cor(df, method="pearson")

#############################################
## Load spatial permutation-based p-values ##
#############################################
load('/ncf/hcp/data/analyses/myelin/code/myelin_spin_test/myelinSpinTest.p.matrix.Rdata')   #permuted p-value matrix
load('/ncf/hcp/data/analyses/myelin/code/myelin_spin_test/myelinSpinTest.spearmans.correlation.matrix.Rdata')  #spearman correlation matrix

dim(permutedp.matrix)
dim(correlation.matrix)

#############################################
## Holmes correction of Spin test p-values ##
#############################################
upper_pmat <- permutedp.matrix[upper.tri(permutedp.matrix)]
p_vec <- as.vector(upper_pmat)
p_vec_corr <- p.adjust(p_vec, method="holm")

## Check significance after correction
length(which(p_vec_corr < 0.05)) # how many p-values are significant? (lower than 0.05)
length(which(p_vec_corr > 0.05)) # how many p-values are not significant? (higher than 0.05)

mixed_mat <- correlation.matrix
mixed_mat[upper.tri(correlation.matrix)] <- permutedp.matrix[upper.tri(permutedp.matrix)]
colnames(mixed_mat) <- rownames(mixed_mat) <- colnames(df)

corrplot(mixed_mat,  method="number", col = "black", p.mat = permutedp.matrix)

## Export image
ggsave("/ncf/hcp/data/analyses/myelin/figures/regional_feature_corrplot.png", width = 10, height = 10, unit="in", dpi=500)


##################
## Scatterplots ##
##################

## Note: region index from spline_vis.R

#####################
## Figure 3-D (R2) ##
#####################
region_index <- c(4,244)
highlight_df <- df[region_index, ]
SA_R2_scatter <- ggplot(data=df, aes(x=SA_rank, y=bayes_R2))
SA_R2_scatter + geom_point(alpha = 0.6, size=2, col="#5C5C5C") +  geom_smooth(method="lm", size=2, col="black") + gbtheme +
geom_point(data=highlight_df, aes(x=SA_rank, y=bayes_R2, col=factor(SA_rank)),size=5, alpha=0.9) + 
  scale_color_manual(values = c("#0057b7",  "#ffd700")) +
  theme(legend.position="none")

ggsave('/ncf/hcp/data/analyses/myelin/figures/SA_scatterplots/bayesR2_SA_scatter.png', width = 6, height = 7 , unit="in", dpi=500)

#################################
## Figure 4-D (Rate of change) ##
#################################
region_index <- c(8,249)
highlight_df <- df[region_index, ]
SA_aroc_scatter <- ggplot(data=df, aes(x=SA_rank, y=rate_of_change))
SA_aroc_scatter + geom_point(alpha = 0.6, size=2, col="#5C5C5C") +  geom_smooth(method="lm", size=2, col="black") + gbtheme +
  geom_point(data=highlight_df, aes(x=SA_rank, y=rate_of_change, col=factor(SA_rank)),size=5, alpha=0.9) + 
  scale_color_manual(values = c("#0057b7",  "#ffd700")) +
  theme(legend.position="none")

ggsave('/ncf/hcp/data/analyses/myelin/figures/SA_scatterplots/aroc_SA_scatter.png', width = 6, height = 7, unit="in", dpi=500)

###############################
## Figure 5-D (Nonlinearity) ##
###############################
region_index <- c(111,184)
highlight_df <- df[region_index, ]
SA_nonlinearity_scatter <- ggplot(data=df, aes(x=SA_rank, y=nonlinearity))
SA_nonlinearity_scatter + geom_point(alpha = 0.6, size=2, col="#5C5C5C") +  geom_smooth(method="lm", size=2, col="black") + gbtheme +
  geom_point(data=highlight_df, aes(x=SA_rank, y=nonlinearity, col=factor(SA_rank)),size=5, alpha=0.9) + 
  scale_color_manual(values = c("#0057b7",  "#ffd700")) +
  scale_y_continuous(breaks=c (0.002, 0.004, 0.006, 0.008), limits=c(0.001, 0.010)) +
  theme(legend.position="none")

ggsave('/ncf/hcp/data/analyses/myelin/figures/SA_scatterplots/nonlinearity_SA_scatter.png', width = 6, height = 7, unit="in", dpi=500)


###################################
## Figure 6-D Age of peak growth ##
###################################
region_index <- c(73, 181)
highlight_df <- df[region_index, ]
SA_peak_scatter <- ggplot(data=df, aes(x=SA_rank, y=median_age_of_peak_growth))
SA_peak_scatter + geom_point(alpha = 0.6, size=2, col="#5C5C5C") +  geom_smooth(method="lm", size=2, col="black") + gbtheme + scale_y_continuous(breaks=c(8, 12, 16, 20), limits = c(8,21)) + 
  geom_point(data=highlight_df, aes(x=SA_rank, y=median_age_of_peak_growth, col=factor(SA_rank)),size=5, alpha=0.9) + 
  scale_color_manual(values = c("#0057b7",  "#ffd700")) +
  theme(legend.position="none")

ggsave('/ncf/hcp/data/analyses/myelin/figures/SA_scatterplots/median_agePeakGrowth_SA_scatter.png', width = 6, height = 7, unit="in", dpi=500)



###################################
## S-A rank ~ func cluster anova ##
###################################
df$func_clusters <- as.factor(df$func_clusters)
df$SA_rank <- as.numeric(df$SA_rank)

oneway_anova <- aov(SA_rank ~ func_clusters, data = df)
summary(oneway_anova)
