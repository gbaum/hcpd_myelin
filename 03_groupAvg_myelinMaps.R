rm(list = ls())

# R Setup ----
library(data.table)
library(ggplot2)
library(readr)

# Load HCP data ----
brain_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv")

####################################
## Age Group-Averaged Myelin Maps ##
####################################
young_df <- subset(brain_data, brain_data$age_floor < 11)
mid_df <- subset(brain_data, brain_data$age_floor > 13 & brain_data$age_floor < 17)
old_df <- subset(brain_data, brain_data$age_floor > 18)

# Mean Regional T1w/T2w across subjects in each Age group
myelin_col_index <- grep("myelin_glasser_v", colnames(brain_data))

## Mean myelin for each group
mean_young_myelin <- as.data.frame(lapply(young_df[myelin_col_index], mean))
mean_young_myelin <- t(mean_young_myelin)
colnames(mean_young_myelin) <- "mean_young_myelin"
mean_young_myelin_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_young_regional_T1wT2w_GlasserMMP_8to10.txt"
mean_young_myelin_cifti_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_young_regional_T1wT2w_GlasserMMP_8to10.pscalar.nii"
write.table(mean_young_myelin, mean_young_myelin_outpath , row.names=FALSE, col.names=FALSE)

mean_mid_myelin <- as.data.frame(lapply(mid_df[myelin_col_index], mean))
mean_mid_myelin <- t(mean_mid_myelin)
colnames(mean_mid_myelin) <- "mean_mid_myelin"
mean_mid_myelin_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_mid_regional_T1wT2w_GlasserMMP_14to16.txt"
mean_mid_myelin_cifti_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_mid_regional_T1wT2w_GlasserMMP_14to16.pscalar.nii"
write.table(mean_mid_myelin, mean_mid_myelin_outpath, row.names=FALSE, col.names=FALSE)

mean_old_myelin <- as.data.frame(lapply(old_df[myelin_col_index], mean))
mean_old_myelin <- t(mean_old_myelin)
colnames(mean_old_myelin) <- "mean_old_myelin"
mean_old_myelin_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_old_regional_T1wT2w_GlasserMMP_19to21.txt"
mean_old_myelin_cifti_outpath <- "/ncf/hcp/data/analyses/myelin/figures/groupAvg_myelin/mean_old_regional_T1wT2w_GlasserMMP_19to21.pscalar.nii"
write.table(mean_old_myelin, mean_old_myelin_outpath, row.names=FALSE, col.names=FALSE)


# Write Cifi Output in Connectome Workbench ----
glasser_template_path <- "/ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_template.pscalar.nii"
schaefer_template_path <- "/ncf/hcp/data/analyses/myelin/parcellations/S1200_Schaefer400x7_template.pscalar.nii"

# cifti.output.path <- "/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model/brm_analysis/n628_mean_posterior_deriv2_Glasser_regional_myelin.pscalar.nii" 
# brainVar.path <- "/ncf/hcp/data/analyses/myelin/brm_output/Glasser-MMP/regional/full_model/brm_analysis/n628_mean_posterior_deriv2_Glasser_regional_myelin.txt"

## Create function that can write out cifti files
write_cifti <- function(template_path, brainVar_path, cifti_output_path) {
  system(paste0("module load connectome-workbench/1.3.2-fasrc01; ", "wb_command -cifti-convert -from-text ", brainVar_path, " ", template_path, " ", cifti_output_path))
  return(print(cifti_output_path))
}

## Export cifti
write_cifti(glasser_template_path, mean_young_myelin_outpath, mean_young_myelin_cifti_outpath)
write_cifti(glasser_template_path, mean_mid_myelin_outpath, mean_mid_myelin_cifti_outpath)
write_cifti(glasser_template_path, mean_old_myelin_outpath, mean_old_myelin_cifti_outpath)
