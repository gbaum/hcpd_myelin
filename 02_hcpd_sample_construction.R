rm(list = ls())

##############################
# Load relevant Libraries ----
##############################
require(ggplot2)
require(Hmisc)
require(MASS)
require(mgcv)
require(stringr)
require(visreg)
require(parallel)
require(multilevel)
require(stats)
require(Formula)
require(reshape2)
require(gratia)

####################################################
## Read in demographics and data quality measures ##
####################################################
harms_QA <- read.csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/freesurfer_qc_20191210.csv", header=TRUE)
hcpd_qa <- subset(harms_QA, harms_QA$Project =="HCD")

demographics_df <- read.csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/GB_FS_QC.csv")

## Merge with demographics
hcpd_data <- merge(hcpd_qa, demographics_df, by="Subject")

## Set variables
hcpd_data$Sex<- as.factor(hcpd_data$Sex)
hcpd_data$Sex<- droplevels(hcpd_data$Sex) # drop empty level of factor 
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)
hcpd_data$Scanner  <- droplevels(hcpd_data$Scanner) # drop empty level of factor 
hcpd_data$SES_RLVL <- as.factor(hcpd_data$SES_RLVL)

# Create age rounded down to integer as factor
hcpd_data$age_int <- as.factor(round(hcpd_data$Age))
hcpd_data$age_floor <- as.factor(floor(hcpd_data$Age))

## Define mean and sum of numNavs measure
hcpd_data$numNavs_mean <- (hcpd_data$numNavs_T1w + hcpd_data$numNavs_T2w) / 2
hcpd_data$numNavs_sum <- (hcpd_data$numNavs_T1w + hcpd_data$numNavs_T2w)

## Subject identifiers
hcpd_data$subject_id <- hcpd_data$Subject
hcpd_data$PIN <- as.character(hcpd_data$subject_id)
nsub <- dim(hcpd_data)[1]

for(i in 1:nsub) {
  tmp_id <- hcpd_data$subject_id[i]
  hcpd_data$PIN[i] <- as.character(paste0(tmp_id,'_V1'))
}

## Merge with Motion data
vnav_data <- read.csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n655_hcpd_vNav_motion_scores.csv")

nda_hcpd_data <- merge(hcpd_data, vnav_data, by="subject_id")

gb_hcpd_data <- subset(nda_hcpd_data, nda_hcpd_data$Age >= 8)

## Look at age distribution
hist(hcpd_data$Age, col="slategray3")
table(hcpd_data$Sex)

# hcpd_data$Subj.ID <- hcpd_data$Subject

######################################
## READ IN GLASSER'S CORRECTED MAPS ##
######################################
library(cifti)
library(data.table)
HCD_myelin_maps_cifti <- read_cifti("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/HCDMyelin_Aug2021/HCD628_Winter2021.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR.dscalar.nii")

myelin_maps <- HCD_myelin_maps_cifti$data
myelin_maps <- transpose(as.data.frame(myelin_maps))
mean_myelin_map <- unlist(lapply(myelin_maps, mean))

map_names <- HCD_myelin_maps_cifti$NamedMap
subject_id <- substr(map_names$map_names, 1, 10)

## Glasser-MMP corrected maps
GlasserMMP_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/parcellated_maps/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Glasser-MMP.txt")
GlasserMMP_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(GlasserMMP_transmitCorr_myelinMaps)))

## Schaefer 400x7 corrected Maps
Schaefer_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/parcellated_maps/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Schaefer_400x7.txt")
Schaefer_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(Schaefer_transmitCorr_myelinMaps)))

## Re-name columns for regional myelin
for(i in 1:400) {
  colnames(Schaefer_transmitCorr_myelinMaps)[i] <- paste("myelin_schaefer_v", i, sep = "")
}

Schaefer_transmitCorr_myelinMaps$subject_id <- subject_id
Schaefer_transmitCorr_myelinMaps$subject_id <- noquote(Schaefer_transmitCorr_myelinMaps$subject_id)

hcpd_data <- merge(Schaefer_transmitCorr_myelinMaps, hcpd_data, by="subject_id")
dim(hcpd_data)

## Regressed maps
# GlasserMMP_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_March2021/Partial.All.MyelinMap_IndPseudoCorr_Reg_MSMAll.32k_fs_LR_Glasser-MMP.txt")


## Re-name columns for regional myelin
for(i in 1:360) {
  colnames(GlasserMMP_transmitCorr_myelinMaps)[i] <- paste("myelin_glasser_v", i, sep = "")
}

GlasserMMP_transmitCorr_myelinMaps$subject_id <- subject_id

GlasserMMP_transmitCorr_myelinMaps$subject_id <- noquote(GlasserMMP_transmitCorr_myelinMaps$subject_id)


hcpd_data <- merge(GlasserMMP_transmitCorr_myelinMaps, hcpd_data, by="subject_id")
dim(hcpd_data)


###############################
## Read in Yeo-7 Myelin Maps ##
###############################
Yeo7_Network_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/parcellated_maps/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Yeo7-Network.txt")
# Yeo7_Network_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_March2021/regressed_maps/Partial.All.MyelinMap_IndPseudoCorr_Reg_MSMAll.32k_fs_LR_Yeo7-Network.txt")

Yeo7_Network_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(Yeo7_Network_transmitCorr_myelinMaps)))
dim(Yeo7_Network_transmitCorr_myelinMaps)

## Re-name columns for Yeo7 network myelin
colnames(Yeo7_Network_transmitCorr_myelinMaps) <- c("Medial_Wall", "Yeo3_myelin", "Yeo6_myelin", "Yeo7_myelin", "Yeo1_myelin", "Yeo5_myelin", "Yeo2_myelin", "Yeo4_myelin")

## Observe strange, unordered label conventions for Yeo7 dlabel file
# Yeo7_dlabel_cifti <- read_cifti("/Users/localadmin/Downloads/RSN-networks.32k_fs_LR.dlabel.nii")
# Yeo7_dlabel_cifti$NamedMap

Yeo7_Network_transmitCorr_myelinMaps$subject_id <- GlasserMMP_transmitCorr_myelinMaps$subject_id

hcpd_data <- merge(hcpd_data, Yeo7_Network_transmitCorr_myelinMaps, by="subject_id")
dim(hcpd_data)

########################################
## Read in Cole-Anticevic Myelin Maps ##
########################################
ColeAnt_12network_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/parcellated_maps/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_ColeAnt.txt")
# ColeAnt_12network_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/data/myelinMaps_March2021/regressed_maps/Partial.All.MyelinMap_IndPseudoCorr_Reg_MSMAll.32k_fs_LR_ColeAnt12-Network.txt")

ColeAnt_12network_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(ColeAnt_12network_transmitCorr_myelinMaps)))
dim(ColeAnt_12network_transmitCorr_myelinMaps)

## Re-name columns for Yeo7 network myelin
colnames(ColeAnt_12network_transmitCorr_myelinMaps) <- c("Visual1_myelin", "Visual2_myelin", "Somatomotor_myelin", "Cingulo_Opercular_myelin", "Dorsal_attention_myelin", "Language_myelin", "Frontoparietal_myelin", "Auditory_myelin", "Default_myelin", "Posterior_Multimodal_myelin", "Ventral_Multimodal_myelin", "Orbito_Affective_myelin")

ColeAnt_12network_transmitCorr_myelinMaps$subject_id <- subject_id

hcpd_data <- merge(hcpd_data, ColeAnt_12network_transmitCorr_myelinMaps, by="subject_id")
dim(hcpd_data)

## Setup model covariates
hcpd_data$Sex <- as.factor(hcpd_data$Sex)
hcpd_data$Scanner <- as.factor(hcpd_data$Scanner)
hcpd_data$scaled_numNavs_T1w <- scale(hcpd_data$numNavs_T1w, center = TRUE, scale = TRUE)
hcpd_data$scaled_numNavs_T2w <- scale(hcpd_data$numNavs_T2w, center = TRUE, scale = TRUE)
# hcpd_data$scaled_HeadSize <- scale(hcpd_data$HeadSize, center = TRUE, scale = TRUE)


##################################
## Read in Nuissance Covariates ##
##################################
headSize_df <- read.csv("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Dec2020/HeadSize.csv")
headSize_df$subject_id <- substr(headSize_df$Subject, 1, 10)

# glasser_covars <- read.csv("/ncf/hcp/data/analyses/myelin/data/myelinMaps_March2021/HCDMyelin/Covariates.csv")

glasser_covars <- read.csv("/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/HCDMyelin_Aug2021/Aug_Covariates.csv")
# glasser_covars$subject_id <- substr(glasser_covars$Subject, 1, 10)
glasser_covars$subject_id <- subject_id

## Read vNav motion metrics
vnav_data <- readr::read_csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/n655_hcpd_vNav_motion_scores.csv")

## Merge with brain data
hcpd_data <- merge(hcpd_data, vnav_data, by="subject_id")
hcpd_data <- merge(hcpd_data, headSize_df, by="subject_id")
hcpd_data <- merge(hcpd_data, glasser_covars, by="subject_id")

names(hcpd_data)
dim(hcpd_data)

############################
## Mean Wholebrain Myelin ##
############################ 
nreg <- 360
myelin_col_index <- grep("myelin_glasser_v", colnames(hcpd_data))

# across regions (within subjects)
tmp_myelin <-as.data.frame(hcpd_data[myelin_col_index])
mean_wholebrain_T1wT2w <- rowMeans(tmp_myelin)
hcpd_data$mean_wholebrain_T1wT2w <- mean_wholebrain_T1wT2w
hist(hcpd_data$mean_wholebrain_T1wT2w, col="slategray3")

## Create Threshold Based on Wholebrain T1wT2w ##
high_wholebrain_thresh <- mean(hcpd_data$mean_wholebrain_T1wT2w) + (4*sd(hcpd_data$mean_wholebrain_T1wT2w))
low_wholebrain_thresh <- mean(hcpd_data$mean_wholebrain_T1wT2w) - (4*sd(hcpd_data$mean_wholebrain_T1wT2w))

## Look at subjects with outlier data
myelin_thresh_idx <- which(hcpd_data$mean_wholebrain_T1wT2w > high_wholebrain_thresh | hcpd_data$mean_wholebrain_T1wT2w < low_wholebrain_thresh)
length(myelin_thresh_idx)
hcpd_data$subject_id[myelin_thresh_idx]

#######################################################
## REMOVE SUBEJECTS WITH WHOLEBRAIN T1w/T2w OUTLIERS ##
#######################################################
# hcpd_data <- hcpd_data[-myelin_thresh_idx,]
# dim(hcpd_data)

####################################
## Read in cortical thickness maps ##
#####################################
nsub <- dim(hcpd_data)[1]
nreg <- 360
glasser_thickness_maps <- as.data.frame(array(rep(NA, nsub*nreg), dim=c(nsub, nreg)))

for(i in 1:nsub) {
  sub_id <- hcpd_data$subject_id[i]
  glasser_thickness_path <- paste0("/users/gbaum/myelin_project/intradb_output/n632_thickness/", sub_id, "_thickness.32k_fs_LR_Glasser-MMP.txt")
  tmp_thickness <-  read.table(glasser_thickness_path, header=FALSE) 
  glasser_thickness_maps[i,] <- t(tmp_thickness)
}

## Re-name columns for regional myelin
for(i in 1:360) {
  colnames(glasser_thickness_maps)[i] <- paste("thickness_glasser_v", i, sep = "")
}

## Merge with hcpd_data
hcpd_data <- cbind(hcpd_data, glasser_thickness_maps)

#######################################
# Visualize Sample Characteristics ----
#######################################

## Graphics setup
apal <- paste0('#', c('000000', 'EAE3D8', 'FFFFFF', 'FFB838', '486D87'))

gbtheme <- theme_classic() +  
  theme(text = element_text(size = 36),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text =  element_text(color = apal[[1]]), axis.title = element_text(color = apal[[1]]),
        axis.ticks.length=unit(0.5, "cm"))

## Age by Sex: Proportion
age_sex_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=Sex, col=Sex)) + geom_bar(position = "fill", alpha=0.7) + scale_y_continuous(labels = scales::percent)
age_sex_plot + gbtheme + scale_fill_manual(values = c("#424C6D",  "#CE7B5B")) + scale_colour_manual(values = c("#424C6D",  "#CE7B5B"))  

## Age by Sex: Count
age_sex_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=Sex, col=Sex)) + geom_bar(alpha=0.7)
age_sex_plot + gbtheme + scale_fill_manual(values = c("#424C6D",  "#CE7B5B")) + scale_colour_manual(values = c("#424C6D",  "#CE7B5B"))  

## Age by SES: Proportion
## Age by SES
rm_idx <- which(is.na(hcpd_data$SES_RLVL)==TRUE)
ses_data <- hcpd_data[-rm_idx,]

age_ses_plot <- ggplot(ses_data, aes(x=age_floor, fill=SES_RLVL, col=SES_RLVL)) + geom_bar(position = "fill", alpha=0.7) + scale_y_continuous(labels = scales::percent)
age_ses_plot + gbtheme + scale_fill_manual(values = c("#424C6D", "#CE7B5B", "#FFD378")) + scale_colour_manual(values = c("#424C6D", "#CE7B5B", "#FFD378"))  

## Age by SES: Count
age_ses_plot <- ggplot(ses_data, aes(x=age_floor, fill=SES_RLVL, col=SES_RLVL)) + geom_bar(alpha=0.7)
age_ses_plot + gbtheme + scale_fill_manual(values = c("#424C6D", "#CE7B5B", "#FFD378")) + scale_colour_manual(values = c("#424C6D", "#CE7B5B", "#FFD378"))  

## Export data as csv
write.csv(hcpd_data, "/ncf/hcp/data/analyses/myelin/data/subject_demographics/n628_hcpd_newCorr_myelin_Aug2021.csv", row.names=FALSE)
write.table(hcpd_data$subject_id, "/ncf/hcp/data/analyses/myelin/subject_list/n628_hcpd_newCorr_subject_list.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
