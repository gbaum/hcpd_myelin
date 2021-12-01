#!/bin/sh

## Load Connectome Workbench module
module load connectome-workbench/1.3.2-fasrc01

## Define input data
# cifti_in=/ncf/hcp/data/analyses/myelin/data/myelinMaps_Dec2020/HCDMyelinMaps/HCD655_Summer2020.All.MyelinMap_Corr_MSMAll.32k_fs_LR.dscalar.nii

cifti_in=/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/HCDMyelin_Aug2021/HCD628_Winter2021.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR.dscalar.nii

cifti_in=/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360.pscalar.nii

# cifti_in=/ncf/hcp/data/analyses/myelin/data/myelinMaps_March2021/HCDMyelin/Partial.All.MyelinMap_IndPseudoCorr_Reg_MSMAll.32k_fs_LR.dscalar.nii

# Parcellations
glasser_atlas=/ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii
yeo7_atlas=/ncf/hcp/data/analyses/myelin/parcellations/RSN-networks.32k_fs_LR.dlabel.nii
cole_atlas=/ncf/hcp/data/analyses/myelin/parcellations/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii
schaefer_atlas=/ncf/hcp/data/analyses/myelin/parcellations/Schaefer2018_400Parcels_7Networks_order.dlabel.nii

# Define Output
outdir=/ncf/hcp/data/analyses/myelin/data/myelinMaps_Aug2021/parcellated_maps

cifti_out=${outdir}/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Glasser-MMP.pscalar.nii
text_out=${outdir}/Partial.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Glasser-MMP.txt

## Extract mean myelin values for each parcel in brain atlas
wb_command -cifti-parcellate ${cifti_in} ${glasser_atlas} COLUMN ${cifti_out} -method MEAN

## Output parcel values as text file
wb_command -cifti-convert -to-text ${cifti_out} ${text_out} 

###########################################################
## Create Text File for Sensorimotor-Association Axis in ##
##             Glasser-360 Parcellation                  ##
###########################################################
cifti_out=/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360.pscalar.nii

text_out=/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360.txt

## Output parcel values as text file
wb_command -cifti-convert -to-text ${cifti_out} ${text_out} 


