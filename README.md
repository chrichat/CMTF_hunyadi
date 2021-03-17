# CMTF_hunyadi
EEG-fMRI fusion of epilepsy data via CMTF

## Description
EEG data (average interictal discharges per patient) are represented in a tensor of patients x channels x time. fMRI data are represented as a matrix of patient x voxels. Voxels are vectorized images, obtained by traditional GLM-based EEG-fMRI images, then masked to reduce the number of voxels. Then, joint blind source separation is performed by 4 different approaches: jointICA, temporal-jointICA, coupled tensor-matrix factorization (CMTF) and restricted CMTF, where the patient signatures are fixed according to the amplitude of spikes and slow waves in the EEG.

## Reference
The code executes the method and reproduces the results and figures described in: 
[1] Hunyadi B., Van Paesschen W., De Vos M., Van Huffel S., ``Coupled tensor-matrix factorization of EEG and fMRI to explore epileptic network activity'', Proc. of the 23rd European Signal Processing Conference 2016 (EUSIPCO), Budapest, Hungary, Sept. 2016.

## Key words/tags
EEG, fMRI Epilepsy, tensor decomposition, blind source separation, multimodal fusion

# Describing the Code 

jica_vs_jmtf_group_tle:

Loads EEG-fMRI data from 10 temporal lobe epilepsy patients. Average EEG waveforms and GLM-based activation maps are taken per patient. The activation maps are masked based on the thresholded average activation map.
 
JointICA, t-jointICA, CMTF and CMTF with fixed patient signatures are compared. For the latter the patient signatures are fixed based on the spike and slow wave amplitude of each patient.
 
Saves signatures, plots EEG temporal signatures and plots the distribution of the fMRI voxel signatures.

jica_vs_jmtf_group_tle_testpat:

Reproducability analysis: the same code above is run multiple times by taking a subset of patients at a time (6-9 patients, all possible subsets). The resulting signatures are saved in .mat files

fusion_resproducability:

Loads the signatures computed with subsets of patients and quantifies their reproducibility as average pairwise correlation. The results are plotted.

Toolbox:
This folder contains all necessary tools for running the previous scripts: 
Spm8
Tensorlab v3
GroupICA
FastICA
Plotting tools


## Running the Code 

After downloading all files and preserving the directory structure, jica_vs_jmtf_group_tle and jica_vs_jmtf_group_tle_testpat can be run as scripts either from command line or pushing the run button in Matlab editor. Setting of the path and loading data are done automatically. fusion_resproducability is called from within jica_vs_jmtf_group_tle_testpat.


## Program Output 
The outputs of the scripts are visualised in Matlab figures or saved in .mat files and will be used in subsequent scripts. In order to interpret the figures, please see the above mentioned paper.


