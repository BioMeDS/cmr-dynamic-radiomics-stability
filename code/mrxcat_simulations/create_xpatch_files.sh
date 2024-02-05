#!/usr/bin/env bash

for i in $(seq 1 24); 
    do
        gzip -d ../../data/mrxcat_simulations/mask_data/nifti_myc_cleared/$i.nii.gz;
        gzip -d ../../data/mrxcat_simulations/mask_data/nifti_myc/cine_act_$i.nii.gz;
        xdelta3 -s ../../data/mrxcat_simulations/mask_data/nifti_myc/cine_act_$i.nii ../../data/mrxcat_simulations/mask_data/nifti_myc_cleared/$i.nii ../../data/mrxcat_simulations/mask_data/nifti_myc_cleared/$i.xdelta;
        rm -f ../../data/mrxcat_simulations/mask_data/nifti_myc_cleared/$i.nii;
        rm -f ../../data/mrxcat_simulations/mask_data/nifti_myc/cine_act_$i.nii;
    done
