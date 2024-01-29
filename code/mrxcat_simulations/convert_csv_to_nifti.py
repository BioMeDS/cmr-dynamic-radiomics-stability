#!/usr/bin/env python
# coding: utf-8
import numpy as np
from numpy import genfromtxt
import nibabel as nib


def convert_csv_to_nifti(csv_file, nifti_file):
    my_data = genfromtxt(f"{csv_file}",  delimiter=',')
    my_data = np.expand_dims(my_data, axis=2)
    nib.save(nib.Nifti1Image(my_data, np.eye(4)), f"{nifti_file}")

# for i in range(0, len(snakemake.input)):
#     convert_csv_to_nifti(snakemake.input[i], snakemake.output[i])
convert_csv_to_nifti(snakemake.input[0], snakemake.output[0])
