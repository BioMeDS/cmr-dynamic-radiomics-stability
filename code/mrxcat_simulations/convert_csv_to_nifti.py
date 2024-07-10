#!/usr/bin/env python
# coding: utf-8
import numpy as np
import nibabel as nib


def convert_csv_to_nifti(csv_file, nifti_file):
    """
    Convert a CSV file to a NIfTI file.

    Args:
        csv_file (str): Path to the CSV file.
        nifti_file (str): Path to the output NIfTI file.

    Returns:
        None
    """
    my_data = np.genfromtxt(f"{csv_file}",  delimiter=',')
    my_data = np.expand_dims(my_data, axis=2)
    nib.save(nib.Nifti1Image(my_data, np.eye(4)), f"{nifti_file}")

convert_csv_to_nifti(snakemake.input[0], snakemake.output[0])
