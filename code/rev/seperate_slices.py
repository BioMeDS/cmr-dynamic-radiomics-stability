import nibabel as nib
import numpy as np
import re
import os
from pathlib import Path

def get_specific_slice(image_path, mask_path, save_loc, given_slice):
    """
    Extracts a specific slice from a 4D image volume and saves it as a 3D NIfTI file.

    Args:
        image_path (str): The path to the 4D image volume.
        mask_path (str): The path to the 4D mask volume.
        save_loc (str): The directory where the extracted slices will be saved.
        given_slice (int, optional): A number giving the slice to be extracted
            If not provided, the function will determine the middle slice of the image volume.

    Raises:
        AssertionError: If the slice is not an integer.

    """
    slice = False
    image_path = Path(image_path)
    mask_path = Path(mask_path)
    if type(given_slice) != int:
        slice = int(given_slice)
    if type(slice) != int:
        assert("Slice isn't an int")
    image_nim = nib.load(image_path)
    mask_nim = nib.load(mask_path)
    image_data = image_nim.get_fdata()
    mask_data = mask_nim.get_fdata()
    mask_data = (mask_data == 2).astype('<f8')
    _, _, _, t1 = image_data.shape
    _, _, _, t2 = mask_data.shape
    for t in range(t1):
        image_data_mod = image_data[:, :, slice, t]
        image_data_mod = np.expand_dims(image_data_mod, axis=2)
        image_data_mod = nib.Nifti1Image(image_data_mod, np.eye(4))
        if not os.path.isdir(f"{save_loc}/{image_path.parent.name}_{slice}/{t}"):
            os.makedirs(f"{save_loc}/{image_path.parent.name}_{slice}/{t}")
        nib.save(image_data_mod, f"{save_loc}/{image_path.parent.name}_{slice}/{t}/image.nii.gz")
    for t in range(t2):
        mask_data_mod = mask_data[:, :, slice, t]
        mask_data_mod = np.expand_dims(mask_data_mod, axis=2)
        mask_data_mod = nib.Nifti1Image(mask_data_mod, np.eye(4))
        nib.save(mask_data_mod, f"{save_loc}/{image_path.parent.name}_{slice}/{t}/segmentation.nii.gz")

get_specific_slice(snakemake.input[0], snakemake.input[1], Path(snakemake.output[0]).parent, snakemake.params[0])
