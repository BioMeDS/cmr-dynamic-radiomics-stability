import numpy as np
import nibabel as nib

def segment_myocard(org_file, segmented_file):
    """
    Crop the original mask and segment the myocardium by choosing the right segmentation marker and save the segmented result.

    Parameters:
    org_file (str): The path to the original file.
    segmented_file (str): The path to save the segmented file.

    Returns:
    None
    """
    mask = nib.load(org_file)
    mask_data = mask.get_fdata()
    mask_data = mask_data[128:640, 256:768, :]
    mask_data = mask_data == 1
    mask_data = mask_data.astype('<f8')
    nib.save(nib.Nifti1Image(mask_data, np.eye(4)), segmented_file)

segment_myocard(snakemake.input[0], snakemake.output[0])
