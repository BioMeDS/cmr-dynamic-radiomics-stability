import numpy as np
import nibabel as nib

def convert_mask_to_myc(org_file, converted_file):
        mask = nib.load(org_file)
        mask_data = mask=mask.get_fdata()
        mask_data = mask_data[128:640,256:768,:]
        mask_data = mask_data == 1
        mask_data = mask_data.astype('<f8')
        nib.save(nib.Nifti1Image(mask_data, np.eye(4)), converted_file)

convert_mask_to_myc(snakemake.input[0], snakemake.output[0])
