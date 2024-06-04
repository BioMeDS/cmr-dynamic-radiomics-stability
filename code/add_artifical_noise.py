import torchio as tio
import nibabel as nib
import torch
import numpy as np

def generate_noise(image_path, seed, noise, output_path):
    torch.manual_seed(seed)
    image_data = nib.load(image_path).get_fdata()
    file_max = image_data.max()
    std = float(noise) * file_max
    image_data = tio.ScalarImage(tensor=image_data)
    subject = tio.Subject(image=image_data)
    transforms = tio.RandomNoise(std=std)
    transformed_image = transforms(subject)
    nib.save(nib.Nifti1Image(np.array(transformed_image['image'].tensor), np.eye(4)), output_path)


generate_noise(snakemake.input[0], snakemake.params[0], snakemake.params[1], snakemake.output[0])