import torchio as tio
import nibabel as nib
import torch
import numpy as np

def generate_noise(image_path, seed, noise, output_path):
    """
    Generate noisy image from the input image.

    Args:
        image_path (str): Path to the input image file.
        seed (int): Seed for random number generation.
        noise (str): Noise level to be applied to the image. Must be convertible to float.
        output_path (str): Path to save the generated noisy image.

    Returns:
        None
    """
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