import skimage as ski
import nibabel as nib
import sys

i1 = sys.argv[1]
i2 = sys.argv[2]

if i1 == i2:  # psnr between identical images is not meaningful
    exit(0)

im1 = nib.load(i1).get_fdata()
im2 = nib.load(i2).get_fdata()

psnr = ski.metrics.peak_signal_noise_ratio(im1, im2, data_range=im1.max())

print(f"{i1}\t{i2}\t{psnr}")
