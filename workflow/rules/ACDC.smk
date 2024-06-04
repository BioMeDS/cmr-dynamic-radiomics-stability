DIGITS = range(61, 81)
TIME = range(0, 30)

rule SplitSlices:
    input:
        # "data/ACDC_dataset/data/patient{digits}/patient0{digits}_4d.nii.gz", 
        # "data/ACDC_dataset/masks/patient{digits}/pat0{digits}_masks_4d.nii.gz"
        "data/ACDC_dataset/data/patient061/patient061_4d.nii.gz", 
        "data/ACDC_dataset/masks/patient061/pat061_masks_4d.nii.gz"
    output:
        # expand("analysis/feature_extraction/ACDC/patient0{digits}/{time}/image.nii.gz",
        # digits = DIGITS, time = TIME),
        # expand("analysis/feature_extraction/ACDC/patient0{digits}/{time}/segmentation.nii.gz",
        # digits = DIGITS, time = TIME),
        "analysis/feature_extraction/ACDC/patient061/1/image.nii.gz",
        "analysis/feature_extraction/ACDC/patient061/1/segmentation.nii.gz",
    params:
        "0"
    conda:
        "../split.yaml"
    script:
        "../../code/seperate_slices.py"