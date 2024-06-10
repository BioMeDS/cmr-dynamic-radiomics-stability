DIGITS = range(61, 81)
SEEDS = ('042', '069', '151', '401', '404')
NOISE = ("0.010", "0.020", "0.030", "0.040")
from glob import glob


wildcard_constraints:
    seeds="^(?!.*0.000)",
    noise="^(?!0)"

def aggregate_input_Acdc(wildcards):
    dir = checkpoints.SplitSlicesAcdc.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

def aggregate_input_noise_Acdc(wildcards):
    dir = checkpoints.SplitSlicesNoiseAcdc.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

checkpoint SplitSlicesAcdc:
    input:
        "data/ACDC_dataset/data/patient0{digits}/patient0{digits}_4d.nii.gz", 
        "data/ACDC_dataset/masks/patient0{digits}/pat0{digits}_masks_4d.nii.gz"
    output:
        directory("analysis/feature_extraction/ACDC/patient0{digits}"),
    params:
        "0"
    conda:
        "../split.yaml"
    script:
        "../../code/seperate_slices.py"


checkpoint SplitSlicesNoiseAcdc:
    input:
        "data/ACDC_dataset/data_noise/patient0{digits}_{noise}_{seed}/patient0{digits}_4d.nii.gz", 
        "data/ACDC_dataset/masks/patient0{digits}/pat0{digits}_masks_4d.nii.gz"
    output:
        directory("analysis/feature_extraction/ACDC_noise/patient0{digits}_{noise}_{seed}")
    params:
        "0"
    conda:
        "../split.yaml"
    script:
        "../../code/seperate_slices.py"

rule GenerateNoiseImagesAcdc:
    input:
        "data/ACDC_dataset/data/patient0{digits}/patient0{digits}_4d.nii.gz"
    output:
        "data/ACDC_dataset/data_noise/patient0{digits}_{noise}_{seed}/patient0{digits}_4d.nii.gz"
    conda:
        "../noise.yaml"
    params:
        "{seed}",
        "{noise}"
    script:
        "../../code/add_artifical_noise.py"
        
rule FeatureExtractionAcdc:
    input:
        aggregate_input_Acdc
    output:
        "analysis/features/ACDC_noise/patient0{digits}_0.000_0.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"
        
rule FeatureExtractionNoiseAcdc:
    input:
        aggregate_input_noise_Acdc
    output:
        "analysis/features/ACDC_noise/patient0{digits}_{noise}_{seed}.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"


rule AcdcFeatureNormalization:
    input:
        expand("analysis/features/ACDC_noise/patient0{{digits}}_{noise}_{seed}.csv", noise = NOISE, seed = SEEDS),
        "analysis/features/ACDC_noise/patient0{digits}_0.000_0.csv"
    output:
        "analysis/features_normalized/ACDC/patient0{digits}.csv"
    conda:
        "../tidyverse.yaml"
    params:
        "patient0{digits}_0.000_0.csv",
        "patient0{digits}"
    script:
        "../../code/mrxcat_simulations/feature_normalization.R"

rule AcdcCalculateMae:
    input:
        "analysis/features_normalized/ACDC/patient0{digits}.csv"
    output:
        "analysis/calculated_mae/ADDC/mae_patient0{digits}.csv"
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"