valid_slices_sub = {"Proband X1": [8],
                    "Proband X2": [8],
                    "Proband X4": [8],
                    "Proband X5": [8],
                    "Proband X6": [8],
                    "Proband X7": [8],
                    "Proband X8": [8],
                    "Proband X9": [9],
                    "Proband X10": [8],
                    "Proband X11": [8],
                    "Proband X12": [8],
                    "Proband X13": [9],
                    "Proband X14": [11],
                    "Proband X15": [7]}

def aggregate_input_Sub(wildcards):
    dir = checkpoints.SplitSlicesSub.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

def aggregate_input_noise_Sub(wildcards):
    dir = checkpoints.SplitSlicesNoiseSub.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

checkpoint SplitSlicesSub:
    input:
        lambda wc: glob("data/subject_data/{folder}/*[!mask].nii.gz".format(folder=wc.folder)), 
        lambda wc: glob("data/subject_data/{folder}/*_mask.nii.gz".format(folder=wc.folder)),
    output:
        directory("analysis/feature_extraction/subject/{folder}"),
    params:
        valid_slices_sub
    conda:
        "../split.yaml"
    script:
        "../../code/seperate_slices.py"


checkpoint SplitSlicesNoiseSub:
    input:
        "data/subject_data_noise/{folder}_{noise}_{seed}/{folder}.nii.gz", 
        lambda wc: glob("data/subject_data/{folder}/*_mask.nii.gz".format(folder=wc.folder)),
    output:
        directory("analysis/feature_extraction/subject_noise/{folder}_{noise}_{seed}"),
    params:
        valid_slices_sub
    conda:
        "../split.yaml"
    script:
        "../../code/seperate_slices.py"

rule GenerateNoiseImagesSub:
    input:
        lambda wc: glob("data/subject_data/{folder}/*[!mask].nii.gz".format(folder=wc.folder)),
    output:
        "data/subject_data_noise/{folder}_{noise}_{seed}/{folder}.nii.gz"
    conda:
        "../noise.yaml"
    params:
        "{seed}",
        "{noise}"
    script:
        "../../code/add_artifical_noise.py"
        
rule FeatureExtractionSub:
    input:
        aggregate_input_Sub
    output:
        "analysis/features/subject_noise/{folder}_0.000_0.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"
        
rule FeatureExtractionNoiseSub:
    input:
        aggregate_input_noise_Sub
    output:
        "analysis/features/subject_noise/{folder}_{noise}_{seed,\\d{3}}.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"


rule SubFeatureNormalization:
    input:
        expand("analysis/features/subject_noise/{{folder}}_{noise}_{seed}.csv", noise = NOISE, seed = SEEDS),
        "analysis/features/subject_noise/{folder}_0.000_0.csv"
    output:
        "analysis/features_normalized/subject/{folder}.csv"
    conda:
        "../tidyverse.yaml"
    params:
        "{folder}_0.000_0.csv",
        "{folder}_"
    script:
        "../../code/mrxcat_simulations/feature_normalization.R"

rule SubCalculateMae:
    input:
        "analysis/features_normalized/subject/{folder}.csv"
    output:
        "analysis/calculated_mae/subject/mae_{folder}.csv"
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"