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
        "data/ACDC_dataset/data_noise/patient0{digits}_{noise}_{seed,\\d{3}}/patient0{digits}_4d.nii.gz"
    conda:
        "../noise.yaml"
    params:
        "{seed,\\d{3}}",
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
        "analysis/features/ACDC_noise/patient0{digits}_{noise}_{seed,\\d{3}}.csv"
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
        "patient0{digits}_"
    script:
        "../../code/mrxcat_simulations/feature_normalization.R"

rule AcdcCalculateMae:
    input:
        "analysis/features_normalized/ACDC/patient0{digits}.csv"
    output:
        "analysis/calculated_mae/ACDC/mae_patient0{digits}.csv",
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"

rule AcdcGeneratePlots:
    input:
        "analysis/features_normalized/ACDC/patient0{digits}.csv"
    output:
        "analysis/plots/ACDC/features_curves/patient0{digits}/top12_features.png"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/feature_plots.R"

rule AcdcGenerateMaePlots:
    input:
        expand("analysis/calculated_mae/ACDC/mae_patient0{digits}.csv", digits = DIGITS)
    output:
        "analysis/tables/rank_table_ACDC.csv",
        "analysis/plots/ACDC/total_mae_vs_snr_mae.png",
        "analysis/plots/ACDC/total_mae_vs_snr_mae_1.png",
        "analysis/plots/ACDC/rank_barcode.png",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/mae_plots.R"

rule AcdcTableRank:
    input:
        "analysis/calculated_mae/ACDC/mae_patient0{digits}.csv"
    output:
        "analysis/tables/ACDC/ranks_patient0{digits}.csv"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/single_ranks.R"