VALID_SLICES = [6, 7, 8]
subject_data_rev = glob("data/subject_data/Proband X8/*[!mask].nii.gz")

FOLDER_REV = [os.path.basename(Path(x).parent) for x in subject_data]

def rev_aggregate_input_Sub(wildcards):
    dir = checkpoints.RevSplitSlicesSub.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

def rev_aggregate_input_noise_Sub(wildcards):
    dir = checkpoints.RevSplitSlicesNoiseSub.get(**wildcards).output[0]
    return glob(f"{dir}/*/*")

checkpoint RevSplitSlicesSub:
    input:
        lambda wc: glob("data/subject_data/{folder}/*[!mask].nii.gz".format(folder=wc.folder)), 
        lambda wc: glob("data/subject_data/{folder}/*_mask.nii.gz".format(folder=wc.folder)),
    output:
        directory("analysis/multiple_slices/feature_extraction/subject/{folder}_{valid_slices}"),
    params:
        "{valid_slices}"
    conda:
        "../split.yaml"
    script:
        "../../code/rev/seperate_slices.py"


checkpoint RevSplitSlicesNoiseSub:
    input:
        "data/subject_data_noise/{folder}_{noise}_{seed}/{folder}.nii.gz", 
        lambda wc: glob("data/subject_data/{folder}/*_mask.nii.gz".format(folder=wc.folder)),
    output:
        directory("analysis/multiple_slices/feature_extraction/subject_noise/{folder}_{noise}_{seed}_{valid_slices}"),
    params:
        "{valid_slices}"
    conda:
        "../split.yaml"
    script:
        "../../code/rev/seperate_slices.py"
        
rule RevFeatureExtractionSub:
    input:
        rev_aggregate_input_Sub
    output:
        "analysis/multiple_slices/features/subject_noise/{valid_slices}/{folder}_0.000_0.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"
        
rule RevFeatureExtractionNoiseSub:
    input:
        rev_aggregate_input_noise_Sub
    output:
        "analysis/multiple_slices/features/subject_noise/{valid_slices}/{folder}_{noise}_{seed,\\d{3}}.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"


rule RevFeatureNormalization:
    input:
        expand("analysis/multiple_slices/features/subject_noise/{{valid_slices}}/{{folder}}_{noise}_{seed}.csv", noise = NOISE, seed = SEEDS),
        "analysis/multiple_slices/features/subject_noise/{valid_slices}/{folder}_0.000_0.csv"
    output:
        "analysis/multiple_slices/features_normalized/subject/{valid_slices}/{folder}.csv"
    conda:
        "../tidyverse.yaml"
    params:
        "{folder}_0.000_0.csv",
        "{folder}_"
    script:
        "../../code/mrxcat_simulations/feature_normalization.R"

rule RevCalculateMae:
    input:
        "analysis/multiple_slices/features_normalized/subject/{valid_slices}/{folder}.csv"
    output:
        "analysis/multiple_slices/calculated_mae/subject/{valid_slices}/mae_{folder}.csv"
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"

rule RevGenerateComp:
    input:
        expand("analysis/multiple_slices/calculated_mae/subject/{valid_slices}/mae_{folder}.csv", folder = FOLDER_REV, valid_slices = VALID_SLICES),
    output:
        "figures/rev/slice6vs7.svg",
        "figures/rev/slice6vs8.svg",
        "figures/rev/slice7vs8.svg",
        "figures/rev/spearman.svg"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/rev/plots.R"

rule RevTableRank:
    input:
        "analysis/multiple_slices/calculated_mae/subject/{valid_slices}/mae_{folder}.csv"
    output:
        "analysis/multiple_slices/tables/subject/{valid_slices}/ranks_{folder}.csv"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/single_ranks.R"

rule Rev:
    input:
        expand("analysis/multiple_slices/tables/subject/{valid_slices}/ranks_{folder}.csv", folder = FOLDER_REV, valid_slices = VALID_SLICES),
        "figures/rev/slice6vs7.svg",
        "figures/rev/slice6vs8.svg",
        "figures/rev/slice7vs8.svg",
        "figures/rev/spearman.svg"


# Rules for resampled image analysis

subject_data = glob("data/subject_data/*/*[!mask].nii.gz")
FOLDER = [os.path.basename(Path(x).parent) for x in subject_data]

rule SubResampleAndNormalize:
    input:
        expand("analysis/features/subject_noise/{{folder}}_{noise}_{seed}.csv", noise = NOISE, seed = SEEDS),
        "analysis/features/subject_noise/{folder}_0.000_0.csv"
    output:
        "analysis/resampled_features/subject/{folder}.csv"
    conda:
        "../tidyverse.yaml"
    params:
        "{folder}_0.000_0.csv",
        "{folder}_"
    script:
        "../../code/rev/feature_resample_normalization.R"


rule SubResampledCalculateMae:
    input:
        "analysis/resampled_features/subject/{folder}.csv"
    output:
        "analysis/resampled_features/calculated_mae/subject/mae_{folder}.csv"
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"

rule SubResampledComp:
    input:
        expand("analysis/resampled_features/calculated_mae/subject/mae_{folder}.csv", folder = FOLDER),
        expand("analysis/calculated_mae/subject/mae_{folder}.csv", folder = FOLDER),
    output:
        "figures/rev/resampled_vs_unresampled.svg",
        "figures/rev/resampled_unresampled_rank_comparison.svg"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/rev/resampled.R"

rule SubResampledTableRank:
    input:
        "analysis/resampled_features/calculated_mae/subject/mae_{folder}.csv"
    output:
        "analysis/resampled_features/tables/subject/ranks_{folder}.csv"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/single_ranks.R"

rule SubResampled:
    input:
        expand("analysis/resampled_features/tables/subject/ranks_{folder}.csv", folder = FOLDER),
        "figures/rev/resampled_vs_unresampled.svg",
        "figures/rev/resampled_unresampled_rank_comparison.svg"
