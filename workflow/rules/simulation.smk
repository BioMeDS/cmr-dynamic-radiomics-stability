SNR = (5, 10, 20, 30)
REPLICATE = range(1, 6)
NUMBER = range(1, 25)
FILE = ("image", "segmentation")

out = expand("data/mrxcat_simulations/snr{snr}_{replicate}/csvs/cine_1x1x1mm_512x512x1x24x4_snr{snr}_fa90_bh{number}.csv", 
               snr=SNR, replicate=REPLICATE, number=NUMBER)
out_50 = expand("data/mrxcat_simulations/snr50_1/csvs/cine_1x1x1mm_512x512x1x24x4_snr50_fa90_bh{number}.csv", 
               number=NUMBER)

rule SimSaveAsCsv:
    input:
        "data/mrxcat_simulations"
    output:
        out, out_50
    conda:
        "../octave_env.yaml"
    shell:
        "cd code/mrxcat_simulations;"
        "octave SaveAsCsv.m --no-gui;"

rule SimConvertCsvToNifti:
    input:
        "data/mrxcat_simulations/snr{snr}_{replicate}/csvs/cine_1x1x1mm_512x512x1x24x4_snr{snr}_fa90_bh{number}.csv"
    output:
        "data/mrxcat_simulations/snr{snr}_{replicate}/nifti/cine_1x1x1mm_512x512x1x24x4_snr{snr}_fa90_bh{number}.nii.gz"
    conda:
        "../conversion_sim.yaml"
    script:
        "../../code/mrxcat_simulations/convert_csv_to_nifti.py"

rule SimConvertMaskBinToNifti:
    input:
        "data/mrxcat_simulations/mask_data/bin_files/cine_act_{number}.bin"
    output:
        "data/mrxcat_simulations/mask_data/nifti_after_conv/cine_act_{number}.nii.gz"
    conda:
        "../conversion_sim.yaml"
    script:
        "../../code/mrxcat_simulations/bin_convert.py"

rule SimMaskSegMyc:
    input:
        "data/mrxcat_simulations/mask_data/nifti_after_conv/cine_act_{number}.nii.gz"
    output:
        "data/mrxcat_simulations/mask_data/nifti_myc/cine_act_{number}.nii.gz"
    conda:
        "../conversion_sim.yaml"
    script:
        "../../code/mrxcat_simulations/segment_myocard.py"

rule UnpackNifti:
    input:
        "data/mrxcat_simulations/mask_data/nifti_myc/cine_act_{number}.nii.gz"
    output:
        "data/mrxcat_simulations/mask_data/nifti_myc/cine_act_{number}.nii"
    shell:
        "gzip -dk {input}"

rule RemovePapillaryMusclesFromMask:
    input:
        "data/mrxcat_simulations/mask_data/nifti_myc/cine_act_{number}.nii",
        "data/mrxcat_simulations/mask_data/nifti_myc_cleared/{number}.xdelta"
    output:
        "data/mrxcat_simulations/mask_data/nifti_myc_cleared/{number}.nii"
    script:
        "../xpatcher.sh"

rule RepackNifti:
    input:
        "data/mrxcat_simulations/mask_data/nifti_myc_cleared/{number}.nii"
    output:
        "data/mrxcat_simulations/mask_data/nifti_myc_cleared/{number}.nii.gz"
    shell:
        "gzip -k {input}"

rule CopySimFilesForExtraction:
    input:
        "data/mrxcat_simulations/snr{snr}_{replicate}/nifti/cine_1x1x1mm_512x512x1x24x4_snr{snr}_fa90_bh{number}.nii.gz",
        "data/mrxcat_simulations/mask_data/nifti_myc_cleared/{number}.nii.gz"
    output:
        "analysis/feature_extraction/mrxcat_simulation/snr{snr}_{replicate}/{number}/image.nii.gz",
        "analysis/feature_extraction/mrxcat_simulation/snr{snr}_{replicate}/{number}/segmentation.nii.gz",
    shell:
        "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}"

rule FeatureExtractionSimulation:
    input:
        expand("analysis/feature_extraction/mrxcat_simulation/snr{{snr}}_{{replicate}}/{number}/{file}.nii.gz",
               number=NUMBER, file=FILE)
    output:
        "analysis/features/mrxcat_simulation/snr{snr}_{replicate}.csv"
    conda:
        "../feature_extraction.yaml"
    script:
        "../../code/mrxcat_simulations/feature_extraction.py"

rule SimulationFeatureNormalization:
    input:
        expand("analysis/features/mrxcat_simulation/snr{snr}_{replicate}.csv", 
               snr=SNR, replicate=REPLICATE),
        "analysis/features/mrxcat_simulation/snr50_1.csv"
    output:
        "analysis/features_normalized/mrxcat_simulation/features.csv"
    conda:
        "../tidyverse.yaml"
    params:
        "snr50_1.csv",
        ""
    script:
        "../../code/mrxcat_simulations/feature_normalization.R"

rule SimulationCalculateMae:
    input:
        "analysis/features_normalized/mrxcat_simulation/features.csv"
    output:
        "analysis/calculated_mae/mrxcat_simulation/mae.csv"
    conda:
        "../calculate_mae.yaml"
    script:
        "../../code/mrxcat_simulations/calculate_mae.py"
 
rule SimTableRank:
    input:
        "analysis/calculated_mae/mrxcat_simulation/mae.csv"
    output:
        "analysis/tables/mrxcat_simulation/ranks.csv"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/single_ranks.R"

rule GeneratePlots:
    input:
        "analysis/features_normalized/mrxcat_simulation/features.csv", 
        "analysis/tables/mrxcat_simulation/ranks.csv"
    output:
        "analysis/plots/mrxcat_simulation/features_curves/top12_features.png",
        "analysis/tables/rank_table_simulation.csv",
        "analysis/plots/mrxcat_simulation/total_mae_vs_snr_mae.png",
        "analysis/plots/mrxcat_simulation/total_mae_vs_snr_mae_1.png",
        "analysis/plots/mrxcat_simulation/rank_barcode.png",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/mrxcat_simulations/feature_plots.R"