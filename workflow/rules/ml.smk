ACDC_PATIENTS = [f"patient{i:03d}" for i in range(1,101)]

rule ML_train:
    input:
        expand("analysis/features_normalized/ACDC/{patient}.csv", patient=ACDC_PATIENTS),
        "figures/supp_tab1.tsv",
        "data/ACDC_dataset/group.tsv",
    output:
        "figures/supp_tab2.tsv"
    conda:
        "../tidymodels.yaml"
    script:
        "../../code/rev/machine_learning_train.R"

rule ML_test:
    input:
        "figures/supp_tab2.tsv",
    output:
        "figures/supp_fig5.svg",
        "analysis/ml/correlation.tsv",
        "analysis/ml/performance.tsv"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/rev/machine_learning_test.R"