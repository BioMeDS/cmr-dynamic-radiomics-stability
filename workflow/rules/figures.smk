rule GenerateFigure2:
    input:
        "analysis/calculated_mae/mrxcat_simulation/mae.csv",
        "analysis/features_normalized/mrxcat_simulation/features.csv"
    output:
        "figures/fig2.svg"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig2.R"