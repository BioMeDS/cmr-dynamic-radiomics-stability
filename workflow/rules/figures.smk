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

rule GenerateFigure3:
    input:
        "analysis/combined/features.tsv",
        "analysis/combined/ranks.tsv"
    output:
        "figures/fig3.svg",
        "figures/fig3_supplement.pdf"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig3.R"