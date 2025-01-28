rule Figure2:
    input:
        "analysis/calculated_mae/mrxcat_simulation/mae.csv",
        "analysis/features_normalized/mrxcat_simulation/features.csv"
    output:
        "figures/fig2.svg"
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig2.R"

rule Figure3:
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

rule Figure4:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/fig4.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig4.R"

rule Figure5:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/fig5.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig5.R"
        
rule SuppFigure1:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/supp_fig1.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/supp_fig1.R"

rule SuppFigure2:
    input:
        "figures/supp_tab1.tsv"
    output:
        "figures/supp_fig2.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/supp_fig2.R"