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

rule GenerateFigure4:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/fig4.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig4.R"

rule GenerateFigure5:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/fig5.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/fig5.R"
        
rule GenerateSuppFigure1:
    input:
        "analysis/combined/ranks.tsv"
    output:
        "figures/supp_fig1.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/supp_fig1.R"

rule GenerateSuppFigure2:
    input:
        "analysis/calculated_mae/summary.csv"
    output:
        "figures/supp_fig2.svg",
    conda:
        "../tidyverse.yaml"
    script:
        "../../code/figures/supp_fig2.R"