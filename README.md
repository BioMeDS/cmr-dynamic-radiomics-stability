# Stability analysis of dynamic radiomics features for cardiac magnetic resonance

## Introduction

This repository contains the source code for the paper

[*Stability of Dynamic Radiomics Features in Cardiac MRI* by Mike D. Klaus, Fabian Laqua, Bettina Baeßler, and Markus J. Ankenbrand. **medRxiv** 2025.04.18.25326051](https://doi.org/10.1101/2025.04.18.25326051)

## General workflow

![workflow](figures/fig1.png)

## Data sets

The large data in this repository is managed via [dvc](https://dvc.org/). That's why you see a lot of `.dvc` files in the folder structure. The content of the corresponding files is not directly included. For many of the MRXCAT and ACDC files and all the final results, these files are included in the [`dvcstore` deposited at Zenodo](https://doi.org/10.5281/zenodo.15239917). In order to use it, download and unzip it into your repository root. Then you can get files from the store with `dvc pull -r zenodo`. For example:

```bash
dvc pull -r zenodo analysis/features_normalized/ACDC/patient001.csv.dvc
```

### MRXCAT

Simulations with different SNR (5,10,20,30) using [MRXCAT](https://biomed.ee.ethz.ch/mrxcat.html) [1] version v1.4 and the breathhold phantom.
Code is in [code/mrxcat_simulations](code/mrxcat_simulations) and simulation results in [data/mrxcat_simulations](data/mrxcat_simulations).
The simulations were performed with Matlab version R2021b. Simulation results have to be converted from the cpx format to the nifti format for further analysis. Steps include:
1. .cpx to .csv convertion
2. .csv to .nifti convertion
The masks were then extracted from the breathhold phantom used for the simulation in MRXCAT.

### ACDC

The ACDC dataset was published as part of the [Automatic Cardiac Detection Challenge](https://www.creatis.insa-lyon.fr/Challenge/acdc/index.html) [2].

### BAE

Data from the publication *A systematic evaluation of three different cardiac T2-mapping sequences at 1.5 and 3T in healthy volunteers* [3].

## Running the analysis

The analysis is implemented as a snakemake workflow (see `Snakemake` and `.smk` files in `workflow/`). If you want to re-run it, you can call this command from the root of the repository:

```bash
snakemake --software-deployment-method conda --cores 8
```

The default target are all the figures and tables in `figures/` (including supplementary ones).
So this command will likely complain about missing input files (the BAE data).
Either edit the inputs for `rule all` manually in `workflow/Snakefile` or call snakemake with the outputs of interest, e.g.:

```bash
snakemake --software-deployment-method conda --cores 8 figures/fig1.svg
```


## References
1. Wissmann, L., Santelli, C., Segars, W.P. et al. MRXCAT: Realistic numerical phantoms for cardiovascular magnetic resonance. J Cardiovasc Magn Reson 16, 63 (2014). https://doi.org/10.1186/s12968-014-0063-3
2. Bernard, O., Lalande, A., Zotti, C., Cervenansky, F., et al. Deep Learning Techniques for Automatic MRI Cardiac Multi-structures Segmentation and Diagnosis: Is the Problem Solved? IEEE Transactions on Medical Imaging, 37, 11, 2514-2525 (2018). https://doi.org/10.1109/TMI.2018.2837502
3. Baeßler, B. et al. A systematic evaluation of three different cardiac T2-mapping sequences at 1.5 and 3T in healthy volunteers. Eur. J. Radiol. 84, 2161–2170 (2015). https://doi.org/10.1016/j.ejrad.2015.08.002