# The Genetic Architecture of Savant Abilities Differs Across Cognitive Domains in Autism

This repository contains the code used to generate the analyses, tables and figures reported in the manuscript *The Genetic Architecture of Savant Abilities Differs Across Cognitive Domains in Autism*.

## Repository contents

| Directory | Contents |
|---|---|
| `data_preparation/` | Construction of the analysis table and cohort comparison tables |
| `reported_analysis/` | Analyses reported in the main and supplementary tables |
| `figure_pipeline/` | Figure-generation scripts and an independent implementation of the reported models |
| `figure_pipeline/figure_data/` | Aggregate inputs for Figures 1–3 |

The files in `figure_pipeline/figure_data/` contain aggregate model results and latent class summaries. Individual-level data are not included.

## Reproduction

Figures 1–3 can be reproduced directly from the included aggregate inputs:

```bash
Rscript figure_pipeline/Figure_v3.0.R \
  --analysis-source figure_pipeline/figure_data \
  --fig-dir figures
```

The full analysis requires the controlled-access phenotype, genotype and clinical source files described in the paper. The main entry points are:

```bash
python data_preparation/build_analysis_master_v1.py \
  --subject-table PATH --ssc-m13 PATH --ssc-m4 PATH \
  --ssc-pheno PATH --korean-pheno PATH --korean-clinical PATH \
  --supertable PATH --output-dir PATH

Rscript reported_analysis/savant_v3.R PATH/analysis_master_subject_level_v3.tsv
Rscript figure_pipeline/reanalyse_savant_v2.9.R \
  --master PATH --out PATH --reference-dir PATH
Rscript figure_pipeline/run_pooled_release_v3.0.R \
  --master PATH --legacy-runner PATH --reference-dir PATH --out PATH
```

The two comparison scripts in `reported_analysis/` read their input tables from the working directory. The preparation scripts expect the source files under a project-level `Data/` directory and write to `Tables/Source_Data/`.

## Implementation notes

The latent class model uses 300 random restarts, a maximum of 20,000 expectation-maximization iterations, a convergence tolerance of 1 × 10⁻¹² and seed 2900. After fitting, confirm that the four class labels match their item-response profiles because class numbering depends on initialization.

The SSC comparison requires calibrated severity scores covering ADOS modules 1–4. The module 1–3 file alone excludes 48 module 4 participants; `build_analysis_master_v1.py` combines the two sources and checks their overlap.

Primary models adjust for age, sex, cohort, ADOS calibrated severity score, Vineland composite score and ancestry principal components 1–5. The sensitivity analysis substitutes the SRS total T-score for ADOS severity. De novo burden models omit the principal components. False discovery rate correction is applied within each reported table using the Benjamini–Hochberg method.

## Software

The R scripts use `data.table`, `dplyr`, `broom`, `poLCA`, `nnet`, `ggplot2`, `patchwork`, `scales`, `openxlsx`, `readr`, `readxl`, `tidyr`, `ragg`, `svglite`, `showtext` and `sysfonts`. The Python script uses `pandas` and `numpy`.

## Data availability

Individual-level participant data are not included in this repository. Simons Simplex Collection data are available to qualified researchers through [SFARI Base](https://base.sfari.org). Whole-genome sequencing data from the Korean families are available through the [MSSNG database](https://research.mss.ng). Additional restrictions and access conditions are described in the paper.

This repository corresponds to the version used for the submitted manuscript.
