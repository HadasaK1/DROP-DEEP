DRIP: Dimensionality Reduction Informs Polygenic Scoring

DRIP is a two-stage framework for polygenic risk prediction that combines:

Phenotype-agnostic dimensionality reduction
Phenotype-specific machine learning prediction

Unlike conventional GWAS-based PRS pipelines, DRIP avoids supervised SNP selection and instead learns compact genomic representations directly from genotype data.

The method was evaluated on UK Biobank data across multiple binary and continuous phenotypes and demonstrated competitive predictive performance with substantially reduced computational runtime.

Paper

DRIP: Dimensionality Reduction Informs Polygenic Scoring
Hadasa Kaufman, Yarden Hochenberg, Michal Linial, Nadav Rappoport

Preprint / publication link: TBD

Overview

The DRIP framework consists of two stages:

Stage 1 — Dimensionality Reduction

Genome-wide SNP data are compressed using:

Principal Component Analysis (PCA)
Autoencoders

The dimensionality reduction step is phenotype-independent and can therefore be reused across multiple traits.

Stage 2 — Prediction Models

Reduced genomic representations are used as input for:

Logistic / Linear Regression
XGBoost
Deep Neural Networks
Repository Structure
data_processing/         Quality control and phenotype creation
dimensionality_reduction/ PCA and autoencoder pipelines
prediction_models/       ML prediction models
baseline_methods/        PRSice-2, PRS-CS, lassosum comparisons
figures/                 Paper figures
notebooks/               Reproducible notebooks
scripts/                 SLURM/bash execution scripts
configs/                 Hyperparameter and path configuration
Data Availability

UK Biobank data are available through application to the
UK Biobank.

Due to UK Biobank restrictions, genotype and phenotype data are not distributed in this repository.

Application ID used in this study: 26664.

Pretrained PCA Transformers

Pretrained chromosome-specific PCA transformers used in the manuscript are available at:

[Google Drive / Zenodo link]

Example structure:

rep1/
  PCA_transformer_1.pkl
  PCA_transformer_2.pkl
  ...
Installation

Clone the repository:

git clone https://github.com/nadavlab/DRIP.git
cd DRIP

Install dependencies:

pip install -r requirements.txt
Running DRIP
1. Quality control
bash scripts/run_qc.sh
2. Train PCA transformers
bash scripts/run_pca.sh
3. Generate reduced genomic representations
bash scripts/run_projection.sh
4. Train prediction models
bash scripts/run_lr.sh
Hardware Requirements

Experiments were conducted using:

UK Biobank genotype data
~342k individuals
~467k SNPs
RTX6000 GPU (58GB RAM) for autoencoders
High-memory CPU nodes for PCA and GWAS analyses
Main Findings
DRIP achieves predictive performance comparable to state-of-the-art PRS methods
PCA-based reduction outperformed more complex autoencoder approaches
DRIP substantially reduces computational runtime
The framework generalizes across diverse phenotypes
Citation

If you use this repository, please cite:

@article{kaufman2026drip,
  title={DRIP: Dimensionality Reduction Informs Polygenic Scoring},
  author={Kaufman, Hadasa and Hochenberg, Yarden and Linial, Michal and Rappoport, Nadav},
  journal={},
  year={2026}
}
Contact

Hadasa Kaufman
The Hebrew University of Jerusalem

For questions or collaborations:
hadasa.kaufman@mail.huji.ac.il
