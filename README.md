# Thermothelomyces TnSeq Project

**Genome-wide functional genomics and temperature adaptation in the thermophilic fungus *Thermothelomyces thermophilus***

This repository contains all Python and Jupyter-based code used to process, analyze, and visualize transposon mutagenesis (TnSeq / BarSeq) and phenotypic data.  
It provides a reproducible framework for computing barcode abundance, gene-level fitness, temperature-dependent germination efficiency, and downstream figure generation.

---

## 🧬 Overview

Thermophilic fungi represent the upper limits of eukaryotic life and serve as ideal models for dissecting mechanisms of high-temperature growth, stress adaptation, and industrial enzyme production.  
This project uses **transposon mutagenesis** combined with **barcode sequencing (BarSeq)** and **phenotypic assays** to identify genes influencing thermotolerance, germination, and metal-nutrient crosstalk in *T. thermophilus*.

The analytical workflow supports:
- Barcode counting and normalization from FASTQ inputs  
- Gene-level fitness estimation using statistical tests (Mann–Whitney U, Wilcoxon signed-rank, FDR correction)  
- Correlation and reproducibility analysis across biological replicates  
- Integration of BarSeq and RNA-seq datasets  
- Visualization (volcano plots, heatmaps, scatter plots, PCA)

---

## 🧰 Repository Structure

Thermothelomyces-Tnseq-project/
│
├── notebooks/ # Jupyter notebooks for data analysis & figure generation
├── scripts/ # Python scripts for data processing and statistics
├── config/ # Optional configs, metadata tables
├── requirements.txt # List of dependencies
├── README.md # Project documentation (this file)
└── .gitignore # Exclusion rules for large data files

> Raw sequencing data (FASTQ, BAM, etc.) and large intermediate outputs are excluded to keep this repository lightweight and reproducible.  

---

## ⚙️ Environment Setup

### Using `conda` (recommended)
```bash
conda create -n tnseq python=3.10
conda activate tnseq
pip install -r requirements.txt
Using pip
python3 -m venv tnseq_env
source tnseq_env/bin/activate
pip install -r requirements.txt
Dependencies (as listed in requirements.txt):
pandas
numpy
scipy
statsmodels
matplotlib
seaborn
plotly
Statistical analysis & plots

Run Jupyter notebooks inside notebooks/ to generate:

Volcano plots of gene-level effects

Replicate correlation plots (R² values)

Germination efficiency heatmaps
Reproducibility Notes

All figures in the associated manuscript were generated from the final notebooks in this repository.
Scripts are deterministic and versioned, ensuring reproducibility across computational environments.

🧑‍💻 Contact

Dr. Olusola A. Ogunyewo
Postdoctoral Researcher – Rachel Brem Lab, UC Berkeley

📧 olusolaogunyewo@gmail.com

🔗 https://github.com/oaogunyewo-git
