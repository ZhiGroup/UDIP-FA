# UDIP-FA: Unsupervised Deep Representation Learning of Fractional Anisotropy Maps

[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.07.04.25330856-blue)](https://doi.org/10.1101/2025.07.04.25330856)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

This repository contains the complete analysis pipeline for the study **"Unveiling genetic architecture of white matter microstructure through unsupervised deep representation learning of fractional anisotropy maps"**.

## 📋 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [UDIP-FA Model Usage](#udip-fa-model-usage)
  - [Data Preparation](#data-preparation)
  - [Training](#training)
  - [Inference](#inference)
  - [Analysis](#analysis)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [Contact](#contact)

## 🔬 Overview

This study introduces **UDIP-FA** (Unsupervised Deep Image Phenotyping of Fractional Anisotropy), a novel deep learning approach for analyzing white matter microstructure in brain imaging data. The pipeline includes:

- **Deep representation learning** of FA maps using customized 3D AutoEncoders.
- **Genome-wide association studies (GWAS)** on learned endophenotypes.
- **Polygenic risk score (PRS)** associations with brain disorders.
- **Network-based drug targeting** analysis.

## 🛠 Installation

### Prerequisites

- Python 3.8 or higher
- R 4.0 or higher
- Git

### Python Dependencies

We recommend using a virtual environment (conda or venv).

```bash
# Create and activate environment
conda create -n udip-fa python=3.8
conda activate udip-fa

# Install dependencies from requirements.txt
pip install -r requirements.txt
```

*Note: Ensure you have a compatible PyTorch version for your CUDA driver installed.*

### R Dependencies

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "data.table", 
                   "ComplexHeatmap", "circlize", "RColorBrewer",
                   "cowplot", "ggpubr", "pheatmap"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE"))
```

## 📁 Repository Structure

```
UDIP-FA/
├── Model/                     # Deep Learning Model & Scripts
│   ├── model.py               # AutoEncoder Architecture (PyTorch)
│   ├── dataset.py             # Dataset Loading Logic
│   ├── Train.py               # Training Script (PyTorch Lightning)
│   ├── inference.py           # Inference Script for generating embeddings
│   └── model_compare.py       # Analysis & Visualization scripts
├── FA_GWAS_all.ipynb          # Main GWAS Analysis Notebook
├── FA_all.R                   # Post-GWAS Analysis (R)
├── FA_network_drug_analysis.R # Network & Drug Analysis (R)
├── requirements.txt           # Python Project Dependencies
└── README.md                  # Project Documentation
```

## 🧠 UDIP-FA Model Usage

The deep learning model is located in the `Model/` directory.

### Data Preparation

Input data should be Affine registered MRI images (NIfTI format).
Prepare a CSV file containing the paths to your images under a column named `mri_names` (or specify your column name during inference).

### Training

To train the AutoEncoder from scratch:

```bash
python Model/Train.py
```

*Note: `Model/Train.py` is configured to use PyTorch Lightning. Adjust hyperparameters (learning rate, batch size, GPUs) directly in the file or by modifying the `LitAutoEncoder` class.*

### Inference

To generate latent representation (endophenotypes) from trained models:

```bash
python Model/inference.py --input_csv /path/to/data.csv \
                          --checkpoint /path/to/model.ckpt \
                          --output_dir /path/to/results
```

**Common Arguments:**
- `--input_csv`: Path to CSV file with image paths.
- `--checkpoint`: Path to the `.ckpt` model file.
- `--output_dir`: Folder to save the output pickle files.
- `--device`: `cuda:0` or `cpu`.

### Analysis

For performing analysis on significant SNPs and feature correlations:

```bash
python Model/model_compare.py
```

This script includes functions to:
1. Plot significant SNPs across different thresholds.
2. Compute and visualize pairwise correlations (CCA, Pearson) between feature sets.

## 🔄 Reproducibility

### Pre-trained Models
The pretrained model can be accessed at this [Google Drive Link](https://drive.google.com/file/d/1wPO-DoaXAD-kil6FZCNOGuWg9cOX9eql/view?usp=drive_link).

### Random Seeds
- Python: `np.random.seed(42)`
- R: `set.seed(42)`

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{zhao2025udip,
  title={Unveiling genetic architecture of white matter microstructure through unsupervised deep representation learning of fractional anisotropy maps},
  author={Zhao, Xingzhong and Xie, Ziqian and He, Wei and Fornage, Myriam and Zhi, Degui},
  journal={medRxiv},
  year={2025},
  doi={10.1101/2025.07.04.25330856}
}
```

## 💬 Contact

- **Xingzhong Zhao** - [xingzhong.zhao@uth.tmc.edu]

---

**Keywords**: white matter, fractional anisotropy, deep learning, GWAS, neuroimaging, brain imaging, genetics, biomarker
