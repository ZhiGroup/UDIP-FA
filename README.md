# UDIP-FA: Unsupervised Deep Representation Learning of Fractional Anisotropy Maps

[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.07.04.25330856-blue)](https://doi.org/10.1101/2025.07.04.25330856)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

This repository contains the complete analysis pipeline for the study "Unveiling genetic architecture of white matter microstructure through unsupervised deep representation learning of fractional anisotropy maps".

<img width="449" height="648" alt="UDIP-FA Overview" src="https://github.com/user-attachments/assets/5c982686-4043-4448-b350-b54dd50b2dc9" />

## 📋 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## 🔬 Overview

This study introduces UDIP-FA (Unsupervised Deep Image Phenotyping of Fractional Anisotropy), a novel deep learning approach for analyzing white matter microstructure in brain imaging data. The pipeline includes:

- **Deep representation learning** of FA maps using unsupervised methods
- **Genome-wide association studies (GWAS)** on learned features
- **Polygenic risk score (PRS)** associations with brain disorders
- **Network-based drug targeting** analysis
- **Comprehensive visualization** and interpretation tools

## 🛠 Installation

### Prerequisites

- Python 3.9 or higher
- R 4.0 or higher
- Git

### Python Dependencies

```bash
pip install pandas numpy matplotlib seaborn
pip install scikit-learn lightgbm statsmodels scipy
pip install nilearn plotly umap-learn tqdm
pip install jupyter notebook ipython
```

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

### Clone Repository

```bash
git clone https://github.com/yourusername/UDIP-FA.git
cd UDIP-FA
```

## 🚀 Usage

### Quick Start

1. **Main Analysis Pipeline**:
   ```bash
   jupyter notebook FA_GWAS_all.ipynb
   ```

2. **Post-GWAS Analysis**:
   ```r
    FA_all.R
   ```

3. **Network and Drug Analysis**:
   ```r
    FA_network_drug_analysis.R
   ```

### Step-by-Step Analysis

The analysis follows a structured workflow:

1. **Data preprocessing and UDIP-FA feature extraction**
2. **Age and sex prediction using FA features**
3. **GWAS analysis on UDIP-FA components**
4. **PRS association studies**
5. **Biological pathway enrichment**
6. **Drug-target network analysis**

## 📁 Repository Structure

```
UDIP-FA/
├── FA_GWAS_all.ipynb          # Main analysis notebook
├── FA_all.R                   # Post-GWAS analysis script
├── FA_network_drug_analysis.R # Network and drug analysis
├── Model/                     # Pre-trained models and weights
├── README.md                  # This file
└── requirements.txt           # Python dependencies
```

### File Descriptions

#### `FA_GWAS_all.ipynb`
Comprehensive Jupyter notebook containing:
- UDIP-FA feature association analyses
- Polygenic risk score (PRS) associations with brain disorders
- Model explainability and interpretability assessments
- Comparative analysis with previous white matter studies

#### `FA_all.R`
R script for post-GWAS processing including:
- Statistical result aggregation and filtering
- Publication-ready figure generation
- Effect size calculations and meta-analysis

#### `FA_network_drug_analysis.R`
Advanced network analysis script featuring:
- Gene-drug interaction network construction
- Therapeutic target identification
- Mechanism of action pathway analysis

### Model
The pretrained model can access at this google driver[https://drive.google.com/file/d/1wPO-DoaXAD-kil6FZCNOGuWg9cOX9eql/view?usp=drive_link]
## 🔄 Reproducibility

### Environment Setup

We recommend using conda for environment management:

```bash
conda create -n udip-fa python=3.8
conda activate udip-fa
pip install -r requirements.txt
```

### Random Seeds

All analyses use fixed random seeds for reproducibility:
- Python: `np.random.seed(42)`
- R: `set.seed(42)`

### Version Information

The analysis was performed using:
- Python 3.9.10
- R 4.2.0
- Key package versions listed in `requirements.txt`

## 📈 Expected Outputs

Running the complete pipeline will generate:

- **GWAS Summary Statistics**: Tab-separated files with association results
- **PRS Association Results**: CSV files with polygenic score analyses
- **Visualization Figures**: High-resolution plots in PNG/PDF format
- **Network Analysis Results**: Gene-drug interaction networks
- **Statistical Reports**: Comprehensive analysis summaries

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
- **Project Link**: [https://github.com/yourusername/UDIP-FA]

## 🙏 Acknowledgments

- UK Biobank for providing the neuroimaging and genetic data
- The research computing facilities for computational resources
- Open-source software communities for tools and libraries used

---

**Keywords**: white matter, fractional anisotropy, deep learning, GWAS, neuroimaging, brain imaging, genetics, biomarker
