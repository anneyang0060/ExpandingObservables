# Reproduction Repository for Online Sparse Regression with Expanding Observables

This repository contains all the code and data necessary to reproduce the results and figures from the paper **"Online Sparse Regression with Expanding Observables"**. It includes the simulation studies (Section 4), real-data analyses (Sections 5.1 and 5.2), and the supplementary material.

## 📁 Repository Structure
├── FNS.R                          # Core functions used across analyses

├── simu/                          # Main simulation study (Section 4) and the supplementary simulation

│   ├── sim.R                      # Master script to run the simulation

│   ├── R-AVAS_FM-ID.R             # Main simulation algorithm

│   ├── Setting_FM-ID.R            # Configuration of simulation parameters

│   ├── organize.R                 # Generates Figures 3 for the paper

│   └── organize_AVAS.R            # Generates Figures in Section 3 of the supplement

├── realdata/                      # Real-data analyses (Section 5)

│   ├── PM2.5/                     # Analysis for Section 5.1 (PM2.5 Data)

│   │   ├── preprocess_for_CMFD.R       # Preprocessing for CMFD data

│   │   ├── preprocess_for_CHAP.R       # Preprocessing for CHAP data

│   │   ├── preprocess_for_PM2.5.R     # Preprocessing for PM2.5 data

│   │   └── RAVAS_for_PM.5.R           # Main analysis producing Figures 4-5

│   └── P2P/                       # Analysis for Section 5.2 (P2P Data)

│       ├── preprocess.R            # Data preprocessing

│       ├── RAVAS_P2P.R             # Main analysis script

│       └── organize.R              # Produces Table 1 for the paper

└── supp/                          # Supplementary material


## 🔧 Detailed Reproduction Guide

### 1. Main Simulation (Section 4)

To reproduce the simulation study from Section 4 of the paper, located in the `simu/` directory:

1.  **Configure Parameters (Optional):** Adjust simulation settings in `simu/Setting_FM-ID.R` if needed.
2.  **Run the Simulation:** Execute the master script `simu/sim.R`. This will call the main algorithm `R-AVAS_FM-ID.R`.
3.  **Generate Figures:** Run `simu/organize.R` to produce the final figures (Figures 3 from the paper).

> **Note:** The scripts in the `simu/` directory can also generate Figures 1 and 2 from the supplementary material by modifying the relevant parameters in `Setting_FM-ID.R` [1](@ref).

### 2. Real-Data Analyses (Section 5)

#### 2.1 PM2.5 Analysis (Section 5.1)

To reproduce the analysis and figures for the PM2.5 data (Section 5.1) in `realdata/PM2.5/`:

0.  **Download Data** from https://data.tpdc.ac.cn/zh-hans/data/d3878788-7bef-4249-9973-07d4b2e46cb8 and https://data.tpdc.ac.cn/home.
1.  **Data Preprocessing:** Run the preprocessing scripts in the following order. Ensure the paths to raw data files are correct within each script.
*   `preprocess_for_CMFD.R`
*   `preprocess_for_CHAP.R`
*   `preprocess_for_PM2.5.R`
2.  **Run Main Analysis:** Execute `RAVAS_for_PM.5.R` to perform the primary analysis.
3.  **Output:** This will generate the results for **Figures 4 and 5** of the paper.

#### 2.2 P2P Analysis (Section 5.2)

To reproduce the analysis for the P2P data (Section 5.2) in `realdata/P2P/`:

0.  **Download Data** from https://www.kaggle.com/datasets/sid321axn/bondora-peer-to-peer-lending-loan-data/data.
1.  **Data Preprocessing:** Run `preprocess.R`.
2.  **Run Main Analysis:** Execute `RAVAS_P2P.R`.
3.  **Generate Table:** Run `organize.R` to compile the results into **Table 1** of the paper.

### 3. Supplementary Material

To reproduce the additional simulation from the supplementary material, located in the `supp/` directory:

1.  **Run Supplementary Simulation:** Execute `supp/sim.R`.
2.  **Generate Supplementary Figure:** Run `supp/organize.R` to produce **Figure 3** from the supplement.

## 📊 Output Summary

| Code Path | Corresponds to Paper Section | Output Generated |
| :--- | :--- | :--- |
| `simu/organize.R` | Section 4 (Main Simulation) | Figures 3 |
| `realdata/PM2.5/RAVAS_for_PM.5.R` | Section 5.1 (PM2.5 Analysis) | Figures 4, 5 |
| `realdata/P2P/organize.R` | Section 5.2 (P2P Analysis) | Table 1 |
| `simu/organize_AVAS.R` | Supplement | Figure 3 |

