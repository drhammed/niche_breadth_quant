# Niche Breadth Quantification - Modular Codebase

This repository contains a modular codebase for comparing niche breadth quantification methods across different ecological scenarios.

## Quick Start

```r
# 1. Set working directory to whatever folder you're working on
setwd("path/toyourfolder")

# 2. Run simulations (generates all 12 scenarios + empirical validation)
rmarkdown::render("01_simulations.Rmd")

# 3. Run metrics agreement analysis
rmarkdown::render("02_metrics_agreement.Rmd")

# 4. Run rank RMSE analysis (supplementary evaluation)
source("03_rmse_analysis.R")
```

## Project Structure

```
Ch1-refactored/
├── README.md                     # This file
├── config.R                      # All configurable params
│
├── R/                            # Function library
│   ├── helper_functions.R        # Utility functions
│   ├── simulation_functions.R    # Unified Gaussian response functions
│   ├── theta_functions.R         # Co-occurrence theta metrics
│   └── niche_breadth_metrics.R   # All 9 niche breadth metrics
│
├── 01_simulations.Rmd            # Main simulation script
├── 02_metrics_agreement.Rmd      # Metrics comparison analysis
├── 03_rmse_analysis.R            # Rank RMSE supplementary analysis
│
├── data/                         # Input data
│   ├── spp_mod.csv              # BBS species occurrence data
│   └── bioclim/                 # WorldClim data
│
└── results/                      # Output files
    ├── simulations/              # RDS files per scenario
    ├── empirical/                # Empirical validation results
    └── figures/                  # figures
```

## Setup Instructions

### 1. Copy Required Data Files

Ensure you have the species and climate data

### 2. Install Required Packages

```r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  doParallel, foreach, data.table, ade4, tidyverse,
  vegan, dplyr, reshape2, ggplot2, corrplot,
  nicheROVER, hypervolume, mFD, ecoCopula, mgcv,
  terra, sf, sp, geodata, dismo,
  FactoMineR, ComplexHeatmap, MatrixCorrelation,
  ggcorrplot, pheatmap, jsonlite, viridis
)
```

## Configuration

All simulation parameters are centralized in `config.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_SPECIES` | 30 | Number of species per simulation |
| `N_SITES` | 500 | Number of sites (communities) |
| `N_REPS` | 30 | Simulation iterations per scenario |
| `NUM_CORES` | 6 | Parallel processing cores |
| `RANDOM_SEED` | 9999 | Reproducibility seed |

Modify `config.R` to adjust parameters.

## Simulation Scenarios

The analysis covers 12 simulation scenarios:

| # | Response | Breadth Dist | Environments | Name |
|---|----------|--------------|--------------|------|
| 1 | Symmetric | Uniform | 1 | Sym_uniform_env1 |
| 2 | Symmetric | Normal | 1 | Sym_normal_env1 |
| 3 | Symmetric | Gamma | 1 | Sym_gamma_env1 |
| 4 | Asymmetric | Uniform | 1 | Asy_uniform_env1 |
| 5 | Asymmetric | Normal | 1 | Asy_normal_env1 |
| 6 | Asymmetric | Gamma | 1 | Asy_gamma_env1 |
| 7 | Symmetric | Uniform | 2 | Sym_uniform_env2 |
| 8 | Symmetric | Normal | 2 | Sym_normal_env2 |
| 9 | Symmetric | Gamma | 2 | Sym_gamma_env2 |
| 10 | Asymmetric | Uniform | 2 | Asy_uniform_env2 |
| 11 | Asymmetric | Normal | 2 | Asy_normal_env2 |
| 12 | Asymmetric | Gamma | 2 | Asy_gamma_env2 |

## Niche Breadth Metrics

Nine metrics are calculated for each species:

| Metric | Function | Package | Description |
|--------|----------|---------|-------------|
| SimpSSI | `co_occur()` | base | Multi-site Simpson index |
| beta.a | `co_occur()` | base | Additive beta partitioning |
| beta.w | `co_occur()` | base | Whittaker's multiplicative beta |
| om_tol | `omi_params_fun()` | ade4 | OMI tolerance |
| nr_hv | `nr_hypervolume_fun()` | nicheROVER | Bayesian hypervolume |
| hv_blond | `hypervolume_blond_fun()` | hypervolume | Kernel density hypervolume |
| nb_Gam | `estimate_nicheBreadth_Gam()` | mgcv | GAM-based niche breadth |
| nb_latent | `estimate_nicheBreadth_Latents()` | ecoCopula | Latent variable model |
| nb_dist | `estimate_nicheBreadth_avg.Dist()` | vegan | Average Hellinger distance |

## Key Functions

### Unified Simulation Function


```r

generate_community(
  n_species = 30,
  n_sites = 500,
  n_env = 2,
  response_type = "symmetric",
  breadth_distribution = "normal"
)
```

### Calculate All Metrics

```r
# Wrapper function calculates all 9 metrics at once
results <- calculate_all_metrics(
  sim.com = presence_absence_matrix,
  env_vars = environmental_variables
)
```

## Output Files

After running both Rmd files:

```
results/
├── simulations/
│   ├── result_df_Sym_uniform_env1.rds
│   ├── result_df_Sym_normal_env1.rds
│   └── ... (12 files total)
├── empirical/
│   └── result_df_empirical.rds
├── matrices_list.json
├── oracle_correlations_summary.csv
├── rank_rmse_summary.csv
└── spearman_vs_rmse_comparison.csv
```

## Authors

- Hammed Akande
- Pedro Peres-Neto

