# Phosphoproteomics Data Normalization and Analysis Pipeline

## Overview

This R script provides a comprehensive pipeline for processing, normalizing, and analyzing phosphoproteomics data, with the goal of identifying differentially expressed phosphosites between two experimental conditions (e.g., Light vs. Dark). The script includes steps for data loading, cleaning, normalization, imputation, linear modeling, and visualization.

## Author
Indrani Bera, January 9, 2025

## Dependencies
The script relies on the following R packages:

*   `tidyverse`
*   `preprocessCore`
*   `data.table`
*   `limma`
*   `ggplot2`
*   `ggrepel`
*   `pheatmap`
*   `FactoMineR`
*   `factoextra`
*   `ggbiplot`
*   `EnhancedVolcano`

To install missing packages, use the following commands:

install.packages(c("tidyverse", "preprocessCore", "data.table", "ggplot2", "ggrepel", "pheatmap","FactoMineR","factoextra"))
BiocManager::install(c("limma","ggbiplot","EnhancedVolcano"))


## Input Files

The script expects two tab-separated values (TSV) files as input:

1.  `report.pg_matrix(10).tsv`: Proteomics data file.
2.  `report.phosphosites_90.tsv`: Phosphoproteomics data file.

These files should contain quantitative data from mass spectrometry experiments.  Specifically, the script expects the following columns:

*   **Proteomics Data**: `Protein.Group`, `Protein.Names`, `Genes`, and numerical columns representing protein intensities across different samples.
*   **Phosphoproteomics Data**: `Protein`, `Protein.Names`, `Residue`, `Site`, `Sequence`, and numerical columns representing phosphosite intensities across different samples.

**Important**: Ensure that the input files are correctly formatted and that column names match the script's expectations. Numerical columns should represent the intensity values for each sample.

## Output Files

The script generates several output files at each stage of the process:

*   `proteome_initial.csv`: Initial proteomics data as CSV.
*   `phosphoproteome_initial.csv`: Initial phosphoproteomics data as CSV.
*   `phosphoproteome_selected.csv`: Phosphoproteomics data with selected columns.
*   `proteome_selected.csv`: Proteomics data with selected columns.
*   `phosphoproteome_bio.csv`: Phosphoproteomics data with technical replicates averaged.
*   `proteome_bio.csv`: Proteomics data with technical replicates averaged.
*   `phosphoproteome_log2.csv`: Log2-transformed phosphoproteomics data.
*   `proteome_log2.csv`: Log2-transformed proteomics data.
*   `phosphoproteome_imputed.csv`: Imputed phosphoproteomics data.
*   `proteome_imputed.csv`: Imputed proteomics data.
*   `phosphoproteome_final.csv`: Final phosphoproteomics data with categorical and numerical columns.
*   `proteome_final.csv`: Final proteomics data with categorical and numerical columns.
*   `filtered_combined_phospho_data.csv`: phosphoproteome data after normalization using linear regression.
*   `filtered_phospho_data_norm.csv`: Combined, filtered, and normalized phosphoproteomics data.
*   `differential_expression_results_light_vs_dark.csv`: Results of differential expression analysis.
*   `volcano_plot_light_vs_dark2.jpeg`: Volcano plot visualizing differential expression.
*   `heatmap_numeric_data.jpeg`: Heatmap of the numeric data.

## Script Steps

1.  **Load Data**: Reads proteomics and phosphoproteomics data from TSV files.
2.  **Data Selection**: Selects relevant columns from both datasets.
3.  **Data Type Assignment**: Assigns appropriate data types to columns (factors for categorical variables and numeric for quantitative data).
4.  **Missing Value Check**: Calculates the percentage of missing values in numerical columns and halts execution if the percentage exceeds 50% in either dataset.
5.  **Technical Replicate Processing**: Averages technical replicates, if present, to generate biological means. Requires an even number of numerical columns.
6.  **Log2 Transformation**: Applies a log2 transformation to the numerical columns.
7.  **Missing Value Imputation**: Imputes missing values (NA, NaN, -Inf) with a small value (0.05).
8.  **Quantile Normalization:** Performs quantile normalization.
9.  **Filtering**: filters the data to retain rows where the "Protein.Names" column exists in both files.
10. **Normalization**: Normalizes the phosphoproteomics data based on the whole proteome data using linear regression.
11. **Differential Expression Analysis**: Performs differential expression analysis using the `limma` package.
12. **Volcano Plot**: Generates a volcano plot to visualize differentially expressed phosphosites.
13. **Heatmap**: Generates a heatmap of the numeric data.
14. **PCA analysis**: Performs PCA analysis on numeric data for dimensionality reduction and visualization.

## Usage

1.  Ensure that all required R packages are installed.
2.  Place the R script and the input data files in the same directory.
3.  Modify the input file names in the script if necessary.
4.  Run the script in an R environment.

