# Chromatin Accessibility and PD-1 Response in Gastric Cancer

This project replicates and improves upon the study **"Chromatin accessibility of circulating CD8+ T cells predicts treatment response to PD-1 blockade in patients with gastric cancer" (Shin et al., 2021)**. The goal is to identify genomic regulatory regions linked to immunotherapy response using ATAC-seq data.

## 🔍 Overview

We analyzed publicly available ATAC-seq data from circulating CD8+ T cells in gastric cancer patients undergoing anti-PD-1 therapy. The original paper identified chromatin regions that distinguish responders from non-responders. We extended this analysis using:

- **Normalization comparisons** (TMM, DESeq2, quantile, and control-peak based)
- **Differential peak analysis**
- **Immune-relevant peak filtering**
- **Peak annotation and visualization**
- **Custom scoring systems**
- **Validation using independent patient data**

## 📁 Project Structure

- `424_Final_Report.Rmd` — Full analysis pipeline (R Markdown)
- `424_Final_Report.pdf` — Final formatted report
- `gastric-atac-analysis.Rproj` — RStudio project
- `.gitignore` — Tracks non-versioned files

## 🧪 Methods Used

- **ATAC-seq pre-processing**
- **DESeq2 / edgeR / quantile normalization**
- **LASSO logistic regression**
- **Immune gene enrichment**
- **PCA, volcano plots, and ROC curves**
- **Snakemake-based pipeline (not included in repo)**

## 📊 Technologies

- R / RStudio
- Tidyverse, edgeR, DESeq2, ChIPseeker, ComplexHeatmap, glmnet
- Snakemake (external pipeline)
- Git & GitHub

## 👩‍🔬 Authors

- Ellie McGregor, Mekail Khattak, Caden Roberts  
*UBC Biomedical Engineering, BMEG 424/524 Final Project*


