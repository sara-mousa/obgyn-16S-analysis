# Female Reproductive Tract Microbiome Analysis - 16S rRNA Workflow

## 🔬 Project Overview
This repository contains a reproducible workflow for analyzing publicly available 16S rRNA sequencing data to investigate differences in samples from the female reproductive tract (vaginal swab and endometrial brushing).

The project aims to characterize microbial composition and diversity, identify differentially abundant taxa, and their predicted functions.

## ❓ Research Question
How do microbial composition and diversity differ between two different female reproductive tract sites (Vaginal vs Endometrial brushing)?

## 🎯 Objectives 
 - Characterize differences in microbial composition and diversity between Vaginal vs Endometrial brushing samples
 - Identify differentially abundant ASVs
 - Explore differentially abundant MetaCyc Pathways
 - Generate reproducible visualizations and statistical analyses 


## 🧬 Dataset Information

 **Study accession:** PRJNA1247240
 
 **Data type:** Paired-end 16S rRNA amplicon sequencing (V4 region) 
 
 **Population:** 
  - Females of reproductive age undergoing Assisted Reproductive Therapy 

 **Source:** NCBI Sequence Read Archive (SRA)

 **Associated publication:** 
 Sola-Leyva, A., Pérez-Prieto, I., Canha-Gouveia, A., Salas-Espejo, E., Molina, N. M., Vargas, E., ... & Altmäe, S. (2026). Comprehensive 16S rRNA gene sequencing and meta-transcriptomic analyses of the female reproductive tract microbiota: two molecular profiles with different messages. Human Reproduction Open, 2026(1), hoag001.

---

## ⚙️ Workflow
  1. Download sequencing data
  2. Quality assessment using FastQC and MultiQC
  3. Read filtering and trimming
  4. ASV inference using DADA2
  5. Chimera removal
  6. Taxonomic assignment
  7. Diversity analysis
  8. Differential abundance analysis using ANCOM-BC2
  9. Functional Prediction using PICRUSt2 
  10.Differential MetaCyc Pathways analysis
  11. Data visualization

---

## 📁 Repository Structure 


---

## 🛠️ Tools Required

### Software
- Bash 
- R
- SRA Toolkit (v3.4.1)
- FastQC (v0.12.1)
- MultiQC (v1.35)
- PICRUSt2 (2.4.1)

### R Packages
- All required R packages are listed at the beginning of each R script.

---

## 📊 Results 
*Analysis in progress*

Final results, tables, and visualizations will be organized in the 3-results/ directory.
---

## 🔁 Reproducibility 
Clone the repository and run the scripts in numerical order.




