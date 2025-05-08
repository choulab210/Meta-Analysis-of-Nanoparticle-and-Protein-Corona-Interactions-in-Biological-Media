# PC-DB-2025
# -*- coding: utf-8 -*-


# Meta-Analysis of Nanoparticle and Protein Corona Interactions in Biological Media
This repository contains all raw data and Python code files for the project entitled  
**"Meta-Analysis of Nanoparticle and Protein Corona Interactions in Biological Media"**

## Respository Structure
Within the **`Code`** folder there are 2 Python scripts: `Prediction_Models` and `Uniprot_Extraction`. Their respective datasets can be found in the **`Datasets`** subfolder.

  ### 1. `Prediction Models`
  This file contains the compiled Python script for development and evaluation of five machine learning classification models
  

 - K-Nearest Neighbors (k-NN)
 - Light Gradient Boosting Machine (LightGBM)
 - Random Forest (RF)
 - Support Vector Machine (SVM)
 - eXtreme Gradient Boosting (XGBoost) 

These models are trained and evaluated using the curated dataset available at [wb.ai2tox.com](https://wb.ai2tox.com). Performance metrics including *accuracy*,  *precision*, *recall*, *F1-Score*, and *ROC-AUC* are reported for comparative analysis. 

  ### 2. `Uniprot Extraction`
This subfolder includes scripts for retrieving protein-specific physicochemical properties from the UniProt database, including Molecular Weight (MW) and Amino Acid Sequences. The [BioPython](https://biopython.org/) (Bio.SeqUtils) package was used to calculate isoelectric point (pI) and Grand Average of Hydropathy (GRAVY) for protein analyses detailed in the manuscript

## Contributors

- **Author:** Alexa Canchola, University of California, Riverside
- **Advisor:** Wei-Chun Chou, University of California, Riverside  


*The manuscript describing the files in this repository is currently under review.*
