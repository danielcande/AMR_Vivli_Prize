# Beyond One-Size-Fits-All: Clustered Resistance Data for Surveillance and Personalised Empiric Therapy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Interactive_Dashboard-brightgreen.svg)](https://shiny.rstudio.com/)

## 🏆 2025 AMR Surveillance Data Challenge Submission

**Team:** Daniel Cande, Gwen Knight, Esther van Kleef, Simon Procter  
**Institution:** London School of Hygiene & Tropical Medicine, University of Oxford  
**Live Dashboard:** [GUARDIAN - Global AMR Dashboard](http://guardian.danielcande.co.uk/AMR_Vivli_Prize/)

## 📋 Overview

This repository contains the code and analysis for our innovative approach to antimicrobial resistance (AMR) surveillance using **latent class analysis (LCA)** to identify clusters of co-occurring resistances in ESKAPE pathogens. We move beyond traditional single antibiotic-pathogen monitoring to reveal hidden multidrug resistance patterns that can inform empiric prescribing decisions.

### 🎯 Key Innovation
- **Cluster-based AMR surveillance** using unsupervised machine learning
- **Clinical decision support** for empiric antibiotic prescribing
- **Forecasting models** for resistance burden prediction through 2028
- **Interactive dashboard** for real-time resistance pattern exploration

## 🧬 Scientific Impact

**Published Results:** Analysis of **352,197 isolates** from **82 countries** (2004-2023) across six ESKAPE pathogens:
- *Enterococcus faecium*
- *Staphylococcus aureus* 
- *Klebsiella pneumoniae*
- ***Acinetobacter baumannii*** (featured exemplar)
- *Pseudomonas aeruginosa*
- *Enterobacter* spp.

**Key Findings:**
- Identified **4 distinct resistance clusters** for *A. baumannii*
- XDR-CRAB phenotype affects **54.5%** of isolates globally
- Regional variation: 90% XDR-CRAB prevalence predicted in South Asia by 2028
- ICU patients have **7.53x higher odds** of XDR-CRAB infection


## 🚀 Quick Start

### Prerequisites
- R version 4.3.3 or higher
- Required R packages (see `requirements.R`)

### Installation

```r
# Clone the repository
git clone https://github.com/danielcande/AMR_Vivli_Prize.git
cd AMR_Vivli_Prize

# Install required packages
source("requirements.R")
```

### Running the Analysis

```r
# 1. Data preprocessing (requires ATLAS dataset)
source("src/01_data_preprocessing.R")

# 2. Latent class analysis
source("src/02_latent_class_analysis.R")

# 3. Epidemiological analysis  
source("src/03_epidemiological_analysis.R")

# 4. Forecasting models
source("src/04_forecasting_models.R")

# 5. Clinical prediction tool
source("src/05_clinical_prediction.R")
```

### Launch Interactive Dashboard

```r
# Local deployment
shiny::runApp("dashboard/")

# Or visit the live version
browseURL("http://guardian.danielcande.co.uk/AMR_Vivli_Prize/")
```

## 📊 Key Methods

### 1. Latent Class Analysis (LCA)
- **Package:** `poLCA`
- **Model selection:** AIC, BIC, Entropy, Clinical validity
- **Validation:** Bootstrap resampling (n=1,000)
- **Missing data:** Full Information Maximum Likelihood (FIML)

### 2. Epidemiological Analysis
- **Method:** Multinomial logistic regression
- **Covariates:** Geography, age, sex, healthcare setting, year
- **Reference:** Pan-susceptible cluster

### 3. Forecasting
- **Algorithm:** XGBoost
- **Horizon:** 2023-2028
- **Validation:** Bootstrap cross-validation (n=1,000)
- **Metrics:** RMSE, MAE

### 4. Clinical Risk Prediction
- **Variables:** Age, sex, location (point-of-care available)
- **Output:** Probability of resistance cluster membership
- **Uncertainty:** 95% bootstrap confidence intervals

## 📈 Dashboard Features

### 🌍 GUARDIAN Dashboard Modules:

1. **Resistance Profiles:** Explore antibiogram clusters by pathogen
2. **Global Overview:** Interactive world maps of resistance patterns  
3. **Phenotype Trends:** Temporal resistance evolution
4. **Forecast Tool:** Regional burden predictions through 2028
5. **Patient Risk Predictor:** Bedside clinical decision support

**Access:** [guardian.danielcande.co.uk/AMR_Vivli_Prize](http://guardian.danielcande.co.uk/AMR_Vivli_Prize/)

## 📋 Data Access

This analysis uses data from the **Vivli AMR Register** (ATLAS dataset):
- **Data Request ID:** 00011468
- **Coverage:** 352,197 isolates, 82 countries, 2004-2023
- **Access:** [Vivli AMR Register](https://amr.vivli.org/)

⚠️ **Note:** Raw data is not included in this repository due to data use agreements. Researchers can request access through Vivli.

## 🎯 Clinical Applications

### Empiric Prescribing Support
- **Risk stratification** by patient characteristics
- **Regional adaptation** of treatment guidelines  
- **Real-time resistance** probability assessment

### Public Health Surveillance
- **Cluster-based monitoring** beyond single drug-bug combinations
- **Early detection** of emerging resistance patterns
- **Targeted interventions** based on resistance phenotypes

## 📚 Key Dependencies

```r
# Core analysis packages
library(poLCA)          # Latent class analysis
library(nnet)           # Multinomial regression  
library(xgboost)        # Forecasting models
library(boot)           # Bootstrap validation

# Data manipulation
library(tidyverse)      # Data wrangling
library(data.table)     # Large dataset handling

# Visualization
library(ggplot2)        # Static plots
library(plotly)         # Interactive plots
library(leaflet)        # Interactive maps

# Dashboard
library(shiny)          # Web application framework
library(shinydashboard) # Dashboard UI
library(DT)             # Interactive tables
```

## 🏆 Competition Results

**2025 AMR Vivli Prize Submission**
- **Challenge:** Innovative re-use of AMR surveillance data
- **Funding:** GARDP, Paratek, Pfizer, Vivli
- **Submission Date:** August 7, 2025

## 👥 Team

- **Daniel Cande** - Lead Developer, Statistical Analysis (LSHTM)
- **Dr. Gwen Knight** - Senior Researcher (LSHTM) 
- **Dr. Esther van Kleef** - Senior Researcher (University of Oxford)
- **Simon Procter** - Senior Researcher (LSHTM)

## 📄 Citation

If you use this work, please cite:

```
Cande, D., Knight, G., van Kleef, E., & Procter, S. (2025). 
Beyond One-Size-Fits-All: Clustered Resistance Data for Surveillance 
and Personalised Empiric Therapy. 2025 AMR Surveillance Data Challenge.
GitHub: https://github.com/danielcande/AMR_Vivli_Prize
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## 📧 Contact

- **Daniel Cande:** lsh2400825@student.lshtm.ac.uk
- **Project Website:** [guardian.danielcande.co.uk](http://guardian.danielcande.co.uk/AMR_Vivli_Prize/)

## 🙏 Acknowledgments

- **Vivli** for providing access to the AMR Register
- **ATLAS Program** contributors for surveillance data
- **Competition sponsors:** GARDP, Paratek, Pfizer, Vivli
- **London School of Hygiene & Tropical Medicine**
- **University of Oxford**

---

<div align="center">

**🔬 Transforming AMR Surveillance Through Data Science 🔬**

*Cluster-based resistance monitoring for the genomic era*

[![Dashboard](https://img.shields.io/badge/🌐_Live_Dashboard-Visit_GUARDIAN-success)](http://guardian.danielcande.co.uk/AMR_Vivli_Prize/)

</div>
