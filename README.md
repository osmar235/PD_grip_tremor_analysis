# PD Grip Tremor Analysis

This repository contains the R analysis code associated with the manuscript:

**Influence of Visual Feedback, Hand Dominance, and Disease Severity on Grip Force Control and Tremor Dynamics in Parkinson's Disease: A Cross-Sectional Secondary Analysis**

## Overview

The scripts in this repository reproduce the statistical analyses used to examine grip force control, tremor dynamics, visual feedback, hand dominance, disease severity, age, sex, and frequency-domain wavelet outcomes in individuals with Parkinson's disease.

The repository includes code for:

1. Linear mixed-effects models for time-domain grip force and tremor outcomes.
2. Linear mixed-effects models for frequency-domain and wavelet-based outcomes.
3. Additional reviewer-requested analyses, including UPDRS-III correlations, sex-stratified summaries, confidence intervals, model diagnostics, and sensitivity analyses.

## Repository structure

```text
PD_grip_tremor_analysis/
│
├── README.md
├── LICENSE
│
├── scripts/
│   ├── 01_MLM_3factor_time_domain_outcomes.R
│   ├── 02_MLM_4factor_frequency_domain_outcomes.R
│   └── 03_reviewer_additional_analyses.R
│
├── data/
│   └── README_data_access.txt
│
└── results/
    └── README_results.txt
```

## Data availability

The participant-level datasets generated and analyzed during the current study are not publicly available due to participant consent and privacy restrictions. De-identified data may be made available from the corresponding author upon reasonable request and subject to applicable ethical, institutional, and data-sharing restrictions.

Because the datasets are not publicly shared, the code is provided for transparency and reproducibility by qualified researchers who obtain access to the data.

## Required data files

After obtaining the data, users should place the required files in the `data/` folder using the filenames expected by the scripts.

The main scripts expect files corresponding to:

* Time-domain/no-frequency outcomes.
* Frequency-domain/wavelet outcomes.
* UPDRS-III and reviewer-requested supplementary analyses.

## Software

The analyses were conducted in R using packages including:

* `lme4`
* `lmerTest`
* `emmeans`
* `dplyr`
* `tidyverse`
* `ggplot2`
* `broom.mixed`
* `readxl`
* `car`
* `ggeffects`
* `sjPlot`
* `boot`

## Scripts

### `01_MLM_3factor_time_domain_outcomes.R`

Runs mixed-effects models for time-domain outcomes, including grip force and tremor measures without the frequency-band factor.

### `02_MLM_4factor_frequency_domain_outcomes.R`

Runs mixed-effects models for frequency-domain and wavelet-based outcomes, including the frequency-band factor.

### `03_reviewer_additional_analyses.R`

Runs additional analyses requested during peer review, including UPDRS-III correlations, sex-stratified summaries, confidence intervals, model diagnostics, and sensitivity analyses.

## Citation

If you use this code, please cite the associated article and this repository DOI.
