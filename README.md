# Code and synthetic data for the paper "Case Fatality Risk of Norovirus in England During a Period of Strain Replacement, 2022/23 - 2024/25 Seasons"

This repository contains synthetic data and the code for norovirus case fatality models from the paper:

"Case Fatality Risk of Norovirus in England During a Period of Strain Replacement, 2022/23 - 2024/25 Seasons"
Authors: Maria L. Tang, Amy Douglas, Cristina Celma, Roberto Vivancos, Gauri Godbole, Thomas Ward, Jonathon Mellor

We include synthetic data and code for the main results for cases and deaths, but do not include hospitalisations as they were not modelled.

Synthetic data for cases and deaths was generated using the [`synthpop`](https://cran.r-project.org/package=synthpop) R package (Nowok et al., 2016). Synthetic data were created to support the open-source running of the code and demonstration of the method. The synthetic data approximates key joint distributions in the source data, but does not exactly reproduce results obtained from the original data.


## /

- `depends.R` - package dependencies

## /data

Contains synthetic data for MOLIS-deaths, MOLIS-SGSS-deaths and SGSS-deaths, as referred to in the paper. Columns in the data include:

The individual-level case data requires the following columns:
- `uid` - unique identifier per patient (in the synthetic data, each patient has exactly one row)
- `specimen_date` - date specimen was taken for norovirus test
- `age` - age of patient at specimen date
- `region` - English region of laboratory where MOLIS test result was processed
- `sgss_lab_region` - English region of laboratory where SGSS test result was processed
- `sgss_lab_broad_region` - English broad region ("North", "Midlands", "South") of laboratory where SGSS result was processed
- `requesting_organisation` - healthcare level where the test was requested ("Primary care", "Secondary care", "Other", "Unknown")
- `genotype` - genotype of norovirus specimen from MOLIS test result
- `death_date` - date of death of patient, NA if no linked death

## /scripts

Contains scripts to run case fatality models and generate descriptive results.

Death threshold can be amended through the parameter `death_threshold` in the model run scripts to run the models using different linkage thresholds.

## /outputs

Descriptive outputs and model outputs are saved in here. Folders are created from models with different death thresholds.

# References

Nowok, B., Raab, G. M., & Dibben, C. (2016).
*Synthpop: Bespoke creation of synthetic data in R.*
Journal of Statistical Software, 74(11), 1–26.
https://doi.org/10.18637/jss.v074.i11
