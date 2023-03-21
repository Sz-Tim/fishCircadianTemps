# Circadian rhythm of preferred temperature in fish: Behavioural thermoregulation linked to daily photocycles in zebrafish and Nile tilapia
**Luisa M. Vera, Gonzalo de Alba, Silvere Santos, Tim M. Szewczyk, Simon A. Mackenzie, Francisco J. Sánchez-Vázquez, Sonia Rey-Planellas. _Journal of Thermal Biology_, 2023**


[![DOI](https://zenodo.org/badge/494023087.svg)](https://zenodo.org/badge/latestdoi/494023087)

The final code and datasets associated with this manuscript are available from the journal website and data dryad. This repository hosts all files, in addition to the development history. 


## Description of the data and file structure

The **R files** in `/code` can be run from 1 to 5:

-   *1_preferredTemp_runModel.R*: Loads raw data, calculates observed preferred temperatures, and fits a Bayesian cosinor model with hourly preferred temperature as the response variable.
-   *2_preferredTemp_processOutput.R*: Reshapes, summarises, and stores model output.
-   *3_chamberComp_runModel.R*: Loads raw data, calculates observed preferred temperatures, and fits a Bayesian cosinor model with hourly proportional composition of fish across chambers as the response variable.
-   *4_chamberComp_processOutput.R*: Reshapes, summarises, and stores model output.
-   *5_figuresAndSummaries.R*: Loads raw data and fitted models, summarises outputs, and generates figures.
-   *00_fn.R*: Helper functions for processing output.

The finalised, clean **data files** in `/data` are:

-   *Tilapia_RawData.csv*: Raw hourly counts per chamber for Nile tilapia; -1 = NA
-   *ZF_Data_RawData.csv:* Raw hourly counts per chamber for Zebrafish; -1 = NA
-   *temperature_Tilapia.csv*: Mean daily temperatures per chamber
-   *temperature_ZF.csv*: Mean daily temperatures per chamber

The output files are:

-   *Table_preferredTemp_summary.csv*: Summarised output from the model in *1_preferredTemp_runModel.R* giving the posterior mean and 95% credible intervals (HPDIs) for preferred temperature and each cosinor parameter. 
-   *Table_chamberComp_summary.csv*: Summarised output from the model in *3_chamberComp_runModel.R* giving the posterior mean and 95% credible intervals (HPDIs) for proportion of fish in each chamber and each cosinor parameter. 

