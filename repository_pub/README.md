# Title of Dataset: Circadian rhythm of body temperature in fish: Behavioural thermoregulation linked to daily photocycles in the zebrafish (Danio rerio) and the Nile tilapia (Oreochromis nilo)

This dataset contains the data files and R scripts required to run all models presented in the article, as well as producing the figures and summarised output.

## Description of the Data and file structure

The R files require the accompanying data files in the same directory, and can be run in numerical order:

-   *1_model_prefTemp.R*: Loads raw data, calculates observed preferred temperatures, and fits a Bayesian Cosinor model with hourly preferred temperature as the response variable and treatment period (Acclimation, Experiment) as predictor.
-   *2_model_nFishPerChamber.R*: Loads raw data and fits a Bayesian Cosinor model with hourly fish per chamber as the response variable and treatment period (Control CTE, Acclimation, Experiment) as predictor.
-   *3_model_timeSmooth.R*: Loads raw data, calculates observed preferred temperatures, and fits a Bayesian Cosinor model with hourly preferred temperature as the response variable and elapsed time as predictor with smoothed effects.
-   *4_output.R*: Loads raw data and fitted models, summarises outputs, and generates figures.

The data files are:

-   *Tilapia_RawData.csv*: Raw hourly counts per chamber for Nile tilapia; -1 = NA
-   *ZF_Data_RawData.csv:* Raw hourly counts per chamber for Zebrafish; -1 = NA
-   *temperature_Tilapia.csv*: Mean daily temperatures per chamber
-   *temperature_ZF.csv*: Mean daily temperatures per chamber

The output files are:

-   *Table_prefTemp_Cosinor_summary.csv*: Summarised output from the model in *1_model_prefTemp.R* giving the posterior mean, median, and 80%, 90%, and 95% credible intervals for each Cosinor parameter during each phase (Acclimation, Experiment)
-   *Table_prefChamber_Cosinor_summary.csv*: Summarised output from the model in *2_model_nFishPerChamber.R* giving the posterior mean and 95% credible interval for each Cosinor parameter for each chamber during each phase (Control CTE, Acclimation, Experiment)
-   *Table_chamberAcrophase.xlsx*: Summarised output from the model in *2_model_nFishPerChamber.R* giving the posterior mean and 95% credible interval for the acrophase of each chamber during each phase (Control CTE, Acclimation, Experiment)

## Sharing/access Information

Links to other publicly accessible locations of the data: None

Was data derived from another source? If yes, list source(s): No
