# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk
# 2. Preferred temperature: process output



# switches ----------------------------------------------------------------

species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")[2]



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
source("code/00_fn.R")

data.df <- read_csv(glue("out/prefTemp_data_{species}.csv"))
time_sc <- readRDS(glue("out/prefTemp_timeScale_{species}.rds"))

out <- readRDS(glue("out/prefTemp_mod_{species}.rds"))



# initialise dataframes ---------------------------------------------------

global.df <- expand_grid(Day=1:13, ZT=0:23) %>%
  mutate(ElapsedTime=ZT+24*(Day-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24) 
tank.df <- expand_grid(ZT=0:23, Day=1:13, Tank=1:3) %>%
  mutate(ElapsedTime=ZT+24*(Day-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24) 



# extract smoothers -------------------------------------------------------

# Global
s_global <- list(
  M=posterior_smooths(out, newdata=global.df, smooth="s(ElapsedTime_sc,m=2)", nlpar="M") +
    c(as_draws_matrix(out, variable="b_M_Intercept")),
  A=posterior_smooths(out, newdata=global.df, smooth="s(ElapsedTime_sc,m=2)", nlpar="A") +
    c(as_draws_matrix(out, variable="b_A_Intercept")),
  phi=posterior_smooths(out, newdata=global.df, smooth="s(ElapsedTime_sc,m=2)", nlpar="phi") +
    c(as_draws_matrix(out, variable="b_phi_Intercept"))
)
s_global$amplitude <- exp(s_global$A)
s_global$acrophase <- phi_to_ZT(s_global$phi) 
s_global$prefTemp <- calc_prefTemp(s_global$M, s_global$A, s_global$phi, 
                                   matrix(global.df$ZT, nrow=nrow(s_global$M), ncol=nrow(global.df), byrow=T))
# Tank
s_tank  <- list(
  prefTemp=posterior_epred(out, newdata=tank.df, re_formula=NULL),
  M=posterior_epred(out, newdata=tank.df, re_formula=NULL, nlpar="M"),
  A=posterior_epred(out, newdata=tank.df, re_formula=NULL, nlpar="A"),
  phi=posterior_epred(out, newdata=tank.df, re_formula=NULL, nlpar="phi")
)
s_tank$amplitude <- exp(s_tank$A)
s_tank$acrophase <- phi_to_ZT(s_tank$phi)



# summarise posteriors ----------------------------------------------------

# Global
hdi_global <- list(prefTemp=HDInterval::hdi(s_global$prefTemp),
                   M=HDInterval::hdi(s_global$M),
                   amplitude=HDInterval::hdi(s_global$amplitude),
                   acrophase=HDInterval::hdi(s_global$acrophase))
global.df <- global.df %>%
  mutate(prefTemp=colMeans(s_global$prefTemp),
         prefTemp_lo=hdi_global$prefTemp[1,],
         prefTemp_hi=hdi_global$prefTemp[2,],
         M=colMeans(s_global$M),
         M_lo=hdi_global$M[1,],
         M_hi=hdi_global$M[2,],
         amplitude=colMeans(s_global$amplitude),
         amplitude_lo=hdi_global$amplitude[1,],
         amplitude_hi=hdi_global$amplitude[2,],
         acrophase=colMeans(s_global$acrophase),
         acrophase_lo=hdi_global$acrophase[1,],
         acrophase_hi=hdi_global$acrophase[2,])
# Tank
hdi_tank <- list(prefTemp=HDInterval::hdi(s_tank$prefTemp),
                  M=HDInterval::hdi(s_tank$M),
                  amplitude=HDInterval::hdi(s_tank$amplitude),
                  acrophase=HDInterval::hdi(s_tank$acrophase))
tank.df <- tank.df %>%
  mutate(prefTemp=colMeans(s_tank$prefTemp),
         prefTemp_lo=hdi_tank$prefTemp[1,],
         prefTemp_hi=hdi_tank$prefTemp[2,],
         M=colMeans(s_tank$M),
         M_lo=hdi_tank$M[1,],
         M_hi=hdi_tank$M[2,],
         amplitude=colMeans(s_tank$amplitude),
         amplitude_lo=hdi_tank$amplitude[1,],
         amplitude_hi=hdi_tank$amplitude[2,],
         acrophase=colMeans(s_tank$acrophase),
         acrophase_lo=hdi_tank$acrophase[1,],
         acrophase_hi=hdi_tank$acrophase[2,])



# store output ------------------------------------------------------------

global.df %>% 
  mutate(Species=species) %>%
  write_csv(glue("out/prefTemp_predGlobal_{species}.csv"))
tank.df %>% 
  mutate(Species=species) %>%
  write_csv(glue("out/prefTemp_predTank_{species}.csv"))
saveRDS(s_global, glue("out/prefTemp_posteriorGlobal_{species}.rds"))
saveRDS(s_tank, glue("out/prefTemp_posteriorTank_{species}.rds"))
