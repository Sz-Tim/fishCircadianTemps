# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk
# 1. Preferred temperature: run model



dir.create("out")
# switches ----------------------------------------------------------------

species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")[1]
iter <- 4000
warmup <- 2000
chains <- 4
stan_args <- list(adapt_delta=0.95, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
source("code/00_fn.R")

data.df <- read_csv(glue("data/{species}_RawData.csv")) %>%
  mutate(FishCount=na_if(FishCount, -1),
         Chamber=factor(Chamber),
         Tank=factor(Tank)) %>%
  left_join(read_csv(glue("data/temperature_{species}.csv")) %>%
              filter(Group=="Experiment") %>%
              pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
              mutate(Chamber=factor(str_sub(Chamber, -1, -1)))) %>%
  group_by(ZT, Group, Tank, Day) %>%
  summarise(prefTemp=sum(FishCount*Temp)/(sum(FishCount))) %>%
  ungroup %>%
  mutate(ElapsedTime=ZT+24*(Day-1))
time_sc <- scale(data.df$ElapsedTime)
data.df <- data.df %>%
  mutate(ElapsedTime_sc=c(time_sc)) %>% 
  filter(complete.cases(.))
write_csv(data.df, glue("out/prefTemp_data_{species}.csv"))
saveRDS(time_sc, glue("out/prefTemp_timeScale_{species}.rds"))



# set priors --------------------------------------------------------------

prior <- c(prior(normal(0, 1), class="b", nlpar="A"),
           prior(normal(0, 1), class="b", nlpar="M"),
           prior(von_mises(0, 1), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
           prior(normal(0, 0.1), class="sigma", lb=0),
           prior(normal(0, 3), class="sds", nlpar="A", lb=0),
           prior(normal(0, 3), class="sds", nlpar="M", lb=0),
           prior(normal(0, 3), class="sds", nlpar="phi", lb=0),
           prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="A"),
           prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="M"),
           prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="phi"),
           prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="A"),
           prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="M"),
           prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="phi"))

if(species=="ZF") {
  prior <- c(prior, 
             prior(normal(27.5, 2), class="b", coef="Intercept", nlpar="M"))
} else if(species=="Tilapia") {
  prior <- c(prior, 
             prior(normal(30, 2), class="b", coef="Intercept", nlpar="M"))
}



# run model ---------------------------------------------------------------

out <- brm(bf(prefTemp ~ M + exp(A) * cos(3.141593*(ZT)/12 + phi),
              M ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
              A ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
              phi ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
              nl=TRUE),
           prior=prior, data=data.df, 
           control=stan_args, iter=iter, warmup=warmup, init=0,
           cores=chains, chains=chains, refresh=100,
           save_model=glue("out/prefTemp_mod_{species}.stan"),
           file=glue("out/prefTemp_mod_{species}"))
