# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk
# 3. Chamber composition: run model




setwd("repository_pub")
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
source("00_fn.R")

data.df <- read_csv(glue("{species}_RawData.csv")) %>%
  mutate(FishCount=na_if(FishCount, -1),
         Chamber=factor(Chamber),
         Tank=factor(Tank)) %>%
  mutate(ElapsedTime=ZT+24*(Day-1))
time_sc <- scale(data.df$ElapsedTime)
data.df <- data.df %>%
  mutate(ElapsedTime_sc=c(time_sc)) %>% 
  filter(complete.cases(.)) %>%
  group_by(Group, Tank, ElapsedTime_sc) %>%
  mutate(FishCount=FishCount+0.001,
         propFish=FishCount/sum(FishCount),
         Chamber=paste0("Ch_", Chamber)) %>%
  ungroup %>%
  select(Group, Tank, ElapsedTime_sc, ZT, Chamber, propFish) %>%
  pivot_wider(names_from=Chamber, values_from=propFish) %>%
  mutate(Y=cbind(Ch_1=Ch_1, Ch_2=Ch_2, Ch_3=Ch_3, Ch_4=Ch_4, Ch_5=Ch_5))
saveRDS(data.df, glue("out/chmbrComp_data_{species}.rds"))
saveRDS(time_sc, glue("out/chmbrComp_timeScale_{species}.rds"))




# set priors --------------------------------------------------------------

prior <- c(prior(normal(0, 1), class="b", nlpar="A1"),
           prior(normal(0, 1), class="b", nlpar="A2"),
           prior(normal(0, 1), class="b", nlpar="A4"),
           prior(normal(0, 1), class="b", nlpar="A5"),
           prior(normal(0, 1), class="b", nlpar="M1"),
           prior(normal(0, 1), class="b", nlpar="M2"),
           prior(normal(0, 1), class="b", nlpar="M4"),
           prior(normal(0, 1), class="b", nlpar="M5"),
           prior(von_mises(0, 1), class="b", nlpar="pphi1", lb=-3.141593, ub=3.141593),
           prior(von_mises(0, 1), class="b", nlpar="pphi2", lb=-3.141593, ub=3.141593),
           prior(von_mises(0, 1), class="b", nlpar="pphi4", lb=-3.141593, ub=3.141593),
           prior(von_mises(0, 1), class="b", nlpar="pphi5", lb=-3.141593, ub=3.141593),
           prior(normal(0, 2), class="sds", nlpar="A1", lb=0),
           prior(normal(0, 2), class="sds", nlpar="A2", lb=0),
           prior(normal(0, 2), class="sds", nlpar="A4", lb=0),
           prior(normal(0, 2), class="sds", nlpar="A5", lb=0),
           prior(normal(0, 2), class="sds", nlpar="M1", lb=0),
           prior(normal(0, 2), class="sds", nlpar="M2", lb=0),
           prior(normal(0, 2), class="sds", nlpar="M4", lb=0),
           prior(normal(0, 2), class="sds", nlpar="M5", lb=0),
           prior(normal(0, 2), class="sds", nlpar="pphi1", lb=0),
           prior(normal(0, 2), class="sds", nlpar="pphi2", lb=0),
           prior(normal(0, 2), class="sds", nlpar="pphi4", lb=0),
           prior(normal(0, 2), class="sds", nlpar="pphi5", lb=0))



# run model ---------------------------------------------------------------

out.exp <- brm(bf(Y ~ Ch) +
                 nlf(muCh1 ~ M1 + exp(A1)*cos(3.141593*(ZT)/12 + pphi1)) +
                 lf(M1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 lf(A1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                 lf(pphi1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 nlf(muCh2 ~ M2 + exp(A2)*cos(3.141593*(ZT)/12 + pphi2)) +
                 lf(M2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 lf(A2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                 lf(pphi2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 nlf(muCh4 ~ M4 + exp(A4)*cos(3.141593*(ZT)/12 + pphi4)) +
                 lf(M4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 lf(A4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                 lf(pphi4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 nlf(muCh5 ~ M5 + exp(A5)*cos(3.141593*(ZT)/12 + pphi5)) +
                 lf(M5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                 lf(A5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                 lf(pphi5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)),
               family=dirichlet(refcat="Ch_3"), prior=prior, data=filter(data.df, Group=="Experiment"), 
               control=stan_args, iter=iter, warmup=warmup, init=0,
               cores=chains, chains=chains, refresh=10,
               save_model=glue("out/chmbrComp_mod_exp_{species}.stan"),
               file=glue("out/chmbrComp_mod_exp_{species}"))

out.ctrl <- brm(bf(Y ~ Ch) +
                      nlf(muCh1 ~ M1 + exp(A1)*cos(3.141593*(ZT)/12 + pphi1)) +
                      lf(M1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      lf(A1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                      lf(pphi1 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      nlf(muCh2 ~ M2 + exp(A2)*cos(3.141593*(ZT)/12 + pphi2)) +
                      lf(M2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      lf(A2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                      lf(pphi2 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      nlf(muCh4 ~ M4 + exp(A4)*cos(3.141593*(ZT)/12 + pphi4)) +
                      lf(M4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      lf(A4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                      lf(pphi4 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      nlf(muCh5 ~ M5 + exp(A5)*cos(3.141593*(ZT)/12 + pphi5)) +
                      lf(M5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) +
                      lf(A5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)) + 
                      lf(pphi5 ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1)),
                    family=dirichlet(refcat="Ch_3"), prior=prior, data=filter(data.df, Group=="Control"), 
                    control=stan_args, iter=iter, warmup=warmup, init=0,
                    cores=chains, chains=chains, refresh=10,
                save_model=glue("out/chmbrComp_mod_ctrl_{species}.stan"),
                file=glue("out/chmbrComp_mod_ctrl_{species}"))
