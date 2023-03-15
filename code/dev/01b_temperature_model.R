# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[2]
iter <- 4000
warmup <- 3000
chains <- 4
stan_args <- list(adapt_delta=0.99, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)

data.df <- readxl::read_xlsx(dir("data", glue("{species}_.*RawData2"), full.names=T),
                             1, col_types=c(rep("numeric", 6), "skip", "skip")) %>%
  mutate(ln_FishCount=log(FishCount+1),
         Chamber=factor(Chamber),
         Group=factor(Group, 
                      levels=c(3, 1, 2), 
                      labels=c("Control", "Acclimation", "Experiment")),
         Tank=factor(Tank)) %>%
  left_join(read_csv(glue("data/temp_{species}.csv")) %>%
              pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
              mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                     Group=factor(Group, levels=c("Control", "Acclimation", "Experiment")))) %>%
  filter(Group != "Control") %>%
  group_by(ZT, Group, Tank, Days) %>%
  summarise(prefTemp=sum(FishCount*Temp)/(sum(FishCount))) %>%
  ungroup %>%
  mutate(ZT_ang=2*pi*ZT/24)



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.))

if(species=="ZF") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(27.5, 2), class="b", nlpar="M"),
                prior(von_mises(0, 0.001), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                prior(normal(0, 0.1), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.1), dpar="sigma", lb=0))
} else if(species=="Tilapia") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(30, 2), class="b", nlpar="M"),
                prior(von_mises(0, 0.001), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                prior(normal(0, 0.1), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.1), dpar="sigma", lb=0))
}


# fit model ---------------------------------------------------------------

out.nl <- brm(bf(prefTemp ~ M + exp(A) * cos(3.141593*(ZT)/12 + phi),
                 M ~ 0 + Group + (0+Group|Tank), 
                 A ~ 0 + Group + (0+Group|Tank), 
                 phi ~ 0 + Group + (0+Group|Tank),
                 sigma ~ Group,
                 nl=TRUE),
              prior=prior.nl, 
              control=stan_args,
              iter=iter, warmup=warmup, init=0,
              data=data.noNA, cores=chains, chains=chains, refresh=50,
              save_model=glue("models/cosinor/mod_temperature_vm_expA_{species}.stan"),
              file=glue("models/cosinor/out_temperature_vm_expA_{species}"))
