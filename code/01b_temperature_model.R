# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[1]
mod_type <- c("RI", "RS")[2]
iter <- 3000
warmup <- 2000
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
  group_by(ZT, Group, Tank, Days) %>%
  summarise(prefTemp=sum(FishCount*Temp)/(sum(FishCount))) %>%
  ungroup %>%
  mutate(ZT_ang=2*pi*ZT/24)



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.))

mod.terms <- c("I(cos(6.283185*ZT/24))", 
               "I(sin(6.283185*ZT/24))",
               "Group")
mod.rand <- ifelse(mod_type=="RI", 
                   "(1|Tank)",
                   paste0("(1+", paste0(mod.terms, collapse="*"), "|Tank)"))

if(species=="ZF") {
  priors <- c(prior(normal(0,2), "b"),
              prior(normal(27.5, 2), "Intercept"),
              prior(cauchy(0,1), "sd"))
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A", lb=0),
                prior(normal(27.5, 2), class="b", nlpar="M"),
                prior(uniform(0, 24), class="b", nlpar="phi", lb=0, ub=24),
                prior(normal(0, 0.5), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.5), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.5), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.5), dpar="sigma", lb=0))
} else {
  priors <- c(prior(normal(0,2), "b"),
              prior(normal(30, 2), "Intercept"),
              prior(cauchy(0,1), "sd"))
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A", lb=0),
                prior(normal(30, 2), class="b", nlpar="M"),
                prior(uniform(0, 24), class="b", nlpar="phi", lb=0, ub=24),
                prior(normal(0, 0.5), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.5), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.5), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.5), dpar="sigma", lb=0))
}




# fit model ---------------------------------------------------------------

out.nl <- brm(bf(prefTemp ~ M + A * cos(3.141593*(ZT + phi)/12),
                 M ~ 0 + Group + (0+Group|Tank), 
                 A ~ 0 + Group + (0+Group|Tank), 
                 phi ~ 0 + Group + (0+Group|Tank),
                 sigma ~ Group,
                 nl=TRUE),
              prior=prior.nl, 
              control=stan_args,
              iter=iter, warmup=warmup, init=0,
              data=data.noNA, cores=chains, chains=chains, refresh=50,
              save_model=glue("models/nl/mod_temperature_wide_{species}.stan"),
              file=glue("models/nl/out_temperature_wide_{species}"))

# out <- brm(bf(paste("prefTemp ~", 
#                     paste0(mod.terms, collapse="*"),
#                     "+", mod.rand), 
#               sigma ~ Group),
#            prior=priors, 
#            control=stan_args,
#            iter=iter, warmup=warmup, init=0,
#            data=data.noNA, cores=4, refresh=50,
#            save_model=glue("models/RF_cor/mod_temperature_{mod_type}_{species}.stan"),
#            file=glue("models/RF_cor/out_temperature_{mod_type}_{species}"))
