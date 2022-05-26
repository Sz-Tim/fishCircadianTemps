# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[1]
mod_type <- c("RI", "RS")[2]
iter <- 2000
warmup <- 1000
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
  ungroup



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
} else {
  priors <- c(prior(normal(0,2), "b"),
              prior(normal(30, 2), "Intercept"),
              prior(cauchy(0,1), "sd"))
}




# fit model ---------------------------------------------------------------

out <- brm(bf(paste("prefTemp ~", 
                    paste0(mod.terms, collapse="*"),
                    "+", mod.rand), 
              sigma ~ Group),
           prior=priors, 
           control=list(adapt_delta=0.99, max_treedepth=20),
           iter=3000, warmup=2000, init=0,
           data=data.noNA, cores=4, refresh=50,
           save_model=glue("models/mod_temperature_{mod_type}.stan"),
           file=glue("models/out_temperature_{species}"))
