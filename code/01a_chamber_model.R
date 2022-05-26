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
         cos_ZT=cos(2*pi*ZT/24),
         sin_ZT=sin(2*pi*ZT/24),
         Chamber=factor(Chamber),
         Group=factor(Group, 
                      levels=c(3, 1, 2), 
                      labels=c("Control", "Acclimation", "Experiment")),
         Tank=factor(Tank)) %>%
  left_join(read_csv(glue("data/temp_{species}.csv")) %>%
              pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
              mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                     Group=factor(Group, levels=c("Control", "Acclimation", "Experiment")))) %>%
  mutate(GrpDay=factor(paste(as.numeric(Group), str_pad(Days, 2, "l", "0"), sep="_")))



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.))

mod.terms <- c("I(cos(6.283185*ZT/24))", 
               "I(sin(6.283185*ZT/24))",
               "Chamber", 
               "Group")
mod.rand <- ifelse(mod_type=="RI", 
                   "(1|Tank) + (1|GrpDay)",
                   paste0("(1+", paste0(mod.terms, collapse="*"), "|Tank)", 
                          " + (1|GrpDay)"))

priors <- c(prior(normal(0, 2), class="b"),
            prior(normal(0, 2), class="Intercept"),
            prior(cauchy(0, 2), class="sd"))


# fit model ---------------------------------------------------------------

out <- brm(bf(paste("ln_FishCount ~", 
                    paste0(mod.terms, collapse="*"),
                    "+", mod.rand), 
              sigma ~ Group),
           prior=priors, 
           control=stan_args,
           iter=iter, warmup=warmup, init=0,
           data=data.noNA, cores=chains, refresh=50,
           save_model=glue("models/mod_{mod_type}.stan"),
           file=glue("models/out_{mod_type}_{species}"))
