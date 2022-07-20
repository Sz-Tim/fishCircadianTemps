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
                   "(1|Tank)",
                   paste0("(1+", paste0(mod.terms, collapse="*"), "|Tank)"))

priors <- c(prior(normal(0, 2), class="b"),
            prior(normal(0, 2), class="Intercept"),
            prior(cauchy(0, 2), class="sd"))
prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A", lb=0),
              prior(normal(1, 1), class="b", nlpar="M", lb=0),
              prior(uniform(0, 24), class="b", nlpar="phi", lb=0, ub=24),
              prior(normal(0, 0.1), class="sd", nlpar="A", lb=0),
              prior(normal(0, 0.1), class="sd", nlpar="M", lb=0),
              prior(normal(0, 0.1), class="sd", nlpar="phi", lb=0),
              prior(normal(0, 0.1), dpar="sigma", lb=0))


# fit model ---------------------------------------------------------------

out.nl <- brm(bf(ln_FishCount ~ M + A * cos(3.141593*(ZT + phi)/12),
                 M ~ 0 + Group*Chamber + (0+Group*Chamber|Tank), 
                 A ~ 0 + Group*Chamber + (0+Group*Chamber|Tank), 
                 phi ~ 0 + Group*Chamber + (0+Group*Chamber|Tank),
                 sigma ~ Group,
                 nl=TRUE),
              prior=prior.nl, 
              control=stan_args,
              iter=iter, warmup=warmup, init=0,
              data=data.noNA, cores=chains, chains=chains, refresh=50,
              save_model=glue("models/nl/mod_count_{species}.stan"),
              file=glue("models/nl/out_count_{species}"))

# out <- brm(bf(paste("ln_FishCount ~", 
#                     paste0(mod.terms, collapse="*"),
#                     "+", mod.rand), 
#               sigma ~ Group),
#            prior=priors, 
#            control=stan_args,
#            iter=iter, warmup=warmup, init=0,
#            data=data.noNA, cores=chains, refresh=50,
#            save_model=glue("models/mod_count_{mod_type}_{species}.stan"),
#            file=glue("models/out_count_{mod_type}_{species}"))
