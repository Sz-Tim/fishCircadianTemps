# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[1]
mod_type <- c("RI", "RS_ZTonly", "RS_noInteractions", "RS_all")[2]
iter <- 3000
warmup <- 2000
chains <- 4
stan_args <- list(adapt_delta=0.99, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)

data.df <- dir("data", glue("{species}_Data_matrix"), full.names=T) %>%
  readxl::read_xlsx(., 1, col_types=c(rep("numeric", 6), "skip", "skip")) %>%
  mutate(ln_FishCount=log(FishCount+1),
         cos_ZT=cos(2*pi*ZT/24),
         sin_ZT=sin(2*pi*ZT/24),
         Chamber=factor(Chamber),
         Group=factor(Group, 
                      levels=c(3, 1, 2), 
                      labels=c("Control", "Acclimation", "Experiment")),
         Tank=factor(Tank))



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.))

form.fixed <- bf(ln_FishCount ~ I(cos(6.283185*ZT/24))*Chamber*Group + 
                   I(sin(6.283185*ZT/24))*Chamber*Group)
mod.form <- switch(
  mod_type,
  RI=update(form.fixed,
            ~. + (1|Tank)),
  RS_ZTonly=update(form.fixed, 
                   ~. + (1 + I(cos(6.283185*ZT/24)) +
                           I(sin(6.283185*ZT/24))|Tank)),
  RS_noInteractions=update(form.fixed, 
                           ~. + (1 + I(cos(6.283185*ZT/24)) +
                                   I(sin(6.283185*ZT/24)) +
                                   Chamber + Group|Tank)),
  RS_all=update(form.fixed, 
                ~. + (1 + I(cos(6.283185*ZT/24))*Chamber*Group + 
                        I(sin(6.283185*ZT/24))*Chamber*Group|Tank))) 

priors <- c(prior(normal(0, 2), class="b"),
            prior(normal(0, 2), class="Intercept"),
            prior(cauchy(0, 2), class="sd"))


# fit model ---------------------------------------------------------------

out <- brm(mod.form,
           prior=priors, 
           control=stan_args,
           iter=iter, warmup=warmup, init=0,
           data=data.noNA, cores=chains, refresh=50,
           save_model=glue("models/mod_{mod_type}.stan"),
           file=glue("models/out_{mod_type}_{species}"))
