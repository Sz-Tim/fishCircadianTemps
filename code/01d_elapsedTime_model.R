# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[2]
iter <- 2000
warmup <- 1000
chains <- 4
stan_args <- list(adapt_delta=0.9, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
theme_set(theme_classic())

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
  mutate(ElapsedTime=ZT+24*(Days-1)+24*3*(Group=="Experiment"),
         ElapsedTime_sc=c(scale(ElapsedTime))) %>%
  filter(Group!="Control")



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.)) %>% 
  mutate(Days=factor(Days))

if(species=="ZF") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(27.5, 2), class="b", coef="Intercept", nlpar="M"),
                prior(normal(0, 1), class="b", nlpar="M"),
                prior(von_mises(0, 1), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                prior(normal(0, 0.1), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.1), class="sigma", lb=0))
} else if(species=="Tilapia") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(30, 2), class="b", coef="Intercept", nlpar="M"),
                prior(normal(0, 1), class="b", nlpar="M"),
                prior(von_mises(0, 1), nlpar="phi", lb=-3.141593, ub=3.141593),
                prior(von_mises(0, 1), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                prior(normal(0, 0.1), class="sd", nlpar="A", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="M", lb=0),
                prior(normal(0, 0.1), class="sd", nlpar="phi", lb=0),
                prior(normal(0, 0.1), class="sigma", lb=0))
}


# fit model ---------------------------------------------------------------

out.nl <- brm(bf(prefTemp ~ M + exp(A) * cos(3.141593*(ZT)/12 + phi),
                 M ~ 1 + s(ElapsedTime_sc) + (1|Tank),
                 A ~ 1 + s(ElapsedTime_sc) + (1|Tank),
                 phi ~ 1 + s(ElapsedTime_sc) + (1|Tank),
                 nl=TRUE),
              prior=prior.nl,
              control=stan_args,
              iter=iter, warmup=warmup, init=0,
              data=data.noNA, cores=chains, chains=chains, refresh=10,
              save_model=glue("models/cosinor/acclimation_s_elapsedTime_{species}.stan"),
              file=glue("models/cosinor/acclimation_s_elapsedTime_{species}"))

pred.dat <- expand_grid(ZT=0:24, Days=1:13, Tank=1:3) %>%
  mutate(ElapsedTime=ZT+24*(Days-1),
         ElapsedTime_sc=c(scale(ElapsedTime))) 
pred.ls  <- list(
  prefTemp=posterior_epred(out.nl, newdata=pred.dat, re_formula=NA),
  M=posterior_epred(out.nl, newdata=pred.dat, re_formula=NA, nlpar="M"),
  A=posterior_epred(out.nl, newdata=pred.dat, re_formula=NA, nlpar="A"),
  phi=posterior_epred(out.nl, newdata=pred.dat, re_formula=NA, nlpar="phi")
)
pred.ls$A <- exp(pred.ls$A)
pred.ls$phi <- pred.ls$phi + 2*pi*(pred.ls$phi < -pi) - 2*pi*(pred.ls$phi > pi)
pred.ls$phi <- (-pred.ls$phi+pi)*12/pi-12
pred.ls$phi <- pred.ls$phi + 24*(pred.ls$phi < 0)
pred.df <- pred.dat %>%
  mutate(pred=colMeans(pred.ls$prefTemp),
         pred_lo=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.025)),
         pred_hi=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.975)),
         pred_M=colMeans(pred.ls$M),
         pred_M_lo=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.025)),
         pred_M_hi=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.975)),
         pred_A=colMeans(pred.ls$A),
         pred_A_lo=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.025)),
         pred_A_hi=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.975)),
         pred_phi=colMeans(pred.ls$phi),
         pred_phi_lo=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.025)),
         pred_phi_hi=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.975)))

temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))

ggplot(pred.df, aes(ElapsedTime, pred)) +
  geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
             linetype=2, colour="grey30") +
  geom_point(data=data.noNA, aes(y=prefTemp, group=Tank, colour=Tank), size=0.25) + 
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.3, colour=NA) +
  geom_line() +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
  ylab("Preferred temperature (ºC)") + ggtitle(species)
ggsave(paste0("figs/sTime_prefTemp_", species, ".jpg"), width=6, height=3)

ggplot(pred.df %>% filter(Tank==1), aes(ElapsedTime, pred_M)) +
  geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
             linetype=2, colour="grey30") +
  geom_ribbon(aes(ymin=pred_M_lo, ymax=pred_M_hi), alpha=0.3, colour=NA) +
  geom_line() +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
  ylab("MESOR (ºC)") + ggtitle(species)
ggsave(paste0("figs/sTime_M_", species, ".jpg"), width=6, height=3)

ggplot(pred.df %>% filter(Tank==1), aes(ElapsedTime, pred_A)) +
  geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
             linetype=2, colour="grey30") +
  geom_ribbon(aes(ymin=pred_A_lo, ymax=pred_A_hi), alpha=0.3, colour=NA) +
  geom_line() +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.2)) + 
  ggtitle(species)
ggsave(paste0("figs/sTime_A_", species, ".jpg"), width=6, height=3)

ggplot(pred.df %>% filter(Tank==1), aes(ElapsedTime, pred_phi)) +
  geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
             linetype=2, colour="grey30") +
  geom_rect(aes(xmin=-6, xmax=-3, ymax=12, ymin=24), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=-6, xmax=-3, ymax=0, ymin=12), 
            colour="grey30", fill="grey30", size=0.25) +
  geom_ribbon(aes(ymin=pred_phi_lo, ymax=pred_phi_hi), alpha=0.3, colour=NA) +
  geom_line() + 
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) + 
  ggtitle(species)
ggsave(paste0("figs/sTime_phi_", species, ".jpg"), width=6, height=3)
 



