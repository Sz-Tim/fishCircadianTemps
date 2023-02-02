# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")[2]
iter <- 4000
warmup <- 2000
chains <- 4
stan_args <- list(adapt_delta=0.9, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
theme_set(theme_classic())
source("code/00_fn.R")

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
  mutate(ElapsedTime=ZT+24*(Days-1)+24*3*(Group=="Experiment"))
time_sc <- scale(data.df$ElapsedTime)
data.df <- data.df %>%
  mutate(ElapsedTime_sc=c(time_sc)) %>%
  filter(Group!="Control")



# model prep --------------------------------------------------------------

data.noNA <- data.df %>% filter(complete.cases(.)) %>% 
  mutate(Days=factor(Days),
         Tank=factor(Tank),
         ElapsedDays=ElapsedTime/24)

if(species=="ZF") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(27.5, 2), class="b", coef="Intercept", nlpar="M"),
                prior(normal(0, 1), class="b", nlpar="M"),
                prior(von_mises(0, 1), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                # prior(normal(0, 0.2), class="sd", nlpar="A", lb=0),
                # prior(normal(0, 0.2), class="sd", nlpar="M", lb=0),
                # prior(normal(0, 0.2), class="sd", nlpar="phi", lb=0),
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
} else if(species=="Tilapia") {
  prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
                prior(normal(30, 2), class="b", coef="Intercept", nlpar="M"),
                prior(normal(0, 1), class="b", nlpar="M"),
                prior(von_mises(0, 1), class="b", nlpar="phi", lb=-3.141593, ub=3.141593),
                # prior(normal(0, 0.2), class="sd", nlpar="A", lb=0),
                # prior(normal(0, 0.2), class="sd", nlpar="M", lb=0),
                # prior(normal(0, 0.2), class="sd", nlpar="phi", lb=0),
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
}


# fit model ---------------------------------------------------------------

# Following Pedersen et al 2019
# out_G <- brm(bf(prefTemp ~ M + exp(A) * cos(3.141593*(ZT)/12 + phi),
#                 M ~ s(ElapsedTime_sc) + (1|Tank),
#                 A ~ s(ElapsedTime_sc) + (1|Tank),
#                 phi ~ s(ElapsedTime_sc) + (1|Tank),
#                 nl=TRUE),
#              prior=prior.nl, 
#              control=stan_args, iter=iter, warmup=warmup, init=0,
#              data=data.noNA, cores=chains, chains=chains, refresh=100,
#              save_model=glue("models/cosinor/s_G_{species}.stan"),
#              file=glue("models/cosinor/s_G_{species}"))

out_GS <- brm(bf(prefTemp ~ M + exp(A) * cos(3.141593*(ZT)/12 + phi),
                 M ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
                 A ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
                 phi ~ s(ElapsedTime_sc, m=2) + s(ElapsedTime_sc, Tank, bs="fs", m=1),
                 nl=TRUE),
              prior=prior.nl, 
              control=stan_args, iter=iter, warmup=warmup, init=0,
              data=data.noNA, cores=chains, chains=chains, refresh=100,
              save_model=glue("models/cosinor/prefTemp_GS_{species}.stan"),
              file=glue("models/cosinor/prefTemp_GS_{species}"))


pred_mu.dat <- expand_grid(Days=1:13, ZT=0:24) %>%
  mutate(ElapsedTime=ZT+24*(Days-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24) 
smooths_GS <- list(
  M=posterior_smooths(out_GS, newdata=pred_mu.dat, smooth="s(ElapsedTime_sc,m=2)", nlpar="M") +
    c(as_draws_matrix(out_GS, variable="b_M_Intercept")),
  A=posterior_smooths(out_GS, newdata=pred_mu.dat, smooth="s(ElapsedTime_sc,m=2)", nlpar="A") +
    c(as_draws_matrix(out_GS, variable="b_A_Intercept")),
  phi=posterior_smooths(out_GS, newdata=pred_mu.dat, smooth="s(ElapsedTime_sc,m=2)", nlpar="phi") +
    c(as_draws_matrix(out_GS, variable="b_phi_Intercept"))
)
smooths_GS$amplitude <- exp(smooths_GS$A)
smooths_GS$acrophase <- phi_to_ZT(smooths_GS$phi) 
smooths_GS$prefTemp <- calc_prefTemp(smooths_GS$M, smooths_GS$A, smooths_GS$phi, 
                                     matrix(pred_mu.dat$ZT, nrow=nrow(smooths_GS$M), ncol=nrow(pred_mu.dat), byrow=T))
hdis <- list(prefTemp=HDInterval::hdi(smooths_GS$prefTemp),
             M=HDInterval::hdi(smooths_GS$M),
             amplitude=HDInterval::hdi(smooths_GS$amplitude),
             acrophase=HDInterval::hdi(smooths_GS$acrophase))
pred_mu.df <- pred_mu.dat %>% 
  mutate(pred=colMeans(smooths_GS$prefTemp),
         pred_lo=hdis$prefTemp[1,],
         pred_hi=hdis$prefTemp[2,],
         M=colMeans(smooths_GS$M),
         M_lo=hdis$M[1,],
         M_hi=hdis$M[2,],
         amplitude=colMeans(smooths_GS$amplitude),
         amplitude_lo=hdis$amplitude[1,],
         amplitude_hi=hdis$amplitude[2,],
         acrophase=colMeans(smooths_GS$acrophase),
         acrophase_lo=hdis$acrophase[1,],
         acrophase_hi=hdis$acrophase[2,])



pred.dat <- expand_grid(ZT=0:24, Days=1:13, Tank=1:3) %>%
  mutate(ElapsedTime=ZT+24*(Days-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24) 
pred.ls  <- list(
  prefTemp=posterior_epred(out_GS, newdata=pred.dat, re_formula=NULL),
  M=posterior_epred(out_GS, newdata=pred.dat, re_formula=NULL, nlpar="M"),
  A=posterior_epred(out_GS, newdata=pred.dat, re_formula=NULL, nlpar="A"),
  phi=posterior_epred(out_GS, newdata=pred.dat, re_formula=NULL, nlpar="phi")
)
pred.ls$A <- exp(pred.ls$A)
pred.ls$phi <- phi_to_ZT(pred.ls$phi)
pred.df <- pred.dat %>%
  mutate(pred=colMeans(pred.ls$prefTemp),
         pred_lo=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.025)),
         pred_hi=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.975)),
         M=colMeans(pred.ls$M),
         M_lo=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.025)),
         M_hi=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.975)),
         amplitude=colMeans(pred.ls$A),
         amplitude_lo=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.025)),
         amplitude_hi=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.975)),
         acrophase=colMeans(pred.ls$phi),
         acrophase_lo=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.025)),
         acrophase_hi=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.975))) %>%
  mutate(Tank=factor(Tank))

pred.df %>% 
  mutate(Species=names(species)) %>%
  write_csv(glue("out/GS_predTank_{species}.csv"))
pred_mu.df %>% 
  mutate(Species=names(species)) %>%
  write_csv(glue("out/GS_predGlobal_{species}.csv"))

# 
# 
# temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))
# 
# ggplot(pred.df, aes(ElapsedDays, pred)) +
#   geom_point(data=data.noNA, aes(y=prefTemp, colour=Tank), size=0.5, shape=1) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=pred_lo, ymax=pred_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("Preferred temperature (ºC)") + ggtitle(species)
# ggsave(paste0("figs/s_prefTemp_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, M)) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=M_lo, ymax=M_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("MESOR (ºC)") + ggtitle(species)
# ggsave(paste0("figs/s_M_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, amplitude)) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=amplitude_lo, ymax=amplitude_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.4)) +
#   ggtitle(species)
# ggsave(paste0("figs/s_A_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, acrophase)) +
#   geom_rect(aes(xmin=-6/24, xmax=-3/24, ymax=12, ymin=24),
#             colour="grey30", fill="white", size=0.25) +
#   geom_rect(aes(xmin=-6/24, xmax=-3/24, ymax=0, ymin=12),
#             colour="grey30", fill="grey30", size=0.25) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=acrophase_lo, ymax=acrophase_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) +
#   ggtitle(species)
# ggsave(paste0("figs/s_phi_", species, ".jpg"), width=6, height=3)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# plot(NA, NA, xlim=c(0,325), ylim=c(28, 35), xlab="Hours", ylab="M")
# for(i in 1:200) {
#   lines(pred_mu.dat$ElapsedTime, smooths_GS$M[i,], col=rgb(0,0,0,0.25))
# }
# 
# plot(NA, NA, xlim=c(0,325), ylim=c(0, 2), xlab="Hours", ylab="Amplitude")
# for(i in 1:200) {
#   lines(pred_mu.dat$ElapsedTime, smooths_GS$amplitude[i,], col=rgb(0,0,0,0.25))
# }
# 
# plot(NA, NA, xlim=c(0,325), ylim=c(-4, 3), xlab="Hours", ylab="A")
# for(i in 1:200) {
#   lines(pred_mu.dat$ElapsedTime, smooths_GS$A[i,], col=rgb(0,0,0,0.25))
# }
# 
# plot(NA, NA, xlim=c(0,325), ylim=c(0, 24), xlab="Hours", ylab="Acrophase")
# for(i in 1:200) {
#   lines(pred_mu.dat$ElapsedTime, smooths_GS$acrophase[i,], col=rgb(0,0,0,0.25))
# }
# 
# plot(NA, NA, xlim=c(0,325), ylim=c(25.5, 35), xlab="Hours", ylab="prefTemp")
# for(i in 1:200) {
#   lines(pred_mu.dat$ElapsedTime, smooths_GS$prefTemp[i,], col=rgb(0,0,0,0.25))
# }
# 
# 
# 
# 
# pred_mu.dat <- expand_grid(ZT=0:24, Days=1:13) %>%
#   mutate(ElapsedTime=ZT+24*(Days-1),
#          ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"))
# 
# pred_mu.ls  <- list(
#   prefTemp=posterior_epred(out, newdata=pred_mu.dat, re_formula=NA),
#   M=posterior_epred(out, newdata=pred_mu.dat, re_formula=NA, nlpar="M"),
#   A=posterior_epred(out, newdata=pred_mu.dat, re_formula=NA, nlpar="A"),
#   phi=posterior_epred(out, newdata=pred_mu.dat, re_formula=NA, nlpar="phi")
# )
# pred_mu.ls  <- list(
#   prefTemp=posterior_predict(out, newdata=pred_mu.dat, re_formula=NA),
#   M=posterior_predict(out, newdata=pred_mu.dat, re_formula=NA, nlpar="M"),
#   A=posterior_predict(out, newdata=pred_mu.dat, re_formula=NA, nlpar="A"),
#   phi=posterior_predict(out, newdata=pred_mu.dat, re_formula=NA, nlpar="phi")
# )
# pred_mu.ls$A <- exp(pred_mu.ls$A)
# pred_mu.ls$phi <- phi_to_ZT(pred_mu.ls$phi)
# pred_mu.df <- pred_mu.dat %>%
#   mutate(pred=colMeans(pred_mu.ls$prefTemp),
#          pred_lo=apply(pred_mu.ls$prefTemp, 2, function(x) quantile(x, probs=0.025)),
#          pred_hi=apply(pred_mu.ls$prefTemp, 2, function(x) quantile(x, probs=0.975)),
#          pred_M=colMeans(pred_mu.ls$M),
#          pred_M_lo=apply(pred_mu.ls$M, 2, function(x) quantile(x, probs=0.025)),
#          pred_M_hi=apply(pred_mu.ls$M, 2, function(x) quantile(x, probs=0.975)),
#          pred_A=colMeans(pred_mu.ls$A),
#          pred_A_lo=apply(pred_mu.ls$A, 2, function(x) quantile(x, probs=0.025)),
#          pred_A_hi=apply(pred_mu.ls$A, 2, function(x) quantile(x, probs=0.975)),
#          pred_phi=colMeans(pred_mu.ls$phi),
#          pred_phi_lo=apply(pred_mu.ls$phi, 2, function(x) quantile(x, probs=0.025)),
#          pred_phi_hi=apply(pred_mu.ls$phi, 2, function(x) quantile(x, probs=0.975)))
# 
# temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))
# 
# ggplot(pred_mu.df, aes(ElapsedTime, pred)) +
#   geom_point(data=data.noNA, aes(y=prefTemp, group=Tank, colour=Tank), size=0.25) +
#   geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("Preferred temperature (ºC)") + ggtitle(species)
# ggsave(paste0("figs/sTime_prefTemp_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred_mu.df, aes(ElapsedTime, pred_M)) +
#   geom_ribbon(aes(ymin=pred_M_lo, ymax=pred_M_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("MESOR (ºC)") + ggtitle(species)
# ggsave(paste0("figs/sTime_M_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred_mu.df, aes(ElapsedTime, pred_A)) +
#   geom_ribbon(aes(ymin=pred_A_lo, ymax=pred_A_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.2)) +
#   ggtitle(species)
# ggsave(paste0("figs/sTime_A_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred_mu.df, aes(ElapsedTime, pred_phi)) +
#   geom_rect(aes(xmin=-6, xmax=-3, ymax=12, ymin=24),
#             colour="grey30", fill="white", size=0.25) +
#   geom_rect(aes(xmin=-6, xmax=-3, ymax=0, ymin=12),
#             colour="grey30", fill="grey30", size=0.25) +
#   geom_ribbon(aes(ymin=pred_phi_lo, ymax=pred_phi_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) +
#   ggtitle(species)
# ggsave(paste0("figs/sTime_phi_", species, ".jpg"), width=6, height=3)
# 
# 
# 
# 
# 
# 
# 
# 
# pred.dat <- expand_grid(ZT=0:24, Days=1:13, Tank=1:3) %>%
#   mutate(ElapsedTime=ZT+24*(Days-1),
#          ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"))
# pred.ls  <- list(
#   prefTemp=posterior_epred(out, newdata=pred.dat, re_formula=NULL),
#   M=posterior_epred(out, newdata=pred.dat, re_formula=NULL, nlpar="M"),
#   A=posterior_epred(out, newdata=pred.dat, re_formula=NULL, nlpar="A"),
#   phi=posterior_epred(out, newdata=pred.dat, re_formula=NULL, nlpar="phi")
# )
# pred.ls$A <- exp(pred.ls$A)
# pred.ls$phi <- phi_to_ZT(pred.ls$phi)
# pred.df <- pred.dat %>%
#   mutate(pred=colMeans(pred.ls$prefTemp),
#          pred_lo=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.025)),
#          pred_hi=apply(pred.ls$prefTemp, 2, function(x) quantile(x, probs=0.975)),
#          pred_M=colMeans(pred.ls$M),
#          pred_M_lo=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.025)),
#          pred_M_hi=apply(pred.ls$M, 2, function(x) quantile(x, probs=0.975)),
#          pred_A=colMeans(pred.ls$A),
#          pred_A_lo=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.025)),
#          pred_A_hi=apply(pred.ls$A, 2, function(x) quantile(x, probs=0.975)),
#          pred_phi=colMeans(pred.ls$phi),
#          pred_phi_lo=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.025)),
#          pred_phi_hi=apply(pred.ls$phi, 2, function(x) quantile(x, probs=0.975))) %>%
#   mutate(Tank=factor(Tank))
# 
# temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))
# 
# ggplot(pred.df, aes(ElapsedTime, pred, group=Tank, colour=Tank, fill=Tank)) +
#   geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
#              linetype=2, colour="grey30") +
#   geom_point(data=data.noNA, aes(y=prefTemp), size=0.25) +
#   geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("Preferred temperature (ºC)") + ggtitle(species)
# ggsave(paste0("figs/s_prefTemp_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime, pred_M, group=Tank, colour=Tank, fill=Tank)) +
#   geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
#              linetype=2, colour="grey30") +
#   geom_ribbon(aes(ymin=pred_M_lo, ymax=pred_M_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("MESOR (ºC)") + ggtitle(species)
# ggsave(paste0("figs/sTime_M_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime, pred_A, group=Tank, colour=Tank, fill=Tank)) +
#   geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
#              linetype=2, colour="grey30") +
#   geom_ribbon(aes(ymin=pred_A_lo, ymax=pred_A_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1, 1.5), limits=c(0, 1.7)) +
#   ggtitle(species)
# ggsave(paste0("figs/sTime_A_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime, pred_phi, group=Tank, colour=Tank, fill=Tank)) +
#   geom_vline(xintercept=max(filter(data.noNA, Group=="Acclimation")$ElapsedTime),
#              linetype=2, colour="grey30") +
#   geom_rect(aes(xmin=-6, xmax=-3, ymax=12, ymin=24),
#             colour="grey30", fill="white", size=0.25) +
#   geom_rect(aes(xmin=-6, xmax=-3, ymax=0, ymin=12),
#             colour="grey30", fill="grey30", size=0.25) +
#   geom_ribbon(aes(ymin=pred_phi_lo, ymax=pred_phi_hi), alpha=0.3, colour=NA) +
#   geom_line() +
#   scale_x_continuous("Elapsed hours", breaks=24*(0:13)) +
#   scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) +
#   ggtitle(species)
# ggsave(paste0("figs/sTime_phi_", species, ".jpg"), width=6, height=3)
# 
# 
# 
# ggplot(pred.df, aes(ElapsedTime/24, pred)) +
#   geom_point(data=data.noNA, aes(y=prefTemp, colour=Tank), size=0.5, shape=1) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=pred_lo, ymax=pred_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("Preferred temperature (ºC)") + ggtitle(species)
# ggsave(paste0("figs/s_prefTemp_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, pred_M)) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=pred_M_lo, ymax=pred_M_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   ylim(temp_rng[[species]][1], temp_rng[[species]][2]) +
#   ylab("MESOR (ºC)") + ggtitle(species)
# ggsave(paste0("figs/s_M_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, pred_A)) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=pred_A_lo, ymax=pred_A_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.4)) +
#   ggtitle(species)
# ggsave(paste0("figs/s_A_", species, ".jpg"), width=6, height=3)
# 
# ggplot(pred.df, aes(ElapsedTime/24, pred_phi)) +
#   geom_rect(aes(xmin=-6/24, xmax=-3/24, ymax=12, ymin=24),
#             colour="grey30", fill="white", size=0.25) +
#   geom_rect(aes(xmin=-6/24, xmax=-3/24, ymax=0, ymin=12),
#             colour="grey30", fill="grey30", size=0.25) +
#   geom_line(aes(group=Tank, colour=Tank)) +
#   geom_ribbon(data=pred_mu.df, aes(ymin=pred_phi_lo, ymax=pred_phi_hi), colour=NA, alpha=0.2) +
#   geom_line(data=pred_mu.df, size=1) +
#   scale_x_continuous("Time (days)", breaks=0:13) +
#   scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) +
#   ggtitle(species)
# ggsave(paste0("figs/s_phi_", species, ".jpg"), width=6, height=3)
