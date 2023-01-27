# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[1]
iter <- 2000
warmup <- 1000
chains <- 4
stan_args <- list(adapt_delta=0.95, max_treedepth=20)



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
source("code/00_fn.R")
chmb_col <- c("1"="#0571b0", "2"="#92c5de", "3"="grey", "4"="#f4a582", "5"="#ca0020")

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
  mutate(GrpDay=factor(paste(as.numeric(Group), str_pad(Days, 2, "l", "0"), sep="_"))) %>%
  mutate(ElapsedTime=ZT+24*(Days-1)+24*3*(Group=="Experiment"))
time_sc <- scale(data.df$ElapsedTime)
data.df <- data.df %>%
  mutate(ElapsedTime_sc=c(time_sc))

data.noNA <- data.df %>% filter(complete.cases(.))

data.noNA %>%
  mutate(Group=if_else(Group=="Control", "Control", "Experiment")) %>%
  group_by(Group, Tank, ElapsedTime) %>%
  mutate(propFish=FishCount/sum(FishCount)) %>%
  ungroup %>%
  select(Group, Tank, ElapsedTime, ZT, Chamber, propFish) %>%
  ggplot(aes(ElapsedTime/24, propFish, fill=Chamber)) + 
  geom_bar(stat="identity", colour=NA) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  facet_grid(Tank~Group) +
  theme_bw()


# dirichlet trial ---------------------------------------------------------

data.wide <- data.noNA %>%
  mutate(Group=if_else(Group=="Control", "Control", "Experiment")) %>%
  group_by(Group, Tank, ElapsedTime_sc) %>%
  mutate(FishCount=FishCount+0.001) %>%
  mutate(propFish=FishCount/sum(FishCount),
         Chamber=paste0("Ch_", Chamber)) %>%
  ungroup %>%
  select(Group, Tank, ElapsedTime_sc, ZT, Chamber, propFish) %>%
  pivot_wider(names_from=Chamber, values_from=propFish) %>%
  filter(Group!="Control") %>%
  # filter(Tank==1) %>%
  mutate(Y=cbind(Ch_1=Ch_1, Ch_2=Ch_2, Ch_3=Ch_3, Ch_4=Ch_4, Ch_5=Ch_5))

data.wide.ctrl <- data.noNA %>%
  mutate(Group=if_else(Group=="Control", "Control", "Experiment")) %>%
  group_by(Group, Tank, ElapsedTime_sc) %>%
  mutate(FishCount=FishCount+0.001) %>%
  mutate(propFish=FishCount/sum(FishCount),
         Chamber=paste0("Ch_", Chamber)) %>%
  ungroup %>%
  select(Group, Tank, ElapsedTime_sc, ZT, Chamber, propFish) %>%
  pivot_wider(names_from=Chamber, values_from=propFish) %>%
  filter(Group=="Control") %>%
  # filter(Tank==1) %>%
  mutate(Y=cbind(Ch_1=Ch_1, Ch_2=Ch_2, Ch_3=Ch_3, Ch_4=Ch_4, Ch_5=Ch_5))




# d_0 <- brm(bf(Y ~ s(ElapsedTime_sc)), 
#            family=dirichlet(), data=data.wide, cores=4, chains=4)
# d_1 <- brm(bf(Y ~ s(ElapsedTime_sc) + s(ZT, bs="cc")), 
#            family=dirichlet(), data=data.wide, cores=4, chains=4)



# dirichlet HGAM ----------------------------------------------------------

# prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A"),
#               prior(normal(0, 2), class="b", coef="Intercept", nlpar="M"),
#               prior(normal(0, 1), class="b", nlpar="M"),
#               prior(von_mises(0, 1), class="b", nlpar="pphi", lb=-3.141593, ub=3.141593),
#               prior(normal(0, 3), class="sds", nlpar="A", lb=0),
#               prior(normal(0, 3), class="sds", nlpar="M", lb=0),
#               prior(normal(0, 3), class="sds", nlpar="pphi", lb=0),
#               prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="A"),
#               prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="M"),
#               prior(normal(0, 3), class="sds", coef="s(ElapsedTime_sc, m = 2)", nlpar="pphi"),
#               prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="A"),
#               prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="M"),
#               prior(normal(0, 0.5), class="sds", coef='s(ElapsedTime_sc, Tank, bs = "fs", m = 1)', nlpar="pphi"))
# 
# out_G <- brm(bf(Y ~ M + exp(A) * cos(3.141593*(ZT)/12 + pphi),
#                 M ~ s(ElapsedTime_sc, m=2),
#                 A ~ s(ElapsedTime_sc, m=2),
#                 pphi ~ s(ElapsedTime_sc, m=2),
#                 nl=TRUE),
#              family=dirichlet(refcat="Ch_3"), prior=prior.nl, 
#              control=stan_args, iter=iter, warmup=warmup, init=0,
#              data=data.wide, cores=1, chains=1, refresh=100,
#              save_model=glue("models/cosinor/ctrl_chProp_G_{species}.stan"),
#              file=glue("models/cosinor/ctrl_chProp_G_{species}"))
# 
# out_GS1 <- brm(bf(Y ~ s(ElapsedTime_sc)),
#               family=dirichlet(refcat="Ch_3"), #prior=prior.nl, 
#               control=stan_args, iter=iter, warmup=warmup, init=0,
#               data=data.wide, cores=1, chains=1, refresh=100,
#               save_model=glue("models/cosinor/ctrl_chProp_GS_ALL_{species}.stan"),
#               file=glue("models/cosinor/ctrl_chProp_GS_ALL_{species}"))
# 
# out_GS2 <- brm(bf(Y ~ M + exp(A)*cos(3.141593*(ZT)/12 + pphi),
#                   M ~ s(ElapsedTime_sc),
#                   A ~ s(ElapsedTime_sc),
#                   pphi ~ s(ElapsedTime_sc),
#                   nl=TRUE),
#                family=dirichlet(refcat="Ch_3"), prior=prior.nl, 
#                control=stan_args, iter=100, warmup=50, init=0,
#                data=data.wide, cores=1, chains=1, refresh=5,
#                save_model=glue("models/cosinor/ctrl_chProp_GS2_{species}.stan"),
#                file=glue("models/cosinor/ctrl_chProp_GS2_{species}"))


prior.nl <- c(prior(normal(0, 1), class="b", nlpar="A1"),
              prior(normal(0, 1), class="b", nlpar="M1"),
              prior(von_mises(0, 1), class="b", nlpar="pphi1", lb=-3.141593, ub=3.141593),
              prior(normal(0, 1), class="b", nlpar="A2"),
              prior(normal(0, 1), class="b", nlpar="M2"),
              prior(von_mises(0, 1), class="b", nlpar="pphi2", lb=-3.141593, ub=3.141593),
              prior(normal(0, 1), class="b", nlpar="A4"),
              prior(normal(0, 1), class="b", nlpar="M4"),
              prior(von_mises(0, 1), class="b", nlpar="pphi4", lb=-3.141593, ub=3.141593),
              prior(normal(0, 1), class="b", nlpar="A5"),
              prior(normal(0, 1), class="b", nlpar="M5"),
              prior(von_mises(0, 1), class="b", nlpar="pphi5", lb=-3.141593, ub=3.141593),
              prior(normal(0, 2), class="sds", nlpar="A1", lb=0),
              prior(normal(0, 2), class="sds", nlpar="M1", lb=0),
              prior(normal(0, 2), class="sds", nlpar="pphi1", lb=0),
              prior(normal(0, 2), class="sds", nlpar="A2", lb=0),
              prior(normal(0, 2), class="sds", nlpar="M2", lb=0),
              prior(normal(0, 2), class="sds", nlpar="pphi2", lb=0),
              prior(normal(0, 2), class="sds", nlpar="A4", lb=0),
              prior(normal(0, 2), class="sds", nlpar="M4", lb=0),
              prior(normal(0, 2), class="sds", nlpar="pphi4", lb=0),
              prior(normal(0, 2), class="sds", nlpar="A5", lb=0),
              prior(normal(0, 2), class="sds", nlpar="M5", lb=0),
              prior(normal(0, 2), class="sds", nlpar="pphi5", lb=0))
              
# out_GS3 <- brm(bf(Y ~ Ch) +
#                  nlf(muCh1 ~ M1 + exp(A1)*cos(3.141593*(ZT)/12 + pphi1)) +
#                  lf(M1 ~ s(ElapsedTime_sc)) +
#                  lf(A1 ~ s(ElapsedTime_sc)) + 
#                  nlf(muCh2 ~ M2 + exp(A2)*cos(3.141593*(ZT)/12 + pphi2)) +
#                  lf(M2 ~ s(ElapsedTime_sc)) +
#                  lf(A2 ~ s(ElapsedTime_sc)) + 
#                  lf(pphi2 ~ s(ElapsedTime_sc)) +
#                  lf(pphi1 ~ s(ElapsedTime_sc)) +
#                  nlf(muCh4 ~ M4 + exp(A4)*cos(3.141593*(ZT)/12 + pphi4)) +
#                  lf(M4 ~ s(ElapsedTime_sc)) +
#                  lf(A4 ~ s(ElapsedTime_sc)) + 
#                  lf(pphi4 ~ s(ElapsedTime_sc)) +
#                  nlf(muCh5 ~ M5 + exp(A5)*cos(3.141593*(ZT)/12 + pphi5)) +
#                  lf(M5 ~ s(ElapsedTime_sc)) +
#                  lf(A5 ~ s(ElapsedTime_sc)) + 
#                  lf(pphi5 ~ s(ElapsedTime_sc)),
#                family=dirichlet(refcat="Ch_1"), prior=prior.nl, 
#                control=stan_args, iter=iter, warmup=warmup, init=0,
#                data=data.wide, cores=1, chains=1, refresh=5,
#                save_model=glue("models/cosinor/ctrl_chProp_GS3_{species}.stan"),
#                file=glue("models/cosinor/ctrl_chProp_GS3_{species}"))

out_GS4 <- brm(bf(Y ~ Ch) +
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
               family=dirichlet(refcat="Ch_3"), prior=prior.nl, 
               control=stan_args, iter=iter, warmup=warmup, init=0,
               data=data.wide, cores=chains, chains=chains, refresh=10,
               save_model=glue("models/cosinor/chProp_GS4_{species}.stan"),
               file=glue("models/cosinor/chProp_GS4_{species}"))

out_GS4.ctrl <- brm(bf(Y ~ Ch) +
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
               family=dirichlet(refcat="Ch_3"), prior=prior.nl, 
               control=stan_args, iter=iter, warmup=warmup, init=0,
               data=data.wide.ctrl, cores=chains, chains=chains, refresh=10,
               save_model=glue("models/cosinor/chProp_GS4_ctrl_{species}.stan"),
               file=glue("models/cosinor/chProp_GS4_ctrl_{species}"))



# 
# 
# 
# pred_mu.dat <- expand_grid(Days=1:13, ZT=0:24) %>%
#   mutate(ElapsedTime=ZT+24*(Days-1),
#          ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
#          ElapsedDays=ElapsedTime/24) 
# 
# post_pred <- posterior_epred(out_GS3, newdata=pred_mu.dat)
# post_mn <- apply(post_pred, 2:3, mean)
# post_hdi <- apply(post_pred, 2:3, HDInterval::hdi)
# 
# post_M <- map(2:5, ~(posterior_smooths(out_GS3, newdata=pred_mu.dat, 
#                                        smooth="s(ElapsedTime_sc)", nlpar=paste0("M", .x)) +
#                        c(as_draws_matrix(out_GS3, variable=paste0("b_M", .x, "_Intercept")))))
# post_M.ar <- array(0, dim=dim(post_pred))
# post_M.ar[,,2] <- exp(post_M[[1]])
# post_M.ar[,,3] <- exp(post_M[[2]])
# post_M.ar[,,4] <- exp(post_M[[3]])
# post_M.ar[,,5] <- exp(post_M[[4]])
# 
# for(i in 1:dim(post_M.ar)[1]) {
#   for(j in 1:dim(post_M.ar)[2]) {
#     post_M.ar[i,j,] <- post_M.ar[i,j,]/sum(post_M.ar[i,j,])
#   }
# }
# post_M_mn <- apply(post_M.ar, 2:3, mean)
# post_M_hdi <- apply(post_M.ar, 2:3, HDInterval::hdi)
# 
# fit.df <- tibble(Chamber=factor(rep(1:5, each=nrow(pred_mu.dat))),
#                  ElapsedDays=rep(pred_mu.dat$ElapsedDays, 5)) %>%
#   mutate(pr_mn=c(post_mn),
#          pr_lo=c(post_hdi[1,,]),
#          pr_hi=c(post_hdi[2,,]),
#          M_mn=c(post_M_mn),
#          M_lo=c(post_M_hdi[1,,]),
#          M_hi=c(post_M_hdi[2,,]))
# 
# 
# 
# ggplot(fit.df, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi, 
#                    group=Chamber, colour=Chamber, fill=Chamber)) +
#   geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
#   scale_colour_manual(values=chmb_col) +
#   scale_fill_manual(values=chmb_col) +
#   scale_x_continuous("Elapsed time (days)", breaks=0:13) +
#   scale_y_continuous("Predicted proportion of fish", limits=c(0, 0.5)) +
#   ggtitle(species) +
#   theme_classic()
# 
# ggplot(fit.df, aes(ElapsedDays, pr_mn, fill=Chamber, group=paste(ElapsedDays, Chamber))) +
#   geom_bar(stat="identity", position="fill", colour=NA) +
#   scale_fill_manual(values=chmb_col) +
#   scale_x_continuous("Elapsed time (days)", breaks=0:13) +
#   scale_y_continuous("Predicted proportion of fish", limits=c(0, 1)) +
#   ggtitle(species) +
#   theme_classic()
# 
# 
# ggplot(fit.df, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi, 
#                    group=Chamber, colour=Chamber, fill=Chamber)) +
#   geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
#   scale_colour_manual(values=chmb_col) +
#   scale_fill_manual(values=chmb_col) +
#   scale_x_continuous("Elapsed time (days)", breaks=0:13) +
#   # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
#   ggtitle(species) +
#   theme_classic()
