# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")[2]
iter <- 4000
warmup <- 2000
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
  # filter(ElapsedTime < (24*6)) %>%
  ggplot(aes(ElapsedTime/24, propFish, fill=Chamber)) + 
  geom_bar(stat="identity", colour="grey30") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  facet_grid(Tank~Group, scales="free_x", space="free_x") +
  theme_bw()

data.noNA %>%
  mutate(Group=if_else(Group=="Control", "Control", "Experiment"),
         Tank=paste("Tank", Tank)) %>%
  group_by(Group, Tank, ElapsedTime) %>%
  mutate(propFish=FishCount/sum(FishCount)) %>%
  ungroup %>%
  select(Group, Tank, ElapsedTime, ZT, Chamber, propFish) %>%
  # filter(ElapsedTime < (24*6)) %>%
  ggplot(aes(ElapsedTime/24, propFish, fill=Chamber)) + 
  geom_area() +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Observed fish distribution", breaks=seq(0,1,by=0.2)) +
  ggtitle(names(species)) +
  facet_grid(Tank~Group, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")


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



fit.ls <- list(
  exp=summarise_N_posteriors(out_GS4, 
                             expand_grid(Days=1:13, ZT=0:23, Tank=3) %>%
                               mutate(ElapsedTime=ZT+24*(Days-1),
                                      ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
                                      ElapsedDays=ElapsedTime/24,
                                      Species=species,
                                      Group="Experimental")),
  ctrl=summarise_N_posteriors(out_GS4.ctrl, 
                              expand_grid(Days=1:6, ZT=0:23, Tank=3) %>%
                                mutate(ElapsedTime=ZT+24*(Days-1),
                                       ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
                                       ElapsedDays=ElapsedTime/24,
                                       Species=species,
                                       Group="Control"))
)
fit.df <- rbind(fit.ls$exp$summaries,
                fit.ls$ctrl$summaries)
pr.post <- list(exp=fit.ls$exp$pr.post,
                ctrl=fit.ls$ctrl$pr.post)
M.post <- list(exp=fit.ls$exp$M.post,
               ctrl=fit.ls$ctrl$M.post)
rm(fit.ls)

write_csv(fit.df, glue("out/GS_predGlobal_NChmbr_{species}.csv"))



fit.dat <- bind_rows(
  expand_grid(Days=1:13, ZT=0:23, Tank=3) %>%
    mutate(ElapsedTime=ZT+24*(Days-1),
           ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
           ElapsedDays=ElapsedTime/24,
           Species=species,
           Group="Experimental"),
  expand_grid(Days=1:6, ZT=0:23, Tank=3) %>%
    mutate(ElapsedTime=ZT+24*(Days-1),
           ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
           ElapsedDays=ElapsedTime/24,
           Species=species,
           Group="Control"))

pr.edge <- array(dim=c(dim(pr.post$exp)[1], 
                       dim(pr.post$exp)[2]+dim(pr.post$ctrl)[2],
                       2))
pr.edge[,,1] <- cbind(apply(pr.post$exp[,,c(1,5)], 1:2, sum), 
                      apply(pr.post$ctrl[,,c(1,5)], 1:2, sum))
pr.edge[,,2] <- cbind(apply(pr.post$exp[,,2:4], 1:2, sum), 
                      apply(pr.post$ctrl[,,2:4], 1:2, sum))

M.edge <- array(dim=c(dim(M.post$exp)[1], 
                      dim(M.post$exp)[2]+dim(M.post$ctrl)[2],
                      2))
M.edge[,,1] <- cbind(apply(M.post$exp[,,c(1,5)], 1:2, sum), 
                     apply(M.post$ctrl[,,c(1,5)], 1:2, sum))
M.edge[,,2] <- cbind(apply(M.post$exp[,,2:4], 1:2, sum), 
                     apply(M.post$ctrl[,,2:4], 1:2, sum))
edge.df <- fit.dat %>%
  mutate(edge_mn=colMeans(pr.edge[,,1]),
         edge_lo=HDInterval::hdi(pr.edge[,,1])[1,],
         edge_hi=HDInterval::hdi(pr.edge[,,1])[2,],
         M_mn=colMeans(M.edge[,,1]),
         M_lo=HDInterval::hdi(M.edge[,,1])[1,],
         M_hi=HDInterval::hdi(M.edge[,,1])[2,])
write_csv(edge.df, glue("out/GS_predEdge_NChmbr_{species}.csv"))


pr.cw <- array(dim=c(dim(pr.post$exp)[1], 
                     dim(pr.post$exp)[2]+dim(pr.post$ctrl)[2],
                     2))
pr.cw[,,1] <- cbind(apply(pr.post$exp[,,1:2], 1:2, sum), 
                    apply(pr.post$ctrl[,,1:2], 1:2, sum))
pr.cw[,,2] <- cbind(apply(pr.post$exp[,,4:5], 1:2, sum), 
                    apply(pr.post$ctrl[,,4:5], 1:2, sum))
M.cw <- array(dim=c(dim(M.post$exp)[1], 
                     dim(M.post$exp)[2]+dim(M.post$ctrl)[2],
                     2))
M.cw[,,1] <- cbind(apply(M.post$exp[,,1:2], 1:2, sum), 
                    apply(M.post$ctrl[,,1:2], 1:2, sum))
M.cw[,,2] <- cbind(apply(M.post$exp[,,4:5], 1:2, sum), 
                    apply(M.post$ctrl[,,4:5], 1:2, sum))
cw.df <-  bind_rows(
  fit.dat %>%
    mutate(Chambers="cold",
           pr_mn=colMeans(pr.cw[,,1]),
           pr_lo=HDInterval::hdi(pr.cw[,,1])[1,],
           pr_hi=HDInterval::hdi(pr.cw[,,1])[2,],
           M_mn=colMeans(M.cw[,,1]),
           M_lo=HDInterval::hdi(M.cw[,,1])[1,],
           M_hi=HDInterval::hdi(M.cw[,,1])[2,]),
  fit.dat %>%
    mutate(Chambers="warm",
           pr_mn=colMeans(pr.cw[,,2]),
           pr_lo=HDInterval::hdi(pr.cw[,,2])[1,],
           pr_hi=HDInterval::hdi(pr.cw[,,2])[2,],
           M_mn=colMeans(M.cw[,,2]),
           M_lo=HDInterval::hdi(M.cw[,,2])[1,],
           M_hi=HDInterval::hdi(M.cw[,,2])[2,]))
write_csv(cw.df, glue("out/GS_predColdWarm_NChmbr_{species}.csv"))




ggplot(edge.df, aes(ElapsedDays, edge_mn, ymin=edge_lo, ymax=edge_hi, colour=Group, fill=Group)) + 
  geom_hline(yintercept=0.4, colour="grey") +
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  scale_y_continuous("Pr(edge chamber)", limits=c(0,1)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2))
ggplot(edge.df, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi, colour=Group, fill=Group)) + 
  geom_hline(yintercept=0.4, colour="grey") +
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  scale_y_continuous("M: Pr(edge chamber)", limits=c(0,1)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2))


ggplot(cw.df, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi, colour=Group, fill=Group)) + 
  geom_hline(yintercept=0.4, colour="grey") +
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  scale_y_continuous("Pr(cold chamber)", limits=c(0,1)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  facet_grid(.~Chambers)
ggplot(cw.df, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi, colour=Group, fill=Group)) + 
  geom_hline(yintercept=0.4, colour="grey") +
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  scale_y_continuous("Pr(cold chamber)", limits=c(0,1)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  facet_grid(.~Chambers)




ggplot(fit.df, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                   group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() + 
  facet_grid(Group~.)

ggplot(fit.df, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() + 
  facet_grid(Chamber~.)

ggplot(fit.df, aes(ElapsedDays, pr_mn, fill=Chamber)) +
  geom_area(colour="grey30") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, 1)) +
  ggtitle(species) +
  theme_classic() + 
  facet_grid(Group~.)

ggplot(fit.df, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                   group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() +
  facet_grid(Group~.)

ggplot(fit.df, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() +
  facet_grid(Chamber~.)

ggplot(fit.df, aes(ElapsedDays, M_mn, fill=Chamber)) +
  geom_area(colour="black") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish (MESOR)", limits=c(0, 1)) +
  ggtitle(species) +
  theme_classic() +
  facet_grid(Group~.)

ggplot(fit.df, aes(ElapsedDays, A_mn, ymin=A_lo, ymax=A_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() +
  facet_grid(Chamber~.)

ggplot(fit.df, aes(ElapsedDays, phi_mn, ymin=phi_lo, ymax=phi_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  ggtitle(species) +
  theme_classic() +
  facet_grid(Chamber~.)
