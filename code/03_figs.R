# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
theme_set(theme_classic())
group_col <- c(`Control CTE`="#767676", Acclimation="#800000", Experiment="#155F83")
cmr.ls <- readRDS("figs/cmr_cmaps.RDS")

species <- c("Nile tilapia"="Tilapia", "Zebrafish"="ZF")
temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))

data.df <- map(species, 
               ~readxl::read_xlsx(dir("data", glue("{.x}_.*RawData2"), full.names=T),
                                  1, col_types=c(rep("numeric", 6), "skip", "skip")) %>%
                 mutate(ln_FishCount=log(FishCount+1),
                        Chamber=factor(Chamber),
                        Group=factor(Group, 
                                     levels=c(3, 1, 2), 
                                     labels=c("Control", "Acclimation", "Experiment")),
                        Tank=factor(Tank)) %>%
                 left_join(read_csv(glue("data/temp_{.x}.csv")) %>%
                             pivot_longer(starts_with("chamber_"), 
                                          names_to="Chamber", values_to="Temp") %>%
                             mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                                    Group=factor(Group, 
                                                 levels=c("Control", 
                                                          "Acclimation", 
                                                          "Experiment")))) %>%
  group_by(ZT, Group, Tank, Days) %>%
  summarise(prefTemp=sum(FishCount*Temp)/(sum(FishCount))) %>%
  ungroup)

data.sum <- map(data.df, 
                ~.x %>%
                  group_by(ZT, Group, Tank) %>%
                  summarise(mn=mean(prefTemp, na.rm=T),
                            se=sd(prefTemp, na.rm=T)/sqrt(sum(!is.na(prefTemp)))))

data.c <- map(species, 
              ~readxl::read_xlsx(dir("data", glue("{.x}_.*RawData2"), full.names=T),
                            1, col_types=c(rep("numeric", 6), "skip", "skip")) %>%
  mutate(ln_FishCount=log(FishCount+1),
         cos_ZT=cos(2*pi*ZT/24),
         sin_ZT=sin(2*pi*ZT/24),
         Chamber=factor(Chamber),
         Group=factor(Group, 
                      levels=c(3, 1, 2), 
                      labels=c("Control", "Acclimation", "Experiment")),
         Tank=factor(Tank)) %>%
  left_join(read_csv(glue("data/temp_{.x}.csv")) %>%
              pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
              mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                     Group=factor(Group, levels=c("Control", "Acclimation", "Experiment")))) %>%
  mutate(GrpDay=factor(paste(as.numeric(Group), str_pad(Days, 2, "l", "0"), sep="_"))))

out <- map(species, 
           ~readRDS(glue("models/noCorr_randEff/out_temperature_RS_{.x}.rds")))
out.c <- map(species, 
           ~readRDS(glue("models/noCorr_randEff/out_count_RS_{.x}.rds")))

pred.df <- map(data.df, 
               ~expand_grid(ZT=unique(.x$ZT),
                            Group=unique(.x$Group)))
preds <- map2(out, pred.df, 
              ~posterior_epred(.x, newdata=.y, re.form=NA))
pred.df <- map2(pred.df, preds, 
                ~.x %>%
                  mutate(pred.mn=apply(.y, 2, mean),
                         CI95_l=apply(.y, 2, function(x) quantile(x, 0.025)),
                         CI90_l=apply(.y, 2, function(x) quantile(x, 0.05)),
                         CI80_l=apply(.y, 2, function(x) quantile(x, 0.1)),
                         CI80_h=apply(.y, 2, function(x) quantile(x, 0.9)),
                         CI90_h=apply(.y, 2, function(x) quantile(x, 0.95)),
                         CI95_h=apply(.y, 2, function(x) quantile(x, 0.975))))



# pref temp by ZT ---------------------------------------------------------

fig1.ls <- map(2:1, 
               ~pred.df[[.x]] %>%
                 mutate(Group=factor(Group, 
                                     labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_jitter(data=data.df[[.x]] %>%
                mutate(Group=factor(Group, 
                                    labels=c("Control CTE", "Acclimation", "Experiment"))), 
              aes(y=prefTemp, shape=Group),
              alpha=0.8, size=0.95, width=0.4, height=0) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.4) +
  geom_rect(aes(xmin=0, xmax=12, ymax=temp_rng[[.x]][2],
                ymin=temp_rng[[.x]][2]-diff(temp_rng[[.x]])*0.02), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymax=temp_rng[[.x]][2],
                ymin=temp_rng[[.x]][2]-diff(temp_rng[[.x]])*0.02), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  scale_x_continuous("ZT", breaks=seq(0, 24, by=6)) +
  scale_y_continuous("Preferred temperature (ºC)", limits=temp_rng[[.x]],
                     breaks=seq(ceiling(temp_rng[[.x]][1]),
                                floor(temp_rng[[.x]][2]),
                                by=2)) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank(), 
        legend.background=element_blank()))

fig.1 <- ggpubr::ggarrange(plotlist=fig1.ls, ncol=1, common.legend=T, 
                           legend="bottom", labels=c("a.", "b."))
ggsave("figs/pub/pref_temp_by_ZT.png", width=4, height=7, units="in")




# effects -----------------------------------------------------------------

sum.t <- readRDS("models/noCorr_randEff/out_temperature_RS_Tilapia.rds") %>%
  as_draws_df(variable="b_.*:Group|b_Group", regex=T) %>% 
  rename_with(~str_replace(.x, "Icos6.283185MUZTD24", "cosZT")) %>% 
  rename_with(~str_replace(.x, "Isin6.283185MUZTD24", "sinZT")) %>% 
  rename_with(~str_remove(.x, "b_")) %>% 
  rename_with(~str_remove(.x, "Group")) %>%
  select(-starts_with(".")) %>%
  pivot_longer(everything(), names_to="param", values_to="draws") %>%
  group_by(param) %>%
  summarise(mn=mean(draws),
            md=median(draws),
            CI80_l=quantile(draws, 0.1),
            CI80_h=quantile(draws, 0.9),
            CI90_l=quantile(draws, 0.05),
            CI90_h=quantile(draws, 0.95),
            CI95_l=quantile(draws, 0.025),
            CI95_h=quantile(draws, 0.975)) %>%
  ungroup %>% 
  mutate(species="Nile tilapia")

sum.z <- readRDS("models/noCorr_randEff/out_temperature_RS_ZF.rds") %>%
  as_draws_df(variable="b_.*:Group|b_Group", regex=T) %>% 
  rename_with(~str_replace(.x, "Icos6.283185MUZTD24", "cosZT")) %>% 
  rename_with(~str_replace(.x, "Isin6.283185MUZTD24", "sinZT")) %>% 
  rename_with(~str_remove(.x, "b_")) %>% 
  rename_with(~str_remove(.x, "Group")) %>%
  select(-starts_with(".")) %>%
  pivot_longer(everything(), names_to="param", values_to="draws") %>%
  group_by(param) %>%
  summarise(mn=mean(draws),
            md=median(draws),
            CI80_l=quantile(draws, 0.1),
            CI80_h=quantile(draws, 0.9),
            CI90_l=quantile(draws, 0.05),
            CI90_h=quantile(draws, 0.95),
            CI95_l=quantile(draws, 0.025),
            CI95_h=quantile(draws, 0.975)) %>%
  ungroup %>% 
  mutate(species="Zebrafish")

sum.df <- bind_rows(sum.t, sum.z) %>%
  mutate(param=factor(param, 
                      levels=c("Acclimation", "Experiment", 
                               "sinZT:Acclimation", "sinZT:Experiment", 
                               "cosZT:Acclimation", "cosZT:Experiment",
                               "cosZT:sinZT:Acclimation", "cosZT:sinZT:Experiment")),
         Group=if_else(grepl("Accl", param), "Acclimation", "Experiment"),
         timeEff=str_remove(str_remove(param, ":Acclimation"), ":Experiment"),
         timeEff=str_replace(timeEff, "Experiment|Acclimation", "Intercept")) %>%
  mutate(Group=factor(Group, levels=c("Experiment", "Acclimation")),
         timeEff=factor(timeEff, levels=c("cosZT:sinZT", "sinZT", "cosZT", "Intercept"),
                        labels=c("cos(ZT):sin(ZT)", "sin(ZT)\n[6-18h cycle]",
                                 "cos(ZT)\n[0-12h cycle]", "Intercept\n[mean temp.]"))) %>%
  mutate(species=factor(species, levels=c("Zebrafish", "Nile tilapia")))
ggplot(sum.df, aes(y=timeEff, x=md, group=Group, colour=Group)) +
  geom_vline(xintercept=0, colour="grey", size=0.5) +
  geom_point(size=2, position=position_dodge(width=0.4)) +
  geom_linerange(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, 
                 position=position_dodge(width=0.4)) +
  geom_linerange(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, 
                 position=position_dodge(width=0.4)) +
  geom_linerange(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, 
                 position=position_dodge(width=0.4)) +
  theme_classic() + 
  scale_colour_manual(values=group_col[2:3]) +
  facet_wrap(~species) +
  labs(x="Effect size") +
  theme(axis.title.y=element_blank(),
        legend.position=c(0.9, 0.13),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key.height=unit(0.4, "cm"),
        legend.key.width=unit(0.3, "cm"))
ggsave("figs/pub/effects.png", width=6, height=3, dpi=300)








# heatmaps ----------------------------------------------------------------
pred.ls <- map2(data.c, out.c,
               ~.x %>% filter(Tank==1) %>%
                 mutate(pred.mn=colMeans(pmax(
                   exp(posterior_predict(.y, 
                                       newdata=.x %>% filter(Tank==1),
                                       re_formula=NA))-1, 
                   0))) %>%
                 mutate(Group=factor(Group, 
                                     labels=c("Control CTE", "Acclimation", "Experiment")),
                        Temp=Temp + rnorm(n(), 0, 0.2)) %>%
                 group_by(Group) %>%
                 group_split())
nBins <- 20
temp.heat.ZF <- map(pred.ls[[2]], 
                    ~.x %>%
                      select(pred.mn, ZT, Group, Temp) %>% 
                      mutate(FishCount=round(pred.mn)) %>% 
                      uncount(FishCount) %>% 
                      ggplot(aes(ZT, Temp)) +
                      geom_density2d_filled(colour=NA, bins=nBins) + 
                      scale_fill_manual(values=lisa::lisa_palette("KatsushikaHokusai", nBins, "continuous"),
                                        guide="none") +
                      scale_x_continuous("ZT", breaks=seq(0, 24, by=6)) +
                      scale_y_continuous(limits=temp_rng[[2]],
                                         breaks=seq(ceiling(temp_rng[[2]][1]),
                                                    floor(temp_rng[[2]][2]),
                                                    by=2)) +
                      facet_wrap(~Group) + 
                      labs(y=""))
temp.heat.ZF[[1]] <- temp.heat.ZF[[1]] + labs(y="Temperature (ºC)")
temp.heat.T <- map(pred.ls[[1]], 
                    ~.x %>%
                      select(pred.mn, ZT, Group, Temp) %>% 
                      mutate(FishCount=round(pred.mn)) %>% 
                      uncount(FishCount) %>% 
                      ggplot(aes(ZT, Temp)) +
                      geom_density2d_filled(colour=NA, bins=nBins) + 
                      scale_fill_manual(values=lisa::lisa_palette("KatsushikaHokusai", nBins, "continuous"),
                                        guide="none") +
                     scale_x_continuous("ZT", breaks=seq(0, 24, by=6)) +
                     scale_y_continuous(limits=temp_rng[[1]],
                                        breaks=seq(ceiling(temp_rng[[1]][1]),
                                                   floor(temp_rng[[1]][2]),
                                                   by=2)) +
                      facet_wrap(~Group) + 
                      labs(y=""))
temp.heat.T[[1]] <- temp.heat.T[[1]] + labs(y="Temperature (ºC)")

ggpubr::ggarrange(plotlist=c(temp.heat.ZF, temp.heat.T), nrow=2, ncol=3, 
                  labels=c("a.", "", "", "b.", "", ""))
ggsave("figs/pub/heatmap_fitted.png", width=7, height=5, units="in")








# tube plot ---------------------------------------------------------------

tube.base <- map(data.c, 
               ~expand_grid(ZT=unique(.x$ZT),
                            Chamber=unique(.x$Chamber),
                            Group=unique(.x$Group)))
tube.df <- map2_dfr(tube.base, out.c,
                ~.x %>%
                  mutate(pred.mn=apply(posterior_epred(.y, newdata=.x, re.form=NA), 
                                       2, mean)),
                .id="Species") %>%
  mutate(Group=factor(Group, 
                      labels=c("Control CTE", "Acclimation", "Experiment")),
         Species=factor(Species, levels=rev(names(species))))
fig.tube <- tube.df %>%
      group_by(Species, Chamber, ZT, Group) %>%
      summarise(mn=mean(exp(pred.mn)-1)) %>%
      ggplot(aes(as.numeric(Chamber), mn, colour=ZT, group=ZT)) + 
      geom_line() +
      facet_grid(Species~Group) +
      scale_colour_gradientn(colours=cmr.ls$seasons, 
                             breaks=c(0,6, 12, 18, 24), limits=c(0, 24)) +
      ylim(0, NA) +
      labs(x="Chamber", y="Number of fish\n(posterior mean)") +
      theme(legend.key.height=unit(0.25, "cm"),
            legend.position="bottom",
            panel.border=element_rect(colour="grey10", size=0.5, fill=NA)) 
ggsave("figs/pub/tube_plot.png", width=5, height=4, units="in")






# acrophase ---------------------------------------------------------------

pred.ZTmax <- map2(out.c, tube.base, 
               ~posterior_epred(.x, newdata=.y, re.form=NA) %>%
                 t %>% as_tibble() %>%
                 cbind(.y) %>%
                 group_by(Chamber, Group) %>%
                 arrange(ZT) %>%
                 mutate(across(starts_with("V"), which.max)) %>%
                 select(-ZT) %>%
                 pivot_longer(starts_with("V"), names_to="iter", values_to="ZT_max"))

map_dfr(pred.ZTmax, 
        ~.x %>% 
          group_by(Chamber, Group, ZT_max) %>%
          summarise(N=n()) %>%
          group_by(Chamber, Group) %>%
          mutate(pr=N/sum(N)),
        .id="Species") %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT_max-1, y=pr, fill=Chamber)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=0.58, ymin=0.55), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=0.58, ymin=0.55), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) + 
  scale_fill_brewer(type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Acrophase probability")
ggsave("figs/pub/radial_acrophase_distr_fitted.png", width=7, height=4)

map_dfr(data.c,
        ~.x%>%
          group_by(Group, Chamber, Tank, Days) %>%
          arrange(ZT) %>%
          summarise(ZT_max=which.max(FishCount)) %>%
          group_by(Group, Chamber, ZT_max) %>%
          summarise(N=n()) %>%
          group_by(Group, Chamber) %>%
          mutate(pr=N/sum(N)),
        .id="Species") %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT_max-1, y=pr, fill=Chamber)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=0.58, ymin=0.55), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=0.58, ymin=0.55), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) + 
  scale_fill_brewer(type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Acrophase probability")
ggsave("figs/pub/radial_acrophase_distr_obs.png", width=7, height=4)

pred.CHmax <- map2(out.c, tube.base, 
                   ~posterior_epred(.x, newdata=.y, re.form=NA) %>%
                     t %>% as_tibble() %>%
                     cbind(.y) %>%
                     group_by(ZT, Group) %>%
                     arrange(Chamber) %>%
                     mutate(across(starts_with("V"), which.max)) %>%
                     select(-Chamber) %>%
                     pivot_longer(starts_with("V"), names_to="iter", values_to="Chamber_max"))

map_dfr(pred.CHmax, 
        ~.x %>% 
          group_by(ZT, Group, Chamber_max) %>%
          summarise(N=n()) %>%
          group_by(ZT, Group) %>%
          mutate(pr=N/sum(N)),
        .id="Species") %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Chamber_max=factor(Chamber_max),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT, y=pr, fill=Chamber_max)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) +
  scale_fill_brewer("Chamber", type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Probability of highest fish density")
ggsave("figs/pub/radial_chamber_distr_fitted.png", width=7, height=4)

map_dfr(data.c,
        ~.x%>%
          group_by(ZT, Group, Tank, Days) %>%
          arrange(Chamber) %>%
          summarise(Chamber_max=which.max(FishCount)) %>%
          group_by(ZT, Group, Chamber_max) %>%
          summarise(N=n()) %>%
          group_by(ZT, Group) %>%
          mutate(pr=N/sum(N)),
        .id="Species") %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Chamber_max=factor(Chamber_max),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT, y=pr, fill=Chamber_max)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) + 
  scale_fill_brewer("Chamber", type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Probability of highest fish density")
ggsave("figs/pub/radial_chamber_distr_obs.png", width=7, height=4)


map_dfr(data.c,
        ~.x%>%
          group_by(ZT, Group, Chamber) %>%
          summarise(meanFish=mean(FishCount, na.rm=T)) %>%
          group_by(ZT, Group) %>%
          mutate(pr=meanFish/sum(meanFish)),
        .id="Species") %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Chamber=factor(Chamber),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT, y=pr, fill=Chamber)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) + 
  scale_fill_brewer("Chamber", type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Proportion of fish (observed)")
ggsave("figs/pub/radial_distr_obs.png", width=7, height=4)


tube.df %>%
  group_by(Species, ZT, Group) %>%
  mutate(pr=pred.mn/sum(pred.mn)) %>%
  ungroup %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Chamber=factor(Chamber),
         Group=factor(Group, labels=c("Control CTE", "Acclimation", "Experiment"))) %>%
  ggplot(aes(ZT, y=pr, fill=Chamber)) + 
  geom_bar(stat="identity", colour="grey30", size=0.1) +
  geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=11.5, xmax=23.5, ymax=1.1, ymin=1.05), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) + 
  scale_fill_brewer("Chamber", type="div", palette="RdBu", direction=-1) +
  facet_grid(Species~Group) +   
  coord_polar(theta="x", start=(-0.5*pi*2)/24) +
  theme(panel.grid.major=element_line(size=0.1, colour="grey")) +
  labs(x="ZT (h)", y="Distribution of fish (fitted)")
ggsave("figs/pub/radial_distr_fit.png", width=7, height=4)


# comparisons -------------------------------------------------------------

# Mean temp
map(out, ~hypothesis(.x, "GroupAcclimation = GroupExperiment"))
map(out, ~hypothesis(.x, "GroupAcclimation = 0"))
map(out, ~hypothesis(.x, "GroupExperiment = 0"))

hyp_b <- c(cos_Ctrl="Icos6.283185MUZTD24 = 0",
           cos_Acc="Icos6.283185MUZTD24:GroupAcclimation = 0",
           cos_Exp="Icos6.283185MUZTD24:GroupExperiment = 0",
           sin_Ctrl="Isin6.283185MUZTD24 = 0",
           sin_Acc="Isin6.283185MUZTD24:GroupAcclimation = 0",
           sin_Exp="Isin6.283185MUZTD24:GroupExperiment = 0")

hypothesis(out[[1]], hyp_b)
hypothesis(out[[2]], hyp_b)


hyp_b <- c(cos_Ctrl="Icos6.283185MUZTD24:Chamber1 = 0",
           cos_Acc="Icos6.283185MUZTD24:GroupAcclimation = 0",
           cos_Exp="Icos6.283185MUZTD24:GroupExperiment = 0",
           sin_Ctrl="Isin6.283185MUZTD24 = 0",
           sin_Acc="Isin6.283185MUZTD24:GroupAcclimation = 0",
           sin_Exp="Isin6.283185MUZTD24:GroupExperiment = 0")
hypothesis(out.c[[1]], hyp_b)
hypothesis(out.c[[2]], hyp_b)

