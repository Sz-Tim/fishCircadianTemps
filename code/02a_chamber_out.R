# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[2]
mod_type <- c("RI", "RS")[2]



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
library(lisa)
theme_set(theme_classic())
group_col <- c(Control="#B68E52", Acclimation="#697A55", Experiment="#8399B3")
cmr.ls <- readRDS("figs/cmr_cmaps.RDS")
temp_rng <- list(ZF=c(24, 31.5), Tilapia=c(25.5, 34.5))[[species]]

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

wtChamber.df <- data.df %>% group_by(ZT, Days, Tank, Group) %>%
  mutate(Chamber=as.numeric(Chamber)) %>%
  summarise(wtChamber=sum(Chamber*FishCount, na.rm=T)/sum(FishCount, na.rm=T)) %>%
  ungroup
wtChamber.sum <- wtChamber.df %>%
  group_by(ZT, Group, Tank) %>%
  summarise(mn=mean(wtChamber, na.rm=T),
            se=sd(wtChamber, na.rm=T)/sqrt(sum(!is.na(wtChamber))))

data.noNA <- data.df %>% filter(complete.cases(.))

out <- readRDS(glue("models/out_count_{mod_type}_{species}.rds"))



# checks ------------------------------------------------------------------

# pp_check(out)
# summary(out)
# conditional_effects(out, conditions=tibble(Chamber=unique(data.df$Chamber)))
# conditional_effects(out, conditions=tibble(Group=unique(data.df$Group)))





# patterns ----------------------------------------------------------------
pred.df <- data.df %>% filter(Tank==1)
pred.lFish <- posterior_epred(out, newdata=pred.df, re_formula=NA)
pred.Fish <- pmax(exp(pred.lFish)-1, 0)


# preferred temperature
prefTemp.df_epred <- t(pred.Fish) %>% as_tibble() %>%
  mutate(id=row_number()) %>%
  pivot_longer(1:4000, names_to="iter", values_to="predFish") %>%
  full_join(pred.df %>% mutate(id=row_number()), by="id") %>%
  group_by(ZT, Group, Days, iter) %>%
  summarise(wtTemp=sum(Temp*predFish)/sum(predFish))
prefTemp.sum_epred <- prefTemp.df_epred %>%
  group_by(ZT, Group) %>%
  summarise(pred.mn=mean(wtTemp, na.rm=T),
            CI95_l=quantile(wtTemp, 0.025, na.rm=T),
            CI90_l=quantile(wtTemp, 0.05, na.rm=T),
            CI80_l=quantile(wtTemp, 0.1, na.rm=T),
            CI80_h=quantile(wtTemp, 0.9, na.rm=T),
            CI90_h=quantile(wtTemp, 0.95, na.rm=T),
            CI95_h=quantile(wtTemp, 0.975, na.rm=T))
prefTemp.sum_epred %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_point(aes(y=pred.mn), shape=1, size=1.75,
             position=position_dodge(width=0.75)) + 
  geom_linerange(aes(ymin=CI80_l, ymax=CI80_h), size=1,
                 position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=CI90_l, ymax=CI90_h), size=0.5,
                 position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=CI95_l, ymax=CI95_h), size=0.25, width=0.2,
                position=position_dodge(width=0.75)) +
  labs(x="ZT", y="Preferred temperature\n(mean + 80%, 90%, 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  ggtitle(species) +
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank(), 
        legend.background=element_blank())
ggsave(paste0("figs/preferred_temp_from_nFish_", species, "_", mod_type,".png"), 
       height=5, width=7, dpi=300)


# heatmap
data.df %>% 
  group_by(ZT, Group, Chamber) %>%
  summarise(mnFish=mean(FishCount, na.rm=T)) %>%
  ggplot(aes(ZT, Chamber, colour=mnFish, fill=mnFish)) + 
  geom_raster(interpolate=F) +
  facet_wrap(~Group) +
  scale_fill_gradientn("Mean observed count", 
                       colours=lisa_palette("KatsushikaHokusai", 1e3, "continuous")) +
  guides(fill=guide_colourbar(title.position="top", title.hjust=0.5)) +
  theme(legend.position="bottom",
        legend.key.width=unit(1, "cm"), 
        legend.background=element_blank()) + 
  ggtitle(species)
ggsave(glue("figs/heatmap_observed_mean_{species}.png"), 
       height=4, width=8, dpi=300)

pred.ls <- pred.df %>% 
  filter(!is.na(FishCount)) %>% 
  group_by(Group) %>%
  group_split()
temp.heat <- map(pred.ls, ~.x %>%
                   select(FishCount, ZT, Group, Temp) %>% 
                   mutate(FishCount=round(FishCount)) %>% 
                   uncount(FishCount) %>% 
                   ggplot(aes(ZT, Temp)) +
                   geom_density2d_filled(colour=NA, bins=20) + 
                   scale_fill_manual(values=lisa_palette("KatsushikaHokusai", 20, "continuous"),
                                     guide="none") +
                   ylim(temp_rng[1], temp_rng[2]) + 
                   facet_wrap(~Group) + 
                   labs(y="", title=" "))
temp.heat[[1]] <- temp.heat[[1]] + 
  labs(title=paste0(species, ": Observed"), y="Temperature")
ggpubr::ggarrange(plotlist=temp.heat, nrow=1)
ggsave(glue("figs/heatmap_observed_temp_{species}.png"), 
       height=3, width=8, dpi=300)

pred.ls <- pred.df %>% 
  mutate(pred.mn=colMeans(pred.Fish)) %>%
  group_by(Group) %>%
  group_split()
temp.heat <- map(pred.ls, ~.x %>%
                   select(pred.mn, ZT, Group, Temp) %>% 
                   mutate(FishCount=round(pred.mn)) %>% 
                   uncount(FishCount) %>% 
                   ggplot(aes(ZT, Temp)) +
                   geom_density2d_filled(colour=NA, bins=20) + 
                   scale_fill_manual(values=lisa_palette("KatsushikaHokusai", 20, "continuous"),
                                     guide="none") +
                   ylim(temp_rng[1], temp_rng[2]) + 
                   facet_wrap(~Group) + 
                   labs(y="", title=" "))
temp.heat[[1]] <- temp.heat[[1]] + 
  labs(title=paste0(species, ": Fitted"), y="Temperature")
ggpubr::ggarrange(plotlist=temp.heat, nrow=1)
ggsave(glue("figs/heatmap_predicted_temp_{species}.png"), 
       height=3, width=8, dpi=300)

heatmap.df <- expand_grid(ZT=unique(data.df$ZT),
                          Chamber=unique(data.df$Chamber),
                          Group=unique(data.df$Group)) 
heatmap.pred <- posterior_epred(out, newdata=heatmap.df, re.form=NA)
heatmap.df <- heatmap.df %>%
  mutate(pred.mn=apply(pmax(exp(heatmap.pred)-1, 0), 2, mean),
         pred.lmn=apply(heatmap.pred, 2, mean))

ggplot(heatmap.df, aes(ZT, Chamber, fill=exp(pred.lmn)-1)) + 
  geom_raster(interpolate=F) +
  facet_wrap(~Group) +
  scale_fill_gradientn("Mean predicted count", 
                       colours=lisa_palette("KatsushikaHokusai", 1e3, "continuous")) +
  guides(fill=guide_colourbar(title.position="top", title.hjust=0.5)) +
  theme(legend.position="bottom",
        legend.key.width=unit(1, "cm"), 
        legend.background=element_blank()) + 
  ggtitle(species)
ggsave(paste0("figs/heatmap_posterior_mean_", species, "_", mod_type,".png"), 
       height=4, width=8, dpi=300)



# preferred chamber
prefChamber.df_epred <- t(heatmap.pred) %>% as_tibble() %>%
  mutate(id=row_number()) %>%
  pivot_longer(1:4000, names_to="iter", values_to="predFish") %>%
  mutate(predFish=pmax(exp(predFish)-1, 0)) %>%
  full_join(heatmap.df %>% mutate(id=row_number()), by="id") %>%
  group_by(ZT, Group, iter) %>%
  mutate(Chamber=as.numeric(Chamber)) %>%
  summarise(wtChamber=sum(Chamber*predFish)/sum(predFish)) %>%
  group_by(ZT, Group) %>%
  summarise(pred.mn=mean(wtChamber, na.rm=T),
            CI95_l=quantile(wtChamber, 0.025, na.rm=T),
            CI90_l=quantile(wtChamber, 0.05, na.rm=T),
            CI80_l=quantile(wtChamber, 0.1, na.rm=T),
            CI80_h=quantile(wtChamber, 0.9, na.rm=T),
            CI90_h=quantile(wtChamber, 0.95, na.rm=T),
            CI95_h=quantile(wtChamber, 0.975, na.rm=T))
write_csv(prefChamber.df_epred, glue("out/predicted_prefChamber_{species}.csv"))

prefChamber.df_epred %>%
  ggplot(aes(ZT, colour=Group)) + 
  geom_point(aes(y=pred.mn), shape=1, size=1.75,
             position=position_dodge(width=0.75)) + 
  geom_linerange(aes(ymin=CI80_l, ymax=CI80_h), size=1,
                 position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=CI90_l, ymax=CI90_h), size=0.5,
                 position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=CI95_l, ymax=CI95_h), size=0.25, width=0.2,
                position=position_dodge(width=0.75)) +
  labs(x="ZT", y="Preferred chamber\n(mean + 80%, 90%, 95% CIs)") +
  ylim(1, 5) +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  theme(legend.position=c(0.85, 0.15), 
        legend.title=element_blank(), 
        legend.background=element_blank()) +
  ggtitle(species)
ggsave(paste0("figs/preferred_chamber_mnOnly_", species, "_", mod_type, ".png"), 
       height=5, width=7, dpi=300)

prefChamber.df_epred %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_jitter(data=wtChamber.df, aes(y=wtChamber, shape=Group),
              alpha=0.8, size=0.95) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  ylim(1, 5) +
  labs(x="ZT", y="Preferred chamber\n(mean + 95% CIs") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  theme(legend.position=c(0.85, 0.15), 
        legend.title=element_blank(), 
        legend.background=element_blank()) +
  ggtitle(species)
ggsave(paste0("figs/preferred_chamber_", species, "_", mod_type,".png"), 
       height=5, width=7, dpi=300)

prefChamber.df_epred %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_point(data=wtChamber.sum, aes(y=mn, shape=Group), alpha=0.8, size=0.95, 
             position=position_dodge(width=0.45)) +
  geom_errorbar(data=wtChamber.sum, aes(ymin=mn-se, ymax=mn+se), size=0.25,
                width=0.3, position=position_dodge(width=0.45)) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  ylim(1, 5) +
  labs(x="ZT", y="Preferred chamber\n(mean + 95% CIs") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  theme(legend.position=c(0.85, 0.15), 
        legend.title=element_blank(), 
        legend.background=element_blank()) +
  ggtitle(species)
ggsave(paste0("figs/preferred_chamber_tankMeans_", species, "_", mod_type,".png"), 
       height=5, width=7, dpi=300)

heatmap.df %>% 
  group_by(Chamber, ZT, Group) %>%
  summarise(mn=mean(exp(pred.lmn)-1)) %>%
  ggplot(aes(as.numeric(Chamber), mn, colour=ZT, group=ZT)) + 
  geom_line() +
  facet_wrap(~Group) +
  scale_colour_gradientn(colours=cmr.ls$seasons) +
  ylim(0, NA) +
  labs(title=species, x="Chamber", y="Number of fish (fitted mean)")
ggsave(paste0("figs/tube_plots_", species, "_", mod_type,".png"), 
       width=8, height=4, dpi=300)

heatmap.df %>% 
  group_by(Chamber, ZT, Group) %>%
  summarise(mn=mean(exp(pred.lmn)-1)) %>%
  ggplot(aes(as.numeric(Chamber), mn, colour=ZT, group=ZT)) + 
  geom_jitter(data=data.df, aes(y=FishCount), alpha=0.5, shape=1, size=0.5) + 
  geom_line() +
  facet_wrap(~Group) +
  scale_colour_gradientn(colours=cmr.ls$seasons) +
  ylim(0, NA) +
  labs(title=species, x="Chamber", y="Number of fish (fitted mean)")
ggsave(paste0("figs/tube_plots_", species, "_", mod_type,"_pts.png"), 
       width=8, height=4, dpi=300)
