# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("Zebra fish"="ZF", "Tilapia"="Tilapia")[2]
mod_type <- c("RI", "RS")[2]



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
library(lisa)
theme_set(theme_classic())
group_col <- c(Control="#B68E52", Acclimation="#697A55", Experiment="#8399B3")
group_col <- lisa_palette("MarcChagall", 5, "discrete")[c(1,2,3)] %>%
  setNames(c("Control", "Acclimation", "Experiment"))
cmr.ls <- readRDS("figs/cmr_cmaps.RDS")
temp_rng <- list(ZF=c(24, 31.5), Tilapia=c(25.5, 34.5))[[species]]

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

data.sum <- data.df %>%
  group_by(ZT, Group, Tank) %>%
  summarise(mn=mean(prefTemp, na.rm=T),
            se=sd(prefTemp, na.rm=T)/sqrt(sum(!is.na(prefTemp))))

out <- readRDS(glue("models/noCorr_randEff/out_temperature_{mod_type}_{species}.rds"))



# checks ------------------------------------------------------------------

# pp_check(out)
# summary(out)
# conditional_effects(out, conditions=tibble(Group=unique(data.df$Group)))




# patterns ----------------------------------------------------------------

pred.df <- expand_grid(ZT=unique(data.df$ZT),
                       Group=unique(data.df$Group)) 
preds <- posterior_epred(out, newdata=pred.df, re.form=NA)
pred.df <- pred.df %>%
  mutate(pred.mn=apply(preds, 2, mean),
         CI95_l=apply(preds, 2, function(x) quantile(x, 0.025)),
         CI90_l=apply(preds, 2, function(x) quantile(x, 0.05)),
         CI80_l=apply(preds, 2, function(x) quantile(x, 0.1)),
         CI80_h=apply(preds, 2, function(x) quantile(x, 0.9)),
         CI90_h=apply(preds, 2, function(x) quantile(x, 0.95)),
         CI95_h=apply(preds, 2, function(x) quantile(x, 0.975)))

write_csv(pred.df, glue("out/tempMod/predicted_prefTemp_{species}.csv"))

pred.df %>%
  ggplot(aes(ZT, colour=Group)) + 
  geom_point(aes(y=pred.mn), shape=1, size=1.75,
             position=position_dodge(width=0.75)) + 
  geom_linerange(aes(ymin=CI80_l, ymax=CI80_h), size=1,
                 position=position_dodge(width=0.75)) +
  geom_linerange(aes(ymin=CI90_l, ymax=CI90_h), size=0.5,
                 position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=CI95_l, ymax=CI95_h), size=0.25, width=0.2,
                position=position_dodge(width=0.75)) +
  geom_rect(aes(xmin=0, xmax=12, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="grey30", size=0.25) +
  labs(x="ZT", y="Preferred temperature\n(mean + 80%, 90%, 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  ggtitle(names(species)) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank(), 
        legend.background=element_blank())
ggsave(paste0("figs/tempMod/preferred_temperature_mnOnly_", species, "_", mod_type, ".png"), 
       height=5, width=7, dpi=300)

pred.df %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_jitter(data=data.df, aes(y=prefTemp, shape=Group),
              alpha=0.8, size=0.95) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  geom_rect(aes(xmin=0, xmax=12, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="grey30", size=0.25) +
  labs(x="ZT", y="Preferred temperature\n(mean + 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  ggtitle(names(species)) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank(), 
        legend.background=element_blank())
ggsave(paste0("figs/tempMod/preferred_temperature_", species, "_", mod_type, ".png"), 
       height=5, width=7, dpi=300)

pred.df %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_point(data=data.sum, aes(y=mn, shape=Group), alpha=0.8, size=0.95, 
             position=position_dodge(width=0.45)) +
  geom_errorbar(data=data.sum, aes(ymin=mn-2*se, ymax=mn+2*se), size=0.25,
                width=0.3, position=position_dodge(width=0.45)) +
  geom_rect(aes(xmin=0, xmax=12, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymin=temp_rng[2]-0.1, ymax=temp_rng[2]), 
            colour="grey30", fill="grey30", size=0.25) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  labs(x="ZT", y="Preferred temperature\n(mean + 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  ggtitle(names(species)) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank(), 
        legend.background=element_blank())
ggsave(paste0("figs/tempMod/preferred_temperature_tankMeans_", species, "_", mod_type, ".png"), 
       height=5, width=7, dpi=300)





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
  mutate(species="Tilapia")

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
  mutate(species="Zebra fish")

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
                                 "cos(ZT)\n[0-12h cycle]", "Intercept\n[mean temp.]")))
ggplot(sum.df, aes(y=timeEff, x=md, group=Group, colour=Group)) +
  geom_vline(xintercept=0, colour="grey", size=0.5) +
  geom_point(size=2, position=position_dodge(width=0.25)) +
  geom_linerange(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, 
                 position=position_dodge(width=0.25)) +
  geom_linerange(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, 
                 position=position_dodge(width=0.25)) +
  geom_linerange(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, 
                 position=position_dodge(width=0.25)) +
  theme_classic() + 
  scale_colour_manual(values=group_col[2:3]) +
  facet_wrap(~species) +
  labs(x="Effect size (median, 80%, 90%, 95% CIs)",
       title="Effects relative to control") +
  theme(axis.title.y=element_blank(),
        legend.position=c(0.92, 0.15),
        legend.title=element_blank(),
        legend.background=element_blank())
ggsave("figs/tempMod/effects.png", width=7, height=4, dpi=300)


sum.df %>% select(species, Group, timeEff, mn, md, starts_with("CI")) %>%
  arrange(species, desc(timeEff), desc(Group)) %>%
  rename(Effect=timeEff, Mean=mn, Median=md) %>%
  write_csv("out/tempMod/effects.csv")


# 95% Confidence:
# T Acc Int
# T Acc cosZT
# Z Acc Int
# Z Exp Int
# Z Acc sinZT

# 90% Confidence:
# T Acc Int
# T Acc cosZT
# T Exp cosZT
# Z Acc Int
# Z Exp Int
# Z Acc sinZT
# Z Exp sinZT


hyp1 <- hypothesis(out, c("Icos6.283185MUZTD24:GroupExperiment = 0",
                          "Icos6.283185MUZTD24:GroupAcclimation = 0",
                          "Isin6.283185MUZTD24:GroupExperiment = 0",
                          "Isin6.283185MUZTD24:GroupAcclimation = 0",
                          "GroupAcclimation = GroupExperiment"))
hypothesis(out, "GroupAcclimation = GroupExperiment", group="Tank", scope="coef")
hypothesis(out, "GroupAcclimation > GroupExperiment")





# Heatmap
tempPost.df <- posterior_predict(out, newdata=pred.df, re.form=NA, ndraws=4e3) %>% 
  as_draws_df() %>% 
  pivot_longer(starts_with("..."), names_to="pred_row", values_to="post") %>%
  mutate(pred_row=as.numeric(str_remove(pred_row, "..."))) %>%
  full_join(., pred.df %>% select(ZT, Group) %>% mutate(pred_row=row_number()))

pred.ls <- tempPost.df %>% 
  group_by(Group) %>%
  group_split()
temp.heat <- map(pred.ls, ~.x %>%
                   ggplot(aes(ZT, post)) +
                   geom_density2d_filled(colour=NA, bins=20) + 
                   scale_fill_manual(values=lisa_palette("KatsushikaHokusai", 20, "continuous"),
                                     guide="none") +
                   ylim(temp_rng[1], temp_rng[2]) + 
                   facet_wrap(~Group) + 
                   labs(y="", title=" "))
temp.heat[[1]] <- temp.heat[[1]] + 
  labs(title=paste0(species, ": Fitted"), y="Temperature")
ggpubr::ggarrange(plotlist=temp.heat, nrow=1)
ggsave(glue("figs/heatmap_fitted_prefTemp_posterior_{species}.png"), 
       height=3, width=8, dpi=300)
