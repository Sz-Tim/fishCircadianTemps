# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk





# switches ----------------------------------------------------------------

species <- c("ZF", "Tilapia")[1]
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

out <- readRDS(glue("models/out_temperature_{mod_type}_{species}.rds"))



# checks ------------------------------------------------------------------

# pp_check(out)
# summary(out)
# conditional_effects(out, conditions=tibble(Group=unique(data.df$Group)))




# patterns ----------------------------------------------------------------

pred.df <- expand_grid(ZT=unique(data.df$ZT),
                       Group=unique(data.df$Group)) 
preds <- posterior_epred(out, newdata=pred.df, re.form=NA)
pred.df <- pred.df %>%
  mutate(pred.mn=apply(preds, 2, median),
         CI95_l=apply(preds, 2, function(x) quantile(x, 0.025)),
         CI90_l=apply(preds, 2, function(x) quantile(x, 0.05)),
         CI80_l=apply(preds, 2, function(x) quantile(x, 0.1)),
         CI80_h=apply(preds, 2, function(x) quantile(x, 0.9)),
         CI90_h=apply(preds, 2, function(x) quantile(x, 0.95)),
         CI95_h=apply(preds, 2, function(x) quantile(x, 0.975)))

pred.df %>%
  ggplot(aes(ZT)) + 
  geom_point(aes(y=pred.mn), size=2) + 
  geom_linerange(aes(ymin=CI80_l, ymax=CI80_h), size=1.5) +
  geom_linerange(aes(ymin=CI90_l, ymax=CI90_h), size=1) +
  geom_linerange(aes(ymin=CI95_l, ymax=CI95_h), size=0.5) +
  facet_wrap(~Group) +
  labs(x="ZT", y="Preferred temperature\n(mean + 80%, 90%, 95% CIs)")

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
  labs(x="ZT", y="Preferred temperature\n(mean + 80%, 90%, 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  ggtitle(species) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank())
ggsave(paste0("figs/preferred_temperature_mnOnly_", species, "_", mod_type, ".png"), 
       width=5, height=5, dpi=300)

pred.df %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_jitter(data=data.df, aes(y=prefTemp, shape=Group),
              alpha=0.8, size=0.95) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  labs(x="ZT", y="Preferred temperature\n(mean + 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  ggtitle(species) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank())
ggsave(paste0("figs/preferred_temperature_", species, "_", mod_type, ".png"), 
       width=5, height=5, dpi=300)

pred.df %>%
  ggplot(aes(ZT, colour=Group, fill=Group)) + 
  geom_point(data=data.sum, aes(y=mn, shape=Group), alpha=0.8, size=0.95, 
             position=position_dodge(width=0.45)) +
  geom_errorbar(data=data.sum, aes(ymin=mn-2*se, ymax=mn+2*se), size=0.25,
                width=0.3, position=position_dodge(width=0.45)) +
  geom_line(aes(y=pred.mn)) + 
  geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h), colour=NA, alpha=0.5) +
  labs(x="ZT", y="Preferred temperature\n(mean + 95% CIs)") +
  scale_fill_manual(values=group_col) +
  scale_colour_manual(values=group_col) + 
  scale_shape_manual(values=1:3) +
  ggtitle(species) + 
  ylim(temp_rng[1], temp_rng[2]) +
  theme(legend.position=c(0.8, 0.15),
        legend.title=element_blank())
ggsave(paste0("figs/preferred_temperature_tankMeans_", species, "_", mod_type, ".png"), 
       width=5, height=5, dpi=300)
