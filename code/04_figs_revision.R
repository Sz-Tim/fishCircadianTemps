# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk



# set up ------------------------------------------------------------------

library(tidyverse)
library(lisa)
library(glue)
theme_set(theme_classic())
# group_col <- lisa::lisa_palette("ClaudeMonet_1", 5)[c(4,5)] %>% 
group_col <- lisa_palette("GustavKlimt", 5)[c(2,4)] %>%
  setNames(c("Control CTE", "Experiment"))
# group_col <- c(`Control CTE`="#1b9e77", Experiment="#7570b3")
group_raw <- c("Control", "Acclimation", "Experiment")
chmb_col <- c("1"="#0571b0", "2"="#92c5de", "3"="grey", "4"="#f4a582", "5"="#ca0020")
cmr.ls <- readRDS("figs/cmr_cmaps.RDS")
sp_col <- c("Nile tilapia"="#1a3431", "Zebrafish"="#6283c8")
species <- c("Nile tilapia"="Tilapia", "Zebrafish"="ZF")
temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))
temp_ctrl <- c(Tilapia=30.1, ZF=26.4)

lower_ZT.df <- tibble(light_1=0:12,
                      light_2=c(0:12)+0.5,
                      dark_1=c(0:12)+0.5,
                      dark_2=1:13,
                      temp_1=min(unlist(temp_rng)),
                      temp_2=min(unlist(temp_rng))+0.2,
                      amp_1=0, amp_2=0.025,
                      acro_1=0, acro_2=0.4,
                      prop_1=-0.05, prop_2=-0.02,
                      ElapsedDays=0, pred=30, M=30, amplitude=0, acrophase=12, pr_mn=0)
lower_ZT.Grp <- bind_rows(
  tibble(Group="Experiment",
         light_1=0:12,
         light_2=c(0:12)+0.5,
         dark_1=c(0:12)+0.5,
         dark_2=1:13,
         temp_1=min(unlist(temp_rng)),
         temp_2=min(unlist(temp_rng))+0.2,
         amp_1=0, amp_2=0.025,
         acro_1=0, acro_2=0.4,
         prop_1=-0.05, prop_2=-0.02,
         ElapsedDays=0, pred=30, M=30, M_mn=0, amplitude=0, acrophase=12, pr_mn=0),
  tibble(Group="Control CTE",
         light_1=0:5,
         light_2=c(0:5)+0.5,
         dark_1=c(0:5)+0.5,
         dark_2=1:6,
         temp_1=min(unlist(temp_rng)),
         temp_2=min(unlist(temp_rng))+0.2,
         amp_1=0, amp_2=0.025,
         acro_1=0, acro_2=0.4,
         prop_1=-0.05, prop_2=-0.02,
         ElapsedDays=0, pred=30, M=30, M_mn=0, amplitude=0, acrophase=12, pr_mn=0))

data.df <- map(species, 
               ~readxl::read_xlsx(dir("data", glue("{.x}_.*RawData2"), full.names=T),
                                  1, col_types=c(rep("numeric", 6), "skip", "skip")) %>%
                 mutate(ln_FishCount=log(FishCount+1),
                        Chamber=factor(Chamber),
                        Group=factor(Group, levels=c(3, 1, 2), labels=group_raw),
                        Tank=factor(Tank)) %>%
                 left_join(read_csv(glue("data/temp_{.x}.csv")) %>%
                             pivot_longer(starts_with("chamber_"), 
                                          names_to="Chamber", values_to="Temp") %>%
                             mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                                    Group=factor(Group, levels=group_raw))) %>%
                 group_by(ZT, Group, Tank, Days) %>%
                 summarise(prefTemp=sum(FishCount*Temp)/(sum(FishCount))) %>%
                 ungroup %>%
                 mutate(ElapsedTime=ZT+24*(Days-1)+24*3*(Group=="Experiment"),
                        ElapsedDays=ElapsedTime/24))

fit.N <- map_dfr(species, ~read_csv(glue("out/GS_predGlobal_NChmbr_{.x}.csv"))) %>%
  mutate(Chamber=factor(Chamber),
         Group=factor(Group, levels=c("Control", "Experimental"), labels=c("Control CTE", "Experiment")),
         SpeciesFull=factor(Species, levels=species, labels=names(species)))
fit.edge <- map_dfr(species, ~read_csv(glue("out/GS_predEdge_NChmbr_{.x}.csv"))) %>%
  mutate(Group=factor(Group, levels=c("Control", "Experimental"), labels=c("Control CTE", "Experiment")),
         SpeciesFull=factor(Species, levels=species, labels=names(species))) %>%
  rename(pr_mn=edge_mn, pr_lo=edge_lo, pr_hi=edge_hi)
fit.cw <- map_dfr(species, ~read_csv(glue("out/GS_predColdWarm_NChmbr_{.x}.csv"))) %>%
  mutate(Chambers=factor(Chambers),
         Group=factor(Group, levels=c("Control", "Experimental"), labels=c("Control CTE", "Experiment")),
         SpeciesFull=factor(Species, levels=species, labels=names(species)))
fit.temp_global <- map_dfr(species, ~read_csv(glue("out/GS_predGlobal_{.x}.csv"))) %>%
  mutate(Spline=factor("Global", levels=c("Global", "Tank")))
fit.temp_tank <- map_dfr(species, ~read_csv(glue("out/GS_predTank_{.x}.csv"))) %>%
  mutate(Tank=factor(Tank))



# Chamber figs
p.N_prMn <- ggplot(fit.N, aes(ElapsedDays, pr_mn, fill=Chamber)) +
  geom_hline(yintercept=c(0,1), colour="grey70", size=0.1) +
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), colour="grey70", size=0.1) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_area(colour="grey30", size=0.2) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=1)) +
  scale_y_continuous("Proportion of fish (posterior mean)", 
                     limits=c(-0.05, 1), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Group, scales="free_x", space="free_x") +
  theme_classic() +
  theme(legend.position="bottom")
ggsave("figs/pub/1_revision1/propFishByChamber_means.png", p.N_prMn, width=7, height=4, dpi=300)

p.N_prMnCI <- fit.N %>%
  mutate(Chamber=paste("Chamber", Chamber)) %>%
  ggplot(aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                  group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey70", size=0.1) +
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), colour="grey70", size=0.1) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Proportion of fish\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom")
ggsave("figs/pub/1_revision1/propFishByChamber_means+CIs.png", p.N_prMnCI, width=9, height=4, dpi=300)


p.Edge_prMnCI <- fit.edge %>%
  mutate(Chamber="Edge chambers") %>%
  ggplot(aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
             group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey50", size=0.25) +
  geom_hline(yintercept=c(0.4), colour="grey50", size=0.25) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Pr(edge chamber)\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom", 
        axis.title.y=element_blank())
ggsave("figs/pub/1_revision1/prEdge_means+CIs.png", p.Edge_prMnCI, width=7, height=4, dpi=300)

p.ColdWarm_prMnCI <- fit.cw %>%
  mutate(Chambers=factor(Chambers, levels=c("cold", "warm"), labels=paste(c("Cold", "Warm"), "chambers"))) %>%
  ggplot(aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
             group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey50", size=0.25) +
  geom_hline(yintercept=c(0.4), colour="grey50", size=0.25) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Proportion of fish\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chambers) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom")
ggsave("figs/pub/1_revision1/prColdWarm_means+CIs.png", p.ColdWarm_MMnCI, width=7, height=5, dpi=300)

ggpubr::ggarrange(p.ColdWarm_prMnCI, p.Edge_prMnCI, ncol=2, legend="bottom", 
                  common.legend=T, widths=c(1, 0.52), labels="auto")
ggsave("figs/pub/1_revision1/prColdWarmEdge_pr_means+CIs.png", width=9, height=5.25, dpi=300)


p.Edge_MMnCI <- fit.edge %>%
  mutate(Chamber="Edge chambers") %>%
  ggplot(aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
             group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey50", size=0.25) +
  geom_hline(yintercept=c(0.4), colour="grey50", size=0.25) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Pr(edge chamber): MESOR\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom", 
        axis.title.y=element_blank())
ggsave("figs/pub/1_revision1/prEdge_M_means+CIs.png", p.Edge_MMnCI, width=7, height=4, dpi=300)

p.ColdWarm_MMnCI <- fit.cw %>%
  mutate(Chambers=factor(Chambers, levels=c("cold", "warm"), labels=paste(c("Cold", "Warm"), "chambers"))) %>%
  ggplot(aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
             group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey50", size=0.25) +
  geom_hline(yintercept=c(0.4), colour="grey50", size=0.25) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Proportion of fish: MESOR\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chambers) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom")
ggsave("figs/pub/1_revision1/prColdWarm_M_means+CIs.png", p.ColdWarm_MMnCI, width=7, height=5, dpi=300)

ggpubr::ggarrange(p.ColdWarm_MMnCI, p.Edge_MMnCI, ncol=2, legend="bottom", 
                  common.legend=T, widths=c(1, 0.52), labels="auto")
ggsave("figs/pub/1_revision1/prColdWarmEdge_M_means+CIs.png", width=9, height=5.25, dpi=300)





p.N_M <- fit.N %>%
  mutate(Chamber=paste("Chamber", Chamber)) %>%
  ggplot(aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                  group=Group, colour=Group, fill=Group)) +
  geom_hline(yintercept=c(0,1), colour="grey70", size=0.1) +
  geom_hline(yintercept=c(0.2,0.4,0.6,0.8), colour="grey70", size=0.1) +
  geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual("", values=group_col) +
  scale_fill_manual("", values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("MESOR: mean proportion of fish\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(SpeciesFull~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom", 
        legend.title=element_blank())
ggsave("figs/pub/1_revision1/propFishByChamber_M_means+CIs.png", p.N_M, width=9, height=4, dpi=300)









# preferred temperature ---------------------------------------------------

p.Temp <- ggplot(fit.temp_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp), size=0.4, alpha=0.25) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.2, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2), limits=c(0,13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F) + 
  geom_rug(data=tibble(ElapsedDays=0, 
                       pred=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=dark_1, xmax=dark_2)) +
  ylab("Preferred temperature (ºC)")
# ggsave("figs/pub/1_revision1/prefTemp.png", p.Temp, width=8, height=4, dpi=300)

p.M <- ggplot(fit.temp_global, aes(ElapsedDays, M, colour=Species, fill=Species)) +
  geom_hline(yintercept=temp_ctrl[1], colour=sp_col[1], linetype=5, size=0.1) +
  geom_hline(yintercept=temp_ctrl[2], colour=sp_col[2], linetype=5, size=0.1) +
  geom_ribbon(aes(ymin=M_lo, ymax=M_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  geom_rug(data=tibble(ElapsedDays=0, 
                       M=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=dark_1, xmax=dark_2)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2), limits=c(0,13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F) + 
  ylab("MESOR (ºC)")
# ggsave("figs/pub/smooth_M.png", p.M, width=8, height=4)

p.A <- ggplot(fit.temp_global, aes(ElapsedDays, amplitude, colour=Species, fill=Species)) +
  geom_ribbon(aes(ymin=amplitude_lo, ymax=amplitude_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=amp_1, ymax=amp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=amp_1, ymax=amp_2, xmin=dark_1, xmax=dark_2)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2), limits=c(0,13)) + 
  scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.5)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F)
# ggsave("figs/pub/smooth_A.png", p.A, width=8, height=4)

p.phi <- ggplot(fit.temp_global, aes(ElapsedDays, acrophase, colour=Species, fill=Species)) +
  geom_rect(aes(xmin=-0.25, xmax=-0.125, ymax=12, ymin=24), 
            colour="grey30", fill="grey30", size=0.25) +
  geom_rect(aes(xmin=-0.25, xmax=-0.125, ymax=0, ymin=12), 
            colour="grey30", fill="white", size=0.25) +
  geom_ribbon(aes(ymin=acrophase_lo, ymax=acrophase_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=acro_1, ymax=acro_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=acro_1, ymax=acro_2, xmin=dark_1, xmax=dark_2)) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2), limits=c(-0.25,13)) + 
  scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F)
# ggsave("figs/pub/smooth_phi.png", p.phi, width=8, height=4)

ggpubr::ggarrange(p.Temp, p.M, p.A, p.phi, ncol=2, nrow=2, common.legend=T, 
                  labels="auto", legend="bottom")
ggsave("figs/pub/1_revision1/prefTemp_all.png", width=9, height=7, dpi=300)


p.Temp2 <- ggplot(fit.temp_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  # geom_hline(yintercept=temp_ctrl[1], colour=sp_col[1], linetype=5, size=0.1) +
  # geom_hline(yintercept=temp_ctrl[2], colour=sp_col[2], linetype=5, size=0.1) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp, shape=Tank), alpha=0.5, size=0.8) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2), limits=c(0,13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.5), drop=F) + 
  scale_shape_manual(values=1:3) +
  guides(shape=guide_legend(override.aes=list(alpha=1, size=1), order=3),
         colour=guide_legend(order=1),
         fill=guide_legend(order=1),
         size=guide_legend(order=2)) +
  geom_rug(data=tibble(ElapsedDays=0, 
                       pred=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=dark_1, xmax=dark_2)) +
  ylab("Preferred temperature (ºC)")
ggsave("figs/pub/1_revision1/prefTemp_ALT_mn+CIs.png", p.Temp2, width=5.5, height=4, dpi=300)

p.legend <- p.M + 
  theme(legend.position="bottom") +
  guides(colour=guide_legend(order=1, ncol=1),
         fill=guide_legend(order=1, ncol=1),
         size=guide_legend(order=2, ncol=1))
ggpubr::ggarrange(p.M, p.A, p.phi, ncol=1, nrow=3, 
                  legend.grob=ggpubr::get_legend(p.legend), legend="bottom", labels="auto")
ggsave("figs/pub/1_revision1/prefTemp_ALT_cosinorPars.png", width=4, height=10, dpi=300)









# values ------------------------------------------------------------------

# Acrophase
fit.temp_global %>% filter(ElapsedDays %in% c(0,7)) %>% 
  select(Species, ElapsedDays, starts_with("acrophase"))

# Average daily max/min preferred temp
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species, Days) %>%
  filter(pred==max(pred)) %>%
  group_by(Species) %>%
  summarise(across(starts_with("pred"), mean))
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species, Days) %>%
  filter(pred==min(pred)) %>%
  group_by(Species) %>%
  summarise(across(starts_with("pred"), mean))
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species, Days) %>%
  arrange(pred) %>%
  mutate(diff=last(pred)-first(pred)) %>%
  slice_head(n=1) %>%
  group_by(Species) %>%
  summarise(diff=mean(diff))

# MESOR
fit.temp_global %>% filter(ElapsedDays==0) %>% 
  group_by(Species) %>%
  summarise(across(starts_with("pred"), mean))
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species) %>%
  summarise(across(starts_with("pred"), mean))

tau <- map(dir("out", "GS_predGlobal_posterior", full.names=T),
           ~readRDS(.x)$prefTemp)
day0 <- which(fit.temp_global$ElapsedDays==0)
day4 <- which(fit.temp_global$ElapsedDays==4)
acc_diff <- map(1:2, ~tau[[.x]][,day0[1]] - tau[[.x]][,day4[1]])
map(acc_diff, mean)
map(acc_diff, HDInterval::hdi)



# Fish proportions --------------------------------------------------------

fit.N %>%
  mutate(Days=floor(ElapsedDays)) %>%
  filter(ElapsedDays >= 7) %>%
  group_by(Group, Species, Chamber, Days) %>% 
  arrange(desc(pr_mn)) %>%
  slice_head(n=1) %>%
  group_by(Group, Species, Chamber) %>%
  summarise(across(starts_with("pr"), mean), ZT=mean((ElapsedDays %% 1)*24))

fit.N %>%
  mutate(Days=floor(ElapsedDays)) %>%
  filter(ElapsedDays < 1) %>%
  group_by(Group, Species, Chamber, Days) %>% 
  arrange(desc(pr_mn)) %>%
  slice_head(n=1) %>%
  group_by(Group, Species, Chamber) %>%
  summarise(across(starts_with("pr"), mean), ZT=mean((ElapsedDays %% 1)*24))
  

fit.N %>%
  filter(Group=="Control CTE") %>%
  group_by(Species, Chamber) %>%
  summarise(across(starts_with("pr"), mean), across(starts_with("M"), mean))



fit.temp_global %>% filter(ElapsedDays %in% c(0,7)) %>% 
  select(Species, ElapsedDays, starts_with("M"))
fit.temp_global %>% filter(ElapsedDays %in% c(0,7)) %>% 
  select(Species, ElapsedDays, starts_with("amplitude"))













# surplus -----------------------------------------------------------------

ggplot(fit.N, aes(ElapsedDays, A_mn, ymin=A_lo, ymax=A_hi,
                  group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  facet_grid(Species~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA))

ggplot(fit.N, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                  group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line(size=1) +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() + 
  facet_grid(Species~Group, scales="free_x", space="free_x")

ggplot(fit.N, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                  group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Group~Species)

ggplot(fit.N, aes(ElapsedDays, M_mn, fill=Chamber)) +
  geom_area(colour="black") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish (MESOR)", limits=c(0, 1)) +
  theme_classic() +
  facet_grid(Group~Species)

ggplot(fit.N, aes(ElapsedDays, phi_mn, ymin=phi_lo, ymax=phi_hi,
                  group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Chamber~Species)
