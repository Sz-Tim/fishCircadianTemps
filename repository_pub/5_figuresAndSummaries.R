# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk
# 5. Create figures and summarise results



dir.create("figs")
# set up ------------------------------------------------------------------

library(tidyverse)
library(glue)
theme_set(theme_classic())
# colours from lisa::lisa_palette("GustavKlimt", 5)
group_col <- c("Control CTE"="#609F5C", "Experiment"="#A27CBA") 
chmb_col <- c("1"="#0571b0", "2"="#92c5de", "3"="grey", "4"="#f4a582", "5"="#ca0020")
sp_col <- c("Zebrafish"="#6283c8", "Nile tilapia"="#1a3431")
species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")



# load data & output ------------------------------------------------------

dataT.df <- map(species, ~read_csv(glue("out/prefTemp_data_{.x}.csv")) %>%
                  mutate(ElapsedDays=ElapsedTime/24,
                         Tank=factor(Tank)))
dataN.df <- map(species, 
                ~readRDS(glue("out/chmbrComp_data_{.x}.rds")) %>%
                  select(-Y) %>% 
                  pivot_longer(starts_with("Ch_"), names_to="Chamber", values_to="propFish") %>%
                  mutate(Group=factor(Group, levels=c("Control", "Experiment"), labels=names(group_col)),
                         Chamber=str_remove(Chamber, "Ch_")) %>%
                  filter(ElapsedTime < (24*13)))


c_temp.sum <- map_dfr(species, 
                      ~read_csv(glue("temperature_{.x}.csv")) %>%
                        pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
                        group_by(Group, Chamber) %>% 
                        summarise(Temp=mean(Temp, na.rm=T)), 
                      .id="Species") %>%
  mutate(M_lo=0, M_mn=0.85, M_hi=0, ElapsedDays=11.5,
         TempRnd=paste0(format(Temp, digits=3), "\u00B0C"),
         Chamber=str_replace(Chamber, "chamber_", "Chamber ")) %>%
  mutate(M_mn=M_mn+(Group=="Control")*0.1,
         Species=factor(Species, levels=names(species)),
         pr_mn=M_mn, pr_lo=M_lo, pr_hi=M_hi,
         Group=factor(Group, levels=c("Control", "Experiment"), labels=names(group_col)))

fit.temp_global <- map_dfr(species, ~read_csv(glue("out/prefTemp_predGlobal_{.x}.csv"))) %>%
  mutate(Spline=factor("Global", levels=c("Global", "Tank")),
         Species=factor(Species, levels=species, labels=names(species)))
fit.temp_tank <- map_dfr(species, ~read_csv(glue("out/prefTemp_predTank_{.x}.csv"))) %>%
  mutate(Tank=factor(Tank),
         Species=factor(Species, levels=species, labels=names(species)))
fit.N <- map_dfr(species, ~read_csv(glue("out/chmbrComp_predGlobal_{.x}.csv"))) %>%
  mutate(Chamber=factor(Chamber),
         Group=factor(Group, levels=c("Control", "Experimental"), labels=names(group_col)),
         Species=factor(Species, levels=species, labels=names(species)))
fit.edge <- map_dfr(species, ~read_csv(glue("out/chmbrComp_predEdge_{.x}.csv"))) %>%
  mutate(Group=factor(Group, levels=c("Control", "Experimental"), labels=names(group_col)),
         Species=factor(Species, levels=species, labels=names(species))) %>%
  rename(pr_mn=edge_mn, pr_lo=edge_lo, pr_hi=edge_hi)
fit.cw <- map_dfr(species, ~read_csv(glue("out/chmbrComp_predColdWarm_{.x}.csv"))) %>%
  mutate(Chambers=factor(Chambers),
         Group=factor(Group, levels=c("Control", "Experimental"), labels=names(group_col)),
         Species=factor(Species, levels=species, labels=names(species)))



# plot elements -----------------------------------------------------------

temp_rng <- list(ZF=c(24, 31.5), Tilapia=c(25.5, 34.5))
temp_ctrl <- c(ZF=26.4, Tilapia=30.1)

lower_ZT.df <- tibble(light_1=0:12,
                      light_2=c(0:12)+0.5,
                      dark_1=c(0:12)+0.5,
                      dark_2=1:13,
                      temp_1=min(unlist(temp_rng)),
                      temp_2=min(unlist(temp_rng))+0.2,
                      amp_1=0, amp_2=0.025,
                      acro_1=0, acro_2=0.4,
                      prop_1=-0.05, prop_2=-0.02,
                      ElapsedDays=0, prefTemp=30, M=30, amplitude=0, acrophase=12, pr_mn=0)
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
         ElapsedDays=0, ElapsedTime=0, pred=30, propFish=0,
         M=30, M_mn=0, amplitude=0, acrophase=12, pr_mn=0),
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
         ElapsedDays=0, ElapsedTime=0, pred=30, propFish=0,
         M=30, M_mn=0, amplitude=0, acrophase=12, pr_mn=0))



# Figure 2 ----------------------------------------------------------------

p.Temp2 <- ggplot(fit.temp_global, aes(ElapsedDays, prefTemp, colour=Species, fill=Species)) +
  geom_point(data=map_dfr(dataT.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp, shape=Tank), alpha=0.5, size=0.8) +
  geom_ribbon(aes(ymin=prefTemp_lo, ymax=prefTemp_hi), alpha=0.25, colour=NA) +
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
                       prefTemp=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=temp_1, ymax=temp_2, xmin=dark_1, xmax=dark_2)) +
  ylab("Preferred temperature (ºC)")
ggsave("figs/Figure_2.png", p.Temp2, width=5.5, height=4, dpi=300)



# Figure 3 ----------------------------------------------------------------

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

p.legend <- p.M + 
  theme(legend.position="bottom") +
  guides(colour=guide_legend(order=1, ncol=1),
         fill=guide_legend(order=1, ncol=1),
         size=guide_legend(order=2, ncol=1))
ggpubr::ggarrange(p.phi, p.M, p.A, ncol=1, nrow=3, 
                  legend.grob=ggpubr::get_legend(p.legend), legend="bottom", labels="auto")
ggsave("figs/Figure_3.png", width=4, height=10, dpi=300)



# Figure 4 ----------------------------------------------------------------

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
  facet_grid(Species~Group, scales="free_x", space="free_x") +
  theme_classic() +
  theme(legend.position="bottom")
ggsave("figs/Figure_4.png", p.N_prMn, width=7, height=4, dpi=300)



# Figure A1-2 -------------------------------------------------------------

obs.prop.ls <- imap(dataN.df,
                    ~ggplot(.x, aes(ElapsedTime/24, propFish, fill=Chamber)) + 
                      geom_hline(yintercept=c(0,1), colour="grey50", size=0.25) +
                      geom_hline(yintercept=c(0.2,0.4,0.6,0.8), colour="grey50", size=0.25) +
                      geom_rect(data=lower_ZT.Grp, fill="white", colour="grey30", size=0.1,
                                aes(ymin=prop_1, ymax=prop_2, xmin=light_1, xmax=light_2)) +
                      geom_rect(data=lower_ZT.Grp, fill="grey30", colour="grey30", size=0.1,
                                aes(ymin=prop_1, ymax=prop_2, xmin=dark_1, xmax=dark_2)) +
                      geom_area(colour="grey30", size=0.2, aes(group=paste(Chamber, Chunk))) +
                      scale_fill_manual(values=chmb_col) +
                      scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
                      scale_y_continuous("Observed fish distribution", breaks=seq(0,1,by=0.2)) +
                      ggtitle(.y) +
                      facet_grid(Tank~Group, scales="free_x", space="free_x") +
                      theme_bw() +
                      theme(panel.grid=element_blank(),
                            legend.position="bottom")) %>%
  setNames(c("A1", "A2"))
iwalk(obs.prop.ls, 
      ~ggsave(glue("figs/Figure_{.y}.png"), .x, width=8, height=5, dpi=300))



# Figure A3 ---------------------------------------------------------------

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
  geom_text(data=c_temp.sum, size=2.5, aes(label=TempRnd), show.legend=F) +
  scale_colour_manual(values=group_col) +
  scale_fill_manual(values=group_col) +
  scale_x_continuous("Elapsed time (days)", breaks=seq(0,13,by=2)) +
  scale_y_continuous("Proportion of fish\n(posterior mean + 95% HDI)", 
                     limits=c(-0.05, NA), breaks=seq(0,1,by=0.2)) +
  facet_grid(Species~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom")
ggsave("figs/Figure_A3.png", p.N_prMnCI, width=9, height=4, dpi=300)



# Figure A4 ---------------------------------------------------------------

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
  facet_grid(Species~Chamber) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom", 
        axis.title.y=element_blank())

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
  facet_grid(Species~Chambers) +
  theme_classic() +
  theme(panel.border=element_rect(fill=NA),
        legend.position="bottom")

ggpubr::ggarrange(p.ColdWarm_MMnCI, p.Edge_MMnCI, ncol=2, legend="bottom", 
                  common.legend=T, widths=c(1, 0.52), labels="auto")
ggsave("figs/Figure_A4.png", width=9, height=5.25, dpi=300)



# Summaries: Preferred temperature ----------------------------------------

# * Acrophase -------------------------------------------------------------
# Day 0 vs. Day 7
fit.temp_global %>% filter(ElapsedDays %in% c(0,7)) %>% 
  select(Species, ElapsedDays, starts_with("acrophase"))
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species) %>%
  summarise(across(starts_with("acrophase"), mean))

# * Daily max/min temperatures --------------------------------------------
# Mean maximum across Days 7-13
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species, Day) %>%
  filter(prefTemp==max(prefTemp)) %>%
  group_by(Species) %>%
  summarise(across(starts_with("prefTemp"), mean))
# Mean minimum across Days 7-13
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species, Day) %>%
  filter(prefTemp==min(prefTemp)) %>%
  group_by(Species) %>%
  summarise(across(starts_with("prefTemp"), mean))

# * MESOR -----------------------------------------------------------------
# Day 0
fit.temp_global %>% filter(ElapsedDays==0) %>% 
  group_by(Species) %>%
  summarise(across(starts_with("M"), mean))
# Mean across Days 7-13
fit.temp_global %>% filter(ElapsedDays>=7) %>% 
  group_by(Species) %>%
  summarise(across(starts_with("M"), mean))
# Difference from Day 0 to Day 7
M <- imap(species, ~readRDS(glue("out/prefTemp_posteriorGlobal_{.x}.rds"))$M)
day0 <- which(fit.temp_global$ElapsedDays==0)[1]
day7 <- which(fit.temp_global$ElapsedDays==7)[1]
acclimation_diff <- map(tau, ~.x[,day7] - .x[,day0])
map(acclimation_diff, mean)
map(acclimation_diff, HDInterval::hdi)



# Summaries: Chamber composition ------------------------------------------
# Mean peaks across Days 7-13
fit.N %>%
  mutate(Day=floor(ElapsedDays)) %>%
  filter(ElapsedDays >= 7) %>%
  group_by(Group, Species, Chamber, Day) %>% 
  arrange(desc(pr_mn)) %>%
  slice_head(n=1) %>%
  group_by(Group, Species, Chamber) %>%
  summarise(across(starts_with("pr"), mean), ZT=mean((ElapsedDays %% 1)*24))
# Mean peak for Day 0
fit.N %>%
  mutate(Day=floor(ElapsedDays)) %>%
  filter(ElapsedDays < 1) %>%
  group_by(Group, Species, Chamber, Day) %>% 
  arrange(desc(pr_mn)) %>%
  slice_head(n=1) %>%
  group_by(Group, Species, Chamber) %>%
  summarise(across(starts_with("pr"), mean), ZT=mean((ElapsedDays %% 1)*24))
# Edge preference during control
chPost <- imap(species, ~readRDS(glue("out/chmbrComp_posteriorGlobal_{.x}.rds")))
M.ctrl <- imap(chPost, ~.x$ctrl$M.post)
M.exp <- imap(chPost, ~.x$exp$M.post)
imap(M.ctrl, ~mean(c(.x[,,1]+.x[,,5])))
imap(M.ctrl, ~HDInterval::hdi(c(.x[,,1]+.x[,,5]), credMass=0.9))
imap(M.ctrl, ~HDInterval::hdi(c(.x[,,1]+.x[,,5]), credMass=0.85))
imap(M.exp, ~mean(c(.x[,,1]+.x[,,5])))
imap(M.exp, ~HDInterval::hdi(c(.x[,,1]+.x[,,5]), credMass=0.9))
imap(M.exp, ~HDInterval::hdi(c(.x[,,1]+.x[,,5]), credMass=0.85))
# Experiment vs. mean(control)
fit.N %>%
  filter(Group=="Experiment") %>%
  full_join(fit.N %>% 
              filter(Group=="Control CTE") %>%
              group_by(Chamber, Species) %>%
              summarise(pr_ctrl=mean(pr_mn))) %>%
  ggplot(aes(ElapsedDays, pr_mn-pr_ctrl, ymin=pr_lo-pr_ctrl, ymax=pr_hi-pr_ctrl)) + 
  geom_hline(yintercept=c(0)) +
  geom_line() + geom_ribbon(alpha=0.25, colour=NA) +
  scale_y_continuous("Pr_exp - Pr_ctrl", breaks=c(-0.2, 0, 0.2, 0.4)) +
  facet_grid(Species~Chamber) +
  theme_bw()
fit.N %>%
  filter(Group=="Experiment") %>%
  full_join(fit.N %>% 
              filter(Group=="Control CTE") %>%
              group_by(Chamber, Species) %>%
              summarise(M_ctrl=mean(M_mn))) %>%
  ggplot(aes(ElapsedDays, M_mn-M_ctrl, ymin=M_lo-M_ctrl, ymax=M_hi-M_ctrl)) + 
  geom_hline(yintercept=c(0)) +
  geom_line() + geom_ribbon(alpha=0.25, colour=NA) +
  scale_y_continuous("M_exp - M_ctrl", breaks=c(-0.2, 0, 0.2, 0.4)) +
  facet_grid(Species~Chamber) +
  theme_bw()
