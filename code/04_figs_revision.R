# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk



# set up ------------------------------------------------------------------

library(tidyverse)
library(lisa)
library(glue)
theme_set(theme_classic())
group_col <- c(`Control CTE`="#1b9e77", Acclimation="#7570b3", Experiment="#d95f02")
group_raw <- c("Control", "Acclimation", "Experiment")
chmb_col <- c("1"="#0571b0", "2"="#92c5de", "3"="grey", "4"="#f4a582", "5"="#ca0020")
cmr.ls <- readRDS("figs/cmr_cmaps.RDS")
sp_col <- c("Nile tilapia"="#1a3431", "Zebrafish"="#6283c8")
species <- c("Nile tilapia"="Tilapia", "Zebrafish"="ZF")
temp_rng <- list(Tilapia=c(25.5, 34.5), ZF=c(24, 31.5))
temp_ctrl <- c(Tilapia=30.1, ZF=26.4)

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
  mutate(Chamber=factor(Chamber))
fit.temp_global <- map_dfr(species, ~read_csv(glue("out/GS_predGlobal_{.x}.csv"))) %>%
  mutate(Spline=factor("Global", levels=c("Global", "Tank")))
fit.temp_tank <- map_dfr(species, ~read_csv(glue("out/GS_predTank_{.x}.csv"))) %>%
  mutate(Tank=factor(Tank))



ggplot(fit.N, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                   group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() + 
  facet_grid(Species~Group)

ggplot(fit.N, aes(ElapsedDays, pr_mn, ymin=pr_lo, ymax=pr_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() + 
  facet_grid(Chamber~Species)

ggplot(fit.N, aes(ElapsedDays, pr_mn, fill=Chamber)) +
  geom_area(colour="grey30") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish", limits=c(0, 1)) +
  theme_classic() + 
  facet_grid(Group~Species)

ggplot(fit.N, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                   group=Chamber, colour=Chamber, fill=Chamber)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_colour_manual(values=chmb_col) +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Group~Species)

ggplot(fit.N, aes(ElapsedDays, M_mn, ymin=M_lo, ymax=M_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Chamber~Species)

ggplot(fit.N, aes(ElapsedDays, M_mn, fill=Chamber)) +
  geom_area(colour="black") +
  scale_fill_manual(values=chmb_col) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  scale_y_continuous("Predicted proportion of fish (MESOR)", limits=c(0, 1)) +
  theme_classic() +
  facet_grid(Group~Species)

ggplot(fit.N, aes(ElapsedDays, A_mn, ymin=A_lo, ymax=A_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Chamber~Species)

ggplot(fit.N, aes(ElapsedDays, phi_mn, ymin=phi_lo, ymax=phi_hi,
                   group=Group, colour=Group, fill=Group)) +
  geom_ribbon(alpha=0.5, colour=NA) + geom_line() +
  scale_x_continuous("Elapsed time (days)", breaks=0:13) +
  # scale_y_continuous("Predicted proportion of fish", limits=c(0, NA)) +
  theme_classic() +
  facet_grid(Chamber~Species)






# preferred temperature ---------------------------------------------------

lower_ZT.df <- tibble(light_1=0:12,
                      light_2=c(0:12)+0.5,
                      dark_1=c(0:12)+0.5,
                      dark_2=1:13,
                      temp_1=min(unlist(temp_rng)),
                      temp_2=min(unlist(temp_rng))+0.2,
                      amp_1=0, amp_2=0.025,
                      acro_1=0, acro_2=0.4,
                      ElapsedDays=0, pred=30, M=30, amplitude=0, acrophase=12)

p.Temp <- ggplot(fit.temp_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp), size=0.25, alpha=0.25) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.2, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13, limits=c(0,13)) + 
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
ggsave("figs/pub/smooth_prefTemp.png", p.Temp, width=8, height=4, dpi=300)

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
  scale_x_continuous("Elapsed time (days)", breaks=0:13, limits=c(0,13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F) + 
  ylab("MESOR (ºC)")
ggsave("figs/pub/smooth_M.png", p.M, width=8, height=4)

p.A <- ggplot(fit.temp_global, aes(ElapsedDays, amplitude, colour=Species, fill=Species)) +
  geom_ribbon(aes(ymin=amplitude_lo, ymax=amplitude_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank)), size=0.3) +
  geom_rect(data=lower_ZT.df, fill="white", colour="grey30", size=0.1,
            aes(ymin=amp_1, ymax=amp_2, xmin=light_1, xmax=light_2)) +
  geom_rect(data=lower_ZT.df, fill="grey30", colour="grey30", size=0.1,
            aes(ymin=amp_1, ymax=amp_2, xmin=dark_1, xmax=dark_2)) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13, limits=c(0,13)) + 
  scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.5)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F)
ggsave("figs/pub/smooth_A.png", p.A, width=8, height=4)

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
  scale_x_continuous("Elapsed time (days)", breaks=0:13, limits=c(-0.25,13)) + 
  scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.3), drop=F)
ggsave("figs/pub/smooth_phi.png", p.phi, width=8, height=4)

ggpubr::ggarrange(p.Temp, p.M, p.A, p.phi, ncol=2, nrow=2, common.legend=T, 
                  labels="auto")
ggsave("figs/pub/smooth_all_new.png", width=10, height=6, dpi=300)


p.Temp2 <- ggplot(fit.temp_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  # geom_hline(yintercept=temp_ctrl[1], colour=sp_col[1], linetype=5, size=0.1) +
  # geom_hline(yintercept=temp_ctrl[2], colour=sp_col[2], linetype=5, size=0.1) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp, shape=Tank), alpha=0.25) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.2, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=fit.temp_tank, aes(group=paste(Species, Tank), linetype=Tank), size=0.3) +
  scale_x_continuous("Elapsed time (days)", breaks=0:13, limits=c(0,13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  scale_size_manual(values=c(1, 0.5), drop=F) + 
  scale_shape_manual(values=1:3) +
  guides(shape=guide_legend(override.aes=list(alpha=1)),
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
ggsave("figs/pub/ALT_smooth_prefTemp.png", p.Temp2, width=8, height=6, dpi=300)

ggpubr::ggarrange(p.M, p.A, p.phi, ncol=1, nrow=3, 
                  common.legend=T, legend="bottom", labels="auto")
ggsave("figs/pub/ALT_smooth_cosinor_params.png", width=5, height=10, dpi=300)


fit.temp_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("M"))
fit.temp_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("A"))
fit.temp_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("phi"))
