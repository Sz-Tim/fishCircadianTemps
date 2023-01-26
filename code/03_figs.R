# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
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
         Group=factor(Group, levels=c(3, 1, 2), labels=group_raw),
         Tank=factor(Tank)) %>%
  left_join(read_csv(glue("data/temp_{.x}.csv")) %>%
              pivot_longer(starts_with("chamber_"), names_to="Chamber", values_to="Temp") %>%
              mutate(Chamber=factor(str_sub(Chamber, -1, -1)),
                     Group=factor(Group, 
                                  levels=group_raw))) %>%
  mutate(GrpDay=factor(paste(as.numeric(Group), str_pad(Days, 2, "l", "0"), sep="_"))))

c_temp.sum <- map(data.c, ~.x %>% 
                    filter(Group != "Control") %>%
                    group_by(Chamber) %>% 
                    summarise(Temp=mean(Temp, na.rm=T)))

out <- map(species, ~readRDS(glue("models/cosinor/out_temperature_vm_expA_{.x}.rds")))
out.c <- map(species, ~readRDS(glue("models/cosinor/out_count_0CI_noLB_Cor_{.x}.rds")))
out.GS_global <- map_dfr(species, ~read_csv(glue("out/GS_predGlobal_{.x}.csv")))
out.GS_tank <- map_dfr(species, ~read_csv(glue("out/GS_predTank_{.x}.csv")))

pred.df <- map(data.df, ~filter(.x, Group != "Control") %>% droplevels) %>%
  map(~expand_grid(ZT=seq(0, 24, length.out=100),
                   Group=unique(.x$Group)))
preds <- map2(out, pred.df, 
              ~posterior_epred(.x, newdata=.y, re.form=NA))
pred.df <- map2(pred.df, preds, 
                ~.x %>%
                  mutate(pred.mn=apply(.y, 2, mean),
                         CI80_l=apply(.y, 2, function(x) hdci(x, 0.8)[,1]),
                         CI80_h=apply(.y, 2, function(x) hdci(x, 0.8)[,2]),
                         CI90_l=apply(.y, 2, function(x) hdci(x, 0.9)[,1]),
                         CI90_h=apply(.y, 2, function(x) hdci(x, 0.9)[,2]),
                         CI95_l=apply(.y, 2, function(x) hdci(x, 0.95)[,1]),
                         CI95_h=apply(.y, 2, function(x) hdci(x, 0.95)[,2])))



# pref temp by ZT ---------------------------------------------------------

fig1.ls <- map(2:1, 
               ~pred.df[[.x]] %>%
                 filter(Group != "Control") %>%
                 mutate(Group=factor(Group, labels=names(group_col)[2:3])) %>%
                 ggplot(aes(ZT)) + 
                 geom_jitter(data=data.df[[.x]] %>%
                               filter(Group != "Control") %>%
                               mutate(Group=factor(Group, labels=names(group_col)[2:3])), 
                             aes(y=prefTemp, shape=Group, colour=Group),
                             alpha=0.8, size=0.95, width=0.4, height=0) +
                 geom_line(aes(y=pred.mn, colour=Group)) + 
                 geom_ribbon(aes(ymin=CI95_l, ymax=CI95_h, fill=Group), colour=NA, alpha=0.4) +
                 geom_rect(aes(xmin=0, xmax=12, ymax=temp_rng[[.x]][2],
                               ymin=temp_rng[[.x]][2]-diff(temp_rng[[.x]])*0.02), 
                           colour="grey30", fill="white", size=0.25) +
                 geom_rect(aes(xmin=12, xmax=24, ymax=temp_rng[[.x]][2],
                               ymin=temp_rng[[.x]][2]-diff(temp_rng[[.x]])*0.02), 
                           colour="grey30", fill="grey30", size=0.25) +
                 geom_text(data=c_temp.sum[[.x]], x=-1.5, colour=chmb_col,
                           size=3.5, aes(label=Chamber, y=Temp)) +
                 scale_fill_manual(values=group_col[2:3]) +
                 scale_colour_manual(values=group_col[2:3]) + 
                 scale_shape_manual(values=1:2) +
                 scale_x_continuous("ZT (h)", breaks=seq(0, 24, by=6), limits=c(-1.25, 24)) +
                 scale_y_continuous("Preferred temperature (ºC)", limits=temp_rng[[.x]],
                                    breaks=seq(ceiling(temp_rng[[.x]][1]),
                                               floor(temp_rng[[.x]][2]),
                                               by=2)) +
                 ggtitle(names(pred.df)[.x]) +
                 theme(legend.position=c(0.8, 0.15),
                       legend.title=element_blank(), 
                       legend.background=element_blank()))

fig.1 <- ggpubr::ggarrange(plotlist=fig1.ls, ncol=1, common.legend=T, 
                           legend="bottom", labels=c("a.", "b."))
ggsave("figs/pub/pref_temp_by_ZT_expA.png", width=4, height=7, units="in")




# effects -----------------------------------------------------------------

paramSum.ls <- map_dfr(
  out,
  ~.x %>%
    as_draws_df(variable="^b_(M|A|phi)", regex=T) %>%
    rename_with(~str_remove(.x, "^b_")) %>%
    rename_with(~str_remove(.x, "Group")) %>%
    rename_with(~str_replace(.x, "Intercept", "Control")) %>%
    select(-starts_with(".")) %>%
    pivot_longer(everything(), names_to="param_OG", values_to="draws") %>%
    filter(param_OG != "phi_Control"),
  .id="Species") %>%
  mutate(Param=str_split_fixed(param_OG, "_", 2)[,1],
         Group=str_split_fixed(param_OG, "_", 2)[,2]) %>%
  # mutate(draws=case_when(param_OG=="phi_Acclimation" & Species=="Nile tilapia" ~ (-draws+pi)*12/pi-12,
                         # param_OG=="phi_Experiment" & Species=="Nile tilapia" ~ (-draws+pi)*12/pi+12,
                         # Param=="phi" & Species=="Zebrafish" ~ (-draws+pi)*12/pi-12,
                         # Param!="phi" ~ draws)) %>%
  mutate(draws=case_when(Param=="phi" ~ draws + 2*pi*(draws< -pi) - 2*pi*(draws > pi),
                         Param!="phi" ~ draws)) %>%
  mutate(draws=case_when(Param=="phi" ~ (-draws + pi)*12/pi-12,
                         Param!="phi" ~ draws)) %>%
  mutate(draws=case_when(Param=="phi" ~ draws + 24*(draws < 0),
                         Param!="phi" ~ draws)) %>%
  mutate(draws=case_when(Param=="A" ~ exp(draws),
                         Param!="A" ~ draws)) %>%
  group_by(Species, Param, Group, param_OG) %>%
  summarise(mn=mean(draws),
            md=median(draws),
            CI80_l=hdci(draws, 0.8)[,1],
            CI80_h=hdci(draws, 0.8)[,2],
            CI90_l=hdci(draws, 0.9)[,1],
            CI90_h=hdci(draws, 0.9)[,2],
            CI95_l=hdci(draws, 0.95)[,1],
            CI95_h=hdci(draws, 0.95)[,2]) %>%
  ungroup %>%
  mutate(Group=factor(Group, levels=group_raw, labels=names(group_col)),
         Param=factor(Param, levels=c("M", "A", "phi"),
                      labels=c("MESOR", "Amplitude", "Acrophase"))) %>%
  filter(Group != "Control CTE") %>% droplevels %>%
  group_by(Param) %>%
  group_split()

param_i <- tibble(Param=c("MESOR", "Amplitude", "Acrophase"),
                  labs=c("MESOR (ºC)", "Amplitude", "Acrophase (h)"),
                  lim=list(range(unlist(temp_rng)), 
                           c(0, 1), 
                           c(0, 24)),
                  breaks=list(seq(min(unlist(temp_rng)), max(unlist(temp_rng)), by=2),
                              c(0, 0.5, 1),
                              c(0, 6, 12, 18, 24)))
temp_rng.df <- tibble(Species=names(species),
                      mn=map_dbl(temp_rng, mean),
                      tmin=map_dbl(temp_rng, ~.x[1]),
                      tmax=map_dbl(temp_rng, ~.x[2]))
                
pt_width <- 0.5
fig5.ls <- map(1:3, 
    ~ggplot(paramSum.ls[[.x]], 
            aes(y=Species, x=mn, colour=Group)) +
      geom_point(size=2, position=position_dodge(width=pt_width)) +
      geom_errorbarh(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, height=0,
                     position=position_dodge(width=pt_width)) +
      geom_errorbarh(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, height=0,
                     position=position_dodge(width=pt_width)) +
      geom_errorbarh(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, height=0.2,
                     position=position_dodge(width=pt_width)) +
      scale_colour_manual(values=group_col[2:3]) + 
      scale_x_continuous(param_i$labs[.x], 
                         limits=param_i$lim[[.x]],
                         breaks=param_i$breaks[[.x]]) +
      {if(.x>1) theme(axis.text.y=element_blank())} +
      {if(.x==1) theme(axis.title.y=element_blank())} +
      {if(.x>1) ylab("")} +
      {if(param_i$Param[.x]=="MESOR") 
        geom_errorbarh(data=temp_rng.df, aes(xmin=tmin, xmax=tmax),
                       # position=position_nudge(y=pt_width/6),
                       colour="grey", size=0.25, height=0.25)} +
      {if(param_i$Param[.x]=="Acrophase") 
        geom_rect(aes(xmin=0, xmax=12, ymax=2.5, ymin=2.58), 
                  colour="grey30", fill="white", size=0.25)} +
      {if(param_i$Param[.x]=="Acrophase") 
        geom_rect(aes(xmin=12, xmax=24, ymax=2.5, ymin=2.58), 
                  colour="grey30", fill="grey30", size=0.25)} +
      theme(legend.title=element_blank()))
p5.abc <- ggpubr::ggarrange(plotlist=fig5.ls, nrow=1, common.legend=T, legend="right", 
                  labels=paste0(letters[1:length(fig5.ls)], "."), 
                  widths=c(1.1, 1, 1))
ggsave("figs/pub/effects_expA.png", width=8, height=2, dpi=300)








# heatmaps ----------------------------------------------------------------
pred.ls <- map2(data.c, out.c,
               ~.x %>% filter(Tank==1) %>%
                 mutate(pred.mn=colMeans(pmax(
                   exp(posterior_predict(.y, 
                                       newdata=.x %>% filter(Tank==1),
                                       re_formula=NA))-1, 
                   0))) %>%
                 mutate(Group=factor(Group, labels=names(group_col)),
                        Temp=Temp + rnorm(n(), 0, 0.2)) %>%
                 group_by(Group) %>%
                 group_split())
nBins <- 20
temp.heat.ZF <- map(2:3, 
                    ~pred.ls[[2]][[.x]] %>%
                      select(pred.mn, ZT, Group, Temp) %>% 
                      mutate(FishCount=round(pred.mn)) %>% 
                      uncount(FishCount) %>% 
                      ggplot(aes(ZT, Temp)) +
                      geom_density2d_filled(colour=NA, bins=nBins) + 
                      geom_rect(aes(xmin=0, xmax=12, ymin=temp_rng[[2]][2],
                                    ymax=temp_rng[[2]][2]+diff(temp_rng[[2]])*0.04), 
                                colour="grey30", fill="white", size=0.25) +
                      geom_rect(aes(xmin=12, xmax=24, ymin=temp_rng[[2]][2],
                                    ymax=temp_rng[[2]][2]+diff(temp_rng[[2]])*0.04), 
                                colour="grey30", fill="grey30", size=0.25) +
                      scale_fill_manual(values=lisa_palette("KatsushikaHokusai", nBins, "continuous"),
                                        guide="none") +
                      scale_x_continuous("ZT (h)", breaks=seq(0, 24, by=6)) +
                      scale_y_continuous(limits=c(temp_rng[[2]][1], 
                                                  temp_rng[[2]][2]+diff(temp_rng[[2]])*0.04),
                                         breaks=seq(ceiling(temp_rng[[2]][1]),
                                                    floor(temp_rng[[2]][2]),
                                                    by=2)) +
                      {if(.x==2) ylab("Temperature (ºC)")} +
                      {if(.x==2) ggtitle("Zebrafish")} +
                      {if(.x>2) ylab("")} +
                      {if(.x>2) ggtitle("")} +
                      {if(.x>2) theme(axis.text.y=element_blank())} +
                      facet_wrap(~Group))
temp.heat.T <- map(2:3, 
                   ~pred.ls[[1]][[.x]] %>%
                     select(pred.mn, ZT, Group, Temp) %>% 
                     mutate(FishCount=round(pred.mn)) %>% 
                     uncount(FishCount) %>% 
                     ggplot(aes(ZT, Temp)) +
                     geom_density2d_filled(colour=NA, bins=nBins) + 
                     geom_rect(aes(xmin=0, xmax=12, ymin=temp_rng[[1]][2],
                                   ymax=temp_rng[[1]][2]+diff(temp_rng[[1]])*0.04), 
                               colour="grey30", fill="white", size=0.25) +
                     geom_rect(aes(xmin=12, xmax=24, ymin=temp_rng[[1]][2],
                                   ymax=temp_rng[[1]][2]+diff(temp_rng[[1]])*0.04), 
                               colour="grey30", fill="grey30", size=0.25) +
                     scale_fill_manual(values=lisa_palette("KatsushikaHokusai", nBins, "continuous"),
                                       guide="none") +
                     scale_x_continuous("ZT (h)", breaks=seq(0, 24, by=6)) +
                     scale_y_continuous(limits=c(temp_rng[[1]][1], 
                                                 temp_rng[[1]][2]+diff(temp_rng[[1]])*0.04),
                                        breaks=seq(ceiling(temp_rng[[1]][1]),
                                                   floor(temp_rng[[1]][2]),
                                                   by=2)) +
                     {if(.x==2) ylab("Temperature (ºC)")} +
                     {if(.x==2) ggtitle("Nile tilapia")} +
                     {if(.x>2) ylab("")} +
                     {if(.x>2) ggtitle("")} +
                     {if(.x>2) theme(axis.text.y=element_blank())} +
                     facet_wrap(~Group))
hm.plot <- ggpubr::ggarrange(plotlist=c(temp.heat.ZF, temp.heat.T), nrow=2, ncol=2, 
                  labels=c("a.", "", "b.", ""), 
                  widths=c(1.05, 1))
ggsave("figs/pub/heatmap_cosinor.png", hm.plot, width=5, height=5.25, units="in")








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
  mutate(Group=factor(Group, labels=names(group_col)),
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
ggsave("figs/pub/tube_plot_cosinor.png", fig.tube, width=5, height=4, units="in")






# acrophase ---------------------------------------------------------------

# WARNING: If the model was re-run, it might use exp(A), which may affect this.
# I think I fixed it, but if anything looks weird that's probably it
acro.df <- expand_grid(ZT=seq(0, 24, length.out=100),
                       Chamber=1:5,
                       Group=levels(data.c[[1]]$Group))
acro.pred <- map_dfr(names(species), 
                 ~posterior_epred(out.c[[.x]], newdata=acro.df, re.form=NA) %>%
                   t() %>%
                   as_tibble %>%
                   bind_cols(acro.df, .) %>%
                   pivot_longer(starts_with("V"),
                                names_to="iter",
                                values_to="pred") %>%
                   mutate(iter=as.numeric(str_sub(iter, 2, -1))) %>%
                   group_by(Chamber, Group, iter) %>%
                   slice_max(order_by=pred), 
                 .id="Species")

acro.sum <- acro.pred %>% ungroup %>%
  mutate(ZT=case_when(Species!=1 | Chamber!=2 | Group!="Experiment" | ZT<12 ~ ZT,
                      Species==1 & Chamber==2 & Group=="Experiment" & ZT>12 ~ ZT-24)) %>%
  group_by(Chamber, Group, Species) %>%
  summarise(mn=mean(ZT),
            md=median(ZT),
            CI80_l=quantile(ZT, 0.1),
            CI80_h=quantile(ZT, 0.9),
            CI90_l=quantile(ZT, 0.05),
            CI90_h=quantile(ZT, 0.95),
            CI95_l=quantile(ZT, 0.025),
            CI95_h=quantile(ZT, 0.975)) %>%
  ungroup %>%
  mutate(Group=factor(Group, levels=group_raw, labels=names(group_col)),
         Species=factor(Species, levels=1:2, labels=names(species)),
         Chamber=factor(Chamber, levels=5:1)) 
acro.sum <- bind_rows(acro.sum %>%
                        mutate(across(starts_with("CI"), ~.x*(.x>0))), 
                      acro.sum %>% 
                        filter(Species=="Nile tilapia" & Chamber==2 & Group=="Experiment") %>%
                        mutate(mn=NA, md=NA) %>%
                        mutate(across(starts_with("CI"), ~.x*(.x<0)+24)))
A_nonZero <- map_dfr(out.c, 
                     ~fixef(.x) %>%
                       as_tibble(rownames="variable") %>%
                       filter(grepl("^A_", variable)) %>%
                       mutate(nonZero=exp(Q2.5) > 0.01),
                     .id="Species") %>%
  mutate(Chamber=factor(str_sub(variable, 10, 10), levels=1:5),
         Group=case_when(grepl("Acclim", variable) ~ "Acclimation",
                         grepl("Experi", variable) ~ "Experiment",
                         !grepl("Acc|Exp", variable) ~ "Control")) %>%
  mutate(Group=factor(Group, levels=group_raw, labels=names(group_col))) %>%
  select(Species, Chamber, Group, nonZero)

acro.sum <- acro.sum %>% 
  full_join(., A_nonZero, by=c("Species", "Chamber", "Group")) %>%
  mutate(CI80_l=if_else(nonZero, CI80_l, md),
         CI80_h=if_else(nonZero, CI80_h, md),
         CI90_l=if_else(nonZero, CI90_l, md),
         CI90_h=if_else(nonZero, CI90_h, md)) %>%
  mutate(Species=factor(Species, levels=rev(names(species))))

pt_width <- 0.5
p5.d <- acro.sum %>%
  mutate(Group=factor(Group, levels=rev(levels(Group)))) %>%
  ggplot(aes(y=Group, x=md, colour=Chamber, shape=nonZero)) +
  geom_point(size=2, position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, height=0.2,
                 position=position_dodge(width=pt_width)) + 
  scale_x_continuous("Acrophase (h)", limits=c(0,24), breaks=c(0,6,12,18,24)) +
  geom_rect(aes(xmin=0, xmax=12, ymax=3.5, ymin=3.58), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymax=3.5, ymin=3.58), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_colour_manual(values=chmb_col) +
  scale_shape_manual("Amplitude > 0", values=c(1, 19)) +
  facet_wrap(~Species) +
  theme(axis.title.y=element_blank(),
        legend.title=element_text(size=10),
        legend.key.height=unit(0.25, "cm"))
p5.d_ALT <- acro.sum %>%
  filter(nonZero) %>%
  mutate(Group=factor(Group, levels=rev(levels(Group)))) %>%
  ggplot(aes(y=Group, x=md, colour=Chamber)) +
  geom_point(size=2, position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, height=0,
                 position=position_dodge(width=pt_width)) + 
  geom_point(data=acro.sum %>%
               filter(nonZero) %>%
               mutate(Group=factor(Group, levels=rev(levels(Group))),
                      CI95_l=replace(CI95_l, CI95_l==0, NA)),
             aes(x=CI95_l), size=1.75, shape="|",
             position=position_dodge(width=pt_width)) + 
  geom_point(data=acro.sum %>%
               filter(nonZero) %>%
               mutate(Group=factor(Group, levels=rev(levels(Group))),
                      CI95_h=replace(CI95_h, CI95_h==24, NA)),
             aes(x=CI95_h), size=1.75, shape="|",
             position=position_dodge(width=pt_width)) + 
  scale_x_continuous("Acrophase (h)", limits=c(0,24), breaks=c(0,6,12,18,24)) +
  geom_rect(aes(xmin=0, xmax=12, ymax=2.5, ymin=2.58), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymax=2.5, ymin=2.58), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_colour_manual(values=chmb_col) +
  facet_wrap(~Species) +
  theme(axis.title.y=element_blank(),
        legend.title=element_text(size=10),
        legend.key.height=unit(0.25, "cm"))
p5.d <- ggpubr::ggarrange(p5.d_ALT, labels="d.")

p5 <- ggpubr::ggarrange(p5.abc, p5.d, common.legend=F, nrow=2, ncol=1, 
                        heights=c(1, 1))
ggsave("figs/pub/effects_all_ALT.png", p5, width=8, height=4, dpi=300)

pt_width <- 0.5
p5.d_alt <- ggplot(acro.sum, aes(x=md, y=Chamber, colour=Group, shape=nonZero)) +
  geom_point(size=2, position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI80_l, xmax=CI80_h), size=1.25, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI90_l, xmax=CI90_h), size=0.75, height=0,
                 position=position_dodge(width=pt_width)) +
  geom_errorbarh(aes(xmin=CI95_l, xmax=CI95_h), size=0.25, height=0.2,
                 position=position_dodge(width=pt_width)) + 
  scale_x_continuous("Acrophase (h)", limits=c(0,24), breaks=c(0,6,12,18,24)) +
  geom_rect(aes(xmin=0, xmax=12, ymax=5.5, ymin=5.58), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=12, xmax=24, ymax=5.5, ymin=5.58), 
            colour="grey30", fill="grey30", size=0.25) +
  scale_colour_manual("", values=group_col) +
  scale_shape_manual("Amplitude > 0", values=c(1, 19)) +
  facet_wrap(~Species)
p5.d_alt <- ggpubr::ggarrange(p5.d_alt, labels="d.")
p5.abc_alt <- ggpubr::ggarrange(plotlist=fig5.ls, nrow=1, common.legend=T, legend="none", 
                                labels=paste0(letters[1:length(fig5.ls)], "."), 
                                widths=c(1.1, 1, 1))
p5_alt <- ggpubr::ggarrange(p5.abc_alt, p5.d_alt, heights=c(0.6, 1),
                            common.legend=T, legend="bottom",
                            nrow=2, ncol=1)
ggsave("figs/pub/effects_all_alt.png", p5_alt, width=8, height=4.5, dpi=300)



acro.sum %>% 
  arrange(Species, Group, desc(Chamber)) %>% 
  select(Species, Group, Chamber, mn, CI95_l, CI95_h, nonZero) %>% 
  mutate(across(where(is.numeric), ~signif(.x, 3))) %>% 
  mutate(entry=if_else(nonZero, paste0(mn, " [", CI95_l, "-", CI95_h, "]"), "---")) %>%
  write_csv("out/chamber_acrophase.csv")




# max temperature ---------------------------------------------------------

pred.df <- map(data.df,
               ~expand_grid(ZT=seq(0, 24, length.out=24*60),
                            Group=unique(.x$Group)))
preds <- map2(out, pred.df,
              ~posterior_epred(.x, newdata=.y, re.form=NA))
pred.df <- map2(pred.df, preds,
                ~.x %>%
                  mutate(pred.mn=apply(.y, 2, mean),
                         CI80_l=apply(.y, 2, function(x) hdci(x, 0.8)[,1]),
                         CI80_h=apply(.y, 2, function(x) hdci(x, 0.8)[,2]),
                         CI90_l=apply(.y, 2, function(x) hdci(x, 0.9)[,1]),
                         CI90_h=apply(.y, 2, function(x) hdci(x, 0.9)[,2]),
                         CI95_l=apply(.y, 2, function(x) hdci(x, 0.95)[,1]),
                         CI95_h=apply(.y, 2, function(x) hdci(x, 0.95)[,2])))
map(pred.df, ~.x %>% group_by(Group) %>% slice_max(pred.mn))
map(pred.df, ~.x %>% group_by(Group) %>% slice_min(pred.mn))


# radial plots ------------------------------------------------------------

pred.CHmax <- map2(out.c, tube.base, 
                   ~posterior_epred(.x, newdata=.y, re.form=NA) %>%
                     t %>% as_tibble() %>%
                     cbind(.y) %>%
                     group_by(ZT, Group) %>%
                     arrange(Chamber) %>%
                     mutate(across(starts_with("V"), which.max)) %>%
                     select(-Chamber) %>%
                     pivot_longer(starts_with("V"), names_to="iter", values_to="Chamber_max"))

chPr.ls <- bind_rows(
  map_dfr(pred.CHmax, 
          ~.x %>% 
            group_by(ZT, Group, Chamber_max) %>%
            summarise(N=n()) %>%
            group_by(ZT, Group) %>%
            mutate(pr=N/sum(N)),
          .id="Species") %>%
    ungroup %>%
    mutate(source="Fitted"),
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
    mutate(source="Observed")) %>%
  mutate(Species=factor(Species, levels=rev(names(species))),
         Chamber_max=factor(Chamber_max),
         Group=factor(Group, labels=names(group_col)),
         source=factor(source, levels=c("Observed", "Fitted"))) %>%
  group_by(Species) %>%
  group_split()

fig7 <- map(chPr.ls, 
    ~.x %>%
      ggplot(aes(ZT, y=pr, fill=Chamber_max)) + 
      geom_bar(stat="identity", colour="grey30", size=0.1) +
      geom_rect(aes(xmin=-0.5, xmax=11.5, ymax=1.1, ymin=1.05), 
                colour="grey30", fill="white", size=0.25) +
      geom_rect(aes(xmin=11.5, xmax=23.5, ymax=1.1, ymin=1.05), 
                colour="grey30", fill="grey30", size=0.25) +
      scale_x_continuous(breaks=seq(0,23,by=2), limits=c(-0.5, 23.5)) +
      scale_y_continuous(breaks=c(0, 0.5, 1)) +
      scale_fill_brewer("Chamber", type="div", palette="RdBu", direction=-1) +
      facet_grid(source~Group) +   
      coord_polar(theta="x", start=(-0.5*pi*2)/24) +
      theme(panel.grid.major=element_line(size=0.1, colour="grey"),
            legend.key.size=unit(0.2, "cm")) +
      labs(x="ZT (h)", y="Probability of highest fish density",
           title=first(.x$Species)))
p <- ggpubr::ggarrange(plotlist=fig7, ncol=1, labels=c("a.", "b."), 
                  common.legend=T, legend="bottom")
ggsave("figs/pub/radial_chamber_distr_cosinor.png", p, width=5.5, height=8)




 # comparisons -------------------------------------------------------------

# Preferred temperature
hyp_M <- c(M_AccExp="M_GroupAcclimation = M_GroupExperiment")
map(out, ~hypothesis(.x, hyp_M))

hyp_A <- c(A_Acc="exp(A_GroupAcclimation) > 0.01",
           A_Exp="exp(A_GroupExperiment) > 0.01",
           A_AccExp="A_GroupAcclimation = A_GroupExperiment")
map(out, ~hypothesis(.x, hyp_A))

hyp_phi <- c(phi_AccExp="phi_GroupAcclimation = phi_GroupExperiment",
             phi_AccExp2="phi_GroupAcclimation - phi_GroupExperiment = -6.283185")
map(out, ~hypothesis(.x, hyp_phi))


# Preferred chamber
hyp_Ctrl_M.df <- expand_grid(a=1:5, b=1:5) %>%
  filter(a!=b) %>%
  mutate(hyp=paste0("M_Chamber", a, " = M_Chamber", b),
         name=paste0("M_", a, b))
hyp_M <- setNames(hyp_Ctrl_M.df$hyp, hyp_Ctrl_M.df$name)
map(out.c, ~hypothesis(.x, hyp_M))

# Tilapia prefer 1 vs. 2, 3, or 4; and 5 vs. 4 with 95% confidence





# chamber cosinor sum -----------------------------------------------------

chamber_param.sum <- map_dfr(
  out.c,
  ~.x %>%
    as_draws_df(variable="^b_(M|A|phi)", regex=T) %>%
    rename_with(~str_remove(.x, "^b_")) %>%
    rename_with(~str_remove(.x, "Group")) %>%
    rename_with(~str_replace(.x, "Intercept", "Control")) %>%
    select(-starts_with(".")) %>%
    pivot_longer(everything(), names_to="param_OG", values_to="draws") %>%
    filter(param_OG != "phi_Control"),
  .id="Species") %>%
  mutate(Param=str_split_fixed(param_OG, "_", 2)[,1],
         Group=case_when(grepl("Acc", param_OG) ~ "Acclimation",
                         grepl("Exp", param_OG) ~ "Experiment",
                         !grepl("Acc", param_OG) & !grepl("Exp", param_OG) ~ "Control"),
         Chamber=str_sub(str_split_fixed(param_OG, "Chamber", 2)[,2], 1, 1)) %>%
  filter(Param != "phi") %>%
  group_by(Species, Param, Group, Chamber) %>%
  mutate(iter=row_number()) %>%
  group_by(Species, Param, Chamber, iter)  %>%
  mutate(draws=case_when(Group=="Control" ~ draws,
                         Group!="Control" ~ first(draws)+draws)) %>%
  mutate(draws=case_when(Param=="M" ~ expm1(draws),
                         Param!="M" ~ draws)) %>%
  mutate(draws=case_when(Param=="A" ~ exp(draws),
                         Param!="A" ~ draws))
  group_by(Species, Param, Group, Chamber) %>%
  summarise(mn=mean(draws),
            md=median(draws),
            CI95_l=hdci(draws, 0.95)[,1],
            CI95_h=hdci(draws, 0.95)[,2]) %>%
  ungroup %>%
  mutate(Group=factor(Group, levels=group_raw, labels=names(group_col)),
         Chamber=factor(Chamber),
         Param=factor(Param, levels=c("M", "A"),
                      labels=c("MESOR", "Amplitude")))




# smooth acclimation ------------------------------------------------------

pred_s.dat <- expand_grid(ZT=0:23, Days=1:13, Tank=1) %>%
  mutate(ElapsedTime=ZT+24*(Days-1),
         ElapsedTime_sc=c(scale(ElapsedTime))) 
pred_s.ls <- map(out.s,
                 ~list(
                   prefTemp=posterior_epred(.x, newdata=pred_s.dat, re_formula=NA),
                   M=posterior_epred(.x, newdata=pred_s.dat, re_formula=NA, nlpar="M"),
                   A=posterior_epred(.x, newdata=pred_s.dat, re_formula=NA, nlpar="A"),
                   phi=posterior_epred(.x, newdata=pred_s.dat, re_formula=NA, nlpar="phi")
                 ))
for(i in seq_along(pred_s.ls)) {
  pred_s.ls[[i]]$A <- exp(pred_s.ls[[i]]$A)
  pred_s.ls[[i]]$phi <- phi_to_ZT(pred_s.ls[[i]]$phi)
}
pred_s.df <- map_dfr(
  pred_s.ls, 
  ~pred_s.dat %>%
    mutate(prefTemp=colMeans(.x$prefTemp),
           prefTemp_lo=apply(.x$prefTemp, 2, function(x) quantile(x, probs=0.025)),
           prefTemp_hi=apply(.x$prefTemp, 2, function(x) quantile(x, probs=0.975)),
           M=colMeans(.x$M),
           M_lo=apply(.x$M, 2, function(x) quantile(x, probs=0.025)),
           M_hi=apply(.x$M, 2, function(x) quantile(x, probs=0.975)),
           A=colMeans(.x$A),
           A_lo=apply(.x$A, 2, function(x) quantile(x, probs=0.025)),
           A_hi=apply(.x$A, 2, function(x) quantile(x, probs=0.975)),
           phi=colMeans(.x$phi),
           phi_lo=apply(.x$phi, 2, function(x) quantile(x, probs=0.025)),
           phi_hi=apply(.x$phi, 2, function(x) quantile(x, probs=0.975))), 
  .id="Species")

p.Temp <- ggplot(pred_s.df, aes(ElapsedTime, prefTemp, colour=Species, fill=Species)) +
  geom_vline(xintercept=72, linetype=2, colour="grey30") +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             size=0.25, alpha=0.25) +
  geom_ribbon(aes(ymin=prefTemp_lo, ymax=prefTemp_hi), alpha=0.3, colour=NA) +
  geom_line() +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  geom_rug(data=tibble(ElapsedTime=0, 
                       prefTemp=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  ylab("Preferred temperature (ºC)")
ggsave("figs/pub/smooth_prefTemp.png", p.Temp, width=8, height=4, dpi=300)

p.M <- ggplot(pred_s.df, aes(ElapsedTime, M, colour=Species, fill=Species)) +
  geom_vline(xintercept=72, linetype=2, colour="grey30") +
  geom_ribbon(aes(ymin=M_lo, ymax=M_hi), alpha=0.3, colour=NA) +
  geom_line() +
  geom_rug(data=tibble(ElapsedTime=0, 
                       M=unlist(temp_rng), 
                       Species=rep(names(sp_col), each=2)), sides="l", size=1) +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) +
  ylab("MESOR (ºC)")
ggsave("figs/pub/smooth_M.png", p.M, width=8, height=4)

p.A <- ggplot(pred_s.df, aes(ElapsedTime, A, colour=Species, fill=Species)) +
  geom_vline(xintercept=72, linetype=2, colour="grey30") +
  geom_ribbon(aes(ymin=A_lo, ymax=A_hi), alpha=0.3, colour=NA) +
  geom_line() +
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_y_continuous("Amplitude", breaks=c(0, 0.5, 1), limits=c(0, 1.2)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) 
ggsave("figs/pub/smooth_A.png", p.A, width=8, height=4)

p.phi <- ggplot(pred_s.df, aes(ElapsedTime, phi, colour=Species, fill=Species)) +
  geom_vline(xintercept=72, linetype=2, colour="grey30") +
  geom_rect(aes(xmin=-6, xmax=-3, ymax=12, ymin=24), 
            colour="grey30", fill="white", size=0.25) +
  geom_rect(aes(xmin=-6, xmax=-3, ymax=0, ymin=12), 
            colour="grey30", fill="grey30", size=0.25) +
  geom_ribbon(aes(ymin=phi_lo, ymax=phi_hi), alpha=0.3, colour=NA) +
  geom_line() + 
  scale_x_continuous("Elapsed hours", breaks=24*(0:13)) + 
  scale_y_continuous("Acrophase (ZT h)", breaks=c(0, 6, 12, 18, 24), limits=c(0,24)) + 
  scale_colour_manual(values=sp_col) +
  scale_fill_manual(values=sp_col) 
ggsave("figs/pub/smooth_phi.png", p.phi, width=8, height=4)

ggpubr::ggarrange(p.Temp, p.M, p.A, p.phi, ncol=2, nrow=2, common.legend=T, 
                  labels="auto")
ggsave("figs/pub/smooth_all.png", width=10, height=6, dpi=300)


pred_s.df %>% filter(ElapsedTime %in% (24*c(0, 3))) %>% 
  select(Species, ElapsedTime, starts_with("M"))
pred_s.df %>% filter(ElapsedTime %in% (24*c(0, 3))) %>% 
  select(Species, ElapsedTime, starts_with("A"))
pred_s.df %>% filter(ElapsedTime %in% (24*c(0, 3))) %>% 
  select(Species, ElapsedTime, starts_with("phi"))







# HGAM --------------------------------------------------------------------

out.GS_global <- out.GS_global %>% mutate(Spline=factor("Global", levels=c("Global", "Tank")))
out.GS_tank <- out.GS_tank %>% mutate(Tank=factor(Tank))
lower_ZT.df <- tibble(light_1=0:12,
                      light_2=c(0:12)+0.5,
                      dark_1=c(0:12)+0.5,
                      dark_2=1:13,
                      temp_1=min(unlist(temp_rng)),
                      temp_2=min(unlist(temp_rng))+0.2,
                      amp_1=0, amp_2=0.025,
                      acro_1=0, acro_2=0.4,
                      ElapsedDays=0, pred=30, M=30, amplitude=0, acrophase=12)

p.Temp <- ggplot(out.GS_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp), size=0.25, alpha=0.25) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.2, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=out.GS_tank, aes(group=paste(Species, Tank)), size=0.3) +
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

p.M <- ggplot(out.GS_global, aes(ElapsedDays, M, colour=Species, fill=Species)) +
  geom_ribbon(aes(ymin=M_lo, ymax=M_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=out.GS_tank, aes(group=paste(Species, Tank)), size=0.3) +
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

p.A <- ggplot(out.GS_global, aes(ElapsedDays, amplitude, colour=Species, fill=Species)) +
  geom_ribbon(aes(ymin=amplitude_lo, ymax=amplitude_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=out.GS_tank, aes(group=paste(Species, Tank)), size=0.3) +
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

p.phi <- ggplot(out.GS_global, aes(ElapsedDays, acrophase, colour=Species, fill=Species)) +
  geom_rect(aes(xmin=-0.25, xmax=-0.125, ymax=12, ymin=24), 
            colour="grey30", fill="grey30", size=0.25) +
  geom_rect(aes(xmin=-0.25, xmax=-0.125, ymax=0, ymin=12), 
            colour="grey30", fill="white", size=0.25) +
  geom_ribbon(aes(ymin=acrophase_lo, ymax=acrophase_hi), alpha=0.25, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=out.GS_tank, aes(group=paste(Species, Tank)), size=0.3) +
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


p.Temp2 <- ggplot(out.GS_global, aes(ElapsedDays, pred, colour=Species, fill=Species)) +
  geom_point(data=map_dfr(data.df, ~filter(.x, Group!="Control"), .id="Species"), 
             aes(y=prefTemp, shape=Tank), alpha=0.25) +
  geom_ribbon(aes(ymin=pred_lo, ymax=pred_hi), alpha=0.2, colour=NA) +
  geom_line(aes(size=Spline)) +
  geom_line(data=out.GS_tank, aes(group=paste(Species, Tank), linetype=Tank), size=0.3) +
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


out.GS_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("M"))
out.GS_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("A"))
out.GS_global %>% filter(ElapsedDays %in% c(0,3)) %>% 
  select(Species, ElapsedDays, starts_with("phi"))

