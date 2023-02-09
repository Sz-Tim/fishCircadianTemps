# Tim Szewczyk
# Circadian rhythms in temperature preferences of fish
# tim.szewczyk@sams.ac.uk
# 4. Chamber composition: process output




setwd("repository_pub")
# switches ----------------------------------------------------------------

species <- c("Zebrafish"="ZF", "Nile tilapia"="Tilapia")[1]



# set up ------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
source("00_fn.R")

data.df <- read_csv(glue("out/chmbrComp_data_{species}.csv"))
time_sc <- readRDS(glue("out/chmbrComp_timeScale_{species}.csv"))

out.exp <- readRDS(glue("out/chmbrComp_mod_exp_{species}"))
out.ctrl <- readRDS(glue("out/chmbrComp_mod_ctrl_{species}"))





# initialise dataframes ---------------------------------------------------

fit.exp <- expand_grid(Day=1:13, ZT=0:23, Tank=3) %>%
  mutate(ElapsedTime=ZT+24*(Day-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24,
         Species=species,
         Group="Experimental")
fit.ctrl <- expand_grid(Day=1:6, ZT=0:23, Tank=3) %>%
  mutate(ElapsedTime=ZT+24*(Day-1),
         ElapsedTime_sc=(ElapsedTime-attr(time_sc, "scaled:center"))/attr(time_sc, "scaled:scale"),
         ElapsedDays=ElapsedTime/24,
         Species=species,
         Group="Control")



# summarise posteriors ----------------------------------------------------

fit.ls <- list(
  exp=summarise_N_posteriors(out.exp, fit.exp),
  ctrl=summarise_N_posteriors(out.ctrl, fit.ctrl)
)
fit.df <- rbind(fit.ls$exp$summaries,
                fit.ls$ctrl$summaries)
pr.post <- list(exp=fit.ls$exp$pr.post,
                ctrl=fit.ls$ctrl$pr.post)
M.post <- list(exp=fit.ls$exp$M.post,
               ctrl=fit.ls$ctrl$M.post)
rm(fit.ls)

write_csv(fit.df, glue("out/chmbrComp_predGlobal_{species}.csv"))



# aggregate chambers for comparisons --------------------------------------

fit.dat <- bind_rows(fit.exp, fit.ctrl)

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
write_csv(edge.df, glue("out/chmbrComp_predEdge_{species}.csv"))


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
write_csv(cw.df, glue("out/chmbrComp_predColdWarm_{species}.csv"))
