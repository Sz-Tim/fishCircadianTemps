

phi_to_ZT <- function(phi) {
  phi <- phi + 2*pi*(phi < -pi) - 2*pi*(phi > pi)
  phi <- (-phi+pi)*12/pi-12
  phi <- phi + 24*(phi < 0)
}


calc_prefTemp <- function(M, A, phi, ZT) {
  M + exp(A) * cos(3.141593*ZT/12 + phi)
}


summarise_N_posteriors <- function(out, newdata) {
  
  # Posterior global smooths for M
  post_M <- map(c(1,2,4,5), ~(posterior_smooths(out, newdata=newdata,
                                                smooth="s(ElapsedTime_sc,m=2)", nlpar=paste0("M", .x)) +
                                c(as_draws_matrix(out, variable=paste0("b_M", .x, "_Intercept")))))
  post_M.ar <- array(0, dim=c(nrow(post_M[[1]]), ncol(post_M[[1]]), 5))
  post_M.ar[,,1] <- exp(post_M[[1]])
  post_M.ar[,,2] <- exp(post_M[[2]])
  post_M.ar[,,3] <- exp(0)
  post_M.ar[,,4] <- exp(post_M[[3]])
  post_M.ar[,,5] <- exp(post_M[[4]])
  for(i in 1:dim(post_M.ar)[1]) {
    for(j in 1:dim(post_M.ar)[2]) {
      post_M.ar[i,j,] <- post_M.ar[i,j,]/sum(post_M.ar[i,j,])
    }
  }
  post_M_mn <- apply(post_M.ar, 2:3, mean)
  post_M_hdi <- apply(post_M.ar, 2:3, HDInterval::hdi)
  
  # Posterior global smooths for A
  post_A <- map(c(1,2,4,5), ~(posterior_smooths(out, newdata=newdata, 
                                                smooth="s(ElapsedTime_sc,m=2)", nlpar=paste0("A", .x)) +
                                c(as_draws_matrix(out, variable=paste0("b_A", .x, "_Intercept")))))
  post_A.ar <- array(0, dim=c(nrow(post_M[[1]]), ncol(post_M[[1]]), 5))
  post_A.ar[,,1] <- exp(post_A[[1]])
  post_A.ar[,,2] <- exp(post_A[[2]])
  post_A.ar[,,3] <- exp(0)
  post_A.ar[,,4] <- exp(post_A[[3]])
  post_A.ar[,,5] <- exp(post_A[[4]])
  post_A_mn <- apply(post_A.ar, 2:3, mean)
  post_A_hdi <- apply(post_A.ar, 2:3, HDInterval::hdi)
  
  # Posterior global smooths for phi
  post_phi <- map(c(1,2,4,5), ~(posterior_smooths(out, newdata=newdata, 
                                                  smooth="s(ElapsedTime_sc,m=2)", nlpar=paste0("pphi", .x)) +
                                  c(as_draws_matrix(out, variable=paste0("b_pphi", .x, "_Intercept")))))
  post_phi.ar <- array(0, dim=c(nrow(post_M[[1]]), ncol(post_M[[1]]), 5))
  post_phi.ar[,,1] <- post_phi[[1]]
  post_phi.ar[,,2] <- post_phi[[2]]
  post_phi.ar[,,3] <- 0
  post_phi.ar[,,4] <- post_phi[[3]]
  post_phi.ar[,,5] <- post_phi[[4]]
  post_phi_mn <- apply(post_phi.ar, 2:3, mean) %>% phi_to_ZT()
  post_phi_hdi <- apply(post_phi.ar, 2:3, HDInterval::hdi) %>% phi_to_ZT()
  
  post_pr <- map(1:4,
                 ~calc_prefTemp(post_M[[.x]], post_A[[.x]], post_phi[[.x]],
                                matrix(newdata$ZT, nrow=nrow(post_M[[.x]]), ncol=nrow(newdata), byrow=T)))
  post_pr.ar <- array(0, dim=c(nrow(post_M[[1]]), ncol(post_M[[1]]), 5))
  post_pr.ar[,,1] <- exp(post_pr[[1]])
  post_pr.ar[,,2] <- exp(post_pr[[2]])
  post_pr.ar[,,3] <- exp(0)
  post_pr.ar[,,4] <- exp(post_pr[[3]])
  post_pr.ar[,,5] <- exp(post_pr[[4]])
  for(i in 1:dim(post_pr.ar)[1]) {
    for(j in 1:dim(post_pr.ar)[2]) {
      post_pr.ar[i,j,] <- post_pr.ar[i,j,]/sum(post_pr.ar[i,j,])
    }
  }
  post_pr_mn <- apply(post_pr.ar, 2:3, mean)
  post_pr_hdi <- apply(post_pr.ar, 2:3, HDInterval::hdi)
  
  
  fit.df <- tibble(Chamber=factor(rep(1:5, each=nrow(newdata))),
                   ElapsedDays=rep(newdata$ElapsedDays, 5)) %>%
    mutate(pr_mn=c(post_pr_mn),
           pr_lo=c(post_pr_hdi[1,,]),
           pr_hi=c(post_pr_hdi[2,,]),
           M_mn=c(post_M_mn),
           M_lo=c(post_M_hdi[1,,]),
           M_hi=c(post_M_hdi[2,,]),
           A_mn=c(post_A_mn),
           A_lo=c(post_A_hdi[1,,]),
           A_hi=c(post_A_hdi[2,,]),
           phi_mn=c(post_phi_mn),
           phi_lo=c(post_phi_hdi[1,,]),
           phi_hi=c(post_phi_hdi[2,,]),
           Group=newdata$Group[1],
           Species=newdata$Species[1])
  return(list(summaries=fit.df, 
              pr.post=post_pr.ar, 
              M.post=post_M.ar, 
              A.post=post_A.ar, 
              phi.post=post_phi.ar))
}
