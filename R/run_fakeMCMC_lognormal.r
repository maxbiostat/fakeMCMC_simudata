generate_and_export <- function(Mu, Sigma,
                                Phi, eff,
                                target,
                                M = 1E4,
                                Nrep = 500,
                                Ncores = 10){
  
  Simus <- parallel::mclapply(1:Nrep,
                              function(i){
                                generate_fake_MCMC(
                                  N = M,
                                  lmean = Mu,
                                  lsd = Sigma,
                                  phi = Phi,
                                  LN = TRUE
                                )
                              }, mc.cores = Ncores)
  
  simu.mat <- matrix(NA, ncol = M, nrow = Nrep)
  for(i in 1:Nrep){
    simu.mat[i, ] <- as.numeric(Simus[[i]])
  }
  
  write.csv(simu.mat,
            file = paste0("../saved_data/LN",
                          target,
                          "_LNAR",
                          "_eff=", eff,
                          "_M=", M,
                          ".csv" ),
            row.names = FALSE)
}
######

library(fakeMCMC)
load("../saved_data/lognormal_targets.RData")

Target <- "moderate"
Mu <- target.info[which(target.info$target == Target), ]$m
Sigma <- target.info[which(target.info$target == Target), ]$s
### IID \approx 0.99999999
Eff <- .99999999
if(Eff > .99){
  Phi <- 0
}else{
  Phi <- fakeMCMC::find_phi_LNAR1(eff = Eff, lmean = Mu, lsd = Sigma)
}

generate_and_export(Mu = Mu,
                    Sigma = Sigma,
                    Phi = Phi,
                    eff = Eff,
                    target = Target)
