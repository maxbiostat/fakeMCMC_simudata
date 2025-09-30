generate_and_export <- function(Alpha, Beta,
                                target,
                                eff,
                                M = 1E4,
                                Nrep = 500,
                                Ncores = 10){
  
  Simus <- parallel::mclapply(1:Nrep,
                              function(i){
                                rgamma(
                                  n = M,
                                  shape = Alpha,
                                  rate = Beta
                                )
                              }, mc.cores = Ncores)
  
  simu.mat <- matrix(NA, ncol = M, nrow = Nrep)
  for(i in 1:Nrep){
    simu.mat[i, ] <- as.numeric(Simus[[i]])
  }
  
  write.csv(simu.mat,
            file = paste0("../saved_data/Gamma",
                          target,
                          "_IID",
                          "_eff=", eff,
                          "_M=", M,
                          ".csv" ),
            row.names = FALSE)
}
######

library(fakeMCMC)
load("../saved_data/gamma_targets.RData")

Target <- "moderate"
Alpha <- target.info[which(target.info$target == Target), ]$gamma_a_KL
Beta <- target.info[which(target.info$target == Target), ]$gamma_b_KL
Eff <- 1

generate_and_export(Alpha =  Alpha,
                    Beta = Beta,
                    eff = Eff,
                    target = Target)