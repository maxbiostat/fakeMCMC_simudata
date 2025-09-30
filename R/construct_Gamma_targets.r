approxGamma_mom <- function(mu, sigma){
  ## approximate a log-normal with parameters mu and sigma by a Gamma
  ### using method of moments
  s2 <- sigma^2
  m <- exp(mu + s2/2)
  v <- (exp(s2) - 1)*m^2
  alpha.star <- m^2/v
  beta.star <-  m/v
  return(
    c(alpha.star, beta.star)
  )
}
approxGamma_KL_bf <- function(mu, sigma){
  the_KL <- function(par){
    
    alpha <- exp(par[1])
    beta <- exp(par[2])
    
    integrand <- function(x){
      logf <- dlnorm(x, meanlog = mu, sdlog = sigma, log = TRUE)
      exp(logf) * 
        (logf - dgamma(x, shape = alpha, rate = beta, log = TRUE) )
    }
    KL <- integrate(integrand, 0, Inf)
    return(KL$value)
  }
  
  Opt <- optim(par = log(approxGamma_mom(mu, sigma)),
               fn =  the_KL)
  
  return(exp(Opt$par))
}
approxGamma_KL_anal <- function(mu, sigma){
  the_exact_KL <- function(par){
    a <- exp(par[1])
    b <- exp(par[2])
    
    res <- -a * log(b) + lgamma(a) - a*mu + b * exp(mu + sigma^2/2)
    
    return(res)
  }
  Opt <- optim(par = log(approxGamma_mom(mu, sigma)),
               fn =  the_exact_KL)
  return(exp(Opt$par))
} 


#######
approximate_and_plot <- function(Mu, Sigma, title = "",
                                 eps = 1E-3, plot = TRUE, leg = TRUE) {
  
  mom <- approxGamma_mom(mu = Mu, sigma = Sigma)
  # klbf <- approxGamma_KL_bf(mu = Mu, sigma = Sigma)
  klexact <- approxGamma_KL_anal(mu = Mu, sigma = Sigma)
  
  plms <- qlnorm(p = c(eps, 1 - eps), meanlog = Mu, sdlog = Sigma)
  
  if(plot){
    curve(dlnorm(x, meanlog = Mu, sdlog = Sigma),
          xlim = plms, lwd = 4, main = paste(title),
          xlab = expression(X), ylab = "Density")
    curve(dgamma(x, shape = mom[1], rate = mom[2]),
          lwd = 3, lty = 3, col = "purple", add = TRUE)
    curve(dgamma(x, shape = klexact[1], rate = klexact[2]),
          lwd = 3, lty = 4, col = "grey50", add = TRUE)
    if(leg){
      legend(x = "topright", 
             legend = c ("Method of moments", "Min KL"),
             col = c("purple", "grey50"), lwd = 3, bty  = 'n')
    }
  }
  return(
    list(par_mom = mom, 
         par_kl = klexact)
  )
}

load("../saved_data/lognormal_targets.RData")

mus <- target.info$m
sigmas <- target.info$sigma

par(mfrow = c(1, 3))
AppAsym <- approximate_and_plot(Mu = mus[3], Sigma = sigmas[3],
                     title = "Asymmetric")
AppMod <- approximate_and_plot(Mu = mus[2], Sigma = sigmas[2], 
                     title = "Moderate", leg = FALSE)
AppSym <- approximate_and_plot(Mu = mus[1], Sigma = sigmas[1],
                     title = "Symmetric", leg = FALSE)


target.info$gamma_a_mom <- c(AppSym$par_mom[1],
                             AppMod$par_mom[1],
                             AppAsym$par_mom[1])
target.info$gamma_b_mom <- c(AppSym$par_mom[2],
                             AppMod$par_mom[2],
                             AppAsym$par_mom[2])
target.info$gamma_a_KL <- c(AppSym$par_kl[1],
                            AppMod$par_kl[1],
                            AppAsym$par_kl[1])
target.info$gamma_b_KL <-  c(AppSym$par_kl[2],
                             AppMod$par_kl[2],
                             AppAsym$par_kl[2])

save(target.info,
     file = "../saved_data/gamma_targets.RData")