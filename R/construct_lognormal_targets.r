## This is the old construction in terms of means and variances
# means <- c(1.64, 1.64)
# vars <- c(.5^2, 2.164285^2)
# mus <- log(means^2/sqrt(means^2 + vars))
# sigmas <- sqrt(log(1 + vars/means^2))

## Now we'll build stuff in terms of median and skewness 
ln_skew <- function(v) (exp(v) + 2)*sqrt(exp(v) - 1)
find_v <- function(s){
  cand_fun <- function(logv){
    ans <- ln_skew(exp(logv))
    return((s - ans)^2)
  }
  Opt <- optimise(cand_fun, interval = c(-10, 10))
  return(exp(Opt$minimum))
}
medians <- c(1, 1, 1)
skews <- c(ln_skew(1)/10, ln_skew(1)/5, ln_skew(1))


mus <- log(medians)
sigmas <- sqrt(sapply(skews, find_v))
cbind(mus, sigmas)

target.info <- data.frame(
  target = c("symmetric", "moderate", "asymmetric"),
  m = mus,
  sigma = sigmas,
  v = sigmas^2
)

save(target.info,
     file = "../saved_data/lognormal_targets.RData")