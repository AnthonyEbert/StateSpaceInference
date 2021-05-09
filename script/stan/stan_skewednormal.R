
options(mc.cores = parallel::detectCores())

library(StateSpaceInference)
library(parallel)
library(rstan)

sessionInfo()

seed <- 2
set.seed(seed)

#cl <- makeCluster(parallel::detectCores())
#cl = "mclapply"
#cl <- NULL

TT <- 20
x <- rnorm(TT)
x <- cumsum(x)

z <- generate_stan_skew(TT, x, c(0.25, 2))

y <- matrix(0, nrow = TT, ncol = 100)
for (j in 1:TT) {
  y[j, ] <- z[[j]]
}

datastan <- list(TT = TT, y = y)
fit <- stan(
  file = "../../../script/stan/skewednormal.stan",
  data = datastan,
  iter = 1000,
  chains = 2*parallel::detectCores(),
  cores = parallel::detectCores(),
  pars = c("sigma", "gamma")
)

print(fit)
theta_stan = extract(fit, pars = c("sigma", "gamma"))

save.image()
save(theta_stan, file = "theta_stan.RData")

stan_df <- data.frame(
  sigma = theta_stan$sigma,
  gamma = theta_stan$gamma,
  weight = 1/length(theta_stan$sigma),
  seed = seed,
  type = "stan"
)

save(stan_df, file = "stan_df.RData")

saveRDS(stan_df, file = paste0("theta_stan_", seed,".RData"))

#library(shinystan)
#launch_shinystan(fit)
