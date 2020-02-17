
options(mc.cores = parallel::detectCores())

library(StateSpaceInference)
library(parallel)
library(rstan)

sessionInfo()

seed <- 10
set.seed(seed)

cl <- makeCluster(parallel::detectCores())
#cl = "mclapply"
#cl <- NULL

# length of the time series
TT <- 20
# parameters
alpha <- 2; beta <- 0; gamma <- 0.1 * sqrt(1/2); mu <- 1; phi <- 0.80; sh <- 0.6; s_v <- 1
# simulating the hidden states
h <- rep(0, TT)
h[1] <- rnorm(1, mu/(1-phi), sd = sqrt(sh^2/(1-phi^2)))
for (t in 2:TT) {
  h[t] <- mu + phi * (h[t - 1]) + sh * rnorm(1)
}

# emission of the observations
yobs <- exp(h/2) * stabledist::rstable(TT, alpha, beta, gamma, s_v)


dat <- list(TT = TT, y = yobs, mu = mu, sh = sh)
fit <- stan(
  file = "../../../script/stan/stochvol.stan",
  model_name = "example",
  data = dat,
  iter = 100000,
  chains = parallel::detectCores(),
  cores = parallel::detectCores(),
  init = rep(list(list(theta = phi, x = h)), parallel::detectCores()),
  control = list(adapt_delta = 0.99),
  pars = "theta"
)
print(fit)
theta_stan = extract(fit, pars = "theta")$theta

save.image()
save(theta_stan, file = "theta_stan.RData")

stan_df <- data.frame(
  value = theta_stan,
  weight = 1/length(theta_stan),
  seed = seed,
  type = "stan"
)

saveRDS(stan_df, file = paste0("theta_stan_", seed,".RData"))



#init = rep(list(list(theta = phi)), parallel::detectCores())
