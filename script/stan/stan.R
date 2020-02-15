
library(StateSpaceInference)
library(parallel)
library(rstan)

sessionInfo()
set.seed(1)

# length of the time series
TT <- 40
# parameters
alpha <- 2; beta <- 0; gamma <- sqrt(2 * 0.1)/2; mu <- -0.2; phi <- 0.95; sh <- 0.6; s_v <- 1
# simulating the hidden states
h <- rep(0, TT)
h[1] <- rnorm(1, mu, sh/(sqrt(1-phi^2)))
for (t in 2:TT) {
  h[t] <- mu + phi * (h[t - 1]) + sh * rnorm(1)
}

# emission of the observations
yobs <- exp(h/2) * rnorm(TT, s_v, 0.1)

dat <- list(TT = TT, y = yobs, mu = mu, sh = sh)
fit <- stan(
  file = "../../../script/stan/stochvol.stan",
  model_name = "example",
  data = dat,
  iter = 500000,
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


#init = rep(list(list(theta = phi)), parallel::detectCores())
