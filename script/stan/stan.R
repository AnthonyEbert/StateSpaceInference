
library(StateSpaceInference)
library(parallel)
library(rstan)

sessionInfo()

#cl <- makeCluster(parallel::detectCores() - 1)
#cl = "mclapply"
cl <- NULL

# length of the time series
TT <- 40
# parameters
mu <- -0.2; phi <- 0.7; sh <- 0.6
# simulating the hidden states
h <- rep(0, TT)
h[1] <- rnorm(1)
for (t in 2:TT) {
  h[t] <- mu + phi * h[t - 1] + sh * rnorm(1)
}

# emission of the observations
yobs <- exp(h/2) * rnorm(TT, 0, 1)

dat <- list(TT = TT, y = yobs, mu = mu, sh = sh)
fit <- stan(file = "script/stan/stochvol.stan", model_name = "example",
            data = dat, iter = 100, chains = 3)
print(fit)
