library(sn)
library(rstan)
options(mc.cores = 2)
rstan_options(auto_write = TRUE)

# Simulating data
TT <- 40
x <- rnorm(TT)
x <- cumsum(x)

y <- matrix(0, nrow = TT, ncol = 10)
for (j in 1:10) {
  y[, j] <- x + rsn(TT, 0, .25, .5)
}

datastan <- list(TT = 40, y = y)
fit <- stan(
  file = "script/stan/skewednormal.stan",
  model_name = "sn",
  data = datastan,
  iter = 100000,
  chains = 2,
  cores = 2,
  pars = c("sigma", "gamma")
)

library(shinystan)
launch_shinystan(fit)
