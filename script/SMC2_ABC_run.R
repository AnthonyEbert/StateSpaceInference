
library(parallel)
library(StateSpaceModel)
sessionInfo()

cl <- makeCluster(parallel::detectCores() - 1)

lfunc <- function(...){parLapply(cl, ...)}
#lfunc <- function(...){mclapply(...)}
#lfunc <- function(...){lapply(...)}

library(StateSpaceModel)

# IBIS Algorithm
set.seed(3)

# Prior :
# sig ~ exp(1)
# tau ~ exp(1)
TT <- 30
true_theta <- c(0.25, 0.5)
#y <- data_simulator(true_theta, TT)


lower <- 0
upper <- 3.5
sd_t <- 1
init <- min(rgamma(1, 100, 100), upper - 1)
a_logit <- 0.9
dist_coef <- 0.5
true_states <- generate_state(init, TT, lower, upper, sd_t, a = a_logit)

#true_states <- c(0.5, 1, 2, 1, 3, 0.5, 1, 1, 2, 0.5, 0.5, 0.5, 0.5, 2, 2, 2, 4, 0.5, 1, 1, 1)

lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = true_states)
kern <- function(x){return(decay_func(x, alpha = true_theta[1], delta = true_theta[2]))}

y_history <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)
y <- hist(y_history, breaks = seq(0, TT*10, by = 10), plot = FALSE)$counts

hist(y_history, breaks = TT * 10)

plot(lambda_fun, add = TRUE, col = "red")

simulator <- StateSpaceModel::hawkes_simulator

inp <- list(
  lower = lower,
  upper = upper,
  sd_t = sd_t,
  a_logit = a_logit,
  history = y_history,
  simulator = simulator
)

library(StateSpaceInference)

rstate <- function(n, theta, inp){
  sapply(rep(NA, n), generate_statey, n = 1, lower = inp$lower, upper = inp$upper, sd = inp$sd_t, a = inp$a_logit)
}

loss = loss_hawkes

control <- list(
  Ntheta = 50,
  Nx = 200,
  TT = TT,
  pacc = 0.01
)

prior_sample <- rgamma(control$Ntheta, 10, 40)

x_list <- SMC2_ABC(prior_sample, rstate, loss, inp, control = control, cl = cl, dt = 10)

plotrix::weighted.hist(sapply(x_list, function(x){x$theta}), w = sapply(x_list, function(x){x$omega})/sum(sapply(x_list, function(x){x$omega})), breaks = 100)

