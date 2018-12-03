
library(StateSpaceModel)

# IBIS Algorithm
set.seed(3)

# Prior :
# sig ~ exp(1)
# tau ~ exp(1)
TT <- 40
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







