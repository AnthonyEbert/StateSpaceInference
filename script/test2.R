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

library(StateSpaceInference)

Ntheta <- 200
Nx <- 500
acceptance_threshold <- 0.02
x_list <- as.list(rep(NA, Ntheta))
x_list <- lapply(x_list, function(i){list()})

eps <- rep(NA, 40)

for(m in 1:Ntheta){
  x_list[[m]]$theta <- rgamma(1, 1*10, 4*10)
  x_list[[m]]$x <- sapply(rep(NA, Nx), generate_statey, n = 1, lower = lower, upper = upper, sd = sd_t, a = a_logit)
  x_list[[m]]$w <- rep(1, Nx)
  x_list[[m]]$p <- 1
  x_list[[m]]$omega <- 1
}


for(t in 2:40){
  x_list <- lfunc(x_list, theta_filter, time1 = t*10 - 10, time2 = t*10, loss_args = inp)

  distances <- sapply(x_list, function(x){x$distance})
  eps[t] <- quantile(as.numeric(distances), probs = acceptance_threshold)

  for(m in 1:Ntheta){
    x_list[[m]]$w <- (x_list[[m]]$distance <= eps[t])*1
    x_list[[m]]$p <- mean(x_list[[m]]$w)
    if (x_list[[m]]$p == 0) {
      x_list[[m]]$w <- rep(1, Nx)
    }
    x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
    x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p
  }
}

