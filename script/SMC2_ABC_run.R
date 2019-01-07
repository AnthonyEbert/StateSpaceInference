
library(parallel)
library(StateSpaceInference)
library(ggplot2)
library(ggalt)
sessionInfo()

#cl <- makeCluster(parallel::detectCores())
cl = "mclapply"
#cl <- NULL

set.seed(3)

TT <- 60
true_theta <- c(0.5, 0.5)
lower <- 0
upper <- 3.5
sd_t <- 1
init <- min(rgamma(1, 100, 100), upper - 1)
a_logit <- 0.9
dist_coef <- 0.5
true_states <- generate_state(init, TT, lower, upper, sd_t, a = a_logit)

y <- NULL

#y <- hawkes_simulator(true_states[1], true_theta, NULL, 0, 10)
for(tp in 1:TT){
  y <- hawkes_simulator(true_states[tp], true_theta, y$history, tp * 10 - 10, tp * 10)
}

y_history <- y$history

lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = true_states)
kern <- function(x){return(decay_func(x, alpha = true_theta[1], delta = true_theta[2]))}
#
#y_history <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)
#y <- hist(y_history, breaks = seq(0, TT*10, by = 10), plot = FALSE)$counts

hist(y_history, breaks = TT * 10)

plot(lambda_fun, add = TRUE, col = "red")

simulator <- hawkes_simulator

inp <- list(
  lower = lower,
  upper = upper,
  sd_t = sd_t,
  a_logit = a_logit,
  history = y_history,
  simulator = simulator
)

loss = loss_hawkes


Ntheta = 2000
Nx = 1000
pacc = 0.5

lower_theta <- c(0.3, 0.3)
upper_theta <- c(0.7, 0.7)

prior_sample <- data.frame(theta1 = runif(Ntheta, lower_theta[1], upper_theta[1]), theta2 = runif(Ntheta, lower_theta[2], upper_theta[2]))

prior_sample <- as.matrix(prior_sample, ncol = 2)

trans <- function(x, trans_args){
  theta1 <- log(x[,1])
  theta2 <- log(x[,2])
  return(cbind(theta1, theta2))
}

invtrans <- function(x, trans_args){
  theta1 <- exp(x[,1])
  theta2 <- exp(x[,2])
  return(cbind(theta1, theta2))
}

full_list <- SMC2_ABC(prior_sample, dprior = dHawkes, loss, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 10, ESS_threshold = 0.5, TT = TT, trans = trans, invtrans = invtrans, cov_coef = 0.25^2)

state_df <- get_state(full_list, probs = c(0.25, 0.5, 0.75))

state_df$state <- true_states

theta_df <- get_parameter(full_list)


save.image()	
save(state_df, file = "state_df.RData")	
save(theta_df, file = "theta_df.RData")



ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")
