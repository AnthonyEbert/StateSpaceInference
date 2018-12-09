
library(parallel)
library(StateSpaceInference)
library(ggplot2)
library(ggalt)
sessionInfo()

#cl <- makeCluster(parallel::detectCores() - 1)
cl = "mclapply"


set.seed(3)

TT <- 40
true_theta <- c(0.25, 0.5)
lower <- 0
upper <- 3.5
sd_t <- 1
init <- min(rgamma(1, 100, 100), upper - 1)
a_logit <- 0.9
dist_coef <- 0.5
true_states <- generate_state(init, TT, lower, upper, sd_t, a = a_logit)

lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = true_states)
kern <- function(x){return(decay_func(x, alpha = true_theta[1], delta = true_theta[2]))}

y_history <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)
y <- hist(y_history, breaks = seq(0, TT*10, by = 10), plot = FALSE)$counts

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


Ntheta = 4000
Nx = 500
pacc = 0.005

prior_sample <- data.frame(theta1 = rgamma(Ntheta, 10, 40), theta2 = rgamma(Ntheta, 10, 20))

prior_sample <- as.matrix(prior_sample, ncol = 2)

full_list <- SMC2_ABC(prior_sample, dprior = dHawkes, loss, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 10, ESS_threshold = 0.1, TT = TT, trans = log, invtrans = exp, resample_times = 40)

state_df <- get_state(full_list)

state_df$state <- true_states

theta_df <- get_parameter(full_list)


save.image()
save(state_df, file = "state_df.RData")
save(theta_df, file = "theta_df.RData")



ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df) + aes(x = Value, weights = Weight, col = Time) + geom_density() + facet_wrap(~Parameter)
