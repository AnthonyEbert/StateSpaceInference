
library(StateSpaceInference)
sessionInfo()
set.seed(3)

# TT <- 60
# true_theta <- c(0.25, 0.5)
# lower <- 0
# upper <- 3.5
# sd_t <- 1
# init <- min(rgamma(1, 100, 100), upper - 1)
# a_logit <- 0.9
# dist_coef <- 0.5
# true_states <- generate_state(init, TT, lower, upper, sd_t, a = a_logit)
#
# lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = true_states)
# kern <- function(x){return(decay_func(x, alpha = true_theta[1], delta = true_theta[2]))}
#
# y_history <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)

y_fun <- function(theta1, theta2, states = NULL){
  if(is.null(states)){
    states <- generate_state(NULL, 4, 0, 3.5, sd = 1, 0.9)
  }

  theta <- c(theta1, theta2)
  TT <- length(states)
  lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = states)
  kern <- function(x){return(decay_func(x, alpha = theta[1], delta = theta[2]))}
  output <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)

  output <- output[which(output >= TT * 10 - 10)]

  n_events <- length(output)

  skew   <- e1071::skewness(output - (TT*10 - 10))
  if(is.nan(skew)){skew <- 0}
  output <- c(TT*10 - 10, output, TT * 10)
  diffs2    <- sum((diff(output))^2)
  diffs3     <- sum((diff(output))^3)
  min_diff <- min(diff(output))

  return(c(state = states[TT], n = n_events, diffs2 = diffs2, diffs3 = diffs3, min_diff = min_diff, skew = skew))
}

## Theta1

t1 <- seq(0.1, 0.4, by = 1e-4)
t1_mat <- t(mapply(y_fun, theta1 = t1, theta2 = 0.5))
t1_df <- cbind.data.frame(theta1 = t1, theta2 = 0.5, t1_mat)
t1_lm <- lm(theta1 ~ poly(n,6, raw = TRUE) + diffs2 + diffs3 + min_diff, data = t1_df)
summary(t1_lm)
attr(t1_lm$terms, "predvars")
as.numeric(coef(t1_lm))

state1_lm <- lm(state ~ poly(n,6, raw = TRUE) + diffs2 + diffs3 + min_diff, data = t1_df)
summary(state1_lm)

t2 <- seq(0.2, 1, by = 1e-4)
t2_mat <- t(mapply(y_fun, theta1 = 0.25, theta2 = t2))
t2_df <- cbind.data.frame(theta1 = 0.25, theta2 = t2, t2_mat)
t2_lm <- lm(theta2 ~ poly(n,6, raw = TRUE) + diffs2 + diffs3 + min_diff, data = t2_df)
summary(t2_lm)
attr(t2_lm$terms, "predvars")
as.numeric(coef(t2_lm))

state2_lm <- lm(state ~ poly(n,6, raw = TRUE) + diffs2 + diffs3 + min_diff, data = t2_df)
summary(state2_lm)
attr(state2_lm$terms, "predvars")
as.numeric(coef(state2_lm))


theta_mat <- expand.grid(seq(0.1, 0.4, by = .0025), seq(0.3, 0.8, by = .0025))
theta_mat2 <- t(mapply(y_fun, theta1 = theta_mat[,1], theta2 = theta_mat[,2]))
theta_df <- cbind.data.frame(theta1 = theta_mat[,1], theta2 = theta_mat[,2], theta_mat2)
theta_lm <- lm(cbind(theta1, theta2, state) ~ log(n + 1) + log(diffs2 + 1e-6) + log(diffs3 + 1e-6)  + skew, data = theta_df)

theta1_lm <-
