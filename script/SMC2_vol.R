
library(StateSpaceInference)
library(parallel)
library(ggplot2)
library(ggalt)

sessionInfo()
set.seed(1)

#cl <- makeCluster(parallel::detectCores() - 1)
cl = "mclapply"
#cl <- NULL

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


true_states <- h

inp <- list(
  alpha = alpha,
  beta = beta,
  gamma = gamma,
  mu = mu,
  s_h = sh,
  s_v = s_v,
  y = yobs
)

Ntheta <- 80
Nx <- 50000
pacc = 0.005

prior_sample <- data.frame(rprior_vol(Ntheta))

prior_sample <- as.matrix(prior_sample, ncol = 1)

trans <- function(x, trans_args){
  theta1 <- qnorm((x + 1)/2)
  return(theta1)
}

invtrans <- function(x, trans_args){
  theta1 <- 2*pnorm(x) - 1
  return(theta1)
}

acceptance_correction <- function(x){
  0.5/(dnorm(qnorm((x+1)/2)))
}

full_list <- SMC2_ABC(prior_sample, dprior = dprior_vol, loss = loss_volatility, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 1, ESS_threshold = 0.5, TT = TT, trans = trans, invtrans = invtrans, cov_coef = 0.5, acceptance_correction = acceptance_correction)


state_df <- get_state(full_list, probs = c(0.25, 0.5, 0.75))

state_df$state <- true_states

theta_df <- get_parameter(full_list)

save.image()
save(state_df, file = "state_df.RData")
save(theta_df, file = "theta_df.RData")

ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")
