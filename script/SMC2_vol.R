
library(StateSpaceInference)
library(parallel)
library(ggplot2)
library(ggalt)

sessionInfo()

cl <- makeCluster(parallel::detectCores() - 1)

# length of the time series
TT <- 20
# parameters
alpha <- 1.75; beta <- 0.1; mu <- -0.2; phi <- 0.95; s_h <- 0.6; s_v <- 0.8
# simulating the hidden states
h <- rep(0, TT)
h[1] <- rnorm(1)
for (t in 2:TT) {
  h[t] <- mu + phi * h[t - 1] + s_h * rnorm(1)
}

true_states <- h

# emission of the observations
yobs <- exp(h/2) * stable(TT, alpha, beta, 0, s_v)


inp <- list(
  alpha = alpha,
  beta = beta,
  mu = mu,
  s_h = s_h,
  s_v = s_v,
  y = yobs
)

Ntheta <- 80
Nx <- 5000
pacc = 0.005

prior_sample <- data.frame(rprior_vol(Ntheta))

prior_sample <- as.matrix(prior_sample, ncol = 1)

trans <- function(x, trans_args){
  theta1 <- qnorm((x + 1)/2)
  return(theta1)
}

invtrans <- function(x, trans_args){
  theta1 <- pnorm(2*x - 1)
  return(theta1)
}

acceptance_correction <- function(x){
  0.5/(dnorm(qnorm((x+1)/2)))
}

full_list <- SMC2_ABC(prior_sample, dprior = dprior_vol, loss = loss_volatility, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 1, ESS_threshold = 0.5, TT = TT, trans = trans, invtrans = invtrans, cov_coef = 0.5, acceptance_correction = acceptance_correction)


state_df <- get_state(full_list, probs = c(0.25, 0.5, 0.75))

state_df$state <- true_states

theta_df <- get_parameter(full_list)

ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")
