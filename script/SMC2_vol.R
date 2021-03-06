
library(StateSpaceInference)
library(parallel)
library(ggplot2)
library(ggalt)
library(dplyr)

sessionInfo()

seed <- 14
set.seed(seed)

cl <- makeCluster(parallel::detectCores())
#cl = "mclapply"
#cl <- NULL

# length of the time series
TT <- 20
# parameters
alpha <- 2; beta <- 0; gamma <- 0.1 * sqrt(1/2); mu <- 1; phi <- 0.80; sh <- 0.6; s_v <- 1
# simulating the hidden states
h <- rep(0, TT)
h[1] <- rnorm(1, mu/(1-phi), sd = sqrt(sh^2/(1-phi^2)))
for (t in 2:TT) {
  h[t] <- mu + phi * (h[t - 1]) + sh * rnorm(1)
}

# emission of the observations
yobs <- exp(h/2) * stabledist::rstable(TT, alpha, beta, gamma, s_v)


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

Ntheta <- 20
Nx <- 2500
pacc = 0.05

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

full_list <- SMC2_ABC(prior_sample, dprior = dprior_vol, loss = loss_volatility, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 1, ESS_threshold = 0.5, TT = TT, trans = trans, invtrans = invtrans, cov_coef = 1.5, acceptance_correction = acceptance_correction)


state_df <- get_state(full_list, probs = c(0.025, 0.5, 0.975))

state_df$state <- true_states

theta_df <- get_parameter(full_list)

ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")

theta_df <- theta_df %>%
  filter(time == TT) %>%
  mutate(seed = seed, type = "ABC") %>%
  select(-parameter, -time)

save.image()
#save(state_df, file = "state_df.RData")
saveRDS(theta_df, file = paste0("theta_df_", seed,".RData"))
