
library(StateSpaceInference)
library(parallel)
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

Ntheta <- 50
Nx <- 100
pacc = 0.05

prior_sample <- data.frame(rprior_vol(Ntheta))

prior_sample <- as.matrix(prior_sample, ncol = 1)

x_list <- SMC2_ABC(prior_sample, dprior = dprior_vol, loss = loss_volatility, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 1, ESS_threshold = 0.1, TT = TT, trans = I, invtrans = I)

plotrix::weighted.hist(sapply(x_list, function(x){x$theta}), w = sapply(x_list, function(x){x$omega})/sum(sapply(x_list, function(x){x$omega})))

plot(density(sapply(x_list, function(x){x$theta}), weights = sapply(x_list, function(x){x$omega})/sum(sapply(x_list, function(x){x$omega}))))

q_mat <- array(unlist(attr(x_list, "q_l"), recursive = FALSE), dim = c(3, 20))

q_df <- data.frame(time = seq(1, TT), lower = q_mat[1,], med = q_mat[2,], upper = q_mat[3,], state = h)

library(ggalt)

ggplot(q_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
