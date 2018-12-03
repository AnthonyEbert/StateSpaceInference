
theta_df <- data.frame(
  theta = rep(c(0.5, 0.25), each = 2),
  init = runif(4, 0.5, 3)
)

Ntheta <- 4
Nx <- 10
x_list <- as.list(rep(NA, Ntheta))
x_list <- lapply(x_list, function(i){list()})

eps <- rep(NA, 5)

inp <- list(
  lower = lower,
  upper = upper,
  sd_t = sd_t,
  a_logit = a_logit,
  history = history,
  simulator = simulator,
  time1 = 10,
  time2 = 20
)

for(m in 1:4){
  x_list[[m]]$theta <- theta_df$theta[m]
  x_list[[m]]$x <- sapply(rep(NA, 10), generate_statey, n = 1, lower = lower, upper = upper, sd = sd_t, a = a_logit)
  x_list[[m]]$w <- 1
}

for(m in 1:4){
  output <- particle_filter(x_list[[m]]$theta, x_list[[m]]$x, w = x_list$w, loss = loss_hawkes, loss_args = inp, threshold = 2, lfunc = lfunc, time1 = 10, time2 = 20)

  x_list[[m]]$distance <- sapply(output, function(x){x$distance})

  x_list[[m]]$x <- sapply(output, function(x){x$x})
}

  distances <- sapply(x_list, function(x){x$distance})

  eps <- quantile(as.numeric(distances), probs = 0.2)
