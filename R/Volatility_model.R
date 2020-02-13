# Volatility model

#' @export
loss_volatility <- function(x, theta, time1, time2, inp){
  phi <- theta[1]

  if(gtools::invalid(x)){
    x <- rnorm(1, inp$mu/(1-phi), inp$s_h^2/(1-phi^2))
  } else {
    x <- inp$mu + phi * x + inp$s_h * rnorm(1)
  }

  if(is.null(inp$gamma)){
    inp$gamma <- 0
  }

  y_prop <- exp(x/2) * stable(1, inp$alpha, inp$beta, inp$gamma, inp$s_v)
  output <- abs(inp$y[time2] - y_prop)

  return(list(distance = output, x = x))
}

#' stable distribution
#' @export
#' @references
#' Vankov, E., & Ensor, K. B. (2014). Stochastic Volatility Filtering with Intractable Likelihoods. arXiv:1405.4323.
stable <- function(n, alpha, beta, gamma, delta) {
  U <- runif(n, min = -.5 * pi, max = .5 * pi)
  W <- rexp(n)
  if (dplyr::near(alpha, 1.)) {
    X <- 2 / pi * ((0.5 * pi + beta * U) * tan(U) -
                     beta * log(pi * W * cos(U) / (0.5*pi + beta * W)))
    Y <- gamma * X + (delta + beta * 2 / pi * gamma * log(gamma))
  } else {
    psi <- atan(beta * tan(pi * alpha * .5))
    X <- sin(alpha * (psi + U)) / (cos(alpha * psi) * cos(U))^(1/alpha) *
      (cos(alpha * psi + (1 - alpha) * U) / W) ^ ((1 - alpha) / alpha)
    Y <- gamma * X + delta
  }
  return(Y)
}

#' @export
rprior_vol <- function(N) {
  phi <- runif(N, min = -1, max = +1)
  return(phi)
}

#' @export
dprior_vol <- function(theta){
  dunif(theta[1], -1, 1)
}
