#' @export
dHawkes <- function(theta){
  dunif(theta[1], 0.1, 0.5) * dunif(theta[2], 0.2, 1)
}

#' @export
loss_hawkes <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- rgamma(1, 10, 10)
  } else {
    x <- generate_state(x, 1, lower = inp$lower, upper = inp$upper, sd = inp$sd, a = inp$a)
  }

  output <- dist_h(x, theta, inp$history, time1 = time1, time2 = time2, simulator = inp$simulator)

  return(list(distance = output, x = x))
}

#' @export
dist_h <- function(x, theta, history, time1, time2, simulator){

  if(gtools::invalid(which(history <= time1))){
    history <- NULL
  }


  sim_out <- simulator(x, theta, history[which(history <= time1 & history >= time1 - 400)], time1, time2, Ni = Inf)

  dist_out <- (sim_out$n - length(which(history > time1 & history < time2)))^2 / max(1, length(which(history > time1 & history < time2)))

  return(dist_out)
}


#' Generate the state
#' @param init initial value
#' @param n number of states to generate
#' @param lower lower bound
#' @param upper upper bound
#' @param sd standard deviation of normal
#' @param a coefficient to apply
#' @export
#' @examples
#' init <- 2.5
#' x <- generate_state(init, 10, 0, 3.5, 1, 0.9)
generate_state <- function(init, n, lower, upper, sd, a = 1){

  if(gtools::invalid(init)){
    init <- rgamma(1, 10, 10)
    if(n == 1){
      return(init)
    }
  }

  one_step <- FALSE

  if(n == 1){
    one_step <- TRUE
    n <- 2
  }

  x <- rep(NA, n)
  x[1] <- init

  for(t in 2:n){
    x[t] <- gtools::inv.logit(a * gtools::logit(x[t-1], min = lower, max = upper) + rnorm(1, sd = sd), min = lower, max = upper)
  }

  if(one_step){
    x <- x[2]
  }

  return(x)
}
