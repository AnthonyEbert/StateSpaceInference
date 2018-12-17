#' @export
dHawkes <- function(theta){
  dunif(theta[1], 0.1, 0.5) * dunif(theta[2], 0.2, 1)
}

#' @export
loss_hawkes <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- gtools::inv.logit(inp$a * gtools::logit(inp$lower + 0.5 * (inp$upper - inp$lower), min = inp$lower, max = inp$upper) + rnorm(1, sd = inp$sd), min = inp$lower, max = inp$upper)
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

  sim_out <- simulator(x, theta, history[which(history <= time1 & history >= time1 - 100)], time1, time2, Ni = Inf)

  cotemp_history <- history[which(history > time1 & history < time2)]

  n_out <- (sim_out$n - length(cotemp_history))^2 / max(1, length(cotemp_history))

  true_input <- sort(c(cotemp_history, time2, time1))
  sim_input <- sort(c(sim_out$new_history, time2, time1))

  diff_out <- mean(diff(true_input[which(!is.na(true_input))])) - mean(diff(sim_input[which(!is.na(sim_input))]))

  if(is.nan(diff_out)){
    kol <- 2+2
  }

  dist_out <- n_out + abs(diff_out)


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
