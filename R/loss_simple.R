#' @export
loss_simple <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- rgamma(1, 10, 10)
  } else {
    x <- generate_state(x, 1, lower = inp$lower, upper = inp$upper, sd = inp$sd, a = inp$a)
  }

  y <- rnorm(1000, x, theta[1])
  output <- sum(
    (mean(y) - mean(inp$y[[time2]]))^2 +
    (4*(sd(y) - sd(inp$y[[time2]])))^2
  )

  #output <- dist_h(x, theta, inp$history, time1 = time1, time2 = time2, simulator = inp$simulator)

  return(list(distance = output, x = x))
}

#' @export
generate_simple <- function(TT, true_states, theta){
  y <- list()
  for(i in 1:TT){
    y[[i]] <- rnorm(1000, true_states[i], theta)
  }
  return(y)
}

