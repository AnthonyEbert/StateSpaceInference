#' @export
loss_simple <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- rgamma(1, 10, 10)
  } else {
    x <- generate_state(x, 1, lower = inp$lower, upper = inp$upper, sd = inp$sd, a = inp$a)
  }

  y <- sn::rsn(1e1, dp = sn::cp2dp(cp = c(x, theta), "SN"))

  ss_obs <- c(mean(inp$y[[time2]]), sd(inp$y[[time2]]), e1071::skewness(inp$y[[time2]]))
  ss_sim <- c(mean(y), sd(y), e1071::skewness(y))

  # ss_obs <- sn::summary(sn::selm(inp$y[[time2]]~1, opt.method = "BFGS"))@param.table[,1]
  # ss_sim <- sn::summary(sn::selm(y~1), opt.method = "BFGS")@param.table[,1]

  return(list(distance = sqrt(sum((ss_obs - ss_sim)^2)), x = x))
}

#' @export
generate_simple <- function(TT, true_states, theta){
  y <- list()
  for(i in 1:TT){
    y[[i]] <- as.numeric(sn::rsn(1e1, dp = sn::cp2dp(cp = c(true_states[i], theta), "SN")))
  }
  return(y)
}

