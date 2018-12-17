
#' @export
particle_filter <- function(theta, x, w, loss, loss_args, threshold, time1, time2, lfunc){
  a <- sample(1:I(dim(x)[1]), prob = w, replace = TRUE)
  dist_out <- apply(x[a, , drop = FALSE], MARGIN = 1, loss, theta = theta, inp = loss_args, time1 = time1, time2 = time2)
  return(dist_out)
}

#' @export
particle_filter2 <- function(theta, x, w, loss, loss_args, threshold, time1, time2, lfunc){

  prior <- function(n){
    a <- sample(1:I(dim(x)[1]), size = n, prob = w, replace = TRUE)
    x_out <- apply(x[a, , drop = FALSE], MARGIN = 1, generate_state_simple, theta = theta, inp = loss_args, time1 = time1, time2 = time2)
    out <- data.frame(x_out)
    return(out)
  }

  orig_value <- x

  prior_eval <- function(z, inp = loss_args){
    z_logit <- gtools::logit(z, min = inp$lower, max = inp$upper)
    if(is.na(orig_value[1])){
      orig_value <- apply(orig_value, 1, generate_state, n = 1, lower = inp$lower, upper = inp$upper, sd = inp$sd_t)
    }
    orig_logit <- gtools::logit(orig_value, min = inp$lower, max = inp$upper)


    output <- sum(w * dnorm(z_logit, as.numeric(orig_logit), sd = inp$sd_t))/sum(w)
    return(output)
  }


  loss_fun <- function(x, inp){
    dist_out <- loss_state_simple(x, theta = theta, time1 = time1, time2 = time2, inp = inp)
    return(dist_out)
  }

  output <- protoABC::abc_start(prior, loss_fun, distance_args = loss_args, method = "RABC", control = list(n = dim(x)[1], pacc_final = 0.2, prior_eval = prior_eval), output_control = list(include_dist = TRUE))

  out_list <- as.list(rep(NA, dim(output)[1]))
  out_list <- lapply(out_list, function(i){list()})

  for(i in 1:dim(output)[1]){
    out_list[[i]]$distance <- output[i,2]
    out_list[[i]]$x        <- output[i,1]
  }

  return(out_list)
}


