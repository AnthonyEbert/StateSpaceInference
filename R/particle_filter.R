
#' @export
particle_filter <- function(theta, x, w, loss, loss_args, threshold, time1, time2, lfunc){
  a <- sample(1:I(dim(x)[1]), prob = w, replace = TRUE)
  dist_out <- apply(x[a, , drop = FALSE], MARGIN = 1, loss, theta = theta, inp = loss_args, time1 = time1, time2 = time2)
  return(dist_out)
}

#' @export
particle_filter2 <- function(theta, x, w, loss, loss_args, threshold, time1, time2, lfunc){

  prior <- function(n, inp = loss_args){
    a <- sample(1:I(dim(x)[1]), size = n, prob = w, replace = TRUE)
    x_out <- apply(x[a, , drop = FALSE], MARGIN = 1, generate_state, lower = inp$lower, upper = inp$upper, sd = inp$sd_t, n = 1)

    z_out <- gtools::logit(x_out, min = inp$lower, max = inp$upper)
    out <- data.frame(z_out)
    return(out)
  }

  orig_value <- x

  prior_eval <- function(z, inp = loss_args){
    # if(z < inp$lower | z > inp$upper){
    #   return(0)
    # }
    #
    # z_logit <- gtools::logit(z, min = inp$lower, max = inp$upper)
    #
    # if(inp$kern){
    #   z_logit <- x
    # }

    z_logit <- z

    if(is.na(orig_value[1])){
      #orig_value <- apply(orig_value, 1, generate_state, n = 1, lower = inp$lower, upper = inp$upper, sd = inp$sd_t)
      orig_value <- rep(1, length(orig_value))
    }
    orig_logit <- gtools::logit(orig_value, min = inp$lower, max = inp$upper)


    output <- sum(w * dnorm(z_logit, as.numeric(orig_logit), sd = inp$sd_t))/sum(w)

    if(gtools::invalid(output)){
      output <- 0
    }
    return(output)
  }


  loss_fun <- function(x, inp = loss_args){
    z <- gtools::inv.logit(x, min = inp$lower, max = inp$upper)
    dist_out <- loss_hawkes_simple(z, theta = theta, time1 = time1, time2 = time2, inp = inp)
    return(dist_out)
  }

  kern_loss <- function(x, inp = loss_args){
    # z <- gtools::inv.logit(x, min = inp$lower, max = inp$upper)
    dist_out <- loss_fun(x, inp = inp)
    dnorm(dist_out, 0, 1)
  }

  output <- protoABC::abc_start(prior, loss_fun, distance_args = loss_args, method = "RABC", control = list(n = dim(x)[1], pacc_final = 0.05, eps_final = 1, prior_eval = prior_eval), output_control = list(include_dist = TRUE))

  loss_args$kern <- TRUE
  if(time1 > 2){
    kol <- 2+2
  }

  evidence <- protoABC::evidence(prior_eval, kern_loss, output[,-I(dim(output)[2]), drop = FALSE], kern_args = loss_args, control = list(cov_func = var), rw_cov = matrix(0.1))

  out_list <- as.list(rep(NA, dim(output)[1]))
  out_list <- lapply(out_list, function(i){list()})

  for(i in 1:dim(output)[1]){
    out_list[[i]]$distance <- output[i,2]
    out_list[[i]]$x        <- gtools::inv.logit(output[i,1], min = loss_args$lower, max = loss_args$upper)
    out_list[[i]]$evidence <- mean(evidence)
  }

  return(out_list)
}


