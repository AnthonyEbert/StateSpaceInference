
#' @export
SMC2_ABC <- function(prior_sample, dprior, rstate, loss, loss_args, control, dt = 1, ESS_threshold = 0.1, eps = NULL, cl = NULL, lfunc = NULL, TT){

  if(is.null(lfunc)){
    parallel <- ifelse(is.null(cl), 1,
                       ifelse("cluster" %in% class(cl), 2,
                              ifelse(cl == "mclapply", 3,
                                     ifelse(cl == "test", 4, 5)))
    )

    lfunc <- make_listfunc(parallel, cl)
  }

  n_params <- ifelse(gtools::invalid(dim(prior_sample)[2]), 1, dim(prior_sample)[2])

  Ntheta <- control$Ntheta
  Nx <- control$Nx

  x_list <- as.list(rep(NA, Ntheta))
  x_list <- lapply(x_list, function(i){list()})

  if(is.null(eps)){
    eps <- rep(NA, TT)
  }

  for(m in 1:Ntheta){
    x_list[[m]]$theta <- prior_sample[m,, drop = FALSE]
    x_list[[m]]$x <- rstate(Nx, prior_sample[m,], loss_args)
    x_list[[m]]$w <- rep(1, Nx)
    x_list[[m]]$p <- NULL
    x_list[[m]]$omega <- 1
  }

  q_l <- list()


  for(t in 1:TT){
    x_list <- lfunc(x_list, theta_filter, time1 = t*dt - dt, time2 = t*dt, loss = loss, loss_args = loss_args)

    distances <- sapply(x_list, function(x){x$distance})
    eps[t] <- quantile(as.numeric(distances), probs = control$pacc)

    for(m in 1:Ntheta){
      x_list[[m]]$w <- (x_list[[m]]$distance <= eps[t])*1
      x_list[[m]]$p <- c(x_list[[m]]$p, mean(x_list[[m]]$w))
    }

    size_x <- ifelse(gtools::invalid(dim(x_list[[1]]$x)[2]), 1, dim(x_list[[1]]$x)[2])

    x_mat <- array(as.numeric(unlist(lapply(x_list, function(x){x$x}))), dim = c(Nx, Ntheta, size_x))
    w_mat <- array(as.numeric(unlist(lapply(x_list, function(x){x$w}))), dim = c(Nx, Ntheta))
    w_mat <- array(w_mat, dim = c(dim(w_mat), size_x))

    x_l <- lapply(1:size_x, function(i){x_mat[,,i]})
    w_l <- lapply(1:size_x, function(i){w_mat[,,i]})

    q_l[[t]] <- mapply(Hmisc::wtd.quantile, x_l, w_l, MoreArgs = list(probs = c(0.025, 0.5, 0.975), normwt = TRUE))

    for(m in 1:Ntheta){
      if (x_list[[m]]$p[t] == 0) {
        x_list[[m]]$w <- rep(1, Nx)
      }
      x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
      x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p[t]
    }

    omegas <- sapply(x_list, function(x){x$omega})
    thetas <- matrix(sapply(x_list, function(x){x$theta}), ncol = n_params)

    ESS <- sum(omegas) ^ 2 / sum(omegas^2)
    ifelse(ESS_threshold > 0, print(ESS), print(paste("SMC: ", ESS)))

    if(ESS < Ntheta * ESS_threshold){
      print("resample")
      post_cov <- cov.wt(log(thetas), wt=omegas)$cov
      nb <- 0
      aa <- sample(1:Ntheta, Ntheta, prob = omegas/sum(omegas), replace = TRUE)
      proposed_log_theta <- log(thetas) + mvtnorm::rmvnorm(Ntheta, sigma = post_cov)
      proposed_thetas <- exp(proposed_log_theta)

      x_list_prop <- SMC2_ABC(proposed_thetas, dprior = dprior, rstate = rstate, loss = loss, loss_args = loss_args, control = control, lfunc = lfunc, dt = dt, ESS_threshold = 0, eps = eps[1:t], TT = t)

      for(m in 1:Ntheta){
        proposed_Z_hat <- prod(x_list_prop[[m]]$p)
        old_Z_hat      <- prod(x_list[[aa[m]]]$p)

        MH_ratio <- dprior(proposed_thetas[m,])  * proposed_Z_hat /
          (dprior(thetas[aa[m],]) * old_Z_hat)

        un <- runif(1)

        if(un < MH_ratio){
          nb <- nb + 1
          x_list[[m]] <- x_list_prop[[m]]
        } else {
          x_list[[m]] <- x_list[[aa[m]]]
        }


        x_list[[m]]$omega <- 1
      }
      cat("acceptance rate: ", nb/Ntheta, "\n")

    }

  }

  x_list$q_l <- q_l

  return(x_list)
}


make_listfunc <- function(parallel, cl){
  output <- switch(parallel,
                   lapply,
                   function(X, FUN, ...){
                     return(parallel::parLapply(cl, X, FUN, ...))
                   },
                   parallel::mclapply
  )

  return(output)
}
