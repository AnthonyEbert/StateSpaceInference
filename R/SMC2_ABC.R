
#' @export
SMC2_ABC <- function(prior_sample, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp = 1, ESS_threshold = 0.1, eps = NULL, cl = NULL, TT, trans = I, invtrans = I, size_x = 1){


  parallel <- ifelse(is.null(cl), 1,
                     ifelse("cluster" %in% class(cl), 2,
                            ifelse(cl == "mclapply", 3,
                                   ifelse(cl == "test", 4, 5)))
  )

  lfunc <- make_listfunc(parallel, cl)

  n_params <- ifelse(gtools::invalid(dim(prior_sample)[2]), 1, dim(prior_sample)[2])

  x_list <- as.list(rep(NA, Ntheta))
  x_list <- lapply(x_list, function(i){list()})

  if(is.null(eps)){
    eps <- rep(NA, TT)
  }

  for(m in 1:Ntheta){
    x_list[[m]]$theta <- prior_sample[m, , drop = FALSE]
    x_list[[m]]$x <- matrix(nrow = Nx, ncol = size_x)
    x_list[[m]]$w <- rep(1, Nx)
    x_list[[m]]$p <- NULL
    x_list[[m]]$omega <- 1
  }

  q_l <- list()


  for(tp in 1:TT){
    x_list <- lfunc(x_list, theta_filter, time1 = tp*dtp - dtp, time2 = tp*dtp, loss = loss, loss_args = loss_args)

    distances <- sapply(x_list, function(x){x$distance})
    if(is.na(eps[tp])){
      eps[tp] <- quantile(as.numeric(distances), probs = pacc)
    }

    for(m in 1:Ntheta){
      x_list[[m]]$w <- (x_list[[m]]$distance <= eps[tp])*1
      x_list[[m]]$p <- c(x_list[[m]]$p, mean(x_list[[m]]$w))
    }

    # Saving quantiles of x for plotting only -------
    size_x <- ifelse(gtools::invalid(dim(x_list[[1]]$x)[2]), 1, dim(x_list[[1]]$x)[2])

    x_mat <- aperm(abind::abind(lapply(x_list, function(x){x$x}), along = 3), c(3,1,2))

    #x_mat <- t(sapply(lapply(x_list, function(x){x$x}), I))
    w_mat <- t(sapply(lapply(x_list, function(x){x$w}), I))
    w_mat <- array(w_mat, dim = c(dim(w_mat), size_x))

    x_l <- lapply(1:Ntheta, function(i){abind::adrop(x_mat[i,,, drop = FALSE], drop = 1)})
    w_l <- lapply(1:Ntheta, function(i){w_mat[i,]})

    q_mean <- mapply(function(x,y){apply(x, 2, weighted.mean, w = y)}, x_l, w_l)

    q_l[[tp]] <- mapply(Hmisc::wtd.quantile, x_l, w_l, MoreArgs = list(probs = c(0.025, 0.5, 0.975), normwt = TRUE))

    # End state quantile saving ------------

    for(m in 1:Ntheta){
      if (x_list[[m]]$p[tp] == 0) {
        x_list[[m]]$w <- rep(1, Nx)
      }
      x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
      x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p[tp]
    }

    omegas <- sapply(x_list, function(x){x$omega})
    thetas <- matrix(sapply(x_list, function(x){x$theta}), ncol = n_params)

    ESS <- sum(omegas) ^ 2 / sum(omegas^2)
    ifelse(ESS_threshold > 0, print(paste0(tp, ". ", ESS)), print(paste("SMC: ", ESS)))

    if(ESS < Ntheta * ESS_threshold){
      print("resample")
      post_cov <- cov.wt(trans(thetas), wt=omegas)$cov
      nb <- 0
      aa <- sample(1:Ntheta, Ntheta, prob = omegas/sum(omegas), replace = TRUE)
      proposed_log_theta <- trans(thetas[aa, , drop = FALSE]) + mvtnorm::rmvnorm(Ntheta, sigma = post_cov)
      probs <- apply(invtrans(proposed_log_theta), 1, dprior)
      probs[is.na(probs)] <- 0
      bb <- sample(1:Ntheta, Ntheta, replace = TRUE, prob = probs)
      proposed_thetas <- invtrans(proposed_log_theta[bb, , drop = FALSE])

      x_list_prop <- SMC2_ABC(proposed_thetas, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp, ESS_threshold = 0, eps = eps[1:tp], cl,  TT = tp, trans = trans, invtrans = invtrans)

      for(m in 1:Ntheta){
        proposed_Z_hat <- prod(x_list_prop[[m]]$p)
        old_Z_hat      <- prod(x_list[[aa[m]]]$p)

        MH_ratio <- proposed_Z_hat / old_Z_hat

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

  attr(x_list, "q_l") <- q_l

  return(x_list)
}

#' @export
make_listfunc <- function(parallel, cl){
  output <- switch(parallel,
                   lapply,
                   function(X, fun, ...){
                     return(parallel::parLapply(cl, X, fun, ...))
                   },
                   parallel::mclapply
  )

  return(output)
}
