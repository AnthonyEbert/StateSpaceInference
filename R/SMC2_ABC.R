
#' @export
SMC2_ABC <- function(prior_sample, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp = 1, ESS_threshold = 0.1, eps = NULL, cl = NULL, TT, trans = I, invtrans = I, resample_times = NA, trans_args = list(), cov_coef = 1){


  parallel <- ifelse(is.null(cl), 1,
                     ifelse("cluster" %in% class(cl), 2,
                            ifelse(cl == "mclapply", 3,
                                   ifelse(cl == "test", 4, 5)))
  )

  lfunc <- make_listfunc(parallel, cl)

  n_params <- ifelse(gtools::invalid(dim(prior_sample)[2]), 1, dim(prior_sample)[2])

  x_list <- as.list(rep(NA, Ntheta))
  x_list <- lapply(x_list, function(i){list()})

  full_list <- list()
  if(is.null(eps)){
    eps <- rep(NA, TT)
  }

  for(m in 1:Ntheta){
    x_list[[m]]$theta <- prior_sample[m, , drop = FALSE]
    x_list[[m]]$x <- matrix(NA, nrow = Nx)
    x_list[[m]]$w <- rep(1, Nx)
    x_list[[m]]$p <- NULL
    x_list[[m]]$pprod <- NULL
    x_list[[m]]$omega <- 1
  }

  q_l <- list()


  for(tp in 1:TT){
    x_list <- lfunc(x_list, theta_filter, time1 = tp*dtp - dtp, time2 = tp*dtp, loss = loss, loss_args = loss_args)

    distances <- sapply(x_list, function(x){x$distance})
    if(tp > 1){
      d_weights <- sapply(full_list[[tp - 1]], function(x){rep(x$omega, Nx)})
    } else {
      d_weights <- rep(1, length(distances))
    }

    if(is.na(eps[tp])){
      eps[tp] <- Hmisc::wtd.quantile(distances, weights = d_weights/sum(d_weights), probs = pacc, normwt = TRUE)
      #eps[tp] <- quantile(as.numeric(distances), probs = pacc)
    }

    for(m in 1:Ntheta){
      x_list[[m]]$w     <- (x_list[[m]]$distance <= eps[tp])*1
      x_list[[m]]$pprod <- prod(x_list[[m]]$pprod, mean(x_list[[m]]$w))
      x_list[[m]]$p     <- mean(x_list[[m]]$w)

      if (x_list[[m]]$p == 0) {
        x_list[[m]]$w <- rep(1, Nx)
      }
      x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
      x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p
    }

    omegas <- sapply(x_list, function(x){x$omega})
    thetas <- matrix(t(sapply(x_list, function(x){x$theta})), ncol = n_params)

    ESS <- sum(omegas) ^ 2 / sum(omegas^2)
    for(i in 1:dim(thetas)[2]){
      print(as.numeric(Hmisc::wtd.quantile(thetas[,i], weights = omegas/sum(omegas), normwt = TRUE)))
    }
    ifelse(ESS_threshold > 0, print(paste0(tp, ". ", ESS)), print(paste("SMC: ", ESS)))

    full_list[[tp]] <- x_list

    if(ESS < Ntheta * ESS_threshold | tp %in% resample_times){
      print("resample")
      post_cov <- cov_coef * cov.wt(trans(thetas, trans_args = trans_args), wt=omegas)$cov
      nb <- 0
      aa <- sample(1:Ntheta, Ntheta, prob = omegas/sum(omegas), replace = TRUE)
      proposed_log_theta <- trans(thetas[aa, , drop = FALSE], trans_args = trans_args) + mvtnorm::rmvnorm(Ntheta, sigma = post_cov)
      proposed_thetas <- invtrans(proposed_log_theta, trans_args = trans_args)

      full_list_prop <- SMC2_ABC(proposed_thetas, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp, ESS_threshold = 0, eps = eps[1:tp], cl,  TT = tp, trans = trans, invtrans = invtrans)

      x_list_prop <- full_list_prop[[tp]]

      proposed_omegas <- sapply(x_list_prop, function(x){x$omega})
      p_aa <- sample(1:Ntheta, Ntheta, prob = proposed_omegas/sum(proposed_omegas), replace = TRUE)

      for(m in 1:Ntheta){
        proposed_Z_hat <- x_list_prop[[p_aa[m]]]$pprod
        old_Z_hat      <- x_list[[aa[m]]]$pprod

        MH_ratio <- proposed_Z_hat * dprior(proposed_thetas[p_aa[m],]) / (old_Z_hat * dprior(thetas[aa[m],]))

        un <- runif(1)

        if(un < MH_ratio){
          nb <- nb + 1
          for(time_star in 1:tp){
            full_list[[time_star]][[m]] <- full_list_prop[[time_star]][[p_aa[m]]]
          }
        } else {
          for(time_star in 1:tp){
            full_list[[time_star]][[m]] <- full_list[[time_star]][[aa[m]]]
          }
        }

        x_list[[m]] <- full_list[[tp]][[m]]

        x_list[[m]]$omega <- 1
      }
      cat("acceptance rate: ", nb/Ntheta, "\n")

    }


    full_list[[tp]] <- x_list

  }

  return(full_list)
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
