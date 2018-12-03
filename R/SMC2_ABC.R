
#' @export
SMC2_ABC <- function(prior_sample, rstate, loss, loss_args, control, dt = 1, cl = NULL){

  parallel <- ifelse(is.null(cl), 1,
                     ifelse("cluster" %in% class(cl), 2,
                            ifelse(cl == "mclapply", 3,
                                   ifelse(cl == "test", 4, 5)))
  )

  lfunc <- make_listfunc(parallel, cl)

  Ntheta <- control$Ntheta
  Nx <- control$Nx

  x_list <- as.list(rep(NA, Ntheta))
  x_list <- lapply(x_list, function(i){list()})

  eps <- rep(NA, control$TT)

  for(m in 1:Ntheta){
    x_list[[m]]$theta <- prior_sample[m]
    x_list[[m]]$x <- rstate(Nx, prior_sample[m], loss_args)
    x_list[[m]]$w <- rep(1, Nx)
    x_list[[m]]$p <- 1
    x_list[[m]]$omega <- 1
  }


  for(t in 2:control$TT){
    x_list <- lfunc(x_list, theta_filter, time1 = t*dt - dt, time2 = t*dt, loss = loss, loss_args = inp)

    distances <- sapply(x_list, function(x){x$distance})
    eps[t] <- quantile(as.numeric(distances), probs = control$pacc)

    for(m in 1:Ntheta){
      x_list[[m]]$w <- (x_list[[m]]$distance <= eps[t])*1
      x_list[[m]]$p <- mean(x_list[[m]]$w)
      if (x_list[[m]]$p == 0) {
        x_list[[m]]$w <- rep(1, Nx)
      }
      x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
      x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p
    }
  }

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
