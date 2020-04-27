
#' @export
#' @importFrom dplyr %>%
get_state <- function(full_list, dim = 1, probs = c(0.025, 0.5, 0.975), cl = NULL){


  TT <- length(full_list)

  Ntheta <- length(full_list[[2]])
  Nx <- dim(full_list[[2]][[1]]$x)[1]

  q_mat <- list()

  if(is.na(probs[1])){
    x_array <- array(dim = c(TT, Ntheta, Nx))
    w_array <- array(dim = c(TT, Ntheta, Nx))
  }

  output <- parallel::parLapply(cl, full_list, function(x_list){

    x_mat <- t(abind::abind(lapply(x_list, function(x){abind::abind(x$x, along = 2)})))

    w_mat <- t(sapply(x_list, function(x){x$w * x$omega}))

    if(is.na(probs[1])){
      x_array[tp,,] <- x_mat
      w_array[tp,,] <- w_mat
    }

    if(!is.na(probs[1])){
      #x <- as.numeric(x_mat)
      w <- as.numeric(w_mat)

      # state_sample <- x_mat[indicies,]

      q_mat <- apply(x_mat, 2, Hmisc::wtd.quantile, weights = w/sum(w) * 1e6, probs = probs, normwt = FALSE)

      #q_mat[tp,] <- Hmisc::wtd.quantile(x, weights = w / sum(w), probs = probs, normwt = FALSE)
    }

    return(q_mat)
  })

  # if(!is.na(probs[1])){
  #   output <- data.frame(time = seq(1, TT), lower = q_mat[,1], med = q_mat[,2], upper = q_mat[,3])
  # } else {
  #   output <- list(x_array = x_array, w_array = w_array)
  # }

  return(output)
}
