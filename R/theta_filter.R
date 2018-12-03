
#' @export
theta_filter <- function(input_list, time1, time2, loss_args){

  output <- particle_filter(input_list$theta, input_list$x, input_list$w, loss = loss_hawkes, loss_args = loss_args, time1 = time1, time2 = time2)

  output_list <- input_list

  output_list$distance <- sapply(output, function(x){x$distance})

  output_list$x <- sapply(output, function(x){x$x})

  return(output_list)
}

#' @export
weighter_filter <- function(input_list, eps){

  output <- lapply(input_list$distances)

  distances <- input_list$distance

  output <- (distances <= eps) * 1

  input_list$w <- output

  return(input_list)

}
