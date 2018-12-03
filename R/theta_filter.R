
theta_filter <- function(input_list, time1, time2, loss_args){
  output <- particle_filter(input$theta, input$x, loss = loss_hawkes, loss_args = loss_args, time1 = time1, time2 = time2)

  input_list$distance <- sapply(output, function(x){x$distance})

  input_list$x <- sapply(output, function(x){x$x})

  return(input_list)
}
