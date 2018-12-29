library(ggtern)
input <- get_parameter(full_list, spread = TRUE)
theta2_df <- input[which(input$time == 15),]
dens <- kde2d.weighted(theta2_df$`1`, theta2_df$`2`, w = theta2_df$weight)
dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
ggplot(theta2_df, aes(x = `1`, y = `2`)) +
  geom_point() +
  geom_contour(aes(x=x, y=y, z=z), data= dfdens)
