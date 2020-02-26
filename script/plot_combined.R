
library(ggplot2)
library(ggalt)
library(dplyr)

x = readRDS("combined_df.rds")

true_theta <- data.frame(value = c(0.8))

ggplot(x) +
  aes(x = value, weight = weight, col = type) +
  geom_density(adjust = 1) +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab(expression(Parameter~value:~theta)) +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  ylab("ABC posterior density") +
  xlim(c(0.5, 1))
