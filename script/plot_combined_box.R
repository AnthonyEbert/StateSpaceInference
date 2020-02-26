
library(ggplot2)
library(ggalt)
library(dplyr)

x = readRDS("runs/Volatility/18/combined_df.rds")

true_theta <- data.frame(value = c(0.8))

x1 <- x[-1,] %>%
  group_by(type, seed) %>%
  do(data.frame(matrix(Hmisc::wtd.quantile(.$value, .$weight * 1e4, normwt = FALSE, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), nrow = 1)))

ggplot(x) +
  aes(y = value, weight = weight, col = type, group = seed) +
  geom_boxplot() +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab(expression(Parameter~value:~theta)) +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  ylab("ABC posterior density") +
  xlim(c(0.5, 1))
