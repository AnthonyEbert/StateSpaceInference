library(ggplot2)
library(ggalt)
library(dplyr)

load("state_df.RData")
load("theta_df.RData")

fullpath = getwd()
directoryname = basename(fullpath)

nr <- dim(state_df)[1]

state_df <- state_df %>%
  bind_rows(state_df[nr,])

state_df$time[nr] <- state_df$time[nr] + 1

state_plot <- ggplot(state_df) +
  aes(x = time, y = state, ymin = lower, ymax = upper) +
  geom_point(shape = 16, size = 1, col = "red") +
  geom_pointrange(mapping = aes(x = time, y = med), shape = 95, size = 0.5) +
  ggthemes::theme_base() +
  xlab("Epoch: j") +
  ylab(expression(State:~x[j]))

ggsave(paste0("Vol_", directoryname, "_state.pdf"), height = 10, width = 15, units = "cm", plot = state_plot)

theta_of_interest <- theta_df[which(theta_df$time == max(theta_df$time)),] %>% select(-time, -parameter)

limits <- data.frame(value = c(-1,1), weight = 0)

theta_of_interest <- bind_rows(theta_of_interest, limits)

true_theta <- data.frame(value = c(0.95))

parameter_plot <- ggplot(theta_of_interest) +
  aes(x = value, weights = weight) +
  geom_density(adjust = 1) +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab(expression(Parameter~value:~phi)) +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  ylab("ABC posterior density")

ggsave(paste0("Vol_", directoryname, "_parameter.pdf"), height = 7, width = 10, units = "cm", plot = parameter_plot)
