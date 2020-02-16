library(ggplot2)
library(ggalt)
library(dplyr)

load("state_df.RData")
load("theta_df.RData")
load("theta_stan.RData")

fullpath = getwd()
directoryname = basename(fullpath)

nr <- dim(state_df)[1]

state_df <- state_df %>%
  bind_rows(state_df[nr,])

state_df$time[nr] <- state_df$time[nr] + 1

state_plot <- ggplot(state_df) +
  aes(x = time, y = state, ymin = lower, ymax = upper) +
  geom_point(size = 1, col = "red") +
  geom_errorbar(mapping = aes(x = time, y = med), shape = 95, size = 0.5) +
  ggthemes::theme_base() +
  xlab("Time: t") +
  ylab(expression(State:~x[t]))

ggsave(paste0("Vol_", directoryname, "_state.pdf"), height = 10, width = 15, units = "cm", plot = state_plot)

theta_of_interest <- theta_df[which(theta_df$time == max(theta_df$time)),] %>% select(-time, -parameter)

limits <- data.frame(value = c(0.5,1), weight = 0)

true_theta <- data.frame(value = c(0.8))

theta_of_interest <- bind_rows(theta_of_interest, limits) %>%
  mutate(type = "Approximate posterior")

theta_stan <- data.frame(value = theta_stan, weight = 1/length(theta_stan), type = "True posterior")

theta_of_interest <- bind_rows(theta_of_interest, theta_stan)

parameter_plot <- ggplot(theta_of_interest) +
  aes(x = value, weights = weight, col = type) +
  geom_density(adjust = 1) +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab(expression(Parameter~value:~theta)) +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "bottom", legend.title=element_blank()) +
  ylab("ABC posterior density")

ggsave(paste0("Vol_", directoryname, "_parameter.pdf"), height = 8.4, width = 12, units = "cm", plot = parameter_plot)
