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
  aes(x = time * 10 - 10, y = med / 3.5, ymin = lower/3.5, ymax = upper/3.5) +
  geom_step() +
  geom_ribbon(alpha = 0.2, stat = "stepribbon") +
  geom_step(mapping = aes(x = time * 10 - 10, y = state/3.5), col = "red") +
  ggthemes::theme_base() +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(Continuous~time:~tau)) +
  ylab(expression(Transformed~state:~logit^{-1}~(x[t])))

ggsave(paste0("Hawkes_", directoryname, "_state.pdf"), height = 10, width = 15, units = "cm", plot = state_plot)

theta_of_interest <- theta_df[which(theta_df$time == max(theta_df$time)),] %>% select(-time)

limits <- data.frame(parameter = factor(rep(c(1,2), each = 2)), value = rep(c(0.3, 0.7), 2), weight = 0)

theta_of_interest <- bind_rows(theta_of_interest, limits)

levels(theta_of_interest$parameter) <- c("theta[1]", "theta[2]")

true_theta <- data.frame(parameter = factor(c("theta[1]", "theta[2]")), value = c(0.5, 0.5))

parameter_plot <- ggplot(theta_of_interest) +
  aes(x = value, weights = weight) +
  geom_density(bw = 0.025) +
  facet_wrap(~parameter, labeller = label_parsed) +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab("Parameter value") +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  ylab("ABC posterior density")

ggsave(paste0("Hawkes_", directoryname, "_parameter.pdf"), height = 7, width = 10, units = "cm", plot = parameter_plot)
