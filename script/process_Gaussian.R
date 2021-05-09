# Process Gaussian skewed

library(dplyr)
library(forcats)
library(ggplot2)
library(ggalt)

load("theta_df.RData")
load("stan_df.RData")

theta_df_2 <- theta_df %>% filter(time == max(time)) %>% select(-time) %>% mutate(type = "ABC")

theta_df_2$parameter <- forcats::fct_recode(theta_df_2$parameter, sigma = "1", gamma = "2")

stan_df_2 <- stan_df %>%
  select(-seed) %>%
  tidyr::gather(key = parameter, value = value, -type, -weight)

combined_df <- bind_rows(theta_df_2, stan_df_2)

theta_of_interest <- combined_df

limits <- data.frame(parameter = rep(c("sigma","gamma"), each = 2), value = c(c(0.1, 0.5), c(0.2, 0.8)), weight = 0, type = "limits")

#theta_of_interest <- bind_rows(theta_of_interest, limits)

true_theta <- data.frame(parameter = factor(c("sigma", "gamma")), value = c(0.25, 0.5))

parameter_plot <- ggplot(theta_of_interest) +
  aes(x = value, weights = weight, color = factor(type)) +
  geom_density(bw = 0.01) +
  facet_wrap(~parameter, labeller = label_parsed, scales = "free_x") +
  geom_vline(data = true_theta, mapping = aes(xintercept = value), col = "red") +
  xlab("Parameter value") +
  scale_y_continuous(expand = expand_scale(mult = c(0,0.05))) +
  ggthemes::theme_few() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  ylab("ABC posterior density")

parameter_plot
