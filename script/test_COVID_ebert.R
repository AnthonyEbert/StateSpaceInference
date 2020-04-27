library(dplyr)
library(StateSpaceInference)
library(StateSpaceCOVID)


#italy_data = COVID19::covid19() %>% filter(country == "Italy") %>% process_data()

x = c(5e7,1000,0,0,0,0)
names(x) <- names_david

theta <- c(beta = 0.01, delta = 0.025, kappa = 0.7, eta = 0.45, gamma = 0.3)
log_alpha <- exp(rnorm(1))

pop <- sum(australia_data[1,])

inp <- list(transitions = transitions_david, lvrates = lvrates_ebert, data = australia_data, pop = pop)

output = adaptivetau::ssa.adaptivetau(x, transitions_david, inp$lvrates, params = as.list(c(theta, alpha = exp(log_alpha))), tf = 1)

## Set up SMCABC

Ntheta = 100
Nx = 400
pacc = 0.25
cl = parallel::makeCluster(parallel::detectCores())
TT = 48

prior.p1 = c(0, 0, 0, 0, 0);
prior.p2 = c(0.15, 0.15, 2, 1, 2);

prior = protoABC::prior_unif(prior.p1, prior.p2, var_names = names(theta))
dprior = protoABC::prior_unif(prior.p1, prior.p2, var_names = names(theta), eval = TRUE)

prior_sample <- prior(Ntheta)

trans <- function(x, trans_args){
  return(log(x))
}

invtrans <- function(x, trans_args){
  return(exp(x))
}

full_list <- StateSpaceInference::SMC2_ABC(
  prior_sample,
  dprior = dprior,
  loss = loss_ebert,
  loss_args = inp,
  Ntheta = Ntheta,
  Nx = Nx,
  pacc = pacc,
  cl = cl,
  ESS_threshold = 0.1,
  TT = TT,
  cov_coef = 0.5^2,
  trans = trans,
  invtrans = invtrans
)

library(ggplot2)

state_df <- get_state(full_list, probs = c(0.25, 0.5, 0.75))
theta_df <- get_parameter(full_list)

state_list <-
  lapply(state_df, function(i) {
    as.data.frame(i) %>% mutate(quant = rownames(.)) %>% tidyr::gather(, ,-quant)
  }) %>% bind_rows(.id = "column_label") %>%
  tidyr::spread(key = quant, value = value) %>%
  rename(lower = `25%`, med = `50%`, upper = `75%`)

ggplot(state_list) +
  aes(x = as.numeric(column_label), ymin = lower, ymax = upper, y = med) +
  facet_wrap(~key, scales = "free") +
  geom_crossbar()


ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper, y = ) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")
