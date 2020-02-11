data {
  int<lower=0> TT;   // # time points (equally spaced)
  vector[TT] y;      // mean corrected return at time t
  real mu;
  real sh;
}
parameters {
  real<lower=-1> theta;         // white noise shock scale
  vector[TT] x;                 // log volatility at time t
}
model {
  theta ~ uniform(-1, 1);
  x[1] ~ normal(mu/(1-theta), sh^2/(1-theta^2));
  for (t in 2:TT)
    x[t] ~ normal(mu + theta * x[t-1], sh);
  y ~ normal(0, exp(0.5 * x));
}
