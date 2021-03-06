data {
  int<lower=0> TT;      // # time points (equally spaced)
  matrix[TT,100] y;      // mean corrected return at time t
}
parameters {
  real<lower=.1, upper=.5> sigma;
  real<lower=0.2, upper=4> gamma;
  vector[TT] x;                 // latent states
}

model {
  // x_std ~ std_normal();
  sigma ~ uniform(0.1, 0.5);
  gamma ~ uniform(0.2, 4);
  x[1] ~ normal(0, 1);
  for (t in 2:TT)
    x[t] ~ normal(x[t-1], 1);
  for (t in 1:TT)
    for (j in 1:100)
      y[t,j] ~ skew_normal(x[t], sigma, gamma);
}
