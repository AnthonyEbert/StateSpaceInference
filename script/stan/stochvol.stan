data {
  int<lower=0> TT;   // # time points (equally spaced)
  vector[TT] y;      // mean corrected return at time t
  real mu;
  real sh;
}
parameters {
  real<lower=-1> theta;         // white noise shock scale
  vector[TT] x_std;                 // log volatility at time t
}
transformed parameters {
  vector[TT] x = x_std * sh;  // now h ~ normal(0, sigma)
  x[1] /= sqrt(1 - theta * theta);  // rescale h[1]
  x += mu;
  for (t in 2:TT)
    x[t] += theta * (x[t-1]);
}
model {
  theta ~ uniform(-1, 1);
  x_std ~ std_normal();
  // x[1] ~ normal(mu/(1-theta), sh^2/(1-theta^2));
  // for (t in 2:TT)
  //   x[t] ~ normal(mu + theta * x[t-1], sh);
  y ~ normal(0, exp(0.5 * x));
}
