data {
  int<lower=0> TT;      // # time points (equally spaced)
  matrix[TT,10] y;      // mean corrected return at time t
}
parameters {
  real<lower=.1, upper=.5> sigma;
  real<lower=.2, upper=.8> gamma;
  vector[TT] x;                 // latent states
}
// transformed parameters {
//   vector[TT] x = x_std * sh;  // now h ~ normal(0, sigma)
//   x[1] /= sqrt(1 - theta * theta);  // rescale h[1]
//   x += mu;
//   for (t in 2:TT)
//     x[t] += theta * (x[t-1]);
// }
model {
  // x_std ~ std_normal();
  x[1] ~ normal(0, 1);
  for (t in 2:TT)
    x[t] ~ normal(x[t-1], 1);
  for (t in 1:TT)
    for (j in 1:10)
      y[t,j] ~ skew_normal(x[t], sigma, gamma);
}
