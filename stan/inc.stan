data {
  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_t]; // population size
  real<lower=0> Years[n_t]; // years of the notification data

  real YearSurveyed;

  real<lower=0> Inc_mu[n_t];
  real<lower=0> Inc_sig[n_t];

  // Prior knowledge
  // real<lower=0> scale_dur;
}
parameters {
  real<lower=0, upper=1> inc0;
  real<lower=0, upper=0.1> adr;
}
transformed parameters {
  vector<lower=0>[n_t] inc;

  for (i in 1:n_t) {
    inc[i] = inc0 * exp(- adr * (Years[i] - YearSurveyed));
  }
}
model {
  inc0 ~ uniform(0, 1);
  
  for (i in 1:n_t) {
    target += normal_lpdf(Inc_mu[i] | inc[i], Inc_sig[i]);
  }
}
