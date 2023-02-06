data {
  // Data from a prevalence survey
  int<lower=0> N;
  int<lower=0> Asym;
  int<lower=0> Sym;
  int<lower=0> CS;
  real YearSurveyed; // timing of the survey

  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  real<lower=0> Years[n_t]; // years of the notification data  
  
  int<lower=0> Pop[n_t]; // population size
  
  real<lower=0> Inc_mu[n_t];
  real<lower=0> Inc_sig[n_t];
  real<lower=0> Mor_mu[n_t];
  real<lower=0> Mor_sig[n_t];
  int<lower=0> Case[n_t]; // notification counts

  real<lower=0> Amp_age[n_t];
  real<lower=0> Amp_ep[n_t];

  // Prior knowledge
  real<lower=0> scale_dur;

  // Exogenous variables
  real<lower=0> r_death_s;
  real<lower=0> r_death_bg;
  real<lower=0, upper=1> p_tx_die;
}
parameters {
  real<lower=0, upper=1> prv0;
  real<lower=-0.2, upper=0.2> adr;
  real<lower=0> r_sym;
  real<lower=0> r_aware;
  real<lower=0> r_det;
  real<lower=0.1, upper = 0.3> r_sc;
  real<lower=0, upper=0.4> p_under;
  real<lower=0, upper=1> rr_die_a;
  real<lower=0.5, upper=1> ppv;
}
transformed parameters {
  real<lower=0> inc0;
  vector<lower=0>[n_t] inc;
  
  real<lower=0> mor0;
  vector<lower=0>[n_t] mor;
  
  real<lower=0> r_death_a;
  real<lower=0> r_death_tx;
  
  real<lower=0> ra;
  real<lower=0> rs;
  real<lower=0> rc;

  real<lower=0> a0;
  real<lower=0> c0;
  
  real<lower=0, upper=1> pr_a;
  real<lower=0, upper=1> pr_s;
  real<lower=0, upper=1> pr_c;
  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;

  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] tp;
  vector<lower=0>[n_t] nr;
  
  r_death_tx = 2 * p_tx_die / (1 - p_tx_die);
  r_death_a = r_death_s * rr_die_a;
  
  ra = r_sc + r_death_a + r_death_bg;
  rs = r_sc + r_death_s + r_death_bg;
  rc = r_sc + r_death_s + r_death_bg;
  
  a0 = (rs + r_aware - adr) / r_sym;
  c0 = r_aware / (rc + r_det - adr);
    
  pr_a = a0 / (a0 + 1 + c0);
  pr_s = 1 / (a0 + 1 + c0);
  pr_c = c0 / (a0 + 1 + c0);
  
  prv_a = prv0 * pr_a;
  prv_s = prv0 * pr_s;
  prv_c = prv0 * pr_c;
  
  inc0 = (ra + r_sym - adr) * prv_a;
  
  mor0 = r_death_a * prv_a + r_death_s * (prv_s + prv_c) + r_death_tx * prv_c * r_det / (2 + r_death_tx - adr);
  
  
  for (i in 1:n_t) {
    inc[i] = inc0 * exp(- adr * (Years[i] - YearSurveyed));
    mor[i] = mor0 * exp(- adr * (Years[i] - YearSurveyed));
    prv[i] = prv0 * exp(- adr * (Years[i] - YearSurveyed)) * (Amp_ep[i] / Amp_age[i]);
    
    tp[i] = prv[i] * pr_c * r_det / ppv;
    nr[i] = tp[i] / ppv * (1 - p_under);
  }
}
model {
  prv0 ~ uniform(0, 1);

  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_aware ~ inv_gamma(scale_dur, scale_dur);
  r_det ~ inv_gamma(scale_dur, scale_dur);
  r_sc ~ uniform(0.1, 0.3);

  p_under ~ beta(1, 5);

  adr ~ uniform(-0.2, 0.2);

  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sym | N, prv_s);
  target += binomial_lpmf(CS | N, prv_c);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(Case[i] | nr[i] * Pop[i]);
    target += normal_lpdf(Inc_mu[i] | inc[i], Inc_sig[i]);
    target += normal_lpdf(Mor_mu[i] | mor[i], Mor_sig[i]);
  }
}
