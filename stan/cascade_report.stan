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
  
  real<lower=0> inc0;
  real adr;
  
  real<lower=0> Mor_mu[n_t];
  real<lower=0> Mor_sig[n_t];
  int<lower=0> Case[n_t]; // notification counts

  real<lower=0> Amp;

  // Prior knowledge
  real<lower=0> scale_dur;

  // Exogenous variables
  real<lower=0> r_death_s;
  real<lower=0> r_death_bg;
  real<lower=0, upper=1> p_tx_die;
  
  real<lower=0.5, upper=1> cap_report;
}
parameters {
  real<lower=0> r_sym;
  real<lower=0> r_aware;
  real<lower=0> r_det;
  real<lower=0.1, upper = 0.3> r_sc;
  real<lower=0, upper=1> rr_die_a;
  
  real<lower=0, upper=0.5> report0;
  real<lower=0> rt_report;
  
  real<lower=0.5, upper=0.85> ppv;
}
transformed parameters {
  vector<lower=0>[n_t] inc;
  
  real<lower=0> mor0;
  vector<lower=0>[n_t] mor;
  
  real<lower=0> r_death_a;
  real<lower=0> r_death_tx;
  
  real<lower=0> ra;
  real<lower=0> rs;
  real<lower=0> rc;

  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;
  real<lower=0, upper=1> prv0;
  
  vector<lower=0>[n_t] prv;  
  vector<lower=0>[n_t] prv_t_a;
  vector<lower=0>[n_t] prv_t_s;
  vector<lower=0>[n_t] prv_t_c;

  vector<lower=0>[n_t] tp;
  vector<lower=0>[n_t] nr;
  vector<lower=0, upper=1>[n_t] p_under;
  
  r_death_tx = 2 * p_tx_die / (1 - p_tx_die);
  r_death_a = r_death_s * rr_die_a;
  
  ra = r_sc + r_death_a + r_death_bg;
  rs = r_sc + r_death_s + r_death_bg;
  rc = r_sc + r_death_s + r_death_bg;
  
  prv_a = inc0 / (r_sym + ra - adr);
  prv_s = prv_a * r_sym / (r_aware + rs - adr);
  prv_c = prv_s * r_aware / (r_det + rc - adr);
  prv0 = prv_a + prv_s + prv_c;
  
  mor0 = r_death_a * prv_a + r_death_s * (prv_s + prv_c) + r_death_tx * prv_c * r_det / (2 + r_death_tx - adr);
  
  
  for (i in 1:n_t) {
    p_under[i] = 1 - (report0 + (cap_report - report0) / (1 + exp(- rt_report * (Years[i] - 2020))));
    
    inc[i] = inc0 * exp(- adr * (Years[i] - YearSurveyed));
    mor[i] = mor0 * exp(- adr * (Years[i] - YearSurveyed));
    prv[i] = prv0 * exp(- adr * (Years[i] - YearSurveyed));
    prv_t_a[i] = prv_a * exp(- adr * (Years[i] - YearSurveyed));
    prv_t_s[i] = prv_s * exp(- adr * (Years[i] - YearSurveyed));
    prv_t_c[i] = prv_c * exp(- adr * (Years[i] - YearSurveyed));
    
    tp[i] = prv_t_c[i] * r_det;
    nr[i] = tp[i] / ppv * (1 - p_under[i]);
  }
}
model {
  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_aware ~ inv_gamma(scale_dur, scale_dur);
  r_det ~ inv_gamma(scale_dur, scale_dur);
  r_sc ~ uniform(0.1, 0.3);

  target += binomial_lpmf(Asym | N, prv_a / Amp);
  target += binomial_lpmf(Sym | N, prv_s / Amp);
  target += binomial_lpmf(CS | N, prv_c / Amp);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(Case[i] | nr[i] * Pop[i]);
    target += normal_lpdf(Mor_mu[i] | mor[i], Mor_sig[i]);
  }
}
