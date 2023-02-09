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
}
parameters {
  real<lower=0, upper=1> prv0;
  real<lower=0> r_sym;
  real<lower=0> r_aware0;
  real<lower=0> r_det0;
  real<lower=1, upper=1.5> rr_det_t;
  
  real<lower=0.1, upper = 0.3> r_sc;
  real<lower=0, upper=0.3> p_under;
  real<lower=0, upper=1> rr_die_a;
  real<lower=0.6, upper=0.85> ppv;
}
transformed parameters {
  vector<lower=0>[n_t] inc;
  vector<lower=0>[n_t] mor;
  
  real<lower=0> r_death_a;
  real<lower=0> r_death_tx;
  
  vector<lower=0>[n_t] r_aware;
  vector<lower=0>[n_t] r_det;
  
  real<lower=0> ra;
  real<lower=0> rs;
  real<lower=0> rc;

  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;
  
  vector<lower=0>[n_t] prv;  
  vector<lower=0>[n_t] prv_t_a;
  vector<lower=0>[n_t] prv_t_s;
  vector<lower=0>[n_t] prv_t_c;
  
  
  vector<lower=0>[n_t] tp;
  vector<lower=0>[n_t] nr;
  
  r_death_tx = 2 * p_tx_die / (1 - p_tx_die);
  r_death_a = r_death_s * rr_die_a;
  
  ra = r_sc + r_death_a + r_death_bg;
  rs = r_sc + r_death_s + r_death_bg;
  rc = r_sc + r_death_s + r_death_bg;
  
  prv_a = inc0 / (r_sym + ra - adr);
  prv_s = prv_a * r_sym / (r_aware0 + rs - adr);
  prv_c = prv_s * r_aware0 / (r_det0 + rc - adr);
  
  
  for (i in 1:n_t) {
    inc[i] = inc0 * exp(- adr * (Years[i] - YearSurveyed));
    r_aware[i] = r_aware0 * rr_det_t ^ (Years[i] - YearSurveyed);
    r_det[i] = r_det0 * rr_det_t ^ (Years[i] - YearSurveyed);
    
    prv_t_a[i] = inc[i] / (r_sym + ra - adr);
    prv_t_s[i] = r_sym * prv_t_a[i] / (r_sym + rs - adr);
    prv_t_c[i] = r_aware[i] * prv_t_s[i] / (r_det[i] + rc - adr);
    
    mor[i] = r_death_a * prv_t_a[i];
    mor[i] += r_death_s * (prv_t_s[i] + prv_t_c[i]);
    mor[i] += r_death_tx * prv_t_c[i] * r_det[i] / (2 + r_death_tx - adr);;
    
    prv[i] = prv_t_a[i] + prv_t_s[i] + prv_t_c[i];
    
    tp[i] = prv_t_c[i] * r_det[i] / ppv;
    nr[i] = tp[i] / ppv * (1 - p_under);
  }
}
model {
  prv0 ~ uniform(0, 1);

  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_aware0 ~ inv_gamma(scale_dur, scale_dur);
  r_det0 ~ inv_gamma(scale_dur, scale_dur);
  r_sc ~ uniform(0.1, 0.3);

  p_under ~ beta(1, 5);

  target += binomial_lpmf(Asym | N, prv_a / Amp);
  target += binomial_lpmf(Sym | N, prv_s / Amp);
  target += binomial_lpmf(CS | N, prv_c / Amp);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(Case[i] | nr[i] * Pop[i]);
    target += normal_lpdf(Mor_mu[i] | mor[i], Mor_sig[i]);
  }
}
