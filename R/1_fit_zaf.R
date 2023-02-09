library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


load("data/ZAF/d_burden.rdata")
load("data/ZAF/d_cases.rdata")
load("data/ZAF/d_prev.rdata")
load("data/ZAF/d_pop.rdata")


d_pop <- d_pop_all %>% 
  filter(Year %in% d_burden$Year) %>% 
  mutate(dea = N_Pop * R_Die)


year0 <- d_prev$Year[1]
adr <- -mean(diff(log(d_burden$Inc_M)))

d_case_all <- d_case_all %>% 
  filter(Year >= 2014)

d_burden <- d_burden %>% 
  mutate(
    k = exp(-adr * (Year - year0)),
    inc0 = mean(Inc_M / k),
    inc = inc0 * k
  )


dat <- list(
  Years = d_burden$Year,
  n_t = nrow(d_burden),
  adr = adr,
  inc0 = d_burden$inc0[1],
  Mor_mu = d_burden$Mor_M,
  Mor_sig = (d_burden$Mor_U - d_burden$Mor_L) / 2 / 1.96,
  Pop = d_pop$N_Pop,
  Case = d_case_all$N_Case,
  Amp = 0.9090663,
  YearSurveyed = year0
)



prev <- local({
  prev <- d_prev %>% filter(Tag == "All")
  n_subject = prev$N_Subject[1]
  n_prev <- set_names(prev$N_Prev, prev$State)
  
  res <- as.list(n_prev[c("Asym", "Sym", "CS")])
  res$N = n_subject
  res
})


exo <- list(
  r_death_s = 0.2,
  r_death_bg = d_pop %>% summarise(r = weighted.mean(R_Die, N_Pop)) %>% pull(r),
  scale_dur = 1,
  p_tx_die = 0.09,
  cap_report = 0.8,
  ppv = 0.7
)


dat <- c(dat, prev, exo)


m_cas <- rstan::stan_model("stan/cascade_report.stan")

post <- rstan::sampling(m_cas, data = dat, iter = 5000, warmup = 4500)

summary(post)$summary

summary(post, pars='nr')$summary
ts.plot(data.frame(mu = dat$Case / dat$Pop, e = summary(post, pars='nr')$summary[, 'mean']))

ts.plot(data.frame(e = 1 - summary(post, pars='p_under')$summary[, 'mean']))

tab <- data.frame(rstan::extract(post, pars = c("prv0", "r_death_a", "r_death_tx", 
                                                "r_sym", "r_aware", "r_det", "r_sc", "rt_report"))) %>% 
  as_tibble() %>% 
  bind_cols(exo) %>% 
  mutate(
    Year0 = dat$YearSurveyed,
    inc0 = dat$inc0,
    adr = dat$adr
  )



write_csv(tab, file = here::here("results", "pars_ZAF.csv"))
