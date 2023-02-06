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

d_case <- d_case_all %>% 
  filter(Year %in% d_burden$Year)

dat <- list(
  Years = d_burden$Year,
  n_t = nrow(d_burden),
  Inc_mu = d_burden$Inc_M,
  Inc_sig = (d_burden$Inc_U - d_burden$Inc_L) / 2 / 1.96,
  Pop = d_pop$N_Pop,
  Case = d_case$N_Case,
  Amp_age = d_case$Amp_age,
  Amp_ep = d_case$Amp_ep,
  YearSurveyed = d_prev$Year[1]
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
  r_death_a = 0,
  r_death_s = 0.12,
  r_death_bg = d_pop %>% summarise(r = weighted.mean(R_Die, N_Pop)) %>% pull(r),
  scale_dur = 1,
  ppv = 0.7
)


dat <- c(dat, prev, exo)



m_cas <- rstan::stan_model("stan/cascade.stan")

post <- rstan::sampling(m_cas, data = dat, iter = 5000, warmup = 4500)

summary(post)$summary


tab <- data.frame(rstan::extract(post, pars = c("prv0", "adr", "r_sym", "r_aware", "r_det", "r_sc", "p_under"))) %>% 
  as_tibble() %>% 
  bind_cols(exo) %>% 
  mutate(
    Year0 = dat$YearSurveyed,
    prv = prv0 * exp(- adr * (2023 - Year0))
  )


tab

write_csv(tab, file = here::here("results", "pars_ZAF.csv"))



