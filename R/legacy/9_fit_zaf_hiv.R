library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


load("data/ZAF/d_burden.rdata")
load("data/ZAF/d_cases.rdata")
load("data/ZAF/d_prev.rdata")
load("data/ZAF/d_pop.rdata")
load("data/ZAF/d_hiv.rdata")


d_case <- d_case_hiv %>% 
  group_by(Tag) %>% 
  mutate(
    Tag = ifelse(Tag == "PLHIV", "HIV", "NonHIV"),
    Amp = mean(Amp_all)
  ) %>% 
  select(Year, Tag, Case, Amp) %>% 
  pivot_wider(values_from = c(Case, Amp), names_from = Tag)


d_prev <- d_prev %>% filter(Tag != "All") %>% mutate(Tag = ifelse(Tag == "PLHIV", "HIV", Tag))
year0 <- d_prev$Year[1]


d_gp <- d_hiv %>% 
  mutate(
    Year, 
    Pop_HIV = PrvHIV, 
    Pop_NonHIV = N_Pop - Pop_HIV,
    morbg = (N_Pop * R_Die - MorHIV) / N_Pop,
    Mor_HIV = MorHIV / Pop_HIV + morbg,
    Mor_NonHIV = morbg,
    MorTB_HIV = MorTB_HIV / Pop_HIV
  ) %>% 
  select(
    Year, starts_with(c("Pop_", "Mor_", "MorTB_")), IncTB_HIV
  ) %>% 
  inner_join(d_burden %>% select(-Country), "Year") %>% 
  inner_join(d_case) %>% 
  mutate(
    Inc_HIV = IncTB_HIV,
    Inc_NonHIV = (Pop_HIV + Pop_NonHIV) * Inc_M - IncTB_HIV,
    Inc_HIV = Inc_HIV / Pop_HIV,
    Inc_NonHIV = Inc_NonHIV / Pop_NonHIV,
    Mor_mu_HIV = Mor_M,
    Mor_mu_NonHIV = ((Pop_HIV + Pop_NonHIV) * Mor_M - Mor_mu_HIV) / Pop_NonHIV,
    Mor_L_HIV = Mor_L * Mor_mu_HIV / Mor_M,
    Mor_U_HIV = Mor_U * Mor_mu_HIV / Mor_M,
    Mor_L_NonHIV = Mor_L * Mor_mu_HIV / Mor_M,
    Mor_U_NonHIV = Mor_U * Mor_mu_HIV / Mor_M,
    Mor_sig_HIV = (Mor_U_HIV - Mor_L_HIV) / 2 / 1.96,
    Mor_sig_NonHIV = (Mor_U_NonHIV - Mor_L_NonHIV) / 2 / 1.96
  ) %>% 
  select(Year, starts_with(c("Pop_", "Mor_", "Inc_", "Mor_mu_", "Mor_sig", "Case", "Amp"))) %>% 
  pivot_longer(-Year) %>% 
  extract(name, c("Index", "Group"), "(\\S+)_(HIV|NonHIV)") %>% 
  filter(!is.na(Index)) %>% 
  pivot_wider(names_from = Index) %>% 
  group_by(Group) %>% 
  mutate(
    adr = - mean(diff(log(Inc))),
    k = exp(-adr * (Year - year0)),
    inc0 = mean(Inc / k),
    inc = inc0 * k
  ) %>% 
  select(- k) %>% 
  rename(Years = Year) %>% 
  ungroup()



# Model fitting
m_cas <- rstan::stan_model("stan/cascade.stan")


posts <- lapply(c(HIV = "HIV", NonHIV = "NonHIV"), function(gp) {
  ds <- d_gp %>% 
    filter(Group == gp)
  
  d0 <- ds %>% 
    select(Years, Mor_mu, Mor_sig, Pop, Case) %>% 
    as.list()
  
  d0$n_t <- nrow(ds)
  d0$adr <- ds$adr[1]
  d0$inc0 <- ds$inc0[1]
  d0$Amp <- ds$Amp[1]
  d0$YearSurveyed = year0
  
  prev <- d_prev %>% filter(Tag == gp)
  n_subject = prev$N_Subject[1]
  n_prev <- set_names(prev$N_Prev, prev$State)
  
  d1 <- as.list(n_prev[c("Asym", "Sym", "CS")])
  d1$N <- n_subject
  
  exo <- list(
    r_death_s = 0.12,
    r_death_bg = ds %>% summarise(r = weighted.mean(Mor, Pop)) %>% pull(r),
    scale_dur = 1,
    p_tx_die = 0.09
  )
  
  dat <- c(d0, d1, exo)
  
  post <- rstan::sampling(m_cas, data = dat, iter = 5000, warmup = 4500)

  list(
    post = post,
    exo = exo,
    dat = dat
  )
}) 
  





print(summary(posts[[1]]$post)$summary)

print(summary(posts[[2]]$post)$summary)


summary(posts[[1]]$post, pars='nr')$summary
summary(posts[[2]]$post, pars='nr')$summary



tabs <- lapply(posts, function(res) {
  post <- res$post
  exo <- res$exo
  dat <- res$dat
  
  data.frame(rstan::extract(post, pars = c("prv0", "r_death_a", "r_death_tx", 
                                           "r_sym", "r_aware", "r_det", "r_sc", "p_under", "ppv"))) %>% 
    as_tibble() %>% 
    bind_cols(exo) %>% 
    mutate(
      Year0 = dat$YearSurveyed,
      inc0 = dat$inc0,
      adr = dat$adr,
      cap_report = 1 - p_under,
      rt_report = 0
    )
}) 
  

write_csv(tabs$HIV, file = here::here("results", "pars_ZAF_HIV.csv"))
write_csv(tabs$NonHIV, file = here::here("results", "pars_ZAF_NonHIV.csv"))

