library(tidyverse)



load(here::here("data", "IND", "d_cases.rdata"))
load(here::here("data", "IND", "d_report.rdata"))
load(here::here("data", "IND", "d_burden.rdata"))
load(here::here("data", "IND", "d_prev.rdata"))
load(here::here("data", "IND", "d_pop.rdata"))




amp <- mean(d_case_all$Amp_all)
year0 <- d_prev$Year[1]


d_ts <- d_case_all %>% 
  left_join(d_report %>% select(Year, starts_with("Pr_Case_"))) %>% 
  filter(Year >= 2017) %>% 
  mutate(
    x = Pr_Case_Public + Pr_Case_Private,
    N_Case_Pu = N_Case * (Pr_Case_Public) / x,
    N_Case_Pr = N_Case * Pr_Case_Private /x 
  ) %>% 
  select(Year, N_Case_Pu, N_Case_Pr) %>% 

  left_join(d_burden) %>% 
  mutate(
    adr = -mean(diff(log(d_burden$Inc_M))),
    k = exp(-adr * (Year - year0)),
    inc0 = mean(Inc_M / k),
    inc = inc0 * k,
    Mor_sig = (Mor_U - Mor_L) / 2 / 1.96
  ) %>% 
  left_join(d_pop_all)
  


d_ts


dat <- with(d_ts, {
  list(
    Years = Year,
    n_t = length(Year),
    adr = adr,
    inc0 = inc0[1],
    
    Mor_mu = Mor_M,
    Mor_sig = Mor_sig,
    Pop = N_Pop,
    Case_Pu = round(N_Case_Pu),
    Case_Pr = round(N_Case_Pr),
    Amp = amp,
    YearSurveyed = year0
  )
})


prev <- local({
  prev <- d_prev %>% filter(Tag == "All")
  n_subject = prev$N_Subject[1]
  n_prev <- set_names(prev$N_Prev, prev$State)
  
  res <- as.list(n_prev[c("Asym", "Sym", "CS")])
  res$N = n_subject
  res
})


exo <- list(
  r_death_s = 0.12,
  r_death_bg = d_ts %>% summarise(r = weighted.mean(R_Die, N_Pop)) %>% pull(r),
  scale_dur = 1,
  p_tx_die = 0.05,
  cap_report = 0.8
)


dat <- c(dat, prev, exo)



d_ts %>% 
  select(Year, starts_with("N_Case"), inc0, adr, N_Pop, Amp) %>% 
  mutate(
    k = exp(-adr * (Year - year0)),
    IncR = inc0 * k,
    CNR_Pu = N_Case_Pu / N_Pop,
    CNR_Pr = N_Case_Pr / N_Pop,
    PrevA = prev$Asym / prev$N * k * Amp,
    PrevS = prev$Sym / prev$N * k * Amp,
    PrevC = prev$CS / prev$N * k * Amp,
    r_onset = IncR / PrevA - 0.1,
    Onset = r_onset * PrevA,
    
    Onset - (r_det_pu + r_det_pr + 0.3),
    
    
    r_aware = Onset / PrevS - 0.3,
    Aware = r_aware * PrevC,
    
    PPV = 0.85,
    TP_Pu = CNR_Pu * PPV,
    r_det_pu = TP_Pu / PrevC,
    
    r_det_pr = r_det_pu * 0.1,
    ar_det = mean(diff(log(r_det_pu))),

    
    r_det = r_det_pu + r_det_pr,
    
    ps = Onset / (r_det + 0.3 - adr),
    Det_Pu = ps * r_det_pu,
    Det_Pr = ps * r_det_pr,
    
    PPV_Pu = Det_Pu / CNR_Pu,
    
    TP_Pu = CNR_Pu * 1,
    TP_Pr = CNR_Pr * 1,
    TP = TP_Pu + TP_Pr,
    Det_Pr = PrevC * r_det - TP_Pu,
    Under_Pr = Det_Pr / TP_Pr
  ) %>% 
  select(-inc0, -adr, -N_Pop, -starts_with("N_Case"))# %>% 
  ggplot() + 
  geom_line(aes(x = Year, y = r_det_pu))

