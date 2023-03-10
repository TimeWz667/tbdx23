library(tidyverse)



for (iso in c("ZAF", "IND")) {
  load(here::here("data", iso, "d_pop.rdata"))
  load(here::here("data", iso, "d_tx.rdata"))
  load(here::here("data", iso, "d_cases.rdata"))
  load(here::here("data", iso, "d_burden.rdata"))
  
  
  targets <- d_case_all %>% 
    left_join(d_pop_all) %>% 
    left_join(d_burden) %>% 
    mutate(
      CNR_mu = N_Case / N_Pop,
      CNR_eps = sqrt(CNR_mu * (1 - CNR_mu) / N_Pop),
      IncR_mu = Inc_M,
      IncR_eps = (Inc_U - Inc_L) / 2,
      MorR_mu = Mor_M,
      MorR_eps = (Mor_U - Mor_L) / 2
    ) %>% 
    select(Year, N = N_Pop, ends_with("_mu"), ends_with("_eps")) %>% 
    pivot_longer(-c(Year, N)) %>% 
    separate(name, c("Index", "name")) %>% 
    pivot_wider() %>% 
    rename(m = mu) %>% 
    mutate(Tag = "All", error = eps / m)

  
  targets
  
  write_csv(targets, file = here::here("pars", iso, "targets.csv"))
  
}



# India Targets with public/private ----
iso <- "IND"

targets <- read_csv(here::here("pars", iso, "targets.csv"))

load(here::here("data", iso, "d_report.rdata"))



targets_pupr <- targets %>% 
  filter(Index == "CNR") %>% 
  left_join(
    d_report %>% 
      mutate(
        Pr_Pub = Pr_Case_Public + Pr_Case_ACF,
        Pr_Pri = Pr_Case_Private
      ) %>% select(Year, Pr_Pub, Pr_Pri)
  ) %>% 
  mutate(
    m_Pub = m * Pr_Pub,
    m_Pri = m * Pr_Pri
  ) %>% 
  select(Year, N, Index, starts_with("m_")) %>% 
  pivot_longer(c(m_Pub, m_Pri), values_to = "m") %>% 
  extract(name, "Tag", "m_(\\S+)") %>% 
  mutate(error = 0.1, eps = m * error)

  
write_csv(targets_pupr, file = here::here("pars", iso, "targets_pupr.csv"))


# South Africa targets with HIV ----
iso <- "ZAF"

load(here::here("data", iso, "d_pop.rdata"))
load(here::here("data", iso, "d_cases.rdata"))
load(here::here("data", iso, "d_burden.rdata"))
load(here::here("data", iso, "d_hiv.rdata"))


d_case <- d_case_hiv %>% 
  group_by(Tag) %>% 
  mutate(
    Tag = ifelse(Tag == "PLHIV", "HIV", "NonHIV"),
    Amp = mean(Amp_all)
  ) %>% 
  select(Year, Tag, Case) %>% 
  pivot_wider(values_from = Case, names_from = Tag, names_prefix = "Case_")


targets_hiv <- d_hiv %>% 
  mutate(
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
    IncR_HIV = Inc_HIV / Pop_HIV,
    IncR_NonHIV = Inc_NonHIV / Pop_NonHIV,
    
    IncR_L_HIV = IncR_HIV * Inc_L / Inc_M,
    IncR_U_HIV = IncR_HIV * Inc_U / Inc_M,
    IncR_L_NonHIV = IncR_NonHIV * Inc_L / Inc_M,
    IncR_U_NonHIV = IncR_NonHIV * Inc_U / Inc_M,
    
    IncR_eps_HIV = (IncR_U_HIV - IncR_L_HIV) / 2 / 1.96,
    IncR_eps_NonHIV = (IncR_U_NonHIV - IncR_L_NonHIV) / 2 / 1.96,
    
    Mor_mu_HIV = Mor_M,
    Mor_mu_NonHIV = ((Pop_HIV + Pop_NonHIV) * Mor_M - Mor_mu_HIV) / Pop_NonHIV,
    Mor_L_HIV = Mor_L * Mor_mu_HIV / Mor_M,
    Mor_U_HIV = Mor_U * Mor_mu_HIV / Mor_M,
    Mor_L_NonHIV = Mor_L * Mor_mu_HIV / Mor_M,
    Mor_U_NonHIV = Mor_U * Mor_mu_HIV / Mor_M,
    MorR_mu_HIV = Mor_mu_HIV,
    MorR_mu_NonHIV = Mor_mu_NonHIV,
    MorR_eps_HIV = (Mor_U_HIV - Mor_L_HIV) / 2,
    MorR_eps_NonHIV = (Mor_U_NonHIV - Mor_L_NonHIV) / 2,
  ) %>% 
  select(Year, starts_with(c("Pop_", "IncR_", "MorR_mu_", "MorR_eps_", "Case", "Amp"))) %>% 
  pivot_longer(-Year) %>% 
  extract(name, c("Index", "Tag"), "(\\S+)_(HIV|NonHIV)") %>% 
  filter(!is.na(Index)) %>% 
  pivot_wider(names_from = Index) %>% 
  mutate(
    CNR_mu = Case / Pop,
    CNR_eps = sqrt(CNR_mu * (1 - CNR_mu) / Pop) * 1.96
  ) %>% 
  rename(IncR_mu = IncR) %>% 
  select(Year, Tag, N = Pop, ends_with("_mu"), ends_with("_eps")) %>% 
  pivot_longer(-c(Year, Tag, N)) %>% 
  separate(name, c("Index", "name")) %>% 
  pivot_wider() %>% 
  rename(m = mu) %>% 
  mutate(error = eps / m)
  arrange(Tag, Year)


  write_csv(targets_hiv, file = here::here("pars", iso, "targets_hiv.csv"))

