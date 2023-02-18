library(tidyverse)
library(jsonlite)

source(here::here("r", "fn_fit.R"))
source(here::here("r", "fn_cas.R"))


year0 <- 2010


for (country in c("ZAF", "IND")) {
  # Demography
  iso = glue::as_glue(country)
  
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
      IncR_eps = (Inc_U - Inc_L) / 2 / 1.96,
      MorR_mu = Mor_M,
      MorR_eps = (Mor_U - Mor_L) / 2 / 1.96
    ) %>% 
    select(Year, N_Pop, ends_with("_mu"), ends_with("_eps"))
  
  write_csv(targets, file = here::here("data", "pars", iso, "targets.csv"))
  
  
  # Demography
  pars_pop <- get_pars_pop(d_pop_all, year0)
  
  jsonlite::write_json(pars_pop, here::here("data", "pars", iso, "pars_pop.json"), digits = 10, auto_unbox = T)
  
  
  # Treatment outcome
  ptx <- d_tx %>% 
    group_by(Country, Tag, Stat) %>% 
    summarise(across(starts_with("value_"), weighted.mean, w = N)) %>% 
    extract(Stat, "Type", "TxOut(\\S+)")
  
  tags <- unique(ptx$Tag)
  tags <- set_names(tags, tags)
  
  ptx <- lapply(tags, function(tag) {
    p <- ptx %>% filter(Tag == tag)
    
    as.list(set_names(p$value_M, p$Type))
  })
  
  
  jsonlite::write_json(ptx, here::here("data", "pars", iso, "pars_tx.json"), digits = 10, auto_unbox = T)
  
}



iso = glue::as_glue("ZAF")

load(here::here("data", iso, "d_pop.rdata"))
load(here::here("data", iso, "d_tx.rdata"))
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
    MorR_eps_HIV = (Mor_U_HIV - Mor_L_HIV) / 2 / 1.96,
    MorR_eps_NonHIV = (Mor_U_NonHIV - Mor_L_NonHIV) / 2 / 1.96,
  ) %>% 
  select(Year, starts_with(c("Pop_", "IncR_", "MorR_mu_", "MorR_eps_", "Case", "Amp"))) %>% 
  pivot_longer(-Year) %>% 
  extract(name, c("Index", "Group"), "(\\S+)_(HIV|NonHIV)") %>% 
  filter(!is.na(Index)) %>% 
  pivot_wider(names_from = Index) %>% 
  mutate(
    CNR_mu = Case / Pop,
    CNR_eps = sqrt(CNR_mu * (1 - CNR_mu) / Pop)
  )


d_gp



