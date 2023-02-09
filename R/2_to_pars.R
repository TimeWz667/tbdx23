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
