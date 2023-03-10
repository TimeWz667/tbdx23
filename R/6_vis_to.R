library(tidyverse)


to <- read_csv(here::here("data", "TB_outcomes.csv"))


to %>% 
  filter(iso3 %in% c("IND", "ZAF")) %>% 
  filter(year >= 2015) %>% 
  select(country, year, newrel_succ, newrel_coh, newrel_died) %>% 
  group_by(country) %>%
  summarise(
    p_succ = sum(newrel_succ) / sum(newrel_coh),
    p_died = sum(newrel_died) / sum(newrel_coh)
  )





to <- read_csv(here::here("data", "ITR_Tx.csv"))



to %>% 
  filter(State == "India") %>% 
  mutate(
    across(starts_with("N_Tx"), function(x) as.numeric(gsub("\\s+\\S+", "", x)))
  ) %>%
  summarise(
    p_succ_pub = sum(N_Tx_Succ_Pub) / sum(N_Tx_Ini_Pub),
    p_die_pub = sum(N_Tx_Die_Pub) / sum(N_Tx_Ini_Pub),
    p_succ_pri = sum(N_Tx_Succ_Pri) / sum(N_Tx_Ini_Pri),
    p_die_pri = sum(N_Tx_Die_Pri) / sum(N_Tx_Ini_Pri)
  )
