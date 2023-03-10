library(tidyverse)



d_tx_itr <- read_csv(here::here("data", "ITR_Tx.csv")) %>% 
  filter(State == "India") %>% 
  mutate(
    across(starts_with("N_Tx"), function(x) as.numeric(gsub("\\s+\\S+", "", x)))
  ) %>%
  summarise(
    p_succ_pub = sum(N_Tx_Succ_Pub) / sum(N_Tx_Ini_Pub),
    p_die_pub = sum(N_Tx_Die_Pub) / sum(N_Tx_Ini_Pub),
    p_succ_pri = sum(N_Tx_Succ_Pri) / sum(N_Tx_Ini_Pri),
    p_die_pri = sum(N_Tx_Die_Pri) / sum(N_Tx_Ini_Pri),
    p_succ_all = sum(N_Tx_Succ_Pub + N_Tx_Succ_Pri) / sum(N_Tx_Ini_Pub + N_Tx_Ini_Pri),
    p_die_all = sum(N_Tx_Die_Pub + N_Tx_Die_Pri) / sum(N_Tx_Ini_Pub + N_Tx_Ini_Pri)
  )


d_tx_itr


save(d_tx_itr, file = here::here("data", iso, "d_tx_itr.rdata"))

