library(tidyverse)



db <- read_csv(here::here("data", "TB_burden_countries.csv"))

d_burden_all <- db %>% 
  filter(iso3 %in% c("IND", "ZAF")) %>% 
  select(
    Year = year, Country = country,
    Inc_M = e_inc_100k, Inc_L = e_inc_100k_lo, Inc_U = e_inc_100k_hi,
    Mor_M = e_mort_100k, Mor_L = e_mort_100k_lo, Mor_U = e_mort_100k_hi
  ) %>% 
  mutate(
    across(Inc_M:Mor_U, function(x) x * 1e-5)
  )


d_burden <- d_burden_all %>% filter(Country == "India")
save(d_burden, file = here::here("data", "IND", "d_burden.rdata"))
