library(tidyverse)




pars_cs <- read_csv(here::here("pars", "IND", "pars_cas.csv"))

summary(pars_cs %>% select(r_sym, r_aware))


pars_cs <- read_csv(here::here("pars", "ZAF", "pars_cas.csv"))

summary(pars_cs %>% select(r_sym, r_aware))

