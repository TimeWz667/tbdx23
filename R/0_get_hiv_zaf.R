library(tidyverse)


reform <- function(df, index) {
  df %>% 
    filter(Country == "South Africa") %>% 
    pivot_longer(-Country) %>% 
    extract(name, c("Year", "Stat"), "(\\d+)(\\S*)") %>% 
    mutate(
      Year = as.numeric(Year),
      value = as.numeric(gsub(" ", "", value)),
      Stat = case_when(
        Stat == "" ~ "M",
        Stat == "_lower" ~ "L",
        Stat == "_upper" ~ "U"
      ),
      Index = index
    ) %>% 
    filter(!is.na(Stat)) %>% 
    pivot_wider(names_from = Stat) %>% 
    select(Country:U)
}


pr_plhiv <- read_csv(here::here("data", "AIDS", "PLHIV.csv")) %>% 
  reform(index = "PrvHIV")

mor_hiv <- read_csv(here::here("data", "AIDS", "MorHIV.csv")) %>% 
  reform(index = "MorHIV")

pr_art <- read_csv(here::here("data", "AIDS", "PrART.csv")) %>% 
  reform(index = "PrART_HIV")

inc_tb <- read_csv(here::here("data", "AIDS", "IncTB_PLHIV.csv")) %>% 
  reform(index = "IncTB_HIV")

mor_tb_hiv <- read_csv(here::here("data", "AIDS", "MorTB_PLHIV.csv")) %>% 
  reform(index = "MorTB_HIV")

load(here::here("data", "ZAF", "d_pop.rdata"))


d_hiv <- bind_rows(pr_plhiv, mor_hiv, pr_art, inc_tb, mor_tb_hiv) %>% 
  select(Country, Year, Index, M) %>% 
  pivot_wider(names_from = Index, values_from = M) %>% 
  filter(Year >= 2009 & Year < 2020) %>% 
  left_join(d_pop_all) %>% 
  mutate(PrART_HIV = PrART_HIV / 100) %>% 
  select(-Amp, -Tag)

save(d_hiv, file = here::here("data", "ZAF", "d_hiv.rdata"))


