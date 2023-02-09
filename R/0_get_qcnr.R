library(tidyverse)



dat <- read_csv(here::here("data", "TB_notifications.csv"))



ind <- dat %>% 
  filter(iso3 == "IND") %>% 
  select(Country = country, iso = iso3, Year = year, starts_with("m_")) %>% 
  pivot_longer(starts_with("m_")) %>% 
  extract(name, "month", "m_(\\d+)") %>% 
  mutate(
    q = ceiling(as.numeric(month) / 3),
    Time = Year + (q - 1) / 4 + 1 / 8
  ) %>% 
  group_by(Country, iso, Time, Year, q) %>% 
  summarise(value = sum(value)) %>% 
  filter(!is.na(value))


zaf <- dat %>% 
  filter(iso3 == "ZAF") %>% 
  select(Country = country, iso = iso3, Year = year, starts_with("q_")) %>% 
  pivot_longer(starts_with("q_")) %>% 
  extract(name, "q", "q_(\\d+)") %>% 
  mutate(
    q = as.numeric(q),
    Time = Year + (q - 1) / 4 + 1 / 8
  ) %>% 
  group_by(Country, iso, Time, Year, q) %>% 
  summarise(value = sum(value))


pop <- bind_rows(list(
  local({
    load(here::here("data", "IND", "d_pop.rdata"))
    
    as_tibble(approx(d_pop_all$Year, d_pop_all$N_Pop, xout = ind$Time)) %>% 
      mutate(
        iso = "IND"
      ) %>% 
      rename(Time = x, Pop = y)
    
  }),
  local({
    load(here::here("data", "ZAF", "d_pop.rdata"))
    
    as_tibble(approx(d_pop_all$Year, d_pop_all$N_Pop, xout = zaf$Time)) %>% 
      mutate(
        iso = "ZAF"
      ) %>% 
      rename(Time = x, Pop = y)
    
  })
))


qcnr <- bind_rows(ind, zaf) %>% 
  left_join(pop) %>% 
  mutate(QCNR = value / Pop, CNR = QCNR * 4)

qcnr %>% 
  ggplot() +
  geom_line(aes(x = Time, y = CNR, colour = iso))



qcnr %>% 
  group_by(Country, iso, Year) %>% 
  summarise(CNR = sum(value) / mean(Pop) * (4 / n()))

qcnr %>% filter(iso == "IND") %>% 
  select(Country, Time, Year, q, QCNR) %>% 
  write_csv(file = here::here("data", "pars", "IND", "targets_q.csv"))


qcnr %>% filter(iso == "ZAF") %>% 
  select(Country, Time, Year, q, QCNR) %>% 
  write_csv(file = here::here("data", "pars", "ZAF", "targets_q.csv"))

