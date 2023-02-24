

dat <- read_csv(here::here("data", "TB_provisional_notifications.csv"))



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





cnr <- bind_rows(lapply(c("ZAF", "IND"), function(iso) {
  tar <- read_csv(here::here("data", "pars", iso, "targets.csv"))
  tarq  <- read_csv(here::here("data", "pars", iso, "targets_q.csv"))
  
  bind_rows(
    tar %>% 
      mutate(
        Time = Year + 0.5,
        QCNR = CNR_mu / 4,
        q = 2.5, iso = iso
      ) %>% 
      select(iso, Year, q, Time, CNR = CNR_mu, QCNR),
    
    tarq %>% 
      mutate(
        CNR = QCNR * 4
      ) %>% 
      select(iso, Year, q, Time, CNR, QCNR)
  )
}))




cnr %>% 
  ggplot() + 
  geom_point(aes(x = Time, y = CNR, colour = iso)) +
  geom_line(aes(x = Time, y = CNR, group = iso))  + 
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  expand_limits(y=0)






