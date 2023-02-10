library(tidyverse)


theme_set(theme_bw())



country = "IND"
iso = glue::as_glue(country)




mss0 <- read_csv(here::here("results", "I_" + iso, "RunBaseline.csv"))
mss1 <- read_csv(here::here("results", "I_" + iso, "RunIntv.csv"))


mss = bind_rows(mss0, mss1) %>% 
  select(Time, IncR, MorR, CNR, Scenario) %>%
  pivot_longer(-c(Time, Scenario), names_to = "Index") %>% 
  group_by(Time, Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  )


avt <- mss1 %>% 
  select(Time, CumInc, CumMor, Key, Scenario) %>% 
  group_by(Key, Scenario) %>% 
  mutate(
    CumInc = CumInc - CumInc[1],
    CumMor = CumMor - CumMor[1]
  ) %>% 
  ungroup()


avt <- avt %>% 
  left_join(avt %>% 
              filter(Scenario == "Baseline") %>% 
              select(Time, CumInc0 = CumInc, CumMor0 = CumMor, Key)) %>% 
  mutate(
    AvtInc = (1 - CumInc / CumInc0),
    AvtMor = (1 - CumMor / CumMor0),
    AvtInc = ifelse(is.na(AvtInc), 0, AvtInc),
    AvtMor = ifelse(is.na(AvtMor), 0, AvtMor)
  ) %>% 
  select(Time, AvtInc, AvtMor, Scenario) %>%
  pivot_longer(-c(Time, Scenario), names_to = "Index") %>% 
  group_by(Time, Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  )



mss %>% 
  ggplot() +
  #geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  facet_wrap(Index~., scales = "free_y", nrow=3) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  expand_limits(y = 0)


avt %>% 
  ggplot() +
  #geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  facet_wrap(Index~., scales = "free_y", nrow=3) +
  scale_y_continuous("%", labels = scales::percent_format()) + 
  expand_limits(y = 0)
  




