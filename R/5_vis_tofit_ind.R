library(tidyverse)


theme_set(theme_bw())



country = "IND"
iso = glue::as_glue(country)


d2plot <- local({
  load("data/" + iso + "/d_burden.rdata")
  
  d_burden %>% 
    select(-Country)  %>% 
    mutate(
      Inc_M = Inc_M,
      Inc_L = Inc_L,
      Inc_U = Inc_U
    )%>% 
    filter(Year >= 2014) %>% 
    pivot_longer(- Year) %>% 
    separate(name, c("Index", "name")) %>% 
    pivot_wider() %>% 
    mutate(Index = paste0(Index, "R"))
})


tar <- read_csv(here::here("data", "pars", iso, "targets.csv"))
tarq  <- read_csv(here::here("data", "pars", iso, "targets_q.csv"))



mss0 <- read_csv(here::here("results", "C_" + iso, "RunPost.csv"))
mss1 <- read_csv(here::here("results", "D_" + iso, "RunPost.csv")) 

mss0 %>% 
  select(Year = Time, IncR, MorR, CNR, Key) %>% 
  pivot_longer(c(IncR, MorR, CNR), names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ungroup() %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U)) + 
  geom_point(data = tar %>% select(Year, CNR = CNR_mu) %>% mutate(Index = "CNR"), 
             aes(x = Year, y = CNR)) +
  geom_point(data = tarq %>% mutate(CNR = QCNR * 4, Index = "CNR"), 
             aes(x = Time, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)


mss1 %>% 
  select(Year = Time, IncR, MorR, CNR, Key) %>% 
  pivot_longer(c(IncR, MorR, CNR), names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ungroup() %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_point(data = tarq %>% mutate(CNR = QCNR, Index = "CNR"), 
             aes(x = Time, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)

bind_rows(
  mss0 %>% filter(Time <= 2020),
  mss1 %>% mutate(CNR = CNR * 4)
) %>% 
  select(Year = Time, CNR = CNR, Key) %>% 
  pivot_longer(c(CNR), names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ungroup() %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_point(data = tar %>% select(Year, CNR = CNR_mu) %>% mutate(Index = "CNR"), 
             aes(x = Year, y = CNR)) +
  geom_point(data = tarq %>% mutate(CNR = QCNR * 4, Index = "CNR"), 
             aes(x = Time, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)



bind_rows(
  mss0 %>% filter(Time <= 2020),
  mss1 %>% mutate(CNR = CNR * 4)
) %>% 
  select(Year = Time, IncR = IncR_apx, MorR = MorR_apx, Key) %>% 
  pivot_longer(c(IncR, MorR), names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ungroup() %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)



