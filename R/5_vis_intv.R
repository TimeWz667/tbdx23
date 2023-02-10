library(tidyverse)


theme_set(theme_bw())



country = "IND"
iso = glue::as_glue(country)


scs <- c(
  Baseline = "Baseline",
  Pub90 = "90% dx in publice",
  PubEng90 = "90% dx in public and engaged private",
  RedDelay2 = "0.5% care-seeking delays",
  ACF2 = "2 times per year ACF reached"
)


tar <- read_csv(here::here("data", "pars", iso, "targets.csv"))


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
  ) %>% 
  ungroup() %>% 
  mutate(
    Scenario = factor(Scenario, names(scs))
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
  ) %>% 
  ungroup() %>% 
  mutate(
    Scenario = factor(Scenario, names(scs))
  )



g_intv <- mss %>% 
  #filter(Time > 2018) %>% 
  ggplot() +
  #geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U)) + 
  geom_point(data = tar %>% select(Year, CNR = CNR_mu) %>% mutate(Index = "CNR"), 
             aes(x = Year, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  scale_color_discrete(labels = scs) +
  expand_limits(y = 0) +
  facet_wrap(Index~., scales = "free_y",  
             labeller = labeller(Index=c(CNR="Case notification rate", IncR = "Incidence", MorR = "Mortality"))) +
  theme(legend.position = "bottom", legend.direction = "vertical")


g_intv


g_avt <- avt %>% 
  ggplot() +
  #geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  scale_y_continuous("%", labels = scales::percent_format()) +
  scale_color_discrete(labels = scs) +
  facet_wrap(Index~., scales = "free_y", labeller = labeller(Index=c(AvtInc = "Incidence", AvtMor = "Mortality"))) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom", legend.direction = "vertical")

g_avt  



ggsave(g_intv, filename = here::here("results", "figs", "g_intv_ind.png"), width = 8, height = 4.5)
ggsave(g_avt, filename = here::here("results", "figs", "g_avt_ind.png"), width = 6, height = 4.5)



mss0 %>% 
  group_by(Key, Scenario) %>% 
  summarise(adr = mean(diff(log(IncR)))) %>% 
  pull(adr) %>% hist()







