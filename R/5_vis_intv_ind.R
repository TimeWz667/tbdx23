library(tidyverse)


theme_set(theme_bw())



country = "IND"
iso = glue::as_glue(country)


t_end <- 2030

scs <- c(
  Baseline = "Baseline",
  Dx90 = "Increased diagnostic uptake in healthcare settings",
  PPM90 = "Increased diagnostic uptake, with PPM",
  #"Dx+PPM90" = "90% dx and 90% PPM",
  RedDelay2 = "half care-seeking delays",
  ACF2 = "Proactive case-finding, symptomatic TB",
  PPM90_RedDelay2 = "Combined: 90% ppm + half cs",
  AsymACF = "ACF, 10% of asymptomatic TB in vulnerable popn",
  Dx90_ACF2 = "Combined: 90% dx + ACF 2 times per year",
  PPM90_ACF2 = "All measures combined"
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


mss1 %>% 
  filter(Scenario == "ACF2") %>% 
  select(Time, CumACF, Key) %>% 
  filter(Time %in% c(2030, 2035)) %>% 
  group_by(Time) %>% 
  summarise(
    M = median(CumACF),
    L = quantile(CumACF, 0.25),
    U = quantile(CumACF, 0.75)
  ) %>% 
  mutate(
    across(c(M, L, U), function(x) x / (Time - 2023))
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
  filter(Time <= t_end) %>% 
  #filter(Time > 2018) %>% 
  ggplot() +
  #geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U)) + 
  geom_point(data = tar %>% select(Year, CNR = CNR_mu) %>% mutate(Index = "CNR"), 
             aes(x = Year, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = c(2015, 2023, 2025, 2030, 2035)) +
  scale_fill_brewer(palette="Set1", labels = scs) +
  expand_limits(y = 0) +
  facet_wrap(Index~., scales = "free_y",  
             labeller = labeller(Index=c(CNR="Case notification rate", IncR = "Incidence", MorR = "Mortality"))) +
  theme(legend.position = "bottom", legend.direction = "vertical")


g_intv


g_intv0 <- mss %>% 
  filter(Index %in% c("IncR", "MorR")) %>% 
  filter(Time > 2021) %>% 
  filter(Time <= t_end) %>% 
  filter(!Scenario %in% c("RedDelay2", "PPM90_RedDelay2", "Dx90_ACF2")) %>% 
  ggplot() +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = c(2015, 2023, 2025, 2030, 2035)) +
  scale_colour_brewer(palette="Dark2", labels = scs) +
  expand_limits(y = 0) +
  facet_wrap(Index~., scales = "free_y",  
             labeller = labeller(Index=c(CNR="Case notification rate", IncR = "Incidence", MorR = "Mortality"))) +
  theme(legend.position = "bottom", legend.direction = "vertical")


g_intv0



g_avt <- avt %>% 
  filter(!Scenario %in% c("RedDelay2", "PPM90_RedDelay2", "Dx90_ACF2")) %>% 
  filter(Time <= t_end) %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  scale_y_continuous("%", labels = scales::percent_format()) +
  scale_x_continuous("Year", breaks = c(2023, 2025, 2030, 2035)) +
  scale_fill_brewer(palette="Dark2", labels = scs) +
  scale_colour_brewer(palette="Dark2", labels = scs) +
  facet_wrap(Index~., scales = "free_y", labeller = labeller(Index=c(AvtInc = "Incidence", AvtMor = "Mortality"))) +
  expand_limits(y = 0) +
  theme(legend.position = "bottom", legend.direction = "vertical")

g_avt  



ggsave(g_intv, filename = here::here("results", "figs", "g_intv_ind.png"), width = 8, height = 4.5)
ggsave(g_intv0, filename = here::here("results", "figs", "g_intv0_ind.png"), width = 6, height = 4.5)
ggsave(g_avt, filename = here::here("results", "figs", "g_avt_ind.png"), width = 6, height = 4.5)
