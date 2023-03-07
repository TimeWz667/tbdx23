library(tidyverse)


theme_set(theme_bw())


ci_range <- 0.95
ru <- c((1 - ci_range) / 2, (1 + ci_range) / 2)
t_end <- 2030


scs <- c(
  Baseline = "Baseline",
  Dx90 = "(A) Increased diagnostic uptake in healthcare settings",
  PPM90 = "(B) Increased diagnostic uptake, with PPM",
  #"Dx+PPM90" = "90% dx and 90% PPM",
  #RedDelay2 = "(C) half care-seeking delays",
  ACF2 = "(C) Proactive case-finding, symptomatic TB",
  AsymACF = "(D) ACF, 10% of asymptomatic TB in vulnerable popn",
  PPM90_ACF2 = "Combined measures: (B) + (C) + (D)",
  Dx90_ACF2 = "Combined measures: (A) + (C) + (D)",
  Combined = "All measures combined"
)


sel_scs <- list(
  IND = c(
    Baseline = "Baseline",
    Dx90 = "Increased diagnostic uptake in healthcare settings",
    PPM90 = "Increased diagnostic uptake, with PPM",
    ACF2 = "Proactive case-finding, symptomatic TB",
    AsymACF = "ACF, 10% of asymptomatic TB in vulnerable popn",
    PPM90_ACF2 = "All measures combined"
  ),
  ZAF = c(
    Baseline = "Baseline",
    Dx90 = "Increased diagnostic uptake in healthcare settings",
    #RedDelay2 = "(C) half care-seeking delays",
    ACF2 = "Proactive case-finding, symptomatic TB",
    AsymACF = "ACF, 10% of asymptomatic TB in vulnerable popn",
    Dx90_ACF2 = "All measures combined"
  )
)


d_ts <- bind_rows(lapply(c(IND = "IND", ZAF = "ZAF"), function(country) {
  iso <- glue::as_glue(country)
  
  
  mss0 <- read_csv(here::here("results", "I_" + iso, "RunBaseline.csv"))
  mss1 <- read_csv(here::here("results", "I_" + iso, "RunIntv.csv"))
  
  mss <- bind_rows(mss0, mss1) %>% 
    select(Time, IncR, MorR, CNR, Scenario) %>%
    pivot_longer(-c(Time, Scenario), names_to = "Index") %>% 
    group_by(Time, Scenario, Index) %>% 
    summarise(
      M = median(value),
      L = quantile(value, ru[1]),
      U = quantile(value, ru[2])
    ) %>% 
    ungroup() %>% 
    filter(Scenario %in% names(sel_scs[[country]])) %>% 
    mutate(
      Scenario = ifelse(Scenario %in% c("PPM90_ACF2", "Dx90_ACF2"), "Combined", Scenario),
      Scenario = factor(Scenario, names(scs)),
      Country = country
    ) %>% 
    filter(Time <= t_end)
  
}))


tar <- bind_rows(lapply(c(IND = "IND", ZAF = "ZAF"), function(country) {
  read_csv(here::here("data", "pars", country, "targets.csv")) %>% 
    mutate(Country = country)
  
})) %>% select(Year, Country, CNR = CNR_mu, IncR = IncR_mu, MorR = MorR_mu) %>% 
  pivot_longer(c(CNR, IncR, MorR), names_to = "Index") %>% 
  filter(Year > 2015)


g_ts <- d_ts %>% 
  filter(Time > 2015) %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  # geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U)) + 
  geom_point(data = tar, 
             aes(x = Year, y = value)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  scale_x_continuous("Year", breaks = c(2015, 2023, 2025, 2030, 2035)) +
  scale_fill_brewer(palette="Dark2", labels = scs) +
  scale_colour_brewer(palette="Dark2", labels = scs) +
  expand_limits(y = 0) +
  facet_grid(Index~Country, scales = "free_y", 
             labeller = labeller(Index=c(CNR="Case notification rate", IncR = "Incidence", MorR = "Mortality"),
                                 Country=c(IND = "India", ZAF = "South Africa"))) +
  theme(legend.position = "bottom", legend.direction = "vertical")



ggsave(g_ts, filename = here::here("results", "figs", "g_ts.png"), width = 8, height = 7.5)


