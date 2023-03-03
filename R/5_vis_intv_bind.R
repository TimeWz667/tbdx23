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
  Dx90_ACF2 = "Combined measures: (A) + (C) + (D)"
)


sel_scs <- list(
  IND = c(
    Dx90 = "Increased diagnostic uptake in healthcare settings",
    PPM90 = "Increased diagnostic uptake, with PPM",
    ACF2 = "Proactive case-finding, symptomatic TB",
    AsymACF = "ACF, 10% of asymptomatic TB in vulnerable popn",
    PPM90_ACF2 = "All measures combined"
  ),
  ZAF = c(
    Dx90 = "Increased diagnostic uptake in healthcare settings",
    #RedDelay2 = "(C) half care-seeking delays",
    ACF2 = "Proactive case-finding, symptomatic TB",
    AsymACF = "ACF, 10% of asymptomatic TB in vulnerable popn",
    Dx90_ACF2 = "All measures combined"
  )
)


d_avt <- bind_rows(lapply(c(IND = "IND", ZAF = "ZAF"), function(country) {
  iso <- glue::as_glue(country)
  
  
  mss0 <- read_csv(here::here("results", "I_" + iso, "RunBaseline.csv"))
  mss1 <- read_csv(here::here("results", "I_" + iso, "RunIntv.csv"))
  
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
      L = quantile(value, ru[1]),
      U = quantile(value, ru[2])
    ) %>% 
    ungroup() %>% 
    mutate(
      Scenario = factor(Scenario, names(scs)),
      Country = country
    ) %>% 
    filter(Scenario %in% names(sel_scs[[country]])) %>% 
    filter(Time <= t_end)
  
}))


g_avt <- d_avt %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(x = 2030, y = 0, label = "No intervention", vjust = 1.2, hjust = 1) + 
  scale_y_continuous("", labels = scales::percent_format()) +
  scale_x_continuous("Year", breaks = c(2023, 2025, 2030, 2035)) +
  scale_fill_brewer(palette="Dark2", labels = scs) +
  scale_colour_brewer(palette="Dark2", labels = scs) +
  facet_grid(Index~Country, scales = "free_y", 
             labeller = labeller(Index=c(AvtInc = "Cumulative Cases Averted", AvtMor = "Cumulative Deaths Averted"),
                                 Country=c(IND = "India", ZAF = "South Africa"))) +
  expand_limits(y = -0.03) +
  theme(legend.position = "bottom", legend.direction = "vertical")



ggsave(g_avt, filename = here::here("results", "figs", "g_avt.png"), width = 8, height = 7.5)








