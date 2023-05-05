library(tidyverse)


theme_set(theme_bw())


ci_range <- 0.95
ru <- c((1 - ci_range) / 2, (1 + ci_range) / 2)
t_end <- 2030


scs <- c(
  Baseline = "Baseline",
  ACF2_b = "Sens = 50%, Sy > SSM",
  ACF2_cxr = "Sens = 80%, CXR > Xpert",
  #"Dx+PPM90" = "90% dx and 90% PPM",
  #RedDelay2 = "(C) half care-seeking delays",
  ACF2_xpert = "Sens = 95%, (CXR, Sy) > Xpert"
)


sel_scs <- list(
  IND = c(
    ACF2_b = "SSM",
    ACF2_cxr = "xpert",
    ACF2_xpert = "cxr + xpert"
  ),
  ZAF = c(
    ACF2_b = "SSM",
    ACF2_cxr = "xpert",
    ACF2_xpert = "cxr + xpert"
  )
)



d_avt <- bind_rows(lapply(c(IND = "IND", ZAF = "ZAF"), function(country) {
  iso <- glue::as_glue(country)
  
  mss1 <- read_csv(here::here("results", "I_" + iso, "RunSens.csv"))
  
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
    filter(Scenario %in% names(sel_scs[[country]])) %>% 
    mutate(
      Scenario = ifelse(Scenario %in% c("PPM90_ACF2", "Dx90_ACF2"), "Combined", Scenario),
      Scenario = factor(Scenario, names(scs)),
      Country = country
    ) %>% 
    
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
  theme(legend.position = "bottom", legend.direction = "vertical")# + 
#labs(caption = "India: (B) + (C) + (D)\nSouth Africa: (A) + (C) + (D)")

g_avt



