library(tidyverse)


theme_set(theme_bw())


ci_range <- 0.95
ru <- c((1 - ci_range) / 2, (1 + ci_range) / 2)
t_end <- 2030


scs <- c(
  Baseline = "Baseline",
  Dx90 = "New diagnostics in healthcare settings",
  PPM90 = "New diagnostics in healthcare settings, with PPM",
  #"Dx+PPM90" = "90% dx and 90% PPM",
  #RedDelay2 = "(C) half care-seeking delays",
  ACF2 = "Proactive case-finding, symptomatic TB",
  AsymACF = "Detecting 10% of asymptomatic TB in vulnerable population annually",
  PPM90_ACF2 = "Combined measures: (B) + (C) + (D)",
  Dx90_ACF2 = "Combined measures: (A) + (C) + (D)",
  Combined = "All measures combined"
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
  mss1 <- read_csv(here::here("results", "I_" + iso, "RunIntv80.csv"))
  
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


tab_avt <- d_avt %>% 
  filter(Time == 2030) %>% 
  mutate(
    Scenario = scs[Scenario],
    across(c(M, L, U), scales::percent_format(accuracy = 0.1)),
    mlu = paste0(M, " (", L, " - ", U, ")")
  ) %>% 
  arrange(Country, Index, Scenario, Index)


 
# write_csv(tab_avt, here::here("results", "avt.csv"))



d_prev <- bind_rows(local({
  load("data/IND/d_prev.rdata")
  d_prev
}), local({
  load("data/ZAF/d_prev.rdata")
  d_prev
})) %>% filter(State != "PreTx" & Tag == "All")


g_tbps <- d_prev %>% 
  mutate(
    Prev = N_Prev / N_Subject,
    State = factor(State, rev(c("Asym", "Sym", "CS")))
  ) %>% 
  ggplot() + 
  geom_bar(aes(x = Prev, y = 1, fill = State), stat = "identity", position = "fill") +
  scale_x_continuous("Proportion of untreated TB, %", labels = scales::percent_format()) +
  scale_y_discrete("", labels = " ") +
  scale_fill_discrete("State", 
                      labels=c(Asym = "Subclinical", 
                               Sym = "Pre care-seeking", 
                               CS = "Sought care"), guide = guide_legend(reverse = TRUE)) +
  facet_grid(.~Country) +
  theme(legend.position = "bottom")

g_tbps


g_tbpsv <- d_prev %>% 
  mutate(
    Country = factor(Country, c("South Africa", "India")),
    Prev = N_Prev / N_Subject,
    State = factor(State, rev(c("Asym", "Sym", "CS")))
  ) %>% 
  ggplot() + 
  geom_bar(aes(x = Prev, y = Country, fill = State), stat = "identity", position = "fill") +
  scale_x_continuous("Proportion of untreated TB, %", labels = scales::percent_format()) +
  scale_y_discrete("") +
  scale_fill_discrete("State", labels=c(Asym = "Subclinical", Sym = "Pre care-seeking", CS = "Sought care"), guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom")





g_bind2 <- ggpubr::ggarrange(g_tbps, g_avt, nrow = 2, heights = c(1, 4))
g_bind3 <- ggpubr::ggarrange(g_tbpsv, g_avt, nrow = 2, heights = c(1, 4))

# ggsave(g_avt, filename = here::here("results", "figs", "g_avt.png"), width = 8, height = 7.5)
# ggsave(g_bind2, filename = here::here("results", "figs", "g_avt_v2.png"), width = 8, height = 10.5)
# ggsave(g_bind3, filename = here::here("results", "figs", "g_avt_v3.png"), width = 8, height = 10.5)
# 



ps <- lapply(c("India", "South Africa"), function(cnt) {
  d_prev %>% 
    mutate(
      Country = factor(Country, c("South Africa", "India")),
      Prev = N_Prev / N_Subject,
      State = factor(State, rev(c("Asym", "Sym", "CS")))
    ) %>% 
    filter(Country == cnt) %>% 
    ggplot() + 
    geom_bar(aes(x = Prev, y = 1, fill = State), stat = "identity", position = "fill") +
    scale_x_continuous("Proportion of untreated TB, %", labels = scales::percent_format()) +
    scale_y_discrete("", labels = " ") +
    scale_fill_discrete("Cascade", 
                        labels=c(Asym = "Subclinical", 
                                 Sym = "Pre care-seeking", 
                                 CS = "Sought care"), guide = guide_legend(reverse = TRUE)) +
    facet_grid(.~Country) +
    theme(legend.position = "bottom") +
    labs(subtitle = ifelse(cnt == "India", "(A)", "(B)"))
})


g_tbps <- ggpubr::ggarrange(plotlist = ps, common.legend = T, legend = "bottom")



label.panel <- tibble(
  Index = c("AvtInc", "AvtInc", "AvtMor", "AvtMor"),
  Country = c("IND", "ZAF", "IND", "ZAF"),
  Lab = c("(C)", "(D)", "(E)", "(F)"),
  y = rep(c(0.22, 0.33), each = 2)
)



g_avt <- d_avt %>% 
  ggplot() +
  geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Time, y = M, colour = Scenario)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(x = 2030, y = 0, label = "No intervention", vjust = 1.2, hjust = 1) +
  geom_text(data = label.panel, aes(x = 2022.5, y = y, label = Lab)) +
  scale_y_continuous("", labels = scales::percent_format()) +
  scale_x_continuous("Year", breaks = c(2023, 2025, 2030, 2035)) +
  scale_fill_brewer(palette="Dark2", labels = scs) +
  scale_colour_brewer(palette="Dark2", labels = scs) +
  facet_grid(Index~Country, scales = "free_y", 
             labeller = labeller(Index=c(AvtInc = "Cumulative Cases Averted", AvtMor = "Cumulative Deaths Averted"),
                                 Country=c(IND = "India", ZAF = "South Africa"))) +
  expand_limits(y = -0.03) +
  theme(legend.position = "bottom", legend.direction = "vertical")

g_bind4 <- ggpubr::ggarrange(g_tbps, g_avt, nrow = 2, heights = c(1, 4))


as <- apply(label.panel, 1, function(gp) {
  d_avt %>% 
    filter(Country == gp["Country"] & Index == gp["Index"]) %>% 
    ggplot() +
    geom_ribbon(aes(x = Time, ymin = L, ymax = U, fill = Scenario), alpha = 0.2) +
    geom_line(aes(x = Time, y = M, colour = Scenario)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_text(x = 2030, y = 0, label = "No intervention", vjust = 1.2, hjust = 1) +
    scale_y_continuous("", labels = scales::percent_format()) +
    scale_x_continuous("Year", breaks = c(2023, 2025, 2030, 2035)) +
    scale_fill_brewer(palette="Dark2", labels = scs) +
    scale_colour_brewer(palette="Dark2", labels = scs) +
    facet_wrap(Country ~ Index, scales = "free_y", 
               labeller = labeller(Index=c(AvtInc = "Cumulative Cases Averted", AvtMor = "Cumulative Deaths Averted"),
                                   Country=c(IND = "India", ZAF = "South Africa"))) +
    expand_limits(y = -0.03) +
    theme(legend.position = "bottom", legend.direction = "vertical") +
    labs(subtitle = gp["Lab"])
})

g_avt <- ggpubr::ggarrange(plotlist = as, common.legend = T, legend = "bottom")

g_bind5 <- ggpubr::ggarrange(g_tbps, g_avt, nrow = 2, heights = c(1, 4))

ggsave(g_bind5, filename = here::here("results", "figs", "g_avt_v5.png"), width = 8, height = 10.5)
