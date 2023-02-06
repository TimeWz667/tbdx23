library(tidyverse)
library(jsonlite)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


source(here::here("r", "fn_cas.R"))



pars_ind <- read_csv(here::here("results", "pars_IND.csv"))
pars_zaf <- read_csv(here::here("results", "pars_ZAF.csv"))

js <- list(
  IND = pars_ind[1:300, ],
  ZAF = pars_zaf[1:300, ]
)


# jsonlite::write_json(js, here::here("results/pars.json"), digits = 10)


pars <- bind_rows(
  pars_ind %>% mutate(Country = "IND"),
  pars_zaf %>% mutate(Country = "ZAF")
) %>% 
  calc_prev()


pars0 <- pars %>% 
  calc_reform(pdx0 = 0.4, pdx1 = 0.7)


p_cas0 <- pars0 %>% 
  calc_cascade(sc = "Baseline")


p_cas1 <- pars0 %>% 
  calc_intv(or_pdx0 = 13, or_pdx1 = 4) %>% 
  calc_cascade(sc = "90% dx")


p_cas2 <- pars0 %>% 
  calc_intv(rr_csi = 2, rr_recsi = 2) %>% 
  calc_cascade(sc = "2x care-seeking")


p_cas3 <- pars0 %>% 
  calc_intv(or_pdx0 = 13, or_pdx1 = 4, rr_csi = 2, rr_recsi = 2) %>% 
  calc_cascade(sc = "Both")



cas <- bind_rows(
  p_cas0, p_cas1, p_cas2, p_cas3
) %>% 
  select(Country, Scenario, starts_with("Delay"), starts_with("Pr_")) %>% 
  group_by(Country, Scenario) %>% 
  summarise(across(everything(), list(
    M = function(x) quantile(x, 0.5, na.rm = T),
    L = function(x) quantile(x, 0.025, na.rm = T),
    U = function(x) quantile(x, 0.975, na.rm = T)
  ))) %>% 
  pivot_longer(-c(Country, Scenario)) %>% 
  separate(name, c("Type", "Index", "name")) %>% 
  pivot_wider() %>% 
  mutate(
    Scenario = factor(Scenario, c("Baseline", "2x care-seeking", "90% dx", "Both"))
  )



g_delay <- cas %>%
  filter(Type == "Delay") %>% 
  mutate(Index = factor(Index, c("Tot", "Sys", "Pat"))) %>% 
  ggplot() +
  geom_bar(aes(y = Index, x = M, fill = Scenario), stat = "identity", alpha = 0.7, position = "dodge") +
  geom_pointrange(aes(y = Index, x = M, xmin = L, xmax = U, group = Scenario), position = position_dodge2(1)) + 
  facet_grid(. ~ Country) +
  scale_x_continuous("Delay, months", labels = scales::number_format(scale = 12)) +
  scale_y_discrete("Type", labels = c(Pat = "Patient", Sys = "System", Tot = "Total")) +
  expand_limits(x = 0)


g_cascade <- cas %>%
  filter(Type == "Pr") %>% 
  ggplot() +
  geom_bar(aes(x = Index, y = M, fill = Scenario), stat = "identity", alpha = 0.7, position = "dodge") +
  geom_pointrange(aes(x = Index, y = M, ymin = L, ymax = U, group = Scenario), position = position_dodge2(1)) + 
  facet_grid(. ~ Country) +
  scale_y_continuous("Cascade, % incidence", labels = scales::percent_format()) +
  scale_x_discrete("Type", labels = c(Det = "Detected", Notif = "Notified")) +
  expand_limits(y = c(0, 1))


g_delay
g_cascade


ggsave(g_delay, filename = here::here("results", "g_delay.png"), width = 6, height = 4)
ggsave(g_cascade, filename = here::here("results", "g_cascade.png"), width = 6, height = 4)


