library(tidyverse)
library(jsonlite)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))



calc_prev <- function(df) {
  df %>% 
    mutate(
      prv = prv0 * exp(- adr * (2023 - Year0)),
      ra = r_sc + r_death_a + r_death_bg,
      rs = r_sc + r_death_s + r_death_bg,
      rc = r_sc + r_death_s + r_death_bg,
      
      a0 = (rs + r_aware - adr) / r_sym,
      c0 = r_aware / (rc + r_det - adr),
      
      pr_a = a0 / (a0 + 1 + c0),
      pr_s = 1 / (a0 + 1 + c0),
      pr_c = c0 / (a0 + 1 + c0),
      
      prv_a = prv * pr_a,
      prv_s = prv * pr_s,
      prv_c = prv * pr_c,
    ) %>% 
    select(-a0, -c0)
}


calc_reform <- function(df, pdx0 = 0.4, pdx1 = 0.7) {
  df %>% 
    mutate(
      det = r_det * prv_c,
      pdx0 = pdx0,
      pdx1 = pdx1,
      r_csi = r_aware,
      det0 = r_csi * pdx0 * prv_s,
      det1 = det - det0,
      fn0 = r_csi * (1 - pdx0) * prv_s,
      r_recsi = det1 / (pdx1 * prv_c)
    ) %>% 
    select(Country, prv_a, prv_s, prv_c, r_sym, 
           r_csi, r_recsi, pdx0, pdx1, ra, rs, rc, adr, p_under)
}


calc_intv <- function(df, or_pdx0 = 1, or_pdx1 = 1, 
                      rr_csi = 1, rr_recsi = 1) {
  df %>% 
    mutate(
    odd = or_pdx0 * pdx0 / (1 - pdx0),
    pdx0 = odd / (1 + odd),
    odd = or_pdx1 * pdx1 / (1 - pdx1),
    pdx1 = odd / (1 + odd),
    r_csi = rr_csi * r_csi,
    r_recsi = rr_recsi * r_recsi
  ) %>% 
    select(-odd)
}


calc_cascade <- function(df, sc = "baseline") {
  df %>% 
    mutate(
      inc = (r_sym + ra - adr) * prv_a,
      det0 = r_csi * pdx0 * prv_s,
      det1 = r_recsi * pdx1 * prv_c,
      dur_a = 1 / (r_sym + ra),
      dur_s = 1 / (r_csi + rs),
      dur_c = 1 / (r_recsi * pdx1 + rc),
      drop_a = ra * dur_a,
      drop_s = rs * dur_s,
      drop_c = rc * dur_c,
      
      pd0 = r_csi * pdx0 * dur_s,
      pd1 = r_recsi * pdx1 * dur_c,
      k = pd0 + pd1,
      
      Delay_Pat = dur_s,
      Delay_Sys = pd1 / k * dur_c,
      Delay_Tot = Delay_Pat + Delay_Sys,
      
      Pr_Det = r_sym * dur_a * pd0,
      Pr_Det = Pr_Det + r_sym * dur_a * (r_csi * (1 - pdx0) * dur_s) * pd1,
      Pr_Notif = Pr_Det * (1 - p_under),
      Scenario = sc
    )
}




pars_ind <- read_csv(here::here("results", "pars_IND.csv"))
pars_zaf <- read_csv(here::here("results", "pars_ZAF.csv"))

js <- list(
  IND = pars_ind[1:300, ],
  ZAF = pars_zaf[1:300, ]
)


# jsonlite::write_json(js, here::here("results/pars.json"), digits = 10)


orr <- function(p0, or) {
  odd = p0 / (1 - p0) * or
  odd / (1 + odd)
}

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
    M = median,
    L = function(x) quantile(x, 0.025),
    U = function(x) quantile(x, 0.975)
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
  expand_limits(x = c(0, 1))


g_delay
g_cascade


ggsave(g_delay, filename = here::here("results", "g_delay.png"), width = 6, height = 4)
ggsave(g_cascade, filename = here::here("results", "g_cascade.png"), width = 6, height = 4)


