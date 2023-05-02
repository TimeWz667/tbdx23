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


load(here::here("data", "IND", "d_prev.rdata"))


mss0 <- read_csv(here::here("results", "C_" + iso, "RunPost.csv"))
mss1 <- read_csv(here::here("results", "D_" + iso, "RunPost.csv")) 


## Pre-COVID epi
g_pre <- mss0 %>% 
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
  geom_pointrange(data = d2plot %>% filter(Year <= 2019), aes(x = Year, y = M, ymin = L, ymax = U)) + 
  geom_point(data = tar %>% select(Year, CNR = CNR_mu) %>% mutate(Index = "CNR"), 
             aes(x = Year, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y", labeller = labeller(Index=c(IncR="Incidence", MorR="Mortality", CNR="Case notification"))) +
  expand_limits(y = 0)

g_pre


## COVID disruption

g_covid <- bind_rows(
  mss1 %>% mutate(CNR = CNR)
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
  geom_point(data = tarq %>% mutate(CNR = QCNR, Index = "CNR"), 
             aes(x = Time, y = CNR)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y", labeller = labeller(Index=c(CNR="Case notification, quarterly"))) +
  expand_limits(y = 0)


g_covid


## TBPS


prev_pr <- d_prev %>% 
  ungroup() %>% 
  filter(Tag == "All" & State != "PreTx") %>% 
  mutate(
    Index = State,
    N = sum(N_Prev),
    M = N_Prev / N,
    L = qbinom(0.025, N, M) / N,
    U = qbinom(0.975, N, M) / N
  ) %>%
  mutate(Index = factor(Index, levels = c("Asym", "Sym", "CS")))


g_prev <- bind_rows(mss0, mss1) %>% 
  filter(Time >= 2018 & Time <= 2021) %>%
  mutate(
    Asym = PrevA / Prev,
    Sym = PrevS / Prev,
    CS = PrevC / Prev
  ) %>% 
  select(Year = Time, Asym, Sym, CS, Key) %>% 
  pivot_longer(c(Asym, Sym, CS), names_to = "Index") %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  ungroup() %>%
  mutate(Index = factor(Index, c("Asym", "Sym", "CS"))) %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_pointrange(data = prev_pr, aes(x = 2019, y = M, ymin = L, ymax = U)) + 
  scale_y_continuous("percent", labels = scales::percent) + 
  facet_wrap(.~Index, scales = "free_y", labeller = labeller(Index=c(Asym="Asymptomatic", Sym="Before care-seeking", CS="Sought care"))) +
  expand_limits(y = c(0, 1))

g_prev



g_gof <- ggpubr::ggarrange(
  g_pre + labs(subtitle = "Pre-COVID-19 TB epidemiology"),
  g_prev + labs(subtitle = "TB Prevalence by states"),
  nrow = 2
)



ggsave(g_gof, filename = here::here("results", "figs", "g_gof_ind.png"), width = 8, height = 7)





