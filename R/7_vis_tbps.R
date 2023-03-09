library(tidyverse)



theme_set(theme_bw())





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
  geom_bar(aes(x = Prev, y = Country, fill = State), stat = "identity") +
  scale_x_continuous("Prevalence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  scale_fill_discrete("State", labels=c(Asym = "Subclinical", Sym = "Pre care-seeking", CS = "Sought care"), guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom") +
  labs(subtitle = "Untreated TB from TBPS")



g_tbps_cs <- d_prev %>% 
  mutate(
    Prev = N_Prev / N_Subject,
    State = factor(State, rev(c("Asym", "Sym", "CS")))
  ) %>% 
  ggplot() + 
  geom_bar(aes(x = Prev, y = Country, fill = State), stat = "identity", position = "fill") +
  scale_x_continuous("Proportion, %", labels = scales::percent_format()) +
  scale_fill_discrete("State", labels=c(Asym = "Subclinical", Sym = "Pre care-seeking", CS = "Sought care"), guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "bottom") +
  labs(subtitle = "Untreated TB from TBPS")



ggsave(g_tbps, filename = here::here("results", "g_tbps.png"), width = 5, height = 3)
ggsave(g_tbps_cs, filename = here::here("results", "g_tbps_fill.png"), width = 5, height = 3)




1.149369


d_prev %>% 
  mutate(
    Pr = N_Prev / N_Subject * rep(c(1.149369, 0.9090663), each = 3),
    L = qbinom(0.025, size = N_Subject, prob = Pr) / N_Subject,
    U = qbinom(0.975, size = N_Subject, prob = Pr) / N_Subject,                           
    
    across(c(Pr, L, U), scales::number_format(scale = 1e5))
  )

0.9090663


