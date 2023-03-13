library(tidyverse)


theme_set(theme_bw())



country = "IND"
iso = glue::as_glue(country)

rci <- 0.95
ci_l <- (1 - rci) / 2
ci_u <- 1 - ci_l


tar_all <- read_csv(here::here("pars", iso, "targets.csv"))
tar_pupr <- read_csv(here::here("pars", iso, "targets_pupr.csv")) %>% filter(Year >= 2017)
tar_prev <- read_csv(here::here("pars", iso, "targets_prev.csv"))
tar_q  <- read_csv(here::here("pars", iso, "targets_q.csv"))



mss0 <- read_csv(here::here("out", "A_" + iso, "RunPost.csv"))

mss1 <- read_csv(here::here("results", "B_" + iso, "RunPost.csv")) 


# mss <- bind_rows(
#   mss0 %>% filter(Time <= 2020),
#   mss1 %>% mutate(CNR = CNR * 4)
# )

mss <- mss0

stats <- mss %>% 
  pivot_longer(-c(Year, Key), names_to = "Index") %>% 
  filter(!is.na(value)) %>% 
  group_by(Year, Index) %>% 
  summarise(
    M = mean(value),
    L = quantile(value, ci_l),
    U = quantile(value, ci_u)
  ) %>% 
  ungroup()


stats_ts <- stats %>% 
  filter(Index %in% c("IncR", "MorR", "CNR"))
  

stats_prev <- stats %>% 
  filter(Index %in% c("PrevA", "PrevS", "PrevC"))


stats_pupr <- stats %>% 
  filter(Index %in% c("CNR_Pub", "CNR_Pri")) %>% 
  separate(Index, c("Index", "Tag"), "_")



stats_ts %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_pointrange(data = tar_all, aes(x = Year, y = m, ymin = m - eps, ymax = m + eps)) + 
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)


stats_prev %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_pointrange(data = tar_prev, aes(x = Year, y = m, ymin = m - eps, ymax = m + eps)) + 
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(.~Index, scales = "free_y") +
  expand_limits(y = 0)


stats_pupr %>% 
  ggplot() +
  geom_ribbon(aes(x = Year, ymin = L, ymax = U), alpha = 0.2) +
  geom_line(aes(x = Year, y = M)) + 
  geom_pointrange(data = tar_pupr, aes(x = Year, y = m, ymin = m - eps, ymax = m + eps)) + 
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
  facet_wrap(Tag~Index, scales = "free_y") +
  expand_limits(y = 0)

