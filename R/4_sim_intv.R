library(tidyverse)
library(odin)

theme_set(theme_bw())


source(here::here("r", "fn_cas.R"))
source(here::here("r", "fn_fit.R"))



m <- odin::odin("odin/d_tb.R")
year0 <- 2010
pars <- m$private_fields$user


for (country in c("IND", "ZAF")) {
  iso <- glue::as_glue(country)
  
  
  # Load data
  load("data/" + iso + "/d_burden.rdata")
  load("data/" + iso + "/d_pop.rdata")
  
  d2plot <- d_burden %>% 
    select(-Country) %>% 
    pivot_longer(- Year) %>% 
    separate(name, c("Index", "name")) %>% 
    pivot_wider() %>% 
    mutate(Index = paste0(Index, "R"))
  
  
  pars_pop <- get_pars_pop(d_pop_all, year0)
  
  post_pars <-read_csv("results/post_pars_" + iso + ".csv") %>% 
    select(-Year0)
  
  
  # Initialise model
  p <- c(post_pars[1, ])
  p <- p[names(p) %in% pars]
  cm <- m$new(user = c(p, pars_pop))
  
  
  # Setup scenarios
  post0 <- post_pars %>% 
    calc_reform(pdx0 = 0.4, pdx1 = 0.7) %>% 
    mutate(sc = "baseline")
  
  
  post1 <- post0 %>% 
    calc_intv(or_pdx0 = 13, or_pdx1 = 4) %>% 
    mutate(sc = "90% dx")
  
  
  post2 <- post0 %>% 
    calc_intv(rr_csi = 2, rr_recsi = 2) %>% 
    mutate(sc = "2x care-seeking")
  
  
  post3 <- post0 %>% 
    calc_intv(or_pdx0 = 13, or_pdx1 = 4, rr_csi = 2, rr_recsi = 2) %>% 
    mutate(sc = "Both")
  
  
  
  
  # Simulate baseline
  ys0 <- bind_rows(lapply(1:nrow(post0), function(i) {
    p <- c(post0[i, ])
    p <- p[names(p) %in% pars]
    
    cm$set_user(user = p)
    ys <- cm$run(seq(1700, 2036, 1)) %>% 
      data.frame() %>% as_tibble() %>% 
      # select(-starts_with("Y.")) %>% 
      filter(t >= 2010) %>% 
      rename(Year = t) %>% 
      mutate(Key = i)
    
  })) %>% mutate(sc = "Baseline")
  
  
  y0s <- lapply(ys0$Key %>% unique(), function(i) {
    ys0 %>% 
      filter(Year == 2023) %>% 
      filter(Key == i) %>% 
      select(starts_with("Y.")) %>% 
      unlist() %>% 
      unname()
  })
  
  
  ys.start <- ys0 %>% 
    group_by(Key) %>% 
    select(-starts_with("Y.")) %>% 
    mutate(
      IncR = c((IncR[-n()] + IncR[-1]) / 2, 0),
      MorR = c((MorR[-n()] + MorR[-1]) / 2, 0),
      Inc = c(Inc[-1], 0),
      Mor = c(Mor[-1], 0)
    ) %>% 
    ungroup() %>% 
    filter(Year == 2022) %>% 
    select(-sc)
  
  
  sim_intv <- function(ps, ys0, y0s) {
    sc <- ps$sc[1]
    ks <- ys0$Key %>% unique()
    
    
    sim <- bind_rows(lapply(ks, function(i) {
      p <- c(ps[i, ])
      p <- p[names(p) %in% pars]
      p$Y0 = y0s[[i]]
      
      cm$set_user(user = p)
      ys <- cm$run(seq(2023, 2036, 1)) %>% 
        data.frame() %>% as_tibble() %>% 
        select(-starts_with("Y.")) %>% 
        rename(Year = t) %>% 
        mutate(Key = i)
    }))
    
    bind_rows(sim, ys.start) %>% 
      mutate(sc = sc) %>% 
      group_by(Key) %>% 
      mutate(
        IncR = c((IncR[-n()] + IncR[-1]) / 2, 0),
        MorR = c((MorR[-n()] + MorR[-1]) / 2, 0),
        Inc = c(Inc[-1], 0),
        Mor = c(Mor[-1], 0)
      ) %>% 
      filter(Year <= 2035) %>% 
      arrange(Key, Year)
  }
  
  
  si <- sim_intv(post0, ys0, y0s)
  
  
  res <- bind_rows(
    sim_intv(post0, ys0, y0s),
    sim_intv(post1, ys0, y0s),
    sim_intv(post2, ys0, y0s),
    sim_intv(post3, ys0, y0s)
  ) 
  
  stats <- res %>% 
    select(Year, IncR, MorR, Key, sc) %>% 
    pivot_longer(c(MorR, IncR), names_to = "Index") %>% 
    group_by(Year, sc, Index) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.25),
      U = quantile(value, 0.75)
    )
  
  
  avt <- res %>% 
    filter(sc != "baseline") %>% 
    left_join(res %>% filter(sc == "baseline") %>% select(Year, Key, Inc0 = Inc, Mor0 = Mor)) %>% 
    mutate(
      AInc = 1 - Inc / Inc0,
      AMor = 1 - Mor / Mor0,
      AInc = ifelse(is.na(AInc), 0, AInc),
      AMor = ifelse(is.na(AMor), 0, AMor)
    ) %>% 
    select(Year, AInc, AMor, Key, sc) %>% 
    pivot_longer(c(AInc, AMor), names_to = "Index") %>% 
    group_by(Year, sc, Index) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.25),
      U = quantile(value, 0.75)
    )
  
  
  
  
  ys_baseline <- ys0 %>% 
    group_by(Key) %>% 
    mutate(
      IncR = c((IncR[-n()] + IncR[-1]) / 2, 0),
      MorR = c((MorR[-n()] + MorR[-1]) / 2, 0),
      Inc = c(Inc[-1], 0),
      Mor = c(Mor[-1], 0)
    ) %>% 
    ungroup() %>% 
    select(Year, MorR, IncR) %>% 
    pivot_longer(-Year, names_to = "Index") %>% 
    group_by(Year, Index) %>% 
    filter(Year <= 2023) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.25),
      U = quantile(value, 0.75)
    )
  
  
  g_intv0 <- stats %>% 
    rename(Scenario = sc) %>% 
    filter(Year > 2022) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) + 
    scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
    facet_wrap(.~Index, scales = "free_y", 
               labeller = labeller(Index = c(IncR = "Incidence", MorR = "Mortality"))) +
    expand_limits(y = 0)
  
  g_intv <- g_intv0 + 
    geom_ribbon(data = ys_baseline, aes(x = Year, ymin = L, ymax = U), alpha = 0.1) +
    geom_line(data = ys_baseline, aes(x = Year, y = M)) + 
    geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U))
  
  
  g_intv
  g_intv0
  
  
  g_avt <- avt %>% 
    rename(Scenario = sc) %>% 
    ggplot() +
    geom_ribbon(aes(x = Year, ymin = L, ymax = U, fill = Scenario), alpha = 0.1) +
    geom_line(aes(x = Year, y = M, colour = Scenario)) + 
    scale_y_continuous("Case averted", labels = scales::percent_format()) + 
    facet_wrap(.~Index, scales = "free_y", 
               labeller = labeller(Index = c(AInc = "Incident cases", AMor = "TB-related deaths"))) +
    expand_limits(y = 0)
    
  
  g_avt
  
  ggsave(g_intv0, filename = "results/g_intv0_" + iso + ".png", width = 7, height = 4)
  ggsave(g_intv, filename = "results/g_intv_" + iso + ".png", width = 7, height = 4)
  ggsave(g_avt, filename = "results/g_avt_" + iso + ".png", width = 7, height = 4)
}




