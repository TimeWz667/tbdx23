library(tidyverse)
library(tidybayes)
library(odin)

theme_set(theme_bw())


source(here::here("r", "fn_cas.R"))
source(here::here("r", "fn_fit.R"))


m <- odin::odin("odin/d_tb.R")


year0 <- 2010
p_eps <- 0.02
n_post <- 300

for (country in c("IND", "ZAF")) {
  iso <- glue::glue(country)
  
  
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
  
  pars_cs0 <- read_csv(here::here("results", "pars_" + iso + ".csv")) %>% 
    mutate(Country = "" + iso) %>% 
    calc_prev()
  
  pars_cs <- pars_cs0 %>% 
    calc_reform(pdx0 = 0, pdx1 = 1) %>% 
    select(r_sym, r_sc, r_death_a, r_death_s, r_death_tx, r_csi, r_recsi, pdx0, pdx1, p_under)
  
  
  r_pars <- function() {
    p <- r_prior(country)
    p <- c(p, pars_cs[sample(1:nrow(pars_cs), 1), ])
    p
  }
  
  
  # Initialise Model
  p <- r_prior(country)
  p <- c(p, pars_cs[sample(1:nrow(pars_cs), 1), ])
  cm <- m$new(user = c(pars_pop, p))
  
  
  # Initialise ABC
  dis <- sapply(1:200, function(i) {
    p <- r_pars()
    cm$set_user(user = p)
    calc_distance(cm, p, d_burden)
  })
  
  
  eps <- quantile(dis, p_eps)
  print(eps)
  
  # Gather posterior
  posterior <- lapply(1:n_post, function(key) {
    dist <- Inf
    while (dist > eps) {
      p <- r_pars()
      cm$set_user(user = p)
      dist <- calc_distance(cm, p, d_burden)
    }
    
    ys <- cm$run(seq(1700, 2035, 1)) %>% 
      data.frame() %>% as_tibble() %>% 
      select(-starts_with("Y.")) %>% 
      rename(Year = t) %>% 
      mutate(Key = key)
    
    list(
      Pars = p,
      Ys = ys
    )
  })
  
  
  pss <- bind_rows(lapply(posterior, function(x) x$Pars)) %>% 
    left_join(pars_cs0 %>% select(-adr))
  
  
  yss <- bind_rows(lapply(posterior, function(x) x$Ys))
  
  g_fit <- yss %>% 
    select(Year, Key, MorR, IncR) %>% 
    pivot_longer(-c(Year, Key), names_to = "Index") %>% 
    filter(Year >= year0) %>% 
    ggplot() +
    geom_line(aes(x = Year, y = value, group = Key), alpha = 0.1) +
    geom_pointrange(data = d2plot, aes(x = Year, y = M, ymin = L, ymax = U)) + 
    scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) + 
    facet_wrap(.~Index, scales = "free_y") +
    expand_limits(y = 0)
  
  write_csv(yss, "results/post_run_" + iso + ".csv")
  write_csv(pss, "results/post_pars_" + iso + ".csv")
  ggsave(g_fit, filename = "results/g_fit_" + iso + ".png", width = 7, height = 4)
}
