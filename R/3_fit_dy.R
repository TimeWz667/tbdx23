library(tidyverse)
library(odin)


source(here::here("r", "fn_cas.R"))




year0 <- 2010

load("data/IND/d_pop.rdata")


pars_pop <- local({
  pars_pop <-  tibble(
    years = move_avg2(d_pop_all$Year),
    gr = diff(log(d_pop_all$N_Pop)),
    dr = move_avg2(d_pop_all$N_Pop * d_pop_all$R_Die) / move_avg2(d_pop_all$N_Pop)
  ) %>% 
    mutate(
      br = gr + dr
    ) %>% 
    select(-gr) %>% 
    filter(years < 2050) %>% 
    as.list()
  
  pars_pop$years <- c(0, pars_pop$years)
  pars_pop$dr <- c(pars_pop$dr[1], pars_pop$dr)
  pars_pop$br <- c(pars_pop$br[1], pars_pop$br)
  
  y0 <- d_pop_all %>% filter(Year == year0) %>% pull(N_Pop)
  
  pars_pop <- c(
    pars_pop,
    list(
      Y0 = c(0.95, 0.01, 0, 0, 0, 0.03, 0.01, 0) * y0,
      Year0 = year0
    )
  )
})



pars_cs <- read_csv(here::here("results", "pars_IND.csv")) %>% 
  mutate(Country = "IND") %>% 
  calc_prev() %>% 
  calc_reform(pdx0 = 0, pdx1 = 1) %>% 
  select(r_sym, r_sc, r_death_a, r_death_s, r_csi, r_recsi, pdx0, pdx1, p_under)






m <- odin::odin("odin/d_tb.R")

p_cs <- c(pars_cs[1, ])

p <- c(pars_pop, p_cs)
p$beta <- 15
p$adr <- 0.02

cm <- m$new(user = p)




ys <- cm$run(seq(1700, 2035, 1)) %>% 
  data.frame() %>% as_tibble() %>% 
  select(-starts_with("Y.")) %>% 
  filter(t > 2010) %>% 
  mutate(
    Mor = Mor - Mor[1],
    Inc = Inc - Inc[1]
  ) %>% 
  rename(Year = t)


ys %>% 
  select(Year, MorR, IncR) %>% 
  pivot_longer(-Year, names_to = "Index") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = value)) +
  facet_wrap(.~Index, scales = "free_y") + 
  scale_y_continuous("", labels = scales::number_format(scale = 1e5)) + 
  expand_limits(y = 0)


d_pop_all %>% 
  ggplot() +
  geom_line(aes(x = Year, y = N_Pop)) +
  geom_point(data = ys, aes(x = Year, y = N)) +
  scale_y_continuous("Population, B", labels = scales::number_format(scale = 1e-9)) + 
  expand_limits(y = 0)



load("data/IND/d_burden.rdata")



ys %>% inner_join(d_burden) %>% 
  mutate(
    d_inc = (IncR - Inc_M) ^ 2 / (Inc_U - Inc_L) * 2 * 1.96,
    d_mor = (MorR - Mor_M) ^ 2 / (Mor_U - Mor_L) * 2 * 1.96,
    dist = d_inc + d_mor
  ) %>% 
  summarise(dist = sum(dist))



