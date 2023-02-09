
r_prior <- function(iso = "IND") {
  list(
    p_primary = runif(1, 0.09, 0.14),
    r_react = runif(1, 0.0005, 0.0015),
    r_relapse = runif(1, 0.0011, 0.002),
    r_clear = runif(1, 0.02, 0.04),
    beta = runif(1, 1, ifelse(iso == "IND", 20, 40)),
    p_im = runif(1, 0.69, 0.86),
    adr = runif(1, 0, ifelse(iso == "IND", 0.05, 0.2))
  )
}


get_pars_pop <- function(d_pop_all, year0) {
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
        N0 = y0,
        Year0 = year0
      )
    )
  })
}


calc_distance <- function(cm, p, d_burden) {
  cm$set_user(user = p)
  cm$run(seq(1700, 2035, 1)) %>% 
    data.frame() %>% as_tibble() %>% 
    select(-starts_with("Y.")) %>% 
    rename(Year = t) %>% 
    inner_join(d_burden, by = "Year") %>% 
    mutate(
      d_inc = (IncR - Inc_M) ^ 2 / (Inc_U - Inc_L) * 2 * 1.96,
      d_mor = (MorR - Mor_M) ^ 2 / (Mor_U - Mor_L) * 2 * 1.96,
      dist = d_inc + d_mor
    ) %>% 
    summarise(dist = sum(dist)) %>% 
    pull(dist)
}



