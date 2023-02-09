library(tidyverse)
library(odin)

theme_set(theme_bw())


source(here::here("r", "fn_cas.R"))
source(here::here("r", "fn_fit.R"))


year0 <- 2010


m <- odin::odin("odin/d_cas.R")


country = "IND"
iso = glue::as_glue(country)




# Load data
load("data/" + iso + "/d_burden.rdata")
load("data/" + iso + "/d_pop.rdata")
load("data/" + iso + "/d_cases.rdata")
load("data/" + iso + "/d_prev.rdata")



d_prev %>% filter(Tag == "All") %>% 
  left_join(d_pop_all) %>% 
  left_join(d_case_all)





d_burden

year0 <- 2020
adr <- - mean(diff(log(d_burden$Inc_M)))
k <- exp(- adr *(d_burden$Year - year0))
inc0 <- mean(d_burden$Inc_M / k)

d <- tibble(
  Year = 2014:2022,
  year0 = 2020,
  inc0 = inc0,
  adr = adr
) %>% 
  mutate(
    k = exp(- adr * (Year - year0)),
    Inc = k * inc0
  ) %>% 
  bind_cols(
    d_prev %>% 
      mutate(
        prev = N_Prev / N_Subject * 1.13
      ) %>% 
      filter(Tag == "All") %>% 
      ungroup() %>% 
      select(State, prev) %>% 
      pivot_wider(names_from = State, values_from = prev)
  ) %>% 
  left_join(d_case_all %>% 
              left_join(d_pop_all) %>% 
              mutate(CNR = N_Case / N_Pop) %>% 
              select(Year, CNR))


d2 <- d %>% 
  mutate(
    ra = 0.12 * 0.4 + 0.2,
    rs = 0.12 + 0.2,
    rc = 0.12 + 0.2,
    r_sym = Inc / Asym - ra + adr,
    r_aware = r_sym * Asym / Sym - rs + adr,
    r_det = r_aware * Sym / CS - rc + adr,
    Det = r_det * CS
  ) %>% 
  left_join(d_case_all) %>% 
  mutate(
    r_sym = r_sym[Year == year0],
    r_aware = r_aware[Year == year0],
    r_det = r_det[Year == year0],
    Asym = Asym * k,
    Sym = Sym * k
  )



fn <- function(x, d2) {
  rr_det_t = x[1]
  al = x[2]
  
  y = d2 %>% 
    mutate(
      r_aware_t = r_aware * rr_det_t ^ (Year - year0),
      r_det_t = r_det * rr_det_t ^ (Year - year0),
      Sym = r_sym * Asym / (r_aware_t + rs - adr),
      CS = r_aware_t * Sym / (r_det_t + rc - adr),
      CNR_hat = r_det_t * CS * al
    )
  
  sum((y$CNR / y$CNR_hat - 1) ^ 2, na.rm = T)
}

x <- c(2, 1)

opt <- nlminb(x, fn, lower = 0.1, upper = 8, d2=d2)

opt

x = opt$par

d2 %>% 
  mutate(
    r_det0 = opt$par[1],
    rr_det_t = opt$par[2],
    al = opt$par[3],
    r_det = r_det0 * rr_det_t ^ (Year - year0),
    CS = r_aware * Sym / (r_det + rc - adr),
    CNR_hat = r_det * CS * al
  ) %>% 
  select(Year, Asym, Sym, CS, CNR, CNR_hat, r_det)



r_aware * Sym = r_det_t * CS + (rc - adr) * CS

d_case_all %>% 
  left_join(d_pop_all) %>% 
  mutate(CNR = N_Case / N_Pop) %>% 
  select(Year, CNR)


pars_cs <- read_csv(here::here("results", "pars_" + iso + ".csv")) %>% 
  mutate(Country = "" + iso) %>% 
  calc_prev() %>% 
  calc_reform(pdx0 = 0, pdx1 = 1)


dd = d_case_all %>% 
  left_join(d_pop_all) %>% 
  select(Year, N_Case, N_Pop) %>% 
  left_join(d_burden) %>% 
  mutate(
    CNR = N_Case / N_Pop
  ) %>% 
  select(Year, IncR = Inc_M, CNR) %>% 
  mutate(
    r_onset = 0.924,
    r_det = CNR / IncR
  )


library(TSA)

ts.plot(dd$r_det)


pars_cs %>% 
  select(starts_with("r_")) %>% 
  summarise(across(everything(), mean))


r_pars <- function() {
  p <- r_prior(country)
  p <- c(p, pars_cs[sample(1:nrow(pars_cs), 1), 
                    c('r_sym', 'r_sc', 'r_death_a', 'r_death_s', 'r_death_tx', 
                      'r_csi', 'r_recsi', 'pdx0', 'pdx1', 'p_under')])
  p
}

pcs <- pars_cs[sample(1:nrow(pars_cs), 1), ]



p <- c(pcs)[c('r_sym', 'r_sc', 'r_death_a', 'r_death_s', 'r_death_tx', 
         'r_csi', 'r_recsi', 'pdx0', 'pdx1', 'p_under')]

p$Y0 <- c(
  (1 - pcs$prv),
   pcs$prv_a,
  pcs$prv_s,
  pcs$prv_c,
  0
) * pars_pop$N0


# Initialise Model

cm <- m$new(user = c(pars_pop, p))

ys <- cm$run(seq(2000, 2030, 1))

ys %>% data.frame() %>% tibble()





