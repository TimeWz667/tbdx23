library(tidyverse)
library(odin)

theme_set(theme_bw())


source(here::here("r", "fn_cas.R"))
source(here::here("r", "fn_fit.R"))


year0 <- 2010



country = "IND"
iso = glue::as_glue(country)




# Load data
load("data/" + iso + "/d_burden.rdata")
load("data/" + iso + "/d_pop.rdata")
load("data/" + iso + "/d_cases.rdata")
load("data/" + iso + "/d_prev.rdata")


cnr_q <- read_csv(here::here("data", "pars", iso, "targets_q.csv"))


cnr <- bind_rows(
  d_case_all %>% 
    left_join(d_pop_all) %>% 
    mutate(CNR = N_Case / N_Pop) %>% 
    select(Year, CNR),
  cnr_q %>% 
    group_by(Year) %>% 
    summarise(CNR = mean(QCNR * 4))
) %>% 
  mutate(
    Time = Year + 0.5,
    COVID = (Year >= 2020 & Year < 2022) + 0,
    dt = Time - 2010
  )

cnr


cnrq <- bind_rows(
  d_case_all %>% 
    left_join(d_pop_all) %>% 
    mutate(CNR = N_Case / N_Pop, Time = Year + 0.5) %>% 
    select(Year, Time, CNR),
  cnr_q %>% 
    mutate(CNR = QCNR * 4)  %>% 
    select(Year, Time, CNR)
)

cnrq %>% 
  ggplot() +
  geom_line(aes(x = Time, y = CNR))


res <- lm(log(CNR)~ dt + COVID, data = cnr %>% mutate(dt = Year - 2010))
summary(res)


cnr <- cnr %>% 
  mutate(
    CNR_mu = exp(fitted(res))
  )


cnr %>% 
  ggplot() +
  geom_line(aes(x = Time, y = CNR_mu)) +
  geom_point(aes(x = Time, y = CNR)) +
  expand_limits(y = 0)





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
  year0 = year0,
  inc0 = inc0 * 1,
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
  left_join(cnr)


 d %>% 
  mutate(
    ra = 0.12 * 1 + 0.2,
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
    Asym = Asym * k
  ) %>% 
  mutate(
    Sym = r_sym * Asym / (r_aware + rs - adr),
    CS = r_aware * Sym / (r_det + rc - adr),
    Det_hat = r_det * CS,
    gap = CNR_mu * 0.5 / Det,
    p_report = 0.2 + 0.8 / (1 + exp(- 0.22 * (Year - 2020)))
  ) %>% 
  select(Year, Det_hat, Det, CNR_mu, gap, p_report) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = gap)) +
  geom_point(aes(x = Year, y = p_report))



fn <- function(x, d2, year_ref=2010) {
  y <- d2 %>% 
    mutate(
      year_ref = year_ref,
      k_det = x[1],
      rr_det_t = x[2],
      al = 1.4,
      k = k_det * exp(rr_det_t * (Year - year_ref) ^ 1),
      r_det_t = r_det * k,
      r_aware_t = r_aware * k,
      Sym = r_sym * Asym / (r_aware_t + rs - adr),
      CS = r_aware_t * Sym / (r_det_t + rc - adr),
      CNR_hat = r_det_t * CS * al
    )
  
  sum((y$CNR / y$CNR_hat - 1) ^ 2, na.rm = T)
}

opt <- nlminb(c(2, 1), fn, lower = 0, upper = 4, d2 = d2)
opt

d2 %>% 
  mutate(
    year_ref = 2010,
    k_det = opt$par[1],
    rr_det_t = opt$par[2],
    al = 1.4,
    k = k_det * exp(rr_det_t * (Year - year_ref) ^ 1),
    r_det_t = r_det * k,
    r_aware_t = r_aware * k,
    Sym = r_sym * Asym / (r_aware_t + rs - adr),
    CS = r_aware_t * Sym / (r_det_t + rc - adr),
    CNR_hat = r_det_t * CS * al
  ) %>% 
  select(Year, Asym, Sym, CS, CNR, CNR_hat, k, r_aware_t, r_det_t) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = CNR_hat)) +
  geom_point(aes(x = Year, y = CNR)) +
  expand_limits(y = 0)


