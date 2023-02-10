library(tidyverse)


source(here::here("r", "fn_fit.R"))


load(here::here("data", "ZAF", "d_pop.rdata"))
load(here::here("data", "ZAF", "d_hiv.rdata"))


year0 <- 2010

pars_pop <- get_pars_pop(d_pop_all, year0)



dd <- d_hiv %>% 
  mutate(
    U = N_Pop - PrvHIV,
    H = PrvHIV * (1 - PrART_HIV),
    A = PrvHIV * PrART_HIV,
  ) %>% 
  select(-Country) %>% 
  mutate(
    mu_h = coef(lm(MorHIV ~ H - 1)),
    mu = (R_Die * N_Pop - mu_h * H) / N_Pop
  )

d_hiv %>% 
  mutate(
    U = N_Pop - PrvHIV,
    H = PrvHIV * (1 - PrART_HIV),
    A = PrvHIV * PrART_HIV,
    MorHIV / H
  )


d1 <- dd %>% select(Year, U:mu) %>% 
  mutate(
    dA = (A - lag(A)) / 2,
    A = (A + lag(A)) / 2,
    dH = (H - lag(H)) / 2,
    H = (H + lag(H)) / 2,
    dU = (U - lag(U)) / 2,
    U = (U + lag(U)) / 2,
    mu = (mu + lag(mu)) / 2,
    mu_h = mu_h[1],
    r_art = (dA + mu * A) / H,
    r_hiv = (dH + (r_art + mu + mu_h) * H) / U, 
    dt = Year - 2010,
    kk = log(0.045 / (r_art - 0.07) - 1)
  )
d1 <- d1[-1, ]


coh <- lm(log(r_hiv) ~ dt, data = d1) %>% coef() %>% unname()
coa = lm(kk ~ dt, data = d1) %>% coef() %>% unname()

pars_hiv <- list(
  r_hiv0 = exp(coh[1]),
  rt_hiv = coh[2],
  r_die_hiv = dd$mu_h[1] %>% unname(),
  t0_art = 2010 - coa[1] / coa[2],
  rt_art = coa[2]
)

pars_hiv


pars <- c(pars_pop, pars_hiv)

names(pars)

y0 <- d_hiv %>% 
  mutate(
    U = N_Pop - PrvHIV,
    H = PrvHIV * (1 - PrART_HIV),
    A = PrvHIV * PrART_HIV
  )

pars$Y0 <- unlist(y0[1, c("U", "H", "A")])

model <- odin::odin("odin/d_hiv.R")

cm <- model$new(user = pars)


cm$set_user(
  r_hiv0 = x[1],
  rt_hiv = - x[2],
  rt_art = - x[3],
  r_die_hiv = x[4]
)

ys <- cm$run(seq(2010, 2019, 1)) %>% data.frame() %>% tibble()

sims <- ys %>% 
  select(Year = t, PrvHIV, PrART_HIV=PrART, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index")


d2plot <- d_hiv %>% 
  mutate(
    PrvHIV = PrvHIV / N_Pop,
    MorHIV = MorHIV / N_Pop
  ) %>% 
  select(Year, PrvHIV, PrART_HIV, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index")


sims %>% 
  ggplot() +
  geom_line(aes(x = Year, y = value)) +
  geom_point(data = d2plot, aes(x = Year, y = value)) +
  facet_wrap(.~Index, scales = "free_y")






d2fit <- d_hiv %>% 
  mutate(
    PrvHIV = PrvHIV / N_Pop,
    MorHIV = MorHIV / N_Pop
  ) %>% 
  select(Year, PrvHIV, PrART_HIV, MorHIV)


# x.start <- unlist(pars_hiv[c("r_hiv0", "rt_hiv", "rt_art", "r_die_hiv")])
# 
# 
# opt <- nlminb(x.start, function(x) {
#   cm$set_user(
#     r_hiv0 = x[1],
#     rt_hiv = - x[2],
#     rt_art = - x[3],
#     r_die_hiv = x[4]
#   )
#   
#   ys <- cm$run(seq(2010, 2019, 1)) %>% data.frame() %>% tibble()
#   
#   dis <- (d2fit$PrvHIV / ys$PrvHIV - 1) ^ 2
#   dis <- dis + (d2fit$PrART_HIV / ys$PrART - 1) ^ 2
#   dis <- dis + (d2fit$MorHIV / ys$MorHIV - 1) ^ 2
#   sum(dis)
# }, lower = 0, upper = 2)
# 
# x = opt$par
# 
# opt




sims <- ys %>% 
  select(Year = t, PrvHIV, PrART_HIV=PrART, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index")


d2plot <- d_hiv %>% 
  mutate(
    PrvHIV = PrvHIV / N_Pop,
    MorHIV = MorHIV / N_Pop
  ) %>% 
  select(Year, PrvHIV, PrART_HIV, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index")


sims %>% 
  ggplot() +
  geom_line(aes(x = Year, y = value)) +
  geom_point(data = d2plot, aes(x = Year, y = value)) +
  facet_wrap(.~Index, scales = "free_y")



sims <- ys %>% 
  select(Year = t, r_art=R_ART, r_hiv=R_HIV) %>% 
  pivot_longer(-Year, names_to = "Index")


d2plot <- d1 %>% 
  select(Year, r_art, r_hiv) %>% 
  pivot_longer(-Year, names_to = "Index")

sims %>% 
  ggplot() +
  geom_line(aes(x = Year, y = value)) +
  geom_point(data = d2plot, aes(x = Year, y = value)) +
  facet_wrap(.~Index, scales = "free_y")








