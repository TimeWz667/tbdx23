library(tidyverse)


source(here::here("r", "fn_fit.R"))


load(here::here("data", "ZAF", "d_pop.rdata"))
load(here::here("data", "ZAF", "d_hiv.rdata"))


year0 <- 2010


pars_pop <- get_pars_pop(d_pop_all, year0)



dd <- tibble(
  years = seq(2009.5, 2019.5, 0.01)
) %>% 
  mutate(
    N_Pop = spline(d_hiv$Year + 0.5, d_hiv$N_Pop, xout = years, method = "fmm")$y,
    PrvHIV = spline(d_hiv$Year + 0.5, d_hiv$PrvHIV / d_hiv$N_Pop, xout = years, method = "fmm")$y,
    PrART = spline(d_hiv$Year + 0.5, d_hiv$PrART_HIV, xout = years, method = "fmm")$y,
    mu_h = spline(d_hiv$Year + 0.5, d_hiv$MorHIV / d_hiv$PrvHIV, xout = years, method = "fmm")$y,
    U = N_Pop * (1 - PrvHIV),
    H = N_Pop * PrvHIV * (1 - PrART),
    A = N_Pop * PrvHIV * PrART,
    dr = spline(pars_pop$years, pars_pop$dr, xout = years, method = "fmm")$y,
    br = spline(pars_pop$years, pars_pop$br, xout = years, method = "fmm")$y
  )


hiv0 <- dd[1, ] %>% select(U, H, A) %>% unlist()
hiv0


rates <- dd %>% 
  select(years, ue=U, he=H, ae=A, mu_h, dr, br)
  
  

p <- as.list(rates)
p$Y0 <- hiv0
p$Year0 <- 2010


model <- odin::odin("odin/d_hiv_pre.R")

cm <- model$new(user = p)


ys <- cm$run(seq(2009.5, 2019.5, 0.1)) %>% data.frame() %>% tibble()

ys %>% 
  ggplot() +
  geom_line(aes(x = t, y = R_ART)) +
  geom_line(aes(x = t, y = R_HIV))



r_art <- function(x, ys) {
  r_art0 <- x[1]
  r_art1 <- x[1] + x[2]
  rt_art <- x[3]
  t0_art <- x[4]
  
  years <- ys$t
  
  r_art0 + (r_art1 - r_art0) / (1 + exp(- rt_art * (years - t0_art)))
}


opt_art <- nlminb(c(0.13, 0.07, 3, 2015.5), objective = function(x) {
  sum((r_art(x, ys) - ys$R_ART) ^ 2)
}, lower = c(0.1, 0, 1, 2010), upper = c(0.2, 0.1, 5, 2020))


plot(ys$t, ys$R_ART, type = 'l')
points(ys$t, r_art(opt_art$par, ys))


r_hiv <- function(x, ys) {
  r_hiv0 <- x[1]
  rt_hiv <- x[2]
  
  years <- ys$t
  
  r_hiv0 * exp(- rt_hiv * (years - 2010))
}


opt_hiv <- nlminb(c(0.01, 0.08), objective = function(x) {
  sum((r_hiv(x, ys) - ys$R_HIV) ^ 2)
}, lower = c(0.005, 0.01), upper = c(0.015, 0.05))


plot(ys$t, ys$R_HIV, type = 'l')
points(ys$t, r_hiv(opt_hiv$par, ys))


dr_hiv <- function(x, ys) {
  dr_hiv0 <- x[1]
  drt_hiv <- x[2]
  dr_hiv1 <- x[3]
  years <- ys$t
  
  dr_hiv0 * exp(- drt_hiv * (years - 2010)) + dr_hiv1
}


opt_drhiv <- nlminb(c(0.01, 0.08, 0.0004), objective = function(x) {
  sum((dr_hiv(x, ys) - ys$mu_h_t) ^ 2)
}, lower = c(0.02, 0.01, 0), upper = c(0.1, 1, 0.01))


plot(ys$t, ys$mu_h_t, type = 'l')
points(ys$t, dr_hiv(opt_drhiv$par, ys))



pars_hiv <- list(
  r_art0 = opt_art$par[1],
  r_art1 = opt_art$par[1] + opt_art$par[2],
  rt_art = opt_art$par[3],
  t0_art = opt_art$par[4],
  r_hiv0 = opt_hiv$par[1],
  rt_hiv = opt_hiv$par[2],
  dr_hiv0 = opt_drhiv$par[1],
  drt_hiv = opt_drhiv$par[2],
  dr_hiv1 = opt_drhiv$par[3]
)

pars_hiv


pars <- c(pars_pop, pars_hiv)

pars$Y0 <- hiv0


model <- odin::odin("odin/d_hiv.R")

cm <- model$new(user = pars)



d2fit <- d_hiv %>% 
  filter(Year >= 2010) %>% 
  mutate(
    t = Year + 0.5,
    PrvHIV = PrvHIV / N_Pop,
    PrART_HIV = PrART_HIV,
    MorHIV = MorHIV / N_Pop
  ) %>% 
  select(Year = t, PrvHIV, PrART_HIV, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index", values_to = "Data")
  

x.start <- c(opt_art$par, opt_hiv$par) #, opt_drhiv$par)
x.lower <- c(c(0.1, 0, 1, 2010),       c(0.005, 0.01)) #, c(0.02, 0.01, 0))
x.upper <- c(c(0.5, 0.5, 5, 2020),     c(0.015, 0.05)) #, c(0.1, 1, 0.01))


opt <- nlminb(x.start, objective = function(x) {
  cm$set_user(
    r_art0 = x[1],
    r_art1 = x[1] + x[2],
    rt_art = x[3],
    t0_art = x[4],
    r_hiv0 = x[5],
    rt_hiv = x[6]
  )
  
  ys <- cm$run(seq(2009.5, 2019.5, 0.5)) %>% data.frame() %>% tibble()
  d <- ys %>% 
    select(Year = t, PrvHIV, PrART_HIV = PrART) %>% 
    pivot_longer(-Year, names_to = "Index") %>% 
    inner_join(d2fit, c("Year", "Index")) %>% 
    summarise(dsq = sum((Data - value) ^ 2)) %>% pull(dsq)
  d
}, lower = x.lower, upper = x.upper)



x <- opt$par
#x <- x.start

cm$set_user(
  r_art0 = x[1],
  r_art1 = x[1] + x[2],
  rt_art = x[3],
  t0_art = x[4],
  r_hiv0 = x[5],
  rt_hiv = x[6]
)

ys <- cm$run(seq(2009.5, 2036, 0.5)) %>% data.frame() %>% tibble()

ys %>% 
  select(Year = t, PrvHIV, PrART_HIV=PrART, MorHIV) %>% 
  pivot_longer(-Year, names_to = "Index") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = value)) +
  geom_point(data = d2fit, aes(x = Year, y = Data)) + 
  facet_wrap(.~Index, scales = "free_y")


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




