## ODEs -----
deriv(Y[]) <- d_pop[i] + d_hiv[i]
dim(Y) <- n_hiv


## Initial values -----
initial(Y[]) <- Y0[i]



Y0[] <- user()
dim(Y0) <- n_hiv

Year0 <- user(2010)



## Output -----
output(N) <- n

output(PrvHIV) <- (H + A) / n
output(PrART) <- A / (H + A)
output(MorHIV) <- mor / n

output(R_HIV) <- r_hiv
output(R_ART) <- r_art

## Summary -----
n <- sum(Y)


## dims -----
years[] <- user() # data times for AIM by 1
dim(years) <- user()


## lengths -----
n_hiv <- 3 # U, H, A

n_years <- length(years)


## demography -----
br[] <- user()
dim(br) <- n_years

br_t <- interpolate(years, br, "linear")

dr[] <- user()
dim(dr) <- n_years

dr_t <- interpolate(years, dr, "linear")

r_die_hiv <- user(0.03)

mu <- (dr_t * n - r_die_hiv * Y[2]) / n

mor <- r_die_hiv * Y[2]

d_pop[1] <- br_t * n - mu * U
d_pop[2] <- - (r_die_hiv + mu) * H
d_pop[3] <- - mu * A
dim(d_pop) <- n_hiv



## HIV
r_hiv0 <- user()
rt_hiv <- user()

t0_art <- user()
rt_art <- user()


r_hiv <- r_hiv0 * exp(rt_hiv * (t - Year0))
r_art <- 0.07 + 0.045 / (1 + exp(rt_art * (t - t0_art)))


U <- Y[1]
H <- Y[2]
A <- Y[3]


d_hiv[1] <- - r_hiv * U
d_hiv[2] <- r_hiv * U - r_art * H
d_hiv[3] <- r_art * H
dim(d_hiv) <- 3

