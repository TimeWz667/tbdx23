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
output(mu_h_t) <- TRUE

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

mu_h[] <- user()
dim(mu_h) <- n_years
mu_h_t <- interpolate(years, mu_h, "linear")


mu <- (dr_t * n - mu_h_t * Y[2]) / n

mor <- mu_h_t * Y[2]

d_pop[1] <- br_t * n - mu * U
d_pop[2] <- - (mu_h_t + mu) * H
d_pop[3] <- - mu * A
dim(d_pop) <- n_hiv



## HIV
ue[] <- user()
dim(ue) <- n_years
ue_t <- interpolate(years, ue, "linear")

he[] <- user()
dim(he) <- n_years
he_t <- interpolate(years, he, "linear")

ae[] <- user()
dim(ae) <- n_years
ae_t <- interpolate(years, ae, "linear")


U <- Y[1]
H <- Y[2]
A <- Y[3]

r_art <- 10 * (ae_t - A) / A #(ae_t - A - d_pop[3]) / H
r_hiv <- 10 * (he_t - H) / H #(he_t - H - d_pop[2] + r_art * H) / U



d_hiv[1] <- - r_hiv * U
d_hiv[2] <- r_hiv * U - r_art * H
d_hiv[3] <- r_art * H
dim(d_hiv) <- 3
