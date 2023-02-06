## ODEs -----
deriv(Y[]) <- d_pop[i] + d_tb[i] - adj * Y[i]
dim(Y) <- n_tb


deriv(Inc) <- inc
deriv(Mor) <- mor

adj <- if (t > Year0) 0 else sum(d_pop) / sum(Y)



## Initial values -----
initial(Y[]) <- Y0[i]

initial(Inc) <- 0
initial(Mor) <- 0

Y0[] <- user()
dim(Y0) <- n_tb

Year0 <- user(2010)



## Output -----
output(N) <- n
output(MorR) <- mor / n
output(IncR) <- inc / n


## Summary -----
n <- sum(Y)


## dims -----
years[] <- user() # data times for AIM by 1
dim(years) <- user()


## lengths -----
n_tb <- 8 # U, A, S, C, Tx, FL, SL, R

n_years <- length(years)


## demography -----
br[] <- user()
dim(br) <- n_years

br_t <- interpolate(years, br, "linear")

dr[] <- user()
dim(dr) <- n_years

dr_t <- interpolate(years, dr, "linear")


r_death_a <- user() 
r_death_s <- user() 

mu_t[] <- (dr_t * sum(Y) - r_death_a * Y[2] - r_death_s * (Y[3] + Y[4])) / sum(Y)
dim(mu_t) <- n_tb 


mor <- - (d_pop[2] + d_pop[3] + d_pop[4])

d_pop[1] <- br_t * sum(Y) - mu_t[1] * Y[1]
d_pop[2] <- - (r_death_a + mu_t[i]) * Y[i]
d_pop[3] <- - (r_death_s + mu_t[i]) * Y[i]
d_pop[4] <- - (r_death_s + mu_t[i]) * Y[i]
d_pop[5:n_tb] <- - mu_t[i] * Y[i]
dim(d_pop) <- c(n_tb)



## TB
r_sym <- user() 
r_sc <- user() 

r_csi <- user() 
r_recsi <- user() 
pdx0 <- user() 
pdx1 <- user()
p_under <- user()

dur_tx <- 0.5


r_act <- r_lat * p_primary / (1 - p_primary)
p_primary <- user(0.1)
r_lat <- user(0.5)

r_react <- user(0.0001)
r_relapse <- user(0.0001)



beta <- user(10)
p_im <- user(0.3)
adr <- user(0.01)

# U, A, S, C, Tx, FL, SL, R
U <- Y[1]
Asym <- Y[2]
Sym <- Y[3]
ExCS <- Y[4]
Tx <- Y[5]

FL <- Y[6]
SL <- Y[7]
Rec <- Y[8]


det0 <- r_csi * pdx0 * Sym
fn0 <- r_csi * (1 - pdx0) * Sym
det1 <- r_recsi * pdx1 * ExCS
det <- det0 + det1

# inc <- 2e-3 * n
inc <- r_act * FL + r_react * SL + r_relapse * Rec

k_beta <- if (t > Year0) exp(- adr * sqrt(t - Year0)) else 1
foi <- beta * k_beta * (Asym + Sym + ExCS) / n
rfoi <- foi * (1 - p_im)

d_tb[1] <- - foi * U
# A, S, C, Tx
d_tb[2] <- inc - (r_sym + r_sc) * Asym
d_tb[3] <- r_sym * Asym - (r_csi + r_sc) * Sym
d_tb[4] <- fn0 - (r_recsi * pdx1 + r_sc) * ExCS
d_tb[5] <- det - 1 / dur_tx * Tx
# FL, SL, R
d_tb[6] <- foi * U + rfoi * (SL + Rec) - (r_act + r_lat) * FL
d_tb[7] <- r_lat * FL + r_sc * (Asym + Sym + ExCS) - (rfoi + r_react) * SL
d_tb[8] <- 1 / dur_tx * Tx - (rfoi + r_relapse) * Rec
dim(d_tb) <- n_tb
