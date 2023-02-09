## ODEs -----
deriv(Y[]) <- d_pop[i] + d_tb[i] - adj * Y[i]
dim(Y) <- n_tb


deriv(Inc) <- inc
deriv(Mor) <- mor_tb

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
output(MorR) <- mor_tb / n
output(IncR) <- inc / n


## Summary -----
n <- sum(Y)


## dims -----
years[] <- user() # data times for AIM by 1
dim(years) <- user()


## lengths -----
n_tb <- 5 # U, A, S, C, Tx

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
r_death_tx <- user()

mor_tb <- - r_death_a * Y[2] + r_death_s * (Y[3] + Y[4]) + r_death_tx * Y[5]

d_pop[1] <- br_t * sum(Y) - dr_t * Y[1] + mor_tb
d_pop[2] <- - (r_death_a + dr_t) * Y[i]
d_pop[3] <- - (r_death_s + dr_t) * Y[i]
d_pop[4] <- - (r_death_s + dr_t) * Y[i]
d_pop[5] <- - (r_death_tx + dr_t) * Y[i]
dim(d_pop) <- n_tb



## TB
r_sym <- user() 
r_sc <- user() 

r_csi <- user() 
r_recsi <- user() 
pdx0 <- user() 
pdx1 <- user()
p_under <- user()


dur_tx <- 0.5

r_tx_succ <- 1 / dur_tx

adr <- user(0.01)

# U, A, S, C, Tx, FL, SL, R
U <- Y[1]
Asym <- Y[2]
Sym <- Y[3]
ExCS <- Y[4]
Tx <- Y[5]


det0 <- r_csi * pdx0 * Sym
fn0 <- r_csi * (1 - pdx0) * Sym
det1 <- r_recsi * pdx1 * ExCS
det <- det0 + det1

# inc <- 2e-3 * n

adr_t <- if (t > Year0) adr else 0
inc <- (r_sym + r_death_a + dr_t - adr_t) * Asym 

d_tb[1] <- - inc + r_sc * (Asym + Sym + ExCS) + r_tx_succ * Tx
# A, S, C, Tx
d_tb[2] <- inc - (r_sym + r_sc) * Asym
d_tb[3] <- r_sym * Asym - r_sc * Sym - det0 - fn0
d_tb[4] <- fn0 - (r_recsi * pdx1 + r_sc) * ExCS
d_tb[5] <- det - r_tx_succ * Tx
# FL, SL, R
dim(d_tb) <- n_tb
