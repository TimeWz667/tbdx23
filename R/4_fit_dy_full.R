

### Load the model file
f <- system.file("models/SIR.txt", package = "odin2data")
test_m <- odin::odin(f, verbose=F)


### Generate a prior and set up it as a list
r_prior <- function() {
  list(beta = runif(1, 0.1, 5), gamma = runif(1, 0.1, 0.3))
}

d_prior <- function(pars) {
  dunif(pars$beta, 0.1, 5, log = T) + dunif(pars$gamma, 0.1, 0.3, log = T)
}

times = seq(0, 10, 0.2)
y0 <- c(995, 5, 0)

### Compile all elements as a simulation model
sim <- odin2data::compile_model(d_prior, r_prior, y0, ts_sim = times, m_sim = test_m)


test_data <- data.frame(
  t = 1:5,
  incidence = c(20, 49, 109, 184, 206) / 1000
)

### Compile the model with data
lf <- odin2data::compile_model_likefree(test_data, sim)



post <- odin2data::fit(lf, 1000, method = "abcsmc", max_round = 20, alpha = 0.80)


summary(post)


