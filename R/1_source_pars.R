library(tidyverse)


if (!file.exists("pars")) {
  file.create("pars", showWarnings = F)
  file.create("pars/IND", showWarnings = F)
  file.create("pars/ZAF", showWarnings = F)
}



# Demographics -----

iso = "IND"

load(here::here("data", iso, "d_pop.rdata"))


move_avg2 <- function(x) {
  (x[-1] + x[-length(x)]) / 2
}


for (iso in c("IND", "ZAF")) {
  pars_pop <- local({
    pars_pop <-  tibble(
      years = move_avg2(d_pop_all$Year),
      N = move_avg2(d_pop_all$N_Pop),
      gr = diff(log(d_pop_all$N_Pop)),
      dr = move_avg2(d_pop_all$N_Pop * d_pop_all$R_Die) / N
    ) %>% 
      mutate(
        br = gr + dr
      ) %>% 
      select(-gr) %>% 
      filter(years <= 2050) %>% 
      as.list()
    
    pars_pop$years <- c(0, pars_pop$years)
    pars_pop$dr <- c(pars_pop$dr[1], pars_pop$dr)
    pars_pop$br <- c(pars_pop$br[1], pars_pop$br)
    
    pars_pop$N <- c(pars_pop$N[1], pars_pop$N)
    
    pars_pop
  })
  
  
  jsonlite::write_json(pars_pop, here::here("pars", iso, "pars_pop.json"), digits = 10, auto_unbox = T)
}


# Treatment outcomes
iso <- "IND"
load(here::here("data", iso, "d_tx_itr.rdata"))

d_tx_itr <- as.list(d_tx_itr)


ptx <- list(
  All = list(
    Succ = d_tx_itr$p_succ_all,
    Die = d_tx_itr$p_die_all
  ),
  Pub = list(
    Succ = d_tx_itr$p_succ_pub,
    Death = d_tx_itr$p_die_pub
  ),
  Pri = list(
    Succ = d_tx_itr$p_succ_pri,
    Death = d_tx_itr$p_die_pri
  )
)

jsonlite::write_json(ptx, here::here("pars", iso, "pars_tx.json"), digits = 10, auto_unbox = T)



iso <- "ZAF"
load(here::here("data", iso, "d_tx.rdata"))

ptx <- d_tx %>% 
  group_by(Country, Tag, Stat) %>% 
  summarise(across(starts_with("value_"), weighted.mean, w = N)) %>% 
  extract(Stat, "Type", "TxOut(\\S+)")

tags <- unique(ptx$Tag)
tags <- set_names(tags, tags)

ptx <- lapply(tags, function(tag) {
  p <- ptx %>% filter(Tag == tag)
  
  as.list(set_names(p$value_M, p$Type))
})


jsonlite::write_json(ptx, here::here("pars", iso, "pars_tx.json"), digits = 10, auto_unbox = T)
