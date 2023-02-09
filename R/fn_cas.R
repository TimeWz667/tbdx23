
calc_prev <- function(df) {
  df %>% 
    mutate(
      inc = inc0 * exp(- adr * (2023 - Year0)),

      ra = r_sc + r_death_a + r_death_bg,
      rs = r_sc + r_death_s + r_death_bg,
      rc = r_sc + r_death_s + r_death_bg,
      
      r_aware = r_aware0 * rr_det_t ^ (2023 - Year0),
      r_det = r_det0 * rr_det_t ^ (2023 - Year0),
      
      prv_a = inc / (r_sym + ra - adr),
      prv_s = prv_a * r_sym / (r_aware + rs - adr),
      prv_c = prv_s * r_aware / (r_det + rc - adr),
      pdx0 = 0, 
      pdx1 = 1
    )
}


calc_reform <- function(df, pdx0 = 0.4, pdx1 = 0.7) {
  df %>% 
    select(- pdx0, - pdx1) %>% 
    mutate(
      det = r_det * prv_c,
      pdx0 = pdx0,
      pdx1 = pdx1,
      r_csi = r_aware,
      det0 = r_csi * pdx0 * prv_s,
      det1 = det - det0,
      fn0 = r_csi * (1 - pdx0) * prv_s,
      r_recsi = det1 / (pdx1 * prv_c)
    ) %>% 
    select(-det, -det0, -det1, -fn0)
}


calc_intv <- function(df, or_pdx0 = 1, or_pdx1 = 1, 
                      rr_csi = 1, rr_recsi = 1) {
  df %>% 
    mutate(
      odd = or_pdx0 * pdx0 / (1 - pdx0),
      pdx0 = odd / (1 + odd),
      odd = or_pdx1 * pdx1 / (1 - pdx1),
      pdx1 = odd / (1 + odd),
      r_csi = rr_csi * r_csi,
      r_recsi = rr_recsi * r_recsi
    ) %>% 
    select(-odd)
}


calc_cascade <- function(df, sc = "baseline") {
  df %>% 
    mutate(
      inc = (r_sym + ra - adr) * prv_a,
      det0 = r_csi * pdx0 * prv_s,
      det1 = r_recsi * pdx1 * prv_c,
      dur_a = 1 / (r_sym + ra),
      dur_s = 1 / (r_csi + rs),
      dur_c = 1 / (r_recsi * pdx1 + rc),
      drop_a = ra * dur_a,
      drop_s = rs * dur_s,
      drop_c = rc * dur_c,
      
      pd0 = r_csi * pdx0 * dur_s,
      pd1 = r_recsi * pdx1 * dur_c,
      k = pd0 + pd1,
      
      Delay_Pat = dur_s,
      Delay_Sys = pd1 / k * dur_c,
      Delay_Tot = Delay_Pat + Delay_Sys,
      
      Pr_Det = r_sym * dur_a * pd0,
      Pr_Det = Pr_Det + r_sym * dur_a * (r_csi * (1 - pdx0) * dur_s) * pd1,
      Pr_Notif = Pr_Det * (1 - p_under),
      Scenario = sc
    )
}


move_avg2 <- function(x) {
  (x[-1] + x[-length(x)]) / 2
}

