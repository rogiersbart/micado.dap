# gamma data ----

set.seed(123456)
rename_gamma <- function(df) dplyr::rename(
    df,
    nuclide = Nuclide,
    energy = Energy,
    time = meas_time,
    detector = Detector,
    net = NetCounts,
    background = BckgCounts
  )
sample_poisson <- function(df) dplyr::mutate(
  df,
  net = rpois(nrow(df), net),
  background = rpois(nrow(df), background)
)
vc1_gamma <- readRDS(system.file("model/data/gamma/vc1/Counts.rds", package = "micado.dap")) |>
  rename_gamma() |>
  sample_poisson()
usethis::use_data(vc1_gamma, overwrite = TRUE)
vc3_gamma <- readRDS(system.file("model/data/gamma/vc3/Counts.rds", package = "micado.dap")) |>
  rename_gamma() |>
  sample_poisson()
usethis::use_data(vc3_gamma, overwrite = TRUE)

# pncc data ----

pncc_obs <- function(true_mu_alph, true_tot_pu, true_am241r, case) {
  obs_CCpassive <- readRDS(system.file("model/data/neutron/", case, "/CCpassive.rds", package = "micado.dap"))[1,2] |>
    as.numeric()
  sel_nucl_pu<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241","Pu-242")
  true_mass<-c(true_am241r*true_tot_pu,true_mu_alph*true_tot_pu/100)
  comp_pn<-true_mass/sum(true_mass)
  # load nuclear data qiqr for sel_nucl_pu
  qiqr_pu<-readRDS(system.file("model/data/neutron/Qi_Qref_full.rds", package = "micado.dap")) |>
    dplyr::arrange(radionuclide) |>
    dplyr::filter(radionuclide %in% sel_nucl_pu) |>
    dplyr::select(`Q/Qref`) |>
    tidyr::replace_na(list(`Q/Qref`=0)) |>
    purrr::as_vector() |>
    unname()
  w_sum_sf_pair_rate<-sum(comp_pn*qiqr_pu)
  m240Pu_eq<-sum(true_mass)*w_sum_sf_pair_rate
  obs_Rexp_m <- m240Pu_eq*obs_CCpassive
  alfa_rexp<-0.01
  sigma_Rexp_m<-alfa_rexp*obs_Rexp_m
  c(reals_per_second = m240Pu_eq*obs_CCpassive + rnorm(1)*sigma_Rexp_m)
}
true_mu_alph<-c(2.029,	58.475,	26.127,	6.210,	7.159)
true_tot_pu<-0.374
true_am241r<-0.0736
vc1_pncc <- pncc_obs(true_mu_alph, true_tot_pu, true_am241r, "vc1")
usethis::use_data(vc1_pncc, overwrite = TRUE)
true_mu_alph<-c(2.356,	55.506,	25.960,	8.975,	7.203)
true_tot_pu<-0.149
true_am241r<-0.00442
vc3_pncc <- pncc_obs(true_mu_alph, true_tot_pu, true_am241r, "vc3")
usethis::use_data(vc3_pncc, overwrite = TRUE)

# ancc data ----

obs_CCactive_239Pu<-readRDS(system.file("model/data/neutron/vc1/CCactive.rds", package = "micado.dap"))[1,2] |>
  as.numeric()
obs_CCactive_241Pu<-readRDS(system.file("model/data/neutron/vc1/CCactive.rds", package = "micado.dap"))[1,3] |>
  as.numeric()
sel_nucl_pu<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241","Pu-242")
true_mass<-c(true_am241r*true_tot_pu,true_mu_alph*true_tot_pu/100)
true_tot_mass <- sum(true_mass)
obs_csa <- 1/true_tot_mass*(true_mass[3]*obs_CCactive_239Pu + true_mass[5]*obs_CCactive_241Pu)
alfa_csa <- 0.01
sigma_csa <- alfa_csa*obs_csa
vc1_ancc <- c(counts_per_second = obs_csa + rnorm(1)*sigma_csa)
usethis::use_data(vc1_ancc, overwrite = TRUE)
