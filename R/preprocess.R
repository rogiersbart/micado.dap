prep <- function(
  gamma = NULL, # NULL or gamma data
  pncc = NULL, # NULL or pncc data
  ancc = NULL, # NULL or ancc data
  vc, # 1 or 3
  case = "realistic" # or "mysterious"
) {
  if (length(gamma) == 0) gamma <- NULL
  if (length(pncc) == 0) pncc <- NULL
  if (length(ancc) == 0) ancc <- NULL
  test_case<-paste0("vc", vc)
  meas_case <- c("gamma", "pncc", "ancc")[! sapply(list(gamma = gamma, pncc = pncc, ancc = ancc), is.null)]
  DoCorruptCounts<-FALSE
  homoscedastic_error<-FALSE
  InformativeSetup<-case
  #### load data ####
  if (identical(meas_case, c("gamma"))) idx_case <- 1
  if (identical(meas_case, c("gamma","pncc"))) idx_case <- 2
  if (identical(meas_case, c("pncc"))) idx_case <- 3
  if (identical(meas_case, c("gamma","ancc"))) idx_case <- 4
  if (identical(meas_case, c("ancc"))) idx_case <- 5
  if (identical(meas_case, c("gamma","pncc","ancc"))) idx_case <- 6
  #### general setup and uncertainty level ####
  if(test_case=="vc1"){
    nseg<-5
    npu<-5
    sel_nucl<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241","Pu-242")
    sel_nucl_gamma<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241")
    idx_gamma<-c(1,2,3,4,5)
    n_nucl<-length(sel_nucl)
    n_nucl_gam<-length(idx_gamma)
    true_mu_alph<-c(2.029,	58.475,	26.127,	6.210,	7.159)
    neff<-4
    vc_case<-0
    true_tot_pu<-0.374
    true_am241r<-0.0736
    e_bias<-0

  }else if(test_case=="vc3"){
    nseg<-6
    npu<-5

    sel_nucl<-c("Am-241", "Co-60", "Cs-134", "Cs-137","Eu-154",
                "Pu-238", "Pu-239","Pu-240","Pu-241",
                "Pu-242","Ru-106", "Sb-125")
    idx_gamma<-c(seq(7),11,12)
    sel_nucl_gamma<-sel_nucl[idx_gamma]
    n_nucl<-length(sel_nucl)
    n_nucl_gam<-length(idx_gamma)
    true_mu_alph<-c(2.356,	55.506,	25.960,	8.975,	7.203)
    neff<-6+1-1
    vc_case<-1
    nbg<-6
    true_tot_pu<-0.149
    true_am241r<-0.00442
    true_co60<-3.061E-06
    true_ru106<-1.745E-05
    true_sb125<-8.537E-06
    true_cs134<- 1.268E-05
    true_cs137<- 1.396E-03
    true_eu154<- 5.722E-06

    e_bias<-0
  }
  if (homoscedastic_error){
    rel_err<-numeric(npu)+0.05
  }else{rel_err<-c(0.02,0.02,0.02,0.02,0.1)}
  noise_v<-1+rnorm(npu)*rel_err

  mu_alph<-true_mu_alph*noise_v
  mu_alph<-100*mu_alph/sum(mu_alph)

  if(InformativeSetup=="full"){
    alph_puv<-true_mu_alph*10000
    if(test_case=="vc1"){
      pu_bounds<-c(0.372,0.376) # truth is 0.374
      alph_puf<-c(0.400, 0.400, 0.075, 0.075, 0.050)*10000
      alph_lam<-c(1,1,1,1)
      am241r_bounds<-c(0.0734,0.0738) # truth is 0.0736
    }else if(test_case=="vc3"){
      pu_bounds<-c(0.147,0.151) # truth is 0.149
      alph_puf <-c(
        0.05775195, 0.04960252,
        0.07858283, 0.71058296,
        0.10347975, 0.00011)
      alph_puf<-alph_puf/sum(alph_puf)
      alph_puf<-alph_puf*10000

      alph_lam<-numeric(neff)+1.0

      am241r_bounds<-c(0.00440,0.00444) # truth is 0.00442

      # truths are:
      # co60<-3.061E-06
      # ru106<-1.745E-05
      # sb125<-8.537E-06
      # cs134<- 1.268E-05
      # cs137<- 1.396E-03
      # eu154<- 5.722E-06

      co60_bounds<-c(3.060E-06,3.062E-06)
      cs134_bounds<-c(1.267E-05,1.269E-05)
      cs137_bounds<-c(1.395E-03,1.397E-03)
      eu154_bounds<-c(5.721E-06,5.723E-06)
      ru106_bounds<-c(1.744E-05,1.746E-05)
      sb125_bounds<-c(8.536E-06,8.538E-06)

      alph_bgf<-c(
        0.08499708, 0.47051953,
        0.00011, 0.18343759,
        0.26104580, 0.00011
      )
      alph_bgf<-alph_bgf/sum(alph_bgf)
      alph_bgf<-alph_bgf*10000


    }
    warmup <- 2000
    max_treedepth <- 10
  }else if(InformativeSetup=="realistic"){
    alph_puv<-mu_alph
    if(test_case=="vc1"){
      pu_bounds<-c(0.01,1)
      alph_puf<-c(1,1,1,1,1)
      alph_lam<-c(1,1,1,1)
      am241r_bounds<-c(0.07,0.08)

    }else if(test_case=="vc3"){

      pu_bounds<-c(0.01,1)
      alph_puf <-c(1,1,1,1,1,1)

      alph_lam<-numeric(neff)+1.0

      am241r_bounds<-c(0.004,0.005)

      co60_bounds<-c(1E-06,1E-05)
      cs134_bounds<-c(0.5E-05,5E-05)
      cs137_bounds<-c(0.5E-03,5E-03)
      eu154_bounds<-c(1E-06,1E-05)
      ru106_bounds<-c(0.5E-05,5E-05)
      sb125_bounds<-c(1E-06,1E-05)

      alph_bgf<-c(1,1,1,1,1,1)

    }
    warmup <- 2000
    max_treedepth <- 10
  }else if(InformativeSetup=="mysterious"){
    alph_puv<-c(1,1,1,1,1)
    if(test_case=="vc1"){
      pu_bounds<-c(0.01,1)
      alph_puf<-c(1,1,1,1,1)
      alph_lam<-c(1,1,1,1)
      am241r_bounds<-c(0.01,0.2)
    }else if(test_case=="vc3"){

      pu_bounds<-c(0.01,1)
      alph_puf <-c(1,1,1,1,1,1)

      alph_lam<-numeric(neff)+1.0

      am241r_bounds<-c(0.001,0.1)

      co60_bounds<-c(1E-06,1E-05)
      cs134_bounds<-c(1E-06,1E-04)
      cs137_bounds<-c(1E-04,1E-02)
      eu154_bounds<-c(1E-06,1E-04)
      ru106_bounds<-c(1E-06,1E-04)
      sb125_bounds<-c(1E-06,1E-04)

      alph_bgf<-c(1,1,1,1,1,1)

    }
    warmup <- 4000
    max_treedepth <- 20
  }
  #### gamma ####
  # considered nuclides and level of prior info

  # efficiencies
  eff <- readRDS(system.file("model/data/gamma/",test_case,"/eff.rds", package = "micado.dap"))
  n_eff<-length(eff) #n_eff is not used
  # each 3D eff array (elem of the eff list) is structured as
  # source location * energy * detector location
  ndet<-dim(eff[[1]])[3]
  nseg<-dim(eff[[1]])[1]
  if(test_case=="vc1"){
    e1<-eff[[1]]
    e2<-eff[[2]]
    e3<-eff[[3]]
    e4<-eff[[4]]
  }else if(test_case=="vc3"){
    e1<-eff[[1]]
    e2<-eff[[2]]
    e3<-eff[[3]]
    e4<-eff[[4]]
    e5<-eff[[5]]
    e6<-eff[[6]]

  }
  # count data,test_case
  dfc <- if (is.null(gamma)) {
    readRDS(system.file("model/data/gamma/",test_case,"/Counts.rds", package = "micado.dap"))
  } else {
    gamma |>
      dplyr::rename(
        Nuclide = nuclide,
        Energy = energy,
        meas_time = time,
        Detector = detector,
        NetCounts = net,
        BckgCounts = background
      )
  }
  npeaks<-length(unique(dfc$Energy))

  # dfc must be sorted first by Detector, then Nuclide and then Energy
  dfc <- dfc |>
    dplyr::arrange(Detector,Nuclide,Energy) |>
    dplyr::mutate(id=rep(seq(1:npeaks),ndet))

  dt<-dfc$meas_time[1]

  # get number of peaks per nuclide
  idx_p<-readRDS(system.file("model/data/gamma/",test_case,"/Counts.rds", package = "micado.dap")) |>
    dplyr::filter(Detector=="D1") |>
    dplyr::arrange(Nuclide, Energy) |>
    dplyr::select(Nuclide) |>
    table()

  idx_p<-as.vector(idx_p)

  obs_net<-dfc$NetCounts
  #sd_net<-dfc$Unc_NetCounts
  obs_bckg<-dfc$BckgCounts
  #sd_back<-dfc$Unc_BckgCounts
  obs_gross_ref<-obs_net+obs_bckg
  #sd_gross<-sqrt(dfc$Unc_NetCounts^2+dfc$Unc_BkgrCounts^2)

  # seed<-1978
  # set.seed(seed)
  if(DoCorruptCounts){ # resample the observed counts from their poisson means
    obs_gross<-numeric(length(obs_gross_ref))
    for (i in 1:length(obs_gross_ref)){
      obs_gross[i]<-rpois(1,obs_gross_ref[i])
    }
  }else{obs_gross<-obs_gross_ref}

  # round and convert to integer for poisson likelihood in stan
  obs_gross<-round(obs_gross)

  # specific activities (Bq/g)
  spec_a<-readRDS(system.file("model/data/gamma/specific_activity.rds", package = "micado.dap")) |>
    dplyr::filter(Nuclide %in% sel_nucl) |>
    dplyr::arrange(Nuclide) |>
    dplyr::select(Value)
  spa<-spec_a$Value

  # REb #
  dfe <- readRDS(system.file("model/data/gamma/",test_case,"/REb.rds", package = "micado.dap")) |>
    dplyr::mutate(id=seq(1:npeaks))
  # dfe is sorted get the indices of the energies in the 3D eff arrays
  # Dimension 2 the 3D eff arrays are the energies peaks sorted in ascending order
  # So we need to get the energy indices when sorted first by nuclide and then by energy
  # = the ide vector

  ide <- dfe |>
    dplyr::arrange(Energy,Nuclide) |>
    dplyr::mutate(id=seq(npeaks)) |>
    dplyr::arrange(Nuclide,Energy) |>
    dplyr::select(id) |>
    purrr::as_vector()

  Ig<-dfe$b

  #### pncc ####
  # load "true" CCpassive
  if (is.null(pncc)) {
    obs_CCpassive <- readRDS(system.file("model/data/neutron/",test_case,"/CCpassive.rds", package = "micado.dap"))[1,2] |>
      as.numeric()
    # Compute observed Rexp_m
    # Among the considered nuclides in both vc1 and vc3, pncc is only
    # sensitive to the even Pu-isotopes, Pu-238, Pu240 and Pu-242
    # here true_mass only contains the Pu-related isotopes:
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

    #observed number of reals per second (what is actually measured)
    obs_Rexp_m <- m240Pu_eq*obs_CCpassive

    alfa_rexp<-0.01
    sigma_Rexp_m<-alfa_rexp*obs_Rexp_m

    if(DoCorruptCounts){
      obs_Rexp_m <- m240Pu_eq*obs_CCpassive + rnorm(1)*sigma_Rexp_m
    }
  } else {
    obs_Rexp_m <- pncc[[1]]
    sigma_Rexp_m <- 0.01 * obs_Rexp_m
  }

  # nuclear data
  Qref<-readRDS(system.file("model/data/neutron/Qref.rds", package = "micado.dap"))
  Qi_Qref<-readRDS(system.file("model/data/neutron/Qi_Qref_full.rds", package = "micado.dap"))

  if(test_case=="vc3"){
    Qi_Qref<- Qi_Qref |>
      tibble::add_row(radionuclide = "Co-60", Q = NA_real_, "Q/Qref"=NA_real_) |>
      tibble::add_row(radionuclide = "Cs-134", Q = NA_real_, "Q/Qref"=NA_real_) |>
      tibble::add_row(radionuclide = "Cs-137", Q = NA_real_, "Q/Qref"=NA_real_) |>
      tibble::add_row(radionuclide = "Eu-154", Q = NA_real_, "Q/Qref"=NA_real_) |>
      tibble::add_row(radionuclide = "Sb-125", Q = NA_real_, "Q/Qref"=NA_real_) |>
      tibble::add_row(radionuclide = "Ru-106", Q = NA_real_, "Q/Qref"=NA_real_) |>
      dplyr::arrange(radionuclide)
  }

  qiqr<-Qi_Qref |>
    dplyr::arrange(radionuclide) |>
    dplyr::filter(radionuclide %in% sel_nucl) |>
    dplyr::select(`Q/Qref`) |>
    tidyr::replace_na(list(`Q/Qref`=0)) |>
    purrr::as_vector() |>
    unname()

  qiqr<-array(rep(qiqr,each=nseg),c(nseg,n_nucl))

  #Fs_prime_240Pu<-readRDS(system.file("model/data/neutron/Fs_prime_240Pu.rds", package = "micado.dap"))
  # for consistency with forward calculations:
  Fs_prime_240Pu <- 478 # instead of 479.5558

  eff<-readRDS(system.file("model/data/neutron/",test_case,"/eff_passive.rds", package = "micado.dap"))

  e1_pncc<-eff[[1]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()
  e2_pncc<-eff[[2]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()
  e3_pncc<-eff[[3]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()
  e4_pncc<-eff[[4]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()

  # from CCpassive to Rsim (efficiency)
  e1_pncc<-e1_pncc*nseg/Fs_prime_240Pu
  e2_pncc<-e2_pncc*nseg/Fs_prime_240Pu
  e3_pncc<-e3_pncc*nseg/Fs_prime_240Pu
  e4_pncc<-e4_pncc*nseg/Fs_prime_240Pu

  if(test_case=="vc3"){
    e5_pncc<-eff[[5]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()
    e6_pncc<-eff[[6]] |> dplyr::select(`CCpassive [c/s/g240Pu]`) |> unlist() |> as.numeric()

    # from CCpassive to Rsim (efficiency)
    e5_pncc<-e5_pncc*nseg/Fs_prime_240Pu
    e6_pncc<-e6_pncc*nseg/Fs_prime_240Pu

  }
  #### ancc ####
  if(test_case=="vc1"){
    if (is.null(ancc)) {
      obs_CCactive_239Pu<-readRDS(system.file("model/data/neutron/",test_case,"/CCactive.rds", package = "micado.dap"))[1,2] |>
        as.numeric()
      obs_CCactive_241Pu<-readRDS(system.file("model/data/neutron/",test_case,"/CCactive.rds", package = "micado.dap"))[1,3] |>
        as.numeric()
      # Among the considered nuclides in vc1, ancc is only
      # sensitive to the fissile materials, Pu239 and Pu-241
      # for vc3, no ancc signal measured
      # here true_mass only contains the Pu-related isotopes:
      sel_nucl_pu<-c("Am-241","Pu-238","Pu-239","Pu-240","Pu-241","Pu-242")
      true_mass<-c(true_am241r*true_tot_pu,true_mu_alph*true_tot_pu/100)
      true_tot_mass <- sum(true_mass)
      # observed number of count per seconds (what is actually measured)
      # obs_csa mass_239Pu*CCactive_239Pu+mass_241Pu*CCactive_241Pu
      obs_csa <- 1/true_tot_mass*(true_mass[3]*obs_CCactive_239Pu + true_mass[5]*obs_CCactive_241Pu)

      alfa_csa <- 0.01
      sigma_csa <- alfa_csa*obs_csa

      if(DoCorruptCounts){
        obs_csa <- obs_csa + rnorm(1)*sigma_csa
      }
    } else {
      obs_csa <- ancc[[1]]
      sigma_csa <- 0.01 * obs_csa
    }

    # to go from inferred CCactive_239Pu, CCactive_241Pu and mass vector to obs_csa:

    # m239Pu_eq = obs_csa / CCactive_239Pu
    # m239Pu_eq = 1/CCactive_239Pu *(CCactive_239Pu*mass_239Pu/tot_mass + CCactive_241Pu*mass_241Pu/tot_mass)
    # sim_csa = m239Pu_eq * CCactive_239Pu

    eff<-readRDS(system.file("model/data/neutron/",test_case,"/eff_active.rds", package = "micado.dap"))

    e11_ancc<-eff[[1]] |> dplyr::select(`CCactive [c/s/g239Pu]`) |> unlist() |> as.numeric()
    e12_ancc<-eff[[2]] |> dplyr::select(`CCactive [c/s/g239Pu]`) |> unlist() |> as.numeric()
    e13_ancc<-eff[[3]] |> dplyr::select(`CCactive [c/s/g239Pu]`) |> unlist() |> as.numeric()
    e14_ancc<-eff[[4]] |> dplyr::select(`CCactive [c/s/g239Pu]`) |> unlist() |> as.numeric()

    e21_ancc<-eff[[1]] |> dplyr::select(`CCactive [c/s/g241Pu]`) |> unlist() |> as.numeric()
    e22_ancc<-eff[[2]] |> dplyr::select(`CCactive [c/s/g241Pu]`) |> unlist() |> as.numeric()
    e23_ancc<-eff[[3]] |> dplyr::select(`CCactive [c/s/g241Pu]`) |> unlist() |> as.numeric()
    e24_ancc<-eff[[4]] |> dplyr::select(`CCactive [c/s/g241Pu]`) |> unlist() |> as.numeric()

    e11_ancc<-e11_ancc*nseg
    e12_ancc<-e12_ancc*nseg
    e13_ancc<-e13_ancc*nseg
    e14_ancc<-e14_ancc*nseg

    e21_ancc<-e21_ancc*nseg
    e22_ancc<-e22_ancc*nseg
    e23_ancc<-e23_ancc*nseg
    e24_ancc<-e24_ancc*nseg
  }
  if (test_case == "vc1") return(
    list(
      Nseg=nseg,
      Npu=npu,
      Nm=ndet*npeaks,
      Neff=neff,
      Nnucl=n_nucl,
      Nnucl_gam=n_nucl_gam,
      Npeaks=npeaks,
      Ndet=ndet,
      alph_puv=alph_puv,
      alph_puf=alph_puf,
      alph_lam=alph_lam,
      obs_net=obs_net,
      obs_bckg=obs_bckg,
      obs_gross=obs_gross,
      idx_gamma=idx_gamma,
      spa=spa,
      e1=e1,
      e2=e2,
      e3=e3,
      e4=e4,
      idx_p=idx_p,
      Ig=Ig,
      ide=ide,
      dt=dt,
      pu_bounds=pu_bounds,
      am241r_bounds=am241r_bounds,
      vc_case=vc_case,
      obs_rexp_m=obs_Rexp_m,
      sigma_rexp_m = sigma_Rexp_m,
      e1_pncc = e1_pncc,
      e2_pncc = e2_pncc,
      e3_pncc = e3_pncc,
      e4_pncc = e4_pncc,
      qiqr = qiqr,
      fs_prime_240Pu=Fs_prime_240Pu,
      obs_csa = obs_csa,
      sigma_csa = sigma_csa,
      e11_ancc = e11_ancc,
      e12_ancc = e12_ancc,
      e13_ancc = e13_ancc,
      e14_ancc = e14_ancc,
      e21_ancc = e21_ancc,
      e22_ancc = e22_ancc,
      e23_ancc = e23_ancc,
      e24_ancc = e24_ancc,
      idx_case=idx_case
    )
  )
  if (test_case == "vc3") return(
    list(
      Nseg=nseg,
      Npu=npu,
      Nm=ndet*npeaks,
      Neff=neff,
      Nnucl=n_nucl,
      Nnucl_gam=n_nucl_gam,
      Npeaks=npeaks,
      Ndet=ndet,
      alph_puv=alph_puv,
      alph_puf=alph_puf,
      alph_lam=alph_lam,
      obs_net=obs_net,
      obs_bckg=obs_bckg,
      obs_gross=obs_gross,
      idx_gamma=idx_gamma,
      spa=spa,
      e1=e1,
      e2=e2,
      e3=e3,
      e4=e4,
      idx_p=idx_p,
      Ig=Ig,
      ide=ide,
      dt=dt,
      pu_bounds=pu_bounds,
      am241r_bounds=am241r_bounds,
      vc_case=vc_case,
      e5=e5,
      e6=e6,
      co60_bounds=co60_bounds,
      cs134_bounds=cs134_bounds,
      cs137_bounds=cs137_bounds,
      eu154_bounds=eu154_bounds,
      ru106_bounds=ru106_bounds,
      sb125_bounds=sb125_bounds,
      Nbg=nbg,
      alph_bgf=alph_bgf,
      obs_rexp_m=obs_Rexp_m,
      sigma_rexp_m = sigma_Rexp_m,
      e1_pncc = e1_pncc,
      e2_pncc = e2_pncc,
      e3_pncc = e3_pncc,
      e4_pncc = e4_pncc,
      e5_pncc = e5_pncc,
      e6_pncc = e6_pncc,
      qiqr = qiqr,
      fs_prime_240Pu=Fs_prime_240Pu,
      idx_case=idx_case,
      e_bias = e_bias
    )
  )
}

