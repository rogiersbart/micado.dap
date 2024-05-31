// MICADO - VC3

functions {
  row_vector col_sums(matrix X) {
    row_vector[cols(X)] s;
    s = rep_row_vector(1, rows(X)) * X;
    return s;
  }

  vector row_sums(matrix X) {
    vector[rows(X)] s;
    s = X * rep_vector(1, cols(X));
    return s;
  }
}
data {
  int<lower=0> Nseg;
  int<lower=0> Npu;
  int<lower=0> Nm;
  int<lower=0> Neff;
  int<lower=0> Nnucl;
  int<lower=0> Nnucl_gam;
  int<lower=0> Npeaks;
  int<lower=0> Ndet;
  vector<lower=0>[Npu] alph_puv;
  vector<lower=0>[Nseg] alph_puf;
  vector<lower=0>[Neff] alph_lam;
  vector[Nm] obs_net;
  vector[Nm] obs_bckg;
  array[Nm] int obs_gross;
  array[Nnucl_gam] int idx_gamma;
  vector[Nnucl] spa;
  array[Nseg] matrix[Npeaks, Ndet] e1;
  array[Nseg] matrix[Npeaks, Ndet] e2;
  array[Nseg] matrix[Npeaks, Ndet] e3;
  array[Nseg] matrix[Npeaks, Ndet] e4;
  array[Nnucl_gam] int idx_p;
  row_vector[Npeaks] Ig;
  array[Npeaks] int<lower=0> ide;
  real dt;
  vector[2] pu_bounds;
  vector[2] am241r_bounds;
  // vc_case=0 means VC1,vc_case=1 means VC3
  int<lower=0, upper=1> vc_case;

  // VC3 variables
  array[Nseg] matrix[Npeaks, Ndet] e5;
  array[Nseg] matrix[Npeaks, Ndet] e6;

  vector[2] co60_bounds;
  vector[2] cs134_bounds;
  vector[2] cs137_bounds;
  vector[2] eu154_bounds;
  vector[2] ru106_bounds;
  vector[2] sb125_bounds;
  int<lower=0> Nbg;
  vector<lower=0>[Nbg] alph_bgf;

  // pncc
  real<lower=0> obs_rexp_m;
  real<lower=0> sigma_rexp_m;
  vector[Nseg] e1_pncc;
  vector[Nseg] e2_pncc;
  vector[Nseg] e3_pncc;
  vector[Nseg] e4_pncc;
  vector[Nseg] e5_pncc;
  vector[Nseg] e6_pncc;
  matrix[Nseg, Nnucl] qiqr;
  real fs_prime_240Pu;

  //all
  int idx_case;
  int e_bias;
}
transformed data {
  real<lower=0> upper_pu;
  real<lower=0> lower_pu;
  real<lower=0> upper_am241r;
  real<lower=0> lower_am241r;
  vector[Nm] lower_bckg_st;
  matrix[Nseg, Nnucl] spa_arr;
  vector[Nm] mu_bckg;
  vector[Nm] sd_bckg;
  vector[Nm] mu_bckg_st;
  vector[Nm] sd_bckg_st;

  // VC3 variables
  real<lower=0> upper_co60;
  real<lower=0> lower_co60;
  real<lower=0> upper_cs134;
  real<lower=0> lower_cs134;
  real<lower=0> upper_cs137;
  real<lower=0> lower_cs137;
  real<lower=0> upper_eu154;
  real<lower=0> lower_eu154;
  real<lower=0> upper_ru106;
  real<lower=0> lower_ru106;
  real<lower=0> upper_sb125;
  real<lower=0> lower_sb125;

  // Ne=num_elements(obs_net);
  upper_pu = pu_bounds[2];
  lower_pu = pu_bounds[1];
  upper_am241r = am241r_bounds[2];
  lower_am241r = am241r_bounds[1];
  lower_bckg_st = -1.0 * obs_bckg ./ sqrt(obs_bckg);
  mu_bckg = obs_bckg;
  sd_bckg = sqrt(mu_bckg);
  spa_arr = rep_matrix(to_row_vector(spa), Nseg);

  // VC3
  upper_co60 = co60_bounds[2];
  lower_co60 = co60_bounds[1];
  upper_cs134 = cs134_bounds[2];
  lower_cs134 = cs134_bounds[1];
  upper_cs137 = cs137_bounds[2];
  lower_cs137 = cs137_bounds[1];
  upper_eu154 = eu154_bounds[2];
  lower_eu154 = eu154_bounds[1];
  upper_ru106 = ru106_bounds[2];
  lower_ru106 = ru106_bounds[1];
  upper_sb125 = sb125_bounds[2];
  lower_sb125 = sb125_bounds[1];
}
parameters {
  simplex[Npu] puv;
  simplex[Nseg] puf;
  //simplex[Neff] lam_e;
  array[Nseg] simplex[Neff] lam_e;

  real<lower=0, upper=1> totpu_st;
  //real totpu_st;

  real<lower=0, upper=1> am241r_st;
  vector<lower=lower_bckg_st, upper=100>[Nm] bckg_st;

  // VC3
  real<lower=0, upper=1> co60_st;
  real<lower=0, upper=1> cs134_st;
  real<lower=0, upper=1> cs137_st;
  real<lower=0, upper=1> eu154_st;
  real<lower=0, upper=1> ru106_st;
  real<lower=0, upper=1> sb125_st;
  simplex[Nseg] co60f;
  simplex[Nseg] cs134f;
  simplex[Nseg] cs137f;
  simplex[Nseg] eu154f;
  simplex[Nseg] ru106f;
  simplex[Nseg] sb125f;

  vector[Nseg] ecb;
}
transformed parameters {
  real totpu;
  real am241r;

  // VC3
  real co60;
  real cs134;
  real cs137;
  real eu154;
  real ru106;
  real sb125;

  totpu = log(lower_pu) + totpu_st * (log(upper_pu) - log(lower_pu));
  totpu = exp(totpu);

  am241r = lower_am241r + am241r_st * (upper_am241r - lower_am241r);

  // VC3
  co60 = log(lower_co60) + co60_st * (log(upper_co60) - log(lower_co60));
  co60 = exp(co60);
  cs134 = log(lower_cs134) + cs134_st * (log(upper_cs134) - log(lower_cs134));
  cs134 = exp(cs134);
  cs137 = log(lower_cs137) + cs137_st * (log(upper_cs137) - log(lower_cs137));
  cs137 = exp(cs137);
  eu154 = log(lower_eu154) + eu154_st * (log(upper_eu154) - log(lower_eu154));
  eu154 = exp(eu154);
  ru106 = log(lower_ru106) + ru106_st * (log(upper_ru106) - log(lower_ru106));
  ru106 = exp(ru106);
  sb125 = log(lower_sb125) + sb125_st * (log(upper_sb125) - log(lower_sb125));
  sb125 = exp(sb125);

  // those for the posterior predictive
  real sim_rexp_m;
  vector[Nm] sim_gross;

  { // local?
    // variables //
    matrix[Nseg, Npu] puv_arr;
    vector[Nseg] pu_z;
    vector[Nseg] am241_z;
    matrix[Nseg, Npu] pu_iso_z;
    matrix[Nseg, Nnucl] mass;
    array[Nseg] matrix[Npeaks, Ndet] e0;
    matrix[Nseg, Nnucl_gam] a0;

    matrix[Nseg, Npeaks] ap;
    array[Nseg] matrix[Npeaks, Ndet] lam_e1;
    array[Nseg] matrix[Npeaks, Ndet] lam_e2;
    array[Nseg] matrix[Npeaks, Ndet] lam_e3;
    array[Nseg] matrix[Npeaks, Ndet] lam_e4;
    int istart;
    int iend;
    row_vector[Nm] cps;
    vector[Nm] sim_net;
    vector[Nm] bckg;
    array[Ndet] matrix[Npeaks, Nseg] e00;

    // VC3
    array[Nseg] matrix[Npeaks, Ndet] lam_e5;
    array[Nseg] matrix[Npeaks, Ndet] lam_e6;
    vector[Nseg] co60_z;
    vector[Nseg] cs134_z;
    vector[Nseg] cs137_z;
    vector[Nseg] eu154_z;
    vector[Nseg] ru106_z;
    vector[Nseg] sb125_z;

    // pncc
    matrix[Nseg, Nnucl] comp_pncc;
    vector[Nseg] w_sum_sf_pair_rate;
    vector[Nseg] m240Pu_eq;
    vector[Nseg] e0_pncc;

    // forward model //

    pu_z = totpu * puf;
    am241_z = am241r * pu_z;
    puv_arr = rep_matrix(to_row_vector(puv), Nseg); // Nseg x Npu
    pu_iso_z = puv_arr .* rep_matrix(to_vector(pu_z), Npu); // Nseg x Npu

    co60_z = co60 * co60f;
    ru106_z = ru106 * ru106f;
    sb125_z = sb125 * sb125f;
    cs134_z = cs134 * cs134f;
    cs137_z = cs137 * cs137f;
    eu154_z = eu154 * eu154f;

    mass[ : , 1] = am241_z;
    mass[ : , 2] = co60_z;
    mass[ : , 3] = cs134_z;
    mass[ : , 4] = cs137_z;
    mass[ : , 5] = eu154_z;
    mass[ : , 6 : 10] = pu_iso_z;
    mass[ : , 11] = ru106_z;
    mass[ : , 12] = sb125_z;

    // ** gamma **

    if (idx_case < 3) {
      a0 = spa_arr[ : , idx_gamma] .* mass[ : , idx_gamma];

      for (i in 1 : Nseg) {
        lam_e1[i,  : ,  : ] = rep_matrix(lam_e[i, 1], Npeaks, Ndet);
        lam_e2[i,  : ,  : ] = rep_matrix(lam_e[i, 2], Npeaks, Ndet);
        lam_e3[i,  : ,  : ] = rep_matrix(lam_e[i, 3], Npeaks, Ndet);
        lam_e4[i,  : ,  : ] = rep_matrix(lam_e[i, 4], Npeaks, Ndet);

        lam_e5[i,  : ,  : ] = rep_matrix(lam_e[i, 5], Npeaks, Ndet);
        lam_e6[i,  : ,  : ] = rep_matrix(lam_e[i, 6], Npeaks, Ndet);

        e0[i,  : ,  : ] = lam_e1[i,  : ,  : ] .* e1[i,  : ,  : ]
                          + lam_e2[i,  : ,  : ] .* e2[i,  : ,  : ]
                          + lam_e3[i,  : ,  : ] .* e3[i,  : ,  : ]
                          + lam_e4[i,  : ,  : ] .* e4[i,  : ,  : ]
                          + lam_e5[i,  : ,  : ] .* e5[i,  : ,  : ]
                          + lam_e6[i,  : ,  : ] .* e6[i,  : ,  : ];

        if (e_bias == 1) {
          e0[i,  : ,  : ] = e0[i,  : ,  : ] * ecb[i];
        }

        istart = 0;
        iend = 0;
        for (j in 1 : Nnucl_gam) {
          istart = iend + 1;
          iend = iend + idx_p[j];
          ap[i, istart : iend] = rep_row_vector(a0[i, j], idx_p[j]);
        }
        ap[i,  : ] = ap[i,  : ] .* Ig;
      }

      // this ugly hack is necessary as ap .*e0[,,i] doesn't work
      // switchf from e0 (Nseg,Npeaks,Ndet) to e00 (Ndet,Npeaks,Nseg)
      for (i in 1 : Nseg) {
        for (j in 1 : Npeaks) {
          for (k in 1 : Ndet) {
            e00[k, j, i] = e0[i, j, k];
          }
        }
      }

      for (i in 1 : Ndet) {
        istart = 1 + (i - 1) * Npeaks;
        iend = i * Npeaks;
        cps[istart : iend] = col_sums(ap .* e00[i,  : ,  : ][ide,  : ]');
      }

      sim_net = to_vector(cps * dt);

      bckg = bckg_st .* sd_bckg + mu_bckg;
      sim_gross = sim_net + bckg;

    }

    if (idx_case >= 2) {
      comp_pncc = mass ./ rep_matrix(row_sums(mass), Nnucl);
      w_sum_sf_pair_rate = row_sums(comp_pncc .* qiqr);
      m240Pu_eq = row_sums(mass) .* w_sum_sf_pair_rate;

      for (i in 1 : Nseg) {
        e0_pncc[i] = lam_e[i, 1] .* e1_pncc[i] + lam_e[i, 2] .* e2_pncc[i]
                     + lam_e[i, 3] .* e3_pncc[i] + lam_e[i, 4] .* e4_pncc[i]
                     + lam_e[i, 5] .* e5_pncc[i] + lam_e[i, 6] .* e6_pncc[i];
      }
      sim_rexp_m = fs_prime_240Pu * sum(m240Pu_eq .* e0_pncc);

    }
  } // local?
}
model {

  // (non-uniform) priors
  target += dirichlet_lpdf(puv | alph_puv);
  target += dirichlet_lpdf(puf | alph_puf);
  for (i in 1 : Nseg) {
    target += dirichlet_lpdf(lam_e[i] | alph_lam);
  }
  target += std_normal_lpdf(bckg_st);
  target += dirichlet_lpdf(co60f | alph_bgf);
  target += dirichlet_lpdf(cs134f | alph_bgf);
  target += dirichlet_lpdf(cs137f | alph_bgf);
  target += dirichlet_lpdf(eu154f | alph_bgf);
  target += dirichlet_lpdf(ru106f | alph_bgf);
  target += dirichlet_lpdf(sb125f | alph_bgf);
  target += normal_lpdf(ecb | 1, 0.01);

  // likelihoods
  // ** gamma **
  if (idx_case < 3) {

    // likelihood //
    target += poisson_lupmf(obs_gross | sim_gross);
  }

  if (idx_case >= 2) {

    target += normal_lpdf(obs_rexp_m | sim_rexp_m, sigma_rexp_m);
  }
}
generated quantities {
  // gamma
  array[Nm] int ppr_gross;
  if (idx_case < 3) {
    for (i in 1 : Nm)
      ppr_gross[i] = poisson_rng(sim_gross[i]) - obs_gross[i];
  }
  // pncc
  real ppr_rexp_m;
  if (idx_case >= 2)
    ppr_rexp_m = normal_rng(sim_rexp_m, sigma_rexp_m) - obs_rexp_m;
}
