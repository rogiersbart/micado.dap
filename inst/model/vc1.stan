// MICADO - VC1

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

  // pncc
  real<lower=0> obs_rexp_m;
  real<lower=0> sigma_rexp_m;
  vector[Nseg] e1_pncc;
  vector[Nseg] e2_pncc;
  vector[Nseg] e3_pncc;
  vector[Nseg] e4_pncc;
  matrix[Nseg, Nnucl] qiqr;
  real fs_prime_240Pu;

  //ancc
  real<lower=0> obs_csa;
  real<lower=0> sigma_csa;
  vector[Nseg] e11_ancc;
  vector[Nseg] e12_ancc;
  vector[Nseg] e13_ancc;
  vector[Nseg] e14_ancc;
  vector[Nseg] e21_ancc;
  vector[Nseg] e22_ancc;
  vector[Nseg] e23_ancc;
  vector[Nseg] e24_ancc;

  //all
  int idx_case;
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

  upper_pu = pu_bounds[2];
  lower_pu = pu_bounds[1];
  upper_am241r = am241r_bounds[2];
  lower_am241r = am241r_bounds[1];
  lower_bckg_st = -1.0 * obs_bckg ./ sqrt(obs_bckg);
  mu_bckg = obs_bckg;
  sd_bckg = sqrt(mu_bckg);
  spa_arr = rep_matrix(to_row_vector(spa), Nseg);
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
}
transformed parameters {
  real totpu;
  real am241r;

  totpu = log(lower_pu) + totpu_st * (log(upper_pu) - log(lower_pu));
  totpu = exp(totpu);

  am241r = lower_am241r + am241r_st * (upper_am241r - lower_am241r);

  // those for the posterior predictive
  vector[Nm] sim_gross;
  real sim_rexp_m;
  real sim_csa;

  {
    // variables //
    //vector[Npu] puv_arr[Nseg];

    matrix[Nseg, Npu] puv_arr;
    //matrix[Nseg,Neff] lam_e_arr;
    vector[Nseg] pu_z;
    vector[Nseg] am241_z;
    //vector[Nnucl] mass[Nseg];
    matrix[Nseg, Npu] pu_iso_z;
    matrix[Nseg, Nnucl] mass;
    array[Nseg] matrix[Npeaks, Ndet] e0;
    //real e0[Nseg,Npeaks,Ndet];
    matrix[Nseg, Nnucl_gam] a0;

    matrix[Nseg, Npeaks] ap;
    //real lam_e1[Nseg,Npeaks,Ndet];
    array[Nseg] matrix[Npeaks, Ndet] lam_e1;
    array[Nseg] matrix[Npeaks, Ndet] lam_e2;
    array[Nseg] matrix[Npeaks, Ndet] lam_e3;
    array[Nseg] matrix[Npeaks, Ndet] lam_e4;
    int istart;
    int iend;
    row_vector[Nm] cps;
    array[Ndet] matrix[Npeaks, Nseg] e00;
    vector[Nm] sim_net;
    vector[Nm] bckg;

    // pncc
    matrix[Nseg, Nnucl] comp_pncc;
    vector[Nseg] w_sum_sf_pair_rate;
    vector[Nseg] m240Pu_eq;
    vector[Nseg] e0_pncc;

    //ancc
    matrix[Nseg, Nnucl] comp_ancc;
    vector[Nseg] m239Pu_eq; // m_ref
    vector[Nseg] e01_ancc; // CCactive_239Pu
    vector[Nseg] e02_ancc; // CCactive_241Pu

    // forward model //

    pu_z = totpu * puf;
    am241_z = am241r * pu_z;
    puv_arr = rep_matrix(to_row_vector(puv), Nseg); // Nseg x Npu
    pu_iso_z = puv_arr .* rep_matrix(to_vector(pu_z), Npu); // Nseg x Npu

    mass = append_col(am241_z, pu_iso_z);

    // ** gamma **
    if (idx_case < 3 || idx_case == 4 || idx_case == 6) {
      a0 = spa_arr[ : , idx_gamma] .* mass[ : , idx_gamma];

      for (i in 1 : Nseg) {
        lam_e1[i,  : ,  : ] = rep_matrix(lam_e[i, 1], Npeaks, Ndet);
        lam_e2[i,  : ,  : ] = rep_matrix(lam_e[i, 2], Npeaks, Ndet);
        lam_e3[i,  : ,  : ] = rep_matrix(lam_e[i, 3], Npeaks, Ndet);
        lam_e4[i,  : ,  : ] = rep_matrix(lam_e[i, 4], Npeaks, Ndet);

        e0[i,  : ,  : ] = lam_e1[i,  : ,  : ] .* e1[i,  : ,  : ]
                          + lam_e2[i,  : ,  : ] .* e2[i,  : ,  : ]
                          + lam_e3[i,  : ,  : ] .* e3[i,  : ,  : ]
                          + lam_e4[i,  : ,  : ] .* e4[i,  : ,  : ];

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
      // switch from e0 (Nseg,Npeaks,Ndet) to e00 (Ndet,Npeaks,Nseg)
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
    // ** pncc **
    if (idx_case == 2 || idx_case == 3 || idx_case == 6) {
      comp_pncc = mass ./ rep_matrix(row_sums(mass), Nnucl);
      w_sum_sf_pair_rate = row_sums(comp_pncc .* qiqr);
      m240Pu_eq = row_sums(mass) .* w_sum_sf_pair_rate;

      for (i in 1 : Nseg) {
        e0_pncc[i] = lam_e[i, 1] .* e1_pncc[i] + lam_e[i, 2] .* e2_pncc[i]
                     + lam_e[i, 3] .* e3_pncc[i] + lam_e[i, 4] .* e4_pncc[i];
      }
      sim_rexp_m = fs_prime_240Pu * sum(m240Pu_eq .* e0_pncc);
    }
    // ** ancc **
    if (idx_case >= 4) {
      comp_ancc = mass ./ rep_matrix(row_sums(mass), Nnucl);

      for (i in 1 : Nseg) {
        e01_ancc[i] = lam_e[i, 1] .* e11_ancc[i] + lam_e[i, 2] .* e12_ancc[i]
                      + lam_e[i, 3] .* e13_ancc[i]
                      + lam_e[i, 4] .* e14_ancc[i];
        e02_ancc[i] = lam_e[i, 1] .* e21_ancc[i] + lam_e[i, 2] .* e22_ancc[i]
                      + lam_e[i, 3] .* e23_ancc[i]
                      + lam_e[i, 4] .* e24_ancc[i];
      }

      m239Pu_eq = 1.0 ./ e01_ancc
                  .* (e01_ancc .* comp_ancc[ : , 3]
                      + e02_ancc .* comp_ancc[ : , 5]);
      sim_csa = sum(m239Pu_eq .* e01_ancc) * 1.0 ./ Nseg;
    }
  }
}
model {
  // (non-uniform) priors //

  target += dirichlet_lpdf(puv | alph_puv);
  target += dirichlet_lpdf(puf | alph_puf);
  for (i in 1 : Nseg) {
    //lam_e[i] ~ dirichlet(alph_lam);
    target += dirichlet_lpdf(lam_e[i] | alph_lam);
  }
  target += std_normal_lpdf(bckg_st);

  if (idx_case < 3 || idx_case == 4 || idx_case == 6) {
    // likelihood //
    target += poisson_lupmf(obs_gross | sim_gross);
  }
  if (idx_case == 2 || idx_case == 3 || idx_case == 6) {
    target += normal_lpdf(obs_rexp_m | sim_rexp_m, sigma_rexp_m);
  }
  if (idx_case >= 4) {
    target += normal_lpdf(obs_csa | sim_csa, sigma_csa);
  }
}
generated quantities {
  // gamma
  array[Nm] int ppr_gross;
  if (idx_case < 3 || idx_case == 4 || idx_case == 6) {
    for (i in 1 : Nm) {
      ppr_gross[i] = poisson_rng(sim_gross[i]) - obs_gross[i];
    }
  }
  // pncc
  real ppr_rexp_m;
  if (idx_case == 2 || idx_case == 3 || idx_case == 6) {
    ppr_rexp_m = normal_rng(sim_rexp_m, sigma_rexp_m) - obs_rexp_m;
  }
  // ancc
  real ppr_csa;
  if (idx_case >= 4) {
    ppr_csa = normal_rng(sim_csa, sigma_csa) - obs_csa;
  }
}
