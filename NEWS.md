
# micado.dap 0.2.1

* Upgraded {rui} version requirement

# micado.dap 0.2.0

* Corrected the PNCC and ANCC "true" measurement data.
* Added on load checking of CmdStan availability, and automatic installation if
  it is missing (still depends on availability of Rtools).
* Perturbation of the "true" gamma, pncc and ancc data was added.
* New chi-squared discrepancy p-values introduced as additional model
  diagnostics, using a non-parametric rank-based approach.
* Moved part of the Stan model block to transformed parameters, to avoid
  coding the forward model twice.
* Replaced generated quantities in Stan models by posterior predictive
  residuals.
* Fixed correct identification of measurement combination case in internal
  `prep()`.

# micado.dap 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* New `ip()` helper for determining where the HTTP API can be accessed.
* New `dap_start()`, `dap_stop()`, `dap_on()`, `dap_url()` and `dap_request()`
  for starting the HTTP API server, and making remote requests from an R client.
* New `vc1_model()` and `vc3_model()` interfaces to the probabilistic models.
* Added the "true" measurement data sets `vc1_gamma`, `vc1_pncc`, `vc1_ancc`,
  `vc3_gamma` and `vc3_ancc`.
