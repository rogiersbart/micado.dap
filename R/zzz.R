.onAttach <- function(libname, pkgname) {
  rui::console("# {.href [MICADO](https://www.micado-project.eu/)} Data Analysis Pipeline")
  rui::console("i R package version {packageVersion(\"micado.dap\")}.")
  rui::console("i Copyright 2022 {.href [SCK CEN](https://www.sckcen.be/en)}.")
  rui::console("i All rights reserved.")
  invisible()
}
.onLoad <- function(libname, pkgname) {
  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    rui::console("! CmdStan is not found.")
    rui::console("~ Installing CmdStan")
    cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
    cmdstanr::install_cmdstan(quiet = TRUE)
    if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
      rui::console("x")
    } else {
      rui::console("v")
    }
  }
}
