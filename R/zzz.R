.onAttach <- function(libname, pkgname) {
  std::err("# {.href [MICADO](https://www.micado-project.eu/)} Data Analysis Pipeline")
  std::err("i R package version {packageVersion(\"micado.dap\")}.")
  std::err("i Copyright 2022 {.href [SCK CEN](https://www.sckcen.be/en)}.")
  std::err("i All rights reserved.")
  invisible()
}
.onLoad <- function(libname, pkgname) {
  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE) |> suppressMessages())) {
    std::err("! CmdStan is not found.")
  # FIXME issue with cmdstan installation?
  #   std::err("~ Installing CmdStan")
  #   cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  #   cmdstanr::install_cmdstan(quiet = TRUE)
  #   if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  #     std::err("x")
  #   } else {
  #     std::err("v")
  #   }
  }
}
