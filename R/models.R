#' Probabilistic model for Virtual Case 1
#'
#' @param gamma Data frame with the GAMMA measurement data, including columns
#'   nuclide, energy, time, detector, net and background. See [`vc1_gamma`] for
#'   an example.
#' @param pncc PNCC reals per second. See [`vc1_pncc`] for an example.
#' @param ancc ANCC counts per second. See [`vc1_ancc`] for an example.
#' @param case The virtual case to consider. Either `"vc1"` or `"vc3"`.
#' @param legacy Logical flag to indicate whether the Plutonium isotopic vector
#'   should be considered unknown or not (defaults to `FALSE`).
#' @param cores The number of cores to use for the probabilistic models. Can be
#'   1 (default), 2 or 4.
#' @param quiet Whether to suppress the Stan messages. Defaults to `TRUE`.
#'
#' @return A list with a `masses` data frame, and a `diagnostics` list.
#' @export
#' @seealso [`vc3_model()`], [`get_fit()`], [`get_json()`]
vc1_model <- function(gamma = NULL, pncc = NULL, ancc = NULL, legacy = FALSE, cores = 1, quiet = TRUE) {
  on.exit(std::err("x"))
  std::err("~ Preparing data")
  if (! Sys.getenv("MICADO.DAP.CORES") == "") cores <- as.numeric(Sys.getenv("MICADO.DAP.CORES"))
  model_data <- prep(gamma, pncc, ancc, vc = 1, case = if (legacy) "mysterious" else "realistic")
  std::err("v")
  std::err("~ Loading model")
  model_vc1 <- cmdstanr::cmdstan_model(
    # stan_file = system.file("model/vc1.stan", package = "micado.dap"),
    exe_file = system.file("model/vc1.exe", package = "micado.dap"),
    compile = FALSE
  )
  std::err("v")
  std::err("~ Sampling")
  capture <- if (quiet) capture.output else I
  capture(
    model_fit <- model_vc1$sample(
      iter_warmup = if (legacy) 4000 else 2000,
      iter_sampling = 1000,
      data = model_data,
      seed = 0,
      chains = 4,
      parallel_chains = cores,
      max_treedepth = if (legacy) 20 else 10,
      show_messages = !quiet,
      show_exceptions = !quiet,
      refresh = 0
    )
  )
  std::err("v")
  std::err("~ Postprocessing")
  masses <- model_fit$draws(
    c("totpu", "am241r", paste0("puv[", 1:5, "]")),
    format = "df"
  ) |>
    tibble::as_tibble() |>
    dplyr::select(-dplyr::starts_with(".")) |>
    dplyr::transmute(
      pu238 = totpu * `puv[1]`,
      pu239 = totpu * `puv[2]`,
      pu240 = totpu * `puv[3]`,
      pu241 = totpu * `puv[4]`,
      pu242 = totpu * `puv[5]`,
      am241 = totpu * am241r
    ) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "nuclide") |>
    dplyr::group_by(nuclide) |>
    dplyr::summarise(
      mean = mean(value),
      sd = sd(value),
      q025 = quantile(value, 0.025, names = FALSE),
      q975 = quantile(value, 0.975, names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(nuclide = nuclide_pretty(nuclide))
  p_gamma <- NULL
  if (!is.null(gamma)) {
    ppr <- model_fit$draws("ppr_gross", format = "df") |>
        dplyr::as_tibble() |>
        dplyr::select(dplyr::starts_with("ppr_"))
    ppr_scores_squared <- ppr |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) score_rank(x)^2))
    obs_scores <- ppr |>
      dplyr::summarise(dplyr::across(dplyr::everything(), score_zero)) |>
      as.numeric()
    ssq_obs <- sum(obs_scores^2)
    ssq_pp <- rowSums(ppr_scores_squared)
    p_gamma <- sum(ssq_pp > ssq_obs) / length(ssq_pp)
  }
  p_pncc <- NULL
  if (!is.null(pncc)) {
    ppr <- model_fit$draws("ppr_rexp_m", format = "df") |>
      dplyr::as_tibble() |>
      dplyr::select(dplyr::starts_with("ppr_"))
    ppr_scores_squared <- ppr |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) score_rank(x)^2))
    obs_scores <- ppr |>
      dplyr::summarise(dplyr::across(dplyr::everything(), score_zero)) |>
      as.numeric()
    ssq_obs <- sum(obs_scores^2)
    ssq_pp <- rowSums(ppr_scores_squared)
    p_pncc <- sum(ssq_pp > ssq_obs) / length(ssq_pp)
  }
  p_ancc <- NULL
  if (!is.null(ancc)) {
    ppr <- model_fit$draws("ppr_csa", format = "df") |>
      dplyr::as_tibble() |>
      dplyr::select(dplyr::starts_with("ppr_"))
    ppr_scores_squared <- ppr |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) score_rank(x)^2))
    obs_scores <- ppr |>
      dplyr::summarise(dplyr::across(dplyr::everything(), score_zero)) |>
      as.numeric()
    ssq_obs <- sum(obs_scores^2)
    ssq_pp <- rowSums(ppr_scores_squared)
    p_ancc <- sum(ssq_pp > ssq_obs) / length(ssq_pp)
  }
  output <- list(
    masses = masses,
    diagnostics = list(
      converged = all(model_fit$summary(NULL, "rhat")$rhat < 1.02, na.rm = TRUE),
      p_gamma = p_gamma,
      p_pncc = p_pncc,
      p_ancc = p_ancc
    )
  )
  attr(output, "fit") <- model_fit
  attr(output, "json") <- to_json(output)
  std::err("v")
  output
}

#' Probabilistic model for Virtual Case 3
#'
#' @inheritParams vc1_model
#'
#' @return A list with a `masses` data frame, and a `diagnostics` list.
#' @export
#' @seealso [`vc1_model()`], [`get_fit()`], [`get_json()`]
vc3_model <- function(gamma = NULL, pncc = NULL, legacy = FALSE, cores = 1, quiet = TRUE) {
  on.exit(std::err("x"))
  std::err("~ Preparing data")
  if (! Sys.getenv("MICADO.DAP.CORES") == "") cores <- as.numeric(Sys.getenv("MICADO.DAP.CORES"))
  model_data <- prep(gamma, pncc, vc = 3, case = if (legacy) "mysterious" else "realistic")
  std::err("v")
  std::err("~ Loading model")
  model_vc3 <- cmdstanr::cmdstan_model(
    # stan_file = system.file("model/vc3.stan", package = "micado.dap"),
    exe_file = system.file("model/vc3.exe", package = "micado.dap"),
    compile = FALSE
  )
  std::err("v")
  std::err("~ Sampling")
  capture <- if (quiet) capture.output else I
  capture(
    model_fit <- model_vc3$sample(
      iter_warmup = if (legacy) 4000 else 2000,
      iter_sampling = 1000,
      data = model_data,
      seed = 0,
      chains = 4,
      parallel_chains = cores,
      max_treedepth = if (legacy) 20 else 10,
      show_messages = !quiet,
      show_exceptions = !quiet,
      refresh = 0
    )
  )
  std::err("v")
  std::err("~ Postprocessing")
  masses <- model_fit$draws(
    c("totpu", "am241r", paste0("puv[", 1:5, "]"), "co60", "ru106", "sb125", "cs134", "cs137", "eu154"),
    format = "df"
  ) |>
    tibble::as_tibble() |>
    dplyr::select(-dplyr::starts_with(".")) |>
    dplyr::transmute(
      pu238 = totpu * `puv[1]`,
      pu239 = totpu * `puv[2]`,
      pu240 = totpu * `puv[3]`,
      pu241 = totpu * `puv[4]`,
      pu242 = totpu * `puv[5]`,
      am241 = totpu * am241r,
      co60, ru106, sb125, cs134, cs137, eu154
    ) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "nuclide") |>
    dplyr::group_by(nuclide) |>
    dplyr::summarise(
      mean = mean(value),
      sd = sd(value),
      q025 = quantile(value, 0.025, names = FALSE),
      q975 = quantile(value, 0.975, names = FALSE),
      .groups = "drop"
    ) |>
    dplyr::mutate(nuclide = nuclide_pretty(nuclide))
  p_gamma <- NULL
  if (!is.null(gamma)) {
    ppr <- model_fit$draws("ppr_gross", format = "df") |>
      dplyr::as_tibble() |>
      dplyr::select(dplyr::starts_with("ppr_"))
    ppr_scores_squared <- ppr |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) score_rank(x)^2))
    obs_scores <- ppr |>
      dplyr::summarise(dplyr::across(dplyr::everything(), score_zero)) |>
      as.numeric()
    ssq_obs <- sum(obs_scores^2)
    ssq_pp <- rowSums(ppr_scores_squared)
    p_gamma <- sum(ssq_pp > ssq_obs) / length(ssq_pp)
  }
  p_pncc <- NULL
  if (!is.null(pncc)) {
    ppr <- model_fit$draws("ppr_rexp_m", format = "df") |>
      dplyr::as_tibble() |>
      dplyr::select(dplyr::starts_with("ppr_"))
    ppr_scores_squared <- ppr |>
      dplyr::mutate(dplyr::across(dplyr::everything(), \(x) score_rank(x)^2))
    obs_scores <- ppr |>
      dplyr::summarise(dplyr::across(dplyr::everything(), score_zero)) |>
      as.numeric()
    ssq_obs <- sum(obs_scores^2)
    ssq_pp <- rowSums(ppr_scores_squared)
    p_pncc <- sum(ssq_pp > ssq_obs) / length(ssq_pp)
  }
  p_ancc <- NULL
  output <- list(
    masses = masses,
    diagnostics = list(
      converged = all(model_fit$summary(NULL, "rhat")$rhat < 1.02, na.rm = TRUE),
      p_gamma = p_gamma,
      p_pncc = p_pncc,
      p_ancc = p_ancc
    )
  )
  attr(output, "fit") <- model_fit
  attr(output, "json") <- to_json(output)
  std::err("v")
  output
}

#' Access the underlying CmdStanMCMC object
#'
#' @param model_result Output of a call to [`vc1_model()`] or [`vc3_model()`]
#'
#' @return The underlying `CmdStanMCMC` object
#' @export
#' @seealso [`vc1_model()`], [`vc3_model()`], [`get_json()`]
get_fit <- function(model_result) {
  attr(model_result, "fit")
}

#' Access the JSON for DigiWaste Bridge
#'
#' @param model_result Output of a call to [`vc1_model()`] or [`vc3_model()`]
#'
#' @return The JSON for DigiWaste Bridge
#' @export
#' @seealso [`vc1_model()`], [`vc3_model()`], [`get_fit()`]
get_json <- function(model_result) {
  attr(model_result, "json")
}

nuclide_pretty <- function(x) {
  first <- substr(x, 1, 1) |> toupper()
  second <- substr(x, 2, 2)
  number <- substr(x, 3, nchar(x))
  paste0(first, second, "-", number)
}
nuclide_json <- function(x) {
  c(
    `Pu-238` = "Plutonium238",
    `Pu-239` = "Plutonium239",
    `Pu-240` = "Plutonium240",
    `Pu-241` = "Plutonium241",
    `Pu-242` = "Plutonium242",
    `Am-241` = "Americium241",
    `Co-60`  = "Cobalt60",
    `Ru-106` = "Ruthenium106",
    `Sb-125` = "Antimony125",
    `Cs-134` = "Cesium134",
    `Cs-137` = "Cesium137",
    `Eu-154` = "Europium154"
  )[x]
}
to_json <- function(output) {
  masses <- output$masses |>
    dplyr::rename(val = mean, error = sd) |>
    dplyr::rowwise() |>
    dplyr::mutate(`CI95%` = list(c(q025, q975)),
                  nuclide = nuclide_json(nuclide)) |>
    dplyr::select(-q025, -q975) |>
    dplyr::ungroup()
  masses <- split(masses |> dplyr::select(-nuclide), masses$nuclide) |>
    purrr::map(~list(
      type = "mass",
      unit = "g",
      val = .$val,
      error = .$error,
      `CI95%` = .$`CI95%`[[1]] |> unname()
    ))
  diagnostics <- list(
    converged = output$diagnostics$converged,
    p_gamma = output$diagnostics$p_gamma,
    p_pncc = output$diagnostics$p_pncc,
    p_ancc = output$diagnostics$p_ancc
  )
  jsonlite::toJSON(
    masses |> append(list(Diagnostics = diagnostics)),
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = I(4)
  )
}

recompile <- function(model = c("vc1", "vc3")) {
  if ("vc1" %in% model) {
    model_vc1 <- cmdstanr::cmdstan_model(
      stan_file = system.file("model/vc1.stan", package = "micado.dap"),
      exe_file = system.file("model/vc1.exe", package = "micado.dap"),
      compile = FALSE
    )
    model_vc1$compile()
  }
  if ("vc3" %in% model) {
    model_vc3 <- cmdstanr::cmdstan_model(
      stan_file = system.file("model/vc3.stan", package = "micado.dap"),
      exe_file = system.file("model/vc3.exe", package = "micado.dap"),
      compile = FALSE
    )
    model_vc3$compile()
  }
}
# percent_rank <- function(x) rank(x) / length(x) - 1 / length(x) / 2
# score_rank <- function(x) qnorm(percent_rank(x))
# percent_zero <- function(x) ecdf(x)(0)
# score_zero <- function(x) qnorm(percent_zero(x))
score_rank <- function(x) (x - mean(x))/sd(x)
score_zero <- function(x) (0 - mean(x))/sd(x)
