server <- NULL

#' Start the MICADO DAP server
#'
#' For accessing the MICADO DAP through the HTTP API, the server has to be
#' started from R first.
#'
#' @inheritParams vc1_model
#' @export
#' @seealso [`dap_stop()`], [`dap_on()`], [`dap_url()`], [`dap_request()`]
dap_start <- function(cores = 1) {
  assignInNamespace(
    "server",
    callr::r_bg(
      function() {
        plumber::plumb(system.file("api/plumber.R", package = "micado.dap"))$run(host = "0.0.0.0", port = 8000, swagger = FALSE)
      },
      env = c(callr::rcmd_safe_env(), "MICADO.DAP.CORES" = cores),
      supervise = TRUE
    ),
    "micado.dap"
  )
  std::err("v The MICADO DAP server was successfully started.")
}

#' Stop the MICADO DAP server
#'
#' The MICADO DAP server can always be stopped explicitly if needed. Normally,
#' however, the background R process will be killed whenever R exits.
#'
#' @export
#' @seealso [`dap_start()`], [`dap_on()`], [`dap_url()`], [`dap_request()`]
dap_stop <- function() {
  kill <- micado.dap:::server$kill()
  if (kill) {
    assignInNamespace("server", NULL, "micado.dap")
    std::err("v The MIDACO DAP server was successfully terminated.")
  } else {
    std::err("x The MIDACO DAP server seems unavailable.")
  }
}

#' Check if the MICADO DAP server is running
#'
#' For testing purposes, one can check whether the background R process that is
#' still serving the HTTP API is still active.
#'
#' @return Logical
#' @export
#' @seealso [`dap_start()`], [`dap_stop()`], [`dap_url()`], [`dap_request()`]
dap_on <- function() {
  ! is.null(micado.dap:::server)
}

#' Get the URL at which the HTTP API can be accessed
#'
#' This is useful to check where to access the HTTP API from a remote machine.
#'
#' @return The MICADO DAP URL
#' @export
#' @seealso [`dap_start()`], [`dap_stop()`], [`dap_on()`], [`dap_request()`], [`ip()`]
dap_url <- function() {
  paste0("http://", ip(), ":8000/dap")
}

#' Submit an HTTP request to the MICADO DAP
#'
#' @inheritParams vc1_model
#' @param ip IP address of the machine that is running the MICADO DAP server.
#'   Can be `"localhost"` (default) for local testing purposes.
#'
#' @return The probabilistic interpretation results, as a list with a `masses`
#'   data frame and a `diagnostics` list.
#' @export
#' @seealso [`dap_start()`], [`dap_stop()`], [`dap_on()`], [`dap_url()`]
dap_request <- function(gamma = NULL, pncc = NULL, ancc = NULL, case = NULL, legacy = FALSE, ip = "localhost") {
  request_json <- jsonlite::toJSON(
    list(
      gamma = gamma,
      pncc = pncc,
      ancc = ancc,
      case = case,
      legacy = legacy,
      client = "R"
    ),
    digits = NA
  )
  url <- httr::parse_url(paste0("http://", ip, ":8000"))
  url$path <- "dap"
  httr::build_url(url) |>
    httr::POST(body = request_json) |>
    httr::content(as = "text", encoding = "UTF-8") |>
    jsonlite::fromJSON()
}

#' Find a machine's IP address through system commands.
#'
#' For the current machine, we use the first IPv4 listed by `ipconfig`, for
#' a remote machine, the address returned by `nslookup` is used.
#'
#' @param machine Name of the machine to get the IP address for. If `NULL`
#'   (default), the current machine is used.
#'
#' @return An IP address.
#' @export
#' @seealso [`dap_url()`]
ip <- function(machine = NULL) {
  if (is.null(machine)) return(ip_this_machine())
  cmd <- paste0("nslookup ", machine)
  txt <- system(cmd, intern = TRUE)
  gsub(".*? ([[:digit:]])", "\\1", tail(txt[grep("Address", txt)], 1))
}

ip_this_machine <- function() {
  gsub(
    ".*? ([[:digit:]])",
    "\\1",
    system("ipconfig", intern=T)[grep("IPv4", system("ipconfig", intern = T))]
  )[1]
}
