#* MICADO DAP Pipeline Web API
#* @post /dap
#* @serializer text
function(req) {
  # req <- jsonlite::fromJSON(req$body)
  # return(jsonlite::toJSON(list(case = req$body$case)))
  library(micado.dap)
  if (req$body$case == "vc1") {
    output <- vc1_model(
      req$body$gamma,
      req$body$pncc,
      req$body$ancc,
      req$body$legacy
    )
  }
  if (req$body$case == "vc3") {
    output <- vc3_model(
      req$body$gamma,
      req$body$pncc,
      req$body$legacy
    )
  }
  if ("client" %in% names(req$body)) {
    if (req$body$client == "R") {
      return(jsonlite::toJSON(list(masses = output$masses, diagnostics = output$diagnostics)))
    }
  }
  get_json(output)
}
