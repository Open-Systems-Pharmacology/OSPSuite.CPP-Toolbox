#' Compile a simulation to a shared library
#'
#' @example
#'
#' sim = SimulationCompiled$new('Aciclovir')
#' parameter = sim$getP()
#' sim$run(parameter)
#' tgrid = sim$getT()
#' sim$setT(c(tgrid, tgrid[length(tgrid)]+1.0))
#' paramter[1] = paramter[1] + 0.1
#' sim$run(parameter)
#'
#' @export
SimulationCompiled <- R6::R6Class("SimulationCompiled",
  private = list(
    sim = NULL
  ),
  public = list(
    #' Constructor for SimulationCompiled class
    #'
    #' @param name Name of the shared library of a compiled model created by exportSimulationCpp and compileSimulationCpp.
    #'
    initialize = function(name) {
      #stopifnot(is.character(name), length(name) == 1) # type checks should be done by Rcpp
      self$sim <- .CompiledSimulation_load(name)
    },
    run = function(parameter) {
      return(.CompiledSimulation_run(self$sim,p))
    },
    getP = function() {
      return(.CompiledSimulation_getP(self$sim))
    },
    getT = function() {
      return(.CompiledSimulation_getT(self$sim))
    },
    setT = function(timeGrid) {
      .CompiledSimulation_setT(self$sim, timeGrid)
    }
  )
)
