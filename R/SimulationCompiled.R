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
      private$sim <- .CompiledSimulation_load(name)
    },
    #' Run the associated simulation with a set of given parameters.
    #'
    #' @param parameter A vector of all free parameters.
    #'
    #' @returns A matrix of observations.
    run = function(parameter) {
      return(.CompiledSimulation_run(private$sim, parameter))
    },
    #' Get the current vector of free p parameters.
    #'
    #' @returns A vector of all free parameters.
    getP = function() {
      return(.CompiledSimulation_getP(private$sim))
    },
    #' Get the current vector integration output time points
    #'
    #' @param parameter A vector of all free parameters.
    #'
    #' @returns A vector of all current output time points.
    getT = function() {
      return(.CompiledSimulation_getT(private$sim))
    },
    #' Set a new time grid for the integration.
    #'
    #' @param timeGrid The vector of the new time grid for the integration.
    setT = function(timeGrid) {
      .CompiledSimulation_setT(private$sim, timeGrid)
    }
  )
)
