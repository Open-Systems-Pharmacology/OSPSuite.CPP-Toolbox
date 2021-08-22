#' Export an individual simulation to C++
#'
#' @param simulation ospsuite simulation to be exported.
#' @param simulationBatchOption ospsuite SimulationBatchOption object with free parameter and molecule information.
#' @param outputName Both filename (without .cpp extension), as well as simulation identifier. Has to be valid C++ identifier and unique among multiple simulations.
#' @param outputPath Directory for C++ file output.
#'
#' @return File name if successful, NULL otherwise.
#'
#' @export
exportSimulationCpp <- function(simulation, simulationBatchOption, outputName, outputPath = getwd()) {
  batchFactory <- rClr::clrCallStatic("OSPSuite.R.Api", "GetSimulationBatchFactory")
  batchSim <- clrCall(batchFactory, "Create", simulation$ref, simulationBatchOption$ref)
  clrCall(batchSim, "ExportToCPPCode", outputPath, FALSE, outputName)
  filename <- paste0(file.path(outputPath, outputName), ".cpp")
  if(file.exists(filename)) {
    return(filename)
  } else {
    return(NULL)
  }
}

#' Compile a simulation to a shared library
#'
#' @param filename Simulation to be compiled.
#'
#' @return library name if successful, NULL otherwise.
#'
#' @export
compileSimulationCpp <- function(filename) {
  cppInc <- paste("-I", system.file("include", package = "ospsuiteCpp"))
  cppModel <- paste0(system.file("include", package = "ospsuiteCpp"), "/model.cpp")
  libPath <- system.file("libs", package = "ospsuiteCpp")
  libName <- list.files(path = libPath, pattern = "ospsuite.*", recursive = TRUE, full.names = TRUE)
  cppLib <- paste0("-L",dirname(libName), " -l", tools::file_path_sans_ext(basename(libName)))

  libExt <- tools::file_ext(basename(libName))
  system(paste("g++ -fPIC -O2 -Wall  -mfpmath=sse -msse2", cppInc, cppLib, cppModel, filename, "-shared -o ",
               paste0(tools::file_path_sans_ext(basename(filename)), ".", libExt)))

  if(file.exists(libName)) {
    return(libName)
  } else {
    return(NULL)
  }
}
