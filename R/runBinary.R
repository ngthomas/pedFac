
#' return the path where pedigraph should be in the R system paths
#'
pedigraph_binary_path <- function() {
  #bin_name <- paste("pedigraph", Sys.info()["sysname"], sep = "-")
  bin_name <- "pedigraph"
  if(Sys.info()["sysname"] == "Windows") {
    bin_name <- paste(bin_name, ".exe", sep = "")
  }
  file.path(system.file(package = "pedfac"), "exec", bin_name)
}

#' return TRUE if pedigraph exists in described path
#'
pedigraph_exists <- function() {
  file.exists(pedigraph_binary_path())
}

#' file path to be used in a call to pedigraph
#'
#' This version checks to make sure it is there and throws an
#' error with a suggestion of how to get it if it is not there.
#' @export
pedigraph_binary <- function() {
  if(!pedigraph_exists()) stop("Can't find the pedigraph executable where it was expected
                             at ", pedigraph_binary_path(), ".")

  pedigraph_binary_path()
}
