.onAttach <- function(libname, pkgname) {
  packageStartupMessage("For guided examples, see vignette('bkmrGSOverview'). This package is based on the bkmr package version 0.2.2. For a BKMR analysis without effect modification, use the latest bkmr package.")
}

release_questions <- function() {
  c(
    "Have you updated the vignette and posted to GitHub?"
  )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables("ex_data")
}