.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::red(
      "Note: This package 'concordance' is deprecated and no longer maintained.\nPlease use the new 'coactivity' package going forward.\n\ndevtools::install_github(\"knrumsey/coactivity\")"
    )
  )
}
