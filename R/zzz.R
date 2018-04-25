.onLoad <- function(libname = find.package("MaizeGeneFamilies"), pkgname = "MaizeGeneFamilies") {
  data( "maizeUstilagoDeGenes", package="MaizeGeneFamilies" )
  data( "maizeGenome", package="MaizeGeneFamilies" )
  data( "maizeCellCycleHouseKeeping", package="MaizeGeneFamilies" )
  data( "Zmays_pfam_go", package="MaizeGeneFamilies" )
  data( "PFam_to_InterPro_mapping", package="MaizeGeneFamilies" )
  data( "MaizeGeneSets", package="MaizeGeneFamilies" )
}
#'####################
#' Define constants: #
#'####################
