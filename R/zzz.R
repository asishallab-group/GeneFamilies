.onLoad <- function(libname = find.package("GeneFamilies"), pkgname = "MaizeGeneFamilies") {
  data( "maizeUstilagoDeGenes", package="GeneFamilies" )
  data( "maizeGenome", package="GeneFamilies" )
  data( "maizeCellCycleHouseKeeping", package="GeneFamilies" )
  data( "Zmays_pfam_go", package="GeneFamilies" )
  data( "PFam_to_InterPro_mapping", package="GeneFamilies" )
  data( "MaizeGeneSets", package="GeneFamilies" )
  data( "maize_de_mapManBin_overrep", package="GeneFamilies" )
  data( "maizeGeneAnnotationsAndMappings", package="GeneFamilies" )
}
#'####################
#' Define constants: #
#'####################
