.onLoad <- function(libname = find.package("GeneFamilies"), pkgname = "GeneFamilies") {
    data("codingSequences", package = "GeneFamilies")
    data("families", package = "GeneFamilies")
    data("interPro", package = "GeneFamilies")
    data("orthologsTandems", package = "GeneFamilies")
    data("cafe_result", package = "GeneFamilies")
    data("RNA_Seq_RPKM_and_profiles", package = "GeneFamilies")
    data("pairwiseKaKs", package = "GeneFamilies")
    # data('orthologsKaKs', package = 'GeneFamilies')
    data("fubar_results", package = "GeneFamilies")
    data("BUSTED_Results", package = "GeneFamilies")
    data("GeneGroups", package = "GeneFamilies")
    data("IprBasedEntropies", package = "GeneFamilies")
    data("ExpressionProfileDistances", package = "GeneFamilies")
    data("ExpressionProfileDistanceDistributions", package = "GeneFamilies")
    data("ExpressionProfileDistancesPerTissueDistributions", package = "GeneFamilies")
    # After loading the species' CDS define this constant:
    if (exists("ath.cds") && exists("cla.cds") && exists("cme.cds") && exists("csa.cds") && 
        exists("fve.cds")) {
        assign("spec.gene.ids", list(ath = names(ath.cds), cla = names(cla.cds), 
            cme = names(cme.cds), csa = names(csa.cds), fve = names(fve.cds)), 
            envir = globalenv())
    }
}
#'####################
#' Define constants: #
#'####################
#' 'spec.regex.gene.ids' - Regular expressions to identify the species a gene
#' belongs to. They are supposed to be used with the perl=TRUE flag.
#' @export
spec.regex.gene.ids <- list(aet = "^aet_\\d+$", aly = "^((scaffold_)|(fgenesh2_kg\\.)|(Al_scaffold_))$", 
    ath = "^AT[0-9CM]G\\d+(\\.\\d)?$", bra = "^Bra\\d+$", cru = "^PAC:\\d+$", esa = "^Thhalv\\d+m$", 
    tpa = "^(Tp([0-9mc]g)|(_un)\\d+(_\\d+)?)$", chi = "^CARHR\\d+(\\.\\d)?$")

