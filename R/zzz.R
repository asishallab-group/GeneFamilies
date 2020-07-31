.onLoad <- function(libname = find.package("GeneFamilies"), pkgname = "GeneFamilies") {
    data("codingSequences", package = "GeneFamilies")
    data("families", package = "GeneFamilies")
    data("interPro", package = "GeneFamilies")
    data("familyHumanReadableDescriptions", package = "GeneFamilies")
    data("orthologsTandems", package = "GeneFamilies")
    data("cafe_result", package = "GeneFamilies")
    data("RNA_Seq_RPKM_and_profiles", package = "GeneFamilies")
    data("pairwiseKaKs", package = "GeneFamilies")
    # data("orthologsKaKs", package = "GeneFamilies")
    data("fubar_results", package = "GeneFamilies")
    data("BUSTED_Results", package = "GeneFamilies")
    data("GeneGroups", package = "GeneFamilies")
    data("IprBasedEntropies", package = "GeneFamilies")
    data("ExpressionProfileDistances", package = "GeneFamilies")
    data("ExpressionProfileDistanceDistributions", package = "GeneFamilies")
    data("ExpressionProfileDistancesPerTissueDistributions", package = "GeneFamilies")
    data("geneCopyNumbers", package = "GeneFamilies")
    data("correlationExpressionCopyNumber", package = "GeneFamilies")
    data("mapManRootBinAnnotations", package = "GeneFamilies")
    data("rpkmExpressionProfiles", package = "GeneFamilies")
    data("rpkmNormalizedExpressionProfiles", package = "GeneFamilies")
    data("differentiallyExpressedGenes", package = "GeneFamilies")
    data("domainBasedSubfunctionalization", package = "GeneFamilies")
    data("tandemsDf", package = "GeneFamilies")
    data("DuplicatedKsStats", package = "GeneFamilies")
    # After loading the species' CDS define this constant:
    if (exists("aet.cds") && exists("aly.cds") && exists("ath.cds") && exists("bra.cds") && 
        exists("cru.cds") && exists("esa.cds") && exists("tpa.cds") && exists("chi.cds")) {
        assign("spec.gene.ids", list(aet = names(aet.cds), aly = names(aly.cds), 
            ath = names(ath.cds), bra = names(bra.cds), cru = names(cru.cds), esa = names(esa.cds), 
            tpa = names(tpa.cds), chi = names(chi.cds)), envir = globalenv())
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
