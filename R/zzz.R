.onLoad <- function(libname = find.package("GeneFamilies"), pkgname = "GeneFamilies") {
    data("codingSequences", package = "GeneFamilies")
    data("families", package = "GeneFamilies")
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

