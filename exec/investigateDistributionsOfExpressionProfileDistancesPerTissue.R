require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfExpressionProfileDistancesPerTissue.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute statistics about Expression Profile Distances within different
#' classes of gene groups. In this, consider orthologs to be the 'background'.
gene.classes <- list(Orthologs = sub("\\.\\d+$", "", orthologs.genes), `Non-Orthologs` = sub("\\.\\d+$", 
    "", setdiff(names(c(ath.cds, chi.cds)), orthologs.genes)), all = sub("\\.\\d+$", 
    "", names(c(ath.cds, chi.cds))))
#' - Tandems:
tandems.exp.prof.dists.tissue.orth.dist <- Reduce(rbind, mclapply(names(tandems.exp.prof.dists.tissue), 
    function(x) expressionProfileDistStatsPerTissues(tandems.exp.prof.dists.tissue[[x]], 
        gene.classes, x)))
#' - Families:
families.exp.prof.dists.tissue.orth.dist <- Reduce(rbind, mclapply(names(families.exp.prof.dists.tissue), 
    function(x) expressionProfileDistStatsPerTissues(families.exp.prof.dists.tissue[[x]], 
        gene.classes, x)))


#' Save results:
save(tandems.exp.prof.dists.tissue.orth.dist, families.exp.prof.dists.tissue.orth.dist, 
    file = file.path(input.args[[1]], "data", "ExpressionProfileDistancesPerTissueDistributions.RData"))


message("DONE")
