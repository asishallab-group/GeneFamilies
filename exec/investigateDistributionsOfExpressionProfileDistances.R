require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/investigateDistributionsOfExpressionProfileDistances.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

#' Compute statistics about Expression Profile Distances within different
#' classes of gene groups. In this, consider orthologs to be the 'background'.
gene.classes <- list(Orthologs = sub("\\.\\d+$", "", orthologs.genes), `Non-Orthologs` = sub("\\.\\d+$", 
    "", setdiff(names(all.cds), orthologs.genes)))
#' - Tandems:
tandems.exp.prof.dists.orth.dist <- setNames(mclapply(names(tandems.exp.prof.dists), 
    function(x) classSpecificExpressionProfileDists(tandems.exp.prof.dists[[x]], 
        gene.classes = gene.classes, fam.name = x)), names(tandems.exp.prof.dists))
tandems.exp.prof.dists.orth.dist.df <- Reduce(rbind, mclapply(names(tandems.exp.prof.dists.orth.dist), 
    function(x) tandems.exp.prof.dists.orth.dist[[x]][["stats"]]))
#' - Families:
families.exp.prof.dists.orth.dist <- setNames(mclapply(names(families.exp.prof.dists), 
    function(x) classSpecificExpressionProfileDists(families.exp.prof.dists[[x]], 
        gene.classes = gene.classes, fam.name = x)), names(families.exp.prof.dists))
families.exp.prof.dists.orth.dist.df <- Reduce(rbind, mclapply(names(families.exp.prof.dists.orth.dist), 
    function(x) families.exp.prof.dists.orth.dist[[x]][["stats"]]))
#' - Orthologs:
#' (Of course no distinction can be made between Non-Orthologs and Orthologs
#' within _Orthologs_ themselves, but none the less we need the statistics for
#' each cluster of orthologs.)
orthologs.exp.prof.dists.stats <- setNames(mclapply(names(orthologs.exp.prof.dists), 
    function(x) classSpecificExpressionProfileDists(orthologs.exp.prof.dists[[x]], 
        gene.classes = gene.classes["Orthologs"], fam.name = x)), names(orthologs.exp.prof.dists))
orthologs.exp.prof.dists.stats.df <- Reduce(rbind, mclapply(names(orthologs.exp.prof.dists.stats), 
    function(x) orthologs.exp.prof.dists.stats[[x]][["stats"]]))


#' Save results:
save(tandems.exp.prof.dists.orth.dist, tandems.exp.prof.dists.orth.dist.df, families.exp.prof.dists.orth.dist, 
    families.exp.prof.dists.orth.dist.df, orthologs.exp.prof.dists.stats.df, file = file.path(input.args[[1]], 
        "data", "ExpressionProfileDistanceDistributions.RData"))


message("DONE")
