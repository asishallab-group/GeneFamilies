require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/computeExpressionProfileDistances.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

#' Compute for each group pairwise gene expression profile distances:
#' - Tandems:
tandems.exp.prof.dists <- mclapply(tandems.lst, expressionProfilesDists)
tandems.exp.prof.dists.tissue <- mclapply(tandems.lst, expressionProfilesDists, per.tissue = TRUE)
#' - Orthologs:
orthologs.exp.prof.dists <- mclapply(orthologs.lst, expressionProfilesDists)
orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, expressionProfilesDists, 
    per.tissue = TRUE)
#' - Gene Families:
non.singleton.fams <- families.df$id[which(families.df$size > 1)]
families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists)
families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists, 
    per.tissue = TRUE)


#' Save results:
save(tandems.exp.prof.dists, tandems.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
    orthologs.exp.prof.dists.tissue, families.exp.prof.dists, families.exp.prof.dists.tissue, 
    file = file.path(input.args[[1]], "data", "ExpressionProfileDistances.RData"))


message("DONE")
