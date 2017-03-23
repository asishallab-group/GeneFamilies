require(GeneFamilies)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/computeAnnotationBasedShannonEntropies.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute annotation based Shannon Entropy for the following gene groups.
#' - Tandems:
tandems.ipr.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(tandems.lst), 
    tandems.lst)
#' - Orthologs:
orthologs.ipr.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(orthologs.lst), 
    orthologs.lst)
#' - Gene Families:
families.ipr.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(families.lst), 
    families.lst)


#' Make a distinction between Orthologs and Non-Orthologs within each gene
#' group.
#' - Tandems:
tandems.ipr.entropies.orth.dist <- computeOrthologSpecificAnnotationDiversityPerGroup(tandems, 
    group.col = "Family", gene.col = "Gene")
#' - Gene Families:
families.ipr.entropies.orth.dist <- computeOrthologSpecificAnnotationDiversityPerGroup(families.genes.df, 
    group.col = "Family", gene.col = "Gene")


#' Save results:
save(tandems.ipr.entropies, orthologs.ipr.entropies, families.ipr.entropies, tandems.ipr.entropies.orth.dist, 
    families.ipr.entropies.orth.dist, file = file.path(input.args[[1]], "data", "IprBasedEntropies.RData"))


message("DONE")
