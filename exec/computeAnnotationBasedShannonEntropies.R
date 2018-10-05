require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/computeAnnotationBasedShannonEntropies.R path/2/MaizeGeneFamilies")

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


#' In the following use normalized Shannon Entropy see
#' en.wikipedia.org/wiki/Entropy_(information_theory)#Efficiency
options(GeneFamilies.entropy.function = MaizeGeneFamilies::shannonEntropy)

#' Compute NORMALIZED annotation based Shannon Entropy for the following gene
#' groups.  - Tandems:
tandems.ipr.norm.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(tandems.lst), 
    tandems.lst)
#' - Orthologs:
orthologs.ipr.norm.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(orthologs.lst), 
    orthologs.lst)
#' - Gene Families:
families.ipr.norm.entropies <- geneGroupsAnnotationBasedShannonEntropy(names(families.lst), 
    families.lst)


#' Make a NORMALIZED distinction between Orthologs and Non-Orthologs within
#' each gene group.
#' - Tandems:
tandems.ipr.norm.entropies.orth.dist <- computeOrthologSpecificAnnotationDiversityPerGroup(tandems, 
    group.col = "Family", gene.col = "Gene")
#' - Gene Families:
families.ipr.norm.entropies.orth.dist <- computeOrthologSpecificAnnotationDiversityPerGroup(families.genes.df, 
    group.col = "Family", gene.col = "Gene")


#' Save results:
save(tandems.ipr.entropies, orthologs.ipr.entropies, families.ipr.entropies, 
    tandems.ipr.entropies.orth.dist, families.ipr.entropies.orth.dist, 
    tandems.ipr.norm.entropies, orthologs.ipr.norm.entropies, families.ipr.norm.entropies, 
    tandems.ipr.norm.entropies.orth.dist, families.ipr.norm.entropies.orth.dist, 
    file = file.path(input.args[[1]], "data", "IprBasedEntropies.RData"))


message("DONE")
