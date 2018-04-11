require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/plotDistributionsOfExpressionProfileDistancesPerTissue.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


tissues <- c("seedling", "cotyledon", "developing leaf", "flower stage 9", "flower stage 16")
gene.classes <- c("Orthologs", "Non-Orthologs")
#' Plot median expression profile distances per tissue and gene class, i.e.
#' orthologs and non-orthologs.
#' - Tandems:
plot.df <- tandems.exp.prof.dists.tissue.orth.dist[with(tandems.exp.prof.dists.tissue.orth.dist, 
    which(Statistic == "median" & Gene.Class != "all")), ]
plotTissueSpecificExpressionProfileDistanceDistribs(plot.df, tissues, gene.classes, 
    file.path(input.args[[1]], "inst", "tandemsTissueSpecificExpressionProfileDistancesBoxplot.pdf"), 
    "Tandems")
#' - Expanded Families:
plot.df <- families.exp.prof.dists.tissue.orth.dist[with(families.exp.prof.dists.tissue.orth.dist, 
    which(Statistic == "median" & Gene.Class != "all" & Family %in% families.exp)), 
    ]
plotTissueSpecificExpressionProfileDistanceDistribs(plot.df, tissues, gene.classes, 
    file.path(input.args[[1]], "inst", "expandedFamsTissueSpecificExpressionProfileDistancesBoxplot.pdf"), 
    "Expanded Families")
#' - Positively Selected Families:
plot.df <- families.exp.prof.dists.tissue.orth.dist[with(families.exp.prof.dists.tissue.orth.dist, 
    which(Statistic == "median" & Gene.Class != "all" & Family %in% fams.pos.sel)), 
    ]
plotTissueSpecificExpressionProfileDistanceDistribs(plot.df, tissues, gene.classes, 
    file.path(input.args[[1]], "inst", "posSelFamsTissueSpecificExpressionProfileDistancesBoxplot.pdf"), 
    "Positively Selected Families")


#' T-Tests in order to investigate significance and extend of differences
#' between the distributions of per gene-group and tissue median expression
#' profile distances.
#' - Tandems:
tandems.exp.prof.dists.tissue.orth.dist.t.tests <- testHourglassModel(tandems.exp.prof.dists.tissue.orth.dist)
#' - Expanded Families:
families.expanded.exp.prof.dists.tissue.orth.dist.t.tests <- testHourglassModel(families.exp.prof.dists.tissue.orth.dist[which(families.exp.prof.dists.tissue.orth.dist$Family %in% 
    families.exp), ])
#' - Positively Selected Families:
families.psel.exp.prof.dists.tissue.orth.dist.t.tests <- testHourglassModel(families.exp.prof.dists.tissue.orth.dist[which(families.exp.prof.dists.tissue.orth.dist$Family %in% 
    fams.pos.sel), ])


#' Save results:
write.table(tandems.exp.prof.dists.tissue.orth.dist.t.tests, file.path(input.args[[1]], 
    "inst", "hourglassModel_T_Tests_tandems.tsv"), sep = "\t", row.names = FALSE, 
    quote = FALSE)
write.table(families.expanded.exp.prof.dists.tissue.orth.dist.t.tests, file.path(input.args[[1]], 
    "inst", "hourglassModel_T_Tests_expandedFamilies.tsv"), sep = "\t", row.names = FALSE, 
    quote = FALSE)
write.table(families.psel.exp.prof.dists.tissue.orth.dist.t.tests, file.path(input.args[[1]], 
    "inst", "hourglassModel_T_Tests_positivelySelectedFamilies.tsv"), sep = "\t", 
    row.names = FALSE, quote = FALSE)


message("DONE")
