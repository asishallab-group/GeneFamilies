require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/loadFubarResults.R path/2/families/working_directory path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Find all FUBAR result files and load them as tables:
fubar.fls <- system(paste("find", input.args[[1]], "-name 'cluster_*_ml_tree_no_node_labels.newick.fubar.csv' -type f"), 
    intern = TRUE)
names(fubar.fls) <- regmatches(fubar.fls, regexpr("cluster_\\d+", fubar.fls))

fubar.tbl <- Reduce(rbind, mclapply(names(fubar.fls), function(fam) {
    readFubarTable(fubar.fls[[fam]], fam)
}))


#' Identify families with decisive evidence for positive selection at at least
#' a single codon [https://en.wikipedia.org/wiki/Bayes_factor]:
fubar.fams.decisive.evidence <- sort(unique(fubar.tbl[which(fubar.tbl$BayesFactor > 
    100), "family"]))


#' Save results:
save(fubar.tbl, fubar.fams.decisive.evidence, file = file.path(input.args[[2]], "data", 
    "fubar_results.RData"))

message("DONE")
