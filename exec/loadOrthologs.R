require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadOrthologs.R path/2/GeneFamilies/data pairwiseReciprBestHitsConcat.txt")
message("To know how to create the input file see this package's Vignette 'GeneFamilies'")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load data:
orthologs <- orthologsFromPairwiseReciprocalBestHits(read.table(input.args[[2]], 
    sep = "\t", comment.char = "", quote = "", na.strings = "", colClasses = rep("character", 
        2), stringsAsFactors = FALSE))

#' Save loaded data:
save(orthologs, file = file.path(input.args[[1]], "orthologs.RData"))

message("DONE")
