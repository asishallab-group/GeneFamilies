require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/loadMaizeGenomeAndMapManAnnos.R maizeProteins_AA.fa mapManBinAnnos.txt path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Load Maize Proteins:
maize.aas <- read.fasta(input.args[[1]], seqtype = "AA", strip.desc = TRUE, 
    as.string = TRUE)


#' Load MapMan-Bin Annotations:
maize.mapMan <- fread(input.args[[2]], data.table = FALSE, stringsAsFactors = FALSE, 
    header = TRUE, sep = "\t")
for (i in colnames(maize.mapMan)) {
    maize.mapMan[, i] <- sub("^'", "", sub("'$", "", maize.mapMan[, i]))
}


#' Save results:
save( maize.aas, maize.mapMan, file=file.path( input.args[[3]], 'data', 'maizeGenome.RData' ) )


message("DONE")
