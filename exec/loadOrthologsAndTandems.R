require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadOrthologsAndParalogs.R path/2/MaizeGeneFamilies/data all_vs_all_tabular_blast_out.txt eight_brassicaceae_orthologs.txt eight_brassicaceae_tandems.txt")
message("EXPEXTED FORMATS:")
message("- all_vs_all_tabular_blast_out.txt is the result of calling BLAT or BLAST on all coding sequences in an 'all vs all' approach. The result file is expected to be in the tabular Blast output format ")
message("- eight_brassicaceae_orthologs.txt is a TAB separated table holding orthologous gene clusters. The table is expected to have the following header:\naly\tath\tcru\tchi\taet\tbra\tesa\ttpa\n")
message("- eight_brassicaceae_tandems.txt is a TAb separated table with header 'Family Gene' and rows in the form of 'tandem_cluster_1 CARHR1234.1'.")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load data:
all.vs.all.sim <- fread(input.args[[2]], data.table = FALSE, header = FALSE, stringsAsFactors = FALSE, 
    sep = "\t", na.strings = "", colClasses = c(rep("character", 2), rep("numeric", 
        10)))
orthologs <- read.table(input.args[[3]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 8))
tandems <- read.table(input.args[[4]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 2))

#' Save loaded data:
save(orthologs, tandems, file = file.path(input.args[[1]], "orthologsTandems.RData"))
save(all.vs.all.sim, file = file.path(input.args[[1]], "pairwiseSequenceSimilarities.RData"))

message("DONE")
