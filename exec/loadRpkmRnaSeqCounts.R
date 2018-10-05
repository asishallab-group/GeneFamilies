require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/loadRpkmRnaSeqCounts.R RPKM_counts_table.tsv path/2/MaizeGeneFamilies")
message("Note, that the table in argument file 'RPKM_counts_table.tsv' is expected to be TAB-Delimited and have a header line:\n", 
    "id tissue rank expression variance")

input.args <- commandArgs(trailingOnly = TRUE)

#' Read RPKM normalized counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    comment.char = "", quote = "", na.strings = "", colClasses = c(rep("character", 
        3), rep("numeric", "2")))
#' Transform into expression profiles:
genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))
rna.seq.exp.profils <- do.call("rbind", mclapply(genes, function(x) {
    y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
    x.df <- as.data.frame(t(setNames(y[, "expression"]/sum(y[, "expression"], na.rm = TRUE), 
        y$tissue)), stringsAsFactors = FALSE)
    x.df$gene <- x
    x.df
}))
#' Add matching gene-names _with_ expression variants:
rna.seq.exp.profils$gene.exp.var <- as.character(unlist(mclapply(rna.seq.exp.profils$gene, 
    function(x) {
        x.exp.var <- names(all.cds)[grepl(paste("^", x, sep = ""), names(all.cds))]
        if (length(x.exp.var) == 1) {
            x.exp.var
        } else {
            warning("Found != 1 matching expression variants for '", x, "': ", paste(x.exp.var, 
                collapse = ", "), " !")
            NA
        }
    })))

#' Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(input.args[[2]], 
    "data", "RNA_Seq_RPKM_and_profiles.RData"))
write.table(rna.seq.exp.profils, file.path(input.args[[2]], "inst", "RNA_Seq_RPKM_and_profiles.tsv"), 
    sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
