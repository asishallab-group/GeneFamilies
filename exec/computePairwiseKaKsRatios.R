require(MaizeGeneFamilies)
data("pairwiseSequenceSimilarities", package = "MaizeGeneFamilies")

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/computePairwiseKaKsRatios.R output_table.tsv [start-index stop-index]")

input.args <- commandArgs(trailingOnly = TRUE)

#' Which genes are not members of singleton families?
families.non.singleton.inds <- families.df[which(families.df$size > 1), "id"]
#' Compute Ka/Ks for them using their closest homologs:
genes.2.process <- unlist(families.lst[families.non.singleton.inds])
#' Do we work in parallel?
i.start <- if (length(input.args) > 1) as.numeric(input.args[[2]]) else 1
i.stop <- if (length(input.args) > 1) {
    if (as.numeric(input.args[[3]]) > length(genes.2.process)) 
        length(genes.2.process) else as.numeric(input.args[[3]])
} else length(genes.2.process)
#' Process batch:
genes.2.process <- genes.2.process[i.start:i.stop]

#' Write results into this data.frame:
pairwise.Ka.Ks <- data.frame(gene.a = c(), gene.b = c(), Ka = c(), Ks = c(), w = c(), 
    stringsAsFactors = FALSE)

#' Compute Ka-Ks-Ratios for closest homologous genes of non-singleton families:
while (length(genes.2.process) > 0) {
    gene.a <- genes.2.process[[1]]
    genes.2.process <- setdiff(genes.2.process, gene.a)
    if (!gene.a %in% pairwise.Ka.Ks$gene.a) {
        gene.b <- closestHomolog(gene.a)
        if (!is.na(gene.b) && "NA" != gene.b) {
            genes.2.process <- setdiff(genes.2.process, gene.b)
            genes.Ka.Ks <- computeKsPipeline(gene.a, gene.b)
            if (!is.null(genes.Ka.Ks) && !is.na(genes.Ka.Ks)) {
                i <- nrow(pairwise.Ka.Ks) + 1
                pairwise.Ka.Ks[i, "gene.a"] <- gene.a
                pairwise.Ka.Ks[i, "gene.b"] <- gene.b
                pairwise.Ka.Ks[i, "Ka"] <- genes.Ka.Ks[["Ka"]]
                pairwise.Ka.Ks[i, "Ks"] <- genes.Ka.Ks[["Ks"]]
                pairwise.Ka.Ks[i, "w"] <- genes.Ka.Ks[["w"]]
            }
        }
    }
}

#' Save results:
write.table(pairwise.Ka.Ks, input.args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)

message("DONE")
