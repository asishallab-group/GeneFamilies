require(MaizeGeneFamilies)
data("pairwiseSequenceSimilarities", package = "MaizeGeneFamilies")

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/computeOrthologsKaKsRatios.R orthologs_Ka_Ks_Batch.tsv [ start-index stop-index ]")
message("NOTE: Use parameters 'start-index' and 'stop-index' to execute this computation in parallel.")

input.args <- commandArgs(trailingOnly = TRUE)

#' Which genes are not members of singleton families?
families.non.singleton.inds <- families.df[which(families.df$size > 1), "id"]

#' And which of them are Orthologs and do still lack pairwise Ka/Ks-values:
all.orths <- as.character(unlist(orthologs))
genes.2.process <- setdiff(all.orths, unlist(pairwise.Ka.Ks[, c("gene.a", "gene.b")]))
#' Do we work in parallel?
i.start <- if (length(input.args) > 1) as.numeric(input.args[[2]]) else 1
i.stop <- if (length(input.args) > 1) {
    if (as.numeric(input.args[[3]]) > length(genes.2.process)) 
        length(genes.2.process) else as.numeric(input.args[[3]])
} else length(genes.2.process)
#' Process batch:
genes.2.process <- genes.2.process[i.start:i.stop]
#' Reduced sequence similarity result table for the above orthologs:
orths.seq.sim <- all.vs.all.sim[which(all.vs.all.sim$V1 %in% all.orths | all.vs.all.sim$V2 %in% 
    all.orths), ]

#' Compute Ka-Ks-Ratios for closest orthologous genes:
orths.Ka.Ks <- data.frame(gene.a = c(), gene.b = c(), Ka = c(), Ks = c(), w = c(), 
    stringsAsFactors = FALSE)

while (length(genes.2.process) > 0) {
    gene.a <- genes.2.process[[1]]
    genes.2.process <- setdiff(genes.2.process, gene.a)
    if (!gene.a %in% orths.Ka.Ks$gene.a) {
        gene.b <- closestHomolog(gene.a, seq.sim.tbl = orths.seq.sim)
        if (!is.na(gene.b) && "NA" != gene.b) {
            genes.2.process <- setdiff(genes.2.process, gene.b)
            genes.Ka.Ks <- computeKsPipeline(gene.a, gene.b)
            if (!is.null(genes.Ka.Ks) && !is.na(genes.Ka.Ks)) {
                i <- nrow(orths.Ka.Ks) + 1
                orths.Ka.Ks[i, "gene.a"] <- gene.a
                orths.Ka.Ks[i, "gene.b"] <- gene.b
                orths.Ka.Ks[i, "Ka"] <- genes.Ka.Ks[["Ka"]]
                orths.Ka.Ks[i, "Ks"] <- genes.Ka.Ks[["Ks"]]
                orths.Ka.Ks[i, "w"] <- genes.Ka.Ks[["w"]]
            }
        }
    }
}

#' Save results:
write.table(orths.Ka.Ks, input.args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)

message("DONE")

repl <- function(x, y = aly.nm) {
    for (x.col in c("gene.a", "gene.b")) {
        i <- which(x[, x.col] %in% y[, 2])
        x[i, x.col] <- as.character(mclapply(x[i, x.col], function(z) y[which(y[, 
            2] == z), 1]))
    }
    x
}
