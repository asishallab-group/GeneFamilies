require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/exec/assessDosageEffect.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

s.g.i <- spec.gene.ids[c("ath", "chi")]
gene.copy.number.df <- data.frame(Gene = unique(unlist(s.g.i)), tandem.no = 1, copy.no = 1, 
    stringsAsFactors = FALSE)
gene.copy.number.df$tandem.no <- unlist(mclapply(gene.copy.number.df$Gene, function(gene) {
    if (gene %in% tandems.genes) {
        gene.spec <- speciesForGeneId(gene, spec.genes = s.g.i)
        # Make use of the fact that each gene can only be found in a single tandem
        # cluster:
        gene.tand.clst <- tandems[which(tandems$Gene == gene), "Family"][[1]]
        length(intersect(s.g.i[[gene.spec]], tandems.lst[[gene.tand.clst]]))
    } else 1
}))
gene.copy.number.df$copy.no <- unlist(mclapply(gene.copy.number.df$Gene, function(gene) {
    if (!gene %in% orthologs.genes && gene %in% families.genes.df$Gene) {
        gene.spec <- speciesForGeneId(gene, spec.genes = s.g.i)
        gene.fam <- families.genes.df[[which(families.genes.df$Gene == gene), "Family"]]
        length(setdiff(intersect(s.g.i[[gene.spec]], families.lst[[gene.fam]]), orthologs.genes))
    } else 1
}))
gene.copy.number.df$Gene.no.expr.var <- sub("\\.\\d+$", "", gene.copy.number.df$Gene)


# Save results:
save(gene.copy.number.df, file = file.path(input.args[[1]], "data", "geneCopyNumbers.RData"))

message("DONE")
