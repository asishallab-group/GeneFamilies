require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/generateRpkmExpressionProfiles.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Transform into expression profiles:
genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))
rpkm.expr.profiles.df <- do.call("rbind", mclapply(genes, function(x) {
    y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
    x.df <- as.data.frame(t(setNames(y[, "expression"], y$tissue)), stringsAsFactors = FALSE)
    x.df$gene <- x
    x.df
}))
#' Add matching gene-names _with_ expression variants:
rpkm.expr.profiles.df$gene.exp.var <- as.character(unlist(mclapply(rpkm.expr.profiles.df$gene, 
    function(x) {
        x.exp.var <- names(all.cds)[grepl(paste("^", x, sep = ""), names(all.cds))]
        if (length(x.exp.var) == 1) {
            x.exp.var
        } else {
            warning("Found != 1 matching expression variants for '", x, 
                "': ", paste(x.exp.var, collapse = ", "), " !")
            NA
        }
    })))


#' Save results:
save(rpkm.expr.profiles.df, file = file.path(input.args[[1]], "data", "rpkmExpressionProfiles.RData"))

message("DONE")
