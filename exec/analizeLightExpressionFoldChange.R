require(GeneFamilies)

message("USAGE: Rscript path/2/exec/analizeLightExpression.R 1_UP_in_Ch_only.txt 2_UP_in_At_only.txt 3_DOWN_in_Ch_only.txt 4_DOWN_in_At_only.txt 5_UP_in_Ch+At_shared.txt 6_DOWN_in_Ch+At_shared.txt path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Ignore expression variants in MapMan-Bin Annotations:
mercator.df$IDENTIFIER.no.expr.var <- sub("\\.\\d+$", "", mercator.df$IDENTIFIER)


#' Read lists of differentially expressed (DE) genes:
chi.light.up <- tolower(readLines(input.args[[1]]))
chi.light.up.mapMan <- intersect(chi.light.up, mercator.df$IDENTIFIER.no.expr.var)
chi.light.up.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], chi.light.up)
ath.light.up <- tolower(readLines(input.args[[2]]))
ath.light.up.mapMan <- intersect(ath.light.up, mercator.df$IDENTIFIER.no.expr.var)
ath.light.up.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], ath.light.up)
chi.light.down <- tolower(readLines(input.args[[3]]))
chi.light.down.mapMan <- intersect(chi.light.down, mercator.df$IDENTIFIER.no.expr.var)
chi.light.down.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], chi.light.down)
ath.light.down <- tolower(readLines(input.args[[4]]))
ath.light.down.mapMan <- intersect(ath.light.down, mercator.df$IDENTIFIER.no.expr.var)
ath.light.down.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], ath.light.down)
ath.chi.light.up <- tolower(readLines(input.args[[5]]))
ath.chi.light.up.mapMan <- intersect(ath.chi.light.up, mercator.df$IDENTIFIER.no.expr.var)
ath.chi.light.up.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], ath.chi.light.up)
ath.chi.light.down <- tolower(readLines(input.args[[6]]))
ath.chi.light.down.mapMan <- intersect(ath.chi.light.down, mercator.df$IDENTIFIER.no.expr.var)
ath.chi.light.down.mapMan.not <- setdiff(mercator.df$IDENTIFIER.no.expr.var[grepl("^at", 
    mercator.df$IDENTIFIER.no.expr.var)], ath.chi.light.down)



#' Function to carry out exact Fischer test:
lightFischerTest <- function(genes.de, genes.not.de, univ.genes = union(genes.de, 
    genes.not.de)) {
    res.df <- Reduce(rbind, lapply(names(mercator.lst), function(mmBin) {
        x <- mercator.lst[[mmBin]]
        mmBin.light.genes <- intersect(tolower(sub("\\.\\d+$", "", x)), 
            univ.genes)
        if (length(mmBin.light.genes) > 0) {
            cont.tbl <- generateContingencyTable(genes.de, genes.not.de, 
                mmBin.light.genes, setdiff(univ.genes, mmBin.light.genes), 
                "DE", "MapManBin")
            p.val <- fisher.test(cont.tbl, alternative = "greater")$p.value
            data.frame(BINCODE = mmBin, p.value = p.val, stringsAsFactors = FALSE)
        } else NULL
    }))
    # Correct for multiple hypothesis testing:
    res.df$p.adjusted <- p.adjust(res.df$p.value, method = "BY")
    # Retain significant ones only:
    res.df.sign <- res.df[which(res.df$p.adjusted <= 0.05), ]
    # Add MapMan-Bin Names and Descriptions:
    res.df.sign$NAME <- sapply(res.df.sign$BINCODE, function(mmBin) mercator.df[which(mercator.df$BINCODE == 
        mmBin), "NAME"][[1]])
    res.df.sign$DESCRIPTION <- sapply(res.df.sign$BINCODE, function(mmBin) mercator.df[which(mercator.df$BINCODE == 
        mmBin), "DESCRIPTION"][[1]])
    # Done:
    res.df.sign
}


#' Identify overrepresented functions wihtin the above DE genes:
chi.light.up.fish <- lightFischerTest(chi.light.up.mapMan, chi.light.up.mapMan.not)
ath.light.up.fish <- lightFischerTest(ath.light.up.mapMan, ath.light.up.mapMan.not)
chi.light.down.fish <- lightFischerTest(chi.light.down.mapMan, chi.light.down.mapMan.not)
ath.light.down.fish <- lightFischerTest(ath.light.down.mapMan, ath.light.down.mapMan.not)
ath.chi.light.up.fish <- lightFischerTest(ath.chi.light.up.mapMan, ath.chi.light.up.mapMan.not)
ath.chi.light.down.fish <- lightFischerTest(ath.chi.light.down.mapMan, 
    ath.chi.light.down.mapMan.not)


#' Write output files:
write.table(chi.light.up.fish, file.path(input.args[[7]], "inst", "light_genes_de_foldChange_chi_up_overrep_MapManBins.txt"), 
    row.names = FALSE, sep = "\t")
write.table(ath.light.up.fish, file.path(input.args[[7]], "inst", "light_genes_de_foldChange_ath_up_overrep_MapManBins.txt"), 
    row.names = FALSE, sep = "\t")
write.table(chi.light.down.fish, file.path(input.args[[7]], "inst", "light_genes_de_foldChange_chi_down_overrep_MapManBins.txt"), 
    row.names = FALSE, sep = "\t")
write.table(ath.light.down.fish, file.path(input.args[[7]], "inst", "light_genes_de_foldChange_ath_down_overrep_MapManBins.txt"), 
    row.names = FALSE, sep = "\t")
write.table(ath.chi.light.up.fish, file.path(input.args[[7]], "inst", "light_genes_de_foldChange_ath_chi_up_overrep_MapManBins.txt"), 
    row.names = FALSE, sep = "\t")
write.table(ath.chi.light.down.fish, file.path(input.args[[7]], "inst", 
    "light_genes_de_foldChange_ath_chi_down_overrep_MapManBins.txt"), row.names = FALSE, sep = "\t")


message("DONE")
