require(MaizeGeneFamilies)


message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/analyzeMaizeDeGenes.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' ********
#' MapMan *
#' ********

#' Prepare exact Fischer tests:
maize.w.mapMan <- intersect(toLowerCutTail(names(maize.aas)), toLowerCutTail(maize.mapMan$IDENTIFIER))
maize.genes.expr <- Reduce(union, list(mmbs.v.mms$target_id, mmbs.v.sgstc$target_id, 
    mms.v.sgltc$target_id, sgstc.v.sgltc$target_id, seetc.v.mms$target_id, 
    seetc.v.sgltc$target_id))
maize.expr.w.mapMan <- intersect(maize.w.mapMan, toLowerCutTail(maize.genes.expr))


#' Venn Diagramm has revealed, that there are only four sets of Diff Expr Genes
#' (DEG):
sgltc.inter.sgstc.de.genes.w.mmb <- intersect(intersect(toLowerCutTail(mms.v.sgltc.de.genes), 
    toLowerCutTail(mmbs.v.sgstc.de.genes)), maize.expr.w.mapMan)
mms.v.sgltc.de.genes.uniq.w.mmb <- setdiff(intersect(toLowerCutTail(mms.v.sgltc.de.genes), 
    maize.expr.w.mapMan), sgltc.inter.sgstc.de.genes.w.mmb)
mmbs.v.sgstc.de.genes.uniq.w.mmb <- setdiff(intersect(toLowerCutTail(mmbs.v.sgstc.de.genes), 
    maize.expr.w.mapMan), sgltc.inter.sgstc.de.genes.w.mmb)
sgltc.union.sgstc.de.genes.w.mmb <- intersect(union(toLowerCutTail(mms.v.sgltc.de.genes), 
    toLowerCutTail(mmbs.v.sgstc.de.genes)), maize.expr.w.mapMan)


#' exact Fischer test to find over represented MapMan-Bin annotations:
maize.mapMan$IDENTIFIER.san <- toLowerCutTail(maize.mapMan$IDENTIFIER)
mmBins <- unique(maize.mapMan$BINCODE)


#' Function to execute Fischer test:
maizeFischerTest <- function(genes.de, genes.not.de, univ.genes = union(genes.de, 
    genes.not.de)) {
    res.df <- Reduce(rbind, lapply(mmBins, function(mmBin) {
        mmBin.expr.genes <- intersect(maize.mapMan[which(maize.mapMan$BINCODE == 
            mmBin), "IDENTIFIER.san"], maize.expr.w.mapMan)
        if (length(mmBin.expr.genes) > 0) {
            cont.tbl <- generateContingencyTable(genes.de, genes.not.de, 
                mmBin.expr.genes, setdiff(univ.genes, mmBin.expr.genes), 
                "DE", "MapManBin")
            p.val <- fisher.test(cont.tbl, alternative = "greater")$p.value
            data.frame(BINCODE = mmBin, p.value = p.val, stringsAsFactors = FALSE)
        } else NULL
    }))
    if (!is.null(res.df)) {
        # Correct for multiple hypothesis testing:
        res.df$p.adjusted <- p.adjust(res.df$p.value, method = "fdr")
        # Retain significant ones only:
        res.df.sign <- res.df[which(res.df$p.adjusted <= 0.05), ]
        if (nrow(res.df.sign) > 0) {
            # Add MapMan-Bin Names and Descriptions:
            res.df.sign$NAME <- sapply(res.df.sign$BINCODE, function(mmBin) maize.mapMan[which(maize.mapMan$BINCODE == 
                mmBin), "NAME"][[1]])
            res.df.sign$DESCRIPTION <- sapply(res.df.sign$BINCODE, function(mmBin) maize.mapMan[which(maize.mapMan$BINCODE == 
                mmBin), "DESCRIPTION"][[1]])
            # Done:
            res.df.sign
        } else {
            warning("No significant P-Values found.")
            res.df
        }
    } else NULL
}


#' Fisher tests:
mms.v.sgltc.de.genes.uniq.fish <- maizeFischerTest(mms.v.sgltc.de.genes.uniq.w.mmb, 
    setdiff(maize.expr.w.mapMan, mms.v.sgltc.de.genes.uniq.w.mmb))
sgltc.inter.sgstc.de.genes.fish <- maizeFischerTest(sgltc.inter.sgstc.de.genes.w.mmb, 
    setdiff(maize.expr.w.mapMan, sgltc.inter.sgstc.de.genes.w.mmb))
sgltc.union.sgstc.de.genes.fish <- maizeFischerTest(sgltc.union.sgstc.de.genes.w.mmb, 
    setdiff(maize.expr.w.mapMan, sgltc.union.sgstc.de.genes.w.mmb))
mmbs.v.sgstc.de.genes.uniq.fish <- maizeFischerTest(mmbs.v.sgstc.de.genes.uniq.w.mmb, 
    setdiff(maize.expr.w.mapMan, mmbs.v.sgstc.de.genes.uniq.w.mmb))


#' Analysis showed that only 'sgltc.inter.sgstc.de.genes.fish' and
#' 'sgltc.union.sgstc.de.genes.fish' contain significant results:
overrep.mapManBins <- unique(union(sgltc.inter.sgstc.de.genes.fish$BINCODE, 
    sgltc.union.sgstc.de.genes.fish$BINCODE))
overrep.mapManBins.genes.df <- Reduce(rbind, lapply(overrep.mapManBins, 
    function(mmBin) {
        binGenes.df <- maize.mapMan[(grepl(paste("^", mmBin, "[0-9.]*", 
            sep = ""), maize.mapMan$BINCODE) & maize.mapMan$TYPE == "TRUE"), 
            ]
        mmbs.v.sgstc.de.genes.mmBin <- binGenes.df[which(binGenes.df$IDENTIFIER.san %in% 
            toLowerCutTail(mmbs.v.sgstc.de.genes)), c("IDENTIFIER.san", 
            "BINCODE", "NAME")]
        mmbs.v.sgstc.de.genes.mmBin$DEG.set <- "mmbs.v.sgstc"
        mms.v.sgltc.de.genes.mmBin <- binGenes.df[which(binGenes.df$IDENTIFIER.san %in% 
            toLowerCutTail(mms.v.sgltc.de.genes)), c("IDENTIFIER.san", 
            "BINCODE", "NAME")]
        mms.v.sgltc.de.genes.mmBin$DEG.set <- "mms.v.sgltc"
        # Generate a Venn Diagramm:
        mmBin.venn.lst <- list(mmbs.v.sgstc = unique(mmbs.v.sgstc.de.genes.mmBin$IDENTIFIER.san), 
            mms.v.sgltc = unique(mms.v.sgltc.de.genes.mmBin$IDENTIFIER.san))
        mmBin.venn <- venn(mmBin.venn.lst)
        pdf(file.path(input.args[[1]], "inst", paste("mapManBin_DEG_sgst_sgltc_", 
            mmBin, "_venn.pdf", sep = "")))
        plot(mmBin.venn)
        dev.off()
        # Generate result data.frame:
        res.df <- rbind(mmbs.v.sgstc.de.genes.mmBin, mms.v.sgltc.de.genes.mmBin)
        res.df$overrep.mapManBin <- mmBin
        unique(res.df)
    }))


#' Write table
write.table(overrep.mapManBins.genes.df, file.path(input.args[[1]], "inst", 
    "overrep_mapManBins_DEG_tbl.txt"), sep = "\t", row.names = FALSE)


#' Save results:
save(mms.v.sgltc.de.genes.uniq.fish, sgltc.inter.sgstc.de.genes.fish, mmbs.v.sgstc.de.genes.uniq.fish, 
    overrep.mapManBins.genes.df, file = file.path(input.args[[1]], "data", 
        "maize_de_mapManBin_overrep.RData"))


message("DONE")
