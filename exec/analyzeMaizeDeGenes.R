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
mmbs.v.sgstc.de.genes.uniq.fish <- maizeFischerTest(mmbs.v.sgstc.de.genes.uniq.w.mmb, 
    setdiff(maize.expr.w.mapMan, mmbs.v.sgstc.de.genes.uniq.w.mmb))


#' ******
#' PFam *
#' ******
maize.w.pfam <- intersect(toLowerCutTail(names(maize.aas)), maize.pfam$gene.san.id)
maize.expr.w.pfam <- intersect(maize.w.pfam, toLowerCutTail(maize.genes.expr))
maize.pfam.compound <- data.frame(gene = maize.w.pfam, compound.Pfam = sapply(maize.w.pfam, 
    function(x) {
        paste(sort(unique(maize.pfam[which(maize.pfam$gene.san.id == x), 
            "Pfam"])), collapse = ",")
    }), stringsAsFactors = FALSE)


#' Venn Diagramm has revealed, that there are only four sets of Diff Expr Genes
#' (DEG):
sgltc.inter.sgstc.de.genes.w.pfam <- intersect(intersect(toLowerCutTail(mms.v.sgltc.de.genes), 
    toLowerCutTail(mmbs.v.sgstc.de.genes)), maize.expr.w.pfam)
mms.v.sgltc.de.genes.uniq.w.pfam <- setdiff(intersect(toLowerCutTail(mms.v.sgltc.de.genes), 
    maize.expr.w.pfam), sgltc.inter.sgstc.de.genes.w.pfam)
mmbs.v.sgstc.de.genes.uniq.w.pfam <- setdiff(intersect(toLowerCutTail(mmbs.v.sgstc.de.genes), 
    maize.expr.w.pfam), sgltc.inter.sgstc.de.genes.w.pfam)


#' Function to execute Fischer test:
maizeFischerTestGen <- function(genes.de, genes.not.de, univ.genes = union(genes.de, 
    genes.not.de), anno.tbl, at.gene.col = 3, at.anno.col = 2) {
    annos <- unique(anno.tbl[anno.tbl[, at.gene.col] %in% univ.genes, at.anno.col])
    res.df <- Reduce(rbind, lapply(annos, function(anno.i) {
        expr.genes.w.anno <- intersect(anno.tbl[which(anno.tbl[, at.anno.col] == 
            anno.i), at.gene.col], genes.de)
        if (length(expr.genes.w.anno) > 0) {
            cont.tbl <- generateContingencyTable(genes.de, genes.not.de, 
                expr.genes.w.anno, setdiff(univ.genes, expr.genes.w.anno), 
                "DE", "Annotation")
            p.val <- fisher.test(cont.tbl, alternative = "greater")$p.value
            data.frame(Annotation = anno.i, p.value = p.val, stringsAsFactors = FALSE)
        } else NULL
    }))
    if (!is.null(res.df)) {
        # Correct for multiple hypothesis testing:
        res.df$p.adjusted <- p.adjust(res.df$p.value, method = "fdr")
        # Retain significant ones only:
        res.df.sign <- res.df[which(res.df$p.adjusted <= 0.05), ]
        if (nrow(res.df.sign) > 0) {
            res.df.sign
        } else {
            warning("No significant P-Values found.")
            res.df
        }
    } else NULL
}


#' Function to lookup InterPro NAME for given PFam IDs:
iprNamesForPFamIds <- function(pfam.ids) {
    Reduce(rbind, lapply(pfam.ids, function(pfam.id) {
        ipr.id <- pfam.2.ipr[which(pfam.2.ipr$V1 == pfam.id), 2]
        if (!is.null(ipr.id) && !is.na(ipr.id) && length(ipr.id) > 0 && 
            ipr.id %in% names(ipr.db)) {
            ipr.entry <- ipr.db[[ipr.id]]
            data.frame(IPR.ID = ipr.id, IPR.NAME = ipr.entry[["NAME"]], 
                stringsAsFactors = FALSE)
        } else {
            data.frame(IPR.ID = NA, IPR.NAME = NA, stringsAsFactors = FALSE)
        }
    }))
}


#' Execute Fischer tests to find enriched PFams among DE genes:
mms.v.sgltc.de.genes.uniq.fish.pfam <- maizeFischerTestGen(mms.v.sgltc.de.genes.uniq.w.pfam, 
    setdiff(maize.expr.w.pfam, mms.v.sgltc.de.genes.uniq.w.pfam), anno.tbl = maize.pfam)
mms.v.sgltc.de.genes.uniq.fish.pfam <- cbind(mms.v.sgltc.de.genes.uniq.fish.pfam, 
    iprNamesForPFamIds(mms.v.sgltc.de.genes.uniq.fish.pfam$Annotation))

sgltc.inter.sgstc.de.genes.fish.pfam <- maizeFischerTestGen(sgltc.inter.sgstc.de.genes.w.pfam, 
    setdiff(maize.expr.w.pfam, sgltc.inter.sgstc.de.genes.w.pfam), anno.tbl = maize.pfam)
sgltc.inter.sgstc.de.genes.fish.pfam <- cbind(sgltc.inter.sgstc.de.genes.fish.pfam, 
    iprNamesForPFamIds(sgltc.inter.sgstc.de.genes.fish.pfam$Annotation))

mmbs.v.sgstc.de.genes.uniq.fish.pfam <- maizeFischerTestGen(mmbs.v.sgstc.de.genes.uniq.w.pfam, 
    setdiff(maize.expr.w.pfam, mmbs.v.sgstc.de.genes.uniq.w.pfam), anno.tbl = maize.pfam)
mmbs.v.sgstc.de.genes.uniq.fish.pfam <- cbind(mmbs.v.sgstc.de.genes.uniq.fish.pfam, 
    iprNamesForPFamIds(mmbs.v.sgstc.de.genes.uniq.fish.pfam$Annotation))



message("DONE")
