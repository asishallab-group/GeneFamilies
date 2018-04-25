require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/loadMaizeExpressionData.R MBS_vs_MMS.tsv MBS_vs_SGsTC.tsv MMS_vs_SGlTC.tsv SGsTC_vs_SGlTC.tsv SeeTC_vs_MMS.tsv SeeTC_vs_SGlTC.tsv path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Read expression tables:
mmbs.v.mms <- read.table(input.args[[1]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)
mmbs.v.sgstc <- read.table(input.args[[2]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)
mms.v.sgltc <- read.table(input.args[[3]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)
sgstc.v.sgltc <- read.table(input.args[[4]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)
seetc.v.mms <- read.table(input.args[[5]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)
seetc.v.sgltc <- read.table(input.args[[6]], sep = "\t", stringsAsFactors = FALSE, 
    header = TRUE)


#' Adjust P-Values for multiple hypothesis testing:
mmbs.v.mms$padj <- p.adjust(mmbs.v.mms$pval, method = "fdr")
mmbs.v.sgstc$padj <- p.adjust(mmbs.v.sgstc$pval, method = "fdr")
mms.v.sgltc$padj <- p.adjust(mms.v.sgltc$pval, method = "fdr")
sgstc.v.sgltc$padj <- p.adjust(sgstc.v.sgltc$pval, method = "fdr")
seetc.v.mms$padj <- p.adjust(seetc.v.mms$pval, method = "fdr")
seetc.v.sgltc$padj <- p.adjust(seetc.v.sgltc$pval, method = "fdr")


#' Add log2-fold-changes:
log2Fc <- function(counts.tbl, untreated.name, treated.name) {
    unt.col <- colnames(counts.tbl)[grepl(paste("^", untreated.name, "\\d\\.est_counts$", 
        sep = ""), colnames(counts.tbl))]
    trt.col <- colnames(counts.tbl)[grepl(paste("^", treated.name, "\\d\\.est_counts$", 
        sep = ""), colnames(counts.tbl))]
    apply(counts.tbl, 1, function(x) {
        log2(mean(as.numeric(x[trt.col]), na.rm = TRUE)/mean(as.numeric(x[unt.col]), 
            na.rm = TRUE))
    })
}
mmbs.v.mms$log2.folg.change <- log2Fc( mmbs.v.mms, "MBS", "MMS" )
mmbs.v.sgstc$log2.folg.change <- log2Fc( mmbs.v.sgstc, "MBS", "SGsTC" )
mms.v.sgltc$log2.folg.change <- log2Fc( mms.v.sgltc, "MBS", "SGsTC" )
sgstc.v.sgltc$log2.folg.change <- log2Fc( sgstc.v.sgltc, "SGsTC", "SGlTC" )
seetc.v.mms$log2.folg.change <- log2Fc( seetc.v.mms, "MMS", "SeeTC" )
seetc.v.sgltc$log2.folg.change <- log2Fc( seetc.v.sgltc, "SGlTC", "SeeTC" )


#' Identify differentially expressed genes:
mmbs.v.mms.de.genes <- mmbs.v.mms[which(mmbs.v.mms$padj <= 0.05), "target_id"]
mmbs.v.sgstc.de.genes <- mmbs.v.sgstc[which(mmbs.v.sgstc$padj <= 0.05), 
    "target_id"]
mms.v.sgltc.de.genes <- mms.v.sgltc[which(mms.v.sgltc$padj <= 0.05), "target_id"]
sgstc.v.sgltc.de.genes <- sgstc.v.sgltc[which(sgstc.v.sgltc$padj <= 0.05), 
    "target_id"]
seetc.v.mms.de.genes <- seetc.v.mms[which(seetc.v.mms$padj <= 0.05), "target_id"]
seetc.v.sgltc.de.genes <- seetc.v.sgltc[which(seetc.v.sgltc$padj <= 0.05), 
    "target_id"]


#' Infer intersections:
venn.lst <- list(mmbs.v.mms = mmbs.v.mms.de.genes, mmbs.v.sgstc = mmbs.v.sgstc.de.genes, 
    mms.v.sgltc = mms.v.sgltc.de.genes, sgstc.v.sgltc = sgstc.v.sgltc.de.genes, 
    seetc.v.mms = seetc.v.mms.de.genes, seetc.v.sgltc = seetc.v.sgltc.de.genes)
maize.ustilago.de.venn <- venn(venn.lst[sapply(venn.lst, function(x) length(x) > 
    0)])
pdf(file.path(input.args[[7]], "inst", "maize_ustilago_de_venn_diagr.pdf"))
plot(maize.ustilago.de.venn)
dev.off()


#' Above Venn Diagram reveals only two sets are of interest. Venn with
#' house-keeping and cell-cycle genes:
venn.house.cell.cyc.lst <- list(house = toLowerCutTail(house), cell.cycle = toLowerCutTail(cell.cycle$V1), 
    mmbs.v.sgstc = mmbs.v.sgstc.de.genes, mms.v.sgltc = mms.v.sgltc.de.genes, 
    seetc.v.mms = seetc.v.sgltc.de.genes)
maize.ustilago.cell.house.de.venn <- venn(venn.house.cell.cyc.lst)
pdf(file.path(input.args[[7]], "inst", "maize_ustilago_cell_cycle_housekeeping_de_venn_diagr.pdf"))
plot(maize.ustilago.cell.house.de.venn)
dev.off()


#' Investigate intersections:
maize.gene.lst <- list(mmbs.v.sgstc = toLowerCutTail(mmbs.v.sgstc.de.genes), 
    mms.v.sgltc = toLowerCutTail(mms.v.sgltc.de.genes), seetc.v.mms = toLowerCutTail(seetc.v.sgltc.de.genes))
for (gene.class in unique(maize.gene.sets$V3)) {
    lst <- setNames(list(toLowerCutTail(maize.gene.sets[which(maize.gene.sets$V3 == 
        gene.class), 1])), gene.class)
    venn.lst <- venn(append(maize.gene.lst, lst))
    # Plot:
    pdf(file.path(input.args[[7]], "inst", paste("maize_ustilago_de_", 
        gene.class, "_venn_diagr.pdf", sep = "")))
    plot(venn.lst)
    dev.off()
    # Write Gene Lists:
    for (de.class in names(maize.gene.lst)) {
        de.inter.class <- intersect(unlist(lst), maize.gene.lst[[de.class]])
        fl.pth <- file.path(input.args[[7]], "inst", paste("intersect", 
            de.class, gene.class, ".txt", sep = "_"))
        writeLines(de.inter.class, fl.pth)
    }
}


#' Write lists of DE genes:
writeLines(sort(mmbs.v.mms.de.genes), file.path(input.args[[7]], "inst", 
    "mmbs.v.mms_de_genes.txt"))
writeLines(sort(mmbs.v.sgstc.de.genes), file.path(input.args[[7]], "inst", 
    "mmbs.v.sgstc_de_genes.txt"))
writeLines(sort(mms.v.sgltc.de.genes), file.path(input.args[[7]], "inst", 
    "mms.v.sgltc_de_genes.txt"))
writeLines(sort(sgstc.v.sgltc.de.genes), file.path(input.args[[7]], "inst", 
    "sgstc.v.sgltc_de_genes.txt"))
writeLines(sort(seetc.v.mms.de.genes), file.path(input.args[[7]], "inst", 
    "seetc.v.mms_de_genes.txt"))
writeLines(sort(seetc.v.sgltc.de.genes), file.path(input.args[[7]], "inst", 
    "seetc.v.sgltc_de_genes.txt"))


#' Save results:
save(mmbs.v.mms, mmbs.v.sgstc, mms.v.sgltc, sgstc.v.sgltc, seetc.v.mms, 
    seetc.v.sgltc, mmbs.v.mms.de.genes, mmbs.v.sgstc.de.genes, mms.v.sgltc.de.genes, 
    sgstc.v.sgltc.de.genes, seetc.v.mms.de.genes, seetc.v.sgltc.de.genes, 
    maize.ustilago.de.venn, file = file.path(input.args[[7]], "data", "maizeUstilagoDeGenes.RData"))


message("DONE")
