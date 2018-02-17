require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/exec/assessDosageEffect.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

s.g.i <- spec.gene.ids[c("ath", "chi")]
message("Computing correlation between gene copy number and their RPKM values")
message("Species are: ", paste(names(s.g.i), collapse = ", "))

#'######################
#' Regression analyses #
#'######################

# Families that have paralogs in one or both of the investigated
# species (see s.g.i):
genes.paralogs <- gene.copy.number.df[which(gene.copy.number.df$copy.no > 
    1), "Gene"]
fams.paralogs <- unique(families.genes.df[which(families.genes.df$Gene %in% 
    genes.paralogs), "Family"])
# RPKM counts for genes within the above families:
paralogs.no.expr.var <- sub("\\.\\d+$", "", genes.paralogs)
rpkm.rna.seq.counts.paralogs <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id %in% 
    paralogs.no.expr.var), ]
rm(genes.paralogs)  # Clean up!


# For mean RPKM of all tissues:
message("Computing correlation between expression levels and paralog compy numbers..")
gene.copy.no.rpkm.corr.df <- correlationRpkmCopyNo(lapply.funk = mclapply)$data
message("Computing correlation between expression levels and tandem compy numbers..")
gene.tand.no.rpkm.corr.df <- correlationRpkmCopyNo(c.copy.no.col = "tandem.no", 
    lapply.funk = mclapply)$data
gc()  # Clean up!
# For each family with paralogs:
message("Computing correlation between expression levels and paralog compy numbers - for each family..")
fams.paralogs.copy.no.rpkm.corr <- setNames(mclapply(fams.paralogs, function(fam) {
    fam.genes.no.expr.var <- sub("\\.\\d+$", "", families.lst[[fam]])
    fam.rpkm.df <- rpkm.rna.seq.counts.paralogs[which(with(rpkm.rna.seq.counts.paralogs, 
        id %in% fam.genes.no.expr.var)), ]
    if (!is.null(fam.rpkm.df) && any(fam.rpkm.df$expression > 0)) {
        correlationRpkmCopyNo(rpkm.df = fam.rpkm.df)
    } else NA
}), fams.paralogs)
fams.paralogs.copy.no.rpkm.corr.df <- data.frame(Family = names(fams.paralogs.copy.no.rpkm.corr), 
    r.squared = as.numeric(lapply(fams.paralogs.copy.no.rpkm.corr, function(x) if ("r.squared" %in% 
        names(x) && !is.na(x$r.squared)) x$r.squared else NA)), stringsAsFactors = FALSE)
fams.paralogs.copy.no.rpkm.corr.df$description <- as.character(lapply(fams.paralogs.copy.no.rpkm.corr.df$Family, 
    function(fam) {
        families.HRD.df[which(families.HRD.df$id == fam), "description"]
    }))
fams.paralogs.copy.no.rpkm.corr.df$size <- as.character(lapply(fams.paralogs.copy.no.rpkm.corr.df$Family, 
    function(fam) {
        families.HRD.df[which(families.HRD.df$id == fam), "size"]
    }))
rm(fams.paralogs.copy.no.rpkm.corr)  # Clean up!
gc()
# For each tandem cluster:
tand.cls.of.interest <- tandems.lst[unique(tandems[which(tandems$Gene %in% 
    unlist(s.g.i)), "Family"])]
message("Computing correlation between expression levels and tandem compy numbers - for each tandem cluster..")
tandems.rpkm.corr <- setNames(mclapply(tand.cls.of.interest, function(tand.cl) {
    tand.cl.genes.no.expr.var <- sub("\\.\\d+$", "", tand.cl)
    tand.cl.rpkm.df <- rpkm.rna.seq.counts[which(with(rpkm.rna.seq.counts, 
        id %in% tand.cl.genes.no.expr.var)), ]
    if (!is.null(tand.cl.rpkm.df) && any(tand.cl.rpkm.df$expression > 0)) {
        correlationRpkmCopyNo(rpkm.df = tand.cl.rpkm.df, c.copy.no.col = "tandem.no")
    } else NA
}), names(tand.cls.of.interest))
tandems.rpkm.corr.df <- data.frame(Family = names(tandems.rpkm.corr), r.squared = as.numeric(lapply(tandems.rpkm.corr, 
    function(x) if ("r.squared" %in% names(x) && !is.na(x$r.squared)) x$r.squared else NA)), 
    stringsAsFactors = FALSE)
tandems.rpkm.corr.df$size <- as.numeric(lapply(tandems.rpkm.corr.df$Family, 
    function(x) {
        length(tandems.lst[[x]])
    }))
rm(tandems.rpkm.corr)  # Clean up!
gc()


# Save results:
save(gene.copy.no.rpkm.corr.df, gene.tand.no.rpkm.corr.df, fams.paralogs.copy.no.rpkm.corr.df, 
    tandems.rpkm.corr.df, file = file.path(input.args[[1]], "data", "correlationExpressionCopyNumber.RData"))


message("DONE")
