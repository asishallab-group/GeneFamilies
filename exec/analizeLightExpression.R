require(GeneFamilies)

message("USAGE: Rscript path/2/exec/analizeLightExpression.R path/2/ath.light.allgenes.txt path/chi.light.allgenes.txt path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Read light expression tables for the two species:
ath.light <- read.table(input.args[[1]], header = TRUE, sep = " ", stringsAsFactors = FALSE)
ath.light$id.lower.case <- tolower(ath.light$id)
chi.light <- read.table(input.args[[2]], header = TRUE, sep = " ", stringsAsFactors = FALSE)
chi.light$id.lower.case <- tolower(chi.light$id)
ath.chi.light <- rbind(ath.light, chi.light)


#' Ignore expression variants in MapMan-Bin Annotations:
mercator.df$IDENTIFIER.no.expr.var <- sub("\\.\\d+$", "", mercator.df$IDENTIFIER)


#' Using orthology, identify genes differentially expressed in both, and
#' uniquely in each of the species:
ath.light.de <- ath.light[which(ath.light$padj <= 0.05), ]
chi.light.de <- chi.light[which(chi.light$padj <= 0.05), ]
chi.light.de.orths <- orthologs[which(orthologs$chi %in% chi.light$id), 
    "ath"]
ath.light.de.orths <- orthologs[which(orthologs$ath %in% ath.light$id), 
    "chi"]
ath.light.de.uniq <- ath.light.de[which(!ath.light.de$id %in% chi.light.de.orths), 
    ]
chi.light.de.uniq <- chi.light.de[which(!chi.light.de$id %in% ath.light.de.orths), 
    ]
ath.chi.light.orth <- ath.chi.light[which(ath.chi.light$id %in% ath.light.de.orths | 
    ath.chi.light$id %in% chi.light.de.orths), ]


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
            log2FoldChange <- ath.chi.light[which(ath.chi.light$id.lower.case %in% 
                intersect(genes.de, mmBin.light.genes)), "log2FoldChange"]
            percent.up <- length(which(log2FoldChange > 0))/length(log2FoldChange)
            data.frame(BINCODE = mmBin, p.value = p.val, percent.up = percent.up, 
                stringsAsFactors = FALSE)
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


#' Fischer tests to identify enriched MapMan-Bins among diff expr genes -
#' Genes diff. expr. only in A.thaliana:.
ath.light.de.uniq.lc <- intersect(ath.light.de.uniq$id.lower.case, mercator.df$IDENTIFIER.no.expr.var)
ath.light.not.de.uniq.lc <- setdiff(intersect(ath.light$id.lower.case, 
    mercator.df$IDENTIFIER.no.expr.var), ath.light.de.uniq.lc)
ath.light.fish <- lightFischerTest(ath.light.de.uniq.lc, ath.light.not.de.uniq.lc)
#' Genes diff. expr. only in C.hirsuta:
chi.light.de.uniq.lc <- intersect(chi.light.de.uniq$id.lower.case, mercator.df$IDENTIFIER.no.expr.var)
chi.light.not.de.uniq.lc <- setdiff(intersect(chi.light$id.lower.case, 
    mercator.df$IDENTIFIER.no.expr.var), chi.light.de.uniq.lc)
chi.light.fish <- lightFischerTest(chi.light.de.uniq.lc, chi.light.not.de.uniq.lc)
#' Genes diff. expr. in both species based on orthology:
ath.chi.light.de.orth <- ath.chi.light.orth[which(ath.chi.light.orth$padj <= 
    0.05 & !ath.chi.light.orth$id %in% union(ath.light.de.uniq$id, chi.light.de.uniq$id) & 
    ath.chi.light.orth$id.lower.case %in% mercator.df$IDENTIFIER.no.expr.var), 
    "id.lower.case"]
ath.chi.light.not.de.orth <- ath.chi.light.orth[which(!ath.chi.light.orth$id.lower.case %in% 
    ath.chi.light.de.orth & ath.chi.light.orth$id.lower.case %in% mercator.df$IDENTIFIER.no.expr.var), 
    "id.lower.case"]
ath.chi.light.fish <- lightFischerTest(ath.chi.light.de.orth, ath.chi.light.not.de.orth)


#' Find DE genes of the above significant MapMan-Bins:
lightMapManSignGenes <- function(light.fisher.mapMan.sign, light.genes.de) {
    Reduce(rbind, lapply(light.fisher.mapMan.sign$BINCODE, function(mmBin) {
        mmBin.genes <- mercator.lst[[mmBin]]
        mmBin.genes.de <- intersect(tolower(sub("\\.\\d+$", "", mmBin.genes)), 
            light.genes.de$id.lower.case)
        data.frame(BINCODE = mmBin, NAME = light.fisher.mapMan.sign[which(light.fisher.mapMan.sign$BINCODE == 
            mmBin), "NAME"][[1]], gene = mmBin.genes.de, log2FoldChange = sapply(mmBin.genes.de, 
            function(x) light.genes.de[which(light.genes.de$id.lower.case == 
                x), "log2FoldChange"]), stringsAsFactors = FALSE)
    }))
}


#' Assess diff. expr. genes involved with significantly overrepresented MapMan-Bins - 
#' - in A.thaliana:
ath.light.de.sign.MapMan.genes <- lightMapManSignGenes(ath.light.fish, 
    ath.light.de)
#' - in C.hirsuta:
chi.light.de.sign.MapMan.genes <- lightMapManSignGenes(chi.light.fish, 
    chi.light.de)
#' - in both species:
ath.chi.light.de.sign.MapMan.genes <- lightMapManSignGenes(ath.chi.light.fish, 
    ath.chi.light.orth[which(ath.chi.light.orth$padj <= 0.05 & !ath.chi.light.orth$id %in% 
        union(ath.light.de.uniq$id, chi.light.de.uniq$id) & ath.chi.light.orth$id.lower.case %in% 
        mercator.df$IDENTIFIER.no.expr.var), ])


#' Write output files -
#' for A.thaliana: 
write.table(ath.light.fish[, setdiff(colnames(ath.light.fish), "p.value")], 
    file.path(input.args[[3]], "inst", "light_ath_fisher_tests_sign_mapManBins.txt"), 
    row.names = FALSE, sep = "\t", quote = FALSE)
write.table(ath.light.de.sign.MapMan.genes, file.path(input.args[[3]], 
    "inst", "light_ath_fisher_tests_sign_mapManBins_genes.txt"), row.names = FALSE, 
    sep = "\t", quote = FALSE)
#' for C.hirsuta:
write.table(chi.light.fish[, setdiff(colnames(chi.light.fish), "p.value")], 
    file.path(input.args[[3]], "inst", "light_chi_fisher_tests_sign_mapManBins.txt"), 
    row.names = FALSE, sep = "\t", quote = FALSE)
write.table(chi.light.de.sign.MapMan.genes, file.path(input.args[[3]], 
    "inst", "light_chi_fisher_tests_sign_mapManBins_genes.txt"), row.names = FALSE, 
    sep = "\t", quote = FALSE)
#' for both species:
write.table(ath.chi.light.fish[, setdiff(colnames(ath.chi.light.fish), 
    "p.value")], file.path(input.args[[3]], "inst", "light_ath_chi_fisher_tests_sign_mapManBins.txt"), 
    row.names = FALSE, sep = "\t", quote = FALSE)
write.table(ath.chi.light.de.sign.MapMan.genes, file.path(input.args[[3]], 
    "inst", "light_ath_chi_fisher_tests_sign_mapManBins_genes.txt"), row.names = FALSE, 
    sep = "\t", quote = FALSE)


#' Save results:
save(ath.light, chi.light, ath.light.fish, chi.light.fish, ath.chi.light.fish, 
    file = file.path(input.args[[3]], "data", "lightExpressionData.RData"))

message("DONE")
