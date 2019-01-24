require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/expressionAngles.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' **************************************
#' Analysis of expression vector space. *
#' **************************************


#' Define expression vector space:
expr.cols <- c("cotyledon", "seedling", "developing leaf", "flower stage 9", 
    "flower stage 16")
n.dims <- length(expr.cols)
message("Defining expression vector space with the following axes: ", paste(expr.cols, 
    collapse = ", "))


#' Infer angle to diagonal as measure of tissue specificity:
tands.expr.angle.diag.df <- data.frame(gene = tands.expr, angle.diag = as.numeric(mclapply(tands.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$gene.exp.var == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)
orths.expr.angle.diag.df <- data.frame(gene = orths.expr, angle.diag = as.numeric(mclapply(orths.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)
dupl.expr.angle.diag.df <- data.frame(gene = dupl.expr, angle.diag = as.numeric(mclapply(dupl.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$gene.exp.var == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)
psel.expr.angle.diag.df <- data.frame(gene = psel.expr, angle.diag = as.numeric(mclapply(psel.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$gene.exp.var == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)


#' Plot results:
p.lst <- list(tandem = tands.expr.angle.diag.df, ortholog = orths.expr.angle.diag.df, 
    duplicated = dupl.expr.angle.diag.df, pos.selected = psel.expr.angle.diag.df)
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))
plot.df <- p.df[!is.nan(p.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = c("ortholog", "duplicated", 
    "tandem"))
pdf(file.path(input.args[[1]], "inst", "expressionAngleToDiagonalBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(angle.diag ~ gene.type, data = plot.df, xlab = "Type of Gene", 
    ylab = "relative tissue specificity", border = pushed.colors, col = addAlpha(pushed.colors))
dev.off()
plot.df$rel.vers <- 1 - plot.df$angle.diag
pdf(file.path(input.args[[1]], "inst", "relativeExpressionVersatilityBoxplot.pdf"))
boxplot(rel.vers ~ gene.type, data = plot.df, xlab = "Type of Gene", ylab = "relative tissue versatility", 
    border = pushed.colors, col = addAlpha(pushed.colors), outline = FALSE)
dev.off()


#' Investigate what happened in each gene group after duplication:
#' - Comparing tandems and orthologs
tands.w.orths.angles.df <- Reduce(rbind, mclapply(names(tands.w.orths), 
    function(fam.name) {
        genes <- intersect(tands.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, tand.classifier)
    }))
#' - Comparing tandems, pos.selected, and orthologs
tands.psel.w.orths.angles.df <- Reduce(rbind, mclapply(names(tands.w.orths), 
    function(fam.name) {
        genes <- intersect(tands.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, tand.psel.classifier)
    }))
#' - Comparing duplicated and orthologs
dupl.w.orths.angles.df <- Reduce(rbind, mclapply(names(dupl.w.orths), function(fam.name) {
    genes <- intersect(dupl.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
    exprVecSpaceEvolutionAfterDupl(genes, fam.name, dupl.classifier)
}))
#' - Comparing duplicated, pos.selected, and orthologs
dupl.psel.w.orths.angles.df <- Reduce(rbind, mclapply(names(dupl.w.orths), 
    function(fam.name) {
        genes <- intersect(dupl.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, dupl.psel.classifier)
    }))


#' Plot the above:
p.lst <- list(tandem = tands.w.orths.angles.df[which(tands.w.orths.angles.df$tandem.dist.vec.clouds > 
    0), ], duplicated = dupl.w.orths.angles.df[which(dupl.w.orths.angles.df$duplicated.dist.vec.clouds > 
    0), ])
p.df <- Reduce(rbind, lapply(names(p.lst), function(x) {
    tryCatch({
        y <- p.lst[[x]]
        data.frame(group.type = x, diff.tissue.versatility = (y[, paste(x, 
            ".tiss.vers", sep = "")] - y$mean.base.tiss.vers), angle.between.orth.2.diag.vecs = y[, 
            paste(x, ".tiss.change", sep = "")], stringsAsFactors = FALSE)
    }, error = function(e) browser())
}))
pdf(file.path(input.args[[1]], "inst", "meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf"))
boxplot(diff.tissue.versatility ~ group.type, data = p.df, col = addAlpha(colors), 
    outline = FALSE, border = colors, xlab = "Type of gene group", ylab = "Mean difference of tissue versatility after duplication", 
    pch = "-")
stripchart(diff.tissue.versatility ~ group.type, vertical = TRUE, data = p.df, 
    method = "jitter", add = TRUE, pch = ".", col = colors)
dev.off()
pdf(file.path(input.args[[1]], "inst", "meanTissueVersatilityDiffsAfterDuplicationHistograms.pdf"))
old.par <- par(mfrow = c(2, 1))
xlim <- c(min(p.df$diff.tissue.versatility, na.rm = TRUE), max(p.df$diff.tissue.versatility, 
    na.rm = TRUE))
hist(p.df[which(p.df$group.type == "duplicated"), "diff.tissue.versatility"], 
    breaks = 30, main = "duplicated", col = addAlpha(colors[[1]]), border = colors[[1]], 
    xlab = "", xlim = xlim)
hist(p.df[which(p.df$group.type == "tandem"), "diff.tissue.versatility"], 
    breaks = 30, main = "tandem", col = addAlpha(colors[[2]]), border = colors[[2]], 
    xlab = "Mean difference of tissue versatility after duplication", xlim = xlim)
dev.off()
par(old.par)


pdf(file.path(input.args[[1]], "inst", "afterDuplicationAngleBetweenOrth2DiagVecsBoxplot.pdf"))
boxplot(angle.between.orth.2.diag.vecs ~ group.type, data = p.df, col = addAlpha(colors), 
    outline = FALSE, border = colors, xlab = "Type of gene group", ylab = "Angle between orthogonals on expression space diagional", 
    pch = "-")
stripchart(angle.between.orth.2.diag.vecs ~ group.type, vertical = TRUE, 
    data = p.df, method = "jitter", add = TRUE, pch = ".", col = colors)
dev.off()
pdf(file.path(input.args[[1]], "inst", "afterDuplicationAngleBetweenOrth2DiagVecsHistograms.pdf"))
old.par <- par(mfrow = c(2, 1))
xlim <- c(min(p.df$angle.between.orth.2.diag.vecs, na.rm = TRUE), max(p.df$angle.between.orth.2.diag.vecs, 
    na.rm = TRUE))
hist(p.df[which(p.df$group.type == "duplicated"), "angle.between.orth.2.diag.vecs"], 
    breaks = 30, main = "duplicated", col = addAlpha(colors[[1]]), border = colors[[1]], 
    xlab = "", xlim = xlim)
hist(p.df[which(p.df$group.type == "tandem"), "angle.between.orth.2.diag.vecs"], 
    breaks = 30, main = "tandem", col = addAlpha(colors[[2]]), border = colors[[2]], 
    xlab = "Angle between orthogonals on expression space diagional", xlim = xlim)
dev.off()
par(old.par)

#' Plot both changes in tissue specificity:
p.df.1 <- dupl.w.orths.angles.df[with(dupl.w.orths.angles.df, which(!is.na(duplicated.tiss.change) & 
    !is.nan(duplicated.tiss.change))), ]
p.df.1$class <- "Duplicated"
p.df.2 <- tands.w.orths.angles.df[with(tands.w.orths.angles.df, which(!is.na(tandem.tiss.change) & 
    !is.nan(tandem.tiss.change))), ]
colnames(p.df.2) <- sub("tandem", "duplicated", colnames(p.df.2))
p.df.2$class <- "Tandem"
cols.i <- c("class")
p.df <- rbind(p.df.1, p.df.2)

#' Plot all
pdf(file.path(input.args[[1]], "inst", "meanChangeInAngleBetweenOrthOnDiagsBoxplot.pdf"))
pushed.colors <- brewer.pal(3, "Dark2")
boxplot(duplicated.tiss.change ~ class, data = p.df, xlab = "Type of Gene", 
    ylab = "angle between diagonal orthologs", border = pushed.colors, 
    col = addAlpha(pushed.colors))
dev.off()


#' Boxplot log2 fold change of ortholog, tandem, and duplicated expressed
#' genes. Gene Expression is inferred on the intra-species level during fruit
#' development and shade reaction.
p.df <- Reduce(rbind, mclapply(list(fruit.ath.intra, fruit.chi.intra, light.ath.intra, 
    light.chi.intra), function(x) {
    y <- x[which(!is.na(x$log2FoldChange) & !is.infinite(x$log2FoldChange) & 
        !is.null(x$log2FoldChange)), c("id", "log2FoldChange")]
    y$gene.class <- NA
    y$absLog2FoldChange <- abs(y$log2FoldChange)
    gene.sets <- list(ortholog = orths.expr, tandem = sub("\\.\\d+$", "", 
        tands.expr), duplicated = sub("\\.\\d+$", "", dupl.expr))
    for (gene.set.name in names(gene.sets)) {
        gene.set <- gene.sets[[gene.set.name]]
        i <- which(x$id %in% gene.set)
        y[i, "gene.class"] <- gene.set.name
    }
    y
}))
p.df$gene.class <- factor(p.df$gene.class, levels = c("ortholog", "duplicated", 
    "tandem"))
p.df.sel <- p.df[which(!is.na(p.df$gene.class) & p.df$absLog2FoldChange >= 
    2), ]

colors <- brewer.pal(length(unique(p.df.sel$gene.class)), "Dark2")
pdf(file.path(input.args[[1]], "inst", "absLog2FoldChangeIntraSpeciesBoxplots.pdf"))
boxplot(absLog2FoldChange ~ gene.class, data = p.df.sel, col = addAlpha(colors), 
    outline = FALSE, border = colors, xlab = "Type of gene group", ylab = "absolute Log2FoldChange", 
    main = "Fold-Change within intra-species\ndifferentially expressed genes", 
    pch = "-")
stripchart(absLog2FoldChange ~ gene.class, vertical = TRUE, data = p.df.sel, 
    method = "jitter", add = TRUE, pch = ".", col = colors)
dev.off()


message("DONE")
