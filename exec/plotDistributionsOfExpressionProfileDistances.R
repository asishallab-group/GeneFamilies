require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/plotDistributionsOfExpressionProfileDistances.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Plot distributions
#' - Compare Group.Class:
p.df <- Reduce(rbind, list(data.frame(Group.Class = "Orth", median.exp.prof.dist = orthologs.exp.prof.dists.stats.df[with(orthologs.exp.prof.dists.stats.df, 
    which(Statistic == "median")), "Value"]), data.frame(Group.Class = "Tand", 
    median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[with(tandems.exp.prof.dists.orth.dist.df, 
        which(Statistic == "median")), "Value"]), data.frame(Group.Class = "Exp", 
    median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[with(families.exp.prof.dists.orth.dist.df, 
        which(Statistic == "median" & Family %in% families.exp)), "Value"]), 
    data.frame(Group.Class = "PSel", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[with(families.exp.prof.dists.orth.dist.df, 
        which(Statistic == "median" & Family %in% fams.pos.sel)), "Value"])))
b.cols <- brewer.pal(length(unique(p.df$Group.Class)), "Dark2")
p.df$Group.Class <- factor(p.df$Group.Class, levels = c("Orth", "Tand", 
    "Exp", "PSel"))
#' Plot twice, once with and once without 'PSel':
pl.df <- p.df
for (w.psel in c(TRUE, FALSE)) {
    if (w.psel) {
        pl.fl <- "medianExpressionProfileDistancesPerGeneGroupClassesBoxplot.pdf"
        pl.df$Group.Class <- factor(pl.df$Group.Class, levels = c("Orth", 
            "Tand", "Exp", "PSel"))
    } else {
        pl.fl <- "medianExpressionProfileDistancesPerGeneGroupClassesBoxplot_NoPosSel.pdf"
        pl.df <- pl.df[which(pl.df$Group.Class != "PSel"), ]
        pl.df$Group.Class <- factor(pl.df$Group.Class, levels = c("Orth", 
            "Tand", "Exp"))
    }
    pdf(file.path(input.args[[1]], "inst", pl.fl))
    boxplot(median.exp.prof.dist ~ Group.Class, data = pl.df, outpch = "-", 
        cex = 0.5, border = b.cols, col = addAlpha(b.cols), lwd = 2.5, 
        ylab = "median Euclidean distance between expression profiles")
    dev.off()
}


#' Test for significant difference among the distributions:
inter.group.tests <- cbind(data.frame(alternative.hypothesis = c("Dist of median expr. dists of Tand is above Exp", 
    "Dist of median expr. dists of Exp is above Orth", "Dist of median expr. dists of PSel is above Orth", 
    "Dist of median expr. dists of Exp is above PSel"), stringsAsFactors = FALSE), 
    Reduce(rbind, lapply(list(c("Tand", "Exp"), c("Exp", "Orth"), c("PSel", 
        "Orth"), c("Exp", "PSel")), function(dist.pair) {
        pair.t <- t.test(p.df[which(p.df$Group.Class == dist.pair[[2]]), 
            "median.exp.prof.dist"], p.df[which(p.df$Group.Class == dist.pair[[1]]), 
            "median.exp.prof.dist"], alternative = "less")$p.value
        pair.w <- wilcox.test(p.df[which(p.df$Group.Class == dist.pair[[2]]), 
            "median.exp.prof.dist"], p.df[which(p.df$Group.Class == dist.pair[[1]]), 
            "median.exp.prof.dist"], alternative = "less")$p.value
        data.frame(t.test.p.val = pair.t, wilcox.test.p.val = pair.w, stringsAsFactors = FALSE)
    })))
#' Adjust for multiple hypothesis testing:
inter.group.tests$t.test.p.adj <- p.adjust(inter.group.tests$t.test.p.val, 
    method = "fdr")
inter.group.tests$wilcox.test.p.adj <- p.adjust(inter.group.tests$wilcox.test.p.val, 
    method = "fdr")
write.table(inter.group.tests, file.path(input.args[[1]], "inst", "medianExpressionProfileDistancesPerGeneGroupClassesBoxplot_h-tests.txt"), 
    row.names = FALSE, sep = "\t")


#' - Comparing withing each Group.Class the Gene.Class:
plot.df <- Reduce(rbind, list(data.frame(Group.Class = "Tand", Gene.Class = "Orthologs", 
    median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[which(with(tandems.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Orthologs")), "Value"], 
    stringsAsFactors = FALSE), data.frame(Group.Class = "Tand", Gene.Class = "Non-Orthologs", 
    median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[which(with(tandems.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Non-Orthologs")), "Value"], 
    stringsAsFactors = FALSE), data.frame(Group.Class = "Exp", Gene.Class = "Orthologs", 
    median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Orthologs" & Family %in% 
            families.exp)), "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "Exp", 
    Gene.Class = "Non-Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Non-Orthologs" & Family %in% 
            families.exp)), "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "PSel", 
    Gene.Class = "Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Orthologs" & Family %in% 
            fams.pos.sel)), "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "PSel", 
    Gene.Class = "Non-Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Non-Orthologs" & Family %in% 
            fams.pos.sel)), "Value"], stringsAsFactors = FALSE)))
plot.df$Group.Class <- factor(plot.df$Group.Class, levels = c("Tand", "Exp", 
    "PSel"))
#' Plot twice, once with and once without 'PSel':
pl.df <- plot.df
for (w.psel in c(TRUE, FALSE)) {
    if (w.psel) {
        pl.fl <- "medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot.pdf"
    } else {
        pl.fl <- "medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot_NoPosSel.pdf"
        pl.df <- pl.df[which(pl.df$Group.Class != "PSel"), ]
    }
    g.pl <- ggplot(pl.df, aes(y = median.exp.prof.dist, x = Group.Class, 
        fill = Gene.Class)) + ylab("median Euclidean distance between expression profiles") + 
        geom_boxplot(outlier.shape = "-", outlier.size = 1) + theme_bw() + 
        theme(plot.background = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.border = element_blank(), 
            axis.line.x = element_line(color = "black"), axis.text.x = element_text(size = 12), 
            axis.line.y = element_line(color = "black"), legend.position = "bottom")
    ggsave(plot = g.pl, filename = file.path(input.args[[1]], "inst", pl.fl))
}
#' Test for significant difference among the distributions:
intra.group.tests <- rbind(cbind(data.frame(alternative.hypothesis = c("Dist of median expr. dists of Tand\\Orth is above Tand^Orth", 
    "Dist of median expr. dists of Exp\\Orth is above Exp^Orth", "Dist of median expr. dists of PSel\\Orth is above PSel^Orth"), 
    stringsAsFactors = FALSE), Reduce(rbind, lapply(list("Tand", "Exp", 
    "PSel"), function(gr.cl) {
    pair.t <- t.test(plot.df[which(plot.df$Group.Class == gr.cl & plot.df$Gene.Class == 
        "Non-Orthologs"), "median.exp.prof.dist"], plot.df[which(plot.df$Group.Class == 
        gr.cl & plot.df$Gene.Class == "Orthologs"), "median.exp.prof.dist"], 
        alternative = "less")$p.value
    pair.w <- wilcox.test(plot.df[which(plot.df$Group.Class == gr.cl & 
        plot.df$Gene.Class == "Orthologs"), "median.exp.prof.dist"], plot.df[which(plot.df$Group.Class == 
        gr.cl & plot.df$Gene.Class == "Non-Orthologs"), "median.exp.prof.dist"], 
        alternative = "less")$p.value
    data.frame(t.test.p.val = pair.t, wilcox.test.p.val = pair.w, stringsAsFactors = FALSE)
}))), data.frame(alternative.hypothesis = "Dist of median expr. dists of Tand\\Orth is above Exp\\Orth", 
    t.test.p.val = t.test(plot.df[which(plot.df$Group.Class == "Tand" & 
        plot.df$Gene.Class == "Non-Orthologs"), "median.exp.prof.dist"], 
        plot.df[which(plot.df$Group.Class == "Exp" & plot.df$Gene.Class == 
            "Non-Orthologs"), "median.exp.prof.dist"], alternative = "less")$p.value, 
    wilcox.test.p.val = wilcox.test(plot.df[which(plot.df$Group.Class == 
        "Tand" & plot.df$Gene.Class == "Non-Orthologs"), "median.exp.prof.dist"], 
        plot.df[which(plot.df$Group.Class == "Exp" & plot.df$Gene.Class == 
            "Non-Orthologs"), "median.exp.prof.dist"], alternative = "less")$p.value, 
    stringsAsFactors = FALSE))
#' Adjust for multiple hypothesis testing:
intra.group.tests$t.test.p.adj <- p.adjust(intra.group.tests$t.test.p.val, 
    method = "fdr")
intra.group.tests$wilcox.test.p.adj <- p.adjust(intra.group.tests$wilcox.test.p.val, 
    method = "fdr")
write.table(intra.group.tests, file.path(input.args[[1]], "inst", "medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot_h-tests.txt"), 
    row.names = FALSE, sep = "\t")


#' Investigate further the distributions of non-ortholog Tandems and non-ortholog Expanded:
my.df <- plot.df[with(plot.df, which(Group.Class %in% c("Tand", "Exp") & 
    Gene.Class == "Non-Orthologs" & !is.na(median.exp.prof.dist))), ]
g.p <- ggplot(my.df, aes(median.exp.prof.dist, fill = Group.Class)) + geom_histogram(alpha = 0.5, 
    aes(y = ..density..), position = "identity", binwidth = .01)
ggsave(g.p, filename = "./inst/medianExpressionDistsTandemNonOrthsAndExpandedNonOrthsHist.pdf")


message("DONE")
