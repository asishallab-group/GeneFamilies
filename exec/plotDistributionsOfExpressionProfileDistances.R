require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/plotDistributionsOfExpressionProfileDistances.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Plot distributions
#' - Compare Group.Class:
p.df <- Reduce(rbind, list(data.frame(Group.Class = "Orth", median.exp.prof.dist = orthologs.exp.prof.dists.stats.df[with(orthologs.exp.prof.dists.stats.df, 
    which(Statistic == "median")), "Value"]), data.frame(Group.Class = "Tand", median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[with(tandems.exp.prof.dists.orth.dist.df, 
    which(Statistic == "median")), "Value"]), data.frame(Group.Class = "Exp", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[with(families.exp.prof.dists.orth.dist.df, 
    which(Statistic == "median" & Family %in% families.exp)), "Value"]), data.frame(Group.Class = "PSel", 
    median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[with(families.exp.prof.dists.orth.dist.df, 
        which(Statistic == "median" & Family %in% fams.pos.sel)), "Value"])))
p.df$Group.Class <- factor(p.df$Group.Class, levels = c("Orth", "Tand", "Exp", "PSel"))
b.cols <- brewer.pal(length(unique(p.df$Group.Class)), "Dark2")
pdf(file.path(input.args[[1]], "inst", "medianExpressionProfileDistancesPerGeneGroupClassesBoxplot.pdf"))
boxplot(median.exp.prof.dist ~ Group.Class, data = p.df, outpch = "-", cex = 0.5, 
    border = b.cols, col = addAlpha(b.cols), lwd = 2.5, ylab = "median Euclidean distance between expression profiles")
dev.off()


#' - Comparing withing each Group.Class the Gene.Class:
plot.df <- Reduce(rbind, list(data.frame(Group.Class = "Tand", Gene.Class = "Orthologs", 
    median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[which(with(tandems.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Orthologs")), "Value"], stringsAsFactors = FALSE), 
    data.frame(Group.Class = "Tand", Gene.Class = "Non-Orthologs", median.exp.prof.dist = tandems.exp.prof.dists.orth.dist.df[which(with(tandems.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Non-Orthologs")), "Value"], stringsAsFactors = FALSE), 
    data.frame(Group.Class = "Exp", Gene.Class = "Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
        Statistic == "median" & Gene.Class == "Orthologs" & Family %in% families.exp)), 
        "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "Exp", Gene.Class = "Non-Orthologs", 
        median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
            Statistic == "median" & Gene.Class == "Non-Orthologs" & Family %in% families.exp)), 
            "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "PSel", 
        Gene.Class = "Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
            Statistic == "median" & Gene.Class == "Orthologs" & Family %in% fams.pos.sel)), 
            "Value"], stringsAsFactors = FALSE), data.frame(Group.Class = "PSel", 
        Gene.Class = "Non-Orthologs", median.exp.prof.dist = families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
            Statistic == "median" & Gene.Class == "Non-Orthologs" & Family %in% fams.pos.sel)), 
            "Value"], stringsAsFactors = FALSE)))
plot.df$Group.Class <- factor(plot.df$Group.Class, levels = c("Tand", "Exp", "PSel"))
g.pl <- ggplot(plot.df, aes(y = median.exp.prof.dist, x = Group.Class, fill = Gene.Class)) + 
    ylab("median Euclidean distance between expression profiles") + geom_boxplot(outlier.shape = "-", 
    outlier.size = 1) + theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.x = element_line(color = "black"), 
    axis.text.x = element_text(size = 12), axis.line.y = element_line(color = "black"), 
    legend.position = "bottom")
ggsave(plot = g.pl, filename = file.path(input.args[[1]], "inst", "medianExpressionProfileDistancesPerGeneGroupClassesAndOrthologDistinctionBoxplot.pdf"))


message("DONE")
