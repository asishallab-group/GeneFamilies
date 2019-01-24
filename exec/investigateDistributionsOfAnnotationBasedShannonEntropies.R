require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfAnnotationBasedShannonEntropies.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

#' Plot the distributions of function diversities for the different types of
#' gene-groups:
entropies.lst <- list(absolute = list(Orth = orthologs.ipr.entropies, Tand = tandems.ipr.entropies, 
    Exp = families.ipr.entropies[families.exp], PSel = families.ipr.entropies[fams.pos.sel]), 
    normalized = list(Orth = orthologs.ipr.norm.entropies, Tand = tandems.ipr.norm.entropies, 
        Exp = families.ipr.norm.entropies[families.exp], PSel = families.ipr.norm.entropies[fams.pos.sel]))
for (entr.type in names(entropies.lst)) {
    for (with.pos.sel in c(TRUE, FALSE)) {
        if (with.pos.sel) {
            pdf(file.path(input.args[[1]], "inst", paste(entr.type, "FunctionDiversitiesPerGeneGroupClassesBoxplot.pdf", 
                sep = "")))
            plot.lst <- entropies.lst[[entr.type]]
            lvls <- c("Orth", "Tand", "Exp", "PSel")
        } else {
            pdf(file.path(input.args[[1]], "inst", paste(entr.type, "FunctionDiversitiesPerGeneGroupClassesBoxplot_NoPosSel.pdf", 
                sep = "")))
            plot.lst <- entropies.lst[[entr.type]]
            plot.lst <- plot.lst[setdiff(names(plot.lst), "PSel")]
            lvls <- c("Orth", "Tand", "Exp")
        }
        plot.df <- Reduce(rbind, lapply(names(plot.lst), function(x) {
            data.frame(Entropy = unlist(plot.lst[[x]]), Group.Class = x, 
                stringsAsFactors = FALSE)
        }))
        plot.df$Group.Class <- factor(plot.df$Group.Class, levels = lvls)
        b.cols <- brewer.pal(length(plot.lst), "Dark2")
        boxplot(Entropy ~ Group.Class, data = plot.df, outpch = "-", cex = 0.5, 
            border = b.cols, col = addAlpha(b.cols), lwd = 2.5)
        dev.off()
    }
}


#' Plot the same, but comparing within each group 'background' orthologs with
#' non-orthologs:
entropies.comp.lst <- list(absolute = list(Tand = tandems.ipr.entropies.orth.dist, 
    Exp = families.ipr.entropies.orth.dist[which(families.ipr.entropies.orth.dist$Family %in% 
        families.exp), ], PSel = families.ipr.entropies.orth.dist[which(families.ipr.entropies.orth.dist$Family %in% 
        fams.pos.sel), ]), normalized = list(Tand = tandems.ipr.norm.entropies.orth.dist, 
    Exp = families.ipr.norm.entropies.orth.dist[which(families.ipr.norm.entropies.orth.dist$Family %in% 
        families.exp), ], PSel = families.ipr.norm.entropies.orth.dist[which(families.ipr.norm.entropies.orth.dist$Family %in% 
        fams.pos.sel), ]))
for (entr.type in names(entropies.comp.lst)) {
    for (with.pos.sel in c(TRUE, FALSE)) {
        if (with.pos.sel) {
            oplot.lst <- entropies.comp.lst[[entr.type]]
            plot.file <- "FunctionDiversitiesPerGeneGroupClassesAndOrthologDistinctionBoxplot.pdf"
        } else {
            oplot.lst <- entropies.comp.lst[[entr.type]]
            oplot.lst <- oplot.lst[setdiff(names(oplot.lst), "PSel")]
            plot.file <- "FunctionDiversitiesPerGeneGroupClassesAndOrthologDistinctionBoxplot_NoPosSel.pdf"
        }
        oplot.df <- Reduce(rbind, lapply(names(oplot.lst), function(x) {
            x.df <- oplot.lst[[x]]
            y.df <- rbind(data.frame(Entropy = x.df$Orthologs, Gene.Class = "Orthologs", 
                stringsAsFactors = FALSE), data.frame(Entropy = x.df$"Non-Orthologs", 
                Gene.Class = "Non-Orthologs", stringsAsFactors = FALSE))
            y.df$Group.Class <- x
            y.df
        }))
        oplot.df$Group.Class <- factor(oplot.df$Group.Class, levels = names(oplot.lst))
        g.pl <- ggplot(oplot.df, aes(y = Entropy, x = Group.Class, fill = Gene.Class)) + 
            geom_boxplot(outlier.shape = "-", outlier.size = 1) + theme_bw() + 
            theme(plot.background = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.border = element_blank(), 
                axis.line.x = element_line(color = "black"), axis.text.x = element_text(size = 12), 
                axis.line.y = element_line(color = "black"), legend.position = "bottom")
        
        ggsave(plot = g.pl, filename = file.path(input.args[[1]], "inst", 
            paste(entr.type, plot.file, sep = "")))
        if (with.pos.sel) {
            #' Do t-Tests for each group, testing the alternative hypothesis that orthologs
            #' have lower diversities than non-orthologs:
            ipr.entropies.orth.dist.t.tests <- Reduce(rbind, lapply(names(oplot.lst), 
                function(x) {
                  x.df <- oplot.lst[[x]]
                  t.h.test <- t.test(x.df[, "Orthologs"], x.df[, "Non-Orthologs"], 
                    alternative = "less")
                  w.h.test <- wilcox.test(x.df[, "Orthologs"], x.df[, "Non-Orthologs"], 
                    alternative = "less")
                  data.frame(Group.Class = x, Alternative.Hypothesis = "orthologs have lower mean (t-test) or diversity distribution (wilcox test) than non-orthologs", 
                    t.p.value = t.h.test$p.value, effect.size.t = as.numeric(t.h.test$statistic), 
                    w.p.value = w.h.test$p.value, stringsAsFactors = FALSE)
                }))
            ipr.entropies.orth.dist.t.tests$t.p.adj <- p.adjust(ipr.entropies.orth.dist.t.tests$t.p.value, 
                method = "fdr")
            ipr.entropies.orth.dist.t.tests$w.p.adj <- p.adjust(ipr.entropies.orth.dist.t.tests$w.p.value, 
                method = "fdr")
            ipr.entropies.orth.inter.group.class.t.tests <- data.frame(Alternative.Hypothesis = c("Expanded families' non-orthologs have higher mean (t.test) diversity (wilcox.test) than positively selected families", 
                "Positively selected families' non-orthologs have higher mean (t.test) diversity (wilcox.test) than tandem clusters", 
                "Expanded families non-orthologs have higher mean (t.test) diversity (wilcox.test) than tandem clusters"), 
                t.p.value = rep(NA, 3), effect.size.t = rep(NA, 3), w.p.value = rep(NA, 
                  3), stringsAsFactors = FALSE)
            t.1 <- t.test(oplot.lst[["PSel"]]$"Non-Orthologs", oplot.lst[["Exp"]]$"Non-Orthologs", 
                alternative = "less")
            w.1 <- wilcox.test(oplot.lst[["PSel"]]$"Non-Orthologs", oplot.lst[["Exp"]]$"Non-Orthologs", 
                alternative = "less")
            ipr.entropies.orth.inter.group.class.t.tests[1, "t.p.value"] <- t.1$p.value
            ipr.entropies.orth.inter.group.class.t.tests[1, "effect.size.t"] <- as.numeric(t.1$statistic)
            ipr.entropies.orth.inter.group.class.t.tests[1, "w.p.value"] <- w.1$p.value
            t.2 <- t.test(oplot.lst[["Tand"]]$"Non-Orthologs", oplot.lst[["PSel"]]$"Non-Orthologs", 
                alternative = "less")
            w.2 <- wilcox.test(oplot.lst[["Tand"]]$"Non-Orthologs", oplot.lst[["PSel"]]$"Non-Orthologs", 
                alternative = "less")
            ipr.entropies.orth.inter.group.class.t.tests[2, "t.p.value"] <- t.2$p.value
            ipr.entropies.orth.inter.group.class.t.tests[2, "effect.size.t"] <- as.numeric(t.2$statistic)
            ipr.entropies.orth.inter.group.class.t.tests[2, "w.p.value"] <- w.2$p.value
            t.3 <- t.test(oplot.lst[["Tand"]]$"Non-Orthologs", oplot.lst[["Exp"]]$"Non-Orthologs", 
                alternative = "less")
            w.3 <- wilcox.test(oplot.lst[["Tand"]]$"Non-Orthologs", oplot.lst[["Exp"]]$"Non-Orthologs", 
                alternative = "less")
            ipr.entropies.orth.inter.group.class.t.tests[3, "t.p.value"] <- t.3$p.value
            ipr.entropies.orth.inter.group.class.t.tests[3, "effect.size.t"] <- as.numeric(t.3$statistic)
            ipr.entropies.orth.inter.group.class.t.tests[3, "w.p.value"] <- w.3$p.value
            ipr.entropies.orth.inter.group.class.t.tests$t.p.adj <- p.adjust(ipr.entropies.orth.inter.group.class.t.tests$t.p.value, 
                method = "fdr")
            ipr.entropies.orth.inter.group.class.t.tests$w.p.adj <- p.adjust(ipr.entropies.orth.inter.group.class.t.tests$w.p.value, 
                method = "fdr")
            #' Save results:
            write.table(ipr.entropies.orth.dist.t.tests, file = file.path(input.args[[1]], 
                "inst", paste(entr.type, "FunctionDiversitiesPerGeneGroupClasses_h-tests.tsv", 
                  sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
            write.table(ipr.entropies.orth.inter.group.class.t.tests, file = file.path(input.args[[1]], 
                "inst", paste(entr.type, "FunctionDiversitiesInterGeneGroupClasses_h-tests.tsv", 
                  sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
}


message("DONE")
