require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/plotGeneGroupSizes.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

tandems.df <- data.frame(Family = unique(tandems$Family), stringsAsFactors = FALSE)
tandems.df$size <- as.integer(unlist(mclapply(tandems.df$Family, function(x) length(which(tandems$Family == 
    x)))))


#' Plot
p.df <- rbind(rbind(data.frame(size = tandems.df$size, group.type = "Tandem", 
    stringsAsFactors = FALSE), data.frame(size = families.df$size, group.type = "Family", 
    stringsAsFactors = FALSE)), data.frame(size = families.df[which(families.df$id %in% 
    families.exp.cnt), "size"], group.type = "Trans Duplicated", stringsAsFactors = FALSE))

colors <- brewer.pal(3, "Set1")
pdf(file.path(input.args[[1]], "inst", "geneFamilyAndTandemClusterSizesBoxplot.pdf"))
boxplot(size ~ group.type, data = p.df, col = addAlpha(colors), outline = FALSE, 
    border = colors, xlab = "Type of gene group", ylab = "Number of genes")
stripchart(size ~ group.type, vertical = TRUE, data = p.df, method = "jitter", 
    add = TRUE, pch = ".", col = colors)
dev.off()


#' Save results:
save(tandems.df, file = file.path(input.args[[1]], "data", "tandemsDf.RData"))

message("DONE")
