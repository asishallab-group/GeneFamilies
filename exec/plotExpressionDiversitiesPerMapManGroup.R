require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/plotExpressionDiversitiesPerMapManGroup.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Generate paired boxplots for each MapMan root-bin comparing function
#' diversity (in terms of domain-architecture-entropy) between orthologs and
#' families (that might intersect with orthologs!):
mm.bins <- setdiff(unique(families.mapMan.root.bins.df$description.mapMan), 
    "not assigned")
fams.plot.df <- Reduce(rbind, mclapply(mm.bins, function(mm.b) {
    mm.fams <- families.mapMan.root.bins.df[which(families.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    mm.orth.cls <- orthologs.mapMan.root.bins.df[which(orthologs.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    if (!is.null(mm.fams) && length(mm.fams) > 0) {
        Reduce(rbind, lapply(mm.fams, function(fam.nm) {
            my.df <- families.exp.prof.dists.orth.dist.df[which(with(families.exp.prof.dists.orth.dist.df, 
                Family == fam.nm & Statistic == "median")), ]
            if (!is.null(my.df) && nrow(my.df) > 0) {
                my.df$MapManBin <- mm.b
                my.df[, c("MapManBin", "Gene.Class", "Value")]
            } else NULL
        }))
    } else {
        NULL
    }
}))
#' Generate Boxplots for each of the MapMan root-bins:
mm.bins <- unique(fams.plot.df$MapManBin)
n.bins <- length(mm.bins)
pdf <- pdf(file.path(input.args[[1]], "inst", "expressionDiversityPerMapManRootBinFamiliesBoxplots.pdf"), 
    width = 21, height = 21)
old.par <- par(mfrow = c(ceiling(sqrt(n.bins)), floor(sqrt(n.bins))))
for (b in mm.bins) {
    plot.df <- fams.plot.df[which(fams.plot.df$MapManBin == b), c("Value", 
        "Gene.Class")]
    boxplot(Value ~ Gene.Class, data = plot.df, main = b, sub = paste("No Families =", 
        nrow(plot.df)))
}
dev.off()
par(old.par)


#' Generate paired boxplots for each MapMan root-bin comparing function
#' diversity (in terms of domain-architecture-entropy) between orthologs and
#' tandems (that might intersect with orthologs!):
mm.bins <- setdiff(unique(tandems.mapMan.root.bins.df$description.mapMan), 
    "not assigned")
tands.plot.df <- Reduce(rbind, mclapply(mm.bins, function(mm.b) {
    mm.tands <- tandems.mapMan.root.bins.df[which(tandems.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    mm.orth.cls <- orthologs.mapMan.root.bins.df[which(orthologs.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    if (!is.null(mm.tands) && length(mm.tands) > 0) {
        Reduce(rbind, lapply(mm.tands, function(fam.nm) {
            my.df <- tandems.exp.prof.dists.orth.dist.df[which(with(tandems.exp.prof.dists.orth.dist.df, 
                Family == fam.nm & Statistic == "median")), ]
            if (!is.null(my.df) && nrow(my.df) > 0) {
                my.df$MapManBin <- mm.b
                my.df[, c("MapManBin", "Gene.Class", "Value")]
            } else NULL
        }))
    } else {
        NULL
    }
}))
#' Generate Boxplots for each of the MapMan root-bins:
mm.bins <- unique(tands.plot.df$MapManBin)
n.bins <- length(mm.bins)
pdf <- pdf(file.path(input.args[[1]], "inst", "expressionDiversityPerMapManRootBinTandemsBoxplots.pdf"), 
    width = 21, height = 21)
old.par <- par(mfrow = c(ceiling(sqrt(n.bins)), floor(sqrt(n.bins))))
for (b in mm.bins) {
    plot.df <- tands.plot.df[which(tands.plot.df$MapManBin == b), c("Value", 
        "Gene.Class")]
    boxplot(Value ~ Gene.Class, data = plot.df, main = b, sub = paste("No Tandem-Clusters =", 
        nrow(plot.df)))
}
dev.off()
par(old.par)


message("DONE")
