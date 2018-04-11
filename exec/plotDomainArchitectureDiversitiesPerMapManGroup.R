require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/plotDomainArchitectureDiversitiesPerMapManGroup.R path/2/MaizeGeneFamilies")

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
    rbind(if (!is.null(mm.fams) && length(mm.fams) > 0) {
        Reduce(rbind, lapply(mm.fams, function(fam.nm) {
            data.frame(MapManBin = mm.b, Gene.Group = "Family", Value = families.ipr.entropies[[fam.nm]], 
                stringsAsFactors = FALSE)
        }))
    } else {
        NULL
    }, if (!is.null(mm.orth.cls) && length(mm.orth.cls) > 0) {
        Reduce(rbind, lapply(mm.orth.cls, function(o.cl.nm) {
            data.frame(MapManBin = mm.b, Gene.Group = "Ortholog-Cluster", 
                Value = orthologs.ipr.entropies[[o.cl.nm]], stringsAsFactors = FALSE)
        }))
    } else {
        NULL
    })
}))
#' Generate Boxplots for each of the MapMan root-bins:
mm.bins <- unique(fams.plot.df$MapManBin)
n.bins <- length(mm.bins)
pdf <- pdf(file.path(input.args[[1]], "inst", "domainArchitectureDiversityPerMapManRootBinFamiliesBoxplots.pdf"), 
    width = 21, height = 21)
old.par <- par(mfrow = c(ceiling(sqrt(n.bins)), floor(sqrt(n.bins))))
for (b in mm.bins) {
    plot.df <- fams.plot.df[which(fams.plot.df$MapManBin == b), c("Value", 
        "Gene.Group")]
    s.t <- paste("No Families=", length(which(plot.df$Gene.Group == "Family")), 
        " , No Ortholog-Clusters=", length(which(plot.df$Gene.Group == 
            "Ortholog-Cluster")), sep = "")
    boxplot(Value ~ Gene.Group, data = plot.df, main = b, sub = s.t)
}
dev.off()
par(old.par)


#' Generate paired boxplots for each MapMan root-bin comparing function
#' diversity (in terms of domain-architecture-entropy) between orthologs and
#' tandems:
mm.bins <- setdiff(unique(tandems.mapMan.root.bins.df$description.mapMan), 
    "not assigned")
tands.plot.df <- Reduce(rbind, mclapply(mm.bins, function(mm.b) {
    mm.tands <- tandems.mapMan.root.bins.df[which(tandems.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    mm.orth.cls <- orthologs.mapMan.root.bins.df[which(orthologs.mapMan.root.bins.df$description.mapMan == 
        mm.b), "id"]
    rbind(if (!is.null(mm.tands) && length(mm.tands) > 0) {
        Reduce(rbind, lapply(mm.tands, function(fam.nm) {
            data.frame(MapManBin = mm.b, Gene.Group = "Tandem-Cluster", 
                Value = tandems.ipr.entropies[[fam.nm]], stringsAsFactors = FALSE)
        }))
    } else {
        NULL
    }, if (!is.null(mm.orth.cls) && length(mm.orth.cls) > 0) {
        Reduce(rbind, lapply(mm.orth.cls, function(o.cl.nm) {
            data.frame(MapManBin = mm.b, Gene.Group = "Ortholog-Cluster", 
                Value = orthologs.ipr.entropies[[o.cl.nm]], stringsAsFactors = FALSE)
        }))
    } else {
        NULL
    })
}))
#' Generate Boxplots for each of the MapMan root-bins:
mm.bins <- unique(tands.plot.df$MapManBin)
n.bins <- length(mm.bins)
pdf <- pdf(file.path(input.args[[1]], "inst", "domainArchitectureDiversityPerMapManRootBinTandemsBoxplots.pdf"), 
    width = 21, height = 21)
old.par <- par(mfrow = c(ceiling(sqrt(n.bins)), floor(sqrt(n.bins))))
for (b in mm.bins) {
    plot.df <- tands.plot.df[which(tands.plot.df$MapManBin == b), c("Value", 
        "Gene.Group")]
    #' Inverse alphabetical order on X-Axis, please:
    plot.df$Gene.Group <- factor(plot.df$Gene.Group, levels = c("Tandem-Cluster", 
        "Ortholog-Cluster"))
    s.t <- paste("No Tandem-Clstr=", length(which(plot.df$Gene.Group == 
        "Tandem-Cluster")), " , No Ortholog-Clstr=", length(which(plot.df$Gene.Group == 
        "Ortholog-Cluster")), sep = "")
    boxplot(Value ~ Gene.Group, data = plot.df, main = b, sub = s.t)
}
dev.off()
par(old.par)


message("DONE")
