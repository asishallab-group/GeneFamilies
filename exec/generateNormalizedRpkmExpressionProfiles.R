require(MaizeGeneFamilies)
#' Load pairwise sequence similarity results (all.vs.all.sim):
data("pairwiseSequenceSimilarities", package = "MaizeGeneFamilies")

options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/exec/generateNormalizedRpkmExpressionProfiles.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Normalize the RPKM expression vectors, maintaining their direction but
#' scaling them up to length 1.0:
expr.cols <- c("cotyledon", "developing leaf", "seedling", "flower stage 9", 
    "flower stage 16")
message("Expression Vector Space defined by the following axes:")
message(paste(expr.cols, collapse = ", "))
rpkm.expr.profiles.norm.df <- merge(rpkm.expr.profiles.df[, c("gene", "gene.exp.var")], 
    Reduce(rbind, mclapply(1:nrow(rpkm.expr.profiles.df), function(i) {
        x <- rpkm.expr.profiles.df[i, ]
        x[, expr.cols] <- normalizeVector(x[, expr.cols])
        x[, c("gene", expr.cols)]
    })), by = "gene")


#' Measure each expression profile's angle to the vector space diagonal:
rpkm.expr.profiles.norm.df$angleDiag <- as.numeric(apply(rpkm.expr.profiles.df[, 
    expr.cols], 1, function(x) rad2deg(acos(cosDiag(x)))))


#' Measure each expression profile's distance to the diagonal vector:
norm.diag <- normalizeVector(rep(1, length(expr.cols)))
rpkm.expr.profiles.norm.df$distDiag <- as.numeric(apply(rpkm.expr.profiles.norm.df[, 
    expr.cols], 1, function(x) euclNorm(x - norm.diag)))


#' Compute distances between normalized expression vectors for each gene and
#' its closest ortholog, if present:
tands.expr <- intersect(rpkm.expr.profiles.norm.df$gene, sub("\\.\\d+$", 
    "", tandems.genes))
orths.expr <- intersect(rpkm.expr.profiles.norm.df$gene, orthologs.genes)
a.v.a.s <- all.vs.all.sim
a.v.a.s$V1 <- sub("\\.\\d+$", "", a.v.a.s$V1)
a.v.a.s$V2 <- sub("\\.\\d+$", "", a.v.a.s$V2)
a.v.a.s <- a.v.a.s[which(a.v.a.s$V1 %in% setdiff(rpkm.expr.profiles.norm.df$gene, 
    orths.expr) & a.v.a.s$V2 %in% orths.expr), ]


#' Infer each gene's closest ortholog in terms of sequence similarity:
rpkm.expr.profiles.norm.df$closest.ortholog <- as.character(mclapply(rpkm.expr.profiles.norm.df$gene, 
    function(x) {
        closestHomolog(x, seq.sim.tbl = a.v.a.s)
    }))
#'Compute difference in angles to diagonal and distance between a gene and its
#'closest ortholog, where possible:
rpkm.expr.profiles.norm.df <- merge(rpkm.expr.profiles.norm.df, Reduce(rbind, 
    mclapply(1:nrow(rpkm.expr.profiles.norm.df), function(i) {
        x <- rpkm.expr.profiles.norm.df[i, ]
        if (x$closest.ortholog %in% rpkm.expr.profiles.norm.df$gene) {
            rpkm.g <- rpkm.expr.profiles.norm.df[which(rpkm.expr.profiles.norm.df$gene == 
                x$gene), ]
            rpkm.o <- rpkm.expr.profiles.norm.df[which(rpkm.expr.profiles.norm.df$gene == 
                x$closest.ortholog), ]
            x$dist.closest.orth <- euclNorm(rpkm.g[, expr.cols] - rpkm.o[, 
                expr.cols])
            x$diff.angle.diag <- rpkm.g[, "angleDiag"] - rpkm.o[, "angleDiag"]
            x[, c("gene", "dist.closest.orth", "diff.angle.diag")]
        } else NULL
    })), by = "gene", all.x = TRUE)


#' Save results:
save(rpkm.expr.profiles.norm.df, file = file.path(input.args[[1]], "data", 
    "rpkmNormalizedExpressionProfiles.RData"))


message("DONE")
