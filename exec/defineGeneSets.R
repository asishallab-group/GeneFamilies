require(GeneFamilies)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/defineGeneSets.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Define sets of genes and gene families showing signs of positive selection:
fubar.fams.decisive.evidence.genes <- unique(unlist(families.lst[fubar.fams.decisive.evidence]))
fubar.fams.decisive.evidence.genes.nev <- removeExpressionVariant(fubar.fams.decisive.evidence.genes)
busted.fams <- unique(busted.df[which(busted.df$p.value <= 0.05), "family"])
busted.genes <- unique(busted.df[which(busted.df$p.value <= 0.05), "gene"])
busted.genes.nev <- removeExpressionVariant(busted.genes)
fams.pos.sel <- union(fubar.fams.decisive.evidence, busted.fams)
genes.pos.sel <- union(fubar.fams.decisive.evidence.genes, busted.genes)
genes.pos.sel.nev <- removeExpressionVariant(genes.pos.sel)


#' Define sets of genes and gene families showing species specific sign of
#' significant expansion or contraction:
families.exp.cnt <- rownames(cafe.result.df)[which(apply(cafe.result.df, 
    1, function(x) any(unlist(x) <= 0.05)))]
families.exp.cnt.genes <- unique(unlist(families.lst[families.exp.cnt]))
families.exp.cnt.genes.nev <- removeExpressionVariant(families.exp.cnt.genes)
families.spec.cols <- c("A.thaliana", "A.lyrata", "C.rubella", "C.hirsuta", 
    "A.arabicum", "B.rapa", "E.salsugineum", "S.parvula")
families.unq <- families.df$id[which(apply(families.df[, families.spec.cols], 
    1, function(x) {
        length(x[which(x > 0)]) == 1 && sum(x) > 1
    }))]
col.maps <- setNames(sub("CAFE.", "", colnames(cafe.result.df), fixed = TRUE), 
    colnames(cafe.result.df))
families.exp <- rownames(cafe.result.df)[as.logical(unlist(mclapply(rownames(cafe.result.df), 
    function(x) {
        y <- cafe.result.df[x, ]
        i <- which(y <= 0.05)
        res <- FALSE
        if (length(i) > 0) {
            i.spec <- col.maps[i]
            z <- families.df[which(families.df$id == x), ]
            res <- any(z[, i.spec] > median(unlist(z[, families.spec.cols])))
        }
        res
    })))]
families.unq.genes <- unlist(families.lst[families.unq])
families.unq.genes.nev <- removeExpressionVariant(families.unq.genes)
families.exp.cnt.unq <- union(families.exp.cnt, families.unq)
families.exp.cnt.unq.genes <- unlist(families.lst[families.exp.cnt.unq])
families.exp.cnt.unq.genes.nev <- removeExpressionVariant(families.exp.cnt.unq.genes)
families.exp.genes <- unlist(families.lst[families.exp])
families.exp.genes.nev <- removeExpressionVariant(families.exp.genes)


#' Conserved families are those that have an identical number of genes within
#' each species:
families.conserved <- families.df$id[which(families.df$size%%length(families.spec.cols) == 
    0)]
families.conserved.genes <- unlist(families.lst[families.conserved])
families.conserved.genes.nev <- removeExpressionVariant(families.conserved.genes)


#' Generate data containers required for the computation of annotation based
#' function diversity (Shannon-Entropies).
#' - Tandem Duplicates:
tandem.nms <- unique(tandems$Family)
tandems.lst <- setNames(mclapply(tandem.nms, function(x) tandems[which(tandems$Family == 
    x), "Gene"]), tandem.nms)
tandems.genes <- unlist(tandems.lst)
tandems.genes.nev <- removeExpressionVariant(tandems.genes)
#' - Orthologs:
orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) unlist(orthologs[x, 
    ])), orths.nms)
orthologs.genes <- unlist(orthologs.lst)


#' Define Duplicated Gene Sets with Orthologs:
t.i <- which(as.logical(lapply(tandems.lst, function(t.c) any(removeExpressionVariant(t.c) %in% 
    orthologs.genes))))
tands.w.orths <- lapply(tandems.lst[t.i], function(t.c) removeExpressionVariant(t.c))
d.i <- which(as.logical(lapply(families.lst[families.exp.cnt.unq], function(f.e) any(removeExpressionVariant(f.e) %in% 
    orthologs.genes))))
dupl.w.orths <- lapply(families.lst[families.exp.cnt.unq][d.i], function(f.e) removeExpressionVariant(f.e))
p.i <- which(as.logical(lapply(families.lst[fams.pos.sel], function(f.p) any(removeExpressionVariant(f.p) %in% 
    orthologs.genes))))
psel.w.orths <- lapply(families.lst[fams.pos.sel][p.i], function(f.p) removeExpressionVariant(f.p))

#' Define lists that can be used to classify above groups into subgroups:
tand.classifier <- list(ortholog = orthologs.genes, tandem = setdiff(unlist(tands.w.orths), 
    orthologs.genes))
tand.psel.classifier <- list(ortholog = orthologs.genes, tandem = setdiff(unlist(tands.w.orths), 
    union(orthologs.genes, genes.pos.sel.nev)), pos.selected = setdiff(intersect(unlist(tands.w.orths), 
    genes.pos.sel.nev), orthologs.genes))
dupl.classifier <- list(ortholog = orthologs.genes, duplicated = setdiff(unlist(dupl.w.orths), 
    orthologs.genes))
dupl.psel.classifier <- list(ortholog = orthologs.genes, duplicated = setdiff(unlist(dupl.w.orths), 
    union(genes.pos.sel.nev, orthologs.genes)), pos.selected = setdiff(intersect(unlist(dupl.w.orths), 
    genes.pos.sel.nev), orthologs.genes))
psel.classifier <- list(ortholog = orthologs.genes, pos.selected = setdiff(unlist(psel.w.orths), 
    orthologs.genes))


#' Define sets of expressed genes:
tands.expr <- intersect(tandems.genes, rpkm.expr.profiles.df$gene.exp.var)
tands.expr.nev <- removeExpressionVariant(tands.expr)
orths.expr <- setdiff(intersect(orthologs.genes, rpkm.expr.profiles.df$gene), 
    removeExpressionVariant(tands.expr))
dupl.expr <- setdiff(intersect(families.exp.genes, rpkm.expr.profiles.df$gene.exp.var), 
    union(tands.expr, rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$gene %in% 
        orths.expr), "gene.exp.var"]))
dupl.expr.nev <- removeExpressionVariant(dupl.expr)
psel.expr <- Reduce(intersect, list(families.exp.genes, rpkm.expr.profiles.df$gene.exp.var, 
    genes.pos.sel))
psel.expr.nev <- removeExpressionVariant(psel.expr)


#' Plot results:
#' - Venn Diagram of the five basic gene groups:
pdf(file.path(input.args, "inst", "GeneGroupsVenn.pdf"))
lst <- list(Dupl = families.exp.genes.nev, Orth = orthologs.genes, Tand = tandems.genes.nev, 
    PSel = genes.pos.sel.nev)
venn(lst)
dev.off()

#' - Venn Diagram of groups of expressed genes:
pdf(file.path(input.args, "inst", "ExpressedGeneGroupsVenn.pdf"))
lst.exp <- list(Dupl = dupl.expr.nev, Orth = orths.expr, Tand = tands.expr.nev, 
    PSel = psel.expr.nev)
venn(lst.exp)
dev.off()


#' Save results:
save(fubar.fams.decisive.evidence.genes, fubar.fams.decisive.evidence.genes.nev, 
    busted.fams, busted.genes, busted.genes.nev, fams.pos.sel, genes.pos.sel, 
    genes.pos.sel.nev, families.exp.cnt, families.exp.cnt.genes, families.exp.cnt.genes.nev, 
    families.spec.cols, families.unq, families.unq.genes, families.unq.genes.nev, 
    families.exp.cnt.unq, families.exp.cnt.unq.genes, families.exp.cnt.unq.genes.nev, 
    families.exp, families.exp.genes, families.exp.genes.nev, families.conserved, 
    families.conserved.genes, families.conserved.genes.nev, tandems.lst, 
    tandems.genes, tandems.genes.nev, orthologs.lst, orthologs.genes, tands.w.orths, 
    dupl.w.orths, psel.w.orths, tand.classifier, tand.psel.classifier, 
    dupl.classifier, dupl.psel.classifier, psel.classifier, file = file.path(input.args[[1]], 
        "data", "GeneGroups.RData"))


message("DONE")
