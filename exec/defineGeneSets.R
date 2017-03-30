require(GeneFamilies)
options( mc.cores=detectCores() )

message("USAGE: Rscript path/2/GeneFamilies/exec/defineGeneSets.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Define sets of genes and gene families showing signs of positive selection:
fubar.fams.decisive.evidence.genes <- unique(unlist(families.lst[fubar.fams.decisive.evidence]))
busted.fams <- unique(busted.df[which(busted.df$p.value <= 0.05), "family"])
busted.genes <- unique(busted.df[which(busted.df$p.value <= 0.05), "gene"])
fams.pos.sel <- union(fubar.fams.decisive.evidence, busted.fams)
genes.pos.sel <- union(fubar.fams.decisive.evidence.genes, busted.genes)


#' Define sets of genes and gene families showing species specific sign of
#' significant expansion or contraction:
families.exp.cnt <- rownames(cafe.result.df)[which(apply(cafe.result.df, 1, function(x) any(unlist(x) <= 
    0.05)))]
families.exp.cnt.genes <- unique(unlist(families.lst[families.exp.cnt]))
families.spec.cols <- c("A.thaliana", "A.lyrata", "C.rubella", "C.hirsuta", "A.arabicum", 
    "B.rapa", "E.salsugineum", "S.parvula")
families.unq <- families.df$id[which(apply(families.df[, families.spec.cols], 1, 
    function(x) {
        length(x[which(x > 0)]) == 1 && sum(x) > 1
    }))]
col.maps <- setNames(sub("CAFE.", "", colnames(cafe.result.df), fixed = TRUE), colnames(cafe.result.df))
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
families.exp.cnt.unq <- union(families.exp.cnt, families.unq)
families.exp.cnt.unq.genes <- unlist(families.lst[families.exp.cnt.unq])
families.exp.genes <- unlist( families.lst[ families.exp ] )


#' Conserved families are those that have an identical number of genes within
#' each species:
families.conserved <- families.df$id[which(families.df$size%%length(families.spec.cols) == 
    0)]
families.conserved.genes <- unlist(families.lst[families.conserved])


#' Generate data containers required for the computation of annotation based
#' function diversity (Shannon-Entropies).
#' - Tandem Duplicates:
tandem.nms <- unique(tandems$Family)
tandems.lst <- setNames(mclapply(tandem.nms, function(x) tandems[which(tandems$Family == 
    x), "Gene"]), tandem.nms)
tandems.genes <- unlist(tandems.lst)
#' - Orthologs:
orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) unlist(orthologs[x, 
    ])), orths.nms)
orthologs.genes <- unlist(orthologs.lst)


#' Plot results:
#' - Venn Diagram of the five basic gene groups:
pdf(file.path(input.args, "inst", "GeneGroupsVenn.pdf"))
lst <- list(Exp = families.exp.genes, Orth = orthologs.genes, Tand = tandems.genes, 
    PSel = genes.pos.sel)
venn(lst)
dev.off()


#' Save results:
save(fubar.fams.decisive.evidence.genes, busted.fams, busted.genes, fams.pos.sel, 
    genes.pos.sel, families.exp.cnt, families.exp.cnt.genes, families.spec.cols, 
    families.unq, families.unq.genes, families.exp.cnt.unq, families.exp.cnt.unq.genes, families.exp, families.exp.genes,
    families.conserved, families.conserved.genes, tandems.lst, tandems.genes, orthologs.lst, 
    orthologs.genes, file = file.path(input.args[[1]], "data", "GeneGroups.RData"))


message("DONE")
