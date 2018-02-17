#' For the argument \code{genes} construct a vector space based on the distinct
#' annotations found for the respective genes in arument \code{annot.df}. Each
#' distinct annotation will be an axis, each gene will be assigned a vector,
#' with coordinates the discrete number of times an annotation has been
#' assigned to the respective gene.
#'
#' @param genes A character vector of gene IDs
#' @param annot.df A data.frame or matrix holding the gene annotations. Default
#' is \code{all.ipr}.
#' @param gene.col The column name or number of \code{annot.df} in which to
#' look up the gene IDs found in \code{genes}. Default is \code{1}.
#' @param annot.col The column name or number of \code{annot.df} in which to
#' look up the genes annotations. Default is \code{2}.
#'
#' @export
#' @return \code{NA} if no annotations are found for any gene in \code{genes},
#' otherwise a data.frame with a column \code{gene} in which the gene IDs are
#' held and one column for each distinct annotation or the respective vector
#' space axes.
constructAnnotationVectorSpace <- function(genes, annot.df = all.ipr, gene.col = 1, 
    annot.col = 2) {
    genes.w.annots <- sort(intersect(genes, annot.df[, gene.col]))
    if (length(genes.w.annots) > 0) {
        genes.a.df <- annot.df[which(annot.df[, gene.col] %in% genes.w.annots), 
            ]
        uniq.annots <- sort(unique(genes.a.df[, annot.col]))
        cbind(data.frame(gene = genes.w.annots, stringsAsFactors = FALSE), 
            Reduce(cbind, lapply(uniq.annots, function(a) {
                g.a <- genes.a.df[which(genes.a.df[, annot.col] == a), 
                  gene.col]
                g.a.col <- rep(0, length(genes.w.annots))
                if (length(g.a) > 0) 
                  g.a.col[which(genes.w.annots %in% g.a)] <- 1
                as.data.frame(matrix(g.a.col, ncol = 1, dimnames = list(c(), 
                  a)), stringsAsFactors = FALSE)
            })))
    } else NA
}


#' Invokes \code{GeneFamilies::statVectorCloud} on each subset of genes in
#' \code{annot.vec.space} classified by the gene-sets in \code{gene.sets}.
#'
#' @param annot.vec.space A data.frame result of calling function
#' \code{constructAnnotationVectorSpace}
#' @param gene.sets A named list, where names indicate gene-sets (or classes),
#' and values are gene IDs including the ones references in
#' \code{annot.vec.space}
#'
#' @export
#' @return A named list, with one entry per gene-class in \code{gene.sets}.
#' Values are the results of invoking \code{GeneFamilies::statVectorCloud} on
#' each respective subset.
geneClassesVectorClouds <- function(annot.vec.space, gene.sets) {
    # Validate a
    ument
    if (is.null(annot.vec.space) || is.na(annot.vec.space)) 
        return(NA)
    
    Reduce(append, lapply(names(gene.sets), function(g.s.n) {
        set.genes <- which(annot.vec.space[, "gene"] %in% gene.sets[[g.s.n]])
        if (length(set.genes) > 0) {
            stat.vec.cl <- statVectorCloud(annot.vec.space[set.genes, setdiff(colnames(annot.vec.space), 
                "gene"), drop = FALSE])
            setNames(list(stat.vec.cl), g.s.n)
        } else NULL
    }))
}


#' Identifies instances of subfunctionalization occurred in the evolution from
#' the base genes to the evolved genes. Subfunc. is identified based on the
#' gene annotations. All possible combinations of evolved genes are tested with
#' each base gene. Subfunc conditions are 1. The tested evolved genes
#' annotations (EA) are not a superset of the base gene annotations (BGA), 2.
#' the intersection of EA is smaller than the union of EA, and 3. The union of
#' EA is identical to the BGA, or, if not strict, the union of the EA is a
#' subset of the BGA.
#'
#' @param base.genes Character vector of gene IDs, identifying the base genes.
#' @param evolved.genes Character vector of gene IDs, identifying the evolved
#' genes.
#' @param annot.df A data.frame of gene annotations. Default is
#' \code{GeneFamilies::all.ipr}
#' @param gene.col The column name or index in \code{annot.df} in which to
#' lookup the gene IDs. Default is \code{1}
#' @param anno.col The column name or index in \code{annot.df} in which to
#' lookup the gene annotations. Default is \code{2}.
#' @param strict Boolean indicating whether to apply the strict criterium in
#' the second condition, or not (see above for details). Default is
#' \code{TRUE}.
#'
#' @export
#' @return \code{NA} if there are no annotations for at least one base gene and
#' at least two evolved.genes. \code{NULL} is returned if the previous
#' conditions are met, but no subfunctionilizations could be identified, and
#' finally a data.frame is returned with the identified subfunctionilizations.
#' In this, columns are: base.gene, base.annotations, evolved.genes, and
#' evolved.annotations.
subfunctionalizationAnnotationBased <- function(base.genes, evolved.genes, 
    annot.df = all.ipr, gene.col = 1, annot.col = 2, strict = TRUE) {
    # At least a single base genes and at least two evolved genes need to
    # have annotations:
    if (any(base.genes %in% annot.df[, gene.col]) && length(which(evolved.genes %in% 
        annot.df[, gene.col])) > 1) {
        # Compound annotations for base and evolved genes, and identification
        # of duplicated annotations:
        base.genes.annots <- setNames(lapply(base.genes, function(x) sort(annot.df[which(annot.df[, 
            gene.col] == x), annot.col])), base.genes)
        b.g.a.duplicated <- duplicated(base.genes.annots)
        evolved.genes.annots <- setNames(lapply(evolved.genes, function(x) sort(annot.df[which(annot.df[, 
            gene.col] == x), annot.col])), evolved.genes)
        e.g.a.duplicated <- duplicated(evolved.genes.annots)
        # For each base gene, test different combinations of evolved genes for
        # possible subfunctionalization:
        Reduce(rbind, lapply(names(base.genes.annots[!b.g.a.duplicated]), 
            function(b.g) {
                b.g.a <- annot.df[which(annot.df[, gene.col] == b.g), annot.col]
                # Return:
                res.df <- NULL
                # A single annotation can not be subfunctionalized:
                if (length(b.g.a) > 1) {
                  e.g <- names(evolved.genes.annots[!e.g.a.duplicated][sapply(evolved.genes.annots[!e.g.a.duplicated], 
                    function(x) {
                      # Subfunctionalization can only have occurred, if evolved genes have a
                      # TRUE subset of the base gene's annotations:
                      cond.1 <- length(intersect(b.g.a, x)) > 0
                      cond.2 <- length(setdiff(b.g.a, x)) > 0
                      cond.1 && cond.2
                    })])
                  # We need at least two evolved genes to have subfunctionilization:
                  if (length(e.g) > 1) {
                    # Init while loop:
                    keep.testing <- TRUE
                    m <- length(e.g)
                    e.g.i.s <- combn(e.g, m, simplify = FALSE)
                    i <- 1
                    # Test combinations of evolved genes for subfunctionalization:
                    while (keep.testing) {
                      # Init current iteration:
                      e.g.i <- e.g.i.s[[i]]
                      e.g.a <- evolved.genes.annots[!e.g.a.duplicated][e.g.i]
                      e.g.a.intersect <- Reduce(intersect, e.g.a)
                      e.g.a.union <- Reduce(union, e.g.a)
                      # TEST two conditions for subfunctionalization: 1.  Evolved genes do
                      # not have identical annotations,
                      cond.e.g <- length(e.g.a.intersect) < length(e.g.a.union)
                      # and 2. The union of the evolved genes annotations is identical to the
                      # base genes annotations (strict), or at least a subset of the base
                      # genes annotations (not strict).
                      cond.b.g <- if (strict) 
                        identical(b.g.a, e.g.a.union) else all(b.g.a %in% e.g.a.union)
                      if (cond.e.g && cond.b.g) {
                        # Found subfunctionalization, construct result data.frame:
                        df.b.g <- if (!any(b.g.a.duplicated)) {
                          b.g
                        } else {
                          paste(sort(names(base.genes.annots[which(sapply(base.genes.annots, 
                            function(x) identical(b.g.a, x)))])), collapse = ",")
                        }
                        df.b.a <- paste(b.g.a, collapse = ",")
                        if (!any(e.g.a.duplicated)) {
                          df.e.g <- paste(sort(e.g.i), collapse = ",")
                          df.e.a <- paste(sapply(sort(e.g.i), function(x) paste(x, 
                            "(", paste(evolved.genes.annots[!e.g.a.duplicated][[x]], 
                              collapse = ","), ")", sep = "")), collapse = ", ")
                        } else {
                          e.g.compl <- sort(union(e.g.i, names(evolved.genes.annots[e.g.a.duplicated])[which(sapply(evolved.genes.annots[e.g.a.duplicated], 
                            function(e.g.dupl.annot) {
                              any(sapply(e.g.a, function(e.g.annot) {
                                identical(e.g.annot, e.g.dupl.annot)
                              }))
                            }))]))
                          df.e.g <- paste(e.g.compl, collapse = ",")
                          df.e.a <- paste(sapply(e.g.compl, function(x) paste(x, 
                            "(", paste(evolved.genes.annots[[x]], collapse = ","), 
                            ")", sep = "")), collapse = ", ")
                        }
                        res.df <- data.frame(base.genes = df.b.g, base.annotations = df.b.a, 
                          evolved.genes = df.e.g, evolved.annotations = df.e.a, 
                          stringsAsFactors = FALSE)
                        # Do not look for a smaller set of evolved genes, that fullfill the
                        # subfunctionalization criteria:
                        keep.testing <- FALSE
                      }
                      if (keep.testing) {
                        if (m == 2 && i == length(e.g.i.s)) {
                          # Tested all possible combinations of evolved genes.
                          keep.testing <- FALSE
                        } else {
                          # Prepare next iteration
                          if (i < length(e.g.i.s)) {
                            i <- i + 1
                          } else {
                            m <- m - 1
                            i <- 1
                            e.g.i.s <- combn(e.g, m, simplify = FALSE)
                          }
                        }
                      }
                    }
                  }
                  res.df  # return value
                }
            }))
    } else NA
}


#' Testing \code{subfunctionalizationAnnotationBased}
#'
#' @export
#' @return \code{TRUE} if and only of all tests are passed successfully.
testSubfunctionalizationAnnotationBased <- function() {
    base.genes <- LETTERS[1:3]
    evolved.genes <- letters[1:5]
    annot.df <- read.table(stringsAsFactors = FALSE, text = "
A dom_1
A dom_2
B dom_2
B dom_3
C dom_1
C dom_2
a dom_1
b dom_2
c dom_2
c dom_3
d dom_4
d dom_5
e dom_1
e dom_3")
    # Test non strict:
    res.df.1 <- read.table(stringsAsFactors = FALSE, header = TRUE, quote = "\"", 
        colClasses = rep("character", 4), sep = ";", text = "\"base.genes\";\"base.annotations\";\"evolved.genes\";\"evolved.annotations\"
\"A,C\";\"dom_1,dom_2\";\"a,b,c,e\";\"a(dom_1), b(dom_2), c(dom_2,dom_3), e(dom_1,dom_3)\"
\"B\";\"dom_2,dom_3\";\"b,e\";\"b(dom_2), e(dom_1,dom_3)\"")
    y.1 <- subfunctionalizationAnnotationBased(base.genes, evolved.genes, 
        annot.df, strict = FALSE)
    test.1 <- identical(res.df.1, y.1)
    # Test strict:
    res.df.2 <- read.table(stringsAsFactors = FALSE, header = TRUE, quote = "\"", 
        colClasses = rep("character", 4), sep = ";", text = "\"base.genes\";\"base.annotations\";\"evolved.genes\";\"evolved.annotations\"
\"A,C\";\"dom_1,dom_2\";\"a,b\";\"a(dom_1), b(dom_2)\"")
    test.2 <- identical(res.df.2, subfunctionalizationAnnotationBased(base.genes, 
        evolved.genes, annot.df, strict = TRUE))
    # Test strict with evolved genes with identical annotations:
    annot.df.2 <- rbind(annot.df, data.frame(V1 = "f", V2 = "dom_1", stringsAsFactors = FALSE))
    res.df.3 <- read.table(stringsAsFactors = FALSE, header = TRUE, quote = "\"", 
        colClasses = rep("character", 4), sep = ";", text = "\"base.genes\";\"base.annotations\";\"evolved.genes\";\"evolved.annotations\"
\"A,C\";\"dom_1,dom_2\";\"a,b,f\";\"a(dom_1), b(dom_2), f(dom_1)\"")
    evolved.genes.2 <- c(evolved.genes, "f")
    y.3 <- subfunctionalizationAnnotationBased(base.genes, evolved.genes.2, 
        annot.df.2, strict = TRUE)
    test.3 <- identical(res.df.3, y.3)
    # Return:
    all(c(test.1, test.2, test.3))
}


#' Invokes \code{subfunctionalizationAnnotationBased} for all gene.groups
#' (families) contained in argument \code{gene.groups.lst}, using argument
#' \code{gene.classes} to separate 'base genes' from 'evolved genes'.
#'
#' @param gene.groups.lst A named list of character vectors holding gene IDs.
#' @param gene.classes A named list of gene identifiers separated into a base
#' class and a class of evolved genes.
#'
#' @export
#' @return NULL, if no subfunctionalization could be identified within any of
#' the gene families, or a data frame with the identified subfunctionalizations
#' (see \code{subfunctionalizationAnnotationBased} for details). An additional
#' column holds the name of the gene family.
geneGroupSubfunctionalization <- function(gene.groups.lst, gene.classes) {
    Reduce(rbind, mclapply(names(gene.groups.lst), function(gene.grp.nm) {
        genes <- gene.groups.lst[[gene.grp.nm]]
        base.genes <- intersect(genes, tand.classifier$ortholog)
        evolved.genes <- intersect(genes, tand.classifier$tandem)
        res.df <- subfunctionalizationAnnotationBased(base.genes, evolved.genes)
        if (!is.null(res.df) && !is.na(res.df)) {
            res.df$family <- gene.grp.nm
            res.df
        } else NULL
    }))
}


#' Measures annotations gained, lost, and shared between base genes and evolved
#' genes.
#'
#' @param genes A group of genes in the form of a character vector of gene IDs
#' @param gene.sets A named list defining base gene and evolved gene
#' identifiers.
#' @param base.gene.set The name of list-entry in \code{gene.sets} indicating
#' which is to be interpreted as base genes. Default is \code{'ortholog'}.
#' @param annot.df A data.frame of gene annotations. Default is
#' \code{GeneFamilies::all.ipr}
#' @param gene.col The column name or index in \code{annot.df} in which to
#' lookup the gene IDs. Default is \code{1}
#' @param anno.col The column name or index in \code{annot.df} in which to
#' lookup the gene annotations. Default is \code{2}.
#'
#' @export
#' @return A data.frame with the following columns: n.gained, n.lost, n.shared,
#' and total.
annotationEvolution <- function(genes, gene.sets, base.gene.set = "ortholog", 
    annot.df = all.ipr, gene.col = 1, annot.col = 2) {
    g.i <- lapply(gene.sets, function(x) which(annot.df[, gene.col] %in% 
        intersect(genes, x)))
    if (all(sapply(g.i, function(x) length(x) > 0))) {
        gene.sets.annot <- setNames(lapply(names(g.i), function(gene.class) {
            sort(unique(annot.df[g.i[[gene.class]], annot.col]))
        }), names(g.i))
        base.annot <- gene.sets.annot[[base.gene.set]]
        Reduce(rbind, lapply(setdiff(names(gene.sets.annot), base.gene.set), 
            function(gene.set.name) {
                g.s.a <- gene.sets.annot[[gene.set.name]]
                n.gained <- length(setdiff(g.s.a, base.annot))
                n.lost <- length(setdiff(base.annot, g.s.a))
                n.shared <- length(intersect(base.annot, g.s.a))
                total <- length(union(base.annot, g.s.a))
                data.frame(n.gained = n.gained, n.lost = n.lost, n.shared = n.shared, 
                  total = total, stringsAsFactors = FALSE)
            }))
    } else NA
}
