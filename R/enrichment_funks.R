#' Generates a matrix that can directly be used as input for hypothesis
#' testing, e.g. with base::fisher.test(...). Validates the resulting table
#' based on the input data. Validation checks whether the intersections of
#' elements having attribute A and B and those having attribute A and NOT B
#' truly form the whole set of all elements with attribute A. The same
#' validation is done for the other rows and columns. Error-Codes are:
#' 1 - row.1 does not fully separate into col.1 and col.2
#' 2 - row.2 does not fully separate into col.1 and col.2
#' 3 - col.1 does not fully separate into row.1 and row.2
#' 4 - col.2 does not fully separate into row.1 and row.2
#' 5 - the sum of all cells sizes in the contingency table is not identical to
#' the size of the mathematical set of all elements contained in both rows and
#' columns (row.1, row.2, col.1, col.2).
#'
#' @param row.1 a vector of elements that have attribute A
#' @param row.2 a vector of elements that do not have attribute A
#' @param col.1 a vector of elements that have attribute B
#' @param col.2 a vector of elements that do not have attribute B
#' @param row.category the name of attribute A true for elements in row.1 and
#' false for elements in row.2
#' @param col.category the name of attribute B true for elements in col.1 and
#' false for elements in col.2
#' @param categories default is 'T' and 'F' indicating whether row.1 and col.1
#' indicate elements that have attribute A, and elements in row.2 and col.2 do
#' not have the respective attribute.
#' @param validate set to FALSE, if validation is not wanted. Default is TRUE
#'
#' @return A 2 x 2 contingency table.
#' @export
generateContingencyTable <- function(row.1, row.2, col.1, col.2, row.category, col.category, 
    categories = c("T", "F"), validate = TRUE) {
    # Generates a contingency table of the given argument rows and columns. The
    # resulting table can be used for e.g. exact Fisher tests and the like.
    cont.tbl <- matrix(c(length(intersect(row.1, col.1)), length(intersect(row.2, 
        col.1)), length(intersect(row.1, col.2)), length(intersect(row.2, col.2))), 
        nrow = 2, ncol = 2, dimnames = setNames(list(categories, categories), c(row.category, 
            col.category)))
    if (validate) {
        val.res <- c(sum(cont.tbl[1, ]) == length(row.1), sum(cont.tbl[2, ]) == length(row.2), 
            sum(cont.tbl[, 1]) == length(col.1), sum(cont.tbl[, 2]) == length(col.2), 
            sum(cont.tbl) == length(unique(c(row.1, row.2, col.1, col.2))))
        if (!all(val.res)) 
            stop("Validation of contingency table failed. Error-Code(s): ", paste(which(!val.res), 
                collapse = ","))
    }
    cont.tbl
}

#' Removes trailing splice variant identifier from gene IDs; e.g. CARHR213450.1
#' will result in CARHR213450 
#'
#' @param gene.ids a character vector of gene identifiers
#'
#' @return The modified 'gene.ids' without the trailing splice variant IDs.
#' @export
removeSpliceVariant <- function(gene.ids) {
    sub("\\.\\d+$", "", gene.ids)
}

#' Function appends a row to the argument data.frame 'df', which is supposed to
#' hold two columns 'tested.hypothesis' and 'p.val'.
#'
#' @param df an instance of class base::data.frame with two columns
#' @param hypo string describing the hypothesis tested
#' @param p.val numerical result of testing 'hypo' as alternative hypothesis.
#'
#' @return The extended data.frame 'df' with a new row holding the arguments
#' 'hypo' and 'p.val'.
#' @export
appendToTestedHypotheses <- function(df, hypo, p.val) {
    a.df <- data.frame(tested.hypothesis = hypo, p.value = p.val, stringsAsFactors = FALSE)
    rbind(df, a.df)
}

#' Counts the number of times an annotation, e.g. a Pfam or InterPro entry, is
#' annotated in argument annotation table 'anno.tbl'.
#'
#' @param anno string contained in column 2 of 'anno.tbl'
#' @param anno.tbl a data.frame of two columns, column 1 holds gene identifiers
#' and column 2 annotations, e.g. InterPro or Pfam entries.
#'
#' @return integer - the number of genes 'anno' has been found annotated to.
#' @export
numberAnnotations <- function(anno, anno.tbl) {
    nrow(anno.tbl[which(anno.tbl$V2 == anno), ])
}

#' Generate a 2x2 contingency table for statistical hypothesis testing of over-
#' or underrepresentation of 'anno' in the case 'anno.tbl' when compared with
#' the whole universe 'universe.tbl'. Please note, that 'anno.tbl' is expected
#' to be a true subset of 'universe.tbl.'
#'
#' @param anno string representing the annotation, e.g. InterPro entry
#' @param anno.tbl annotation table: column 1 holds gene identifiers, and
#' column 2 annotations, e.g. InterPro entries. 'anno.tbl' represents the
#' 'case'
#' @param universe.tbl annotation table, true superset of 'anno.tbl'.
#' Represents the 'universe' to test alternative hypotheses against.
#' @param anno.tbl.anno.col The column of 'anno.tbl' in which to lookup the
#' case gene identifiers. Default is 1.
#' @param universe.tbl.anno.col The column of 'universe.tbl' in which to lookup
#' the universe gene identifiers. Default is 1.
#' @param dim.names used to set the row and column names of the generated
#' contingency table. See base::matrix(...,dimnames) for more details. Default
#' value is 'dim.names=list( c( 'annotated', 'not.annotated' ), c( 'case',
#' 'universe' ) )'
#'
#' @return An instance of base::matrix representing the contingency table 
#' @export
annotationMatrix <- function(anno, anno.tbl, universe.tbl, anno.tbl.anno.col = 1, 
    universe.tbl.anno.col = 1, dim.names = list(c("annotated", "not.annotated"), 
        c("case", "universe"))) {
    univ.min.case <- universe.tbl[which(!(universe.tbl[, universe.tbl.anno.col] %in% 
        anno.tbl[, anno.tbl.anno.col])), ]
    no.annos.case <- numberAnnotations(anno, anno.tbl)
    no.annos.universe <- numberAnnotations(anno, univ.min.case)
    matrix(c(no.annos.case, (nrow(anno.tbl) - no.annos.case), no.annos.universe, 
        (nrow(univ.min.case) - no.annos.universe)), ncol = 2, dimnames = dim.names)
}

#' Tests for every function annotation in 'case.anno.tbl' the alternative
#' hypothesis provided with 'test.dir'. Resulting p.values are adjusted for
#' multiple hypothesis testing, before being returned. Note, if applied on
#' large data sets, run in parallel on multiple cores and set
#' options('mc.cores'=...) accordingly.
#'
#' @param case.anno.tbl the annotations for which the alternative hypothesis is
#' to be tested
#' @param universe.anno.tbl the 'universe' annotations to compare the case with
#' @param test.dir the direction of the alternative hypothesis. Default defines
#' it as 'greater'
#' @param p.adjust.method the method used to adjust resulting p.values for
#' multiple hypothesis testing. Default is 'BY', see base::p.adjust for more
#' details.
#'
#' @return A data.frame holding the annotation in column 'PFam', the adjusted
#' and un-adjusted P.Value in two additional columns.
#' @export
enrichmentTests <- function(case.anno.tbl, universe.anno.tbl, test.dir = "greater", 
    p.adjust.method = "BY") {
    case.uniq.annos <- unique(case.anno.tbl$V2)
    res <- as.data.frame(cbind(PFam = case.uniq.annos, P.Value = as.numeric(mclapply(case.uniq.annos, 
        function(anno) {
            fisher.test(annotationMatrix(anno, case.anno.tbl, universe.anno.tbl), 
                alternative = test.dir)$p.value
        }))), stringsAsFactors = FALSE)
    res$P.Value <- as.numeric(res$P.Value)
    cbind(res, P.Value.adjusted = p.adjust(res$P.Value, method = p.adjust.method))
}

#' Filters 'enrichment.tbl' for those annotations that have a significant
#' adjusted p.value.
#'
#' @param enrichment.tbl a data.frame returned by function 'enrichmentTests'
#' @param significance.level the p.value cutoff to apply, default is '<= 0.05'
#'
#' @return A subset of 'enrichment.tbl' holding only the significant results.
#' @export
selectSignificantEnrichments <- function(enrichment.tbl, significance.level = 0.05) {
    enrichment.tbl[which(enrichment.tbl$P.Value.adjusted <= significance.level), 
        ]
}

#' Adds a column to 'enrichment.tbl' holding human readable descriptions for
#' the annotations found o be enriched. This function has been written
#' primarily for use with InterPro annotations.
#'
#' @param enrichment.tbl a data.frame result of function 'enrichmentTests' or
#' 'selectSignificantEnrichments'
#' @param annotation.db a list in which the names are function annotation
#' identifiers, e.g. IPR666666
#' @param ano.desc.field the name of the list entry to select in order to
#' obtain the description from members of 'annotation.db'
#'
#' @return A copy of 'enrichment.tbl' with an additional column holding the
#' decriptions of the respective annotations.
#' @export
cbindAnnotationDescription <- function(enrichment.tbl, annotation.db = ipr.db, ano.desc.field = "NAME") {
    enrichment.tbl$description <- as.character(lapply(enrichment.tbl[, 1], function(x) {
        y <- annotation.db[[x]]
        if (!is.null(y) && ano.desc.field %in% names(y)) y[[ano.desc.field]] else NA
    }))
    enrichment.tbl
}

#' Computes a Gene Ontology (GO) enrichment test using the functions from #'
#' package GOstats and GSEABase.
#'
#' @param gsc An instance of gene set collection [package GEABase]
#' @param gene.ids A character vector of gene IDs among which to look for
#'        enriched GO terms
#' @param univ.gene.ids A character vector of 'background' gene IDs to be used
#'        as reference GO annotations when looking for overrepresentations in
#'        'gene.ids'
#' @param ontologies A character vector of the Gene Ontology categories to
#'        compute enrichment tests for, default is 'BP', 'CC', and 'MF'
#' @param pvalue.cutoff The significance level, default to 0.01
#' @param cond If set to TRUE already found to be enriched GO annotations wont
#'        be counted when testing enrichment for parental terms, default is
#'        FALSE
#' @param test.dir Set to 'over' to infer overrepresentation, set to 'under'
#'        if the opposite is wanted.
#' @param p.adjust.method The method used to adjust p-values for multiple
#'        hypothesis testing, default is 'fdr'. See function p.adjust for
#'        details.
#'
#' @return A named list with enriched GO terms in the argument GO categories.
#' @export
goEnrichTest <- function(gsc, gene.ids, univ.gene.ids, ontologies = c("BP", "CC", 
    "MF"), pvalue.cutoff = 0.01, cond = FALSE, test.dir = "over", p.adjust.method = "fdr") {
    setNames(mclapply(ontologies, function(go.ont) {
        tryCatch({
            ghgr <- GOstats::hyperGTest(Category::GSEAGOHyperGParams(name = "", geneSetCollection = gsc, 
                geneIds = gene.ids, universeGeneIds = univ.gene.ids, ontology = go.ont, 
                pvalueCutoff = pvalue.cutoff, conditional = cond, testDirection = test.dir))
            procGOHyperGResult(ghgr, pvalue.cutoff)
        }, error = function(e) {
            warning("Error in goEnrichTest(...) when computing GO-Ontology ", go.ont, 
                " ", e)
            NA
        })
    }), ontologies)
}


#' Converts an instance of GOHyperGResult into a table containing only the
#' significantly enriched Gene Ontology terms.
#'
#' @param ghgr The instance of GOHyperGResult to process
#' @param pvalue.cutoff The significance level to be applied
#'
#' @return A data.frame with three columns: 1. 'GO.term', 2. 'p.value', and 3.
#'         'GO.category' - or NULL if no significantly enriched GO terms were
#'         found.
#' @export
procGOHyperGResult <- function(ghgr, pvalue.cutoff = 0.01) {
    go.ont <- strsplit(ghgr@testName, " ")[[2]]
    x.pv <- pvalues(ghgr)
    x.res <- x.pv[which(x.pv <= pvalue.cutoff)]
    if (length(x.res) > 0) {
        x.go.names <- as.character(lapply(names(x.res), function(go.acc) {
            go.term <- GO.db::GOTERM[[go.acc]]
            if (!is.null(go.term)) go.term@Term
        }))
        data.frame(GO.term = names(x.res), GO.category = go.ont, GO.name = x.go.names, 
            p.value = as.numeric(x.res), stringsAsFactors = FALSE)
    } else NULL
}

#' Method extracts all GO ontology specific enrichment result and combines them
#' into a single data frame, re-adjusting the Pvalues for multiple hypothesis
#' testing.
#'
#' @param go.enrich.lst The result of calling method goEnrichTest(...).
#' @param p.adjust.method The method to use when adjusting p-values for
#'        multiple hypothesis testing. Default is 'fdr' [see p.adjust for
#'        details]
#'
#' @return A data.frame with the merged enriched GO terms and adjusted
#' p-values.
#' @export
joinGOEnrichResults <- function(go.enrich.lst, p.adjust.method = "fdr") {
    res <- do.call("rbind", go.enrich.lst)
    if (!is.null(res)) 
        res$p.value.adjusted <- p.adjust(res$p.value, method = p.adjust.method)
    res
}

#' Identifies those composite annotations that are less frequent than the most
#' frequent ones found for the argument genes and based upon argument gene
#' (function) annotations. For more details about composite annotations see
#' function 'compositeAnnotations(...)'
#'
#' @param gene.accs The identifiers or accessions of the genes to compute the
#' entropy for
#' @param gene.annos The data.frame holding the annotations for the genes in
#' 'gene.accs'. Default is all available InterPro annotations expected to be
#' found in 'all.ipr'
#' @param gene.col the column of 'gene.annos' in which to lookup the gene
#' identifiers or gene accessions. Default is 1
#' @param anno.col the column of 'gene.annos' in which to lookup the function
#' annotation for the genes in 'gene.accs'. Default is 2
#'
#' @export
#' @return An instance of base::table
nonConservedCompositeAnnotations <- function(gene.accs, gene.annos = all.ipr, gene.col = 1, 
    anno.col = 2) {
    gene.comp.annos <- compositeAnnotations(gene.accs, gene.annos, gene.col, anno.col)
    gene.annos.freqs <- table(gene.comp.annos)
    max.freq <- max(gene.annos.freqs)
    non.cons.funks <- names(gene.annos.freqs[which(gene.annos.freqs < max.freq)])
    gene.comp.annos[which(gene.comp.annos %in% non.cons.funks)]
}

#' Test for function nonConservedCompositeAnnotations.
#'
#' @export
#' @return TRUE if and only if all tests pass
testNonConservedCompositeAnnotations <- function() {
    gene.accs <- c("A", "B", "C", "D")
    gene.annos <- read.table(stringsAsFactors = FALSE, text = "A funk_1
A funk_2
A funk_3
B funk_1
B funk_2
B funk_3
C funk_1
C funk_3
D funk_4
D funk_5
D funk_6
E funk_1
E funk_3
F funk_1
F funk_3")
    t1 <- identical(setNames(c("funk_1,funk_3", "funk_4,funk_5,funk_6"), c("C", "D")), 
        nonConservedCompositeAnnotations(gene.accs, gene.annos))
    t2 <- length(nonConservedCompositeAnnotations(gene.accs[1:2], gene.annos)) == 
        0
    t3 <- length(nonConservedCompositeAnnotations(c("X", "Y"), gene.annos)) == 0
    t4 <- identical(setNames(c("funk_1,funk_2,funk_3", "funk_1,funk_2,funk_3", "funk_4,funk_5,funk_6", 
        ""), c("A", "B", "D", "G")), nonConservedCompositeAnnotations(c(gene.accs, 
        "E", "F", "G"), gene.annos))
    all(t1, t2, t3, t4)
}

#' Identifies those gene annotations that are non conserved within the argument
#' gene-group. Non conserved is inferred as those annotations that do not
#' appear or miss from the most frequent composite annotation.
#'
#' @param gene.accs Character-Vector of gene-IDs.
#' @param gene.annos A data.frame of gene annotations
#' @param gene.col The column of \code{gene.annos} in which to lookup the gene
#' IDs. Default is 1
#' @param anno.col The column of \code{gene.annos} in which to lookup the gene
#' annotations. Default is 2
#' @param rm.genes.without.annos If set to TRUE genes without any annotation
#' are ignored, otherwise this will cause all conserved annotations to appear
#' in the result. Default is TRUE
#'
#' @export
#' @return A data.frame with two columns: 'Gene' holds the gene identifiers of
#' genes that have non conserved annotations, and column
#' 'Non.Conserved.Annotation' holds the respective annotations.
nonConservedAnnotations <- function(gene.accs, gene.annos = all.ipr, gene.col = 1, 
    anno.col = 2, rm.genes.without.annos = TRUE) {
    if (rm.genes.without.annos) {
        gene.accs <- gene.accs[which(gene.accs %in% gene.annos[, gene.col])]
    }
    if (!is.null(gene.accs) && !is.na(gene.accs) && length(gene.accs) > 0) {
        gene.comp.annos <- compositeAnnotations(gene.accs, gene.annos, gene.col, 
            anno.col)
        gene.annos.freqs <- table(gene.comp.annos)
        max.freq <- max(gene.annos.freqs)
        cons.annos <- sort(unique(unlist(strsplit(names(gene.annos.freqs[which(gene.annos.freqs == 
            max.freq)]), ","))))
        non.cons.comp.annos <- names(gene.annos.freqs[which(gene.annos.freqs < max.freq)])
        genes.nca <- names(gene.comp.annos[which(gene.comp.annos %in% non.cons.comp.annos)])
        do.call("rbind", lapply(genes.nca, function(x) {
            x.annos <- gene.annos[which(gene.annos[, gene.col] == x), anno.col]
            nc.annos <- setdiff(union(x.annos, cons.annos), intersect(x.annos, cons.annos))
            if (!is.null(nc.annos) && !is.na(nc.annos) && length(nc.annos) > 0) {
                data.frame(Gene = x, Non.Conserved.Annotation = nc.annos, stringsAsFactors = FALSE)
            } else NULL
        }))
    } else NULL
}

#' Test function \code{nonConservedAnnotations}
#'
#' @export
#' @return TRUE if and only if all tests pass
testNonConservedAnnotations <- function() {
    gene.accs <- 1:3
    gene.annos <- read.table(stringsAsFactors = FALSE, text = "1 A
1 B
2 A
2 B
3 A
3 B
3 C
4 A
4 D")
    t.1 <- identical(data.frame(Gene = "3", Non.Conserved.Annotation = "C", stringsAsFactors = FALSE), 
        nonConservedAnnotations(gene.accs, gene.annos))
    res <- nonConservedAnnotations(c(gene.accs, "4"), gene.annos)
    t.2 <- all(c("B", "C", "D") %in% res$Non.Conserved.Annotation) && all(c("3", 
        "4") %in% res$Gene)
    t.3 <- is.null(nonConservedAnnotations(1:2, gene.annos))
    t.4 <- is.null(nonConservedAnnotations(c(1:2, "666"), gene.annos, rm.genes.without.annos = TRUE))
    t.5 <- is.null(nonConservedAnnotations(c("666", "999"), gene.annos, rm.genes.without.annos = TRUE))
    all(c(t.1, t.2, t.3, t.4, t.5))
}

#' For each species identifies those annotations that are non conserved within
#' those gene groups (e.g. families) that are indicated by \code{group.ids}.
#' Non conservation is inferred by
#' \code{nonConservedCompositeAnnotations(...)}.
#'
#' @param genes.df A data.frame with three columns, two of which need to be
#' named 'Gene' and 'Species'. The first column needs to hold the gene-group
#' IDs.
#' @param group.col Identifies the column of \code{genes.df} in which to lookup
#' gene-group-identifiers. Default is \code{'Tandem.Cluster'}.
#' @param group.ids Vector holding indicating in which gene-groups to look for
#' non conserved composite annotations. Default is \code{unique(genes.df[,
#' group.col])}.
#' @param species The species for which to identify non conserved composite
#' annotations. Default is \code{unique(genes.df$Species)}.
#' @param gene.annos The data.frame holding the annotations for the genes in
#' 'gene.accs'. Default is all available InterPro annotations expected to be
#' found in 'all.ipr'
#' @param annos.gene.col The column of \code{gene.annos} in which to lookup the gene
#' identifiers or gene accessions. Default is 1
#' @param anno.col the column of \code{gene.annos} in which to lookup the function
#' annotation for the genes in 'gene.accs'. Default is 2
#'
#' @export
#' @return A data.frame with four columns: 1. The names of the gene-groups, 2.
#' Species, 3.  genes, and 4. the non conserved annotations (now
#' non-composite!).
collectNonConservedAnnotations <- function(genes.df, group.col = "Tandem.Cluster", 
    group.ids = unique(genes.df[, group.col]), gene.annos = all.ipr, annos.gene.col = 1, 
    annos.anno.col = 2) {
    do.call("rbind", mclapply(group.ids, function(x) {
        clstr <- genes.df[which(genes.df[, group.col] == x), ]
        gene.accs <- clstr$Gene
        z <- nonConservedAnnotations(gene.accs, gene.annos = gene.annos, gene.col = annos.gene.col, 
            anno.col = annos.anno.col)
        if (!is.null(z) && !is.na(z) && nrow(z) > 0) {
            z$Gene.Group <- x
            z$Species <- as.character(unlist(lapply(z$Gene, function(i) clstr[which(clstr$Gene == 
                i), "Species"][[1]])))
            # Re-order columns for convenience:
            z[, c("Gene.Group", "Species", "Gene", "Non.Conserved.Annotation")]
        } else NULL
    }))
}

#' Tests the NULL Hypothesis for each gene (function) annotation in
#' \code{annos.2.test} that its annotation-frequency in the \code{case.genes}
#' is explicable with the background distribution observed in
#' \code{universe.annos}. The default alternative hypothesis is that the actual
#' case frequency is 'greater'.
#'
#' @param case.genes A character vector of gene identifiers defining the gene
#' group of interest for which to infer enriched annotations.
#' @param universe.annos A data.frame with at least two columns: 1. of gene
#' identifiers and another of (function) annotations. Default is this package's
#' data \code{all.ipr}.
#' @param univ.gene.col The column of \code{universe.annos} in which to lookup
#' gene identifier. Default is 1.
#' @param univ.anno.col The column of \code{universe.annos} in which to lookup
#' annotations. Default is 2.
#' @param annos.2.test The annotations for which to test enrichment. Default is
#' all annotations found for the \code{case.genes}.
#' @param alt.hypothesis The alternative hypotheses to test. Default is that
#' the case frequencies are significantly greater than the one observed in the
#' 'universe'.
#' @param p.adjust.method The method used to adjust P-Values for multiple
#' hypothesis testing. Default is 'BY'. See \code{?p.adjust} for more details.
#'
#' @export
#' @return A numeric vector of P-Values adjusted to multiple hypothesis
#' testing. Names are the \code{annos.2.test}.
enrichedAnnotations <- function(case.genes, universe.annos = all.ipr, univ.gene.col = 1, 
    univ.anno.col = 2, annos.2.test = unique(universe.annos[which(universe.annos[, 
        univ.gene.col] %in% case.genes), univ.anno.col]), alt.hypothesis = "greater", 
    p.adjust.method = "BY") {
    univ.genes <- universe.annos[, univ.gene.col]
    case.genes.with.annos <- intersect(univ.genes, case.genes)
    if (length(case.genes.with.annos) > 0) {
        setNames(p.adjust(as.numeric(unlist(mclapply(annos.2.test, function(i.anno) {
            genes.with.anno <- universe.annos[which(universe.annos[, univ.anno.col] == 
                i.anno), univ.gene.col]
            genes.without.anno <- setdiff(univ.genes, genes.with.anno)
            non.case.genes <- setdiff(univ.genes, case.genes.with.annos)
            cont.tbl <- generateContingencyTable(genes.with.anno, genes.without.anno, 
                case.genes.with.annos, non.case.genes, i.anno, "Gene.Group")
            fisher.test(cont.tbl, alternative = alt.hypothesis)$p.value
        }))), method = p.adjust.method), annos.2.test)
    } else {
        warning("Function 'enrichedAnnotations': None of the genes in argument", 
            " 'case.genes' had annotations in 'universe.annos'. Returning NA.")
        NA
    }
}

#' Selects those results obtained from calling \code{enrichedAnnotations} that
#' are significant and adds informative information to the annotations
#' extracted from the \code{anno.db}.
#'
#' @param enriched.annotations A numeric vector with values the P-Values and
#' names the annotations. See \code{enrichedAnnotations(...)} for details.
#' @param p.val.cutoff The significance cuto-off to be applied. Default is 0.01
#' @param anno.db A list with names being the annotations as in the names of
#' \code{enriched.annotations} and values are lists with a slot holding human
#' readable information about the annotation.
#' @param anno.name The slot name of the lists held in \code{anno.db} to
#' extract for human readable information about the annotations.
#'
#' @export
#' @return A data.frame with three columns: 1. Annotation, 2. P.Value.Adjusted,
#' and 3. Annotation.Name
selectSignificantlyEnrichedAnnotations <- function(enriched.annotations, p.val.cutoff = 0.01, 
    anno.db = ipr.db, anno.name = "NAME") {
    sign.enr.annos <- enriched.annotations[which(enriched.annotations < p.val.cutoff)]
    data.frame(Annotation = names(sign.enr.annos), P.Value.Adjusted = sign.enr.annos, 
        Annotation.Name = as.character(unlist(lapply(names(sign.enr.annos), function(x) {
            y <- anno.db[[x]]
            if (!is.null(y) && !is.na(y) && length(y) > 0) {
                y[[anno.name]]
            } else NA
        }))), stringsAsFactors = FALSE)
}

#' Infer enriched annotations specific to the argument type of gene-groups, add
#' human readable information to it and return it as a data.frame.
#'
#' @param gene.group.type The type of the gene.groups for which to identify
#' specific enrichments.
#' @param gene.group.type.enriched.annos A list of enriched annotations per
#' gene group. Names must include \code{gene.group.type}. Default is the
#' enriched annotations found in this packages data.
#' @param anno.db A list with names being the annotations as in the names of
#' \code{enriched.annotations} and values are lists with a slot holding human
#' readable information about the annotation.
#' @param anno.name The slot name of the lists held in \code{anno.db} to
#' extract for human readable information about the annotations.
#'
#' @export
#' @return A data.frame with two columns: 1. Annotation and 2. Annotation.Name
selectGeneGroupSpecificEnrichedAnnotations <- function(gene.group.type, gene.group.type.enriched.annos = list(positively.selected = fams.pos.sel.non.conserved.iprs.enrich.df$Annotation, 
    expanded.contracted = exp.fams.df.non.conserved.iprs.enrich.df$Annotation, tandems = bra8.tandem.clusters.non.conserved.iprs.enrich.df$Annotation, 
    conserved = fams.conserved.non.conserved.iprs.enrich.df$Annotation), anno.db = ipr.db, 
    anno.name = "NAME") {
    uniq.enrich.annos <- setdiff(gene.group.type.enriched.annos[[gene.group.type]], 
        unique(unlist(gene.group.type.enriched.annos[setdiff(names(gene.group.type.enriched.annos), 
            gene.group.type)])))
    data.frame(Annotation = uniq.enrich.annos, Annotation.Name = as.character(unlist(lapply(uniq.enrich.annos, 
        function(x) {
            y <- anno.db[[x]]
            if (!is.null(y) && !is.na(y) && length(y) > 0) {
                y[[anno.name]]
            } else NA
        }))), stringsAsFactors = FALSE)
}
