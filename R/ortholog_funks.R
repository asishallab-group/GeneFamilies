#' Identify homologous gene pairs in a table of significantly similar genes.
#'
#' @param seq.sim.tbl a data frame of two character columns and a third score
#' column. It is highly important that this table represents only gene pairs of
#' TWO species, not more or less.
#'
#' @return A subset of the argument data frame, holding only the two character
#' columns, containing the candidate homologs. 
#' @export
homologs <- function(seq.sim.tbl) {
    # Exclude self matches:
    x <- seq.sim.tbl[which(seq.sim.tbl[, "V1"] != seq.sim.tbl[, "V2"]), ]
    # Order by score, such that best scoring hits are listed first for each query:
    x <- x[with(x, order(V1, -V3)), c("V1", "V2")]
    # Retain only the best hit for each query:
    x <- x[which(!duplicated(x$V1)), ]
    # Identify reciprocal best hits:
    y <- data.frame(V1 = x$V2, V2 = x$V1, stringsAsFactors = FALSE)
    dbl.x <- rbind(x, y)
    y <- dbl.x[duplicated(dbl.x[, 1:2]), ]
    rownames(y) <- 1:nrow(y)
    y
}

#' Test of function 'homologs(...)'
#'
#' @return TRUE if and only if the test succeeded.
#' @export
testHomologs <- function() {
    df <- data.frame(V1 = c(1, 1, 1:5), V2 = c(1, 4, 4, 3, 2, 1, 8), V11 = c(0, 
        0, 1, rep(0, 4)), stringsAsFactors = FALSE)
    expected.res <- data.frame(V1 = 4:1, V2 = 1:4, V11 = rep(0, 2), stringsAsFactors = FALSE)
    res <- homologs(df)
    # Ignore increased row-names in 'identical' comparison:
    all(as.logical(lapply(1:ncol(res), function(i) identical(as.numeric(res[, i]), 
        as.numeric(expected.res[, i])))))
}

#' Identifies orthologous gene pairs as genes being significantly similar in
#' sequence and belonging to different species.
#'
#' @param homologs.tbl a data.frame of two character columns - possibly the
#' result obtained from invoking homologous(...)
#' @param gene.regex.1 a string encoding a PERL regular expression to match
#' gene identifiers to species A
#' @param gene.regex.2 a string encoding a PERL regular expression to match
#' gene identifiers to species B
#'
#' @return A subset of argument 'homologs.tbl' containg the orthologous gene
#' pairs. Note that if pair (A,B) is contained the reciprocal pair (B,A) is,
#' too. This should ease lookup.
#' @export
orthologs <- function(homologs.tbl, gene.regex.1, gene.regex.2) {
    paralogs.1 <- grepl(gene.regex.1, homologs.tbl$V1, perl = TRUE) & grepl(gene.regex.1, 
        homologs.tbl$V2, perl = TRUE)
    paralogs.2 <- grepl(gene.regex.2, homologs.tbl$V1, perl = TRUE) & grepl(gene.regex.2, 
        homologs.tbl$V2, perl = TRUE)
    homologs.tbl[which(!(paralogs.1 | paralogs.2)), ]
}

#' Parses 'psl' default BLAT output formatted files to infer orthologs between
#' species. Uses PERL regular expressions to identify which gene identifiers
#' belongs to which species. Hence by using two such regular expressions
#' orthologs between two species can be extracted, even though the BLAT output
#' file contains similarity pairs from more than two species. Note that
#' orthology is inferred based on reciprocal best search results.
#'
#' @param blat.df a data frame holding the output of BLAT in psl format. Se
#' function readBlatPSLoutput(...) for more details
#' @param query.col the number of the column in which to lookup the query
#' sequence's name
#' @param target.col the number of the column in which to lookup the target
#' sequence's name
#' @param matching.bases.col the number of the column in which to lookup the
#' number of nucleotide bases matching between query and target
#' @param qSize.col the number of the column in which to lookup the query
#' sequence's size
#' @param tSize.col the number of the column in which to lookup the target
#' sequence's size
#' @param spec.regexs a character vector containg two regular expressions
#' matching gene identifiers from each of the two species orthologies are to be
#' found for
#'
#' @return A data frame of two columns holding the gene identifiers of the
#' orthologous gene pairs. Note, that for reasons of easing lookup inverse
#' pairs are also contained.
#' @export
orthologsFromBLATpslTable <- function(blat.df, query.col = 10, target.col = 14, 
    matching.bases.col = 1, qSize.col = 11, tSize.col = 15, spec.regexs = c("CAHR\\d+(\\.\\d)?", 
        "AT[0-9MC]G\\d+(\\.\\d)?")) {
    # Retain only those gene pairs, where each gene matches one of the species'
    # regular expressions:
    b.df <- blat.df[which((grepl(spec.regexs[[1]], blat.df[, query.col], perl = TRUE) | 
        grepl(spec.regexs[[1]], blat.df[, target.col], perl = TRUE)) & (grepl(spec.regexs[[2]], 
        blat.df[, query.col], perl = TRUE) | grepl(spec.regexs[[2]], blat.df[, 
        target.col], perl = TRUE))), ]
    b.df <- b.df[which((b.df[, matching.bases.col] >= ((b.df[, qSize.col] + b.df[, 
        tSize.col])/4))), c(query.col, target.col, matching.bases.col)]
    colnames(b.df) <- paste("V", 1:ncol(b.df), sep = "")
    # Identify homologous candidates and retain orthologs as reciprocal best hits
    # of genes belonging to different species:
    orthologs(homologs(b.df), spec.regexs[[1]], spec.regexs[[2]])
}

#' Wraps function base::read.table with the proper arguments needed to
#' efficiently and correctly read in a BLAT output table in PSL format.
#'
#' @param blat.psl.path valid file path to the respective BLAT psl output file
#'
#' @return The read in table as data frame
#' @export
readBlatPSLoutput <- function(blat.psl.path) {
    read.table(blat.psl.path, sep = "\t", skip = 5, stringsAsFactors = FALSE, colClasses = c(rep("numeric", 
        8), rep("character", 2), rep("numeric", 3), "character", rep("numeric", 
        4), rep("character", 3)), comment.char = "", quote = "", na.strings = "")
}

#' Identifies those pairs of orthologous genes where one member appears in all
#' pairwise species comparisons.
#'
#' @param pairwise.orths a list of two column data frames holding orthologous
#' gene pairs identified in pairwise species comparisons. The names of this
#' list reflect the species the anchoring species was compared to iteratively
#' @param intersect.col the name or index of the column of each of the data
#' frames in 'pairwise.orths' to use as a source of gene identifiers
#' @param ortholog.col the name or index of the column to extract ortholog
#' genes from when they match anchor species orthologs conserved over all
#' species comparisons
#' @param anchor.species the name of the species which was used in all pairwise
#' comparisons.
#'
#' @return a two column data frame of orthologous gene groups where an ortholog
#' was found in each species.
#' @export
conservedOrthologs <- function(pairwise.orths, intersect.col = "V1", ortholog.col = "V2", 
    anchor.species = "chi") {
    orths.genes <- lapply(pairwise.orths, function(x) unique(as.character(unlist(x[, 
        intersect.col]))))
    conserved.genes <- base::Reduce(intersect, orths.genes)
    conserved.orths <- data.frame(V1 = conserved.genes, stringsAsFactors = FALSE)
    colnames(conserved.orths) <- anchor.species
    cbind(conserved.orths, do.call("cbind", lapply(pairwise.orths, function(x) x[which(x[, 
        intersect.col] %in% conserved.genes), ortholog.col])), stringsAsFactors = FALSE)
}

#' Assuming 'gene.ids' holds gene identifiers that lack splice variant suffixes
#' like '.1' but refer to genes contained in this packages coding sequences
#' (data) this function aims to identify the matching splice variants.
#'
#' @param 'gene.ids' character vector of gene identifiers to be matched to
#' coding sequence identifiers held in this packages data.
#' @param 'genes' character vector of gene identifiers with splice variant
#' suffixes to be matched against the argument gene identifiers.
#'
#' @return Character vector of the original gene identifiers or the matching
#' splice variants if the original argument identifiers were not found in
#' 'genes'. NA wherever the argument in 'gene.ids' was not found in 'genes' and
#' could also not be matched to an entry of 'genes'.
#' @export
matchingSpliceVariant <- function(gene.ids, genes = as.character(unlist(spec.gene.ids))) {
    i <- which(!gene.ids %in% genes)
    gene.ids[i] <- as.character(unlist(mclapply(gene.ids[i], function(x) {
        y <- which(grepl(paste("\\b", x, "\\.\\d+\\b", sep = ""), genes, perl = TRUE))
        if (length(y) != 1) {
            warning("Could not identify unique matching splicing variant for gene identifier '", 
                x, "'!")
            NA
        } else genes[[y]]
    })))
    gene.ids
}

#' Extracts clusters of highly conserved orthologous genes from pairwise
#' reciprocal best hits obtained from sequence similarity searches. The
#' searches are explained in detail in this package's Vignette
#' \code{GeneFamilies}, section \code{Ortholog Identification}. Only ortholog
#' clusters that have the exact same size as the numer of investigated species
#' are returned.
#'
#' @param pair.best.hits.tbl An instance of \code{base::data.frame} holding the
#' final result of the concatonated filtered pairwise sequence similarity
#' searches. See Vignette \code{GeneFamilies} for details on how to generate
#' this input table.
#' @param spec.gene.ids.arg An instance of \code{base::list} with names being
#' the investigated species and values all protein accessions belonging to the
#' respective species. Default is \code{spec.gene.ids} (see file \code{zzz.R}
#' for more details).
#' @param anchor.gene.col The column-name or column-index of
#' \code{pair.best.hits.tbl} in which to find the anchor genes, i.e. those
#' genes of the species all pairwise sequence similarity searches were
#' conducted with. Default is \code{1}.
#' @param candidate.gene.col The column-name or column-index of
#' \code{pair.best.hits.tbl} in which to find the candidate genes, that were
#' the best hits in the pairwise sequence similarity searches. Default is
#' \code{2}.
#'
#' @return An instance of \code{base::data.frame} with one ortholog cluster per
#' row and one column for each of the investigated species. An additional
#' column holds the ortholog cluster name, composed as
#' \code{'ortholog_cluster_<number>'}.
#' @export
orthologsFromPairwiseReciprocalBestHits <- function(pair.best.hits.tbl, spec.gene.ids.arg = spec.gene.ids, 
    anchor.gene.col = 1, candidate.gene.col = 2) {
    anchor.genes <- unique(pair.best.hits.tbl[, anchor.gene.col])
    orths.df <- as.data.frame(matrix(vector(), ncol = length(spec.gene.ids), nrow = 0, 
        dimnames = list(c(), sort(names(spec.gene.ids)))), stringsAsFactors = FALSE)
    for (x in anchor.genes) {
        orth.cands <- pair.best.hits.tbl[which(pair.best.hits.tbl[, anchor.gene.col] == 
            x), candidate.gene.col]
        if (length(orth.cands) == length(spec.gene.ids.arg) - 1) {
            orth.genes <- c(orth.cands, x)
            orth.genes.named <- setNames(orth.genes, unlist(lapply(orth.genes, 
                function(y) speciesForGeneId(y, spec.gene.ids.arg))))
            o.df.row <- nrow(orths.df) + 1
            orths.df[o.df.row, ] <- orth.genes.named[sort(names(spec.gene.ids))]
        }
    }
    o.cl.names <- paste("ortholog_cluster_", 1:nrow(orths.df), sep = "")
    orths.df$Cluster <- o.cl.names
    orths.df[, c("Cluster", sort(names(spec.gene.ids)))]
}

#' Informs whether two genes are located on the same chromosome and are within
#' a maximum distance. Distance is measured NOT in nucleotides, but in number
#' of separating genes.
#'
#' @param gene.a The gene identifier for the first gene in question
#' @param gene.b The gene identifier for the second gene in question
#' @param neighborhood.tbl an instance of \code{base::data.frame} with at least
#' two columns. One holding the gene identifiers, and the other the chromosome
#' name. This table MUST be sorted by chromose and gene position, i.e. genes
#' separated by three rows and being located on the same chromosome are
#' expected to be separated by three genes on the respective chromosome.
#' @param gene.col Column name or number of \code{neighborhood.tbl} indicating
#' where to find the gene identifier. Default is
#' \code{getOption('GeneFamilies.neighborhood.gene.col', 2)}.
#' @param chromosome.col Column name or number of \code{neighborhood.tbl}
#' indicating where to find the chromosome names. Default is
#' \code{getOption('GeneFamilies.neighborhood.chromosome.col', 1)}.
#' @param max.dist An integer defining the maximum number of genes allowed to
#' separate two genes to be still considered 'neighbors'. Default is
#' \code{getOption('GeneFamilies.neighborhood.max.dist', 9)}
#'
#' @return A boolean or NA, if the two genes are not found in
#' \code{neighborhood.tbl}.
#' @export
areNeighbors <- function(gene.a, gene.b, neighborhood.tbl, gene.col = getOption("GeneFamilies.neighborhood.gene.col", 
    2), chromosome.col = getOption("GeneFamilies.neighborhood.chromosome.col", 
    1), max.dist = getOption("GeneFamilies.neighborhood.max.dist", 9)) {
    genes <- c(gene.a, gene.b)
    i <- sort(which(neighborhood.tbl[, gene.col] %in% genes))
    if (length(i) != 2) {
        warning("GeneFamilies::areNeighbors - Not all supplied argument genes are in argument 'neighborhood.tbl'")
        return(NA)
    }
    n.chr <- unique(neighborhood.tbl[i, chromosome.col])
    if (length(n.chr) == 1) {
        sum(i) - i[1] <= max.dist
    } else FALSE
}

#' Identifies tandem duplicated genes within a gene family. Conditions are that
#' the tandems belonging to the same species are within a gene neighborhod of
#' (default) nine genes on the same chromosome. See
#' \code{GeneFamilies::areNeighbors} for more details.
#' _Note_ that this approach does NOT exclude orthologous genes. If you want to
#' distinguish those do a subsequent \code{setdiff} with the result of this
#' function.
#'
#' @param gene.family A character vector holding the gene identifiers of the
#' genes comprising a single gene family.
#' @param spec.neighborhoods A list with names representing the species and
#' values being the tables indicating gene neighborhood. See the Vignette
#' \code{GeneFamilies} on how to generate these neighborhood-tables. Also see
#' function \code{GeneFamilies::areNeighbors} for more details on these tables'
#' format.
#' @param spec.gene.ids.arg An instance of \code{base::list} with names being
#' the investigated species and values all protein accessions belonging to the
#' respective species. Default is \code{spec.gene.ids} (see file \code{zzz.R}
#' for more details).
#'
#' @return A character vector holding all Tandem Duplicated Genes found within
#' the argument gene family.
#' @export
tandemsFromFamiliesAndNeighborhood <- function(gene.family, spec.neighborhoods, 
    spec.gene.ids.arg = spec.gene.ids) {
    if (length(gene.family) >= 2) {
        genes.df <- data.frame(gene = gene.family, species = unlist(lapply(gene.family, 
            speciesForGeneId, spec.gene.ids.arg)), stringsAsFactors = FALSE)
        tands <- c()
        for (spec in names(spec.gene.ids.arg)) {
            spec.genes <- genes.df[which(genes.df$species == spec), "gene"]
            if (length(spec.genes) > 1) {
                all.pairs <- combn(spec.genes, 2)
                spec.genes.neighbrs <- c()
                for (i in 1:ncol(all.pairs)) {
                  x <- all.pairs[, i]
                  if (areNeighbors(x[1], x[2], spec.neighborhoods[[spec]])) {
                    spec.genes.neighbrs <- union(spec.genes.neighbrs, x)
                  }
                }
                if (length(spec.genes.neighbrs) > 0) 
                  tands <- union(tands, spec.genes.neighbrs)
            }
        }
        return(tands)
    } else c()
}

#' Function testing \code{GeneFamilies::tandemsFromFamiliesAndNeighborhood}.
#'
#' @param \code{NULL}
#'
#' @return TRUE if and only if all tests pass
#' @export
testTandemsFromFamiliesAndNeighborhood <- function() {
    spec.nghbr.tbl <- list(ath = data.frame(V1 = c(rep("chr1", 3), rep("chr2", 
        3)), V2 = LETTERS[1:6], stringsAsFactors = FALSE), chi = data.frame(V1 = paste("chr", 
        c(1:3, 3), sep = ""), V2 = letters[1:4], stringsAsFactors = FALSE))
    gene.fam <- c(LETTERS[1:6], letters[1:4])
    s.g.i <- list(ath = LETTERS[1:6], chi = letters[1:4])
    t.1 <- identical(c(LETTERS[1:6], letters[3:4]), tandemsFromFamiliesAndNeighborhood(gene.fam, 
        spec.nghbr.tbl, s.g.i))
    options(GeneFamilies.neighborhood.max.dist = 0)
    t.2 <- identical(NULL, tandemsFromFamiliesAndNeighborhood(gene.fam, spec.nghbr.tbl, 
        s.g.i))
    options(GeneFamilies.neighborhood.max.dist = NULL)
    all(c(t.1, t.2))
}
