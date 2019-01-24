#' Computes the Euclidean norm (length) of vector \code{x}.
#'
#' @param x numeric vector
#'
#' @return A scalar numeric value, the length of the vector \code{x}.
#' @export
euclNorm <- function(x) {
    if (hasInvalidVecSpaceComponents(x)) 
        return(NaN)
    sqrt(sum(x^2))
}

#' Calculates Reads Per Kilobase Million (RPKM)
#' RPKM = (Rg * 10^9) / (T * Lg)
#' where
#' Rg: number of reads mapped to a particular transcript g = count
#' T = total number of transcripts sampled in run
#' Lg: length of transcript g (bases) 
#' Note, that because we use gene length in bases and not in kilobases the
#' scaling factor is 10^9 and not 10^6.
#'
#' @param raw.count.tbl A data.frame with at least three columns: One must
#' contain gene IDs (accessions), another the genes' lengths, and the other
#' columns are expected to hold the read counts.
#' @param gene.col The column name or number in which to find the gene
#' identifier (accessions). Default is 'gene'.
#' @param gene.length.col The column name or number in which to find the gene
#' lengths. Default is 'gene.length'.
#' @param replicate.cols The names or numbers of the columns in which to find
#' the raw read counts. Default is \code{setdiff(colnames(raw.count.tbl),
#' c(gene.col, gene.length.col))}
#'
#' @export
#' @return An extended version of \code{raw.count.tbl} with added columns of
#' the RPKM values, one column for each column of raw reads in argument
#' \code{replicate.cols}.
transcriptsPerKilobaseMillion <- function(raw.count.tbl, gene.col = "gene", 
    gene.length.col = "gene.length", replicate.cols = setdiff(colnames(raw.count.tbl), 
        c(gene.col, gene.length.col))) {
    for (repl.col in replicate.cols) {
        total.counts <- sum(raw.count.tbl[, repl.col], na.rm = TRUE)
        raw.count.tbl[, paste(repl.col, "_RPKM", sep = "")] <- 1e+09 * 
            (raw.count.tbl[, repl.col])/(total.counts * raw.count.tbl[, 
            gene.length.col])
    }
    raw.count.tbl
}

#' Computes the distance (\code{base::dist}) between the expression profiles of
#' the genes in argument \code{gene.accessions}. 
#'
#' @param gene.accessions A character vector of gene identifiers some of which
#' are expected to be present in \code{expression.profiles}
#' @param expression.profiles A data.frame with gene expression profiles.
#' Default is this package's \code{rna.seq.exp.profils} data.
#' @param expr.prof.gene.col A number or string identifying the column of
#' \code{expression.profiles} in which to lookup the gene accessions. Default
#' is \code{'gene'}.
#' @param expr.prof.gene.variant.col A number or string identifying the column
#' of \code{expression.profiles} in which to lookup the gene accessions
#' including their expression variant identifier (e.g. 'AT1G12345.3'). Default
#' is \code{gene.exp.var}.
#' @param tissues A character or integer indexing those columns of
#' \code{expression.profiles} in which the actual relative expression per
#' tissue is stored. In short these columns are the expression profiles.
#' Default is \code{setdiff(colnames(expression.profiles),
#' c(expr.prof.gene.col, expr.prof.gene.variant.col))}.
#' @param dist.method A string to pass as argument \code{method} to
#' \code{base::dist}. Default is \code{'euclidean'}.
#' @param per.tissue Logical indicating wether to compute distances in the
#' euclidean space created by the dimensions of all tissues, or if set to
#' \code{TRUE} computes the one-dimensional distances per tissue. Default is
#' \code{FALSE}.
#'
#' @export
#' @return An instance of base::dist holding pairwise distances between the
#' argument genes expression profiles. If not at least two expression profiles
#' are found for the argument \code{gene.accessions} \code{NA} is returned. In
#' case the argument \code{per.tissue} is set to \code{TRUE} a list of
#' one-dimensional distances is returned; names of the list are the respective
#' tissues.
expressionProfilesDists <- function(gene.accessions, expression.profiles = rna.seq.exp.profils, 
    expr.prof.gene.col = "gene", expr.prof.gene.variant.col = "gene.exp.var", 
    tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, 
        expr.prof.gene.variant.col)), dist.method = "euclidean", per.tissue = FALSE) {
    inds <- which(expression.profiles[, expr.prof.gene.col] %in% gene.accessions | 
        expression.profiles[, expr.prof.gene.variant.col] %in% gene.accessions)
    if (length(inds) > 1) {
        exp.profs <- expression.profiles[inds, ]
        rownames(exp.profs) <- exp.profs[, expr.prof.gene.col]
        if (per.tissue) {
            setNames(lapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)
    } else NA
}

#' Extracts the inter- and intra-species distances from argument object
#' \code{expr.prof.dists} which is expected to be an instance of class
#' \code{base::dist} (See \code{expressionProfilesDists} for more details).
#'
#' @param expr.prof.dists an instance of \code{base::dist} a single result from
#' calling \code{expressionProfilesDists}.
#' @param gene.ids.to.specs.regex A named character where the names are species
#' and the values are regular expressions capable of uniquely identifying which
#' species a gene accession belongs to. Default is \code{list(
#' 'ath'='^AT[0-9CM]G\\d+$', 'chi'='^CARHR\\d+$' )} as within this package's
#' application we only have expression data for genes of Arabidopsis and C.
#' hirsuta and in this function's case we are dealing with genes without
#' expression variant identifiers.
#'
#' @export
#' @return A list with either data.frame or numeric values. Data.frames have
#' three columns, two indicating gene pairs and the third the extracted
#' distances, while numeric values point to maximum distances. Result list
#' names are the species found as names in \code{gene.ids.to.specs.regex},
#' \code{inter.species}, \code{max.inter.species}, and for each species
#' \code{spec_i} one \code{max.spec_i}.
interVsIntraSpeciesExpressionProfileDists <- function(expr.prof.dists, 
    gene.ids.to.specs.regex = list(ath = "^AT[0-9CM]G\\d+$", chi = "^CARHR\\d+$")) {
    #' Validate input:
    if (is.null(expr.prof.dists) || all(is.na(expr.prof.dists))) 
        return(NA)
    #' Initialize:
    genes <- attr(expr.prof.dists, "Label")
    genes.specs.df <- data.frame(gene = genes, species = as.character(lapply(genes, 
        speciesForGeneId_Regex, spec.regexs = gene.ids.to.specs.regex)), 
        stringsAsFactors = FALSE)
    species <- intersect(names(gene.ids.to.specs.regex), genes.specs.df$species)
    #' Validate species:
    if (length(species) == 0) 
        stop("Found ZERO species matching gene accessions in argument 'expr.prof.dists' and being present in argument 'gene.ids.to.specs.regex'.")
    #' Extract the requested intra-species distances:
    intra.spec.dists <- setNames(lapply(species, function(spec) {
        spec.genes <- genes.specs.df[which(genes.specs.df$species == spec), 
            "gene"]
        if (length(spec.genes) > 0) {
            x <- expand.grid(spec.genes, spec.genes)
            #' Do not extract meaningless distances for pairs of identical
            #' genes, e.g. (gene_a, gene_a):
            x <- x[which(x[, 1] != x[, 2]), ]
            if (nrow(x) > 0) {
                x$dist <- as.numeric(apply(x, 1, function(g.p) {
                  as.matrix(expr.prof.dists)[[g.p[[1]], g.p[[2]]]]
                }))
                x
            } else NA
        } else NA
    }), species)
    #' Extract the requested inter-species distances:
    inter.spec.dists <- if (length(species) > 1) {
        dist.df <- do.call("rbind", apply(combn(species, 2), 2, function(spec.pair) {
            x.spec.genes <- genes.specs.df[which(genes.specs.df$species == 
                spec.pair[[1]]), "gene"]
            y.spec.genes <- genes.specs.df[which(genes.specs.df$species == 
                spec.pair[[2]]), "gene"]
            expand.grid(x.spec.genes, y.spec.genes)
        }))
        dist.df$dist <- as.numeric(apply(dist.df, 1, function(g.p) {
            as.matrix(expr.prof.dists)[[g.p[[1]], g.p[[2]]]]
        }))
        dist.df
    } else NA
    #' Return result as list:
    res.lst <- intra.spec.dists
    res.lst$inter.species <- inter.spec.dists
    for (stat in c("max", "median", "mean")) {
        res.lst <- append(res.lst, setNames(lapply(intra.spec.dists, function(x) {
            if (!is.null(x) && !is.na(x) && nrow(x) > 0) 
                eval(call(stat, x$dist, na.rm = TRUE)) else NA
        }), paste(stat, names(intra.spec.dists), sep = ".")))
        inter.spec.stat.nm <- paste(stat, ".inter.species", sep = "")
        res.lst[[inter.spec.stat.nm]] <- if (length(species) > 1) 
            eval(call(stat, inter.spec.dists$dist, na.rm = TRUE)) else NA
        abs.stat.nm <- paste(stat, ".abs", sep = "")
        res.lst[[abs.stat.nm]] <- eval(call(stat, expr.prof.dists, na.rm = TRUE))
    }
    res.lst
}

#' Generates a data.frame with one row for each family with expression profile
#' distance measurements and columns the maximum distances found for intra-,
#' inter-species, and all pairwise distances.
#'
#' @param fams.inter.intra.exp.prof.dists The result of invoking function
#' \code{interVsIntraSpeciesExpressionProfileDists}.
#' @param family.names The names of those families to generate rows in the
#' result data frame. Set to all families that actually have non NA results in
#' \code{fams.inter.intra.exp.prof.dists}. Default is
#' \code{names(fams.inter.intra.exp.prof.dists)}.
#'
#' @export
#' @return An instance of \code{data.frame} with one row for each family in
#' argument \code{family.names} and columns holding the maximum values
#' explained above.
maxValuesFromInterVsIntraSpeciesExpressionProfileDists <- function(fams.inter.intra.exp.prof.dists, 
    family.names = names(fams.inter.intra.exp.prof.dists)) {
    fams.inter.intra.exp.prof.dists.df <- data.frame(family = family.names, 
        max.ath = as.numeric(NA), max.chi = as.numeric(NA), max.inter.species = as.numeric(NA), 
        max.abs = as.numeric(NA), median.ath = as.numeric(NA), median.chi = as.numeric(NA), 
        median.inter.species = as.numeric(NA), median.abs = as.numeric(NA), 
        mean.ath = as.numeric(NA), mean.chi = as.numeric(NA), mean.inter.species = as.numeric(NA), 
        mean.abs = as.numeric(NA), stringsAsFactors = FALSE)
    for (fam.name in family.names) {
        fam.dists <- fams.inter.intra.exp.prof.dists[[fam.name]]
        if (!is.null(fam.dists) && !all(is.na(fam.dists)) && length(fam.dists) > 
            0) {
            for (stat in c("max", "median", "mean")) {
                stat.vals <- names(fam.dists)[grepl(stat, names(fam.dists), 
                  fixed = TRUE)]
                fams.inter.intra.exp.prof.dists.df[which(fams.inter.intra.exp.prof.dists.df$family == 
                  fam.name), stat.vals] <- fam.dists[stat.vals]
            }
        }
    }
    fams.inter.intra.exp.prof.dists.df
}

#' Extract different \code{stats} per tissue specific expression profile distances.
#'
#' @param expr.prof.dists.per.tissue An instance of list with names identifying
#' tissues and values being instances of \code{base::dist}. See function
#' \code{expressionProfilesDists(..., per.tissue = TRUE)}.
#' @param gene.classes An instance of \code{list} where the names are the
#' different gene classes and the values character vectors intersecting with
#' the argument genes and thus segregating them into the respective
#' gene-classes.
#' @param fam.name The name of the group of genes the distances were
#' computed for. Will be the row-name of the returned numeric matrix. Default
#' is \code{NA}.
#' @param stats A character vector holding the names of the functions to apply
#' on the the tissue specific distances. Default is \code{stats = c('max',
#' 'median', 'mean', 'maxMinusMin')}.
#'
#' @export
#' @return An instance of \code{matrix} with numeric values holding the stats
#' inferred for each group of tissue specific distances as found in argument
#' \code{expr.prof.dists.per.tissue}. Returns \code{NULL} if
#' \code{all(is.na(expr.prof.dists.per.tissue))}.
expressionProfileDistStatsPerTissues <- function(expr.prof.dists.per.tissue, 
    gene.classes, fam.name = NA, stats = c("max", "median", "mean", "maxMinusMin")) {
    res.df <- Reduce(rbind, lapply(names(expr.prof.dists.per.tissue), function(tissue) {
        epd.tissue <- expr.prof.dists.per.tissue[[tissue]]
        stats.df <- classSpecificExpressionProfileDists(epd.tissue, gene.classes, 
            fam.name, stats)$stats
        stats.df$Tissue <- tissue
        stats.df
    }))
    res.df
}


#' Computes the absolute distance between maximum and minimum, optionally
#' removing \code{NA}.
#'
#' @param x a numeric
#' @param na.rm If set to \code{TRUE}, the default, \code{NA} values in
#' \code{x} are removed before computing \code{max} and \code{min}.
#'
#' @export
#' @return A single numeric value 
maxMinusMin <- function(x, na.rm = TRUE) {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
}


#' Computes the within \code{stats} of expression profile distances for classes
#' of gene groups separating the argument genes into sub-groups.
#'
#' @param e.p.d An instance of \code{stats::dist} holding the pairwise
#' expression profile distances for the argument genes.
#' @param gene.classes An instance of \code{list} where the names are the
#' different gene classes and the values character vectors intersecting with
#' the argument genes and thus segregating them into the respective
#' gene-classes.
#' @param fam.name A string to be assigned as \code{Family} column in the
#' \code{base::data.frame} in slot \code{stats} of the returned object. Default
#' is \code{NA}.
#' @param stats A character vector holding the names of the functions to apply
#' on the the gene.class specific distances. Default is \code{stats = c('max',
#' 'median', 'mean', 'maxMinusMin')}.
#'
#' @export
#' @return An instance of \code{base::list} with one \code{base::dist} per
#' Gene.Class in \code{names(gene.classes)} and an additional
#' \code{base::data.frame} named \code{'stats'}. The latter holds the values
#' returned by the functions in \code{stats} applied to each separate class of
#' argument genes.
classSpecificExpressionProfileDists <- function(e.p.d, gene.classes, fam.name = NA, 
    stats = c("max", "median", "mean", "maxMinusMin")) {
    epd.m <- as.matrix(e.p.d)
    genes <- colnames(epd.m)
    genes.classified <- lapply(gene.classes, function(x) intersect(genes, 
        x))
    res <- lapply(genes.classified, function(x) as.dist(epd.m[rownames(epd.m) %in% 
        x, colnames(epd.m) %in% x]))
    res$stats <- Reduce(rbind, lapply(stats, function(i.stat) {
        Reduce(rbind, lapply(names(gene.classes), function(k.class) {
            genes.dist <- res[[k.class]]
            dist.stat <- if (length(genes.dist) > 0 && !all(is.na(genes.dist) | 
                is.nan(genes.dist) | is.infinite(genes.dist))) {
                eval(call(i.stat, unlist(res[[k.class]]), na.rm = TRUE))
            } else NA
            data.frame(Family = fam.name, Statistic = i.stat, Gene.Class = k.class, 
                Value = dist.stat, stringsAsFactors = FALSE)
        }))
    }))
    res
}


#' Identifies the value that corresponds to the upper whisker end in an
#' boxplot. \code{NA} are ignored
#'
#' @param x a numerical vector
#'
#' @export
#' @return A single number
upperWhisker <- function(x) {
    q.s <- quantile(x, seq(0.25, 1, 0.25), na.rm = TRUE)
    min(max(x, na.rm = TRUE), q.s[[3]] + 1.5 * (q.s[[3]] - q.s[[1]]), na.rm = TRUE)
}


#' Adds an alpha value to colors.
#'
#' @param A character vector of colors.
#' @param numeric alpha value between 0 and 1. Default is \code{.1}
#'
#' @export
#' @return A vector of colors adjusted with the argument alpha value.
addAlpha <- function(col, alpha = 0.1) {
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], 
        alpha = alpha))
}

#' Plots horizontal boxplots of per tissue distributions of expression profile
#' distances. Draws a line from each upper Whisker to the next in order to
#' highlight the hourglass shape. Tissues must be given in the order according
#' to a plant's life.
#'
#' @param plot.df A data.frame with the following columns: \code{'Tissue',
#' 'Gene.Class', 'Value'}.
#' @param tissues A character vector of tissues held in \code{plot.df$Tissue}.
#' The order found in this argument will be the order on the Y-Axis of the
#' horizontal plot.
#' @param gene.classes A character vector of gene.classes to group paired
#' boxplots by. Must be of length 2.
#' @param pdf.path The valid file path to the PDF in which the horizontal
#' boxplot is to be stored.
#' @param plot.title A string defining the headline printedabove the plot. See
#' argument \code{main} in function \code{base::boxplot} for more details.
#' Default is \code{''}.
#'
#' @export
#' @return The output of \code{dev.off()}
plotTissueSpecificExpressionProfileDistanceDistribs <- function(plot.df, 
    tissues, gene.classes, pdf.path, plot.title = "") {
    plot.df$interaction <- interaction(plot.df$Tissue, plot.df$Gene.Class)
    plot.df$interaction <- factor(plot.df$interaction, levels = unlist(lapply(tissues, 
        paste, gene.classes, sep = ".")))
    line.df <- rbind(data.frame(Tissue = tissues, Value = as.numeric(lapply(tissues, 
        function(x) upperWhisker(plot.df[which(plot.df$Tissue == x & plot.df$Gene.Class == 
            gene.classes[[1]]), "Value"]))), Gene.Class = gene.classes[[1]], 
        stringsAsFactors = FALSE), data.frame(Tissue = tissues, Value = as.numeric(lapply(tissues, 
        function(x) upperWhisker(plot.df[which(plot.df$Tissue == x & plot.df$Gene.Class == 
            gene.classes[[2]]), "Value"]))), Gene.Class = gene.classes[[2]], 
        stringsAsFactors = FALSE))
    br.cols <- brewer.pal(length(levels(plot.df$interaction)), "Paired")
    pdf(pdf.path)
    par(mar = c(5, 15, 4, 2) + 0.1)
    ticks.at <- sort(c(seq(1, length(tissues) * 2 + 3, by = 3), seq(2, 
        length(tissues) * 2 + 4, by = 3)))
    boxplot(Value ~ interaction, data = plot.df, outline = FALSE, horizontal = TRUE, 
        las = 1, border = br.cols, col = addAlpha(br.cols), lwd = 2.5, 
        xlab = "median of Euclidean distances between expression profiles", 
        main = plot.title, names = NA, at = ticks.at)
    axis(2, at = ticks.at, labels = unlist(lapply(tissues, function(x) paste(x, 
        " (", gene.classes, ")", sep = ""))), las = 1)
    line.cols <- brewer.pal(3, "Greys")[2:3]
    lines(x = line.df[which(line.df$Gene.Class == gene.classes[[1]]), ]$Value, 
        y = seq(1, length(tissues) * 2 + 3, by = 3), col = line.cols[[1]], 
        lwd = 3)
    lines(x = line.df[which(line.df$Gene.Class == gene.classes[[2]]), ]$Value, 
        y = seq(2, length(tissues) * 2 + 4, by = 3), col = line.cols[[2]], 
        lwd = 3)
    dev.off()
}


#' Test the hypothesis of the developmental hourglass by carrying out various
#' t-tests. In these, compare the distributions of expression profile distances
#' within a given tissue between gene classes. This should show which gene
#' class contributes to the hourglass model. Also compare for a given gene
#' class early and late developmental tissues with intermediate samples. This
#' should show wether early and late stages are more diverse than intermediate
#' ones, and thus give the hourglass shape.
#'
#' @param x an instance of \code{base::data.frame} with columns: 'Family',
#' 'Statistic', 'Gene.Class', 'Value', 'Tissue'. For example this package's
#' data \code{tandems.exp.prof.dists.tissue.orth.dist}.
#' @param tissues a character vector identifying the tissues to be
#' investigated. Default is \code{c('seedling', 'cotyledon', 'developing leaf',
#' 'flower stage 9', 'flower stage 16')}.
#' @param gene.classes a character vector identifying the gene classes to be
#' investigated. Default is \code{c('Orthologs', 'Non-Orthologs')}.
#' @param early.late.indx A numeric vector identifying which entries of
#' \code{tissues} are considered early and late developmental stages. Default
#' is \code{c(2,4)}.
#' @param intermediate.indx A single integer identifying which entry in
#' \code{tissues} is considered to be the intermediate developmental stage.
#' Default is \code{3}.
#' @param statistic a string identifying which statistic of \code{x} to be
#' investigated. Default is \code{'median'}.
#'
#' @export
#' @return A \code{base::data.frame} with columns: 'Tissue',
#' 'Alternative.Hypothesis', 'p.value', 'effect.size', 'p.value.adjusted'
testHourglassModel <- function(x, tissues = c("seedling", "cotyledon", 
    "developing leaf", "flower stage 9", "flower stage 16"), gene.classes = c("Orthologs", 
    "Non-Orthologs"), early.late.indx = c(2, 4), intermediate.indx = 3, 
    statistic = "median") {
    res.df <- rbind(Reduce(rbind, lapply(tissues, function(tissue) {
        y.ort <- x[which(x$Tissue == tissue & x$Statistic == statistic & 
            x$Gene.Class == "Orthologs" & !is.na(x$Value)), "Value"]
        y.no <- x[which(x$Tissue == tissue & x$Statistic == statistic & 
            x$Gene.Class == "Non-Orthologs" & !is.na(x$Value)), "Value"]
        h.test <- t.test(y.ort, y.no, alternative = "less")
        data.frame(Tissue = tissue, Alternative.Hypothesis = "Non-Orthologs have higher mean expression profile distances than Orthologs", 
            p.value = h.test$p.value, effect.size = h.test$statistic, stringsAsFactors = FALSE)
    })), Reduce(rbind, unlist(lapply(gene.classes, function(gene.class) {
        lapply(tissues[early.late.indx], function(tissue) {
            y.dev.leaf <- x[which(x$Tissue == tissues[intermediate.indx] & 
                x$Statistic == statistic & x$Gene.Class == gene.class & 
                !is.na(x$Value)), "Value"]
            y.foreground <- x[which(x$Tissue == tissue & x$Statistic == 
                statistic & x$Gene.Class == gene.class & !is.na(x$Value)), 
                "Value"]
            h.test <- t.test(y.dev.leaf, y.foreground, alternative = "less")
            hypo <- paste(gene.class, "have higher median expression profile distances in", 
                tissue, "than in developing leaf")
            data.frame(Tissue = tissue, Alternative.Hypothesis = hypo, 
                p.value = h.test$p.value, effect.size = h.test$statistic, 
                stringsAsFactors = FALSE)
        })
    }), recursive = FALSE)))
    res.df$p.value.adjusted <- p.adjust(res.df$p.value, method = "BY")
    res.df
}

#' Generates a data.frame of the Genes' names, copy numbers, and expression
#' values (RPKM). Furthermore, computes a lm model for the observed data. If
#' expression levels are extracted for more than a single tissue the
#' \code{mean} value is used. Note, that only genes are considered that have
#' some expression in any tissue, never expressed genes are discarded from the
#' analysis.
#'
#' @param rpkm.df An instance of data frame holding the genes expression counts
#' (RPKM). Default is \code{rpkm.rna.seq.counts}.
#' @param r.gene.col A string or integer identifying the column of
#' \code{rpkm.df} in which to find the gene accessions (IDs). Default is
#' \code{'id'}.
#' @param r.rpkm.col A string or integer identifying the column of
#' \code{rpkm.df} in which to find the genes' expression counts (RPKM). Default
#' is \code{'expression'}.
#' @param r.tissue.col A string or integer identifying the column of
#' \code{rpkm.df} in which to find the tissue in which the expression was
#' measured. Default is \code{'tissue'}.
#' @param copy.no.df An instance of \code{base:data.frame} in which to find the
#' genes' copy numbers. Default is \code{'gene.copy.number.df'}.
#' @param c.gene.col A string or integer identifying the column of
#' \code{copy.no.df} in which the gene accessions (IDs) are stored. Default is
#' \code{'Gene.no.expr.var'}.
#' @param c.copy.no.col A string or integer identifying the column of
#' \code{copy.no.df} in which the genes' copy numbers are held. Default is
#' \code{'copy.no'}.
#' @param tissues A character of minimum length one in which the tissues to
#' extract expression levels for are stored. Default is \code{unique(rpkm.df[,
#' r.tissue.col])}
#' @param stat.funk The statistical measure applied on the expression values
#' found for a respective gene within the argument \code{tissues}. Default is
#' \code{base::mean}.
#' @param lapply.funk One of \code{base::lapply} or \code{parallel::mclapply}
#' to indicate which function shall be used to iterate over the respective
#' genes. Default is \code{lapply}.
#'
#' @export
#' @return A list with entries: \code{data} the above data.frame, \code{lm} the
#' above mentioned computed generalized linear model, and the \code{R^2} (R
#' squared) value.
correlationRpkmCopyNo <- function(rpkm.df = rpkm.rna.seq.counts, r.gene.col = "id", 
    r.rpkm.col = "expression", r.tissue.col = "tissue", copy.no.df = gene.copy.number.df, 
    c.gene.col = "Gene.no.expr.var", c.copy.no.col = "copy.no", tissues = unique(rpkm.df[, 
        r.tissue.col]), stat.funk = mean, lapply.funk = lapply) {
    r.df <- rpkm.df[which(rpkm.df[, r.tissue.col] %in% tissues), ]
    genes <- intersect(r.df[which(r.df[, r.rpkm.col] > 0), r.gene.col], 
        copy.no.df[, c.gene.col])
    corr.df <- cbind(data.frame(Gene = genes, row.names = genes, stringsAsFactors = FALSE), 
        Reduce(rbind, lapply.funk(genes, function(gene) {
            g.rpkm <- stat.funk(r.df[which(r.df[, r.gene.col] == gene), 
                r.rpkm.col])
            g.copy.no <- copy.no.df[which(copy.no.df[, c.gene.col] == gene), 
                c.copy.no.col]
            data.frame(row.names = gene, RPKM = g.rpkm, Copies = g.copy.no, 
                stringsAsFactors = FALSE)
        })))
    lm.model <- lm(formula = RPKM ~ Copies, data = corr.df)
    list(data = corr.df, lm = lm.model, r.squared = summary(lm.model)$r.squared)
}

#' Compute cosinus of the angle between argument vector \code{x} and the vector
#' space diagonal.
#'
#' @param x numeric representing a vector in n-dimensional space
#' @param d.v numeric the diagonal vector in n-dimensional space. Default is
#' \code{rep(1, length(x))}.
#'
#' @export
#' @return Returns the cosinus of the angle between argument vector \code{x}
#' and the vector space diagonal. If any value in argument \code{x} is infinite
#' or not a number or if \code{x} is \code{NULL} \code{NaN} is returned.
cosDiag <- function(x, d.v = rep(1, length(x))) {
    cosAngleVec(x, d.v)
}

#' Computes the cosine of the angle between the two n-dimensional argument
#' vectors \code{a} and \code{b}.
#'
#' @param a numerical representing first vector
#' @param b numerical representing second vector
#'
#' @export
#' @return cosine of angle between the two argument vectors or NaN if any entry
#' in the vectors is not a number.
cosAngleVec <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'cosAngleVec(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    sum(a * b)/(euclNorm(a) * euclNorm(b))
}

#' Generic function to check wether any vector component is non-numerical.
#'
#' @param x generic object
#'
#' @export
#' @return TRUE if and only if any component in any vector in \code{x} is
#' non-numerical.
hasInvalidVecSpaceComponents <- function(x) {
    UseMethod("hasInvalidVecSpaceComponents", x)
}

#' Generic function to check wether any vector component is non-numerical.
#'
#' @param x numeric vector
#'
#' @export
#' @return TRUE if and only if any component in any vector in \code{x} is
#' non-numerical.
hasInvalidVecSpaceComponents.numeric <- function(x) {
    any(is.na(x) | is.nan(x) | is.null(x) | is.infinite(x))
}

#' Generic function to check wether any vector component is non-numerical.
#'
#' @param x A numeric matrix in which each row is interpreted as a single
#' vector.
#'
#' @export
#' @return TRUE if and only if any component in any vector in \code{x} is
#' non-numerical.
hasInvalidVecSpaceComponents.matrix <- function(x) {
    any(as.logical(apply(x, 1, hasInvalidVecSpaceComponents)))
}

#' Generic function to check wether any vector component is non-numerical.
#'
#' @param x A data.frame in which each row is interpreted as a single vector.
#'
#' @export
#' @return TRUE if and only if any component in any vector in \code{x} is
#' non-numerical.
hasInvalidVecSpaceComponents.data.frame <- function(x) {
    any(as.logical(apply(x, 1, hasInvalidVecSpaceComponents)))
}

#' Computes the scalar projection of vector \code{a} onto vector \code{b}, i.e.
#' the length of the segment of \code{b} 'overshadowed' by \code{a}. See
#' en.wikipedia.org/wiki/Vector_projection
#'
#' @param a numerical representing first vector
#' @param b numerical representing second vector
#'
#' @export
#' @return A scalar value, the length of the segment of \code{b} overshadowed
#' by \code{a}.
scalarProjection <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'scalarProjection(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    
    euclNorm(a) * cosAngleVec(a, b)
}


#' Computes the vector projection of vector \code{a} onto vector \code{b}, i.e.
#' the vector-segment of \code{b} 'overshadowed' by \code{a}. See
#' en.wikipedia.org/wiki/Vector_projection
#'
#' @param a numerical representing first vector
#' @param b numerical representing second vector
#'
#' @export
#' @return The vector resulting from projecting vector \code{a} upon vector
#' \code{b}.
vectorProjection <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'vectorProjection(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    scalarProjection(a, b) * b/euclNorm(b)
}

#' Computes the statistic \code{stat} (default is \code{base::mean}) vector of
#' the cloud of vectors. Also infers the deviance in terms of
#' \code{deviance.funk} (default is \code{base::sd}) of each axis of the vector
#' space.
#'
#' @param vecs.df A matrix of data.frame of row vectors of the same vector
#' space
#' @param stat A function used to infer the vector representing the argument
#' cloud. The function is applied on each column seperately in order to infer
#' the resulting vector's components. Default is
#' \code{getOption('GeneFamilies.vector.cloud.stat', mean)}. 
#' @param deviance.funk Function applied on each column of the cloud
#' \code{vecs.df} in order to infer the cloud's spread around it's representing
#' vector. Default is \code{getOption('GeneFamilies.vector.cloud.deviance',
#' sd)}.
#'
#' @export
#' @return A list with three named entries: \code{stat.vec} the vector
#' representing the cloud, \code{deviance.vec} the vector representing the
#' spread around \code{stat.vec}, and the vector orthogonal on the vector space
#' diagonal pointing to the \code{stat.vec}.
statVectorCloud <- function(vecs.df, stat = getOption("GeneFamilies.vector.cloud.stat", 
    mean), deviance.funk = getOption("GeneFamilies.vector.cloud.deviance", 
    sd)) {
    if ((class(vecs.df) == "data.frame" || class(vecs.df) == "matrix") && 
        nrow(vecs.df) > 0) {
        i <- which(!as.logical(apply(vecs.df, 1, hasInvalidVecSpaceComponents)))
        if (length(i) == 0) 
            return(NA)
        v.df <- vecs.df[i, , drop = FALSE]
        stat.vec <- as.numeric(apply(v.df, 2, stat))
        dev.vec <- as.numeric(apply(v.df, 2, deviance.funk))
        orth.on.diag.2.stat.vec <- stat.vec - vectorProjection(stat.vec, 
            rep(1, length(stat.vec)))
        list(stat.vec = stat.vec, deviance.vec = dev.vec, orth.on.diag.2.stat.vec = orth.on.diag.2.stat.vec)
    }
}


#' Helper function that converts argument \code{x} to the number zero if and
#' only if \code{is.na(x)}.
#'
#' @param x any scalar
#'
#' @export
#' @return \code{0.0} if \code{is.na(x)} or \code{x} otherwise
naAsZero <- function(x) {
    if (is.na(x)) 
        0 else x
}

#' Measures the distance between two vector clouds (see
#' \code{GeneFamilies::statVectorCloud}) within the same vector space. The
#' vector between the two \code{stat.vec}s representing \code{cl.a} and
#' \code{cl.b} is computed. The deviance spaces of the two respective clouds
#' are projected onto this distance vector, and the remaining 'un-shadowed'
#' length returned. If this length is positive the two clouds do not overlap
#' considering their representatives and the irespective deviance spreads
#' around them.
#'
#' @param cl.a Result of invoking \code{GeneFamilies::statVectorCloud} on a set
#' of vectors.
#' @param cl.b Result of invoking \code{GeneFamilies::statVectorCloud} on a set
#' of vectors.
#'
#' @export
#' @return Numeric - The length of the difference vector being outside the
#' respective clouds' deviances.
distVectorClouds <- function(cl.a, cl.b) {
    if (is.na(cl.a) || is.na(cl.b)) 
        return(NA)
    cl.diff <- cl.a$stat.vec - cl.b$stat.vec
    dev.proj.a <- naAsZero(scalarProjection(cl.a$deviance.vec, cl.diff))
    dev.proj.b <- naAsZero(scalarProjection(cl.b$deviance.vec, cl.diff))
    euclNorm(cl.diff) - (abs(dev.proj.a) + abs(dev.proj.b))
}

#' Drops a perpendicular from vector \code{x} to the vector space's diagonal.
#' The perpendicular points towards \code{x}.
#'
#' @param x numeric representing a vector in n-dimensional space
#' @param d.v numeric the diagonal vector in n-dimensional space. Default is
#' \code{rep(1, length(x))}.
#'
#' @export
#' @return A vector in the same space of \code{x}, the perpendicular.
perpVecToDiagonal <- function(x, d.v = rep(1, length(x))) {
    if (hasInvalidVecSpaceComponents(x)) 
        return(NaN)
    x - vectorProjection(x, d.v)
}

#' Normalizes the vector \code{x} to length 1.0.
#'
#' @param numeric vector 
#'
#' @export
#' @return \code{x / euclNorm(x)} or \code{NaN} if
#' \code{hasInvalidVecSpaceComponents(x)} is \code{TRUE}.
normalizeVector <- function(x) {
    if (hasInvalidVecSpaceComponents(x)) 
        return(NaN)
    x/euclNorm(x)
}

#' Returns the vector resulting from subtracting \code{e.v.b} from
#' \code{e.v.a}. Both argument vectors are normalized before subtraction.
#'
#' @param e.v.a numeric
#' @param e.v.b numeric
#'
#' @export
#' @return \code{NA} if any of the arguments contains components that are non
#' legal coordinates (see \code{hasInvalidVecSpaceComponents}). Otherwise
#' \code{normalizeVector(e.v.a) - normalizeVector(e.v.b)} is returned.
distNormExprVectrs <- function(e.v.a, e.v.b) {
    if (hasInvalidVecSpaceComponents(rbind(e.v.a, e.v.b))) 
        return(NA)
    normalizeVector(e.v.a) - normalizeVector(e.v.b)
}

#' Converts an angle's size given in \code{radians} to \code{degrees}.
#'
#' @param rad a numeric
#'
#' @export
#' @return The converted \code{rad} value in degrees.
rad2deg <- function(rad) {
    (rad * 180)/(pi)
}

#' Based on an expression vector space infers the evolution of gene function
#' after duplication by comparing base genes and groups of evolved genes.
#'
#' @param genes character vector of gene identifiers
#' @param family.name The name of the gene group comprising argument \code{genes} 
#' @param classifier.lst A named list of at least two entries. Entries should
#' be lists of gene IDs belonging to their respective (named) classes.
#' @param base.class A name of argument \code{classifier.lst} to be used as
#' base gene class. Default is \code{'ortholog'}.
#' @param expr.vecs A set of vectors in the respective expression vector space.
#' Each row holds the expression vector for its respective gene. Default is
#' \code{rpkm.expr.profiles.df}.
#' @param gene.col The name or index of a column of argument \code{expr.vecs}
#' indicating in which column to lookup gene IDs.
#' @param vec.space.axes A character or integer vector indicating the columns
#' of argument \code{expr.vecs} which correspond to the respective vector space
#' axes. Default is \code{c('cotyledon', 'seedling', 'developing leaf', 'flower
#' stage 9', 'flower stage 16')}.
#'
#' @export
#' @return A data.frame with the following columns: Family (set to argument
#' \code{family.name}, mean.base.tiss.vers, and for each class 'EvCl' of
#' evolved genes EvCl.tiss.vers, EvCl.tiss.change, and EvCl.dist.vec.clouds
exprVecSpaceEvolutionAfterDupl <- function(genes, family.name, classifier.lst, 
    base.class = "ortholog", expr.vecs = rpkm.expr.profiles.df, gene.col = "gene", 
    vec.space.axes = c("cotyledon", "seedling", "developing leaf", "flower stage 9", 
        "flower stage 16")) {
    res <- NULL
    if (length(genes) > 0) {
        genes.classed <- lapply(classifier.lst, function(x) intersect(intersect(genes, 
            x), expr.vecs[, gene.col]))
        g.c.i <- sapply(genes.classed, function(x) length(x) > 0)
        evolved.classes <- setdiff(names(genes.classed), base.class)
        # If base class and at least one other class has genes with expression
        # values:
        if (length(genes.classed[[base.class]]) > 0 && any(g.c.i[evolved.classes])) {
            base.class.cloud <- statVectorCloud(expr.vecs[which(expr.vecs[, 
                gene.col] %in% genes.classed[[base.class]]), vec.space.axes])
            res <- cbind(data.frame(Family = family.name, mean.base.tiss.vers = (1 - 
                cosDiag(base.class.cloud$stat.vec)/sqrt(2)), stringsAsFactors = FALSE), 
                Reduce(cbind, lapply(names(g.c.i[evolved.classes]), function(evol.class) {
                  # Check whether evol.class has expression values for its genes:
                  if (g.c.i[[evol.class]]) {
                    e.g <- genes.classed[[evol.class]]
                    e.g.cloud <- statVectorCloud(expr.vecs[which(expr.vecs[, 
                      gene.col] %in% e.g), vec.space.axes])
                    evol.df <- data.frame(V1 = (1 - cosDiag(e.g.cloud$stat.vec)/sqrt(2)), 
                      V2 = rad2deg(acos(cosAngleVec(base.class.cloud$orth.on.diag.2.stat.vec, 
                        e.g.cloud$orth.on.diag.2.stat.vec))), V3 = distVectorClouds(base.class.cloud, 
                        e.g.cloud), stringsAsFactors = FALSE)
                    colnames(evol.df) <- paste(evol.class, c("tiss.vers", 
                      "tiss.change", "dist.vec.clouds"), sep = ".")
                    evol.df
                  } else {
                    evol.df <- data.frame(V1 = NA, V2 = NA, V3 = NA, stringsAsFactors = FALSE)
                    colnames(evol.df) <- paste(evol.class, c("tiss.vers", 
                      "tiss.change", "dist.vec.clouds"), sep = ".")
                    evol.df
                  }
                })))
        }
    }
    res
}


#' Testing function \code{exprVecSpaceEvolutionAfterDupl}.
#'
#' @export
#' @return \code{TRUE} if and only if all tests pass
testExprVecSpaceEvolutionAfterDupl <- function() {
    genes <- LETTERS[1:6]
    fam <- "Test_Family"
    g.class <- list(ortholog = LETTERS[1:2], tandem = LETTERS[3:4], pos.selected = LETTERS[5:6])
    tissues <- paste("tiss", LETTERS[1:3], sep = ".")
    vec.space <- read.table(stringsAsFactors = FALSE, header = TRUE, text = "gene tiss.A tiss.B tiss.C
A 10 2 2
B 9 2 3
C 2 10 2
D 3 9 2
E 2 2 10
F 3 2 9")
    y.1 <- exprVecSpaceEvolutionAfterDupl(genes, fam, g.class, expr.vecs = vec.space, 
        vec.space.axes = tissues)
    test.1.a <- identical(colnames(y.1), c("Family", "mean.base.tiss.vers", 
        "tandem.tiss.vers", "tandem.tiss.change", "tandem.dist.vec.clouds", 
        "pos.selected.tiss.vers", "pos.selected.tiss.change", "pos.selected.dist.vec.clouds"))
    test.1.b <- identical(y.1$Family, fam)
    res.df.1 <- read.table(stringsAsFactors = FALSE, header = TRUE, text = "
Family mean.base.tiss.vers tandem.tiss.vers tandem.tiss.change tandem.dist.vec.clouds pos.selected.tiss.vers pos.selected.tiss.change pos.selected.dist.vec.clouds
Test_Family 0.4298759 0.4298759 120 9.720577 0.4298759 113.164 9.899495")
    test.1.c <- identical(round(y.1[, 2:8], digits = 2), round(res.df.1[, 
        2:8], digits = 2))
    test.1 <- all(c(test.1.a, test.1.b, test.1.c))
    all(test.1)
}
