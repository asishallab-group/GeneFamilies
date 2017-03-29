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
transcriptsPerKilobaseMillion <- function(raw.count.tbl, gene.col = "gene", gene.length.col = "gene.length", 
    replicate.cols = setdiff(colnames(raw.count.tbl), c(gene.col, gene.length.col))) {
    for (repl.col in replicate.cols) {
        total.counts <- sum(raw.count.tbl[, repl.col], na.rm = TRUE)
        raw.count.tbl[, paste(repl.col, "_RPKM", sep = "")] <- 1e+09 * (raw.count.tbl[, 
            repl.col])/(total.counts * raw.count.tbl[, gene.length.col])
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
    expr.prof.gene.col = "gene", expr.prof.gene.variant.col = "gene.exp.var", tissues = setdiff(colnames(expression.profiles), 
        c(expr.prof.gene.col, expr.prof.gene.variant.col)), dist.method = "euclidean", 
    per.tissue = FALSE) {
    inds <- which(expression.profiles[, expr.prof.gene.col] %in% gene.accessions | 
        expression.profiles[, expr.prof.gene.variant.col] %in% gene.accessions)
    if (length(inds) > 1) {
        exp.profs <- expression.profiles[inds, ]
        rownames(exp.profs) <- exp.profs[, expr.prof.gene.col]
        if (per.tissue) {
            setNames(lapply(tissues, function(tissue) dist(setNames(exp.profs[, tissue], 
                rownames(exp.profs)), method = dist.method)), tissues)
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
interVsIntraSpeciesExpressionProfileDists <- function(expr.prof.dists, gene.ids.to.specs.regex = list(ath = "^AT[0-9CM]G\\d+$", 
    chi = "^CARHR\\d+$")) {
    #' Validate input:
    if (is.null(expr.prof.dists) || all(is.na(expr.prof.dists))) 
        return(NA)
    #' Initialize:
    genes <- attr(expr.prof.dists, "Label")
    genes.specs.df <- data.frame(gene = genes, species = as.character(lapply(genes, 
        speciesForGeneId_Regex, spec.regexs = gene.ids.to.specs.regex)), stringsAsFactors = FALSE)
    species <- intersect(names(gene.ids.to.specs.regex), genes.specs.df$species)
    #' Validate species:
    if (length(species) == 0) 
        stop("Found ZERO species matching gene accessions in argument 'expr.prof.dists' and being present in argument 'gene.ids.to.specs.regex'.")
    #' Extract the requested intra-species distances:
    intra.spec.dists <- setNames(lapply(species, function(spec) {
        spec.genes <- genes.specs.df[which(genes.specs.df$species == spec), "gene"]
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
            x.spec.genes <- genes.specs.df[which(genes.specs.df$species == spec.pair[[1]]), 
                "gene"]
            y.spec.genes <- genes.specs.df[which(genes.specs.df$species == spec.pair[[2]]), 
                "gene"]
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
    fams.inter.intra.exp.prof.dists.df <- data.frame(family = family.names, max.ath = as.numeric(NA), 
        max.chi = as.numeric(NA), max.inter.species = as.numeric(NA), max.abs = as.numeric(NA), 
        median.ath = as.numeric(NA), median.chi = as.numeric(NA), median.inter.species = as.numeric(NA), 
        median.abs = as.numeric(NA), mean.ath = as.numeric(NA), mean.chi = as.numeric(NA), 
        mean.inter.species = as.numeric(NA), mean.abs = as.numeric(NA), stringsAsFactors = FALSE)
    for (fam.name in family.names) {
        fam.dists <- fams.inter.intra.exp.prof.dists[[fam.name]]
        if (!is.null(fam.dists) && !all(is.na(fam.dists)) && length(fam.dists) > 
            0) {
            for (stat in c("max", "median", "mean")) {
                stat.vals <- names(fam.dists)[grepl(stat, names(fam.dists), fixed = TRUE)]
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
expressionProfileDistStatsPerTissues <- function(expr.prof.dists.per.tissue, gene.classes, 
    fam.name = NA, stats = c("max", "median", "mean", "maxMinusMin")) {
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
    genes.classified <- lapply(gene.classes, function(x) intersect(genes, x))
    res <- lapply(genes.classified, function(x) as.dist(epd.m[rownames(epd.m) %in% 
        x, colnames(epd.m) %in% x]))
    res$stats <- Reduce(rbind, lapply(stats, function(i.stat) {
        Reduce(rbind, lapply(names(gene.classes), function(k.class) {
            genes.dist <- res[[k.class]]
            dist.stat <- if (length(genes.dist) > 0 && !all(is.na(genes.dist) | is.nan(genes.dist) | 
                is.infinite(genes.dist))) {
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
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
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
plotTissueSpecificExpressionProfileDistanceDistribs <- function(plot.df, tissues, 
    gene.classes, pdf.path, plot.title = "") {
    plot.df$interaction <- interaction(plot.df$Tissue, plot.df$Gene.Class)
    plot.df$interaction <- factor(plot.df$interaction, levels = unlist(lapply(tissues, 
        paste, gene.classes, sep = ".")))
    line.df <- rbind(data.frame(Tissue = tissues, Value = as.numeric(lapply(tissues, 
        function(x) upperWhisker(plot.df[which(plot.df$Tissue == x & plot.df$Gene.Class == 
            gene.classes[[1]]), "Value"]))), Gene.Class = gene.classes[[1]], stringsAsFactors = FALSE), 
        data.frame(Tissue = tissues, Value = as.numeric(lapply(tissues, function(x) upperWhisker(plot.df[which(plot.df$Tissue == 
            x & plot.df$Gene.Class == gene.classes[[2]]), "Value"]))), Gene.Class = gene.classes[[2]], 
            stringsAsFactors = FALSE))
    br.cols <- brewer.pal(length(levels(plot.df$interaction)), "Paired")
    pdf(pdf.path)
    par(mar = c(5, 15, 4, 2) + 0.1)
    ticks.at <- sort(c(seq(1, length(tissues) * 2 + 3, by = 3), seq(2, length(tissues) * 
        2 + 4, by = 3)))
    boxplot(Value ~ interaction, data = plot.df, outline = FALSE, horizontal = TRUE, 
        las = 1, border = br.cols, col = addAlpha(br.cols), lwd = 2.5, xlab = "median of Euclidean distances between expression profiles", 
        main = plot.title, names = NA, at = ticks.at)
    axis(2, at = ticks.at, labels = unlist(lapply(tissues, function(x) paste(x, " (", 
        gene.classes, ")", sep = ""))), las = 1)
    line.cols <- brewer.pal(3, "Greys")[2:3]
    lines(x = line.df[which(line.df$Gene.Class == gene.classes[[1]]), ]$Value, y = seq(1, 
        length(tissues) * 2 + 3, by = 3), col = line.cols[[1]], lwd = 3)
    lines(x = line.df[which(line.df$Gene.Class == gene.classes[[2]]), ]$Value, y = seq(2, 
        length(tissues) * 2 + 4, by = 3), col = line.cols[[2]], lwd = 3)
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
#' @param x an instance of \code{base::data.frame} with columns: "Family",
#' "Statistic", "Gene.Class", "Value", "Tissue". For example this package's
#' data \code{tandems.exp.prof.dists.tissue.orth.dist}.
#' @param tissues a character vector identifying the tissues to be
#' investigated. Default is \code{c("seedling", "cotyledon", "developing leaf",
#' "flower stage 9", "flower stage 16")}.
#' @param gene.classes a character vector identifying the gene classes to be
#' investigated. Default is \code{c("Orthologs", "Non-Orthologs")}.
#' @param early.late.indx A numeric vector identifying which entries of
#' \code{tissues} are considered early and late developmental stages. Default
#' is \code{c(2,4)}.
#' @param intermediate.indx A single integer identifying which entry in
#' \code{tissues} is considered to be the intermediate developmental stage.
#' Default is \code{3}.
#' @param statistic a string identifying which statistic of \code{x} to be
#' investigated. Default is \code{"median"}.
#'
#' @export
#' @return A \code{base::data.frame} with columns: "Tissue",
#' "Alternative.Hypothesis", "p.value", "effect.size", "p.value.adjusted"
testHourglassModel <- function(x, tissues = c("seedling", "cotyledon", "developing leaf", 
    "flower stage 9", "flower stage 16"), gene.classes = c("Orthologs", "Non-Orthologs"), 
    early.late.indx = c(2, 4), intermediate.indx = 3, statistic = "median") {
    res.df <- rbind(Reduce(rbind, lapply(tissues, function(tissue) {
        y.ort <- x[which(x$Tissue == tissue & x$Statistic == statistic & x$Gene.Class == 
            "Orthologs" & !is.na(x$Value)), "Value"]
        y.no <- x[which(x$Tissue == tissue & x$Statistic == statistic & x$Gene.Class == 
            "Non-Orthologs" & !is.na(x$Value)), "Value"]
        h.test <- t.test(y.ort, y.no, alternative = "less")
        data.frame(Tissue = tissue, Alternative.Hypothesis = "Non-Orthologs have higher mean expression profile distances than Orthologs", 
            p.value = h.test$p.value, effect.size = h.test$statistic, stringsAsFactors = FALSE)
    })), Reduce(rbind, unlist(lapply(gene.classes, function(gene.class) {
        lapply(tissues[early.late.indx], function(tissue) {
            y.dev.leaf <- x[which(x$Tissue == tissues[intermediate.indx] & x$Statistic == 
                statistic & x$Gene.Class == gene.class & !is.na(x$Value)), "Value"]
            y.foreground <- x[which(x$Tissue == tissue & x$Statistic == statistic & 
                x$Gene.Class == gene.class & !is.na(x$Value)), "Value"]
            h.test <- t.test(y.dev.leaf, y.foreground, alternative = "less")
            hypo <- paste(gene.class, "have higher median expression profile distances in", 
                tissue, "than in developing leaf")
            data.frame(Tissue = tissue, Alternative.Hypothesis = hypo, p.value = h.test$p.value, 
                effect.size = h.test$statistic, stringsAsFactors = FALSE)
        })
    }), recursive = FALSE)))
    res.df$p.value.adjusted <- p.adjust(res.df$p.value, method = "BY")
    res.df
}
