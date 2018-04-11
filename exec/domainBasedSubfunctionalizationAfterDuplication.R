require(MaizeGeneFamilies)
options(mc.cores = getMcCores())


message("USAGE:\n[MCCORES=<Int>] Rscript path/2/MaizeGeneFamilies/exec/domainBasedSubfunctionalizationAfterDuplication.R duplicated|tandem start_index end_index output_table.txt")

input.args <- commandArgs(trailingOnly = TRUE)

#' Parse input arguments:
if (input.args[[1]] == "tandem") {
    gene.set <- tands.w.orths
    g.classifier <- tand.classifier
} else if (input.args[[1]] == "duplicated") {
    gene.set <- dupl.w.orths
    g.classifier <- dupl.classifier
} else {
    stop("Please specify gene set as one of 'tandem' or 'duplicated'.")
}

s.i <- as.integer(input.args[[2]])
start.i <- if (s.i >= 0 && s.i < length(gene.set)) s.i else stop("start_index not in range of '", 
    input.args[[1]], "'.")
rm(s.i)

e.i <- as.integer(input.args[[3]])
end.i <- if (e.i > 0 && e.i <= length(gene.set)) e.i else stop("end_index not in range of of '", 
    input.args[[1]], "'.")
rm(e.i)

out.tbl <- if (!is.null(input.args[[4]])) file.path(input.args[[4]]) else stop("Please specify a valid output table path. ", 
    input.args[[4]], " is not a valid path.")
#' END: parse input arguments


#' Investigate possible subfunctionilization events between base ortholog genes
#' and evolved duplicated/tandem genes:
gene.set.i <- gene.set[start.i:end.i]
gene.set.subfunc <- geneGroupSubfunctionalization(gene.set.i, g.classifier)

#' Save results:
if( !is.null(gene.set.subfunc) && ! is.na(gene.set.subfunc) ) {
  write.table( gene.set.subfunc, out.tbl, row.names=FALSE, sep="\t", quote=TRUE )
} else {
  message("NOTIFICATION: No subfunctionaliziation could be identified within the argument gene set.")
  message("NOT writing an output table.")
}


message("DONE")
