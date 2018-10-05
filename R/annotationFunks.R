#' Generate a profile of gene annotations. Count all annotations that are
#' groupable by a regular expression. This function has been designed for
#' MapMan-Bins.
#'
#' @param anno.df An instance of \code{base::data.frame} with at least one
#' column, holding the annotation codes. Default is \code{maize.mapMan}.
#' @param ad.i logical or integer vector indicating which subset of argument
#' \code{anno.df} to consider. Default is to consider it completely,
#' \code{1:nrow(anno.df)}.
#' @param anno.regex A regular expression to identify groups of annotations.
#' Default is \code{'^\\d+\\.\\d+$'}.
#' @param ad.anno.col The column index or name of argument \code{anno.df} in
#' which to find the annotation codes. Default is \code{'BINCODE'}.
#' @param ad.anno.col The column index or name of argument \code{anno.df} in
#' which to find the annotation name. Default is \code{'NAME'}.
#'
#' @export
#' @return An instance of \code{base::data.frame} with the following columns:
#' 'Anno.Group' holds the group of annotation, 'Count' the number of matching
#' annotations, and 'Name' holds the annotations' names.
annotationProfile <- function(anno.df = maize.mapMan, ad.i = 1:nrow(anno.df), 
    anno.regex = "^\\d+", ad.anno.col = "BINCODE", ad.anno.name = "NAME") {
    anno.groups <- sort(unique(unlist(regmatches(anno.df[ad.i, ad.anno.col], 
        regexec(anno.regex, anno.df[ad.i, ad.anno.col])))))
    Reduce(rbind, lapply(anno.groups, function(mm.b) {
        anno.count <- length(which(grepl(paste("^", mm.b, sep = ""), anno.df[ad.i, 
            ad.anno.col])))
        anno.name <- anno.df[which(anno.df[, ad.anno.col] == mm.b), ad.anno.name][[1]]
        data.frame(Anno.Group = mm.b, Count = anno.count, Name = anno.name, 
            stringsAsFactors = FALSE)
    }))
}
