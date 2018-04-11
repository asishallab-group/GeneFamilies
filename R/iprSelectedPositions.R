#' Parses the result file of InterProScan into a data.frame.
#'
#' @param path.2.iprscan.res The file path to the argument InterProScan result
#' table.
#'
#' @return A data.frame with 14 columns.
#' @export
#' 
readInterProScanResultTable <- function(path.2.iprscan.res) {
    read.table(path.2.iprscan.res, sep = "\t", stringsAsFactors = FALSE, 
        fill = TRUE, comment.char = "", quote = "", na.string = "")
}


#' Uses function 'readInterProScanResultTable(path.2.iprscan.res)' and retains
#' only two columns: 1. Gene Identifier, and 2. InterPro Identifier.
#'
#' @param path.2.iprscan.res The file path to the argument InterProScan result
#' table.
#' @param gene.col The column identifier for the column in which to find the
#' gene identifiers.
#' @param ipr.col The column identifier for the column in which to find the
#' InterPro identifiers.
#'
#' @export
#' @return An instance of base::data.frame
readMinimumInterProScanResultTable <- function(path.2.iprscan.res, gene.col = 1, 
    ipr.col = 12) {
    ipr.df <- readInterProScanResultTable(path.2.iprscan.res)
    if (!is.null(ipr.df) && !is.na(ipr.df) && nrow(ipr.df) > 0) {
        ipr.df <- ipr.df[which(!is.na(ipr.df[, gene.col]) & !is.na(ipr.df[, 
            ipr.col]) & ipr.df[, ipr.col] != "NA"), c(gene.col, ipr.col)]
        colnames(ipr.df) <- paste("V", 1:2, sep = "")
    }
    unique(ipr.df)
}


#' Identifies conserved protein domains (InterPro) overlapping with the
#' conserved homologous amino acid (codon) 'sel.aa'.
#'
#' @param sel.aa The index of the conserved amino acid under selection in the
#' corresponding multiple sequence alignment (MSA)
#' @param gene The gene accession / ID as used in 'msa.fasta'
#' @param msa.fasta The result of
#' \code{seqinr::read.fasta(path_2_AAs_MSA.fasta, seqtype='AA', as.string=TRUE,
#' strip.desc=TRUE)}.
#' @param iprscan.tbl The result of calling
#' readInterProScanResultTable(path_2_interproscan_result_table.tsv)
#'
#' @return  A character vector of matching InterPro entries.
#' @export
#' 
iprWithSelectedAA <- function(sel.aa, gene, msa.fasta, iprscan.tbl) {
    sel.aa.unaligned.pos <- unalignedAAforAlignedAAPos(gene, sel.aa, msa.fasta)
    if (!is.na(sel.aa.unaligned.pos)) 
        domainsForPos(gene, sel.aa.unaligned.pos, iprscan.tbl) else NA
}


#' Computes the position in the unaligned sequence of a character indicated by
#' \code{sel.aa}.
#'
#' @param  gene The gene accession / ID as used in 'msa.fasta'
#' @param  sel.aa The index of the homologous amino acid subject to selection
#' (integer coordinate)
#' @param msa.fasta The result of
#' \code{seqinr::read.fasta(path_2_AAs_MSA.fasta, seqtype='AA', as.string=TRUE,
#' strip.desc=TRUE)}.
#' @param gap.char The character or regular expression used to identify non
#' sequence characters in the aligned sequences. Default is
#' \code{getOption('MaizeGeneFamilies.gap.char', '-')}.
#' @param gap.char.matching.fixed boolean indicating the value to pass to
#' \code{grepl} and \code{gsub} argument \code{fixed=}. Set to \code{TRUE} if
#' \code{gap.char} is literal and atomic. Default is
#' getOption('MaizeGeneFamilies.gap.char.matching.fixed', TRUE).
#'
#' @return  An integer; either NA if the position is a non sequence character,
#' or the corresponding un-aligned position.
#' @export
unalignedAAforAlignedAAPos <- function(gene, sel.aa, msa.fasta, gap.char = getOption("MaizeGeneFamilies.gap.char", 
    "-"), gap.char.matching.fixed = getOption("MaizeGeneFamilies.gap.char.matching.fixed", 
    TRUE)) {
    tryCatch({
        aa.algn.seq <- msa.fasta[[gene]][[1]]
        aa.algn.char <- substr(aa.algn.seq, sel.aa, sel.aa)
        if (grepl(gap.char, aa.algn.char, fixed = gap.char.matching.fixed)) {
            NA
        } else {
            nchar(gsub(gap.char, "", substr(aa.algn.seq, 1, sel.aa), fixed = gap.char.matching.fixed))
        }
    }, error = function(e) {
        browser()
    })
}


#' Looks up the conserved protein domains (InterPro) that have been annotated
#' to gene and overlap with amino acid position 'aa.pos'.
#'
#' @param gene The gene accession or ID as used in iprscan.tbl
#' @param aa.pos The amino acid position (non aligned) to be overlapped by any
#' InterPros
#' @param gene.col The column index or name of iprscan.tbl in which to look up
#' the gene accessions
#' @param start.col The column index or name of iprscan.tbl in which to look up
#' the domain start positions
#' @param end.col The column index or name of iprscan.tbl in which to look up
#' the domain end positions
#' @param ipr.col The column index or name of iprscan.tbl in which to look up
#' the InterPro identifier
#'
#' @return  A character vector of matching InterPro domains.
#' @export
#' 
domainsForPos <- function(gene, aa.pos, iprscan.tbl, gene.col = "V1", start.col = "V7", 
    end.col = "V8", ipr.col = "V12") {
    x <- iprscan.tbl[which(iprscan.tbl[, gene.col] == gene), ]
    sort(unique(x[which(x[, start.col] <= aa.pos & x[, end.col] >= aa.pos), 
        ipr.col]), na.last = NA)
}


#' Replaces sanitized protein identifiers in XStringSet 'xstring.set' with the
#' original ones held in 'name.maps' table.
#'
#' @param xstring.set A vector or list of strings
#' @param name.maps A data.frame with at least two columns one holding the
#' sanitized and the other the original protein IDs
#' @param san.col The column index or name of 'name.maps' holding the sanitized
#' protein IDs
#' @param orig.col The column index or name of 'name.maps' holding the original
#' protein IDs
#'
#' @return  A copy of argument 'xstring.set' with the sanitized protein IDs
#' replaced with their originals.
#' @export
#' 
replaceSanitizedWithOriginalIDs <- function(xstring.set, name.maps, san.col = "sanitized", 
    orig.col = "original") {
    names(xstring.set) <- as.character(lapply(names(xstring.set), function(x) {
        name.maps[which(name.maps[, san.col] == x), orig.col]
    }))
    xstring.set
}
