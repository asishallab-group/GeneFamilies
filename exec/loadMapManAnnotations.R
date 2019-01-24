require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamiliesexec/exec/loadMapManAnnotations.R mercatorMapManResults.txt path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


mercator.df <- fread(input.args[[1]], header = TRUE, stringsAsFactors = FALSE, 
    colClasses = c(rep("character", 4), "logical"), na.strings = "", sep = "\t", 
    data.table = FALSE)
#' Strip quotes from columns:
for (c.nm in colnames(mercator.df)[1:4]) {
    mercator.df[, c.nm] <- sub("^'", "", sub("'$", "", mercator.df[, c.nm]))
}


#' List of Bins (names) and their genes (values):
mapManBins <- sort(unique(mercator.df$BINCODE))
mercator.lst <- setNames(lapply(mapManBins, function(mmBin) {
    if (is.na(mercator.df[which(mercator.df$BINCODE == mmBin), "TYPE"][[1]])) {
        mercator.df[which(grepl(paste("^", mmBin, sep = ""), mercator.df$BINCODE) & 
            mercator.df$TYPE), "IDENTIFIER"]
    } else {
        mercator.df[which(mercator.df$BINCODE == mmBin), "IDENTIFIER"]
    }
}), mapManBins)


#' Number of genes annotated to each MapMan-Bin:
mercator.annos.stat <- data.frame(BINCODE = mapManBins, no.annos = sapply(mapManBins, 
    function(mmBin) length(mercator.lst[[mmBin]])), stringsAsFactors = FALSE)


#' Save results:
save(mercator.df, mercator.lst, mercator.annos.stat, file = file.path(input.args[[2]], 
    "data", "mapManBinAnnotations.RData"))


message("DONE")
