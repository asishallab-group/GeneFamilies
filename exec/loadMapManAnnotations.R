require(GeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamiliesexec/exec/loadMapManAnnotations.R mercatorMapManResults.txt path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


mercator.df <- fread(input.args[[1]], header = TRUE, stringsAsFactors = FALSE, 
    colClasses = rep("character", 5), na.strings = "", sep = "\t", data.table = FALSE)
#' Strip quotes from columns:
for ( c.nm in colnames(mercator.df) ) {
  mercator.df[,c.nm] <- sub("^'", "", sub( "'$", "", mercator.df[,c.nm] ) )
}


save(mercator.df, file = file.path(input.args[[2]], "data", "mapManBinAnnotations.RData"))


message("DONE")
