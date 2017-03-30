require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadInterProData.R path/2/GeneFamilies/data interpro.xml ath_ipr.txt fvesca_v1.0_ipr.txt watermelon_v1_ipr.txt chi_iprscan_out.tsv cru_iprscan_out.tsv esa_iprscan_out.tsv tpa_iprscan_out.tsv")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load InterPro annotations per genome:
ath.ipr <- unique(read.table(input.args[[3]], sep = "\t", comment.char = "", na.string = "", 
    quote = "", colClasses = rep("character", 2), stringsAsFactors = FALSE))
fve.ipr <- unique(read.table(input.args[[4]], sep = "\t", comment.char = "", na.string = "", 
    quote = "", colClasses = rep("character", 2), stringsAsFactors = FALSE))
cla.ipr <- unique(read.table(input.args[[5]], sep = "\t", comment.char = "", na.string = "", 
    quote = "", colClasses = rep("character", 2), stringsAsFactors = FALSE))
# bra.ipr <- read.table(input.args[[5]], sep='\t', comment.char='',
# na.string='', quote='', colClasses=rep('character',2),
# stringsAsFactors=FALSE) chi.ipr <- read.table(input.args[[6]], sep='\t',
# comment.char='', na.string='', quote='', colClasses=rep('character',2),
# stringsAsFactors=FALSE) cru.ipr <- read.table(input.args[[7]], sep='\t',
# comment.char='', na.string='', quote='', colClasses=rep('character',2),
# stringsAsFactors=FALSE) esa.ipr <- read.table(input.args[[8]], sep='\t',
# comment.char='', na.string='', quote='', colClasses=rep('character',2),
# stringsAsFactors=FALSE) tpa.ipr <- read.table(input.args[[9]], sep='\t',
# comment.char='', na.string='', quote='', colClasses=rep('character',2),
# stringsAsFactors=FALSE)
#' Concatonate them into a single base::data.frame:
all.ipr <- unique(Reduce(rbind, list(ath.ipr, fve.ipr, cla.ipr)))

#' Load the InterPro-Database from the provided XML file:
ipr.db <- parseInterProXML(input.args[[2]])

#' Save loaded data:
save(ath.ipr, fve.ipr, cla.ipr, all.ipr, ipr.db, file = file.path(input.args[[1]], 
    "interPro.RData"))

message("DONE")
