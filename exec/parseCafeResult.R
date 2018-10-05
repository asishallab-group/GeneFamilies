require(GeneFamilies)
options(mc.cores=detectCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/parseCafeResult.R path/2/cafe_result.txt.cafe path/2/MaizeGeneFamilies/data")

input.args <- commandArgs(trailingOnly = TRUE)

spec.names <- c("B.rapa","S.parvula","E.salsugineum","C.hirsuta","C.rubella","A.thaliana","A.lyrata","A.arabicum")

#' Infer how to obtain the species specific P-Values from CAFE results:
node.maps <- sub("^.+:", "", readLines(input.args[[1]], n = 3)[[3]], perl = TRUE)
spec.nodes <- setNames(lapply(spec.names, function(x) getSpeciesNodeNumberFromCafe(node.maps, 
    x)), spec.names)
node.pos <- strsplit(gsub("[()\\s]", "", gsub("\\) \\(", ",", sub("^.+:", "", readLines(input.args[[1]], 
    n = 4)[[4]], perl = TRUE), perl = TRUE), perl = TRUE), ",")[[1]]
spec.pos <- unlist(lapply(spec.nodes, function(x) which(node.pos == x)))

#' Extract CAFE species specific P-Values, correct for multiple hypothesis
#' testing, and add them to the families.df data.frame:
cf.tbl <- read.table(input.args[[1]], stringsAsFactors = FALSE, skip = 9, header = TRUE)
cafe.result.df <- do.call("rbind", mclapply(cf.tbl$Viterbi.P.values, function(x) {
    parseCafeBranchSpecificPValues(x, spec.pos = spec.pos)
}))
rownames(cafe.result.df) <- cf.tbl$ID
colnames(cafe.result.df) <- paste("CAFE", colnames(cafe.result.df), sep = ".")
#' Adjust for multiple hypothesis testing, do this not wining about NA values:
suppressWarnings({
    for (i in 1:ncol(cafe.result.df)) {
        cafe.result.df[, i] <- p.adjust(cafe.result.df[, i], method = "BY")
    }
})

#' Save results:
save(cafe.result.df, file = file.path(input.args[[2]], "cafe_result.RData")) 
