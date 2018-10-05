require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/loadBustedResults.R path/2/Families/Working-Dir path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' This script needs data on which genes have been found to have a Ka/Ks ratio
#' significantly greater than 1: 
all.KaKs <- unique(unlist(pairwise.Ka.Ks[which(pairwise.Ka.Ks$w > 1), c("gene.a", 
    "gene.b")]))


#' Load BUSTED result files:
busted.out.fls <- system(paste("find", input.args[[1]], "-type f -name 'cluster_*_CDS_MSA.nex.BUSTED.json'"), 
    intern = TRUE)
names(busted.out.fls) <- regmatches(busted.out.fls, regexec("cluster_\\d+", busted.out.fls))
busted.p.vals <- setNames(as.numeric(unlist(mclapply(busted.out.fls, function(x) as.numeric(system(paste("tail -1", 
    x, "| sed -e 's/[^0-9.e-]//g'"), intern = TRUE))))), names(busted.out.fls))


#' Infer which genes have strong evidence for positive selection:
busted.df <- data.frame( family=as.character(NA), gene=all.KaKs, p.value=NA, stringsAsFactors = FALSE )
for ( fam in names( busted.p.vals )) {
  fam.genes <- families.lst[[fam]]
  fam.foreground.genes <- intersect( fam.genes, all.KaKs )
  i <- which( busted.df$gene %in% fam.foreground.genes )
  busted.df[ i, "p.value" ] <- busted.p.vals[[fam]]
  busted.df[ i, "family" ] <- fam
}


#' Save results:
save( busted.df, file=file.path( input.args[[2]], "data", "BUSTED_Results.RData" ) )


message("DONE")
