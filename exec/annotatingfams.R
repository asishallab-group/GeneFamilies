require(GeneFamilies)
require( parallel )
options('mc.cores'=detectCores())
require(AHRD.on.gene.clusters)

message("USAGE: Rscript path/2/GeneFamilies/exec/annotatingfams.R path path/2/GeneFamilies/data)

# annotate families
mcl.fams.hrds <- mclapply(families.lst, FUN=annotateCluster, ipr.annos=all.ipr,      interpro.database=ipr.db )

# get only name and bind to df
families.df.desc <- families.df
families.df.desc$desc <- unlist(lapply(mcl.fams.hrds, function(x) { if( length(x) > 0) { y <- x[["most.frequent.IPRs"]][[1]]$NAME; if (! is.null(y) && length(y) > 0) { y } else { NA } } else NA }))

# save result
save(families.df.desc, file=file.path(input.args[[1]], "AHDR_annotated_fams.RData"))

message("DENME UN <3")
