require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/mapManRootBinAnnotations.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' No blacklisting nor filtering using MapMan annotations:
options(AHRD.prot.desc.blacklist = list())
options(AHRD.token.blacklist = list())
options(AHRD.prot.desc.filter = list())
options(AHRD.token.filter = list())


#' Extract MapMan root-bins for the anntotated genes:
root.bins <- mercator.df[grepl("^\\d+$", mercator.df$BINCODE), c("BINCODE", 
    "NAME")]
rownames(root.bins) <- root.bins$BINCODE
i <- which(mercator.df$IDENTIFIER != "")
mapMan.root.bins <- data.frame(SHORT.ID = mercator.df$IDENTIFIER[i], filtered.description = as.character(apply(mercator.df[i, 
    ], 1, function(m.row) {
    root.bins[sub("\\..+$", "", m.row[[1]]), "NAME"]
})), stringsAsFactors = FALSE)


#' For gene families:
families.hrd.mapMan <- setNames(mclapply(families.lst, function(cl) {
    AhrdOnGeneClusters::annotateClusterWithProteinDescription(unlist(cl), 
        prot.desc.db = mapMan.root.bins)
}), names(families.lst))
families.mapMan.root.bins.df <- Reduce(rbind, mclapply(names(families.hrd.mapMan), 
    function(cl.nm) {
        fam.hrd.mm <- families.hrd.mapMan[[cl.nm]]
        if (!is.null(fam.hrd.mm) && !is.na(fam.hrd.mm)) {
            data.frame(id = cl.nm, description.mapMan = fam.hrd.mm$descriptions, 
                entropy.mm.desc = fam.hrd.mm$entropy, stringsAsFactors = FALSE)
        } else NULL
    }))
rm(families.hrd.mapMan)  # Clean up!
gc()



#' For tandem clusters:
tandems.hrd.mapMan <- setNames(mclapply(tandems.lst, function(cl) {
    AhrdOnGeneClusters::annotateClusterWithProteinDescription(unlist(cl), 
        prot.desc.db = mapMan.root.bins)
}), names(tandems.lst))
tandems.mapMan.root.bins.df <- Reduce(rbind, mclapply(names(tandems.hrd.mapMan), 
    function(cl.nm) {
        fam.hrd.mm <- tandems.hrd.mapMan[[cl.nm]]
        if (!is.null(fam.hrd.mm) && !is.na(fam.hrd.mm)) {
            data.frame(id = cl.nm, description.mapMan = fam.hrd.mm$descriptions, 
                entropy.mm.desc = fam.hrd.mm$entropy, stringsAsFactors = FALSE)
        } else NULL
    }))
rm(tandems.hrd.mapMan)  # Clean up!
gc()


#' For ortholog clusters:
orthologs.hrd.mapMan <- setNames(mclapply(orthologs.lst, function(cl) {
    AhrdOnGeneClusters::annotateClusterWithProteinDescription(unlist(cl), 
        prot.desc.db = mapMan.root.bins)
}), names(orthologs.lst))
orthologs.mapMan.root.bins.df <- Reduce(rbind, mclapply(names(orthologs.hrd.mapMan), 
    function(cl.nm) {
        fam.hrd.mm <- orthologs.hrd.mapMan[[cl.nm]]
        if (!is.null(fam.hrd.mm) && !is.na(fam.hrd.mm)) {
            data.frame(id = cl.nm, description.mapMan = fam.hrd.mm$descriptions, 
                entropy.mm.desc = fam.hrd.mm$entropy, stringsAsFactors = FALSE)
        } else NULL
    }))
rm(orthologs.hrd.mapMan)  # Clean up!
gc()


#'Save results:
save(mapMan.root.bins, families.mapMan.root.bins.df, tandems.mapMan.root.bins.df, 
    orthologs.mapMan.root.bins.df, file = file.path(input.args[[1]], "data", 
        "mapManRootBinAnnotations.RData"))

message("DONE")
