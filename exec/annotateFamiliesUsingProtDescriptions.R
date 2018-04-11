require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/annotateFamiliesUsingProtDescriptions.R protein_description_database.txt path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

# Load the protein descrition database:
gf.prot.desc.db <- fread(input.args[[1]], sep = "\t", na.strings = "", 
    header = TRUE, stringsAsFactors = FALSE, colClasses = rep("character", 
        4), data.table = FALSE)


# Annotate Gene Families based on Protein Descriptions found for their
# respective members. See AHRD.on.gene.clusters for more details.
families.hrd.pd <- setNames(mclapply(families.lst, function(fam) {
    AhrdOnGeneClusters::annotateClusterWithProteinDescription(unlist(fam), 
        prot.desc.db = gf.prot.desc.db, prot.id.col = "ID")
}), names(families.lst))
families.hrd.pd.df <- Reduce(rbind, mclapply(names(families.hrd.pd), function(fam.nm) {
    fam.hrd.prot <- families.hrd.pd[[fam.nm]]
    if (!is.null(fam.hrd.prot) && !is.na(fam.hrd.prot)) {
        data.frame(id = fam.nm, description.prot = paste(fam.hrd.prot$descriptions, 
            collapse = ","), entropy.prot.desc = fam.hrd.prot$entropy, 
            stringsAsFactors = FALSE)
    } else NULL
}))
families.HRD.df <- merge(families.HRD.df, families.hrd.pd.df, by = "id")



# Annotate Gene Families based on MapMan-Bin root-bin annotations found
# for their respective members. See AHRD.on.gene.clusters for more
# details.
if (exists("mercator.df")) {
    #' No blacklisting nor filtering using MapMan annotations:
    options(AHRD.prot.desc.blacklist = c())
    options(AHRD.token.blacklist = c())
    options(AHRD.prot.desc.filter = c())
    root.bins <- mercator.df[grepl("^\\d+$", mercator.df$BINCODE), c("BINCODE", 
        "NAME")]
    rownames(root.bins) <- root.bins$BINCODE
    i <- which(mercator.df$IDENTIFIER != "")
    mapMan.root.bins <- data.frame(SHORT.ID = mercator.df$IDENTIFIER[i], 
        filtered.description = as.character(apply(mercator.df[i, ], 1, 
            function(m.row) {
                root.bins[sub("\\..+$", "", m.row[[1]]), "NAME"]
            })), stringsAsFactors = FALSE)
    #' For families:
    families.hrd.mapMan <- setNames(mclapply(families.lst, function(fam) {
        AhrdOnGeneClusters::annotateClusterWithProteinDescription(unlist(fam), 
            prot.desc.db = mapMan.root.bins)
    }), names(families.lst))
    families.hrd.mapMan.df <- Reduce(rbind, mclapply(names(families.hrd.mapMan), 
        function(fam.nm) {
            fam.hrd.mm <- families.hrd.mapMan[[fam.nm]]
            if (!is.null(fam.hrd.mm) && !is.na(fam.hrd.mm)) {
                data.frame(id = fam.nm, description.mapMan = paste(fam.hrd.mm$descriptions, 
                  collapse = ","), entropy.mm.desc = fam.hrd.mm$entropy, 
                  stringsAsFactors = FALSE)
            } else NULL
        }))
    families.HRD.df <- merge(families.HRD.df, families.hrd.mapMan.df, by = "id")
    rm(families.hrd.mapMan.df, families.hrd.mapMan)  # Clean up!
    gc()
} else warning("No MapMan-Bin annotations found. If you want MapMan based gene family descriptions run ./exec/loadMapManAnnotations.R first.")


# Save results:
save(families.hrd, families.HRD.df, file = file.path(input.args[[2]], "data", 
    "familyHumanReadableDescriptions.RData"))

message("DONE")
