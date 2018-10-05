require(GeneFamilies)
options(mc.cores = getMcCores())


message("USAGE: Rscript path/2/GeneFamilies/exec/annotateFamilies.R path/2/MaizeGeneFamilies/data")

input.args <- commandArgs(trailingOnly = TRUE)


#' Compute the human readable descritions for the Gene-Families:
families.hrd <- mclapply(families.lst, annotateCluster, ipr.annos = all.ipr, 
    interpro.database = ipr.db)

#' Create a new data.frame extended version of 'families.df' with an extra
#' column holding the HRDs:
families.HRD.df <- families.df
families.HRD.df$description <- as.character(unlist(mclapply(families.df[, "id"], 
    function(fam.name) {
        if (fam.name %in% names(families.hrd)) {
            fam.hrd <- families.hrd[[fam.name]]
            if (!is.null(fam.hrd) && !is.na(fam.hrd) && length(fam.hrd) > 0) {
                fam.freq.iprs <- fam.hrd$most.frequent.IPRs
                if (!is.null(fam.freq.iprs) && !is.na(fam.freq.iprs) && length(fam.freq.iprs) > 
                  0) {
                  paste(unlist(lapply(fam.freq.iprs, function(x) x$SHORT.NAME)), 
                    collapse = ",")
                } else NA
            } else NA
        } else NA
    })))


#' Save results:
save(families.hrd, families.HRD.df, file = file.path(input.args[[1]], "familyHumanReadableDescriptions.RData"))

message("DONE")
