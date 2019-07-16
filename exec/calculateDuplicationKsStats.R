require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/calculateFamiliesKsStats.R families-working-dir data-dir")

input.args <- commandArgs(trailingOnly = TRUE)

#' Use absolut paths:
input.args <- lapply(input.args, normalizePath)

#' Find paths to family working dirs:
fams.dirs <- system(paste("find", input.args[[1]], "-maxdepth 1 -type d -name 'cluster_*'"), 
    intern = TRUE)
names(fams.dirs) <- sub("^.*/", "", fams.dirs)

#' Calculate Ks statistics between trans duplicated and ancestral orthologs
dupl.w.orths.dirs <- fams.dirs[names(dupl.w.orths)]

dupl.w.orths.Ks.df <- do.call(rbind, mclapply(names(dupl.w.orths.dirs), 
    function(fam.name) {
        tryCatch({
            fam.dir <- dupl.w.orths.dirs[[fam.name]]
            fam.cds <- read.alignment(file.path(fam.dir, paste0(fam.name, 
                "_CDS_MSA.fasta")), format = "fasta")
            fam.nm <- read.table(file.path(fam.dir, paste0(fam.name, 
                "_name_mappings_table.txt")), sep = "\t", header = TRUE, 
                comment.char = "", na.strings = "", stringsAsFactors = FALSE, 
                quote = "")
            fam.orths <- fam.nm[which(sub("\\.\\d+$", "", fam.nm$original, 
                perl = TRUE) %in% orthologs.genes), "sanitized"]
            fam.ka.ks.tbl <- KaKsStatistics(fam.cds, fam.name, 
                fam.orths, fam.dir, header = TRUE)
            if (!is.null(fam.ka.ks.tbl) && nrow(fam.ka.ks.tbl) > 
                0) {
                ks.tbl <- fam.ka.ks.tbl[, c(1, 4)]
                ks.tbl$Gene.Group <- fam.name
                ks.tbl
            } else NULL
        }, error = function(e) {
            message("An ERROR occurred when calculating the KaKsStatistics for family '", 
                fam.name, "' in directory '", fam.dir, ":\n", 
                e)
            return(NULL)
        })
    }))


#' Calculate Ks statistics for Tandem Clusters with Orthologs:
tands.w.orths.Ks.df <- do.call(rbind, mclapply(names(tands.w.orths), 
    function(tan.cl.name) {
        tryCatch({
            tan.cl.dir <- file.path(input.args[[1]], tan.cl.name)
            dir.create(tan.cl.dir)
            alignCodingSequencesPipeline(getSpliceVariantSeqs(all.cds, 
                tands.w.orths[[tan.cl.name]]), tan.cl.dir, tan.cl.name)
            tan.cl.cds <- read.alignment(file.path(tan.cl.dir, 
                paste0(tan.cl.name, "_CDS_MSA.fasta")), format = "fasta")
            tan.cl.nm <- read.table(file.path(tan.cl.dir, paste0(tan.cl.name, 
                "_name_mappings_table.txt")), sep = "\t", header = TRUE, 
                comment.char = "", na.strings = "", stringsAsFactors = FALSE, 
                quote = "")
            tan.cl.orths <- tan.cl.nm[which(tan.cl.nm$original %in% 
                orthologs.genes), "sanitized"]
            tan.cl.ka.ks.tbl <- KaKsStatistics(tan.cl.cds, tan.cl.name, 
                tan.cl.orths, tan.cl.dir, header = TRUE)
            if (!is.null(tan.cl.ka.ks.tbl) && nrow(tan.cl.ka.ks.tbl) > 
                0) {
                ks.tbl <- tan.cl.ka.ks.tbl[, c(1, 4)]
                ks.tbl$Gene.Group <- tan.cl.name
                ks.tbl
            } else NULL
        }, error = function(e) {
            message("An ERROR occurred when calculating the KaKsStatistics for tandem cluster '", 
                tan.cl.name, "':\n", e)
            return(NULL)
        })
    }))


# Save results:
save(dupl.w.orths.Ks.df, tands.w.orths.Ks.df, file = file.path(input.args[[2]], 
    "DuplicatedKsStats.RData"))

message("DONE")
