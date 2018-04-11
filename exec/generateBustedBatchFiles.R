require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/generateBustedBatchFiles.R path/2/families_dir")

#' This script needs data on which genes have been found to have a Ka/Ks ratio
#' significantly greater than 1: 
all.KaKs <- unique(unlist(pairwise.Ka.Ks[which(pairwise.Ka.Ks$w > 1), c("gene.a", 
    "gene.b")]))


#' Read Input:
input.args <- commandArgs(trailingOnly = TRUE)


#' Find families for which the phylogenetic pipeline has produced results:
fams <- system(paste("find", input.args[[1]], "-maxdepth 1 -type d -regex '.*cluster_[0-9]+'"), 
    intern = TRUE)
fam.nms <- regmatches(fams, regexpr("cluster_\\d+", fams))
names(fams) <- fam.nms


#' Generate BUSTED input files:
fams.no.foreground.genes <- c()
for (fam.name in names(fams)) {
    tryCatch({
        fam.dir <- fams[[fam.name]]
        fam.tree <- read.tree(file.path(fam.dir, paste(fam.name, "_ml_tree.newick", 
            sep = "")))
        fam.tree.orig <- read.tree(file.path(fam.dir, paste(fam.name, "_ml_tree_original_identifiers.newick", 
            sep = "")))
        fam.foreground.genes <- intersect(all.KaKs, fam.tree.orig$tip.label)
        if ( length( fam.foreground.genes ) > 0 ) {
        fam.tree.4.busted <- tree4Busted(fam.tree, foregroundNodes(fam.tree.orig, 
            fam.foreground.genes))
        fam.hyphy.busted.tree.path <- file.path(fam.dir, paste(fam.name, "_ml_tree_4_BUSTED.newick", 
            sep = ""))
        write.tree(fam.tree.4.busted, fam.hyphy.busted.tree.path)
        fam.cds.msa.nexus.path <- file.path(fam.dir, paste(fam.name, "_CDS_MSA.nex", 
            sep = ""))
        write.nexus.data(loadDnaFasta(file.path(fam.dir, paste(fam.name, "_CDS_MSA.fasta", 
            sep = "")), sanitize.names = FALSE), file = fam.cds.msa.nexus.path, datablock = FALSE, 
            interleaved = FALSE)
        fam.hyphy.busted.batch.file.path <- file.path(fam.dir, paste(fam.name, "_hyphy_busted_input.bf", 
            sep = ""))
        brew(text = hyphy.busted.bf, output = fam.hyphy.busted.batch.file.path)
        } else {
          fams.no.foreground.genes <- append( fams.no.foreground.genes, fam.name )
        }
    }, error = function(e) {
        message("Could not generate BUSTED HyPhy Batch-File for family '", fam.name, 
            "'. An error has occurred:")
        message(e)  #' Continue...
    })
}
#' Inform about families without foreground genes (pairwise Ka/Ks with closest
#' Homolog was greater than 1.0):
message(length(fams.no.foreground.genes), " families had no foreground genes based on pairwise Ka.Ks measurements with closest homologs:")
message(paste(sort(fams.no.foreground.genes), collapse=", "))


message("DONE")
