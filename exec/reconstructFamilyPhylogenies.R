require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/reconstructFamilyPhylogenies.R working_dir n_cores [Batch-Start Batch-Stop]")
message("NOTE: If you want to reconstruct the Gene families.lst' Phylogenies in parallel, you can provide both Batch-Start and Batch-Stop as indices of the non-singleton families.lst for which to run the Phylogenetic Pipeline.")

#' Working folder:
input.args <- commandArgs(trailingOnly = TRUE)
o.dir <- input.args[[1]]
if (!file.exists(o.dir)) dir.create(o.dir, recursive = TRUE)
#' Cores to use:
n.threads <- if (length(input.args) > 1) as.integer(input.args[[2]]) else detectCores()
options(mc.cores = n.threads)
#' Batch size:
families.non.singleton.inds <- families.df[ which( families.df$size > 1 ), "id" ]
i.start <- if (length(input.args) > 2) as.integer(input.args[[3]]) else 1
i.stop <- if (length(input.args) > 3) as.integer(input.args[[4]]) else length(families.non.singleton.inds)


#'#########################################
#' Process each non singleton gene family #
#'#########################################
for (fam.name in names(families.lst[families.non.singleton.inds[i.start:i.stop]])) {
    message("Generating phylogeny of family ", fam.name)
    tryCatch({
        fam <- as.character(unlist(families.lst[[fam.name]]))
        # Family's directory:
        fam.dir <- file.path(o.dir, fam.name)
        if (!file.exists(fam.dir)) 
            dir.create(fam.dir, recursive = TRUE)
        # Store coding sequences:
        fam.cds.orig <- getFamilyDNAStringSet(fam)
        # Sanitize the gene identifiers:
        fam.name.maps <- sanitizeNames(fam.cds.orig)
        fam.cds <- fam.cds.orig
        names(fam.cds) <- fam.name.maps$sanitized
        fam.cds.path <- file.path(fam.dir, paste(fam.name, "_CDS.fasta", sep = ""))
        write.fasta(sequences=fam.cds, names=names(fam.cds), file.out=fam.cds.path)
        fam.name.maps.path <- file.path(fam.dir, paste(fam.name, "_name_mappings_table.txt", 
            sep = ""))
        write.table(fam.name.maps, fam.name.maps.path, row.names = FALSE, sep = "\t", 
            quote = FALSE)
        # Convert to AA and align the AA-sequences:
        translate2AASeqs(fam.cds.path)
        # Remove invalid AA-Sequences, i.e. AA-Seqs with premature stop-codons:
        fam.aa.path <- file.path(fam.dir, paste(fam.name, "_CDS_macse_AA.fasta", 
            sep = ""))
        fam.aas <- read.fasta(fam.aa.path, seqtype="AA", as.string=TRUE, strip.desc=TRUE)
        fam.aas.san <- fam.aas[validateAAStringSet(fam.aas)]
        # Warn about removed AA-Seqs:
        fam.cds.san <- if (length(fam.aas.san) < length(fam.aas)) {
            len.diff <- length(fam.aas) - length(fam.aas.san)
            warning(len.diff, " amino-acid-sequences had a premature stop codon and were removed from further analysis.")
            fam.cds[names(fam.aas.san)]
        } else fam.cds
        # If only a single sequence is left, continue with next family:
        if (length(fam.cds.san) < 2) {
            warning("Family ", fam.name, " has less than two sequnces left after all quality checks. Stopping phylognentic reconstruction.")
            next
        }
        # Write out the sanitized amino acid seqs:
        fam.aas.san.path <- file.path(fam.dir, paste(fam.name, "_AA_sanitized.fasta", 
            sep = ""))
        write.fasta(sequences=fam.aas.san, names=names(fam.aas.san), file.out=fam.aas.san.path)
        # Generate a multiple sequence alignment:
        fam.aas.msa.path <- file.path(fam.dir, paste(fam.name, "_AA_sanitized_MSA.fasta", 
            sep = ""))
        alignAAStringSet(fam.aas.san.path, fam.aas.msa.path)
        fam.aas.san.msa <- read.fasta(fam.aas.msa.path, seqtype="AA", as.string=TRUE, strip.desc=TRUE)
        # Use the aligned AA-Seqs as quide to align the CDS Sequences:
        fam.cds.msa <- alignCDSSetWithAlignedAAsAsGuide(fam.cds.san, fam.aas.san.msa)
        fam.cds.msa.path <- file.path(fam.dir, paste(fam.name, "_CDS_MSA.fasta", 
            sep = ""))
        write.fasta(fam.cds.msa, names=names(fam.cds.msa), file.out=fam.cds.msa.path)
        # Generate Phylogenetic maximum likelihood Tree:
        fam.tree.path <- file.path(fam.dir, paste(fam.name, "_ml_tree.newick", sep = ""))
        fam.tree.orig.path <- file.path(fam.dir, paste(fam.name, "_ml_tree_original_identifiers.newick", 
            sep = ""))
        buildPhylogeneticTree(fam.cds.msa.path, fam.tree.path)
        fam.tree <- read.tree(fam.tree.path)
        fam.tree.orig <- fam.tree
        fam.tree.orig$tip.label <- replaceWithOriginal(fam.tree$tip.label, fam.name.maps)
        write.tree(fam.tree.orig, fam.tree.orig.path)
        fam.tree.no.node.labels <- fam.tree
        fam.tree.no.node.labels$node.label <- NULL
        fam.tree.no.node.labels.path <- file.path(fam.dir, paste(fam.name, "_ml_tree_no_node_labels.newick", 
            sep = ""))
        write.tree(fam.tree.no.node.labels, fam.tree.no.node.labels.path)
    }, error = function(e) {
        warning("\n")
        warning("Encountered an error while reconstructing phylogeny of family ", 
            fam.name)
        warning(e)
        warning("Continuing with next family...")
        warning("\n")
    })
}

message("Successfully completed") 
