#' Brew template to generate a family's BranchREL input file
#' @export
hyphy.branch.site.bf <- "inputRedirect = {};\ninputRedirect[\"01\"]=\"Universal\";\ninputRedirect[\"02\"]=\"Yes\";\ninputRedirect[\"03\"]=\"Yes\";\ninputRedirect[\"04\"]=\"<%= fam.cds.msa.path %>\";\ninputRedirect[\"05\"]=\"<%= fam.tree.4.paml.path %>\";\ninputRedirect[\"06\"]=\"<%= fam.hyphy.branch.site.output.path %>\";\n\nExecuteAFile (\"<%= hyphy.batch.files.dir %>/BranchSiteREL.bf\", inputRedirect);"

#' Brew template to generate a family's MEME input file
#' @export
hyphy.meme.bf <- "inputRedirect = {};\ninputRedirect[\"01\"]=\"Universal\";\ninputRedirect[\"02\"]=\"New Analysis\";\ninputRedirect[\"03\"]=\"<%= fam.cds.msa.path %>\";\ninputRedirect[\"04\"]=\"Custom\";\ninputRedirect[\"05\"]=\"110240\";\ninputRedirect[\"06\"]=\"<%= fam.tree.no.node.labels.path %>\";\ninputRedirect[\"07\"]=\"<%= fam.hyphy.meme.log.path %>\";\ninputRedirect[\"08\"]=\"Estimate dN/dS only\";\ninputRedirect[\"09\"]=\"MEME\";\ninputRedirect[\"10\"]=\"0.1\";\ninputRedirect[\"11\"]=\"N\";\ninputRedirect[\"12\"]=\"<%= fam.hyphy.meme.output.path %>\";\n\nExecuteAFile (\"<%= hyphy.batch.files.dir %>/QuickSelectionDetection.bf\", inputRedirect);"

#' Brew template to generate a family's FUBAR input file
#' @export
hyphy.fubar.bf <- "inputRedirect = {};\ninputRedirect[\"01\"]=\"Universal\";\ninputRedirect[\"02\"]=\"1\";\ninputRedirect[\"03\"]=\"<%= fam.cds.msa.path %>\";\ninputRedirect[\"04\"]=\"<%= fam.tree.no.node.lables.path %>\";\ninputRedirect[\"05\"]=\"20\";\ninputRedirect[\"06\"]=\"5\";\ninputRedirect[\"07\"]=\"2000000\";\ninputRedirect[\"08\"]=\"1000000\";\ninputRedirect[\"09\"]=\"100\";\ninputRedirect[\"10\"]=\"0.5\";\ninputRedirect[\"11\"]=\"<%= fam.hyphy.fubar.output.path %>\";\n\nExecuteAFile (\"<%= hyphy.batch.files.dir %>/FUBAR.bf\", inputRedirect);"

#' Brew template to generate a family's BUSTED input file
#' @export
hyphy.busted.bf <- "inputRedirect = {};\ninputRedirect[\"01\"]=\"Universal\";\ninputRedirect[\"02\"]=\"<%= fam.cds.msa.nexus.path %>\";\ninputRedirect[\"03\"]=\"<%= fam.hyphy.busted.tree.path %>\";\ninputRedirect[\"04\"]=\"Set TEST\";\ninputRedirect[\"05\"]=\"\";\n\nExecuteAFile (\"<%= hyphy.batch.files.dir %>/BUSTED.bf\", inputRedirect);"

#' Brew template to generate a protein group's Mr Bayes input file
#' @export
mr.bayes.prot.tree <- "\nbegin mrbayes;\nlog start filename=<%= mr.bayes.log.file %> replace;\nprset aamodelpr=fixed(wag);\nlset rates=invgamma Ngammacat=4;\nset autoclose=yes;\nmcmc ngen=20000 printfreq=500 samplefreq=10\nnchains=<%= mr.bayes.n.chains %> savebrlens=yes starttree=random\nfilename=<%= mr.bayes.out.file %>;\nsumt filename=<%= mr.bayes.out.file %> burnin=1000 showtreeprobs=yes\ncontype=allcompat Conformat=Simple;\nquit;\nend;\n"

#' Brew template to generate a protein group's Mr Bayes input file using
#' additional e.g. morphological knowledge
#' @export
mr.bayes.mixed.prot.tree <- "\nbegin mrbayes;\n<% for (x in names(mr.bayes.char.sets)) { %>\n<% y <- mr.bayes.char.sets[[x]] %>\ncharset <%= x %> = <%= y[[1]] %> - <%= y[[2]] %>;\n<% } %>\npartition favored = <%= length(mr.bayes.char.sets) %>: <%= paste( names( mr.bayes.char.sets ), collapse=\", \" ) %>;\nset partition = favored;\n\nlog start filename=<%= mr.bayes.log.file %> replace;\nprset aamodelpr=fixed(wag);\nlset rates=invgamma Ngammacat=4;\nset autoclose=yes;\nmcmc ngen=20000 printfreq=500 samplefreq=10\nnchains=<%= mr.bayes.n.chains %> savebrlens=yes starttree=random\nfilename=<%= mr.bayes.out.file %>;\nsumt filename=<%= mr.bayes.out.file %> burnin=1000 showtreeprobs=yes\ncontype=allcompat Conformat=Simple;\nquit;\nend;\n"

#' Brew template to generate a coding sequences group's Mr Bayes input file
#' using additional e.g. morphological knowledge
#' @export
mr.bayes.mixed.cds.tree <- "\nbegin mrbayes;\n<% for (x in names(mr.bayes.char.sets)) { %>\n<% y <- mr.bayes.char.sets[[x]] %>\ncharset <%= x %> = <%= y[[1]] %> - <%= y[[2]] %>;\n<% } %>\npartition favored = <%= length(mr.bayes.char.sets) %>: <%= paste( names( mr.bayes.char.sets ), collapse=\", \" ) %>;\nset partition = favored;\n\nlog start filename=<%= mr.bayes.log.file %> replace;\nlset nst=6 rates=invgamma Nucmodel=Codon;\nset autoclose=yes;\nmcmc ngen=20000 printfreq=500 samplefreq=10\nnchains=<%= mr.bayes.n.chains %> savebrlens=yes starttree=random\nfilename=<%= mr.bayes.out.file %>;\nsumt filename=<%= mr.bayes.out.file %> burnin=1000 showtreeprobs=yes\ncontype=allcompat Conformat=Simple;\nquit;\nend;\n"

#' Brew template to generate an Mr Bayes input file for a multiple alignment of
#' Coding Sequences (CDS MSA).
#' @export
mr.bayes.cds.tree <- "\nbegin mrbayes;\nlog start filename=<%= mr.bayes.log.file %> replace;\nlset nst=6 rates=invgamma Nucmodel=Codon;\nset autoclose=yes;\nmcmc ngen=20000 printfreq=500 samplefreq=10\nnchains=<%= mr.bayes.n.chains %> savebrlens=yes starttree=random\nfilename=<%= mr.bayes.out.file %>;\nsumt filename=<%= mr.bayes.out.file %> burnin=1000 showtreeprobs=yes\ncontype=allcompat Conformat=Simple;\nquit;\nend;\n"

#' Removes node labels and branch lengths from instance of ape::phylo
#'
#' @param phyl.tree instance of class ape::phylo
#'
#' @return A copy of 'phyl.tree' without branch length and without node labels
#' @export
removeNodeLabelsAndBranchLengths <- function(phyl.tree) {
    phyl.tree$node.label <- NULL
    phyl.tree$edge.length <- NULL
    phyl.tree
}

#' Annotates a phylogenetic tree in such amanner that foreground branches are
#' recognized by HYPHY's BUSTED analysis.
#'
#' @param phylo.tree an instance of class ape::phylo
#' @param foreground.nodes integer vector holding the nodes that are selected
#' as foreground
#'
#' @export
#' @return A BUSTED annotated copy of 'phylo.tree'
tree4Busted <- function(phylo.tree, foreground.nodes) {
    if (!is.null(foreground.nodes) && length(foreground.nodes) > 0) {
        #' Mark foreground tips / leaves:
        fg.tips <- intersect(1:length(phylo.tree$tip.label), foreground.nodes)
        if (!is.null(fg.tips) && length(fg.tips) > 0) 
            phylo.tree$tip.label[fg.tips] <- paste(phylo.tree$tip.label[fg.tips], 
                "{TEST}", sep = "")
        #' Mark foreground inner nodes:
        fg.in.nds <- which(1:phylo.tree$Nnode %in% (foreground.nodes - length(phylo.tree$tip.label)))
        phylo.tree$node.label <- character(length(phylo.tree$node.label))
        if (!is.null(fg.in.nds) && length(fg.in.nds) > 0) 
            phylo.tree$node.label[fg.in.nds] <- "{TEST}"
    }
    phylo.tree
}

#' Extracts pairs of genes from a CDS multiple sequence alignment. The ones
#' extracted are identified by computing all pairwise sequence based distances
#' from the multiple sequence alignment. Then those pairs representing the
#' argument statistics, e.g. 'min', 'mean', 'median', and 'max' are used to
#' generate the axt file for subsequent analysis with KaKs_Calculator
#' (https://sourceforge.net/projects/kakscalculator2/).
#'
#' @param cds.msa The result of invoking \code{seqinr::read.fasta} representing
#' the multiple alignment of coding sequences.
#' @param dist.stats The statistics that are to be matched with the respective
#' pairwise sequence based distances. Default is \code{c('min', 'median',
#' 'mean', 'max')}.
#'
#' @return A named list. Names are those of argument \code{dist.stats}. These
#' names can be concatonated with the letter 'n' if a pair of genes fulfills
#' the criterion of various statistics. Values are the gene pairs themselfes.
#' @export
axtForSelectedDistances <- function(cds.msa, dist.stats=c('min', 'median', 'mean', 'max')) {
    anchor.cds <- toString(cds.msa[[anchor.gene]])
    other.genes <- setdiff(names(cds.msa), anchor.gene)
    setNames(lapply(other.genes, function(x) {
        list(anchor.cds, toString(cds.msa[[x]]))
    }), other.genes)
}

#' Extracts pairs of genes from a CDS multiple sequence alignment.
#'
#' @param cds.msa an instance of Biostrings::DNAStringSet representing the
#' multiple alignment of coding sequences.
#' @param anchor.gene the gene identifier for which to build pairs with the
#' other genes present in 'cds.msa'.
#'
#' @return A named list. Names are the non anchor gene identifier that are
#' paired with the anchor gene. Values are the gene pairs themselfes.
#' @export
axt <- function(cds.msa, anchor.gene) {
    anchor.cds <- toString(cds.msa[[anchor.gene]])
    other.genes <- setdiff(names(cds.msa), anchor.gene)
    setNames(lapply(other.genes, function(x) {
        list(anchor.cds, toString(cds.msa[[x]]))
    }), other.genes)
}

#' Writes pairwise aligned sequences to a file of AXT format for Ka/Ks
#' calculation.
#'
#' @param axt a named list as returned from function 'axt(...)'
#' @param axt.path the valid file path indicating the output AXT file to
#' generate
#'
#' @return the result of calling readLines 
#' @export
writeAxt <- function(axt, axt.path) {
    axt.lines <- paste(as.character(unlist(lapply(names(axt), function(nm) {
        gene.pair <- axt[[nm]]
        paste(nm, gene.pair[[1]], gene.pair[[2]], sep = "\n")
    }))), collapse = "\n\n")
    writeLines(axt.lines, con = axt.path)
}


#' Identifies duplicated inverse pairs and discards them.
#'
#' @param clst.hmlgs.list a named list in which for each pair member A is the
#' name and member B is the value
#' @param ret.as.list indicates whether to return the result as a list, subset
#' of the original 'clst.hmlgs.list', or as a character matrix of two columns one
#' column per member. For the latter option set 'ret.as.list' to FALSE.
#'
#' @return The subset of pairs that are not inverse duplicates, either as list
#' or as character matrix.
#' @export
removeDuplicatedInversePairs <- function(clst.hmlgs.list, ret.as.list = TRUE) {
    srtd.pairs <- do.call("rbind", lapply(names(clst.hmlgs.list), function(x) {
        sort(c(x, clst.hmlgs.list[[x]]))
    }))
    non.dupl.pairs <- srtd.pairs[!duplicated(srtd.pairs), , drop = FALSE]
    if (ret.as.list) 
        setNames(as.list(non.dupl.pairs[, 2]), non.dupl.pairs[, 1]) else non.dupl.pairs
}

#' Identifies those tip and inner nodes of the phylogenetic tree that span
#' subtrees whose tips are true subsets of the argument 'foreground.nodes'.
#'
#' @param phylo.tr instance of ape::phylo representing a phylogenetic tree
#' @param foreground.nodes a vector of tip indices or tip labels that represent
#' candidate foreground tips.
#'
#' @export
#' @return integer vector of inner node indices that span subtrees containing
#' only foreground tips merged with the argument 'foreground.nodes'
foregroundNodes <- function(phylo.tr, foreground.nodes) {
    #' If tip labels are passed as argument infer their index in the tree's nodes
    #' slot:
    if (inherits(foreground.nodes, "character")) 
        foreground.nodes <- which(phylo.tr$tip.label %in% foreground.nodes)
    #' Find those inner nodes that spawn subtrees whose tips are true subsets of
    #' 'foreground.nodes'
    in.nds <- setdiff(phylo.tr$edge[, 1], 1:length(phylo.tr$tip.label))
    in.nds.spawned <- setNames(lapply(in.nds, function(x) Descendants(phylo.tr, x, 
        type = "tips")[[1]]), in.nds)
    union(foreground.nodes, as.integer(names(in.nds.spawned[as.logical(lapply(in.nds.spawned, 
        function(x) {
            length(setdiff(x, foreground.nodes)) == 0
        }))])))
}

#' Uses base::read.table to parse the output of 'KaKs_Calculator'.
#'
#' @param path.2.KaKs.Calc.output the valid file path to the output table to be
#' parsed
#'
#' @return An instance of base::data.frame with three columns: 'Sequence',
#' 'Ka.Ks', 'P.Value.Fisher.'
#' @export
readKaKsCalculatorOutput <- function(path.2.KaKs.Calc.output) {
    read.table(path.2.KaKs.Calc.output, header = TRUE, sep = "\t", comment.char = "", 
        quote = "", colClasses = c("character", rep("NULL", 3), rep("numeric", 2), 
            rep("NULL", 16)))
}

#' Extracts the gene pairs' identifiers from the 'Sequence' column of
#' 'hom.kaks.tbl' and replaces this column with the parse result.
#'
#' @param hom.kaks.tbl the result of calling 'readKaKsCalculatorOutput(...)' on
#' the result file of KaKs_Calculator used on homologs, that were NOT conserved
#' orthologs.
#'
#' @return A modified version of 'hom.kaks.tbl' in which the column 'Sequence'
#' is replaced with two columns 'gene.1' and 'gene.2'.
#' @export
addGeneIdsCols2HomologsKaKsTable <- function(hom.kaks.tbl) {
    x <- do.call("rbind", strsplit(hom.kaks.tbl$Sequence, "~"))
    colnames(x) <- c("gene.1", "gene.2")
    cbind(x, hom.kaks.tbl[, 2:3], stringsAsFactors = FALSE)
}

#' Uses 'gene.1' as the first member of the orthologous gene pair and extracts
#' the gene identifier of the pairs second member from the 'Sequence' column of
#' 'orth.kaks.tbl' and replaces this column with the result.
#'
#' @param orth.kaks.tbl the result of calling 'readKaKsCalculatorOutput(...)'
#' on the result file of KaKs_Calculator used on conserved orthologs.
#' @param gene.1 the identifier of the anchoring gene used within the currently
#' processed group of conserved orthologs.
#'
#' @return A modified version of 'orth.kaks.tbl' in which the column 'Sequence'
#' is replaced with two columns 'gene.1' and 'gene.2'.
#' @export
addGeneIdsCols2OrthologsKaKsTable <- function(orth.kaks.tbl, gene.1) {
    x <- data.frame(gene.1 = gene.1, gene.2 = orth.kaks.tbl$Sequence, stringsAsFactors = FALSE)
    cbind(x, orth.kaks.tbl[, 2:3], stringsAsFactors = FALSE)
}

#' Parses a FUBAR output table for significant posterior probabilities of a
#' codon being subject to positive selection. If any codons fit the criteria, a
#' respective table is returned, otherwise \code{NULL} is returned.
#'
#' @param path.2.fubar.csv a valid file path to the FUBAR csv file to parse
#' @param family.name A string indicating the name of the gene cluster (family)
#' subjected to the FUBAR analysis.
#' @param sign.post.prob The cutoff value the posterior probabilities have to
#' meet (\code{>=}) in order to be considered significant. Default is
#' \code{.9}.
#'
#' @return An instance of base::data.frame holding the tabular FUBAR output
#' with just the following columns: 1. Codon, 2. Post.Prob, 3. BayesFactor, and
#' 4. family. \code{NULL} is returned if no codon matches the significance
#' criteria.
#' @export
readFubarTable <- function(path.2.fubar.csv, family.name, sign.post.prob = 0.9) {
    fubar.txt <- system(paste("awk -F \",\" '{if ( $5 >=", sign.post.prob, ") print $1 \"\\t\" $5 \"\\t\" $7 }'", 
        path.2.fubar.csv), intern = TRUE)
    if (length(fubar.txt) > 1) {
        x <- read.table(text = fubar.txt, skip = 1, stringsAsFactors = FALSE, sep = "\t", 
            comment.char = "", quote = "", colClasses = rep("numeric", 2))
        colnames(x) <- c("Codon", "Post.Prob", "BayesFactor")
        x$family <- family.name
        x
    } else NULL
}

#' Translates the FUBAR reported expected absolute NUMBER of false positives
#' into a relative false positive rate for each family within which positively
#' selected codons have been identified.
#'
#' @param fubar.tbl.res a named vector in which names are families and values
#' are data.frames of FUBAR results
#' @param fubar.expct.fls.pos.n named vector in which names are families and
#' values are the expected absolute number of false positives as reported by
#' FUBAR.
#'
#' @return A named vector in which names are families and values are relative
#' false positive rates per family. 
#' @export
inferExpectedFalsePositiveRate <- function(fubar.tbl.res, fubar.expct.fls.pos.n) {
    setNames(as.numeric(lapply(names(fubar.expct.fls.pos.n), function(fam) {
        as.numeric(fubar.expct.fls.pos.n[[fam]])/nrow(fubar.tbl.res[[fam]])
    })), names(fubar.expct.fls.pos.n))
}

#' Generates a valid input file for Mr Bayes, which will generate a Bayesian
#' phylogeny from a multiple amino acid sequence alignment.
#'
#' @param aa.msa an instance of Biostrings::AAStringSet representing the MSA
#' @param mr.bayes.dir a valid file path to the output directory in which to
#' store the Mr Bayes analysis related files
#' @param gene.group.name the name of the group of proteins to generate the
#' Bayesian tree for
#' @param mr.bayes.n.chains number of chains to use in Mr. Bayes, default is
#' 'mc.cores' or 4, if 'mc.cores' is not set.
#'
#' @return NULL
#' @export
mrBayesProteinTree <- function(aa.msa, mr.bayes.dir, gene.group.name, mr.bayes.n.chains = getOption("mc.cores", 
    4)) {
    mr.bayes.log.file <- paste(gene.group.name, "_mr_bayes_out.log", sep = "")
    mr.bayes.out.file <- paste(gene.group.name, "_mr_bayes_out", sep = "")
    gene.group.msa.nexus.path <- file.path(mr.bayes.dir, paste(gene.group.name, "_AA_MSA.nex", 
        sep = ""))
    write.nexus.data(as.list(aa.msa), gene.group.msa.nexus.path, format = "PROTEIN")
    mr.bayes.inp.file <- file(gene.group.msa.nexus.path, "a")
    brew(text = mr.bayes.prot.tree, output = mr.bayes.inp.file)
    close(mr.bayes.inp.file)
}

#' Generates a valid input file for Mr Bayes, which will generate a Bayesian
#' phylogeny from a multiple coding sequence (codons) alignment.
#'
#' @param cds.msa an instance of Biostrings::DNAStringSet representing the MSA
#' of coding sequences (CODONS)
#' @param mr.bayes.dir a valid file path to the output directory in which to
#' store the Mr Bayes analysis related files
#' @param gene.group.name the name of the group of proteins to generate the
#' Bayesian tree for
#' @param mr.bayes.n.chains number of chains to use in Mr. Bayes, default is
#' 'mc.cores' or 4, if 'mc.cores' is not set.
#'
#' @return NULL
#' @export
mrBayesCdsTree <- function(cds.msa, mr.bayes.dir, gene.group.name, mr.bayes.n.chains = getOption("mc.cores", 
    4)) {
    mr.bayes.log.file <- paste(gene.group.name, "_mr_bayes_out.log", sep = "")
    mr.bayes.out.file <- paste(gene.group.name, "_mr_bayes_out", sep = "")
    gene.group.msa.nexus.path <- file.path(mr.bayes.dir, paste(gene.group.name, "_CDS_MSA.nex", 
        sep = ""))
    write.nexus.data(as.list(cds.msa), gene.group.msa.nexus.path, format = "dna")
    mr.bayes.inp.file <- file(gene.group.msa.nexus.path, "a")
    brew(text = mr.bayes.cds.tree, output = mr.bayes.inp.file)
    close(mr.bayes.inp.file)
}


#' Generates a valid input file for Mr Bayes, which will generate a Bayesian
#' phylogeny from a mixed multiple sequence alignment, e.g. CDS AND an alignment of
#' additional e.g. morphological data.
#'
#' @param msa.lst A list of multiple sequence alignments, a mix of DNA and
#' Standard data
#' @param msa.types a character vector of the types of each MSA: DNA, Codon, or
#' Standard
#' @param mr.bayes.dir The file directory in which to store the output file
#' @param the name of the analysis. Will be the head of the file name holding
#' the input for Mr Bayes
#'
#' @return NULL
#' @export
mrBayesTreeFromMixedMSA <- function(msa.lst, msa.types, mr.bayes.dir, gene.group.name, 
    mr.bayes.n.chains = getOption("mc.cores", 4)) {
    mr.bayes.log.file <- paste(gene.group.name, "_mr_bayes_out.log", sep = "")
    mr.bayes.out.file <- paste(gene.group.name, "_mr_bayes_out", sep = "")
    gene.group.msa.nexus.path <- file.path(mr.bayes.dir, paste(gene.group.name, ".nex", 
        sep = ""))
    gene.ids <- names(msa.lst[[1]])
    mixed.msa <- setNames(lapply(gene.ids, function(x) {
        paste(lapply(msa.lst, function(i.msa) toString(i.msa[[x]])), collapse = "")
    }), gene.ids)
    write.nexus.data(mixed.msa, gene.group.msa.nexus.path, format = "dna")
    nchar.mix <- sum(unlist(lapply(msa.lst, function(x) nchar(x)[[1]])))
    mixed.nex <- readChar(gene.group.msa.nexus.path, file.info(gene.group.msa.nexus.path)$size)
    mixed.nex <- sub("NCHAR=1", paste("NCHAR=", nchar.mix, sep = ""), mixed.nex, 
        fixed = TRUE)
    n <- 1
    mixed.nex <- sub("DATATYPE=DNA", paste("DATATYPE=MIXED(", paste(lapply(1:length(msa.types), 
        function(i) {
            msa.len <- nchar(msa.lst[[i]])[[1]]
            msa.str <- paste(msa.types[[i]], ":", n, "-", msa.len + (n - 1), sep = "")
            n <<- n + msa.len
            msa.str
        }), collapse = ","), ")", sep = ""), mixed.nex, fixed = TRUE)
    writeLines(toString(mixed.nex), gene.group.msa.nexus.path)
    n <- 1
    mr.bayes.char.sets <- setNames(lapply(msa.lst, function(x) {
        msa.lst <- list(n, nchar(x)[[1]] + (n - 1))
        n <<- n + nchar(x)[[1]]
        msa.lst
    }), names(msa.lst))
    gene.group.msa.nexus.file <- file(gene.group.msa.nexus.path, "a")
    brew(text = mr.bayes.mixed.cds.tree, output = gene.group.msa.nexus.file)
    close(gene.group.msa.nexus.file)
}

#' Generates a valid input file for Mr Bayes, which will generate a Bayesian
#' phylogeny from a multiple amino acid sequence alignment AND an alignment of
#' additional e.g. morphological data.
#'
#' @param aa.msa an instance of Biostrings::AAStringSet representing the MSA
#' @param extra.msa a mamed list of additional e.g. morphological data,
#' aligned. Names have to be identical with names of 'aa.msa'.
#' @param mr.bayes.dir a valid file path to the output directory in which to
#' store the Mr Bayes analysis related files
#' @param gene.group.name the name of the group of proteins to generate the
#' Bayesian tree for
#' @param mr.bayes.n.chains number of chains to use in Mr. Bayes, default is
#' 'mc.cores' or 4, if 'mc.cores' is not set.
#'
#' @return NULL
#' @export
mrBayesProteinTreeFromMixedMSA <- function(aa.msa, extra.msa, mr.bayes.dir, gene.group.name, 
    mr.bayes.n.chains = getOption("mc.cores", 4)) {
    mr.bayes.log.file <- paste(gene.group.name, "_mr_bayes_out.log", sep = "")
    mr.bayes.out.file <- paste(gene.group.name, "_mr_bayes_out", sep = "")
    gene.group.msa.nexus.path <- file.path(mr.bayes.dir, paste(gene.group.name, "_mixed_AA_extra_MSA.nex", 
        sep = ""))
    mixed.msa <- setNames(lapply(names(aa.msa), function(x) {
        paste(toString(aa.msa[[x]]), extra.msa[[x]], sep = "")
    }), names(aa.msa))
    write.nexus.data(mixed.msa, gene.group.msa.nexus.path, format = "PROTEIN")
    nchar.aa <- unique(nchar(aa.msa))[[1]]
    nchar.extra <- unique(nchar(extra.msa))[[1]]
    mixed.nex <- readChar(gene.group.msa.nexus.path, file.info(gene.group.msa.nexus.path)$size)
    mixed.nex <- sub("NCHAR=1", paste("NCHAR=", (nchar.aa + nchar.extra), sep = ""), 
        mixed.nex, fixed = TRUE)
    mixed.nex <- sub("DATATYPE=PROTEIN", paste("DATATYPE=MIXED(PROTEIN:1-", nchar.aa, 
        ",Standard:", (nchar.aa + 1), "-", (nchar.aa + nchar.extra), ")", sep = ""), 
        mixed.nex, fixed = TRUE)
    writeLines(toString(mixed.nex), gene.group.msa.nexus.path)
    mr.bayes.char.sets <- list(aaseqs = list(1, nchar.aa), extra = list((1 + nchar.aa), 
        (nchar.aa + nchar.extra)))
    gene.group.msa.nexus.file <- file(gene.group.msa.nexus.path, "a")
    brew(text = mr.bayes.mixed.prot.tree, output = gene.group.msa.nexus.file)
    close(gene.group.msa.nexus.file)
}

#' Uses GNU sed to replace the sanitized gene identifiers with the original
#' ones in the Fig-Tree file - consensus tree generated by Mr. Bayes.
#'
#' @param name.mappings data.frame with two columns 'original', and
#' 'sanitized', holding the pairs of original gene ID to sanitized gene ID.
#' @param path.2.fig.tree the valid file path to the tree file to process
#' @param path.2.fig.tree.orig the valid file path to the output file in which
#' the original gene identifiers are used. Must not be identical to
#' 'path.2.fig.tree'. Default is 'path.2.fig.tree' with its suffix '.con.tre'
#' modified to '_original_identifiers.con.tre'.
#'
#' @return NULL
#' @export
replaceWithOriginalInFigTree <- function(name.mappings, path.2.fig.tree, path.2.fig.tree.orig = sub(".con.tre", 
    "_original_identifiers.con.tre", path.2.fig.tree, fixed = TRUE)) {
    system(paste(c("sed", paste("-e 's/\\b", name.mappings$sanitized, "\\b/", name.mappings$original, 
        "/'", sep = ""), path.2.fig.tree, ">", path.2.fig.tree.orig), collapse = " "))
}

#' Parses a HyPhy MEME branch result file to extract significant LEAVES with
#' codons showing signs of positive selection.
#'
#' @param meme.file a valid file path to a MEME branch result file, e.g.
#' 'fam8119_hyphy_meme_output.txt.branches'
#' @param bayes.factor.sign.cutoff The numeric value above which to consider a
#' Bayes-Factor to be significant. Default is highly conservative and set to
#' 100 aka 'decisive'
#' @param fam.name.reg.exp The regular expression used within
#' \code{base::regexec(...)} to extract the gene family's name from the
#' \code{meme.file} argument
#'
#' @return A data.frame with four columns: Family, Site, Branch, and
#' BayesianFactor, where the rows hold the significant results. NULL is
#' returned if no significant results were contained in \code{meme.file}.
#' @export
parseMemeBranchResults <- function(meme.file, bayes.factor.sign.cutoff = 100, fam.name.reg.exp = "/(fam\\d+)/") {
    fam.name <- regmatches(meme.file, regexec(fam.name.reg.exp, meme.file))[[1]][[2]]
    fam.dir <- sub(paste(fam.name, "_hyphy_meme_output.txt.branches", sep = ""), 
        "", meme.file, fixed = TRUE)
    meme.signif.branches <- system(paste("sed -e '1,2d'", meme.file, "| awk -F \",\" '{ if($1 ~ /Site/ || $4 >", 
        bayes.factor.sign.cutoff, "&& $2 ~ /PROT/) { print $1 \"\\t\" $2 \"\\t\" $4 } }'"), 
        intern = TRUE)
    if (!is.null(meme.signif.branches) && !is.na(meme.signif.branches) && length(meme.signif.branches) > 
        1) {
        fam.name.maps <- read.table(file.path(fam.dir, paste(fam.name, "_name_mappings_table.txt", 
            sep = "")), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", 
            comment.char = "", na.strings = "", colClasses = rep("character", 2))
        x <- cbind(Family = fam.name, read.table(text = meme.signif.branches, sep = "\t", 
            quote = "", comment.char = "", na.strings = "", colClasses = c("integer", 
                "character", "numeric"), header = TRUE, stringsAsFactors = FALSE))
        x$Branch <- replaceWithOriginal(x$Branch, fam.name.maps)
        x
    } else NULL
}
