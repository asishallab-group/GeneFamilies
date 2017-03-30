require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadCodingSequences.R path/2/GeneFamilies/data ath.fa cla.fa cme.fa csa.fa fve.fa")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load coding sequences in DNA alphabetically:
ath.cds <- loadDnaFasta(input.args[[2]])
cla.cds <- loadDnaFasta(input.args[[3]])
cme.cds <- loadDnaFasta(input.args[[4]])
csa.cds <- loadDnaFasta(input.args[[5]])
fve.cds <- loadDnaFasta(input.args[[6]])

#' Compile a single list holding ALL CDS of all eight Brassicaceaen genomes:
all.cds <- Reduce( append, list( ath.cds,cla.cds,cme.cds,csa.cds,fve.cds))

#' Save results:
save(ath.cds,cla.cds,cme.cds,fve.cds,csa.cds, all.cds,
    file = file.path(input.args[[1]], "codingSequences.RData"))
message("DONE")
