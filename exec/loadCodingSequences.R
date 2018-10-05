require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadCodingSequences.R path/2/MaizeGeneFamilies/data aet.fa aly.fa ath.fa bra.fa chi.fa cru.fa esa.fa tpa.fa")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load coding sequences in DNA alphabetically:
aet.cds <- loadDnaFasta(input.args[[2]])
aly.cds <- loadDnaFasta(input.args[[3]])
ath.cds <- loadDnaFasta(input.args[[4]])
bra.cds <- loadDnaFasta(input.args[[5]])
chi.cds <- loadDnaFasta(input.args[[6]])
cru.cds <- loadDnaFasta(input.args[[7]])
esa.cds <- loadDnaFasta(input.args[[8]])
tpa.cds <- loadDnaFasta(input.args[[9]])

#' Compile a single list holding ALL CDS of all eight Brassicaceaen genomes:
all.cds <- Reduce(append, list(aet.cds, aly.cds, ath.cds, bra.cds, chi.cds, cru.cds, 
    esa.cds, tpa.cds))

#' Save results:
save(aet.cds, aly.cds, ath.cds, bra.cds, chi.cds, cru.cds, esa.cds, tpa.cds, all.cds, 
    file = file.path(input.args[[1]], "codingSequences.RData"))
message("DONE")
