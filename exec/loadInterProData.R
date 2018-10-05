require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/loadInterProData.R path/2/MaizeGeneFamilies/data aet_iprscan_out.tsv aly_iprscan_out.tsv ath_iprscan_out.tsv bra_iprscan_out.tsv chi_iprscan_out.tsv cru_iprscan_out.tsv esa_iprscan_out.tsv tpa_iprscan_out.tsv interpro.xml")

input.args <- commandArgs(trailingOnly = TRUE)

#' Load InterPro annotations per genome:
aet.ipr <- readMinimumInterProScanResultTable(input.args[[2]])
aly.ipr <- readMinimumInterProScanResultTable(input.args[[3]])
ath.ipr <- readMinimumInterProScanResultTable(input.args[[4]])
bra.ipr <- readMinimumInterProScanResultTable(input.args[[5]])
chi.ipr <- readMinimumInterProScanResultTable(input.args[[6]])
cru.ipr <- readMinimumInterProScanResultTable(input.args[[7]])
esa.ipr <- readMinimumInterProScanResultTable(input.args[[8]])
tpa.ipr <- readMinimumInterProScanResultTable(input.args[[9]])
#' Concatonate them into a single base::data.frame:
all.ipr <- unique(Reduce(rbind, list(aet.ipr, aly.ipr, ath.ipr, bra.ipr, chi.ipr, 
    cru.ipr, esa.ipr, tpa.ipr)))

#' Load the InterPro-Database from the provided XML file:
ipr.db <- parseInterProXML(input.args[[10]])

#' Save loaded data:
save(aet.ipr, aly.ipr, ath.ipr, bra.ipr, chi.ipr, cru.ipr, esa.ipr, tpa.ipr, all.ipr, 
    ipr.db, file = file.path(input.args[[1]], "interPro.RData"))

message("DONE")
