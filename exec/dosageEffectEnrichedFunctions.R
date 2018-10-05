require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/dosageEffectEnrichedFunctions.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)
message("For those families and tandem clusters showing dossage effect (correlation between copy number and RPKM expression) we identify the enriched protein functions.")

# Investigated species:
s.g.i <- spec.gene.ids[c("ath", "chi")]

# Families with dossage effect:
i.f <- which(fams.paralogs.copy.no.rpkm.corr.df$r.squared >= 0.9)
fams.paralogs.dosage.efct <- fams.paralogs.copy.no.rpkm.corr.df[i.f, "Family"]
genes.f <- intersect(unlist(s.g.i), unlist(families.lst[fams.paralogs.dosage.efct]))
fams.paralogs.copy.no.rpkm.corr.ipr.enrich <- cbindAnnotationDescription(selectSignificantEnrichments(enrichmentTests(all.ipr[which(all.ipr$V1 %in% 
    genes.f), ], all.ipr)))


# Tandem Clusters with dossage effect:
i.t <- which(tandems.rpkm.corr.df$r.squared >= 0.9)
tandem.clstrs.dosage.efct <- tandems.rpkm.corr.df[i.t, "Family"]
genes.t <- intersect(unlist(s.g.i), unlist(tandems.lst[tandem.clstrs.dosage.efct]))
tandems.rpkm.corr.ipr.enrich <- cbindAnnotationDescription(selectSignificantEnrichments(enrichmentTests(all.ipr[which(all.ipr$V1 %in% 
    genes.t), ], all.ipr)))


# Save results:
save(fams.paralogs.copy.no.rpkm.corr.ipr.enrich, tandems.rpkm.corr.ipr.enrich, 
    fams.paralogs.dosage.efct, tandem.clstrs.dosage.efct, file = file.path(input.args[[1]], 
        "data", "dosageEffectEnrichedFunctions.RData"))

message("DONE")
