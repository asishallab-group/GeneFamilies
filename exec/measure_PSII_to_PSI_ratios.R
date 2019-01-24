require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/measure_PSII_to_PSI_ratios.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Identify MapMan-Bins indicating proteins of the Photosystem I:
ps.I.bins <- unique(mercator.df[grepl("^1.1.4.2[0-9.]*$", mercator.df$BINCODE), 
    "BINCODE"])
ps.I.chi <- unique(mercator.df[which(mercator.df$BINCODE %in% ps.I.bins & 
    grepl("^carhr", mercator.df$IDENTIFIER.no.expr.var)), "IDENTIFIER.no.expr.var"])
ps.I.ath <- unique(mercator.df[which(mercator.df$BINCODE %in% ps.I.bins & 
    grepl("^at", mercator.df$IDENTIFIER.no.expr.var)), "IDENTIFIER.no.expr.var"])


#' Identify MapMan-Bins indicating proteins of the Photosystem II:
ps.II.bins <- unique(mercator.df[grepl("^1.1.1.2[0-9.]*$", mercator.df$BINCODE), 
    "BINCODE"])
ps.II.chi <- unique(mercator.df[which(mercator.df$BINCODE %in% ps.II.bins & 
    grepl("^carhr", mercator.df$IDENTIFIER.no.expr.var)), "IDENTIFIER.no.expr.var"])
ps.II.ath <- unique(mercator.df[which(mercator.df$BINCODE %in% ps.II.bins & 
    grepl("^at", mercator.df$IDENTIFIER.no.expr.var)), "IDENTIFIER.no.expr.var"])


#' Measure ratio of PS-II over PS-I;
#' - before artificial shade treatment:
ps.I.chi.median.baseMeanA <- median(chi.light[which(chi.light$id.lower.case %in% 
    ps.I.chi), "baseMeanA"], na.rm = TRUE)
ps.I.ath.median.baseMeanA <- median(ath.light[which(ath.light$id.lower.case %in% 
    ps.I.ath), "baseMeanA"], na.rm = TRUE)
ps.II.chi.median.baseMeanA <- median(chi.light[which(chi.light$id.lower.case %in% 
    ps.II.chi), "baseMeanA"], na.rm = TRUE)
ps.II.ath.median.baseMeanA <- median(ath.light[which(ath.light$id.lower.case %in% 
    ps.II.ath), "baseMeanA"], na.rm = TRUE)
#' - after artificial shade treatment:
ps.I.chi.median.baseMeanB <- median(chi.light[which(chi.light$id.lower.case %in% 
    ps.I.chi), "baseMeanB"], na.rm = TRUE)
ps.I.ath.median.baseMeanB <- median(ath.light[which(ath.light$id.lower.case %in% 
    ps.I.ath), "baseMeanB"], na.rm = TRUE)
ps.II.chi.median.baseMeanB <- median(chi.light[which(chi.light$id.lower.case %in% 
    ps.II.chi), "baseMeanB"], na.rm = TRUE)
ps.II.ath.median.baseMeanB <- median(ath.light[which(ath.light$id.lower.case %in% 
    ps.II.ath), "baseMeanB"], na.rm = TRUE)

#' Store results in a table:
ps.expr.levels <- data.frame(species = c(rep(c("ath", "ath", "chi", "chi"), 
    2)), photosystem = c(rep("PS-I", 4), rep("PS-II", 4)), condition = c(rep(c("untreated", 
    "artificial shade"), 4)), median.expr.level = c(ps.I.ath.median.baseMeanA, 
    ps.I.ath.median.baseMeanB, ps.I.chi.median.baseMeanA, ps.I.chi.median.baseMeanB, 
    ps.II.ath.median.baseMeanA, ps.II.ath.median.baseMeanB, ps.II.chi.median.baseMeanA, 
    ps.II.chi.median.baseMeanB), stringsAsFactors = FALSE)

#' Compute ratios:
ps.II.over.ps.I.ratios <- data.frame(ps.II.over.ps.I = c((ps.II.chi.median.baseMeanA/ps.I.chi.median.baseMeanA), 
    (ps.II.chi.median.baseMeanB/ps.I.chi.median.baseMeanB), (ps.II.ath.median.baseMeanA/ps.I.ath.median.baseMeanA), 
    (ps.II.ath.median.baseMeanB/ps.I.ath.median.baseMeanB)), condition = c("untreated", 
    "artificial shade", "untreated", "artificial shade"), species = c(rep("chi", 
    2), rep("ath", 2)), stringsAsFactors = FALSE)


#' Save results:
write.table(ps.II.over.ps.I.ratios, file.path(input.args[[1]], "inst", 
    "PS_II_over_PS_I_ratios.txt"), row.names = FALSE, quote = FALSE, sep = "\t")


message("DONE")
