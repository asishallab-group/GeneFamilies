require(MaizeGeneFamilies)


message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/maizePlots.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

#' MapMan Profiles:
maize.mapMan.profile <- annotationProfile()
mms.v.sgltc.mapMan.profile <- annotationProfile(ad.i = which(maize.mapMan$IDENTIFIER.san %in% 
    toLowerCutTail(mms.v.sgltc.de.genes)))
mmbs.v.sgstc.mapMan.profile <- annotationProfile(ad.i = which(maize.mapMan$IDENTIFIER.san %in% 
    toLowerCutTail(mmbs.v.sgstc.de.genes)))


#' Plot the above profiles:
extractCounts <- function(prfl.df, a.grps) {
    total.cnt <- sum(prfl.df$Count)
    sapply(a.grps, function(mm.grp) {
        if (mm.grp %in% prfl.df$Anno.Group) {
            prfl.df[which(prfl.df$Anno.Group == mm.grp), "Count"]/total.cnt
        } else 0
    })
}

p.df <- maize.mapMan.profile[order(maize.mapMan.profile$Count), ]
anno.groups <- p.df$Anno.Group
p.mtrx <- matrix(extractCounts(maize.mapMan.profile, anno.groups), ncol = 1)
p.mtrx <- cbind(p.mtrx, extractCounts(mms.v.sgltc.mapMan.profile, anno.groups))
p.mtrx <- cbind(p.mtrx, extractCounts(mmbs.v.sgstc.mapMan.profile, anno.groups))
rownames(p.mtrx) <- p.df$Name
overrep.mapManBins.groups <- sort(unique(unlist(regmatches(overrep.mapManBins.genes.df$BINCODE, 
    regexec("^\\d+", overrep.mapManBins.genes.df$BINCODE)))))
p.names <- unlist(lapply(p.df$Anno.Group, function(x) {
    y <- p.df[which(p.df$Anno.Group == x), "Name"]
    if (x %in% overrep.mapManBins.groups) 
        paste(y, "*") else paste(y, " ")
}))
pdf(file.path(input.args[[1]], "inst", "mapManBinProfiles.pdf"))
old.par <- par(mar = c(5, 17, 4, 2) + 0.1)
p.col <- brewer.pal(3, "Dark2")
barplot(t(p.mtrx), beside = TRUE, horiz = TRUE, las = 1, names.arg = p.names, 
    col = p.col, border = p.col, xlab="annotation frequency")
legend("right", legend = c("Maize Genome", "Mesophyll", "Bundle Sheath"), 
    fill = p.col, border = p.col, bty = "n")
par(old.par)
dev.off()
